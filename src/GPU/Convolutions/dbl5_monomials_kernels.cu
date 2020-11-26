// The file dbl5_monomials_kernels.cu defines the kernels specified
// in dbl5_monomials_kernels.h.

/* The algorithm to compute forward, backward, and cross products
 * (denoted respectively by arrays f, b, and c)
 * for a monomial cff*x[0]*x[1]* .. *x[n-1] goes as follows:
 *
 * f[0] := cff*x[0]
 * for i from 1 to n-1 do f[i] := f[i-1]*x[i]
 * if n > 2 then
 *    b[0] := x[n-1]*x[n-2]
 *    for i from 1 to n-3 do b[i] := b[i-1]*x[n-2-i]
 *    b[n-3] := b[n-3]*cff
 *    if n = 3 then
 *       c[0] = f[0]*x[2]
 *    else
 *       for i from 0 to n-4 do c[i] := f[i]*b[n-4-i]
 *       c[n-3] := f[n-3]*x[n-1]
 *
 * Compared to the evaluation and differentiation of a product of variables,
 * (without coefficient cff), two extra multiplications must be done,
 * but this is better than n+1 multiplications with cff afterwards. */

#include "dbl5_convolutions_kernels.cu"
#include "dbl5_monomials_kernels.h"

void GPU_dbl5_speel
 ( int BS, int nvr, int deg, int *idx,
   double *cfftb, double *cffix, double *cffmi, double *cffrg, double *cffpk,
   double *inputtb, double *inputix, double *inputmi,
   double *inputrg, double *inputpk,
   double *forwardtb, double *forwardix, double *forwardmi,
   double *forwardrg, double *forwardpk,
   double *backwardtb, double *backwardix, double *backwardmi,
   double *backwardrg, double *backwardpk,
   double *crosstb, double *crossix, double *crossmi,
   double *crossrg, double *crosspk )
{
   const int deg1 = deg+1;
   int ix1,ix2,ix3;

   ix1 = idx[0]*deg1;                                     // f[0] = cff*x[0]
   dbl5_padded_convolute<<<1,BS>>>
      (cfftb,cffix,cffmi,cffrg,cffpk,
       &inputtb[ix1],&inputix[ix1],&inputmi[ix1],&inputrg[ix1],&inputpk[ix1],
       forwardtb,forwardix,forwardmi,forwardrg,forwardpk,deg1);

   for(int i=1; i<nvr; i++)                            // f[i] = f[i-1]*x[i]
   {
      ix2 = idx[i]*deg1; ix3 = i*deg1; ix1 = ix3 - deg1;
      dbl5_padded_convolute<<<1,BS>>>
         (&forwardtb[ix1],&forwardix[ix1],&forwardmi[ix1],
          &forwardrg[ix1],&forwardpk[ix1],
          &inputtb[ix2],&inputix[ix2],&inputmi[ix2],
          &inputrg[ix2],&inputpk[ix2],
          &forwardtb[ix3],&forwardix[ix3],&forwardmi[ix3],
          &forwardrg[ix3],&forwardpk[ix3],deg1);
   }
   if(nvr > 2)
   {
      ix1 = idx[nvr-1]*deg1; ix2 = idx[nvr-2]*deg1;  // b[0] = x[n-1]*x[n-2]
      dbl5_padded_convolute<<<1,BS>>>
         (&inputtb[ix1],&inputix[ix1],&inputmi[ix1],
          &inputrg[ix1],&inputpk[ix1],
          &inputtb[ix2],&inputix[ix2],&inputmi[ix2],
          &inputrg[ix2],&inputpk[ix2],
          backwardtb,backwardix,backwardmi,backwardrg,backwardpk,deg1);

      for(int i=1; i<nvr-2; i++)                   // b[i] = b[i-1]*x[n-2-i]
      {
         ix2 = idx[nvr-2-i]*deg1; ix3 = i*deg1; ix1 = ix3 - deg1;
         dbl5_padded_convolute<<<1,BS>>>
            (&backwardtb[ix1],&backwardix[ix1],&backwardmi[ix1],
             &backwardrg[ix1],&backwardpk[ix1],
             &inputtb[ix2],&inputix[ix2],&inputmi[ix2],
             &inputrg[ix2],&inputpk[ix2],
             &backwardtb[ix3],&backwardix[ix3],&backwardmi[ix3],
             &backwardrg[ix3],&backwardpk[ix3],deg1);
      }
      ix3 = (nvr-3)*deg1; ix2 = (nvr-2)*deg1;         // b[n-2] = b[n-3]*cff
      dbl5_padded_convolute<<<1,BS>>>
         (&backwardtb[ix3],&backwardix[ix3],&backwardmi[ix3],
          &backwardrg[ix3],&backwardpk[ix3],
          cfftb,cffix,cffmi,cffrg,cffpk,
          &backwardtb[ix2],&backwardix[ix2],&backwardmi[ix2],
          &backwardrg[ix2],&backwardpk[ix2],deg1);

      if(nvr == 3)                                       // c[0] = f[0]*x[2]
      {
         ix2 = idx[2]*deg1;
         dbl5_padded_convolute<<<1,BS>>>
            (forwardtb,forwardix,forwardmi,forwardrg,forwardpk,
             &inputtb[ix2],&inputix[ix2],&inputmi[ix2],
             &inputrg[ix2],&inputpk[ix2],
             crosstb,crossix,crossmi,crossrg,crosspk,deg1);
      }
      else
      {
         for(int i=0; i<nvr-3; i++)                  // c[i] = f[i]*b[n-4-i]
         {
            ix1 = i*deg1; ix2 = (nvr-4-i)*deg1;
            dbl5_padded_convolute<<<1,BS>>>
               (&forwardtb[ix1],&forwardix[ix1],&forwardmi[ix1],
                &forwardrg[ix1],&forwardpk[ix1],
                &backwardtb[ix2],&backwardix[ix2],&backwardmi[ix2],
                &backwardrg[ix2],&backwardpk[ix2],
                &crosstb[ix1],&crossix[ix1],&crossmi[ix1],
                &crossrg[ix1],&crosspk[ix1],deg1);
         }
         ix1 = (nvr-3)*deg1; ix2 = idx[nvr-1]*deg1; // c[n-3] = f[n-3]*x[n-1]
         dbl5_padded_convolute<<<1,BS>>>
            (&forwardtb[ix1],&forwardix[ix1],&forwardmi[ix1],
             &forwardrg[ix1],&forwardpk[ix1],
             &inputtb[ix2],&inputix[ix2],&inputmi[ix2],
             &inputrg[ix2],&inputpk[ix2],
             &crosstb[ix1],&crossix[ix1],&crossmi[ix1],
             &crossrg[ix1],&crosspk[ix1],deg1);
      }
   }
}

void GPU_cmplx5_speel
 ( int BS, int nvr, int deg, int *idx,
   double *cffretb, double *cffreix, double *cffremi,
   double *cffrerg, double *cffrepk,
   double *cffimtb, double *cffimix, double *cffimmi,
   double *cffimrg, double *cffimpk,
   double *inputretb, double *inputreix, double *inputremi,
   double *inputrerg, double *inputrepk,
   double *inputimtb, double *inputimix, double *inputimmi,
   double *inputimrg, double *inputimpk,
   double *forwardretb, double *forwardreix, double *forwardremi,
   double *forwardrerg, double *forwardrepk,
   double *forwardimtb, double *forwardimix, double *forwardimmi,
   double *forwardimrg, double *forwardimpk,
   double *backwardretb, double *backwardreix, double *backwardremi,
   double *backwardrerg, double *backwardrepk,
   double *backwardimtb, double *backwardimix, double *backwardimmi,
   double *backwardimrg, double *backwardimpk,
   double *crossretb, double *crossreix, double *crossremi,
   double *crossrerg, double *crossrepk,
   double *crossimtb, double *crossimix, double *crossimmi,
   double *crossimrg, double *crossimpk )
{
   const int deg1 = deg+1;
   int ix1,ix2,ix3;

   ix1 = idx[0]*deg1;                                     // f[0] = cff*x[0]
   cmplx5_padded_convolute<<<1,BS>>>
      (cffretb,cffreix,cffremi,cffrerg,cffrepk,
       cffimtb,cffimix,cffimmi,cffimrg,cffimpk,
       &inputretb[ix1],&inputreix[ix1],&inputremi[ix1],
       &inputrerg[ix1],&inputrepk[ix1],
       &inputimtb[ix1],&inputimix[ix1],&inputimmi[ix1],
       &inputimrg[ix1],&inputimpk[ix1],
       forwardretb,forwardreix,forwardremi,forwardrerg,forwardrepk,
       forwardimtb,forwardimix,forwardimmi,forwardimrg,forwardimpk,deg1); 

   for(int i=1; i<nvr; i++)                            // f[i] = f[i-i]*x[i]
   {
      ix2 = idx[i]*deg1; ix3 = i*deg1; ix1 = ix3 - deg1;
      cmplx5_padded_convolute<<<1,BS>>>
         (&forwardretb[ix1],&forwardreix[ix1],&forwardremi[ix1],
          &forwardrerg[ix1],&forwardrepk[ix1],
          &forwardimtb[ix1],&forwardimix[ix1],&forwardimmi[ix1],
          &forwardimrg[ix1],&forwardimpk[ix1],
          &inputretb[ix2],&inputreix[ix2],&inputremi[ix2],
          &inputrerg[ix2],&inputrepk[ix2],
          &inputimtb[ix2],&inputimix[ix2],&inputimmi[ix2],
          &inputimrg[ix2],&inputimpk[ix2],
          &forwardretb[ix3],&forwardreix[ix3],&forwardremi[ix3],
          &forwardrerg[ix3],&forwardrepk[ix3],
          &forwardimtb[ix3],&forwardimix[ix3],&forwardimmi[ix3],
          &forwardimrg[ix3],&forwardimpk[ix3],deg1);
   }
   if(nvr > 2)
   {
      ix1 = idx[nvr-1]*deg1; ix2 = idx[nvr-2]*deg1;  // b[0] = x[n-1]*x[n-2]
      cmplx5_padded_convolute<<<1,BS>>>
         (&inputretb[ix1],&inputreix[ix1],&inputremi[ix1],
          &inputrerg[ix1],&inputrepk[ix1],
          &inputimtb[ix1],&inputimix[ix1],&inputimmi[ix1],
          &inputimrg[ix1],&inputimpk[ix1],
          &inputretb[ix2],&inputreix[ix2],&inputremi[ix2],
          &inputrerg[ix2],&inputrepk[ix2],
          &inputimtb[ix2],&inputimix[ix2],&inputimmi[ix2],
          &inputimrg[ix2],&inputimpk[ix2],
          backwardretb,backwardreix,backwardremi,backwardrerg,backwardrepk,
          backwardimtb,backwardimix,backwardimmi,backwardimrg,backwardimpk,
          deg1);

      for(int i=1; i<nvr-2; i++)                   // b[i] = b[i-1]*x[n-2-i]
      {
         ix2 = idx[nvr-2-i]*deg1; ix3 = i*deg1; ix1 = ix3 - deg1;
         cmplx5_padded_convolute<<<1,BS>>>
            (&backwardretb[ix1],&backwardreix[ix1],&backwardremi[ix1],
             &backwardrerg[ix1],&backwardrepk[ix1],
             &backwardimtb[ix1],&backwardimix[ix1],&backwardimmi[ix1],
             &backwardimrg[ix1],&backwardimpk[ix1],
             &inputretb[ix2],&inputreix[ix2],&inputremi[ix2],
             &inputrerg[ix2],&inputrepk[ix2],
             &inputimtb[ix2],&inputimix[ix2],&inputimmi[ix2],
             &inputimrg[ix2],&inputimpk[ix2],
             &backwardretb[ix3],&backwardreix[ix3],&backwardremi[ix3],
             &backwardrerg[ix3],&backwardrepk[ix3],
             &backwardimtb[ix3],&backwardimix[ix3],&backwardimmi[ix3],
             &backwardimrg[ix3],&backwardimpk[ix3],deg1);
      }
      ix3 = (nvr-3)*deg1; ix2 = (nvr-2)*deg1;         // b[n-2] = b[n-3]*cff
      cmplx5_padded_convolute<<<1,BS>>>
         (&backwardretb[ix3],&backwardreix[ix3],&backwardremi[ix3],
          &backwardrerg[ix3],&backwardrepk[ix3],
          &backwardimtb[ix3],&backwardimix[ix3],&backwardimmi[ix3],
          &backwardimrg[ix3],&backwardimpk[ix3],
          cffretb,cffreix,cffremi,cffrerg,cffrepk,
          cffimtb,cffimix,cffimmi,cffimrg,cffimpk,
          &backwardretb[ix2],&backwardreix[ix2],&backwardremi[ix2],
          &backwardrerg[ix2],&backwardrepk[ix2],
          &backwardimtb[ix2],&backwardimix[ix2],&backwardimmi[ix2],
          &backwardimrg[ix2],&backwardimpk[ix2],deg1);

      if(nvr == 3)                                       // c[0] = f[0]*x[2]
      {
         ix2 = idx[2]*deg1;
         cmplx5_padded_convolute<<<1,BS>>>
            (forwardretb,forwardreix,forwardremi,forwardrerg,forwardrepk,
             forwardimtb,forwardimix,forwardimmi,forwardimrg,forwardimpk,
             &inputretb[ix2],&inputreix[ix2],&inputremi[ix2],
             &inputrerg[ix2],&inputrepk[ix2],
             &inputimtb[ix2],&inputimix[ix2],&inputimmi[ix2],
             &inputimrg[ix2],&inputimpk[ix2],
             crossretb,crossreix,crossremi,crossrerg,crossrepk,
             crossimtb,crossimix,crossimmi,crossimrg,crossimpk,deg1);
      }
      else
      {
         for(int i=0; i<nvr-3; i++)                  // c[i] = f[i]*b[n-4-i]
         {
            ix1 = i*deg1; ix2 = (nvr-4-i)*deg1;
            cmplx5_padded_convolute<<<1,BS>>>
               (&forwardretb[ix1],&forwardreix[ix1],&forwardremi[ix1],
                &forwardrerg[ix1],&forwardrepk[ix1],
                &forwardimtb[ix1],&forwardimix[ix1],&forwardimmi[ix1],
                &forwardimrg[ix1],&forwardimpk[ix1],
                &backwardretb[ix2],&backwardreix[ix2],&backwardremi[ix2],
                &backwardrerg[ix2],&backwardrepk[ix2],
                &backwardimtb[ix2],&backwardimix[ix2],&backwardimmi[ix2],
                &backwardimrg[ix2],&backwardimpk[ix2],
                &crossretb[ix1],&crossreix[ix1],&crossremi[ix1],
                &crossrerg[ix1],&crossrepk[ix1],
                &crossimtb[ix1],&crossimix[ix1],&crossimmi[ix1],
                &crossimrg[ix1],&crossimpk[ix1],deg1);
         }
         ix1 = (nvr-3)*deg1; ix2 = idx[nvr-1]*deg1; // c[n-3] = f[n-3]*x[n-1]
         cmplx5_padded_convolute<<<1,BS>>>
            (&forwardretb[ix1],&forwardreix[ix1],&forwardremi[ix1],
             &forwardrerg[ix1],&forwardrepk[ix1],
             &forwardimtb[ix1],&forwardimix[ix1],&forwardimmi[ix1],
             &forwardimrg[ix1],&forwardimpk[ix1],
             &inputretb[ix2],&inputreix[ix2],&inputremi[ix2],
             &inputrerg[ix2],&inputrepk[ix2],
             &inputimtb[ix2],&inputimix[ix2],&inputimmi[ix2],
             &inputimrg[ix2],&inputimpk[ix2],
             &crossretb[ix1],&crossreix[ix1],&crossremi[ix1],
             &crossrerg[ix1],&crossrepk[ix1],
             &crossimtb[ix1],&crossimix[ix1],&crossimmi[ix1],
             &crossimrg[ix1],&crossimpk[ix1],deg1);
      }
   }
}

void GPU_dbl5_evaldiff
 ( int BS, int dim, int nvr, int deg, int *idx,
   double *cfftb, double *cffix, double *cffmi, double *cffrg, double *cffpk,
   double **inputtb, double **inputix, double **inputmi,
   double **inputrg, double **inputpk,
   double **outputtb, double **outputix, double **outputmi,
   double **outputrg, double **outputpk )
{
   const int deg1 = deg+1;      // length of all vectors
   double *inputtb_d;           // inputtb_d is inputtb on the device
   double *inputix_d;           // inputix_d is inputix on the device
   double *inputmi_d;           // inputmi_d is inputmi on the device
   double *inputrg_d;           // inputrg_d is inputrg on the device
   double *inputpk_d;           // inputpk_d is inputpk on the device
   double *forwardtb_d;         // highest forward products on the device
   double *forwardix_d;         // second highest forward products
   double *forwardmi_d;         // middle forward products
   double *forwardrg_d;         // second lowest forward products
   double *forwardpk_d;         // lowest forward products
   double *backwardtb_d;        // highest backward products on the device
   double *backwardix_d;        // second highest backward products
   double *backwardmi_d;        // middle backward products
   double *backwardrg_d;        // second lowest backward products
   double *backwardpk_d;        // lowest backward products
   double *crosstb_d;           // highest cross products on the device
   double *crossix_d;           // second highest cross products
   double *crossmi_d;           // middle cross products
   double *crossrg_d;           // second lowest cross products
   double *crosspk_d;           // lowest cross products
   double *cfftb_d;             // cfftb_d is cfftb on device
   double *cffix_d;             // cffix_d is cffix on device
   double *cffmi_d;             // cffmi_d is cffmi on device
   double *cffrg_d;             // cffrg_d is cffrg on device
   double *cffpk_d;             // cffpk_d is cffpk on device

   size_t szcff = deg1*sizeof(double);
   size_t szdim = dim*(deg1)*sizeof(double);
   size_t sznvr = nvr*(deg1)*sizeof(double);
   size_t sznvr1 = (nvr-1)*(deg1)*sizeof(double);
   size_t sznvr2 = (nvr-2)*(deg1)*sizeof(double);

   cudaMalloc((void**)&cfftb_d,szcff);
   cudaMalloc((void**)&cffix_d,szcff);
   cudaMalloc((void**)&cffmi_d,szcff);
   cudaMalloc((void**)&cffrg_d,szcff);
   cudaMalloc((void**)&cffpk_d,szcff);
   cudaMalloc((void**)&inputtb_d,szdim);
   cudaMalloc((void**)&inputix_d,szdim);
   cudaMalloc((void**)&inputmi_d,szdim);
   cudaMalloc((void**)&inputrg_d,szdim);
   cudaMalloc((void**)&inputpk_d,szdim);
   cudaMalloc((void**)&forwardtb_d,sznvr);
   cudaMalloc((void**)&forwardix_d,sznvr);
   cudaMalloc((void**)&forwardmi_d,sznvr);
   cudaMalloc((void**)&forwardrg_d,sznvr);
   cudaMalloc((void**)&forwardpk_d,sznvr);
   cudaMalloc((void**)&backwardtb_d,sznvr1);
   cudaMalloc((void**)&backwardix_d,sznvr1);
   cudaMalloc((void**)&backwardmi_d,sznvr1);
   cudaMalloc((void**)&backwardrg_d,sznvr1);
   cudaMalloc((void**)&backwardpk_d,sznvr1);
   cudaMalloc((void**)&crosstb_d,sznvr2);
   cudaMalloc((void**)&crossix_d,sznvr2);
   cudaMalloc((void**)&crossmi_d,sznvr2);
   cudaMalloc((void**)&crossrg_d,sznvr2);
   cudaMalloc((void**)&crosspk_d,sznvr2);

   double *inputtb_h = new double[dim*(deg1)];
   double *inputix_h = new double[dim*(deg1)];
   double *inputmi_h = new double[dim*(deg1)];
   double *inputrg_h = new double[dim*(deg1)];
   double *inputpk_h = new double[dim*(deg1)];
   int ix = 0;
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++)
      {
         inputtb_h[ix] = inputtb[i][j];
         inputix_h[ix] = inputix[i][j];
         inputmi_h[ix] = inputmi[i][j];
         inputrg_h[ix] = inputrg[i][j];
         inputpk_h[ix++] = inputpk[i][j];
      }

   cudaMemcpy(cfftb_d,cfftb,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffix_d,cffix,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffmi_d,cffmi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffrg_d,cffrg,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffpk_d,cffpk,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(inputtb_d,inputtb_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputix_d,inputix_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputmi_d,inputmi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputrg_d,inputrg_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputpk_d,inputpk_h,szdim,cudaMemcpyHostToDevice);

   if(BS == deg1)
   {
      GPU_dbl5_speel
         (BS,nvr,deg,idx,cfftb_d,cffix_d,cffmi_d,cffrg_d,cffpk_d,
          inputtb_d,inputix_d,inputmi_d,inputrg_d,inputpk_d,
          forwardtb_d,forwardix_d,forwardmi_d,forwardrg_d,forwardpk_d,
          backwardtb_d,backwardix_d,backwardmi_d,backwardrg_d,backwardpk_d,
          crosstb_d,crossix_d,crossmi_d,crossrg_d,crosspk_d);
   }
   double *forwardtb_h = new double[(deg1)*nvr];
   double *forwardix_h = new double[(deg1)*nvr];
   double *forwardmi_h = new double[(deg1)*nvr];
   double *forwardrg_h = new double[(deg1)*nvr];
   double *forwardpk_h = new double[(deg1)*nvr];
   double *backwardtb_h = new double[(deg1)*(nvr-1)];
   double *backwardix_h = new double[(deg1)*(nvr-1)];
   double *backwardmi_h = new double[(deg1)*(nvr-1)];
   double *backwardrg_h = new double[(deg1)*(nvr-1)];
   double *backwardpk_h = new double[(deg1)*(nvr-1)];
   double *crosstb_h = new double[(deg1)*(nvr-2)];
   double *crossix_h = new double[(deg1)*(nvr-2)];
   double *crossmi_h = new double[(deg1)*(nvr-2)];
   double *crossrg_h = new double[(deg1)*(nvr-2)];
   double *crosspk_h = new double[(deg1)*(nvr-2)];
  
   cudaMemcpy(forwardtb_h,forwardtb_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardix_h,forwardix_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardmi_h,forwardmi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardrg_h,forwardrg_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardpk_h,forwardpk_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardtb_h,backwardtb_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardix_h,backwardix_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardmi_h,backwardmi_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardrg_h,backwardrg_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardpk_h,backwardpk_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(crosstb_h,crosstb_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossix_h,crossix_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossmi_h,crossmi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossrg_h,crossrg_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crosspk_h,crosspk_d,sznvr2,cudaMemcpyDeviceToHost);

   int offset = (nvr-1)*deg1;            // assign value of the monomial
   for(int i=0; i<deg1; i++)
   {
      outputtb[dim][i] = forwardtb_h[offset+i];
      outputix[dim][i] = forwardix_h[offset+i];
      outputmi[dim][i] = forwardmi_h[offset+i];
      outputrg[dim][i] = forwardrg_h[offset+i];
      outputpk[dim][i] = forwardpk_h[offset+i];
   }
   ix = idx[nvr-1];                      // derivative with respect to x[n-1]
   offset = (nvr-2)*deg1;
   for(int i=0; i<deg1; i++)
   {
      outputtb[ix][i] = forwardtb_h[offset+i];
      outputix[ix][i] = forwardix_h[offset+i];
      outputmi[ix][i] = forwardmi_h[offset+i];
      outputrg[ix][i] = forwardrg_h[offset+i];
      outputpk[ix][i] = forwardpk_h[offset+i];
   }
   ix = idx[0];                          // derivative with respect to x[0]
   offset = (nvr-2)*deg1;
   for(int i=0; i<deg1; i++)
   {
      outputtb[ix][i] = backwardtb_h[offset+i];
      outputix[ix][i] = backwardix_h[offset+i];
      outputmi[ix][i] = backwardmi_h[offset+i];
      outputrg[ix][i] = backwardrg_h[offset+i];
      outputpk[ix][i] = backwardpk_h[offset+i];
   }
   for(int k=1; k<nvr-1; k++)            // derivative with respect to x[k]
   {
      ix = idx[k]; offset = (k-1)*deg1;
      for(int i=0; i<deg1; i++)
      {
         outputtb[ix][i] = crosstb_h[offset+i];
         outputix[ix][i] = crossix_h[offset+i];
         outputmi[ix][i] = crossmi_h[offset+i];
         outputrg[ix][i] = crossrg_h[offset+i];
         outputpk[ix][i] = crosspk_h[offset+i];
      }
   }
}

void GPU_cmplx5_evaldiff
 ( int BS, int dim, int nvr, int deg, int *idx,
   double *cffretb, double *cffreix, double *cffremi,
   double *cffrerg, double *cffrepk,
   double *cffimtb, double *cffimix, double *cffimmi,
   double *cffimrg, double *cffimpk,
   double **inputretb, double **inputreix, double **inputremi,
   double **inputrerg, double **inputrepk,
   double **inputimtb, double **inputimix, double **inputimmi,
   double **inputimrg, double **inputimpk,
   double **outputretb, double **outputreix, double **outputremi,
   double **outputrerg, double **outputrepk,
   double **outputimtb, double **outputimix, double **outputimmi,
   double **outputimrg, double **outputimpk )
{
   const int deg1 = deg+1;          // length of all vectors
   double *inputretb_d;             // inputretb_d is inputretb on the device
   double *inputreix_d;             // inputreix_d is inputreix on the device
   double *inputremi_d;             // inputremi_d is inputremi on the device
   double *inputrerg_d;             // inputrerg_d is inputrerg on the device
   double *inputrepk_d;             // inputrepk_d is inputrepk on the device
   double *inputimtb_d;             // inputimtb_d is inputretb on the device
   double *inputimix_d;             // inputimix_d is inputreix on the device
   double *inputimmi_d;             // inputimmi_d is inputremi on the device
   double *inputimrg_d;             // inputimrg_d is inputrerg on the device
   double *inputimpk_d;             // inputimpk_d is inputrepk on the device
   double *forwardretb_d;
   double *forwardreix_d;
   double *forwardremi_d;
   double *forwardrerg_d;
   double *forwardrepk_d;
   double *forwardimtb_d;
   double *forwardimix_d;
   double *forwardimmi_d;
   double *forwardimrg_d;
   double *forwardimpk_d;
   double *backwardretb_d;
   double *backwardreix_d;
   double *backwardremi_d;
   double *backwardrerg_d;
   double *backwardrepk_d;
   double *backwardimtb_d;
   double *backwardimix_d;
   double *backwardimmi_d;
   double *backwardimrg_d;
   double *backwardimpk_d;
   double *crossretb_d;
   double *crossreix_d;
   double *crossremi_d;
   double *crossrerg_d;
   double *crossrepk_d;
   double *crossimtb_d;
   double *crossimix_d;
   double *crossimmi_d;
   double *crossimrg_d;
   double *crossimpk_d;
   double *cffretb_d;               // cffretb_d is cffretb on the device
   double *cffreix_d;               // cffreix_d is cffreix on the device
   double *cffremi_d;               // cffremi_d is cffremi on the device
   double *cffrerg_d;               // cffremi_d is cffrerg on the device
   double *cffrepk_d;               // cffrepk_d is cffrepk on the device
   double *cffimtb_d;               // cffimtb_d is cffimtb on the device
   double *cffimix_d;               // cffimix_d is cffimix on the device
   double *cffimmi_d;               // cffimmi_d is cffimmi on the device
   double *cffimrg_d;               // cffimrg_d is cffimrg on the device
   double *cffimpk_d;               // cffimpk_d is cffimpk on the device

   size_t szdim = dim*(deg1)*sizeof(double);
   size_t sznvr = nvr*(deg1)*sizeof(double);
   size_t sznvr1 = (nvr-1)*(deg1)*sizeof(double);
   size_t sznvr2 = (nvr-2)*(deg1)*sizeof(double);
   size_t szcff = deg1*sizeof(double);

   cudaMalloc((void**)&cffretb_d,szcff);
   cudaMalloc((void**)&cffreix_d,szcff);
   cudaMalloc((void**)&cffremi_d,szcff);
   cudaMalloc((void**)&cffrerg_d,szcff);
   cudaMalloc((void**)&cffrepk_d,szcff);
   cudaMalloc((void**)&cffimtb_d,szcff);
   cudaMalloc((void**)&cffimix_d,szcff);
   cudaMalloc((void**)&cffimmi_d,szcff);
   cudaMalloc((void**)&cffimrg_d,szcff);
   cudaMalloc((void**)&cffimpk_d,szcff);
   cudaMalloc((void**)&inputretb_d,szdim);
   cudaMalloc((void**)&inputreix_d,szdim);
   cudaMalloc((void**)&inputremi_d,szdim);
   cudaMalloc((void**)&inputrerg_d,szdim);
   cudaMalloc((void**)&inputrepk_d,szdim);
   cudaMalloc((void**)&inputimtb_d,szdim);
   cudaMalloc((void**)&inputimix_d,szdim);
   cudaMalloc((void**)&inputimmi_d,szdim);
   cudaMalloc((void**)&inputimrg_d,szdim);
   cudaMalloc((void**)&inputimpk_d,szdim);
   cudaMalloc((void**)&forwardretb_d,sznvr);
   cudaMalloc((void**)&forwardreix_d,sznvr);
   cudaMalloc((void**)&forwardremi_d,sznvr);
   cudaMalloc((void**)&forwardrerg_d,sznvr);
   cudaMalloc((void**)&forwardrepk_d,sznvr);
   cudaMalloc((void**)&forwardimtb_d,sznvr);
   cudaMalloc((void**)&forwardimix_d,sznvr);
   cudaMalloc((void**)&forwardimmi_d,sznvr);
   cudaMalloc((void**)&forwardimrg_d,sznvr);
   cudaMalloc((void**)&forwardimpk_d,sznvr);
   cudaMalloc((void**)&backwardretb_d,sznvr1);
   cudaMalloc((void**)&backwardreix_d,sznvr1);
   cudaMalloc((void**)&backwardremi_d,sznvr1);
   cudaMalloc((void**)&backwardrerg_d,sznvr1);
   cudaMalloc((void**)&backwardrepk_d,sznvr1);
   cudaMalloc((void**)&backwardimtb_d,sznvr1);
   cudaMalloc((void**)&backwardimix_d,sznvr1);
   cudaMalloc((void**)&backwardimmi_d,sznvr1);
   cudaMalloc((void**)&backwardimrg_d,sznvr1);
   cudaMalloc((void**)&backwardimpk_d,sznvr1);
   cudaMalloc((void**)&crossretb_d,sznvr2);
   cudaMalloc((void**)&crossreix_d,sznvr2);
   cudaMalloc((void**)&crossremi_d,sznvr2);
   cudaMalloc((void**)&crossrerg_d,sznvr2);
   cudaMalloc((void**)&crossrepk_d,sznvr2);
   cudaMalloc((void**)&crossimtb_d,sznvr2);
   cudaMalloc((void**)&crossimix_d,sznvr2);
   cudaMalloc((void**)&crossimmi_d,sznvr2);
   cudaMalloc((void**)&crossimrg_d,sznvr2);
   cudaMalloc((void**)&crossimpk_d,sznvr2);

   double *inputretb_h = new double[dim*(deg1)];
   double *inputreix_h = new double[dim*(deg1)];
   double *inputremi_h = new double[dim*(deg1)];
   double *inputrerg_h = new double[dim*(deg1)];
   double *inputrepk_h = new double[dim*(deg1)];
   double *inputimtb_h = new double[dim*(deg1)];
   double *inputimix_h = new double[dim*(deg1)];
   double *inputimmi_h = new double[dim*(deg1)];
   double *inputimrg_h = new double[dim*(deg1)];
   double *inputimpk_h = new double[dim*(deg1)];
   int ix = 0;
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++)
      {
         inputretb_h[ix] = inputretb[i][j];
         inputreix_h[ix] = inputreix[i][j];
         inputremi_h[ix] = inputremi[i][j];
         inputrerg_h[ix] = inputrerg[i][j];
         inputrepk_h[ix] = inputrepk[i][j];
         inputimtb_h[ix] = inputimtb[i][j];
         inputimix_h[ix] = inputimix[i][j];
         inputimmi_h[ix] = inputimmi[i][j];
         inputimrg_h[ix] = inputimrg[i][j];
         inputimpk_h[ix++] = inputimpk[i][j];
      }

   cudaMemcpy(cffretb_d,cffretb,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffreix_d,cffreix,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffremi_d,cffremi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffrerg_d,cffrerg,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffrepk_d,cffrepk,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffimtb_d,cffimtb,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffimix_d,cffimix,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffimmi_d,cffimmi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffimrg_d,cffimrg,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffimpk_d,cffimpk,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(inputretb_d,inputretb_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputreix_d,inputreix_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputremi_d,inputremi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputrerg_d,inputrerg_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputrepk_d,inputrepk_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputimtb_d,inputimtb_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputimix_d,inputimix_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputimmi_d,inputimmi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputimrg_d,inputimrg_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputimpk_d,inputimpk_h,szdim,cudaMemcpyHostToDevice);

   if(BS == deg1)
   {
      GPU_cmplx5_speel(BS,nvr,deg,idx,
         cffretb_d,cffreix_d,cffremi_d,cffrerg_d,cffrepk_d,
         cffimtb_d,cffimix_d,cffimmi_d,cffimrg_d,cffimpk_d,
         inputretb_d,inputreix_d,inputremi_d,inputrerg_d,inputrepk_d,
         inputimtb_d,inputimix_d,inputimmi_d,inputimrg_d,inputimpk_d,
         forwardretb_d,forwardreix_d,forwardremi_d,
         forwardrerg_d,forwardrepk_d,
         forwardimtb_d,forwardimix_d,forwardimmi_d,
         forwardimrg_d,forwardimpk_d,
         backwardretb_d,backwardreix_d,backwardremi_d,
         backwardrerg_d,backwardrepk_d,
         backwardimtb_d,backwardimix_d,backwardimmi_d,
         backwardimrg_d,backwardimpk_d,
         crossretb_d,crossreix_d,crossremi_d,crossrerg_d,crossrepk_d,
         crossimtb_d,crossimix_d,crossimmi_d,crossimrg_d,crossimpk_d);
   }
   double *forwardretb_h = new double[(deg1)*nvr];
   double *forwardreix_h = new double[(deg1)*nvr];
   double *forwardremi_h = new double[(deg1)*nvr];
   double *forwardrerg_h = new double[(deg1)*nvr];
   double *forwardrepk_h = new double[(deg1)*nvr];
   double *forwardimtb_h = new double[(deg1)*nvr];
   double *forwardimix_h = new double[(deg1)*nvr];
   double *forwardimmi_h = new double[(deg1)*nvr];
   double *forwardimrg_h = new double[(deg1)*nvr];
   double *forwardimpk_h = new double[(deg1)*nvr];
   double *backwardretb_h = new double[(deg1)*(nvr-1)];
   double *backwardreix_h = new double[(deg1)*(nvr-1)];
   double *backwardremi_h = new double[(deg1)*(nvr-1)];
   double *backwardrerg_h = new double[(deg1)*(nvr-1)];
   double *backwardrepk_h = new double[(deg1)*(nvr-1)];
   double *backwardimtb_h = new double[(deg1)*(nvr-1)];
   double *backwardimix_h = new double[(deg1)*(nvr-1)];
   double *backwardimmi_h = new double[(deg1)*(nvr-1)];
   double *backwardimrg_h = new double[(deg1)*(nvr-1)];
   double *backwardimpk_h = new double[(deg1)*(nvr-1)];
   double *crossretb_h = new double[(deg1)*(nvr-2)];
   double *crossreix_h = new double[(deg1)*(nvr-2)];
   double *crossremi_h = new double[(deg1)*(nvr-2)];
   double *crossrerg_h = new double[(deg1)*(nvr-2)];
   double *crossrepk_h = new double[(deg1)*(nvr-2)];
   double *crossimtb_h = new double[(deg1)*(nvr-2)];
   double *crossimix_h = new double[(deg1)*(nvr-2)];
   double *crossimmi_h = new double[(deg1)*(nvr-2)];
   double *crossimrg_h = new double[(deg1)*(nvr-2)];
   double *crossimpk_h = new double[(deg1)*(nvr-2)];
  
   cudaMemcpy(forwardretb_h,forwardretb_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardreix_h,forwardreix_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardremi_h,forwardremi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardrerg_h,forwardrerg_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardrepk_h,forwardrepk_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardimtb_h,forwardimtb_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardimix_h,forwardimix_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardimmi_h,forwardimmi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardimrg_h,forwardimrg_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardimpk_h,forwardimpk_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardretb_h,backwardretb_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardreix_h,backwardreix_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardremi_h,backwardremi_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardrerg_h,backwardrerg_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardrepk_h,backwardrepk_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardimtb_h,backwardimtb_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardimix_h,backwardimix_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardimmi_h,backwardimmi_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardimrg_h,backwardimrg_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardimpk_h,backwardimpk_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossretb_h,crossretb_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossreix_h,crossreix_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossremi_h,crossremi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossrerg_h,crossrerg_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossrepk_h,crossrepk_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossimtb_h,crossimtb_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossimix_h,crossimix_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossimmi_h,crossimmi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossimrg_h,crossimrg_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossimpk_h,crossimpk_d,sznvr2,cudaMemcpyDeviceToHost);

   int offset = (nvr-1)*deg1;
   for(int i=0; i<deg1; i++)   // assign value of the monomial
   {
      outputretb[dim][i] = forwardretb_h[offset+i];
      outputreix[dim][i] = forwardreix_h[offset+i];
      outputremi[dim][i] = forwardremi_h[offset+i];
      outputrerg[dim][i] = forwardrerg_h[offset+i];
      outputrepk[dim][i] = forwardrepk_h[offset+i];
      outputimtb[dim][i] = forwardimtb_h[offset+i];
      outputimix[dim][i] = forwardimix_h[offset+i];
      outputimmi[dim][i] = forwardimmi_h[offset+i];
      outputimrg[dim][i] = forwardimrg_h[offset+i];
      outputimpk[dim][i] = forwardimpk_h[offset+i];
   }
   ix = idx[nvr-1];
   offset = (nvr-2)*deg1;
   for(int i=0; i<deg1; i++)  // derivative with respect to x[n-1]
   {
      outputretb[ix][i] = forwardretb_h[offset+i];
      outputreix[ix][i] = forwardreix_h[offset+i];
      outputremi[ix][i] = forwardremi_h[offset+i];
      outputrerg[ix][i] = forwardrerg_h[offset+i];
      outputrepk[ix][i] = forwardrepk_h[offset+i];
      outputimtb[ix][i] = forwardimtb_h[offset+i];
      outputimix[ix][i] = forwardimix_h[offset+i];
      outputimmi[ix][i] = forwardimmi_h[offset+i];
      outputimrg[ix][i] = forwardimrg_h[offset+i];
      outputimpk[ix][i] = forwardimpk_h[offset+i];
   }
   ix = idx[0]; 
   offset = (nvr-2)*deg1;
   for(int i=0; i<deg1; i++)   // derivative with respect to x[0]
   {
      outputretb[ix][i] = backwardretb_h[offset+i];
      outputreix[ix][i] = backwardreix_h[offset+i];
      outputremi[ix][i] = backwardremi_h[offset+i];
      outputrerg[ix][i] = backwardrerg_h[offset+i];
      outputrepk[ix][i] = backwardrepk_h[offset+i];
      outputimtb[ix][i] = backwardimtb_h[offset+i];
      outputimix[ix][i] = backwardimix_h[offset+i];
      outputimmi[ix][i] = backwardimmi_h[offset+i];
      outputimrg[ix][i] = backwardimrg_h[offset+i];
      outputimpk[ix][i] = backwardimpk_h[offset+i];
   }
   for(int k=1; k<nvr-1; k++)  // derivative with respect to x[k]
   {
      ix = idx[k]; offset = (k-1)*deg1;
      for(int i=0; i<deg1; i++)
      {
         outputretb[ix][i] = crossretb_h[offset+i];
         outputreix[ix][i] = crossreix_h[offset+i];
         outputremi[ix][i] = crossremi_h[offset+i];
         outputrerg[ix][i] = crossrerg_h[offset+i];
         outputrepk[ix][i] = crossrepk_h[offset+i];
         outputimtb[ix][i] = crossimtb_h[offset+i];
         outputimix[ix][i] = crossimix_h[offset+i];
         outputimmi[ix][i] = crossimmi_h[offset+i];
         outputimrg[ix][i] = crossimrg_h[offset+i];
         outputimpk[ix][i] = crossimpk_h[offset+i];
      }
   }
}
