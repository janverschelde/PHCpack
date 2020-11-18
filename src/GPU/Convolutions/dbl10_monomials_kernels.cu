// The file dbl10_monomials_kernels.cu defines the kernels specified
// in dbl10_monomials_kernels.h.

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

#include "dbl10_convolutions_kernels.cu"
#include "dbl10_monomials_kernels.h"

void GPU_dbl10_speel
 ( int BS, int nvr, int deg, int *idx,
   double *cffrtb, double *cffrix, double *cffrmi, double *cffrrg,
   double *cffrpk, double *cffltb, double *cfflix, double *cfflmi,
   double *cfflrg, double *cfflpk,
   double *inputrtb, double *inputrix, double *inputrmi, double *inputrrg,
   double *inputrpk, double *inputltb, double *inputlix, double *inputlmi,
   double *inputlrg, double *inputlpk,
   double *forwardrtb, double *forwardrix, double *forwardrmi,
   double *forwardrrg, double *forwardrpk,
   double *forwardltb, double *forwardlix, double *forwardlmi,
   double *forwardlrg, double *forwardlpk,
   double *backwardrtb, double *backwardrix, double *backwardrmi,
   double *backwardrrg, double *backwardrpk,
   double *backwardltb, double *backwardlix, double *backwardlmi,
   double *backwardlrg, double *backwardlpk,
   double *crossrtb, double *crossrix, double *crossrmi, double *crossrrg,
   double *crossrpk, double *crossltb, double *crosslix, double *crosslmi,
   double *crosslrg, double *crosslpk )
{
   const int deg1 = deg+1;
   int ix1,ix2,ix3;

   ix1 = idx[0]*deg1;                                     // f[0] = cff*x[0]
   dbl10_convolute<<<1,BS>>>
      (cffrtb,cffrix,cffrmi,cffrrg,cffrpk,cffltb,cfflix,cfflmi,cfflrg,cfflpk,
       &inputrtb[ix1],&inputrix[ix1],&inputrmi[ix1],
       &inputrrg[ix1],&inputrpk[ix1],
       &inputltb[ix1],&inputlix[ix1],&inputlmi[ix1],
       &inputlrg[ix1],&inputlpk[ix1],
       forwardrtb,forwardrix,forwardrmi,forwardrrg,forwardrpk,
       forwardltb,forwardlix,forwardlmi,forwardlrg,forwardlpk,deg1);

   for(int i=1; i<nvr; i++)                            // f[i] = f[i-1]*x[i]
   {
      ix2 = idx[i]*deg1; ix3 = i*deg1; ix1 = ix3 - deg1;
      dbl10_convolute<<<1,BS>>>
         (&forwardrtb[ix1],&forwardrix[ix1],&forwardrmi[ix1],
          &forwardrrg[ix1],&forwardrpk[ix1],
          &forwardltb[ix1],&forwardlix[ix1],&forwardlmi[ix1],
          &forwardlrg[ix1],&forwardlpk[ix1],
          &inputrtb[ix2],&inputrix[ix2],&inputrmi[ix2],
          &inputrrg[ix2],&inputrpk[ix2],
          &inputltb[ix2],&inputlix[ix2],&inputlmi[ix2],
          &inputlrg[ix2],&inputlpk[ix2],
          &forwardrtb[ix3],&forwardrix[ix3],&forwardrmi[ix3],
          &forwardrrg[ix3],&forwardrpk[ix3],
          &forwardltb[ix3],&forwardlix[ix3],&forwardlmi[ix3],
          &forwardlrg[ix3],&forwardlpk[ix3],deg1);
   }
   if(nvr > 2)
   {
      ix1 = idx[nvr-1]*deg1; ix2 = idx[nvr-2]*deg1;  // b[0] = x[n-1]*x[n-2]
      dbl10_convolute<<<1,BS>>>
         (&inputrtb[ix1],&inputrix[ix1],&inputrmi[ix1],
          &inputrrg[ix1],&inputrpk[ix1],
          &inputltb[ix1],&inputlix[ix1],&inputlmi[ix1],
          &inputlrg[ix1],&inputlpk[ix1],
          &inputrtb[ix2],&inputrix[ix2],&inputrmi[ix2],
          &inputrrg[ix2],&inputrpk[ix2],
          &inputltb[ix2],&inputlix[ix2],&inputlmi[ix2],
          &inputlrg[ix2],&inputlpk[ix2],
          backwardrtb,backwardrix,backwardrmi,backwardrrg,backwardrpk,
          backwardltb,backwardlix,backwardlmi,backwardlrg,backwardlpk,deg1);

      for(int i=1; i<nvr-2; i++)                   // b[i] = b[i-1]*x[n-2-i]
      {
         ix2 = idx[nvr-2-i]*deg1; ix3 = i*deg1; ix1 = ix3 - deg1;
         dbl10_convolute<<<1,BS>>>
            (&backwardrtb[ix1],&backwardrix[ix1],&backwardrmi[ix1],
             &backwardrrg[ix1],&backwardrpk[ix1],
             &backwardltb[ix1],&backwardlix[ix1],&backwardlmi[ix1],
             &backwardlrg[ix1],&backwardlpk[ix1],
             &inputrtb[ix2],&inputrix[ix2],&inputrmi[ix2],
             &inputrrg[ix2],&inputrpk[ix2],
             &inputltb[ix2],&inputlix[ix2],&inputlmi[ix2],
             &inputlrg[ix2],&inputlpk[ix2],
             &backwardrtb[ix3],&backwardrix[ix3],&backwardrmi[ix3],
             &backwardrrg[ix3],&backwardrpk[ix3],
             &backwardltb[ix3],&backwardlix[ix3],&backwardlmi[ix3],
             &backwardlrg[ix3],&backwardlpk[ix3],deg1);
      }
      ix3 = (nvr-3)*deg1; ix2 = (nvr-2)*deg1;         // b[n-2] = b[n-3]*cff
      dbl10_convolute<<<1,BS>>>
         (&backwardrtb[ix3],&backwardrix[ix3],&backwardrmi[ix3],
          &backwardrrg[ix3],&backwardrpk[ix3],
          &backwardltb[ix3],&backwardlix[ix3],&backwardlmi[ix3],
          &backwardlrg[ix3],&backwardlpk[ix3],
          cffrtb,cffrix,cffrmi,cffrrg,cffrpk,
          cffltb,cfflix,cfflmi,cfflrg,cfflpk,
          &backwardrtb[ix2],&backwardrix[ix2],&backwardrmi[ix2],
          &backwardrrg[ix2],&backwardrpk[ix2],
          &backwardltb[ix2],&backwardlix[ix2],&backwardlmi[ix2],
          &backwardlrg[ix2],&backwardlpk[ix2],deg1);

      if(nvr == 3)                                       // c[0] = f[0]*x[2]
      {
         ix2 = idx[2]*deg1;
         dbl10_convolute<<<1,BS>>>
            (forwardrtb,forwardrix,forwardrmi,forwardrrg,forwardrpk,
             forwardltb,forwardlix,forwardlmi,forwardlrg,forwardlpk,
             &inputrtb[ix2],&inputrix[ix2],&inputrmi[ix2],
             &inputrrg[ix2],&inputrpk[ix2],
             &inputltb[ix2],&inputlix[ix2],&inputlmi[ix2],
             &inputlrg[ix2],&inputlpk[ix2],
             crossrtb,crossrix,crossrmi,crossrrg,crossrpk,
             crossltb,crosslix,crosslmi,crosslrg,crosslpk,deg1);
      }
      else
      {
         for(int i=0; i<nvr-3; i++)                  // c[i] = f[i]*b[n-4-i]
         {
            ix1 = i*deg1; ix2 = (nvr-4-i)*deg1;
            dbl10_convolute<<<1,BS>>>
               (&forwardrtb[ix1],&forwardrix[ix1],&forwardrmi[ix1],
                &forwardrrg[ix1],&forwardrpk[ix1],
                &forwardltb[ix1],&forwardlix[ix1],&forwardlmi[ix1],
                &forwardlrg[ix1],&forwardlpk[ix1],
                &backwardrtb[ix2],&backwardrix[ix2],&backwardrmi[ix2],
                &backwardrrg[ix2],&backwardrpk[ix2],
                &backwardltb[ix2],&backwardlix[ix2],&backwardlmi[ix2],
                &backwardlrg[ix2],&backwardlpk[ix2],
                &crossrtb[ix1],&crossrix[ix1],&crossrmi[ix1],
                &crossrrg[ix1],&crossrpk[ix1],
                &crossltb[ix1],&crosslix[ix1],&crosslmi[ix1],
                &crosslrg[ix1],&crosslpk[ix1],deg1);
         }
         ix1 = (nvr-3)*deg1; ix2 = idx[nvr-1]*deg1; // c[n-3] = f[n-3]*x[n-1]
         dbl10_convolute<<<1,BS>>>
            (&forwardrtb[ix1],&forwardrix[ix1],&forwardrmi[ix1],
             &forwardrrg[ix1],&forwardrpk[ix1],
             &forwardltb[ix1],&forwardlix[ix1],&forwardlmi[ix1],
             &forwardlrg[ix1],&forwardlpk[ix1],
             &inputrtb[ix2],&inputrix[ix2],&inputrmi[ix2],
             &inputrrg[ix2],&inputrpk[ix2],
             &inputltb[ix2],&inputlix[ix2],&inputlmi[ix2],
             &inputlrg[ix2],&inputlpk[ix2],
             &crossrtb[ix1],&crossrix[ix1],&crossrmi[ix1],
             &crossrrg[ix1],&crossrpk[ix1],
             &crossltb[ix1],&crosslix[ix1],&crosslmi[ix1],
             &crosslrg[ix1],&crosslpk[ix1],deg1);
      }
   }
}

void GPU_cmplx10_speel
 ( int BS, int nvr, int deg, int *idx,
   double *cffrertb, double *cffrerix, double *cffrermi, double *cffrerrg,
   double *cffrerpk, double *cffreltb, double *cffrelix, double *cffrelmi,
   double *cffrelrg, double *cffrelpk,
   double *cffimrtb, double *cffimrix, double *cffimrmi, double *cffimrrg,
   double *cffimrpk, double *cffimltb, double *cffimlix, double *cffimlmi,
   double *cffimlrg, double *cffimlpk,
   double *inputrertb, double *inputrerix, double *inputrermi,
   double *inputrerrg, double *inputrerpk,
   double *inputreltb, double *inputrelix, double *inputrelmi,
   double *inputrelrg, double *inputrelpk,
   double *inputimrtb, double *inputimrix, double *inputimrmi,
   double *inputimrrg, double *inputimrpk,
   double *inputimltb, double *inputimlix, double *inputimlmi,
   double *inputimlrg, double *inputimlpk,
   double *forwardrertb, double *forwardrerix, double *forwardrermi,
   double *forwardrerrg, double *forwardrerpk,
   double *forwardreltb, double *forwardrelix, double *forwardrelmi,
   double *forwardrelrg, double *forwardrelpk,
   double *forwardimrtb, double *forwardimrix, double *forwardimrmi,
   double *forwardimrrg, double *forwardimrpk,
   double *forwardimltb, double *forwardimlix, double *forwardimlmi,
   double *forwardimlrg, double *forwardimlpk,
   double *backwardrertb, double *backwardrerix, double *backwardrermi,
   double *backwardrerrg, double *backwardrerpk,
   double *backwardreltb, double *backwardrelix, double *backwardrelmi,
   double *backwardrelrg, double *backwardrelpk,
   double *backwardimrtb, double *backwardimrix, double *backwardimrmi,
   double *backwardimrrg, double *backwardimrpk,
   double *backwardimltb, double *backwardimlix, double *backwardimlmi,
   double *backwardimlrg, double *backwardimlpk,
   double *crossrertb, double *crossrerix, double *crossrermi,
   double *crossrerrg, double *crossrerpk,
   double *crossreltb, double *crossrelix, double *crossrelmi,
   double *crossrelrg, double *crossrelpk,
   double *crossimrtb, double *crossimrix, double *crossimrmi,
   double *crossimrrg, double *crossimrpk,
   double *crossimltb, double *crossimlix, double *crossimlmi,
   double *crossimlrg, double *crossimlpk )
{
}

void GPU_dbl10_evaldiff
 ( int BS, int dim, int nvr, int deg, int *idx,
   double *cffrtb, double *cffrix, double *cffrmi, double *cffrrg,
   double *cffrpk, double *cffltb, double *cfflix, double *cfflmi,
   double *cfflrg, double *cfflpk,
   double **inputrtb, double **inputrix, double **inputrmi,
   double **inputrrg, double **inputrpk,
   double **inputltb, double **inputlix, double **inputlmi,
   double **inputlrg, double **inputlpk,
   double **outputrtb, double **outputrix, double **outputrmi,
   double **outputrrg, double **outputrpk,
   double **outputltb, double **outputlix, double **outputlmi,
   double **outputlrg, double **outputlpk )
{
   const int deg1 = deg+1;      // length of all vectors
   double *inputrtb_d;          // inputrtb_d is inputrtb on the device
   double *inputrix_d;          // inputrix_d is inputrix on the device
   double *inputrmi_d;          // inputrmi_d is inputrmi on the device
   double *inputrrg_d;          // inputrrg_d is inputrrg on the device
   double *inputrpk_d;          // inputrpk_d is inputrpk on the device
   double *inputltb_d;          // inputltb_d is inputltb on the device
   double *inputlix_d;          // inputlix_d is inputlix on the device
   double *inputlmi_d;          // inputlmi_d is inputlmi on the device
   double *inputlrg_d;          // inputlrg_d is inputlrg on the device
   double *inputlpk_d;          // inputlpk_d is inputlpk on the device
   double *forwardrtb_d;        // highest forward products on the device
   double *forwardrix_d;        // second highest forward products
   double *forwardrmi_d;        // third highest forward products
   double *forwardrrg_d;        // fourth highest forward products
   double *forwardrpk_d;        // fifth highest forward products
   double *forwardltb_d;        // fifth lowest forward products
   double *forwardlix_d;        // fourth lowest forward products
   double *forwardlmi_d;        // third lowest forward products
   double *forwardlrg_d;        // second lowest forward products
   double *forwardlpk_d;        // lowest forward products
   double *backwardrtb_d;       // highest backward products on the device
   double *backwardrix_d;       // second highest backward products
   double *backwardrmi_d;       // third highest backward products
   double *backwardrrg_d;       // fourth highest backward products
   double *backwardrpk_d;       // fifth highest backward products
   double *backwardltb_d;       // fifth lowest backward products
   double *backwardlix_d;       // fourth lowest backward products
   double *backwardlmi_d;       // third lowest backward products
   double *backwardlrg_d;       // second lowest backward products
   double *backwardlpk_d;       // lowest backward products
   double *crossrtb_d;          // highest cross products on the device
   double *crossrix_d;          // second highest cross products
   double *crossrmi_d;          // third highest cross products
   double *crossrrg_d;          // fourth highest cross products
   double *crossrpk_d;          // fifth highest cross products
   double *crossltb_d;          // fifth lowest cross products
   double *crosslix_d;          // fourth lowest cross products
   double *crosslmi_d;          // third lowest cross products
   double *crosslrg_d;          // second lowest cross products
   double *crosslpk_d;          // lowest cross products
   double *cffrtb_d;            // cffrtb_d is cffrtb on device
   double *cffrix_d;            // cffrix_d is cffrix on device
   double *cffrmi_d;            // cffrmi_d is cffrmi on device
   double *cffrrg_d;            // cffrrg_d is cffrrg on device
   double *cffrpk_d;            // cffrpk_d is cffrpk on device
   double *cffltb_d;            // cffltb_d is cffltb on device
   double *cfflix_d;            // cfflix_d is cfflix on device
   double *cfflmi_d;            // cfflmi_d is cfflmi on device
   double *cfflrg_d;            // cfflrg_d is cfflrg on device
   double *cfflpk_d;            // cfflpk_d is cfflpk on device

   size_t szcff = deg1*sizeof(double);
   size_t szdim = dim*(deg1)*sizeof(double);
   size_t sznvr = nvr*(deg1)*sizeof(double);
   size_t sznvr1 = (nvr-1)*(deg1)*sizeof(double);
   size_t sznvr2 = (nvr-2)*(deg1)*sizeof(double);
   size_t szidx = nvr*sizeof(int);

   cudaMalloc((void**)&cffrtb_d,szcff);
   cudaMalloc((void**)&cffrix_d,szcff);
   cudaMalloc((void**)&cffrmi_d,szcff);
   cudaMalloc((void**)&cffrrg_d,szcff);
   cudaMalloc((void**)&cffrpk_d,szcff);
   cudaMalloc((void**)&cffltb_d,szcff);
   cudaMalloc((void**)&cfflix_d,szcff);
   cudaMalloc((void**)&cfflmi_d,szcff);
   cudaMalloc((void**)&cfflrg_d,szcff);
   cudaMalloc((void**)&cfflpk_d,szcff);
   cudaMalloc((void**)&inputrtb_d,szdim);
   cudaMalloc((void**)&inputrix_d,szdim);
   cudaMalloc((void**)&inputrmi_d,szdim);
   cudaMalloc((void**)&inputrrg_d,szdim);
   cudaMalloc((void**)&inputrpk_d,szdim);
   cudaMalloc((void**)&inputltb_d,szdim);
   cudaMalloc((void**)&inputlix_d,szdim);
   cudaMalloc((void**)&inputlmi_d,szdim);
   cudaMalloc((void**)&inputlrg_d,szdim);
   cudaMalloc((void**)&inputlpk_d,szdim);
   cudaMalloc((void**)&forwardrtb_d,sznvr);
   cudaMalloc((void**)&forwardrix_d,sznvr);
   cudaMalloc((void**)&forwardrmi_d,sznvr);
   cudaMalloc((void**)&forwardrrg_d,sznvr);
   cudaMalloc((void**)&forwardrpk_d,sznvr);
   cudaMalloc((void**)&forwardltb_d,sznvr);
   cudaMalloc((void**)&forwardlix_d,sznvr);
   cudaMalloc((void**)&forwardlmi_d,sznvr);
   cudaMalloc((void**)&forwardlrg_d,sznvr);
   cudaMalloc((void**)&forwardlpk_d,sznvr);
   cudaMalloc((void**)&backwardrtb_d,sznvr1);
   cudaMalloc((void**)&backwardrix_d,sznvr1);
   cudaMalloc((void**)&backwardrmi_d,sznvr1);
   cudaMalloc((void**)&backwardrrg_d,sznvr1);
   cudaMalloc((void**)&backwardrpk_d,sznvr1);
   cudaMalloc((void**)&backwardltb_d,sznvr1);
   cudaMalloc((void**)&backwardlix_d,sznvr1);
   cudaMalloc((void**)&backwardlmi_d,sznvr1);
   cudaMalloc((void**)&backwardlrg_d,sznvr1);
   cudaMalloc((void**)&backwardlpk_d,sznvr1);
   cudaMalloc((void**)&crossrtb_d,sznvr2);
   cudaMalloc((void**)&crossrix_d,sznvr2);
   cudaMalloc((void**)&crossrmi_d,sznvr2);
   cudaMalloc((void**)&crossrrg_d,sznvr2);
   cudaMalloc((void**)&crossrpk_d,sznvr2);
   cudaMalloc((void**)&crossltb_d,sznvr2);
   cudaMalloc((void**)&crosslix_d,sznvr2);
   cudaMalloc((void**)&crosslmi_d,sznvr2);
   cudaMalloc((void**)&crosslrg_d,sznvr2);
   cudaMalloc((void**)&crosslpk_d,sznvr2);

   double *inputrtb_h = new double[dim*(deg1)];
   double *inputrix_h = new double[dim*(deg1)];
   double *inputrmi_h = new double[dim*(deg1)];
   double *inputrrg_h = new double[dim*(deg1)];
   double *inputrpk_h = new double[dim*(deg1)];
   double *inputltb_h = new double[dim*(deg1)];
   double *inputlix_h = new double[dim*(deg1)];
   double *inputlmi_h = new double[dim*(deg1)];
   double *inputlrg_h = new double[dim*(deg1)];
   double *inputlpk_h = new double[dim*(deg1)];
   int ix = 0;
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++)
      {
         inputrtb_h[ix] = inputrtb[i][j];
         inputrix_h[ix] = inputrix[i][j];
         inputrmi_h[ix] = inputrmi[i][j];
         inputrrg_h[ix] = inputrrg[i][j];
         inputrpk_h[ix] = inputrpk[i][j];
         inputltb_h[ix] = inputltb[i][j];
         inputlix_h[ix] = inputlix[i][j];
         inputlmi_h[ix] = inputlmi[i][j];
         inputlrg_h[ix] = inputlrg[i][j];
         inputlpk_h[ix++] = inputlpk[i][j];
      }

   cudaMemcpy(cffrtb_d,cffrtb,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffrix_d,cffrix,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffrmi_d,cffrmi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffrrg_d,cffrrg,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffrpk_d,cffrpk,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffltb_d,cffltb,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cfflix_d,cfflix,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cfflmi_d,cfflmi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cfflrg_d,cfflrg,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cfflpk_d,cfflpk,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(inputrtb_d,inputrtb_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputrix_d,inputrix_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputrmi_d,inputrmi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputrrg_d,inputrrg_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputrpk_d,inputrpk_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputltb_d,inputltb_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputlix_d,inputlix_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputlmi_d,inputlmi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputlrg_d,inputlrg_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputlpk_d,inputlpk_h,szdim,cudaMemcpyHostToDevice);

   if(BS == deg1)
   {
      GPU_dbl10_speel
         (BS,nvr,deg,idx,
          cffrtb_d,cffrix_d,cffrmi_d,cffrrg_d,cffrpk_d,
          cffltb_d,cfflix_d,cfflmi_d,cfflrg_d,cfflpk_d,
          inputrtb_d,inputrix_d,inputrmi_d,inputrrg_d,inputrpk_d,
          inputltb_d,inputlix_d,inputlmi_d,inputlrg_d,inputlpk_d,
          forwardrtb_d,forwardrix_d,forwardrmi_d,forwardrrg_d,forwardrpk_d,
          forwardltb_d,forwardlix_d,forwardlmi_d,forwardlrg_d,forwardlpk_d,
          backwardrtb_d,backwardrix_d,backwardrmi_d,backwardrrg_d,
          backwardrpk_d,backwardltb_d,backwardlix_d,backwardlmi_d,
          backwardlrg_d,backwardlpk_d,
          crossrtb_d,crossrix_d,crossrmi_d,crossrrg_d,crossrpk_d,
          crossltb_d,crosslix_d,crosslmi_d,crosslrg_d,crosslpk_d);
   }
   double *forwardrtb_h = new double[(deg1)*nvr];
   double *forwardrix_h = new double[(deg1)*nvr];
   double *forwardrmi_h = new double[(deg1)*nvr];
   double *forwardrrg_h = new double[(deg1)*nvr];
   double *forwardrpk_h = new double[(deg1)*nvr];
   double *forwardltb_h = new double[(deg1)*nvr];
   double *forwardlix_h = new double[(deg1)*nvr];
   double *forwardlmi_h = new double[(deg1)*nvr];
   double *forwardlrg_h = new double[(deg1)*nvr];
   double *forwardlpk_h = new double[(deg1)*nvr];
   double *backwardrtb_h = new double[(deg1)*(nvr-1)];
   double *backwardrix_h = new double[(deg1)*(nvr-1)];
   double *backwardrmi_h = new double[(deg1)*(nvr-1)];
   double *backwardrrg_h = new double[(deg1)*(nvr-1)];
   double *backwardrpk_h = new double[(deg1)*(nvr-1)];
   double *backwardltb_h = new double[(deg1)*(nvr-1)];
   double *backwardlix_h = new double[(deg1)*(nvr-1)];
   double *backwardlmi_h = new double[(deg1)*(nvr-1)];
   double *backwardlrg_h = new double[(deg1)*(nvr-1)];
   double *backwardlpk_h = new double[(deg1)*(nvr-1)];
   double *crossrtb_h = new double[(deg1)*(nvr-2)];
   double *crossrix_h = new double[(deg1)*(nvr-2)];
   double *crossrmi_h = new double[(deg1)*(nvr-2)];
   double *crossrrg_h = new double[(deg1)*(nvr-2)];
   double *crossrpk_h = new double[(deg1)*(nvr-2)];
   double *crossltb_h = new double[(deg1)*(nvr-2)];
   double *crosslix_h = new double[(deg1)*(nvr-2)];
   double *crosslmi_h = new double[(deg1)*(nvr-2)];
   double *crosslrg_h = new double[(deg1)*(nvr-2)];
   double *crosslpk_h = new double[(deg1)*(nvr-2)];
  
   cudaMemcpy(forwardrtb_h,forwardrtb_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardrix_h,forwardrix_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardrmi_h,forwardrmi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardrrg_h,forwardrrg_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardrpk_h,forwardrpk_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardltb_h,forwardltb_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardlix_h,forwardlix_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardlmi_h,forwardlmi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardlrg_h,forwardlrg_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardlpk_h,forwardlpk_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardrtb_h,backwardrtb_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardrix_h,backwardrix_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardrmi_h,backwardrmi_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardrrg_h,backwardrrg_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardrpk_h,backwardrpk_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardltb_h,backwardltb_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardlix_h,backwardlix_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardlmi_h,backwardlmi_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardlrg_h,backwardlrg_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardlpk_h,backwardlpk_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossrtb_h,crossrtb_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossrix_h,crossrix_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossrmi_h,crossrmi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossrrg_h,crossrrg_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossrpk_h,crossrpk_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossltb_h,crossltb_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crosslix_h,crosslix_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crosslmi_h,crosslmi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crosslrg_h,crosslrg_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crosslpk_h,crosslpk_d,sznvr2,cudaMemcpyDeviceToHost);

   int offset = (nvr-1)*deg1;            // assign value of the monomial
   for(int i=0; i<deg1; i++)
   {
      outputrtb[dim][i] = forwardrtb_h[offset+i];
      outputrix[dim][i] = forwardrix_h[offset+i];
      outputrmi[dim][i] = forwardrmi_h[offset+i];
      outputrrg[dim][i] = forwardrrg_h[offset+i];
      outputrpk[dim][i] = forwardrpk_h[offset+i];
      outputltb[dim][i] = forwardltb_h[offset+i];
      outputlix[dim][i] = forwardlix_h[offset+i];
      outputlmi[dim][i] = forwardlmi_h[offset+i];
      outputlrg[dim][i] = forwardlrg_h[offset+i];
      outputlpk[dim][i] = forwardlpk_h[offset+i];
   }
   ix = idx[nvr-1];                      // derivative with respect to x[n-1]
   offset = (nvr-2)*deg1;
   for(int i=0; i<deg1; i++)
   {
      outputrtb[ix][i] = forwardrtb_h[offset+i];
      outputrix[ix][i] = forwardrix_h[offset+i];
      outputrmi[ix][i] = forwardrmi_h[offset+i];
      outputrrg[ix][i] = forwardrrg_h[offset+i];
      outputrpk[ix][i] = forwardrpk_h[offset+i];
      outputltb[ix][i] = forwardltb_h[offset+i];
      outputlix[ix][i] = forwardlix_h[offset+i];
      outputlmi[ix][i] = forwardlmi_h[offset+i];
      outputlrg[ix][i] = forwardlrg_h[offset+i];
      outputlpk[ix][i] = forwardlpk_h[offset+i];
   }
   ix = idx[0];                          // derivative with respect to x[0]
   offset = (nvr-2)*deg1;
   for(int i=0; i<deg1; i++)
   {
      outputrtb[ix][i] = backwardrtb_h[offset+i];
      outputrix[ix][i] = backwardrix_h[offset+i];
      outputrmi[ix][i] = backwardrmi_h[offset+i];
      outputrrg[ix][i] = backwardrrg_h[offset+i];
      outputrpk[ix][i] = backwardrpk_h[offset+i];
      outputltb[ix][i] = backwardltb_h[offset+i];
      outputlix[ix][i] = backwardlix_h[offset+i];
      outputlmi[ix][i] = backwardlmi_h[offset+i];
      outputlrg[ix][i] = backwardlrg_h[offset+i];
      outputlpk[ix][i] = backwardlpk_h[offset+i];
   }
   for(int k=1; k<nvr-1; k++)            // derivative with respect to x[k]
   {
      ix = idx[k]; offset = (k-1)*deg1;
      for(int i=0; i<deg1; i++)
      {
         outputrtb[ix][i] = crossrtb_h[offset+i];
         outputrix[ix][i] = crossrix_h[offset+i];
         outputrmi[ix][i] = crossrmi_h[offset+i];
         outputrrg[ix][i] = crossrrg_h[offset+i];
         outputrpk[ix][i] = crossrpk_h[offset+i];
         outputltb[ix][i] = crossltb_h[offset+i];
         outputlix[ix][i] = crosslix_h[offset+i];
         outputlmi[ix][i] = crosslmi_h[offset+i];
         outputlrg[ix][i] = crosslrg_h[offset+i];
         outputlpk[ix][i] = crosslpk_h[offset+i];
      }
   }
}

void GPU_cmplx10_evaldiff
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
}
