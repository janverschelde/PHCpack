// The file dbl3_monomials_kernels.cu defines the kernels specified
// in dbl3_monomials_kernels.h.

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

#include "dbl3_convolutions_kernels.cu"
#include "dbl3_monomials_kernels.h"

void GPU_dbl3_speel
 ( int BS, int nvr, int deg, int *idx,
   double *cffhi, double *cffmi, double *cfflo,
   double *inputhi, double *inputmi, double *inputlo,
   double *forwardhi, double *forwardmi, double *forwardlo,
   double *backwardhi, double *backwardmi, double *backwardlo,
   double *crosshi, double *crossmi, double *crosslo )
{
   const int deg1 = deg+1;
   int ix1,ix2,ix3;

   ix1 = idx[0]*deg1;                                     // f[0] = cff*x[0]
   dbl3_padded_convolute<<<1,BS>>>
      (cffhi,cffmi,cfflo,&inputhi[ix1],&inputmi[ix1],&inputlo[ix1],
       forwardhi,forwardmi,forwardlo,deg1);

   for(int i=1; i<nvr; i++)                            // f[i] = f[i-1]*x[i]
   {
      ix2 = idx[i]*deg1; ix3 = i*deg1; ix1 = ix3 - deg1;
      dbl3_padded_convolute<<<1,BS>>>
         (&forwardhi[ix1],&forwardmi[ix1],&forwardlo[ix1],
          &inputhi[ix2],&inputmi[ix2],&inputlo[ix2],
          &forwardhi[ix3],&forwardmi[ix3],&forwardlo[ix3],deg1);
   }
   if(nvr > 2)
   {
      ix1 = idx[nvr-1]*deg1; ix2 = idx[nvr-2]*deg1;  // b[0] = x[n-1]*x[n-2]
      dbl3_padded_convolute<<<1,BS>>>
         (&inputhi[ix1],&inputmi[ix1],&inputlo[ix1],
          &inputhi[ix2],&inputmi[ix2],&inputlo[ix2],
          backwardhi,backwardmi,backwardlo,deg1);

      for(int i=1; i<nvr-2; i++)                   // b[i] = b[i-1]*x[n-2-i]
      {
         ix2 = idx[nvr-2-i]*deg1; ix3 = i*deg1; ix1 = ix3 - deg1;
         dbl3_padded_convolute<<<1,BS>>>
            (&backwardhi[ix1],&backwardmi[ix1],&backwardlo[ix1],
             &inputhi[ix2],&inputmi[ix2],&inputlo[ix2],
             &backwardhi[ix3],&backwardmi[ix3],&backwardlo[ix3],deg1);
      }
      ix3 = (nvr-3)*deg1; ix2 = (nvr-2)*deg1;         // b[n-2] = b[n-3]*cff
      dbl3_padded_convolute<<<1,BS>>>
         (&backwardhi[ix3],&backwardmi[ix3],&backwardlo[ix3],
          cffhi,cffmi,cfflo,
          &backwardhi[ix2],&backwardmi[ix2],&backwardlo[ix2],deg1);

      if(nvr == 3)                                       // c[0] = f[0]*x[2]
      {
         ix2 = idx[2]*deg1;
         dbl3_padded_convolute<<<1,BS>>>
            (forwardhi,forwardmi,forwardlo,
             &inputhi[ix2],&inputmi[ix2],&inputlo[ix2],
             crosshi,crossmi,crosslo,deg1);
      }
      else
      {
         for(int i=0; i<nvr-3; i++)                  // c[i] = f[i]*b[n-4-i]
         {
            ix1 = i*deg1; ix2 = (nvr-4-i)*deg1;
            dbl3_padded_convolute<<<1,BS>>>
               (&forwardhi[ix1],&forwardmi[ix1],&forwardlo[ix1],
                &backwardhi[ix2],&backwardmi[ix2],&backwardlo[ix2],
                &crosshi[ix1],&crossmi[ix1],&crosslo[ix1],deg1);
         }
         ix1 = (nvr-3)*deg1; ix2 = idx[nvr-1]*deg1; // c[n-3] = f[n-3]*x[n-1]
         dbl3_padded_convolute<<<1,BS>>>
            (&forwardhi[ix1],&forwardmi[ix1],&forwardlo[ix1],
             &inputhi[ix2],&inputmi[ix2],&inputlo[ix2],
             &crosshi[ix1],&crossmi[ix1],&crosslo[ix1],deg1);
      }
   }
}

void GPU_cmplx3_speel
 ( int BS, int nvr, int deg, int *idx,
   double *cffrehi, double *cffremi, double *cffrelo,
   double *cffimhi, double *cffimmi, double *cffimlo,
   double *inputrehi, double *inputremi, double *inputrelo,
   double *inputimhi, double *inputimmi, double *inputimlo,
   double *forwardrehi, double *forwardremi, double *forwardrelo,
   double *forwardimhi, double *forwardimmi, double *forwardimlo,
   double *backwardrehi, double *backwardremi, double *backwardrelo,
   double *backwardimhi, double *backwardimmi, double *backwardimlo,
   double *crossrehi, double *crossremi, double *crossrelo,
   double *crossimhi, double *crossimmi, double *crossimlo )
{
   const int deg1 = deg+1;
   int ix1,ix2,ix3;

   ix1 = idx[0]*deg1;                                     // f[0] = cff*x[0]
   cmplx3_padded_convolute<<<1,BS>>>
      (cffrehi,cffremi,cffrelo,cffimhi,cffimmi,cffimlo,
       &inputrehi[ix1],&inputremi[ix1],&inputrelo[ix1],
       &inputimhi[ix1],&inputimmi[ix1],&inputimlo[ix1],
       forwardrehi,forwardremi,forwardrelo,
       forwardimhi,forwardimmi,forwardimlo,deg1); 

   for(int i=1; i<nvr; i++)                            // f[i] = f[i-i]*x[i]
   {
      ix2 = idx[i]*deg1; ix3 = i*deg1; ix1 = ix3 - deg1;
      cmplx3_padded_convolute<<<1,BS>>>
         (&forwardrehi[ix1],&forwardremi[ix1],&forwardrelo[ix1],
          &forwardimhi[ix1],&forwardimmi[ix1],&forwardimlo[ix1],
          &inputrehi[ix2],&inputremi[ix2],&inputrelo[ix2],
          &inputimhi[ix2],&inputimmi[ix2],&inputimlo[ix2],
          &forwardrehi[ix3],&forwardremi[ix3],&forwardrelo[ix3],
          &forwardimhi[ix3],&forwardimmi[ix3],&forwardimlo[ix3],deg1);
   }
   if(nvr > 2)
   {
      ix1 = idx[nvr-1]*deg1; ix2 = idx[nvr-2]*deg1;  // b[0] = x[n-1]*x[n-2]
      cmplx3_padded_convolute<<<1,BS>>>
         (&inputrehi[ix1],&inputremi[ix1],&inputrelo[ix1],
          &inputimhi[ix1],&inputimmi[ix1],&inputimlo[ix1],
          &inputrehi[ix2],&inputremi[ix2],&inputrelo[ix2],
          &inputimhi[ix2],&inputimmi[ix2],&inputimlo[ix2],
          backwardrehi,backwardremi,backwardrelo,
          backwardimhi,backwardimmi,backwardimlo,deg1);

      for(int i=1; i<nvr-2; i++)                   // b[i] = b[i-1]*x[n-2-i]
      {
         ix2 = idx[nvr-2-i]*deg1; ix3 = i*deg1; ix1 = ix3 - deg1;
         cmplx3_padded_convolute<<<1,BS>>>
            (&backwardrehi[ix1],&backwardremi[ix1],&backwardrelo[ix1],
             &backwardimhi[ix1],&backwardimmi[ix1],&backwardimlo[ix1],
             &inputrehi[ix2],&inputremi[ix2],&inputrelo[ix2],
             &inputimhi[ix2],&inputimmi[ix2],&inputimlo[ix2],
             &backwardrehi[ix3],&backwardremi[ix3],&backwardrelo[ix3],
             &backwardimhi[ix3],&backwardimmi[ix3],&backwardimlo[ix3],deg1);
      }
      ix3 = (nvr-3)*deg1; ix2 = (nvr-2)*deg1;         // b[n-2] = b[n-3]*cff
      cmplx3_padded_convolute<<<1,BS>>>
         (&backwardrehi[ix3],&backwardremi[ix3],&backwardrelo[ix3],
          &backwardimhi[ix3],&backwardimmi[ix3],&backwardimlo[ix3],
          cffrehi,cffremi,cffrelo,cffimhi,cffimmi,cffimlo,
          &backwardrehi[ix2],&backwardremi[ix2],&backwardrelo[ix2],
          &backwardimhi[ix2],&backwardimmi[ix2],&backwardimlo[ix2],deg1);

      if(nvr == 3)                                       // c[0] = f[0]*x[2]
      {
         ix2 = idx[2]*deg1;
         cmplx3_padded_convolute<<<1,BS>>>
            (forwardrehi,forwardremi,forwardrelo,
             forwardimhi,forwardimmi,forwardimlo,
             &inputrehi[ix2],&inputremi[ix2],&inputrelo[ix2],
             &inputimhi[ix2],&inputimmi[ix2],&inputimlo[ix2],
             crossrehi,crossremi,crossrelo,
             crossimhi,crossimmi,crossimlo,deg1);
      }
      else
      {
         for(int i=0; i<nvr-3; i++)                  // c[i] = f[i]*b[n-4-i]
         {
            ix1 = i*deg1; ix2 = (nvr-4-i)*deg1;
            cmplx3_padded_convolute<<<1,BS>>>
               (&forwardrehi[ix1],&forwardremi[ix1],&forwardrelo[ix1],
                &forwardimhi[ix1],&forwardimmi[ix1],&forwardimlo[ix1],
                &backwardrehi[ix2],&backwardremi[ix2],&backwardrelo[ix2],
                &backwardimhi[ix2],&backwardimmi[ix2],&backwardimlo[ix2],
                &crossrehi[ix1],&crossremi[ix1],&crossrelo[ix1],
                &crossimhi[ix1],&crossimmi[ix1],&crossimlo[ix1],deg1);
         }
         ix1 = (nvr-3)*deg1; ix2 = idx[nvr-1]*deg1; // c[n-3] = f[n-3]*x[n-1]
         cmplx3_padded_convolute<<<1,BS>>>
            (&forwardrehi[ix1],&forwardremi[ix1],&forwardrelo[ix1],
             &forwardimhi[ix1],&forwardimmi[ix1],&forwardimlo[ix1],
             &inputrehi[ix2],&inputremi[ix2],&inputrelo[ix2],
             &inputimhi[ix2],&inputimmi[ix2],&inputimlo[ix2],
             &crossrehi[ix1],&crossremi[ix1],&crossrelo[ix1],
             &crossimhi[ix1],&crossimmi[ix1],&crossimlo[ix1],deg1);
      }
   }
}

void GPU_dbl3_evaldiff
 ( int BS, int dim, int nvr, int deg, int *idx,
   double *cffhi, double *cffmi, double *cfflo,
   double **inputhi, double **inputmi, double **inputlo,
   double **outputhi, double **outputmi, double **outputlo )
{
   const int deg1 = deg+1;            // length of all vectors
   double *inputhi_d;                 // inputhi_d is input on the device
   double *inputmi_d;                 // inputmi_d is input on the device
   double *inputlo_d;                 // inputlo_d is input on the device
   double *forwardhi_d;               // high forward products on the device
   double *forwardmi_d;               // middle forward products on the device
   double *forwardlo_d;               // low forward products on the device
   double *backwardhi_d;              // high backward products on the device
   double *backwardmi_d;              // middle backward products on the device
   double *backwardlo_d;              // low backward products on the device
   double *crosshi_d;                 // high cross products on the device
   double *crossmi_d;                 // middle cross products on the device
   double *crosslo_d;                 // low cross products on the device
   double *cffhi_d;                   // cffhi_d is cffhi on device
   double *cffmi_d;                   // cffmi_d is cffmi on device
   double *cfflo_d;                   // cfflo_d is cfflo on device

   size_t szcff = deg1*sizeof(double);
   size_t szdim = dim*(deg1)*sizeof(double);
   size_t sznvr = nvr*(deg1)*sizeof(double);
   size_t sznvr1 = (nvr-1)*(deg1)*sizeof(double);
   size_t sznvr2 = (nvr-2)*(deg1)*sizeof(double);

   cudaMalloc((void**)&cffhi_d,szcff);
   cudaMalloc((void**)&cffmi_d,szcff);
   cudaMalloc((void**)&cfflo_d,szcff);
   cudaMalloc((void**)&inputhi_d,szdim);
   cudaMalloc((void**)&inputmi_d,szdim);
   cudaMalloc((void**)&inputlo_d,szdim);
   cudaMalloc((void**)&forwardhi_d,sznvr);
   cudaMalloc((void**)&forwardmi_d,sznvr);
   cudaMalloc((void**)&forwardlo_d,sznvr);
   cudaMalloc((void**)&backwardhi_d,sznvr1);
   cudaMalloc((void**)&backwardmi_d,sznvr1);
   cudaMalloc((void**)&backwardlo_d,sznvr1);
   cudaMalloc((void**)&crosshi_d,sznvr2);
   cudaMalloc((void**)&crossmi_d,sznvr2);
   cudaMalloc((void**)&crosslo_d,sznvr2);

   double *inputhi_h = new double[dim*(deg1)];
   double *inputmi_h = new double[dim*(deg1)];
   double *inputlo_h = new double[dim*(deg1)];
   int ix = 0;
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++)
      {
         inputhi_h[ix] = inputhi[i][j];
         inputmi_h[ix] = inputmi[i][j];
         inputlo_h[ix++] = inputlo[i][j];
      }

   cudaMemcpy(cffhi_d,cffhi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffmi_d,cffmi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cfflo_d,cfflo,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(inputhi_d,inputhi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputmi_d,inputmi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputlo_d,inputlo_h,szdim,cudaMemcpyHostToDevice);

   if(BS == deg1)
   {
      GPU_dbl3_speel
         (BS,nvr,deg,idx,cffhi_d,cffmi_d,cfflo_d,inputhi_d,inputmi_d,
          inputlo_d,forwardhi_d,forwardmi_d,forwardlo_d,backwardhi_d,
          backwardmi_d,backwardlo_d,crosshi_d,crossmi_d,crosslo_d);
   }
   double *forwardhi_h = new double[(deg1)*nvr];
   double *forwardmi_h = new double[(deg1)*nvr];
   double *forwardlo_h = new double[(deg1)*nvr];
   double *backwardhi_h = new double[(deg1)*(nvr-1)];
   double *backwardmi_h = new double[(deg1)*(nvr-1)];
   double *backwardlo_h = new double[(deg1)*(nvr-1)];
   double *crosshi_h = new double[(deg1)*(nvr-2)];
   double *crossmi_h = new double[(deg1)*(nvr-2)];
   double *crosslo_h = new double[(deg1)*(nvr-2)];
  
   cudaMemcpy(forwardhi_h,forwardhi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardmi_h,forwardmi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardlo_h,forwardlo_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardhi_h,backwardhi_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardmi_h,backwardmi_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardlo_h,backwardlo_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(crosshi_h,crosshi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossmi_h,crossmi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crosslo_h,crosslo_d,sznvr2,cudaMemcpyDeviceToHost);

   int offset = (nvr-1)*deg1;            // assign value of the monomial
   for(int i=0; i<deg1; i++)
   {
      outputhi[dim][i] = forwardhi_h[offset+i];
      outputmi[dim][i] = forwardmi_h[offset+i];
      outputlo[dim][i] = forwardlo_h[offset+i];
   }
   ix = idx[nvr-1];                      // derivative with respect to x[n-1]
   offset = (nvr-2)*deg1;
   for(int i=0; i<deg1; i++)
   {
      outputhi[ix][i] = forwardhi_h[offset+i];
      outputmi[ix][i] = forwardmi_h[offset+i];
      outputlo[ix][i] = forwardlo_h[offset+i];
   }
   ix = idx[0];                          // derivative with respect to x[0]
   offset = (nvr-2)*deg1;
   for(int i=0; i<deg1; i++)
   {
      outputhi[ix][i] = backwardhi_h[offset+i];
      outputmi[ix][i] = backwardmi_h[offset+i];
      outputlo[ix][i] = backwardlo_h[offset+i];
   }
   for(int k=1; k<nvr-1; k++)            // derivative with respect to x[k]
   {
      ix = idx[k]; offset = (k-1)*deg1;
      for(int i=0; i<deg1; i++)
      {
         outputhi[ix][i] = crosshi_h[offset+i];
         outputmi[ix][i] = crossmi_h[offset+i];
         outputlo[ix][i] = crosslo_h[offset+i];
      }
   }
}

void GPU_cmplx3_evaldiff
 ( int BS, int dim, int nvr, int deg, int *idx,
   double *cffrehi, double *cffremi, double *cffrelo,
   double *cffimhi, double *cffimmi, double *cffimlo,
   double **inputrehi, double **inputremi, double **inputrelo,
   double **inputimhi, double **inputimmi, double **inputimlo,
   double **outputrehi, double **outputremi, double **outputrelo,
   double **outputimhi, double **outputimmi, double **outputimlo )
{
   const int deg1 = deg+1;          // length of all vectors
   double *inputrehi_d;             // inputrehi_d is inputrehi on the device
   double *inputremi_d;             // inputremi_d is inputremi on the device
   double *inputrelo_d;             // inputrelo_d is inputrelo on the device
   double *inputimhi_d;             // inputimhi_d is inputrehi on the device
   double *inputimmi_d;             // inputimmi_d is inputremi on the device
   double *inputimlo_d;             // inputimlo_d is inputrelo on the device
   double *forwardrehi_d;
   double *forwardremi_d;
   double *forwardrelo_d;
   double *forwardimhi_d;
   double *forwardimmi_d;
   double *forwardimlo_d;
   double *backwardrehi_d;
   double *backwardremi_d;
   double *backwardrelo_d;
   double *backwardimhi_d;
   double *backwardimmi_d;
   double *backwardimlo_d;
   double *crossrehi_d;
   double *crossremi_d;
   double *crossrelo_d;
   double *crossimhi_d;
   double *crossimmi_d;
   double *crossimlo_d;
   double *cffrehi_d;               // cffrehi_d is cffrehi on the device
   double *cffremi_d;               // cffremi_d is cffremi on the device
   double *cffrelo_d;               // cffrelo_d is cffrelo on the device
   double *cffimhi_d;               // cffimhi_d is cffimhi on the device
   double *cffimmi_d;               // cffimmi_d is cffimmi on the device
   double *cffimlo_d;               // cffimlo_d is cffimlo on the device

   size_t szdim = dim*(deg1)*sizeof(double);
   size_t sznvr = nvr*(deg1)*sizeof(double);
   size_t sznvr1 = (nvr-1)*(deg1)*sizeof(double);
   size_t sznvr2 = (nvr-2)*(deg1)*sizeof(double);
   size_t szcff = deg1*sizeof(double);

   cudaMalloc((void**)&cffrehi_d,szcff);
   cudaMalloc((void**)&cffremi_d,szcff);
   cudaMalloc((void**)&cffrelo_d,szcff);
   cudaMalloc((void**)&cffimhi_d,szcff);
   cudaMalloc((void**)&cffimmi_d,szcff);
   cudaMalloc((void**)&cffimlo_d,szcff);
   cudaMalloc((void**)&inputrehi_d,szdim);
   cudaMalloc((void**)&inputremi_d,szdim);
   cudaMalloc((void**)&inputrelo_d,szdim);
   cudaMalloc((void**)&inputimhi_d,szdim);
   cudaMalloc((void**)&inputimmi_d,szdim);
   cudaMalloc((void**)&inputimlo_d,szdim);
   cudaMalloc((void**)&forwardrehi_d,sznvr);
   cudaMalloc((void**)&forwardremi_d,sznvr);
   cudaMalloc((void**)&forwardrelo_d,sznvr);
   cudaMalloc((void**)&forwardimhi_d,sznvr);
   cudaMalloc((void**)&forwardimmi_d,sznvr);
   cudaMalloc((void**)&forwardimlo_d,sznvr);
   cudaMalloc((void**)&backwardrehi_d,sznvr1);
   cudaMalloc((void**)&backwardremi_d,sznvr1);
   cudaMalloc((void**)&backwardrelo_d,sznvr1);
   cudaMalloc((void**)&backwardimhi_d,sznvr1);
   cudaMalloc((void**)&backwardimmi_d,sznvr1);
   cudaMalloc((void**)&backwardimlo_d,sznvr1);
   cudaMalloc((void**)&crossrehi_d,sznvr2);
   cudaMalloc((void**)&crossremi_d,sznvr2);
   cudaMalloc((void**)&crossrelo_d,sznvr2);
   cudaMalloc((void**)&crossimhi_d,sznvr2);
   cudaMalloc((void**)&crossimmi_d,sznvr2);
   cudaMalloc((void**)&crossimlo_d,sznvr2);

   double *inputrehi_h = new double[dim*(deg1)];
   double *inputremi_h = new double[dim*(deg1)];
   double *inputrelo_h = new double[dim*(deg1)];
   double *inputimhi_h = new double[dim*(deg1)];
   double *inputimmi_h = new double[dim*(deg1)];
   double *inputimlo_h = new double[dim*(deg1)];
   int ix = 0;
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++)
      {
         inputrehi_h[ix] = inputrehi[i][j];
         inputremi_h[ix] = inputremi[i][j];
         inputrelo_h[ix] = inputrelo[i][j];
         inputimhi_h[ix] = inputimhi[i][j];
         inputimmi_h[ix] = inputimmi[i][j];
         inputimlo_h[ix++] = inputimlo[i][j];
      }

   cudaMemcpy(cffrehi_d,cffrehi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffremi_d,cffremi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffrelo_d,cffrelo,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffimhi_d,cffimhi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffimmi_d,cffimmi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffimlo_d,cffimlo,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(inputrehi_d,inputrehi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputremi_d,inputremi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputrelo_d,inputrelo_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputimhi_d,inputimhi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputimmi_d,inputimmi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputimlo_d,inputimlo_h,szdim,cudaMemcpyHostToDevice);

   if(BS == deg1)
   {
      GPU_cmplx3_speel(BS,nvr,deg,idx,
         cffrehi_d,cffremi_d,cffrelo_d,cffimhi_d,cffimmi_d,cffimlo_d,
         inputrehi_d,inputremi_d,inputrelo_d,
         inputimhi_d,inputimmi_d,inputimlo_d,
         forwardrehi_d,forwardremi_d,forwardrelo_d,
         forwardimhi_d,forwardimmi_d,forwardimlo_d,
         backwardrehi_d,backwardremi_d,backwardrelo_d,
         backwardimhi_d,backwardimmi_d,backwardimlo_d,
         crossrehi_d,crossremi_d,crossrelo_d,
         crossimhi_d,crossimmi_d,crossimlo_d);
   }
   double *forwardrehi_h = new double[(deg1)*nvr];
   double *forwardremi_h = new double[(deg1)*nvr];
   double *forwardrelo_h = new double[(deg1)*nvr];
   double *forwardimhi_h = new double[(deg1)*nvr];
   double *forwardimmi_h = new double[(deg1)*nvr];
   double *forwardimlo_h = new double[(deg1)*nvr];
   double *backwardrehi_h = new double[(deg1)*(nvr-1)];
   double *backwardremi_h = new double[(deg1)*(nvr-1)];
   double *backwardrelo_h = new double[(deg1)*(nvr-1)];
   double *backwardimhi_h = new double[(deg1)*(nvr-1)];
   double *backwardimmi_h = new double[(deg1)*(nvr-1)];
   double *backwardimlo_h = new double[(deg1)*(nvr-1)];
   double *crossrehi_h = new double[(deg1)*(nvr-2)];
   double *crossremi_h = new double[(deg1)*(nvr-2)];
   double *crossrelo_h = new double[(deg1)*(nvr-2)];
   double *crossimhi_h = new double[(deg1)*(nvr-2)];
   double *crossimmi_h = new double[(deg1)*(nvr-2)];
   double *crossimlo_h = new double[(deg1)*(nvr-2)];
  
   cudaMemcpy(forwardrehi_h,forwardrehi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardremi_h,forwardremi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardrelo_h,forwardrelo_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardimhi_h,forwardimhi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardimmi_h,forwardimmi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardimlo_h,forwardimlo_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardrehi_h,backwardrehi_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardremi_h,backwardremi_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardrelo_h,backwardrelo_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardimhi_h,backwardimhi_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardimmi_h,backwardimmi_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardimlo_h,backwardimlo_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossrehi_h,crossrehi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossremi_h,crossremi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossrelo_h,crossrelo_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossimhi_h,crossimhi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossimmi_h,crossimmi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossimlo_h,crossimlo_d,sznvr2,cudaMemcpyDeviceToHost);

   int offset = (nvr-1)*deg1;
   for(int i=0; i<deg1; i++)   // assign value of the monomial
   {
      outputrehi[dim][i] = forwardrehi_h[offset+i];
      outputremi[dim][i] = forwardremi_h[offset+i];
      outputrelo[dim][i] = forwardrelo_h[offset+i];
      outputimhi[dim][i] = forwardimhi_h[offset+i];
      outputimmi[dim][i] = forwardimmi_h[offset+i];
      outputimlo[dim][i] = forwardimlo_h[offset+i];
   }
   ix = idx[nvr-1];
   offset = (nvr-2)*deg1;
   for(int i=0; i<deg1; i++)  // derivative with respect to x[n-1]
   {
      outputrehi[ix][i] = forwardrehi_h[offset+i];
      outputremi[ix][i] = forwardremi_h[offset+i];
      outputrelo[ix][i] = forwardrelo_h[offset+i];
      outputimhi[ix][i] = forwardimhi_h[offset+i];
      outputimmi[ix][i] = forwardimmi_h[offset+i];
      outputimlo[ix][i] = forwardimlo_h[offset+i];
   }
   ix = idx[0]; 
   offset = (nvr-2)*deg1;
   for(int i=0; i<deg1; i++)   // derivative with respect to x[0]
   {
      outputrehi[ix][i] = backwardrehi_h[offset+i];
      outputremi[ix][i] = backwardremi_h[offset+i];
      outputrelo[ix][i] = backwardrelo_h[offset+i];
      outputimhi[ix][i] = backwardimhi_h[offset+i];
      outputimmi[ix][i] = backwardimmi_h[offset+i];
      outputimlo[ix][i] = backwardimlo_h[offset+i];
   }
   for(int k=1; k<nvr-1; k++)  // derivative with respect to x[k]
   {
      ix = idx[k]; offset = (k-1)*deg1;
      for(int i=0; i<deg1; i++)
      {
         outputrehi[ix][i] = crossrehi_h[offset+i];
         outputremi[ix][i] = crossremi_h[offset+i];
         outputrelo[ix][i] = crossrelo_h[offset+i];
         outputimhi[ix][i] = crossimhi_h[offset+i];
         outputimmi[ix][i] = crossimmi_h[offset+i];
         outputimlo[ix][i] = crossimlo_h[offset+i];
      }
   }
}
