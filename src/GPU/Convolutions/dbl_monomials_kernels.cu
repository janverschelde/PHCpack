// The file dbl_monomials_kernels.cu defines the kernels with prototypes
// in dbl_monomials_kernels.h.

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

#include "dbl_convolutions_kernels.cu"
#include "dbl_monomials_kernels.h"

void GPU_dbl_speel
 ( int BS, int nvr, int deg, int *idx, double *cff, double *input,
   double *forward, double *backward, double *cross )
{
   const int deg1 = deg+1;
   int ix1,ix2,ix3;

   ix1 = idx[0]*deg1;                                     // f[0] = cff*x[0]
   dbl_padded_convolute<<<1,BS>>>(cff,&input[ix1],forward,deg1);

   for(int i=1; i<nvr; i++)                            // f[i] = f[i-1]*x[i]
   {
      ix2 = idx[i]*deg1; ix3 = i*deg1; ix1 = ix3 - deg1;
      dbl_padded_convolute<<<1,BS>>>
         (&forward[ix1],&input[ix2],&forward[ix3],deg1);
   }
   if(nvr > 2)
   {
      ix1 = idx[nvr-1]*deg1; ix2 = idx[nvr-2]*deg1;  // b[0] = x[n-1]*x[n-2]
      dbl_padded_convolute<<<1,BS>>>(&input[ix1],&input[ix2],backward,deg1);

      for(int i=1; i<nvr-2; i++)                   // b[i] = b[i-1]*x[n-2-i]
      {
         ix2 = idx[nvr-2-i]*deg1; ix3 = i*deg1; ix1 = ix3 - deg1;
         dbl_padded_convolute<<<1,BS>>>
            (&backward[ix1],&input[ix2],&backward[ix3],deg1);
      }
      ix3 = (nvr-3)*deg1; ix2 = (nvr-2)*deg1;         // b[n-2] = b[n-3]*cff
      dbl_padded_convolute<<<1,BS>>>(&backward[ix3],cff,&backward[ix2],deg1);

      if(nvr == 3)                                       // c[0] = f[0]*x[2]
      {
         ix2 = idx[2]*deg1;
         dbl_padded_convolute<<<1,BS>>>(forward,&input[ix2],cross,deg1);
      }
      else
      {
         for(int i=0; i<nvr-3; i++)                  // c[i] = f[i]*b[n-4-i]
         {
            ix1 = i*deg1; ix2 = (nvr-4-i)*deg1;
            dbl_padded_convolute<<<1,BS>>>
               (&forward[ix1],&backward[ix2],&cross[ix1],deg1);
         }
         ix1 = (nvr-3)*deg1; ix2 = idx[nvr-1]*deg1; // c[n-3] = f[n-3]*x[n-1]
         dbl_padded_convolute<<<1,BS>>>
            (&forward[ix1],&input[ix2],&cross[ix1],deg1);
      }
   }
}

void GPU_cmplx_speel
 ( int BS, int nvr, int deg, int *idx, double *cffre, double *cffim,
   double *inputre, double *inputim, double *forwardre, double *forwardim,
   double *backwardre, double *backwardim, double *crossre, double *crossim )
{
   const int deg1 = deg+1;
   int ix1,ix2,ix3;

   ix1 = idx[0]*deg1;                                     // f[0] = cff*x[0]
   cmplx_padded_convolute<<<1,BS>>>
      (cffre,cffim,&inputre[ix1],&inputim[ix1],forwardre,forwardim,deg1); 

   for(int i=1; i<nvr; i++)                            // f[i] = f[i-i]*x[i]
   {
      ix2 = idx[i]*deg1; ix3 = i*deg1; ix1 = ix3 - deg1;
      cmplx_padded_convolute<<<1,BS>>>
         (&forwardre[ix1],&forwardim[ix1],&inputre[ix2],&inputim[ix2],
          &forwardre[ix3],&forwardim[ix3],deg1);
   }
   if(nvr > 2)
   {
      ix1 = idx[nvr-1]*deg1; ix2 = idx[nvr-2]*deg1;  // b[0] = x[n-1]*x[n-2]
      cmplx_padded_convolute<<<1,BS>>>
         (&inputre[ix1],&inputim[ix1],&inputre[ix2],&inputim[ix2],
          backwardre,backwardim,deg1);

      for(int i=1; i<nvr-2; i++)                   // b[i] = b[i-1]*x[n-2-i]
      {
         ix2 = idx[nvr-2-i]*deg1; ix3 = i*deg1; ix1 = ix3 - deg1;
         cmplx_padded_convolute<<<1,BS>>>
            (&backwardre[ix1],&backwardim[ix1],&inputre[ix2],&inputim[ix2],
             &backwardre[ix3],&backwardim[ix3],deg1);
      }
      ix3 = (nvr-3)*deg1; ix2 = (nvr-2)*deg1;         // b[n-2] = b[n-3]*cff
      cmplx_padded_convolute<<<1,BS>>>
         (&backwardre[ix3],&backwardim[ix3],cffre,cffim,
          &backwardre[ix2],&backwardim[ix2],deg1);

      if(nvr == 3)                                       // c[0] = f[0]*x[2]
      {
         ix2 = idx[2]*deg1;
         cmplx_padded_convolute<<<1,BS>>>
            (forwardre,forwardim,&inputre[ix2],&inputim[ix2],
             crossre,crossim,deg1);
      }
      else
      {
         for(int i=0; i<nvr-3; i++)                  // c[i] = f[i]*b[n-4-i]
         {
            ix1 = i*deg1; ix2 = (nvr-4-i)*deg1;
            cmplx_padded_convolute<<<1,BS>>>
               (&forwardre[ix1],&forwardim[ix1],&backwardre[ix2],
                &backwardim[ix2],&crossre[ix1],&crossim[ix1],deg1);
         }
         ix1 = (nvr-3)*deg1; ix2 = idx[nvr-1]*deg1; // c[n-3] = f[n-3]*x[n-1]
         cmplx_padded_convolute<<<1,BS>>>
            (&forwardre[ix1],&forwardim[ix1],&inputre[ix2],&inputim[ix2],
             &crossre[ix1],&crossim[ix1],deg1);
      }
   }
}

void GPU_dbl_evaldiff
 ( int BS, int dim, int nvr, int deg, int *idx, double *cff,
   double **input, double **output )
{
   const int deg1 = deg+1;            // length of all vectors
   double *input_d;                   // input_d is input on the device
   double *forward_d;                 // forward products on the device
   double *backward_d;                // backward products on the device
   double *cross_d;                   // cross products on the device
   double *cff_d;                     // cff_d is cff on device

   size_t szcff = deg1*sizeof(double);
   size_t szdim = dim*(deg1)*sizeof(double);
   size_t sznvr = nvr*(deg1)*sizeof(double);
   size_t sznvr1 = (nvr-1)*(deg1)*sizeof(double);
   size_t sznvr2 = (nvr-2)*(deg1)*sizeof(double);

   cudaMalloc((void**)&cff_d,szcff);
   cudaMalloc((void**)&input_d,szdim);
   cudaMalloc((void**)&forward_d,sznvr);
   cudaMalloc((void**)&backward_d,sznvr1);
   cudaMalloc((void**)&cross_d,sznvr2);

   double *input_h = new double[dim*(deg1)];
   int ix = 0;
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++) input_h[ix++] = input[i][j];

   cudaMemcpy(cff_d,cff,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(input_d,input_h,szdim,cudaMemcpyHostToDevice);

   if(BS == deg1)
   {
      GPU_dbl_speel
         (BS,nvr,deg,idx,cff_d,input_d,forward_d,backward_d,cross_d);
   }
   double *forward_h = new double[(deg1)*nvr];
   double *backward_h = new double[(deg1)*(nvr-1)];
   double *cross_h = new double[(deg1)*(nvr-2)];
  
   cudaMemcpy(forward_h,forward_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(backward_h,backward_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(cross_h,cross_d,sznvr2,cudaMemcpyDeviceToHost);

   int offset = (nvr-1)*deg1;            // assign value of the monomial
   for(int i=0; i<deg1; i++)
      output[dim][i] = forward_h[offset+i];

   ix = idx[nvr-1];                      // derivative with respect to x[n-1]
   offset = (nvr-2)*deg1;
   for(int i=0; i<deg1; i++) output[ix][i] = forward_h[offset+i];

   ix = idx[0];                          // derivative with respect to x[0]
   offset = (nvr-2)*deg1;
   for(int i=0; i<deg1; i++) output[ix][i] = backward_h[offset+i];

   for(int k=1; k<nvr-1; k++)            // derivative with respect to x[k]
   {
      ix = idx[k]; offset = (k-1)*deg1;
      for(int i=0; i<deg1; i++) output[ix][i] = cross_h[offset+i];
   }
}

void GPU_cmplx_evaldiff
 ( int BS, int dim, int nvr, int deg, int *idx, double *cffre, double *cffim,
   double **inputre, double **inputim, double **outputre, double **outputim )
{
   const int deg1 = deg+1;            // length of all vectors
   double *inputre_d;                 // inputre_d is inputre on the device
   double *inputim_d;                 // inputim_d is inputre on the device
   double *forwardre_d;
   double *forwardim_d;
   double *backwardre_d;
   double *backwardim_d;
   double *crossre_d;
   double *crossim_d;
   double *cffre_d;                   // cffre_d is cffre on the device
   double *cffim_d;                   // cffim_d is cffim on the device

   size_t szdim = dim*(deg1)*sizeof(double);
   size_t sznvr = nvr*(deg1)*sizeof(double);
   size_t sznvr1 = (nvr-1)*(deg1)*sizeof(double);
   size_t sznvr2 = (nvr-2)*(deg1)*sizeof(double);
   size_t szcff = deg1*sizeof(double);

   cudaMalloc((void**)&cffre_d,szcff);
   cudaMalloc((void**)&cffim_d,szcff);
   cudaMalloc((void**)&inputre_d,szdim);
   cudaMalloc((void**)&inputim_d,szdim);
   cudaMalloc((void**)&forwardre_d,sznvr);
   cudaMalloc((void**)&forwardim_d,sznvr);
   cudaMalloc((void**)&backwardre_d,sznvr1);
   cudaMalloc((void**)&backwardim_d,sznvr1);
   cudaMalloc((void**)&crossre_d,sznvr2);
   cudaMalloc((void**)&crossim_d,sznvr2);

   double *inputre_h = new double[dim*(deg1)];
   double *inputim_h = new double[dim*(deg1)];
   int ix = 0;
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++)
      {
         inputre_h[ix] = inputre[i][j]; inputim_h[ix++] = inputim[i][j];
      }

   cudaMemcpy(cffre_d,cffre,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffim_d,cffim,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(inputre_d,inputre_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputim_d,inputim_h,szdim,cudaMemcpyHostToDevice);

   if(BS == deg1)
   {
      GPU_cmplx_speel
         (BS,nvr,deg,idx,cffre_d,cffim_d,inputre_d,inputim_d,forwardre_d,
          forwardim_d,backwardre_d,backwardim_d,crossre_d,crossim_d);
   }
   double *forwardre_h = new double[(deg1)*nvr];
   double *forwardim_h = new double[(deg1)*nvr];
   double *backwardre_h = new double[(deg1)*(nvr-1)];
   double *backwardim_h = new double[(deg1)*(nvr-1)];
   double *crossre_h = new double[(deg1)*(nvr-2)];
   double *crossim_h = new double[(deg1)*(nvr-2)];
  
   cudaMemcpy(forwardre_h,forwardre_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardim_h,forwardim_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardre_h,backwardre_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardim_h,backwardim_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossre_h,crossre_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossim_h,crossim_d,sznvr2,cudaMemcpyDeviceToHost);

   int offset = (nvr-1)*deg1;
   for(int i=0; i<deg1; i++)   // assign value of the monomial
   {
      outputre[dim][i] = forwardre_h[offset+i];
      outputim[dim][i] = forwardim_h[offset+i];
   }
   ix = idx[nvr-1];
   offset = (nvr-2)*deg1;
   for(int i=0; i<deg1; i++)  // derivative with respect to x[n-1]
   {
      outputre[ix][i] = forwardre_h[offset+i];
      outputim[ix][i] = forwardim_h[offset+i];
   }
   ix = idx[0]; 
   offset = (nvr-2)*deg1;
   for(int i=0; i<deg1; i++)   // derivative with respect to x[0]
   {
      outputre[ix][i] = backwardre_h[offset+i];
      outputim[ix][i] = backwardim_h[offset+i];
   }
   for(int k=1; k<nvr-1; k++)  // derivative with respect to x[k]
   {
      ix = idx[k]; offset = (k-1)*deg1;
      for(int i=0; i<deg1; i++)
      {
         outputre[ix][i] = crossre_h[offset+i];
         outputim[ix][i] = crossim_h[offset+i];
      }
   }
}
