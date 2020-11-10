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

#include "dbl_convolutions_kernels.h"
#include "dbl_monomials_kernels.h"

__device__ void dbl_convolute
 ( double *x, double *y, double *z, int dim, int k )
{
   z[k] = x[0]*y[k];

   for(int i=1; i<=k; i++) z[k] = z[k] + x[i]*y[k-i];
}

__device__ void cmplx_convolute
 ( double *xre, double *xim, double *yre, double *yim,
   double *zre, double *zim, int dim, int k )
{
   double xr,xi,yr,yi,zr,zi;

   xr = xre[0]; xi = xim[0];    // z[k] = x[0]*y[k]
   yr = yre[k]; yi = yim[k];
   zr = xr*yr - xi*yi;
   zi = xr*yi + xi*yr;
   zre[k] = zr; zim[k] = zi;

   for(int i=1; i<=k; i++) // z[k] = z[k] + x[i]*y[k-i]
   {
      xr = xre[i];   xi = xim[i];
      yr = yre[k-i]; yi = yim[k-i];
      zr = xr*yr - xi*yi;
      zi = xr*yi + xi*yr;
      zre[k] += zr;
      zim[k] += zi;
   }
}

__global__ void GPU_dbl_speel
 ( int nvr, int deg, int *idx, double *cff, double *input,
   double *forward, double *backward, double *cross )
{
   const int k = threadIdx.x;
   const int deg1 = deg+1;
   int ix1,ix2;

   __shared__ double xv[d_shmemsize];
   __shared__ double yv[d_shmemsize];
   __shared__ double zv[d_shmemsize];
  
   xv[k] = cff[k]; ix1 = idx[0]*deg1+k; yv[k] = input[ix1]; 
   __syncthreads(); dbl_convolute(xv,yv,zv,deg1,k);
   forward[k] = zv[k];                               // f[0] = cff*x[0]

   for(int i=1; i<nvr; i++)
   {
      xv[k] = zv[k]; ix2 = idx[i]*deg1+k; yv[k] = input[ix2];
      __syncthreads(); dbl_convolute(xv,yv,zv,deg1,k);
      forward[i*deg1+k] = zv[k];                    // f[i] = f[i-1]*x[i]
   }
   if(nvr > 2)
   {
      ix1 = idx[nvr-1]*deg1+k; xv[k] = input[ix1];
      ix2 = idx[nvr-2]*deg1+k; yv[k] = input[ix2];
      __syncthreads(); dbl_convolute(xv,yv,zv,deg1,k);
      backward[k] = zv[k];                           // b[0] = x[n-1]*x[n-2]
      for(int i=1; i<nvr-2; i++)
      {
         xv[k] = zv[k]; ix2 = idx[nvr-2-i]*deg1+k; yv[k] = input[ix2];
         __syncthreads(); dbl_convolute(xv,yv,zv,deg1,k);
         backward[i*deg1+k] = zv[k];                // b[i] = b[i-1]*x[n-2-i]
      }
      xv[k] = zv[k]; yv[k] = cff[k];
      __syncthreads(); dbl_convolute(xv,yv,zv,deg1,k);
      ix2 = (nvr-3)*deg1+k;
      backward[ix2] = zv[k];                         // b[n-3] = b[n-3]*cff

      if(nvr == 3)
      {
         xv[k] = forward[k]; ix2 = idx[2]*deg1+k; yv[k] = input[ix2];
         __syncthreads(); dbl_convolute(xv,yv,zv,deg1,k);
         cross[k] = zv[k];                          // c[0] = f[0]*x[2]
      }
      else
      {
         for(int i=0; i<nvr-3; i++)
         {
            ix1 = i*deg1+k;         xv[k] = forward[ix1];
            ix2 = (nvr-4-i)*deg1+k; yv[k] = backward[ix2];
            __syncthreads(); dbl_convolute(xv,yv,zv,deg1,k);
            cross[i*deg1+k] = zv[k];                // c[i] = f[i]*b[n-4-i]
         }
         ix1 = (nvr-3)*deg1+k;    xv[k] = forward[ix1];
         ix2 = idx[nvr-1]*deg1+k; yv[k] = input[ix2];
         __syncthreads(); dbl_convolute(xv,yv,zv,deg1,k);
         cross[(nvr-3)*deg1+k] = zv[k];             // c[n-3] = f[n-3]*x[n-1]
      }
   }
}

__global__ void GPU_cmplx_speel
 ( int nvr, int deg, int *idx, double *cffre, double *cffim,
   double *inputre, double *inputim, double *forwardre, double *forwardim,
   double *backwardre, double *backwardim, double *crossre, double *crossim )
{
   const int k = threadIdx.x;
   const int deg1 = deg+1;
   int ix1,ix2;

   __shared__ double xvre[d_shmemsize];
   __shared__ double xvim[d_shmemsize];
   __shared__ double yvre[d_shmemsize];
   __shared__ double yvim[d_shmemsize];
   __shared__ double zvre[d_shmemsize];
   __shared__ double zvim[d_shmemsize];

   xvre[k] = cffre[k]; xvim[k] = cffim[k];
   ix1 = idx[0]*deg1+k; yvre[k] = inputre[ix1]; yvim[k] = inputim[ix1];
   __syncthreads(); cmplx_convolute(xvre,xvim,yvre,yvim,zvre,zvim,deg1,k);
   forwardre[k] = zvre[k];
   forwardim[k] = zvim[k];                            // f[0] = cff*x[0]

   for(int i=1; i<nvr; i++)
   {
      xvre[k] = zvre[k]; xvim[k] = zvim[k];
      ix2 = idx[i]*deg1+k; yvre[k] = inputre[ix2]; yvim[k] = inputim[ix2];
      __syncthreads(); cmplx_convolute(xvre,xvim,yvre,yvim,zvre,zvim,deg1,k);
      ix1 = i*deg1+k;
      forwardre[ix1] = zvre[k];
      forwardim[ix1] = zvim[k];                       // f[i] = f[i-i]*x[i]
   }
   if(nvr > 2)
   {
      ix1 = idx[nvr-1]*deg1+k; xvre[k] = inputre[ix1]; xvim[k] = inputim[ix1];
      ix2 = idx[nvr-2]*deg1+k; yvre[k] = inputre[ix2]; yvim[k] = inputim[ix2];
      __syncthreads(); cmplx_convolute(xvre,xvim,yvre,yvim,zvre,zvim,deg1,k);
      backwardre[k] = zvre[k];
      backwardim[k] = zvim[k];                        // b[0] = x[n-1]*x[n-2]
      for(int i=1; i<nvr-2; i++)
      {
         xvre[k] = zvre[k]; xvim[k] = zvim[k];
         ix2 = idx[nvr-2-i]*deg1+k;
         yvre[k] = inputre[ix2]; yvim[k] = inputim[ix2];
         __syncthreads();
         cmplx_convolute(xvre,xvim,yvre,yvim,zvre,zvim,deg1,k);
         ix1 = i*deg1+k;
         backwardre[ix1] = zvre[k];
         backwardim[ix1] = zvim[k];                // b[i] = b[i]*x[n-2-i]
      }
      xvre[k] = zvre[k];  xvim[k] = zvim[k];
      yvre[k] = cffre[k]; yvim[k] = cffim[k];
      __syncthreads();
      cmplx_convolute(xvre,xvim,yvre,yvim,zvre,zvim,deg1,k);
      ix1 = (nvr-3)*deg1+k;
      backwardre[ix1] = zvre[k];
      backwardim[ix1] = zvim[k];                   // b[n-3] = b[n-3]*cff

      if(nvr == 3)
      {
         xvre[k] = forwardre[k]; xvim[k] = forwardim[k];
         ix2 = idx[2]*deg1+k; yvre[k] = inputre[ix2]; yvim[k] = inputim[ix2];
         __syncthreads();
         cmplx_convolute(xvre,xvim,yvre,yvim,zvre,zvim,deg1,k);
         crossre[k] = zvre[k];
         crossim[k] = zvim[k];                     // c[0] = f[0]*x[2]
      }
      else
      {
         for(int i=0; i<nvr-3; i++)
         {
            ix1 = i*deg1+k;   
            xvre[k] = forwardre[ix1]; xvim[k] = forwardim[ix1];
            ix2 = (nvr-4-i)*deg1+k;
            yvre[k] = backwardre[ix2]; yvim[k] = backwardim[ix2];
            __syncthreads();
            cmplx_convolute(xvre,xvim,yvre,yvim,zvre,zvim,deg1,k);
            ix1 = i*deg1+k;
            crossre[ix1] = zvre[k];
            crossim[ix1] = zvim[k];                // c[i] = f[i]*b[n-4-i]
         }
         ix1 = (nvr-3)*deg1+k;
         xvre[k] = forwardre[ix1]; xvim[k] = forwardim[ix1];
         ix2 = idx[nvr-1]*deg1+k;
         yvre[k] = inputre[ix2];   yvim[k] = inputim[ix2];
         __syncthreads();
         cmplx_convolute(xvre,xvim,yvre,yvim,zvre,zvim,deg1,k);
         ix1 = (nvr-3)*deg1+k;
         crossre[ix1] = zvre[k];
         crossim[ix1] = zvim[k];                  // c[n-3] = f[n-3]*x[n-1]
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
   int *idx_d;                        // idx_d is idx on device

   size_t szcff = deg1*sizeof(double);
   size_t szdim = dim*(deg1)*sizeof(double);
   size_t sznvr = nvr*(deg1)*sizeof(double);
   size_t sznvr2 = (nvr-2)*(deg1)*sizeof(double);
   size_t szidx = nvr*sizeof(int);

   cudaMalloc((void**)&idx_d,szidx);
   cudaMalloc((void**)&cff_d,szcff);
   cudaMalloc((void**)&input_d,szdim);
   cudaMalloc((void**)&forward_d,sznvr);
   cudaMalloc((void**)&backward_d,sznvr2);
   cudaMalloc((void**)&cross_d,sznvr2);

   double *input_h = new double[dim*(deg1)];
   int ix = 0;
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++) input_h[ix++] = input[i][j];

   cudaMemcpy(idx_d,idx,szidx,cudaMemcpyHostToDevice);
   cudaMemcpy(cff_d,cff,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(input_d,input_h,szdim,cudaMemcpyHostToDevice);

   if(BS == deg1)
   {
      GPU_dbl_speel<<<1,BS>>>
         (nvr,deg,idx_d,cff_d,input_d,forward_d,backward_d,cross_d);
   }
   double *forward_h = new double[(deg1)*nvr];
   double *backward_h = new double[(deg1)*(nvr-2)];
   double *cross_h = new double[(deg1)*(nvr-2)];
  
   cudaMemcpy(forward_h,forward_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(backward_h,backward_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(cross_h,cross_d,sznvr2,cudaMemcpyDeviceToHost);

   int offset = (nvr-1)*deg1;            // assign value of the monomial
   for(int i=0; i<deg1; i++)
      output[dim][i] = forward_h[offset+i];

   ix = idx[nvr-1];                      // derivative with respect to x[n-1]
   offset = (nvr-2)*deg1;
   for(int i=0; i<deg1; i++) output[ix][i] = forward_h[offset+i];

   ix = idx[0];                          // derivative with respect to x[0]
   offset = (nvr-3)*deg1;
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
   int *idx_d;                        // idx_d is idx on the device

   size_t szdim = dim*(deg1)*sizeof(double);
   size_t sznvr = nvr*(deg1)*sizeof(double);
   size_t sznvr2 = (nvr-2)*(deg1)*sizeof(double);
   size_t szidx = nvr*sizeof(int);
   size_t szcff = deg1*sizeof(double);

   cudaMalloc((void**)&idx_d,szidx);
   cudaMalloc((void**)&cffre_d,szcff);
   cudaMalloc((void**)&cffim_d,szcff);
   cudaMalloc((void**)&inputre_d,szdim);
   cudaMalloc((void**)&inputim_d,szdim);
   cudaMalloc((void**)&forwardre_d,sznvr);
   cudaMalloc((void**)&forwardim_d,sznvr);
   cudaMalloc((void**)&backwardre_d,sznvr2);
   cudaMalloc((void**)&backwardim_d,sznvr2);
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

   cudaMemcpy(idx_d,idx,szidx,cudaMemcpyHostToDevice);
   cudaMemcpy(cffre_d,cffre,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffim_d,cffim,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(inputre_d,inputre_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputim_d,inputim_h,szdim,cudaMemcpyHostToDevice);

   if(BS == deg1)
   {
      GPU_cmplx_speel<<<1,BS>>>
         (nvr,deg,idx_d,cffre_d,cffim_d,inputre_d,inputim_d,forwardre_d,
          forwardim_d,backwardre_d,backwardim_d,crossre_d,crossim_d);
   }
   double *forwardre_h = new double[(deg1)*nvr];
   double *forwardim_h = new double[(deg1)*nvr];
   double *backwardre_h = new double[(deg1)*(nvr-2)];
   double *backwardim_h = new double[(deg1)*(nvr-2)];
   double *crossre_h = new double[(deg1)*(nvr-2)];
   double *crossim_h = new double[(deg1)*(nvr-2)];
  
   cudaMemcpy(forwardre_h,forwardre_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardim_h,forwardim_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardre_h,backwardre_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardim_h,backwardim_d,sznvr2,cudaMemcpyDeviceToHost);
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
   offset = (nvr-3)*deg1;
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
