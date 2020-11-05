// The file dbl_monomials_kernels.cu defines the kernels with prototypes
// in dbl_monomials_kernels.h.

#include "dbl_convolutions_kernels.h"
#include "dbl_monomials_kernels.h"

using namespace std;

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
 ( int nvr, int deg, int *idx, double *input,
   double *forward, double *backward, double *cross )
{
   const int k = threadIdx.x;
   const int deg1 = deg+1;
   int ix1,ix2;

   __shared__ double xv[d_shmemsize];
   __shared__ double yv[d_shmemsize];
   __shared__ double zv[d_shmemsize];

   ix1 = idx[0]*deg1+k; xv[k] = input[ix1]; 
   ix2 = idx[1]*deg1+k; yv[k] = input[ix2];
   __syncthreads(); dbl_convolute(xv,yv,zv,deg1,k);
   forward[k] = zv[k];

   for(int i=2; i<nvr; i++)
   {
      xv[k] = zv[k];
      ix2 = idx[i]*deg1+k; yv[k] = input[ix2];
      __syncthreads(); dbl_convolute(xv,yv,zv,deg1,k);
      forward[(i-1)*deg1+k] = zv[k];
   }
   if(nvr > 2)
   {
      ix1 = idx[nvr-1]*deg1+k; xv[k] = input[ix1];
      ix2 = idx[nvr-2]*deg1+k; yv[k] = input[ix2];
      __syncthreads(); dbl_convolute(xv,yv,zv,deg1,k);
      backward[k] = zv[k];
      for(int i=1; i<nvr-2; i++)
      {
         xv[k] = zv[k]; 
         ix2 = idx[nvr-2-i]*deg1+k; yv[k] = input[ix2];
         __syncthreads(); dbl_convolute(xv,yv,zv,deg1,k);
         backward[i*deg1+k] = zv[k];
      }
      if(nvr == 3)
      {
         ix1 = idx[0]*deg1+k; xv[k] = input[ix1];
         ix2 = idx[2]*deg1+k; yv[k] = input[ix2];
         __syncthreads(); dbl_convolute(xv,yv,zv,deg1,k);
         cross[k] = zv[k];
      }
      else
      {
         ix1 = idx[0]*deg1+k;  xv[k] = input[ix1];
         ix2 = (nvr-3)*deg1+k; yv[k] = backward[ix2];
         __syncthreads(); dbl_convolute(xv,yv,zv,deg1,k);
         cross[k] = zv[k];
         for(int i=1; i<nvr-3; i++)
         {
            ix1 = (i-1)*deg1+k;     xv[k] = forward[ix1];
            ix2 = (nvr-3-i)*deg1+k; yv[k] = backward[ix2];
            __syncthreads(); dbl_convolute(xv,yv,zv,deg1,k);
            cross[i*deg1+k] = zv[k];
         }
         ix1 = (nvr-4)*deg1+k;    xv[k] = forward[ix1];
         ix2 = idx[nvr-1]*deg1+k; yv[k] = input[ix2];
         __syncthreads(); dbl_convolute(xv,yv,zv,deg1,k);
         cross[(nvr-3)*deg1+k] = zv[k];
      }
   }
}

__global__ void GPU_cmplx_speel
 ( int nvr, int deg, int *idx, double *inputre, double *inputim,
   double *forwardre, double *forwardim, double *backwardre,
   double *backwardim, double *crossre, double *crossim )
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

   ix1 = idx[0]*deg1+k; xvre[k] = inputre[ix1]; xvim[k] = inputim[ix1];
   ix2 = idx[1]*deg1+k; yvre[k] = inputre[ix2]; yvim[k] = inputim[ix2];
   __syncthreads(); cmplx_convolute(xvre,xvim,yvre,yvim,zvre,zvim,deg1,k);
   forwardre[k] = zvre[k]; forwardim[k] = zvim[k];

   for(int i=2; i<nvr; i++)
   {
      xvre[k] = zvre[k]; xvim[k] = zvim[k];
      ix2 = idx[i]*deg1+k; yvre[k] = inputre[ix2]; yvim[k] = inputim[ix2];
      __syncthreads(); cmplx_convolute(xvre,xvim,yvre,yvim,zvre,zvim,deg1,k);
      ix1 = (i-1)*deg1+k;
      forwardre[ix1] = zvre[k]; forwardim[ix1] = zvim[k];
   }
   if(nvr > 2)
   {
      ix1 = idx[nvr-1]*deg1+k; xvre[k] = inputre[ix1]; xvim[k] = inputim[ix1];
      ix2 = idx[nvr-2]*deg1+k; yvre[k] = inputre[ix2]; yvim[k] = inputim[ix2];
      __syncthreads(); cmplx_convolute(xvre,xvim,yvre,yvim,zvre,zvim,deg1,k);
      backwardre[k] = zvre[k]; backwardim[k] = zvim[k];
      for(int i=1; i<nvr-2; i++)
      {
         xvre[k] = zvre[k]; xvim[k] = zvim[k];
         ix2 = idx[nvr-2-i]*deg1+k;
         yvre[k] = inputre[ix2]; yvim[k] = inputim[ix2];
         __syncthreads();
         cmplx_convolute(xvre,xvim,yvre,yvim,zvre,zvim,deg1,k);
         ix1 = i*deg1+k;
         backwardre[ix1] = zvre[k]; backwardim[ix1] = zvim[k];
      }
      if(nvr == 3)
      {
         ix1 = idx[0]*deg1+k; xvre[k] = inputre[ix1]; xvim[k] = inputim[ix1];
         ix2 = idx[2]*deg1+k; yvre[k] = inputre[ix2]; yvim[k] = inputim[ix2];
         __syncthreads();
         cmplx_convolute(xvre,xvim,yvre,yvim,zvre,zvim,deg1,k);
         crossre[k] = zvre[k]; crossim[k] = zvim[k];
      }
      else
      {
         ix1 = idx[0]*deg1+k;  xvre[k] = inputre[ix1]; xvim[k] = inputim[ix1];
         ix2 = (nvr-3)*deg1+k;
         yvre[k] = backwardre[ix2]; yvim[k] = backwardim[ix2];
         __syncthreads();
         cmplx_convolute(xvre,xvim,yvre,yvim,zvre,zvim,deg1,k);
         crossre[k] = zvre[k]; crossim[k] = zvim[k];
         for(int i=1; i<nvr-3; i++)
         {
            ix1 = (i-1)*deg1+k;   
            xvre[k] = forwardre[ix1]; xvim[k] = forwardim[ix1];
            ix2 = (nvr-3-i)*deg1+k;
            yvre[k] = backwardre[ix2]; yvim[k] = backwardim[ix2];
            __syncthreads();
            cmplx_convolute(xvre,xvim,yvre,yvim,zvre,zvim,deg1,k);
            ix1 = i*deg1+k;
            crossre[ix1] = zvre[k]; crossim[ix1] = zvim[k];
         }
         ix1 = (nvr-4)*deg1+k;    
         xvre[k] = forwardre[ix1]; xvim[k] = forwardim[ix1];
         ix2 = idx[nvr-1]*deg1+k;
         yvre[k] = inputre[ix2]; yvim[k] = inputim[ix2];
         __syncthreads();
         cmplx_convolute(xvre,xvim,yvre,yvim,zvre,zvim,deg1,k);
         ix1 = (nvr-3)*deg1+k;
         crossre[ix1] = zvre[k]; crossim[ix1] = zvim[k];
      }
   }
}

void GPU_dbl_evaldiff
 ( int BS, int dim, int nvr, int deg, int *idx, double *cff,
   double **input, double **output )
{
   const int deg1 = deg+1;            // length of all vectors
   double *input_d;                   // input_d is input on the device
   double *forward_d;
   double *backward_d;
   double *cross_d;
   int *idx_d;                        // idx_d is idx on device

   size_t szdim = dim*(deg1)*sizeof(double);
   size_t sznvr1 = (nvr-1)*(deg1)*sizeof(double);
   size_t sznvr2 = (nvr-2)*(deg1)*sizeof(double);
   size_t szidx = nvr*sizeof(int);

   cudaMalloc((void**)&idx_d,szidx);
   cudaMalloc((void**)&input_d,szdim);
   cudaMalloc((void**)&forward_d,sznvr1);
   cudaMalloc((void**)&backward_d,sznvr2);
   cudaMalloc((void**)&cross_d,sznvr2);

   double *input_h = new double[dim*(deg1)];
   int ix = 0;
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++) input_h[ix++] = input[i][j];

   cudaMemcpy(idx_d,idx,szidx,cudaMemcpyHostToDevice);
   cudaMemcpy(input_d,input_h,szdim,cudaMemcpyHostToDevice);

   if(BS = deg1)
   {
      GPU_dbl_speel<<<1,BS>>>
         (nvr,deg,idx_d,input_d,forward_d,backward_d,cross_d);
   }
   double *forward_h = new double[(deg1)*(nvr-1)];
   double *backward_h = new double[(deg1)*(nvr-2)];
   double *cross_h = new double[(deg1)*(nvr-2)];
  
   cudaMemcpy(forward_h,forward_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backward_h,backward_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(cross_h,cross_d,sznvr2,cudaMemcpyDeviceToHost);

   int offset = (nvr-2)*deg1;
   for(int i=0; i<deg1; i++) output[dim][i] = forward_h[offset+i];
   ix = idx[nvr-1]; offset = (nvr-3)*deg1;
   for(int i=0; i<deg1; i++) output[ix][i] = forward_h[offset+i];
   ix = idx[0]; 
   for(int i=0; i<deg1; i++) output[ix][i] = backward_h[offset+i];
   for(int k=1; k<nvr-1; k++)
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
   int *idx_d;                        // idx_d is idx on device

   size_t szdim = dim*(deg1)*sizeof(double);
   size_t sznvr1 = (nvr-1)*(deg1)*sizeof(double);
   size_t sznvr2 = (nvr-2)*(deg1)*sizeof(double);
   size_t szidx = nvr*sizeof(int);

   cudaMalloc((void**)&idx_d,szidx);
   cudaMalloc((void**)&inputre_d,szdim);
   cudaMalloc((void**)&inputim_d,szdim);
   cudaMalloc((void**)&forwardre_d,sznvr1);
   cudaMalloc((void**)&forwardim_d,sznvr1);
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
   cudaMemcpy(inputre_d,inputre_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputim_d,inputim_h,szdim,cudaMemcpyHostToDevice);

   if(BS = deg1)
   {
      GPU_cmplx_speel<<<1,BS>>>
         (nvr,deg,idx_d,inputre_d,inputim_d,forwardre_d,forwardim_d,
          backwardre_d,backwardim_d,crossre_d,crossim_d);
   }
   double *forwardre_h = new double[(deg1)*(nvr-1)];
   double *forwardim_h = new double[(deg1)*(nvr-1)];
   double *backwardre_h = new double[(deg1)*(nvr-2)];
   double *backwardim_h = new double[(deg1)*(nvr-2)];
   double *crossre_h = new double[(deg1)*(nvr-2)];
   double *crossim_h = new double[(deg1)*(nvr-2)];
  
   cudaMemcpy(forwardre_h,forwardre_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardim_h,forwardim_d,sznvr1,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardre_h,backwardre_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardim_h,backwardim_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossre_h,crossre_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossim_h,crossim_d,sznvr2,cudaMemcpyDeviceToHost);

   int offset = (nvr-2)*deg1;
   for(int i=0; i<deg1; i++)
   {
      outputre[dim][i] = forwardre_h[offset+i];
      outputim[dim][i] = forwardim_h[offset+i];
   }
   ix = idx[nvr-1]; offset = (nvr-3)*deg1;
   for(int i=0; i<deg1; i++)
   {
      outputre[ix][i] = forwardre_h[offset+i];
      outputim[ix][i] = forwardim_h[offset+i];
   }
   ix = idx[0]; 
   for(int i=0; i<deg1; i++) 
   {
      outputre[ix][i] = backwardre_h[offset+i];
      outputim[ix][i] = backwardim_h[offset+i];
   }
   for(int k=1; k<nvr-1; k++)
   {
      ix = idx[k]; offset = (k-1)*deg1;
      for(int i=0; i<deg1; i++)
      {
         outputre[ix][i] = crossre_h[offset+i];
         outputim[ix][i] = crossim_h[offset+i];
      }
   }
}
