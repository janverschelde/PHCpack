// The file dbl_monomials_kernels.cu defines the kernels with prototypes
// in dbl_monomials_kernels.h.

#include <iostream>

#include "dbl_convolutions_kernels.h"
#include "dbl_monomials_kernels.h"

using namespace std;

__device__ void dbl_convolute
 ( double *x, double *y, double *z, int dim, int k )
{
   z[k] = x[0]*y[k];

   for(int i=1; i<=k; i++) z[k] = z[k] + x[i]*y[k-i];
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
   __syncthreads();
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
