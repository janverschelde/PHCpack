// Defines code of the functions in cmplx2_norm_kernels.h,
// to compute the 2-norm and normalize a complex vector,
// in double double precision,
// for vectors of small, medium, and large size.

#include <iostream>
#include <cmath>
#include <assert.h>
#include <cstdio>
#include "double_double_gpufun.cu"
#include "cmplx2_norm_kernels.h"

using namespace std;

__global__ void small_normalize_vector
 ( double* vrehi, double* vrelo, double* vimhi, double* vimlo,
   int dim, int dimLog2, double* normhi, double* normlo )
{
   int j = threadIdx.x;

   __shared__ double shvrehi[dd_shmemsize];
   __shared__ double shvrelo[dd_shmemsize];
   __shared__ double shvimhi[dd_shmemsize];
   __shared__ double shvimlo[dd_shmemsize];
   __shared__ double prdhi[dd_shmemsize];
   __shared__ double prdlo[dd_shmemsize];
   __shared__ double sumhi[dd_shmemsize];
   __shared__ double sumlo[dd_shmemsize];

   shvrehi[j] = vrehi[j]; // reading real parts into shared memory
   shvrelo[j] = vrelo[j];
   shvimhi[j] = vimhi[j]; // reading imaginary parts into shared memory
   shvimlo[j] = vimlo[j];

   ddg_sqr(shvrehi[j],shvrelo[j],&sumhi[j],&sumlo[j]);
   ddg_sqr(shvimhi[j],shvimlo[j],&prdhi[j],&prdlo[j]);
   ddg_inc(&sumhi[j],&sumlo[j],prdhi[j],prdlo[j]);

   __syncthreads();
   int powTwo = 1;                          // sum reduction
   for(int k=0; k < dimLog2; k++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < dim)
            ddg_inc(&sumhi[j],&sumlo[j],sumhi[j+powTwo],sumlo[j+powTwo]);
      powTwo = powTwo*2;
      __syncthreads();
   }
   // thread 0 computes the sqrt of the inner product, others wait
   if(j == 0) ddg_sqrt(sumhi[0],sumlo[0],&sumhi[0],&sumlo[0]); 
   if(j == 0)
   {
      *normhi = sumhi[0];
      *normlo = sumlo[0];
   }
   __syncthreads();
   ddg_div(shvrehi[j],shvrelo[j],sumhi[0],sumlo[0],&vrehi[j],&vrelo[j]);
   ddg_div(shvimhi[j],shvimlo[j],sumhi[0],sumlo[0],&vimhi[j],&vimlo[j]);
}

void GPU_norm
 ( double* vrehi_h, double* vrelo_h, double* vimhi_h, double* vimlo_h,
   int dim, int freq, int BS, double* normhi, double* normlo, int blocked )
{
   int BSLog2 = ceil(log2((double) BS)); // ceil for sum reduction

   double* vrehi_d;                      // high real parts on device
   double* vrelo_d;                      // low real parts on device
   double* vimhi_d;                      // high imaginary parts on device
   double* vimlo_d;                      // low imaginary parts on device
   size_t size = dim*sizeof(double);
   cudaMalloc((void**)&vrehi_d,size);
   cudaMalloc((void**)&vrelo_d,size);
   cudaMalloc((void**)&vimhi_d,size);
   cudaMalloc((void**)&vimlo_d,size);
   cudaMemcpy(vrehi_d,vrehi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vrelo_d,vrelo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vimhi_d,vimhi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vimlo_d,vimlo_h,size,cudaMemcpyHostToDevice);
   double* normhi_d;
   double* normlo_d;
   cudaMalloc((void**)&normhi_d,sizeof(double));
   cudaMalloc((void**)&normlo_d,sizeof(double));

   if(dim == BS)
   {
      for(int i=0; i<freq; i++)
         small_normalize_vector<<<1,BS>>>
            (vrehi_d,vrelo_d,vimhi_d,vimlo_d,dim,BSLog2,normhi_d,normlo_d);
   }
   cudaMemcpy(vrehi_h,vrehi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vrelo_h,vrelo_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vimhi_h,vimhi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vimlo_h,vimlo_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(normhi,normhi_d,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(normlo,normlo_d,sizeof(double),cudaMemcpyDeviceToHost);
}
