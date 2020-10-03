// Defines code of the functions in cmplx_norm_kernels.h,
// to compute the 2-norm and normalize a complex vector in double precision,
// for small and large vectors.

#include <iostream>
#include <cmath>
#include <assert.h>
#include <cstdio>
#include "cmplx_norm_kernels.h"

using namespace std;

__global__ void small_normalize_vector
 ( double* vre, double* vim, int dim, int dimLog2, double* twonorm )
{
   int j = threadIdx.x;
   __shared__ double shvre[d_shmemsize];
   __shared__ double shvim[d_shmemsize];
   __shared__ double prd[d_shmemsize];
   shvre[j] = vre[j];    // reading real parts into shared memory
   shvim[j] = vim[j];    // reading imaginary parts into shared memory
   prd[j] = shvre[j]*shvre[j] + shvim[j]*shvim[j];
   __syncthreads();
   int powTwo = 1;                          // sum reduction
   for(int k=0; k < dimLog2; k++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < dim) prd[j] = prd[j] + prd[j+powTwo];
      powTwo = powTwo*2;
      __syncthreads();
   }
   // thread 0 computes the sqrt of the inner product, others wait
   if(j == 0) prd[0] = sqrt(prd[0]); 
   if(j == 0) *twonorm = prd[0];
   __syncthreads();
   vre[j] = shvre[j]/prd[0];
   vim[j] = shvim[j]/prd[0];
}

__global__ void medium_normalize_vector
 ( double* vre, double* vim, int dim, int rnd, int rndLog2,
   int BS, int BSLog2, double* twonorm )
{
   int j = threadIdx.x;
   int powTwo;
   int vBSind = 0;

   __shared__ double shvre[d_shmemsize];
   __shared__ double shvim[d_shmemsize];
   __shared__ double prd[d_shmemsize];
   __shared__ double sums[maxrounds];

   for(int i=0; i<rnd; i++)
   {
      if(vBSind+j >= dim)       // at last round some threads are excluded
         prd[j] = 0.0;
      else
      {
         shvre[j] = vre[vBSind+j];  // reading of vector into shared memory
         shvim[j] = vim[vBSind+j]; 
         prd[j] = shvre[j]*shvre[j] + shvim[j]*shvim[j];
      }
      __syncthreads();
      powTwo = 1;                          // sum reduction
      for(int k=0; k < BSLog2; k++)
      {
         if((j%(powTwo*2)) == 0)
            if(j+powTwo < BS)
               prd[j] = prd[j] + prd[j+powTwo];
         powTwo = powTwo*2;
         __syncthreads();
      }
      // thread 0 copies the sum of this round in sums[i], the others wait
      if(j == 0) sums[i] = prd[0]; 
      __syncthreads();
      vBSind = vBSind + BS;
   }
   powTwo = 1;                          // sum reduction
   for(int k=0; k < rndLog2; k++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < rnd)
            sums[j] = sums[j] + sums[j+powTwo];
      powTwo = powTwo*2;
      __syncthreads();
   }
   if(j == 0) sums[0] = sqrt(sums[0]);
   if(j == 0) *twonorm = sums[0];
   __syncthreads();
   vBSind = 0;
   for(int i=0; i<rnd; i++)
   {
      if(vBSind+j < dim)
      {
         shvre[j] = vre[vBSind+j];           // read into shared memory
         shvim[j] = vim[vBSind+j];
         vre[vBSind+j] = shvre[j]/sums[0];   // normalize vector
         vim[vBSind+j] = shvim[j]/sums[0];
      }
      __syncthreads();
      vBSind = vBSind + BS;
   }
}

__global__ void large_sum_the_squares
 ( double* vre, double* vim, int dim, double* sums, int BS, int BSLog2 )
{
   const int i = blockIdx.x;
   const int j = threadIdx.x;
   const int k = i*BS + j;

   __shared__ double shvre[d_shmemsize];
   __shared__ double shvim[d_shmemsize];
   __shared__ double prd[d_shmemsize];

   shvre[j] = vre[k];
   shvim[j] = vim[k];
   prd[j] = shvre[j]*shvre[j] + shvim[j]*shvim[j];

   __syncthreads();

   int powTwo = 1;                          // sum reduction
   for(int L=0; L < BSLog2; L++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < BS) prd[j] = prd[j] + prd[j+powTwo];
      powTwo = powTwo*2;

      __syncthreads();
   }
   if(j == 0) sums[i] = prd[0];     // thread 0 writes the sum
}

__global__ void large_normalize_vector
 ( double* vre, double* vim, int dim, double* sums,
   int nbsums, int nbsumsLog2, int BS, double* twonorm )
{
   const int i = blockIdx.x;
   const int j = threadIdx.x;
   const int k = i*BS + j;

   __shared__ double shvre[d_shmemsize];
   __shared__ double shvim[d_shmemsize];

   if(j < nbsums) shvre[j] = sums[j];

   __syncthreads();

   int powTwo = 1;                          // sum reduction
   for(int L=0; L < nbsumsLog2; L++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < nbsums) shvre[j] = shvre[j] + shvre[j+powTwo];
      powTwo = powTwo*2;

      __syncthreads();
   }
   __syncthreads();                      // every thread 0 of all blocks
   if(j == 0) *twonorm = sqrt(shvre[0]); // compute the 2-norm and assigns
   __syncthreads();                      // to the output parameter

   if(k < dim)
   {
      shvre[j] = vre[k];
      shvim[j] = vim[k];
      shvre[j] = shvre[j]/(*twonorm);
      shvim[j] = shvim[j]/(*twonorm);
      vre[k] = shvre[j];
      vim[k] = shvim[j];
   }
}

void GPU_norm
 ( double* vre_h, double* vim_h, int dim, int freq, int BS,
   double* twonorm, int blocked )
{
   int BSLog2 = ceil(log2((double) BS)); // ceil for sum reduction

   double* vre_d;                   // allocate for real parts on device
   double* vim_d;                   // allocate for imaginary parts on device
   size_t size = dim*sizeof(double);
   cudaMalloc((void**)&vre_d,size);
   cudaMalloc((void**)&vim_d,size);
   cudaMemcpy(vre_d,vre_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vim_d,vim_h,size,cudaMemcpyHostToDevice);
   double* twonorm_d;
   cudaMalloc((void**)&twonorm_d,sizeof(double));

   if(dim == BS)
   {
      for(int i=0; i<freq; i++)
         small_normalize_vector<<<1,BS>>>(vre_d,vim_d,dim,BSLog2,twonorm_d);
   }
   else if(blocked == 0)
   {
      int rf = ceil(((double) dim)/BS);
      int rfLog2 = ceil(log2((double) rf));
      for(int i=0; i<freq; i++)
         medium_normalize_vector<<<1,BS>>>
            (vre_d,vim_d,dim,rf,rfLog2,BS,BSLog2,twonorm_d);
   }
   else
   {
      const int nblocks = dim/BS;
      const int nblocksLog2 = ceil(log2((double) nblocks));
      double* sums_d; // sums of squares for each block
      size_t sums_size = nblocks*sizeof(double);
      cudaMalloc((void**)&sums_d,sums_size);
      for(int i=0; i<freq; i++)
      {
         large_sum_the_squares<<<nblocks,BS>>>
            (vre_d,vim_d,dim,sums_d,BS,BSLog2);
         large_normalize_vector<<<nblocks,BS>>>
            (vre_d,vim_d,dim,sums_d,nblocks,nblocksLog2,BS,twonorm_d);
      }
   }
   cudaMemcpy(vre_h,vre_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vim_h,vim_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(twonorm,twonorm_d,sizeof(double),cudaMemcpyDeviceToHost);
}
