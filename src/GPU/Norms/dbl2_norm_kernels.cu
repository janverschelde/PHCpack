// This file contains the definition for the functions in dbl2_norm_kernels.h,
// to compute the 2-norm and normalize a double double precision vector,
// for small, medium, and large vectors.

#include <iostream>
#include <cmath>
#include <assert.h>
#include <cstdio>
#include "double_double_gpufun.cu"
#include "dbl2_norm_kernels.h"

using namespace std;

__global__ void small_normalize_vector
 ( double* vhi, double* vlo, int dim, int dimLog2,
   double* normhi, double* normlo )
{
   int j = threadIdx.x;
   __shared__ double shvhi[dd_shmemsize];
   __shared__ double shvlo[dd_shmemsize];
   __shared__ double prdhi[dd_shmemsize];
   __shared__ double prdlo[dd_shmemsize];
   shvhi[j] = vhi[j];    // reading of vector into shared memory
   shvlo[j] = vlo[j];
   ddg_sqr(shvhi[j],shvlo[j],&prdhi[j],&prdlo[j]);
   __syncthreads();
   int powTwo = 1;                          // sum reduction
   for(int k=0; k < dimLog2; k++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < dim)
            ddg_inc(&prdhi[j],&prdlo[j],prdhi[j+powTwo],prdlo[j+powTwo]);
      powTwo = powTwo*2;
      __syncthreads();
   }
   // thread 0 computes the sqrt of the inner product, others wait
   if(j == 0) ddg_sqrt(prdhi[0],prdlo[0],&prdhi[0],&prdlo[0]); 
   if(j == 0)
   {
     *normhi = prdhi[0];
     *normlo = prdlo[0];
   }
   __syncthreads();
   ddg_div(shvhi[j],shvlo[j],prdhi[0],prdlo[0],&vhi[j],&vlo[j]);
}

__global__ void medium_normalize_vector
 ( double* vhi, double* vlo, int dim, int rnd, int rndLog2,
   int BS, int BSLog2, double* normhi, double* normlo )
{
   int j = threadIdx.x;
   int powTwo;
   int vBSind = 0;

   __shared__ double shvhi[dd_shmemsize];
   __shared__ double shvlo[dd_shmemsize];
   __shared__ double prdhi[dd_shmemsize];
   __shared__ double prdlo[dd_shmemsize];
   __shared__ double sumshi[maxrounds];
   __shared__ double sumslo[maxrounds];
/*
   for(int i=0; i<rnd; i++)
   {
      if(vBSind+j >= dim)       // at last round some threads are excluded
         prd[j] = 0.0;
      else
      {
         shv[j] = v[vBSind+j];  // reading of vector into shared memory
         prd[j] = shv[j]*shv[j];
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
         shv[j] = v[vBSind+j];           // read into shared memory
         v[vBSind+j] = shv[j]/sums[0];   // normalize vector
      }
      __syncthreads();
      vBSind = vBSind + BS;
   }
*/
}

__global__ void large_sum_the_squares
 ( double* vhi, double* vlo, int dim, double* sumshi, double* sumslo,
   int BS, int BSLog2 )
{
   const int i = blockIdx.x;
   const int j = threadIdx.x;
   const int k = i*BS + j;

   __shared__ double shvhi[dd_shmemsize];
   __shared__ double shvlo[dd_shmemsize];
   __shared__ double prdhi[dd_shmemsize];
   __shared__ double prdlo[dd_shmemsize];
/*
   shv[j] = v[k];
   prd[j] = shv[j]*shv[j];

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
*/
}

__global__ void large_normalize_vector
 ( double* vhi, double* vlo, int dim, double* sumshi, double* sumslo,
   int nbsums, int nbsumsLog2, int BS, double* normhi, double* normlo )
{
   const int i = blockIdx.x;
   const int j = threadIdx.x;
   const int k = i*BS + j;
/*
   __shared__ double shv[dd_shmemsize];

   if(j < nbsums) shv[j] = sums[j];

   __syncthreads();

   int powTwo = 1;                          // sum reduction
   for(int L=0; L < nbsumsLog2; L++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < nbsums) shv[j] = shv[j] + shv[j+powTwo];
      powTwo = powTwo*2;

      __syncthreads();
   }
   __syncthreads();                    // every thread 0 of all blocks
   if(j == 0) *twonorm = sqrt(shv[0]); // compute the 2-norm and assigns
   __syncthreads();                    // to the output parameter

   if(k < dim)
   {
      shv[j] = v[k];
      shv[j] = shv[j]/(*twonorm);
      v[k] = shv[j];
   }
*/
}

void GPU_norm
 ( double* vhi_h, double* vlo_h, int dim, int freq, int BS,
   double* normhi, double* normlo, int blocked )
{
   int BSLog2 = ceil(log2((double) BS)); // ceil for sum reduction

   double* vhi_d;                   // allocate for vector on device
   double* vlo_d;
   size_t size = dim*sizeof(double);
   cudaMalloc((void**)&vhi_d,size);
   cudaMalloc((void**)&vlo_d,size);
   cudaMemcpy(vhi_d,vhi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vlo_d,vlo_h,size,cudaMemcpyHostToDevice);
   double* normhi_d;
   double* normlo_d;
   cudaMalloc((void**)&normhi_d,sizeof(double));
   cudaMalloc((void**)&normlo_d,sizeof(double));

   if(dim == BS)
   {
      for(int i=0; i<freq; i++)
         small_normalize_vector<<<1,BS>>>
            (vhi_d,vlo_d,dim,BSLog2,normhi_d,normlo_d);
   }
   else if(blocked == 0)
   {
      const int rf = ceil(((double) dim)/BS);
      const int rfLog2 = ceil(log2((double) rf));
      for(int i=0; i<freq; i++)
         medium_normalize_vector<<<1,BS>>>
            (vhi_d,vlo_d,dim,rf,rfLog2,BS,BSLog2,normhi_d,normlo_d);
   }
   else
   {
      const int nblocks = dim/BS;
      const int nblocksLog2 = ceil(log2((double) nblocks));
      double* sumshi_d; // sums of squares for each block
      double* sumslo_d; // low parts of sums of squares
      size_t sums_size = nblocks*sizeof(double);
      cudaMalloc((void**)&sumshi_d,sums_size);
      cudaMalloc((void**)&sumslo_d,sums_size);
      for(int i=0; i<freq; i++)
      {
         large_sum_the_squares<<<nblocks,BS>>>
            (vhi_d,vlo_d,dim,sumshi_d,sumslo_d,BS,BSLog2);
         large_normalize_vector<<<nblocks,BS>>>
            (vhi_d,vlo_d,dim,sumshi_d,sumslo_d,nblocks,nblocksLog2,BS,
             normhi_d,normlo_d);
      }
   }
   cudaMemcpy(vhi_h,vhi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vlo_h,vlo_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(normhi,normhi_d,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(normlo,normlo_d,sizeof(double),cudaMemcpyDeviceToHost);
}
