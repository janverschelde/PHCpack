// This file contains the definition for the functions in dbl_norm_kernels.h,
// to compute the 2-norm and normalize a vector of double precision numbers,
// for small and large vectors.

#include <iostream>
#include <cmath>
#include <assert.h>
#include <cstdio>

#define d_shmemsize 512
#define maxrounds 32

using namespace std;

__global__ void small_normalize_vector
 ( double* v, int dim, int dimLog2, double* twonorm )
{
   int j = threadIdx.x;
   __shared__ double shv[d_shmemsize];
   __shared__ double prd[d_shmemsize];
   shv[j] = v[j];    // reading of vector into shared memory
   prd[j] = shv[j]*shv[j] + shv[j]*shv[j];
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
   v[j] = shv[j]/prd[0];
}

__global__ void large_normalize_vector
 ( double* v, int dim, int rnd, int rndLog2, int BS, int BSLog2,
   double* twonorm )
{
   int j = threadIdx.x;
   int powTwo;
   int vBSind = 0;

   __shared__ double shv[d_shmemsize];
   __shared__ double prd[d_shmemsize];
   __shared__ double sums[maxrounds];

   for(int i=0; i<rnd; i++)
   {
      if(vBSind+j >= dim)       // at last round some threads are excluded
         prd[j] = 0.0;
      else
      {
         shv[j] = v[vBSind+j];  // reading of vector into shared memory
         prd[j] = shv[j]*shv[j] + shv[j]*shv[j];
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
}

void GPU_norm
 ( double* v_h, int dim, int freq, int BS, double* twonorm )
{
   int BSLog2 = ceil(log2((double) BS)); // ceil for sum reduction

   double* v_d;                   // allocate for vector on device
   size_t size = dim*sizeof(double);
   cudaMalloc((void**)&v_d,size);
   cudaMemcpy(v_d,v_h,size,cudaMemcpyHostToDevice);
   double* twonorm_d;
   cudaMalloc((void**)&twonorm_d,sizeof(double));

   if(dim == BS)
      for(int i=0; i<freq; i++)
         small_normalize_vector<<<1,BS>>>(v_d,dim,BSLog2,twonorm_d);
   else
   {
      int rf = ceil(((double) dim)/BS);
      int rfLog2 = ceil(log2((double) rf));
      for(int i=0; i<freq; i++)
         large_normalize_vector<<<1,BS>>>
            (v_d,dim,rf,rfLog2,BS,BSLog2,twonorm_d);
   }

   cudaMemcpy(v_h,v_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(twonorm,twonorm_d,sizeof(double),cudaMemcpyDeviceToHost);
}
