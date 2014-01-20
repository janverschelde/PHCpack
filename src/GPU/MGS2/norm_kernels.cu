// kernels for GPU version of the calculation of the norm of a complex vector

#include <iostream>
#include <gqd.cu>
#include <assert.h>
#include <cstdio>
#include <gqd_type.h>
#include "complex.h"

#define d  0 
#define dd 1
#define qd 2

#ifdef precision
#define p precision
#else
#define p 0
#endif

#if(p == 0)
typedef double T;
#define shmemsize 512
#elif(p == 1)
typedef gdd_real T;
#define shmemsize 256
#else
typedef gqd_real T;
#define shmemsize 128
#endif

#define maxrounds 32

using namespace std;

__global__ void small_normalize_vector
 ( complex<T>* v, int dim, int dimLog2, T* twonorm )
{
   int j = threadIdx.x;
   __shared__ complex<T> shv[shmemsize];
   __shared__ T prd[shmemsize];
   shv[j] = v[j];    // reading of vector into shared memory
   prd[j] = shv[j].real*shv[j].real + shv[j].imag*shv[j].imag;
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
 ( complex<T>* v, int dim, int rnd, int rndLog2, int BS, int BSLog2,
   T* twonorm )
{
   int j = threadIdx.x;
   int powTwo;
   int vBSind = 0;

   __shared__ complex<T> shv[shmemsize];
   __shared__ T prd[shmemsize];
   __shared__ T sums[maxrounds];
   __shared__ complex<T> zero;

   zero.init(0.0,0.0);
   for(int i=0; i<rnd; i++)
   {
      if(vBSind+j >= dim)       // at last round some threads are excluded
         prd[j] = zero.real;
      else
      {
         shv[j] = v[vBSind+j];  // reading of vector into shared memory
         prd[j] = shv[j].real*shv[j].real + shv[j].imag*shv[j].imag;
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

void GPU_norm ( complex<T>* v_h, int dim, int freq, int BS, T* twonorm )
{
   int BSLog2 = ceil(log2((double) BS)); // ceil for sum reduction

   complex<T>* v_d;                      // allocate for vector on device
   size_t size = dim*sizeof(complex<T>);
   cudaMalloc((void**)&v_d,size);
   cudaMemcpy(v_d,v_h,size,cudaMemcpyHostToDevice);
   T* twonorm_d;
   cudaMalloc((void**)&twonorm_d,sizeof(T));

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
   cudaMemcpy(twonorm,twonorm_d,sizeof(T),cudaMemcpyDeviceToHost);
}
