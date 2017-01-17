// This file contains the definition for the functions in all_norm_kernels.h.
// All kernels for GPU version of the calculation of the norm of a 
// complex vector, for double, double double, and quad double precision.
// The file is obtained by copying the norm_kernels.cu three times,
// once for double, once for qdd_real, and once for gqd_real,
// each time replacing the generic type T by double, gdd_real, and gqd_real,
// respectively for double, double double, and quad double precision.
// The three different upper bounds on the shared memories
// have each a separate name.

#include <iostream>
#include <cmath>
#include <assert.h>
#include <gqd.cu>
#include <gqd_type.h>
#include <cstdio>
#include "complex.h"
#include "vector_types.h"
#include "vector_functions.h"

#define d_shmemsize 512
#define dd_shmemsize 256
#define qd_shmemsize 128

#define maxrounds 32

using namespace std;

// for double precision :

__global__ void small_normalize_vector
 ( complex<double>* v, int dim, int dimLog2, double* twonorm )
{
   int j = threadIdx.x;
   __shared__ complex<double> shv[d_shmemsize];
   __shared__ double prd[d_shmemsize];
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
 ( complex<double>* v, int dim, int rnd, int rndLog2, int BS, int BSLog2,
   double* twonorm )
{
   int j = threadIdx.x;
   int powTwo;
   int vBSind = 0;

   __shared__ complex<double> shv[d_shmemsize];
   __shared__ double prd[d_shmemsize];
   __shared__ double sums[maxrounds];
   __shared__ complex<double> zero;

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

void GPU_norm
 ( complex<double>* v_h, int dim, int freq, int BS, double* twonorm )
{
   int BSLog2 = ceil(log2((double) BS)); // ceil for sum reduction

   complex<double>* v_d;                   // allocate for vector on device
   size_t size = dim*sizeof(complex<double>);
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

// for double double precision :

__global__ void small_normalize_vector
 ( complex<gdd_real>* v, int dim, int dimLog2, gdd_real* twonorm )
{
   int j = threadIdx.x;
   __shared__ complex<gdd_real> shv[dd_shmemsize];
   __shared__ gdd_real prd[dd_shmemsize];
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
 ( complex<gdd_real>* v, int dim, int rnd, int rndLog2, int BS, int BSLog2,
   gdd_real* twonorm )
{
   int j = threadIdx.x;
   int powTwo;
   int vBSind = 0;

   __shared__ complex<gdd_real> shv[dd_shmemsize];
   __shared__ gdd_real prd[dd_shmemsize];
   __shared__ gdd_real sums[maxrounds];
   __shared__ complex<gdd_real> zero;

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

void GPU_norm
 ( complex<gdd_real>* v_h, int dim, int freq, int BS, gdd_real* twonorm )
{
   int BSLog2 = ceil(log2((double) BS)); // ceil for sum reduction

   complex<gdd_real>* v_d;                 // allocate for vector on device
   size_t size = dim*sizeof(complex<gdd_real>);
   cudaMalloc((void**)&v_d,size);
   cudaMemcpy(v_d,v_h,size,cudaMemcpyHostToDevice);
   gdd_real* twonorm_d;
   cudaMalloc((void**)&twonorm_d,sizeof(gdd_real));

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
   cudaMemcpy(twonorm,twonorm_d,sizeof(gdd_real),cudaMemcpyDeviceToHost);
}

// for quad double precision :

__global__ void small_normalize_vector
 ( complex<gqd_real>* v, int dim, int dimLog2, gqd_real* twonorm )
{
   int j = threadIdx.x;
   __shared__ complex<gqd_real> shv[qd_shmemsize];
   __shared__ gqd_real prd[qd_shmemsize];
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
 ( complex<gqd_real>* v, int dim, int rnd, int rndLog2, int BS, int BSLog2,
   gqd_real* twonorm )
{
   int j = threadIdx.x;
   int powTwo;
   int vBSind = 0;

   __shared__ complex<gqd_real> shv[qd_shmemsize];
   __shared__ gqd_real prd[qd_shmemsize];
   __shared__ gqd_real sums[maxrounds];
   __shared__ complex<gqd_real> zero;

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

void GPU_norm
 ( complex<gqd_real>* v_h, int dim, int freq, int BS, gqd_real* twonorm )
{
   int BSLog2 = ceil(log2((double) BS)); // ceil for sum reduction

   complex<gqd_real>* v_d;                 // allocate for vector on device
   size_t size = dim*sizeof(complex<gqd_real>);
   cudaMalloc((void**)&v_d,size);
   cudaMemcpy(v_d,v_h,size,cudaMemcpyHostToDevice);
   gqd_real* twonorm_d;
   cudaMalloc((void**)&twonorm_d,sizeof(gqd_real));

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
   cudaMemcpy(twonorm,twonorm_d,sizeof(gqd_real),cudaMemcpyDeviceToHost);
}
