// Kernels for GPU version of the solver

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
#define shmemsize 256
#elif(p == 1)
typedef gdd_real T;
#define shmemsize 128
#else
typedef gqd_real T;
#define shmemsize 85
#endif

using namespace std;

// for double double complex, 127 as high we can go with shared memory
// #define shmemsize 127
// for quad double complex, let us try 63 (64 does not work)
// #define shmemsize 63

/********** Back Substitution Kernel *********/

__global__ void BackS_CoalescI
 ( complex<T>* R, complex<T>* sol, int dim, int Rdim )
{

   // The upper triangular matrix R is outputted by MGS
   // column by column starting from the last column.
   // This allows to achieve coalesced reading of R by threads of
   // the working block of threads of the back substitution kernel.

   int i = threadIdx.x;
   __shared__ complex<T> solH[shmemsize];
   __shared__ complex<T> curr_y_update[shmemsize];
   __shared__ complex<T> R_CURR_COL[shmemsize];

   curr_y_update[i] = R[i];
   int curr_ind;

   for (int u=0; u<dim; u++)
   {
      __syncthreads();
      if(i < (dim-u))
      {
         __syncthreads();
         curr_ind = 0.5*(u+1)*(2*(dim+1)-u) + i;
         __syncthreads();
         R_CURR_COL[i]=R[curr_ind];
         __syncthreads();
         if(i == (dim-u-1)) solH[i] = curr_y_update[i]/R_CURR_COL[i];
         __syncthreads();
         if(i < dim-u-1)
            curr_y_update[i] = curr_y_update[i] - solH[dim-u-1]*R_CURR_COL[i];
         __syncthreads();
      }
      __syncthreads();
   }
   sol[i] = solH[i];
}

/******** Modified Gram-Schmidt Kernel (MGS) *********/

// Blocks of threads of MGS kernel remove components of vectors with indexes
// larger than the index Pv in the direction of the vector with the index PV.
// Threads of the block with the ID = blockIdx.x cooperatively compute 
// the projection of the vector with the index (Pv+blockIdx.x)
// on the Pv-th vector, and subsequently subtruct it from 
// the (Pv+blockIdx.)-th vector.  

// Two almost identical versions of MGS kernel --- GSvu and GSuu --- 
// are used to guarantee purity of experiment 
// as measuring speedups of the linear solver. 
// Speedups are established when solver run multiple times (legitimacy 
// of this approach is described in READ_ME).
// If only the inplace MGS kernel GSuu is used, in all subsequent runs of 
// the solver after the very first one, orthogonalizations will be
// performed on already orthonormal set of vectors.
// In the combined path tracker only GSvv will be used.  	
// GSvu differs from GSvv only by outputting the adjusted vectors 
// into different memory locations 

__global__ void GSvu
 ( complex<T>* v, complex<T>* u, complex<T>* R, int dim, int dimR,
   int k, int Pv, int dimLog, int r )
{
   int i = (Pv+blockIdx.x)*blockDim.x + threadIdx.x; 
   // index of the assigned to thread coordinate of pivotal vector Pv
   //  in the global memory array v1,v2,...vk 
   int j = threadIdx.x;
   int L = blockIdx.x;

   // This indexing helps to organize the output of the array
   // representing the upper triangular matrix R in a way
   // that back sustitution kernel reads in coalesced fashion.
   int Rind = (dimR-1) - (Pv*(Pv+1))/2 - (L*(L + 2*Pv + 3))/2;
   complex<T> sum;

   // Coalesced reading of the pivotal vector
   // and of the assigned to the block vector into the shared memory 
   __shared__ complex<T> vectors[3][shmemsize];
   vectors[0][j] = v[Pv*dim+j];
   vectors[1][j] = v[i];
   __syncthreads();
   
   // Computing the norm of the pivotal vector
   // each block of threads redundantly computes it.
   vectors[2][j] = vectors[0][j]*vectors[0][j].adj();
   __syncthreads();

   int powTwo = 1;                          // sum reduction
   for(int kk=0; kk < dimLog; kk++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < dim)
            vectors[2][j] = vectors[2][j] + vectors[2][j+powTwo];
      powTwo = powTwo*2;
      __syncthreads();
   }
   // thread 0 computes the sqrt of the inner product 
   // as the other threads of the block wait 
   if(threadIdx.x==0) vectors[2][0].real = sqrt(vectors[2][0].real); 
   __syncthreads();
   
   // block 0 outputs the norm of the pivotal vector into
   // the array representing the upper triangular matrix R
   // it also normalizes the pivotal vector, and outputs it 
   // into the global memory in a coalesced way  

   if(L == 0)
   {
      vectors[0][j] = vectors[0][j]/vectors[2][0];
      u[Pv*dim+j] = vectors[0][j];
      if(j == 0) R[Rind] = vectors[2][0];    
   }

   // Each block of threads with an index bigger than zero
   // redundantly normalizes the pivotal vector; then it computes
   // the inner product of the assigned to it vector with the
   // normalized pivotal vector; computes the projection
   // of the assigned vector on the pivotal vector, and subtructs
   // the computed projection from the assigned vector;
   // all these computations are done cooperatively by threads of the block
   // finally thread 0 of the block outputs the computed inner product into
   // the array representing R, and also the threads of the block output 
   // the adjusted assigned vector into the global memory in coalesced fashion.

   if(blockIdx.x != 0)
   {
      vectors[0][j] = vectors[0][j]/vectors[2][0];
      vectors[2][j] = vectors[1][j]*vectors[0][j].adj();
      __syncthreads();
      powTwo = 1;                                // sum reduction
      for(int kk=0; kk<dimLog; kk++)
      {
         if((j%(powTwo*2)) == 0)
            if(j+powTwo < dim)
               vectors[2][j] = vectors[2][j] + vectors[2][j+powTwo];
          powTwo = powTwo*2;
          __syncthreads();
      }
      sum = vectors[2][0];
      vectors[2][j] = sum * vectors[0][j];
      vectors[1][j] = vectors[1][j] - vectors[2][j];

      u[i] = vectors[1][j];
      if(j == 0) R[Rind] = sum;
   }
}

__global__ void GSuu
 ( complex<T>* v, complex<T>* R, int dim, int dimR,
   int k, int Pv, int dimLog, int r )
{
   int i = (Pv + blockIdx.x)*blockDim.x + threadIdx.x; 
   // index of the assigned to thread coordinate of pivotal vector Pv 
   // in the global memory array v1,v2,...vk
   int j = threadIdx.x;
   int L = blockIdx.x;

   // This indexing helps to organize the output of the array
   // representing the upper triangular matrix R
   // in a way that back sustitution kernel reads its in coalesced fashion
   int Rind = (dimR-1) - (Pv*(Pv+1))/2 - (L*(L + 2*Pv + 3))/2;
   complex<T> sum;

   // Coalesced reading of the pivotal vector
   // and of the assigned to the block vector into the shared memory
   __shared__ complex<T> vectors[3][shmemsize];
   vectors[0][j] = v[Pv*dim+j];
   vectors[1][j] = v[i];
   __syncthreads();

   // computing the norm of the pivotal vector
   // each block of threads redundantly computes it

   vectors[2][j] = vectors[0][j]*vectors[0][j].adj();
   __syncthreads();
   int powTwo = 1;                            // sum reduction
   for(int u=0; u<dimLog; u++)
   {
      if((j%(powTwo*2))==0)
         if(j+powTwo < dim)
            vectors[2][j] = vectors[2][j] + vectors[2][j+powTwo];
      powTwo = powTwo*2;
      __syncthreads();
   }
   // thread 0 computes the sqrt of the inner product 
   // as the other threads of the block wait
   if(threadIdx.x==0) vectors[2][0].real = sqrt(vectors[2][0].real);
   __syncthreads();

   // block 0 outputs the norm of the pivotal vector into
   // the array representing the upper triangular matrix R
   // it also normalizes the pivotal vector, and outputs it
   // into the global memory in a coalesced way

   if(blockIdx.x == 0)
   {
      vectors[0][j] = vectors[0][j]/vectors[2][0];
      v[Pv*dim+j] = vectors[0][j];
      if(j == 0) R[Rind] = vectors[2][0];
   }
   // each block of threads with an index bigger than zero
   // redundantly normalizes the pivotal vector; then it computes
   // the inner product of the assigned to it vector with the
   // normalized pivotal vector; computes the projection
   // of the assigned vector on the pivotal vector, and subtructs
   // the computed projection from the assigned vector;
   // all these computations are done cooperatively by threads of the block
   // finally thread 0 of the block outputs the computed inner product into
   // the array representing R, and also the threads of the block output 
   // the adjusted assigned vector into the global memory in coalesced fashion.

   // Having one single statement "if(blockIdx == 0)" respectively "!= 0"
   // followed by specific code for the blocks does not work for quad doubles.

   if(blockIdx.x != 0)
   {
      vectors[0][j] = vectors[0][j]/vectors[2][0];
      vectors[2][j] = vectors[1][j]*vectors[0][j].adj();
      __syncthreads();
      powTwo = 1;                        // sum reduction
      for(int u=0; u<dimLog; u++)
      {
         if((j%(powTwo*2)) == 0)
            if(j+powTwo < dim)
               vectors[2][j] = vectors[2][j] + vectors[2][j+powTwo];
          powTwo = powTwo*2;
          __syncthreads();
      }
      sum = vectors[2][0];
      vectors[2][j] = sum * vectors[0][j];
      vectors[1][j] = vectors[1][j] - vectors[2][j];
      v[i] = vectors[1][j];
      if(j==0) R[Rind] = sum;
   }
}

void GPU_GS
 ( complex<T>* v_h, complex<T>* R_h, complex<T>* sol_h, 
   int dim, int dimR, int k, int r, int BS )
{
   int dimLog = ceil(log2((double)dim)); // ceil for sum reduction

   complex<T>* v_d;
   size_t size = dim*k*sizeof(complex<T>);

   cudaMalloc((void**)&v_d,size);

   cudaMemcpy(v_d,v_h,size,cudaMemcpyHostToDevice);

   complex<T>* u_d;
   cudaMalloc((void**)&u_d,size);

   complex<T>* R_d;
   size_t sizeR = dimR*sizeof(complex<T>);
   cudaMalloc((void**)&R_d, sizeR);

   complex<T>* sol_d;
   size_t sizeSol = dim*sizeof(complex<T>);
   cudaMalloc((void**)&sol_d, sizeSol);

   for(int i=0; i<r; i++)
   {
      GSvu<<<k,BS>>>(v_d,u_d,R_d,dim,dimR,k,0,dimLog,r);
      for(int Pv=1; Pv<k; Pv++)
         GSuu<<<k-Pv,BS>>>(u_d,R_d,dim,dimR,k,Pv,dimLog,r);
      // for Table III in the paper we do not do the back substitution
      BackS_CoalescI<<<1,BS>>>(R_d,sol_d,dim,dimR);
   }

   cudaMemcpy(v_h,u_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(R_h,R_d,sizeR,cudaMemcpyDeviceToHost);
   cudaMemcpy(sol_h,sol_d,sizeSol,cudaMemcpyDeviceToHost);
}
