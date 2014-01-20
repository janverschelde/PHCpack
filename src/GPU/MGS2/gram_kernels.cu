// Kernels to compute the Gram matrix of a sequence of complex vectors
// are defined below.  Because this code prepares for the revision of the
// modified Gram-Schmidt method, a double indexed block structure is not
// applied, instead, kernels are launced with a double loop.

#include <gqd.cu>
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

__global__ void small_gram
 ( complex<T>* v, complex<T>* g, int pivot, int dim, int dimLog2 )
{
   int block = blockIdx.x;
   int j = threadIdx.x;
   __shared__ complex<T> x[shmemsize];
   x[j] = v[pivot*dim+j];
   x[j] = x[j].adj()*v[block*dim+j];
   __syncthreads();
   int powTwo = 1;                          // sum reduction
   for(int k=0; k < dimLog2; k++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < dim) x[j] = x[j] + x[j+powTwo];
      powTwo = powTwo*2;
      __syncthreads();
   }
   if(j == 0) g[pivot*dim+block] = x[0];
}

__global__ void large_gram
 ( complex<T>* v, complex<T>* g, int pivot, int dim,
   int rnd, int rndLog2, int BS, int BSLog2 )
{
   int block = blockIdx.x;
   int j = threadIdx.x;
   int powTwo;
   int vBSind = 0;

   __shared__ complex<T> x[shmemsize];
   __shared__ complex<T> sums[maxrounds];

   for(int i=0; i<rnd; i++)
   {
      if(vBSind+j >= dim)      // exclude some threads in the last round
         x[j].init(0.0,0.0);
      else
      {
         x[j] = v[pivot*dim+j+vBSind];
         x[j] = x[j].adj()*v[block*dim+j+vBSind];
      }
      __syncthreads();
      powTwo = 1;                      // sum reduction within a block
      for(int k=0; k < BSLog2; k++)
      {
         if((j%(powTwo*2)) == 0)
            if(j+powTwo < BS) x[j] = x[j] + x[j+powTwo];
         powTwo = powTwo*2;
         __syncthreads();
      }
      if(j == 0) sums[i] = x[0];
      __syncthreads();
      vBSind = vBSind + BS;
   }
   powTwo = 1;                        // sum reduction of the partial sums
   for(int k=0; k < rndLog2; k++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < rnd)
            sums[j] = sums[j] + sums[j+powTwo];
      powTwo = powTwo*2;
      __syncthreads();
   }
   if(j == 0) g[pivot*dim+block] = sums[0];
}

void GPU_gram ( complex<T>* v_h, complex<T> *g_h, int dim, int freq, int BS )
{
   int BSLog2 = ceil(log2((double) BS)); // ceil for sum reduction
   complex<T>* v_d; // memory allocation on device and data transfer
   size_t size = dim*dim*sizeof(complex<T>);
   cudaMalloc((void**)&v_d,size);
   cudaMemcpy(v_d,v_h,size,cudaMemcpyHostToDevice);
   complex<T>* g_d;
   cudaMalloc((void**)&g_d,size);

   if(dim == BS)
      for(int i=0; i<freq; i++)
         for(int piv=0; piv<dim; piv++)
            small_gram<<<BS,BS>>>(v_d,g_d,piv,BS,BSLog2);
   else
   {
      int rf = ceil(((double) dim)/BS);
      int rfLog2 = ceil(log2((double) rf));

      for(int i=0; i<freq; i++)
         for(int piv=0; piv<dim; piv++)
            large_gram<<<dim,BS>>>(v_d,g_d,piv,dim,rf,rfLog2,BS,BSLog2);
   }

   cudaMemcpy(g_h,g_d,size,cudaMemcpyDeviceToHost);
}
