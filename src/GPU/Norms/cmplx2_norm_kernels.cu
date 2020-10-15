// Defines code of the functions in cmplx2_norm_kernels.h,
// to compute the 2-norm and normalize a complex vector,
// in double double precision,
// for vectors of small, medium, and large size.

#include "double_double_gpufun.cu"
#include "cmplx2_norm_kernels.h"

__global__ void small_normalize_vector
 ( double *vrehi, double *vrelo, double *vimhi, double *vimlo,
   int dim, int dimLog2, double *normhi, double *normlo )
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

__global__ void medium_normalize_vector
 ( double *vrehi, double *vrelo, double *vimhi, double *vimlo,
   int dim, int rnd, int rndLog2, int BS, int BSLog2,
   double *normhi, double *normlo )
{
   int j = threadIdx.x;
   int powTwo;
   int vBSind = 0;

   __shared__ double shvrehi[dd_shmemsize];
   __shared__ double shvrelo[dd_shmemsize];
   __shared__ double shvimhi[dd_shmemsize];
   __shared__ double shvimlo[dd_shmemsize];
   __shared__ double prdhi[dd_shmemsize];
   __shared__ double prdlo[dd_shmemsize];
   __shared__ double acchi[dd_shmemsize];
   __shared__ double acclo[dd_shmemsize];
   __shared__ double sumshi[maxrounds];
   __shared__ double sumslo[maxrounds];

   for(int i=0; i<rnd; i++)
   {
      if(vBSind+j >= dim)       // at last round some threads are excluded
      {
         acchi[j] = 0.0;
         acclo[j] = 0.0;
      }
      else
      {
         shvrehi[j] = vrehi[vBSind+j];  // reading into shared memory
         shvrelo[j] = vrelo[vBSind+j];
         shvimhi[j] = vimhi[vBSind+j]; 
         shvimlo[j] = vimlo[vBSind+j]; 
         ddg_sqr(shvrehi[j],shvrelo[j],&acchi[j],&acclo[j]);
         ddg_sqr(shvimhi[j],shvimlo[j],&prdhi[j],&prdlo[j]);
         ddg_inc(&acchi[j],&acclo[j],prdhi[j],prdlo[j]);
      }
      __syncthreads();
      powTwo = 1;                          // sum reduction
      for(int k=0; k < BSLog2; k++)
      {
         if((j%(powTwo*2)) == 0)
            if(j+powTwo < BS)
               ddg_inc(&acchi[j],&acclo[j],acchi[j+powTwo],acclo[j+powTwo]);
         powTwo = powTwo*2;
         __syncthreads();
      }
      // thread 0 copies the sum of this round in sums[i], the others wait
      if(j == 0)
      {
         sumshi[i] = acchi[0]; 
         sumslo[i] = acclo[0]; 
      }
      __syncthreads();
      vBSind = vBSind + BS;
   }
   powTwo = 1;                          // sum reduction
   for(int k=0; k < rndLog2; k++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < rnd)
            ddg_inc(&sumshi[j],&sumslo[j],sumshi[j+powTwo],sumslo[j+powTwo]);
      powTwo = powTwo*2;
      __syncthreads();
   }
   if(j == 0)
   {
      ddg_sqrt(sumshi[0],sumslo[0],&sumshi[0],&sumslo[0]);
      *normhi = sumshi[0];
      *normlo = sumslo[0];
   }
   __syncthreads();
   vBSind = 0;
   for(int i=0; i<rnd; i++)
   {
      if(vBSind+j < dim)
      {
         shvrehi[j] = vrehi[vBSind+j];       // read into shared memory
         shvrelo[j] = vrelo[vBSind+j];
         shvimhi[j] = vimhi[vBSind+j];
         shvimlo[j] = vimlo[vBSind+j];
         // normalize vector
         ddg_div(shvrehi[j],shvrelo[j],sumshi[0],sumslo[0],
                 &vrehi[vBSind+j],&vrelo[vBSind+j]);
         ddg_div(shvimhi[j],shvimlo[j],sumshi[0],sumslo[0],
                 &vimhi[vBSind+j],&vimlo[vBSind+j]);
      }
      __syncthreads();
      vBSind = vBSind + BS;
   }
}

__global__ void large_sum_the_squares
 ( double *vrehi, double *vrelo, double *vimhi, double *vimlo,
   int dim, double *sumshi, double *sumslo, int BS, int BSLog2 )
{
   const int i = blockIdx.x;
   const int j = threadIdx.x;
   const int k = i*BS + j;

   __shared__ double shvrehi[dd_shmemsize];
   __shared__ double shvrelo[dd_shmemsize];
   __shared__ double shvimhi[dd_shmemsize];
   __shared__ double shvimlo[dd_shmemsize];
   __shared__ double prdhi[dd_shmemsize];
   __shared__ double prdlo[dd_shmemsize];
   __shared__ double acchi[dd_shmemsize];
   __shared__ double acclo[dd_shmemsize];

   shvrehi[j] = vrehi[k];
   shvrelo[j] = vrelo[k];
   shvimhi[j] = vimhi[k];
   shvimlo[j] = vimlo[k];

   ddg_sqr(shvrehi[j],shvrelo[j],&acchi[j],&acclo[j]);
   ddg_sqr(shvimhi[j],shvimlo[j],&prdhi[j],&prdlo[j]);
   ddg_inc(&acchi[j],&acclo[j],prdhi[j],prdlo[j]);

   __syncthreads();

   int powTwo = 1;                          // sum reduction
   for(int L=0; L < BSLog2; L++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < BS)
            ddg_inc(&acchi[j],&acclo[j],acchi[j+powTwo],acclo[j+powTwo]);
      powTwo = powTwo*2;

      __syncthreads();
   }
   if(j == 0)                               // thread 0 writes the sum
   {
      sumshi[i] = acchi[0];
      sumslo[i] = acclo[0];
   }
}

__global__ void large_normalize_vector
 ( double *vrehi, double *vrelo, double *vimhi, double *vimlo,
   int dim, double *sumshi, double *sumslo, int nbsums, int nbsumsLog2,
   int BS, double *normhi, double *normlo )
{
   const int i = blockIdx.x;
   const int j = threadIdx.x;
   const int k = i*BS + j;

   __shared__ double shvrehi[dd_shmemsize];
   __shared__ double shvrelo[dd_shmemsize];
   __shared__ double shvimhi[dd_shmemsize];
   __shared__ double shvimlo[dd_shmemsize];

   if(j < nbsums)
   {
      shvrehi[j] = sumshi[j];
      shvrelo[j] = sumslo[j];
   }

   __syncthreads();

   int powTwo = 1;                          // sum reduction
   for(int L=0; L < nbsumsLog2; L++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < nbsums)
            ddg_inc(&shvrehi[j],&shvrelo[j],
                    shvrehi[j+powTwo],shvrelo[j+powTwo]);
      powTwo = powTwo*2;

      __syncthreads();
   }
   __syncthreads();                      // every thread 0 of all blocks
   if(j == 0)                            // compute the 2-norm and assigns
   {                                     // to the output parameter
      ddg_sqrt(shvrehi[0],shvrelo[0],normhi,normlo);
   }
   __syncthreads();                    

   if(k < dim)
   {
      shvrehi[j] = vrehi[k];
      shvrelo[j] = vrelo[k];
      shvimhi[j] = vimhi[k];
      shvimlo[j] = vimlo[k];

      ddg_div(shvrehi[j],shvrelo[j],*normhi,*normlo,&shvrehi[j],&shvrelo[j]);
      ddg_div(shvimhi[j],shvimlo[j],*normhi,*normlo,&shvimhi[j],&shvimlo[j]);

      vrehi[k] = shvrehi[j];
      vrelo[k] = shvrelo[j];
      vimhi[k] = shvimhi[j];
      vimlo[k] = shvimlo[j];
   }
}

void GPU_norm
 ( double *vrehi_h, double *vrelo_h, double *vimhi_h, double *vimlo_h,
   int dim, int freq, int BS, double *normhi, double *normlo, int blocked )
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
   else if(blocked == 0)
   {
      int rf = ceil(((double) dim)/BS);
      int rfLog2 = ceil(log2((double) rf));
      for(int i=0; i<freq; i++)
         medium_normalize_vector<<<1,BS>>>
            (vrehi_d,vrelo_d,vimhi_d,vimlo_d,dim,rf,rfLog2,
             BS,BSLog2,normhi_d,normlo_d);
   }
   else
   {
      const int nblocks = dim/BS;
      const int nblocksLog2 = ceil(log2((double) nblocks));
      double* sumshi_d; // high parts of sums of squares for each block
      double* sumslo_d; // low parts of sums of squares for each block
      size_t sums_size = nblocks*sizeof(double);
      cudaMalloc((void**)&sumshi_d,sums_size);
      cudaMalloc((void**)&sumslo_d,sums_size);
      for(int i=0; i<freq; i++)
      {
         large_sum_the_squares<<<nblocks,BS>>>
            (vrehi_d,vrelo_d,vimhi_d,vimlo_d,dim,sumshi_d,sumslo_d,BS,BSLog2);
         large_normalize_vector<<<nblocks,BS>>>
            (vrehi_d,vrelo_d,vimhi_d,vimlo_d,dim,sumshi_d,sumslo_d,
             nblocks,nblocksLog2,BS,normhi_d,normlo_d);
      }
   }
   cudaMemcpy(vrehi_h,vrehi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vrelo_h,vrelo_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vimhi_h,vimhi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vimlo_h,vimlo_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(normhi,normhi_d,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(normlo,normlo_d,sizeof(double),cudaMemcpyDeviceToHost);
}
