// This file contains the definition for the functions in dbl3_norm_kernels.h,
// to compute the 2-norm and normalize a triple double precision vector,
// for small, medium, and large vectors.

#include "double_double_gpufun.cu"
#include "quad_double_gpufun.cu"
#include "dbl4_norm_kernels.h"

__global__ void small_normalize_vector
 ( double *vhihi, double *vlohi, double *vhilo, double *vlolo,
   int dim, int dimLog2,
   double *normhihi, double *normlohi, double *normhilo, double *normlolo )
{
   int j = threadIdx.x;
   __shared__ double shvhihi[qd_shmemsize];
   __shared__ double shvlohi[qd_shmemsize];
   __shared__ double shvhilo[qd_shmemsize];
   __shared__ double shvlolo[qd_shmemsize];
   __shared__ double prdhihi[qd_shmemsize];
   __shared__ double prdlohi[qd_shmemsize];
   __shared__ double prdhilo[qd_shmemsize];
   __shared__ double prdlolo[qd_shmemsize];
   shvhihi[j] = vhihi[j];    // reading of vector into shared memory
   shvlohi[j] = vlohi[j];
   shvhilo[j] = vhilo[j];
   shvlolo[j] = vlolo[j];
   qdg_sqr(shvhihi[j],shvlohi[j],shvhilo[j],shvlolo[j],
           &prdhihi[j],&prdlohi[j],&prdhilo[j],&prdlolo[j]);
   __syncthreads();
   int powTwo = 1;                          // sum reduction
   for(int k=0; k < dimLog2; k++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < dim)
            qdg_inc(&prdhihi[j],&prdlohi[j],&prdhilo[j],&prdlolo[j],
                    prdhihi[j+powTwo],prdlohi[j+powTwo],
                    prdhilo[j+powTwo],prdlolo[j+powTwo]);
      powTwo = powTwo*2;
      __syncthreads();
   }
   // thread 0 computes the sqrt of the inner product, others wait
   if(j == 0)
      qdg_sqrt(prdhihi[0],prdlohi[0],prdhilo[0],prdlolo[0],
               &prdhihi[0],&prdlohi[0],&prdhilo[0],&prdlolo[0]); 
   if(j == 0)
   {
      *normhihi = prdhihi[0];
      *normlohi = prdlohi[0];
      *normhilo = prdhilo[0];
      *normlolo = prdlolo[0];
   }
   __syncthreads();
   qdg_div(shvhihi[j],shvlohi[j],shvhilo[j],shvlolo[j],
           prdhihi[0],prdlohi[0],prdhilo[0],prdlolo[0],
           &vhihi[j],&vlohi[j],&vhilo[j],&vlolo[j]);
}

__global__ void medium_normalize_vector
 ( double *vhihi, double *vlohi, double *vhilo, double *vlolo,
   int dim, int rnd, int rndLog2, int BS, int BSLog2,
   double *normhihi, double *normlohi, double *normhilo, double *normlolo )
{
   int j = threadIdx.x;
   int powTwo;
   int vBSind = 0;

   __shared__ double shvhihi[qd_shmemsize];
   __shared__ double shvlohi[qd_shmemsize];
   __shared__ double shvhilo[qd_shmemsize];
   __shared__ double shvlolo[qd_shmemsize];
   __shared__ double prdhihi[qd_shmemsize];
   __shared__ double prdlohi[qd_shmemsize];
   __shared__ double prdhilo[qd_shmemsize];
   __shared__ double prdlolo[qd_shmemsize];
   __shared__ double sumshihi[maxrounds];
   __shared__ double sumslohi[maxrounds];
   __shared__ double sumshilo[maxrounds];
   __shared__ double sumslolo[maxrounds];

   for(int i=0; i<rnd; i++)
   {
      if(vBSind+j >= dim)       // at last round some threads are excluded
      {
         prdhihi[j] = 0.0;
         prdlohi[j] = 0.0;
         prdhilo[j] = 0.0;
         prdlolo[j] = 0.0;
      }
      else
      {
         shvhihi[j] = vhihi[vBSind+j];  // reading into shared memory
         shvlohi[j] = vlohi[vBSind+j];
         shvhilo[j] = vhilo[vBSind+j];
         shvlolo[j] = vlolo[vBSind+j];
         qdg_sqr(shvhihi[j],shvlohi[j],shvhilo[j],shvlolo[j],
                 &prdhihi[j],&prdlohi[j],&prdhilo[j],&prdlolo[j]);
      }
      __syncthreads();
      powTwo = 1;                          // sum reduction
      for(int k=0; k < BSLog2; k++)
      {
         if((j%(powTwo*2)) == 0)
            if(j+powTwo < BS)
               qdg_inc(&prdhihi[j],&prdlohi[j],&prdhilo[j],&prdlolo[j],
                       prdhihi[j+powTwo],prdlohi[j+powTwo],
                       prdhilo[j+powTwo],prdlolo[j+powTwo]);
         powTwo = powTwo*2;
         __syncthreads();
      }
      // thread 0 copies the sum of this round in sums[i], the others wait
      if(j == 0)
      {
         sumshihi[i] = prdhihi[0]; 
         sumslohi[i] = prdlohi[0]; 
         sumshilo[i] = prdhilo[0]; 
         sumslolo[i] = prdlolo[0]; 
      }
      __syncthreads();
      vBSind = vBSind + BS;
   }
   powTwo = 1;                          // sum reduction
   for(int k=0; k < rndLog2; k++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < rnd)
            qdg_inc(&sumshihi[j],&sumslohi[j],&sumshilo[j],&sumslolo[j],
                    sumshihi[j+powTwo],sumslohi[j+powTwo],
                    sumshilo[j+powTwo],sumslolo[j+powTwo]);
      powTwo = powTwo*2;
      __syncthreads();
   }
   if(j == 0) qdg_sqrt(sumshihi[0],sumslohi[0],sumshilo[0],sumslolo[0],
                       &sumshihi[0],&sumslohi[0],&sumshilo[0],&sumslolo[0]); 
   if(j == 0)
   {
      *normhihi = sumshihi[0];
      *normlohi = sumslohi[0];
      *normhilo = sumshilo[0];
      *normlolo = sumslolo[0];
   }
   __syncthreads();
   vBSind = 0;
   for(int i=0; i<rnd; i++)
   {
      if(vBSind+j < dim)
      {
         shvhihi[j] = vhihi[vBSind+j];   // read into shared memory
         shvlohi[j] = vlohi[vBSind+j];   // and normalize the vector
         shvhilo[j] = vhilo[vBSind+j];
         shvlolo[j] = vlolo[vBSind+j];
         qdg_div(shvhihi[j],shvlohi[j],shvhilo[j],shvlolo[j],
                 sumshihi[0],sumslohi[0],sumshilo[0],sumslolo[0],
                 &vhihi[vBSind+j],&vlohi[vBSind+j],
                 &vhilo[vBSind+j],&vlolo[vBSind+j]);
      }
      __syncthreads();
      vBSind = vBSind + BS;
   }
}

__global__ void large_sum_the_squares
 ( double *vhihi, double *vlohi, double *vhilo, double *vlolo, int dim,
   double *sumshihi, double *sumslohi, double *sumshilo, double *sumslolo,
   int BS, int BSLog2 )
{
   const int i = blockIdx.x;
   const int j = threadIdx.x;
   const int k = i*BS + j;

   __shared__ double shvhihi[qd_shmemsize];
   __shared__ double shvlohi[qd_shmemsize];
   __shared__ double shvhilo[qd_shmemsize];
   __shared__ double shvlolo[qd_shmemsize];
   __shared__ double prdhihi[qd_shmemsize];
   __shared__ double prdlohi[qd_shmemsize];
   __shared__ double prdhilo[qd_shmemsize];
   __shared__ double prdlolo[qd_shmemsize];

   shvhihi[j] = vhihi[k];
   shvlohi[j] = vlohi[k];
   shvhilo[j] = vhilo[k];
   shvlolo[j] = vlolo[k];
   qdg_sqr(shvhihi[j],shvlohi[j],shvhilo[j],shvlolo[j],
           &prdhihi[j],&prdlohi[j],&prdhilo[j],&prdlolo[j]);

   __syncthreads();

   int powTwo = 1;                          // sum reduction
   for(int L=0; L < BSLog2; L++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < BS) 
            qdg_inc(&prdhihi[j],&prdlohi[j],&prdhilo[j],&prdlolo[j],
                    prdhihi[j+powTwo],prdlohi[j+powTwo],
                    prdhilo[j+powTwo],prdlolo[j+powTwo]);
      powTwo = powTwo*2;

      __syncthreads();
   }
   if(j == 0)                              // thread 0 writes the sum
   {
      sumshihi[i] = prdhihi[0];
      sumslohi[i] = prdlohi[0];
      sumshilo[i] = prdhilo[0];
      sumslolo[i] = prdlolo[0];
   }
}

__global__ void large_normalize_vector
 ( double *vhihi, double *vlohi, double *vhilo, double *vlolo, int dim,
   double *sumshihi, double *sumslohi, double *sumshilo, double *sumslolo,
   int nbsums, int nbsumsLog2, int BS,
   double *normhihi, double *normlohi, double *normhilo, double *normlolo )
{
   const int i = blockIdx.x;
   const int j = threadIdx.x;
   const int k = i*BS + j;

   __shared__ double shvhihi[qd_shmemsize];
   __shared__ double shvlohi[qd_shmemsize];
   __shared__ double shvhilo[qd_shmemsize];
   __shared__ double shvlolo[qd_shmemsize];

   if(j < nbsums)
   {
      shvhihi[j] = sumshihi[j];
      shvlohi[j] = sumslohi[j];
      shvhilo[j] = sumshilo[j];
      shvlolo[j] = sumslolo[j];
   }
   __syncthreads();

   int powTwo = 1;                          // sum reduction
   for(int L=0; L < nbsumsLog2; L++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < nbsums)
            qdg_inc(&shvhihi[j],&shvlohi[j],&shvhilo[j],&shvlolo[j],
                    shvhihi[j+powTwo],shvlohi[j+powTwo],
                    shvhilo[j+powTwo],shvlolo[j+powTwo]);
      powTwo = powTwo*2;

      __syncthreads();
   }
   __syncthreads();                    // every thread 0 of all blocks
   if(j == 0)                          // compute the 2-norm and assigns
   {                                   // to the output parameter
      qdg_sqrt(shvhihi[0],shvlohi[0],shvhilo[0],shvlolo[0],
               normhihi,normlohi,normhilo,normlolo); 
   }
   __syncthreads();                    // to the output parameter

   if(k < dim)
   {
      shvhihi[j] = vhihi[k];
      shvlohi[j] = vlohi[k];
      shvhilo[j] = vhilo[k];
      shvlolo[j] = vlolo[k];
      qdg_div(shvhihi[j],shvlohi[j],shvhilo[j],shvlolo[j],
              *normhihi,*normlohi,*normhilo,*normlolo,
              &shvhihi[j],&shvlohi[j],&shvhilo[j],&shvlolo[j]);
      vhihi[k] = shvhihi[j];
      vlohi[k] = shvlohi[j];
      vhilo[k] = shvhilo[j];
      vlolo[k] = shvlolo[j];
   }
}

void GPU_norm
 ( double *vhihi_h, double *vlohi_h, double *vhilo_h, double *vlolo_h,
   int dim, int freq, int BS,
   double *normhihi, double *normlohi, double *normhilo,double *normlolo,
   int blocked )
{
   int BSLog2 = ceil(log2((double) BS)); // ceil for sum reduction

   double* vhihi_d;                   // allocate for vector on device
   double* vlohi_d;
   double* vhilo_d;
   double* vlolo_d;
   size_t size = dim*sizeof(double);
   cudaMalloc((void**)&vhihi_d,size);
   cudaMalloc((void**)&vlohi_d,size);
   cudaMalloc((void**)&vhilo_d,size);
   cudaMalloc((void**)&vlolo_d,size);
   cudaMemcpy(vhihi_d,vhihi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vlohi_d,vlohi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vhilo_d,vhilo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vlolo_d,vlolo_h,size,cudaMemcpyHostToDevice);
   double* normhihi_d;
   double* normlohi_d;
   double* normhilo_d;
   double* normlolo_d;
   cudaMalloc((void**)&normhihi_d,sizeof(double));
   cudaMalloc((void**)&normlohi_d,sizeof(double));
   cudaMalloc((void**)&normhilo_d,sizeof(double));
   cudaMalloc((void**)&normlolo_d,sizeof(double));

   if(dim == BS)
   {
      for(int i=0; i<freq; i++)
         small_normalize_vector<<<1,BS>>>
            (vhihi_d,vlohi_d,vhilo_d,vlolo_d,dim,BSLog2,
             normhihi_d,normlohi_d,normhilo_d,normlolo_d);
   }
   else if(blocked == 0)
   {
      const int rf = ceil(((double) dim)/BS);
      const int rfLog2 = ceil(log2((double) rf));
      for(int i=0; i<freq; i++)
         medium_normalize_vector<<<1,BS>>>
            (vhihi_d,vlohi_d,vhilo_d,vlolo_d,dim,rf,rfLog2,BS,BSLog2,
             normhihi_d,normlohi_d,normhilo_d,normlolo_d);
   }
   else
   {
      const int nblocks = dim/BS;
      const int nblocksLog2 = ceil(log2((double) nblocks));
      double* sumshihi_d; // sums of squares for each block
      double* sumslohi_d; // second highest parts of sums of squares
      double* sumshilo_d; // second lowest parts of sums of squares
      double* sumslolo_d; // lowest parts of sums of squares
      size_t sums_size = nblocks*sizeof(double);
      cudaMalloc((void**)&sumshihi_d,sums_size);
      cudaMalloc((void**)&sumslohi_d,sums_size);
      cudaMalloc((void**)&sumshilo_d,sums_size);
      cudaMalloc((void**)&sumslolo_d,sums_size);
      for(int i=0; i<freq; i++)
      {
         large_sum_the_squares<<<nblocks,BS>>>
            (vhihi_d,vlohi_d,vhilo_d,vlolo_d,dim,
             sumshihi_d,sumslohi_d,sumshilo_d,sumslolo_d,BS,BSLog2);
         large_normalize_vector<<<nblocks,BS>>>
            (vhihi_d,vlohi_d,vhilo_d,vlolo_d,dim,
             sumshihi_d,sumslohi_d,sumshilo_d,sumslolo_d,
             nblocks,nblocksLog2,BS,
             normhihi_d,normlohi_d,normhilo_d,normlolo_d);
      }
   }
   cudaMemcpy(vhihi_h,vhihi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vlohi_h,vlohi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vhilo_h,vhilo_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vlolo_h,vlolo_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(normhihi,normhihi_d,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(normlohi,normlohi_d,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(normhilo,normhilo_d,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(normlolo,normlolo_d,sizeof(double),cudaMemcpyDeviceToHost);
}
