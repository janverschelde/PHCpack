// This file contains the definition for the functions in dbl8_norm_kernels.h,
// to compute the 2-norm and normalize an octo double precision vector,
// for small, medium, and large vectors.

#include "double_double_gpufun.cu"
#include "quad_double_gpufun.cu"
#include "octo_double_gpufun.cu"
#include "dbl8_norm_kernels.h"

__global__ void small_normalize_vector
 ( double *vhihihi, double *vlohihi, double *vhilohi, double *vlolohi,
   double *vhihilo, double *vlohilo, double *vhilolo, double *vlololo,
   int dim, int dimLog2, double *normhihihi, double *normlohihi,
   double *normhilohi, double *normlolohi, double *normhihilo,
   double *normlohilo, double *normhilolo, double *normlololo )
{
   int j = threadIdx.x;
   __shared__ double shvhihihi[od_shmemsize];
   __shared__ double shvlohihi[od_shmemsize];
   __shared__ double shvhilohi[od_shmemsize];
   __shared__ double shvlolohi[od_shmemsize];
   __shared__ double shvhihilo[od_shmemsize];
   __shared__ double shvlohilo[od_shmemsize];
   __shared__ double shvhilolo[od_shmemsize];
   __shared__ double shvlololo[od_shmemsize];
   __shared__ double prdhihihi[od_shmemsize];
   __shared__ double prdlohihi[od_shmemsize];
   __shared__ double prdhilohi[od_shmemsize];
   __shared__ double prdlolohi[od_shmemsize];
   __shared__ double prdhihilo[od_shmemsize];
   __shared__ double prdlohilo[od_shmemsize];
   __shared__ double prdhilolo[od_shmemsize];
   __shared__ double prdlololo[od_shmemsize];
   shvhihihi[j] = vhihihi[j];    // reading of vector into shared memory
   shvlohihi[j] = vlohihi[j];
   shvhilohi[j] = vhilohi[j];
   shvlolohi[j] = vlolohi[j];
   shvhihilo[j] = vhihilo[j];
   shvlohilo[j] = vlohilo[j];
   shvhilolo[j] = vhilolo[j];
   shvlololo[j] = vlololo[j];
   odg_sqr(shvhihihi[j],shvlohihi[j],shvhilohi[j],shvlolohi[j],
           shvhihilo[j],shvlohilo[j],shvhilolo[j],shvlololo[j],
           &prdhihihi[j],&prdlohihi[j],&prdhilohi[j],&prdlolohi[j],
           &prdhihilo[j],&prdlohilo[j],&prdhilolo[j],&prdlololo[j]);
   __syncthreads();
   int powTwo = 1;                          // sum reduction
   for(int k=0; k < dimLog2; k++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < dim)
            odg_inc(&prdhihihi[j],&prdlohihi[j],&prdhilohi[j],&prdlolohi[j],
                    &prdhihilo[j],&prdlohilo[j],&prdhilolo[j],&prdlololo[j],
                    prdhihihi[j+powTwo],prdlohihi[j+powTwo],
                    prdhilohi[j+powTwo],prdlolohi[j+powTwo],
                    prdhihilo[j+powTwo],prdlohilo[j+powTwo],
                    prdhilolo[j+powTwo],prdlololo[j+powTwo]);
      powTwo = powTwo*2;
      __syncthreads();
   }
   // thread 0 computes the sqrt of the inner product, others wait
   if(j == 0)
      odg_sqrt(prdhihihi[0],prdlohihi[0],prdhilohi[0],prdlolohi[0],
               prdhihilo[0],prdlohilo[0],prdhilolo[0],prdlololo[0],
               &prdhihihi[0],&prdlohihi[0],&prdhilohi[0],&prdlolohi[0],
               &prdhihilo[0],&prdlohilo[0],&prdhilolo[0],&prdlololo[0]); 
   if(j == 0)
   {
      *normhihihi = prdhihihi[0];
      *normlohihi = prdlohihi[0];
      *normhilohi = prdhilohi[0];
      *normlolohi = prdlolohi[0];
      *normhihilo = prdhihilo[0];
      *normlohilo = prdlohilo[0];
      *normhilolo = prdhilolo[0];
      *normlololo = prdlololo[0];
   }
   __syncthreads();
   odg_div(shvhihihi[j],shvlohihi[j],shvhilohi[j],shvlolohi[j],
           shvhihilo[j],shvlohilo[j],shvhilolo[j],shvlololo[j],
           prdhihihi[0],prdlohihi[0],prdhilohi[0],prdlolohi[0],
           prdhihilo[0],prdlohilo[0],prdhilolo[0],prdlololo[0],
           &vhihihi[j],&vlohihi[j],&vhilohi[j],&vlolohi[j],
           &vhihilo[j],&vlohilo[j],&vhilolo[j],&vlololo[j]);
}

__global__ void medium_normalize_vector
 ( double *vhihihi, double *vlohihi, double *vhilohi, double *vlolohi,
   double *vhihilo, double *vlohilo, double *vhilolo, double *vlololo,
   int dim, int rnd, int rndLog2, int BS, int BSLog2,
   double *normhihihi, double *normlohihi, double *normhilohi,
   double *normlolohi, double *normhihilo, double *normlohilo,
   double *normhilolo, double *normlololo )
{
   int j = threadIdx.x;
   int powTwo;
   int vBSind = 0;

   __shared__ double shvhihihi[od_shmemsize];
   __shared__ double shvlohihi[od_shmemsize];
   __shared__ double shvhilohi[od_shmemsize];
   __shared__ double shvlolohi[od_shmemsize];
   __shared__ double shvhihilo[od_shmemsize];
   __shared__ double shvlohilo[od_shmemsize];
   __shared__ double shvhilolo[od_shmemsize];
   __shared__ double shvlololo[od_shmemsize];
   __shared__ double prdhihihi[od_shmemsize];
   __shared__ double prdlohihi[od_shmemsize];
   __shared__ double prdhilohi[od_shmemsize];
   __shared__ double prdlolohi[od_shmemsize];
   __shared__ double prdhihilo[od_shmemsize];
   __shared__ double prdlohilo[od_shmemsize];
   __shared__ double prdhilolo[od_shmemsize];
   __shared__ double prdlololo[od_shmemsize];
   __shared__ double sumshihihi[maxrounds];
   __shared__ double sumslohihi[maxrounds];
   __shared__ double sumshilohi[maxrounds];
   __shared__ double sumslolohi[maxrounds];
   __shared__ double sumshihilo[maxrounds];
   __shared__ double sumslohilo[maxrounds];
   __shared__ double sumshilolo[maxrounds];
   __shared__ double sumslololo[maxrounds];

   for(int i=0; i<rnd; i++)
   {
      if(vBSind+j >= dim)       // at last round some threads are excluded
      {
         prdhihihi[j] = 0.0;
         prdlohihi[j] = 0.0;
         prdhilohi[j] = 0.0;
         prdlolohi[j] = 0.0;
         prdhihilo[j] = 0.0;
         prdlohilo[j] = 0.0;
         prdhilolo[j] = 0.0;
         prdlololo[j] = 0.0;
      }
      else
      {
         shvhihihi[j] = vhihihi[vBSind+j];  // reading into shared memory
         shvlohihi[j] = vlohihi[vBSind+j];
         shvhilohi[j] = vhilohi[vBSind+j];
         shvlolohi[j] = vlolohi[vBSind+j];
         shvhihilo[j] = vhihilo[vBSind+j];
         shvlohilo[j] = vlohilo[vBSind+j];
         shvhilolo[j] = vhilolo[vBSind+j];
         shvlololo[j] = vlololo[vBSind+j];
         odg_sqr(shvhihihi[j],shvlohihi[j],shvhilohi[j],shvlolohi[j],
                 shvhihilo[j],shvlohilo[j],shvhilolo[j],shvlololo[j],
                 &prdhihihi[j],&prdlohihi[j],&prdhilohi[j],&prdlolohi[j],
                 &prdhihilo[j],&prdlohilo[j],&prdhilolo[j],&prdlololo[j]);
      }
      __syncthreads();
      powTwo = 1;                          // sum reduction
      for(int k=0; k < BSLog2; k++)
      {
         if((j%(powTwo*2)) == 0)
            if(j+powTwo < BS)
               odg_inc(&prdhihihi[j],&prdlohihi[j],
                       &prdhilohi[j],&prdlolohi[j],
                       &prdhihilo[j],&prdlohilo[j],
                       &prdhilolo[j],&prdlololo[j],
                       prdhihihi[j+powTwo],prdlohihi[j+powTwo],
                       prdhilohi[j+powTwo],prdlolohi[j+powTwo],
                       prdhihilo[j+powTwo],prdlohilo[j+powTwo],
                       prdhilolo[j+powTwo],prdlololo[j+powTwo]);
         powTwo = powTwo*2;
         __syncthreads();
      }
      // thread 0 copies the sum of this round in sums[i], the others wait
      if(j == 0)
      {
         sumshihihi[i] = prdhihihi[0]; 
         sumslohihi[i] = prdlohihi[0]; 
         sumshilohi[i] = prdhilohi[0]; 
         sumslolohi[i] = prdlolohi[0]; 
         sumshihilo[i] = prdhihilo[0]; 
         sumslohilo[i] = prdlohilo[0]; 
         sumshilolo[i] = prdhilolo[0]; 
         sumslololo[i] = prdlololo[0]; 
      }
      __syncthreads();
      vBSind = vBSind + BS;
   }
   powTwo = 1;                          // sum reduction
   for(int k=0; k < rndLog2; k++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < rnd)
            odg_inc(&sumshihihi[j],&sumslohihi[j],
                    &sumshilohi[j],&sumslolohi[j],
                    &sumshihilo[j],&sumslohilo[j],
                    &sumshilolo[j],&sumslololo[j],
                    sumshihihi[j+powTwo],sumslohihi[j+powTwo],
                    sumshilohi[j+powTwo],sumslolohi[j+powTwo],
                    sumshihilo[j+powTwo],sumslohilo[j+powTwo],
                    sumshilolo[j+powTwo],sumslololo[j+powTwo]);
      powTwo = powTwo*2;
      __syncthreads();
   }
   if(j == 0)
      odg_sqrt(sumshihihi[0],sumslohihi[0],sumshilohi[0],sumslolohi[0],
               sumshihilo[0],sumslohilo[0],sumshilolo[0],sumslololo[0],
               &sumshihihi[0],&sumslohihi[0],&sumshilohi[0],&sumslolohi[0],
               &sumshihilo[0],&sumslohilo[0],&sumshilolo[0],&sumslololo[0]); 
   if(j == 0)
   {
      *normhihihi = sumshihihi[0];
      *normlohihi = sumslohihi[0];
      *normhilohi = sumshilohi[0];
      *normlolohi = sumslolohi[0];
      *normhihilo = sumshihilo[0];
      *normlohilo = sumslohilo[0];
      *normhilolo = sumshilolo[0];
      *normlololo = sumslololo[0];
   }
   __syncthreads();
   vBSind = 0;
   for(int i=0; i<rnd; i++)
   {
      if(vBSind+j < dim)
      {
         shvhihihi[j] = vhihihi[vBSind+j];   // read into shared memory
         shvlohihi[j] = vlohihi[vBSind+j];   // and normalize the vector
         shvhilohi[j] = vhilohi[vBSind+j];
         shvlolohi[j] = vlolohi[vBSind+j];
         shvhihilo[j] = vhihilo[vBSind+j];
         shvlohilo[j] = vlohilo[vBSind+j];
         shvhilolo[j] = vhilolo[vBSind+j];
         shvlololo[j] = vlololo[vBSind+j];
         odg_div(shvhihihi[j],shvlohihi[j],shvhilohi[j],shvlolohi[j],
                 shvhihilo[j],shvlohilo[j],shvhilolo[j],shvlololo[j],
                 sumshihihi[0],sumslohihi[0],sumshilohi[0],sumslolohi[0],
                 sumshihilo[0],sumslohilo[0],sumshilolo[0],sumslololo[0],
                 &vhihihi[vBSind+j],&vlohihi[vBSind+j],
                 &vhilohi[vBSind+j],&vlolohi[vBSind+j],
                 &vhihilo[vBSind+j],&vlohilo[vBSind+j],
                 &vhilolo[vBSind+j],&vlololo[vBSind+j]);
      }
      __syncthreads();
      vBSind = vBSind + BS;
   }
}
__global__ void large_sum_the_squares
 ( double *vhihihi, double *vlohihi, double *vhilohi, double *vlolohi,
   double *vhihilo, double *vlohilo, double *vhilolo, double *vlololo,
   int dim, double *sumshihihi, double *sumslohihi, double *sumshilohi,
   double *sumslolohi, double *sumshihilo, double *sumslohilo,
   double *sumshilolo, double *sumslololo, int BS, int BSLog2 )
{
   const int i = blockIdx.x;
   const int j = threadIdx.x;
   const int k = i*BS + j;

   __shared__ double shvhihihi[od_shmemsize];
   __shared__ double shvlohihi[od_shmemsize];
   __shared__ double shvhilohi[od_shmemsize];
   __shared__ double shvlolohi[od_shmemsize];
   __shared__ double shvhihilo[od_shmemsize];
   __shared__ double shvlohilo[od_shmemsize];
   __shared__ double shvhilolo[od_shmemsize];
   __shared__ double shvlololo[od_shmemsize];
   __shared__ double prdhihihi[od_shmemsize];
   __shared__ double prdlohihi[od_shmemsize];
   __shared__ double prdhilohi[od_shmemsize];
   __shared__ double prdlolohi[od_shmemsize];
   __shared__ double prdhihilo[od_shmemsize];
   __shared__ double prdlohilo[od_shmemsize];
   __shared__ double prdhilolo[od_shmemsize];
   __shared__ double prdlololo[od_shmemsize];

   shvhihihi[j] = vhihihi[k];
   shvlohihi[j] = vlohihi[k];
   shvhilohi[j] = vhilohi[k];
   shvlolohi[j] = vlolohi[k];
   shvhihilo[j] = vhihilo[k];
   shvlohilo[j] = vlohilo[k];
   shvhilolo[j] = vhilolo[k];
   shvlololo[j] = vlololo[k];
   odg_sqr(shvhihihi[j],shvlohihi[j],shvhilohi[j],shvlolohi[j],
           shvhihilo[j],shvlohilo[j],shvhilolo[j],shvlololo[j],
           &prdhihihi[j],&prdlohihi[j],&prdhilohi[j],&prdlolohi[j],
           &prdhihilo[j],&prdlohilo[j],&prdhilolo[j],&prdlololo[j]);

   __syncthreads();

   int powTwo = 1;                          // sum reduction
   for(int L=0; L < BSLog2; L++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < BS) 
            odg_inc(&prdhihihi[j],&prdlohihi[j],&prdhilohi[j],&prdlolohi[j],
                    &prdhihilo[j],&prdlohilo[j],&prdhilolo[j],&prdlololo[j],
                    prdhihihi[j+powTwo],prdlohihi[j+powTwo],
                    prdhilohi[j+powTwo],prdlolohi[j+powTwo],
                    prdhihilo[j+powTwo],prdlohilo[j+powTwo],
                    prdhilolo[j+powTwo],prdlololo[j+powTwo]);
      powTwo = powTwo*2;

      __syncthreads();
   }
   if(j == 0)                              // thread 0 writes the sum
   {
      sumshihihi[i] = prdhihihi[0];
      sumslohihi[i] = prdlohihi[0];
      sumshilohi[i] = prdhilohi[0];
      sumslolohi[i] = prdlolohi[0];
      sumshihilo[i] = prdhihilo[0];
      sumslohilo[i] = prdlohilo[0];
      sumshilolo[i] = prdhilolo[0];
      sumslololo[i] = prdlololo[0];
   }
}

__global__ void large_normalize_vector
 ( double *vhihihi, double *vlohihi, double *vhilohi, double *vlolohi,
   double *vhihilo, double *vlohilo, double *vhilolo, double *vlololo,
   int dim, double *sumshihihi, double *sumslohihi, double *sumshilohi,
   double *sumslolohi, double *sumshihilo, double *sumslohilo,
   double *sumshilolo, double *sumslololo, int nbsums, int nbsumsLog2,
   int BS, double *normhihihi, double *normlohihi, double *normhilohi,
   double *normlolohi, double *normhihilo, double *normlohilo,
   double *normhilolo, double *normlololo )
{
   const int i = blockIdx.x;
   const int j = threadIdx.x;
   const int k = i*BS + j;

   __shared__ double shvhihihi[od_shmemsize];
   __shared__ double shvlohihi[od_shmemsize];
   __shared__ double shvhilohi[od_shmemsize];
   __shared__ double shvlolohi[od_shmemsize];
   __shared__ double shvhihilo[od_shmemsize];
   __shared__ double shvlohilo[od_shmemsize];
   __shared__ double shvhilolo[od_shmemsize];
   __shared__ double shvlololo[od_shmemsize];

   if(j < nbsums)
   {
      shvhihihi[j] = sumshihihi[j];
      shvlohihi[j] = sumslohihi[j];
      shvhilohi[j] = sumshilohi[j];
      shvlolohi[j] = sumslolohi[j];
      shvhihilo[j] = sumshihilo[j];
      shvlohilo[j] = sumslohilo[j];
      shvhilolo[j] = sumshilolo[j];
      shvlololo[j] = sumslololo[j];
   }
   __syncthreads();

   int powTwo = 1;                          // sum reduction
   for(int L=0; L < nbsumsLog2; L++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < nbsums)
            odg_inc(&shvhihihi[j],&shvlohihi[j],&shvhilohi[j],&shvlolohi[j],
                    &shvhihilo[j],&shvlohilo[j],&shvhilolo[j],&shvlololo[j],
                    shvhihihi[j+powTwo],shvlohihi[j+powTwo],
                    shvhilohi[j+powTwo],shvlolohi[j+powTwo],
                    shvhihilo[j+powTwo],shvlohilo[j+powTwo],
                    shvhilolo[j+powTwo],shvlololo[j+powTwo]);
      powTwo = powTwo*2;

      __syncthreads();
   }
   __syncthreads();                    // every thread 0 of all blocks
   if(j == 0)                          // compute the 2-norm and assigns
   {                                   // to the output parameter
      odg_sqrt(shvhihihi[0],shvlohihi[0],shvhilohi[0],shvlolohi[0],
               shvhihilo[0],shvlohilo[0],shvhilolo[0],shvlololo[0],
               normhihihi,normlohihi,normhilohi,normlolohi,
               normhihilo,normlohilo,normhilolo,normlololo); 
   }
   __syncthreads();                    // to the output parameter

   if(k < dim)
   {
      shvhihihi[j] = vhihihi[k];
      shvlohihi[j] = vlohihi[k];
      shvhilohi[j] = vhilohi[k];
      shvlolohi[j] = vlolohi[k];
      shvhihilo[j] = vhihilo[k];
      shvlohilo[j] = vlohilo[k];
      shvhilolo[j] = vhilolo[k];
      shvlololo[j] = vlololo[k];
      odg_div(shvhihihi[j],shvlohihi[j],shvhilohi[j],shvlolohi[j],
              shvhihilo[j],shvlohilo[j],shvhilolo[j],shvlololo[j],
              *normhihihi,*normlohihi,*normhilohi,*normlolohi,
              *normhihilo,*normlohilo,*normhilolo,*normlololo,
              &shvhihihi[j],&shvlohihi[j],&shvhilohi[j],&shvlolohi[j],
              &shvhihilo[j],&shvlohilo[j],&shvhilolo[j],&shvlololo[j]);
      vhihihi[k] = shvhihihi[j];
      vlohihi[k] = shvlohihi[j];
      vhilohi[k] = shvhilohi[j];
      vlolohi[k] = shvlolohi[j];
      vhihilo[k] = shvhihilo[j];
      vlohilo[k] = shvlohilo[j];
      vhilolo[k] = shvhilolo[j];
      vlololo[k] = shvlololo[j];
   }
}

void GPU_norm
 ( double *vhihihi_h, double *vlohihi_h, double *vhilohi_h, double *vlolohi_h,
   double *vhihilo_h, double *vlohilo_h, double *vhilolo_h, double *vlololo_h,
   int dim, int freq, int BS, double *normhihihi, double *normlohihi,
   double *normhilohi, double *normlolohi, double *normhihilo,
   double *normlohilo, double *normhilolo, double *normlololo, int blocked )
{
   int BSLog2 = ceil(log2((double) BS)); // ceil for sum reduction

   double* vhihihi_d;                   // allocate for vector on device
   double* vlohihi_d;
   double* vhilohi_d;
   double* vlolohi_d;
   double* vhihilo_d;
   double* vlohilo_d;
   double* vhilolo_d;
   double* vlololo_d;
   size_t size = dim*sizeof(double);
   cudaMalloc((void**)&vhihihi_d,size);
   cudaMalloc((void**)&vlohihi_d,size);
   cudaMalloc((void**)&vhilohi_d,size);
   cudaMalloc((void**)&vlolohi_d,size);
   cudaMalloc((void**)&vhihilo_d,size);
   cudaMalloc((void**)&vlohilo_d,size);
   cudaMalloc((void**)&vhilolo_d,size);
   cudaMalloc((void**)&vlololo_d,size);
   cudaMemcpy(vhihihi_d,vhihihi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vlohihi_d,vlohihi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vhilohi_d,vhilohi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vlolohi_d,vlolohi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vhihilo_d,vhihilo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vlohilo_d,vlohilo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vhilolo_d,vhilolo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vlololo_d,vlololo_h,size,cudaMemcpyHostToDevice);
   double* normhihihi_d;
   double* normlohihi_d;
   double* normhilohi_d;
   double* normlolohi_d;
   double* normhihilo_d;
   double* normlohilo_d;
   double* normhilolo_d;
   double* normlololo_d;
   cudaMalloc((void**)&normhihihi_d,sizeof(double));
   cudaMalloc((void**)&normlohihi_d,sizeof(double));
   cudaMalloc((void**)&normhilohi_d,sizeof(double));
   cudaMalloc((void**)&normlolohi_d,sizeof(double));
   cudaMalloc((void**)&normhihilo_d,sizeof(double));
   cudaMalloc((void**)&normlohilo_d,sizeof(double));
   cudaMalloc((void**)&normhilolo_d,sizeof(double));
   cudaMalloc((void**)&normlololo_d,sizeof(double));

   if(dim == BS)
   {
      for(int i=0; i<freq; i++)
         small_normalize_vector<<<1,BS>>>
            (vhihihi_d,vlohihi_d,vhilohi_d,vlolohi_d,
             vhihilo_d,vlohilo_d,vhilolo_d,vlololo_d,dim,BSLog2,
             normhihihi_d,normlohihi_d,normhilohi_d,normlolohi_d,
             normhihilo_d,normlohilo_d,normhilolo_d,normlololo_d);
   }
   else if(blocked == 0)
   {
      const int rf = ceil(((double) dim)/BS);
      const int rfLog2 = ceil(log2((double) rf));
      for(int i=0; i<freq; i++)
         medium_normalize_vector<<<1,BS>>>
            (vhihihi_d,vlohihi_d,vhilohi_d,vlolohi_d,
             vhihilo_d,vlohilo_d,vhilolo_d,vlololo_d,
             dim,rf,rfLog2,BS,BSLog2,
             normhihihi_d,normlohihi_d,normhilohi_d,normlolohi_d,
             normhihilo_d,normlohilo_d,normhilolo_d,normlololo_d);
   }
   else
   {
      const int nblocks = dim/BS;
      const int nblocksLog2 = ceil(log2((double) nblocks));
      double* sumshihihi_d; // sums of squares for each block
      double* sumslohihi_d; // second highest parts of sums of squares
      double* sumshilohi_d; // third highest parts of sums of squares
      double* sumslolohi_d; // fourth highest parts of sums of squares
      double* sumshihilo_d; // fourth lowest parts of sums of squares
      double* sumslohilo_d; // third lowest parts of sums of squares
      double* sumshilolo_d; // second lowest parts of sums of squares
      double* sumslololo_d; // lowest parts of sums of squares
      size_t sums_size = nblocks*sizeof(double);
      cudaMalloc((void**)&sumshihihi_d,sums_size);
      cudaMalloc((void**)&sumslohihi_d,sums_size);
      cudaMalloc((void**)&sumshilohi_d,sums_size);
      cudaMalloc((void**)&sumslolohi_d,sums_size);
      cudaMalloc((void**)&sumshihilo_d,sums_size);
      cudaMalloc((void**)&sumslohilo_d,sums_size);
      cudaMalloc((void**)&sumshilolo_d,sums_size);
      cudaMalloc((void**)&sumslololo_d,sums_size);
      for(int i=0; i<freq; i++)
      {
         large_sum_the_squares<<<nblocks,BS>>>
            (vhihihi_d,vlohihi_d,vhilohi_d,vlolohi_d,
             vhihilo_d,vlohilo_d,vhilolo_d,vlololo_d,dim,
             sumshihihi_d,sumslohihi_d,sumshilohi_d,sumslolohi_d,
             sumshihilo_d,sumslohilo_d,sumshilolo_d,sumslololo_d,BS,BSLog2);
         large_normalize_vector<<<nblocks,BS>>>
            (vhihihi_d,vlohihi_d,vhilohi_d,vlolohi_d,
             vhihilo_d,vlohilo_d,vhilolo_d,vlololo_d,dim,
             sumshihihi_d,sumslohihi_d,sumshilohi_d,sumslolohi_d,
             sumshihilo_d,sumslohilo_d,sumshilolo_d,sumslololo_d,
             nblocks,nblocksLog2,BS,
             normhihihi_d,normlohihi_d,normhilohi_d,normlolohi_d,
             normhihilo_d,normlohilo_d,normhilolo_d,normlololo_d);
      }
   }
   cudaMemcpy(vhihihi_h,vhihihi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vlohihi_h,vlohihi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vhilohi_h,vhilohi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vlolohi_h,vlolohi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vhihilo_h,vhihilo_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vlohilo_h,vlohilo_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vhilolo_h,vhilolo_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vlololo_h,vlololo_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(normhihihi,normhihihi_d,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(normlohihi,normlohihi_d,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(normhilohi,normhilohi_d,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(normlolohi,normlolohi_d,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(normhihilo,normhihilo_d,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(normlohilo,normlohilo_d,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(normhilolo,normhilolo_d,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(normlololo,normlololo_d,sizeof(double),cudaMemcpyDeviceToHost);
}
