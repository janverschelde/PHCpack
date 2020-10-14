// This file contains the definition for the functions in dbl5_norm_kernels.h,
// to compute the 2-norm and normalize a penta double precision vector,
// for small, medium, and large vectors.

#include "double_double_gpufun.cu"
#include "penta_double_gpufun.cu"
#include "dbl5_norm_kernels.h"

__global__ void small_normalize_vector
 ( double *vtb, double *vix, double *vmi, double *vrg, double *vpk, int dim,
   int dimLog2, double *normtb, double *normix, double *normmi,
   double *normrg, double *normpk )
{
   int j = threadIdx.x;
   __shared__ double shvtb[pd_shmemsize];
   __shared__ double shvix[pd_shmemsize];
   __shared__ double shvmi[pd_shmemsize];
   __shared__ double shvrg[pd_shmemsize];
   __shared__ double shvpk[pd_shmemsize];
   __shared__ double prdtb[pd_shmemsize];
   __shared__ double prdix[pd_shmemsize];
   __shared__ double prdmi[pd_shmemsize];
   __shared__ double prdrg[pd_shmemsize];
   __shared__ double prdpk[pd_shmemsize];
   shvtb[j] = vtb[j];    // reading of vector into shared memory
   shvix[j] = vix[j];
   shvmi[j] = vmi[j];
   shvrg[j] = vrg[j];
   shvpk[j] = vpk[j];
   pdg_sqr(shvtb[j],shvix[j],shvmi[j],shvrg[j],shvpk[j],
           &prdtb[j],&prdix[j],&prdmi[j],&prdrg[j],&prdpk[j]);
   __syncthreads();
   int powTwo = 1;                          // sum reduction
   for(int k=0; k < dimLog2; k++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < dim)
            pdg_inc(&prdtb[j],&prdix[j],&prdmi[j],&prdrg[j],&prdpk[j],
                    prdtb[j+powTwo],prdix[j+powTwo],prdmi[j+powTwo],
                    prdrg[j+powTwo],prdpk[j+powTwo]);
      powTwo = powTwo*2;
      __syncthreads();
   }
   // thread 0 computes the sqrt of the inner product, others wait
   if(j == 0)
      pdg_sqrt(prdtb[0],prdix[0],prdmi[0],prdrg[0],prdpk[0],
               &prdtb[0],&prdix[0],&prdmi[0],&prdrg[0],&prdpk[0]); 
   if(j == 0)
   {
      *normtb = prdtb[0];
      *normix = prdix[0];
      *normmi = prdmi[0];
      *normrg = prdrg[0];
      *normpk = prdpk[0];
   }
   __syncthreads();
   pdg_div(shvtb[j],shvix[j],shvmi[j],shvrg[j],shvpk[j],
           prdtb[0],prdix[0],prdmi[0],prdrg[0],prdpk[0],
           &vtb[j],&vix[j],&vmi[j],&vrg[j],&vpk[j]);
}

__global__ void medium_normalize_vector
 ( double *vtb, double *vix, double *vmi, double *vrg, double *vpk, int dim,
   int rnd, int rndLog2, int BS, int BSLog2, double *normtb, double *normix,
   double *normmi, double *normrg, double *normpk )
{
   int j = threadIdx.x;
   int powTwo;
   int vBSind = 0;

   __shared__ double shvtb[pd_shmemsize];
   __shared__ double shvix[pd_shmemsize];
   __shared__ double shvmi[pd_shmemsize];
   __shared__ double shvrg[pd_shmemsize];
   __shared__ double shvpk[pd_shmemsize];
   __shared__ double prdtb[pd_shmemsize];
   __shared__ double prdix[pd_shmemsize];
   __shared__ double prdmi[pd_shmemsize];
   __shared__ double prdrg[pd_shmemsize];
   __shared__ double prdpk[pd_shmemsize];
   __shared__ double sumstb[maxrounds];
   __shared__ double sumsix[maxrounds];
   __shared__ double sumsmi[maxrounds];
   __shared__ double sumsrg[maxrounds];
   __shared__ double sumspk[maxrounds];

   for(int i=0; i<rnd; i++)
   {
      if(vBSind+j >= dim)       // at last round some threads are excluded
      {
         prdtb[j] = 0.0;
         prdix[j] = 0.0;
         prdmi[j] = 0.0;
         prdrg[j] = 0.0;
         prdpk[j] = 0.0;
      }
      else
      {
         shvtb[j] = vtb[vBSind+j];  // reading of vector into shared memory
         shvix[j] = vix[vBSind+j];
         shvmi[j] = vmi[vBSind+j];
         shvrg[j] = vrg[vBSind+j];
         shvpk[j] = vpk[vBSind+j];
         pdg_sqr(shvtb[j],shvix[j],shvmi[j],shvrg[j],shvpk[j],
                 &prdtb[j],&prdix[j],&prdmi[j],&prdrg[j],&prdpk[j]);
      }
      __syncthreads();
      powTwo = 1;                          // sum reduction
      for(int k=0; k < BSLog2; k++)
      {
         if((j%(powTwo*2)) == 0)
            if(j+powTwo < BS)
               pdg_inc(&prdtb[j],&prdix[j],&prdmi[j],&prdrg[j],&prdpk[j],
                       prdtb[j+powTwo],prdix[j+powTwo],prdmi[j+powTwo],
                       prdrg[j+powTwo],prdpk[j+powTwo]);
         powTwo = powTwo*2;
         __syncthreads();
      }
      // thread 0 copies the sum of this round in sums[i], the others wait
      if(j == 0)
      {
         sumstb[i] = prdtb[0]; 
         sumsix[i] = prdix[0]; 
         sumsmi[i] = prdmi[0]; 
         sumsrg[i] = prdrg[0]; 
         sumspk[i] = prdpk[0]; 
      }
      __syncthreads();
      vBSind = vBSind + BS;
   }
   powTwo = 1;                          // sum reduction
   for(int k=0; k < rndLog2; k++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < rnd)
            pdg_inc(&sumstb[j],&sumsix[j],&sumsmi[j],&sumsrg[j],&sumspk[j],
                    sumstb[j+powTwo],sumsix[j+powTwo],sumsmi[j+powTwo],
                    sumsrg[j+powTwo],sumspk[j+powTwo]);
      powTwo = powTwo*2;
      __syncthreads();
   }
   if(j == 0)
      pdg_sqrt(sumstb[0],sumsix[0],sumsmi[0],sumsrg[0],sumspk[0],
               &sumstb[0],&sumsix[0],&sumsmi[0],&sumsrg[0],&sumspk[0]); 
   if(j == 0)
   {
      *normtb = sumstb[0];
      *normix = sumsix[0];
      *normmi = sumsmi[0];
      *normrg = sumsrg[0];
      *normpk = sumspk[0];
   }
   __syncthreads();
   vBSind = 0;
   for(int i=0; i<rnd; i++)
   {
      if(vBSind+j < dim)
      {
         shvtb[j] = vtb[vBSind+j];   // read into shared memory
         shvix[j] = vix[vBSind+j];   // and normalize the vector
         shvmi[j] = vmi[vBSind+j];
         shvrg[j] = vrg[vBSind+j];
         shvpk[j] = vpk[vBSind+j];
         pdg_div(shvtb[j],shvix[j],shvmi[j],shvrg[j],shvpk[j],
                 sumstb[0],sumsix[0],sumsmi[0],sumsrg[0],sumspk[0],
                 &vtb[vBSind+j],&vix[vBSind+j],&vmi[vBSind+j],
                 &vrg[vBSind+j],&vpk[vBSind+j]);
      }
      __syncthreads();
      vBSind = vBSind + BS;
   }
}

__global__ void large_sum_the_squares
 ( double *vtb, double *vix, double *vmi, double *vrg, double *vpk, int dim,
   double *sumstb, double *sumsix, double *sumsmi, double *sumsrg,
   double *sumspk, int BS, int BSLog2 )
{
   const int i = blockIdx.x;
   const int j = threadIdx.x;
   const int k = i*BS + j;

   __shared__ double shvtb[pd_shmemsize];
   __shared__ double shvix[pd_shmemsize];
   __shared__ double shvmi[pd_shmemsize];
   __shared__ double shvrg[pd_shmemsize];
   __shared__ double shvpk[pd_shmemsize];
   __shared__ double prdtb[pd_shmemsize];
   __shared__ double prdix[pd_shmemsize];
   __shared__ double prdmi[pd_shmemsize];
   __shared__ double prdrg[pd_shmemsize];
   __shared__ double prdpk[pd_shmemsize];

   shvtb[j] = vtb[k];
   shvix[j] = vix[k];
   shvmi[j] = vmi[k];
   shvrg[j] = vrg[k];
   shvpk[j] = vpk[k];
   pdg_sqr(shvtb[j],shvix[j],shvmi[j],shvrg[j],shvpk[j],
           &prdtb[j],&prdix[j],&prdmi[j],&prdrg[j],&prdpk[j]);

   __syncthreads();

   int powTwo = 1;                          // sum reduction
   for(int L=0; L < BSLog2; L++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < BS) 
            pdg_inc(&prdtb[j],&prdix[j],&prdmi[j],&prdrg[j],&prdpk[j],
                    prdtb[j+powTwo],prdix[j+powTwo],prdmi[j+powTwo],
                    prdrg[j+powTwo],prdpk[j+powTwo]);
      powTwo = powTwo*2;

      __syncthreads();
   }
   if(j == 0)                              // thread 0 writes the sum
   {
      sumstb[i] = prdtb[0];
      sumsix[i] = prdix[0];
      sumsmi[i] = prdmi[0];
      sumsrg[i] = prdrg[0];
      sumspk[i] = prdpk[0];
   }
}

__global__ void large_normalize_vector
 ( double *vtb, double *vix, double *vmi, double *vrg, double *vpk, int dim,
   double *sumstb, double *sumsix, double *sumsmi, double *sumsrg,
   double *sumspk, int nbsums, int nbsumsLog2, int BS, double *normtb,
   double *normix, double *normmi, double *normrg, double *normpk )
{
   const int i = blockIdx.x;
   const int j = threadIdx.x;
   const int k = i*BS + j;

   __shared__ double shvtb[pd_shmemsize];
   __shared__ double shvix[pd_shmemsize];
   __shared__ double shvmi[pd_shmemsize];
   __shared__ double shvrg[pd_shmemsize];
   __shared__ double shvpk[pd_shmemsize];

   if(j < nbsums)
   {
      shvtb[j] = sumstb[j];
      shvix[j] = sumsix[j];
      shvmi[j] = sumsmi[j];
      shvrg[j] = sumsrg[j];
      shvpk[j] = sumspk[j];
   }
   __syncthreads();

   int powTwo = 1;                          // sum reduction
   for(int L=0; L < nbsumsLog2; L++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < nbsums)
            pdg_inc(&shvtb[j],&shvix[j],&shvmi[j],&shvrg[j],&shvpk[j],
                    shvtb[j+powTwo],shvix[j+powTwo],shvmi[j+powTwo],
                    shvrg[j+powTwo],shvpk[j+powTwo]);
      powTwo = powTwo*2;

      __syncthreads();
   }
   __syncthreads();                    // every thread 0 of all blocks
   if(j == 0)                          // compute the 2-norm and assigns
   {                                   // to the output parameter
      pdg_sqrt(shvtb[0],shvix[0],shvmi[0],shvrg[0],shvpk[0],
               normtb,normix,normmi,normrg,normpk); 
   }
   __syncthreads();                    // to the output parameter

   if(k < dim)
   {
      shvtb[j] = vtb[k];
      shvix[j] = vix[k];
      shvmi[j] = vmi[k];
      shvrg[j] = vrg[k];
      shvpk[j] = vpk[k];
      pdg_div(shvtb[j],shvix[j],shvmi[j],shvrg[j],shvpk[j],
              *normtb,*normix,*normmi,*normrg,*normpk,
              &shvtb[j],&shvix[j],&shvmi[j],&shvrg[j],&shvpk[j]);
      vtb[k] = shvtb[j];
      vix[k] = shvix[j];
      vmi[k] = shvmi[j];
      vrg[k] = shvrg[j];
      vpk[k] = shvpk[j];
   }
}

void GPU_norm
 ( double *vtb_h, double *vix_h, double *vmi_h, double *vrg_h, double *vpk_h,
   int dim, int freq, int BS, double *normtb, double *normix, double *normmi,
   double *normrg, double *normpk, int blocked )
{
   int BSLog2 = ceil(log2((double) BS)); // ceil for sum reduction

   double* vtb_d;                   // allocate for vector on device
   double* vix_d;
   double* vmi_d;
   double* vrg_d;
   double* vpk_d;
   size_t size = dim*sizeof(double);
   cudaMalloc((void**)&vtb_d,size);
   cudaMalloc((void**)&vix_d,size);
   cudaMalloc((void**)&vmi_d,size);
   cudaMalloc((void**)&vrg_d,size);
   cudaMalloc((void**)&vpk_d,size);
   cudaMemcpy(vtb_d,vtb_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vix_d,vix_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vmi_d,vmi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vrg_d,vrg_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vpk_d,vpk_h,size,cudaMemcpyHostToDevice);
   double* normtb_d;
   double* normix_d;
   double* normmi_d;
   double* normrg_d;
   double* normpk_d;
   cudaMalloc((void**)&normtb_d,sizeof(double));
   cudaMalloc((void**)&normix_d,sizeof(double));
   cudaMalloc((void**)&normmi_d,sizeof(double));
   cudaMalloc((void**)&normrg_d,sizeof(double));
   cudaMalloc((void**)&normpk_d,sizeof(double));

   if(dim == BS)
   {
      for(int i=0; i<freq; i++)
         small_normalize_vector<<<1,BS>>>
            (vtb_d,vix_d,vmi_d,vrg_d,vpk_d,dim,BSLog2,
             normtb_d,normix_d,normmi_d,normrg_d,normpk_d);
   }
   else if(blocked == 0)
   {
      const int rf = ceil(((double) dim)/BS);
      const int rfLog2 = ceil(log2((double) rf));
      for(int i=0; i<freq; i++)
         medium_normalize_vector<<<1,BS>>>
            (vtb_d,vix_d,vmi_d,vrg_d,vpk_d,dim,rf,rfLog2,BS,BSLog2,
             normtb_d,normix_d,normmi_d,normrg_d,normpk_d);
   }
   else
   {
      const int nblocks = dim/BS;
      const int nblocksLog2 = ceil(log2((double) nblocks));
      double* sumstb_d; // sums of squares for each block
      double* sumsix_d; // second highest parts of sums of squares
      double* sumsmi_d; // middle parts of sums of squares
      double* sumsrg_d; // second lowest parts of sums of squares
      double* sumspk_d; // lowest parts of sums of squares
      size_t sums_size = nblocks*sizeof(double);
      cudaMalloc((void**)&sumstb_d,sums_size);
      cudaMalloc((void**)&sumsix_d,sums_size);
      cudaMalloc((void**)&sumsmi_d,sums_size);
      cudaMalloc((void**)&sumsrg_d,sums_size);
      cudaMalloc((void**)&sumspk_d,sums_size);
      for(int i=0; i<freq; i++)
      {
         large_sum_the_squares<<<nblocks,BS>>>
            (vtb_d,vix_d,vmi_d,vrg_d,vpk_d,dim,
             sumstb_d,sumsix_d,sumsmi_d,sumsrg_d,sumspk_d,BS,BSLog2);
         large_normalize_vector<<<nblocks,BS>>>
            (vtb_d,vix_d,vmi_d,vrg_d,vpk_d,dim,
             sumstb_d,sumsix_d,sumsmi_d,sumsrg_d,sumspk_d,
             nblocks,nblocksLog2,BS,
             normtb_d,normix_d,normmi_d,normrg_d,normpk_d);
      }
   }
   cudaMemcpy(vtb_h,vtb_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vix_h,vix_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vmi_h,vmi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vrg_h,vrg_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vpk_h,vpk_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(normtb,normtb_d,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(normix,normix_d,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(normmi,normmi_d,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(normrg,normrg_d,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(normpk,normpk_d,sizeof(double),cudaMemcpyDeviceToHost);
}
