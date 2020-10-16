// Defines code of the functions in cmplx5_norm_kernels.h,
// to compute the 2-norm and normalize a complex vector,
// in penta double precision,
// for vectors of small, medium, and large size.

#include "double_double_gpufun.cu"
#include "penta_double_gpufun.cu"
#include "cmplx5_norm_kernels.h"

__global__ void small_normalize_vector
 ( double *vretb, double *vreix, double *vremi, double *vrerg, double *vrepk,
   double *vimtb, double *vimix, double *vimmi, double *vimrg, double *vimpk,
   int dim, int dimLog2,
   double *normtb, double *normix, double *normmi, double *normrg,
   double *normpk )
{
   int j = threadIdx.x;

   __shared__ double shvretb[pd_shmemsize];
   __shared__ double shvreix[pd_shmemsize];
   __shared__ double shvremi[pd_shmemsize];
   __shared__ double shvrerg[pd_shmemsize];
   __shared__ double shvrepk[pd_shmemsize];
   __shared__ double shvimtb[pd_shmemsize];
   __shared__ double shvimix[pd_shmemsize];
   __shared__ double shvimmi[pd_shmemsize];
   __shared__ double shvimrg[pd_shmemsize];
   __shared__ double shvimpk[pd_shmemsize];
   __shared__ double prdtb[pd_shmemsize];
   __shared__ double prdix[pd_shmemsize];
   __shared__ double prdmi[pd_shmemsize];
   __shared__ double prdrg[pd_shmemsize];
   __shared__ double prdpk[pd_shmemsize];
   __shared__ double sumtb[pd_shmemsize];
   __shared__ double sumix[pd_shmemsize];
   __shared__ double summi[pd_shmemsize];
   __shared__ double sumrg[pd_shmemsize];
   __shared__ double sumpk[pd_shmemsize];

   shvretb[j] = vretb[j]; // reading real parts into shared memory
   shvreix[j] = vreix[j];
   shvremi[j] = vremi[j];
   shvrerg[j] = vrerg[j];
   shvrepk[j] = vrepk[j];
   shvimtb[j] = vimtb[j]; // reading imaginary parts into shared memory
   shvimix[j] = vimix[j];
   shvimmi[j] = vimmi[j];
   shvimrg[j] = vimrg[j];
   shvimpk[j] = vimpk[j];

   pdg_sqr(shvretb[j],shvreix[j],shvremi[j],shvrerg[j],shvrepk[j],
            &sumtb[j], &sumix[j], &summi[j], &sumrg[j], &sumpk[j]);
   pdg_sqr(shvimtb[j],shvimix[j],shvimmi[j],shvimrg[j],shvimpk[j],
            &prdtb[j], &prdix[j], &prdmi[j], &prdrg[j], &prdpk[j]);
   pdg_inc(&sumtb[j],&sumix[j],&summi[j],&sumrg[j],&sumpk[j],
            prdtb[j], prdix[j], prdmi[j], prdrg[j], prdpk[j]);

   __syncthreads();
   int powTwo = 1;                          // sum reduction
   for(int k=0; k < dimLog2; k++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < dim)
            pdg_inc(&sumtb[j],&sumix[j],&summi[j],&sumrg[j],&sumpk[j],
                     sumtb[j+powTwo],sumix[j+powTwo],summi[j+powTwo],
                     sumrg[j+powTwo],sumpk[j+powTwo]);
      powTwo = powTwo*2;
      __syncthreads();
   }
   // thread 0 computes the sqrt of the inner product, others wait
   if(j == 0)
      pdg_sqrt( sumtb[0], sumix[0], summi[0], sumrg[0], sumpk[0],
               &sumtb[0],&sumix[0],&summi[0],&sumrg[0],&sumpk[0]); 
   if(j == 0)
   {
      *normtb = sumtb[0];
      *normix = sumix[0];
      *normmi = summi[0];
      *normrg = sumrg[0];
      *normpk = sumpk[0];
   }
   __syncthreads();
   pdg_div(shvretb[j],shvreix[j],shvremi[j],shvrerg[j],shvrepk[j],
             sumtb[0],  sumix[0],  summi[0],  sumrg[0],  sumpk[0],
            &vretb[j], &vreix[j], &vremi[j], &vrerg[j], &vrepk[j]);
   pdg_div(shvimtb[j],shvimix[j],shvimmi[j],shvimrg[j],shvimpk[j],
             sumtb[0],  sumix[0],  summi[0],  sumrg[0],  sumpk[0],
            &vimtb[j], &vimix[j], &vimmi[j], &vimrg[j], &vimpk[j]);
}

__global__ void medium_normalize_vector
 ( double *vretb, double *vreix, double *vremi, double *vrerg, double *vrepk,
   double *vimtb, double *vimix, double *vimmi, double *vimrg, double *vimpk,
   int dim, int rnd, int rndLog2, int BS, int BSLog2,
   double *normtb, double *normix, double *normrg, double *normmi,
   double *normpk )
{
   int j = threadIdx.x;
   int powTwo;
   int vBSind = 0;

   __shared__ double shvretb[pd_shmemsize];
   __shared__ double shvreix[pd_shmemsize];
   __shared__ double shvremi[pd_shmemsize];
   __shared__ double shvrerg[pd_shmemsize];
   __shared__ double shvrepk[pd_shmemsize];
   __shared__ double shvimtb[pd_shmemsize];
   __shared__ double shvimix[pd_shmemsize];
   __shared__ double shvimmi[pd_shmemsize];
   __shared__ double shvimrg[pd_shmemsize];
   __shared__ double shvimpk[pd_shmemsize];
   __shared__ double prdtb[pd_shmemsize];
   __shared__ double prdix[pd_shmemsize];
   __shared__ double prdmi[pd_shmemsize];
   __shared__ double prdrg[pd_shmemsize];
   __shared__ double prdpk[pd_shmemsize];
   __shared__ double acctb[pd_shmemsize];
   __shared__ double accix[pd_shmemsize];
   __shared__ double accmi[pd_shmemsize];
   __shared__ double accrg[pd_shmemsize];
   __shared__ double accpk[pd_shmemsize];
   __shared__ double sumstb[maxrounds];
   __shared__ double sumsix[maxrounds];
   __shared__ double sumsmi[maxrounds];
   __shared__ double sumsrg[maxrounds];
   __shared__ double sumspk[maxrounds];

   for(int i=0; i<rnd; i++)
   {
      if(vBSind+j >= dim)       // at last round some threads are excluded
      {
         acctb[j] = 0.0;
         accix[j] = 0.0;
         accmi[j] = 0.0;
         accrg[j] = 0.0;
         accpk[j] = 0.0;
      }
      else
      {
         shvretb[j] = vretb[vBSind+j];  // reading into shared memory
         shvreix[j] = vreix[vBSind+j];
         shvremi[j] = vremi[vBSind+j];
         shvrerg[j] = vrerg[vBSind+j];
         shvrepk[j] = vrepk[vBSind+j];
         shvimtb[j] = vimtb[vBSind+j]; 
         shvimix[j] = vimix[vBSind+j]; 
         shvimmi[j] = vimmi[vBSind+j]; 
         shvimrg[j] = vimrg[vBSind+j]; 
         shvimpk[j] = vimpk[vBSind+j]; 

         pdg_sqr(shvretb[j],shvreix[j],shvremi[j],shvrerg[j],shvrepk[j],
                  &acctb[j], &accix[j], &accmi[j], &accrg[j], &accpk[j]);
         pdg_sqr(shvimtb[j],shvimix[j],shvimmi[j],shvimrg[j],shvimpk[j],
                  &prdtb[j], &prdix[j], &prdmi[j], &prdrg[j], &prdpk[j]);
         pdg_inc(&acctb[j],&accix[j],&accmi[j],&accrg[j],&accpk[j],
                  prdtb[j], prdix[j], prdmi[j], prdrg[j], prdpk[j]);
      }
      __syncthreads();
      powTwo = 1;                          // sum reduction
      for(int k=0; k < BSLog2; k++)
      {
         if((j%(powTwo*2)) == 0)
            if(j+powTwo < BS)
               pdg_inc(&acctb[j],&accix[j],&accmi[j],&accrg[j],&accpk[j],
                       acctb[j+powTwo],accix[j+powTwo],accmi[j+powTwo],
                       accrg[j+powTwo],accpk[j+powTwo]);
         powTwo = powTwo*2;
         __syncthreads();
      }
      // thread 0 copies the sum of this round in sums[i], the others wait
      if(j == 0)
      {
         sumstb[i] = acctb[0]; 
         sumsix[i] = accix[0]; 
         sumsmi[i] = accmi[0]; 
         sumsrg[i] = accrg[0]; 
         sumspk[i] = accpk[0]; 
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
   {
      pdg_sqrt( sumstb[0], sumsix[0], sumsmi[0], sumsrg[0], sumspk[0],
               &sumstb[0],&sumsix[0],&sumsmi[0],&sumsrg[0],&sumspk[0]);
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
         shvretb[j] = vretb[vBSind+j];     // read into shared memory
         shvreix[j] = vreix[vBSind+j];
         shvremi[j] = vremi[vBSind+j];
         shvrerg[j] = vrerg[vBSind+j];
         shvrepk[j] = vrepk[vBSind+j];
         shvimtb[j] = vimtb[vBSind+j];
         shvimix[j] = vimix[vBSind+j];
         shvimmi[j] = vimmi[vBSind+j];
         shvimrg[j] = vimrg[vBSind+j];
         shvimpk[j] = vimpk[vBSind+j];
         // normalize vector
         pdg_div(shvretb[j],shvreix[j],shvremi[j],shvrerg[j],shvrepk[j],
                  sumstb[0], sumsix[0], sumsmi[0], sumsrg[0], sumspk[0],
                  &vretb[vBSind+j],&vreix[vBSind+j],&vremi[vBSind+j],
                  &vrerg[vBSind+j],&vrepk[vBSind+j]);
         pdg_div(shvimtb[j],shvimix[j],shvimmi[j],shvimrg[j],shvimpk[j],
                  sumstb[0], sumsix[0], sumsmi[0], sumsrg[0], sumspk[0],
                  &vimtb[vBSind+j],&vimix[vBSind+j],&vimmi[vBSind+j],
                  &vimrg[vBSind+j],&vimpk[vBSind+j]);
      }
      __syncthreads();
      vBSind = vBSind + BS;
   }
}

__global__ void large_sum_the_squares
 ( double *vretb, double *vreix, double *vremi, double *vrerg, double *vrepk,
   double *vimtb, double *vimix, double *vimmi, double *vimrg, double *vimpk,
   int dim,
   double *sumstb, double *sumsix, double *sumsmi, double *sumsrg,
   double *sumspk, int BS, int BSLog2 )
{
   const int i = blockIdx.x;
   const int j = threadIdx.x;
   const int k = i*BS + j;

   __shared__ double shvretb[pd_shmemsize];
   __shared__ double shvreix[pd_shmemsize];
   __shared__ double shvremi[pd_shmemsize];
   __shared__ double shvrerg[pd_shmemsize];
   __shared__ double shvrepk[pd_shmemsize];
   __shared__ double shvimtb[pd_shmemsize];
   __shared__ double shvimix[pd_shmemsize];
   __shared__ double shvimmi[pd_shmemsize];
   __shared__ double shvimrg[pd_shmemsize];
   __shared__ double shvimpk[pd_shmemsize];
   __shared__ double prdtb[pd_shmemsize];
   __shared__ double prdix[pd_shmemsize];
   __shared__ double prdmi[pd_shmemsize];
   __shared__ double prdrg[pd_shmemsize];
   __shared__ double prdpk[pd_shmemsize];
   __shared__ double acctb[pd_shmemsize];
   __shared__ double accix[pd_shmemsize];
   __shared__ double accmi[pd_shmemsize];
   __shared__ double accrg[pd_shmemsize];
   __shared__ double accpk[pd_shmemsize];

   shvretb[j] = vretb[k];
   shvreix[j] = vreix[k];
   shvremi[j] = vremi[k];
   shvrerg[j] = vrerg[k];
   shvrepk[j] = vrepk[k];
   shvimtb[j] = vimtb[k];
   shvimix[j] = vimix[k];
   shvimmi[j] = vimmi[k];
   shvimrg[j] = vimrg[k];
   shvimpk[j] = vimpk[k];

   pdg_sqr(shvretb[j],shvreix[j],shvremi[j],shvrerg[j],shvrepk[j],
            &acctb[j], &accix[j], &accmi[j], &accrg[j], &accpk[j]);
   pdg_sqr(shvimtb[j],shvimix[j],shvimmi[j],shvimrg[j],shvimpk[j],
            &prdtb[j], &prdix[j], &prdmi[j], &prdrg[j], &prdpk[j]);
   pdg_inc(&acctb[j],&accix[j],&accmi[j],&accrg[j],&accpk[j],
            prdtb[j], prdix[j], prdmi[j], prdrg[j], prdpk[j]);

   __syncthreads();

   int powTwo = 1;                          // sum reduction
   for(int L=0; L < BSLog2; L++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < BS)
            pdg_inc(&acctb[j],&accix[j],&accmi[j],&accrg[j],&accpk[j],
                     acctb[j+powTwo],accix[j+powTwo],accmi[j+powTwo],
                     accrg[j+powTwo],accpk[j+powTwo]);
      powTwo = powTwo*2;

      __syncthreads();
   }
   if(j == 0)                               // thread 0 writes the sum
   {
      sumstb[i] = acctb[0];
      sumsix[i] = accix[0];
      sumsmi[i] = accmi[0];
      sumsrg[i] = accrg[0];
      sumspk[i] = accpk[0];
   }
}

__global__ void large_normalize_vector
 ( double *vretb, double *vreix, double *vremi, double *vrerg, double *vrepk,
   double *vimtb, double *vimix, double *vimmi, double *vimrg, double *vimpk,
   int dim,
   double *sumstb, double *sumsix, double *sumsmi, double *sumsrg,
   double *sumspk, int nbsums, int nbsumsLog2, int BS,
   double *normtb, double *normix, double *normmi, double *normrg,
   double *normpk )
{
   const int i = blockIdx.x;
   const int j = threadIdx.x;
   const int k = i*BS + j;

   __shared__ double shvretb[pd_shmemsize];
   __shared__ double shvreix[pd_shmemsize];
   __shared__ double shvremi[pd_shmemsize];
   __shared__ double shvrerg[pd_shmemsize];
   __shared__ double shvrepk[pd_shmemsize];
   __shared__ double shvimtb[pd_shmemsize];
   __shared__ double shvimix[pd_shmemsize];
   __shared__ double shvimmi[pd_shmemsize];
   __shared__ double shvimrg[pd_shmemsize];
   __shared__ double shvimpk[pd_shmemsize];

   if(j < nbsums)
   {
      shvretb[j] = sumstb[j];
      shvreix[j] = sumsix[j];
      shvremi[j] = sumsmi[j];
      shvrerg[j] = sumsrg[j];
      shvrepk[j] = sumspk[j];
   }

   __syncthreads();

   int powTwo = 1;                          // sum reduction
   for(int L=0; L < nbsumsLog2; L++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < nbsums)
            pdg_inc(&shvretb[j],     &shvreix[j],       &shvremi[j],
                    &shvrerg[j],     &shvrepk[j],
                     shvretb[j+powTwo],shvreix[j+powTwo],shvremi[j+powTwo],
                     shvrerg[j+powTwo],shvrepk[j+powTwo]);
      powTwo = powTwo*2;

      __syncthreads();
   }
   __syncthreads();                      // every thread 0 of all blocks
   if(j == 0)                            // compute the 2-norm and assigns
   {                                     // to the output parameter
      pdg_sqrt(shvretb[0],shvreix[0],shvremi[0],shvrerg[0],shvrepk[0],
                normtb,    normix,    normmi,    normrg,    normpk);
   }
   __syncthreads();                    

   if(k < dim)
   {
      shvretb[j] = vretb[k];
      shvreix[j] = vreix[k];
      shvremi[j] = vremi[k];
      shvrerg[j] = vrerg[k];
      shvrepk[j] = vrepk[k];
      shvimtb[j] = vimtb[k];
      shvimix[j] = vimix[k];
      shvimmi[j] = vimmi[k];
      shvimrg[j] = vimrg[k];
      shvimpk[j] = vimpk[k];

      pdg_div( shvretb[j], shvreix[j], shvremi[j], shvrerg[j], shvrepk[j],
               *normtb,    *normix,    *normmi,    *normrg,    *normpk,
              &shvretb[j],&shvreix[j],&shvremi[j],&shvrerg[j],&shvrepk[j]);
      pdg_div( shvimtb[j], shvimix[j], shvimmi[j], shvimrg[j], shvimpk[j],
               *normtb,    *normix,    *normmi,    *normrg,    *normpk,
              &shvimtb[j],&shvimix[j],&shvimmi[j],&shvimrg[j],&shvimpk[j]);

      vretb[k] = shvretb[j];
      vreix[k] = shvreix[j];
      vremi[k] = shvremi[j];
      vrerg[k] = shvrerg[j];
      vrepk[k] = shvrepk[j];
      vimtb[k] = shvimtb[j];
      vimix[k] = shvimix[j];
      vimmi[k] = shvimmi[j];
      vimrg[k] = shvimrg[j];
      vimpk[k] = shvimpk[j];
   }
}

void GPU_norm
 ( double *vretb_h, double *vreix_h, double *vremi_h, double *vrerg_h,
   double *vrepk_h,
   double *vimtb_h, double *vimix_h, double *vimmi_h, double *vimrg_h,
   double *vimpk_h,
   int dim, int freq, int BS,
   double *normtb, double *normix, double *normmi, double *normrg,
   double *normpk, int blocked )
{
   int BSLog2 = ceil(log2((double) BS)); // ceil for sum reduction

   double* vretb_d;                 // highest real parts on device
   double* vreix_d;                 // 2nd highest real parts on device
   double* vremi_d;                 // middle real parts on device
   double* vrerg_d;                 // 2nd lowest real parts on device
   double* vrepk_d;                 // lowest real parts on device
   double* vimtb_d;                 // highest imaginary parts on device
   double* vimix_d;                 // 2nd highest imaginary parts on device
   double* vimmi_d;                 // middle imaginary parts on device
   double* vimrg_d;                 // 2nd lowest imaginary parts on device
   double* vimpk_d;                 // lowest imaginary parts on device
   size_t size = dim*sizeof(double);
   cudaMalloc((void**)&vretb_d,size);
   cudaMalloc((void**)&vreix_d,size);
   cudaMalloc((void**)&vremi_d,size);
   cudaMalloc((void**)&vrerg_d,size);
   cudaMalloc((void**)&vrepk_d,size);
   cudaMalloc((void**)&vimtb_d,size);
   cudaMalloc((void**)&vimix_d,size);
   cudaMalloc((void**)&vimmi_d,size);
   cudaMalloc((void**)&vimrg_d,size);
   cudaMalloc((void**)&vimpk_d,size);
   cudaMemcpy(vretb_d,vretb_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vreix_d,vreix_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vremi_d,vremi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vrerg_d,vrerg_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vrepk_d,vrepk_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vimtb_d,vimtb_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vimix_d,vimix_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vimmi_d,vimmi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vimrg_d,vimrg_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vimpk_d,vimpk_h,size,cudaMemcpyHostToDevice);
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
            (vretb_d,vreix_d,vremi_d,vrerg_d,vrepk_d,
             vimtb_d,vimix_d,vimmi_d,vimrg_d,vimpk_d,
             dim,BSLog2,normtb_d,normix_d,normmi_d,normrg_d,normpk_d);
   }
   else if(blocked == 0)
   {
      int rf = ceil(((double) dim)/BS);
      int rfLog2 = ceil(log2((double) rf));
      for(int i=0; i<freq; i++)
         medium_normalize_vector<<<1,BS>>>
            (vretb_d,vreix_d,vremi_d,vrerg_d,vrepk_d,
             vimtb_d,vimix_d,vimmi_d,vimrg_d,vimpk_d,dim,
             rf,rfLog2,BS,BSLog2,
             normtb_d,normix_d,normmi_d,normrg_d,normpk_d);
   }
   else
   {
      const int nblocks = dim/BS;
      const int nblocksLog2 = ceil(log2((double) nblocks));
      double* sumstb_d; // highest parts of sums of squares for each block
      double* sumsix_d; // 2nd highest parts of sums of squares
      double* sumsmi_d; // middle parts of sums of squares
      double* sumsrg_d; // 2nd lowest parts of sums of squares for each block
      double* sumspk_d; // lowest parts of sums of squares for each block
      size_t sums_size = nblocks*sizeof(double);
      cudaMalloc((void**)&sumstb_d,sums_size);
      cudaMalloc((void**)&sumsix_d,sums_size);
      cudaMalloc((void**)&sumsmi_d,sums_size);
      cudaMalloc((void**)&sumsrg_d,sums_size);
      cudaMalloc((void**)&sumspk_d,sums_size);
      for(int i=0; i<freq; i++)
      {
         large_sum_the_squares<<<nblocks,BS>>>
            (vretb_d,vreix_d,vremi_d,vrerg_d,vrepk_d,
             vimtb_d,vimix_d,vimmi_d,vimrg_d,vimpk_d,dim,
             sumstb_d,sumsix_d,sumsmi_d,sumsrg_d,sumspk_d,BS,BSLog2);
         large_normalize_vector<<<nblocks,BS>>>
            (vretb_d,vreix_d,vremi_d,vrerg_d,vrepk_d,
             vimtb_d,vimix_d,vimmi_d,vimrg_d,vimpk_d,dim,
             sumstb_d,sumsix_d,sumsmi_d,sumsrg_d,sumspk_d,
             nblocks,nblocksLog2,BS,
             normtb_d,normix_d,normmi_d,normrg_d,normpk_d);
      }
   }
   cudaMemcpy(vretb_h,vretb_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vreix_h,vreix_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vremi_h,vremi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vrerg_h,vrerg_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vrepk_h,vrepk_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vimtb_h,vimtb_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vimix_h,vimix_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vimmi_h,vimmi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vimrg_h,vimrg_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vimpk_h,vimpk_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(normtb,normtb_d,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(normix,normix_d,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(normmi,normmi_d,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(normrg,normrg_d,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(normpk,normpk_d,sizeof(double),cudaMemcpyDeviceToHost);
}
