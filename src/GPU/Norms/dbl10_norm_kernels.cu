// This file contains the definition for the functions in dbl10_norm_kernels.h,
// to compute the 2-norm and normalize a deca double precision vector,
// for small, medium, and large vectors.

#include "double_double_gpufun.cu"
#include "quad_double_gpufun.cu"
#include "octo_double_gpufun.cu"
#include "deca_double_gpufun.cu"
#include "dbl10_norm_kernels.h"

__global__ void small_normalize_vector
 ( double *vrtb, double *vrix, double *vrmi, double *vrrg, double *vrpk, 
   double *vltb, double *vlix, double *vlmi, double *vlrg, double *vlpk,
   int dim, int dimLog2, double *normrtb, double *normrix, double *normrmi,
   double *normrrg, double *normrpk, double *normltb, double *normlix,
   double *normlmi, double *normlrg, double *normlpk )
{
   int j = threadIdx.x;
   __shared__ double shvrtb[da_shmemsize];
   __shared__ double shvrix[da_shmemsize];
   __shared__ double shvrmi[da_shmemsize];
   __shared__ double shvrrg[da_shmemsize];
   __shared__ double shvrpk[da_shmemsize];
   __shared__ double shvltb[da_shmemsize];
   __shared__ double shvlix[da_shmemsize];
   __shared__ double shvlmi[da_shmemsize];
   __shared__ double shvlrg[da_shmemsize];
   __shared__ double shvlpk[da_shmemsize];
   __shared__ double prdrtb[da_shmemsize];
   __shared__ double prdrix[da_shmemsize];
   __shared__ double prdrmi[da_shmemsize];
   __shared__ double prdrrg[da_shmemsize];
   __shared__ double prdrpk[da_shmemsize];
   __shared__ double prdltb[da_shmemsize];
   __shared__ double prdlix[da_shmemsize];
   __shared__ double prdlmi[da_shmemsize];
   __shared__ double prdlrg[da_shmemsize];
   __shared__ double prdlpk[da_shmemsize];
   shvrtb[j] = vrtb[j];    // reading of vector into shared memory
   shvrix[j] = vrix[j];
   shvrmi[j] = vrmi[j];
   shvrrg[j] = vrrg[j];
   shvrpk[j] = vrpk[j];
   shvltb[j] = vltb[j];
   shvlix[j] = vlix[j];
   shvlmi[j] = vlmi[j];
   shvlrg[j] = vlrg[j];
   shvlpk[j] = vlpk[j];
   dag_sqr(shvrtb[j],shvrix[j],shvrmi[j],shvrrg[j],shvrpk[j],
           shvltb[j],shvlix[j],shvlmi[j],shvlrg[j],shvlpk[j],
           &prdrtb[j],&prdrix[j],&prdrmi[j],&prdrrg[j],&prdrpk[j],
           &prdltb[j],&prdlix[j],&prdlmi[j],&prdlrg[j],&prdlpk[j]);
   __syncthreads();
   int powTwo = 1;                          // sum reduction
   for(int k=0; k < dimLog2; k++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < dim)
            dag_inc(&prdrtb[j],&prdrix[j],&prdrmi[j],&prdrrg[j],&prdrpk[j],
                    &prdltb[j],&prdlix[j],&prdlmi[j],&prdlrg[j],&prdlpk[j],
                    prdrtb[j+powTwo],prdrix[j+powTwo],prdrmi[j+powTwo],
                    prdrrg[j+powTwo],prdrpk[j+powTwo],
                    prdltb[j+powTwo],prdlix[j+powTwo],prdlmi[j+powTwo],
                    prdlrg[j+powTwo],prdlpk[j+powTwo]);
      powTwo = powTwo*2;
      __syncthreads();
   }
   // thread 0 computes the sqrt of the inner prdauct, others wait
   if(j == 0)
      dag_sqrt(prdrtb[0],prdrix[0],prdrmi[0],prdrrg[0],prdrpk[0],
               prdltb[0],prdlix[0],prdlmi[0],prdlrg[0],prdlpk[0],
               &prdrtb[0],&prdrix[0],&prdrmi[0],&prdrrg[0],&prdrpk[0],
               &prdltb[0],&prdlix[0],&prdlmi[0],&prdlrg[0],&prdlpk[0]); 
   if(j == 0)
   {
      *normrtb = prdrtb[0];
      *normrix = prdrix[0];
      *normrmi = prdrmi[0];
      *normrrg = prdrrg[0];
      *normrpk = prdrpk[0];
      *normltb = prdltb[0];
      *normlix = prdlix[0];
      *normlmi = prdlmi[0];
      *normlrg = prdlrg[0];
      *normlpk = prdlpk[0];
   }
   __syncthreads();
   dag_div(shvrtb[j],shvrix[j],shvrmi[j],shvrrg[j],shvrpk[j],
           shvltb[j],shvlix[j],shvlmi[j],shvlrg[j],shvlpk[j],
           prdrtb[0],prdrix[0],prdrmi[0],prdrrg[0],prdrpk[0],
           prdltb[0],prdlix[0],prdlmi[0],prdlrg[0],prdlpk[0],
           &vrtb[j],&vrix[j],&vrmi[j],&vrrg[j],&vrpk[j],
           &vltb[j],&vlix[j],&vlmi[j],&vlrg[j],&vlpk[j]);
}

__global__ void medium_normalize_vector
 ( double *vrtb, double *vrix, double *vrmi, double *vrrg, double *vrpk,
   double *vltb, double *vlix, double *vlmi, double *vlrg, double *vlpk,
   int dim, int rnd, int rndLog2, int BS, int BSLog2,
   double *normrtb, double *normrix, double *normrmi, double *normrrg,
   double *normrpk, double *normltb, double *normlix, double *normlmi,
   double *normlrg, double *normlpk )
{
   int j = threadIdx.x;
   int powTwo;
   int vBSind = 0;

   __shared__ double shvrtb[da_shmemsize];
   __shared__ double shvrix[da_shmemsize];
   __shared__ double shvrmi[da_shmemsize];
   __shared__ double shvrrg[da_shmemsize];
   __shared__ double shvrpk[da_shmemsize];
   __shared__ double shvltb[da_shmemsize];
   __shared__ double shvlix[da_shmemsize];
   __shared__ double shvlmi[da_shmemsize];
   __shared__ double shvlrg[da_shmemsize];
   __shared__ double shvlpk[da_shmemsize];
   __shared__ double prdrtb[da_shmemsize];
   __shared__ double prdrix[da_shmemsize];
   __shared__ double prdrmi[da_shmemsize];
   __shared__ double prdrrg[da_shmemsize];
   __shared__ double prdrpk[da_shmemsize];
   __shared__ double prdltb[da_shmemsize];
   __shared__ double prdlix[da_shmemsize];
   __shared__ double prdlmi[da_shmemsize];
   __shared__ double prdlrg[da_shmemsize];
   __shared__ double prdlpk[da_shmemsize];
   __shared__ double sumsrtb[maxrounds];
   __shared__ double sumsrix[maxrounds];
   __shared__ double sumsrmi[maxrounds];
   __shared__ double sumsrrg[maxrounds];
   __shared__ double sumsrpk[maxrounds];
   __shared__ double sumsltb[maxrounds];
   __shared__ double sumslix[maxrounds];
   __shared__ double sumslmi[maxrounds];
   __shared__ double sumslrg[maxrounds];
   __shared__ double sumslpk[maxrounds];

   for(int i=0; i<rnd; i++)
   {
      if(vBSind+j >= dim)       // at last round some threads are excluded
      {
         prdrtb[j] = 0.0;
         prdrix[j] = 0.0;
         prdrmi[j] = 0.0;
         prdrrg[j] = 0.0;
         prdrpk[j] = 0.0;
         prdltb[j] = 0.0;
         prdlix[j] = 0.0;
         prdlmi[j] = 0.0;
         prdlrg[j] = 0.0;
         prdlpk[j] = 0.0;
      }
      else
      {
         shvrtb[j] = vrtb[vBSind+j];  // reading of vector into shared memory
         shvrix[j] = vrix[vBSind+j];
         shvrmi[j] = vrmi[vBSind+j];
         shvrrg[j] = vrrg[vBSind+j];
         shvrpk[j] = vrpk[vBSind+j];
         shvltb[j] = vltb[vBSind+j];
         shvlix[j] = vlix[vBSind+j];
         shvlmi[j] = vlmi[vBSind+j];
         shvlrg[j] = vlrg[vBSind+j];
         shvlpk[j] = vlpk[vBSind+j];
         dag_sqr(shvrtb[j],shvrix[j],shvrmi[j],shvrrg[j],shvrpk[j],
                 shvltb[j],shvlix[j],shvlmi[j],shvlrg[j],shvlpk[j],
                 &prdrtb[j],&prdrix[j],&prdrmi[j],&prdrrg[j],&prdrpk[j],
                 &prdltb[j],&prdlix[j],&prdlmi[j],&prdlrg[j],&prdlpk[j]);
      }
      __syncthreads();
      powTwo = 1;                          // sum reduction
      for(int k=0; k < BSLog2; k++)
      {
         if((j%(powTwo*2)) == 0)
            if(j+powTwo < BS)
               dag_inc(&prdrtb[j],&prdrix[j],&prdrmi[j],&prdrrg[j],&prdrpk[j],
                       &prdltb[j],&prdlix[j],&prdlmi[j],&prdlrg[j],&prdlpk[j],
                       prdrtb[j+powTwo],prdrix[j+powTwo],prdrmi[j+powTwo],
                       prdrrg[j+powTwo],prdrpk[j+powTwo],
                       prdltb[j+powTwo],prdlix[j+powTwo],prdlmi[j+powTwo],
                       prdlrg[j+powTwo],prdlpk[j+powTwo]);
         powTwo = powTwo*2;
         __syncthreads();
      }
      // thread 0 copies the sum of this round in sums[i], the others wait
      if(j == 0)
      {
         sumsrtb[i] = prdrtb[0]; 
         sumsrix[i] = prdrix[0]; 
         sumsrmi[i] = prdrmi[0]; 
         sumsrrg[i] = prdrrg[0]; 
         sumsrpk[i] = prdrpk[0]; 
         sumsltb[i] = prdltb[0]; 
         sumslix[i] = prdlix[0]; 
         sumslmi[i] = prdlmi[0]; 
         sumslrg[i] = prdlrg[0]; 
         sumslpk[i] = prdlpk[0]; 
      }
      __syncthreads();
      vBSind = vBSind + BS;
   }
   powTwo = 1;                          // sum reduction
   for(int k=0; k < rndLog2; k++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < rnd)
            dag_inc
               (&sumsrtb[j],&sumsrix[j],&sumsrmi[j],&sumsrrg[j],&sumsrpk[j],
                &sumsltb[j],&sumslix[j],&sumslmi[j],&sumslrg[j],&sumslpk[j],
                sumsrtb[j+powTwo],sumsrix[j+powTwo],sumsrmi[j+powTwo],
                sumsrrg[j+powTwo],sumsrpk[j+powTwo],
                sumsltb[j+powTwo],sumslix[j+powTwo],sumslmi[j+powTwo],
                sumslrg[j+powTwo],sumslpk[j+powTwo]);
      powTwo = powTwo*2;
      __syncthreads();
   }
   if(j == 0)
      dag_sqrt(sumsrtb[0],sumsrix[0],sumsrmi[0],sumsrrg[0],sumsrpk[0],
               sumsltb[0],sumslix[0],sumslmi[0],sumslrg[0],sumslpk[0],
               &sumsrtb[0],&sumsrix[0],&sumsrmi[0],&sumsrrg[0],&sumsrpk[0],
               &sumsltb[0],&sumslix[0],&sumslmi[0],&sumslrg[0],&sumslpk[0]); 
   if(j == 0)
   {
      *normrtb = sumsrtb[0];
      *normrix = sumsrix[0];
      *normrmi = sumsrmi[0];
      *normrrg = sumsrrg[0];
      *normrpk = sumsrpk[0];
      *normltb = sumsltb[0];
      *normlix = sumslix[0];
      *normlmi = sumslmi[0];
      *normlrg = sumslrg[0];
      *normlpk = sumslpk[0];
   }
   __syncthreads();
   vBSind = 0;
   for(int i=0; i<rnd; i++)
   {
      if(vBSind+j < dim)
      {
         shvrtb[j] = vrtb[vBSind+j];   // read into shared memory
         shvrix[j] = vrix[vBSind+j];   // and normalize the vector
         shvrmi[j] = vrmi[vBSind+j];
         shvrrg[j] = vrrg[vBSind+j];
         shvrpk[j] = vrpk[vBSind+j];
         shvltb[j] = vltb[vBSind+j];
         shvlix[j] = vlix[vBSind+j];
         shvlmi[j] = vlmi[vBSind+j];
         shvlrg[j] = vlrg[vBSind+j];
         shvlpk[j] = vlpk[vBSind+j];
         dag_div(shvrtb[j],shvrix[j],shvrmi[j],shvrrg[j],shvrpk[j],
                 shvltb[j],shvlix[j],shvlmi[j],shvlrg[j],shvlpk[j],
                 sumsrtb[0],sumsrix[0],sumsrmi[0],sumsrrg[0],sumsrpk[0],
                 sumsltb[0],sumslix[0],sumslmi[0],sumslrg[0],sumslpk[0],
                 &vrtb[vBSind+j],&vrix[vBSind+j],&vrmi[vBSind+j],
                 &vrrg[vBSind+j],&vrpk[vBSind+j],
                 &vltb[vBSind+j],&vlix[vBSind+j],&vlmi[vBSind+j],
                 &vlrg[vBSind+j],&vlpk[vBSind+j]);
      }
      __syncthreads();
      vBSind = vBSind + BS;
   }
}

__global__ void large_sum_the_squares
 ( double *vrtb, double *vrix, double *vrmi, double *vrrg, double *vrpk,
   double *vltb, double *vlix, double *vlmi, double *vlrg, double *vlpk,
   int dim, double *sumsrtb, double *sumsrix, double *sumsrmi,
   double *sumsrrg, double *sumsrpk, double *sumsltb, double *sumslix,
   double *sumslmi, double *sumslrg, double *sumslpk, int BS, int BSLog2 )
{
   const int i = blockIdx.x;
   const int j = threadIdx.x;
   const int k = i*BS + j;

   __shared__ double shvrtb[da_shmemsize];
   __shared__ double shvrix[da_shmemsize];
   __shared__ double shvrmi[da_shmemsize];
   __shared__ double shvrrg[da_shmemsize];
   __shared__ double shvrpk[da_shmemsize];
   __shared__ double shvltb[da_shmemsize];
   __shared__ double shvlix[da_shmemsize];
   __shared__ double shvlmi[da_shmemsize];
   __shared__ double shvlrg[da_shmemsize];
   __shared__ double shvlpk[da_shmemsize];
   __shared__ double prdrtb[da_shmemsize];
   __shared__ double prdrix[da_shmemsize];
   __shared__ double prdrmi[da_shmemsize];
   __shared__ double prdrrg[da_shmemsize];
   __shared__ double prdrpk[da_shmemsize];
   __shared__ double prdltb[da_shmemsize];
   __shared__ double prdlix[da_shmemsize];
   __shared__ double prdlmi[da_shmemsize];
   __shared__ double prdlrg[da_shmemsize];
   __shared__ double prdlpk[da_shmemsize];

   shvrtb[j] = vrtb[k];
   shvrix[j] = vrix[k];
   shvrmi[j] = vrmi[k];
   shvrrg[j] = vrrg[k];
   shvrpk[j] = vrpk[k];
   shvltb[j] = vltb[k];
   shvlix[j] = vlix[k];
   shvlmi[j] = vlmi[k];
   shvlrg[j] = vlrg[k];
   shvlpk[j] = vlpk[k];
   dag_sqr(shvrtb[j],shvrix[j],shvrmi[j],shvrrg[j],shvrpk[j],
           shvltb[j],shvlix[j],shvlmi[j],shvlrg[j],shvlpk[j],
           &prdrtb[j],&prdrix[j],&prdrmi[j],&prdrrg[j],&prdrpk[j],
           &prdltb[j],&prdlix[j],&prdlmi[j],&prdlrg[j],&prdlpk[j]);

   __syncthreads();

   int powTwo = 1;                          // sum reduction
   for(int L=0; L < BSLog2; L++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < BS) 
            dag_inc(&prdrtb[j],&prdrix[j],&prdrmi[j],&prdrrg[j],&prdrpk[j],
                    &prdltb[j],&prdlix[j],&prdlmi[j],&prdlrg[j],&prdlpk[j],
                    prdrtb[j+powTwo],prdrix[j+powTwo],prdrmi[j+powTwo],
                    prdrrg[j+powTwo],prdrpk[j+powTwo],
                    prdltb[j+powTwo],prdlix[j+powTwo],prdlmi[j+powTwo],
                    prdlrg[j+powTwo],prdlpk[j+powTwo]);
      powTwo = powTwo*2;

      __syncthreads();
   }
   if(j == 0)                              // thread 0 writes the sum
   {
      sumsrtb[i] = prdrtb[0];
      sumsrix[i] = prdrix[0];
      sumsrmi[i] = prdrmi[0];
      sumsrrg[i] = prdrrg[0];
      sumsrpk[i] = prdrpk[0];
      sumsltb[i] = prdltb[0];
      sumslix[i] = prdlix[0];
      sumslmi[i] = prdlmi[0];
      sumslrg[i] = prdlrg[0];
      sumslpk[i] = prdlpk[0];
   }
}

__global__ void large_normalize_vector
 ( double *vrtb, double *vrix, double *vrmi, double *vrrg, double *vrpk,
   double *vltb, double *vlix, double *vlmi, double *vlrg, double *vlpk,
   int dim, double *sumsrtb, double *sumsrix, double *sumsrmi,
   double *sumsrrg, double *sumsrpk, double *sumsltb, double *sumslix,
   double *sumslmi, double *sumslrg, double *sumslpk,
   int nbsums, int nbsumsLog2, int BS, double *normrtb, double *normrix,
   double *normrmi, double *normrrg, double *normrpk, double *normltb,
   double *normlix, double *normlmi, double *normlrg, double *normlpk )
{
   const int i = blockIdx.x;
   const int j = threadIdx.x;
   const int k = i*BS + j;

   __shared__ double shvrtb[da_shmemsize];
   __shared__ double shvrix[da_shmemsize];
   __shared__ double shvrmi[da_shmemsize];
   __shared__ double shvrrg[da_shmemsize];
   __shared__ double shvrpk[da_shmemsize];
   __shared__ double shvltb[da_shmemsize];
   __shared__ double shvlix[da_shmemsize];
   __shared__ double shvlmi[da_shmemsize];
   __shared__ double shvlrg[da_shmemsize];
   __shared__ double shvlpk[da_shmemsize];

   if(j < nbsums)
   {
      shvrtb[j] = sumsrtb[j];
      shvrix[j] = sumsrix[j];
      shvrmi[j] = sumsrmi[j];
      shvrrg[j] = sumsrrg[j];
      shvrpk[j] = sumsrpk[j];
      shvltb[j] = sumsltb[j];
      shvlix[j] = sumslix[j];
      shvlmi[j] = sumslmi[j];
      shvlrg[j] = sumslrg[j];
      shvlpk[j] = sumslpk[j];
   }
   __syncthreads();

   int powTwo = 1;                          // sum reduction
   for(int L=0; L < nbsumsLog2; L++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < nbsums)
            dag_inc(&shvrtb[j],&shvrix[j],&shvrmi[j],&shvrrg[j],&shvrpk[j],
                    &shvltb[j],&shvlix[j],&shvlmi[j],&shvlrg[j],&shvlpk[j],
                    shvrtb[j+powTwo],shvrix[j+powTwo],shvrmi[j+powTwo],
                    shvrrg[j+powTwo],shvrpk[j+powTwo],
                    shvltb[j+powTwo],shvlix[j+powTwo],shvlmi[j+powTwo],
                    shvlrg[j+powTwo],shvlpk[j+powTwo]);
      powTwo = powTwo*2;

      __syncthreads();
   }
   __syncthreads();                    // every thread 0 of all blocks
   if(j == 0)                          // compute the 2-norm and assigns
   {                                   // to the output parameter
      dag_sqrt(shvrtb[0],shvrix[0],shvrmi[0],shvrrg[0],shvrpk[0],
               shvltb[0],shvlix[0],shvlmi[0],shvlrg[0],shvlpk[0],
               normrtb,normrix,normrmi,normrrg,normrpk,
               normltb,normlix,normlmi,normlrg,normlpk); 
   }
   __syncthreads();                    // to the output parameter

   if(k < dim)
   {
      shvrtb[j] = vrtb[k];
      shvrix[j] = vrix[k];
      shvrmi[j] = vrmi[k];
      shvrrg[j] = vrrg[k];
      shvrpk[j] = vrpk[k];
      shvltb[j] = vltb[k];
      shvlix[j] = vlix[k];
      shvlmi[j] = vlmi[k];
      shvlrg[j] = vlrg[k];
      shvlpk[j] = vlpk[k];
      dag_div(shvrtb[j],shvrix[j],shvrmi[j],shvrrg[j],shvrpk[j],
              shvltb[j],shvlix[j],shvlmi[j],shvlrg[j],shvlpk[j],
              *normrtb,*normrix,*normrmi,*normrrg,*normrpk,
              *normltb,*normlix,*normlmi,*normlrg,*normlpk,
              &shvrtb[j],&shvrix[j],&shvrmi[j],&shvrrg[j],&shvrpk[j],
              &shvltb[j],&shvlix[j],&shvlmi[j],&shvlrg[j],&shvlpk[j]);
      vrtb[k] = shvrtb[j];
      vrix[k] = shvrix[j];
      vrmi[k] = shvrmi[j];
      vrrg[k] = shvrrg[j];
      vrpk[k] = shvrpk[j];
      vltb[k] = shvltb[j];
      vlix[k] = shvlix[j];
      vlmi[k] = shvlmi[j];
      vlrg[k] = shvlrg[j];
      vlpk[k] = shvlpk[j];
   }
}

void GPU_norm
 ( double *vrtb_h, double *vrix_h, double *vrmi_h, double *vrrg_h,
   double *vrpk_h, double *vltb_h, double *vlix_h, double *vlmi_h,
   double *vlrg_h, double *vlpk_h, int dim, int freq, int BS,
   double *normrtb, double *normrix, double *normrmi, double *normrrg,
   double *normrpk, double *normltb, double *normlix, double *normlmi,
   double *normlrg, double *normlpk, int blocked )
{
   int BSLog2 = ceil(log2((double) BS)); // ceil for sum reduction

   double* vrtb_d;                   // allocate for vector on device
   double* vrix_d;
   double* vrmi_d;
   double* vrrg_d;
   double* vrpk_d;
   double* vltb_d;
   double* vlix_d;
   double* vlmi_d;
   double* vlrg_d;
   double* vlpk_d;
   size_t size = dim*sizeof(double);
   cudaMalloc((void**)&vrtb_d,size);
   cudaMalloc((void**)&vrix_d,size);
   cudaMalloc((void**)&vrmi_d,size);
   cudaMalloc((void**)&vrrg_d,size);
   cudaMalloc((void**)&vrpk_d,size);
   cudaMalloc((void**)&vltb_d,size);
   cudaMalloc((void**)&vlix_d,size);
   cudaMalloc((void**)&vlmi_d,size);
   cudaMalloc((void**)&vlrg_d,size);
   cudaMalloc((void**)&vlpk_d,size);
   cudaMemcpy(vrtb_d,vrtb_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vrix_d,vrix_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vrmi_d,vrmi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vrrg_d,vrrg_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vrpk_d,vrpk_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vltb_d,vltb_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vlix_d,vlix_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vlmi_d,vlmi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vlrg_d,vlrg_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vlpk_d,vlpk_h,size,cudaMemcpyHostToDevice);
   double* normrtb_d;
   double* normrix_d;
   double* normrmi_d;
   double* normrrg_d;
   double* normrpk_d;
   double* normltb_d;
   double* normlix_d;
   double* normlmi_d;
   double* normlrg_d;
   double* normlpk_d;
   cudaMalloc((void**)&normrtb_d,sizeof(double));
   cudaMalloc((void**)&normrix_d,sizeof(double));
   cudaMalloc((void**)&normrmi_d,sizeof(double));
   cudaMalloc((void**)&normrrg_d,sizeof(double));
   cudaMalloc((void**)&normrpk_d,sizeof(double));
   cudaMalloc((void**)&normltb_d,sizeof(double));
   cudaMalloc((void**)&normlix_d,sizeof(double));
   cudaMalloc((void**)&normlmi_d,sizeof(double));
   cudaMalloc((void**)&normlrg_d,sizeof(double));
   cudaMalloc((void**)&normlpk_d,sizeof(double));

   if(dim == BS)
   {
      for(int i=0; i<freq; i++)
         small_normalize_vector<<<1,BS>>>
            (vrtb_d,vrix_d,vrmi_d,vrrg_d,vrpk_d,
             vltb_d,vlix_d,vlmi_d,vlrg_d,vlpk_d,dim,BSLog2,
             normrtb_d,normrix_d,normrmi_d,normrrg_d,normrpk_d,
             normltb_d,normlix_d,normlmi_d,normlrg_d,normlpk_d);
   }
   else if(blocked == 0)
   {
      const int rf = ceil(((double) dim)/BS);
      const int rfLog2 = ceil(log2((double) rf));
      for(int i=0; i<freq; i++)
         medium_normalize_vector<<<1,BS>>>
            (vrtb_d,vrix_d,vrmi_d,vrrg_d,vrpk_d,
             vltb_d,vlix_d,vlmi_d,vlrg_d,vlpk_d,dim,rf,rfLog2,BS,BSLog2,
             normrtb_d,normrix_d,normrmi_d,normrrg_d,normrpk_d,
             normltb_d,normlix_d,normlmi_d,normlrg_d,normlpk_d);
   }
   else
   {
      const int nblocks = dim/BS;
      const int nblocksLog2 = ceil(log2((double) nblocks));
      double* sumsrtb_d; // sums of squares for each block
      double* sumsrix_d; // second highest parts of sums of squares
      double* sumsrmi_d; // third highest parts of sums of squares
      double* sumsrrg_d; // fourth highest parts of sums of squares
      double* sumsrpk_d; // fifth highest parts of sums of squares
      double* sumsltb_d; // fifth lowest parts of sum of squares
      double* sumslix_d; // fourth lowest parts of sums of squares
      double* sumslmi_d; // third lowest parts of sums of squares
      double* sumslrg_d; // second lowest parts of sums of squares
      double* sumslpk_d; // lowest parts of sums of squares
      size_t sums_size = nblocks*sizeof(double);
      cudaMalloc((void**)&sumsrtb_d,sums_size);
      cudaMalloc((void**)&sumsrix_d,sums_size);
      cudaMalloc((void**)&sumsrmi_d,sums_size);
      cudaMalloc((void**)&sumsrrg_d,sums_size);
      cudaMalloc((void**)&sumsrpk_d,sums_size);
      cudaMalloc((void**)&sumsltb_d,sums_size);
      cudaMalloc((void**)&sumslix_d,sums_size);
      cudaMalloc((void**)&sumslmi_d,sums_size);
      cudaMalloc((void**)&sumslrg_d,sums_size);
      cudaMalloc((void**)&sumslpk_d,sums_size);
      for(int i=0; i<freq; i++)
      {
         large_sum_the_squares<<<nblocks,BS>>>
            (vrtb_d,vrix_d,vrmi_d,vrrg_d,vrpk_d,
             vltb_d,vlix_d,vlmi_d,vlrg_d,vlpk_d,dim,
             sumsrtb_d,sumsrix_d,sumsrmi_d,sumsrrg_d,sumsrpk_d,
             sumsltb_d,sumslix_d,sumslmi_d,sumslrg_d,sumslpk_d,BS,BSLog2);
         large_normalize_vector<<<nblocks,BS>>>
            (vrtb_d,vrix_d,vrmi_d,vrrg_d,vrpk_d,
             vltb_d,vlix_d,vlmi_d,vlrg_d,vlpk_d,dim,
             sumsrtb_d,sumsrix_d,sumsrmi_d,sumsrrg_d,sumsrpk_d,
             sumsltb_d,sumslix_d,sumslmi_d,sumslrg_d,sumslpk_d,
             nblocks,nblocksLog2,BS,
             normrtb_d,normrix_d,normrmi_d,normrrg_d,normrpk_d,
             normltb_d,normlix_d,normlmi_d,normlrg_d,normlpk_d);
      }
   }
   cudaMemcpy(vrtb_h,vrtb_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vrix_h,vrix_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vrmi_h,vrmi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vrrg_h,vrrg_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vrpk_h,vrpk_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vltb_h,vltb_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vlix_h,vlix_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vlmi_h,vlmi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vlrg_h,vlrg_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vlpk_h,vlpk_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(normrtb,normrtb_d,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(normrix,normrix_d,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(normrmi,normrmi_d,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(normrrg,normrrg_d,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(normrpk,normrpk_d,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(normltb,normltb_d,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(normlix,normlix_d,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(normlmi,normlmi_d,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(normlrg,normlrg_d,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(normlpk,normlpk_d,sizeof(double),cudaMemcpyDeviceToHost);
}
