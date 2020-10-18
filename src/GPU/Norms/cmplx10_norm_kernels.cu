// Defines code of the functions in cmplx10_norm_kernels.h,
// to compute the 2-norm and normalize a complex vector,
// in deca double precision,
// for vectors of small, medium, and large size.

#include "double_double_gpufun.cu"
#include "quad_double_gpufun.cu"
#include "octo_double_gpufun.cu"
#include "deca_double_gpufun.cu"
#include "cmplx10_norm_kernels.h"

__global__ void small_normalize_vector
 ( double *vrertb, double *vrerix, double *vrermi, double *vrerrg,
   double *vrerpk, double *vreltb, double *vrelix, double *vrelmi,
   double *vrelrg, double *vrelpk,
   double *vimrtb, double *vimrix, double *vimrmi, double *vimrrg,
   double *vimrpk, double *vimltb, double *vimlix, double *vimlmi,
   double *vimlrg, double *vimlpk, int dim, int dimLog2,
   double *normrtb, double *normrix, double *normrmi, double *normrrg,
   double *normrpk, double *normltb, double *normlix, double *normlmi,
   double *normlrg, double *normlpk )
{
   int j = threadIdx.x;

   __shared__ double shvrertb[da_shmemsize];
   __shared__ double shvrerix[da_shmemsize];
   __shared__ double shvrermi[da_shmemsize];
   __shared__ double shvrerrg[da_shmemsize];
   __shared__ double shvrerpk[da_shmemsize];
   __shared__ double shvreltb[da_shmemsize];
   __shared__ double shvrelix[da_shmemsize];
   __shared__ double shvrelmi[da_shmemsize];
   __shared__ double shvrelrg[da_shmemsize];
   __shared__ double shvrelpk[da_shmemsize];
   __shared__ double shvimrtb[da_shmemsize];
   __shared__ double shvimrix[da_shmemsize];
   __shared__ double shvimrmi[da_shmemsize];
   __shared__ double shvimrrg[da_shmemsize];
   __shared__ double shvimrpk[da_shmemsize];
   __shared__ double shvimltb[da_shmemsize];
   __shared__ double shvimlix[da_shmemsize];
   __shared__ double shvimlmi[da_shmemsize];
   __shared__ double shvimlrg[da_shmemsize];
   __shared__ double shvimlpk[da_shmemsize];
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
   __shared__ double sumrtb[da_shmemsize];
   __shared__ double sumrix[da_shmemsize];
   __shared__ double sumrmi[da_shmemsize];
   __shared__ double sumrrg[da_shmemsize];
   __shared__ double sumrpk[da_shmemsize];
   __shared__ double sumltb[da_shmemsize];
   __shared__ double sumlix[da_shmemsize];
   __shared__ double sumlmi[da_shmemsize];
   __shared__ double sumlrg[da_shmemsize];
   __shared__ double sumlpk[da_shmemsize];

   shvrertb[j] = vrertb[j]; // reading real parts into shared memory
   shvrerix[j] = vrerix[j];
   shvrermi[j] = vrermi[j];
   shvrerrg[j] = vrerrg[j];
   shvrerpk[j] = vrerpk[j];
   shvreltb[j] = vreltb[j];
   shvrelix[j] = vrelix[j];
   shvrelmi[j] = vrelmi[j];
   shvrelrg[j] = vrelrg[j];
   shvrelpk[j] = vrelpk[j];
   shvimrtb[j] = vimrtb[j]; // reading imaginary parts into shared memory
   shvimrix[j] = vimrix[j];
   shvimrmi[j] = vimrmi[j];
   shvimrrg[j] = vimrrg[j];
   shvimrpk[j] = vimrpk[j];
   shvimltb[j] = vimltb[j];
   shvimlix[j] = vimlix[j];
   shvimlmi[j] = vimlmi[j];
   shvimlrg[j] = vimlrg[j];
   shvimlpk[j] = vimlpk[j];

   dag_sqr(shvrertb[j],shvrerix[j],shvrermi[j],shvrerrg[j],shvrerpk[j],
           shvreltb[j],shvrelix[j],shvrelmi[j],shvrelrg[j],shvrelpk[j],
            &sumrtb[j], &sumrix[j], &sumrmi[j], &sumrrg[j], &sumrpk[j],
            &sumltb[j], &sumlix[j], &sumlmi[j], &sumlrg[j], &sumlpk[j]);
   dag_sqr(shvimrtb[j],shvimrix[j],shvimrmi[j],shvimrrg[j],shvimrpk[j],
           shvimltb[j],shvimlix[j],shvimlmi[j],shvimlrg[j],shvimlpk[j],
            &prdrtb[j], &prdrix[j], &prdrmi[j], &prdrrg[j], &prdrpk[j],
            &prdltb[j], &prdlix[j], &prdlmi[j], &prdlrg[j], &prdlpk[j]);
   dag_inc(&sumrtb[j],&sumrix[j],&sumrmi[j],&sumrrg[j],&sumrpk[j],
           &sumltb[j],&sumlix[j],&sumlmi[j],&sumlrg[j],&sumlpk[j],
            prdrtb[j], prdrix[j], prdrmi[j], prdrrg[j], prdrpk[j],
            prdltb[j], prdlix[j], prdlmi[j], prdlrg[j], prdlpk[j]);

   __syncthreads();
   int powTwo = 1;                          // sum reduction
   for(int k=0; k < dimLog2; k++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < dim)
            dag_inc(&sumrtb[j],&sumrix[j],&sumrmi[j],&sumrrg[j],&sumrpk[j],
                    &sumltb[j],&sumlix[j],&sumlmi[j],&sumlrg[j],&sumlpk[j],
                     sumrtb[j+powTwo],sumrix[j+powTwo],sumrmi[j+powTwo],
                     sumrrg[j+powTwo],sumrpk[j+powTwo],
                     sumltb[j+powTwo],sumlix[j+powTwo],sumlmi[j+powTwo],
                     sumlrg[j+powTwo],sumlpk[j+powTwo]);
      powTwo = powTwo*2;
      __syncthreads();
   }
   // thread 0 computes the sqrt of the inner product, others wait
   if(j == 0)
   {
      dag_sqrt( sumrtb[0], sumrix[0], sumrmi[0], sumrrg[0], sumrpk[0],
                sumltb[0], sumlix[0], sumlmi[0], sumlrg[0], sumlpk[0],
               &sumrtb[0],&sumrix[0],&sumrmi[0],&sumrrg[0],&sumrpk[0],
               &sumltb[0],&sumlix[0],&sumlmi[0],&sumlrg[0],&sumlpk[0]); 
      *normrtb = sumrtb[0];
      *normrix = sumrix[0];
      *normrmi = sumrmi[0];
      *normrrg = sumrrg[0];
      *normrpk = sumrpk[0];
      *normltb = sumltb[0];
      *normlix = sumlix[0];
      *normlmi = sumlmi[0];
      *normlrg = sumlrg[0];
      *normlpk = sumlpk[0];
   }
   __syncthreads();
   dag_div(shvrertb[j],shvrerix[j],shvrermi[j],shvrerrg[j],shvrerpk[j],
           shvreltb[j],shvrelix[j],shvrelmi[j],shvrelrg[j],shvrelpk[j],
             sumrtb[0],  sumrix[0],  sumrmi[0],  sumrrg[0],  sumrpk[0],
             sumltb[0],  sumlix[0],  sumlmi[0],  sumlrg[0],  sumlpk[0],
            &vrertb[j], &vrerix[j], &vrermi[j], &vrerrg[j], &vrerpk[j],
            &vreltb[j], &vrelix[j], &vrelmi[j], &vrelrg[j], &vrelpk[j]);
   dag_div(shvimrtb[j],shvimrix[j],shvimrmi[j],shvimrrg[j],shvimrpk[j],
           shvimltb[j],shvimlix[j],shvimlmi[j],shvimlrg[j],shvimlpk[j],
             sumrtb[0],  sumrix[0],  sumrmi[0],  sumrrg[0],  sumrpk[0],
             sumltb[0],  sumlix[0],  sumlmi[0],  sumlrg[0],  sumlpk[0],
            &vimrtb[j], &vimrix[j], &vimrmi[j], &vimrrg[j], &vimrpk[j],
            &vimltb[j], &vimlix[j], &vimlmi[j], &vimlrg[j], &vimlpk[j]);
}

__global__ void medium_normalize_vector
 ( double *vrertb, double *vrerix, double *vrermi, double *vrerrg,
   double *vrerpk, double *vreltb, double *vrelix, double *vrelmi,
   double *vrelrg, double *vrelpk,
   double *vimrtb, double *vimrix, double *vimrmi, double *vimrrg,
   double *vimrpk, double *vimltb, double *vimlix, double *vimlmi,
   double *vimlrg, double *vimlpk,
   int dim, int rnd, int rndLog2, int BS, int BSLog2,
   double *normrtb, double *normrix, double *normrmi, double *normrrg,
   double *normrpk, double *normltb, double *normlix, double *normlmi,
   double *normlrg, double *normlpk )
{
   int j = threadIdx.x;
   int powTwo;
   int vBSind = 0;

   __shared__ double shvrertb[da_shmemsize];
   __shared__ double shvrerix[da_shmemsize];
   __shared__ double shvrermi[da_shmemsize];
   __shared__ double shvrerrg[da_shmemsize];
   __shared__ double shvrerpk[da_shmemsize];
   __shared__ double shvreltb[da_shmemsize];
   __shared__ double shvrelix[da_shmemsize];
   __shared__ double shvrelmi[da_shmemsize];
   __shared__ double shvrelrg[da_shmemsize];
   __shared__ double shvrelpk[da_shmemsize];
   __shared__ double shvimrtb[da_shmemsize];
   __shared__ double shvimrix[da_shmemsize];
   __shared__ double shvimrmi[da_shmemsize];
   __shared__ double shvimrrg[da_shmemsize];
   __shared__ double shvimrpk[da_shmemsize];
   __shared__ double shvimltb[da_shmemsize];
   __shared__ double shvimlix[da_shmemsize];
   __shared__ double shvimlmi[da_shmemsize];
   __shared__ double shvimlrg[da_shmemsize];
   __shared__ double shvimlpk[da_shmemsize];
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
   __shared__ double accrtb[da_shmemsize];
   __shared__ double accrix[da_shmemsize];
   __shared__ double accrmi[da_shmemsize];
   __shared__ double accrrg[da_shmemsize];
   __shared__ double accrpk[da_shmemsize];
   __shared__ double accltb[da_shmemsize];
   __shared__ double acclix[da_shmemsize];
   __shared__ double acclmi[da_shmemsize];
   __shared__ double acclrg[da_shmemsize];
   __shared__ double acclpk[da_shmemsize];
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
         accrtb[j] = 0.0;
         accrix[j] = 0.0;
         accrmi[j] = 0.0;
         accrrg[j] = 0.0;
         accrpk[j] = 0.0;
         accltb[j] = 0.0;
         acclix[j] = 0.0;
         acclmi[j] = 0.0;
         acclrg[j] = 0.0;
         acclpk[j] = 0.0;
      }
      else
      {
         shvrertb[j] = vrertb[vBSind+j];  // reading into shared memory
         shvrerix[j] = vrerix[vBSind+j];
         shvrermi[j] = vrermi[vBSind+j];
         shvrerrg[j] = vrerrg[vBSind+j];
         shvrerpk[j] = vrerpk[vBSind+j];
         shvreltb[j] = vreltb[vBSind+j];
         shvrelix[j] = vrelix[vBSind+j];
         shvrelmi[j] = vrelmi[vBSind+j];
         shvrelrg[j] = vrelrg[vBSind+j];
         shvrelpk[j] = vrelpk[vBSind+j];
         shvimrtb[j] = vimrtb[vBSind+j]; 
         shvimrix[j] = vimrix[vBSind+j]; 
         shvimrmi[j] = vimrmi[vBSind+j]; 
         shvimrrg[j] = vimrrg[vBSind+j]; 
         shvimrpk[j] = vimrpk[vBSind+j]; 
         shvimltb[j] = vimltb[vBSind+j]; 
         shvimlix[j] = vimlix[vBSind+j]; 
         shvimlmi[j] = vimlmi[vBSind+j]; 
         shvimlrg[j] = vimlrg[vBSind+j]; 
         shvimlpk[j] = vimlpk[vBSind+j]; 

         dag_sqr(shvrertb[j],shvrerix[j],shvrermi[j],shvrerrg[j],shvrerpk[j],
                 shvreltb[j],shvrelix[j],shvrelmi[j],shvrelrg[j],shvrelpk[j],
                  &accrtb[j], &accrix[j], &accrmi[j], &accrrg[j], &accrpk[j],
                  &accltb[j], &acclix[j], &acclmi[j], &acclrg[j], &acclpk[j]);
         dag_sqr(shvimrtb[j],shvimrix[j],shvimrmi[j],shvimrrg[j],shvimrpk[j],
                 shvimltb[j],shvimlix[j],shvimlmi[j],shvimlrg[j],shvimlpk[j],
                  &prdrtb[j], &prdrix[j], &prdrmi[j], &prdrrg[j], &prdrpk[j],
                  &prdltb[j], &prdlix[j], &prdlmi[j], &prdlrg[j], &prdlpk[j]);
         dag_inc(&accrtb[j],&accrix[j],&accrmi[j],&accrrg[j],&accrpk[j],
                 &accltb[j],&acclix[j],&acclmi[j],&acclrg[j],&acclpk[j],
                  prdrtb[j], prdrix[j], prdrmi[j], prdrrg[j], prdrpk[j],
                  prdltb[j], prdlix[j], prdlmi[j], prdlrg[j], prdlpk[j]);
      }
      __syncthreads();
      powTwo = 1;                          // sum reduction
      for(int k=0; k < BSLog2; k++)
      {
         if((j%(powTwo*2)) == 0)
            if(j+powTwo < BS)
               dag_inc(&accrtb[j],&accrix[j],&accrmi[j],&accrrg[j],&accrpk[j],
                       &accltb[j],&acclix[j],&acclmi[j],&acclrg[j],&acclpk[j],
                        accrtb[j+powTwo],accrix[j+powTwo],accrmi[j+powTwo],
                        accrrg[j+powTwo],accrpk[j+powTwo],
                        accltb[j+powTwo],acclix[j+powTwo],acclmi[j+powTwo],
                        acclrg[j+powTwo],acclpk[j+powTwo]);
         powTwo = powTwo*2;
         __syncthreads();
      }
      // thread 0 copies the sum of this round in sums[i], the others wait
      if(j == 0)
      {
         sumsrtb[i] = accrtb[0]; 
         sumsrix[i] = accrix[0]; 
         sumsrmi[i] = accrmi[0]; 
         sumsrrg[i] = accrrg[0]; 
         sumsrpk[i] = accrpk[0]; 
         sumsltb[i] = accltb[0]; 
         sumslix[i] = acclix[0]; 
         sumslmi[i] = acclmi[0]; 
         sumslrg[i] = acclrg[0]; 
         sumslpk[i] = acclpk[0]; 
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
   {
      dag_sqrt( sumsrtb[0], sumsrix[0], sumsrmi[0], sumsrrg[0], sumsrpk[0],
                sumsltb[0], sumslix[0], sumslmi[0], sumslrg[0], sumslpk[0],
               &sumsrtb[0],&sumsrix[0],&sumsrmi[0],&sumsrrg[0],&sumsrpk[0],
               &sumsltb[0],&sumslix[0],&sumslmi[0],&sumslrg[0],&sumslpk[0]);
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
         shvrertb[j] = vrertb[vBSind+j];     // read into shared memory
         shvrerix[j] = vrerix[vBSind+j];
         shvrermi[j] = vrermi[vBSind+j];
         shvrerrg[j] = vrerrg[vBSind+j];
         shvrerpk[j] = vrerpk[vBSind+j];
         shvreltb[j] = vreltb[vBSind+j];
         shvrelix[j] = vrelix[vBSind+j];
         shvrelmi[j] = vrelmi[vBSind+j];
         shvrelrg[j] = vrelrg[vBSind+j];
         shvrelpk[j] = vrelpk[vBSind+j];
         shvimrtb[j] = vimrtb[vBSind+j];
         shvimrix[j] = vimrix[vBSind+j];
         shvimrmi[j] = vimrmi[vBSind+j];
         shvimrrg[j] = vimrrg[vBSind+j];
         shvimrpk[j] = vimlpk[vBSind+j];
         shvimltb[j] = vimltb[vBSind+j];
         shvimlix[j] = vimlix[vBSind+j];
         shvimlmi[j] = vimlmi[vBSind+j];
         shvimlrg[j] = vimlrg[vBSind+j];
         shvimlpk[j] = vimlpk[vBSind+j];
         // normalize vector
         dag_div(shvrertb[j],shvrerix[j],shvrermi[j],shvrerrg[j],shvrerpk[j],
                 shvreltb[j],shvrelix[j],shvrelmi[j],shvrelrg[j],shvrelpk[j],
                  sumsrtb[0], sumsrix[0], sumsrmi[0], sumsrrg[0], sumsrpk[0],
                  sumsltb[0], sumslix[0], sumslmi[0], sumslrg[0], sumslpk[0],
                  &vrertb[vBSind+j],&vrerix[vBSind+j],&vrermi[vBSind+j],
                  &vrerrg[vBSind+j],&vrerpk[vBSind+j],
                  &vreltb[vBSind+j],&vrelix[vBSind+j],&vrelmi[vBSind+j],
                  &vrelrg[vBSind+j],&vrelpk[vBSind+j]);
         dag_div(shvimrtb[j],shvimrix[j],shvimrmi[j],shvimrrg[j],shvimrpk[j],
                 shvimltb[j],shvimlix[j],shvimlmi[j],shvimlrg[j],shvimlpk[j],
                  sumsrtb[0], sumsrix[0], sumsrmi[0], sumsrrg[0], sumsrpk[0],
                  sumsltb[0], sumslix[0], sumslmi[0], sumslrg[0], sumslpk[0],
                  &vimrtb[vBSind+j],&vimrix[vBSind+j],&vimrmi[vBSind+j],
                  &vimrrg[vBSind+j],&vimrpk[vBSind+j],
                  &vimltb[vBSind+j],&vimlix[vBSind+j],&vimlmi[vBSind+j],
                  &vimlrg[vBSind+j],&vimlpk[vBSind+j]);
      }
      __syncthreads();
      vBSind = vBSind + BS;
   }
}

__global__ void large_sum_the_squares
 ( double *vrertb, double *vrerix, double *vrermi, double *vrerrg,
   double *vrerpk, double *vreltb, double *vrelix, double *vrelmi,
   double *vrelrg, double *vrelpk,
   double *vimrtb, double *vimrix, double *vimrmi, double *vimrrg,
   double *vimrpk, double *vimltb, double *vimlix, double *vimlmi,
   double *vimlrg, double *vimlpk, int dim,
   double *sumsrtb, double *sumsrix, double *sumsrmi, double *sumsrrg,
   double *sumsrpk, double *sumsltb, double *sumslix, double *sumslmi,
   double *sumslrg, double *sumslpk, int BS, int BSLog2 )
{
   const int i = blockIdx.x;
   const int j = threadIdx.x;
   const int k = i*BS + j;

   __shared__ double shvrertb[da_shmemsize];
   __shared__ double shvrerix[da_shmemsize];
   __shared__ double shvrermi[da_shmemsize];
   __shared__ double shvrerrg[da_shmemsize];
   __shared__ double shvrerpk[da_shmemsize];
   __shared__ double shvreltb[da_shmemsize];
   __shared__ double shvrelix[da_shmemsize];
   __shared__ double shvrelmi[da_shmemsize];
   __shared__ double shvrelrg[da_shmemsize];
   __shared__ double shvrelpk[da_shmemsize];
   __shared__ double shvimrtb[da_shmemsize];
   __shared__ double shvimrix[da_shmemsize];
   __shared__ double shvimrmi[da_shmemsize];
   __shared__ double shvimrrg[da_shmemsize];
   __shared__ double shvimrpk[da_shmemsize];
   __shared__ double shvimltb[da_shmemsize];
   __shared__ double shvimlix[da_shmemsize];
   __shared__ double shvimlmi[da_shmemsize];
   __shared__ double shvimlrg[da_shmemsize];
   __shared__ double shvimlpk[da_shmemsize];
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
   __shared__ double accrtb[da_shmemsize];
   __shared__ double accrix[da_shmemsize];
   __shared__ double accrmi[da_shmemsize];
   __shared__ double accrrg[da_shmemsize];
   __shared__ double accrpk[da_shmemsize];
   __shared__ double accltb[da_shmemsize];
   __shared__ double acclix[da_shmemsize];
   __shared__ double acclmi[da_shmemsize];
   __shared__ double acclrg[da_shmemsize];
   __shared__ double acclpk[da_shmemsize];

   shvrertb[j] = vrertb[k];
   shvrerix[j] = vrerix[k];
   shvrermi[j] = vrermi[k];
   shvrerrg[j] = vrerrg[k];
   shvrerpk[j] = vrerpk[k];
   shvreltb[j] = vreltb[k];
   shvrelix[j] = vrelix[k];
   shvrelmi[j] = vrelmi[k];
   shvrelrg[j] = vrelrg[k];
   shvrelpk[j] = vrelpk[k];
   shvimrtb[j] = vimrtb[k];
   shvimrix[j] = vimrix[k];
   shvimrmi[j] = vimrmi[k];
   shvimrrg[j] = vimrrg[k];
   shvimrpk[j] = vimrpk[k];
   shvimltb[j] = vimltb[k];
   shvimlix[j] = vimlix[k];
   shvimlmi[j] = vimlmi[k];
   shvimlrg[j] = vimlrg[k];
   shvimlpk[j] = vimlpk[k];

   dag_sqr(shvrertb[j],shvrerix[j],shvrermi[j],shvrerrg[j],shvrerpk[j],
           shvreltb[j],shvrelix[j],shvrelmi[j],shvrelrg[j],shvrelpk[j],
            &accrtb[j], &accrix[j], &accrmi[j], &accrrg[j], &accrpk[j],
            &accltb[j], &acclix[j], &acclmi[j], &acclrg[j], &acclpk[j]);
   dag_sqr(shvimrtb[j],shvimrix[j],shvimrmi[j],shvimrrg[j],shvimrpk[j],
           shvimltb[j],shvimlix[j],shvimlmi[j],shvimlrg[j],shvimlpk[j],
            &prdrtb[j], &prdrix[j], &prdrmi[j], &prdrrg[j], &prdrpk[j],
            &prdltb[j], &prdlix[j], &prdlmi[j], &prdlrg[j], &prdlpk[j]);
   dag_inc(&accrtb[j],&accrix[j],&accrmi[j],&accrrg[j],&accrpk[j],
           &accltb[j],&acclix[j],&acclmi[j],&acclrg[j],&acclpk[j],
            prdrtb[j], prdrix[j], prdrmi[j], prdrrg[j], prdrpk[j],
            prdltb[j], prdlix[j], prdlmi[j], prdlrg[j], prdlpk[j]);

   __syncthreads();

   int powTwo = 1;                          // sum reduction
   for(int L=0; L < BSLog2; L++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < BS)
            dag_inc(&accrtb[j],&accrix[j],&accrmi[j],&accrrg[j],&accrpk[j],
                    &accltb[j],&acclix[j],&acclmi[j],&acclrg[j],&acclpk[j],
                     accrtb[j+powTwo],accrix[j+powTwo],accrmi[j+powTwo],
                     accrrg[j+powTwo],accrpk[j+powTwo],
                     accltb[j+powTwo],acclix[j+powTwo],acclmi[j+powTwo],
                     acclrg[j+powTwo],acclpk[j+powTwo]);
      powTwo = powTwo*2;

      __syncthreads();
   }
   if(j == 0)                               // thread 0 writes the sum
   {
      sumsrtb[i] = accrtb[0];
      sumsrix[i] = accrix[0];
      sumsrmi[i] = accrmi[0];
      sumsrrg[i] = accrrg[0];
      sumsrpk[i] = accrpk[0];
      sumsltb[i] = accltb[0];
      sumslix[i] = acclix[0];
      sumslmi[i] = acclmi[0];
      sumslrg[i] = acclrg[0];
      sumslpk[i] = acclpk[0];
   }
}

__global__ void large_normalize_vector
 ( double *vrertb, double *vrerix, double *vrermi, double *vrerrg,
   double *vrerpk, double *vreltb, double *vrelix, double *vrelmi,
   double *vrelrg, double *vrelpk,
   double *vimrtb, double *vimrix, double *vimrmi, double *vimrrg,
   double *vimrpk, double *vimltb, double *vimlix, double *vimlmi,
   double *vimlrg, double *vimlpk, int dim,
   double *sumsrtb, double *sumsrix, double *sumsrmi, double *sumsrrg,
   double *sumsrpk, double *sumsltb, double *sumslix, double *sumslmi,
   double *sumslrg, double *sumslpk,
   int nbsums, int nbsumsLog2, int BS,
   double *normrtb, double *normrix, double *normrmi, double *normrrg,
   double *normrpk, double *normltb, double *normlix, double *normlmi,
   double *normlrg, double *normlpk )
{
   const int i = blockIdx.x;
   const int j = threadIdx.x;
   const int k = i*BS + j;

   __shared__ double shvrertb[da_shmemsize];
   __shared__ double shvrerix[da_shmemsize];
   __shared__ double shvrermi[da_shmemsize];
   __shared__ double shvrerrg[da_shmemsize];
   __shared__ double shvrerpk[da_shmemsize];
   __shared__ double shvreltb[da_shmemsize];
   __shared__ double shvrelix[da_shmemsize];
   __shared__ double shvrelmi[da_shmemsize];
   __shared__ double shvrelrg[da_shmemsize];
   __shared__ double shvrelpk[da_shmemsize];
   __shared__ double shvimrtb[da_shmemsize];
   __shared__ double shvimrix[da_shmemsize];
   __shared__ double shvimrmi[da_shmemsize];
   __shared__ double shvimrrg[da_shmemsize];
   __shared__ double shvimrpk[da_shmemsize];
   __shared__ double shvimltb[da_shmemsize];
   __shared__ double shvimlix[da_shmemsize];
   __shared__ double shvimlmi[da_shmemsize];
   __shared__ double shvimlrg[da_shmemsize];
   __shared__ double shvimlpk[da_shmemsize];

   if(j < nbsums)
   {
      shvrertb[j] = sumsrtb[j];
      shvrerix[j] = sumsrix[j];
      shvrermi[j] = sumsrmi[j];
      shvrerrg[j] = sumsrrg[j];
      shvrerpk[j] = sumsrpk[j];
      shvreltb[j] = sumsltb[j];
      shvrelix[j] = sumslix[j];
      shvrelmi[j] = sumslmi[j];
      shvrelrg[j] = sumslrg[j];
      shvrelpk[j] = sumslpk[j];
   }

   __syncthreads();

   int powTwo = 1;                          // sum reduction
   for(int L=0; L < nbsumsLog2; L++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < nbsums)
            dag_inc
               (&shvrertb[j],&shvrerix[j],&shvrermi[j],
                &shvrerrg[j],&shvrerpk[j],
                &shvreltb[j],&shvrelix[j],&shvrelmi[j],
                &shvrelrg[j],&shvrelpk[j],
                shvrertb[j+powTwo],shvrerix[j+powTwo],shvrermi[j+powTwo],
                shvrerrg[j+powTwo],shvrerpk[j+powTwo],
                shvreltb[j+powTwo],shvrelix[j+powTwo],shvrelmi[j+powTwo],
                shvrelrg[j+powTwo],shvrelpk[j+powTwo]);
      powTwo = powTwo*2;

      __syncthreads();
   }
   __syncthreads();                      // every thread 0 of all blocks
   if(j == 0)                            // compute the 2-norm and assigns
   {                                     // to the output parameter
      dag_sqrt(shvrertb[0],shvrerix[0],shvrermi[0],shvrerrg[0],shvrerpk[0],
               shvreltb[0],shvrelix[0],shvrelmi[0],shvrelrg[0],shvrelpk[0],
                normrtb,    normrix,    normrmi,    normrrg,    normrpk,
                normltb,    normlix,    normlmi,    normlrg,    normlpk);
   }
   __syncthreads();                    

   if(k < dim)
   {
      shvrertb[j] = vrertb[k];
      shvrerix[j] = vrerix[k];
      shvrermi[j] = vrermi[k];
      shvrerrg[j] = vrerrg[k];
      shvrerpk[j] = vrerpk[k];
      shvreltb[j] = vreltb[k];
      shvrelix[j] = vrelix[k];
      shvrelmi[j] = vrelmi[k];
      shvrelrg[j] = vrelrg[k];
      shvrelpk[j] = vrelpk[k];
      shvimrtb[j] = vimrtb[k];
      shvimrix[j] = vimrix[k];
      shvimrmi[j] = vimrmi[k];
      shvimrrg[j] = vimrrg[k];
      shvimrpk[j] = vimrpk[k];
      shvimltb[j] = vimltb[k];
      shvimlix[j] = vimlix[k];
      shvimlmi[j] = vimlmi[k];
      shvimlrg[j] = vimlrg[k];
      shvimlpk[j] = vimlpk[k];

      dag_div
         ( shvrertb[j], shvrerix[j], shvrermi[j], shvrerrg[j], shvrerpk[j],
           shvreltb[j], shvrelix[j], shvrelmi[j], shvrelrg[j], shvrelpk[j],
           *normrtb,    *normrix,    *normrmi,    *normrrg,    *normrpk,
           *normltb,    *normlix,    *normlmi,    *normlrg,    *normlpk,
          &shvrertb[j],&shvrerix[j],&shvrermi[j],&shvrerrg[j],&shvrerpk[j],
          &shvreltb[j],&shvrelix[j],&shvrelmi[j],&shvrelrg[j],&shvrelpk[j]);
      dag_div
         ( shvimrtb[j], shvimrix[j], shvimrmi[j], shvimrrg[j], shvimrpk[j],
           shvimltb[j], shvimlix[j], shvimlmi[j], shvimlrg[j], shvimlpk[j],
           *normrtb,    *normrix,    *normrmi,    *normrrg,    *normrpk,
           *normltb,    *normlix,    *normlmi,    *normlrg,    *normlpk,
          &shvimrtb[j],&shvimrix[j],&shvimrmi[j],&shvimrrg[j],&shvimrpk[j],
          &shvimltb[j],&shvimlix[j],&shvimlmi[j],&shvimlrg[j],&shvimlpk[j]);

      vrertb[k] = shvrertb[j];
      vrerix[k] = shvrerix[j];
      vrermi[k] = shvrermi[j];
      vrerrg[k] = shvrerrg[j];
      vrerpk[k] = shvrerpk[j];
      vreltb[k] = shvreltb[j];
      vrelix[k] = shvrelix[j];
      vrelmi[k] = shvrelmi[j];
      vrelrg[k] = shvrelrg[j];
      vrelpk[k] = shvrelpk[j];
      vimrtb[k] = shvimrtb[j];
      vimrix[k] = shvimrix[j];
      vimrmi[k] = shvimrmi[j];
      vimrrg[k] = shvimrrg[j];
      vimrpk[k] = shvimrpk[j];
      vimltb[k] = shvimltb[j];
      vimlix[k] = shvimlix[j];
      vimlmi[k] = shvimlmi[j];
      vimlrg[k] = shvimlrg[j];
      vimlpk[k] = shvimlpk[j];
   }
}

void GPU_norm
 ( double *vrertb_h, double *vrerix_h, double *vrermi_h, double *vrerrg_h,
   double *vrerpk_h, double *vreltb_h, double *vrelix_h, double *vrelmi_h,
   double *vrelrg_h, double *vrelpk_h,
   double *vimrtb_h, double *vimrix_h, double *vimrmi_h, double *vimrrg_h,
   double *vimrpk_h, double *vimltb_h, double *vimlix_h, double *vimlmi_h,
   double *vimlrg_h, double *vimlpk_h, int dim, int freq, int BS,
   double *normrtb, double *normrix, double *normrmi, double *normrrg,
   double *normrpk, double *normltb, double *normlix, double *normlmi,
   double *normlrg, double *normlpk, int blocked )
{
   int BSLog2 = ceil(log2((double) BS)); // ceil for sum reduction

   double* vrertb_d;           // highest real parts on device
   double* vrerix_d;           // second highest real parts on device
   double* vrermi_d;           // third highest real parts on device
   double* vrerrg_d;           // fourth highest real parts on device
   double* vrerpk_d;           // fifth highest real parts on device
   double* vreltb_d;           // fifth lowest real parts on device
   double* vrelix_d;           // fourth lowest real parts on device
   double* vrelmi_d;           // third second lowest real parts on device
   double* vrelrg_d;           // second lowest real parts on device
   double* vrelpk_d;           // lowest real parts on device
   double* vimrtb_d;           // highest imaginary parts on device
   double* vimrix_d;           // second highest imaginary parts on device
   double* vimrmi_d;           // third highest imaginary parts on device
   double* vimrrg_d;           // fourth highest imaginary parts on device
   double* vimrpk_d;           // fifth highest imaginary parts on device
   double* vimltb_d;           // fifth lowest imaginary parts on device
   double* vimlix_d;           // fourth lowest imaginary parts on device
   double* vimlmi_d;           // third lowest imaginary parts on device
   double* vimlrg_d;           // second lowest imaginary parts on device
   double* vimlpk_d;           // lowest imaginary parts on device
   size_t size = dim*sizeof(double);
   cudaMalloc((void**)&vrertb_d,size);
   cudaMalloc((void**)&vrerix_d,size);
   cudaMalloc((void**)&vrermi_d,size);
   cudaMalloc((void**)&vrerrg_d,size);
   cudaMalloc((void**)&vrerpk_d,size);
   cudaMalloc((void**)&vreltb_d,size);
   cudaMalloc((void**)&vrelix_d,size);
   cudaMalloc((void**)&vrelmi_d,size);
   cudaMalloc((void**)&vrelrg_d,size);
   cudaMalloc((void**)&vrelpk_d,size);
   cudaMalloc((void**)&vimrtb_d,size);
   cudaMalloc((void**)&vimrix_d,size);
   cudaMalloc((void**)&vimrmi_d,size);
   cudaMalloc((void**)&vimrrg_d,size);
   cudaMalloc((void**)&vimrpk_d,size);
   cudaMalloc((void**)&vimltb_d,size);
   cudaMalloc((void**)&vimlix_d,size);
   cudaMalloc((void**)&vimlmi_d,size);
   cudaMalloc((void**)&vimlrg_d,size);
   cudaMalloc((void**)&vimlpk_d,size);
   cudaMemcpy(vrertb_d,vrertb_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vrerix_d,vrerix_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vrermi_d,vrermi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vrerrg_d,vrerrg_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vrerpk_d,vrerpk_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vreltb_d,vreltb_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vrelix_d,vrelix_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vrelmi_d,vrelmi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vrelrg_d,vrelrg_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vrelpk_d,vrelpk_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vimrtb_d,vimrtb_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vimrix_d,vimrix_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vimrmi_d,vimrmi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vimrrg_d,vimrrg_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vimrpk_d,vimrpk_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vimltb_d,vimltb_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vimlix_d,vimlix_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vimlmi_d,vimlmi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vimlrg_d,vimlrg_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vimlpk_d,vimlpk_h,size,cudaMemcpyHostToDevice);
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
            (vrertb_d,vrerix_d,vrermi_d,vrerrg_d,vrerpk_d,
             vreltb_d,vrelix_d,vrelmi_d,vrelrg_d,vrelpk_d,
             vimrtb_d,vimrix_d,vimrmi_d,vimrrg_d,vimrpk_d,
             vimltb_d,vimlix_d,vimlmi_d,vimlrg_d,vimlpk_d,
             dim,BSLog2,
             normrtb_d,normrix_d,normrmi_d,normrrg_d,normrpk_d,
             normltb_d,normlix_d,normlmi_d,normlrg_d,normlpk_d);
   }
   else if(blocked == 0)
   {
      int rf = ceil(((double) dim)/BS);
      int rfLog2 = ceil(log2((double) rf));
      for(int i=0; i<freq; i++)
         medium_normalize_vector<<<1,BS>>>
            (vrertb_d,vrerix_d,vrermi_d,vrerrg_d,vrerpk_d,
             vreltb_d,vrelix_d,vrelmi_d,vrelrg_d,vrelpk_d,
             vimrtb_d,vimrix_d,vimrmi_d,vimrrg_d,vimrpk_d,
             vimltb_d,vimlix_d,vimlmi_d,vimlrg_d,vimlpk_d,dim,
             rf,rfLog2,BS,BSLog2,
             normrtb_d,normrix_d,normrmi_d,normrrg_d,normrpk_d,
             normltb_d,normlix_d,normlmi_d,normlrg_d,normlpk_d);
   }
   else
   {
      const int nblocks = dim/BS;
      const int nblocksLog2 = ceil(log2((double) nblocks));
      double* sumsrtb_d; // highest parts of sums of squares for each block
      double* sumsrix_d; // 2nd highest parts of sums of squares
      double* sumsrmi_d; // middle parts of sums of squares
      double* sumsrrg_d; // 2nd lowest parts of sums of squares for each block
      double* sumsrpk_d; // lowest parts of sums of squares for each block
      double* sumsltb_d; // highest parts of sums of squares for each block
      double* sumslix_d; // 2nd highest parts of sums of squares
      double* sumslmi_d; // middle parts of sums of squares
      double* sumslrg_d; // 2nd lowest parts of sums of squares for each block
      double* sumslpk_d; // lowest parts of sums of squares for each block
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
            (vrertb_d,vrerix_d,vrermi_d,vrerrg_d,vrerpk_d,
             vreltb_d,vrelix_d,vrelmi_d,vrelrg_d,vrelpk_d,
             vimrtb_d,vimrix_d,vimrmi_d,vimrrg_d,vimrpk_d,
             vimltb_d,vimlix_d,vimlmi_d,vimlrg_d,vimlpk_d,dim,
             sumsrtb_d,sumsrix_d,sumsrmi_d,sumsrrg_d,sumsrpk_d,
             sumsltb_d,sumslix_d,sumslmi_d,sumslrg_d,sumslpk_d,BS,BSLog2);
         large_normalize_vector<<<nblocks,BS>>>
            (vrertb_d,vrerix_d,vrermi_d,vrerrg_d,vrerpk_d,
             vreltb_d,vrelix_d,vrelmi_d,vrelrg_d,vrelpk_d,
             vimrtb_d,vimrix_d,vimrmi_d,vimrrg_d,vimrpk_d,
             vimltb_d,vimlix_d,vimlmi_d,vimlrg_d,vimlpk_d,dim,
             sumsrtb_d,sumsrix_d,sumsrmi_d,sumsrrg_d,sumsrpk_d,
             sumsltb_d,sumslix_d,sumslmi_d,sumslrg_d,sumslpk_d,
             nblocks,nblocksLog2,BS,
             normrtb_d,normrix_d,normrmi_d,normrrg_d,normrpk_d,
             normltb_d,normlix_d,normlmi_d,normlrg_d,normlpk_d);
      }
   }
   cudaMemcpy(vrertb_h,vrertb_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vrerix_h,vrerix_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vrermi_h,vrermi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vrerrg_h,vrerrg_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vrerpk_h,vrerpk_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vreltb_h,vreltb_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vrelix_h,vrelix_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vrelmi_h,vrelmi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vrelrg_h,vrelrg_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vrelpk_h,vrelpk_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vimrtb_h,vimrtb_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vimrix_h,vimrix_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vimrmi_h,vimrmi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vimrrg_h,vimrrg_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vimrpk_h,vimrpk_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vimltb_h,vimltb_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vimlix_h,vimlix_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vimlmi_h,vimlmi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vimlrg_h,vimlrg_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vimlpk_h,vimlpk_d,size,cudaMemcpyDeviceToHost);
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
