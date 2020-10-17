// Defines code of the functions in cmplx8_norm_kernels.h,
// to compute the 2-norm and normalize a complex vector,
// in octo double precision,
// for vectors of small, medium, and large size.

#include "double_double_gpufun.cu"
#include "quad_double_gpufun.cu"
#include "octo_double_gpufun.cu"
#include "cmplx8_norm_kernels.h"

__global__ void small_normalize_vector
 ( double *vrehihihi, double *vrelohihi, double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo, double *vrehilolo, double *vrelololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo, double *vimhilolo, double *vimlololo,
   int dim, int dimLog2, double *normhihihi, double *normlohihi,
   double *normhilohi, double *normlolohi, double *normhihilo, 
   double *normlohilo, double *normhilolo, double *normlololo )
{
   int j = threadIdx.x;

   __shared__ double shvrehihihi[od_shmemsize];
   __shared__ double shvrelohihi[od_shmemsize];
   __shared__ double shvrehilohi[od_shmemsize];
   __shared__ double shvrelolohi[od_shmemsize];
   __shared__ double shvrehihilo[od_shmemsize];
   __shared__ double shvrelohilo[od_shmemsize];
   __shared__ double shvrehilolo[od_shmemsize];
   __shared__ double shvrelololo[od_shmemsize];
   __shared__ double shvimhihihi[od_shmemsize];
   __shared__ double shvimlohihi[od_shmemsize];
   __shared__ double shvimhilohi[od_shmemsize];
   __shared__ double shvimlolohi[od_shmemsize];
   __shared__ double shvimhihilo[od_shmemsize];
   __shared__ double shvimlohilo[od_shmemsize];
   __shared__ double shvimhilolo[od_shmemsize];
   __shared__ double shvimlololo[od_shmemsize];
   __shared__ double prdhihihi[od_shmemsize];
   __shared__ double prdlohihi[od_shmemsize];
   __shared__ double prdhilohi[od_shmemsize];
   __shared__ double prdlolohi[od_shmemsize];
   __shared__ double prdhihilo[od_shmemsize];
   __shared__ double prdlohilo[od_shmemsize];
   __shared__ double prdhilolo[od_shmemsize];
   __shared__ double prdlololo[od_shmemsize];
   __shared__ double sumhihihi[od_shmemsize];
   __shared__ double sumlohihi[od_shmemsize];
   __shared__ double sumhilohi[od_shmemsize];
   __shared__ double sumlolohi[od_shmemsize];
   __shared__ double sumhihilo[od_shmemsize];
   __shared__ double sumlohilo[od_shmemsize];
   __shared__ double sumhilolo[od_shmemsize];
   __shared__ double sumlololo[od_shmemsize];

   shvrehihihi[j] = vrehihihi[j]; // reading real parts into shared memory
   shvrelohihi[j] = vrelohihi[j];
   shvrehilohi[j] = vrehilohi[j];
   shvrelolohi[j] = vrelolohi[j];
   shvrehihilo[j] = vrehihilo[j]; 
   shvrelohilo[j] = vrelohilo[j];
   shvrehilolo[j] = vrehilolo[j];
   shvrelololo[j] = vrelololo[j];
   shvimhihihi[j] = vimhihihi[j]; // reading imag parts into shared memory
   shvimlohihi[j] = vimlohihi[j];
   shvimhilohi[j] = vimhilohi[j];
   shvimlolohi[j] = vimlolohi[j];
   shvimhihilo[j] = vimhihilo[j]; 
   shvimlohilo[j] = vimlohilo[j];
   shvimhilolo[j] = vimhilolo[j];
   shvimlololo[j] = vimlololo[j];

   odg_sqr(shvrehihihi[j],shvrelohihi[j],shvrehilohi[j],shvrelolohi[j],
           shvrehihilo[j],shvrelohilo[j],shvrehilolo[j],shvrelololo[j],
            &sumhihihi[j], &sumlohihi[j], &sumhilohi[j], &sumlolohi[j],
            &sumhihilo[j], &sumlohilo[j], &sumhilolo[j], &sumlololo[j]);
   odg_sqr(shvimhihihi[j],shvimlohihi[j],shvimhilohi[j],shvimlolohi[j],
           shvimhihilo[j],shvimlohilo[j],shvimhilolo[j],shvimlololo[j],
            &prdhihihi[j], &prdlohihi[j], &prdhilohi[j], &prdlolohi[j],
            &prdhihilo[j], &prdlohilo[j], &prdhilolo[j], &prdlololo[j]);
   odg_inc(&sumhihihi[j],&sumlohihi[j],&sumhilohi[j],&sumlolohi[j],
           &sumhihilo[j],&sumlohilo[j],&sumhilolo[j],&sumlololo[j],
            prdhihihi[j], prdlohihi[j], prdhilohi[j], prdlolohi[j],
            prdhihilo[j], prdlohilo[j], prdhilolo[j], prdlololo[j]);

   __syncthreads();
   int powTwo = 1;                          // sum reduction
   for(int k=0; k < dimLog2; k++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < dim)
            odg_inc(&sumhihihi[j],&sumlohihi[j],&sumhilohi[j],&sumlolohi[j],
                    &sumhihilo[j],&sumlohilo[j],&sumhilolo[j],&sumlololo[j],
                     sumhihihi[j+powTwo],sumlohihi[j+powTwo],
                     sumhilohi[j+powTwo],sumlolohi[j+powTwo],
                     sumhihilo[j+powTwo],sumlohilo[j+powTwo],
                     sumhilolo[j+powTwo],sumlololo[j+powTwo]);
      powTwo = powTwo*2;
      __syncthreads();
   }
   // thread 0 computes the sqrt of the inner product, others wait
   if(j == 0)
      odg_sqrt( sumhihihi[0], sumlohihi[0], sumhilohi[0], sumlolohi[0],
                sumhihilo[0], sumlohilo[0], sumhilolo[0], sumlololo[0],
               &sumhihihi[0],&sumlohihi[0],&sumhilohi[0],&sumlolohi[0],
               &sumhihilo[0],&sumlohilo[0],&sumhilolo[0],&sumlololo[0]); 
   if(j == 0)
   {
      *normhihihi = sumhihihi[0];
      *normlohihi = sumlohihi[0];
      *normhilohi = sumhilohi[0];
      *normlolohi = sumlolohi[0];
      *normhihilo = sumhihilo[0];
      *normlohilo = sumlohilo[0];
      *normhilolo = sumhilolo[0];
      *normlololo = sumlololo[0];
   }
   __syncthreads();
   odg_div(shvrehihihi[j],shvrelohihi[j],shvrehilohi[j],shvrelolohi[j],
           shvrehihilo[j],shvrelohilo[j],shvrehilolo[j],shvrelololo[j],
             sumhihihi[0],  sumlohihi[0],  sumhilohi[0],  sumlolohi[0],
             sumhihilo[0],  sumlohilo[0],  sumhilolo[0],  sumlololo[0],
            &vrehihihi[j], &vrelohihi[j], &vrehilohi[j], &vrelolohi[j],
            &vrehihilo[j], &vrelohilo[j], &vrehilolo[j], &vrelololo[j]);
   odg_div(shvimhihihi[j],shvimlohihi[j],shvimhilohi[j],shvimlolohi[j],
           shvimhihilo[j],shvimlohilo[j],shvimhilolo[j],shvimlololo[j],
             sumhihihi[0],  sumlohihi[0],  sumhilohi[0],  sumlolohi[0],
             sumhihilo[0],  sumlohilo[0],  sumhilolo[0],  sumlololo[0],
            &vimhihihi[j], &vimlohihi[j], &vimhilohi[j], &vimlolohi[j],
            &vimhihilo[j], &vimlohilo[j], &vimhilolo[j], &vimlololo[j]);
}

__global__ void medium_normalize_vector
 ( double *vrehihihi, double *vrelohihi, double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo, double *vrehilolo, double *vrelololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo, double *vimhilolo, double *vimlololo,
   int dim, int rnd, int rndLog2, int BS, int BSLog2,
   double *normhihihi, double *normlohihi, double *normhilohi,
   double *normlolohi, double *normhihilo, double *normlohilo,
   double *normhilolo, double *normlololo )
{
   int j = threadIdx.x;
   int powTwo;
   int vBSind = 0;

   __shared__ double shvrehihihi[od_shmemsize];
   __shared__ double shvrelohihi[od_shmemsize];
   __shared__ double shvrehilohi[od_shmemsize];
   __shared__ double shvrelolohi[od_shmemsize];
   __shared__ double shvrehihilo[od_shmemsize];
   __shared__ double shvrelohilo[od_shmemsize];
   __shared__ double shvrehilolo[od_shmemsize];
   __shared__ double shvrelololo[od_shmemsize];
   __shared__ double shvimhihihi[od_shmemsize];
   __shared__ double shvimlohihi[od_shmemsize];
   __shared__ double shvimhilohi[od_shmemsize];
   __shared__ double shvimlolohi[od_shmemsize];
   __shared__ double shvimhihilo[od_shmemsize];
   __shared__ double shvimlohilo[od_shmemsize];
   __shared__ double shvimhilolo[od_shmemsize];
   __shared__ double shvimlololo[od_shmemsize];
   __shared__ double prdhihihi[od_shmemsize];
   __shared__ double prdlohihi[od_shmemsize];
   __shared__ double prdhilohi[od_shmemsize];
   __shared__ double prdlolohi[od_shmemsize];
   __shared__ double prdhihilo[od_shmemsize];
   __shared__ double prdlohilo[od_shmemsize];
   __shared__ double prdhilolo[od_shmemsize];
   __shared__ double prdlololo[od_shmemsize];
   __shared__ double acchihihi[od_shmemsize];
   __shared__ double acclohihi[od_shmemsize];
   __shared__ double acchilohi[od_shmemsize];
   __shared__ double acclolohi[od_shmemsize];
   __shared__ double acchihilo[od_shmemsize];
   __shared__ double acclohilo[od_shmemsize];
   __shared__ double acchilolo[od_shmemsize];
   __shared__ double acclololo[od_shmemsize];
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
         acchihihi[j] = 0.0;
         acclohihi[j] = 0.0;
         acchilohi[j] = 0.0;
         acclolohi[j] = 0.0;
         acchihilo[j] = 0.0;
         acclohilo[j] = 0.0;
         acchilolo[j] = 0.0;
         acclololo[j] = 0.0;
      }
      else
      {
         shvrehihihi[j] = vrehihihi[vBSind+j];  // reading into shared memory
         shvrelohihi[j] = vrelohihi[vBSind+j];
         shvrehilohi[j] = vrehilohi[vBSind+j];
         shvrelolohi[j] = vrelolohi[vBSind+j];
         shvrehihilo[j] = vrehihilo[vBSind+j]; 
         shvrelohilo[j] = vrelohilo[vBSind+j];
         shvrehilolo[j] = vrehilolo[vBSind+j];
         shvrelololo[j] = vrelololo[vBSind+j];
         shvimhihihi[j] = vimhihihi[vBSind+j]; 
         shvimlohihi[j] = vimlohihi[vBSind+j]; 
         shvimhilohi[j] = vimhilohi[vBSind+j]; 
         shvimlolohi[j] = vimlolohi[vBSind+j]; 
         shvimhihilo[j] = vimhihilo[vBSind+j]; 
         shvimlohilo[j] = vimlohilo[vBSind+j]; 
         shvimhilolo[j] = vimhilolo[vBSind+j]; 
         shvimlololo[j] = vimlololo[vBSind+j]; 

         odg_sqr(shvrehihihi[j],shvrelohihi[j],shvrehilohi[j],shvrelolohi[j],
                 shvrehihilo[j],shvrelohilo[j],shvrehilolo[j],shvrelololo[j],
                  &acchihihi[j], &acclohihi[j], &acchilohi[j], &acclolohi[j],
                  &acchihilo[j], &acclohilo[j], &acchilolo[j], &acclololo[j]);
         odg_sqr(shvimhihihi[j],shvimlohihi[j],shvimhilohi[j],shvimlolohi[j],
                 shvimhihilo[j],shvimlohilo[j],shvimhilolo[j],shvimlololo[j],
                  &prdhihihi[j], &prdlohihi[j], &prdhilohi[j], &prdlolohi[j],
                  &prdhihilo[j], &prdlohilo[j], &prdhilolo[j], &prdlololo[j]);
         odg_inc(&acchihihi[j],&acclohihi[j],&acchilohi[j],&acclolohi[j],
                 &acchihilo[j],&acclohilo[j],&acchilolo[j],&acclololo[j],
                  prdhihihi[j], prdlohihi[j], prdhilohi[j], prdlolohi[j],
                  prdhihilo[j], prdlohilo[j], prdhilolo[j], prdlololo[j]);
      }
      __syncthreads();
      powTwo = 1;                          // sum reduction
      for(int k=0; k < BSLog2; k++)
      {
         if((j%(powTwo*2)) == 0)
            if(j+powTwo < BS)
               odg_inc
                  (&acchihihi[j],&acclohihi[j],&acchilohi[j],&acclolohi[j],
                   &acchihilo[j],&acclohilo[j],&acchilolo[j],&acclololo[j],
                    acchihihi[j+powTwo],acclohihi[j+powTwo],
                    acchilohi[j+powTwo],acclolohi[j+powTwo],
                    acchihilo[j+powTwo],acclohilo[j+powTwo],
                    acchilolo[j+powTwo],acclololo[j+powTwo]);
         powTwo = powTwo*2;
         __syncthreads();
      }
      // thread 0 copies the sum of this round in sums[i], the others wait
      if(j == 0)
      {
         sumshihihi[i] = acchihihi[0]; 
         sumslohihi[i] = acclohihi[0]; 
         sumshilohi[i] = acchilohi[0]; 
         sumslolohi[i] = acclolohi[0]; 
         sumshihilo[i] = acchihilo[0]; 
         sumslohilo[i] = acclohilo[0]; 
         sumshilolo[i] = acchilolo[0]; 
         sumslololo[i] = acclololo[0]; 
      }
      __syncthreads();
      vBSind = vBSind + BS;
   }
   powTwo = 1;                          // sum reduction
   for(int k=0; k < rndLog2; k++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < rnd)
            odg_inc
               (&sumshihihi[j],&sumslohihi[j],&sumshilohi[j],&sumslolohi[j],
                &sumshihilo[j],&sumslohilo[j],&sumshilolo[j],&sumslololo[j],
                 sumshihihi[j+powTwo],sumslohihi[j+powTwo],
                 sumshilohi[j+powTwo],sumslolohi[j+powTwo],
                 sumshihilo[j+powTwo],sumslohilo[j+powTwo],
                 sumshilolo[j+powTwo],sumslololo[j+powTwo]);
      powTwo = powTwo*2;
      __syncthreads();
   }
   if(j == 0)
   {
      odg_sqrt( sumshihihi[0], sumslohihi[0], sumshilohi[0], sumslolohi[0],
                sumshihilo[0], sumslohilo[0], sumshilolo[0], sumslololo[0],
               &sumshihihi[0],&sumslohihi[0],&sumshilohi[0],&sumslolohi[0],
               &sumshihilo[0],&sumslohilo[0],&sumshilolo[0],&sumslololo[0]);
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
         shvrehihihi[j] = vrehihihi[vBSind+j];  // read into shared memory
         shvrelohihi[j] = vrelohihi[vBSind+j];
         shvrehilohi[j] = vrehilohi[vBSind+j];
         shvrelolohi[j] = vrelolohi[vBSind+j];
         shvrehihilo[j] = vrehihilo[vBSind+j];
         shvrelohilo[j] = vrelohilo[vBSind+j];
         shvrehilolo[j] = vrehilolo[vBSind+j];
         shvrelololo[j] = vrelololo[vBSind+j];
         shvimhihihi[j] = vimhihihi[vBSind+j];
         shvimlohihi[j] = vimlohihi[vBSind+j];
         shvimhilohi[j] = vimhilohi[vBSind+j];
         shvimlolohi[j] = vimlolohi[vBSind+j];
         shvimhihilo[j] = vimhihilo[vBSind+j];
         shvimlohilo[j] = vimlohilo[vBSind+j];
         shvimhilolo[j] = vimhilolo[vBSind+j];
         shvimlololo[j] = vimlololo[vBSind+j];
         // normalize vector
         odg_div(shvrehihihi[j],shvrelohihi[j],shvrehilohi[j],shvrelolohi[j],
                 shvrehihilo[j],shvrelohilo[j],shvrehilolo[j],shvrelololo[j],
                  sumshihihi[0], sumslohihi[0], sumshilohi[0], sumslolohi[0],
                  sumshihilo[0], sumslohilo[0], sumshilolo[0], sumslololo[0],
                  &vrehihihi[vBSind+j],&vrelohihi[vBSind+j],
                  &vrehilohi[vBSind+j],&vrelolohi[vBSind+j],
                  &vrehihilo[vBSind+j],&vrelohilo[vBSind+j],
                  &vrehilolo[vBSind+j],&vrelololo[vBSind+j]);
         odg_div(shvimhihihi[j],shvimlohihi[j],shvimhilohi[j],shvimlolohi[j],
                 shvimhihilo[j],shvimlohilo[j],shvimhilolo[j],shvimlololo[j],
                  sumshihihi[0], sumslohihi[0], sumshilohi[0], sumslolohi[0],
                  sumshihilo[0], sumslohilo[0], sumshilolo[0], sumslololo[0],
                  &vimhihihi[vBSind+j],&vimlohihi[vBSind+j],
                  &vimhilohi[vBSind+j],&vimlolohi[vBSind+j],
                  &vimhihilo[vBSind+j],&vimlohilo[vBSind+j],
                  &vimhilolo[vBSind+j],&vimlololo[vBSind+j]);
      }
      __syncthreads();
      vBSind = vBSind + BS;
   }
}

__global__ void large_sum_the_squares
 ( double *vrehihihi, double *vrelohihi, double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo, double *vrehilolo, double *vrelololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo, double *vimhilolo, double *vimlololo,
   int dim,
   double *sumshihihi, double *sumslohihi, double *sumshilohi,
   double *sumslolohi, double *sumshihilo, double *sumslohilo,
   double *sumshilolo, double *sumslololo, int BS, int BSLog2 )
{
   const int i = blockIdx.x;
   const int j = threadIdx.x;
   const int k = i*BS + j;

   __shared__ double shvrehihihi[od_shmemsize];
   __shared__ double shvrelohihi[od_shmemsize];
   __shared__ double shvrehilohi[od_shmemsize];
   __shared__ double shvrelolohi[od_shmemsize];
   __shared__ double shvrehihilo[od_shmemsize];
   __shared__ double shvrelohilo[od_shmemsize];
   __shared__ double shvrehilolo[od_shmemsize];
   __shared__ double shvrelololo[od_shmemsize];
   __shared__ double shvimhihihi[od_shmemsize];
   __shared__ double shvimlohihi[od_shmemsize];
   __shared__ double shvimhilohi[od_shmemsize];
   __shared__ double shvimlolohi[od_shmemsize];
   __shared__ double shvimhihilo[od_shmemsize];
   __shared__ double shvimlohilo[od_shmemsize];
   __shared__ double shvimhilolo[od_shmemsize];
   __shared__ double shvimlololo[od_shmemsize];
   __shared__ double prdhihihi[od_shmemsize];
   __shared__ double prdlohihi[od_shmemsize];
   __shared__ double prdhilohi[od_shmemsize];
   __shared__ double prdlolohi[od_shmemsize];
   __shared__ double prdhihilo[od_shmemsize];
   __shared__ double prdlohilo[od_shmemsize];
   __shared__ double prdhilolo[od_shmemsize];
   __shared__ double prdlololo[od_shmemsize];
   __shared__ double acchihihi[od_shmemsize];
   __shared__ double acclohihi[od_shmemsize];
   __shared__ double acchilohi[od_shmemsize];
   __shared__ double acclolohi[od_shmemsize];
   __shared__ double acchihilo[od_shmemsize];
   __shared__ double acclohilo[od_shmemsize];
   __shared__ double acchilolo[od_shmemsize];
   __shared__ double acclololo[od_shmemsize];

   shvrehihihi[j] = vrehihihi[k];
   shvrelohihi[j] = vrelohihi[k];
   shvrehilohi[j] = vrehilohi[k];
   shvrelolohi[j] = vrelolohi[k];
   shvrehihilo[j] = vrehihilo[k];
   shvrelohilo[j] = vrelohilo[k];
   shvrehilolo[j] = vrehilolo[k];
   shvrelololo[j] = vrelololo[k];
   shvimhihihi[j] = vimhihihi[k];
   shvimlohihi[j] = vimlohihi[k];
   shvimhilohi[j] = vimhilohi[k];
   shvimlolohi[j] = vimlolohi[k];
   shvimhihilo[j] = vimhihilo[k];
   shvimlohilo[j] = vimlohilo[k];
   shvimhilolo[j] = vimhilolo[k];
   shvimlololo[j] = vimlololo[k];

   odg_sqr(shvrehihihi[j],shvrelohihi[j],shvrehilohi[j],shvrelolohi[j],
           shvrehihilo[j],shvrelohilo[j],shvrehilolo[j],shvrelololo[j],
            &acchihihi[j], &acclohihi[j], &acchilohi[j], &acclolohi[j],
            &acchihilo[j], &acclohilo[j], &acchilolo[j], &acclololo[j]);
   odg_sqr(shvimhihihi[j],shvimlohihi[j],shvimhilohi[j],shvimlolohi[j],
           shvimhihilo[j],shvimlohilo[j],shvimhilolo[j],shvimlololo[j],
            &prdhihihi[j], &prdlohihi[j], &prdhilohi[j], &prdlolohi[j],
            &prdhihilo[j], &prdlohilo[j], &prdhilolo[j], &prdlololo[j]);
   odg_inc(&acchihihi[j],&acclohihi[j],&acchilohi[j],&acclolohi[j],
           &acchihilo[j],&acclohilo[j],&acchilolo[j],&acclololo[j],
            prdhihihi[j], prdlohihi[j], prdhilohi[j], prdlolohi[j],
            prdhihilo[j], prdlohilo[j], prdhilolo[j], prdlololo[j]);

   __syncthreads();

   int powTwo = 1;                          // sum reduction
   for(int L=0; L < BSLog2; L++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < BS)
            odg_inc(&acchihihi[j],&acclohihi[j],&acchilohi[j],&acclolohi[j],
                    &acchihilo[j],&acclohilo[j],&acchilolo[j],&acclololo[j],
                     acchihihi[j+powTwo],acclohihi[j+powTwo],
                     acchilohi[j+powTwo],acclolohi[j+powTwo],
                     acchihilo[j+powTwo],acclohilo[j+powTwo],
                     acchilolo[j+powTwo],acclololo[j+powTwo]);
      powTwo = powTwo*2;

      __syncthreads();
   }
   if(j == 0)                               // thread 0 writes the sum
   {
      sumshihihi[i] = acchihihi[0];
      sumslohihi[i] = acclohihi[0];
      sumshilohi[i] = acchilohi[0];
      sumslolohi[i] = acclolohi[0];
      sumshihilo[i] = acchihilo[0];
      sumslohilo[i] = acclohilo[0];
      sumshilolo[i] = acchilolo[0];
      sumslololo[i] = acclololo[0];
   }
}

__global__ void large_normalize_vector
 ( double *vrehihihi, double *vrelohihi, double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo, double *vrehilolo, double *vrelololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo, double *vimhilolo, double *vimlololo,
   int dim,
   double *sumshihihi, double *sumslohihi, double *sumshilohi,
   double *sumslolohi, double *sumshihilo, double *sumslohilo,
   double *sumshilolo, double *sumslololo,
   int nbsums, int nbsumsLog2, int BS, double *normhihihi, double *normlohihi,
   double *normhilohi, double *normlolohi, double *normhihilo,
   double *normlohilo, double *normhilolo, double *normlololo )
{
   const int i = blockIdx.x;
   const int j = threadIdx.x;
   const int k = i*BS + j;

   __shared__ double shvrehihihi[od_shmemsize];
   __shared__ double shvrelohihi[od_shmemsize];
   __shared__ double shvrehilohi[od_shmemsize];
   __shared__ double shvrelolohi[od_shmemsize];
   __shared__ double shvrehihilo[od_shmemsize];
   __shared__ double shvrelohilo[od_shmemsize];
   __shared__ double shvrehilolo[od_shmemsize];
   __shared__ double shvrelololo[od_shmemsize];
   __shared__ double shvimhihihi[od_shmemsize];
   __shared__ double shvimlohihi[od_shmemsize];
   __shared__ double shvimhilohi[od_shmemsize];
   __shared__ double shvimlolohi[od_shmemsize];
   __shared__ double shvimhihilo[od_shmemsize];
   __shared__ double shvimlohilo[od_shmemsize];
   __shared__ double shvimhilolo[od_shmemsize];
   __shared__ double shvimlololo[od_shmemsize];

   if(j < nbsums)
   {
      shvrehihihi[j] = sumshihihi[j];
      shvrelohihi[j] = sumslohihi[j];
      shvrehilohi[j] = sumshilohi[j];
      shvrelolohi[j] = sumslolohi[j];
      shvrehihilo[j] = sumshihilo[j];
      shvrelohilo[j] = sumslohilo[j];
      shvrehilolo[j] = sumshilolo[j];
      shvrelololo[j] = sumslololo[j];
   }

   __syncthreads();

   int powTwo = 1;                          // sum reduction
   for(int L=0; L < nbsumsLog2; L++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < nbsums)
            odg_inc(&shvrehihihi[j],&shvrelohihi[j],&shvrehilohi[j],
                    &shvrelolohi[j],&shvrehihilo[j],&shvrelohilo[j],
                    &shvrehilolo[j],&shvrelololo[j],
                    shvrehihihi[j+powTwo],shvrelohihi[j+powTwo],
                    shvrehilohi[j+powTwo],shvrelolohi[j+powTwo],
                    shvrehihilo[j+powTwo],shvrelohilo[j+powTwo],
                    shvrehilolo[j+powTwo],shvrelololo[j+powTwo]);
      powTwo = powTwo*2;

      __syncthreads();
   }
   __syncthreads();                      // every thread 0 of all blocks
   if(j == 0)                            // compute the 2-norm and assigns
   {                                     // to the output parameter
      odg_sqrt(shvrehihihi[0],shvrelohihi[0],shvrehilohi[0],shvrelolohi[0],
               shvrehihilo[0],shvrelohilo[0],shvrehilolo[0],shvrelololo[0],
                normhihihi,    normlohihi,    normhilohi,    normlolohi,
                normhihilo,    normlohilo,    normhilolo,    normlololo);
   }
   __syncthreads();                    

   if(k < dim)
   {
      shvrehihihi[j] = vrehihihi[k];
      shvrelohihi[j] = vrelohihi[k];
      shvrehilohi[j] = vrehilohi[k];
      shvrelolohi[j] = vrelolohi[k];
      shvrehihilo[j] = vrehihilo[k];
      shvrelohilo[j] = vrelohilo[k];
      shvrehilolo[j] = vrehilolo[k];
      shvrelololo[j] = vrelololo[k];
      shvimhihihi[j] = vimhihihi[k];
      shvimlohihi[j] = vimlohihi[k];
      shvimhilohi[j] = vimhilohi[k];
      shvimlolohi[j] = vimlolohi[k];
      shvimhihilo[j] = vimhihilo[k];
      shvimlohilo[j] = vimlohilo[k];
      shvimhilolo[j] = vimhilolo[k];
      shvimlololo[j] = vimlololo[k];

      odg_div
         ( shvrehihihi[j], shvrelohihi[j], shvrehilohi[j], shvrelolohi[j],
           shvrehihilo[j], shvrelohilo[j], shvrehilolo[j], shvrelololo[j],
           *normhihihi,    *normlohihi,    *normhilohi,    *normlolohi,
           *normhihilo,    *normlohilo,    *normhilolo,    *normlololo,
          &shvrehihihi[j],&shvrelohihi[j],&shvrehilohi[j],&shvrelolohi[j],
          &shvrehihilo[j],&shvrelohilo[j],&shvrehilolo[j],&shvrelololo[j]);
      odg_div
         ( shvimhihihi[j], shvimlohihi[j], shvimhilohi[j], shvimlolohi[j],
           shvimhihilo[j], shvimlohilo[j], shvimhilolo[j], shvimlololo[j],
           *normhihihi,    *normlohihi,    *normhilohi,    *normlolohi,
           *normhihilo,    *normlohilo,    *normhilolo,    *normlololo,
           &shvimhihihi[j],&shvimlohihi[j],&shvimhilohi[j],&shvimlolohi[j],
           &shvimhihilo[j],&shvimlohilo[j],&shvimhilolo[j],&shvimlololo[j]);

      vrehihihi[k] = shvrehihihi[j];
      vrelohihi[k] = shvrelohihi[j];
      vrehilohi[k] = shvrehilohi[j];
      vrelolohi[k] = shvrelolohi[j];
      vrehihilo[k] = shvrehihilo[j];
      vrelohilo[k] = shvrelohilo[j];
      vrehilolo[k] = shvrehilolo[j];
      vrelololo[k] = shvrelololo[j];
      vimhihihi[k] = shvimhihihi[j];
      vimlohihi[k] = shvimlohihi[j];
      vimhilohi[k] = shvimhilohi[j];
      vimlolohi[k] = shvimlolohi[j];
      vimhihilo[k] = shvimhihilo[j];
      vimlohilo[k] = shvimlohilo[j];
      vimhilolo[k] = shvimhilolo[j];
      vimlololo[k] = shvimlololo[j];
  }
}

void GPU_norm
 ( double *vrehihihi_h, double *vrelohihi_h, double *vrehilohi_h,
   double *vrelolohi_h, double *vrehihilo_h, double *vrelohilo_h,
   double *vrehilolo_h, double *vrelololo_h,
   double *vimhihihi_h, double *vimlohihi_h, double *vimhilohi_h,
   double *vimlolohi_h, double *vimhihilo_h, double *vimlohilo_h,
   double *vimhilolo_h, double *vimlololo_h,
   int dim, int freq, int BS, double *normhihihi, double *normlohihi,
   double *normhilohi, double *normlolohi, double *normhihilo,
   double *normlohilo, double *normhilolo, double *normlololo, int blocked )
{
   int BSLog2 = ceil(log2((double) BS)); // ceil for sum reduction

   double* vrehihihi_d;       // highest real parts on device
   double* vrelohihi_d;       // second highest real parts on device
   double* vrehilohi_d;       // third highest real parts on device
   double* vrelolohi_d;       // fourth highest real parts on device
   double* vrehihilo_d;       // fourth lowest real parts on device
   double* vrelohilo_d;       // third lowest real parts on device
   double* vrehilolo_d;       // second lowest real parts on device
   double* vrelololo_d;       // lowest real parts on device
   double* vimhihihi_d;       // highest imaginary parts on device
   double* vimlohihi_d;       // second highest imaginary parts on device
   double* vimhilohi_d;       // third highest imaginary parts on device
   double* vimlolohi_d;       // fourth highest real parts on device
   double* vimhihilo_d;       // fourth lowest imaginary parts on device
   double* vimlohilo_d;       // third lowest imaginary parts on device
   double* vimhilolo_d;       // second lowest imaginary parts on device
   double* vimlololo_d;       // lowest imaginary parts on device
   size_t size = dim*sizeof(double);
   cudaMalloc((void**)&vrehihihi_d,size);
   cudaMalloc((void**)&vrelohihi_d,size);
   cudaMalloc((void**)&vrehilohi_d,size);
   cudaMalloc((void**)&vrelolohi_d,size);
   cudaMalloc((void**)&vrehihilo_d,size);
   cudaMalloc((void**)&vrelohilo_d,size);
   cudaMalloc((void**)&vrehilolo_d,size);
   cudaMalloc((void**)&vrelololo_d,size);
   cudaMalloc((void**)&vimhihihi_d,size);
   cudaMalloc((void**)&vimlohihi_d,size);
   cudaMalloc((void**)&vimhilohi_d,size);
   cudaMalloc((void**)&vimlolohi_d,size);
   cudaMalloc((void**)&vimhihilo_d,size);
   cudaMalloc((void**)&vimlohilo_d,size);
   cudaMalloc((void**)&vimhilolo_d,size);
   cudaMalloc((void**)&vimlololo_d,size);
   cudaMemcpy(vrehihihi_d,vrehihihi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vrelohihi_d,vrelohihi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vrehilohi_d,vrehilohi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vrelolohi_d,vrelolohi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vrehihilo_d,vrehihilo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vrelohilo_d,vrelohilo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vrehilolo_d,vrehilolo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vrelololo_d,vrelololo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vimhihihi_d,vimhihihi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vimlohihi_d,vimlohihi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vimhilohi_d,vimhilohi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vimlolohi_d,vimlolohi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vimhihilo_d,vimhihilo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vimlohilo_d,vimlohilo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vimhilolo_d,vimhilolo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(vimlololo_d,vimlololo_h,size,cudaMemcpyHostToDevice);
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
            (vrehihihi_d,vrelohihi_d,vrehilohi_d,vrelolohi_d,
             vrehihilo_d,vrelohilo_d,vrehilolo_d,vrelololo_d,
             vimhihihi_d,vimlohihi_d,vimhilohi_d,vimlolohi_d,
             vimhihilo_d,vimlohilo_d,vimhilolo_d,vimlololo_d,
             dim,BSLog2,
             normhihihi_d,normlohihi_d,normhilohi_d,normlolohi_d,
             normhihilo_d,normlohilo_d,normhilolo_d,normlololo_d);
   }
   else if(blocked == 0)
   {
      int rf = ceil(((double) dim)/BS);
      int rfLog2 = ceil(log2((double) rf));
      for(int i=0; i<freq; i++)
         medium_normalize_vector<<<1,BS>>>
            (vrehihihi_d,vrelohihi_d,vrehilohi_d,vrelolohi_d,
             vrehihilo_d,vrelohilo_d,vrehilolo_d,vrelololo_d,
             vimhihihi_d,vimlohihi_d,vimhilohi_d,vimlolohi_d,
             vimhihilo_d,vimlohilo_d,vimhilolo_d,vimlololo_d,dim,
             rf,rfLog2,BS,BSLog2,
             normhihihi_d,normlohihi_d,normhilohi_d,normlolohi_d,
             normhihilo_d,normlohilo_d,normhilolo_d,normlololo_d);
   }
   else
   {
      const int nblocks = dim/BS;
      const int nblocksLog2 = ceil(log2((double) nblocks));
      double* sumshihihi_d; // highest parts of sums of squares for each block
      double* sumslohihi_d; // 2nd highest parts of sums of squares
      double* sumshilohi_d; // 3rd highest parts of sums of squares
      double* sumslolohi_d; // 4th highest parts of sums of squares
      double* sumshihilo_d; // 4th lowest parts of sums of squares
      double* sumslohilo_d; // 3rd lowest parts of sums of squares
      double* sumshilolo_d; // 2nd lowest parts of sums of squares 
      double* sumslololo_d; // lowest parts of sums of squares for each block
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
            (vrehihihi_d,vrelohihi_d,vrehilohi_d,vrelolohi_d,
             vrehihilo_d,vrelohilo_d,vrehilolo_d,vrelololo_d,
             vimhihihi_d,vimlohihi_d,vimhilohi_d,vimlolohi_d,
             vimhihilo_d,vimlohilo_d,vimhilolo_d,vimlololo_d,dim,
             sumshihihi_d,sumslohihi_d,sumshilohi_d,sumslolohi_d,
             sumshihilo_d,sumslohilo_d,sumshilolo_d,sumslololo_d,BS,BSLog2);
         large_normalize_vector<<<nblocks,BS>>>
            (vrehihihi_d,vrelohihi_d,vrehilohi_d,vrelolohi_d,
             vrehihilo_d,vrelohilo_d,vrehilolo_d,vrelololo_d,
             vimhihihi_d,vimlohihi_d,vimhilohi_d,vimlolohi_d,
             vimhihilo_d,vimlohilo_d,vimhilolo_d,vimlololo_d,dim,
             sumshihihi_d,sumslohihi_d,sumshilohi_d,sumslolohi_d,
             sumshihilo_d,sumslohilo_d,sumshilolo_d,sumslololo_d,
             nblocks,nblocksLog2,BS,
             normhihihi_d,normlohihi_d,normhilohi_d,normlolohi_d,
             normhihilo_d,normlohilo_d,normhilolo_d,normlololo_d);
      }
   }
   cudaMemcpy(vrehihihi_h,vrehihihi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vrelohihi_h,vrelohihi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vrehilohi_h,vrehilohi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vrelolohi_h,vrelolohi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vrehihilo_h,vrehihilo_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vrelohilo_h,vrelohilo_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vrehilolo_h,vrehilolo_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vrelololo_h,vrelololo_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vimhihihi_h,vimhihihi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vimlohihi_h,vimlohihi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vimhilohi_h,vimhilohi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vimlolohi_h,vimlolohi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vimhihilo_h,vimhihilo_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vimlohilo_h,vimlohilo_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vimhilolo_h,vimhilolo_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(vimlololo_h,vimlololo_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(normhihihi,normhihihi_d,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(normlohihi,normlohihi_d,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(normhilohi,normhilohi_d,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(normlolohi,normlolohi_d,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(normhihilo,normhihilo_d,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(normlohilo,normlohilo_d,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(normhilolo,normhilolo_d,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(normlololo,normlololo_d,sizeof(double),cudaMemcpyDeviceToHost);
}
