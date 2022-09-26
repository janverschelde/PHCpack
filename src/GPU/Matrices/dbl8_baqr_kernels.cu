/* The file dbl8_baqr_kernels.cu defines the functions with prototypes in
 * the file dbl8_baqr_kernels.h. */

#include <iostream>
#include <iomanip>
#ifdef winwalltime
#include "gettimeofday4win.h"
#else
#include <sys/time.h>
#endif
#ifdef gpufun
#include "double_double_gpufun.cu"
#include "quad_double_gpufun.cu"
#include "octo_double_gpufun.cu"
#endif
#include "dbl8_baqr_kernels.h"
#include "octo_double_functions.h"
#include "dbl_baqr_flopcounts.h"

using namespace std;

__global__ void dbl8_small_house
 ( double *x0hihihi, double *x0lohihi, double *x0hilohi, double *x0lolohi,
   double *x0hihilo, double *x0lohilo, double *x0hilolo, double *x0lololo,
   double *x1hihihi, double *x1lohihi, double *x1hilohi, double *x1lolohi,
   double *x1hihilo, double *x1lohilo, double *x1hilolo, double *x1lololo,
   int dim, int dimLog2,
   double *vhihihi, double *vlohihi, double *vhilohi, double *vlolohi,
   double *vhihilo, double *vlohilo, double *vhilolo, double *vlololo,
   double *betahihihi, double *betalohihi,
   double *betahilohi, double *betalolohi,
   double *betahihilo, double *betalohilo,
   double *betahilolo, double *betalololo )
{
   const int j = threadIdx.x;

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

   bool stopflag = false;
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;
   double muhihihi,mulohihi,muhilohi,mulolohi;
   double muhihilo,mulohilo,muhilolo,mulololo;
   double v0hihihi,v0lohihi,v0hilohi,v0lolohi;
   double v0hihilo,v0lohilo,v0hilolo,v0lololo;
   double v0p2hihihi,v0p2lohihi,v0p2hilohi,v0p2lolohi;
   double v0p2hihilo,v0p2lohilo,v0p2hilolo,v0p2lololo;

   shvhihihi[j] = x1hihihi[j];    // reading of vector into shared memory
   shvlohihi[j] = x1lohihi[j];
   shvhilohi[j] = x1hilohi[j];
   shvlolohi[j] = x1lolohi[j];
   shvhihilo[j] = x1hihilo[j];
   shvlohilo[j] = x1lohilo[j];
   shvhilolo[j] = x1hilolo[j];
   shvlololo[j] = x1lololo[j];
   // prd[j] = shv[j]*shv[j];   // for the 2-norm computation
   __syncthreads();
   odg_sqr( shvhihihi[j], shvlohihi[j], shvhilohi[j], shvlolohi[j],
            shvhihilo[j], shvlohilo[j], shvhilolo[j], shvlololo[j],
           &prdhihihi[j],&prdlohihi[j],&prdhilohi[j],&prdlolohi[j],
           &prdhihilo[j],&prdlohilo[j],&prdhilolo[j],&prdlololo[j]);

   __syncthreads();
   vhihihi[j+1] = shvhihihi[j];   // copies x to v, in case beta is zero
   vlohihi[j+1] = shvlohihi[j];
   vhilohi[j+1] = shvhilohi[j];
   vlolohi[j+1] = shvlolohi[j];
   vhihilo[j+1] = shvhihilo[j]; 
   vlohilo[j+1] = shvlohilo[j];
   vhilolo[j+1] = shvhilolo[j];
   vlololo[j+1] = shvlololo[j];
   __syncthreads();
   if(j == 0)
   {
      vhihihi[0] = 1.0;
      vlohihi[0] = 0.0;
      vhilohi[0] = 0.0;
      vlolohi[0] = 0.0;
      vhihilo[0] = 0.0;
      vlohilo[0] = 0.0;
      vhilolo[0] = 0.0;
      vlololo[0] = 0.0;
   }
   __syncthreads();
   int powTwo = 1;                          // sum reduction
   for(int k=0; k < dimLog2; k++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < dim) // prd[j] = prd[j] + prd[j+powTwo];
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
   {                                       // prd[0] is sigma of house
      if((prdhihihi[0] == 0.0) && (prdlohihi[0] == 0.0) &&
         (prdhilohi[0] == 0.0) && (prdlolohi[0] == 0.0) &&
         (prdhihilo[0] == 0.0) && (prdlohilo[0] == 0.0) &&
         (prdhilolo[0] == 0.0) && (prdlololo[0] == 0.0)) 
      {
         *betahihihi = 0.0; *betalohihi = 0.0;
         *betahilohi = 0.0; *betalolohi = 0.0;
         *betahihilo = 0.0; *betalohilo = 0.0;
         *betahilolo = 0.0; *betalololo = 0.0;
         stopflag = true;
      }
   }
   __syncthreads();
   if(stopflag) return;                    // case when sigma is zero
   __syncthreads();
   if(j == 0)                              // thread zero sets beta
   {
      // mu = sqrt((*x0)*(*x0) + prd[0]);
      odg_sqr( *x0hihihi, *x0lohihi, *x0hilohi, *x0lolohi,
               *x0hihilo, *x0lohilo, *x0hilolo, *x0lololo,
              &acchihihi,&acclohihi,&acchilohi,&acclolohi,
              &acchihilo,&acclohilo,&acchilolo,&acclololo);
      odg_inc(&acchihihi,  &acclohihi,  &acchilohi,  &acclolohi,
              &acchihilo,  &acclohilo,  &acchilolo,  &acclololo,
               prdhihihi[0],prdlohihi[0],prdhilohi[0],prdlolohi[0],
               prdhihilo[0],prdlohilo[0],prdhilolo[0],prdlololo[0]);
      odg_sqrt(acchihihi,acclohihi,acchilohi,acclolohi,
               acchihilo,acclohilo,acchilolo,acclololo,
               &muhihihi,&mulohihi,&muhilohi,&mulolohi,
               &muhihilo,&mulohilo,&muhilolo,&mulololo);

      if(*x0hihihi <= 0.0)
      {
         // v0 = *x0 - mu;
         odg_sub(*x0hihihi,*x0lohihi,*x0hilohi,*x0lolohi,
                 *x0hihilo,*x0lohilo,*x0hilolo,*x0lololo,
                  muhihihi, mulohihi, muhilohi, mulolohi,
                  muhihilo, mulohilo, muhilolo, mulololo,
                 &v0hihihi,&v0lohihi,&v0hilohi,&v0lolohi,
                 &v0hihilo,&v0lohilo,&v0hilolo,&v0lololo);
      }
      else
      {
         // v0 = -prd[0]/(*x0 + mu);
         odg_add( *x0hihihi, *x0lohihi, *x0hilohi, *x0lolohi,
                  *x0hihilo, *x0lohilo, *x0hilolo, *x0lololo,
                   muhihihi,  mulohihi,  muhilohi,  mulolohi,
                   muhihilo,  mulohilo,  muhilolo,  mulololo,
                 &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                 &acchihilo,&acclohilo,&acchilolo,&acclololo);
         odg_div(prdhihihi[0],prdlohihi[0],prdhilohi[0],prdlolohi[0],
                 prdhihilo[0],prdlohilo[0],prdhilolo[0],prdlololo[0],
                 acchihihi,   acclohihi,   acchilohi,   acclolohi,
                 acchihilo,   acclohilo,   acchilolo,   acclololo,
                 &v0hihihi,   &v0lohihi,   &v0hilohi,   &v0lolohi,
                 &v0hihilo,   &v0lohilo,   &v0hilolo,   &v0lololo);
         odg_minus(&v0hihihi,&v0lohihi,&v0hilohi,&v0lolohi,
                   &v0hihilo,&v0lohilo,&v0hilolo,&v0lololo);
      }
      // v0p2 = v0*v0;
      odg_sqr(   v0hihihi,   v0lohihi,   v0hilohi,   v0lolohi,
                 v0hihilo,   v0lohilo,   v0hilolo,   v0lololo,
              &v0p2hihihi,&v0p2lohihi,&v0p2hilohi,&v0p2lolohi,
              &v0p2hihilo,&v0p2lohilo,&v0p2hilolo,&v0p2lololo);
      // *beta = 2.0*v0p2/(prd[0] + v0p2);
      odg_add( prdhihihi[0],prdlohihi[0],prdhilohi[0],prdlolohi[0],
               prdhihilo[0],prdlohilo[0],prdhilolo[0],prdlololo[0],
              v0p2hihihi,  v0p2lohihi,  v0p2hilohi,  v0p2lolohi,
              v0p2hihilo,  v0p2lohilo,  v0p2hilolo,  v0p2lololo,
              &acchihihi,  &acclohihi,  &acchilohi,  &acclolohi,
              &acchihilo,  &acclohilo,  &acchilolo,  &acclololo);
      odg_div(v0p2hihihi,v0p2lohihi,v0p2hilohi,v0p2lolohi,
              v0p2hihilo,v0p2lohilo,v0p2hilolo,v0p2lololo,
               acchihihi, acclohihi, acchilohi, acclolohi,
               acchihilo, acclohilo, acchilolo, acclololo,
              betahihihi,betalohihi,betahilohi,betalolohi,
              betahihilo,betalohilo,betahilolo,betalololo);
      odg_mlt_d(betahihihi,betalohihi,betahilohi,betalolohi,
                betahihilo,betalohilo,betahilolo,betalololo,2.0);
      prdhihihi[0] = v0hihihi;
      prdlohihi[0] = v0lohihi;
      prdhilohi[0] = v0hilohi;
      prdlolohi[0] = v0lolohi;
      prdhihilo[0] = v0hihilo;
      prdlohilo[0] = v0lohilo;
      prdhilolo[0] = v0hilolo;
      prdlololo[0] = v0lololo;          // v0 needed for normalization
   }
   __syncthreads();
   // shv[j] = shv[j]/prd[0];
   vhihihi[j+1] = 0.0;
   vlohihi[j+1] = 0.0;
   vhilohi[j+1] = 0.0;
   vlolohi[j+1] = 0.0;
   vhihilo[j+1] = 0.0;
   vlohilo[j+1] = 0.0;
   vhilolo[j+1] = 0.0;
   vlololo[j+1] = 0.0;
   double checksum
      = *betahihihi + *betalohihi + *betahilohi + *betalolohi
      + *betahihilo + *betalohilo + *betahilolo + *betalololo;
   __syncthreads();
   if(1.0 + checksum != 1.0)
   {
      odg_div(shvhihihi[j],shvlohihi[j],shvhilohi[j],shvlolohi[j],
              shvhihilo[j],shvlohilo[j],shvhilolo[j],shvlololo[j],
              prdhihihi[0],prdlohihi[0],prdhilohi[0],prdlolohi[0],
              prdhihilo[0],prdlohilo[0],prdhilolo[0],prdlololo[0],
             &acchihihi,  &acclohihi,  &acchilohi,  &acclolohi,
             &acchihilo,  &acclohilo,  &acchilolo,  &acclololo);
      __syncthreads();
      vhihihi[j+1] = acchihihi;
      vlohihi[j+1] = acclohihi;
      vhilohi[j+1] = acchilohi;
      vlolohi[j+1] = acclolohi;
      vhihilo[j+1] = acchihilo;
      vlohilo[j+1] = acclohilo;
      vhilolo[j+1] = acchilolo;
      vlololo[j+1] = acclololo;
   }
   __syncthreads();
   if(j == 0)
   {
      vhihihi[0] = 1.0;
      vlohihi[0] = 0.0;
      vhilohi[0] = 0.0;
      vlolohi[0] = 0.0;
      vhihilo[0] = 0.0;
      vlohilo[0] = 0.0;
      vhilolo[0] = 0.0;
      vlololo[0] = 0.0;
   }
}

__global__ void cmplx8_small_house
 ( double *x0rehihihi, double *x0relohihi,
   double *x0rehilohi, double *x0relolohi,
   double *x0rehihilo, double *x0relohilo,
   double *x0rehilolo, double *x0relololo,
   double *x0imhihihi, double *x0imlohihi,
   double *x0imhilohi, double *x0imlolohi,
   double *x0imhihilo, double *x0imlohilo,
   double *x0imhilolo, double *x0imlololo,
   double *x1rehihihi, double *x1relohihi,
   double *x1rehilohi, double *x1relolohi,
   double *x1rehihilo, double *x1relohilo,
   double *x1rehilolo, double *x1relololo,
   double *x1imhihihi, double *x1imlohihi,
   double *x1imhilohi, double *x1imlolohi,
   double *x1imhihilo, double *x1imlohilo,
   double *x1imhilolo, double *x1imlololo,
   int dim, int dimLog2,
   double *vrehihihi, double *vrelohihi, double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo, double *vrehilolo, double *vrelololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo, double *vimhilolo, double *vimlololo,
   double *betahihihi, double *betalohihi,
   double *betahilohi, double *betalolohi,
   double *betahihilo, double *betalohilo,
   double *betahilolo, double *betalololo )
{
   const int j = threadIdx.x;

   __shared__ double shvrehihihi[cod_shmemsize];
   __shared__ double shvrelohihi[cod_shmemsize];
   __shared__ double shvrehilohi[cod_shmemsize];
   __shared__ double shvrelolohi[cod_shmemsize];
   __shared__ double shvrehihilo[cod_shmemsize];
   __shared__ double shvrelohilo[cod_shmemsize];
   __shared__ double shvrehilolo[cod_shmemsize];
   __shared__ double shvrelololo[cod_shmemsize];
   __shared__ double shvimhihihi[cod_shmemsize];
   __shared__ double shvimlohihi[cod_shmemsize];
   __shared__ double shvimhilohi[cod_shmemsize];
   __shared__ double shvimlolohi[cod_shmemsize];
   __shared__ double shvimhihilo[cod_shmemsize];
   __shared__ double shvimlohilo[cod_shmemsize];
   __shared__ double shvimhilolo[cod_shmemsize];
   __shared__ double shvimlololo[cod_shmemsize];
   __shared__ double prdhihihi[cod_shmemsize];
   __shared__ double prdlohihi[cod_shmemsize];
   __shared__ double prdhilohi[cod_shmemsize];
   __shared__ double prdlolohi[cod_shmemsize];
   __shared__ double prdhihilo[cod_shmemsize];
   __shared__ double prdlohilo[cod_shmemsize];
   __shared__ double prdhilolo[cod_shmemsize];
   __shared__ double prdlololo[cod_shmemsize];
   __shared__ double v0parts[16];

   bool stopflag = false;
   double muhihihi,mulohihi,muhilohi,mulolohi;
   double muhihilo,mulohilo,muhilolo,mulololo;
   double v0rehihihi,v0relohihi,v0rehilohi,v0relolohi;
   double v0rehihilo,v0relohilo,v0rehilolo,v0relololo;
   double v0imhihihi,v0imlohihi,v0imhilohi,v0imlolohi;
   double v0imhihilo,v0imlohilo,v0imhilolo,v0imlololo;
   double x0radhihihi,x0radlohihi,x0radhilohi,x0radlolohi;
   double x0radhihilo,x0radlohilo,x0radhilolo,x0radlololo;
   double sqrx0hihihi,sqrx0lohihi,sqrx0hilohi,sqrx0lolohi;
   double sqrx0hihilo,sqrx0lohilo,sqrx0hilolo,sqrx0lololo;
   double sqrv0hihihi,sqrv0lohihi,sqrv0hilohi,sqrv0lolohi;
   double sqrv0hihilo,sqrv0lohilo,sqrv0hilolo,sqrv0lololo;
   double inv0rehihihi,inv0relohihi,inv0rehilohi,inv0relolohi;
   double inv0rehihilo,inv0relohilo,inv0rehilolo,inv0relololo;
   double inv0imhihihi,inv0imlohihi,inv0imhilohi,inv0imlolohi;
   double inv0imhihilo,inv0imlohilo,inv0imhilolo,inv0imlololo;
   double zrehihihi,zrelohihi,zrehilohi,zrelolohi;
   double zrehihilo,zrelohilo,zrehilolo,zrelololo;
   double zimhihihi,zimlohihi,zimhilohi,zimlolohi;
   double zimhihilo,zimlohilo,zimhilolo,zimlololo;
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   shvrehihihi[j] = x1rehihihi[j];  // reading of vector into shared memory
   shvrelohihi[j] = x1relohihi[j];
   shvrehilohi[j] = x1rehilohi[j];
   shvrelolohi[j] = x1relolohi[j];
   shvrehihilo[j] = x1rehihilo[j];
   shvrelohilo[j] = x1relohilo[j];
   shvrehilolo[j] = x1rehilolo[j];
   shvrelololo[j] = x1relololo[j];
   shvimhihihi[j] = x1imhihihi[j];
   shvimlohihi[j] = x1imlohihi[j];
   shvimhilohi[j] = x1imhilohi[j];
   shvimlolohi[j] = x1imlolohi[j];
   shvimhihilo[j] = x1imhihilo[j];
   shvimlohilo[j] = x1imlohilo[j];
   shvimhilolo[j] = x1imhilolo[j];
   shvimlololo[j] = x1imlololo[j];
   // prd[j] = shv[j]*shv[j];   // for the 2-norm computation
   // prd[j] = shvre[j]*shvre[j] + shvim[j]*shvim[j];
   __syncthreads();
   odg_sqr(shvrehihihi[j],shvrelohihi[j],shvrehilohi[j],shvrelolohi[j],
           shvrehihilo[j],shvrelohilo[j],shvrehilolo[j],shvrelololo[j],
            &prdhihihi[j], &prdlohihi[j], &prdhilohi[j], &prdlolohi[j],
            &prdhihilo[j], &prdlohilo[j], &prdhilolo[j], &prdlololo[j]);
   odg_sqr(shvimhihihi[j],shvimlohihi[j],shvimhilohi[j],shvimlolohi[j],
           shvimhihilo[j],shvimlohilo[j],shvimhilolo[j],shvimlololo[j],
            &acchihihi,    &acclohihi,    &acchilohi,    &acclolohi,
            &acchihilo,    &acclohilo,    &acchilolo,    &acclololo);
   odg_inc(&prdhihihi[j],&prdlohihi[j],&prdhilohi[j],&prdlolohi[j],
           &prdhihilo[j],&prdlohilo[j],&prdhilolo[j],&prdlololo[j],
            acchihihi,    acclohihi,    acchilohi,    acclolohi,
            acchihilo,    acclohilo,    acchilolo,    acclololo);
   __syncthreads();
   vrehihihi[j+1] = shvrehihihi[j];  // copies x to v, in case beta is zero
   vrelohihi[j+1] = shvrelohihi[j];
   vrehilohi[j+1] = shvrehilohi[j];
   vrelolohi[j+1] = shvrelolohi[j];
   vrehihilo[j+1] = shvrehihilo[j];
   vrelohilo[j+1] = shvrelohilo[j];
   vrehilolo[j+1] = shvrehilolo[j];
   vrelololo[j+1] = shvrelololo[j];
   vimhihihi[j+1] = shvimhihihi[j];
   vimlohihi[j+1] = shvimlohihi[j];
   vimhilohi[j+1] = shvimhilohi[j];
   vimlolohi[j+1] = shvimlolohi[j];
   vimhihilo[j+1] = shvimhihilo[j];
   vimlohilo[j+1] = shvimlohilo[j];
   vimhilolo[j+1] = shvimhilolo[j];
   vimlololo[j+1] = shvimlololo[j];
   __syncthreads();
   if(j == 0)
   {
      vrehihihi[0] = 1.0; vrelohihi[0] = 0.0;
      vrehilohi[0] = 0.0; vrelolohi[0] = 0.0;
      vrehihilo[0] = 0.0; vrelohilo[0] = 0.0;
      vrehilolo[0] = 0.0; vrelololo[0] = 0.0;
   }
   __syncthreads();
   if(j == 0)
   {
      vimhihihi[0] = 0.0; vimlohihi[0] = 0.0;
      vimhilohi[0] = 0.0; vimlolohi[0] = 0.0;
      vimhihilo[0] = 0.0; vimlohilo[0] = 0.0;
      vimhilolo[0] = 0.0; vimlololo[0] = 0.0;
   }
   __syncthreads();
   int powTwo = 1;                          // sum reduction
   for(int k=0; k < dimLog2; k++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < dim)     // prd[j] = prd[j] + prd[j+powTwo];
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
   {                                        // prd[0] is sigma of house
      if((prdhihihi[0] == 0.0) && (prdlohihi[0] == 0.0) &&
         (prdhilohi[0] == 0.0) && (prdlolohi[0] == 0.0) &&
         (prdhihilo[0] == 0.0) && (prdlohilo[0] == 0.0) &&
         (prdhilolo[0] == 0.0) && (prdlololo[0] == 0.0))
      {
         *betahihihi = 0.0; *betalohihi = 0.0;
         *betahilohi = 0.0; *betalolohi = 0.0;
         *betahihilo = 0.0; *betalohilo = 0.0;
         *betahilolo = 0.0; *betalololo = 0.0; stopflag = true;
      }
   }
   __syncthreads();
   if(stopflag) return;                    // case when sigma is zero
   __syncthreads();
   if(j == 0)                              // thread zero sets beta
   {
      // sqrx0 = (*x0re)*(*x0re) + (*x0im)*(*x0im);
      odg_sqr( *x0rehihihi, *x0relohihi, *x0rehilohi, *x0relolohi,
               *x0rehihilo, *x0relohilo, *x0rehilolo, *x0relololo,
              &sqrx0hihihi,&sqrx0lohihi,&sqrx0hilohi,&sqrx0lolohi,
              &sqrx0hihilo,&sqrx0lohilo,&sqrx0hilolo,&sqrx0lololo);
      odg_sqr(*x0imhihihi,*x0imlohihi,*x0imhilohi,*x0imlolohi,
              *x0imhihilo,*x0imlohilo,*x0imhilolo,*x0imlololo,
               &acchihihi, &acclohihi, &acchilohi, &acclolohi,
               &acchihilo, &acclohilo, &acchilolo, &acclololo);
      odg_inc(&sqrx0hihihi,&sqrx0lohihi,&sqrx0hilohi,&sqrx0lolohi,
              &sqrx0hihilo,&sqrx0lohilo,&sqrx0hilolo,&sqrx0lololo,
                 acchihihi,   acclohihi,   acchilohi,   acclolohi,
                 acchihilo,   acclohilo,   acchilolo,   acclololo);
      // x0rad = sqrt(sqrx0);
      odg_sqrt( sqrx0hihihi, sqrx0lohihi, sqrx0hilohi, sqrx0lolohi,
                sqrx0hihilo, sqrx0lohilo, sqrx0hilolo, sqrx0lololo,
               &x0radhihihi,&x0radlohihi,&x0radhilohi,&x0radlolohi,
               &x0radhihilo,&x0radlohilo,&x0radhilolo,&x0radlololo);
      // mu = sqrt(sqrx0 + prd[0]);
      odg_inc(&sqrx0hihihi,&sqrx0lohihi,&sqrx0hilohi,&sqrx0lolohi,
              &sqrx0hihilo,&sqrx0lohilo,&sqrx0hilolo,&sqrx0lololo,
                 prdhihihi[0],prdlohihi[0],prdhilohi[0],prdlolohi[0],
                 prdhihilo[0],prdlohilo[0],prdhilolo[0],prdlololo[0]);
      odg_sqrt(sqrx0hihihi,sqrx0lohihi,sqrx0hilohi,sqrx0lolohi,
               sqrx0hihilo,sqrx0lohilo,sqrx0hilolo,sqrx0lololo,
                 &muhihihi,  &mulohihi,  &muhilohi,  &mulolohi,
                 &muhihilo,  &mulohilo,  &muhilolo,  &mulololo);

      if((x0radhihihi == 0.0) && (x0radlohihi == 0.0) &&
         (x0radhilohi == 0.0) && (x0radlolohi == 0.0) &&
         (x0radhihilo == 0.0) && (x0radlohilo == 0.0) &&
         (x0radhilolo == 0.0) && (x0radlololo == 0.0))
      {
         v0rehihihi = muhihihi; v0relohihi = mulohihi;
         v0rehilohi = muhilohi; v0relolohi = mulolohi;
         v0rehihilo = muhihilo; v0relohilo = mulohilo;
         v0rehilolo = muhilolo; v0relololo = mulololo;

         odg_minus(&v0rehihihi,&v0relohihi,&v0rehilohi,&v0relolohi,
                   &v0rehihilo,&v0relohilo,&v0rehilolo,&v0relololo);

         v0imhihihi = 0.0; v0imlohihi = 0.0;
         v0imhilohi = 0.0; v0imlolohi = 0.0;
         v0imhihilo = 0.0; v0imlohilo = 0.0;
         v0imhilolo = 0.0; v0imlololo = 0.0;
      }
      else
      {
         // mu = mu/x0rad;
         odg_div(   muhihihi,   mulohihi,   muhilohi,   mulolohi,
                    muhihilo,   mulohilo,   muhilolo,   mulololo,
                 x0radhihihi,x0radlohihi,x0radhilohi,x0radlolohi,
                 x0radhihilo,x0radlohilo,x0radhilolo,x0radlololo,
                  &acchihihi, &acclohihi, &acchilohi, &acclolohi,
                  &acchihilo, &acclohilo, &acchilolo, &acclololo);
         muhihihi = acchihihi; mulohihi = acclohihi;
         muhilohi = acchilohi; mulolohi = acclolohi;
         muhihilo = acchihilo; mulohilo = acclohilo;
         muhilolo = acchilolo; mulololo = acclololo;
         // v0re = (*x0re) - mu*(*x0re);
         odg_mul(  muhihihi,     mulohihi,     muhilohi,     mulolohi,
                   muhihilo,     mulohilo,     muhilolo,     mulololo,
                 x0rehihihi[0],x0relohihi[0],x0rehilohi[0],x0relolohi[0],
                 x0rehihilo[0],x0relohilo[0],x0rehilolo[0],x0relololo[0],
                 &acchihihi,   &acclohihi,   &acchilohi,   &acclolohi,
                 &acchihilo,   &acclohilo,   &acchilolo,   &acclololo);
         odg_sub(x0rehihihi[0],x0relohihi[0],x0rehilohi[0],x0relolohi[0],
                 x0rehihilo[0],x0relohilo[0],x0rehilolo[0],x0relololo[0],
                  acchihihi,    acclohihi,    acchilohi,    acclolohi,
                  acchihilo,    acclohilo,    acchilolo,    acclololo,
                &v0rehihihi,  &v0relohihi,  &v0rehilohi,  &v0relolohi,
                &v0rehihilo,  &v0relohilo,  &v0rehilolo,  &v0relololo);
         // v0im = (*x0im) - mu*(*x0im);
         odg_mul(  muhihihi,     mulohihi,     muhilohi,     mulolohi,
                   muhihilo,     mulohilo,     muhilolo,     mulololo,
                 x0imhihihi[0],x0imlohihi[0],x0imhilohi[0],x0imlolohi[0],
                 x0imhihilo[0],x0imlohilo[0],x0imhilolo[0],x0imlololo[0],
                 &acchihihi,   &acclohihi,   &acchilohi,   &acclolohi,
                 &acchihilo,   &acclohilo,   &acchilolo,   &acclololo);
         odg_sub(x0imhihihi[0],x0imlohihi[0],x0imhilohi[0],x0imlolohi[0],
                 x0imhihilo[0],x0imlohilo[0],x0imhilolo[0],x0imlololo[0],
                  acchihihi,    acclohihi,    acchilohi,    acclolohi,
                  acchihilo,    acclohilo,    acchilolo,    acclololo,
                &v0imhihihi,  &v0imlohihi,  &v0imhilohi,  &v0imlolohi,
                &v0imhihilo,  &v0imlohilo,  &v0imhilolo,  &v0imlololo);
      }
      // sqrv0 = v0re*v0re + v0im*v0im;
      odg_sqr(  v0rehihihi,  v0relohihi,  v0rehilohi,  v0relolohi,
                v0rehihilo,  v0relohilo,  v0rehilolo,  v0relololo,
              &sqrv0hihihi,&sqrv0lohihi,&sqrv0hilohi,&sqrv0lolohi,
              &sqrv0hihilo,&sqrv0lohilo,&sqrv0hilolo,&sqrv0lololo);
      odg_sqr(v0imhihihi,v0imlohihi,v0imhilohi,v0imlolohi,
              v0imhihilo,v0imlohilo,v0imhilolo,v0imlololo,
              &acchihihi,&acclohihi,&acchilohi,&acclolohi,
              &acchihilo,&acclohilo,&acchilolo,&acclololo);
      odg_inc(&sqrv0hihihi,&sqrv0lohihi,&sqrv0hilohi,&sqrv0lolohi,
              &sqrv0hihilo,&sqrv0lohilo,&sqrv0hilolo,&sqrv0lololo,
                 acchihihi,   acclohihi,   acchilohi,   acclolohi,
                 acchihilo,   acclohilo,   acchilolo,   acclololo);
      // *beta = 2.0*sqrv0/(prd[0] + sqrv0);
      odg_add(  prdhihihi[0],prdlohihi[0],prdhilohi[0],prdlolohi[0],
                prdhihilo[0],prdlohilo[0],prdhilolo[0],prdlololo[0],
              sqrv0hihihi, sqrv0lohihi, sqrv0hilohi, sqrv0lolohi,
              sqrv0hihilo, sqrv0lohilo, sqrv0hilolo, sqrv0lololo,
               &acchihihi,  &acclohihi,  &acchilohi,  &acclolohi,
               &acchihilo,  &acclohilo,  &acchilolo,  &acclololo);
      odg_div(sqrv0hihihi,sqrv0lohihi,sqrv0hilohi,sqrv0lolohi,
              sqrv0hihilo,sqrv0lohilo,sqrv0hilolo,sqrv0lololo,
                acchihihi,  acclohihi,  acchilohi,  acclolohi,
                acchihilo,  acclohilo,  acchilolo,  acclololo,
               betahihihi, betalohihi, betahilohi, betalolohi,
               betahihilo, betalohilo, betahilolo, betalololo);
      odg_mlt_d(betahihihi,betalohihi,betahilohi,betalolohi,
                betahihilo,betalohilo,betahilolo,betalololo,2.0);

      prdhihihi[0] = sqrv0hihihi;     // sqrv0 needed for normalization
      prdlohihi[0] = sqrv0lohihi;
      prdhilohi[0] = sqrv0hilohi;
      prdlolohi[0] = sqrv0lolohi;
      prdhihilo[0] = sqrv0hihilo; 
      prdlohilo[0] = sqrv0lohilo;
      prdhilolo[0] = sqrv0hilolo;
      prdlololo[0] = sqrv0lololo;
      v0parts[0] = v0rehihihi;      // share v0re with all threads
      v0parts[1] = v0relohihi; 
      v0parts[2] = v0rehilohi; 
      v0parts[3] = v0relolohi; 
      v0parts[4] = v0rehihilo; 
      v0parts[5] = v0relohilo; 
      v0parts[6] = v0rehilolo; 
      v0parts[7] = v0relololo; 
      v0parts[8] = v0imhihihi;      // share v0im with all threads
      v0parts[9] = v0imlohihi;
      v0parts[10] = v0imhilohi;
      v0parts[11] = v0imlolohi;
      v0parts[12] = v0imhihilo; 
      v0parts[13] = v0imlohilo;
      v0parts[14] = v0imhilolo;
      v0parts[15] = v0imlololo;
   }
   __syncthreads(); // important synchronization!
   // inv0re = v0parts[0]/prd[0];               // real part of 1/v[0]
   vrehihihi[j+1] = 0.0;
   vrelohihi[j+1] = 0.0;
   vrehilohi[j+1] = 0.0;
   vrelolohi[j+1] = 0.0;
   vrehihilo[j+1] = 0.0;
   vrelohilo[j+1] = 0.0;
   vrehilolo[j+1] = 0.0;
   vrelololo[j+1] = 0.0;
   vimhihihi[j+1] = 0.0;
   vimlohihi[j+1] = 0.0;
   vimhilohi[j+1] = 0.0;
   vimlolohi[j+1] = 0.0;
   vimhihilo[j+1] = 0.0;
   vimlohilo[j+1] = 0.0;
   vimhilolo[j+1] = 0.0;
   vimlololo[j+1] = 0.0;

   double checksum 
      = prdhihihi[0] + prdlohihi[0] + prdhilohi[0] + prdlolohi[0]
      + prdhihilo[0] + prdlohilo[0] + prdhilolo[0] + prdlololo[0];
   __syncthreads();
   if(1.0 + checksum != 1.0)
   {
      odg_div(v0parts[0],v0parts[1],v0parts[2],v0parts[3],
              v0parts[4],v0parts[5],v0parts[6],v0parts[7],
                  prdhihihi[0], prdlohihi[0], prdhilohi[0], prdlolohi[0],
                  prdhihilo[0], prdlohilo[0], prdhilolo[0], prdlololo[0],
              &inv0rehihihi,&inv0relohihi,&inv0rehilohi,&inv0relolohi,
              &inv0rehihilo,&inv0relohilo,&inv0rehilolo,&inv0relololo);
      // inv0im = -v0parts[1]/prd[0];              // imag part of 1/v[0]
      odg_div(v0parts[8],v0parts[9],v0parts[10],v0parts[11],
              v0parts[12],v0parts[13],v0parts[14],v0parts[15],
                  prdhihihi[0], prdlohihi[0], prdhilohi[0], prdlolohi[0],
                  prdhihilo[0], prdlohilo[0], prdhilolo[0], prdlololo[0],
              &inv0imhihihi,&inv0imlohihi,&inv0imhilohi,&inv0imlolohi,
              &inv0imhihilo,&inv0imlohilo,&inv0imhilolo,&inv0imlololo);
      odg_minus(&inv0imhihihi,&inv0imlohihi,&inv0imhilohi,&inv0imlolohi,
                &inv0imhihilo,&inv0imlohilo,&inv0imhilolo,&inv0imlololo);
      // zre = shvre[j]*inv0re - shvim[j]*inv0im;  // real part of v[j]/v[0]
      __syncthreads();
      odg_mul( shvrehihihi[j],shvrelohihi[j],shvrehilohi[j],shvrelolohi[j],
               shvrehihilo[j],shvrelohilo[j],shvrehilolo[j],shvrelololo[j],
              inv0rehihihi,  inv0relohihi,  inv0rehilohi,  inv0relolohi,
              inv0rehihilo,  inv0relohilo,  inv0rehilolo,  inv0relololo,
                &zrehihihi,    &zrelohihi,    &zrehilohi,    &zrelolohi,
                &zrehihilo,    &zrelohilo,    &zrehilolo,    &zrelololo);
      odg_mul( shvimhihihi[j],shvimlohihi[j],shvimhilohi[j],shvimlolohi[j],
               shvimhihilo[j],shvimlohilo[j],shvimhilolo[j],shvimlololo[j],
              inv0imhihihi,  inv0imlohihi,  inv0imhilohi,  inv0imlolohi,
              inv0imhihilo,  inv0imlohilo,  inv0imhilolo,  inv0imlololo,
                &acchihihi,    &acclohihi,    &acchilohi,    &acclolohi,
                &acchihilo,    &acclohilo,    &acchilolo,    &acclololo);
      odg_dec(&zrehihihi,&zrelohihi,&zrehilohi,&zrelolohi,
              &zrehihilo,&zrelohilo,&zrehilolo,&zrelololo,
               acchihihi, acclohihi, acchilohi, acclolohi,
               acchihilo, acclohilo, acchilolo, acclololo);
      // zim = shvim[j]*inv0re + shvre[j]*inv0im;  // imag part of v[j]/v[0]
      odg_mul( shvimhihihi[j],shvimlohihi[j],shvimhilohi[j],shvimlolohi[j],
               shvimhihilo[j],shvimlohilo[j],shvimhilolo[j],shvimlololo[j],
              inv0rehihihi,  inv0relohihi,  inv0rehilohi,  inv0relolohi,
              inv0rehihilo,  inv0relohilo,  inv0rehilolo,  inv0relololo,
                &zimhihihi,    &zimlohihi,    &zimhilohi,    &zimlolohi,
                &zimhihilo,    &zimlohilo,    &zimhilolo,    &zimlololo);
      odg_mul( shvrehihihi[j],shvrelohihi[j],shvrehilohi[j],shvrelolohi[j],
               shvrehihilo[j],shvrelohilo[j],shvrehilolo[j],shvrelololo[j],
              inv0imhihihi,  inv0imlohihi,  inv0imhilohi,  inv0imlolohi,
              inv0imhihilo,  inv0imlohilo,  inv0imhilolo,  inv0imlololo,
                &acchihihi,    &acclohihi,    &acchilohi,    &acclolohi,
                &acchihilo,    &acclohilo,    &acchilolo,    &acclololo);
      odg_inc(&zimhihihi,&zimlohihi,&zimhilohi,&zimlolohi,
              &zimhihilo,&zimlohilo,&zimhilolo,&zimlololo,
               acchihihi, acclohihi, acchilohi, acclolohi,
               acchihilo, acclohilo, acchilolo, acclololo);
      __syncthreads();
      vrehihihi[j+1] = zrehihihi;
      vrelohihi[j+1] = zrelohihi;
      vrehilohi[j+1] = zrehilohi;
      vrelolohi[j+1] = zrelolohi;
      vrehihilo[j+1] = zrehihilo;
      vrelohilo[j+1] = zrelohilo;
      vrehilolo[j+1] = zrehilolo;
      vrelololo[j+1] = zrelololo;
      vimhihihi[j+1] = zimhihihi;
      vimlohihi[j+1] = zimlohihi;
      vimhilohi[j+1] = zimhilohi;
      vimlolohi[j+1] = zimlolohi;
      vimhihilo[j+1] = zimhihilo;
      vimlohilo[j+1] = zimlohilo;
      vimhilolo[j+1] = zimhilolo;
      vimlololo[j+1] = zimlololo;
   }
   __syncthreads();
   if(j == 0)
   {
      vrehihihi[0] = 1.0; vrelohihi[0] = 0.0;
      vrehilohi[0] = 0.0; vrelolohi[0] = 0.0;
      vrehihilo[0] = 0.0; vrelohilo[0] = 0.0;
      vrehilolo[0] = 0.0; vrelololo[0] = 0.0;
   }
   __syncthreads();
   if(j == 0)
   {
      vimhihihi[0] = 0.0; vimlohihi[0] = 0.0;
      vimhilohi[0] = 0.0; vimlolohi[0] = 0.0;
      vimhihilo[0] = 0.0; vimlohilo[0] = 0.0;
      vimhilolo[0] = 0.0; vimlololo[0] = 0.0;
   }
   __syncthreads();
}

__global__ void dbl8_large_sum_of_squares
 ( double *vhihihi, double *vlohihi, double *vhilohi, double *vlolohi,
   double *vhihilo, double *vlohilo, double *vhilolo, double *vlololo,
   double *sumshihihi, double *sumslohihi,
   double *sumshilohi, double *sumslolohi,
   double *sumshihilo, double *sumslohilo,
   double *sumshilolo, double *sumslololo, int dim, int BS, int BSLog2 )
{
   const int i = blockIdx.x;
   const int j = threadIdx.x;
   const int k = i*BS + j;

   __shared__ double shvhihihi[inner_od_shmemsize];
   __shared__ double shvlohihi[inner_od_shmemsize];
   __shared__ double shvhilohi[inner_od_shmemsize];
   __shared__ double shvlolohi[inner_od_shmemsize];
   __shared__ double shvhihilo[inner_od_shmemsize];
   __shared__ double shvlohilo[inner_od_shmemsize];
   __shared__ double shvhilolo[inner_od_shmemsize];
   __shared__ double shvlololo[inner_od_shmemsize];
   __shared__ double prdhihihi[inner_od_shmemsize];
   __shared__ double prdlohihi[inner_od_shmemsize];
   __shared__ double prdhilohi[inner_od_shmemsize];
   __shared__ double prdlolohi[inner_od_shmemsize];
   __shared__ double prdhihilo[inner_od_shmemsize];
   __shared__ double prdlohilo[inner_od_shmemsize];
   __shared__ double prdhilolo[inner_od_shmemsize];
   __shared__ double prdlololo[inner_od_shmemsize];

   shvhihihi[j] = vhihihi[k];
   shvlohihi[j] = vlohihi[k];
   shvhilohi[j] = vhilohi[k];
   shvlolohi[j] = vlolohi[k];
   shvhihilo[j] = vhihilo[k];
   shvlohilo[j] = vlohilo[k];
   shvhilolo[j] = vhilolo[k];
   shvlololo[j] = vlololo[k];
   __syncthreads();
   if(k >= dim)
   {
      shvhihihi[j] = 0.0;
      shvlohihi[j] = 0.0;
      shvhilohi[j] = 0.0;
      shvlolohi[j] = 0.0;
      shvhihilo[j] = 0.0;
      shvlohilo[j] = 0.0;
      shvhilolo[j] = 0.0;
      shvlololo[j] = 0.0;
   }
   __syncthreads();
   odg_sqr( shvhihihi[j], shvlohihi[j], shvhilohi[j], shvlolohi[j],
            shvhihilo[j], shvlohilo[j], shvhilolo[j], shvlololo[j],
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

__global__ void cmplx8_large_sum_of_squares
 ( double *vrehihihi, double *vrelohihi, double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo, double *vrehilolo, double *vrelololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo, double *vimhilolo, double *vimlololo,
   double *sumshihihi, double *sumslohihi,
   double *sumshilohi, double *sumslolohi,
   double *sumshihilo, double *sumslohilo,
   double *sumshilolo, double *sumslololo, int dim, int BS, int BSLog2 )
{
   const int i = blockIdx.x;
   const int j = threadIdx.x;
   const int k = i*BS + j;

   __shared__ double shvrehihihi[inner_od_shmemsize];
   __shared__ double shvrelohihi[inner_od_shmemsize];
   __shared__ double shvrehilohi[inner_od_shmemsize];
   __shared__ double shvrelolohi[inner_od_shmemsize];
   __shared__ double shvrehihilo[inner_od_shmemsize];
   __shared__ double shvrelohilo[inner_od_shmemsize];
   __shared__ double shvrehilolo[inner_od_shmemsize];
   __shared__ double shvrelololo[inner_od_shmemsize];
   __shared__ double shvimhihihi[inner_od_shmemsize];
   __shared__ double shvimlohihi[inner_od_shmemsize];
   __shared__ double shvimhilohi[inner_od_shmemsize];
   __shared__ double shvimlolohi[inner_od_shmemsize];
   __shared__ double shvimhihilo[inner_od_shmemsize];
   __shared__ double shvimlohilo[inner_od_shmemsize];
   __shared__ double shvimhilolo[inner_od_shmemsize];
   __shared__ double shvimlololo[inner_od_shmemsize];
   __shared__ double prdhihihi[inner_od_shmemsize];
   __shared__ double prdlohihi[inner_od_shmemsize];
   __shared__ double prdhilohi[inner_od_shmemsize];
   __shared__ double prdlolohi[inner_od_shmemsize];
   __shared__ double prdhihilo[inner_od_shmemsize];
   __shared__ double prdlohilo[inner_od_shmemsize];
   __shared__ double prdhilolo[inner_od_shmemsize];
   __shared__ double prdlololo[inner_od_shmemsize];

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
   __syncthreads();
   if(k >= dim)
   {
      shvrehihihi[j] = 0.0;
      shvrelohihi[j] = 0.0;
      shvrehilohi[j] = 0.0;
      shvrelolohi[j] = 0.0;
      shvrehihilo[j] = 0.0;
      shvrelohilo[j] = 0.0;
      shvrehilolo[j] = 0.0;
      shvrelololo[j] = 0.0;
      shvimhihihi[j] = 0.0;
      shvimlohihi[j] = 0.0;
      shvimhilohi[j] = 0.0;
      shvimlolohi[j] = 0.0;
      shvimhihilo[j] = 0.0;
      shvimlohilo[j] = 0.0;
      shvimhilolo[j] = 0.0;
      shvimlololo[j] = 0.0;
   }
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;
   
   __syncthreads();
   odg_sqr(shvrehihihi[j],shvrelohihi[j],shvrehilohi[j],shvrelolohi[j],
           shvrehihilo[j],shvrelohilo[j],shvrehilolo[j],shvrelololo[j],
            &prdhihihi[j], &prdlohihi[j], &prdhilohi[j], &prdlolohi[j],
            &prdhihilo[j], &prdlohilo[j], &prdhilolo[j], &prdlololo[j]);
   odg_sqr(shvimhihihi[j],shvimlohihi[j],shvimhilohi[j],shvimlolohi[j],
           shvimhihilo[j],shvimlohilo[j],shvimhilolo[j],shvimlololo[j],
            &acchihihi,    &acclohihi,    &acchilohi,    &acclolohi,
            &acchihilo,    &acclohilo,    &acchilolo,    &acclololo);
   odg_inc(&prdhihihi[j],&prdlohihi[j],&prdhilohi[j],&prdlolohi[j],
           &prdhihilo[j],&prdlohilo[j],&prdhilolo[j],&prdlololo[j],
            acchihihi,    acclohihi,    acchilohi,    acclolohi,
            acchihilo,    acclohilo,    acchilolo,    acclololo);

   __syncthreads();

   int powTwo = 1;                          // sum reduction
   for(int k=0; k < BSLog2; k++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < BS)     // prd[j] = prd[j] + prd[j+powTwo];
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

__global__ void dbl8_sum_accumulator
 ( double *sumshihihi, double *sumslohihi,
   double *sumshilohi, double *sumslolohi,
   double *sumshihilo, double *sumslohilo,
   double *sumshilolo, double *sumslololo,
   int nbsums, int nbsumsLog2,
   double *acchihihi, double *acclohihi,
   double *acchilohi, double *acclolohi,
   double *acchihilo, double *acclohilo,
   double *acchilolo, double *acclololo )
{
   const int j = threadIdx.x;

   __shared__ double shvhihihi[outer_od_shmemsize];
   __shared__ double shvlohihi[outer_od_shmemsize];
   __shared__ double shvhilohi[outer_od_shmemsize];
   __shared__ double shvlolohi[outer_od_shmemsize];
   __shared__ double shvhihilo[outer_od_shmemsize];
   __shared__ double shvlohilo[outer_od_shmemsize];
   __shared__ double shvhilolo[outer_od_shmemsize];
   __shared__ double shvlololo[outer_od_shmemsize];

   shvhihihi[j] = sumshihihi[j];
   shvlohihi[j] = sumslohihi[j];
   shvhilohi[j] = sumshilohi[j];
   shvlolohi[j] = sumslolohi[j];
   shvhihilo[j] = sumshihilo[j];
   shvlohilo[j] = sumslohilo[j];
   shvhilolo[j] = sumshilolo[j];
   shvlololo[j] = sumslololo[j];

   __syncthreads();

   if(j >= nbsums)
   {
      shvhihihi[j] = 0.0;
      shvlohihi[j] = 0.0;
      shvhilohi[j] = 0.0;
      shvlolohi[j] = 0.0;
      shvhihilo[j] = 0.0;
      shvlohilo[j] = 0.0;
      shvhilolo[j] = 0.0;
      shvlololo[j] = 0.0;
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
   __syncthreads();
   if(j == 0)
   {
      *acchihihi = shvhihihi[0];
      *acclohihi = shvlohihi[0];
      *acchilohi = shvhilohi[0];
      *acclolohi = shvlolohi[0];
      *acchihilo = shvhihilo[0];
      *acclohilo = shvlohilo[0];
      *acchilolo = shvhilolo[0];
      *acclololo = shvlololo[0];
   }
}

__global__ void dbl8_normalize
 ( int dim, int szt,
   double *xhihihi, double *xlohihi, double *xhilohi, double *xlolohi,
   double *xhihilo, double *xlohilo, double *xhilolo, double *xlololo,
   double *v0hihihi, double *v0lohihi, double *v0hilohi, double *v0lolohi,
   double *v0hihilo, double *v0lohilo, double *v0hilolo, double *v0lololo,
   double *vhihihi, double *vlohihi, double *vhilohi, double *vlolohi,
   double *vhihilo, double *vlohilo, double *vhilolo, double *vlololo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int idx = bdx*szt + tdx;  // thread tdx scales idx

   __shared__ double shvhihihi[inner_od_shmemsize];
   __shared__ double shvlohihi[inner_od_shmemsize];
   __shared__ double shvhilohi[inner_od_shmemsize];
   __shared__ double shvlolohi[inner_od_shmemsize];
   __shared__ double shvhihilo[inner_od_shmemsize];
   __shared__ double shvlohilo[inner_od_shmemsize];
   __shared__ double shvhilolo[inner_od_shmemsize];
   __shared__ double shvlololo[inner_od_shmemsize];

   shvhihihi[tdx] = xhihihi[idx];
   shvlohihi[tdx] = xlohihi[idx];
   shvhilohi[tdx] = xhilohi[idx];
   shvlolohi[tdx] = xlolohi[idx];
   shvhihilo[tdx] = xhihilo[idx];
   shvlohilo[tdx] = xlohilo[idx];
   shvhilolo[tdx] = xhilolo[idx];
   shvlololo[tdx] = xlololo[idx];
   __syncthreads();

   double resulthihihi,resultlohihi,resulthilohi,resultlolohi;
   double resulthihilo,resultlohilo,resulthilolo,resultlololo;

   // shv[j] = shv[j]/v0;
   odg_div(    shvhihihi[tdx],shvlohihi[tdx],shvhilohi[tdx],shvlolohi[tdx],
               shvhihilo[tdx],shvlohilo[tdx],shvhilolo[tdx],shvlololo[tdx],
                v0hihihi[0],   v0lohihi[0],   v0hilohi[0],   v0lolohi[0],
                v0hihilo[0],   v0lohilo[0],   v0hilolo[0],   v0lololo[0],
           &resulthihihi, &resultlohihi, &resulthilohi, &resultlolohi,
           &resulthihilo, &resultlohilo, &resulthilolo, &resultlololo);

   __syncthreads();
   if(idx < dim)
   {
      vhihihi[idx] = resulthihihi;
      vlohihi[idx] = resultlohihi;
      vhilohi[idx] = resulthilohi;
      vlolohi[idx] = resultlolohi;
      vhihilo[idx] = resulthihilo;
      vlohilo[idx] = resultlohilo;
      vhilolo[idx] = resulthilolo;
      vlololo[idx] = resultlololo;
   }
}

__global__ void cmplx8_normalize
 ( int dim, int szt,
   double *xrehihihi, double *xrelohihi, double *xrehilohi, double *xrelolohi,
   double *xrehihilo, double *xrelohilo, double *xrehilolo, double *xrelololo,
   double *ximhihihi, double *ximlohihi, double *ximhilohi, double *ximlolohi,
   double *ximhihilo, double *ximlohilo, double *ximhilolo, double *ximlololo,
   double *inv0rehihihi, double *inv0relohihi,
   double *inv0rehilohi, double *inv0relolohi,
   double *inv0rehihilo, double *inv0relohilo,
   double *inv0rehilolo, double *inv0relololo,
   double *inv0imhihihi, double *inv0imlohihi,
   double *inv0imhilohi, double *inv0imlolohi,
   double *inv0imhihilo, double *inv0imlohilo,
   double *inv0imhilolo, double *inv0imlololo,
   double *vrehihihi, double *vrelohihi, double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo, double *vrehilolo, double *vrelololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo,
   double *vimhilolo, double *vimlololo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int idx = bdx*szt + tdx;  // thread tdx scales idx

   __shared__ double shvrehihihi[inner_od_shmemsize];
   __shared__ double shvrelohihi[inner_od_shmemsize];
   __shared__ double shvrehilohi[inner_od_shmemsize];
   __shared__ double shvrelolohi[inner_od_shmemsize];
   __shared__ double shvrehihilo[inner_od_shmemsize];
   __shared__ double shvrelohilo[inner_od_shmemsize];
   __shared__ double shvrehilolo[inner_od_shmemsize];
   __shared__ double shvrelololo[inner_od_shmemsize];
   __shared__ double shvimhihihi[inner_od_shmemsize];
   __shared__ double shvimlohihi[inner_od_shmemsize];
   __shared__ double shvimhilohi[inner_od_shmemsize];
   __shared__ double shvimlolohi[inner_od_shmemsize];
   __shared__ double shvimhihilo[inner_od_shmemsize];
   __shared__ double shvimlohilo[inner_od_shmemsize];
   __shared__ double shvimhilolo[inner_od_shmemsize];
   __shared__ double shvimlololo[inner_od_shmemsize];
   __shared__ double acc[8];
   __shared__ double invre[8];
   __shared__ double invim[8];

   shvrehihihi[tdx] = xrehihihi[idx];
   shvrelohihi[tdx] = xrelohihi[idx];
   shvrehilohi[tdx] = xrehilohi[idx];
   shvrelolohi[tdx] = xrelolohi[idx];
   shvrehihilo[tdx] = xrehihilo[idx];
   shvrelohilo[tdx] = xrelohilo[idx];
   shvrehilolo[tdx] = xrehilolo[idx];
   shvrelololo[tdx] = xrelololo[idx];
   shvimhihihi[tdx] = ximhihihi[idx];
   shvimlohihi[tdx] = ximlohihi[idx];
   shvimhilohi[tdx] = ximhilohi[idx];
   shvimlolohi[tdx] = ximlolohi[idx];
   shvimhihilo[tdx] = ximhihilo[idx];
   shvimlohilo[tdx] = ximlohilo[idx];
   shvimhilolo[tdx] = ximhilolo[idx];
   shvimlololo[tdx] = ximlololo[idx];
   __syncthreads();

   double resulthihihi,resultlohihi,resulthilohi,resultlolohi;
   double resulthihilo,resultlohilo,resulthilolo,resultlololo;
/*
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;
   double invhihihi,invlohihi,invhilohi,invlolohi;
   double invhihilo,invlohilo,invhilolo,invlololo;
 */
   
   // shv[j] = shv[j]/v0;

   // resultre = vre[i]*inv0re - vim[i]*inv0im;
/*
   invhihihi = *inv0rehihihi;
   invlohihi = *inv0relohihi;
   invhilohi = *inv0rehilohi;
   invlolohi = *inv0relolohi;
   invhihilo = *inv0rehihilo;
   invlohilo = *inv0relohilo;
   invhilolo = *inv0rehilolo;
   invlololo = *inv0relololo;
 */
   invre[0] = *inv0rehihihi;
   invre[1] = *inv0relohihi;
   invre[2] = *inv0rehilohi;
   invre[3] = *inv0relolohi;
   invre[4] = *inv0rehihilo;
   invre[5] = *inv0relohilo;
   invre[6] = *inv0rehilolo;
   invre[7] = *inv0relololo;
   odg_mul(shvrehihihi[tdx],shvrelohihi[tdx],shvrehilohi[tdx],shvrelolohi[tdx],
           shvrehihilo[tdx],shvrelohilo[tdx],shvrehilolo[tdx],shvrelololo[tdx],
        // *inv0rehihihi,   *inv0relohihi,   *inv0rehilohi,   *inv0relolohi,
        // *inv0rehihilo,   *inv0relohilo,   *inv0rehilolo,   *inv0relololo,
        //     invhihihi,       invlohihi,       invhilohi,       invlolohi,
        //     invhihilo,       invlohilo,       invhilolo,       invlololo,
              invre[0],        invre[1],        invre[2],        invre[3],
              invre[4],        invre[5],        invre[6],        invre[7],
         &resulthihihi,   &resultlohihi,   &resulthilohi,   &resultlolohi,
         &resulthihilo,   &resultlohilo,   &resulthilolo,   &resultlololo);
/*
   invhihihi = *inv0imhihihi;
   invlohihi = *inv0imlohihi;
   invhilohi = *inv0imhilohi;
   invlolohi = *inv0imlolohi;
   invhihilo = *inv0imhihilo;
   invlohilo = *inv0imlohilo;
   invhilolo = *inv0imhilolo;
   invlololo = *inv0imlololo;
 */
   invim[0] = *inv0imhihihi;
   invim[1] = *inv0imlohihi;
   invim[2] = *inv0imhilohi;
   invim[3] = *inv0imlolohi;
   invim[4] = *inv0imhihilo;
   invim[5] = *inv0imlohilo;
   invim[6] = *inv0imhilolo;
   invim[7] = *inv0imlololo;
   odg_mul(shvimhihihi[tdx],shvimlohihi[tdx],shvimhilohi[tdx],shvimlolohi[tdx],
           shvimhihilo[tdx],shvimlohilo[tdx],shvimhilolo[tdx],shvimlololo[tdx],
      //  *inv0imhihihi,   *inv0imlohihi,   *inv0imhilohi,   *inv0imlolohi,
      //  *inv0imhihilo,   *inv0imlohilo,   *inv0imhilolo,   *inv0imlololo,
      //     invhihihi,       invlohihi,       invhilohi,       invlolohi,
      //     invhihilo,       invlohilo,       invhilolo,       invlololo,
            invim[0],        invim[1],        invim[2],        invim[3],
            invim[4],        invim[5],        invim[6],        invim[7],
      //    &acchihihi,      &acclohihi,      &acchilohi,      &acclolohi,
      //    &acchihilo,      &acclohilo,      &acchilolo,      &acclololo);
           &acc[0],&acc[1],&acc[2],&acc[3],&acc[4],&acc[5],&acc[6],&acc[7]);
   odg_dec(&resulthihihi,&resultlohihi,&resulthilohi,&resultlolohi,
           &resulthihilo,&resultlohilo,&resulthilolo,&resultlololo,
           acc[0],acc[1],acc[2],acc[3],acc[4],acc[5],acc[6],acc[7]);
             //  acchihihi,    acclohihi,    acchilohi,    acclolohi,
             //  acchihilo,    acclohilo,    acchilolo,    acclololo);

   __syncthreads();
   if(idx < dim)
   {
      vrehihihi[idx] = resulthihihi;
      vrelohihi[idx] = resultlohihi;    
      vrehilohi[idx] = resulthilohi;
      vrelolohi[idx] = resultlolohi;    
      vrehihilo[idx] = resulthihilo;
      vrelohilo[idx] = resultlohilo;    
      vrehilolo[idx] = resulthilolo;
      vrelololo[idx] = resultlololo;    
   }
   __syncthreads();
   // zim = vim[i]*inv0re + vre[i]*inv0im;
/*
   invhihihi = *inv0rehihihi;
   invlohihi = *inv0relohihi;
   invhilohi = *inv0rehilohi;
   invlolohi = *inv0relolohi;
   invhihilo = *inv0rehihilo;
   invlohilo = *inv0relohilo;
   invhilolo = *inv0rehilolo;
   invlololo = *inv0relololo;
 */
   odg_mul(shvimhihihi[tdx],shvimlohihi[tdx],shvimhilohi[tdx],shvimlolohi[tdx],
           shvimhihilo[tdx],shvimlohilo[tdx],shvimhilolo[tdx],shvimlololo[tdx],
        // *inv0rehihihi,   *inv0relohihi,   *inv0rehilohi,   *inv0relolohi,
        // *inv0rehihilo,   *inv0relohilo,   *inv0rehilolo,   *inv0relololo,
        //     invhihihi,       invlohihi,       invhilohi,       invlolohi,
        //     invhihilo,       invlohilo,       invhilolo,       invlololo,
              invre[0],        invre[1],        invre[2],        invre[3],
              invre[4],        invre[5],        invre[6],        invre[7],
         &resulthihihi,   &resultlohihi,   &resulthilohi,   &resultlolohi,
         &resulthihilo,   &resultlohilo,   &resulthilolo,   &resultlololo);
/*
   invhihihi = *inv0imhihihi;
   invlohihi = *inv0imlohihi;
   invhilohi = *inv0imhilohi;
   invlolohi = *inv0imlolohi;
   invhihilo = *inv0imhihilo;
   invlohilo = *inv0imlohilo;
   invhilolo = *inv0imhilolo;
   invlololo = *inv0imlololo;
 */
   odg_mul(shvrehihihi[tdx],shvrelohihi[tdx],shvrehilohi[tdx],shvrelolohi[tdx],
           shvrehihilo[tdx],shvrelohilo[tdx],shvrehilolo[tdx],shvrelololo[tdx],
      //  *inv0imhihihi,   *inv0imlohihi,   *inv0imhilohi,   *inv0imlolohi,
      //  *inv0imhihilo,   *inv0imlohilo,   *inv0imhilolo,   *inv0imlololo,
      //     invhihihi,       invlohihi,       invhilohi,       invlolohi,
      //     invhihilo,       invlohilo,       invhilolo,       invlololo,
            invim[0],        invim[1],        invim[2],        invim[3],
            invim[4],        invim[5],        invim[6],        invim[7],
           &acc[0],&acc[1],&acc[2],&acc[3],&acc[4],&acc[5],&acc[6],&acc[7]);
        //  &acchihihi,      &acclohihi,      &acchilohi,      &acclolohi,
        //  &acchihilo,      &acclohilo,      &acchilolo,      &acclololo);
   odg_inc(&resulthihihi,&resultlohihi,&resulthilohi,&resultlolohi,
           &resulthihilo,&resultlohilo,&resulthilolo,&resultlololo,
           acc[0],acc[1],acc[2],acc[3],acc[4],acc[5],acc[6],acc[7]);
          //   acchihihi,    acclohihi,    acchilohi,    acclolohi,
          //   acchihilo,    acclohilo,    acchilolo,    acclololo);

   __syncthreads();
   if(idx < dim)
   {
      vimhihihi[idx] = resulthihihi;
      vimlohihi[idx] = resultlohihi;
      vimhilohi[idx] = resulthilohi;
      vimlolohi[idx] = resultlolohi;
      vimhihilo[idx] = resulthihilo;
      vimlohilo[idx] = resultlohilo;
      vimhilolo[idx] = resulthilolo;
      vimlololo[idx] = resultlololo;
   }
}

__global__ void cmplx8_normalize_rere
 ( int dim, int szt,
   double *xrehihihi, double *xrelohihi, double *xrehilohi, double *xrelolohi,
   double *xrehihilo, double *xrelohilo, double *xrehilolo, double *xrelololo,
   double *inv0rehihihi, double *inv0relohihi,
   double *inv0rehilohi, double *inv0relolohi,
   double *inv0rehihilo, double *inv0relohilo,
   double *inv0rehilolo, double *inv0relololo,
   double *vrehihihi, double *vrelohihi,
   double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo,
   double *vrehilolo, double *vrelololo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int idx = bdx*szt + tdx;  // thread tdx scales idx

   __shared__ double shvrehihihi[inner_od_shmemsize];
   __shared__ double shvrelohihi[inner_od_shmemsize];
   __shared__ double shvrehilohi[inner_od_shmemsize];
   __shared__ double shvrelolohi[inner_od_shmemsize];
   __shared__ double shvrehihilo[inner_od_shmemsize];
   __shared__ double shvrelohilo[inner_od_shmemsize];
   __shared__ double shvrehilolo[inner_od_shmemsize];
   __shared__ double shvrelololo[inner_od_shmemsize];

   double resulthihihi,resultlohihi,resulthilohi,resultlolohi;
   double resulthihilo,resultlohilo,resulthilolo,resultlololo;

   shvrehihihi[tdx] = xrehihihi[idx];
   shvrelohihi[tdx] = xrelohihi[idx];
   shvrehilohi[tdx] = xrehilohi[idx];
   shvrelolohi[tdx] = xrelolohi[idx];
   shvrehihilo[tdx] = xrehihilo[idx];
   shvrelohilo[tdx] = xrelohilo[idx];
   shvrehilolo[tdx] = xrehilolo[idx];
   shvrelololo[tdx] = xrelololo[idx];

   __syncthreads();
   
   // shv[j] = shv[j]/v0;

   // resultre = vre[i]*inv0re - vim[i]*inv0im;
   // -> do only the first multiplication

   odg_mul(shvrehihihi[tdx],shvrelohihi[tdx],shvrehilohi[tdx],shvrelolohi[tdx],
           shvrehihilo[tdx],shvrelohilo[tdx],shvrehilolo[tdx],shvrelololo[tdx],
            inv0rehihihi[0], inv0relohihi[0], inv0rehilohi[0], inv0relolohi[0],
            inv0rehihilo[0], inv0relohilo[0], inv0rehilolo[0], inv0relololo[0],
         &resulthihihi,   &resultlohihi,   &resulthilohi,   &resultlolohi,
         &resulthihilo,   &resultlohilo,   &resulthilolo,   &resultlololo);

   __syncthreads();

   if(idx < dim)
   {
      vrehihihi[idx] = resulthihihi;
      vrelohihi[idx] = resultlohihi;    
      vrehilohi[idx] = resulthilohi;
      vrelolohi[idx] = resultlolohi;    
      vrehihilo[idx] = resulthihilo;
      vrelohilo[idx] = resultlohilo;    
      vrehilolo[idx] = resulthilolo;
      vrelololo[idx] = resultlololo;    
   }
}

__global__ void cmplx8_normalize_imim
 ( int dim, int szt,
   double *ximhihihi, double *ximlohihi, double *ximhilohi, double *ximlolohi,
   double *ximhihilo, double *ximlohilo, double *ximhilolo, double *ximlololo,
   double *inv0imhihihi, double *inv0imlohihi,
   double *inv0imhilohi, double *inv0imlolohi,
   double *inv0imhihilo, double *inv0imlohilo,
   double *inv0imhilolo, double *inv0imlololo,
   double *vrehihihi, double *vrelohihi,
   double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo,
   double *vrehilolo, double *vrelololo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int idx = bdx*szt + tdx;  // thread tdx scales idx

   __shared__ double shvimhihihi[inner_od_shmemsize];
   __shared__ double shvimlohihi[inner_od_shmemsize];
   __shared__ double shvimhilohi[inner_od_shmemsize];
   __shared__ double shvimlolohi[inner_od_shmemsize];
   __shared__ double shvimhihilo[inner_od_shmemsize];
   __shared__ double shvimlohilo[inner_od_shmemsize];
   __shared__ double shvimhilolo[inner_od_shmemsize];
   __shared__ double shvimlololo[inner_od_shmemsize];

   double resulthihihi,resultlohihi,resulthilohi,resultlolohi;
   double resulthihilo,resultlohilo,resulthilolo,resultlololo;
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   shvimhihihi[tdx] = ximhihihi[idx];
   shvimlohihi[tdx] = ximlohihi[idx];
   shvimhilohi[tdx] = ximhilohi[idx];
   shvimlolohi[tdx] = ximlolohi[idx];
   shvimhihilo[tdx] = ximhihilo[idx];
   shvimlohilo[tdx] = ximlohilo[idx];
   shvimhilolo[tdx] = ximhilolo[idx];
   shvimlololo[tdx] = ximlololo[idx];

   __syncthreads();
   
   // shv[j] = shv[j]/v0;

   // resultre = vre[i]*inv0re - vim[i]*inv0im;
   // -> update resultre with the subtraction of imaginary parts

   odg_mul(shvimhihihi[tdx],shvimlohihi[tdx],shvimhilohi[tdx],shvimlolohi[tdx],
           shvimhihilo[tdx],shvimlohilo[tdx],shvimhilolo[tdx],shvimlololo[tdx],
            inv0imhihihi[0], inv0imlohihi[0], inv0imhilohi[0], inv0imlolohi[0],
            inv0imhihilo[0], inv0imlohilo[0], inv0imhilolo[0], inv0imlololo[0],
           &acchihihi,      &acclohihi,      &acchilohi,      &acclolohi,
           &acchihilo,      &acclohilo,      &acchilolo,      &acclololo);

   resulthihihi = vrehihihi[idx];
   resultlohihi = vrelohihi[idx];
   resulthilohi = vrehilohi[idx];
   resultlolohi = vrelolohi[idx];
   resulthihilo = vrehihilo[idx];
   resultlohilo = vrelohilo[idx];
   resulthilolo = vrehilolo[idx];
   resultlololo = vrelololo[idx];

   odg_dec(&resulthihihi,&resultlohihi,&resulthilohi,&resultlolohi,
           &resulthihilo,&resultlohilo,&resulthilolo,&resultlololo,
               acchihihi,    acclohihi,    acchilohi,    acclolohi,
               acchihilo,    acclohilo,    acchilolo,    acclololo);

   __syncthreads();

   if(idx < dim)
   {
      vrehihihi[idx] = resulthihihi;
      vrelohihi[idx] = resultlohihi;    
      vrehilohi[idx] = resulthilohi;
      vrelolohi[idx] = resultlolohi;    
      vrehihilo[idx] = resulthihilo;
      vrelohilo[idx] = resultlohilo;    
      vrehilolo[idx] = resulthilolo;
      vrelololo[idx] = resultlololo;    
   }
}

__global__ void cmplx8_normalize_imre
 ( int dim, int szt,
   double *ximhihihi, double *ximlohihi, double *ximhilohi, double *ximlolohi,
   double *ximhihilo, double *ximlohilo, double *ximhilolo, double *ximlololo,
   double *inv0rehihihi, double *inv0relohihi,
   double *inv0rehilohi, double *inv0relolohi,
   double *inv0rehihilo, double *inv0relohilo,
   double *inv0rehilolo, double *inv0relololo,
   double *vimhihihi, double *vimlohihi,
   double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo,
   double *vimhilolo, double *vimlololo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int idx = bdx*szt + tdx;  // thread tdx scales idx

   __shared__ double shvimhihihi[inner_od_shmemsize];
   __shared__ double shvimlohihi[inner_od_shmemsize];
   __shared__ double shvimhilohi[inner_od_shmemsize];
   __shared__ double shvimlolohi[inner_od_shmemsize];
   __shared__ double shvimhihilo[inner_od_shmemsize];
   __shared__ double shvimlohilo[inner_od_shmemsize];
   __shared__ double shvimhilolo[inner_od_shmemsize];
   __shared__ double shvimlololo[inner_od_shmemsize];

   double resulthihihi,resultlohihi,resulthilohi,resultlolohi;
   double resulthihilo,resultlohilo,resulthilolo,resultlololo;

   shvimhihihi[tdx] = ximhihihi[idx];
   shvimlohihi[tdx] = ximlohihi[idx];
   shvimhilohi[tdx] = ximhilohi[idx];
   shvimlolohi[tdx] = ximlolohi[idx];
   shvimhihilo[tdx] = ximhihilo[idx];
   shvimlohilo[tdx] = ximlohilo[idx];
   shvimhilolo[tdx] = ximhilolo[idx];
   shvimlololo[tdx] = ximlololo[idx];

   __syncthreads();
   
   // shv[j] = shv[j]/v0;

   // resultim = vim[i]*inv0re + vre[i]*inv0im;
   // -> do only the first multiplication

   odg_mul(shvimhihihi[tdx],shvimlohihi[tdx],shvimhilohi[tdx],shvimlolohi[tdx],
           shvimhihilo[tdx],shvimlohilo[tdx],shvimhilolo[tdx],shvimlololo[tdx],
            inv0rehihihi[0], inv0relohihi[0], inv0rehilohi[0], inv0relolohi[0],
            inv0rehihilo[0], inv0relohilo[0], inv0rehilolo[0], inv0relololo[0],
           &resulthihihi,   &resultlohihi,   &resulthilohi,   &resultlolohi,
           &resulthihilo,   &resultlohilo,   &resulthilolo,   &resultlololo);

   __syncthreads();

   if(idx < dim)
   {
      vimhihihi[idx] = resulthihihi;
      vimlohihi[idx] = resultlohihi;
      vimhilohi[idx] = resulthilohi;
      vimlolohi[idx] = resultlolohi;
      vimhihilo[idx] = resulthihilo;
      vimlohilo[idx] = resultlohilo;
      vimhilolo[idx] = resulthilolo;
      vimlololo[idx] = resultlololo;
   }
}

__global__ void cmplx8_normalize_reim
 ( int dim, int szt,
   double *xrehihihi, double *xrelohihi, double *xrehilohi, double *xrelolohi,
   double *xrehihilo, double *xrelohilo, double *xrehilolo, double *xrelololo,
   double *inv0imhihihi, double *inv0imlohihi,
   double *inv0imhilohi, double *inv0imlolohi,
   double *inv0imhihilo, double *inv0imlohilo,
   double *inv0imhilolo, double *inv0imlololo,
   double *vimhihihi, double *vimlohihi,
   double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo,
   double *vimhilolo, double *vimlololo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int idx = bdx*szt + tdx;  // thread tdx scales idx

   __shared__ double shvrehihihi[inner_od_shmemsize];
   __shared__ double shvrelohihi[inner_od_shmemsize];
   __shared__ double shvrehilohi[inner_od_shmemsize];
   __shared__ double shvrelolohi[inner_od_shmemsize];
   __shared__ double shvrehihilo[inner_od_shmemsize];
   __shared__ double shvrelohilo[inner_od_shmemsize];
   __shared__ double shvrehilolo[inner_od_shmemsize];
   __shared__ double shvrelololo[inner_od_shmemsize];

   double resulthihihi,resultlohihi,resulthilohi,resultlolohi;
   double resulthihilo,resultlohilo,resulthilolo,resultlololo;
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   shvrehihihi[tdx] = xrehihihi[idx];
   shvrelohihi[tdx] = xrelohihi[idx];
   shvrehilohi[tdx] = xrehilohi[idx];
   shvrelolohi[tdx] = xrelolohi[idx];
   shvrehihilo[tdx] = xrehihilo[idx];
   shvrelohilo[tdx] = xrelohilo[idx];
   shvrehilolo[tdx] = xrehilolo[idx];
   shvrelololo[tdx] = xrelololo[idx];
   __syncthreads();
   
   // shv[j] = shv[j]/v0;

   // resultim = vim[i]*inv0re + vre[i]*inv0im;
   // -> update resultim with the second term

   odg_mul(shvrehihihi[tdx],shvrelohihi[tdx],shvrehilohi[tdx],shvrelolohi[tdx],
           shvrehihilo[tdx],shvrelohilo[tdx],shvrehilolo[tdx],shvrelololo[tdx],
            inv0imhihihi[0], inv0imlohihi[0], inv0imhilohi[0], inv0imlolohi[0],
            inv0imhihilo[0], inv0imlohilo[0], inv0imhilolo[0], inv0imlololo[0],
           &acchihihi,      &acclohihi,      &acchilohi,      &acclolohi,
           &acchihilo,      &acclohilo,      &acchilolo,      &acclololo);

   resulthihihi = vimhihihi[idx];
   resultlohihi = vimlohihi[idx];
   resulthilohi = vimhilohi[idx];
   resultlolohi = vimlolohi[idx];
   resulthihilo = vimhihilo[idx];
   resultlohilo = vimlohilo[idx];
   resulthilolo = vimhilolo[idx];
   resultlololo = vimlololo[idx];

   odg_inc(&resulthihihi,&resultlohihi,&resulthilohi,&resultlolohi,
           &resulthihilo,&resultlohilo,&resulthilolo,&resultlololo,
               acchihihi,    acclohihi,    acchilohi,    acclolohi,
               acchihilo,    acclohilo,    acchilolo,    acclololo);

   __syncthreads();

   if(idx < dim)
   {
      vimhihihi[idx] = resulthihihi;
      vimlohihi[idx] = resultlohihi;
      vimhilohi[idx] = resulthilohi;
      vimlolohi[idx] = resultlolohi;
      vimhihilo[idx] = resulthihilo;
      vimlohilo[idx] = resultlohilo;
      vimhilolo[idx] = resulthilolo;
      vimlololo[idx] = resultlololo;
   }
}

__global__ void dbl8_small_leftRupdate
 ( int nrows, int ncols, int szt, int k,
   double *Rhihihi, double *Rlohihi, double *Rhilohi, double *Rlolohi,
   double *Rhihilo, double *Rlohilo, double *Rhilolo, double *Rlololo,
   double *vhihihi, double *vlohihi, double *vhilohi, double *vlolohi,
   double *vhihilo, double *vlohilo, double *vhilolo, double *vlololo,
   double *betahihihi, double *betalohihi,
   double *betahilohi, double *betalolohi,
   double *betahihilo, double *betalohilo,
   double *betahilolo, double *betalololo )
{
   const int tdx = threadIdx.x;          // index of thread in block
   const int Roffset = k*nrows + k;
   int Rcolidx;
   double whihihi,wlohihi,whilohi,wlolohi;
   double whihilo,wlohilo,whilolo,wlololo;
   double Rtdxhihihi,Rtdxlohihi,Rtdxhilohi,Rtdxlolohi;
   double Rtdxhihilo,Rtdxlohilo,Rtdxhilolo,Rtdxlololo;
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   __shared__ double shvhihihi[od_shmemsize]; // slice of v
   __shared__ double shvlohihi[od_shmemsize]; 
   __shared__ double shvhilohi[od_shmemsize]; 
   __shared__ double shvlolohi[od_shmemsize]; 
   __shared__ double shvhihilo[od_shmemsize];
   __shared__ double shvlohilo[od_shmemsize]; 
   __shared__ double shvhilolo[od_shmemsize]; 
   __shared__ double shvlololo[od_shmemsize]; 

   shvhihihi[tdx] = vhihihi[tdx];
   shvlohihi[tdx] = vlohihi[tdx];
   shvhilohi[tdx] = vhilohi[tdx];
   shvlolohi[tdx] = vlolohi[tdx];
   shvhihilo[tdx] = vhihilo[tdx];
   shvlohilo[tdx] = vlohilo[tdx];
   shvhilolo[tdx] = vhilolo[tdx];
   shvlololo[tdx] = vlololo[tdx];
   __syncthreads();
   whihihi = 0.0;
   wlohihi = 0.0;
   whilohi = 0.0;
   wlolohi = 0.0;
   whihilo = 0.0;
   wlohilo = 0.0;
   whilolo = 0.0;
   wlololo = 0.0;

   for(int i=0; i<nrows-k; i++)   // loop through rows of R
   {
      Rcolidx = Roffset + i + tdx*nrows;
      __syncthreads();
      Rtdxhihihi = Rhihihi[Rcolidx];
      Rtdxlohihi = Rlohihi[Rcolidx];
      Rtdxhilohi = Rhilohi[Rcolidx];
      Rtdxlolohi = Rlolohi[Rcolidx];
      Rtdxhihilo = Rhihilo[Rcolidx];
      Rtdxlohilo = Rlohilo[Rcolidx];
      Rtdxhilolo = Rhilolo[Rcolidx];
      Rtdxlololo = Rlololo[Rcolidx];
      // w = w + Rtdx*shv[i];
      __syncthreads();
      odg_mul(Rtdxhihihi,  Rtdxlohihi,  Rtdxhilohi,  Rtdxlolohi,
              Rtdxhihilo,  Rtdxlohilo,  Rtdxhilolo,  Rtdxlololo,
               shvhihihi[i],shvlohihi[i],shvhilohi[i],shvlolohi[i],
               shvhihilo[i],shvlohilo[i],shvhilolo[i],shvlololo[i],
              &acchihihi,  &acclohihi,  &acchilohi,  &acclolohi,
              &acchihilo,  &acclohilo,  &acchilolo,  &acclololo);
      odg_inc( &whihihi, &wlohihi, &whilohi, &wlolohi,
               &whihilo, &wlohilo, &whilolo, &wlololo,
              acchihihi,acclohihi,acchilohi,acclolohi,
              acchihilo,acclohilo,acchilolo,acclololo);
   }
   // w = (*beta)*w;
   // odg_mlt(&whi,&wlo,*betahi,*betalo); <-- this does not work!
   __syncthreads();
   odg_mul(*betahihihi,*betalohihi,*betahilohi,*betalolohi,
           *betahihilo,*betalohilo,*betahilolo,*betalololo,
               whihihi,    wlohihi,    whilohi,    wlolohi,
               whihilo,    wlohilo,    whilolo,    wlololo,
            &acchihihi, &acclohihi, &acchilohi, &acclolohi,
            &acchihilo, &acclohilo, &acchilolo, &acclololo);
   whihihi = acchihihi;
   wlohihi = acclohihi;
   whilohi = acchilohi;
   wlolohi = acclolohi;
   whihilo = acchihilo;
   wlohilo = acclohilo;
   whilolo = acchilolo;
   wlololo = acclololo;
   __syncthreads();
   for(int i=0; i<nrows-k; i++)   // update i-th row of R
   {
      Rcolidx = Roffset + i + tdx*nrows;
      __syncthreads();
      Rtdxhihihi = Rhihihi[Rcolidx];
      Rtdxlohihi = Rlohihi[Rcolidx];
      Rtdxhilohi = Rhilohi[Rcolidx];
      Rtdxlolohi = Rlolohi[Rcolidx];
      Rtdxhihilo = Rhihilo[Rcolidx];
      Rtdxlohilo = Rlohilo[Rcolidx];
      Rtdxhilolo = Rhilolo[Rcolidx];
      Rtdxlololo = Rlololo[Rcolidx];
      // Rtdx = Rtdx - shv[i]*w;
      __syncthreads();
      odg_mul(shvhihihi[i],shvlohihi[i],shvhilohi[i],shvlolohi[i],
              shvhihilo[i],shvlohilo[i],shvhilolo[i],shvlololo[i],
                whihihi,     wlohihi,     whilohi,     wlolohi,
                whihilo,     wlohilo,     whilolo,     wlololo,
             &acchihihi,  &acclohihi,  &acchilohi,  &acclolohi,
             &acchihilo,  &acclohilo,  &acchilolo,  &acclololo);
      odg_dec(&Rtdxhihihi,&Rtdxlohihi,&Rtdxhilohi,&Rtdxlolohi,
              &Rtdxhihilo,&Rtdxlohilo,&Rtdxhilolo,&Rtdxlololo,
                acchihihi,  acclohihi,  acchilohi,  acclolohi,
                acchihilo,  acclohilo,  acchilolo,  acclololo);
      __syncthreads();
      // changed nrows-k into ncols-k, where ncols = szt
      if(tdx < ncols-k)
      {
         Rhihihi[Rcolidx] = Rtdxhihihi;
         Rlohihi[Rcolidx] = Rtdxlohihi;
         Rhilohi[Rcolidx] = Rtdxhilohi;
         Rlolohi[Rcolidx] = Rtdxlolohi;
         Rhihilo[Rcolidx] = Rtdxhihilo;
         Rlohilo[Rcolidx] = Rtdxlohilo;
         Rhilolo[Rcolidx] = Rtdxhilolo;
         Rlololo[Rcolidx] = Rtdxlololo;
      }
      __syncthreads();
   }
}

__global__ void cmplx8_small_leftRupdate
 ( int nrows, int ncols, int szt, int k,
   double *Rrehihihi, double *Rrelohihi, double *Rrehilohi, double *Rrelolohi,
   double *Rrehihilo, double *Rrelohilo, double *Rrehilolo, double *Rrelololo,
   double *Rimhihihi, double *Rimlohihi, double *Rimhilohi, double *Rimlolohi,
   double *Rimhihilo, double *Rimlohilo, double *Rimhilolo, double *Rimlololo,
   double *vrehihihi, double *vrelohihi, double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo, double *vrehilolo, double *vrelololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo, double *vimhilolo, double *vimlololo,
   double *betahihihi, double *betalohihi, 
   double *betahilohi, double *betalolohi,
   double *betahihilo, double *betalohilo,
   double *betahilolo, double *betalololo )
{
   const int tdx = threadIdx.x;          // index of thread in block
   const int Roffset = k*nrows + k;
   const double bthihihi = *betahihihi;
   const double btlohihi = *betalohihi;
   const double bthilohi = *betahilohi;
   const double btlolohi = *betalolohi;
   const double bthihilo = *betahihilo;
   const double btlohilo = *betalohilo;
   const double bthilolo = *betahilolo;
   const double btlololo = *betalololo;

   int Rcolidx;
   double w_rehihihi,w_relohihi,w_rehilohi,w_relolohi;
   double w_rehihilo,w_relohilo,w_rehilolo,w_relololo;
   double w_imhihihi,w_imlohihi,w_imhilohi,w_imlolohi;
   double w_imhihilo,w_imlohilo,w_imhilolo,w_imlololo;
   double Rtdx_hihihi,Rtdx_lohihi,Rtdx_hilohi,Rtdx_lolohi;
   double Rtdx_hihilo,Rtdx_lohilo,Rtdx_hilolo,Rtdx_lololo;
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   __shared__ double shvrehihihi[cod_shmemsize]; // slice of v
   __shared__ double shvrelohihi[cod_shmemsize];
   __shared__ double shvrehilohi[cod_shmemsize];
   __shared__ double shvrelolohi[cod_shmemsize];
   __shared__ double shvrehihilo[cod_shmemsize];
   __shared__ double shvrelohilo[cod_shmemsize];
   __shared__ double shvrehilolo[cod_shmemsize];
   __shared__ double shvrelololo[cod_shmemsize];
   __shared__ double shvimhihihi[cod_shmemsize];
   __shared__ double shvimlohihi[cod_shmemsize];
   __shared__ double shvimhilohi[cod_shmemsize];
   __shared__ double shvimlolohi[cod_shmemsize];
   __shared__ double shvimhihilo[cod_shmemsize];
   __shared__ double shvimlohilo[cod_shmemsize];
   __shared__ double shvimhilolo[cod_shmemsize];
   __shared__ double shvimlololo[cod_shmemsize];

   shvrehihihi[tdx] = vrehihihi[tdx];
   shvrelohihi[tdx] = vrelohihi[tdx];
   shvrehilohi[tdx] = vrehilohi[tdx];
   shvrelolohi[tdx] = vrelolohi[tdx];
   shvrehihilo[tdx] = vrehihilo[tdx];
   shvrelohilo[tdx] = vrelohilo[tdx];
   shvrehilolo[tdx] = vrehilolo[tdx];
   shvrelololo[tdx] = vrelololo[tdx];
   shvimhihihi[tdx] = vimhihihi[tdx];
   shvimlohihi[tdx] = vimlohihi[tdx];
   shvimhilohi[tdx] = vimhilohi[tdx];
   shvimlolohi[tdx] = vimlolohi[tdx];
   shvimhihilo[tdx] = vimhihilo[tdx];
   shvimlohilo[tdx] = vimlohilo[tdx];
   shvimhilolo[tdx] = vimhilolo[tdx];
   shvimlololo[tdx] = vimlololo[tdx];
   __syncthreads();
   w_rehihihi = 0.0;
   w_relohihi = 0.0;
   w_rehilohi = 0.0;
   w_relolohi = 0.0;
   w_rehihilo = 0.0;
   w_relohilo = 0.0;
   w_rehilolo = 0.0;
   w_relololo = 0.0;
   w_imhihihi = 0.0;
   w_imlohihi = 0.0;
   w_imhilohi = 0.0;
   w_imlolohi = 0.0;
   w_imhihilo = 0.0;
   w_imlohilo = 0.0;
   w_imhilolo = 0.0;
   w_imlololo = 0.0;

   for(int i=0; i<nrows-k; i++)   // loop through rows of R
   {
      Rcolidx = Roffset + i + tdx*nrows;
      __syncthreads();
      Rtdx_hihihi = Rrehihihi[Rcolidx];
      Rtdx_lohihi = Rrelohihi[Rcolidx];
      Rtdx_hilohi = Rrehilohi[Rcolidx];
      Rtdx_lolohi = Rrelolohi[Rcolidx];
      Rtdx_hihilo = Rrehihilo[Rcolidx];
      Rtdx_lohilo = Rrelohilo[Rcolidx];
      Rtdx_hilolo = Rrehilolo[Rcolidx];
      Rtdx_lololo = Rrelololo[Rcolidx];
      // w = w + Rtdx*shv[i]; beware of the Hermitian transpose!
      // w_re = w_re + Rtdx_re*shvre[i] + Rtdx_im*shvim[i];
      __syncthreads();
      odg_mul(Rtdx_hihihi,   Rtdx_lohihi,   Rtdx_hilohi,   Rtdx_lolohi,
              Rtdx_hihilo,   Rtdx_lohilo,   Rtdx_hilolo,   Rtdx_lololo,
              shvrehihihi[i],shvrelohihi[i],shvrehilohi[i],shvrelolohi[i],
              shvrehihilo[i],shvrelohilo[i],shvrehilolo[i],shvrelololo[i],
               &acchihihi,    &acclohihi,    &acchilohi,    &acclolohi,
               &acchihilo,    &acclohilo,    &acchilolo,    &acclololo);
      odg_inc(&w_rehihihi,&w_relohihi,&w_rehilohi,&w_relolohi,
              &w_rehihilo,&w_relohilo,&w_rehilolo,&w_relololo,
                acchihihi,  acclohihi,  acchilohi,  acclolohi,
                acchihilo,  acclohilo,  acchilolo,  acclololo);

      Rtdx_hihihi = Rimhihihi[Rcolidx];
      Rtdx_lohihi = Rimlohihi[Rcolidx];
      Rtdx_hilohi = Rimhilohi[Rcolidx];
      Rtdx_lolohi = Rimlolohi[Rcolidx];
      Rtdx_hihilo = Rimhihilo[Rcolidx];
      Rtdx_lohilo = Rimlohilo[Rcolidx];
      Rtdx_hilolo = Rimhilolo[Rcolidx];
      Rtdx_lololo = Rimlololo[Rcolidx];
      __syncthreads();
      odg_mul(Rtdx_hihihi,   Rtdx_lohihi,   Rtdx_hilohi,   Rtdx_lolohi,
              Rtdx_hihilo,   Rtdx_lohilo,   Rtdx_hilolo,   Rtdx_lololo,
              shvimhihihi[i],shvimlohihi[i],shvimhilohi[i],shvimlolohi[i],
              shvimhihilo[i],shvimlohilo[i],shvimhilolo[i],shvimlololo[i],
               &acchihihi,    &acclohihi,    &acchilohi,    &acclolohi,
               &acchihilo,    &acclohilo,    &acchilolo,    &acclololo);
      odg_inc(&w_rehihihi,&w_relohihi,&w_rehilohi,&w_relolohi,
              &w_rehihilo,&w_relohilo,&w_rehilolo,&w_relololo,
                acchihihi,  acclohihi,  acchilohi,  acclolohi,
                acchihilo,  acclohilo,  acchilolo,  acclololo);
      // w_im = w_im - Rtdx_im*shvre[i] + Rtdx_re*shvim[i];
      odg_mul(Rtdx_hihihi,   Rtdx_lohihi,   Rtdx_hilohi,   Rtdx_lolohi,
              Rtdx_hihilo,   Rtdx_lohilo,   Rtdx_hilolo,   Rtdx_lololo,
              shvrehihihi[i],shvrelohihi[i],shvrehilohi[i],shvrelolohi[i],
              shvrehihilo[i],shvrelohilo[i],shvrehilolo[i],shvrelololo[i],
               &acchihihi,    &acclohihi,    &acchilohi,    &acclolohi,
               &acchihilo,    &acclohilo,    &acchilolo,    &acclololo);
      odg_dec(&w_imhihihi,&w_imlohihi,&w_imhilohi,&w_imlolohi,
              &w_imhihilo,&w_imlohilo,&w_imhilolo,&w_imlololo,
                acchihihi,  acclohihi,  acchilohi,  acclolohi,
                acchihilo,  acclohilo,  acchilolo,  acclololo);

      Rtdx_hihihi = Rrehihihi[Rcolidx];
      Rtdx_lohihi = Rrelohihi[Rcolidx];
      Rtdx_hilohi = Rrehilohi[Rcolidx];
      Rtdx_lolohi = Rrelolohi[Rcolidx];
      Rtdx_hihilo = Rrehihilo[Rcolidx];
      Rtdx_lohilo = Rrelohilo[Rcolidx];
      Rtdx_hilolo = Rrehilolo[Rcolidx];
      Rtdx_lololo = Rrelololo[Rcolidx];
      odg_mul(Rtdx_hihihi,   Rtdx_lohihi,   Rtdx_hilohi,   Rtdx_lolohi,
              Rtdx_hihilo,   Rtdx_lohilo,   Rtdx_hilolo,   Rtdx_lololo,
              shvimhihihi[i],shvimlohihi[i],shvimhilohi[i],shvimlolohi[i],
              shvimhihilo[i],shvimlohilo[i],shvimhilolo[i],shvimlololo[i],
               &acchihihi,    &acclohihi,    &acchilohi,    &acclolohi,
               &acchihilo,    &acclohilo,    &acchilolo,    &acclololo);
      odg_inc(&w_imhihihi,&w_imlohihi,&w_imhilohi,&w_imlolohi,
              &w_imhihilo,&w_imlohilo,&w_imhilolo,&w_imlololo,
                acchihihi,  acclohihi,  acchilohi,  acclolohi,
                acchihilo,  acclohilo,  acchilolo,  acclololo);
   }
   // w_re = acc*w_re;
   __syncthreads();
   odg_mul(w_rehihihi,w_relohihi,w_rehilohi,w_relolohi,
           w_rehihilo,w_relohilo,w_rehilolo,w_relololo,
             bthihihi,  btlohihi,  bthilohi,  btlolohi,
             bthihilo,  btlohilo,  bthilolo,  btlololo,
           &acchihihi,&acclohihi,&acchilohi,&acclolohi,
           &acchihilo,&acclohilo,&acchilolo,&acclololo);
   w_rehihihi = acchihihi;
   w_relohihi = acclohihi;
   w_rehilohi = acchilohi;
   w_relolohi = acclolohi;
   w_rehihilo = acchihilo;
   w_relohilo = acclohilo;
   w_rehilolo = acchilolo;
   w_relololo = acclololo;
   // w_im = acc*w_im;
   __syncthreads();
   odg_mul(w_imhihihi,w_imlohihi,w_imhilohi,w_imlolohi,
           w_imhihilo,w_imlohilo,w_imhilolo,w_imlololo,
             bthihihi,  btlohihi,  bthilohi,  btlolohi,
             bthihilo,  btlohilo,  bthilolo,  btlololo,
           &acchihihi,&acclohihi,&acchilohi,&acclolohi,
           &acchihilo,&acclohilo,&acchilolo,&acclololo);
   w_imhihihi = acchihihi;
   w_imlohihi = acclohihi;
   w_imhilohi = acchilohi;
   w_imlolohi = acclolohi;
   w_imhihilo = acchihilo;
   w_imlohilo = acclohilo;
   w_imhilolo = acchilolo;
   w_imlololo = acclololo;
   __syncthreads();
   for(int i=0; i<nrows-k; i++)   // update i-th row of Rre
   {
      Rcolidx = Roffset + i + tdx*nrows;
      __syncthreads();
      Rtdx_hihihi = Rrehihihi[Rcolidx];
      Rtdx_lohihi = Rrelohihi[Rcolidx];
      Rtdx_hilohi = Rrehilohi[Rcolidx];
      Rtdx_lolohi = Rrelolohi[Rcolidx];
      Rtdx_hihilo = Rrehihilo[Rcolidx];
      Rtdx_lohilo = Rrelohilo[Rcolidx];
      Rtdx_hilolo = Rrehilolo[Rcolidx];
      Rtdx_lololo = Rrelololo[Rcolidx];
      // Rtdx = Rtdx - shv[i]*w; beware of the Hermitian transpose!
      // Rtdx_re = Rtdx_re - (shvre[i]*w_re + shvim[i]*w_im);
      __syncthreads();
      odg_mul(shvrehihihi[i],shvrelohihi[i],shvrehilohi[i],shvrelolohi[i],
              shvrehihilo[i],shvrelohilo[i],shvrehilolo[i],shvrelololo[i],
               w_rehihihi,    w_relohihi,    w_rehilohi,    w_relolohi,
               w_rehihilo,    w_relohilo,    w_rehilolo,    w_relololo,
               &acchihihi,    &acclohihi,    &acchilohi,    &acclolohi,
               &acchihilo,    &acclohilo,    &acchilolo,    &acclololo);
      odg_dec(&Rtdx_hihihi,&Rtdx_lohihi,&Rtdx_hilohi,&Rtdx_lolohi,
              &Rtdx_hihilo,&Rtdx_lohilo,&Rtdx_hilolo,&Rtdx_lololo,
                 acchihihi,   acclohihi,   acchilohi,   acclolohi,
                 acchihilo,   acclohilo,   acchilolo,   acclololo);
      odg_mul(shvimhihihi[i],shvimlohihi[i],shvimhilohi[i],shvimlolohi[i],
              shvimhihilo[i],shvimlohilo[i],shvimhilolo[i],shvimlololo[i],
               w_imhihihi,    w_imlohihi,    w_imhilohi,    w_imlolohi,
               w_imhihilo,    w_imlohilo,    w_imhilolo,    w_imlololo,
               &acchihihi,    &acclohihi,    &acchilohi,    &acclolohi,
               &acchihilo,    &acclohilo,    &acchilolo,    &acclololo);
      odg_dec(&Rtdx_hihihi,&Rtdx_lohihi,&Rtdx_hilohi,&Rtdx_lolohi,
              &Rtdx_hihilo,&Rtdx_lohilo,&Rtdx_hilolo,&Rtdx_lololo,
                 acchihihi,   acclohihi,   acchilohi,   acclolohi,
                 acchihilo,   acclohilo,   acchilolo,   acclololo);
      __syncthreads();
      // changed nrows-k into ncols-k, where ncols = szt
      if(tdx < ncols-k)
      {
         Rrehihihi[Rcolidx] = Rtdx_hihihi;
         Rrelohihi[Rcolidx] = Rtdx_lohihi;
         Rrehilohi[Rcolidx] = Rtdx_hilohi;
         Rrelolohi[Rcolidx] = Rtdx_lolohi;
         Rrehihilo[Rcolidx] = Rtdx_hihilo;
         Rrelohilo[Rcolidx] = Rtdx_lohilo;
         Rrehilolo[Rcolidx] = Rtdx_hilolo;
         Rrelololo[Rcolidx] = Rtdx_lololo;
      }
   }
   __syncthreads();
   for(int i=0; i<nrows-k; i++)   // update i-th row of Rim
   {
      Rcolidx = Roffset + i + tdx*nrows;
      __syncthreads();
      Rtdx_hihihi = Rimhihihi[Rcolidx];
      Rtdx_lohihi = Rimlohihi[Rcolidx];
      Rtdx_hilohi = Rimhilohi[Rcolidx];
      Rtdx_lolohi = Rimlolohi[Rcolidx];
      Rtdx_hihilo = Rimhihilo[Rcolidx];
      Rtdx_lohilo = Rimlohilo[Rcolidx];
      Rtdx_hilolo = Rimhilolo[Rcolidx];
      Rtdx_lololo = Rimlololo[Rcolidx];
      // Rtdx_im = Rtdx_im - (shvim[i]*w_re - shvre[i]*w_im);
      odg_mul(shvimhihihi[i],shvimlohihi[i],shvimhilohi[i],shvimlolohi[i],
              shvimhihilo[i],shvimlohilo[i],shvimhilolo[i],shvimlololo[i],
               w_rehihihi,    w_relohihi,    w_rehilohi,    w_relolohi,
               w_rehihilo,    w_relohilo,    w_rehilolo,    w_relololo,
               &acchihihi,    &acclohihi,    &acchilohi,    &acclolohi,
               &acchihilo,    &acclohilo,    &acchilolo,    &acclololo);
      odg_dec(&Rtdx_hihihi,&Rtdx_lohihi,&Rtdx_hilohi,&Rtdx_lolohi,
              &Rtdx_hihilo,&Rtdx_lohilo,&Rtdx_hilolo,&Rtdx_lololo,
                 acchihihi,   acclohihi,   acchilohi,   acclolohi,
                 acchihilo,   acclohilo,   acchilolo,   acclololo);
      odg_mul(shvrehihihi[i],shvrelohihi[i],shvrehilohi[i],shvrelolohi[i],
              shvrehihilo[i],shvrelohilo[i],shvrehilolo[i],shvrelololo[i],
               w_imhihihi,    w_imlohihi,    w_imhilohi,    w_imlolohi,
               w_imhihilo,    w_imlohilo,    w_imhilolo,    w_imlololo,
               &acchihihi,    &acclohihi,    &acchilohi,    &acclolohi,
               &acchihilo,    &acclohilo,    &acchilolo,    &acclololo);
      odg_inc(&Rtdx_hihihi,&Rtdx_lohihi,&Rtdx_hilohi,&Rtdx_lolohi,
              &Rtdx_hihilo,&Rtdx_lohilo,&Rtdx_hilolo,&Rtdx_lololo,
                 acchihihi,   acclohihi,   acchilohi,   acclolohi,
                 acchihilo,   acclohilo,   acchilolo,   acclololo);
      __syncthreads();
      // changed nrows-k into ncols-k, where ncols = szt
      if(tdx < ncols-k)
      {
         Rimhihihi[Rcolidx] = Rtdx_hihihi;
         Rimlohihi[Rcolidx] = Rtdx_lohihi;
         Rimhilohi[Rcolidx] = Rtdx_hilohi;
         Rimlolohi[Rcolidx] = Rtdx_lolohi;
         Rimhihilo[Rcolidx] = Rtdx_hihilo;
         Rimlohilo[Rcolidx] = Rtdx_lohilo;
         Rimhilolo[Rcolidx] = Rtdx_hilolo;
         Rimlololo[Rcolidx] = Rtdx_lololo;
      }
      __syncthreads();
   }
}

__global__ void dbl8_RTdotv
 ( int nrows, int szt, int colidx, int Roffset, int dim,
   double *Rhihihi, double *Rlohihi, double *Rhilohi, double *Rlolohi,
   double *Rhihilo, double *Rlohilo, double *Rhilolo, double *Rlololo,
   double *vhihihi, double *vlohihi, double *vhilohi, double *vlolohi,
   double *vhihilo, double *vlohilo, double *vhilolo, double *vlololo,
   double *RTdotvhihihi, double *RTdotvlohihi,
   double *RTdotvhilohi, double *RTdotvlolohi,
   double *RTdotvhihilo, double *RTdotvlohilo,
   double *RTdotvhilolo, double *RTdotvlololo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int idx = bdx*szt + tdx;        // thread tdx computes RTv[idx]

   const int vdx = idx % nrows;          // index in v is column in R^T
   const int row = idx / nrows;          // R is stored column-by-column

   const int Rdx = Roffset + idx + (row+1)*colidx;

   const double Vvalhihihi = vhihihi[vdx];
   const double Vvallohihi = vlohihi[vdx];
   const double Vvalhilohi = vhilohi[vdx];
   const double Vvallolohi = vlolohi[vdx];
   const double Vvalhihilo = vhihilo[vdx];
   const double Vvallohilo = vlohilo[vdx];
   const double Vvalhilolo = vhilolo[vdx];
   const double Vvallololo = vlololo[vdx];
   const double Rvalhihihi = Rhihihi[Rdx];
   const double Rvallohihi = Rlohihi[Rdx];
   const double Rvalhilohi = Rhilohi[Rdx];
   const double Rvallolohi = Rlolohi[Rdx];
   const double Rvalhihilo = Rhihilo[Rdx];
   const double Rvallohilo = Rlohilo[Rdx];
   const double Rvalhilolo = Rhilolo[Rdx];
   const double Rvallololo = Rlololo[Rdx];
   // double result = Rval*Vval;
   double resulthihihi,resultlohihi,resulthilohi,resultlolohi;
   double resulthihilo,resultlohilo,resulthilolo,resultlololo;

   __syncthreads();
   odg_mul(   Rvalhihihi,   Rvallohihi,   Rvalhilohi,   Rvallolohi,
              Rvalhihilo,   Rvallohilo,   Rvalhilolo,   Rvallololo,
              Vvalhihihi,   Vvallohihi,   Vvalhilohi,   Vvallolohi,
              Vvalhihilo,   Vvallohilo,   Vvalhilolo,   Vvallololo,
           &resulthihihi,&resultlohihi,&resulthilohi,&resultlolohi,
           &resulthihilo,&resultlohilo,&resulthilolo,&resultlololo);
   __syncthreads();
   RTdotvhihihi[idx] = resulthihihi;
   RTdotvlohihi[idx] = resultlohihi;
   RTdotvhilohi[idx] = resulthilohi;
   RTdotvlolohi[idx] = resultlolohi;
   RTdotvhihilo[idx] = resulthihilo;
   RTdotvlohilo[idx] = resultlohilo;
   RTdotvhilolo[idx] = resulthilolo;
   RTdotvlololo[idx] = resultlololo;
}

__global__ void cmplx8_RHdotv
 ( int nrows, int szt, int colidx, int Roffset, int dim,
   double *Rrehihihi, double *Rrelohihi, double *Rrehilohi, double *Rrelolohi,
   double *Rrehihilo, double *Rrelohilo, double *Rrehilolo, double *Rrelololo,
   double *Rimhihihi, double *Rimlohihi, double *Rimhilohi, double *Rimlolohi,
   double *Rimhihilo, double *Rimlohilo, double *Rimhilolo, double *Rimlololo,
   double *vrehihihi, double *vrelohihi, double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo, double *vrehilolo, double *vrelololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo, double *vimhilolo, double *vimlololo,
   double *RHdotvrehihihi, double *RHdotvrelohihi,
   double *RHdotvrehilohi, double *RHdotvrelolohi,
   double *RHdotvrehihilo, double *RHdotvrelohilo,
   double *RHdotvrehilolo, double *RHdotvrelololo,
   double *RHdotvimhihihi, double *RHdotvimlohihi,
   double *RHdotvimhilohi, double *RHdotvimlolohi,
   double *RHdotvimhihilo, double *RHdotvimlohilo,
   double *RHdotvimhilolo, double *RHdotvimlololo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int idx = bdx*szt + tdx;        // thread tdx computes RTv[idx]

   const int vdx = idx % nrows;          // index in v is column in R^T
   const int row = idx / nrows;          // R is stored column-by-column

   const int Rdx = Roffset + idx + (row+1)*colidx;

   const double Vvalrehihihi = vrehihihi[vdx];
   const double Vvalrelohihi = vrelohihi[vdx];
   const double Vvalrehilohi = vrehilohi[vdx];
   const double Vvalrelolohi = vrelolohi[vdx];
   const double Vvalrehihilo = vrehihilo[vdx];
   const double Vvalrelohilo = vrelohilo[vdx];
   const double Vvalrehilolo = vrehilolo[vdx];
   const double Vvalrelololo = vrelololo[vdx];
   const double Vvalimhihihi = vimhihihi[vdx];
   const double Vvalimlohihi = vimlohihi[vdx];
   const double Vvalimhilohi = vimhilohi[vdx];
   const double Vvalimlolohi = vimlolohi[vdx];
   const double Vvalimhihilo = vimhihilo[vdx];
   const double Vvalimlohilo = vimlohilo[vdx];
   const double Vvalimhilolo = vimhilolo[vdx];
   const double Vvalimlololo = vimlololo[vdx];
   const double Rvalrehihihi = Rrehihihi[Rdx];
   const double Rvalrelohihi = Rrelohihi[Rdx];
   const double Rvalrehilohi = Rrehilohi[Rdx];
   const double Rvalrelolohi = Rrelolohi[Rdx];
   const double Rvalrehihilo = Rrehihilo[Rdx];
   const double Rvalrelohilo = Rrelohilo[Rdx];
   const double Rvalrehilolo = Rrehilolo[Rdx];
   const double Rvalrelololo = Rrelololo[Rdx];
   const double Rvalimhihihi = Rimhihihi[Rdx];
   const double Rvalimlohihi = Rimlohihi[Rdx];
   const double Rvalimhilohi = Rimhilohi[Rdx];
   const double Rvalimlolohi = Rimlolohi[Rdx];
   const double Rvalimhihilo = Rimhihilo[Rdx];
   const double Rvalimlohilo = Rimlohilo[Rdx];
   const double Rvalimhilolo = Rimhilolo[Rdx];
   const double Rvalimlololo = Rimlololo[Rdx];
   // double result = Rval*Vval;
   double resultrehihihi,resultrelohihi,resultrehilohi,resultrelolohi;
   double resultrehihilo,resultrelohilo,resultrehilolo,resultrelololo;
   double resultimhihihi,resultimlohihi,resultimhilohi,resultimlolohi;
   double resultimhihilo,resultimlohilo,resultimhilolo,resultimlololo;
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   __syncthreads();

   odg_mul(   Rvalrehihihi,   Rvalrelohihi,   Rvalrehilohi,   Rvalrelolohi,
              Rvalrehihilo,   Rvalrelohilo,   Rvalrehilolo,   Rvalrelololo,
              Vvalrehihihi,   Vvalrelohihi,   Vvalrehilohi,   Vvalrelolohi,
              Vvalrehihilo,   Vvalrelohilo,   Vvalrehilolo,   Vvalrelololo,
           &resultrehihihi,&resultrelohihi,&resultrehilohi,&resultrelolohi,
           &resultrehihilo,&resultrelohilo,&resultrehilolo,&resultrelololo);
   odg_mul(Rvalimhihihi,Rvalimlohihi,Rvalimhilohi,Rvalimlolohi,
           Rvalimhihilo,Rvalimlohilo,Rvalimhilolo,Rvalimlololo,
           Vvalimhihihi,Vvalimlohihi,Vvalimhilohi,Vvalimlolohi,
           Vvalimhihilo,Vvalimlohilo,Vvalimhilolo,Vvalimlololo,
             &acchihihi,  &acclohihi,  &acchilohi,  &acclolohi,
             &acchihilo,  &acclohilo,  &acchilolo,  &acclololo);
   odg_inc(&resultrehihihi,&resultrelohihi,&resultrehilohi,&resultrelolohi,
           &resultrehihilo,&resultrelohilo,&resultrehilolo,&resultrelololo,
                 acchihihi,      acclohihi,      acchilohi,      acclolohi,
                 acchihilo,      acclohilo,      acchilolo,      acclololo);
   odg_mul(   Rvalrehihihi,   Rvalrelohihi,   Rvalrehilohi,   Rvalrelolohi,
              Rvalrehihilo,   Rvalrelohilo,   Rvalrehilolo,   Rvalrelololo,
              Vvalimhihihi,   Vvalimlohihi,   Vvalimhilohi,   Vvalimlolohi,
              Vvalimhihilo,   Vvalimlohilo,   Vvalimhilolo,   Vvalimlololo,
           &resultimhihihi,&resultimlohihi,&resultimhilohi,&resultimlolohi,
           &resultimhihilo,&resultimlohilo,&resultimhilolo,&resultimlololo);
   odg_mul(Rvalimhihihi,Rvalimlohihi,Rvalimhilohi,Rvalimlolohi,
           Rvalimhihilo,Rvalimlohilo,Rvalimhilolo,Rvalimlololo,
           Vvalrehihihi,Vvalrelohihi,Vvalrehilohi,Vvalrelolohi,
           Vvalrehihilo,Vvalrelohilo,Vvalrehilolo,Vvalrelololo,
             &acchihihi,  &acclohihi,  &acchilohi,  &acclolohi,
             &acchihilo,  &acclohilo,  &acchilolo,  &acclololo);
   odg_dec(&resultimhihihi,&resultimlohihi,&resultimhilohi,&resultimlolohi,
           &resultimhihilo,&resultimlohilo,&resultimhilolo,&resultimlololo,
                 acchihihi,      acclohihi,      acchilohi,      acclolohi,
                 acchihilo,      acclohilo,      acchilolo,      acclololo);

   __syncthreads();
   RHdotvrehihihi[idx] = resultrehihihi;
   RHdotvrelohihi[idx] = resultrelohihi;
   RHdotvrehilohi[idx] = resultrehilohi;
   RHdotvrelolohi[idx] = resultrelolohi;
   RHdotvrehihilo[idx] = resultrehihilo;
   RHdotvrelohilo[idx] = resultrelohilo;
   RHdotvrehilolo[idx] = resultrehilolo;
   RHdotvrelololo[idx] = resultrelololo;
   RHdotvimhihihi[idx] = resultimhihihi;
   RHdotvimlohihi[idx] = resultimlohihi;
   RHdotvimhilohi[idx] = resultimhilohi;
   RHdotvimlolohi[idx] = resultimlolohi;
   RHdotvimhihilo[idx] = resultimhihilo;
   RHdotvimlohilo[idx] = resultimlohilo;
   RHdotvimhilolo[idx] = resultimhilolo;
   RHdotvimlololo[idx] = resultimlololo;
}

__global__ void cmplx8_RHdotvRe
 ( int nrows, int szt, int colidx, int Roffset, int dim,
   double *Rrehihihi, double *Rrelohihi, double *Rrehilohi, double *Rrelolohi,
   double *Rrehihilo, double *Rrelohilo, double *Rrehilolo, double *Rrelololo,
   double *Rimhihihi, double *Rimlohihi, double *Rimhilohi, double *Rimlolohi,
   double *Rimhihilo, double *Rimlohilo, double *Rimhilolo, double *Rimlololo,
   double *vrehihihi, double *vrelohihi, double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo, double *vrehilolo, double *vrelololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo, double *vimhilolo, double *vimlololo,
   double *RHdotvrehihihi, double *RHdotvrelohihi,
   double *RHdotvrehilohi, double *RHdotvrelolohi,
   double *RHdotvrehihilo, double *RHdotvrelohilo,
   double *RHdotvrehilolo, double *RHdotvrelololo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int idx = bdx*szt + tdx;        // thread tdx computes RTv[idx]

   const int vdx = idx % nrows;          // index in v is column in R^T
   const int row = idx / nrows;          // R is stored column-by-column

   const int Rdx = Roffset + idx + (row+1)*colidx;
/*
   double Vvalhihihi,Vvallohihi,Vvalhilohi,Vvallolohi;
   double Vvalhihilo,Vvallohilo,Vvalhilolo,Vvallololo;
   double Rvalhihihi,Rvallohihi,Rvalhilohi,Rvallolohi;
   double Rvalhihilo,Rvallohilo,Rvalhilolo,Rvallololo;
 */
   double resulthihihi,resultlohihi,resulthilohi,resultlolohi;
   double resulthihilo,resultlohilo,resulthilolo,resultlololo;
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   // double result = Rval*Vval;
 /*
   Vvalhihihi = vrehihihi[vdx];
   Vvallohihi = vrelohihi[vdx];
   Vvalhilohi = vrehilohi[vdx];
   Vvallolohi = vrelolohi[vdx];
   Vvalhihilo = vrehihilo[vdx];
   Vvallohilo = vrelohilo[vdx];
   Vvalhilolo = vrehilolo[vdx];
   Vvallololo = vrelololo[vdx];
   Rvalhihihi = Rrehihihi[Rdx];
   Rvallohihi = Rrelohihi[Rdx];
   Rvalhilohi = Rrehilohi[Rdx];
   Rvallolohi = Rrelolohi[Rdx];
   Rvalhihilo = Rrehihilo[Rdx];
   Rvallohilo = Rrelohilo[Rdx];
   Rvalhilolo = Rrehilolo[Rdx];
   Rvallololo = Rrelololo[Rdx];
  */
   __syncthreads();
   odg_mul(/* Rvalhihihi,   Rvallohihi,   Rvalhilohi,   Rvallolohi,
              Rvalhihilo,   Rvallohilo,   Rvalhilolo,   Rvallololo,
              Vvalhihihi,   Vvallohihi,   Vvalhilohi,   Vvallolohi,
              Vvalhihilo,   Vvallohilo,   Vvalhilolo,   Vvallololo, */
           Rrehihihi[Rdx],Rrelohihi[Rdx],Rrehilohi[Rdx],Rrelolohi[Rdx],
           Rrehihilo[Rdx],Rrelohilo[Rdx],Rrehilolo[Rdx],Rrelololo[Rdx],
           vrehihihi[vdx],vrelohihi[vdx],vrehilohi[vdx],vrelolohi[vdx],
           vrehihilo[vdx],vrelohilo[vdx],vrehilolo[vdx],vrelololo[vdx],
            &resulthihihi, &resultlohihi, &resulthilohi,&resultlolohi,
            &resulthihilo, &resultlohilo, &resulthilolo,&resultlololo);
   __syncthreads();
/*
   Vvalhihihi = vimhihihi[vdx];
   Vvallohihi = vimlohihi[vdx];
   Vvalhilohi = vimhilohi[vdx];
   Vvallolohi = vimlolohi[vdx];
   Vvalhihilo = vimhihilo[vdx];
   Vvallohilo = vimlohilo[vdx];
   Vvalhilolo = vimhilolo[vdx];
   Vvallololo = vimlololo[vdx];
   Rvalhihihi = Rimhihihi[Rdx];
   Rvallohihi = Rimlohihi[Rdx];
   Rvalhilohi = Rimhilohi[Rdx];
   Rvallolohi = Rimlolohi[Rdx];
   Rvalhihilo = Rimhihilo[Rdx];
   Rvallohilo = Rimlohilo[Rdx];
   Rvalhilolo = Rimhilolo[Rdx];
   Rvallololo = Rimlololo[Rdx];
 */
   __syncthreads();
   odg_mul(/* Rvalhihihi,Rvallohihi,Rvalhilohi,Rvallolohi,
              Rvalhihilo,Rvallohilo,Rvalhilolo,Rvallololo,
              Vvalhihihi,Vvallohihi,Vvalhilohi,Vvallolohi,
              Vvalhihilo,Vvallohilo,Vvalhilolo,Vvallololo, */
           Rimhihihi[Rdx],Rimlohihi[Rdx],Rimhilohi[Rdx],Rimlolohi[Rdx],
           Rimhihilo[Rdx],Rimlohilo[Rdx],Rimhilolo[Rdx],Rimlololo[Rdx],
           vimhihihi[vdx],vimlohihi[vdx],vimhilohi[vdx],vimlolohi[vdx],
           vimhihilo[vdx],vimlohilo[vdx],vimhilolo[vdx],vimlololo[vdx],
               &acchihihi,    &acclohihi,    &acchilohi,    &acclolohi,
               &acchihilo,    &acclohilo,    &acchilolo,    &acclololo);
   odg_inc(&resulthihihi,&resultlohihi,&resulthilohi,&resultlolohi,
           &resulthihilo,&resultlohilo,&resulthilolo,&resultlololo,
               acchihihi,    acclohihi,    acchilohi,    acclolohi,
               acchihilo,    acclohilo,    acchilolo,    acclololo);
   __syncthreads();
   RHdotvrehihihi[idx] = resulthihihi;
   RHdotvrelohihi[idx] = resultlohihi;
   RHdotvrehilohi[idx] = resulthilohi;
   RHdotvrelolohi[idx] = resultlolohi;
   RHdotvrehihilo[idx] = resulthihilo;
   RHdotvrelohilo[idx] = resultlohilo;
   RHdotvrehilolo[idx] = resulthilolo;
   RHdotvrelololo[idx] = resultlololo;
}

__global__ void cmplx8_RHdotvReRe
 ( int nrows, int szt, int colidx, int Roffset, int dim,
   double *Rrehihihi, double *Rrelohihi, double *Rrehilohi, double *Rrelolohi,
   double *Rrehihilo, double *Rrelohilo, double *Rrehilolo, double *Rrelololo,
   double *vrehihihi, double *vrelohihi, double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo, double *vrehilolo, double *vrelololo,
   double *RHdotvrehihihi, double *RHdotvrelohihi,
   double *RHdotvrehilohi, double *RHdotvrelolohi,
   double *RHdotvrehihilo, double *RHdotvrelohilo,
   double *RHdotvrehilolo, double *RHdotvrelololo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int idx = bdx*szt + tdx;        // thread tdx computes RTv[idx]

   const int vdx = idx % nrows;          // index in v is column in R^T
   const int row = idx / nrows;          // R is stored column-by-column

   const int Rdx = Roffset + idx + (row+1)*colidx;

   double resulthihihi,resultlohihi,resulthilohi,resultlolohi;
   double resulthihilo,resultlohilo,resulthilolo,resultlololo;

   // double result = Rval*Vval;
   __syncthreads();
   odg_mul(Rrehihihi[Rdx],Rrelohihi[Rdx],Rrehilohi[Rdx],Rrelolohi[Rdx],
           Rrehihilo[Rdx],Rrelohilo[Rdx],Rrehilolo[Rdx],Rrelololo[Rdx],
           vrehihihi[vdx],vrelohihi[vdx],vrehilohi[vdx],vrelolohi[vdx],
           vrehihilo[vdx],vrelohilo[vdx],vrehilolo[vdx],vrelololo[vdx],
            &resulthihihi, &resultlohihi, &resulthilohi,&resultlolohi,
            &resulthihilo, &resultlohilo, &resulthilolo,&resultlololo);
   __syncthreads();
   RHdotvrehihihi[idx] = resulthihihi;
   RHdotvrelohihi[idx] = resultlohihi;
   RHdotvrehilohi[idx] = resulthilohi;
   RHdotvrelolohi[idx] = resultlolohi;
   RHdotvrehihilo[idx] = resulthihilo;
   RHdotvrelohilo[idx] = resultlohilo;
   RHdotvrehilolo[idx] = resulthilolo;
   RHdotvrelololo[idx] = resultlololo;
}

__global__ void cmplx8_RHdotvImIm
 ( int nrows, int szt, int colidx, int Roffset, int dim,
   double *Rimhihihi, double *Rimlohihi, double *Rimhilohi, double *Rimlolohi,
   double *Rimhihilo, double *Rimlohilo, double *Rimhilolo, double *Rimlololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo, double *vimhilolo, double *vimlololo,
   double *RHdotvrehihihi, double *RHdotvrelohihi,
   double *RHdotvrehilohi, double *RHdotvrelolohi,
   double *RHdotvrehihilo, double *RHdotvrelohilo,
   double *RHdotvrehilolo, double *RHdotvrelololo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int idx = bdx*szt + tdx;        // thread tdx computes RTv[idx]

   const int vdx = idx % nrows;          // index in v is column in R^T
   const int row = idx / nrows;          // R is stored column-by-column

   const int Rdx = Roffset + idx + (row+1)*colidx;

   double resulthihihi,resultlohihi,resulthilohi,resultlolohi;
   double resulthihilo,resultlohilo,resulthilolo,resultlololo;
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   // double result = Rval*Vval;

   __syncthreads();
   odg_mul(Rimhihihi[Rdx],Rimlohihi[Rdx],Rimhilohi[Rdx],Rimlolohi[Rdx],
           Rimhihilo[Rdx],Rimlohilo[Rdx],Rimhilolo[Rdx],Rimlololo[Rdx],
           vimhihihi[vdx],vimlohihi[vdx],vimhilohi[vdx],vimlolohi[vdx],
           vimhihilo[vdx],vimlohilo[vdx],vimhilolo[vdx],vimlololo[vdx],
               &acchihihi,    &acclohihi,    &acchilohi,    &acclolohi,
               &acchihilo,    &acclohilo,    &acchilolo,    &acclololo);

   resulthihihi = RHdotvrehihihi[idx];
   resultlohihi = RHdotvrelohihi[idx];
   resulthilohi = RHdotvrehilohi[idx];
   resultlolohi = RHdotvrelolohi[idx];
   resulthihilo = RHdotvrehihilo[idx];
   resultlohilo = RHdotvrelohilo[idx];
   resulthilolo = RHdotvrehilolo[idx];
   resultlololo = RHdotvrelololo[idx];

   odg_inc(&resulthihihi,&resultlohihi,&resulthilohi,&resultlolohi,
           &resulthihilo,&resultlohilo,&resulthilolo,&resultlololo,
               acchihihi,    acclohihi,    acchilohi,    acclolohi,
               acchihilo,    acclohilo,    acchilolo,    acclololo);
   __syncthreads();
   RHdotvrehihihi[idx] = resulthihihi;
   RHdotvrelohihi[idx] = resultlohihi;
   RHdotvrehilohi[idx] = resulthilohi;
   RHdotvrelolohi[idx] = resultlolohi;
   RHdotvrehihilo[idx] = resulthihilo;
   RHdotvrelohilo[idx] = resultlohilo;
   RHdotvrehilolo[idx] = resulthilolo;
   RHdotvrelololo[idx] = resultlololo;
}

__global__ void cmplx8_RHdotvIm
 ( int nrows, int szt, int colidx, int Roffset, int dim,
   double *Rrehihihi, double *Rrelohihi, double *Rrehilohi, double *Rrelolohi,
   double *Rrehihilo, double *Rrelohilo, double *Rrehilolo, double *Rrelololo,
   double *Rimhihihi, double *Rimlohihi, double *Rimhilohi, double *Rimlolohi,
   double *Rimhihilo, double *Rimlohilo, double *Rimhilolo, double *Rimlololo,
   double *vrehihihi, double *vrelohihi, double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo, double *vrehilolo, double *vrelololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo, double *vimhilolo, double *vimlololo,
   double *RHdotvimhihihi, double *RHdotvimlohihi,
   double *RHdotvimhilohi, double *RHdotvimlolohi,
   double *RHdotvimhihilo, double *RHdotvimlohilo,
   double *RHdotvimhilolo, double *RHdotvimlololo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int idx = bdx*szt + tdx;        // thread tdx computes RTv[idx]

   const int vdx = idx % nrows;          // index in v is column in R^T
   const int row = idx / nrows;          // R is stored column-by-column

   const int Rdx = Roffset + idx + (row+1)*colidx;

/*
   double Vvalhihihi,Vvallohihi,Vvalhilohi,Vvallolohi;
   double Vvalhihilo,Vvallohilo,Vvalhilolo,Vvallololo;
   double Rvalhihihi,Rvallohihi,Rvalhilohi,Rvallolohi;
   double Rvalhihilo,Rvallohilo,Rvalhilolo,Rvallololo;
 */
   double resulthihihi,resultlohihi,resulthilohi,resultlolohi;
   double resulthihilo,resultlohilo,resulthilolo,resultlololo;
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   // double result = Rval*Vval;
/*
   Vvalhihihi = vimhihihi[vdx];
   Vvallohihi = vimlohihi[vdx];
   Vvalhilohi = vimhilohi[vdx];
   Vvallolohi = vimlolohi[vdx];
   Vvalhihilo = vimhihilo[vdx];
   Vvallohilo = vimlohilo[vdx];
   Vvalhilolo = vimhilolo[vdx];
   Vvallololo = vimlololo[vdx];
   Rvalhihihi = Rrehihihi[Rdx];
   Rvallohihi = Rrelohihi[Rdx];
   Rvalhilohi = Rrehilohi[Rdx];
   Rvallolohi = Rrelolohi[Rdx];
   Rvalhihilo = Rrehihilo[Rdx];
   Rvallohilo = Rrelohilo[Rdx];
   Rvalhilolo = Rrehilolo[Rdx];
   Rvallololo = Rrelololo[Rdx];
 */
   __syncthreads();
   odg_mul(/* Rvalhihihi,   Rvallohihi,   Rvalhilohi,   Rvallolohi,
              Rvalhihilo,   Rvallohilo,   Rvalhilolo,   Rvallololo,
              Vvalhihihi,   Vvallohihi,   Vvalhilohi,   Vvallolohi,
              Vvalhihilo,   Vvallohilo,   Vvalhilolo,   Vvallololo, */
           Rrehihihi[Rdx],Rrelohihi[Rdx],Rrehilohi[Rdx],Rrelolohi[Rdx],
           Rrehihilo[Rdx],Rrelohilo[Rdx],Rrehilolo[Rdx],Rrelololo[Rdx],
           vimhihihi[vdx],vimlohihi[vdx],vimhilohi[vdx],vimlolohi[vdx],
           vimhihilo[vdx],vimlohilo[vdx],vimhilolo[vdx],vimlololo[vdx],
           &resulthihihi,&resultlohihi,&resulthilohi,&resultlolohi,
           &resulthihilo,&resultlohilo,&resulthilolo,&resultlololo);
   __syncthreads();
/*
   Vvalhihihi = vrehihihi[vdx];
   Vvallohihi = vrelohihi[vdx];
   Vvalhilohi = vrehilohi[vdx];
   Vvallolohi = vrelolohi[vdx];
   Vvalhihilo = vrehihilo[vdx];
   Vvallohilo = vrelohilo[vdx];
   Vvalhilolo = vrehilolo[vdx];
   Vvallololo = vrelololo[vdx];
   Rvalhihihi = Rimhihihi[Rdx];
   Rvallohihi = Rimlohihi[Rdx];
   Rvalhilohi = Rimhilohi[Rdx];
   Rvallolohi = Rimlolohi[Rdx];
   Rvalhihilo = Rimhihilo[Rdx];
   Rvallohilo = Rimlohilo[Rdx];
   Rvalhilolo = Rimhilolo[Rdx];
   Rvallololo = Rimlololo[Rdx];
 */
   __syncthreads();
   odg_mul(/* Rvalhihihi,Rvallohihi,Rvalhilohi,Rvallolohi,
              Rvalhihilo,Rvallohilo,Rvalhilolo,Rvallololo,
              Vvalhihihi,Vvallohihi,Vvalhilohi,Vvallolohi,
              Vvalhihilo,Vvallohilo,Vvalhilolo,Vvallololo, */
           Rimhihihi[Rdx],Rimlohihi[Rdx],Rimhilohi[Rdx],Rimlolohi[Rdx],
           Rimhihilo[Rdx],Rimlohilo[Rdx],Rimhilolo[Rdx],Rimlololo[Rdx],
           vrehihihi[vdx],vrelohihi[vdx],vrehilohi[vdx],vrelolohi[vdx],
           vrehihilo[vdx],vrelohilo[vdx],vrehilolo[vdx],vrelololo[vdx],
           &acchihihi,&acclohihi,&acchilohi,&acclolohi,
           &acchihilo,&acclohilo,&acchilolo,&acclololo);
   odg_dec(&resulthihihi,&resultlohihi,&resulthilohi,&resultlolohi,
           &resulthihilo,&resultlohilo,&resulthilolo,&resultlololo,
               acchihihi,    acclohihi,    acchilohi,    acclolohi,
               acchihilo,    acclohilo,    acchilolo,    acclololo);

   __syncthreads();
   RHdotvimhihihi[idx] = resulthihihi;
   RHdotvimlohihi[idx] = resultlohihi;
   RHdotvimhilohi[idx] = resulthilohi;
   RHdotvimlolohi[idx] = resultlolohi;
   RHdotvimhihilo[idx] = resulthihilo;
   RHdotvimlohilo[idx] = resultlohilo;
   RHdotvimhilolo[idx] = resulthilolo;
   RHdotvimlololo[idx] = resultlololo;
}

__global__ void cmplx8_RHdotvReIm
 ( int nrows, int szt, int colidx, int Roffset, int dim,
   double *Rrehihihi, double *Rrelohihi, double *Rrehilohi, double *Rrelolohi,
   double *Rrehihilo, double *Rrelohilo, double *Rrehilolo, double *Rrelololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo, double *vimhilolo, double *vimlololo,
   double *RHdotvimhihihi, double *RHdotvimlohihi,
   double *RHdotvimhilohi, double *RHdotvimlolohi,
   double *RHdotvimhihilo, double *RHdotvimlohilo,
   double *RHdotvimhilolo, double *RHdotvimlololo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int idx = bdx*szt + tdx;        // thread tdx computes RTv[idx]

   const int vdx = idx % nrows;          // index in v is column in R^T
   const int row = idx / nrows;          // R is stored column-by-column

   const int Rdx = Roffset + idx + (row+1)*colidx;

   double resulthihihi,resultlohihi,resulthilohi,resultlolohi;
   double resulthihilo,resultlohilo,resulthilolo,resultlololo;

   // double result = Rval*Vval;
   __syncthreads();
   odg_mul(Rrehihihi[Rdx],Rrelohihi[Rdx],Rrehilohi[Rdx],Rrelolohi[Rdx],
           Rrehihilo[Rdx],Rrelohilo[Rdx],Rrehilolo[Rdx],Rrelololo[Rdx],
           vimhihihi[vdx],vimlohihi[vdx],vimhilohi[vdx],vimlolohi[vdx],
           vimhihilo[vdx],vimlohilo[vdx],vimhilolo[vdx],vimlololo[vdx],
           &resulthihihi,&resultlohihi,&resulthilohi,&resultlolohi,
           &resulthihilo,&resultlohilo,&resulthilolo,&resultlololo);

   __syncthreads();
   RHdotvimhihihi[idx] = resulthihihi;
   RHdotvimlohihi[idx] = resultlohihi;
   RHdotvimhilohi[idx] = resulthilohi;
   RHdotvimlolohi[idx] = resultlolohi;
   RHdotvimhihilo[idx] = resulthihilo;
   RHdotvimlohilo[idx] = resultlohilo;
   RHdotvimhilolo[idx] = resulthilolo;
   RHdotvimlololo[idx] = resultlololo;
}

__global__ void cmplx8_RHdotvImRe
 ( int nrows, int szt, int colidx, int Roffset, int dim,
   double *Rimhihihi, double *Rimlohihi, double *Rimhilohi, double *Rimlolohi,
   double *Rimhihilo, double *Rimlohilo, double *Rimhilolo, double *Rimlololo,
   double *vrehihihi, double *vrelohihi, double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo, double *vrehilolo, double *vrelololo,
   double *RHdotvimhihihi, double *RHdotvimlohihi,
   double *RHdotvimhilohi, double *RHdotvimlolohi,
   double *RHdotvimhihilo, double *RHdotvimlohilo,
   double *RHdotvimhilolo, double *RHdotvimlololo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int idx = bdx*szt + tdx;        // thread tdx computes RTv[idx]

   const int vdx = idx % nrows;          // index in v is column in R^T
   const int row = idx / nrows;          // R is stored column-by-column

   const int Rdx = Roffset + idx + (row+1)*colidx;

   double resulthihihi,resultlohihi,resulthilohi,resultlolohi;
   double resulthihilo,resultlohilo,resulthilolo,resultlololo;
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   // double result = Rval*Vval;
   __syncthreads();
   odg_mul(Rimhihihi[Rdx],Rimlohihi[Rdx],Rimhilohi[Rdx],Rimlolohi[Rdx],
           Rimhihilo[Rdx],Rimlohilo[Rdx],Rimhilolo[Rdx],Rimlololo[Rdx],
           vrehihihi[vdx],vrelohihi[vdx],vrehilohi[vdx],vrelolohi[vdx],
           vrehihilo[vdx],vrelohilo[vdx],vrehilolo[vdx],vrelololo[vdx],
           &acchihihi,&acclohihi,&acchilohi,&acclolohi,
           &acchihilo,&acclohilo,&acchilolo,&acclololo);

   resulthihihi = RHdotvimhihihi[idx];
   resultlohihi = RHdotvimlohihi[idx];
   resulthilohi = RHdotvimhilohi[idx];
   resultlolohi = RHdotvimlolohi[idx];
   resulthihilo = RHdotvimhihilo[idx];
   resultlohilo = RHdotvimlohilo[idx];
   resulthilolo = RHdotvimhilolo[idx];
   resultlololo = RHdotvimlololo[idx];

   odg_dec(&resulthihihi,&resultlohihi,&resulthilohi,&resultlolohi,
           &resulthihilo,&resultlohilo,&resulthilolo,&resultlololo,
               acchihihi,    acclohihi,    acchilohi,    acclolohi,
               acchihilo,    acclohilo,    acchilolo,    acclololo);

   __syncthreads();
   RHdotvimhihihi[idx] = resulthihihi;
   RHdotvimlohihi[idx] = resultlohihi;
   RHdotvimhilohi[idx] = resulthilohi;
   RHdotvimlolohi[idx] = resultlolohi;
   RHdotvimhihilo[idx] = resulthihilo;
   RHdotvimlohilo[idx] = resultlohilo;
   RHdotvimhilolo[idx] = resulthilolo;
   RHdotvimlololo[idx] = resultlololo;
}

__global__ void dbl8_sum_betaRTdotv
 ( int nrows,
   double *betahihihi, double *betalohihi,
   double *betahilohi, double *betalolohi,
   double *betahihilo, double *betalohilo,
   double *betahilolo, double *betalololo,
   double *RTdotvhihihi, double *RTdotvlohihi,
   double *RTdotvhilohi, double *RTdotvlolohi,
   double *RTdotvhihilo, double *RTdotvlohilo,
   double *RTdotvhilolo, double *RTdotvlololo,
   double *whihihi, double *wlohihi, double *whilohi, double *wlolohi,
   double *whihilo, double *wlohilo, double *whilolo, double *wlololo )
{
   const int tdx = threadIdx.x;  // tdx sums elements on row tdx
   const int offset = tdx*nrows; // number of rows before current row
   int idx;

   double resulthihihi = 0.0;
   double resultlohihi = 0.0;
   double resulthilohi = 0.0;
   double resultlolohi = 0.0;
   double resulthihilo = 0.0;
   double resultlohilo = 0.0;
   double resulthilolo = 0.0;
   double resultlololo = 0.0;
   double Rvalhihihi,Rvallohihi,Rvalhilohi,Rvallolohi;
   double Rvalhihilo,Rvallohilo,Rvalhilolo,Rvallololo;

   for(int i=0; i<nrows; i++)
   {
      idx = offset + i;
      __syncthreads();
      Rvalhihihi = RTdotvhihihi[idx];
      Rvallohihi = RTdotvlohihi[idx];
      Rvalhilohi = RTdotvhilohi[idx];
      Rvallolohi = RTdotvlolohi[idx];
      Rvalhihilo = RTdotvhihilo[idx];
      Rvallohilo = RTdotvlohilo[idx];
      Rvalhilolo = RTdotvhilolo[idx];
      Rvallololo = RTdotvlololo[idx];
      // result = result + Rval;
      __syncthreads();
      odg_inc(&resulthihihi,&resultlohihi,&resulthilohi,&resultlolohi,
              &resulthihilo,&resultlohilo,&resulthilolo,&resultlololo,
                 Rvalhihihi,   Rvallohihi,   Rvalhilohi,   Rvallolohi,
                 Rvalhihilo,   Rvallohilo,   Rvalhilolo,   Rvallololo);
   }
   __syncthreads();
   Rvalhihihi = *betahihihi;
   Rvallohihi = *betalohihi;
   Rvalhilohi = *betahilohi;
   Rvallolohi = *betalolohi;
   Rvalhihilo = *betahihilo;
   Rvallohilo = *betalohilo;
   Rvalhilolo = *betahilolo;
   Rvallololo = *betalololo;
   // w[tdx] = Rval*result;
   __syncthreads();
   odg_mul(  Rvalhihihi,   Rvallohihi,   Rvalhilohi,   Rvallolohi,
             Rvalhihilo,   Rvallohilo,   Rvalhilolo,   Rvallololo,
           resulthihihi, resultlohihi, resulthilohi, resultlolohi,
           resulthihilo, resultlohilo, resulthilolo, resultlololo,
               &whihihi[tdx],&wlohihi[tdx],&whilohi[tdx],&wlolohi[tdx],
               &whihilo[tdx],&wlohilo[tdx],&whilolo[tdx],&wlololo[tdx]);
}

__global__ void cmplx8_sum_betaRHdotv
 ( int nrows,
   double *betahihihi, double *betalohihi,
   double *betahilohi, double *betalolohi,
   double *betahihilo, double *betalohilo,
   double *betahilolo, double *betalololo,
   double *RTdotvrehihihi, double *RTdotvrelohihi,
   double *RTdotvrehilohi, double *RTdotvrelolohi,
   double *RTdotvrehihilo, double *RTdotvrelohilo,
   double *RTdotvrehilolo, double *RTdotvrelololo,
   double *RTdotvimhihihi, double *RTdotvimlohihi,
   double *RTdotvimhilohi, double *RTdotvimlolohi,
   double *RTdotvimhihilo, double *RTdotvimlohilo,
   double *RTdotvimhilolo, double *RTdotvimlololo,
   double *wrehihihi, double *wrelohihi, double *wrehilohi, double *wrelolohi,
   double *wrehihilo, double *wrelohilo, double *wrehilolo, double *wrelololo,
   double *wimhihihi, double *wimlohihi,
   double *wimhilohi, double *wimlolohi,
   double *wimhihilo, double *wimlohilo,
   double *wimhilolo, double *wimlololo )
{
   const int tdx = threadIdx.x;  // tdx sums elements on row tdx
   const int offset = tdx*nrows; // number of rows before current row
   int idx;

   double resultrehihihi = 0.0;
   double resultrelohihi = 0.0;
   double resultrehilohi = 0.0;
   double resultrelolohi = 0.0;
   double resultrehihilo = 0.0;
   double resultrelohilo = 0.0;
   double resultrehilolo = 0.0;
   double resultrelololo = 0.0;
   double resultimhihihi = 0.0;
   double resultimlohihi = 0.0;
   double resultimhilohi = 0.0;
   double resultimlolohi = 0.0;
   double resultimhihilo = 0.0;
   double resultimlohilo = 0.0;
   double resultimhilolo = 0.0;
   double resultimlololo = 0.0;
/*
   double Rvalrehihihi,Rvalrelohihi,Rvalrehilohi,Rvalrelolohi;
   double Rvalrehihilo,Rvalrelohilo,Rvalrehilolo,Rvalrelololo;
   double Rvalimhihihi,Rvalimlohihi,Rvalimhilohi,Rvalimlolohi;
   double Rvalimhihilo,Rvalimlohilo,Rvalimhilolo,Rvalimlololo;
 */
   for(int i=0; i<nrows; i++)
   {
      idx = offset + i;
      __syncthreads();
/*
      Rvalrehihihi = RTdotvrehihihi[idx];
      Rvalrelohihi = RTdotvrelohihi[idx];
      Rvalrehilohi = RTdotvrehilohi[idx];
      Rvalrelolohi = RTdotvrelolohi[idx];
      Rvalrehihilo = RTdotvrehihilo[idx];
      Rvalrelohilo = RTdotvrelohilo[idx];
      Rvalrehilolo = RTdotvrehilolo[idx];
      Rvalrelololo = RTdotvrelololo[idx];
      Rvalimhihihi = RTdotvimhihihi[idx];
      Rvalimlohihi = RTdotvimlohihi[idx];
      Rvalimhilohi = RTdotvimhilohi[idx];
      Rvalimlolohi = RTdotvimlolohi[idx];
      Rvalimhihilo = RTdotvimhihilo[idx];
      Rvalimlohilo = RTdotvimlohilo[idx];
      Rvalimhilolo = RTdotvimhilolo[idx];
      Rvalimlololo = RTdotvimlololo[idx];
 */
      // result = result + Rval;
      __syncthreads();
      odg_inc(&resultrehihihi,&resultrelohihi,&resultrehilohi,&resultrelolohi,
              &resultrehihilo,&resultrelohilo,&resultrehilolo,&resultrelololo,
              RTdotvrehihihi[idx],RTdotvrelohihi[idx],
              RTdotvrehilohi[idx],RTdotvrelolohi[idx],
              RTdotvrehihilo[idx],RTdotvrelohilo[idx],
              RTdotvrehilolo[idx],RTdotvrelololo[idx]);
          //    Rvalrehihihi,   Rvalrelohihi,   Rvalrehilohi,   Rvalrelolohi,
          //    Rvalrehihilo,   Rvalrelohilo,   Rvalrehilolo,   Rvalrelololo);
      odg_inc(&resultimhihihi,&resultimlohihi,&resultimhilohi,&resultimlolohi,
              &resultimhihilo,&resultimlohilo,&resultimhilolo,&resultimlololo,
              RTdotvimhihihi[idx],RTdotvimlohihi[idx],
              RTdotvimhilohi[idx],RTdotvimlolohi[idx],
              RTdotvimhihilo[idx],RTdotvimlohilo[idx],
              RTdotvimhilolo[idx],RTdotvimlololo[idx]);
          //    Rvalimhihihi,   Rvalimlohihi,   Rvalimhilohi,   Rvalimlolohi,
          //    Rvalimhihilo,   Rvalimlohilo,   Rvalimhilolo,   Rvalimlololo);
   }
   __syncthreads();
/*
   Rvalrehihihi = *betahihihi;
   Rvalrelohihi = *betalohihi;
   Rvalrehilohi = *betahilohi;
   Rvalrelolohi = *betalolohi;
   Rvalrehihilo = *betahihilo;
   Rvalrelohilo = *betalohilo;
   Rvalrehilolo = *betahilolo;
   Rvalrelololo = *betalololo;
 */
   // w[tdx] = Rval*result;
   __syncthreads();
   odg_mul(/* Rvalrehihihi,   Rvalrelohihi,   Rvalrehilohi,   Rvalrelolohi,
              Rvalrehihilo,   Rvalrelohilo,   Rvalrehilolo,   Rvalrelololo, */
          betahihihi[0],  betalohihi[0],  betahilohi[0], betalolohi[0],
          betahihilo[0],  betalohilo[0],  betahilolo[0], betalololo[0],
         resultrehihihi, resultrelohihi, resultrehilohi, resultrelolohi,
         resultrehihilo, resultrelohilo, resultrehilolo, resultrelololo,
             &wrehihihi[tdx],&wrelohihi[tdx],&wrehilohi[tdx],&wrelolohi[tdx],
             &wrehihilo[tdx],&wrelohilo[tdx],&wrehilolo[tdx],&wrelololo[tdx]);
   odg_mul(/* Rvalrehihihi,   Rvalrelohihi,   Rvalrehilohi,   Rvalrelolohi,
              Rvalrehihilo,   Rvalrelohilo,   Rvalrehilolo,   Rvalrelololo, */
          betahihihi[0],  betalohihi[0],  betahilohi[0], betalolohi[0],
          betahihilo[0],  betalohilo[0],  betahilolo[0], betalololo[0],
         resultimhihihi, resultimlohihi, resultimhilohi, resultimlolohi,
         resultimhihilo, resultimlohilo, resultimhilolo, resultimlololo,
             &wimhihihi[tdx],&wimlohihi[tdx],&wimhilohi[tdx],&wimlolohi[tdx],
             &wimhihilo[tdx],&wimlohilo[tdx],&wimhilolo[tdx],&wimlololo[tdx]);
}

__global__ void dbl8_medium_subvbetaRTv
 ( int nrows, int ncols, int szt, int k,
   double *Rhihihi, double *Rlohihi, double *Rhilohi, double *Rlolohi,
   double *Rhihilo, double *Rlohilo, double *Rhilolo, double *Rlololo,
   double *vhihihi, double *vlohihi, double *vhilohi, double *vlolohi,
   double *vhihilo, double *vlohilo, double *vhilolo, double *vlololo,
   double *betahihihi, double *betalohihi,
   double *betahilohi, double *betalolohi,
   double *betahihilo, double *betalohilo,
   double *betahilolo, double *betalololo,
   double *whihihi, double *wlohihi, double *whilohi, double *wlolohi,
   double *whihilo, double *wlohilo, double *whilolo, double *wlololo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int Roffset = k*nrows + k;    // start in R
   const int widx = bdx*szt + tdx;     // global thread index 

   const int coldim = ncols - k;       // number of columns in R
   const int bound = coldim*(nrows-k); // bound on Ridx
   const int rowidx = widx / coldim;   // row index
   const int colidx = widx % coldim;   // column index

   const int Ridx = Roffset + nrows*colidx + rowidx;

   __shared__ double shwhihihi[od_shmemsize];  // values in beta*R^T*v
   __shared__ double shwlohihi[od_shmemsize];  // are less in number than szt
   __shared__ double shwhilohi[od_shmemsize];
   __shared__ double shwlolohi[od_shmemsize];
   __shared__ double shwhihilo[od_shmemsize];
   __shared__ double shwlohilo[od_shmemsize]; 
   __shared__ double shwhilolo[od_shmemsize];
   __shared__ double shwlololo[od_shmemsize];
   shwhihihi[tdx] = whihihi[tdx];
   shwlohihi[tdx] = wlohihi[tdx];
   shwhilohi[tdx] = whilohi[tdx];
   shwlolohi[tdx] = wlolohi[tdx];
   shwhihilo[tdx] = whihilo[tdx];
   shwlohilo[tdx] = wlohilo[tdx];
   shwhilolo[tdx] = whilolo[tdx];
   shwlololo[tdx] = wlololo[tdx];
   __syncthreads();

   double Rwidxhihihi = Rhihihi[Ridx];     // number that tdx updates
   double Rwidxlohihi = Rlohihi[Ridx];
   double Rwidxhilohi = Rhilohi[Ridx];
   double Rwidxlolohi = Rlolohi[Ridx];
   double Rwidxhihilo = Rhihilo[Ridx];
   double Rwidxlohilo = Rlohilo[Ridx];
   double Rwidxhilolo = Rhilolo[Ridx];
   double Rwidxlololo = Rlololo[Ridx];
   double vValhihihi = vhihihi[rowidx];    // value in Householder vector
   double vVallohihi = vlohihi[rowidx];
   double vValhilohi = vhilohi[rowidx];
   double vVallolohi = vlolohi[rowidx];
   double vValhihilo = vhihilo[rowidx]; 
   double vVallohilo = vlohilo[rowidx];
   double vValhilolo = vhilolo[rowidx];
   double vVallololo = vlololo[rowidx];
   double wValhihihi = shwhihihi[colidx];  // value in beta*R^T*v
   double wVallohihi = shwlohihi[colidx];
   double wValhilohi = shwhilohi[colidx];
   double wVallolohi = shwlolohi[colidx];
   double wValhihilo = shwhihilo[colidx];
   double wVallohilo = shwlohilo[colidx];
   double wValhilolo = shwhilolo[colidx];
   double wVallololo = shwlololo[colidx];
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;
 
   // Rwidx = Rwidx - vValue*wValue;   // update R[rowidx,colidx]
   __syncthreads();
   odg_mul(vValhihihi,vVallohihi,vValhilohi,vVallolohi,
           vValhihilo,vVallohilo,vValhilolo,vVallololo,
           wValhihihi,wVallohihi,wValhilohi,wVallolohi,
           wValhihilo,wVallohilo,wValhilolo,wVallololo,
           &acchihihi,&acclohihi,&acchilohi,&acclolohi,
           &acchihilo,&acclohilo,&acchilolo,&acclololo);
   odg_dec(&Rwidxhihihi,&Rwidxlohihi,&Rwidxhilohi,&Rwidxlolohi,
           &Rwidxhihilo,&Rwidxlohilo,&Rwidxhilolo,&Rwidxlololo,
              acchihihi,   acclohihi,   acchilohi,   acclolohi,
              acchihilo,   acclohilo,   acchilolo,   acclololo);
   __syncthreads();
   if(widx < bound)                    // if() takes care of padding
   {
      Rhihihi[Ridx] = Rwidxhihihi;
      Rlohihi[Ridx] = Rwidxlohihi;
      Rhilohi[Ridx] = Rwidxhilohi;
      Rlolohi[Ridx] = Rwidxlolohi;
      Rhihilo[Ridx] = Rwidxhihilo;
      Rlohilo[Ridx] = Rwidxlohilo;
      Rhilolo[Ridx] = Rwidxhilolo;
      Rlololo[Ridx] = Rwidxlololo;
   }
}

__global__ void cmplx8_medium_subvbetaRHv
 ( int nrows, int ncols, int szt, int k,
   double *Rrehihihi, double *Rrelohihi, double *Rrehilohi, double *Rrelolohi,
   double *Rrehihilo, double *Rrelohilo, double *Rrehilolo, double *Rrelololo,
   double *Rimhihihi, double *Rimlohihi, double *Rimhilohi, double *Rimlolohi,
   double *Rimhihilo, double *Rimlohilo, double *Rimhilolo, double *Rimlololo,
   double *vrehihihi, double *vrelohihi, double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo, double *vrehilolo, double *vrelololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo, double *vimhilolo, double *vimlololo,
   double *betahihihi, double *betalohihi,
   double *betahilohi, double *betalolohi,
   double *betahihilo, double *betalohilo,
   double *betahilolo, double *betalololo,
   double *wrehihihi, double *wrelohihi, double *wrehilohi, double *wrelolohi,
   double *wrehihilo, double *wrelohilo, double *wrehilolo, double *wrelololo,
   double *wimhihihi, double *wimlohihi,
   double *wimhilohi, double *wimlolohi,
   double *wimhihilo, double *wimlohilo,
   double *wimhilolo, double *wimlololo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int Roffset = k*nrows + k;    // start in R
   const int widx = bdx*szt + tdx;     // global thread index 

   const int coldim = ncols - k;       // number of columns in R
   const int bound = coldim*(nrows-k); // bound on Ridx
   const int rowidx = widx / coldim;   // row index
   const int colidx = widx % coldim;   // column index

   const int Ridx = Roffset + nrows*colidx + rowidx;

   __shared__ double shwrehihihi[cod_shmemsize];  // values in beta*R^H*v are
   __shared__ double shwrelohihi[cod_shmemsize];  // less in number than szt
   __shared__ double shwrehilohi[cod_shmemsize]; 
   __shared__ double shwrelolohi[cod_shmemsize]; 
   __shared__ double shwrehihilo[cod_shmemsize]; 
   __shared__ double shwrelohilo[cod_shmemsize];
   __shared__ double shwrehilolo[cod_shmemsize]; 
   __shared__ double shwrelololo[cod_shmemsize]; 
   __shared__ double shwimhihihi[cod_shmemsize]; 
   __shared__ double shwimlohihi[cod_shmemsize];
   __shared__ double shwimhilohi[cod_shmemsize]; 
   __shared__ double shwimlolohi[cod_shmemsize];
   __shared__ double shwimhihilo[cod_shmemsize]; 
   __shared__ double shwimlohilo[cod_shmemsize];
   __shared__ double shwimhilolo[cod_shmemsize]; 
   __shared__ double shwimlololo[cod_shmemsize];

   shwrehihihi[tdx] = wrehihihi[tdx];
   shwrelohihi[tdx] = wrelohihi[tdx];
   shwrehilohi[tdx] = wrehilohi[tdx];
   shwrelolohi[tdx] = wrelolohi[tdx];
   shwrehihilo[tdx] = wrehihilo[tdx];
   shwrelohilo[tdx] = wrelohilo[tdx];
   shwrehilolo[tdx] = wrehilolo[tdx];
   shwrelololo[tdx] = wrelololo[tdx];
   shwimhihihi[tdx] = wimhihihi[tdx];
   shwimlohihi[tdx] = wimlohihi[tdx];
   shwimhilohi[tdx] = wimhilohi[tdx];
   shwimlolohi[tdx] = wimlolohi[tdx];
   shwimhihilo[tdx] = wimhihilo[tdx];
   shwimlohilo[tdx] = wimlohilo[tdx];
   shwimhilolo[tdx] = wimhilolo[tdx];
   shwimlololo[tdx] = wimlololo[tdx];
   __syncthreads();

   double Rwidxrehihihi = Rrehihihi[Ridx];  // number that tdx updates
   double Rwidxrelohihi = Rrelohihi[Ridx];
   double Rwidxrehilohi = Rrehilohi[Ridx];
   double Rwidxrelolohi = Rrelolohi[Ridx];
   double Rwidxrehihilo = Rrehihilo[Ridx];
   double Rwidxrelohilo = Rrelohilo[Ridx];
   double Rwidxrehilolo = Rrehilolo[Ridx];
   double Rwidxrelololo = Rrelololo[Ridx];
   double Rwidximhihihi = Rimhihihi[Ridx];
   double Rwidximlohihi = Rimlohihi[Ridx];
   double Rwidximhilohi = Rimhilohi[Ridx];
   double Rwidximlolohi = Rimlolohi[Ridx];
   double Rwidximhihilo = Rimhihilo[Ridx];
   double Rwidximlohilo = Rimlohilo[Ridx];
   double Rwidximhilolo = Rimhilolo[Ridx];
   double Rwidximlololo = Rimlololo[Ridx];
   double vValrehihihi = vrehihihi[rowidx];  // value in Householder vector
   double vValrelohihi = vrelohihi[rowidx];
   double vValrehilohi = vrehilohi[rowidx];
   double vValrelolohi = vrelolohi[rowidx];
   double vValrehihilo = vrehihilo[rowidx];
   double vValrelohilo = vrelohilo[rowidx];
   double vValrehilolo = vrehilolo[rowidx];
   double vValrelololo = vrelololo[rowidx];
   double vValimhihihi = vimhihihi[rowidx];
   double vValimlohihi = vimlohihi[rowidx];
   double vValimhilohi = vimhilohi[rowidx];
   double vValimlolohi = vimlolohi[rowidx];
   double vValimhihilo = vimhihilo[rowidx];
   double vValimlohilo = vimlohilo[rowidx];
   double vValimhilolo = vimhilolo[rowidx];
   double vValimlololo = vimlololo[rowidx];
   double wValrehihihi = shwrehihihi[colidx];  // value in beta*R^H*v
   double wValrelohihi = shwrelohihi[colidx];
   double wValrehilohi = shwrehilohi[colidx];
   double wValrelolohi = shwrelolohi[colidx];
   double wValrehihilo = shwrehihilo[colidx];
   double wValrelohilo = shwrelohilo[colidx];
   double wValrehilolo = shwrehilolo[colidx];
   double wValrelololo = shwrelololo[colidx];
   double wValimhihihi = shwimhihihi[colidx];
   double wValimlohihi = shwimlohihi[colidx];
   double wValimhilohi = shwimhilohi[colidx];
   double wValimlolohi = shwimlolohi[colidx];
   double wValimhihilo = shwimhihilo[colidx];
   double wValimlohilo = shwimlohilo[colidx];
   double wValimhilolo = shwimhilolo[colidx];
   double wValimlololo = shwimlololo[colidx];
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;
 
   // Rwidx = Rwidx - vValue*wValue;   // update R[rowidx,colidx]
   // take the Hermitian transpose of w
   __syncthreads();
   odg_mul(vValrehihihi,vValrelohihi,vValrehilohi,vValrelolohi,
           vValrehihilo,vValrelohilo,vValrehilolo,vValrelololo,
           wValrehihihi,wValrelohihi,wValrehilohi,wValrelolohi,
           wValrehihilo,wValrelohilo,wValrehilolo,wValrelololo,
             &acchihihi,  &acclohihi,  &acchilohi,  &acclolohi,
             &acchihilo,  &acclohilo,  &acchilolo,  &acclololo);
   odg_dec(&Rwidxrehihihi,&Rwidxrelohihi,&Rwidxrehilohi,&Rwidxrelolohi,
           &Rwidxrehihilo,&Rwidxrelohilo,&Rwidxrehilolo,&Rwidxrelololo,
                acchihihi,     acclohihi,     acchilohi,     acclolohi,
                acchihilo,     acclohilo,     acchilolo,     acclololo);
   odg_mul(vValimhihihi,vValimlohihi,vValimhilohi,vValimlolohi,
           vValimhihilo,vValimlohilo,vValimhilolo,vValimlololo,
           wValimhihihi,wValimlohihi,wValimhilohi,wValimlolohi,
           wValimhihilo,wValimlohilo,wValimhilolo,wValimlololo,
             &acchihihi,  &acclohihi,  &acchilohi,  &acclolohi,
             &acchihilo,  &acclohilo,  &acchilolo,  &acclololo);
   odg_dec(&Rwidxrehihihi,&Rwidxrelohihi,&Rwidxrehilohi,&Rwidxrelolohi,
           &Rwidxrehihilo,&Rwidxrelohilo,&Rwidxrehilolo,&Rwidxrelololo,
                acchihihi,     acclohihi,     acchilohi,     acclolohi,
                acchihilo,     acclohilo,     acchilolo,     acclololo);
   odg_mul(vValimhihihi,vValimlohihi,vValimhilohi,vValimlolohi,
           vValimhihilo,vValimlohilo,vValimhilolo,vValimlololo,
           wValrehihihi,wValrelohihi,wValrehilohi,wValrelolohi,
           wValrehihilo,wValrelohilo,wValrehilolo,wValrelololo,
             &acchihihi,  &acclohihi,  &acchilohi,  &acclolohi,
             &acchihilo,  &acclohilo,  &acchilolo,  &acclololo);
   odg_dec(&Rwidximhihihi,&Rwidximlohihi,&Rwidximhilohi,&Rwidximlolohi,
           &Rwidximhihilo,&Rwidximlohilo,&Rwidximhilolo,&Rwidximlololo,
                acchihihi,     acclohihi,     acchilohi,     acclolohi,
                acchihilo,     acclohilo,     acchilolo,     acclololo);
   odg_mul(vValrehihihi,vValrelohihi,vValrehilohi,vValrelolohi,
           vValrehihilo,vValrelohilo,vValrehilolo,vValrelololo,
           wValimhihihi,wValimlohihi,wValimhilohi,wValimlolohi,
           wValimhihilo,wValimlohilo,wValimhilolo,wValimlololo,
             &acchihihi,  &acclohihi,  &acchilohi,  &acclolohi,
             &acchihilo,  &acclohilo,  &acchilolo,  &acclololo);
   odg_inc(&Rwidximhihihi,&Rwidximlohihi,&Rwidximhilohi,&Rwidximlolohi,
           &Rwidximhihilo,&Rwidximlohilo,&Rwidximhilolo,&Rwidximlololo,
                acchihihi,     acclohihi,     acchilohi,     acclolohi,
                acchihilo,     acclohilo,     acchilolo,     acclololo);

   __syncthreads();
   if(widx < bound)                    // if() takes care of padding
   {
      Rrehihihi[Ridx] = Rwidxrehihihi;
      Rrelohihi[Ridx] = Rwidxrelohihi;
      Rrehilohi[Ridx] = Rwidxrehilohi;
      Rrelolohi[Ridx] = Rwidxrelolohi;
      Rrehihilo[Ridx] = Rwidxrehihilo;
      Rrelohilo[Ridx] = Rwidxrelohilo;
      Rrehilolo[Ridx] = Rwidxrehilolo;
      Rrelololo[Ridx] = Rwidxrelololo;
      Rimhihihi[Ridx] = Rwidximhihihi;
      Rimlohihi[Ridx] = Rwidximlohihi;
      Rimhilohi[Ridx] = Rwidximhilohi;
      Rimlolohi[Ridx] = Rwidximlolohi;
      Rimhihilo[Ridx] = Rwidximhihilo;
      Rimlohilo[Ridx] = Rwidximlohilo;
      Rimhilolo[Ridx] = Rwidximhilolo;
      Rimlololo[Ridx] = Rwidximlololo;
   }
}

__global__ void cmplx8_medium_subvbetaRHvRe
 ( int nrows, int ncols, int szt, int k,
   double *Rrehihihi, double *Rrelohihi, double *Rrehilohi, double *Rrelolohi,
   double *Rrehihilo, double *Rrelohilo, double *Rrehilolo, double *Rrelololo,
   double *vrehihihi, double *vrelohihi, double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo, double *vrehilolo, double *vrelololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo, double *vimhilolo, double *vimlololo,
   double *wrehihihi, double *wrelohihi, double *wrehilohi, double *wrelolohi,
   double *wrehihilo, double *wrelohilo, double *wrehilolo, double *wrelololo,
   double *wimhihihi, double *wimlohihi,
   double *wimhilohi, double *wimlolohi,
   double *wimhihilo, double *wimlohilo,
   double *wimhilolo, double *wimlololo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int Roffset = k*nrows + k;    // start in R
   const int widx = bdx*szt + tdx;     // global thread index 

   const int coldim = ncols - k;       // number of columns in R
   const int bound = coldim*(nrows-k); // bound on Ridx
   const int rowidx = widx / coldim;   // row index
   const int colidx = widx % coldim;   // column index

   const int Ridx = Roffset + nrows*colidx + rowidx;

   __shared__ double shwrehihihi[cod_shmemsize];  // values in beta*R^H*v are
   __shared__ double shwrelohihi[cod_shmemsize];  // less in number than szt
   __shared__ double shwrehilohi[cod_shmemsize]; 
   __shared__ double shwrelolohi[cod_shmemsize]; 
   __shared__ double shwrehihilo[cod_shmemsize]; 
   __shared__ double shwrelohilo[cod_shmemsize];
   __shared__ double shwrehilolo[cod_shmemsize]; 
   __shared__ double shwrelololo[cod_shmemsize]; 
   __shared__ double shwimhihihi[cod_shmemsize]; 
   __shared__ double shwimlohihi[cod_shmemsize];
   __shared__ double shwimhilohi[cod_shmemsize]; 
   __shared__ double shwimlolohi[cod_shmemsize];
   __shared__ double shwimhihilo[cod_shmemsize]; 
   __shared__ double shwimlohilo[cod_shmemsize];
   __shared__ double shwimhilolo[cod_shmemsize]; 
   __shared__ double shwimlololo[cod_shmemsize];

   shwrehihihi[tdx] = wrehihihi[tdx];
   shwrelohihi[tdx] = wrelohihi[tdx];
   shwrehilohi[tdx] = wrehilohi[tdx];
   shwrelolohi[tdx] = wrelolohi[tdx];
   shwrehihilo[tdx] = wrehihilo[tdx];
   shwrelohilo[tdx] = wrelohilo[tdx];
   shwrehilolo[tdx] = wrehilolo[tdx];
   shwrelololo[tdx] = wrelololo[tdx];
   shwimhihihi[tdx] = wimhihihi[tdx];
   shwimlohihi[tdx] = wimlohihi[tdx];
   shwimhilohi[tdx] = wimhilohi[tdx];
   shwimlolohi[tdx] = wimlolohi[tdx];
   shwimhihilo[tdx] = wimhihilo[tdx];
   shwimlohilo[tdx] = wimlohilo[tdx];
   shwimhilolo[tdx] = wimhilolo[tdx];
   shwimlololo[tdx] = wimlololo[tdx];
   __syncthreads();

   double Rwidxhihihi,Rwidxlohihi,Rwidxhilohi,Rwidxlolohi;
   double Rwidxhihilo,Rwidxlohilo,Rwidxhilolo,Rwidxlololo;
   // value in Householder vector
   // double vValhihihi,vVallohihi,vValhilohi,vVallolohi;
   // double vValhihilo,vVallohilo,vValhilolo,vVallololo;
   // value in beta*R^H*v
   // double wValhihihi,wVallohihi,wValhilohi,wVallolohi;
   // double wValhihilo,wVallohilo,wValhilolo,wVallololo;
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   Rwidxhihihi = Rrehihihi[Ridx];  // number that tdx updates
   Rwidxlohihi = Rrelohihi[Ridx];
   Rwidxhilohi = Rrehilohi[Ridx];
   Rwidxlolohi = Rrelolohi[Ridx];
   Rwidxhihilo = Rrehihilo[Ridx];
   Rwidxlohilo = Rrelohilo[Ridx];
   Rwidxhilolo = Rrehilolo[Ridx];
   Rwidxlololo = Rrelololo[Ridx];
 
   // Rwidx = Rwidx - vValue*wValue;   // update R[rowidx,colidx]
   // take the Hermitian transpose of w
   __syncthreads();
/*
   vValhihihi = vrehihihi[rowidx];  // value in Householder vector
   vVallohihi = vrelohihi[rowidx];
   vValhilohi = vrehilohi[rowidx];
   vVallolohi = vrelolohi[rowidx];
   vValhihilo = vrehihilo[rowidx];
   vVallohilo = vrelohilo[rowidx];
   vValhilolo = vrehilolo[rowidx];
   vVallololo = vrelololo[rowidx];
   wValhihihi = shwrehihihi[colidx];  // value in beta*R^H*v
   wVallohihi = shwrelohihi[colidx];
   wValhilohi = shwrehilohi[colidx];
   wVallolohi = shwrelolohi[colidx];
   wValhihilo = shwrehihilo[colidx];
   wVallohilo = shwrelohilo[colidx];
   wValhilolo = shwrehilolo[colidx];
   wVallololo = shwrelololo[colidx];
 */
   __syncthreads();
   odg_mul(/* vValhihihi,vVallohihi,vValhilohi,vVallolohi,
              vValhihilo,vVallohilo,vValhilolo,vVallololo,
              wValhihihi,wVallohihi,wValhilohi,wVallolohi,
              wValhihilo,wVallohilo,wValhilolo,wVallololo, */
           vrehihihi[rowidx],vrelohihi[rowidx],
           vrehilohi[rowidx],vrelolohi[rowidx],
           vrehihilo[rowidx],vrelohilo[rowidx],
           vrehilolo[rowidx],vrelololo[rowidx],
           shwrehihihi[colidx],shwrelohihi[colidx],
           shwrehilohi[colidx],shwrelolohi[colidx],
           shwrehihilo[colidx],shwrelohilo[colidx],
           shwrehilolo[colidx],shwrelololo[colidx],
           &acchihihi,&acclohihi,&acchilohi,&acclolohi,
           &acchihilo,&acclohilo,&acchilolo,&acclololo);
   odg_dec(&Rwidxhihihi,&Rwidxlohihi,&Rwidxhilohi,&Rwidxlolohi,
           &Rwidxhihilo,&Rwidxlohilo,&Rwidxhilolo,&Rwidxlololo,
              acchihihi,   acclohihi,   acchilohi,   acclolohi,
              acchihilo,   acclohilo,   acchilolo,   acclololo);
/*
   vValhihihi = vimhihihi[rowidx];
   vVallohihi = vimlohihi[rowidx];
   vValhilohi = vimhilohi[rowidx];
   vVallolohi = vimlolohi[rowidx];
   vValhihilo = vimhihilo[rowidx];
   vVallohilo = vimlohilo[rowidx];
   vValhilolo = vimhilolo[rowidx];
   vVallololo = vimlololo[rowidx];
   wValhihihi = shwimhihihi[colidx];
   wVallohihi = shwimlohihi[colidx];
   wValhilohi = shwimhilohi[colidx];
   wVallolohi = shwimlolohi[colidx];
   wValhihilo = shwimhihilo[colidx];
   wVallohilo = shwimlohilo[colidx];
   wValhilolo = shwimhilolo[colidx];
   wVallololo = shwimlololo[colidx];
 */
   __syncthreads();
   odg_mul(/* vValhihihi,vVallohihi,vValhilohi,vVallolohi,
              vValhihilo,vVallohilo,vValhilolo,vVallololo,
              wValhihihi,wVallohihi,wValhilohi,wVallolohi,
              wValhihilo,wVallohilo,wValhilolo,wVallololo, */
           vimhihihi[rowidx],vimlohihi[rowidx],
           vimhilohi[rowidx],vimlolohi[rowidx],
           vimhihilo[rowidx],vimlohilo[rowidx],
           vimhilolo[rowidx],vimlololo[rowidx],
           shwimhihihi[colidx],shwimlohihi[colidx],
           shwimhilohi[colidx],shwimlolohi[colidx],
           shwimhihilo[colidx],shwimlohilo[colidx],
           shwimhilolo[colidx],shwimlololo[colidx],
           &acchihihi,&acclohihi,&acchilohi,&acclolohi,
           &acchihilo,&acclohilo,&acchilolo,&acclololo);
   odg_dec(&Rwidxhihihi,&Rwidxlohihi,&Rwidxhilohi,&Rwidxlolohi,
           &Rwidxhihilo,&Rwidxlohilo,&Rwidxhilolo,&Rwidxlololo,
              acchihihi,   acclohihi,   acchilohi,   acclolohi,
              acchihilo,   acclohilo,   acchilolo,   acclololo);

   __syncthreads();
   if(widx < bound)                    // if() takes care of padding
   {
      Rrehihihi[Ridx] = Rwidxhihihi;
      Rrelohihi[Ridx] = Rwidxlohihi;
      Rrehilohi[Ridx] = Rwidxhilohi;
      Rrelolohi[Ridx] = Rwidxlolohi;
      Rrehihilo[Ridx] = Rwidxhihilo;
      Rrelohilo[Ridx] = Rwidxlohilo;
      Rrehilolo[Ridx] = Rwidxhilolo;
      Rrelololo[Ridx] = Rwidxlololo;
   }
}

__global__ void cmplx8_medium_subvbetaRHvReRe
 ( int nrows, int ncols, int szt, int k,
   double *Rrehihihi, double *Rrelohihi, double *Rrehilohi, double *Rrelolohi,
   double *Rrehihilo, double *Rrelohilo, double *Rrehilolo, double *Rrelololo,
   double *vrehihihi, double *vrelohihi, double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo, double *vrehilolo, double *vrelololo,
   double *wrehihihi, double *wrelohihi, double *wrehilohi, double *wrelolohi,
   double *wrehihilo, double *wrelohilo, double *wrehilolo, double *wrelololo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int Roffset = k*nrows + k;    // start in R
   const int widx = bdx*szt + tdx;     // global thread index 

   const int coldim = ncols - k;       // number of columns in R
   const int bound = coldim*(nrows-k); // bound on Ridx
   const int rowidx = widx / coldim;   // row index
   const int colidx = widx % coldim;   // column index

   const int Ridx = Roffset + nrows*colidx + rowidx;

   __shared__ double shwrehihihi[cod_shmemsize];  // values in beta*R^H*v are
   __shared__ double shwrelohihi[cod_shmemsize];  // less in number than szt
   __shared__ double shwrehilohi[cod_shmemsize]; 
   __shared__ double shwrelolohi[cod_shmemsize]; 
   __shared__ double shwrehihilo[cod_shmemsize]; 
   __shared__ double shwrelohilo[cod_shmemsize];
   __shared__ double shwrehilolo[cod_shmemsize]; 
   __shared__ double shwrelololo[cod_shmemsize]; 

   shwrehihihi[tdx] = wrehihihi[tdx];
   shwrelohihi[tdx] = wrelohihi[tdx];
   shwrehilohi[tdx] = wrehilohi[tdx];
   shwrelolohi[tdx] = wrelolohi[tdx];
   shwrehihilo[tdx] = wrehihilo[tdx];
   shwrelohilo[tdx] = wrelohilo[tdx];
   shwrehilolo[tdx] = wrehilolo[tdx];
   shwrelololo[tdx] = wrelololo[tdx];
   __syncthreads();

   double Rwidxhihihi,Rwidxlohihi,Rwidxhilohi,Rwidxlolohi;
   double Rwidxhihilo,Rwidxlohilo,Rwidxhilolo,Rwidxlololo;
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   Rwidxhihihi = Rrehihihi[Ridx];  // number that tdx updates
   Rwidxlohihi = Rrelohihi[Ridx];
   Rwidxhilohi = Rrehilohi[Ridx];
   Rwidxlolohi = Rrelolohi[Ridx];
   Rwidxhihilo = Rrehihilo[Ridx];
   Rwidxlohilo = Rrelohilo[Ridx];
   Rwidxhilolo = Rrehilolo[Ridx];
   Rwidxlololo = Rrelololo[Ridx];
 
   // Rwidx = Rwidx - vValue*wValue;   // update R[rowidx,colidx]
   // take the Hermitian transpose of w
   __syncthreads();
   odg_mul(vrehihihi[rowidx],vrelohihi[rowidx],
           vrehilohi[rowidx],vrelolohi[rowidx],
           vrehihilo[rowidx],vrelohilo[rowidx],
           vrehilolo[rowidx],vrelololo[rowidx],
           shwrehihihi[colidx],shwrelohihi[colidx],
           shwrehilohi[colidx],shwrelolohi[colidx],
           shwrehihilo[colidx],shwrelohilo[colidx],
           shwrehilolo[colidx],shwrelololo[colidx],
           &acchihihi,&acclohihi,&acchilohi,&acclolohi,
           &acchihilo,&acclohilo,&acchilolo,&acclololo);
   odg_dec(&Rwidxhihihi,&Rwidxlohihi,&Rwidxhilohi,&Rwidxlolohi,
           &Rwidxhihilo,&Rwidxlohilo,&Rwidxhilolo,&Rwidxlololo,
              acchihihi,   acclohihi,   acchilohi,   acclolohi,
              acchihilo,   acclohilo,   acchilolo,   acclololo);

   __syncthreads();
   if(widx < bound)                    // if() takes care of padding
   {
      Rrehihihi[Ridx] = Rwidxhihihi;
      Rrelohihi[Ridx] = Rwidxlohihi;
      Rrehilohi[Ridx] = Rwidxhilohi;
      Rrelolohi[Ridx] = Rwidxlolohi;
      Rrehihilo[Ridx] = Rwidxhihilo;
      Rrelohilo[Ridx] = Rwidxlohilo;
      Rrehilolo[Ridx] = Rwidxhilolo;
      Rrelololo[Ridx] = Rwidxlololo;
   }
}

__global__ void cmplx8_medium_subvbetaRHvImIm
 ( int nrows, int ncols, int szt, int k,
   double *Rrehihihi, double *Rrelohihi, double *Rrehilohi, double *Rrelolohi,
   double *Rrehihilo, double *Rrelohilo, double *Rrehilolo, double *Rrelololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo, double *vimhilolo, double *vimlololo,
   double *wimhihihi, double *wimlohihi, double *wimhilohi, double *wimlolohi,
   double *wimhihilo, double *wimlohilo, double *wimhilolo, double *wimlololo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int Roffset = k*nrows + k;    // start in R
   const int widx = bdx*szt + tdx;     // global thread index 

   const int coldim = ncols - k;       // number of columns in R
   const int bound = coldim*(nrows-k); // bound on Ridx
   const int rowidx = widx / coldim;   // row index
   const int colidx = widx % coldim;   // column index

   const int Ridx = Roffset + nrows*colidx + rowidx;

   __shared__ double shwimhihihi[cod_shmemsize];  // values in beta*R^H*v are
   __shared__ double shwimlohihi[cod_shmemsize];  // less in number than szt
   __shared__ double shwimhilohi[cod_shmemsize]; 
   __shared__ double shwimlolohi[cod_shmemsize];
   __shared__ double shwimhihilo[cod_shmemsize]; 
   __shared__ double shwimlohilo[cod_shmemsize];
   __shared__ double shwimhilolo[cod_shmemsize]; 
   __shared__ double shwimlololo[cod_shmemsize];

   shwimhihihi[tdx] = wimhihihi[tdx];
   shwimlohihi[tdx] = wimlohihi[tdx];
   shwimhilohi[tdx] = wimhilohi[tdx];
   shwimlolohi[tdx] = wimlolohi[tdx];
   shwimhihilo[tdx] = wimhihilo[tdx];
   shwimlohilo[tdx] = wimlohilo[tdx];
   shwimhilolo[tdx] = wimhilolo[tdx];
   shwimlololo[tdx] = wimlololo[tdx];
   __syncthreads();

   double Rwidxhihihi,Rwidxlohihi,Rwidxhilohi,Rwidxlolohi;
   double Rwidxhihilo,Rwidxlohilo,Rwidxhilolo,Rwidxlololo;
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   Rwidxhihihi = Rrehihihi[Ridx];  // number that tdx updates
   Rwidxlohihi = Rrelohihi[Ridx];
   Rwidxhilohi = Rrehilohi[Ridx];
   Rwidxlolohi = Rrelolohi[Ridx];
   Rwidxhihilo = Rrehihilo[Ridx];
   Rwidxlohilo = Rrelohilo[Ridx];
   Rwidxhilolo = Rrehilolo[Ridx];
   Rwidxlololo = Rrelololo[Ridx];
 
   // Rwidx = Rwidx - vValue*wValue;   // update R[rowidx,colidx]
   // take the Hermitian transpose of w
   __syncthreads();
   odg_mul(vimhihihi[rowidx],vimlohihi[rowidx],
           vimhilohi[rowidx],vimlolohi[rowidx],
           vimhihilo[rowidx],vimlohilo[rowidx],
           vimhilolo[rowidx],vimlololo[rowidx],
           shwimhihihi[colidx],shwimlohihi[colidx],
           shwimhilohi[colidx],shwimlolohi[colidx],
           shwimhihilo[colidx],shwimlohilo[colidx],
           shwimhilolo[colidx],shwimlololo[colidx],
           &acchihihi,&acclohihi,&acchilohi,&acclolohi,
           &acchihilo,&acclohilo,&acchilolo,&acclololo);
   odg_dec(&Rwidxhihihi,&Rwidxlohihi,&Rwidxhilohi,&Rwidxlolohi,
           &Rwidxhihilo,&Rwidxlohilo,&Rwidxhilolo,&Rwidxlololo,
              acchihihi,   acclohihi,   acchilohi,   acclolohi,
              acchihilo,   acclohilo,   acchilolo,   acclololo);

   __syncthreads();
   if(widx < bound)                    // if() takes care of padding
   {
      Rrehihihi[Ridx] = Rwidxhihihi;
      Rrelohihi[Ridx] = Rwidxlohihi;
      Rrehilohi[Ridx] = Rwidxhilohi;
      Rrelolohi[Ridx] = Rwidxlolohi;
      Rrehihilo[Ridx] = Rwidxhihilo;
      Rrelohilo[Ridx] = Rwidxlohilo;
      Rrehilolo[Ridx] = Rwidxhilolo;
      Rrelololo[Ridx] = Rwidxlololo;
   }
}

__global__ void cmplx8_medium_subvbetaRHvIm
 ( int nrows, int ncols, int szt, int k,
   double *Rimhihihi, double *Rimlohihi, double *Rimhilohi, double *Rimlolohi,
   double *Rimhihilo, double *Rimlohilo, double *Rimhilolo, double *Rimlololo,
   double *vrehihihi, double *vrelohihi, double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo, double *vrehilolo, double *vrelololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo, double *vimhilolo, double *vimlololo,
   double *wrehihihi, double *wrelohihi, double *wrehilohi, double *wrelolohi,
   double *wrehihilo, double *wrelohilo, double *wrehilolo, double *wrelololo,
   double *wimhihihi, double *wimlohihi,
   double *wimhilohi, double *wimlolohi,
   double *wimhihilo, double *wimlohilo,
   double *wimhilolo, double *wimlololo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int Roffset = k*nrows + k;    // start in R
   const int widx = bdx*szt + tdx;     // global thread index 

   const int coldim = ncols - k;       // number of columns in R
   const int bound = coldim*(nrows-k); // bound on Ridx
   const int rowidx = widx / coldim;   // row index
   const int colidx = widx % coldim;   // column index

   const int Ridx = Roffset + nrows*colidx + rowidx;

   __shared__ double shwrehihihi[cod_shmemsize];  // values in beta*R^H*v are
   __shared__ double shwrelohihi[cod_shmemsize];  // less in number than szt
   __shared__ double shwrehilohi[cod_shmemsize]; 
   __shared__ double shwrelolohi[cod_shmemsize]; 
   __shared__ double shwrehihilo[cod_shmemsize]; 
   __shared__ double shwrelohilo[cod_shmemsize];
   __shared__ double shwrehilolo[cod_shmemsize]; 
   __shared__ double shwrelololo[cod_shmemsize]; 
   __shared__ double shwimhihihi[cod_shmemsize]; 
   __shared__ double shwimlohihi[cod_shmemsize];
   __shared__ double shwimhilohi[cod_shmemsize]; 
   __shared__ double shwimlolohi[cod_shmemsize];
   __shared__ double shwimhihilo[cod_shmemsize]; 
   __shared__ double shwimlohilo[cod_shmemsize];
   __shared__ double shwimhilolo[cod_shmemsize]; 
   __shared__ double shwimlololo[cod_shmemsize];

   shwrehihihi[tdx] = wrehihihi[tdx];
   shwrelohihi[tdx] = wrelohihi[tdx];
   shwrehilohi[tdx] = wrehilohi[tdx];
   shwrelolohi[tdx] = wrelolohi[tdx];
   shwrehihilo[tdx] = wrehihilo[tdx];
   shwrelohilo[tdx] = wrelohilo[tdx];
   shwrehilolo[tdx] = wrehilolo[tdx];
   shwrelololo[tdx] = wrelololo[tdx];
   shwimhihihi[tdx] = wimhihihi[tdx];
   shwimlohihi[tdx] = wimlohihi[tdx];
   shwimhilohi[tdx] = wimhilohi[tdx];
   shwimlolohi[tdx] = wimlolohi[tdx];
   shwimhihilo[tdx] = wimhihilo[tdx];
   shwimlohilo[tdx] = wimlohilo[tdx];
   shwimhilolo[tdx] = wimhilolo[tdx];
   shwimlololo[tdx] = wimlololo[tdx];
   __syncthreads();

   double Rwidxhihihi,Rwidxlohihi,Rwidxhilohi,Rwidxlolohi;
   double Rwidxhihilo,Rwidxlohilo,Rwidxhilolo,Rwidxlololo;
   // value in Householder vector
   // double vValhihihi,vVallohihi,vValhilohi,vVallolohi;
   // double vValhihilo,vVallohilo,vValhilolo,vVallololo;
   // value in beta*R^H*v
   // double wValhihihi,wVallohihi,wValhilohi,wVallolohi;
   // double wValhihilo,wVallohilo,wValhilolo,wVallololo;
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   __syncthreads();
   Rwidxhihihi = Rimhihihi[Ridx];
   Rwidxlohihi = Rimlohihi[Ridx];
   Rwidxhilohi = Rimhilohi[Ridx];
   Rwidxlolohi = Rimlolohi[Ridx];
   Rwidxhihilo = Rimhihilo[Ridx];
   Rwidxlohilo = Rimlohilo[Ridx];
   Rwidxhilolo = Rimhilolo[Ridx];
   Rwidxlololo = Rimlololo[Ridx];

   __syncthreads();
/*
   vValhihihi = vimhihihi[rowidx];
   vVallohihi = vimlohihi[rowidx];
   vValhilohi = vimhilohi[rowidx];
   vVallolohi = vimlolohi[rowidx];
   vValhihilo = vimhihilo[rowidx];
   vVallohilo = vimlohilo[rowidx];
   vValhilolo = vimhilolo[rowidx];
   vVallololo = vimlololo[rowidx];
   wValhihihi = shwrehihihi[colidx];  // value in beta*R^H*v
   wVallohihi = shwrelohihi[colidx];
   wValhilohi = shwrehilohi[colidx];
   wVallolohi = shwrelolohi[colidx];
   wValhihilo = shwrehihilo[colidx];
   wVallohilo = shwrelohilo[colidx];
   wValhilolo = shwrehilolo[colidx];
   wVallololo = shwrelololo[colidx];
 */
   __syncthreads();
   odg_mul(/* vValhihihi,vVallohihi,vValhilohi,vVallolohi,
              vValhihilo,vVallohilo,vValhilolo,vVallololo,
              wValhihihi,wVallohihi,wValhilohi,wVallolohi,
              wValhihilo,wVallohilo,wValhilolo,wVallololo, */
           vimhihihi[rowidx],vimlohihi[rowidx],
           vimhilohi[rowidx],vimlolohi[rowidx],
           vimhihilo[rowidx],vimlohilo[rowidx],
           vimhilolo[rowidx],vimlololo[rowidx],
           shwrehihihi[colidx],shwrelohihi[colidx],
           shwrehilohi[colidx],shwrelolohi[colidx],
           shwrehihilo[colidx],shwrelohilo[colidx],
           shwrehilolo[colidx],shwrelololo[colidx],
           &acchihihi,&acclohihi,&acchilohi,&acclolohi,
           &acchihilo,&acclohilo,&acchilolo,&acclololo);
   odg_dec(&Rwidxhihihi,&Rwidxlohihi,&Rwidxhilohi,&Rwidxlolohi,
           &Rwidxhihilo,&Rwidxlohilo,&Rwidxhilolo,&Rwidxlololo,
              acchihihi,   acclohihi,   acchilohi,   acclolohi,
              acchihilo,   acclohilo,   acchilolo,   acclololo);
/*
   vValhihihi = vrehihihi[rowidx];  // value in Householder vector
   vVallohihi = vrelohihi[rowidx];
   vValhilohi = vrehilohi[rowidx];
   vVallolohi = vrelolohi[rowidx];
   vValhihilo = vrehihilo[rowidx];
   vVallohilo = vrelohilo[rowidx];
   vValhilolo = vrehilolo[rowidx];
   vVallololo = vrelololo[rowidx];
   wValhihihi = shwimhihihi[colidx];
   wVallohihi = shwimlohihi[colidx];
   wValhilohi = shwimhilohi[colidx];
   wVallolohi = shwimlolohi[colidx];
   wValhihilo = shwimhihilo[colidx];
   wVallohilo = shwimlohilo[colidx];
   wValhilolo = shwimhilolo[colidx];
   wVallololo = shwimlololo[colidx];
 */
   __syncthreads();
   odg_mul(/* vValhihihi,vVallohihi,vValhilohi,vVallolohi,
              vValhihilo,vVallohilo,vValhilolo,vVallololo,
              wValhihihi,wVallohihi,wValhilohi,wVallolohi,
              wValhihilo,wVallohilo,wValhilolo,wVallololo, */
           vrehihihi[rowidx],vrelohihi[rowidx],
           vrehilohi[rowidx],vrelolohi[rowidx],
           vrehihilo[rowidx],vrelohilo[rowidx],
           vrehilolo[rowidx],vrelololo[rowidx],
           shwimhihihi[colidx],shwimlohihi[colidx],
           shwimhilohi[colidx],shwimlolohi[colidx],
           shwimhihilo[colidx],shwimlohilo[colidx],
           shwimhilolo[colidx],shwimlololo[colidx],
           &acchihihi,&acclohihi,&acchilohi,&acclolohi,
           &acchihilo,&acclohilo,&acchilolo,&acclololo);
   odg_inc(&Rwidxhihihi,&Rwidxlohihi,&Rwidxhilohi,&Rwidxlolohi,
           &Rwidxhihilo,&Rwidxlohilo,&Rwidxhilolo,&Rwidxlololo,
              acchihihi,   acclohihi,   acchilohi,   acclolohi,
              acchihilo,   acclohilo,   acchilolo,   acclololo);

   __syncthreads();
   if(widx < bound)                    // if() takes care of padding
   {
      Rimhihihi[Ridx] = Rwidxhihihi;
      Rimlohihi[Ridx] = Rwidxlohihi;
      Rimhilohi[Ridx] = Rwidxhilohi;
      Rimlolohi[Ridx] = Rwidxlolohi;
      Rimhihilo[Ridx] = Rwidxhihilo;
      Rimlohilo[Ridx] = Rwidxlohilo;
      Rimhilolo[Ridx] = Rwidxhilolo;
      Rimlololo[Ridx] = Rwidxlololo;
   }
}

__global__ void cmplx8_medium_subvbetaRHvImRe
 ( int nrows, int ncols, int szt, int k,
   double *Rimhihihi, double *Rimlohihi, double *Rimhilohi, double *Rimlolohi,
   double *Rimhihilo, double *Rimlohilo, double *Rimhilolo, double *Rimlololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo, double *vimhilolo, double *vimlololo,
   double *wrehihihi, double *wrelohihi, double *wrehilohi, double *wrelolohi,
   double *wrehihilo, double *wrelohilo, double *wrehilolo, double *wrelololo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int Roffset = k*nrows + k;    // start in R
   const int widx = bdx*szt + tdx;     // global thread index 

   const int coldim = ncols - k;       // number of columns in R
   const int bound = coldim*(nrows-k); // bound on Ridx
   const int rowidx = widx / coldim;   // row index
   const int colidx = widx % coldim;   // column index

   const int Ridx = Roffset + nrows*colidx + rowidx;

   __shared__ double shwrehihihi[cod_shmemsize];  // values in beta*R^H*v are
   __shared__ double shwrelohihi[cod_shmemsize];  // less in number than szt
   __shared__ double shwrehilohi[cod_shmemsize]; 
   __shared__ double shwrelolohi[cod_shmemsize]; 
   __shared__ double shwrehihilo[cod_shmemsize]; 
   __shared__ double shwrelohilo[cod_shmemsize];
   __shared__ double shwrehilolo[cod_shmemsize]; 
   __shared__ double shwrelololo[cod_shmemsize]; 

   shwrehihihi[tdx] = wrehihihi[tdx];
   shwrelohihi[tdx] = wrelohihi[tdx];
   shwrehilohi[tdx] = wrehilohi[tdx];
   shwrelolohi[tdx] = wrelolohi[tdx];
   shwrehihilo[tdx] = wrehihilo[tdx];
   shwrelohilo[tdx] = wrelohilo[tdx];
   shwrehilolo[tdx] = wrehilolo[tdx];
   shwrelololo[tdx] = wrelololo[tdx];
   __syncthreads();

   double Rwidxhihihi,Rwidxlohihi,Rwidxhilohi,Rwidxlolohi;
   double Rwidxhihilo,Rwidxlohilo,Rwidxhilolo,Rwidxlololo;
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   __syncthreads();
   Rwidxhihihi = Rimhihihi[Ridx];
   Rwidxlohihi = Rimlohihi[Ridx];
   Rwidxhilohi = Rimhilohi[Ridx];
   Rwidxlolohi = Rimlolohi[Ridx];
   Rwidxhihilo = Rimhihilo[Ridx];
   Rwidxlohilo = Rimlohilo[Ridx];
   Rwidxhilolo = Rimhilolo[Ridx];
   Rwidxlololo = Rimlololo[Ridx];

   __syncthreads();
   odg_mul(vimhihihi[rowidx],vimlohihi[rowidx],
           vimhilohi[rowidx],vimlolohi[rowidx],
           vimhihilo[rowidx],vimlohilo[rowidx],
           vimhilolo[rowidx],vimlololo[rowidx],
           shwrehihihi[colidx],shwrelohihi[colidx],
           shwrehilohi[colidx],shwrelolohi[colidx],
           shwrehihilo[colidx],shwrelohilo[colidx],
           shwrehilolo[colidx],shwrelololo[colidx],
           &acchihihi,&acclohihi,&acchilohi,&acclolohi,
           &acchihilo,&acclohilo,&acchilolo,&acclololo);
   odg_dec(&Rwidxhihihi,&Rwidxlohihi,&Rwidxhilohi,&Rwidxlolohi,
           &Rwidxhihilo,&Rwidxlohilo,&Rwidxhilolo,&Rwidxlololo,
              acchihihi,   acclohihi,   acchilohi,   acclolohi,
              acchihilo,   acclohilo,   acchilolo,   acclololo);
   __syncthreads();
   if(widx < bound)                    // if() takes care of padding
   {
      Rimhihihi[Ridx] = Rwidxhihihi;
      Rimlohihi[Ridx] = Rwidxlohihi;
      Rimhilohi[Ridx] = Rwidxhilohi;
      Rimlolohi[Ridx] = Rwidxlolohi;
      Rimhihilo[Ridx] = Rwidxhihilo;
      Rimlohilo[Ridx] = Rwidxlohilo;
      Rimhilolo[Ridx] = Rwidxhilolo;
      Rimlololo[Ridx] = Rwidxlololo;
   }
}

__global__ void cmplx8_medium_subvbetaRHvReIm
 ( int nrows, int ncols, int szt, int k,
   double *Rimhihihi, double *Rimlohihi, double *Rimhilohi, double *Rimlolohi,
   double *Rimhihilo, double *Rimlohilo, double *Rimhilolo, double *Rimlololo,
   double *vrehihihi, double *vrelohihi, double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo, double *vrehilolo, double *vrelololo,
   double *wimhihihi, double *wimlohihi, double *wimhilohi, double *wimlolohi,
   double *wimhihilo, double *wimlohilo, double *wimhilolo, double *wimlololo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int Roffset = k*nrows + k;    // start in R
   const int widx = bdx*szt + tdx;     // global thread index 

   const int coldim = ncols - k;       // number of columns in R
   const int bound = coldim*(nrows-k); // bound on Ridx
   const int rowidx = widx / coldim;   // row index
   const int colidx = widx % coldim;   // column index

   const int Ridx = Roffset + nrows*colidx + rowidx;

   __shared__ double shwimhihihi[cod_shmemsize];  // values in beta*R^H*v are
   __shared__ double shwimlohihi[cod_shmemsize];  // less in number than szt
   __shared__ double shwimhilohi[cod_shmemsize]; 
   __shared__ double shwimlolohi[cod_shmemsize];
   __shared__ double shwimhihilo[cod_shmemsize]; 
   __shared__ double shwimlohilo[cod_shmemsize];
   __shared__ double shwimhilolo[cod_shmemsize]; 
   __shared__ double shwimlololo[cod_shmemsize];

   shwimhihihi[tdx] = wimhihihi[tdx];
   shwimlohihi[tdx] = wimlohihi[tdx];
   shwimhilohi[tdx] = wimhilohi[tdx];
   shwimlolohi[tdx] = wimlolohi[tdx];
   shwimhihilo[tdx] = wimhihilo[tdx];
   shwimlohilo[tdx] = wimlohilo[tdx];
   shwimhilolo[tdx] = wimhilolo[tdx];
   shwimlololo[tdx] = wimlololo[tdx];
   __syncthreads();

   double Rwidxhihihi,Rwidxlohihi,Rwidxhilohi,Rwidxlolohi;
   double Rwidxhihilo,Rwidxlohilo,Rwidxhilolo,Rwidxlololo;
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   __syncthreads();
   Rwidxhihihi = Rimhihihi[Ridx];
   Rwidxlohihi = Rimlohihi[Ridx];
   Rwidxhilohi = Rimhilohi[Ridx];
   Rwidxlolohi = Rimlolohi[Ridx];
   Rwidxhihilo = Rimhihilo[Ridx];
   Rwidxlohilo = Rimlohilo[Ridx];
   Rwidxhilolo = Rimhilolo[Ridx];
   Rwidxlololo = Rimlololo[Ridx];

   __syncthreads();
   odg_mul(vrehihihi[rowidx],vrelohihi[rowidx],
           vrehilohi[rowidx],vrelolohi[rowidx],
           vrehihilo[rowidx],vrelohilo[rowidx],
           vrehilolo[rowidx],vrelololo[rowidx],
           shwimhihihi[colidx],shwimlohihi[colidx],
           shwimhilohi[colidx],shwimlolohi[colidx],
           shwimhihilo[colidx],shwimlohilo[colidx],
           shwimhilolo[colidx],shwimlololo[colidx],
           &acchihihi,&acclohihi,&acchilohi,&acclolohi,
           &acchihilo,&acclohilo,&acchilolo,&acclololo);
   odg_inc(&Rwidxhihihi,&Rwidxlohihi,&Rwidxhilohi,&Rwidxlolohi,
           &Rwidxhihilo,&Rwidxlohilo,&Rwidxhilolo,&Rwidxlololo,
              acchihihi,   acclohihi,   acchilohi,   acclolohi,
              acchihilo,   acclohilo,   acchilolo,   acclololo);

   __syncthreads();
   if(widx < bound)                    // if() takes care of padding
   {
      Rimhihihi[Ridx] = Rwidxhihihi;
      Rimlohihi[Ridx] = Rwidxlohihi;
      Rimhilohi[Ridx] = Rwidxhilohi;
      Rimlolohi[Ridx] = Rwidxlolohi;
      Rimhihilo[Ridx] = Rwidxhihilo;
      Rimlohilo[Ridx] = Rwidxlohilo;
      Rimhilolo[Ridx] = Rwidxhilolo;
      Rimlololo[Ridx] = Rwidxlololo;
   }
}

__global__ void dbl8_beta_times_V
 ( int nrows, int szt,
   double *Bhihihi, double *Blohihi, double *Bhilohi, double *Blolohi,
   double *Bhihilo, double *Blohilo, double *Bhilolo, double *Blololo,
   double *Vhihihi, double *Vlohihi, double *Vhilohi, double *Vlolohi,
   double *Vhihilo, double *Vlohilo, double *Vhilolo, double *Vlololo,
   double *Whihihi, double *Wlohihi, double *Whilohi, double *Wlolohi,
   double *Whihilo, double *Wlohilo, double *Whilolo, double *Wlololo )
{
   const int bdx = blockIdx.x;        // index of block
   const int tdx = threadIdx.x;       // index of thread in block
   const int idx = bdx*szt + tdx;     // thread tdx computes W[idx]
   double resulthihihi,resultlohihi,resulthilohi,resultlolohi;
   double resulthihilo,resultlohilo,resulthilolo,resultlololo;

   __shared__ double shvhihihi[od_shmemsize]; // to store a slice of V
   __shared__ double shvlohihi[od_shmemsize];
   __shared__ double shvhilohi[od_shmemsize];
   __shared__ double shvlolohi[od_shmemsize];
   __shared__ double shvhihilo[od_shmemsize]; 
   __shared__ double shvlohilo[od_shmemsize];
   __shared__ double shvhilolo[od_shmemsize];
   __shared__ double shvlololo[od_shmemsize];

   shvhihihi[tdx] = Vhihihi[idx]; // thread tdx loads the data 
   shvlohihi[tdx] = Vlohihi[idx]; // at the global index
   shvhilohi[tdx] = Vhilohi[idx];
   shvlolohi[tdx] = Vlolohi[idx];
   shvhihilo[tdx] = Vhihilo[idx];
   shvlohilo[tdx] = Vlohilo[idx];
   shvhilolo[tdx] = Vhilolo[idx];
   shvlololo[tdx] = Vlololo[idx];

   // result = -B[0]*shv[tdx];
   __syncthreads();
   odg_mul(     -Bhihihi[0],   -Blohihi[0],   -Bhilohi[0],   -Blolohi[0],
                -Bhihilo[0],   -Blohilo[0],   -Bhilolo[0],   -Blololo[0],
               shvhihihi[tdx],shvlohihi[tdx],shvhilohi[tdx],shvlolohi[tdx],
               shvhihilo[tdx],shvlohilo[tdx],shvhilolo[tdx],shvlololo[tdx],
           &resulthihihi, &resultlohihi, &resulthilohi, &resultlolohi,
           &resulthihilo, &resultlohilo, &resulthilolo, &resultlololo);

   __syncthreads();
   if(idx < nrows)
   {
      Whihihi[idx] = resulthihihi;
      Wlohihi[idx] = resultlohihi;
      Whilohi[idx] = resulthilohi;
      Wlolohi[idx] = resultlolohi;
      Whihilo[idx] = resulthihilo;
      Wlohilo[idx] = resultlohilo;
      Whilolo[idx] = resulthilolo;
      Wlololo[idx] = resultlololo;
   }
}

__global__ void cmplx8_beta_times_V
 ( int nrows, int szt,
   double *Bhihihi, double *Blohihi, double *Bhilohi, double *Blolohi,
   double *Bhihilo, double *Blohilo, double *Bhilolo, double *Blololo,
   double *Vrehihihi, double *Vrelohihi, double *Vrehilohi, double *Vrelolohi,
   double *Vrehihilo, double *Vrelohilo, double *Vrehilolo, double *Vrelololo,
   double *Vimhihihi, double *Vimlohihi, double *Vimhilohi, double *Vimlolohi,
   double *Vimhihilo, double *Vimlohilo, double *Vimhilolo, double *Vimlololo,
   double *Wrehihihi, double *Wrelohihi, double *Wrehilohi, double *Wrelolohi,
   double *Wrehihilo, double *Wrelohilo, double *Wrehilolo, double *Wrelololo,
   double *Wimhihihi, double *Wimlohihi,
   double *Wimhilohi, double *Wimlolohi,
   double *Wimhihilo, double *Wimlohilo,
   double *Wimhilolo, double *Wimlololo )
{
   const int bdx = blockIdx.x;        // index of block
   const int tdx = threadIdx.x;       // index of thread in block
   const int idx = bdx*szt + tdx;     // thread tdx computes W[idx]
   double resulthihihi,resultlohihi,resulthilohi,resultlolohi;
   double resulthihilo,resultlohilo,resulthilolo,resultlololo;

   const double minBhihihi = -Bhihihi[0];
   const double minBlohihi = -Blohihi[0];
   const double minBhilohi = -Bhilohi[0];
   const double minBlolohi = -Blolohi[0];
   const double minBhihilo = -Bhihilo[0];
   const double minBlohilo = -Blohilo[0];
   const double minBhilolo = -Bhilolo[0];
   const double minBlololo = -Blololo[0];

   __shared__ double shvrehihihi[cod_shmemsize]; // to store a slice of V
   __shared__ double shvrelohihi[cod_shmemsize];
   __shared__ double shvrehilohi[cod_shmemsize];
   __shared__ double shvrelolohi[cod_shmemsize];
   __shared__ double shvrehihilo[cod_shmemsize];
   __shared__ double shvrelohilo[cod_shmemsize];
   __shared__ double shvrehilolo[cod_shmemsize];
   __shared__ double shvrelololo[cod_shmemsize];
   __shared__ double shvimhihihi[cod_shmemsize];
   __shared__ double shvimlohihi[cod_shmemsize];
   __shared__ double shvimhilohi[cod_shmemsize];
   __shared__ double shvimlolohi[cod_shmemsize];
   __shared__ double shvimhihilo[cod_shmemsize];
   __shared__ double shvimlohilo[cod_shmemsize];
   __shared__ double shvimhilolo[cod_shmemsize];
   __shared__ double shvimlololo[cod_shmemsize];

   shvrehihihi[tdx] = Vrehihihi[idx]; // thread tdx loads data
   shvrelohihi[tdx] = Vrelohihi[idx]; // at the global index
   shvrehilohi[tdx] = Vrehilohi[idx];
   shvrelolohi[tdx] = Vrelolohi[idx];
   shvrehihilo[tdx] = Vrehihilo[idx];
   shvrelohilo[tdx] = Vrelohilo[idx];
   shvrehilolo[tdx] = Vrehilolo[idx];
   shvrelololo[tdx] = Vrelololo[idx];

   // resultre = -B[0]*shvre[tdx];
   __syncthreads();
   odg_mul(minBhihihi,      minBlohihi,      minBhilohi,      minBlolohi,
           minBhihilo,      minBlohilo,      minBhilolo,      minBlololo,
          shvrehihihi[tdx],shvrelohihi[tdx],shvrehilohi[tdx],shvrelolohi[tdx],
          shvrehihilo[tdx],shvrelohilo[tdx],shvrehilolo[tdx],shvrelololo[tdx],
        &resulthihihi,   &resultlohihi,   &resulthilohi,   &resultlolohi,
        &resulthihilo,   &resultlohilo,   &resulthilolo,   &resultlololo);
   __syncthreads();
   if(idx < nrows)
   {
      Wrehihihi[idx] = resulthihihi;
      Wrelohihi[idx] = resultlohihi;
      Wrehilohi[idx] = resulthilohi;
      Wrelolohi[idx] = resultlolohi;
      Wrehihilo[idx] = resulthihilo;
      Wrelohilo[idx] = resultlohilo;
      Wrehilolo[idx] = resulthilolo;
      Wrelololo[idx] = resultlololo;
   }
   // resultim = -B[0]*shvim[tdx];
   __syncthreads();
   shvimhihihi[tdx] = Vimhihihi[idx];
   shvimlohihi[tdx] = Vimlohihi[idx];
   shvimhilohi[tdx] = Vimhilohi[idx];
   shvimlolohi[tdx] = Vimlolohi[idx];
   shvimhihilo[tdx] = Vimhihilo[idx];
   shvimlohilo[tdx] = Vimlohilo[idx];
   shvimhilolo[tdx] = Vimhilolo[idx];
   shvimlololo[tdx] = Vimlololo[idx];
   __syncthreads();
   odg_mul(minBhihihi,      minBlohihi,      minBhilohi,      minBlolohi,
           minBhihilo,      minBlohilo,      minBhilolo,      minBlololo,
          shvimhihihi[tdx],shvimlohihi[tdx],shvimhilohi[tdx],shvimlolohi[tdx],
          shvimhihilo[tdx],shvimlohilo[tdx],shvimhilolo[tdx],shvimlololo[tdx],
        &resulthihihi,   &resultlohihi,   &resulthilohi,   &resultlolohi,
        &resulthihilo,   &resultlohilo,   &resulthilolo,   &resultlololo);
   __syncthreads();
   if(idx < nrows)
   {
      Wimhihihi[idx] = resulthihihi;
      Wimlohihi[idx] = resultlohihi;
      Wimhilohi[idx] = resulthilohi;
      Wimlolohi[idx] = resultlolohi;
      Wimhihilo[idx] = resulthihilo;
      Wimlohilo[idx] = resultlohilo;
      Wimhilolo[idx] = resulthilolo;
      Wimlololo[idx] = resultlololo;
   }
}

__global__ void dbl8_initialize_WYT
 ( int dim, int szt,
   double *Vhihihi, double *Vlohihi, double *Vhilohi, double *Vlolohi,
   double *Vhihilo, double *Vlohilo, double *Vhilolo, double *Vlololo,
   double *Whihihi, double *Wlohihi, double *Whilohi, double *Wlolohi,
   double *Whihilo, double *Wlohilo, double *Whilolo, double *Wlololo,
   double *WYThihihi, double *WYTlohihi,
   double *WYThilohi, double *WYTlolohi,
   double *WYThihilo, double *WYTlohilo,
   double *WYThilolo, double *WYTlololo )
{
   const int bdx = blockIdx.x;        // index of block
   const int tdx = threadIdx.x;       // index of thread in block
   const int idx = bdx*szt + tdx;     // global index of the thread
   const int row = idx / dim;         // row index in WYT
   const int col = idx % dim;         // column index in WYT

   const double Vvalhihihi = Vhihihi[col];
   const double Vvallohihi = Vlohihi[col];
   const double Vvalhilohi = Vhilohi[col];
   const double Vvallolohi = Vlolohi[col];
   const double Vvalhihilo = Vhihilo[col];
   const double Vvallohilo = Vlohilo[col];
   const double Vvalhilolo = Vhilolo[col];
   const double Vvallololo = Vlololo[col];
   const double Wvalhihihi = Whihihi[row];
   const double Wvallohihi = Wlohihi[row];
   const double Wvalhilohi = Whilohi[row];
   const double Wvallolohi = Wlolohi[row];
   const double Wvalhihilo = Whihilo[row];
   const double Wvallohilo = Wlohilo[row];
   const double Wvalhilolo = Whilolo[row];
   const double Wvallololo = Wlololo[row];
   // const double result = Vval*Wval;
   double resulthihihi,resultlohihi,resulthilohi,resultlolohi;
   double resulthihilo,resultlohilo,resulthilolo,resultlololo;

   __syncthreads();
   odg_mul(   Vvalhihihi,   Vvallohihi,   Vvalhilohi,   Vvallolohi,
              Vvalhihilo,   Vvallohilo,   Vvalhilolo,   Vvallololo,
              Wvalhihihi,   Wvallohihi,   Wvalhilohi,   Wvallolohi,
              Wvalhihilo,   Wvallohilo,   Wvalhilolo,   Wvallololo,
           &resulthihihi,&resultlohihi,&resulthilohi,&resultlolohi,
           &resulthihilo,&resultlohilo,&resulthilolo,&resultlololo);

   __syncthreads();
   if(idx < dim*dim)
   {
      WYThihihi[idx] = resulthihihi;
      WYTlohihi[idx] = resultlohihi;
      WYThilohi[idx] = resulthilohi;
      WYTlolohi[idx] = resultlolohi;
      WYThihilo[idx] = resulthihilo;
      WYTlohilo[idx] = resultlohilo;
      WYThilolo[idx] = resulthilolo;
      WYTlololo[idx] = resultlololo;
   }
}

__global__ void cmplx8_initialize_WYH
 ( int dim, int szt,
   double *Vrehihihi, double *Vrelohihi, double *Vrehilohi, double *Vrelolohi,
   double *Vrehihilo, double *Vrelohilo, double *Vrehilolo, double *Vrelololo,
   double *Vimhihihi, double *Vimlohihi, double *Vimhilohi, double *Vimlolohi,
   double *Vimhihilo, double *Vimlohilo, double *Vimhilolo, double *Vimlololo,
   double *Wrehihihi, double *Wrelohihi, double *Wrehilohi, double *Wrelolohi,
   double *Wrehihilo, double *Wrelohilo, double *Wrehilolo, double *Wrelololo,
   double *Wimhihihi, double *Wimlohihi, double *Wimhilohi, double *Wimlolohi,
   double *Wimhihilo, double *Wimlohilo, double *Wimhilolo, double *Wimlololo,
   double *WYTrehihihi, double *WYTrelohihi,
   double *WYTrehilohi, double *WYTrelolohi,
   double *WYTrehihilo, double *WYTrelohilo,
   double *WYTrehilolo, double *WYTrelololo,
   double *WYTimhihihi, double *WYTimlohihi,
   double *WYTimhilohi, double *WYTimlolohi,
   double *WYTimhihilo, double *WYTimlohilo,
   double *WYTimhilolo, double *WYTimlololo )
{
   const int bdx = blockIdx.x;        // index of block
   const int tdx = threadIdx.x;       // index of thread in block
   const int idx = bdx*szt + tdx;     // global index of the thread
   const int row = idx / dim;         // row index in WYH
   const int col = idx % dim;         // column index in WYH

   const double Vvalrehihihi = Vrehihihi[col];
   const double Vvalrelohihi = Vrelohihi[col];
   const double Vvalrehilohi = Vrehilohi[col];
   const double Vvalrelolohi = Vrelolohi[col];
   const double Vvalrehihilo = Vrehihilo[col];
   const double Vvalrelohilo = Vrelohilo[col];
   const double Vvalrehilolo = Vrehilolo[col];
   const double Vvalrelololo = Vrelololo[col];
   const double Vvalimhihihi = Vimhihihi[col];
   const double Vvalimlohihi = Vimlohihi[col];
   const double Vvalimhilohi = Vimhilohi[col];
   const double Vvalimlolohi = Vimlolohi[col];
   const double Vvalimhihilo = Vimhihilo[col];
   const double Vvalimlohilo = Vimlohilo[col];
   const double Vvalimhilolo = Vimhilolo[col];
   const double Vvalimlololo = Vimlololo[col];
   const double Wvalrehihihi = Wrehihihi[row];
   const double Wvalrelohihi = Wrelohihi[row];
   const double Wvalrehilohi = Wrehilohi[row];
   const double Wvalrelolohi = Wrelolohi[row];
   const double Wvalrehihilo = Wrehihilo[row];
   const double Wvalrelohilo = Wrelohilo[row];
   const double Wvalrehilolo = Wrehilolo[row];
   const double Wvalrelololo = Wrelololo[row];
   const double Wvalimhihihi = Wimhihihi[row];
   const double Wvalimlohihi = Wimlohihi[row];
   const double Wvalimhilohi = Wimhilohi[row];
   const double Wvalimlolohi = Wimlolohi[row];
   const double Wvalimhihilo = Wimhihilo[row];
   const double Wvalimlohilo = Wimlohilo[row];
   const double Wvalimhilolo = Wimhilolo[row];
   const double Wvalimlololo = Wimlololo[row];
   // const double result = Vval*Wval;
   double resulthihihi,resultlohihi,resulthilohi,resultlolohi;
   double resulthihilo,resultlohilo,resulthilolo,resultlololo;
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   // take the Hermitian transpose of V
   __syncthreads();
   odg_mul( Vvalrehihihi, Vvalrelohihi, Vvalrehilohi, Vvalrelolohi,
            Vvalrehihilo, Vvalrelohilo, Vvalrehilolo, Vvalrelololo,
            Wvalrehihihi, Wvalrelohihi, Wvalrehilohi, Wvalrelolohi,
            Wvalrehihilo, Wvalrelohilo, Wvalrehilolo, Wvalrelololo,
           &resulthihihi,&resultlohihi,&resulthilohi,&resultlolohi,
           &resulthihilo,&resultlohilo,&resulthilolo,&resultlololo);
   odg_mul(Vvalimhihihi,Vvalimlohihi,Vvalimhilohi,Vvalimlolohi,
           Vvalimhihilo,Vvalimlohilo,Vvalimhilolo,Vvalimlololo,
           Wvalimhihihi,Wvalimlohihi,Wvalimhilohi,Wvalimlolohi,
           Wvalimhihilo,Wvalimlohilo,Wvalimhilolo,Wvalimlololo,
             &acchihihi,  &acclohihi,  &acchilohi,  &acclolohi,
             &acchihilo,  &acclohilo,  &acchilolo,  &acclololo);
   odg_inc(&resulthihihi,&resultlohihi,&resulthilohi,&resultlolohi,
           &resulthihilo,&resultlohilo,&resulthilolo,&resultlololo,
               acchihihi,    acclohihi,    acchilohi,    acclolohi,
               acchihilo,    acclohilo,    acchilolo,    acclololo);

   __syncthreads();
   if(idx < dim*dim)
   {
      WYTrehihihi[idx] = resulthihihi;
      WYTrelohihi[idx] = resultlohihi;
      WYTrehilohi[idx] = resulthilohi;
      WYTrelolohi[idx] = resultlolohi;
      WYTrehihilo[idx] = resulthihilo;
      WYTrelohilo[idx] = resultlohilo;
      WYTrehilolo[idx] = resulthilolo;
      WYTrelololo[idx] = resultlololo;
   }
   __syncthreads();
   odg_mul( Vvalrehihihi, Vvalrelohihi, Vvalrehilohi, Vvalrelolohi,
            Vvalrehihilo, Vvalrelohilo, Vvalrehilolo, Vvalrelololo,
            Wvalimhihihi, Wvalimlohihi, Wvalimhilohi, Wvalimlolohi,
            Wvalimhihilo, Wvalimlohilo, Wvalimhilolo, Wvalimlololo,
           &resulthihihi,&resultlohihi,&resulthilohi,&resultlolohi,
           &resulthihilo,&resultlohilo,&resulthilolo,&resultlololo);
   odg_mul(Vvalimhihihi,Vvalimlohihi,Vvalimhilohi,Vvalimlolohi,
           Vvalimhihilo,Vvalimlohilo,Vvalimhilolo,Vvalimlololo,
           Wvalrehihihi,Wvalrelohihi,Wvalrehilohi,Wvalrelolohi,
           Wvalrehihilo,Wvalrelohilo,Wvalrehilolo,Wvalrelololo,
             &acchihihi,  &acclohihi,  &acchilohi,  &acclolohi,
             &acchihilo,  &acclohilo,  &acchilolo,  &acclololo);
   odg_dec(&resulthihihi,&resultlohihi,&resulthilohi,&resultlolohi,
           &resulthihilo,&resultlohilo,&resulthilolo,&resultlololo,
               acchihihi,    acclohihi,    acchilohi,    acclolohi,
               acchihilo,    acclohilo,    acchilolo,    acclololo);

   __syncthreads();
   if(idx < dim*dim)
   {
      WYTimhihihi[idx] = resulthihihi;
      WYTimlohihi[idx] = resultlohihi;
      WYTimhilohi[idx] = resulthilohi;
      WYTimlolohi[idx] = resultlolohi;
      WYTimhihilo[idx] = resulthihilo;
      WYTimlohilo[idx] = resultlohilo;
      WYTimhilolo[idx] = resulthilolo;
      WYTimlololo[idx] = resultlololo;
   }
}

__global__ void dbl8_update_WYT
 ( int dim, int szt,
   double *Vhihihi, double *Vlohihi, double *Vhilohi, double *Vlolohi,
   double *Vhihilo, double *Vlohilo, double *Vhilolo, double *Vlololo,
   double *Whihihi, double *Wlohihi, double *Whilohi, double *Wlolohi,
   double *Whihilo, double *Wlohilo, double *Whilolo, double *Wlololo,
   double *WYThihihi, double *WYTlohihi,
   double *WYThilohi, double *WYTlolohi,
   double *WYThihilo, double *WYTlohilo,
   double *WYThilolo, double *WYTlololo )
{
   const int bdx = blockIdx.x;        // index of block
   const int tdx = threadIdx.x;       // index of thread in block
   const int idx = bdx*szt + tdx;     // global index of the thread
   const int row = idx / dim;         // row index in WYT
   const int col = idx % dim;         // column index in WYT

   const double Vvalhihihi = Vhihihi[col];
   const double Vvallohihi = Vlohihi[col];
   const double Vvalhilohi = Vhilohi[col];
   const double Vvallolohi = Vlolohi[col];
   const double Vvalhihilo = Vhihilo[col];
   const double Vvallohilo = Vlohilo[col];
   const double Vvalhilolo = Vhilolo[col];
   const double Vvallololo = Vlololo[col];
   __syncthreads();
   const double Wvalhihihi = Whihihi[row];
   const double Wvallohihi = Wlohihi[row];
   const double Wvalhilohi = Whilohi[row];
   const double Wvallolohi = Wlolohi[row];
   const double Wvalhihilo = Whihilo[row];
   const double Wvallohilo = Wlohilo[row];
   const double Wvalhilolo = Whilolo[row];
   const double Wvallololo = Wlololo[row];
   __syncthreads();
   double resulthihihi = WYThihihi[idx];
   double resultlohihi = WYTlohihi[idx];
   double resulthilohi = WYThilohi[idx];
   double resultlolohi = WYTlolohi[idx];
   double resulthihilo = WYThihilo[idx];
   double resultlohilo = WYTlohilo[idx];
   double resulthilolo = WYThilolo[idx];
   double resultlololo = WYTlololo[idx];
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   // result = result + Vval*Wval;

   __syncthreads();
   odg_mul(Vvalhihihi,Vvallohihi,Vvalhilohi,Vvallolohi,
           Vvalhihilo,Vvallohilo,Vvalhilolo,Vvallololo,
           Wvalhihihi,Wvallohihi,Wvalhilohi,Wvallolohi,
           Wvalhihilo,Wvallohilo,Wvalhilolo,Wvallololo,
           &acchihihi,&acclohihi,&acchilohi,&acclolohi,
           &acchihilo,&acclohilo,&acchilolo,&acclololo);
   odg_inc(&resulthihihi,&resultlohihi,&resulthilohi,&resultlolohi,
           &resulthihilo,&resultlohilo,&resulthilolo,&resultlololo,
               acchihihi,    acclohihi,    acchilohi,    acclolohi,
               acchihilo,    acclohilo,    acchilolo,    acclololo);
   
   __syncthreads();
   if(idx < dim*dim)
   {
      WYThihihi[idx] = resulthihihi;
      WYTlohihi[idx] = resultlohihi;
      WYThilohi[idx] = resulthilohi;
      WYTlolohi[idx] = resultlolohi;
      WYThihilo[idx] = resulthihilo;
      WYTlohilo[idx] = resultlohilo;
      WYThilolo[idx] = resulthilolo;
      WYTlololo[idx] = resultlololo;
   }
}

__global__ void cmplx8_update_WYH
 ( int dim, int szt,
   double *Vrehihihi, double *Vrelohihi, double *Vrehilohi, double *Vrelolohi,
   double *Vrehihilo, double *Vrelohilo, double *Vrehilolo, double *Vrelololo,
   double *Vimhihihi, double *Vimlohihi, double *Vimhilohi, double *Vimlolohi,
   double *Vimhihilo, double *Vimlohilo, double *Vimhilolo, double *Vimlololo,
   double *Wrehihihi, double *Wrelohihi, double *Wrehilohi, double *Wrelolohi,
   double *Wrehihilo, double *Wrelohilo, double *Wrehilolo, double *Wrelololo,
   double *Wimhihihi, double *Wimlohihi, double *Wimhilohi, double *Wimlolohi,
   double *Wimhihilo, double *Wimlohilo, double *Wimhilolo, double *Wimlololo,
   double *WYHrehihihi, double *WYHrelohihi,
   double *WYHrehilohi, double *WYHrelolohi,
   double *WYHrehihilo, double *WYHrelohilo,
   double *WYHrehilolo, double *WYHrelololo,
   double *WYHimhihihi, double *WYHimlohihi,
   double *WYHimhilohi, double *WYHimlolohi,
   double *WYHimhihilo, double *WYHimlohilo,
   double *WYHimhilolo, double *WYHimlololo )
{
   const int bdx = blockIdx.x;        // index of block
   const int tdx = threadIdx.x;       // index of thread in block
   const int idx = bdx*szt + tdx;     // global index of the thread
   const int row = idx / dim;         // row index in WYH
   const int col = idx % dim;         // column index in WYH

   const double Vvalrehihihi = Vrehihihi[col];
   const double Vvalrelohihi = Vrelohihi[col];
   const double Vvalrehilohi = Vrehilohi[col];
   const double Vvalrelolohi = Vrelolohi[col];
   const double Vvalrehihilo = Vrehihilo[col];
   const double Vvalrelohilo = Vrelohilo[col];
   const double Vvalrehilolo = Vrehilolo[col];
   const double Vvalrelololo = Vrelololo[col];
   const double Vvalimhihihi = Vimhihihi[col];
   const double Vvalimlohihi = Vimlohihi[col];
   const double Vvalimhilohi = Vimhilohi[col];
   const double Vvalimlolohi = Vimlolohi[col];
   const double Vvalimhihilo = Vimhihilo[col];
   const double Vvalimlohilo = Vimlohilo[col];
   const double Vvalimhilolo = Vimhilolo[col];
   const double Vvalimlololo = Vimlololo[col];
   __syncthreads();
   const double Wvalrehihihi = Wrehihihi[row];
   const double Wvalrelohihi = Wrelohihi[row];
   const double Wvalrehilohi = Wrehilohi[row];
   const double Wvalrelolohi = Wrelolohi[row];
   const double Wvalrehihilo = Wrehihilo[row];
   const double Wvalrelohilo = Wrelohilo[row];
   const double Wvalrehilolo = Wrehilolo[row];
   const double Wvalrelololo = Wrelololo[row];
   const double Wvalimhihihi = Wimhihihi[row];
   const double Wvalimlohihi = Wimlohihi[row];
   const double Wvalimhilohi = Wimhilohi[row];
   const double Wvalimlolohi = Wimlolohi[row];
   const double Wvalimhihilo = Wimhihilo[row];
   const double Wvalimlohilo = Wimlohilo[row];
   const double Wvalimhilolo = Wimhilolo[row];
   const double Wvalimlololo = Wimlololo[row];
   // const double result = Vval*Wval;
   __syncthreads();
   double resulthihihi = WYHrehihihi[idx];
   double resultlohihi = WYHrelohihi[idx];
   double resulthilohi = WYHrehilohi[idx];
   double resultlolohi = WYHrelolohi[idx];
   double resulthihilo = WYHrehihilo[idx];
   double resultlohilo = WYHrelohilo[idx];
   double resulthilolo = WYHrehilolo[idx];
   double resultlololo = WYHrelololo[idx];
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   // take the Hermitian transpose of V
   __syncthreads();
   odg_mul(Vvalrehihihi,Vvalrelohihi,Vvalrehilohi,Vvalrelolohi,
           Vvalrehihilo,Vvalrelohilo,Vvalrehilolo,Vvalrelololo,
           Wvalrehihihi,Wvalrelohihi,Wvalrehilohi,Wvalrelolohi,
           Wvalrehihilo,Wvalrelohilo,Wvalrehilolo,Wvalrelololo,
             &acchihihi,  &acclohihi,  &acchilohi,  &acclolohi,
             &acchihilo,  &acclohilo,  &acchilolo,  &acclololo);
   odg_inc(&resulthihihi,&resultlohihi,&resulthilohi,&resultlolohi,
           &resulthihilo,&resultlohilo,&resulthilolo,&resultlololo,
               acchihihi,    acclohihi,    acchilohi,    acclolohi,
               acchihilo,    acclohilo,    acchilolo,    acclololo);
   odg_mul(Vvalimhihihi,Vvalimlohihi,Vvalimhilohi,Vvalimlolohi,
           Vvalimhihilo,Vvalimlohilo,Vvalimhilolo,Vvalimlololo,
           Wvalimhihihi,Wvalimlohihi,Wvalimhilohi,Wvalimlolohi,
           Wvalimhihilo,Wvalimlohilo,Wvalimhilolo,Wvalimlololo,
             &acchihihi,  &acclohihi,  &acchilohi,  &acclolohi,
             &acchihilo,  &acclohilo,  &acchilolo,  &acclololo);
   odg_inc(&resulthihihi,&resultlohihi,&resulthilohi,&resultlolohi,
           &resulthihilo,&resultlohilo,&resulthilolo,&resultlololo,
               acchihihi,    acclohihi,    acchilohi,    acclolohi,
               acchihilo,    acclohilo,    acchilolo,    acclololo);

   if(idx < dim*dim)
   {
      WYHrehihihi[idx] = resulthihihi;
      WYHrelohihi[idx] = resultlohihi;
      WYHrehilohi[idx] = resulthilohi;
      WYHrelolohi[idx] = resultlolohi;
      WYHrehihilo[idx] = resulthihilo;
      WYHrelohilo[idx] = resultlohilo;
      WYHrehilolo[idx] = resulthilolo;
      WYHrelololo[idx] = resultlololo;
   }
   __syncthreads();
   resulthihihi = WYHimhihihi[idx];
   resultlohihi = WYHimlohihi[idx];
   resulthilohi = WYHimhilohi[idx];
   resultlolohi = WYHimlolohi[idx];
   resulthihilo = WYHimhihilo[idx];
   resultlohilo = WYHimlohilo[idx];
   resulthilolo = WYHimhilolo[idx];
   resultlololo = WYHimlololo[idx];

   odg_mul(Vvalrehihihi,Vvalrelohihi,Vvalrehilohi,Vvalrelolohi,
           Vvalrehihilo,Vvalrelohilo,Vvalrehilolo,Vvalrelololo,
           Wvalimhihihi,Wvalimlohihi,Wvalimhilohi,Wvalimlolohi,
           Wvalimhihilo,Wvalimlohilo,Wvalimhilolo,Wvalimlololo,
             &acchihihi,  &acclohihi,  &acchilohi,  &acclolohi,
             &acchihilo,  &acclohilo,  &acchilolo,  &acclololo);
   odg_inc(&resulthihihi,&resultlohihi,&resulthilohi,&resultlolohi,
           &resulthihilo,&resultlohilo,&resulthilolo,&resultlololo,
               acchihihi,    acclohihi,    acchilohi,    acclolohi,
               acchihilo,    acclohilo,    acchilolo,    acclololo);
   odg_mul(Vvalimhihihi,Vvalimlohihi,Vvalimhilohi,Vvalimlolohi,
           Vvalimhihilo,Vvalimlohilo,Vvalimhilolo,Vvalimlololo,
           Wvalrehihihi,Wvalrelohihi,Wvalrehilohi,Wvalrelolohi,
           Wvalrehihilo,Wvalrelohilo,Wvalrehilolo,Wvalrelololo,
             &acchihihi,  &acclohihi,  &acchilohi,  &acclolohi,
             &acchihilo,  &acclohilo,  &acchilolo,  &acclololo);
   odg_dec(&resulthihihi,&resultlohihi,&resulthilohi,&resultlolohi,
           &resulthihilo,&resultlohilo,&resulthilolo,&resultlololo,
               acchihihi,    acclohihi,    acchilohi,    acclolohi,
               acchihilo,    acclohilo,    acchilolo,    acclololo);

   __syncthreads();
   if(idx < dim*dim)
   {
      WYHimhihihi[idx] = resulthihihi;
      WYHimlohihi[idx] = resultlohihi;
      WYHimhilohi[idx] = resulthilohi;
      WYHimlolohi[idx] = resultlolohi;
      WYHimhihilo[idx] = resulthihilo;
      WYHimlohilo[idx] = resultlohilo;
      WYHimhilolo[idx] = resulthilolo;
      WYHimlololo[idx] = resultlololo;
   }
}

__global__ void dbl8_beta_next_W
 ( int nrows, int szt,
   double *Bhihihi, double *Blohihi, double *Bhilohi, double *Blolohi,
   double *Bhihilo, double *Blohilo, double *Bhilolo, double *Blololo,
   double *Vhihihi, double *Vlohihi, double *Vhilohi, double *Vlolohi,
   double *Vhihilo, double *Vlohilo, double *Vhilolo, double *Vlololo,
   double *Whihihi, double *Wlohihi, double *Whilohi, double *Wlolohi,
   double *Whihilo, double *Wlohilo, double *Whilolo, double *Wlololo,
   double *WYThihihi, double *WYTlohihi,
   double *WYThilohi, double *WYTlolohi,
   double *WYThihilo, double *WYTlohilo,
   double *WYThilolo, double *WYTlololo )
{
   const int bdx = blockIdx.x;        // index of block
   const int tdx = threadIdx.x;       // index of thread in block
   const int idx = bdx*szt + tdx;     // global index of the thread
   const int WYToff = idx*nrows;      // start of idx row in YWT
   const double mybetahihihi = Bhihihi[0];
   const double mybetalohihi = Blohihi[0];
   const double mybetahilohi = Bhilohi[0];
   const double mybetalolohi = Blolohi[0];
   const double mybetahihilo = Bhihilo[0];
   const double mybetalohilo = Blohilo[0];
   const double mybetahilolo = Bhilolo[0];
   const double mybetalololo = Blololo[0];
   int vdx,ydx;
   double resulthihihi,resultlohihi,resulthilohi,resultlolohi;
   double resulthihilo,resultlohilo,resulthilolo,resultlololo;
   double WYTvalhihihi,WYTvallohihi,WYTvalhilohi,WYTvallolohi;
   double WYTvalhihilo,WYTvallohilo,WYTvalhilolo,WYTvallololo;
   double Vvalhihihi,Vvallohihi,Vvalhilohi,Vvallolohi;
   double Vvalhihilo,Vvallohilo,Vvalhilolo,Vvallololo;
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   __shared__ double shVhihihi[od_shmemsize];   // to store a slice of V
   __shared__ double shVlohihi[od_shmemsize];
   __shared__ double shVhilohi[od_shmemsize];
   __shared__ double shVlolohi[od_shmemsize];
   __shared__ double shVhihilo[od_shmemsize];
   __shared__ double shVlohilo[od_shmemsize];
   __shared__ double shVhilolo[od_shmemsize];
   __shared__ double shVlololo[od_shmemsize];

   shVhihihi[tdx] = Vhihihi[idx]; // thread tdx loads the data
   shVlohihi[tdx] = Vlohihi[idx]; // at the global index
   shVhilohi[tdx] = Vhilohi[idx];
   shVlolohi[tdx] = Vlolohi[idx];
   shVhihilo[tdx] = Vhihilo[idx];
   shVlohilo[tdx] = Vlohilo[idx];
   shVhilolo[tdx] = Vhilolo[idx];
   shVlololo[tdx] = Vlololo[idx];

   __syncthreads();
   resulthihihi = shVhihihi[tdx]; // thread tdx computes the value
   resultlohihi = shVlohihi[tdx]; // at the index idx
   resulthilohi = shVhilohi[tdx];
   resultlolohi = shVlolohi[tdx];
   resulthihilo = shVhihilo[tdx];
   resultlohilo = shVlohilo[tdx];
   resulthilolo = shVhilolo[tdx];
   resultlololo = shVlololo[tdx];

   for(int i=0; i<nrows/szt; i++)
   {
      vdx = i*szt + tdx;                 // index in V and in YWT
      __syncthreads();
      shVhihihi[tdx] = Vhihihi[vdx];     // threads load next szt values
      shVlohihi[tdx] = Vlohihi[vdx];
      shVhilohi[tdx] = Vhilohi[vdx];
      shVlolohi[tdx] = Vlolohi[vdx];
      shVhihilo[tdx] = Vhihilo[vdx];
      shVlohilo[tdx] = Vlohilo[vdx];
      shVhilolo[tdx] = Vhilolo[vdx];
      shVlololo[tdx] = Vlololo[vdx];

      __syncthreads();
      for(int j=0; j<szt; j++)           // multiply szt values with YWT
      {
         ydx = WYToff + i*szt + j;       // WYT is stored row by row
         __syncthreads();
         WYTvalhihihi = WYThihihi[ydx];
         WYTvallohihi = WYTlohihi[ydx];
         WYTvalhilohi = WYThilohi[ydx];
         WYTvallolohi = WYTlolohi[ydx];
         WYTvalhihilo = WYThihilo[ydx];
         WYTvallohilo = WYTlohilo[ydx];
         WYTvalhilolo = WYThilolo[ydx];
         WYTvallololo = WYTlololo[ydx];
         __syncthreads();
         Vvalhihihi = shVhihihi[j];
         Vvallohihi = shVlohihi[j];
         Vvalhilohi = shVhilohi[j];
         Vvallolohi = shVlolohi[j];
         Vvalhihilo = shVhihilo[j];
         Vvallohilo = shVlohilo[j];
         Vvalhilolo = shVhilolo[j];
         Vvallololo = shVlololo[j];
         // result = result + WYTval*Vvalue;
         __syncthreads();
         odg_mul(  Vvalhihihi,  Vvallohihi,  Vvalhilohi,  Vvallolohi,
                   Vvalhihilo,  Vvallohilo,  Vvalhilolo,  Vvallololo,
                 WYTvalhihihi,WYTvallohihi,WYTvalhilohi,WYTvallolohi,
                 WYTvalhihilo,WYTvallohilo,WYTvalhilolo,WYTvallololo,
                   &acchihihi,  &acclohihi,  &acchilohi,  &acclolohi,
                   &acchihilo,  &acclohilo,  &acchilolo,  &acclololo);
         odg_inc(&resulthihihi,&resultlohihi,&resulthilohi,&resultlolohi,
                 &resulthihilo,&resultlohilo,&resulthilolo,&resultlololo,
                     acchihihi,    acclohihi,    acchilohi,    acclolohi,
                     acchihilo,    acclohilo,    acchilolo,    acclololo);
      }
      __syncthreads();
   }
   int quot = nrows/szt;
   int rest = nrows - quot*szt;          // remainder to compute

   vdx = quot*szt + tdx;                 // next index to compute
   __syncthreads();
   shVhihihi[tdx] = Vhihihi[vdx];
   shVlohihi[tdx] = Vlohihi[vdx];
   shVhilohi[tdx] = Vhilohi[vdx];
   shVlolohi[tdx] = Vlolohi[vdx];
   shVhihilo[tdx] = Vhihilo[vdx];
   shVlohilo[tdx] = Vlohilo[vdx];
   shVhilolo[tdx] = Vhilolo[vdx];
   shVlololo[tdx] = Vlololo[vdx];

   for(int j=0; j<rest; j++)            // rest < szt prevents overflow
   {
      __syncthreads();
      ydx = WYToff + quot*szt + j;
      WYTvalhihihi = WYThihihi[ydx];
      WYTvallohihi = WYTlohihi[ydx];
      WYTvalhilohi = WYThilohi[ydx];
      WYTvallolohi = WYTlolohi[ydx];
      WYTvalhihilo = WYThihilo[ydx];
      WYTvallohilo = WYTlohilo[ydx];
      WYTvalhilolo = WYThilolo[ydx];
      WYTvallololo = WYTlololo[ydx];
      __syncthreads();
      Vvalhihihi = shVhihihi[j];
      Vvallohihi = shVlohihi[j];
      Vvalhilohi = shVhilohi[j];
      Vvallolohi = shVlolohi[j];
      Vvalhihilo = shVhihilo[j];
      Vvallohilo = shVlohilo[j];
      Vvalhilolo = shVhilolo[j];
      Vvallololo = shVlololo[j];
      // result = result + WYTval*Vvalue;
      __syncthreads();
      odg_mul(  Vvalhihihi,  Vvallohihi,  Vvalhilohi,  Vvallolohi,
                Vvalhihilo,  Vvallohilo,  Vvalhilolo,  Vvallololo,
              WYTvalhihihi,WYTvallohihi,WYTvalhilohi,WYTvallolohi,
              WYTvalhihilo,WYTvallohilo,WYTvalhilolo,WYTvallololo,
                &acchihihi,  &acclohihi,  &acchilohi,  &acclolohi,
                &acchihilo,  &acclohilo,  &acchilolo,  &acclololo);
      odg_inc(&resulthihihi,&resultlohihi,&resulthilohi,&resultlolohi,
              &resulthihilo,&resultlohilo,&resulthilolo,&resultlololo,
                  acchihihi,    acclohihi,    acchilohi,    acclolohi,
                  acchihilo,    acclohilo,    acchilolo,    acclololo);
   }
   // result = -mybeta*result;
   __syncthreads();
   odg_mul(-mybetahihihi,-mybetalohihi,-mybetahilohi,-mybetalolohi,
           -mybetahihilo,-mybetalohilo,-mybetahilolo,-mybetalololo,
            resulthihihi, resultlohihi, resulthilohi, resultlolohi,
            resulthihilo, resultlohilo, resulthilolo, resultlololo,
              &acchihihi,   &acclohihi,   &acchilohi,   &acclolohi,
              &acchihilo,   &acclohilo,   &acchilolo,   &acclololo);

   __syncthreads();
   if(idx < nrows) 
   {
      Whihihi[idx] = acchihihi;
      Wlohihi[idx] = acclohihi;
      Whilohi[idx] = acchilohi;
      Wlolohi[idx] = acclolohi;
      Whihilo[idx] = acchihilo;
      Wlohilo[idx] = acclohilo;
      Whilolo[idx] = acchilolo;
      Wlololo[idx] = acclololo;
   }
}

__global__ void cmplx8_beta_next_W
 ( int nrows, int szt,
   double *Bhihihi, double *Blohihi, double *Bhilohi, double *Blolohi,
   double *Bhihilo, double *Blohilo, double *Bhilolo, double *Blololo,
   double *Vrehihihi, double *Vrelohihi, double *Vrehilohi, double *Vrelolohi,
   double *Vrehihilo, double *Vrelohilo, double *Vrehilolo, double *Vrelololo,
   double *Vimhihihi, double *Vimlohihi, double *Vimhilohi, double *Vimlolohi,
   double *Vimhihilo, double *Vimlohilo, double *Vimhilolo, double *Vimlololo,
   double *Wrehihihi, double *Wrelohihi, double *Wrehilohi, double *Wrelolohi,
   double *Wrehihilo, double *Wrelohilo, double *Wrehilolo, double *Wrelololo,
   double *Wimhihihi, double *Wimlohihi, double *Wimhilohi, double *Wimlolohi,
   double *Wimhihilo, double *Wimlohilo, double *Wimhilolo, double *Wimlololo,
   double *WYHrehihihi, double *WYHrelohihi,
   double *WYHrehilohi, double *WYHrelolohi,
   double *WYHrehihilo, double *WYHrelohilo,
   double *WYHrehilolo, double *WYHrelololo,
   double *WYHimhihihi, double *WYHimlohihi,
   double *WYHimhilohi, double *WYHimlolohi,
   double *WYHimhihilo, double *WYHimlohilo,
   double *WYHimhilolo, double *WYHimlololo )
{
   const int bdx = blockIdx.x;        // index of block
   const int tdx = threadIdx.x;       // index of thread in block
   const int idx = bdx*szt + tdx;     // global index of the thread
   const int WYHoff = idx*nrows;      // start of idx row in YWH
   const double mybetahihihi = Bhihihi[0];
   const double mybetalohihi = Blohihi[0];
   const double mybetahilohi = Bhilohi[0];
   const double mybetalolohi = Blolohi[0];
   const double mybetahihilo = Bhihilo[0];
   const double mybetalohilo = Blohilo[0];
   const double mybetahilolo = Bhilolo[0];
   const double mybetalololo = Blololo[0];
   // __syncthreads();
   int vdx,ydx;
   double resultrehihihi,resultrelohihi,resultrehilohi,resultrelolohi;
   double resultrehihilo,resultrelohilo,resultrehilolo,resultrelololo;
   double resultimhihihi,resultimlohihi,resultimhilohi,resultimlolohi;
   double resultimhihilo,resultimlohilo,resultimhilolo,resultimlololo;
   double WYHvalhihihi,WYHvallohihi,WYHvalhilohi,WYHvallolohi;
   double WYHvalhihilo,WYHvallohilo,WYHvalhilolo,WYHvallololo;
   double Vvalhihihi,Vvallohihi,Vvalhilohi,Vvallolohi;
   double Vvalhihilo,Vvallohilo,Vvalhilolo,Vvallololo;
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   __shared__ double shVrehihihi[cod_shmemsize]; // to store a slice of V
   __shared__ double shVrelohihi[cod_shmemsize];
   __shared__ double shVrehilohi[cod_shmemsize];
   __shared__ double shVrelolohi[cod_shmemsize];
   __shared__ double shVrehihilo[cod_shmemsize];
   __shared__ double shVrelohilo[cod_shmemsize];
   __shared__ double shVrehilolo[cod_shmemsize];
   __shared__ double shVrelololo[cod_shmemsize];
   __shared__ double shVimhihihi[cod_shmemsize];
   __shared__ double shVimlohihi[cod_shmemsize];
   __shared__ double shVimhilohi[cod_shmemsize];
   __shared__ double shVimlolohi[cod_shmemsize];
   __shared__ double shVimhihilo[cod_shmemsize];
   __shared__ double shVimlohilo[cod_shmemsize];
   __shared__ double shVimhilolo[cod_shmemsize];
   __shared__ double shVimlololo[cod_shmemsize];

   shVrehihihi[tdx] = Vrehihihi[idx]; // thread tdx loads data
   shVrelohihi[tdx] = Vrelohihi[idx]; // at the global index
   shVrehilohi[tdx] = Vrehilohi[idx];
   shVrelolohi[tdx] = Vrelolohi[idx];
   shVrehihilo[tdx] = Vrehihilo[idx]; 
   shVrelohilo[tdx] = Vrelohilo[idx];
   shVrehilolo[tdx] = Vrehilolo[idx];
   shVrelololo[tdx] = Vrelololo[idx];
   shVimhihihi[tdx] = Vimhihihi[idx];
   shVimlohihi[tdx] = Vimlohihi[idx];
   shVimhilohi[tdx] = Vimhilohi[idx];
   shVimlolohi[tdx] = Vimlolohi[idx];
   shVimhihilo[tdx] = Vimhihilo[idx];
   shVimlohilo[tdx] = Vimlohilo[idx];
   shVimhilolo[tdx] = Vimhilolo[idx];
   shVimlololo[tdx] = Vimlololo[idx];

   __syncthreads();
   resultrehihihi = shVrehihihi[tdx]; // thread tdx computes the value at idx
   resultrelohihi = shVrelohihi[tdx];
   resultrehilohi = shVrehilohi[tdx];
   resultrelolohi = shVrelolohi[tdx];
   resultrehihilo = shVrehihilo[tdx];
   resultrelohilo = shVrelohilo[tdx];
   resultrehilolo = shVrehilolo[tdx];
   resultrelololo = shVrelololo[tdx];
   resultimhihihi = shVimhihihi[tdx];
   resultimlohihi = shVimlohihi[tdx];
   resultimhilohi = shVimhilohi[tdx];
   resultimlolohi = shVimlolohi[tdx];
   resultimhihilo = shVimhihilo[tdx];
   resultimlohilo = shVimlohilo[tdx];
   resultimhilolo = shVimhilolo[tdx];
   resultimlololo = shVimlololo[tdx];

   for(int i=0; i<nrows/szt; i++)
   {
      vdx = i*szt + tdx;                 // index in V and in YWT
      __syncthreads();
      shVrehihihi[tdx] = Vrehihihi[vdx]; // threads load next szt values
      shVrelohihi[tdx] = Vrelohihi[vdx];
      shVrehilohi[tdx] = Vrehilohi[vdx];
      shVrelolohi[tdx] = Vrelolohi[vdx];
      shVrehihilo[tdx] = Vrehihilo[vdx];
      shVrelohilo[tdx] = Vrelohilo[vdx];
      shVrehilolo[tdx] = Vrehilolo[vdx];
      shVrelololo[tdx] = Vrelololo[vdx];
      shVimhihihi[tdx] = Vimhihihi[vdx];
      shVimlohihi[tdx] = Vimlohihi[vdx];
      shVimhilohi[tdx] = Vimhilohi[vdx];
      shVimlolohi[tdx] = Vimlolohi[vdx];
      shVimhihilo[tdx] = Vimhihilo[vdx];
      shVimlohilo[tdx] = Vimlohilo[vdx];
      shVimhilolo[tdx] = Vimhilolo[vdx];
      shVimlololo[tdx] = Vimlololo[vdx];

      __syncthreads();
      for(int j=0; j<szt; j++)           // multiply szt values with YWT
      {
         ydx = WYHoff + i*szt + j;       // WYH is stored row by row
         // take the Hermitian transpose of V

         Vvalhihihi = shVrehihihi[j];
         Vvallohihi = shVrelohihi[j];
         Vvalhilohi = shVrehilohi[j];
         Vvallolohi = shVrelolohi[j];
         Vvalhihilo = shVrehihilo[j];
         Vvallohilo = shVrelohilo[j];
         Vvalhilolo = shVrehilolo[j];
         Vvallololo = shVrelololo[j];
         WYHvalhihihi = WYHrehihihi[ydx];
         WYHvallohihi = WYHrelohihi[ydx];
         WYHvalhilohi = WYHrehilohi[ydx];
         WYHvallolohi = WYHrelolohi[ydx];
         WYHvalhihilo = WYHrehihilo[ydx];
         WYHvallohilo = WYHrelohilo[ydx];
         WYHvalhilolo = WYHrehilolo[ydx];
         WYHvallololo = WYHrelololo[ydx];
         odg_mul(  Vvalhihihi,  Vvallohihi,  Vvalhilohi,  Vvallolohi,
                   Vvalhihilo,  Vvallohilo,  Vvalhilolo,  Vvallololo,
                 WYHvalhihihi,WYHvallohihi,WYHvalhilohi,WYHvallolohi,
                 WYHvalhihilo,WYHvallohilo,WYHvalhilolo,WYHvallololo,
                   &acchihihi,  &acclohihi,  &acchilohi,  &acclolohi,
                   &acchihilo,  &acclohilo,  &acchilolo,  &acclololo);
         odg_inc(&resultrehihihi,&resultrelohihi,
                 &resultrehilohi,&resultrelolohi,
                 &resultrehihilo,&resultrelohilo,
                 &resultrehilolo,&resultrelololo,
                       acchihihi,      acclohihi,  acchilohi,  acclolohi,
                       acchihilo,      acclohilo,  acchilolo,  acclololo);

         Vvalhihihi = shVimhihihi[j];
         Vvallohihi = shVimlohihi[j];
         Vvalhilohi = shVimhilohi[j];
         Vvallolohi = shVimlolohi[j];
         Vvalhihilo = shVimhihilo[j];
         Vvallohilo = shVimlohilo[j];
         Vvalhilolo = shVimhilolo[j];
         Vvallololo = shVimlololo[j];
         WYHvalhihihi = WYHimhihihi[ydx];
         WYHvallohihi = WYHimlohihi[ydx];
         WYHvalhilohi = WYHimhilohi[ydx];
         WYHvallolohi = WYHimlolohi[ydx];
         WYHvalhihilo = WYHimhihilo[ydx];
         WYHvallohilo = WYHimlohilo[ydx];
         WYHvalhilolo = WYHimhilolo[ydx];
         WYHvallololo = WYHimlololo[ydx];
         odg_mul(  Vvalhihihi,  Vvallohihi,  Vvalhilohi,  Vvallolohi,
                   Vvalhihilo,  Vvallohilo,  Vvalhilolo,  Vvallololo,
                 WYHvalhihihi,WYHvallohihi,WYHvalhilohi,WYHvallolohi,
                 WYHvalhihilo,WYHvallohilo,WYHvalhilolo,WYHvallololo,
                   &acchihihi,  &acclohihi,  &acchilohi,  &acclolohi,
                   &acchihilo,  &acclohilo,  &acchilolo,  &acclololo);
         odg_dec(&resultrehihihi,&resultrelohihi,
                 &resultrehilohi,&resultrelolohi,
                 &resultrehihilo,&resultrelohilo,
                 &resultrehilolo,&resultrelololo,
                       acchihihi,      acclohihi,  acchilohi,  acclolohi,
                       acchihilo,      acclohilo,  acchilolo,  acclololo);

         // have already imaginary values for V
         WYHvalhihihi = WYHrehihihi[ydx];
         WYHvallohihi = WYHrelohihi[ydx];
         WYHvalhilohi = WYHrehilohi[ydx];
         WYHvallolohi = WYHrelolohi[ydx];
         WYHvalhihilo = WYHrehihilo[ydx];
         WYHvallohilo = WYHrelohilo[ydx];
         WYHvalhilolo = WYHrehilolo[ydx];
         WYHvallololo = WYHrelololo[ydx];
         odg_mul(  Vvalhihihi,  Vvallohihi,  Vvalhilohi,  Vvallolohi,
                   Vvalhihilo,  Vvallohilo,  Vvalhilolo,  Vvallololo,
                 WYHvalhihihi,WYHvallohihi,WYHvalhilohi,WYHvallolohi,
                 WYHvalhihilo,WYHvallohilo,WYHvalhilolo,WYHvallololo,
                   &acchihihi,  &acclohihi,  &acchilohi,  &acclolohi,
                   &acchihilo,  &acclohilo,  &acchilolo,  &acclololo);
         odg_inc(&resultimhihihi,&resultimlohihi,
                 &resultimhilohi,&resultimlolohi,
                 &resultimhihilo,&resultimlohilo,
                 &resultimhilolo,&resultimlololo,
                       acchihihi,      acclohihi,  acchilohi,  acclolohi,
                       acchihilo,      acclohilo,  acchilolo,  acclololo);

         Vvalhihihi = shVrehihihi[j];
         Vvallohihi = shVrelohihi[j];
         Vvalhilohi = shVrehilohi[j];
         Vvallolohi = shVrelolohi[j];
         Vvalhihilo = shVrehihilo[j];
         Vvallohilo = shVrelohilo[j];
         Vvalhilolo = shVrehilolo[j];
         Vvallololo = shVrelololo[j];
         WYHvalhihihi = WYHimhihihi[ydx];
         WYHvallohihi = WYHimlohihi[ydx];
         WYHvalhilohi = WYHimhilohi[ydx];
         WYHvallolohi = WYHimlolohi[ydx];
         WYHvalhihilo = WYHimhihilo[ydx];
         WYHvallohilo = WYHimlohilo[ydx];
         WYHvalhilolo = WYHimhilolo[ydx];
         WYHvallololo = WYHimlololo[ydx];
         odg_mul(  Vvalhihihi,  Vvallohihi,  Vvalhilohi,  Vvallolohi,
                   Vvalhihilo,  Vvallohilo,  Vvalhilolo,  Vvallololo,
                 WYHvalhihihi,WYHvallohihi,WYHvalhilohi,WYHvallolohi,
                 WYHvalhihilo,WYHvallohilo,WYHvalhilolo,WYHvallololo,
                   &acchihihi,  &acclohihi,  &acchilohi,  &acclolohi,
                   &acchihilo,  &acclohilo,  &acchilolo,  &acclololo);
         odg_inc(&resultimhihihi,&resultimlohihi,
                 &resultimhilohi,&resultimlolohi,
                 &resultimhihilo,&resultimlohilo,
                 &resultimhilolo,&resultimlololo,
                       acchihihi,      acclohihi,  acchilohi,  acclolohi,
                       acchihilo,      acclohilo,  acchilolo,  acclololo);
      }
      __syncthreads();
   }
   int quot = nrows/szt;
   int rest = nrows - quot*szt;          // remainder to compute

   vdx = quot*szt + tdx;                 // next index to compute
   __syncthreads();
   shVrehihihi[tdx] = Vrehihihi[vdx];
   shVrelohihi[tdx] = Vrelohihi[vdx];
   shVrehilohi[tdx] = Vrehilohi[vdx];
   shVrelolohi[tdx] = Vrelolohi[vdx];
   shVrehihilo[tdx] = Vrehihilo[vdx];
   shVrelohilo[tdx] = Vrelohilo[vdx];
   shVrehilolo[tdx] = Vrehilolo[vdx];
   shVrelololo[tdx] = Vrelololo[vdx];
   shVimhihihi[tdx] = Vimhihihi[vdx];
   shVimlohihi[tdx] = Vimlohihi[vdx];
   shVimhilohi[tdx] = Vimhilohi[vdx];
   shVimlolohi[tdx] = Vimlolohi[vdx];
   shVimhihilo[tdx] = Vimhihilo[vdx];
   shVimlohilo[tdx] = Vimlohilo[vdx];
   shVimhilolo[tdx] = Vimhilolo[vdx];
   shVimlololo[tdx] = Vimlololo[vdx];

   for(int j=0; j<rest; j++)            // rest < szt prevents overflow
   {
      // result = result + WYTval*Vvalue;
      // take the Hermitian transpose of V

      Vvalhihihi = shVrehihihi[j];
      Vvallohihi = shVrelohihi[j];
      Vvalhilohi = shVrehilohi[j];
      Vvallolohi = shVrelolohi[j];
      Vvalhihilo = shVrehihilo[j];
      Vvallohilo = shVrelohilo[j];
      Vvalhilolo = shVrehilolo[j];
      Vvallololo = shVrelololo[j];
      WYHvalhihihi = WYHrehihihi[ydx];
      WYHvallohihi = WYHrelohihi[ydx];
      WYHvalhilohi = WYHrehilohi[ydx];
      WYHvallolohi = WYHrelolohi[ydx];
      WYHvalhihilo = WYHrehihilo[ydx];
      WYHvallohilo = WYHrelohilo[ydx];
      WYHvalhilolo = WYHrehilolo[ydx];
      WYHvallololo = WYHrelololo[ydx];
      odg_mul(  Vvalhihihi,  Vvallohihi,  Vvalhilohi,  Vvallolohi,
                Vvalhihilo,  Vvallohilo,  Vvalhilolo,  Vvallololo,
              WYHvalhihihi,WYHvallohihi,WYHvalhilohi,WYHvallolohi,
              WYHvalhihilo,WYHvallohilo,WYHvalhilolo,WYHvallololo,
                &acchihihi,  &acclohihi,  &acchilohi,  &acclolohi,
                &acchihilo,  &acclohilo,  &acchilolo,  &acclololo);
      odg_inc(&resultrehihihi,&resultrelohihi,&resultrehilohi,&resultrelolohi,
              &resultrehihilo,&resultrelohilo,&resultrehilolo,&resultrelololo,
                    acchihihi,      acclohihi,      acchilohi,      acclolohi,
                    acchihilo,      acclohilo,      acchilolo,      acclololo);

      Vvalhihihi = shVimhihihi[j];
      Vvallohihi = shVimlohihi[j];
      Vvalhilohi = shVimhilohi[j];
      Vvallolohi = shVimlolohi[j];
      Vvalhihilo = shVimhihilo[j];
      Vvallohilo = shVimlohilo[j];
      Vvalhilolo = shVimhilolo[j];
      Vvallololo = shVimlololo[j];
      WYHvalhihihi = WYHimhihihi[ydx];
      WYHvallohihi = WYHimlohihi[ydx];
      WYHvalhilohi = WYHimhilohi[ydx];
      WYHvallolohi = WYHimlolohi[ydx];
      WYHvalhihilo = WYHimhihilo[ydx];
      WYHvallohilo = WYHimlohilo[ydx];
      WYHvalhilolo = WYHimhilolo[ydx];
      WYHvallololo = WYHimlololo[ydx];
      odg_mul(  Vvalhihihi,  Vvallohihi,  Vvalhilohi,  Vvallolohi,
                Vvalhihilo,  Vvallohilo,  Vvalhilolo,  Vvallololo,
              WYHvalhihihi,WYHvallohihi,WYHvalhilohi,WYHvallolohi,
              WYHvalhihilo,WYHvallohilo,WYHvalhilolo,WYHvallololo,
                &acchihihi,  &acclohihi,  &acchilohi,  &acclolohi,
                &acchihilo,  &acclohilo,  &acchilolo,  &acclololo);
      odg_dec(&resultrehihihi,&resultrelohihi,&resultrehilohi,&resultrelolohi,
              &resultrehihilo,&resultrelohilo,&resultrehilolo,&resultrelololo,
                    acchihihi,      acclohihi,      acchilohi,      acclolohi,
                    acchihilo,      acclohilo,      acchilolo,      acclololo);

      // already have the imaginary values for V
      WYHvalhihihi = WYHrehihihi[ydx];
      WYHvallohihi = WYHrelohihi[ydx];
      WYHvalhilohi = WYHrehilohi[ydx];
      WYHvallolohi = WYHrelolohi[ydx];
      WYHvalhihilo = WYHrehihilo[ydx];
      WYHvallohilo = WYHrelohilo[ydx];
      WYHvalhilolo = WYHrehilolo[ydx];
      WYHvallololo = WYHrelololo[ydx];
      odg_mul(  Vvalhihihi,  Vvallohihi,  Vvalhilohi,  Vvallolohi,
                Vvalhihilo,  Vvallohilo,  Vvalhilolo,  Vvallololo,
              WYHvalhihihi,WYHvallohihi,WYHvalhilohi,WYHvallolohi,
              WYHvalhihilo,WYHvallohilo,WYHvalhilolo,WYHvallololo,
                &acchihihi,  &acclohihi,  &acchilohi,  &acclolohi,
                &acchihilo,  &acclohilo,  &acchilolo,  &acclololo);
      odg_inc(&resultimhihihi,&resultimlohihi,&resultimhilohi,&resultimlolohi,
              &resultimhihilo,&resultimlohilo,&resultimhilolo,&resultimlololo,
                    acchihihi,      acclohihi,      acchilohi,      acclolohi,
                    acchihilo,      acclohilo,      acchilolo,      acclololo);

      Vvalhihihi = shVrehihihi[j];
      Vvallohihi = shVrelohihi[j];
      Vvalhilohi = shVrehilohi[j];
      Vvallolohi = shVrelolohi[j];
      Vvalhihilo = shVrehihilo[j];
      Vvallohilo = shVrelohilo[j];
      Vvalhilolo = shVrehilolo[j];
      Vvallololo = shVrelololo[j];
      WYHvalhihihi = WYHimhihihi[ydx];
      WYHvallohihi = WYHimlohihi[ydx];
      WYHvalhilohi = WYHimhilohi[ydx];
      WYHvallolohi = WYHimlolohi[ydx];
      WYHvalhihilo = WYHimhihilo[ydx];
      WYHvallohilo = WYHimlohilo[ydx];
      WYHvalhilolo = WYHimhilolo[ydx];
      WYHvallololo = WYHimlololo[ydx];
      odg_mul(  Vvalhihihi,  Vvallohihi,  Vvalhilohi,  Vvallolohi,
                Vvalhihilo,  Vvallohilo,  Vvalhilolo,  Vvallololo,
              WYHvalhihihi,WYHvallohihi,WYHvalhilohi,WYHvallolohi,
              WYHvalhihilo,WYHvallohilo,WYHvalhilolo,WYHvallololo,
                &acchihihi,  &acclohihi,  &acchilohi,  &acclolohi,
                &acchihilo,  &acclohilo,  &acchilolo,  &acclololo);
      odg_inc(&resultimhihihi,&resultimlohihi,&resultimhilohi,&resultimlolohi,
              &resultimhihilo,&resultimlohilo,&resultimhilolo,&resultimlololo,
                    acchihihi,      acclohihi,      acchilohi,      acclolohi,
                    acchihilo,      acclohilo,      acchilolo,      acclololo);
   }
   // result = -mybeta*result;
   __syncthreads();
   odg_mul( -mybetahihihi, -mybetalohihi, -mybetahilohi, -mybetalolohi,
            -mybetahihilo, -mybetalohilo, -mybetahilolo, -mybetalololo,
           resultrehihihi,resultrelohihi,resultrehilohi,resultrelolohi,
           resultrehihilo,resultrelohilo,resultrehilolo,resultrelololo,
               &acchihihi,    &acclohihi,    &acchilohi,    &acclolohi,
               &acchihilo,    &acclohilo,    &acchilolo,    &acclololo);

   // __syncthreads();
   if(idx < nrows) 
   {
      Wrehihihi[idx] = acchihihi;
      Wrelohihi[idx] = acclohihi;
      Wrehilohi[idx] = acchilohi;
      Wrelolohi[idx] = acclolohi;
      Wrehihilo[idx] = acchihilo;
      Wrelohilo[idx] = acclohilo;
      Wrehilolo[idx] = acchilolo;
      Wrelololo[idx] = acclololo;
   }
   __syncthreads();
   odg_mul( -mybetahihihi, -mybetalohihi, -mybetahilohi, -mybetalolohi,
            -mybetahihilo, -mybetalohilo, -mybetahilolo, -mybetalololo,
           resultimhihihi,resultimlohihi,resultimhilohi,resultimlolohi,
           resultimhihilo,resultimlohilo,resultimhilolo,resultimlololo,
               &acchihihi,    &acclohihi,    &acchilohi,    &acclolohi,
               &acchihilo,    &acclohilo,    &acchilolo,    &acclololo);

   __syncthreads();
   if(idx < nrows) 
   {
      Wimhihihi[idx] = acchihihi;
      Wimlohihi[idx] = acclohihi;
      Wimhilohi[idx] = acchilohi;
      Wimlolohi[idx] = acclolohi;
      Wimhihilo[idx] = acchihilo;
      Wimlohilo[idx] = acclohilo;
      Wimhilolo[idx] = acchilolo;
      Wimlololo[idx] = acclololo;
   }
}

__global__ void dbl8_small_WYT
 ( int nrows, int szt,
   double *Whihihi, double *Wlohihi, double *Whilohi, double *Wlolohi,
   double *Whihilo, double *Wlohilo, double *Whilolo, double *Wlololo,
   double *Vhihihi, double *Vlohihi, double *Vhilohi, double *Vlolohi,
   double *Vhihilo, double *Vlohilo, double *Vhilolo, double *Vlololo,
   double *WYThihihi, double *WYTlohihi,
   double *WYThilohi, double *WYTlolohi,
   double *WYThihilo, double *WYTlohilo,
   double *WYThilolo, double *WYTlololo )
{
   const int bdx = blockIdx.x;           // index of block
   const int tdx = threadIdx.x;          // index of thread in block
   const int offset = bdx*szt + tdx;     // offset in result
   const int row = offset / nrows;
   const int col = offset % nrows;       // thread 0 computes WYT[row][col]

   double resulthihihi = 0.0;
   double resultlohihi = 0.0;
   double resulthilohi = 0.0;
   double resultlolohi = 0.0;
   double resulthihilo = 0.0;
   double resultlohilo = 0.0;
   double resulthilolo = 0.0;
   double resultlololo = 0.0;
   double ahihihi,alohihi,ahilohi,alolohi;
   double ahihilo,alohilo,ahilolo,alololo;
   double bhihihi,blohihi,bhilohi,blolohi;
   double bhihilo,blohilo,bhilolo,blololo;
   double chihihi,clohihi,chilohi,clolohi;
   double chihilo,clohilo,chilolo,clololo;

   for(int k=0; k<szt; k++)
   {
      __syncthreads();
      ahihihi = Whihihi[k*nrows + row];   // if(nrows == szt) then row = bdx
      alohihi = Wlohihi[k*nrows + row];
      ahilohi = Whilohi[k*nrows + row];
      alolohi = Wlolohi[k*nrows + row];
      ahihilo = Whihilo[k*nrows + row];
      alohilo = Wlohilo[k*nrows + row];
      ahilolo = Whilolo[k*nrows + row];
      alololo = Wlololo[k*nrows + row];
      __syncthreads();
      bhihihi = Vhihihi[k*nrows + col];   // if(nrows == szt) then col = tdx
      blohihi = Vlohihi[k*nrows + col]; 
      bhilohi = Vhilohi[k*nrows + col]; 
      blolohi = Vlolohi[k*nrows + col]; 
      bhihilo = Vhihilo[k*nrows + col]; 
      blohilo = Vlohilo[k*nrows + col]; 
      bhilolo = Vhilolo[k*nrows + col]; 
      blololo = Vlololo[k*nrows + col]; 
      // result = result + a*b;
      __syncthreads();
      odg_mul( ahihihi, alohihi, ahilohi, alolohi,
               ahihilo, alohilo, ahilolo, alololo,
               bhihihi, blohihi, bhilohi, blolohi,
               bhihilo, blohilo, bhilolo, blololo,
              &chihihi,&clohihi,&chilohi,&clolohi,
              &chihilo,&clohilo,&chilolo,&clololo);
      odg_inc(&resulthihihi,&resultlohihi,&resulthilohi,&resultlolohi,
              &resulthihilo,&resultlohilo,&resulthilolo,&resultlololo,
                    chihihi,      clohihi,      chilohi,      clolohi,
                    chihilo,      clohilo,      chilolo,      clololo);
   }
   __syncthreads();
   WYThihihi[offset] = resulthihihi;
   WYTlohihi[offset] = resultlohihi;
   WYThilohi[offset] = resulthilohi;
   WYTlolohi[offset] = resultlolohi;
   WYThihilo[offset] = resulthihilo;
   WYTlohilo[offset] = resultlohilo;
   WYThilolo[offset] = resulthilolo;
   WYTlololo[offset] = resultlololo;
}

__global__ void cmplx8_small_WYH
 ( int nrows, int szt,
   double *Wrehihihi, double *Wrelohihi, double *Wrehilohi, double *Wrelolohi,
   double *Wrehihilo, double *Wrelohilo, double *Wrehilolo, double *Wrelololo,
   double *Wimhihihi, double *Wimlohihi, double *Wimhilohi, double *Wimlolohi,
   double *Wimhihilo, double *Wimlohilo, double *Wimhilolo, double *Wimlololo,
   double *Yrehihihi, double *Yrelohihi, double *Yrehilohi, double *Yrelolohi,
   double *Yrehihilo, double *Yrelohilo, double *Yrehilolo, double *Yrelololo,
   double *Yimhihihi, double *Yimlohihi, double *Yimhilohi, double *Yimlolohi,
   double *Yimhihilo, double *Yimlohilo, double *Yimhilolo, double *Yimlololo,
   double *WYTrehihihi, double *WYTrelohihi,
   double *WYTrehilohi, double *WYTrelolohi,
   double *WYTrehihilo, double *WYTrelohilo,
   double *WYTrehilolo, double *WYTrelololo,
   double *WYTimhihihi, double *WYTimlohihi,
   double *WYTimhilohi, double *WYTimlolohi,
   double *WYTimhihilo, double *WYTimlohilo,
   double *WYTimhilolo, double *WYTimlololo )
{
   const int bdx = blockIdx.x;           // index of block
   const int tdx = threadIdx.x;          // index of thread in block
   const int offset = bdx*szt + tdx;     // offset in result
   const int row = offset / nrows;
   const int col = offset % nrows;       // thread 0 computes WYT[row][col]

   double resultrehihihi = 0.0;
   double resultrelohihi = 0.0;
   double resultrehilohi = 0.0;
   double resultrelolohi = 0.0;
   double resultrehihilo = 0.0;
   double resultrelohilo = 0.0;
   double resultrehilolo = 0.0;
   double resultrelololo = 0.0;
   double resultimhihihi = 0.0;
   double resultimlohihi = 0.0;
   double resultimhilohi = 0.0;
   double resultimlolohi = 0.0;
   double resultimhihilo = 0.0;
   double resultimlohilo = 0.0;
   double resultimhilolo = 0.0;
   double resultimlololo = 0.0;
   double a_rehihihi,a_relohihi,a_rehilohi,a_relolohi;
   double a_rehihilo,a_relohilo,a_rehilolo,a_relololo;
   double a_imhihihi,a_imlohihi,a_imhilohi,a_imlolohi;
   double a_imhihilo,a_imlohilo,a_imhilolo,a_imlololo;
   double b_rehihihi,b_relohihi,b_rehilohi,b_relolohi;
   double b_rehihilo,b_relohilo,b_rehilolo,b_relololo;
   double b_imhihihi,b_imlohihi,b_imhilohi,b_imlolohi;
   double b_imhihilo,b_imlohilo,b_imhilolo,b_imlololo;
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;
   int Widx,Yidx;

   for(int k=0; k<szt; k++)
   {
      Widx = k*nrows + row;
      __syncthreads();
      a_rehihihi = Wrehihihi[Widx];    // if(nrows == szt) then row = bdx
      a_relohihi = Wrelohihi[Widx];
      a_rehilohi = Wrehilohi[Widx];
      a_relolohi = Wrelolohi[Widx];
      a_rehihilo = Wrehihilo[Widx];
      a_relohilo = Wrelohilo[Widx];
      a_rehilolo = Wrehilolo[Widx];
      a_relololo = Wrelololo[Widx];
      a_imhihihi = Wimhihihi[Widx]; 
      a_imlohihi = Wimlohihi[Widx]; 
      a_imhilohi = Wimhilohi[Widx]; 
      a_imlolohi = Wimlolohi[Widx]; 
      a_imhihilo = Wimhihilo[Widx]; 
      a_imlohilo = Wimlohilo[Widx]; 
      a_imhilolo = Wimhilolo[Widx]; 
      a_imlololo = Wimlololo[Widx]; 
      Yidx = k*nrows + col;
      __syncthreads();
      b_rehihihi = Yrehihihi[Yidx];   // if(nrows == szt) then col = tdx
      b_relohihi = Yrelohihi[Yidx];
      b_rehilohi = Yrehilohi[Yidx];
      b_relolohi = Yrelolohi[Yidx];
      b_rehihilo = Yrehihilo[Yidx];
      b_relohilo = Yrelohilo[Yidx];
      b_rehilolo = Yrehilolo[Yidx];
      b_relololo = Yrelololo[Yidx];
      b_imhihihi = Yimhihihi[Yidx];
      b_imlohihi = Yimlohihi[Yidx];
      b_imhilohi = Yimhilohi[Yidx];
      b_imlolohi = Yimlolohi[Yidx];
      b_imhihilo = Yimhihilo[Yidx];
      b_imlohilo = Yimlohilo[Yidx];
      b_imhilolo = Yimhilolo[Yidx];
      b_imlololo = Yimlololo[Yidx];
      // result = result + a*b; with Hermitian transpose of Y
      // resultre = resultre + a_re*b_re + a_im*b_im;
      __syncthreads();
      odg_mul(a_rehihihi,a_relohihi,a_rehilohi,a_relolohi,
              a_rehihilo,a_relohilo,a_rehilolo,a_relololo,
              b_rehihihi,b_relohihi,b_rehilohi,b_relolohi,
              b_rehihilo,b_relohilo,b_rehilolo,b_relololo,
              &acchihihi,&acclohihi,&acchilohi,&acclolohi,
              &acchihilo,&acclohilo,&acchilolo,&acclololo);
      odg_inc(&resultrehihihi,&resultrelohihi,&resultrehilohi,&resultrelolohi,
              &resultrehihilo,&resultrelohilo,&resultrehilolo,&resultrelololo,
                    acchihihi,      acclohihi,      acchilohi,      acclolohi,
                    acchihilo,      acclohilo,      acchilolo,      acclololo);
      odg_mul(a_imhihihi,a_imlohihi,a_imhilohi,a_imlolohi,
              a_imhihilo,a_imlohilo,a_imhilolo,a_imlololo,
              b_imhihihi,b_imlohihi,b_imhilohi,b_imlolohi,
              b_imhihilo,b_imlohilo,b_imhilolo,b_imlololo,
              &acchihihi,&acclohihi,&acchilohi,&acclolohi,
              &acchihilo,&acclohilo,&acchilolo,&acclololo);
      odg_inc(&resultrehihihi,&resultrelohihi,&resultrehilohi,&resultrelolohi,
              &resultrehihilo,&resultrelohilo,&resultrehilolo,&resultrelololo,
                    acchihihi,      acclohihi,      acchilohi,      acclolohi,
                    acchihilo,      acclohilo,      acchilolo,      acclololo);
      // resultim = resultim + a_im*b_re - a_re*b_im;
      odg_mul(a_imhihihi,a_imlohihi,a_imhilohi,a_imlolohi,
              a_imhihilo,a_imlohilo,a_imhilolo,a_imlololo,
              b_rehihihi,b_relohihi,b_rehilohi,b_relolohi,
              b_rehihilo,b_relohilo,b_rehilolo,b_relololo,
              &acchihihi,&acclohihi,&acchilohi,&acclolohi,
              &acchihilo,&acclohilo,&acchilolo,&acclololo);
      odg_inc(&resultimhihihi,&resultimlohihi,&resultimhilohi,&resultimlolohi,
              &resultimhihilo,&resultimlohilo,&resultimhilolo,&resultimlololo,
                    acchihihi,      acclohihi,      acchilohi,      acclolohi,
                    acchihilo,      acclohilo,      acchilolo,      acclololo);
      odg_mul(a_rehihihi,a_relohihi,a_rehilohi,a_relolohi,
              a_rehihilo,a_relohilo,a_rehilolo,a_relololo,
              b_imhihihi,b_imlohihi,b_imhilohi,b_imlolohi,
              b_imhihilo,b_imlohilo,b_imhilolo,b_imlololo,
              &acchihihi,&acclohihi,&acchilohi,&acclolohi,
              &acchihilo,&acclohilo,&acchilolo,&acclololo);
      odg_dec(&resultimhihihi,&resultimlohihi,&resultimhilohi,&resultimlolohi,
              &resultimhihilo,&resultimlohilo,&resultimhilolo,&resultimlololo,
                    acchihihi,      acclohihi,      acchilohi,      acclolohi,
                    acchihilo,      acclohilo,      acchilolo,      acclololo);
   }
   __syncthreads();
   WYTrehihihi[offset] = resultrehihihi;
   WYTrelohihi[offset] = resultrelohihi;
   WYTrehilohi[offset] = resultrehilohi;
   WYTrelolohi[offset] = resultrelolohi;
   WYTrehihilo[offset] = resultrehihilo;
   WYTrelohilo[offset] = resultrelohilo;
   WYTrehilolo[offset] = resultrehilolo;
   WYTrelololo[offset] = resultrelololo;
   WYTimhihihi[offset] = resultimhihihi;
   WYTimlohihi[offset] = resultimlohihi;
   WYTimhilohi[offset] = resultimhilohi;
   WYTimlolohi[offset] = resultimlolohi;
   WYTimhihilo[offset] = resultimhihilo;
   WYTimlohilo[offset] = resultimlohilo;
   WYTimhilolo[offset] = resultimhilolo;
   WYTimlololo[offset] = resultimlololo;
}

__global__ void dbl8_small_QWYT
 ( int dim, int rowdim, int szt, int coloff,
   double *Qhihihi, double *Qlohihi, double *Qhilohi, double *Qlolohi,
   double *Qhihilo, double *Qlohilo, double *Qhilolo, double *Qlololo,
   double *WYThihihi, double *WYTlohihi, double *WYThilohi, double *WYTlolohi,
   double *WYThihilo, double *WYTlohilo, double *WYThilolo, double *WYTlololo,
   double *QWYThihihi, double *QWYTlohihi,
   double *QWYThilohi, double *QWYTlolohi,
   double *QWYThihilo, double *QWYTlohilo,
   double *QWYThilolo, double *QWYTlololo )
{
   const int bdx = blockIdx.x;         // index of block
   const int tdx = threadIdx.x;        // index of thread in block
   const int offset = bdx*szt + tdx;   // offset in result
   const int row = offset / rowdim;
   const int col = offset % rowdim;    // thread 0 computes QWYT[row][col]

   double resulthihihi = 0.0;
   double resultlohihi = 0.0;
   double resulthilohi = 0.0;
   double resultlolohi = 0.0;
   double resulthihilo = 0.0;
   double resultlohilo = 0.0;
   double resulthilolo = 0.0;
   double resultlololo = 0.0;
   double ahihihi,alohihi,ahilohi,alolohi;
   double ahihilo,alohilo,ahilolo,alololo;
   double bhihihi,blohihi,bhilohi,blolohi;
   double bhihilo,blohilo,bhilolo,blololo;
   double chihihi,clohihi,chilohi,clolohi;
   double chihilo,clohilo,chilolo,clololo;
   int idx;

   for(int k=0; k<rowdim; k++)       // run over rowdim, not just szt
   {                                 // coloff shifts by col*row elements
      idx = row*dim + coloff + k;
      __syncthreads();
      ahihihi = Qhihihi[idx];        // row = bdx,
      alohihi = Qlohihi[idx];
      ahilohi = Qhilohi[idx];        // if dim == szt, coloff == 0
      alolohi = Qlolohi[idx];
      ahihilo = Qhihilo[idx];
      alohilo = Qlohilo[idx];
      ahilolo = Qhilolo[idx];
      alololo = Qlololo[idx];
      idx = k*rowdim + col;
      __syncthreads();
      bhihihi = WYThihihi[idx];      // if(dim == szt) then col = tdx
      blohihi = WYTlohihi[idx];
      bhilohi = WYThilohi[idx];
      blolohi = WYTlolohi[idx]; 
      bhihilo = WYThihilo[idx];
      blohilo = WYTlohilo[idx];
      bhilolo = WYThilolo[idx];
      blololo = WYTlololo[idx]; 
      // result = result + a*b;
      __syncthreads();
      odg_mul( ahihihi, alohihi, ahilohi, alolohi,
               ahihilo, alohilo, ahilolo, alololo,
               bhihihi, blohihi, bhilohi, blolohi,
               bhihilo, blohilo, bhilolo, blololo,
              &chihihi,&clohihi,&chilohi,&clolohi,
              &chihilo,&clohilo,&chilolo,&clololo);
      odg_inc(&resulthihihi,&resultlohihi,&resulthilohi,&resultlolohi,
              &resulthihilo,&resultlohilo,&resulthilolo,&resultlololo,
                    chihihi,      clohihi,      chilohi,      clolohi,
                    chihilo,      clohilo,      chilolo,      clololo);
   }
   __syncthreads();
   QWYThihihi[offset] = resulthihihi;  // no column offset in saving QWYT
   QWYTlohihi[offset] = resultlohihi;
   QWYThilohi[offset] = resulthilohi;
   QWYTlolohi[offset] = resultlolohi;
   QWYThihilo[offset] = resulthihilo;
   QWYTlohilo[offset] = resultlohilo;
   QWYThilolo[offset] = resulthilolo;
   QWYTlololo[offset] = resultlololo;
}

__global__ void cmplx8_small_QWYH
 ( int dim, int rowdim, int szt, int coloff,
   double *Qrehihihi, double *Qrelohihi, double *Qrehilohi, double *Qrelolohi,
   double *Qrehihilo, double *Qrelohilo, double *Qrehilolo, double *Qrelololo,
   double *Qimhihihi, double *Qimlohihi, double *Qimhilohi, double *Qimlolohi,
   double *Qimhihilo, double *Qimlohilo, double *Qimhilolo, double *Qimlololo,
   double *WYTrehihihi, double *WYTrelohihi,
   double *WYTrehilohi, double *WYTrelolohi,
   double *WYTrehihilo, double *WYTrelohilo,
   double *WYTrehilolo, double *WYTrelololo,
   double *WYTimhihihi, double *WYTimlohihi,
   double *WYTimhilohi, double *WYTimlolohi,
   double *WYTimhihilo, double *WYTimlohilo,
   double *WYTimhilolo, double *WYTimlololo,
   double *QWYTrehihihi, double *QWYTrelohihi,
   double *QWYTrehilohi, double *QWYTrelolohi,
   double *QWYTrehihilo, double *QWYTrelohilo,
   double *QWYTrehilolo, double *QWYTrelololo,
   double *QWYTimhihihi, double *QWYTimlohihi,
   double *QWYTimhilohi, double *QWYTimlolohi,
   double *QWYTimhihilo, double *QWYTimlohilo,
   double *QWYTimhilolo, double *QWYTimlololo )
{
   const int bdx = blockIdx.x;         // index of block
   const int tdx = threadIdx.x;        // index of thread in block
   const int offset = bdx*szt + tdx;   // offset in result
   const int row = offset / rowdim;
   const int col = offset % rowdim;    // thread 0 computes QWYT[row][col]

   double resultrehihihi = 0.0;
   double resultrelohihi = 0.0;
   double resultrehilohi = 0.0;
   double resultrelolohi = 0.0;
   double resultrehihilo = 0.0;
   double resultrelohilo = 0.0;
   double resultrehilolo = 0.0;
   double resultrelololo = 0.0;
   double resultimhihihi = 0.0;
   double resultimlohihi = 0.0;
   double resultimhilohi = 0.0;
   double resultimlolohi = 0.0;
   double resultimhihilo = 0.0;
   double resultimlohilo = 0.0;
   double resultimhilolo = 0.0;
   double resultimlololo = 0.0;
   double a_rehihihi,a_relohihi,a_rehilohi,a_relolohi;
   double a_rehihilo,a_relohilo,a_rehilolo,a_relololo;
   double a_imhihihi,a_imlohihi,a_imhilohi,a_imlolohi;
   double a_imhihilo,a_imlohilo,a_imhilolo,a_imlololo;
   double b_rehihihi,b_relohihi,b_rehilohi,b_relolohi;
   double b_rehihilo,b_relohilo,b_rehilolo,b_relololo;
   double b_imhihihi,b_imlohihi,b_imhilohi,b_imlolohi;
   double b_imhihilo,b_imlohilo,b_imhilolo,b_imlololo;
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;
   int Qidx,WYTidx;

   for(int k=0; k<rowdim; k++)          // run over rowdim, not just szt
   {                                    // coloff shifts by col*row elements
      Qidx = row*dim + coloff + k;
      __syncthreads();
      a_rehihihi = Qrehihihi[Qidx];     // row = bdx,
      a_relohihi = Qrelohihi[Qidx];
      a_rehilohi = Qrehilohi[Qidx];
      a_relolohi = Qrelolohi[Qidx];
      a_rehihilo = Qrehihilo[Qidx];
      a_relohilo = Qrelohilo[Qidx];
      a_rehilolo = Qrehilolo[Qidx];
      a_relololo = Qrelololo[Qidx];
      a_imhihihi = Qimhihihi[Qidx];     // if dim == szt, coloff == 0
      a_imlohihi = Qimlohihi[Qidx];
      a_imhilohi = Qimhilohi[Qidx];
      a_imlolohi = Qimlolohi[Qidx];
      a_imhihilo = Qimhihilo[Qidx];
      a_imlohilo = Qimlohilo[Qidx];
      a_imhilolo = Qimhilolo[Qidx];
      a_imlololo = Qimlololo[Qidx];
      WYTidx = k*rowdim + col;
      __syncthreads();
      b_rehihihi = WYTrehihihi[WYTidx];   // if(dim == szt) then col = tdx
      b_relohihi = WYTrelohihi[WYTidx];
      b_rehilohi = WYTrehilohi[WYTidx];
      b_relolohi = WYTrelolohi[WYTidx];
      b_rehihilo = WYTrehihilo[WYTidx];
      b_relohilo = WYTrelohilo[WYTidx];
      b_rehilolo = WYTrehilolo[WYTidx];
      b_relololo = WYTrelololo[WYTidx];
      b_imhihihi = WYTimhihihi[WYTidx];
      b_imlohihi = WYTimlohihi[WYTidx];
      b_imhilohi = WYTimhilohi[WYTidx];
      b_imlolohi = WYTimlolohi[WYTidx];
      b_imhihilo = WYTimhihilo[WYTidx];
      b_imlohilo = WYTimlohilo[WYTidx];
      b_imhilolo = WYTimhilolo[WYTidx];
      b_imlololo = WYTimlololo[WYTidx];
      // result = result + a*b;
      // resultre = resultre + a_re*b_re - a_im*b_im;
      __syncthreads();
      odg_mul(a_rehihihi,a_relohihi,a_rehilohi,a_relolohi,
              a_rehihilo,a_relohilo,a_rehilolo,a_relololo,
              b_rehihihi,b_relohihi,b_rehilohi,b_relolohi,
              b_rehihilo,b_relohilo,b_rehilolo,b_relololo,
              &acchihihi,&acclohihi,&acchilohi,&acclolohi,
              &acchihilo,&acclohilo,&acchilolo,&acclololo);
      odg_inc(&resultrehihihi,&resultrelohihi,&resultrehilohi,&resultrelolohi,
              &resultrehihilo,&resultrelohilo,&resultrehilolo,&resultrelololo,
                    acchihihi,      acclohihi,      acchilohi,      acclolohi,
                    acchihilo,      acclohilo,      acchilolo,      acclololo);
      odg_mul(a_imhihihi,a_imlohihi,a_imhilohi,a_imlolohi,
              a_imhihilo,a_imlohilo,a_imhilolo,a_imlololo,
              b_imhihihi,b_imlohihi,b_imhilohi,b_imlolohi,
              b_imhihilo,b_imlohilo,b_imhilolo,b_imlololo,
              &acchihihi,&acclohihi,&acchilohi,&acclolohi,
              &acchihilo,&acclohilo,&acchilolo,&acclololo);
      odg_dec(&resultrehihihi,&resultrelohihi,&resultrehilohi,&resultrelolohi,
              &resultrehihilo,&resultrelohilo,&resultrehilolo,&resultrelololo,
                    acchihihi,      acclohihi,      acchilohi,      acclolohi,
                    acchihilo,      acclohilo,      acchilolo,      acclololo);
      // resultim = resultim + a_im*b_re + a_re*b_im;
      odg_mul(a_imhihihi,a_imlohihi,a_imhilohi,a_imlolohi,
              a_imhihilo,a_imlohilo,a_imhilolo,a_imlololo,
              b_rehihihi,b_relohihi,b_rehilohi,b_relolohi,
              b_rehihilo,b_relohilo,b_rehilolo,b_relololo,
              &acchihihi,&acclohihi,&acchilohi,&acclolohi,
              &acchihilo,&acclohilo,&acchilolo,&acclololo);
      odg_inc(&resultimhihihi,&resultimlohihi,&resultimhilohi,&resultimlolohi,
              &resultimhihilo,&resultimlohilo,&resultimhilolo,&resultimlololo,
                    acchihihi,      acclohihi,      acchilohi,      acclolohi,
                    acchihilo,      acclohilo,      acchilolo,      acclololo);
      odg_mul(a_rehihihi,a_relohihi,a_rehilohi,a_relolohi,
              a_rehihilo,a_relohilo,a_rehilolo,a_relololo,
              b_imhihihi,b_imlohihi,b_imhilohi,b_imlolohi,
              b_imhihilo,b_imlohilo,b_imhilolo,b_imlololo,
              &acchihihi,&acclohihi,&acchilohi,&acclolohi,
              &acchihilo,&acclohilo,&acchilolo,&acclololo);
      odg_inc(&resultimhihihi,&resultimlohihi,&resultimhilohi,&resultimlolohi,
              &resultimhihilo,&resultimlohilo,&resultimhilolo,&resultimlololo,
                    acchihihi,      acclohihi,      acchilohi,      acclolohi,
                    acchihilo,      acclohilo,      acchilolo,      acclololo);
   }
   __syncthreads();
   QWYTrehihihi[offset] = resultrehihihi;  // no column offset in saving QWYT
   QWYTrelohihi[offset] = resultrelohihi;
   QWYTrehilohi[offset] = resultrehilohi;
   QWYTrelolohi[offset] = resultrelolohi;
   QWYTrehihilo[offset] = resultrehihilo;
   QWYTrelohilo[offset] = resultrelohilo;
   QWYTrehilolo[offset] = resultrehilolo;
   QWYTrelololo[offset] = resultrelololo;
   QWYTimhihihi[offset] = resultimhihihi;
   QWYTimlohihi[offset] = resultimlohihi;
   QWYTimhilohi[offset] = resultimhilohi;
   QWYTimlolohi[offset] = resultimlolohi;
   QWYTimhihilo[offset] = resultimhihilo;
   QWYTimlohilo[offset] = resultimlohilo;
   QWYTimhilolo[offset] = resultimhilolo;
   QWYTimlololo[offset] = resultimlololo;
}

__global__ void dbl8_small_YWTC
 ( int nrows, int ncols, int rowdim, int coldim, int szt,
   int rowoff, int coloff,
   double *YWThihihi, double *YWTlohihi, double *YWThilohi, double *YWTlolohi,
   double *YWThihilo, double *YWTlohilo, double *YWThilolo, double *YWTlololo,
   double *Chihihi, double *Clohihi, double *Chilohi, double *Clolohi,
   double *Chihilo, double *Clohilo, double *Chilolo, double *Clololo,
   double *YWTChihihi, double *YWTClohihi,
   double *YWTChilohi, double *YWTClolohi,
   double *YWTChihilo, double *YWTClohilo,
   double *YWTChilolo, double *YWTClololo )
{
   const int bdx = blockIdx.x;         // bdx*szt done by previous blocks
   const int tdx = threadIdx.x;        // index of thread in block
   const int offset = bdx*szt + tdx;   // offset in result
   const int row = offset / coldim;    // 1st thread does YWTC[row][col]
   const int col = offset % coldim;
   const int colCoff0 = (coloff+col)*nrows + rowoff; // 1st element in C

   double resulthihihi = 0.0;
   double resultlohihi = 0.0;
   double resulthilohi = 0.0;
   double resultlolohi = 0.0;
   double resulthihilo = 0.0;
   double resultlohilo = 0.0;
   double resulthilolo = 0.0;
   double resultlololo = 0.0;
   double ahihihi,alohihi,ahilohi,alolohi;
   double ahihilo,alohilo,ahilolo,alololo;
   double bhihihi,blohihi,bhilohi,blolohi;
   double bhihilo,blohilo,bhilolo,blololo;
   double chihihi,clohihi,chilohi,clolohi;
   double chihilo,clohilo,chilolo,clololo;
   int idx;

   for(int k=0; k<rowdim; k++)         // innermost loop runs over rowdim
   {
      idx = row*rowdim + k;
      __syncthreads();
      ahihihi = YWThihihi[idx];        // YWT is stored row by row
      alohihi = YWTlohihi[idx];
      ahilohi = YWThilohi[idx];
      alolohi = YWTlolohi[idx];
      ahihilo = YWThihilo[idx];
      alohilo = YWTlohilo[idx];
      ahilolo = YWThilolo[idx];
      alololo = YWTlololo[idx];
      idx = colCoff0 + k;
      __syncthreads();
      bhihihi = Chihihi[idx];         // but C is stored column by column
      blohihi = Clohihi[idx];
      bhilohi = Chilohi[idx];
      blolohi = Clolohi[idx];
      bhihilo = Chihilo[idx];
      blohilo = Clohilo[idx];
      bhilolo = Chilolo[idx];
      blololo = Clololo[idx];
      // result = result + a*b;
      __syncthreads();
      odg_mul( ahihihi, alohihi, ahilohi, alolohi,
               ahihilo, alohilo, ahilolo, alololo,
               bhihihi, blohihi, bhilohi, blolohi,
               bhihilo, blohilo, bhilolo, blololo,
              &chihihi,&clohihi,&chilohi,&clolohi,
              &chihilo,&clohilo,&chilolo,&clololo);
      odg_inc(&resulthihihi,&resultlohihi,&resulthilohi,&resultlolohi,
              &resulthihilo,&resultlohilo,&resulthilolo,&resultlololo,
                    chihihi,      clohihi,      chilohi,      clolohi,
                    chihilo,      clohilo,      chilolo,      clololo);
   }
   idx = (coloff + col)*nrows + (rowoff + row);
   __syncthreads();
   YWTChihihi[idx] = resulthihihi;
   YWTClohihi[idx] = resultlohihi;
   YWTChilohi[idx] = resulthilohi;
   YWTClolohi[idx] = resultlolohi;
   YWTChihilo[idx] = resulthihilo;
   YWTClohilo[idx] = resultlohilo;
   YWTChilolo[idx] = resulthilolo;
   YWTClololo[idx] = resultlololo;
}

__global__ void cmplx8_small_YWHC
 ( int nrows, int ncols, int rowdim, int coldim, int szt,
   int rowoff, int coloff,
   double *YWTrehihihi, double *YWTrelohihi,
   double *YWTrehilohi, double *YWTrelolohi,
   double *YWTrehihilo, double *YWTrelohilo,
   double *YWTrehilolo, double *YWTrelololo,
   double *YWTimhihihi, double *YWTimlohihi,
   double *YWTimhilohi, double *YWTimlolohi,
   double *YWTimhihilo, double *YWTimlohilo,
   double *YWTimhilolo, double *YWTimlololo,
   double *Crehihihi, double *Crelohihi, double *Crehilohi, double *Crelolohi,
   double *Crehihilo, double *Crelohilo, double *Crehilolo, double *Crelololo,
   double *Cimhihihi, double *Cimlohihi, double *Cimhilohi, double *Cimlolohi,
   double *Cimhihilo, double *Cimlohilo, double *Cimhilolo, double *Cimlololo,
   double *YWTCrehihihi, double *YWTCrelohihi,
   double *YWTCrehilohi, double *YWTCrelolohi,
   double *YWTCrehihilo, double *YWTCrelohilo,
   double *YWTCrehilolo, double *YWTCrelololo,
   double *YWTCimhihihi, double *YWTCimlohihi,
   double *YWTCimhilohi, double *YWTCimlolohi,
   double *YWTCimhihilo, double *YWTCimlohilo,
   double *YWTCimhilolo, double *YWTCimlololo )
{
   const int bdx = blockIdx.x;         // bdx*szt done by previous blocks
   const int tdx = threadIdx.x;        // index of thread in block
   const int offset = bdx*szt + tdx;   // offset in result
   const int row = offset / coldim;    // 1st thread does YWTC[row][col]
   const int col = offset % coldim;
   const int colCoff0 = (coloff+col)*nrows + rowoff; // 1st element in C
   int idx;

   double resultrehihihi = 0.0;
   double resultrelohihi = 0.0;
   double resultrehilohi = 0.0;
   double resultrelolohi = 0.0;
   double resultrehihilo = 0.0;
   double resultrelohilo = 0.0;
   double resultrehilolo = 0.0;
   double resultrelololo = 0.0;
   double resultimhihihi = 0.0;
   double resultimlohihi = 0.0;
   double resultimhilohi = 0.0;
   double resultimlolohi = 0.0;
   double resultimhihilo = 0.0;
   double resultimlohilo = 0.0;
   double resultimhilolo = 0.0;
   double resultimlololo = 0.0;
   double a_hihihi,a_lohihi,a_hilohi,a_lolohi;
   double a_hihilo,a_lohilo,a_hilolo,a_lololo;
   double b_hihihi,b_lohihi,b_hilohi,b_lolohi;
   double b_hihilo,b_lohilo,b_hilolo,b_lololo;
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;
   int YWTidx,Cidx;

   for(int k=0; k<rowdim; k++)            // innermost loop runs over rowdim
   {
      YWTidx = row*rowdim + k;
      __syncthreads();
      a_hihihi = YWTrehihihi[YWTidx];   // YWT is stored row by row
      a_lohihi = YWTrelohihi[YWTidx];
      a_hilohi = YWTrehilohi[YWTidx];
      a_lolohi = YWTrelolohi[YWTidx];
      a_hihilo = YWTrehihilo[YWTidx];
      a_lohilo = YWTrelohilo[YWTidx];
      a_hilolo = YWTrehilolo[YWTidx];
      a_lololo = YWTrelololo[YWTidx];
      Cidx = colCoff0 + k;
      __syncthreads();
      b_hihihi = Crehihihi[Cidx];    // but C is stored column by column
      b_lohihi = Crelohihi[Cidx];
      b_hilohi = Crehilohi[Cidx];
      b_lolohi = Crelolohi[Cidx];
      b_hihilo = Crehihilo[Cidx];
      b_lohilo = Crelohilo[Cidx];
      b_hilolo = Crehilolo[Cidx];
      b_lololo = Crelololo[Cidx];
      // result = result + a*b;
      // resultre = resultre + a_re*b_re - a_im*b_im;
      __syncthreads();
      odg_mul(a_hihihi,  a_lohihi,  a_hilohi,  a_lolohi,
              a_hihilo,  a_lohilo,  a_hilolo,  a_lololo,
              b_hihihi,  b_lohihi,  b_hilohi,  b_lolohi,
              b_hihilo,  b_lohilo,  b_hilolo,  b_lololo,
            &acchihihi,&acclohihi,&acchilohi,&acclolohi,
            &acchihilo,&acclohilo,&acchilolo,&acclololo);
      odg_inc(&resultrehihihi,&resultrelohihi,&resultrehilohi,&resultrelolohi,
              &resultrehihilo,&resultrelohilo,&resultrehilolo,&resultrelololo,
                    acchihihi,      acclohihi,      acchilohi,      acclolohi,
                    acchihilo,      acclohilo,      acchilolo,      acclololo);

      a_hihihi = YWTimhihihi[YWTidx];
      a_lohihi = YWTimlohihi[YWTidx];
      a_hilohi = YWTimhilohi[YWTidx];
      a_lolohi = YWTimlolohi[YWTidx];
      a_hihilo = YWTimhihilo[YWTidx];
      a_lohilo = YWTimlohilo[YWTidx];
      a_hilolo = YWTimhilolo[YWTidx];
      a_lololo = YWTimlololo[YWTidx];
      b_hihihi = Cimhihihi[Cidx];
      b_lohihi = Cimlohihi[Cidx];
      b_hilohi = Cimhilohi[Cidx];
      b_lolohi = Cimlolohi[Cidx];
      b_hihilo = Cimhihilo[Cidx];
      b_lohilo = Cimlohilo[Cidx];
      b_hilolo = Cimhilolo[Cidx];
      b_lololo = Cimlololo[Cidx];

      odg_mul(a_hihihi,  a_lohihi,  a_hilohi,  a_lolohi,
              a_hihilo,  a_lohilo,  a_hilolo,  a_lololo,
              b_hihihi,  b_lohihi,  b_hilohi,  b_lolohi,
              b_hihilo,  b_lohilo,  b_hilolo,  b_lololo,
            &acchihihi,&acclohihi,&acchilohi,&acclolohi,
            &acchihilo,&acclohilo,&acchilolo,&acclololo);
      odg_dec(&resultrehihihi,&resultrelohihi,&resultrehilohi,&resultrelolohi,
              &resultrehihilo,&resultrelohilo,&resultrehilolo,&resultrelololo,
                    acchihihi,      acclohihi,      acchilohi,      acclolohi,
                    acchihilo,      acclohilo,      acchilolo,      acclololo);
      // resultim = resultim + a_im*b_re + a_re*b_im;
      // a_im is already available
      b_hihihi = Crehihihi[Cidx];
      b_lohihi = Crelohihi[Cidx];
      b_hilohi = Crehilohi[Cidx];
      b_lolohi = Crelolohi[Cidx];
      b_hihilo = Crehihilo[Cidx];
      b_lohilo = Crelohilo[Cidx];
      b_hilolo = Crehilolo[Cidx];
      b_lololo = Crelololo[Cidx];
      odg_mul(a_hihihi,  a_lohihi,  a_hilohi,  a_lolohi,
              a_hihilo,  a_lohilo,  a_hilolo,  a_lololo,
              b_hihihi,  b_lohihi,  b_hilohi,  b_lolohi,
              b_hihilo,  b_lohilo,  b_hilolo,  b_lololo,
            &acchihihi,&acclohihi,&acchilohi,&acclolohi,
            &acchihilo,&acclohilo,&acchilolo,&acclololo);
      odg_inc(&resultimhihihi,&resultimlohihi,&resultimhilohi,&resultimlolohi,
              &resultimhihilo,&resultimlohilo,&resultimhilolo,&resultimlololo,
                    acchihihi,      acclohihi,      acchilohi,      acclolohi,
                    acchihilo,      acclohilo,      acchilolo,      acclololo);

      a_hihihi = YWTrehihihi[YWTidx];
      a_lohihi = YWTrelohihi[YWTidx];
      a_hilohi = YWTrehilohi[YWTidx];
      a_lolohi = YWTrelolohi[YWTidx];
      a_hihilo = YWTrehihilo[YWTidx];
      a_lohilo = YWTrelohilo[YWTidx];
      a_hilolo = YWTrehilolo[YWTidx];
      a_lololo = YWTrelololo[YWTidx];
      b_hihihi = Cimhihihi[Cidx];
      b_lohihi = Cimlohihi[Cidx];
      b_hilohi = Cimhilohi[Cidx];
      b_lolohi = Cimlolohi[Cidx];
      b_hihilo = Cimhihilo[Cidx];
      b_lohilo = Cimlohilo[Cidx];
      b_hilolo = Cimhilolo[Cidx];
      b_lololo = Cimlololo[Cidx];
      odg_mul(a_hihihi,  a_lohihi,  a_hilohi,  a_lolohi,
              a_hihilo,  a_lohilo,  a_hilolo,  a_lololo,
              b_hihihi,  b_lohihi,  b_hilohi,  b_lolohi,
              b_hihilo,  b_lohilo,  b_hilolo,  b_lololo,
            &acchihihi,&acclohihi,&acchilohi,&acclolohi,
            &acchihilo,&acclohilo,&acchilolo,&acclololo);
      odg_inc(&resultimhihihi,&resultimlohihi,&resultimhilohi,&resultimlolohi,
              &resultimhihilo,&resultimlohilo,&resultimhilolo,&resultimlololo,
                    acchihihi,      acclohihi,      acchilohi,      acclolohi,
                    acchihilo,      acclohilo,      acchilolo,      acclololo);
   }
   __syncthreads();
   idx = (coloff + col)*nrows + (rowoff + row);
   YWTCrehihihi[idx] = resultrehihihi;
   YWTCrelohihi[idx] = resultrelohihi;
   YWTCrehilohi[idx] = resultrehilohi;
   YWTCrelolohi[idx] = resultrelolohi;
   YWTCrehihilo[idx] = resultrehihilo;
   YWTCrelohilo[idx] = resultrelohilo;
   YWTCrehilolo[idx] = resultrehilolo;
   YWTCrelololo[idx] = resultrelololo;
   YWTCimhihihi[idx] = resultimhihihi;
   YWTCimlohihi[idx] = resultimlohihi;
   YWTCimhilohi[idx] = resultimhilohi;
   YWTCimlolohi[idx] = resultimlolohi;
   YWTCimhihilo[idx] = resultimhihilo;
   YWTCimlohilo[idx] = resultimlohilo;
   YWTCimhilolo[idx] = resultimhilolo;
   YWTCimlololo[idx] = resultimlololo;
}

__global__ void dbl8_small_Qupdate
 ( int dim, int rowdim, int szt, int coloff,
   double *Qhihihi, double *Qlohihi, double *Qhilohi, double *Qlolohi,
   double *Qhihilo, double *Qlohilo, double *Qhilolo, double *Qlololo,
   double *QWYThihihi, double *QWYTlohihi,
   double *QWYThilohi, double *QWYTlolohi,
   double *QWYThihilo, double *QWYTlohilo,
   double *QWYThilolo, double *QWYTlololo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int offset = bdx*szt + tdx;   // offset in result
   const int row = offset / rowdim;
   const int col = offset % rowdim;
   const int idx1 = row*dim + coloff + col;

   double ahihihi,alohihi,ahilohi,alolohi;
   double ahihilo,alohilo,ahilolo,alololo;
   double bhihihi,blohihi,bhilohi,blolohi;
   double bhihilo,blohilo,bhilolo,blololo;

   ahihihi = Qhihihi[idx1];   // row = bdx, if dim == szt, coloff == 0
   alohihi = Qlohihi[idx1];
   ahilohi = Qhilohi[idx1];
   alolohi = Qlolohi[idx1];
   ahihilo = Qhihilo[idx1];
   alohilo = Qlohilo[idx1];
   ahilolo = Qhilolo[idx1];
   alololo = Qlololo[idx1];
   __syncthreads();
   bhihihi = QWYThihihi[offset];  // if(dim == szt) then col = tdx
   blohihi = QWYTlohihi[offset];
   bhilohi = QWYThilohi[offset];
   blolohi = QWYTlolohi[offset];
   bhihilo = QWYThihilo[offset];
   blohilo = QWYTlohilo[offset];
   bhilolo = QWYThilolo[offset];
   blololo = QWYTlololo[offset];
   // a = a + b;
   __syncthreads();
   odg_inc(&ahihihi,&alohihi,&ahilohi,&alolohi,
           &ahihilo,&alohilo,&ahilolo,&alololo,
            bhihihi, blohihi, bhilohi, blolohi,
            bhihilo, blohilo, bhilolo, blololo);
   __syncthreads();
   Qhihihi[idx1] = ahihihi;
   Qlohihi[idx1] = alohihi;
   Qhilohi[idx1] = ahilohi;
   Qlolohi[idx1] = alolohi;
   Qhihilo[idx1] = ahihilo;
   Qlohilo[idx1] = alohilo;
   Qhilolo[idx1] = ahilolo;
   Qlololo[idx1] = alololo;
}

__global__ void cmplx8_small_Qupdate
 ( int dim, int rowdim, int szt, int coloff,
   double *Qrehihihi, double *Qrelohihi, double *Qrehilohi, double *Qrelolohi,
   double *Qrehihilo, double *Qrelohilo, double *Qrehilolo, double *Qrelololo,
   double *Qimhihihi, double *Qimlohihi, double *Qimhilohi, double *Qimlolohi,
   double *Qimhihilo, double *Qimlohilo, double *Qimhilolo, double *Qimlololo,
   double *QWYTrehihihi, double *QWYTrelohihi,
   double *QWYTrehilohi, double *QWYTrelolohi,
   double *QWYTrehihilo, double *QWYTrelohilo,
   double *QWYTrehilolo, double *QWYTrelololo,
   double *QWYTimhihihi, double *QWYTimlohihi,
   double *QWYTimhilohi, double *QWYTimlolohi,
   double *QWYTimhihilo, double *QWYTimlohilo,
   double *QWYTimhilolo, double *QWYTimlololo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int offset = bdx*szt + tdx;   // offset in result
   const int row = offset / rowdim;
   const int col = offset % rowdim;
   const int idx1 = row*dim + coloff + col;

   double a_rehihihi,a_relohihi,a_rehilohi,a_relolohi;
   double a_rehihilo,a_relohilo,a_rehilolo,a_relololo;
   double a_imhihihi,a_imlohihi,a_imhilohi,a_imlolohi;
   double a_imhihilo,a_imlohilo,a_imhilolo,a_imlololo;
   double b_rehihihi,b_relohihi,b_rehilohi,b_relolohi;
   double b_rehihilo,b_relohilo,b_rehilolo,b_relololo;
   double b_imhihihi,b_imlohihi,b_imhilohi,b_imlolohi;
   double b_imhihilo,b_imlohilo,b_imhilolo,b_imlololo;

   a_rehihihi = Qrehihihi[idx1];  // row = bdx, if dim == szt, coloff == 0
   a_relohihi = Qrelohihi[idx1];
   a_rehilohi = Qrehilohi[idx1];
   a_relolohi = Qrelolohi[idx1];
   a_rehihilo = Qrehihilo[idx1];
   a_relohilo = Qrelohilo[idx1];
   a_rehilolo = Qrehilolo[idx1];
   a_relololo = Qrelololo[idx1];
   a_imhihihi = Qimhihihi[idx1];
   a_imlohihi = Qimlohihi[idx1];
   a_imhilohi = Qimhilohi[idx1];
   a_imlolohi = Qimlolohi[idx1];
   a_imhihilo = Qimhihilo[idx1];
   a_imlohilo = Qimlohilo[idx1];
   a_imhilolo = Qimhilolo[idx1];
   a_imlololo = Qimlololo[idx1];
   __syncthreads();
   b_rehihihi = QWYTrehihihi[offset];  // if(dim == szt) then col = tdx
   b_relohihi = QWYTrelohihi[offset];
   b_rehilohi = QWYTrehilohi[offset];
   b_relolohi = QWYTrelolohi[offset];
   b_rehihilo = QWYTrehihilo[offset];
   b_relohilo = QWYTrelohilo[offset];
   b_rehilolo = QWYTrehilolo[offset];
   b_relololo = QWYTrelololo[offset];
   b_imhihihi = QWYTimhihihi[offset];
   b_imlohihi = QWYTimlohihi[offset];
   b_imhilohi = QWYTimhilohi[offset];
   b_imlolohi = QWYTimlolohi[offset];
   b_imhihilo = QWYTimhihilo[offset];
   b_imlohilo = QWYTimlohilo[offset];
   b_imhilolo = QWYTimhilolo[offset];
   b_imlololo = QWYTimlololo[offset];
   // a_re = a_re + b_re;
   __syncthreads();
   odg_inc(&a_rehihihi,&a_relohihi,&a_rehilohi,&a_relolohi,
           &a_rehihilo,&a_relohilo,&a_rehilolo,&a_relololo,
            b_rehihihi, b_relohihi, b_rehilohi, b_relolohi,
            b_rehihilo, b_relohilo, b_rehilolo, b_relololo);
   // a_im = a_im + b_im;
   odg_inc(&a_imhihihi,&a_imlohihi,&a_imhilohi,&a_imlolohi,
           &a_imhihilo,&a_imlohilo,&a_imhilolo,&a_imlololo,
            b_imhihihi, b_imlohihi, b_imhilohi, b_imlolohi,
            b_imhihilo, b_imlohilo, b_imhilolo, b_imlololo);

   __syncthreads();
   Qrehihihi[idx1] = a_rehihihi;
   Qrelohihi[idx1] = a_relohihi;
   Qrehilohi[idx1] = a_rehilohi;
   Qrelolohi[idx1] = a_relolohi;
   Qrehihilo[idx1] = a_rehihilo;
   Qrelohilo[idx1] = a_relohilo;
   Qrehilolo[idx1] = a_rehilolo;
   Qrelololo[idx1] = a_relololo;
   Qimhihihi[idx1] = a_imhihihi;
   Qimlohihi[idx1] = a_imlohihi;
   Qimhilohi[idx1] = a_imhilohi;
   Qimlolohi[idx1] = a_imlolohi;
   Qimhihilo[idx1] = a_imhihilo;
   Qimlohilo[idx1] = a_imlohilo;
   Qimhilolo[idx1] = a_imhilolo;
   Qimlololo[idx1] = a_imlololo;
}

__global__ void dbl8_small_R_add_YWTC
 ( int nrows, int coldim, int szt, int rowoff, int coloff,
   double *Rhihihi, double *Rlohihi, double *Rhilohi, double *Rlolohi,
   double *Rhihilo, double *Rlohilo, double *Rhilolo, double *Rlololo,
   double *YWTChihihi, double *YWTClohihi,
   double *YWTChilohi, double *YWTClolohi,
   double *YWTChihilo, double *YWTClohilo,
   double *YWTChilolo, double *YWTClololo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int offset = bdx*szt + tdx;   // offset in result
   const int row = offset / coldim;    // thread updates R[row][col]
   const int col = offset % coldim;
   const int idx = (coloff + col)*nrows + (rowoff + row);
 
   double ahihihi,alohihi,ahilohi,alolohi;
   double ahihilo,alohilo,ahilolo,alololo;
   double bhihihi,blohihi,bhilohi,blolohi;
   double bhihilo,blohilo,bhilolo,blololo;
   
   ahihihi = Rhihihi[idx];
   alohihi = Rlohihi[idx];
   ahilohi = Rhilohi[idx];
   alolohi = Rlolohi[idx];
   ahihilo = Rhihilo[idx];
   alohilo = Rlohilo[idx];
   ahilolo = Rhilolo[idx];
   alololo = Rlololo[idx];
   __syncthreads();
   bhihihi = YWTChihihi[idx];
   blohihi = YWTClohihi[idx];
   bhilohi = YWTChilohi[idx];
   blolohi = YWTClolohi[idx];
   bhihilo = YWTChihilo[idx];
   blohilo = YWTClohilo[idx];
   bhilolo = YWTChilolo[idx];
   blololo = YWTClololo[idx];
   // a = a + b;
   __syncthreads();
   odg_inc(&ahihihi,&alohihi,&ahilohi,&alolohi,
           &ahihilo,&alohilo,&ahilolo,&alololo,
            bhihihi, blohihi, bhilohi, blolohi,
            bhihilo, blohilo, bhilolo, blololo);
  
   __syncthreads();
   Rhihihi[idx] = ahihihi;
   Rlohihi[idx] = alohihi;
   Rhilohi[idx] = ahilohi;
   Rlolohi[idx] = alolohi;
   Rhihilo[idx] = ahihilo;
   Rlohilo[idx] = alohilo;
   Rhilolo[idx] = ahilolo;
   Rlololo[idx] = alololo;
}

__global__ void cmplx8_small_R_add_YWHC
 ( int nrows, int coldim, int szt, int rowoff, int coloff,
   double *Rrehihihi, double *Rrelohihi, double *Rrehilohi, double *Rrelolohi,
   double *Rrehihilo, double *Rrelohilo, double *Rrehilolo, double *Rrelololo,
   double *Rimhihihi, double *Rimlohihi, double *Rimhilohi, double *Rimlolohi,
   double *Rimhihilo, double *Rimlohilo, double *Rimhilolo, double *Rimlololo,
   double *YWTCrehihihi, double *YWTCrelohihi,
   double *YWTCrehilohi, double *YWTCrelolohi,
   double *YWTCrehihilo, double *YWTCrelohilo,
   double *YWTCrehilolo, double *YWTCrelololo,
   double *YWTCimhihihi, double *YWTCimlohihi,
   double *YWTCimhilohi, double *YWTCimlolohi,
   double *YWTCimhihilo, double *YWTCimlohilo,
   double *YWTCimhilolo, double *YWTCimlololo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int offset = bdx*szt + tdx;   // offset in result
   const int row = offset / coldim;    // thread updates R[row][col]
   const int col = offset % coldim;
   const int idx = (coloff + col)*nrows + (rowoff + row);
 
   double a_hihihi,a_lohihi,a_hilohi,a_lolohi;
   double a_hihilo,a_lohilo,a_hilolo,a_lololo;
   double b_hihihi,b_lohihi,b_hilohi,b_lolohi;
   double b_hihilo,b_lohilo,b_hilolo,b_lololo;
   
   a_hihihi = Rrehihihi[idx];
   a_lohihi = Rrelohihi[idx];
   a_hilohi = Rrehilohi[idx];
   a_lolohi = Rrelolohi[idx];
   a_hihilo = Rrehihilo[idx];
   a_lohilo = Rrelohilo[idx];
   a_hilolo = Rrehilolo[idx];
   a_lololo = Rrelololo[idx];
   b_hihihi = YWTCrehihihi[idx];
   b_lohihi = YWTCrelohihi[idx];
   b_hilohi = YWTCrehilohi[idx];
   b_lolohi = YWTCrelolohi[idx];
   b_hihilo = YWTCrehihilo[idx];
   b_lohilo = YWTCrelohilo[idx];
   b_hilolo = YWTCrehilolo[idx];
   b_lololo = YWTCrelololo[idx];
   // a_re = a_re + b_re;
   __syncthreads();
   odg_inc(&a_hihihi,&a_lohihi,&a_hilohi,&a_lolohi,
           &a_hihilo,&a_lohilo,&a_hilolo,&a_lololo,
            b_hihihi, b_lohihi, b_hilohi, b_lolohi,
            b_hihilo, b_lohilo, b_hilolo, b_lololo);
   __syncthreads();
   Rrehihihi[idx] = a_hihihi;
   Rrelohihi[idx] = a_lohihi;
   Rrehilohi[idx] = a_hilohi;
   Rrelolohi[idx] = a_lolohi;
   Rrehihilo[idx] = a_hihilo;
   Rrelohilo[idx] = a_lohilo;
   Rrehilolo[idx] = a_hilolo;
   Rrelololo[idx] = a_lololo;

   // a_im = a_im + b_im;
   a_hihihi = Rimhihihi[idx];
   a_lohihi = Rimlohihi[idx];
   a_hilohi = Rimhilohi[idx];
   a_lolohi = Rimlolohi[idx];
   a_hihilo = Rimhihilo[idx];
   a_lohilo = Rimlohilo[idx];
   a_hilolo = Rimhilolo[idx];
   a_lololo = Rimlololo[idx];
   b_hihihi = YWTCimhihihi[idx];
   b_lohihi = YWTCimlohihi[idx];
   b_hilohi = YWTCimhilohi[idx];
   b_lolohi = YWTCimlolohi[idx];
   b_hihilo = YWTCimhihilo[idx];
   b_lohilo = YWTCimlohilo[idx];
   b_hilolo = YWTCimhilolo[idx];
   b_lololo = YWTCimlololo[idx];
   __syncthreads();
   odg_inc(&a_hihihi,&a_lohihi,&a_hilohi,&a_lolohi,
           &a_hihilo,&a_lohilo,&a_hilolo,&a_lololo,
            b_hihihi, b_lohihi, b_hilohi, b_lolohi,
            b_hihilo, b_lohilo, b_hilolo, b_lololo);
   __syncthreads();
   Rimhihihi[idx] = a_hihihi;
   Rimlohihi[idx] = a_lohihi;
   Rimhilohi[idx] = a_hilohi;
   Rimlolohi[idx] = a_lolohi;
   Rimhihilo[idx] = a_hihilo;
   Rimlohilo[idx] = a_lohilo;
   Rimhilolo[idx] = a_hilolo;
   Rimlololo[idx] = a_lololo;
}

void GPU_dbl8_small_house
 ( int nrows, int ncols, int szt, int nbt,
   int colidx, int nrows1, int k, int L,
   double *Ahihihi_h, double *Alohihi_h, double *Ahilohi_h, double *Alolohi_h,
   double *Ahihilo_h, double *Alohilo_h, double *Ahilolo_h, double *Alololo_h,
   double *Ahihihi_d, double *Alohihi_d, double *Ahilohi_d, double *Alolohi_d,
   double *Ahihilo_d, double *Alohilo_d, double *Ahilolo_d, double *Alololo_d,
   double *vhihihi_h, double *vlohihi_h, double *vhilohi_h, double *vlolohi_h,
   double *vhihilo_h, double *vlohilo_h, double *vhilolo_h, double *vlololo_h,
   double *Vhihihi_d, double *Vlohihi_d, double *Vhilohi_d, double *Vlolohi_d,
   double *Vhihilo_d, double *Vlohilo_d, double *Vhilolo_d, double *Vlololo_d,
   double *betahihihi_h, double *betalohihi_h,
   double *betahilohi_h, double *betalolohi_h,
   double *betahihilo_h, double *betalohilo_h,
   double *betahilolo_h, double *betalololo_h,
   double *betahihihi_d, double *betalohihi_d,
   double *betahilohi_d, double *betalolohi_d,
   double *betahihilo_d, double *betalohilo_d,
   double *betahilolo_d, double *betalololo_d,
   double *lapms, long long int *add, long long int *mul, long long int *div,
   long long int *sqrtfun, bool verbose )
{
   const int nrLog2 = ceil(log2((double) nrows1));
   const int rowidx = colidx*(nrows+1);       // start of number in A_h
   const int nVrows = nrows - k*szt;          // dimension of V matrix

   if(verbose)
   {
      cout << "nrows : " << nrows
           << "  nVrows : " << nVrows
           << "  ncols : " << ncols
           << "  szt : " << szt
           << "  nbt : " << nbt << endl;
      cout << "k : " << k 
           << "  L : " << L
           << "  nrows1 : " << nrows1
           << "  colidx : " << colidx
           << "  rowidx : " << rowidx << endl;
   }
   if(L > 0)
   {
      for(int i=0; i<L; i++)             // insert zeros
      {
         vhihihi_h[i] = 0.0;
         vlohihi_h[i] = 0.0;
         vhilohi_h[i] = 0.0;
         vlolohi_h[i] = 0.0;
         vhihilo_h[i] = 0.0;
         vlohilo_h[i] = 0.0;
         vhilolo_h[i] = 0.0;
         vlololo_h[i] = 0.0;
      }
      cudaMemcpy(&Vhihihi_d[L*nVrows],vhihihi_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vlohihi_d[L*nVrows],vlohihi_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vhilohi_d[L*nVrows],vhilohi_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vlolohi_d[L*nVrows],vlolohi_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vhihilo_d[L*nVrows],vhihilo_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vlohilo_d[L*nVrows],vlohilo_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vhilolo_d[L*nVrows],vhilolo_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vlololo_d[L*nVrows],vlololo_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
   }
   if(nrows1 == 0)
   {
      betahihihi_h[L] = 0.0;
      betalohihi_h[L] = 0.0;
      betahilohi_h[L] = 0.0;
      betalolohi_h[L] = 0.0;
      betahihilo_h[L] = 0.0;
      betalohilo_h[L] = 0.0;
      betahilolo_h[L] = 0.0;
      betalololo_h[L] = 0.0;
      vhihihi_h[0] = 1.0;
      vlohihi_h[0] = 0.0;
      vhilohi_h[0] = 0.0;
      vlolohi_h[0] = 0.0;
      vhihilo_h[0] = 1.0;
      vlohilo_h[0] = 0.0;
      vhilolo_h[0] = 0.0;
      vlololo_h[0] = 0.0;
      cudaMemcpy(&betahihihi_d[L],&betahihihi_h[L],sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&betalohihi_d[L],&betalohihi_h[L],sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&betahilohi_d[L],&betahilohi_h[L],sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&betalolohi_d[L],&betalolohi_h[L],sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&betahihilo_d[L],&betahihilo_h[L],sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&betalohilo_d[L],&betalohilo_h[L],sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&betahilolo_d[L],&betahilolo_h[L],sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&betalololo_d[L],&betalololo_h[L],sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vhihihi_d[L*nVrows+L],vhihihi_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vlohihi_d[L*nVrows+L],vlohihi_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vhilohi_d[L*nVrows+L],vhilohi_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vlolohi_d[L*nVrows+L],vlolohi_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vhihilo_d[L*nVrows+L],vhihilo_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vlohilo_d[L*nVrows+L],vlohilo_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vhilolo_d[L*nVrows+L],vhilolo_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vlololo_d[L*nVrows+L],vlololo_h,sizeof(double),
                 cudaMemcpyHostToDevice);
   }
   else
   {
      cudaEvent_t start,stop;           // to measure time spent by kernels 
      cudaEventCreate(&start);
      cudaEventCreate(&stop);
      float milliseconds;

      cudaEventRecord(start);
      dbl8_small_house<<<1,nrows1>>>
         (&Ahihihi_d[rowidx],&Alohihi_d[rowidx],
          &Ahilohi_d[rowidx],&Alolohi_d[rowidx],
          &Ahihilo_d[rowidx],&Alohilo_d[rowidx],
          &Ahilolo_d[rowidx],&Alololo_d[rowidx],
          &Ahihihi_d[rowidx+1],&Alohihi_d[rowidx+1],
          &Ahilohi_d[rowidx+1],&Alolohi_d[rowidx+1],
          &Ahihilo_d[rowidx+1],&Alohilo_d[rowidx+1],
          &Ahilolo_d[rowidx+1],&Alololo_d[rowidx+1],
          nrows1,nrLog2,&Vhihihi_d[L*nVrows+L],&Vlohihi_d[L*nVrows+L],
                        &Vhilohi_d[L*nVrows+L],&Vlolohi_d[L*nVrows+L],
                        &Vhihilo_d[L*nVrows+L],&Vlohilo_d[L*nVrows+L],
                        &Vhilolo_d[L*nVrows+L],&Vlololo_d[L*nVrows+L],
          &betahihihi_d[L],&betalohihi_d[L],&betahilohi_d[L],&betalolohi_d[L],
          &betahihilo_d[L],&betalohilo_d[L],&betahilolo_d[L],&betalololo_d[L]);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds,start,stop);
      *lapms += milliseconds;
      flopcount_dbl_small_house(nrows1,nrLog2,add,mul,div,sqrtfun);
   }
   cudaMemcpy(&betahihihi_h[L],&betahihihi_d[L],sizeof(double),
              cudaMemcpyDeviceToHost);
   cudaMemcpy(&betalohihi_h[L],&betalohihi_d[L],sizeof(double),
              cudaMemcpyDeviceToHost);
   cudaMemcpy(&betahilohi_h[L],&betahilohi_d[L],sizeof(double),
              cudaMemcpyDeviceToHost);
   cudaMemcpy(&betalolohi_h[L],&betalolohi_d[L],sizeof(double),
              cudaMemcpyDeviceToHost);
   cudaMemcpy(&betahihilo_h[L],&betahihilo_d[L],sizeof(double),
              cudaMemcpyDeviceToHost);
   cudaMemcpy(&betalohilo_h[L],&betalohilo_d[L],sizeof(double),
              cudaMemcpyDeviceToHost);
   cudaMemcpy(&betahilolo_h[L],&betahilolo_d[L],sizeof(double),
              cudaMemcpyDeviceToHost);
   cudaMemcpy(&betalololo_h[L],&betalololo_d[L],sizeof(double),
              cudaMemcpyDeviceToHost);
   if(verbose)
   {
      const size_t szhouse = nVrows*sizeof(double);

      cudaMemcpy(vhihihi_h,&Vhihihi_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vlohihi_h,&Vlohihi_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vhilohi_h,&Vhilohi_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vlolohi_h,&Vlolohi_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vhihilo_h,&Vhihilo_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vlohilo_h,&Vlohilo_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vhilolo_h,&Vhilolo_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vlololo_h,&Vlololo_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);

      cout << scientific << setprecision(16)
           << "beta[" << colidx << "] : "
           << betahihihi_h[L] << "  " << betalohihi_h[L] << endl
           << "          "
           << betahilohi_h[L] << "  " << betalolohi_h[L] << endl
           << "          "
           << betahihilo_h[L] << "  " << betalohilo_h[L] << endl
           << "          "
           << betahilolo_h[L] << "  " << betalololo_h[L] << endl;

      for(int i=0; i<nVrows; i++)
         cout << "v[" << i << "] : "
              << vhihihi_h[i] << "  " << vlohihi_h[i] << endl
              << "       "
              << vhilohi_h[i] << "  " << vlolohi_h[i] << endl
              << "       "
              << vhihilo_h[i] << "  " << vlohilo_h[i] << endl
              << "       "
              << vhilolo_h[i] << "  " << vlololo_h[i] << endl;
   }
}

void GPU_cmplx8_small_house
 ( int nrows, int ncols, int szt, int nbt,
   int colidx, int nrows1, int k, int L,
   double *Arehihihi_h, double *Arelohihi_h,
   double *Arehilohi_h, double *Arelolohi_h,
   double *Arehihilo_h, double *Arelohilo_h,
   double *Arehilolo_h, double *Arelololo_h,
   double *Aimhihihi_h, double *Aimlohihi_h,
   double *Aimhilohi_h, double *Aimlolohi_h,
   double *Aimhihilo_h, double *Aimlohilo_h,
   double *Aimhilolo_h, double *Aimlololo_h,
   double *Arehihihi_d, double *Arelohihi_d,
   double *Arehilohi_d, double *Arelolohi_d,
   double *Arehihilo_d, double *Arelohilo_d,
   double *Arehilolo_d, double *Arelololo_d,
   double *Aimhihihi_d, double *Aimlohihi_d,
   double *Aimhilohi_d, double *Aimlolohi_d,
   double *Aimhihilo_d, double *Aimlohilo_d,
   double *Aimhilolo_d, double *Aimlololo_d,
   double *vrehihihi_h, double *vrelohihi_h,
   double *vrehilohi_h, double *vrelolohi_h,
   double *vrehihilo_h, double *vrelohilo_h,
   double *vrehilolo_h, double *vrelololo_h,
   double *vimhihihi_h, double *vimlohihi_h,
   double *vimhilohi_h, double *vimlolohi_h,
   double *vimhihilo_h, double *vimlohilo_h,
   double *vimhilolo_h, double *vimlololo_h,
   double *Vrehihihi_d, double *Vrelohihi_d,
   double *Vrehilohi_d, double *Vrelolohi_d,
   double *Vrehihilo_d, double *Vrelohilo_d,
   double *Vrehilolo_d, double *Vrelololo_d,
   double *Vimhihihi_d, double *Vimlohihi_d,
   double *Vimhilohi_d, double *Vimlolohi_d,
   double *Vimhihilo_d, double *Vimlohilo_d,
   double *Vimhilolo_d, double *Vimlololo_d,
   double *betahihihi_h, double *betalohihi_h,
   double *betahilohi_h, double *betalolohi_h,
   double *betahihilo_h, double *betalohilo_h,
   double *betahilolo_h, double *betalololo_h,
   double *betahihihi_d, double *betalohihi_d,
   double *betahilohi_d, double *betalolohi_d,
   double *betahihilo_d, double *betalohilo_d,
   double *betahilolo_d, double *betalololo_d,
   double *lapms, long long int *add, long long int *mul, long long int *div,
   long long int *sqrtfun, bool verbose )
{
   const int nrLog2 = ceil(log2((double) nrows1));
   const int rowidx = colidx*(nrows+1);       // start of number in A_h
   const int nVrows = nrows - k*szt;          // dimension of V matrix

   if(verbose)
   {
      cout << "nrows : " << nrows
           << "  nVrows : " << nVrows
           << "  ncols : " << ncols
           << "  szt : " << szt
           << "  nbt : " << nbt << endl;
      cout << "k : " << k 
           << "  L : " << L
           << "  nrows1 : " << nrows1
           << "  colidx : " << colidx
           << "  rowidx : " << rowidx << endl;
   }
   if(L > 0)
   {
      for(int i=0; i<L; i++)   // insert zeros
      {
         vrehihihi_h[i] = 0.0; vrelohihi_h[i] = 0.0;
         vrehilohi_h[i] = 0.0; vrelolohi_h[i] = 0.0;
         vrehihilo_h[i] = 0.0; vrelohilo_h[i] = 0.0;
         vrehilolo_h[i] = 0.0; vrelololo_h[i] = 0.0;
         vimhihihi_h[i] = 0.0; vimlohihi_h[i] = 0.0;
         vimhilohi_h[i] = 0.0; vimlolohi_h[i] = 0.0;
         vimhihilo_h[i] = 0.0; vimlohilo_h[i] = 0.0;
         vimhilolo_h[i] = 0.0; vimlololo_h[i] = 0.0;
      }
      cudaMemcpy(&Vrehihihi_d[L*nVrows],vrehihihi_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vrelohihi_d[L*nVrows],vrelohihi_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vrehilohi_d[L*nVrows],vrehilohi_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vrelolohi_d[L*nVrows],vrelolohi_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vrehihilo_d[L*nVrows],vrehihilo_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vrelohilo_d[L*nVrows],vrelohilo_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vrehilolo_d[L*nVrows],vrehilolo_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vrelololo_d[L*nVrows],vrelololo_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vimhihihi_d[L*nVrows],vimhihihi_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vimlohihi_d[L*nVrows],vimlohihi_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vimhilohi_d[L*nVrows],vimhilohi_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vimlolohi_d[L*nVrows],vimlolohi_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vimhihilo_d[L*nVrows],vimhihilo_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vimlohilo_d[L*nVrows],vimlohilo_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vimhilolo_d[L*nVrows],vimhilolo_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vimlololo_d[L*nVrows],vimlololo_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
   }
   if(nrows1 == 0)
   {
      betahihihi_h[L] = 0.0; betalohihi_h[L] = 0.0;
      betahilohi_h[L] = 0.0; betalolohi_h[L] = 0.0;
      betahihilo_h[L] = 0.0; betalohilo_h[L] = 0.0;
      betahilolo_h[L] = 0.0; betalololo_h[L] = 0.0;
      vrehihihi_h[0] = 1.0; vrelohihi_h[0] = 0.0;
      vrehilohi_h[0] = 0.0; vrelolohi_h[0] = 0.0;
      vrehihilo_h[0] = 0.0; vrelohilo_h[0] = 0.0;
      vrehilolo_h[0] = 0.0; vrelololo_h[0] = 0.0;
      vimhihihi_h[0] = 0.0; vimlohihi_h[0] = 0.0;
      vimhilohi_h[0] = 0.0; vimlolohi_h[0] = 0.0;
      vimhihilo_h[0] = 0.0; vimlohilo_h[0] = 0.0;
      vimhilolo_h[0] = 0.0; vimlololo_h[0] = 0.0;
      cudaMemcpy(&betahihihi_d[L],&betahihihi_h[L],sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&betalohihi_d[L],&betalohihi_h[L],sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&betahilohi_d[L],&betahilohi_h[L],sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&betalolohi_d[L],&betalolohi_h[L],sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&betahihilo_d[L],&betahihilo_h[L],sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&betalohilo_d[L],&betalohilo_h[L],sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&betahilolo_d[L],&betahilolo_h[L],sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&betalololo_d[L],&betalololo_h[L],sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vrehihihi_d[L*nVrows+L],vrehihihi_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vrelohihi_d[L*nVrows+L],vrelohihi_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vrehilohi_d[L*nVrows+L],vrehilohi_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vrelolohi_d[L*nVrows+L],vrelolohi_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vrehihilo_d[L*nVrows+L],vrehihilo_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vrelohilo_d[L*nVrows+L],vrelohilo_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vrehilolo_d[L*nVrows+L],vrehilolo_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vrelololo_d[L*nVrows+L],vrelololo_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vimhihihi_d[L*nVrows+L],vimhihihi_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vimlohihi_d[L*nVrows+L],vimlohihi_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vimhilohi_d[L*nVrows+L],vimhilohi_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vimlolohi_d[L*nVrows+L],vimlolohi_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vimhihilo_d[L*nVrows+L],vimhihilo_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vimlohilo_d[L*nVrows+L],vimlohilo_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vimhilolo_d[L*nVrows+L],vimhilolo_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vimlololo_d[L*nVrows+L],vimlololo_h,sizeof(double),
                 cudaMemcpyHostToDevice);
   }
   else
   {
      if(verbose)  // to verify the input column is correct ...
      {
         cout << "The column of A : " << endl;
         cudaMemcpy(&Arehihihi_h[rowidx],&Arehihihi_d[rowidx],
                   (nrows1+1)*sizeof(double),cudaMemcpyDeviceToHost);
         cudaMemcpy(&Arelohihi_h[rowidx],&Arelohihi_d[rowidx],
                   (nrows1+1)*sizeof(double),cudaMemcpyDeviceToHost);
         cudaMemcpy(&Arehilohi_h[rowidx],&Arehilohi_d[rowidx],
                   (nrows1+1)*sizeof(double),cudaMemcpyDeviceToHost);
         cudaMemcpy(&Arelolohi_h[rowidx],&Arelolohi_d[rowidx],
                   (nrows1+1)*sizeof(double),cudaMemcpyDeviceToHost);
         cudaMemcpy(&Arehihilo_h[rowidx],&Arehihilo_d[rowidx],
                   (nrows1+1)*sizeof(double),cudaMemcpyDeviceToHost);
         cudaMemcpy(&Arelohilo_h[rowidx],&Arelohilo_d[rowidx],
                   (nrows1+1)*sizeof(double),cudaMemcpyDeviceToHost);
         cudaMemcpy(&Arehilolo_h[rowidx],&Arehilolo_d[rowidx],
                   (nrows1+1)*sizeof(double),cudaMemcpyDeviceToHost);
         cudaMemcpy(&Arelololo_h[rowidx],&Arelololo_d[rowidx],
                   (nrows1+1)*sizeof(double),cudaMemcpyDeviceToHost);
         cudaMemcpy(&Aimhihihi_h[rowidx],&Aimhihihi_d[rowidx],
                   (nrows1+1)*sizeof(double),cudaMemcpyDeviceToHost);
         cudaMemcpy(&Aimlohihi_h[rowidx],&Aimlohihi_d[rowidx],
                   (nrows1+1)*sizeof(double),cudaMemcpyDeviceToHost);
         cudaMemcpy(&Aimhilohi_h[rowidx],&Aimhilohi_d[rowidx],
                   (nrows1+1)*sizeof(double),cudaMemcpyDeviceToHost);
         cudaMemcpy(&Aimlolohi_h[rowidx],&Aimlolohi_d[rowidx],
                   (nrows1+1)*sizeof(double),cudaMemcpyDeviceToHost);
         cudaMemcpy(&Aimhihilo_h[rowidx],&Aimhihilo_d[rowidx],
                   (nrows1+1)*sizeof(double),cudaMemcpyDeviceToHost);
         cudaMemcpy(&Aimlohilo_h[rowidx],&Aimlohilo_d[rowidx],
                   (nrows1+1)*sizeof(double),cudaMemcpyDeviceToHost);
         cudaMemcpy(&Aimhilolo_h[rowidx],&Aimhilolo_d[rowidx],
                   (nrows1+1)*sizeof(double),cudaMemcpyDeviceToHost);
         cudaMemcpy(&Aimlololo_h[rowidx],&Aimlololo_d[rowidx],
                   (nrows1+1)*sizeof(double),cudaMemcpyDeviceToHost);

         for(int i=0; i<=nrows1; i++)
         {
            cout << "A[" << i << "]re : " 
                 << Arehihihi_h[i] << "  " << Arelohihi_h[i] << endl
                 << "         "
                 << Arehilohi_h[i] << "  " << Arelolohi_h[i] << endl
                 << "         "
                 << Arehihilo_h[i] << "  " << Arelohilo_h[i] << endl
                 << "         "
                 << Arehilolo_h[i] << "  " << Arelololo_h[i] << endl;
            cout << "A[" << i << "]im : " 
                 << Aimhihihi_h[i] << "  " << Aimlohihi_h[i] << endl
                 << "         "
                 << Aimhilohi_h[i] << "  " << Aimlolohi_h[i] << endl
                 << "         "
                 << Aimhihilo_h[i] << "  " << Aimlohilo_h[i] << endl
                 << "         "
                 << Aimhilolo_h[i] << "  " << Aimlololo_h[i] << endl;
         }
      }
      cudaEvent_t start,stop;           // to measure time spent by kernels 
      cudaEventCreate(&start);
      cudaEventCreate(&stop);
      float milliseconds;

      cudaEventRecord(start);
      cmplx8_small_house<<<1,nrows1>>>
         (&Arehihihi_d[rowidx],&Arelohihi_d[rowidx],
          &Arehilohi_d[rowidx],&Arelolohi_d[rowidx],
          &Arehihilo_d[rowidx],&Arelohilo_d[rowidx],
          &Arehilolo_d[rowidx],&Arelololo_d[rowidx],
          &Aimhihihi_d[rowidx],&Aimlohihi_d[rowidx],
          &Aimhilohi_d[rowidx],&Aimlolohi_d[rowidx],
          &Aimhihilo_d[rowidx],&Aimlohilo_d[rowidx],
          &Aimhilolo_d[rowidx],&Aimlololo_d[rowidx],
          &Arehihihi_d[rowidx+1],&Arelohihi_d[rowidx+1],
          &Arehilohi_d[rowidx+1],&Arelolohi_d[rowidx+1],
          &Arehihilo_d[rowidx+1],&Arelohilo_d[rowidx+1],
          &Arehilolo_d[rowidx+1],&Arelololo_d[rowidx+1],
          &Aimhihihi_d[rowidx+1],&Aimlohihi_d[rowidx+1],
          &Aimhilohi_d[rowidx+1],&Aimlolohi_d[rowidx+1],
          &Aimhihilo_d[rowidx+1],&Aimlohilo_d[rowidx+1],
          &Aimhilolo_d[rowidx+1],&Aimlololo_d[rowidx+1],nrows1,nrLog2,
          &Vrehihihi_d[L*nVrows+L],&Vrelohihi_d[L*nVrows+L],
          &Vrehilohi_d[L*nVrows+L],&Vrelolohi_d[L*nVrows+L],
          &Vrehihilo_d[L*nVrows+L],&Vrelohilo_d[L*nVrows+L],
          &Vrehilolo_d[L*nVrows+L],&Vrelololo_d[L*nVrows+L],
          &Vimhihihi_d[L*nVrows+L],&Vimlohihi_d[L*nVrows+L],
          &Vimhilohi_d[L*nVrows+L],&Vimlolohi_d[L*nVrows+L],
          &Vimhihilo_d[L*nVrows+L],&Vimlohilo_d[L*nVrows+L],
          &Vimhilolo_d[L*nVrows+L],&Vimlololo_d[L*nVrows+L],
          &betahihihi_d[L],&betalohihi_d[L],&betahilohi_d[L],&betalolohi_d[L],
          &betahihilo_d[L],&betalohilo_d[L],&betahilolo_d[L],&betalololo_d[L]);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds,start,stop);
      *lapms += milliseconds;
      flopcount_cmplx_small_house(nrows1,nrLog2,add,mul,div,sqrtfun);
   }
   cudaMemcpy(&betahihihi_h[L],&betahihihi_d[L],sizeof(double),
              cudaMemcpyDeviceToHost);
   cudaMemcpy(&betalohihi_h[L],&betalohihi_d[L],sizeof(double),
              cudaMemcpyDeviceToHost);
   cudaMemcpy(&betahilohi_h[L],&betahilohi_d[L],sizeof(double),
              cudaMemcpyDeviceToHost);
   cudaMemcpy(&betalolohi_h[L],&betalolohi_d[L],sizeof(double),
              cudaMemcpyDeviceToHost);
   cudaMemcpy(&betahihilo_h[L],&betahihilo_d[L],sizeof(double),
              cudaMemcpyDeviceToHost);
   cudaMemcpy(&betalohilo_h[L],&betalohilo_d[L],sizeof(double),
              cudaMemcpyDeviceToHost);
   cudaMemcpy(&betahilolo_h[L],&betahilolo_d[L],sizeof(double),
              cudaMemcpyDeviceToHost);
   cudaMemcpy(&betalololo_h[L],&betalololo_d[L],sizeof(double),
              cudaMemcpyDeviceToHost);
   if(verbose)
   {
      const size_t szhouse = nVrows*sizeof(double);

      cudaMemcpy(vrehihihi_h,&Vrehihihi_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vrelohihi_h,&Vrelohihi_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vrehilohi_h,&Vrehilohi_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vrelolohi_h,&Vrelolohi_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vrehihilo_h,&Vrehihilo_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vrelohilo_h,&Vrelohilo_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vrehilolo_h,&Vrehilolo_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vrelololo_h,&Vrelololo_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vimhihihi_h,&Vimhihihi_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vimlohihi_h,&Vimlohihi_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vimhilohi_h,&Vimhilohi_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vimlolohi_h,&Vimlolohi_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vimhihilo_h,&Vimhihilo_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vimlohilo_h,&Vimlohilo_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vimhilolo_h,&Vimhilolo_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vimlololo_h,&Vimlololo_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);

      cout << scientific << setprecision(16)
           << "beta[" << colidx << "] : "
           << betahihihi_h[L] << "  " << betalohihi_h[L] << endl
           << "          "
           << betahilohi_h[L] << "  " << betalolohi_h[L] << endl
           << "          "
           << betahihilo_h[L] << "  " << betalohilo_h[L] << endl
           << "          "
           << betahilolo_h[L] << "  " << betalololo_h[L] << endl;
      for(int i=0; i<nVrows; i++)
      {
         cout << "v[" << i << "]re : "
              << vrehihihi_h[i] << "  " << vrelohihi_h[i] << endl
              << "          "
              << vrehilohi_h[i] << "  " << vrelolohi_h[i] << endl
              << "          "
              << vrehihilo_h[i] << "  " << vrelohilo_h[i] << endl
              << "          "
              << vrehilolo_h[i] << "  " << vrelololo_h[i] << endl;
         cout << "v[" << i << "]im : "
              << vimhihihi_h[i] << "  " << vimlohihi_h[i] << endl
              << "          "
              << vimhilohi_h[i] << "  " << vimlolohi_h[i] << endl
              << "          "
              << vimhihilo_h[i] << "  " << vimlohilo_h[i] << endl
              << "          "
              << vimhilolo_h[i] << "  " << vimlololo_h[i] << endl;
      }
   }
}

void GPU_dbl8_large_house
 ( int nrows, int ncols, int szt, int nbt,
   int colidx, int nrows1, int k, int L,
   double *Ahihihi_h, double *Alohihi_h, double *Ahilohi_h, double *Alolohi_h,
   double *Ahihilo_h, double *Alohilo_h, double *Ahilolo_h, double *Alololo_h,
   double *Ahihihi_d, double *Alohihi_d, double *Ahilohi_d, double *Alolohi_d,
   double *Ahihilo_d, double *Alohilo_d, double *Ahilolo_d, double *Alololo_d,
   double *vhihihi_h, double *vlohihi_h, double *vhilohi_h, double *vlolohi_h,
   double *vhihilo_h, double *vlohilo_h, double *vhilolo_h, double *vlololo_h,
   double *Vhihihi_d, double *Vlohihi_d, double *Vhilohi_d, double *Vlolohi_d,
   double *Vhihilo_d, double *Vlohilo_d, double *Vhilolo_d, double *Vlololo_d,
   double *betahihihi_h, double *betalohihi_h,
   double *betahilohi_h, double *betalolohi_h,
   double *betahihilo_h, double *betalohilo_h,
   double *betahilolo_h, double *betalololo_h,
   double *betahihihi_d, double *betalohihi_d,
   double *betahilohi_d, double *betalolohi_d,
   double *betahihilo_d, double *betalohilo_d,
   double *betahilolo_d, double *betalololo_d,
   double *sumshihihi_h, double *sumslohihi_h,
   double *sumshilohi_h, double *sumslolohi_h,
   double *sumshihilo_h, double *sumslohilo_h,
   double *sumshilolo_h, double *sumslololo_h,
   double *sumshihihi_d, double *sumslohihi_d,
   double *sumshilohi_d, double *sumslolohi_d,
   double *sumshihilo_d, double *sumslohilo_d,
   double *sumshilolo_d, double *sumslololo_d,
   double *sigmahihihi_h, double *sigmalohihi_h, 
   double *sigmahilohi_h, double *sigmalolohi_h, 
   double *sigmahihilo_h, double *sigmalohilo_h, 
   double *sigmahilolo_h, double *sigmalololo_h, 
   double *sigmahihihi_d, double *sigmalohihi_d,
   double *sigmahilohi_d, double *sigmalolohi_d,
   double *sigmahihilo_d, double *sigmalohilo_d,
   double *sigmahilolo_d, double *sigmalololo_d,
   double *lapms, long long int *add, long long int *mul, long long int *div,
   long long int *sqrtfun, bool verbose )
{
   // nrows1 = nrows - colidx - 1 = size of Householder vector
   const int nblocks = ceil(((double) nrows1)/szt); // sufficient threads
   const int nblLog2 = ceil(log2((double) nblocks));
   const int sztLog2 = ceil(log2((double) szt));
   const int rowidx = colidx*(nrows+1);         // start of number in A_h
   const int nVrows = nrows - k*szt;             // dimension of V matrix

   cudaEvent_t start,stop;            // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;

   if(L > 0)
   {
      for(int i=0; i<L; i++)             // insert zeros
      {
         vhihihi_h[i] = 0.0;
         vlohihi_h[i] = 0.0;
         vhilohi_h[i] = 0.0;
         vlolohi_h[i] = 0.0;
         vhihilo_h[i] = 0.0;
         vlohilo_h[i] = 0.0;
         vhilolo_h[i] = 0.0;
         vlololo_h[i] = 0.0;
      }
   }
   vhihihi_h[L] = 1.0;                    // set one on the diagonal
   vlohihi_h[L] = 0.0;
   vhilohi_h[L] = 0.0;
   vlolohi_h[L] = 0.0;
   vhihilo_h[L] = 0.0;
   vlohilo_h[L] = 0.0;
   vhilolo_h[L] = 0.0;
   vlololo_h[L] = 0.0;

   cudaMemcpy(&Vhihihi_d[L*nVrows],vhihihi_h,(L+1)*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&Vlohihi_d[L*nVrows],vlohihi_h,(L+1)*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&Vhilohi_d[L*nVrows],vhilohi_h,(L+1)*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&Vlolohi_d[L*nVrows],vlolohi_h,(L+1)*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&Vhihilo_d[L*nVrows],vhihilo_h,(L+1)*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&Vlohilo_d[L*nVrows],vlohilo_h,(L+1)*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&Vhilolo_d[L*nVrows],vhilolo_h,(L+1)*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&Vlololo_d[L*nVrows],vlololo_h,(L+1)*sizeof(double),
              cudaMemcpyHostToDevice);

   if(verbose)
   {
      cout << "-> launching " << nblocks << " blocks of "
           << szt << " threads to compute the sum of squares ..." << endl;
      cout << "nrows1 : " << nrows1 << "  rowidx : " << rowidx;
      cout << "  ceil(log2(#blocks)) : " << nblLog2;
      cout << "  ceil(log2(szt)) : " << sztLog2 << endl;
   }
   for(int i=0; i<nblocks; i++)
   {
      sumshihihi_h[i] = 0.0;
      sumslohihi_h[i] = 0.0;
      sumshilohi_h[i] = 0.0;
      sumslolohi_h[i] = 0.0;
      sumshihilo_h[i] = 0.0;
      sumslohilo_h[i] = 0.0;
      sumshilolo_h[i] = 0.0;
      sumslololo_h[i] = 0.0;
   }
   cudaMemcpy(sumshihihi_d,sumshihihi_h,nblocks*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(sumslohihi_d,sumslohihi_h,nblocks*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(sumshilohi_d,sumshilohi_h,nblocks*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(sumslolohi_d,sumslolohi_h,nblocks*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(sumshihilo_d,sumshihilo_h,nblocks*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(sumslohilo_d,sumslohilo_h,nblocks*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(sumshilolo_d,sumshilolo_h,nblocks*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(sumslololo_d,sumslololo_h,nblocks*sizeof(double),
              cudaMemcpyHostToDevice);

   cudaEventRecord(start);
   dbl8_large_sum_of_squares<<<nblocks,szt>>>
      (&Ahihihi_d[rowidx+1],&Alohihi_d[rowidx+1],
       &Ahilohi_d[rowidx+1],&Alolohi_d[rowidx+1],
       &Ahihilo_d[rowidx+1],&Alohilo_d[rowidx+1],
       &Ahilolo_d[rowidx+1],&Alololo_d[rowidx+1],
       sumshihihi_d,sumslohihi_d,sumshilohi_d,sumslolohi_d,
       sumshihilo_d,sumslohilo_d,sumshilolo_d,sumslololo_d,
       nrows1,szt,sztLog2);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_dbl_large_sum_of_squares(nblocks,szt,sztLog2,add,mul);

   if(verbose)
   {
      cout << "-> launching 1 block of " << nblocks
           << " threads to accumulate the sums ..." << endl;
   }
   cudaEventRecord(start);
   dbl8_sum_accumulator<<<1,nblocks>>>
      (sumshihihi_d,sumslohihi_d,sumshilohi_d,sumslolohi_d,
       sumshihilo_d,sumslohilo_d,sumshilolo_d,sumslololo_d,
       nblocks,nblLog2,
       sigmahihihi_d,sigmalohihi_d,sigmahilohi_d,sigmalolohi_d,
       sigmahihilo_d,sigmalohilo_d,sigmahilolo_d,sigmalololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_dbl_sum_accumulator(nblocks,nblLog2,add);

   cudaMemcpy(sigmahihihi_h,sigmahihihi_d,sizeof(double),
              cudaMemcpyDeviceToHost);
   cudaMemcpy(sigmalohihi_h,sigmalohihi_d,sizeof(double),
              cudaMemcpyDeviceToHost);
   cudaMemcpy(sigmahilohi_h,sigmahilohi_d,sizeof(double),
              cudaMemcpyDeviceToHost);
   cudaMemcpy(sigmalolohi_h,sigmalolohi_d,sizeof(double),
              cudaMemcpyDeviceToHost);
   cudaMemcpy(sigmahihilo_h,sigmahihilo_d,sizeof(double),
              cudaMemcpyDeviceToHost);
   cudaMemcpy(sigmalohilo_h,sigmalohilo_d,sizeof(double),
              cudaMemcpyDeviceToHost);
   cudaMemcpy(sigmahilolo_h,sigmahilolo_d,sizeof(double),
              cudaMemcpyDeviceToHost);
   cudaMemcpy(sigmalololo_h,sigmalololo_d,sizeof(double),
              cudaMemcpyDeviceToHost);

   bool done = false;

   if((sigmahihihi_h[0] == 0.0) && (sigmalohihi_h[0] == 0.0) &&
      (sigmahilohi_h[0] == 0.0) && (sigmalolohi_h[0] == 0.0) &&
      (sigmahihilo_h[0] == 0.0) && (sigmalohilo_h[0] == 0.0) &&
      (sigmahilolo_h[0] == 0.0) && (sigmalololo_h[0] == 0.0))
   {
      betahihihi_h[L] = 0.0; betalohihi_h[L] = 0.0;
      betahilohi_h[L] = 0.0; betalolohi_h[L] = 0.0;
      betahihilo_h[L] = 0.0; betalohilo_h[L] = 0.0;
      betahilolo_h[L] = 0.0; betalololo_h[L] = 0.0;
      done = true;

      if(verbose)
         cout << "Zero sigma value encountered." << endl;
   }
   else // beta is computed on the host instead of by one GPU thread
   {
      // const double x0hi = Ahi_h[rowidx];
      // const double x0lo = Alo_h[rowidx];
      double acchihihi,acclohihi,acchilohi,acclolohi;
      double acchihilo,acclohilo,acchilolo,acclololo;
      double muhihihi,mulohihi,muhilohi,mulolohi;
      double muhihilo,mulohilo,muhilolo,mulololo;
      double v0hihihi,v0lohihi,v0hilohi,v0lolohi;
      double v0hihilo,v0lohilo,v0hilolo,v0lololo;
      double v0p2hihihi,v0p2lohihi,v0p2hilohi,v0p2lolohi;
      double v0p2hihilo,v0p2lohilo,v0p2hilolo,v0p2lololo;
      double x0hihihi,x0lohihi,x0hilohi,x0lolohi;
      double x0hihilo,x0lohilo,x0hilolo,x0lololo;

      cudaMemcpy(&x0hihihi,&Ahihihi_d[rowidx],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&x0lohihi,&Alohihi_d[rowidx],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&x0hilohi,&Ahilohi_d[rowidx],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&x0lolohi,&Alolohi_d[rowidx],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&x0hihilo,&Ahihilo_d[rowidx],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&x0lohilo,&Alohilo_d[rowidx],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&x0hilolo,&Ahilolo_d[rowidx],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&x0lololo,&Alololo_d[rowidx],sizeof(double),
                 cudaMemcpyDeviceToHost);

      // mu = sqrt((*x0)*(*x0) + sigma[0]);
      odf_sqr(x0hihihi,x0lohihi,x0hilohi,x0lolohi,
              x0hihilo,x0lohilo,x0hilolo,x0lololo,
              &acchihihi,&acclohihi,&acchilohi,&acclolohi,
              &acchihilo,&acclohilo,&acchilolo,&acclololo);
      odf_inc(&acchihihi,&acclohihi,&acchilohi,&acclolohi,
              &acchihilo,&acclohilo,&acchilolo,&acclololo,
              sigmahihihi_h[0],sigmalohihi_h[0],
              sigmahilohi_h[0],sigmalolohi_h[0],
              sigmahihilo_h[0],sigmalohilo_h[0],
              sigmahilolo_h[0],sigmalololo_h[0]);
      odf_sqrt(acchihihi,acclohihi,acchilohi,acclolohi,
               acchihilo,acclohilo,acchilolo,acclololo,
               &muhihihi,&mulohihi,&muhilohi,&mulolohi,
               &muhihilo,&mulohilo,&muhilolo,&mulololo);
      if(x0hihihi <= 0.0)
      {
         // v0 = *x0 - mu;
         odf_sub( x0hihihi, x0lohihi, x0hilohi, x0lolohi,
                  x0hihilo, x0lohilo, x0hilolo, x0lololo,
                  muhihihi, mulohihi, muhilohi, mulolohi,
                  muhihilo, mulohilo, muhilolo, mulololo,
                 &v0hihihi,&v0lohihi,&v0hilohi,&v0lolohi,
                 &v0hihilo,&v0lohilo,&v0hilolo,&v0lololo);
      }
      else
      {
         // v0 = -sigma[0]/(*x0 + mu);
         odf_add(  x0hihihi,  x0lohihi,  x0hilohi,  x0lolohi,
                   x0hihilo,  x0lohilo,  x0hilolo,  x0lololo,
                   muhihihi,  mulohihi,  muhilohi,  mulolohi,
                   muhihilo,  mulohilo,  muhilolo,  mulololo,
                 &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                 &acchihilo,&acclohilo,&acchilolo,&acclololo);
         odf_div(sigmahihihi_h[0],sigmalohihi_h[0],
                 sigmahilohi_h[0],sigmalolohi_h[0],
                 sigmahihilo_h[0],sigmalohilo_h[0],
                 sigmahilolo_h[0],sigmalololo_h[0],
                   acchihihi,       acclohihi,     acchilohi,     acclolohi,
                   acchihilo,       acclohilo,     acchilolo,     acclololo,
                   &v0hihihi,       &v0lohihi,     &v0hilohi,     &v0lolohi,
                   &v0hihilo,       &v0lohilo,     &v0hilolo,     &v0lololo);
         odf_minus(&v0hihihi,&v0lohihi,&v0hilohi,&v0lolohi,
                   &v0hihilo,&v0lohilo,&v0hilolo,&v0lololo);
      }
      // v0p2 = v0*v0;
      odf_sqr(   v0hihihi,   v0lohihi,   v0hilohi,   v0lolohi,
                 v0hihilo,   v0lohilo,   v0hilolo,   v0lololo,
              &v0p2hihihi,&v0p2lohihi,&v0p2hilohi,&v0p2lolohi,
              &v0p2hihilo,&v0p2lohilo,&v0p2hilolo,&v0p2lololo);
      // *beta = 2.0*v0p2/(sigma[0] + v0p2);
      odf_add(sigmahihihi_h[0],sigmalohihi_h[0],
              sigmahilohi_h[0],sigmalolohi_h[0],
              sigmahihilo_h[0],sigmalohilo_h[0],
              sigmahilolo_h[0],sigmalololo_h[0],
               v0p2hihihi,      v0p2lohihi,      v0p2hilohi,      v0p2lolohi,
               v0p2hihilo,      v0p2lohilo,      v0p2hilolo,      v0p2lololo,
               &acchihihi,      &acclohihi,      &acchilohi,      &acclolohi,
               &acchihilo,      &acclohilo,      &acchilolo,      &acclololo);
      odf_div( v0p2hihihi,      v0p2lohihi,      v0p2hilohi,      v0p2lolohi,
               v0p2hihilo,      v0p2lohilo,      v0p2hilolo,      v0p2lololo,
                acchihihi,       acclohihi,       acchilohi,       acclolohi,
                acchihilo,       acclohilo,       acchilolo,       acclololo,
              &betahihihi_h[L],&betalohihi_h[L],
              &betahilohi_h[L],&betalolohi_h[L],
              &betahihilo_h[L],&betalohilo_h[L],
              &betahilolo_h[L],&betalololo_h[L]);
      odf_mlt_d(&betahihihi_h[L],&betalohihi_h[L],
                &betahilohi_h[L],&betalolohi_h[L],
                &betahihilo_h[L],&betalohilo_h[L],
                &betahilolo_h[L],&betalololo_h[L],2.0);
      sigmahihihi_h[0] = v0hihihi;
      sigmalohihi_h[0] = v0lohihi;
      sigmahilohi_h[0] = v0hilohi;
      sigmalolohi_h[0] = v0lolohi;
      sigmahihilo_h[0] = v0hihilo;
      sigmalohilo_h[0] = v0lohilo;
      sigmahilolo_h[0] = v0hilolo;
      sigmalololo_h[0] = v0lololo;           // v0 needed for normalization
      // update the flop counts
      *add += 3;
      *mul += 3;
      *div += 2;
      *sqrtfun += 1;
   }
   if(verbose)
   {
      cout << scientific << setprecision(16)
           << "beta[" << colidx << "] : "
           << betahihihi_h[L] << "  " << betalohihi_h[L] << endl
           << "          "
           << betahilohi_h[L] << "  " << betalolohi_h[L] << endl
           << "          "
           << betahihilo_h[L] << "  " << betalohilo_h[L] << endl
           << "          "
           << betahilolo_h[L] << "  " << betalololo_h[L] << endl;
   }
   cudaMemcpy(&betahihihi_d[L],&betahihihi_h[L],sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&betalohihi_d[L],&betalohihi_h[L],sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&betahilohi_d[L],&betahilohi_h[L],sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&betalolohi_d[L],&betalolohi_h[L],sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&betahihilo_d[L],&betahihilo_h[L],sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&betalohilo_d[L],&betalohilo_h[L],sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&betahilolo_d[L],&betahilolo_h[L],sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&betalololo_d[L],&betalololo_h[L],sizeof(double),
              cudaMemcpyHostToDevice);

   if(!done)  // normalization needed
   {
      // (sigmahi_h, sigmalo_h) has the values for (v0hi, v0lo).
      cudaMemcpy(sigmahihihi_d,sigmahihihi_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(sigmalohihi_d,sigmalohihi_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(sigmahilohi_d,sigmahilohi_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(sigmalolohi_d,sigmalolohi_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(sigmahihilo_d,sigmahihilo_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(sigmalohilo_d,sigmalohilo_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(sigmahilolo_d,sigmahilolo_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(sigmalololo_d,sigmalololo_h,sizeof(double),
                 cudaMemcpyHostToDevice);

      if(verbose)
      {
         cout << "-> launching " << nblocks << " blocks of "
              << szt << " threads to normalize ..." << endl;
         cout << "   nrows1 : " << nrows1
              << "  rowidx : " << rowidx << "  nVrows : " << nVrows << endl;
      }
      cudaEventRecord(start);
      dbl8_normalize<<<nblocks,szt>>>
         (nrows1,szt,
          &Ahihihi_d[rowidx+1],&Alohihi_d[rowidx+1],
          &Ahilohi_d[rowidx+1],&Alolohi_d[rowidx+1],
          &Ahihilo_d[rowidx+1],&Alohilo_d[rowidx+1],
          &Ahilolo_d[rowidx+1],&Alololo_d[rowidx+1],
          sigmahihihi_d,sigmalohihi_d,sigmahilohi_d,sigmalolohi_d,
          sigmahihilo_d,sigmalohilo_d,sigmahilolo_d,sigmalololo_d,
          &Vhihihi_d[L*nVrows+L+1],&Vlohihi_d[L*nVrows+L+1],
          &Vhilohi_d[L*nVrows+L+1],&Vlolohi_d[L*nVrows+L+1],
          &Vhihilo_d[L*nVrows+L+1],&Vlohilo_d[L*nVrows+L+1],
          &Vhilolo_d[L*nVrows+L+1],&Vlololo_d[L*nVrows+L+1]);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds,start,stop);
      *lapms += milliseconds;
      flopcount_dbl_normalize(nblocks,szt,div);
   }
   if(verbose)
   {
      const size_t szhouse = nVrows*sizeof(double);

      cudaMemcpy(&betahihihi_h[L],&betahihihi_d[L],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&betalohihi_h[L],&betalohihi_d[L],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&betahilohi_h[L],&betahilohi_d[L],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&betalolohi_h[L],&betalolohi_d[L],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&betahihilo_h[L],&betahihilo_d[L],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&betalohilo_h[L],&betalohilo_d[L],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&betahilolo_h[L],&betahilolo_d[L],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&betalololo_h[L],&betalololo_d[L],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vhihihi_h,&Vhihihi_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vlohihi_h,&Vlohihi_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vhilohi_h,&Vhilohi_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vlolohi_h,&Vlolohi_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vhihilo_h,&Vhihilo_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vlohilo_h,&Vlohilo_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vhilolo_h,&Vhilolo_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vlololo_h,&Vlololo_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);

      cout << scientific << setprecision(16)
           << "beta[" << colidx << "] : "
           << betahihihi_h[L] << "  " << betalohihi_h[L] << endl
           << "           "
           << betahilohi_h[L] << "  " << betalolohi_h[L] << endl
           << "           "
           << betahihilo_h[L] << "  " << betalohilo_h[L] << endl
           << "           "
           << betahilolo_h[L] << "  " << betalololo_h[L] << endl;
      for(int i=0; i<nVrows; i++)
         cout << "v[" << i << "] : "
              << vhihihi_h[i] << "  " << vlohihi_h[i] << endl
              << "       "
              << vhilohi_h[i] << "  " << vlolohi_h[i] << endl
              << "       "
              << vhihilo_h[i] << "  " << vlohilo_h[i] << endl
              << "       "
              << vhilolo_h[i] << "  " << vlololo_h[i] << endl;
   }
}

void GPU_cmplx8_large_house
 ( int nrows, int ncols, int szt, int nbt,
   int colidx, int nrows1, int k, int L,
   double *Arehihihi_h, double *Arelohihi_h,
   double *Arehilohi_h, double *Arelolohi_h,
   double *Arehihilo_h, double *Arelohilo_h,
   double *Arehilolo_h, double *Arelololo_h,
   double *Aimhihihi_h, double *Aimlohihi_h,
   double *Aimhilohi_h, double *Aimlolohi_h,
   double *Aimhihilo_h, double *Aimlohilo_h,
   double *Aimhilolo_h, double *Aimlololo_h,
   double *Arehihihi_d, double *Arelohihi_d,
   double *Arehilohi_d, double *Arelolohi_d,
   double *Arehihilo_d, double *Arelohilo_d,
   double *Arehilolo_d, double *Arelololo_d,
   double *Aimhihihi_d, double *Aimlohihi_d,
   double *Aimhilohi_d, double *Aimlolohi_d,
   double *Aimhihilo_d, double *Aimlohilo_d,
   double *Aimhilolo_d, double *Aimlololo_d,
   double *vrehihihi_h, double *vrelohihi_h,
   double *vrehilohi_h, double *vrelolohi_h,
   double *vrehihilo_h, double *vrelohilo_h,
   double *vrehilolo_h, double *vrelololo_h,
   double *vimhihihi_h, double *vimlohihi_h,
   double *vimhilohi_h, double *vimlolohi_h,
   double *vimhihilo_h, double *vimlohilo_h,
   double *vimhilolo_h, double *vimlololo_h,
   double *Vrehihihi_d, double *Vrelohihi_d,
   double *Vrehilohi_d, double *Vrelolohi_d,
   double *Vrehihilo_d, double *Vrelohilo_d,
   double *Vrehilolo_d, double *Vrelololo_d,
   double *Vimhihihi_d, double *Vimlohihi_d,
   double *Vimhilohi_d, double *Vimlolohi_d,
   double *Vimhihilo_d, double *Vimlohilo_d,
   double *Vimhilolo_d, double *Vimlololo_d,
   double *betahihihi_h, double *betalohihi_h,
   double *betahilohi_h, double *betalolohi_h,
   double *betahihilo_h, double *betalohilo_h,
   double *betahilolo_h, double *betalololo_h,
   double *betahihihi_d, double *betalohihi_d,
   double *betahilohi_d, double *betalolohi_d,
   double *betahihilo_d, double *betalohilo_d,
   double *betahilolo_d, double *betalololo_d,
   double *sumshihihi_h, double *sumslohihi_h,
   double *sumshilohi_h, double *sumslolohi_h,
   double *sumshihilo_h, double *sumslohilo_h,
   double *sumshilolo_h, double *sumslololo_h,
   double *sumshihihi_d, double *sumslohihi_d,
   double *sumshilohi_d, double *sumslolohi_d,
   double *sumshihilo_d, double *sumslohilo_d,
   double *sumshilolo_d, double *sumslololo_d,
   double *sigmahihihi_h, double *sigmalohihi_h,
   double *sigmahilohi_h, double *sigmalolohi_h,
   double *sigmahihilo_h, double *sigmalohilo_h,
   double *sigmahilolo_h, double *sigmalololo_h,
   double *sigmahihihi_d, double *sigmalohihi_d,
   double *sigmahilohi_d, double *sigmalolohi_d,
   double *sigmahihilo_d, double *sigmalohilo_d,
   double *sigmahilolo_d, double *sigmalololo_d,
   double *lapms, long long int *add, long long int *mul, long long int *div,
   long long int *sqrtfun, bool verbose )
{
   // nrows1 = nrows - colidx - 1 = size of Householder vector
   const int nblocks = ceil(((double) nrows1)/szt); // sufficient threads
   const int nblLog2 = ceil(log2((double) nblocks));
   const int sztLog2 = ceil(log2((double) szt));
   const int rowidx = colidx*(nrows+1);         // start of number in A_h
   const int nVrows = nrows - k*szt;             // dimension of V matrix

   cudaEvent_t start,stop;            // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;

   if(L > 0)
   {
      for(int i=0; i<L; i++)             // insert zeros
      {
         vrehihihi_h[i] = 0.0; vrelohihi_h[i] = 0.0;
         vrehilohi_h[i] = 0.0; vrelolohi_h[i] = 0.0;
         vrehihilo_h[i] = 0.0; vrelohilo_h[i] = 0.0;
         vrehilolo_h[i] = 0.0; vrelololo_h[i] = 0.0;
         vimhihihi_h[i] = 0.0; vimlohihi_h[i] = 0.0;
         vimhilohi_h[i] = 0.0; vimlolohi_h[i] = 0.0;
         vimhihilo_h[i] = 0.0; vimlohilo_h[i] = 0.0;
         vimhilolo_h[i] = 0.0; vimlololo_h[i] = 0.0;
      }
   }
   vrehihihi_h[L] = 1.0; vrelohihi_h[L] = 0.0; // set one on the diagonal
   vrehilohi_h[L] = 0.0; vrelolohi_h[L] = 0.0;
   vrehihilo_h[L] = 0.0; vrelohilo_h[L] = 0.0;
   vrehilolo_h[L] = 0.0; vrelololo_h[L] = 0.0;
   vimhihihi_h[L] = 0.0; vimlohihi_h[L] = 0.0;
   vimhilohi_h[L] = 0.0; vimlolohi_h[L] = 0.0;
   vimhihilo_h[L] = 0.0; vimlohilo_h[L] = 0.0;
   vimhilolo_h[L] = 0.0; vimlololo_h[L] = 0.0;

   cudaMemcpy(&Vrehihihi_d[L*nVrows],vrehihihi_h,(L+1)*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&Vrelohihi_d[L*nVrows],vrelohihi_h,(L+1)*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&Vrehilohi_d[L*nVrows],vrehilohi_h,(L+1)*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&Vrelolohi_d[L*nVrows],vrelolohi_h,(L+1)*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&Vrehihilo_d[L*nVrows],vrehihilo_h,(L+1)*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&Vrelohilo_d[L*nVrows],vrelohilo_h,(L+1)*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&Vrehilolo_d[L*nVrows],vrehilolo_h,(L+1)*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&Vrelololo_d[L*nVrows],vrelololo_h,(L+1)*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&Vimhihihi_d[L*nVrows],vimhihihi_h,(L+1)*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&Vimlohihi_d[L*nVrows],vimlohihi_h,(L+1)*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&Vimhilohi_d[L*nVrows],vimhilohi_h,(L+1)*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&Vimlolohi_d[L*nVrows],vimlolohi_h,(L+1)*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&Vimhihilo_d[L*nVrows],vimhihilo_h,(L+1)*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&Vimlohilo_d[L*nVrows],vimlohilo_h,(L+1)*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&Vimhilolo_d[L*nVrows],vimhilolo_h,(L+1)*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&Vimlololo_d[L*nVrows],vimlololo_h,(L+1)*sizeof(double),
              cudaMemcpyHostToDevice);

   if(verbose)
   {
      cout << "-> launching " << nblocks << " blocks of "
           << szt << " threads to compute the sum of squares ..." << endl;
      cout << "nrows1 : " << nrows1 << "  rowidx : " << rowidx;
      cout << "  ceil(log2(#blocks)) : " << nblLog2;
      cout << "  ceil(log2(szt)) : " << sztLog2 << endl;
   }
   for(int i=0; i<nblocks; i++)
   {
      sumshihihi_h[i] = 0.0; sumslohihi_h[i] = 0.0;
      sumshilohi_h[i] = 0.0; sumslolohi_h[i] = 0.0;
      sumshihilo_h[i] = 0.0; sumslohilo_h[i] = 0.0;
      sumshilolo_h[i] = 0.0; sumslololo_h[i] = 0.0;
   }
   cudaMemcpy(sumshihihi_d,sumshihihi_h,nblocks*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(sumslohihi_d,sumslohihi_h,nblocks*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(sumshilohi_d,sumshilohi_h,nblocks*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(sumslolohi_d,sumslolohi_h,nblocks*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(sumshihilo_d,sumshihilo_h,nblocks*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(sumslohilo_d,sumslohilo_h,nblocks*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(sumshilolo_d,sumshilolo_h,nblocks*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(sumslololo_d,sumslololo_h,nblocks*sizeof(double),
              cudaMemcpyHostToDevice);

   cudaEventRecord(start);
   cmplx8_large_sum_of_squares<<<nblocks,szt>>>
      (&Arehihihi_d[rowidx+1],&Arelohihi_d[rowidx+1],
       &Arehilohi_d[rowidx+1],&Arelolohi_d[rowidx+1],
       &Arehihilo_d[rowidx+1],&Arelohilo_d[rowidx+1],
       &Arehilolo_d[rowidx+1],&Arelololo_d[rowidx+1],
       &Aimhihihi_d[rowidx+1],&Aimlohihi_d[rowidx+1],
       &Aimhilohi_d[rowidx+1],&Aimlolohi_d[rowidx+1],
       &Aimhihilo_d[rowidx+1],&Aimlohilo_d[rowidx+1],
       &Aimhilolo_d[rowidx+1],&Aimlololo_d[rowidx+1],
       sumshihihi_d,sumslohihi_d,sumshilohi_d,sumslolohi_d,
       sumshihilo_d,sumslohilo_d,sumshilolo_d,sumslololo_d,
       nrows1,szt,sztLog2);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_cmplx_large_sum_of_squares(nblocks,szt,sztLog2,add,mul);

   if(verbose)
   {
      cout << "-> launching 1 block of " << nblocks
           << " threads to accumulate the sums ..." << endl;
   }
   cudaEventRecord(start);
   dbl8_sum_accumulator<<<1,nblocks>>>
      (sumshihihi_d,sumslohihi_d,sumshilohi_d,sumslolohi_d,
       sumshihilo_d,sumslohilo_d,sumshilolo_d,sumslololo_d,
       nblocks,nblLog2,
       sigmahihihi_d,sigmalohihi_d,sigmahilohi_d,sigmalolohi_d,
       sigmahihilo_d,sigmalohilo_d,sigmahilolo_d,sigmalololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_dbl_sum_accumulator(nblocks,nblLog2,add);

   cudaMemcpy(sigmahihihi_h,sigmahihihi_d,sizeof(double),
              cudaMemcpyDeviceToHost);
   cudaMemcpy(sigmalohihi_h,sigmalohihi_d,sizeof(double),
              cudaMemcpyDeviceToHost);
   cudaMemcpy(sigmahilohi_h,sigmahilohi_d,sizeof(double),
              cudaMemcpyDeviceToHost);
   cudaMemcpy(sigmalolohi_h,sigmalolohi_d,sizeof(double),
              cudaMemcpyDeviceToHost);
   cudaMemcpy(sigmahihilo_h,sigmahihilo_d,sizeof(double),
              cudaMemcpyDeviceToHost);
   cudaMemcpy(sigmalohilo_h,sigmalohilo_d,sizeof(double),
              cudaMemcpyDeviceToHost);
   cudaMemcpy(sigmahilolo_h,sigmahilolo_d,sizeof(double),
              cudaMemcpyDeviceToHost);
   cudaMemcpy(sigmalololo_h,sigmalololo_d,sizeof(double),
              cudaMemcpyDeviceToHost);

   bool done = false;

   if((sigmahihihi_h[0] == 0.0) && (sigmalohihi_h[0] == 0.0) &&
      (sigmahilohi_h[0] == 0.0) && (sigmalolohi_h[0] == 0.0) &&
      (sigmahihilo_h[0] == 0.0) && (sigmalohilo_h[0] == 0.0) &&
      (sigmahilolo_h[0] == 0.0) && (sigmalololo_h[0] == 0.0))
   {
      betahihihi_h[L] = 0.0; betalohihi_h[L] = 0.0;
      betahilohi_h[L] = 0.0; betalolohi_h[L] = 0.0;
      betahihilo_h[L] = 0.0; betalohilo_h[L] = 0.0;
      betahilolo_h[L] = 0.0; betalololo_h[L] = 0.0; done = true;
      if(verbose)
         cout << "Zero sigma value encountered." << endl;
   }
   else // beta is computed on the host instead of by one GPU thread
   {
      double acchihihi,acclohihi,acchilohi,acclolohi;
      double acchihilo,acclohilo,acchilolo,acclololo;
      double muhihihi,mulohihi,muhilohi,mulolohi;
      double muhihilo,mulohilo,muhilolo,mulololo;
      double sqrv0hihihi,sqrv0lohihi,sqrv0hilohi,sqrv0lolohi;
      double sqrv0hihilo,sqrv0lohilo,sqrv0hilolo,sqrv0lololo;
      double x0rehihihi,x0relohihi,x0rehilohi,x0relolohi;
      double x0rehihilo,x0relohilo,x0rehilolo,x0relololo;
      double x0imhihihi,x0imlohihi,x0imhilohi,x0imlolohi;
      double x0imhihilo,x0imlohilo,x0imhilolo,x0imlololo;
      double sqrx0hihihi,sqrx0lohihi,sqrx0hilohi,sqrx0lolohi;
      double sqrx0hihilo,sqrx0lohilo,sqrx0hilolo,sqrx0lololo;
      double v0rehihihi,v0relohihi,v0rehilohi,v0relolohi;
      double v0rehihilo,v0relohilo,v0rehilolo,v0relololo;
      double v0imhihihi,v0imlohihi,v0imhilohi,v0imlolohi;
      double v0imhihilo,v0imlohilo,v0imhilolo,v0imlololo;
      double x0radhihihi,x0radlohihi,x0radhilohi,x0radlolohi;
      double x0radhihilo,x0radlohilo,x0radhilolo,x0radlololo;
      double inv0rehihihi,inv0relohihi,inv0rehilohi,inv0relolohi;
      double inv0rehihilo,inv0relohilo,inv0rehilolo,inv0relololo;
      double inv0imhihihi,inv0imlohihi,inv0imhilohi,inv0imlolohi;
      double inv0imhihilo,inv0imlohilo,inv0imhilolo,inv0imlololo;

      cudaMemcpy(&x0rehihihi,&Arehihihi_d[rowidx],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&x0relohihi,&Arelohihi_d[rowidx],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&x0rehilohi,&Arehilohi_d[rowidx],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&x0relolohi,&Arelolohi_d[rowidx],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&x0rehihilo,&Arehihilo_d[rowidx],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&x0relohilo,&Arelohilo_d[rowidx],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&x0rehilolo,&Arehilolo_d[rowidx],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&x0relololo,&Arelololo_d[rowidx],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&x0imhihihi,&Aimhihihi_d[rowidx],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&x0imlohihi,&Aimlohihi_d[rowidx],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&x0imhilohi,&Aimhilohi_d[rowidx],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&x0imlolohi,&Aimlolohi_d[rowidx],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&x0imhihilo,&Aimhihilo_d[rowidx],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&x0imlohilo,&Aimlohilo_d[rowidx],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&x0imhilolo,&Aimhilolo_d[rowidx],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&x0imlololo,&Aimlololo_d[rowidx],sizeof(double),
                 cudaMemcpyDeviceToHost);

      // sqrx0 = xre[0]*xre[0] + xim[0]*xim[0];
      odf_sqr(x0rehihihi,x0relohihi,x0rehilohi,x0relolohi,
              x0rehihilo,x0relohilo,x0rehilolo,x0relololo,
              &sqrx0hihihi,&sqrx0lohihi,&sqrx0hilohi,&sqrx0lolohi,
              &sqrx0hihilo,&sqrx0lohilo,&sqrx0hilolo,&sqrx0lololo);
      odf_sqr(x0imhihihi,x0imlohihi,x0imhilohi,x0imlolohi,
              x0imhihilo,x0imlohilo,x0imhilolo,x0imlololo,
              &acchihihi,&acclohihi,&acchilohi,&acclolohi,
              &acchihilo,&acclohilo,&acchilolo,&acclololo);
      odf_inc(&sqrx0hihihi,&sqrx0lohihi,&sqrx0hilohi,&sqrx0lolohi,
              &sqrx0hihilo,&sqrx0lohilo,&sqrx0hilolo,&sqrx0lololo,
              acchihihi,acclohihi,acchilohi,acclolohi,
              acchihilo,acclohilo,acchilolo,acclololo);
      // x0rad = sqrt(sqrx0);
      odf_sqrt(sqrx0hihihi,sqrx0lohihi,sqrx0hilohi,sqrx0lolohi,
               sqrx0hihilo,sqrx0lohilo,sqrx0hilolo,sqrx0lololo,
               &x0radhihihi,&x0radlohihi,&x0radhilohi,&x0radlolohi,
               &x0radhihilo,&x0radlohilo,&x0radhilolo,&x0radlololo);
      // mu = sqrt(sqrx0 + sigma); // norm of the vector x
      odf_inc(&sqrx0hihihi,&sqrx0lohihi,&sqrx0hilohi,&sqrx0lolohi,
              &sqrx0hihilo,&sqrx0lohilo,&sqrx0hilolo,&sqrx0lololo,
              *sigmahihihi_h,*sigmalohihi_h,*sigmahilohi_h,*sigmalolohi_h,
              *sigmahihilo_h,*sigmalohilo_h,*sigmahilolo_h,*sigmalololo_h);
      odf_sqrt(sqrx0hihihi,sqrx0lohihi,sqrx0hilohi,sqrx0lolohi,
               sqrx0hihilo,sqrx0lohilo,sqrx0hilolo,sqrx0lololo,
               &muhihihi,&mulohihi,&muhilohi,&mulolohi,
               &muhihilo,&mulohilo,&muhilolo,&mulololo);

      if((x0radhihihi == 0.0) && (x0radlohihi == 0.0) &&
         (x0radhilohi == 0.0) && (x0radlolohi == 0.0) &&
         (x0radhihilo == 0.0) && (x0radlohilo == 0.0) &&
         (x0radhilolo == 0.0) && (x0radlololo == 0.0))
      {
         v0rehihihi = -muhihihi;
         v0relohihi = -mulohihi; 
         v0rehilohi = -muhilohi;
         v0relolohi = -mulolohi;
         v0rehihilo = -muhihilo;
         v0relohilo = -mulohilo;
         v0rehilolo = -muhilolo;
         v0relololo = -mulololo;
         v0imhihihi = 0.0;
         v0imlohihi = 0.0;
         v0imhilohi = 0.0;
         v0imlolohi = 0.0;
         v0imhihilo = 0.0;
         v0imlohilo = 0.0;
         v0imhilolo = 0.0;
         v0imlololo = 0.0;
      }
      else // if(x0rad /= 0.0)   // xre[0]/xrad = cos(angle)
      {                          // xim[0]/xrad = sin(angle)
         // mu = mu/x0rad;
         odf_div(muhihihi,mulohihi,muhilohi,mulolohi,
                 muhihilo,mulohilo,muhilolo,mulololo,
                 x0radhihihi,x0radlohihi,x0radhilohi,x0radlolohi,
                 x0radhihilo,x0radlohilo,x0radhilolo,x0radlololo,
                 &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                 &acchihilo,&acclohilo,&acchilolo,&acclololo);
         muhihihi = acchihihi; mulohihi = acclohihi;
         muhilohi = acchilohi; mulolohi = acclolohi;
         muhihilo = acchihilo; mulohilo = acclohilo;
         muhilolo = acchilolo; mulololo = acclololo;
         // vre[0] = xre[0] - mu*xre[0];
         odf_mul(muhihihi,mulohihi,muhilohi,mulolohi,
                 muhihilo,mulohilo,muhilolo,mulololo,
                 x0rehihihi,x0relohihi,x0rehilohi,x0relolohi,
                 x0rehihilo,x0relohilo,x0rehilolo,x0relololo,
                 &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                 &acchihilo,&acclohilo,&acchilolo,&acclololo);
         odf_sub(x0rehihihi,x0relohihi,x0rehilohi,x0relolohi,
                 x0rehihilo,x0relohilo,x0rehilolo,x0relololo,
                 acchihihi,acclohihi,acchilohi,acclolohi,
                 acchihilo,acclohilo,acchilolo,acclololo,
                 &v0rehihihi,&v0relohihi,&v0rehilohi,&v0relolohi,
                 &v0rehihilo,&v0relohilo,&v0rehilolo,&v0relololo);
         // vim[0] = xim[0] - mu*xim[0];
         odf_mul(muhihihi,mulohihi,muhilohi,mulolohi,
                 muhihilo,mulohilo,muhilolo,mulololo,
                 x0imhihihi,x0imlohihi,x0imhilohi,x0imlolohi,
                 x0imhihilo,x0imlohilo,x0imhilolo,x0imlololo,
                 &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                 &acchihilo,&acclohilo,&acchilolo,&acclololo);
         odf_sub(x0imhihihi,x0imlohihi,x0imhilohi,x0imlolohi,
                 x0imhihilo,x0imlohilo,x0imhilolo,x0imlololo,
                 acchihihi,acclohihi,acchilohi,acclolohi,
                 acchihilo,acclohilo,acchilolo,acclololo,
                 &v0imhihihi,&v0imlohihi,&v0imhilohi,&v0imlolohi,
                 &v0imhihilo,&v0imlohilo,&v0imhilolo,&v0imlololo);
      }
      // sqrv0 = vre[0]*vre[0] + vim[0]*vim[0];
      odf_sqr(v0rehihihi,v0relohihi,v0rehilohi,v0relolohi,
              v0rehihilo,v0relohilo,v0rehilolo,v0relololo,
              &sqrv0hihihi,&sqrv0lohihi,&sqrv0hilohi,&sqrv0lolohi,
              &sqrv0hihilo,&sqrv0lohilo,&sqrv0hilolo,&sqrv0lololo);
      odf_sqr(v0imhihihi,v0imlohihi,v0imhilohi,v0imlolohi,
              v0imhihilo,v0imlohilo,v0imhilolo,v0imlololo,
              &acchihihi,&acclohihi,&acchilohi,&acclolohi,
              &acchihilo,&acclohilo,&acchilolo,&acclololo);
      odf_inc(&sqrv0hihihi,&sqrv0lohihi,&sqrv0hilohi,&sqrv0lolohi,
              &sqrv0hihilo,&sqrv0lohilo,&sqrv0hilolo,&sqrv0lololo,
              acchihihi,acclohihi,acchilohi,acclolohi,
              acchihilo,acclohilo,acchilolo,acclololo);
      // *beta = 2.0*sqrv0/(sigma + sqrv0);
      odf_inc(sigmahihihi_h,sigmalohihi_h,sigmahilohi_h,sigmalolohi_h,
              sigmahihilo_h,sigmalohilo_h,sigmahilolo_h,sigmalololo_h,
              sqrv0hihihi,sqrv0lohihi,sqrv0hilohi,sqrv0lolohi,
              sqrv0hihilo,sqrv0lohilo,sqrv0hilolo,sqrv0lololo);
      odf_div(sqrv0hihihi,sqrv0lohihi,sqrv0hilohi,sqrv0lolohi,
              sqrv0hihilo,sqrv0lohilo,sqrv0hilolo,sqrv0lololo,
              *sigmahihihi_h,*sigmalohihi_h,*sigmahilohi_h,*sigmalolohi_h,
              *sigmahihilo_h,*sigmalohilo_h,*sigmahilolo_h,*sigmalololo_h,
              &betahihihi_h[L],&betalohihi_h[L],
              &betahilohi_h[L],&betalolohi_h[L],
              &betahihilo_h[L],&betalohilo_h[L],
              &betahilolo_h[L],&betalololo_h[L]);
      odf_mlt_d(&betahihihi_h[L],&betalohihi_h[L],
                &betahilohi_h[L],&betalolohi_h[L],
                &betahihilo_h[L],&betalohilo_h[L],
                &betahilolo_h[L],&betalololo_h[L],2.0);
      // inv0re = vre[0]/sqrv0;  // real part of 1/v[0]
      odf_div(v0rehihihi,v0relohihi,v0rehilohi,v0relolohi,
              v0rehihilo,v0relohilo,v0rehilolo,v0relololo,
              sqrv0hihihi,sqrv0lohihi,sqrv0hilohi,sqrv0lolohi,
              sqrv0hihilo,sqrv0lohilo,sqrv0hilolo,sqrv0lololo,
              &inv0rehihihi,&inv0relohihi,&inv0rehilohi,&inv0relolohi,
              &inv0rehihilo,&inv0relohilo,&inv0rehilolo,&inv0relololo);
      // inv0im = -vim[0]/sqrv0; // imaginary part of 1/v[0]
      odf_div(v0imhihihi,v0imlohihi,v0imhilohi,v0imlolohi,
              v0imhihilo,v0imlohilo,v0imhilolo,v0imlololo,
              sqrv0hihihi,sqrv0lohihi,sqrv0hilohi,sqrv0lolohi,
              sqrv0hihilo,sqrv0lohilo,sqrv0hilolo,sqrv0lololo,
              &inv0imhihihi,&inv0imlohihi,&inv0imhilohi,&inv0imlolohi,
              &inv0imhihilo,&inv0imlohilo,&inv0imhilolo,&inv0imlololo);
      odf_minus(&inv0imhihihi,&inv0imlohihi,&inv0imhilohi,&inv0imlolohi,
                &inv0imhihilo,&inv0imlohilo,&inv0imhilolo,&inv0imlololo);
      *sigmahihihi_h = inv0rehihihi;
      *sigmalohihi_h = inv0relohihi;
      *sigmahilohi_h = inv0rehilohi;
      *sigmalolohi_h = inv0relolohi;
      *sigmahihilo_h = inv0rehihilo;
      *sigmalohilo_h = inv0relohilo;
      *sigmahilolo_h = inv0rehilolo;
      *sigmalololo_h = inv0relololo;
      betahihihi_h[szt] = inv0imhihihi;
      betalohihi_h[szt] = inv0imlohihi;
      betahilohi_h[szt] = inv0imhilohi;
      betalolohi_h[szt] = inv0imlolohi;
      betahihilo_h[szt] = inv0imhihilo;
      betalohilo_h[szt] = inv0imlohilo;
      betahilolo_h[szt] = inv0imhilolo;
      betalololo_h[szt] = inv0imlololo;
      // update the flop counts
      *add += 6;
      *mul += 7;
      *div += 4;
      *sqrtfun += 2;
   }
   if(verbose)
   {
      cout << scientific << setprecision(16)
           << "beta[" << colidx << "] : "
           << betahihihi_h[L] << "  " << betalohihi_h[L] << endl
           << "          "
           << betahilohi_h[L] << "  " << betalolohi_h[L] << endl
           << "          "
           << betahihilo_h[L] << "  " << betalohilo_h[L] << endl
           << "          "
           << betahilolo_h[L] << "  " << betalololo_h[L] << endl;
   }
   cudaMemcpy(&betahihihi_d[L],&betahihihi_h[L],sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&betalohihi_d[L],&betalohihi_h[L],sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&betahilohi_d[L],&betahilohi_h[L],sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&betalolohi_d[L],&betalolohi_h[L],sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&betahihilo_d[L],&betahihilo_h[L],sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&betalohilo_d[L],&betalohilo_h[L],sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&betahilolo_d[L],&betahilolo_h[L],sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&betalololo_d[L],&betalololo_h[L],sizeof(double),
              cudaMemcpyHostToDevice);

   if(!done)  // normalization needed
   {
      // (sigmahi_h, sigmalo_h) has the values for (v0rehi, v0relo)
      cudaMemcpy(sigmahihihi_d,sigmahihihi_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(sigmalohihi_d,sigmalohihi_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(sigmahilohi_d,sigmahilohi_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(sigmalolohi_d,sigmalolohi_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(sigmahihilo_d,sigmahihilo_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(sigmalohilo_d,sigmalohilo_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(sigmahilolo_d,sigmahilolo_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(sigmalololo_d,sigmalololo_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      // (betahi_h[szt], betalo_h[szt]) has the values for (v0imhi, v0rimlo).
      cudaMemcpy(&betahihihi_d[szt],&betahihihi_h[szt],sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&betalohihi_d[szt],&betalohihi_h[szt],sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&betahilohi_d[szt],&betahilohi_h[szt],sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&betalolohi_d[szt],&betalolohi_h[szt],sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&betahihilo_d[szt],&betahihilo_h[szt],sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&betalohilo_d[szt],&betalohilo_h[szt],sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&betahilolo_d[szt],&betahilolo_h[szt],sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&betalololo_d[szt],&betalololo_h[szt],sizeof(double),
                 cudaMemcpyHostToDevice);

      if(verbose)
      {
         cout << "-> launching " << nblocks << " blocks of "
              << szt << " threads to normalize ..." << endl;
         cout << "   nrows1 : " << nrows1
              << "  rowidx : " << rowidx << "  nVrows : " << nVrows << endl;
      }
/*
      cudaEventRecord(start);
      cmplx8_normalize<<<nblocks,szt>>>
         (nrows1,szt,&Arehihihi_d[rowidx+1],&Arelohihi_d[rowidx+1],
                     &Arehilohi_d[rowidx+1],&Arelolohi_d[rowidx+1],
                     &Arehihilo_d[rowidx+1],&Arelohilo_d[rowidx+1],
                     &Arehilolo_d[rowidx+1],&Arelololo_d[rowidx+1],
                     &Aimhihihi_d[rowidx+1],&Aimlohihi_d[rowidx+1],
                     &Aimhilohi_d[rowidx+1],&Aimlolohi_d[rowidx+1],
                     &Aimhihilo_d[rowidx+1],&Aimlohilo_d[rowidx+1],
                     &Aimhilolo_d[rowidx+1],&Aimlololo_d[rowidx+1],
          sigmahihihi_d,sigmalohihi_d,sigmahilohi_d,sigmalolohi_d,
          sigmahihilo_d,sigmalohilo_d,sigmahilolo_d,sigmalololo_d,
          &betahihihi_d[szt],&betalohihi_d[szt],
          &betahilohi_d[szt],&betalolohi_d[szt],
          &betahihilo_d[szt],&betalohilo_d[szt],
          &betahilolo_d[szt],&betalololo_d[szt],
          &Vrehihihi_d[L*nVrows+L+1],&Vrelohihi_d[L*nVrows+L+1],
          &Vrehilohi_d[L*nVrows+L+1],&Vrelolohi_d[L*nVrows+L+1],
          &Vrehihilo_d[L*nVrows+L+1],&Vrelohilo_d[L*nVrows+L+1],
          &Vrehilolo_d[L*nVrows+L+1],&Vrelololo_d[L*nVrows+L+1],
          &Vimhihihi_d[L*nVrows+L+1],&Vimlohihi_d[L*nVrows+L+1],
          &Vimhilohi_d[L*nVrows+L+1],&Vimlolohi_d[L*nVrows+L+1],
          &Vimhihilo_d[L*nVrows+L+1],&Vimlohilo_d[L*nVrows+L+1],
          &Vimhilolo_d[L*nVrows+L+1],&Vimlololo_d[L*nVrows+L+1]);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds,start,stop);
      *lapms += milliseconds;
*/
      cudaEventRecord(start);
      cmplx8_normalize_rere<<<nblocks,szt>>>
         (nrows1,szt,&Arehihihi_d[rowidx+1],&Arelohihi_d[rowidx+1],
                     &Arehilohi_d[rowidx+1],&Arelolohi_d[rowidx+1],
                     &Arehihilo_d[rowidx+1],&Arelohilo_d[rowidx+1],
                     &Arehilolo_d[rowidx+1],&Arelololo_d[rowidx+1],
          sigmahihihi_d,sigmalohihi_d,sigmahilohi_d,sigmalolohi_d,
          sigmahihilo_d,sigmalohilo_d,sigmahilolo_d,sigmalololo_d,
          &Vrehihihi_d[L*nVrows+L+1],&Vrelohihi_d[L*nVrows+L+1],
          &Vrehilohi_d[L*nVrows+L+1],&Vrelolohi_d[L*nVrows+L+1],
          &Vrehihilo_d[L*nVrows+L+1],&Vrelohilo_d[L*nVrows+L+1],
          &Vrehilolo_d[L*nVrows+L+1],&Vrelololo_d[L*nVrows+L+1]);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds,start,stop);
      *lapms += milliseconds;

      cudaEventRecord(start);
      cmplx8_normalize_imim<<<nblocks,szt>>>
         (nrows1,szt, &Aimhihihi_d[rowidx+1],&Aimlohihi_d[rowidx+1],
                     &Aimhilohi_d[rowidx+1],&Aimlolohi_d[rowidx+1],
                     &Aimhihilo_d[rowidx+1],&Aimlohilo_d[rowidx+1],
                     &Aimhilolo_d[rowidx+1],&Aimlololo_d[rowidx+1],
          &betahihihi_d[szt],&betalohihi_d[szt],
          &betahilohi_d[szt],&betalolohi_d[szt],
          &betahihilo_d[szt],&betalohilo_d[szt],
          &betahilolo_d[szt],&betalololo_d[szt],
          &Vrehihihi_d[L*nVrows+L+1],&Vrelohihi_d[L*nVrows+L+1],
          &Vrehilohi_d[L*nVrows+L+1],&Vrelolohi_d[L*nVrows+L+1],
          &Vrehihilo_d[L*nVrows+L+1],&Vrelohilo_d[L*nVrows+L+1],
          &Vrehilolo_d[L*nVrows+L+1],&Vrelololo_d[L*nVrows+L+1]);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds,start,stop);
      *lapms += milliseconds;

      cudaEventRecord(start);
      cmplx8_normalize_imre<<<nblocks,szt>>>
         (nrows1,szt,&Aimhihihi_d[rowidx+1],&Aimlohihi_d[rowidx+1],
                     &Aimhilohi_d[rowidx+1],&Aimlolohi_d[rowidx+1],
                     &Aimhihilo_d[rowidx+1],&Aimlohilo_d[rowidx+1],
                     &Aimhilolo_d[rowidx+1],&Aimlololo_d[rowidx+1],
          sigmahihihi_d,sigmalohihi_d,sigmahilohi_d,sigmalolohi_d,
          sigmahihilo_d,sigmalohilo_d,sigmahilolo_d,sigmalololo_d,
          &Vimhihihi_d[L*nVrows+L+1],&Vimlohihi_d[L*nVrows+L+1],
          &Vimhilohi_d[L*nVrows+L+1],&Vimlolohi_d[L*nVrows+L+1],
          &Vimhihilo_d[L*nVrows+L+1],&Vimlohilo_d[L*nVrows+L+1],
          &Vimhilolo_d[L*nVrows+L+1],&Vimlololo_d[L*nVrows+L+1]);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds,start,stop);
      *lapms += milliseconds;

      cudaEventRecord(start);
      cmplx8_normalize_reim<<<nblocks,szt>>>
         (nrows1,szt,&Arehihihi_d[rowidx+1],&Arelohihi_d[rowidx+1],
                     &Arehilohi_d[rowidx+1],&Arelolohi_d[rowidx+1],
                     &Arehihilo_d[rowidx+1],&Arelohilo_d[rowidx+1],
                     &Arehilolo_d[rowidx+1],&Arelololo_d[rowidx+1],
          &betahihihi_d[szt],&betalohihi_d[szt],
          &betahilohi_d[szt],&betalolohi_d[szt],
          &betahihilo_d[szt],&betalohilo_d[szt],
          &betahilolo_d[szt],&betalololo_d[szt],
          &Vimhihihi_d[L*nVrows+L+1],&Vimlohihi_d[L*nVrows+L+1],
          &Vimhilohi_d[L*nVrows+L+1],&Vimlolohi_d[L*nVrows+L+1],
          &Vimhihilo_d[L*nVrows+L+1],&Vimlohilo_d[L*nVrows+L+1],
          &Vimhilolo_d[L*nVrows+L+1],&Vimlololo_d[L*nVrows+L+1]);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds,start,stop);
      *lapms += milliseconds;

      flopcount_cmplx_normalize(nblocks,szt,add,mul);
   }
   if(verbose)
   {
      const size_t szhouse = nVrows*sizeof(double);

      cudaMemcpy(&betahihihi_h[L],&betahihihi_d[L],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&betalohihi_h[L],&betalohihi_d[L],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&betahilohi_h[L],&betahilohi_d[L],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&betalolohi_h[L],&betalolohi_d[L],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&betahihilo_h[L],&betahihilo_d[L],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&betalohilo_h[L],&betalohilo_d[L],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&betahilolo_h[L],&betahilolo_d[L],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&betalololo_h[L],&betalololo_d[L],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vrehihihi_h,&Vrehihihi_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vrelohihi_h,&Vrelohihi_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vrehilohi_h,&Vrehilohi_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vrelolohi_h,&Vrelolohi_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vrehihilo_h,&Vrehihilo_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vrelohilo_h,&Vrelohilo_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vrehilolo_h,&Vrehilolo_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vrelololo_h,&Vrelololo_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vimhihihi_h,&Vimhihihi_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vimlohihi_h,&Vimlohihi_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vimhilohi_h,&Vimhilohi_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vimlolohi_h,&Vimlolohi_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vimhihilo_h,&Vimhihilo_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vimlohilo_h,&Vimlohilo_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vimhilolo_h,&Vimhilolo_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vimlololo_h,&Vimlololo_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);

      cout << scientific << setprecision(16)
           << "beta[" << colidx << "] : "
           << betahihihi_h[L] << "  " << betalohihi_h[L] << endl
           << "          "
           << betahilohi_h[L] << "  " << betalolohi_h[L] << endl
           << "          "
           << betahihilo_h[L] << "  " << betalohilo_h[L] << endl
           << "          "
           << betahilolo_h[L] << "  " << betalololo_h[L] << endl;

      for(int i=0; i<nVrows; i++)
      {
         cout << "v[" << i << "]re : "
              << vrehihihi_h[i] << "  " << vrelohihi_h[i] << endl
              << "         "
              << vrehilohi_h[i] << "  " << vrelolohi_h[i] << endl
              << "         "
              << vrehihilo_h[i] << "  " << vrelohilo_h[i] << endl
              << "         "
              << vrehilolo_h[i] << "  " << vrelololo_h[i] << endl;
         cout << "v[" << i << "]im : "
              << vimhihihi_h[i] << "  " << vimlohihi_h[i] << endl
              << "         "
              << vimhilohi_h[i] << "  " << vimlolohi_h[i] << endl
              << "         "
              << vimhihilo_h[i] << "  " << vimlohilo_h[i] << endl
              << "         "
              << vimhilolo_h[i] << "  " << vimlololo_h[i] << endl;
      }
   }
}

void GPU_dbl8_small_leftRupdate
 ( int nrows, int ncols, int szt, int colidx, int k, int L,
   double *Ahihihi_h, double *Alohihi_h, double *Ahilohi_h, double *Alolohi_h,
   double *Ahihilo_h, double *Alohilo_h, double *Ahilolo_h, double *Alololo_h,
   double *Ahihihi_d, double *Alohihi_d, double *Ahilohi_d, double *Alolohi_d,
   double *Ahihilo_d, double *Alohilo_d, double *Ahilolo_d, double *Alololo_d,
   double *Vhihihi_d, double *Vlohihi_d, double *Vhilohi_d, double *Vlolohi_d,
   double *Vhihilo_d, double *Vlohilo_d, double *Vhilolo_d, double *Vlololo_d,
   double *betahihihi_h, double *betalohihi_h,
   double *betahilohi_h, double *betalolohi_h,
   double *betahihilo_h, double *betalohilo_h,
   double *betahilolo_h, double *betalololo_h,
   double *betahihihi_d, double *betalohihi_d,
   double *betahilohi_d, double *betalolohi_d,
   double *betahihilo_d, double *betalohilo_d,
   double *betahilolo_d, double *betalololo_d,
   double *lapms, long long int *add, long long int *mul, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int endcol = (k+1)*szt;     // 1 + last column index in tile
   const int nVrows = nrows - k*szt;          // dimension of V matrix

   cudaEventRecord(start);           // 2nd argument: ncols -> szt
   // changed second argument ncols into szt
   // to avoid updating the next tile
   dbl8_small_leftRupdate<<<1,nrows-colidx>>>
      (nrows,endcol,szt,colidx,
       Ahihihi_d,Alohihi_d,Ahilohi_d,Alolohi_d,
       Ahihilo_d,Alohilo_d,Ahilolo_d,Alololo_d,
       &Vhihihi_d[L*nVrows+L],&Vlohihi_d[L*nVrows+L],
       &Vhilohi_d[L*nVrows+L],&Vlolohi_d[L*nVrows+L],
       &Vhihilo_d[L*nVrows+L],&Vlohilo_d[L*nVrows+L],
       &Vhilolo_d[L*nVrows+L],&Vlololo_d[L*nVrows+L],
       &betahihihi_d[L],&betalohihi_d[L],&betahilohi_d[L],&betalolohi_d[L],
       &betahihilo_d[L],&betalohilo_d[L],&betahilolo_d[L],&betalololo_d[L]);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_dbl_small_leftRupdate(nrows,ncols,szt,colidx,add,mul);

   if(verbose)
   {
      const int dim = nrows*ncols;
      const size_t sznum = dim*sizeof(double);

      cudaMemcpy(Ahihihi_h,Ahihihi_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Alohihi_h,Alohihi_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Ahilohi_h,Ahilohi_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Alolohi_h,Alolohi_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Ahihilo_h,Ahihilo_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Alohilo_h,Alohilo_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Ahilolo_h,Ahilolo_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Alololo_h,Alololo_d,sznum,cudaMemcpyDeviceToHost);

      cout << "the matrix after the update :" << endl;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
            cout << "A_d[" << i << "][" << j << "] : "
                 << Ahihihi_h[j*nrows+i] << "  "
                 << Alohihi_h[j*nrows+i] << endl
                 << "            "
                 << Ahilohi_h[j*nrows+i] << "  "
                 << Alolohi_h[j*nrows+i] << endl
                 << "            "
                 << Ahihilo_h[j*nrows+i] << "  "
                 << Alohilo_h[j*nrows+i] << endl
                 << "            "
                 << Ahilolo_h[j*nrows+i] << "  "
                 << Alololo_h[j*nrows+i] << endl;
   }
}

void GPU_cmplx8_small_leftRupdate
 ( int nrows, int ncols, int szt, int colidx, int k, int L,
   double *Arehihihi_h, double *Arelohihi_h,
   double *Arehilohi_h, double *Arelolohi_h,
   double *Arehihilo_h, double *Arelohilo_h,
   double *Arehilolo_h, double *Arelololo_h,
   double *Aimhihihi_h, double *Aimlohihi_h,
   double *Aimhilohi_h, double *Aimlolohi_h,
   double *Aimhihilo_h, double *Aimlohilo_h,
   double *Aimhilolo_h, double *Aimlololo_h,
   double *Arehihihi_d, double *Arelohihi_d,
   double *Arehilohi_d, double *Arelolohi_d,
   double *Arehihilo_d, double *Arelohilo_d,
   double *Arehilolo_d, double *Arelololo_d,
   double *Aimhihihi_d, double *Aimlohihi_d,
   double *Aimhilohi_d, double *Aimlolohi_d,
   double *Aimhihilo_d, double *Aimlohilo_d,
   double *Aimhilolo_d, double *Aimlololo_d,
   double *Vrehihihi_d, double *Vrelohihi_d,
   double *Vrehilohi_d, double *Vrelolohi_d,
   double *Vrehihilo_d, double *Vrelohilo_d,
   double *Vrehilolo_d, double *Vrelololo_d,
   double *Vimhihihi_d, double *Vimlohihi_d,
   double *Vimhilohi_d, double *Vimlolohi_d,
   double *Vimhihilo_d, double *Vimlohilo_d,
   double *Vimhilolo_d, double *Vimlololo_d,
   double *betahihihi_h, double *betalohihi_h,
   double *betahilohi_h, double *betalolohi_h,
   double *betahihilo_h, double *betalohilo_h,
   double *betahilolo_h, double *betalololo_h,
   double *betahihihi_d, double *betalohihi_d,
   double *betahilohi_d, double *betalolohi_d,
   double *betahihilo_d, double *betalohilo_d,
   double *betahilolo_d, double *betalololo_d,
   double *lapms, long long int *add, long long int *mul, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int endcol = (k+1)*szt;     // 1 + last column index in tile
   const int nVrows = nrows - k*szt;          // dimension of V matrix

   cudaEventRecord(start);
   cmplx8_small_leftRupdate<<<1,nrows-colidx>>>
      (nrows,endcol,szt,colidx,
       Arehihihi_d,Arelohihi_d,Arehilohi_d,Arelolohi_d,
       Arehihilo_d,Arelohilo_d,Arehilolo_d,Arelololo_d,
       Aimhihihi_d,Aimlohihi_d,Aimhilohi_d,Aimlolohi_d,
       Aimhihilo_d,Aimlohilo_d,Aimhilolo_d,Aimlololo_d,
       &Vrehihihi_d[L*nVrows+L],&Vrelohihi_d[L*nVrows+L],
       &Vrehilohi_d[L*nVrows+L],&Vrelolohi_d[L*nVrows+L],
       &Vrehihilo_d[L*nVrows+L],&Vrelohilo_d[L*nVrows+L],
       &Vrehilolo_d[L*nVrows+L],&Vrelololo_d[L*nVrows+L],
       &Vimhihihi_d[L*nVrows+L],&Vimlohihi_d[L*nVrows+L],
       &Vimhilohi_d[L*nVrows+L],&Vimlolohi_d[L*nVrows+L],
       &Vimhihilo_d[L*nVrows+L],&Vimlohilo_d[L*nVrows+L],
       &Vimhilolo_d[L*nVrows+L],&Vimlololo_d[L*nVrows+L],
       &betahihihi_d[L],&betalohihi_d[L],&betahilohi_d[L],&betalolohi_d[L],
       &betahihilo_d[L],&betalohilo_d[L],&betahilolo_d[L],&betalololo_d[L]);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_cmplx_small_leftRupdate(nrows,ncols,szt,colidx,add,mul);

   if(verbose)
   {
      const int dim = nrows*ncols;
      const size_t sznum = dim*sizeof(double);

      cudaMemcpy(Arehihihi_h,Arehihihi_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Arelohihi_h,Arelohihi_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Arehilohi_h,Arehilohi_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Arelolohi_h,Arelolohi_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Arehihilo_h,Arehihilo_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Arelohilo_h,Arelohilo_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Arehilolo_h,Arehilolo_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Arelololo_h,Arelololo_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Aimhihihi_h,Aimhihihi_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Aimlohihi_h,Aimlohihi_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Aimhilohi_h,Aimhilohi_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Aimlolohi_h,Aimlolohi_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Aimhihilo_h,Aimhihilo_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Aimlohilo_h,Aimlohilo_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Aimhilolo_h,Aimhilolo_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Aimlololo_h,Aimlololo_d,sznum,cudaMemcpyDeviceToHost);

      cout << "the matrix after the update :" << endl;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
         {
            cout << "A_d[" << i << "][" << j << "]re : "
                 << Arehihihi_h[j*nrows+i] << "  "
                 << Arelohihi_h[j*nrows+i] << endl
                 << "              "
                 << Arehilohi_h[j*nrows+i] << "  "
                 << Arelolohi_h[j*nrows+i] << endl
                 << "              "
                 << Arehihilo_h[j*nrows+i] << "  "
                 << Arelohilo_h[j*nrows+i] << endl
                 << "              "
                 << Arehilolo_h[j*nrows+i] << "  "
                 << Arelololo_h[j*nrows+i] << endl;
            cout << "A_d[" << i << "][" << j << "]im : "
                 << Aimhihihi_h[j*nrows+i] << "  "
                 << Aimlohihi_h[j*nrows+i] << endl
                 << "              "
                 << Aimhilohi_h[j*nrows+i] << "  "
                 << Aimlolohi_h[j*nrows+i] << endl
                 << "              "
                 << Aimhihilo_h[j*nrows+i] << "  "
                 << Aimlohilo_h[j*nrows+i] << endl
                 << "              "
                 << Aimhilolo_h[j*nrows+i] << "  "
                 << Aimlololo_h[j*nrows+i] << endl;
         }
   }
}

void GPU_dbl8_medium_leftRupdate
 ( int nrows, int ncols, int szt, int colidx, int k, int L,
   double *Ahihihi_h, double *Alohihi_h, double *Ahilohi_h, double *Alolohi_h,
   double *Ahihilo_h, double *Alohilo_h, double *Ahilolo_h, double *Alololo_h,
   double *Ahihihi_d, double *Alohihi_d, double *Ahilohi_d, double *Alolohi_d,
   double *Ahihilo_d, double *Alohilo_d, double *Ahilolo_d, double *Alololo_d,
   double *Vhihihi_d, double *Vlohihi_d, double *Vhilohi_d, double *Vlolohi_d,
   double *Vhihilo_d, double *Vlohilo_d, double *Vhilolo_d, double *Vlololo_d,
   double *betahihihi_h, double *betalohihi_h,
   double *betahilohi_h, double *betalolohi_h,
   double *betahihilo_h, double *betalohilo_h,
   double *betahilolo_h, double *betalololo_h,
   double *betahihihi_d, double *betalohihi_d,
   double *betahilohi_d, double *betalolohi_d,
   double *betahihilo_d, double *betalohilo_d,
   double *betahilolo_d, double *betalololo_d,
   double *RTdotvhihihi_h, double *RTdotvlohihi_h,
   double *RTdotvhilohi_h, double *RTdotvlolohi_h,
   double *RTdotvhihilo_h, double *RTdotvlohilo_h,
   double *RTdotvhilolo_h, double *RTdotvlololo_h,
   double *RTdotvhihihi_d, double *RTdotvlohihi_d,
   double *RTdotvhilohi_d, double *RTdotvlolohi_d,
   double *RTdotvhihilo_d, double *RTdotvlohilo_d,
   double *RTdotvhilolo_d, double *RTdotvlololo_d,
   double *whihihi_h, double *wlohihi_h, double *whilohi_h, double *wlolohi_h,
   double *whihilo_h, double *wlohilo_h, double *whilolo_h, double *wlololo_h,
   double *whihihi_d, double *wlohihi_d, double *whilohi_d, double *wlolohi_d,
   double *whihilo_d, double *wlohilo_d, double *whilolo_d, double *wlololo_d,
   double *RTvlapms, double *redlapms,
   long long int *add, long long int *mul, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int endcol = (k+1)*szt;     // 1 + last column index in tile
   const int nVrows = nrows - k*szt;          // dimension of V matrix
   const int nhouse = nrows - colidx;  // length of Householder vector
   // total number of entries in R that will be modified
   const int RToffset = colidx*nrows;
   const int dimRTdotv = endcol - colidx;
   // total number of entries in R that will be modified
   const int sizenum = (nrows - colidx)*dimRTdotv;
   const int nbrblocks = (int) ceil(sizenum/((double) szt));

   // changed second argument ncols into endcol
   // to avoid updating the next tile
   // dbl_medium_betaRTv<<<nbrblocks,szt>>>
   //   (nrows,endcol,szt,colidx,A_d,&V_d[L*nVrows+L],&beta_d[L],w_d);
   // number of threads must be ncols - colidx, not endcol - colidx
   // dbl2_small_betaRTv<<<1,nrows-colidx>>> // nrows ...
   //   (nrows,endcol,szt,colidx,Ahi_d,Alo_d,
   //    &Vhi_d[L*nVrows+L],&Vlo_d[L*nVrows+L],
   //    &betahi_d[L],&betalo_d[L],whi_d,wlo_d);

   if(verbose)
   {
      cout << "-> launching " << nbrblocks << " blocks of " << szt
           << " threads to compute RTdotv ..." << endl;
      cout << "   nhouse : " << nhouse << "  RToffset : " << RToffset
           << "  dimRTdotv : " << dimRTdotv << endl;
   }

   cudaEventRecord(start);
   dbl8_RTdotv<<<nbrblocks,szt>>>
      (nhouse,szt,colidx,RToffset,dimRTdotv,
       Ahihihi_d,Alohihi_d,Ahilohi_d,Alolohi_d,
       Ahihilo_d,Alohilo_d,Ahilolo_d,Alololo_d,
       &Vhihihi_d[L*nVrows+L],&Vlohihi_d[L*nVrows+L],
       &Vhilohi_d[L*nVrows+L],&Vlolohi_d[L*nVrows+L],
       &Vhihilo_d[L*nVrows+L],&Vlohilo_d[L*nVrows+L],
       &Vhilolo_d[L*nVrows+L],&Vlololo_d[L*nVrows+L],
       RTdotvhihihi_d,RTdotvlohihi_d,RTdotvhilohi_d,RTdotvlolohi_d,
       RTdotvhihilo_d,RTdotvlohilo_d,RTdotvhilolo_d,RTdotvlololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *RTvlapms += milliseconds;
   cudaEventRecord(start);
   dbl8_sum_betaRTdotv<<<1,dimRTdotv>>>
      (nhouse,&betahihihi_d[L],&betalohihi_d[L],
              &betahilohi_d[L],&betalolohi_d[L],
              &betahihilo_d[L],&betalohilo_d[L],
              &betahilolo_d[L],&betalololo_d[L],
       RTdotvhihihi_d,RTdotvlohihi_d,RTdotvhilohi_d,RTdotvlolohi_d,
       RTdotvhihilo_d,RTdotvlohilo_d,RTdotvhilolo_d,RTdotvlololo_d,
       whihihi_d,wlohihi_d,whilohi_d,wlolohi_d,
       whihilo_d,wlohilo_d,whilolo_d,wlololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *RTvlapms += milliseconds;
   flopcount_dbl_RTdotv(nhouse,szt,mul);
   flopcount_dbl_sum_betaRTdotv(nhouse,dimRTdotv,add,mul);

   if(verbose)
   {
      cout << "-> launching " << nbrblocks << " blocks of " << szt
           << " threads to update " << sizenum << " numbers ..." << endl;
      cout << "   nrows : " << nrows << "  endcol : " << endcol
           << "  szt : " << szt << "  colidx : " << colidx << endl;
   }
   cudaEventRecord(start);
   dbl8_medium_subvbetaRTv<<<nbrblocks,szt>>>
      (nrows,endcol,szt,colidx,
       Ahihihi_d,Alohihi_d,Ahilohi_d,Alolohi_d,
       Ahihilo_d,Alohilo_d,Ahilolo_d,Alololo_d,
       &Vhihihi_d[L*nVrows+L],&Vlohihi_d[L*nVrows+L],
       &Vhilohi_d[L*nVrows+L],&Vlolohi_d[L*nVrows+L],
       &Vhihilo_d[L*nVrows+L],&Vlohilo_d[L*nVrows+L],
       &Vhilolo_d[L*nVrows+L],&Vlololo_d[L*nVrows+L],
       &betahihihi_d[L],&betalohihi_d[L],&betahilohi_d[L],&betalolohi_d[L],
       &betahihilo_d[L],&betalohilo_d[L],&betahilolo_d[L],&betalololo_d[L],
       whihihi_d,wlohihi_d,whilohi_d,wlolohi_d,
       whihilo_d,wlohilo_d,whilolo_d,wlololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *redlapms += milliseconds;
   flopcount_dbl_medium_subvbetaRTv(nrows,endcol,szt,colidx,add,mul);

   if(verbose)
   {
      const int dim = nrows*ncols;
      const size_t sznum = dim*sizeof(double);
      const size_t szbRTv = dimRTdotv*sizeof(double);
      const size_t szRTdotv = nVrows*szbRTv;

      cudaMemcpy(RTdotvhihihi_h,RTdotvhihihi_d,szRTdotv,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(RTdotvlohihi_h,RTdotvlohihi_d,szRTdotv,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(RTdotvhilohi_h,RTdotvhilohi_d,szRTdotv,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(RTdotvlolohi_h,RTdotvlolohi_d,szRTdotv,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(RTdotvhihilo_h,RTdotvhihilo_d,szRTdotv,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(RTdotvlohilo_h,RTdotvlohilo_d,szRTdotv,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(RTdotvhilolo_h,RTdotvhilolo_d,szRTdotv,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(RTdotvlololo_h,RTdotvlololo_d,szRTdotv,
                 cudaMemcpyDeviceToHost);

      cout << "the matrix R^T dot v : " << endl;
      int ix = 0;
      for(int i=0; i<endcol-colidx; i++)
      {
         for(int j=0; j<nhouse; j++)      // must use nhouse
         {
            cout << "RTdotv[" << i << "][" << j << "] : "
                 << RTdotvhihihi_h[ix] << "  "
                 << RTdotvlohihi_h[ix] << endl
                 << "               "
                 << RTdotvhilohi_h[ix] << "  "
                 << RTdotvlolohi_h[ix] << endl
                 << "               "
                 << RTdotvhihilo_h[ix] << "  "
                 << RTdotvlohilo_h[ix] << endl
                 << "               "
                 << RTdotvhilolo_h[ix] << "  "
                 << RTdotvlololo_h[ix] << endl;
            ix = ix + 1;
         }
      }
      cudaMemcpy(whihihi_h,whihihi_d,szbRTv,cudaMemcpyDeviceToHost);
      cudaMemcpy(wlohihi_h,wlohihi_d,szbRTv,cudaMemcpyDeviceToHost);
      cudaMemcpy(whilohi_h,whilohi_d,szbRTv,cudaMemcpyDeviceToHost);
      cudaMemcpy(wlolohi_h,wlolohi_d,szbRTv,cudaMemcpyDeviceToHost);
      cudaMemcpy(whihilo_h,whihilo_d,szbRTv,cudaMemcpyDeviceToHost);
      cudaMemcpy(wlohilo_h,wlohilo_d,szbRTv,cudaMemcpyDeviceToHost);
      cudaMemcpy(whilolo_h,whilolo_d,szbRTv,cudaMemcpyDeviceToHost);
      cudaMemcpy(wlololo_h,wlololo_d,szbRTv,cudaMemcpyDeviceToHost);

      cout << "the vector w = beta*R^T*v : " << endl;
      for(int i=0; i<endcol-colidx; i++)
         cout << "w[" << i << "] : "
              << whihihi_h[i] << "  " << wlohihi_h[i] << endl
              << "       "
              << whilohi_h[i] << "  " << wlolohi_h[i] << endl
              << "       "
              << whihilo_h[i] << "  " << wlohilo_h[i] << endl
              << "       "
              << whilolo_h[i] << "  " << wlololo_h[i] << endl;

      cudaMemcpy(Ahihihi_h,Ahihihi_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Alohihi_h,Alohihi_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Ahilohi_h,Ahilohi_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Alolohi_h,Alolohi_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Ahihilo_h,Ahihilo_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Alohilo_h,Alohilo_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Ahilolo_h,Ahilolo_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Alololo_h,Alololo_d,sznum,cudaMemcpyDeviceToHost);

      cout << "the matrix after the update :" << endl;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
            cout << "A_d[" << i << "][" << j << "] : "
                 << Ahihihi_h[j*nrows+i] << "  "
                 << Alohihi_h[j*nrows+i] << endl
                 << "            "
                 << Ahilohi_h[j*nrows+i] << "  "
                 << Alolohi_h[j*nrows+i] << endl
                 << "            "
                 << Ahihilo_h[j*nrows+i] << "  "
                 << Alohilo_h[j*nrows+i] << endl
                 << "            "
                 << Ahilolo_h[j*nrows+i] << "  "
                 << Alololo_h[j*nrows+i] << endl;
   }
}

void GPU_cmplx8_medium_leftRupdate
 ( int nrows, int ncols, int szt, int colidx, int k, int L,
   double *Arehihihi_h, double *Arelohihi_h,
   double *Arehilohi_h, double *Arelolohi_h,
   double *Arehihilo_h, double *Arelohilo_h,
   double *Arehilolo_h, double *Arelololo_h,
   double *Aimhihihi_h, double *Aimlohihi_h,
   double *Aimhilohi_h, double *Aimlolohi_h,
   double *Aimhihilo_h, double *Aimlohilo_h,
   double *Aimhilolo_h, double *Aimlololo_h,
   double *Arehihihi_d, double *Arelohihi_d,
   double *Arehilohi_d, double *Arelolohi_d,
   double *Arehihilo_d, double *Arelohilo_d,
   double *Arehilolo_d, double *Arelololo_d,
   double *Aimhihihi_d, double *Aimlohihi_d,
   double *Aimhilohi_d, double *Aimlolohi_d,
   double *Aimhihilo_d, double *Aimlohilo_d,
   double *Aimhilolo_d, double *Aimlololo_d,
   double *Vrehihihi_d, double *Vrelohihi_d,
   double *Vrehilohi_d, double *Vrelolohi_d,
   double *Vrehihilo_d, double *Vrelohilo_d,
   double *Vrehilolo_d, double *Vrelololo_d,
   double *Vimhihihi_d, double *Vimlohihi_d,
   double *Vimhilohi_d, double *Vimlolohi_d,
   double *Vimhihilo_d, double *Vimlohilo_d,
   double *Vimhilolo_d, double *Vimlololo_d,
   double *betahihihi_h, double *betalohihi_h,
   double *betahilohi_h, double *betalolohi_h,
   double *betahihilo_h, double *betalohilo_h,
   double *betahilolo_h, double *betalololo_h,
   double *betahihihi_d, double *betalohihi_d,
   double *betahilohi_d, double *betalolohi_d,
   double *betahihilo_d, double *betalohilo_d,
   double *betahilolo_d, double *betalololo_d,
   double *RHdotvrehihihi_h, double *RHdotvrelohihi_h,
   double *RHdotvrehilohi_h, double *RHdotvrelolohi_h,
   double *RHdotvrehihilo_h, double *RHdotvrelohilo_h,
   double *RHdotvrehilolo_h, double *RHdotvrelololo_h,
   double *RHdotvimhihihi_h, double *RHdotvimlohihi_h,
   double *RHdotvimhilohi_h, double *RHdotvimlolohi_h,
   double *RHdotvimhihilo_h, double *RHdotvimlohilo_h,
   double *RHdotvimhilolo_h, double *RHdotvimlololo_h,
   double *RHdotvrehihihi_d, double *RHdotvrelohihi_d,
   double *RHdotvrehilohi_d, double *RHdotvrelolohi_d,
   double *RHdotvrehihilo_d, double *RHdotvrelohilo_d,
   double *RHdotvrehilolo_d, double *RHdotvrelololo_d,
   double *RHdotvimhihihi_d, double *RHdotvimlohihi_d,
   double *RHdotvimhilohi_d, double *RHdotvimlolohi_d,
   double *RHdotvimhihilo_d, double *RHdotvimlohilo_d,
   double *RHdotvimhilolo_d, double *RHdotvimlololo_d,
   double *wrehihihi_h, double *wrelohihi_h,
   double *wrehilohi_h, double *wrelolohi_h,
   double *wrehihilo_h, double *wrelohilo_h,
   double *wrehilolo_h, double *wrelololo_h,
   double *wimhihihi_h, double *wimlohihi_h,
   double *wimhilohi_h, double *wimlolohi_h,
   double *wimhihilo_h, double *wimlohilo_h,
   double *wimhilolo_h, double *wimlololo_h,
   double *wrehihihi_d, double *wrelohihi_d,
   double *wrehilohi_d, double *wrelolohi_d,
   double *wrehihilo_d, double *wrelohilo_d,
   double *wrehilolo_d, double *wrelololo_d,
   double *wimhihihi_d, double *wimlohihi_d,
   double *wimhilohi_d, double *wimlolohi_d,
   double *wimhihilo_d, double *wimlohilo_d,
   double *wimhilolo_d, double *wimlololo_d,
   double *RHvlapms, double *redlapms,
   long long int *add, long long int *mul, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int endcol = (k+1)*szt;     // 1 + last column index in tile
   const int nVrows = nrows - k*szt;          // dimension of V matrix
   const int nhouse = nrows - colidx;  // length of Householder vector
   // total number of entries in R that will be modified
   const int RHoffset = colidx*nrows;
   const int dimRHdotv = endcol - colidx;
   // total number of entries in R that will be modified
   const int sizenum = (nrows - colidx)*dimRHdotv;
   const int nbrblocks = (int) ceil(sizenum/((double) szt));

   if(verbose)
   {
      cout << "-> launching " << nbrblocks << " blocks of " << szt
           << " threads to compute RHdotv ..." << endl;
      cout << "   nhouse : " << nhouse << "  RHoffset : " << RHoffset
           << "  dimRHdotv : " << dimRHdotv << endl;
   }
/*
   cudaEventRecord(start);
   cmplx8_RHdotv<<<nbrblocks,szt>>>
      (nhouse,szt,colidx,RHoffset,dimRHdotv,
       Arehihihi_d,Arelohihi_d,Arehilohi_d,Arelolohi_d,
       Arehihihi_d,Arelohihi_d,Arehilohi_d,Arelolohi_d,
       Aimhihilo_d,Aimlohilo_d,Aimhilolo_d,Aimlololo_d,
       Aimhihilo_d,Aimlohilo_d,Aimhilolo_d,Aimlololo_d,
       &Vrehihihi_d[L*nVrows+L],&Vrelohihi_d[L*nVrows+L],
       &Vrehilohi_d[L*nVrows+L],&Vrelolohi_d[L*nVrows+L],
       &Vrehihilo_d[L*nVrows+L],&Vrelohilo_d[L*nVrows+L],
       &Vrehilolo_d[L*nVrows+L],&Vrelololo_d[L*nVrows+L],
       &Vimhihihi_d[L*nVrows+L],&Vimlohihi_d[L*nVrows+L],
       &Vimhilohi_d[L*nVrows+L],&Vimlolohi_d[L*nVrows+L],
       &Vimhihilo_d[L*nVrows+L],&Vimlohilo_d[L*nVrows+L],
       &Vimhilolo_d[L*nVrows+L],&Vimlololo_d[L*nVrows+L],
       RHdotvrehihihi_d,RHdotvrelohihi_d,RHdotvrehilohi_d,RHdotvrelolohi_d,
       RHdotvrehihilo_d,RHdotvrelohilo_d,RHdotvrehilolo_d,RHdotvrelololo_d,
       RHdotvimhihihi_d,RHdotvimlohihi_d,RHdotvimhilohi_d,RHdotvimlolohi_d,
       RHdotvimhihilo_d,RHdotvimlohilo_d,RHdotvimhilolo_d,RHdotvimlololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *RHvlapms += milliseconds;
 */

/*
   cudaEventRecord(start);
   cmplx8_RHdotvRe<<<nbrblocks,szt>>>
      (nhouse,szt,colidx,RHoffset,dimRHdotv,
       Arehihihi_d,Arelohihi_d,Arehilohi_d,Arelolohi_d,
       Arehihihi_d,Arelohihi_d,Arehilohi_d,Arelolohi_d,
       Aimhihilo_d,Aimlohilo_d,Aimhilolo_d,Aimlololo_d,
       Aimhihilo_d,Aimlohilo_d,Aimhilolo_d,Aimlololo_d,
       &Vrehihihi_d[L*nVrows+L],&Vrelohihi_d[L*nVrows+L],
       &Vrehilohi_d[L*nVrows+L],&Vrelolohi_d[L*nVrows+L],
       &Vrehihilo_d[L*nVrows+L],&Vrelohilo_d[L*nVrows+L],
       &Vrehilolo_d[L*nVrows+L],&Vrelololo_d[L*nVrows+L],
       &Vimhihihi_d[L*nVrows+L],&Vimlohihi_d[L*nVrows+L],
       &Vimhilohi_d[L*nVrows+L],&Vimlolohi_d[L*nVrows+L],
       &Vimhihilo_d[L*nVrows+L],&Vimlohilo_d[L*nVrows+L],
       &Vimhilolo_d[L*nVrows+L],&Vimlololo_d[L*nVrows+L],
       RHdotvrehihihi_d,RHdotvrelohihi_d,RHdotvrehilohi_d,RHdotvrelolohi_d,
       RHdotvrehihilo_d,RHdotvrelohilo_d,RHdotvrehilolo_d,RHdotvrelololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *RHvlapms += milliseconds;
 */

   cudaEventRecord(start);
   cmplx8_RHdotvReRe<<<nbrblocks,szt>>>
      (nhouse,szt,colidx,RHoffset,dimRHdotv,
       Arehihihi_d,Arelohihi_d,Arehilohi_d,Arelolohi_d,
       Arehihilo_d,Arelohilo_d,Arehilolo_d,Arelololo_d,
       &Vrehihihi_d[L*nVrows+L],&Vrelohihi_d[L*nVrows+L],
       &Vrehilohi_d[L*nVrows+L],&Vrelolohi_d[L*nVrows+L],
       &Vrehihilo_d[L*nVrows+L],&Vrelohilo_d[L*nVrows+L],
       &Vrehilolo_d[L*nVrows+L],&Vrelololo_d[L*nVrows+L],
       RHdotvrehihihi_d,RHdotvrelohihi_d,RHdotvrehilohi_d,RHdotvrelolohi_d,
       RHdotvrehihilo_d,RHdotvrelohilo_d,RHdotvrehilolo_d,RHdotvrelololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *RHvlapms += milliseconds;

   cudaEventRecord(start);
   cmplx8_RHdotvImIm<<<nbrblocks,szt>>>
      (nhouse,szt,colidx,RHoffset,dimRHdotv,
       Aimhihihi_d,Aimlohihi_d,Aimhilohi_d,Aimlolohi_d,
       Aimhihilo_d,Aimlohilo_d,Aimhilolo_d,Aimlololo_d,
       &Vimhihihi_d[L*nVrows+L],&Vimlohihi_d[L*nVrows+L],
       &Vimhilohi_d[L*nVrows+L],&Vimlolohi_d[L*nVrows+L],
       &Vimhihilo_d[L*nVrows+L],&Vimlohilo_d[L*nVrows+L],
       &Vimhilolo_d[L*nVrows+L],&Vimlololo_d[L*nVrows+L],
       RHdotvrehihihi_d,RHdotvrelohihi_d,RHdotvrehilohi_d,RHdotvrelolohi_d,
       RHdotvrehihilo_d,RHdotvrelohilo_d,RHdotvrehilolo_d,RHdotvrelololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *RHvlapms += milliseconds;
/*
   cudaEventRecord(start);
   cmplx8_RHdotvIm<<<nbrblocks,szt>>>
      (nhouse,szt,colidx,RHoffset,dimRHdotv,
       Arehihihi_d,Arelohihi_d,Arehilohi_d,Arelolohi_d,
       Arehihihi_d,Arelohihi_d,Arehilohi_d,Arelolohi_d,
       Aimhihilo_d,Aimlohilo_d,Aimhilolo_d,Aimlololo_d,
       Aimhihilo_d,Aimlohilo_d,Aimhilolo_d,Aimlololo_d,
       &Vrehihihi_d[L*nVrows+L],&Vrelohihi_d[L*nVrows+L],
       &Vrehilohi_d[L*nVrows+L],&Vrelolohi_d[L*nVrows+L],
       &Vrehihilo_d[L*nVrows+L],&Vrelohilo_d[L*nVrows+L],
       &Vrehilolo_d[L*nVrows+L],&Vrelololo_d[L*nVrows+L],
       &Vimhihihi_d[L*nVrows+L],&Vimlohihi_d[L*nVrows+L],
       &Vimhilohi_d[L*nVrows+L],&Vimlolohi_d[L*nVrows+L],
       &Vimhihilo_d[L*nVrows+L],&Vimlohilo_d[L*nVrows+L],
       &Vimhilolo_d[L*nVrows+L],&Vimlololo_d[L*nVrows+L],
       RHdotvimhihihi_d,RHdotvimlohihi_d,RHdotvimhilohi_d,RHdotvimlolohi_d,
       RHdotvimhihilo_d,RHdotvimlohilo_d,RHdotvimhilolo_d,RHdotvimlololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *RHvlapms += milliseconds;
 */

   cudaEventRecord(start);
   cmplx8_RHdotvReIm<<<nbrblocks,szt>>>
      (nhouse,szt,colidx,RHoffset,dimRHdotv,
       Arehihihi_d,Arelohihi_d,Arehilohi_d,Arelolohi_d,
       Arehihilo_d,Arelohilo_d,Arehilolo_d,Arelololo_d,
       &Vimhihihi_d[L*nVrows+L],&Vimlohihi_d[L*nVrows+L],
       &Vimhilohi_d[L*nVrows+L],&Vimlolohi_d[L*nVrows+L],
       &Vimhihilo_d[L*nVrows+L],&Vimlohilo_d[L*nVrows+L],
       &Vimhilolo_d[L*nVrows+L],&Vimlololo_d[L*nVrows+L],
       RHdotvimhihihi_d,RHdotvimlohihi_d,RHdotvimhilohi_d,RHdotvimlolohi_d,
       RHdotvimhihilo_d,RHdotvimlohilo_d,RHdotvimhilolo_d,RHdotvimlololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *RHvlapms += milliseconds;

   cudaEventRecord(start);
   cmplx8_RHdotvImRe<<<nbrblocks,szt>>>
      (nhouse,szt,colidx,RHoffset,dimRHdotv,
       Aimhihihi_d,Aimlohihi_d,Aimhilohi_d,Aimlolohi_d,
       Aimhihilo_d,Aimlohilo_d,Aimhilolo_d,Aimlololo_d,
       &Vrehihihi_d[L*nVrows+L],&Vrelohihi_d[L*nVrows+L],
       &Vrehilohi_d[L*nVrows+L],&Vrelolohi_d[L*nVrows+L],
       &Vrehihilo_d[L*nVrows+L],&Vrelohilo_d[L*nVrows+L],
       &Vrehilolo_d[L*nVrows+L],&Vrelololo_d[L*nVrows+L],
       RHdotvimhihihi_d,RHdotvimlohihi_d,RHdotvimhilohi_d,RHdotvimlolohi_d,
       RHdotvimhihilo_d,RHdotvimlohilo_d,RHdotvimhilolo_d,RHdotvimlololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *RHvlapms += milliseconds;

   cudaEventRecord(start);
   cmplx8_sum_betaRHdotv<<<1,dimRHdotv>>>
      (nhouse,
          &betahihihi_d[L],&betalohihi_d[L],&betahilohi_d[L],&betalolohi_d[L],
          &betahihilo_d[L],&betalohilo_d[L],&betahilolo_d[L],&betalololo_d[L],
       RHdotvrehihihi_d,RHdotvrelohihi_d,RHdotvrehilohi_d,RHdotvrelolohi_d,
       RHdotvrehihilo_d,RHdotvrelohilo_d,RHdotvrehilolo_d,RHdotvrelololo_d,
       RHdotvimhihihi_d,RHdotvimlohihi_d,RHdotvimhilohi_d,RHdotvimlolohi_d,
       RHdotvimhihilo_d,RHdotvimlohilo_d,RHdotvimhilolo_d,RHdotvimlololo_d,
            wrehihihi_d,     wrelohihi_d,     wrehilohi_d,     wrelolohi_d,
            wrehihilo_d,     wrelohilo_d,     wrehilolo_d,     wrelololo_d,
            wimhihihi_d,     wimlohihi_d,     wimhilohi_d,     wimlolohi_d,
            wimhihilo_d,     wimlohilo_d,     wimhilolo_d,     wimlololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *RHvlapms += milliseconds;
   flopcount_cmplx_RHdotv(nhouse,szt,add,mul);
   flopcount_cmplx_sum_betaRHdotv(nhouse,dimRHdotv,add,mul);

   if(verbose)
   {
      cout << "-> launching " << nbrblocks << " blocks of " << szt
           << " threads to update " << sizenum << " numbers ..." << endl;
      cout << "   nrows : " << nrows << "  endcol : " << endcol
           << "  szt : " << szt << "  colidx : " << colidx << endl;
   }
/*
   cudaEventRecord(start);
   cmplx8_medium_subvbetaRHv<<<nbrblocks,szt>>>
      (nrows,endcol,szt,colidx,
       Arehihihi_d,Arelohihi_d,Arehilohi_d,Arelolohi_d,
       Arehihilo_d,Arelohilo_d,Arehilolo_d,Arelololo_d,
       Aimhihihi_d,Aimlohihi_d,Aimhilohi_d,Aimlolohi_d,
       Aimhihilo_d,Aimlohilo_d,Aimhilolo_d,Aimlololo_d,
       &Vrehihihi_d[L*nVrows+L],&Vrelohihi_d[L*nVrows+L],
       &Vrehilohi_d[L*nVrows+L],&Vrelolohi_d[L*nVrows+L],
       &Vrehihilo_d[L*nVrows+L],&Vrelohilo_d[L*nVrows+L],
       &Vrehilolo_d[L*nVrows+L],&Vrelololo_d[L*nVrows+L],
       &Vimhihihi_d[L*nVrows+L],&Vimlohihi_d[L*nVrows+L],
       &Vimhilohi_d[L*nVrows+L],&Vimlolohi_d[L*nVrows+L],
       &Vimhihilo_d[L*nVrows+L],&Vimlohilo_d[L*nVrows+L],
       &Vimhilolo_d[L*nVrows+L],&Vimlololo_d[L*nVrows+L],
       &betahihihi_d[L],&betalohihi_d[L],&betahilohi_d[L],&betalolohi_d[L],
       &betahihilo_d[L],&betalohilo_d[L],&betahilolo_d[L],&betalololo_d[L],
       wrehihihi_d,wrelohihi_d,wrehilohi_d,wrelolohi_d,
       wrehihilo_d,wrelohilo_d,wrehilolo_d,wrelololo_d,
       wimhihihi_d,wimlohihi_d,wimhilohi_d,wimlolohi_d,
       wimhihilo_d,wimlohilo_d,wimhilolo_d,wimlololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *redlapms += milliseconds;
 */
/*
   cudaEventRecord(start);
   cmplx8_medium_subvbetaRHvRe<<<nbrblocks,szt>>>
      (nrows,endcol,szt,colidx,
       Arehihihi_d,Arelohihi_d,Arehilohi_d,Arelolohi_d,
       Arehihilo_d,Arelohilo_d,Arehilolo_d,Arelololo_d,
       &Vrehihihi_d[L*nVrows+L],&Vrelohihi_d[L*nVrows+L],
       &Vrehilohi_d[L*nVrows+L],&Vrelolohi_d[L*nVrows+L],
       &Vrehihilo_d[L*nVrows+L],&Vrelohilo_d[L*nVrows+L],
       &Vrehilolo_d[L*nVrows+L],&Vrelololo_d[L*nVrows+L],
       &Vimhihihi_d[L*nVrows+L],&Vimlohihi_d[L*nVrows+L],
       &Vimhilohi_d[L*nVrows+L],&Vimlolohi_d[L*nVrows+L],
       &Vimhihilo_d[L*nVrows+L],&Vimlohilo_d[L*nVrows+L],
       &Vimhilolo_d[L*nVrows+L],&Vimlololo_d[L*nVrows+L],
       wrehihihi_d,wrelohihi_d,wrehilohi_d,wrelolohi_d,
       wrehihilo_d,wrelohilo_d,wrehilolo_d,wrelololo_d,
       wimhihihi_d,wimlohihi_d,wimhilohi_d,wimlolohi_d,
       wimhihilo_d,wimlohilo_d,wimhilolo_d,wimlololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *redlapms += milliseconds;
 */
   cudaEventRecord(start);
   cmplx8_medium_subvbetaRHvReRe<<<nbrblocks,szt>>>
      (nrows,endcol,szt,colidx,
       Arehihihi_d,Arelohihi_d,Arehilohi_d,Arelolohi_d,
       Arehihilo_d,Arelohilo_d,Arehilolo_d,Arelololo_d,
       &Vrehihihi_d[L*nVrows+L],&Vrelohihi_d[L*nVrows+L],
       &Vrehilohi_d[L*nVrows+L],&Vrelolohi_d[L*nVrows+L],
       &Vrehihilo_d[L*nVrows+L],&Vrelohilo_d[L*nVrows+L],
       &Vrehilolo_d[L*nVrows+L],&Vrelololo_d[L*nVrows+L],
       wrehihihi_d,wrelohihi_d,wrehilohi_d,wrelolohi_d,
       wrehihilo_d,wrelohilo_d,wrehilolo_d,wrelololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *redlapms += milliseconds;

   cudaEventRecord(start);
   cmplx8_medium_subvbetaRHvImIm<<<nbrblocks,szt>>>
      (nrows,endcol,szt,colidx,
       Arehihihi_d,Arelohihi_d,Arehilohi_d,Arelolohi_d,
       Arehihilo_d,Arelohilo_d,Arehilolo_d,Arelololo_d,
       &Vimhihihi_d[L*nVrows+L],&Vimlohihi_d[L*nVrows+L],
       &Vimhilohi_d[L*nVrows+L],&Vimlolohi_d[L*nVrows+L],
       &Vimhihilo_d[L*nVrows+L],&Vimlohilo_d[L*nVrows+L],
       &Vimhilolo_d[L*nVrows+L],&Vimlololo_d[L*nVrows+L],
       wimhihihi_d,wimlohihi_d,wimhilohi_d,wimlolohi_d,
       wimhihilo_d,wimlohilo_d,wimhilolo_d,wimlololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *redlapms += milliseconds;
/*
   cudaEventRecord(start);
   cmplx8_medium_subvbetaRHvIm<<<nbrblocks,szt>>>
      (nrows,endcol,szt,colidx,
       Aimhihihi_d,Aimlohihi_d,Aimhilohi_d,Aimlolohi_d,
       Aimhihilo_d,Aimlohilo_d,Aimhilolo_d,Aimlololo_d,
       &Vrehihihi_d[L*nVrows+L],&Vrelohihi_d[L*nVrows+L],
       &Vrehilohi_d[L*nVrows+L],&Vrelolohi_d[L*nVrows+L],
       &Vrehihilo_d[L*nVrows+L],&Vrelohilo_d[L*nVrows+L],
       &Vrehilolo_d[L*nVrows+L],&Vrelololo_d[L*nVrows+L],
       &Vimhihihi_d[L*nVrows+L],&Vimlohihi_d[L*nVrows+L],
       &Vimhilohi_d[L*nVrows+L],&Vimlolohi_d[L*nVrows+L],
       &Vimhihilo_d[L*nVrows+L],&Vimlohilo_d[L*nVrows+L],
       &Vimhilolo_d[L*nVrows+L],&Vimlololo_d[L*nVrows+L],
       wrehihihi_d,wrelohihi_d,wrehilohi_d,wrelolohi_d,
       wrehihilo_d,wrelohilo_d,wrehilolo_d,wrelololo_d,
       wimhihihi_d,wimlohihi_d,wimhilohi_d,wimlolohi_d,
       wimhihilo_d,wimlohilo_d,wimhilolo_d,wimlololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *redlapms += milliseconds;
 */
   cudaEventRecord(start);
   cmplx8_medium_subvbetaRHvImRe<<<nbrblocks,szt>>>
      (nrows,endcol,szt,colidx,
       Aimhihihi_d,Aimlohihi_d,Aimhilohi_d,Aimlolohi_d,
       Aimhihilo_d,Aimlohilo_d,Aimhilolo_d,Aimlololo_d,
       &Vimhihihi_d[L*nVrows+L],&Vimlohihi_d[L*nVrows+L],
       &Vimhilohi_d[L*nVrows+L],&Vimlolohi_d[L*nVrows+L],
       &Vimhihilo_d[L*nVrows+L],&Vimlohilo_d[L*nVrows+L],
       &Vimhilolo_d[L*nVrows+L],&Vimlololo_d[L*nVrows+L],
       wrehihihi_d,wrelohihi_d,wrehilohi_d,wrelolohi_d,
       wrehihilo_d,wrelohilo_d,wrehilolo_d,wrelololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *redlapms += milliseconds;
   flopcount_cmplx_medium_subvbetaRHv(nrows,endcol,szt,colidx,add,mul);

   cudaEventRecord(start);
   cmplx8_medium_subvbetaRHvReIm<<<nbrblocks,szt>>>
      (nrows,endcol,szt,colidx,
       Aimhihihi_d,Aimlohihi_d,Aimhilohi_d,Aimlolohi_d,
       Aimhihilo_d,Aimlohilo_d,Aimhilolo_d,Aimlololo_d,
       &Vrehihihi_d[L*nVrows+L],&Vrelohihi_d[L*nVrows+L],
       &Vrehilohi_d[L*nVrows+L],&Vrelolohi_d[L*nVrows+L],
       &Vrehihilo_d[L*nVrows+L],&Vrelohilo_d[L*nVrows+L],
       &Vrehilolo_d[L*nVrows+L],&Vrelololo_d[L*nVrows+L],
       wimhihihi_d,wimlohihi_d,wimhilohi_d,wimlolohi_d,
       wimhihilo_d,wimlohilo_d,wimhilolo_d,wimlololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *redlapms += milliseconds;

   if(verbose)
   {
      const int dim = nrows*ncols;
      const size_t sznum = dim*sizeof(double);
      const size_t szbRHv = dimRHdotv*sizeof(double);
      const size_t szRHdotv = nVrows*szbRHv;

      cudaMemcpy(RHdotvrehihihi_h,RHdotvrehihihi_d,szRHdotv,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(RHdotvrelohihi_h,RHdotvrelohihi_d,szRHdotv,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(RHdotvrehilohi_h,RHdotvrehilohi_d,szRHdotv,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(RHdotvrelolohi_h,RHdotvrelolohi_d,szRHdotv,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(RHdotvrehihilo_h,RHdotvrehihilo_d,szRHdotv,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(RHdotvrelohilo_h,RHdotvrelohilo_d,szRHdotv,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(RHdotvrehilolo_h,RHdotvrehilolo_d,szRHdotv,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(RHdotvrelololo_h,RHdotvrelololo_d,szRHdotv,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(RHdotvimhihihi_h,RHdotvimhihihi_d,szRHdotv,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(RHdotvimlohihi_h,RHdotvimlohihi_d,szRHdotv,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(RHdotvimhilohi_h,RHdotvimhilohi_d,szRHdotv,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(RHdotvimlolohi_h,RHdotvimlolohi_d,szRHdotv,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(RHdotvimhihilo_h,RHdotvimhihilo_d,szRHdotv,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(RHdotvimlohilo_h,RHdotvimlohilo_d,szRHdotv,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(RHdotvimhilolo_h,RHdotvimhilolo_d,szRHdotv,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(RHdotvimlololo_h,RHdotvimlololo_d,szRHdotv,
                 cudaMemcpyDeviceToHost);

      cout << "the matrix R^H dot v : " << endl;
      int ix = 0;
      for(int i=0; i<endcol-colidx; i++)
      {
         for(int j=0; j<nhouse; j++)      // must use nhouse
         {
            cout << "RHdotv[" << i << "][" << j << "]re : "
                 << RHdotvrehihihi_h[ix] << "  "
                 << RHdotvrelohihi_h[ix] << endl
                 << "                 "
                 << RHdotvrehilohi_h[ix] << "  "
                 << RHdotvrelolohi_h[ix] << endl
                 << "                 "
                 << RHdotvrehihilo_h[ix] << "  "
                 << RHdotvrelohilo_h[ix] << endl
                 << "                 "
                 << RHdotvrehilolo_h[ix] << "  "
                 << RHdotvrelololo_h[ix] << endl;
            cout << "RHdotv[" << i << "][" << j << "]im : "
                 << RHdotvimhihihi_h[ix] << "  "
                 << RHdotvimlohihi_h[ix] << endl
                 << "                 "
                 << RHdotvimhilohi_h[ix] << "  "
                 << RHdotvimlolohi_h[ix] << endl
                 << "                 "
                 << RHdotvimhihilo_h[ix] << "  "
                 << RHdotvimlohilo_h[ix] << endl
                 << "                 "
                 << RHdotvimhilolo_h[ix] << "  "
                 << RHdotvimlololo_h[ix] << endl;
            ix = ix + 1;
         }
      }
      cudaMemcpy(wrehihihi_h,wrehihihi_d,szbRHv,cudaMemcpyDeviceToHost);
      cudaMemcpy(wrelohihi_h,wrelohihi_d,szbRHv,cudaMemcpyDeviceToHost);
      cudaMemcpy(wrehilohi_h,wrehilohi_d,szbRHv,cudaMemcpyDeviceToHost);
      cudaMemcpy(wrelolohi_h,wrelolohi_d,szbRHv,cudaMemcpyDeviceToHost);
      cudaMemcpy(wrehihilo_h,wrehihilo_d,szbRHv,cudaMemcpyDeviceToHost);
      cudaMemcpy(wrelohilo_h,wrelohilo_d,szbRHv,cudaMemcpyDeviceToHost);
      cudaMemcpy(wrehilolo_h,wrehilolo_d,szbRHv,cudaMemcpyDeviceToHost);
      cudaMemcpy(wrelololo_h,wrelololo_d,szbRHv,cudaMemcpyDeviceToHost);
      cudaMemcpy(wimhihihi_h,wimhihihi_d,szbRHv,cudaMemcpyDeviceToHost);
      cudaMemcpy(wimlohihi_h,wimlohihi_d,szbRHv,cudaMemcpyDeviceToHost);
      cudaMemcpy(wimhilohi_h,wimhilohi_d,szbRHv,cudaMemcpyDeviceToHost);
      cudaMemcpy(wimlolohi_h,wimlolohi_d,szbRHv,cudaMemcpyDeviceToHost);
      cudaMemcpy(wimhihilo_h,wimhihilo_d,szbRHv,cudaMemcpyDeviceToHost);
      cudaMemcpy(wimlohilo_h,wimlohilo_d,szbRHv,cudaMemcpyDeviceToHost);
      cudaMemcpy(wimhilolo_h,wimhilolo_d,szbRHv,cudaMemcpyDeviceToHost);
      cudaMemcpy(wimlololo_h,wimlololo_d,szbRHv,cudaMemcpyDeviceToHost);

      cout << "the vector w = beta*R^H*v : " << endl;
      for(int i=0; i<endcol-colidx; i++)
      {
         cout << "w[" << i << "]re : "
              << wrehihihi_h[i] << "  " << wimlohihi_h[i] << endl
              << "         "
              << wrehilohi_h[i] << "  " << wimlolohi_h[i] << endl
              << "         "
              << wrehihilo_h[i] << "  " << wimlohilo_h[i] << endl
              << "         "
              << wrehilolo_h[i] << "  " << wimlololo_h[i] << endl;
         cout << "w[" << i << "]im : "
              << wimhihihi_h[i] << "  " << wimlohihi_h[i] << endl
              << "         "
              << wimhilohi_h[i] << "  " << wimlolohi_h[i] << endl
              << "         "
              << wimhihilo_h[i] << "  " << wimlohilo_h[i] << endl
              << "         "
              << wimhilolo_h[i] << "  " << wimlololo_h[i] << endl;
      }
      cudaMemcpy(Arehihihi_h,Arehihihi_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Arelohihi_h,Arelohihi_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Arehilohi_h,Arehilohi_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Arelolohi_h,Arelolohi_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Arehihilo_h,Arehihilo_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Arelohilo_h,Arelohilo_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Arehilolo_h,Arehilolo_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Arelololo_h,Arelololo_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Aimhihihi_h,Aimhihihi_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Aimlohihi_h,Aimlohihi_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Aimhilohi_h,Aimhilohi_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Aimlolohi_h,Aimlolohi_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Aimhihilo_h,Aimhihilo_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Aimlohilo_h,Aimlohilo_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Aimhilolo_h,Aimhilolo_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Aimlololo_h,Aimlololo_d,sznum,cudaMemcpyDeviceToHost);

      cout << "the matrix after the update :" << endl;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
         {
            cout << "A_d[" << i << "][" << j << "]re : "
                 << Arehihihi_h[j*nrows+i] << "  "
                 << Arelohihi_h[j*nrows+i] << endl
                 << "              "
                 << Arehilohi_h[j*nrows+i] << "  "
                 << Arelolohi_h[j*nrows+i] << endl
                 << "              "
                 << Arehihilo_h[j*nrows+i] << "  "
                 << Arelohilo_h[j*nrows+i] << endl
                 << "              "
                 << Arehilolo_h[j*nrows+i] << "  "
                 << Arelololo_h[j*nrows+i] << endl;
            cout << "A_d[" << i << "][" << j << "]im : "
                 << Aimhihihi_h[j*nrows+i] << "  "
                 << Aimlohihi_h[j*nrows+i] << endl
                 << "              "
                 << Aimhilohi_h[j*nrows+i] << "  "
                 << Aimlolohi_h[j*nrows+i] << endl
                 << "              "
                 << Aimhihilo_h[j*nrows+i] << "  "
                 << Aimlohilo_h[j*nrows+i] << endl
                 << "              "
                 << Aimhilolo_h[j*nrows+i] << "  "
                 << Aimlololo_h[j*nrows+i] << endl;
         }
   }
}

void GPU_dbl8_medium_VB_to_W
 ( int nrows, int ncols, int szt, int idx,
   double *Vhihihi_h, double *Vlohihi_h, double *Vhilohi_h, double *Vlolohi_h,
   double *Vhihilo_h, double *Vlohilo_h, double *Vhilolo_h, double *Vlololo_h,
   double *Vhihihi_d, double *Vlohihi_d, double *Vhilohi_d, double *Vlolohi_d,
   double *Vhihilo_d, double *Vlohilo_d, double *Vhilolo_d, double *Vlololo_d,
   double *Whihihi_h, double *Wlohihi_h, double *Whilohi_h, double *Wlolohi_h,
   double *Whihilo_h, double *Wlohilo_h, double *Whilolo_h, double *Wlololo_h,
   double *Whihihi_d, double *Wlohihi_d, double *Whilohi_d, double *Wlolohi_d,
   double *Whihilo_d, double *Wlohilo_d, double *Whilolo_d, double *Wlololo_d,
   double *WYThihihi_h, double *WYTlohihi_h,
   double *WYThilohi_h, double *WYTlolohi_h,
   double *WYThihilo_h, double *WYTlohilo_h,
   double *WYThilolo_h, double *WYTlololo_h,
   double *WYThihihi_d, double *WYTlohihi_d,
   double *WYThilohi_d, double *WYTlolohi_d,
   double *WYThihilo_d, double *WYTlohilo_d,
   double *WYThilolo_d, double *WYTlololo_d,
   double *betahihihi_h, double *betalohihi_h,
   double *betahilohi_h, double *betalolohi_h,
   double *betahihilo_h, double *betalohilo_h,
   double *betahilolo_h, double *betalololo_h,
   double *betahihihi_d, double *betalohihi_d,
   double *betahilohi_d, double *betalolohi_d,
   double *betahihilo_d, double *betalohilo_d,
   double *betahilolo_d, double *betalololo_d,
   double *lapms, long long int *add, long long int *mul, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int rowdim = nrows - idx*szt;
   const int nbrblocks1 = (int) ceil(rowdim/((double) szt));

   cudaEventRecord(start);
   dbl8_beta_times_V<<<nbrblocks1,szt>>>
      (rowdim,szt,betahihihi_d,betalohihi_d,betahilohi_d,betalolohi_d,
                  betahihilo_d,betalohilo_d,betahilolo_d,betalololo_d,
                     Vhihihi_d,   Vlohihi_d,   Vhilohi_d,   Vlolohi_d,
                     Vhihilo_d,   Vlohilo_d,   Vhilolo_d,   Vlololo_d,
                     Whihihi_d,   Wlohihi_d,   Whilohi_d,   Wlolohi_d,
                     Whihilo_d,   Wlohilo_d,   Whilolo_d,   Wlololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_dbl_beta_times_V(rowdim,mul);

   const int nbrblocks2 = (int) ceil(rowdim*rowdim/((double) szt));

   cudaEventRecord(start);
   dbl8_initialize_WYT<<<nbrblocks2,szt>>>
      (rowdim,szt,  Vhihihi_d,  Vlohihi_d,  Vhilohi_d,  Vlolohi_d,
                    Vhihilo_d,  Vlohilo_d,  Vhilolo_d,  Vlololo_d,
                    Whihihi_d,  Wlohihi_d,  Whilohi_d,  Wlolohi_d,
                    Whihilo_d,  Wlohilo_d,  Whilolo_d,  Wlololo_d,
                  WYThihihi_d,WYTlohihi_d,WYThilohi_d,WYTlolohi_d,
                  WYThihilo_d,WYTlohilo_d,WYThilolo_d,WYTlololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_dbl_initialize_WYT(rowdim,mul);

   for(int j=1; j<szt; j++)
   {
      cudaEventRecord(start);
      dbl8_beta_next_W<<<nbrblocks1,szt>>>
         (rowdim,szt,
          &betahihihi_d[j],&betalohihi_d[j],&betahilohi_d[j],&betalolohi_d[j],
          &betahihilo_d[j],&betalohilo_d[j],&betahilolo_d[j],&betalololo_d[j],
          &Vhihihi_d[j*rowdim],&Vlohihi_d[j*rowdim],
          &Vhilohi_d[j*rowdim],&Vlolohi_d[j*rowdim],
          &Vhihilo_d[j*rowdim],&Vlohilo_d[j*rowdim],
          &Vhilolo_d[j*rowdim],&Vlololo_d[j*rowdim],
          &Whihihi_d[j*rowdim],&Wlohihi_d[j*rowdim],
          &Whilohi_d[j*rowdim],&Wlolohi_d[j*rowdim],
          &Whihilo_d[j*rowdim],&Wlohilo_d[j*rowdim],
          &Whilolo_d[j*rowdim],&Wlololo_d[j*rowdim],
          WYThihihi_d,WYTlohihi_d,WYThilohi_d,WYTlolohi_d,
          WYThihilo_d,WYTlohilo_d,WYThilolo_d,WYTlololo_d);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds,start,stop);
      *lapms += milliseconds;
      flopcount_dbl_beta_next_W(rowdim,add,mul);

      cudaEventRecord(start);
      dbl8_update_WYT<<<nbrblocks2,szt>>>
         (rowdim,szt,
          &Vhihihi_d[j*rowdim],&Vlohihi_d[j*rowdim],
          &Vhilohi_d[j*rowdim],&Vlolohi_d[j*rowdim],
          &Vhihilo_d[j*rowdim],&Vlohilo_d[j*rowdim],
          &Vhilolo_d[j*rowdim],&Vlololo_d[j*rowdim],
          &Whihihi_d[j*rowdim],&Wlohihi_d[j*rowdim],
          &Whilohi_d[j*rowdim],&Wlolohi_d[j*rowdim],
          &Whihilo_d[j*rowdim],&Wlohilo_d[j*rowdim],
          &Whilolo_d[j*rowdim],&Wlololo_d[j*rowdim],
          WYThihihi_d,WYTlohihi_d,WYThilohi_d,WYTlolohi_d,
          WYThihilo_d,WYTlohilo_d,WYThilolo_d,WYTlololo_d);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds,start,stop);
      *lapms += milliseconds;
      flopcount_dbl_update_WYT(rowdim,add,mul);
   }
   if(verbose)
   {
      const size_t szbeta = szt*sizeof(double);
      const size_t szhouse = rowdim*sizeof(double);
      const size_t szVandW = szt*szhouse;
      const size_t szmat = rowdim*rowdim*sizeof(double);

      cudaMemcpy(betahihihi_h,betahihihi_d,szbeta,cudaMemcpyDeviceToHost);
      cudaMemcpy(betalohihi_h,betalohihi_d,szbeta,cudaMemcpyDeviceToHost);
      cudaMemcpy(betahilohi_h,betahilohi_d,szbeta,cudaMemcpyDeviceToHost);
      cudaMemcpy(betalolohi_h,betalolohi_d,szbeta,cudaMemcpyDeviceToHost);
      cudaMemcpy(betahihilo_h,betahihilo_d,szbeta,cudaMemcpyDeviceToHost);
      cudaMemcpy(betalohilo_h,betalohilo_d,szbeta,cudaMemcpyDeviceToHost);
      cudaMemcpy(betahilolo_h,betahilolo_d,szbeta,cudaMemcpyDeviceToHost);
      cudaMemcpy(betalololo_h,betalololo_d,szbeta,cudaMemcpyDeviceToHost);

      cout << "the betas :" << endl;
      for(int j=0; j<szt; j++)
         cout << "beta[" << j << "] : "
              << betahihihi_h[j] << "  " << betalohihi_h[j] << endl
              << "          "
              << betahilohi_h[j] << "  " << betalolohi_h[j] << endl
              << "          "
              << betahihilo_h[j] << "  " << betalohilo_h[j] << endl
              << "          "
              << betahilolo_h[j] << "  " << betalololo_h[j] << endl;

      cudaMemcpy(Vhihihi_h,Vhihihi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Vlohihi_h,Vlohihi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Vhilohi_h,Vhilohi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Vlolohi_h,Vlolohi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Vhihilo_h,Vhihilo_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Vlohilo_h,Vlohilo_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Vhilolo_h,Vhilolo_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Vlololo_h,Vlololo_d,szVandW,cudaMemcpyDeviceToHost);

      cout << "the columns of the V matrix :" << endl;
      int ix = 0;
      for(int j=0; j<szt; j++) 
         for(int i=0; i<rowdim; i++) 
         {
            cout << "V[" << i << "][" << j << "] : "
                 << Vhihihi_h[ix] << "  " << Vlohihi_h[ix] << endl
                 << "          "
                 << Vhilohi_h[ix] << "  " << Vlolohi_h[ix] << endl
                 << "          "
                 << Vhihilo_h[ix] << "  " << Vlohilo_h[ix] << endl
                 << "          "
                 << Vhilolo_h[ix] << "  " << Vlololo_h[ix] << endl;
            ix = ix + 1;
         }

      cudaMemcpy(Whihihi_h,Whihihi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Wlohihi_h,Wlohihi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Whilohi_h,Whilohi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Wlolohi_h,Wlolohi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Whihilo_h,Whihilo_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Wlohilo_h,Wlohilo_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Whilolo_h,Whilolo_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Wlololo_h,Wlololo_d,szVandW,cudaMemcpyDeviceToHost);

      cout << "the columns of the W matrix :" << endl;
      ix = 0;
      for(int j=0; j<szt; j++) 
         for(int i=0; i<rowdim; i++) 
         {
            cout << "W[" << i << "][" << j << "] : "
                 << Whihihi_h[ix] << "  " << Wlohihi_h[ix] << endl
                 << "          "
                 << Whilohi_h[ix] << "  " << Wlolohi_h[ix] << endl
                 << "          "
                 << Whihilo_h[ix] << "  " << Wlohilo_h[ix] << endl
                 << "          "
                 << Whilolo_h[ix] << "  " << Wlololo_h[ix] << endl;
            ix = ix + 1;
         }

      cudaMemcpy(WYThihihi_h,WYThihihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYTlohihi_h,WYTlohihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYThilohi_h,WYThilohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYTlolohi_h,WYTlolohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYThihilo_h,WYThihilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYTlohilo_h,WYTlohilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYThilolo_h,WYThilolo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYTlololo_h,WYTlololo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the WYT matrix :" << endl;
      ix = 0;
      for(int i=0; i<rowdim; i++) 
         for(int j=0; j<rowdim; j++) 
         {
            cout << "WYT[" << i << "][" << j << "] : "
                 << WYThihihi_h[ix] << "  " << WYTlohihi_h[ix] << endl
                 << "            "
                 << WYThilohi_h[ix] << "  " << WYTlolohi_h[ix] << endl
                 << "            "
                 << WYThihilo_h[ix] << "  " << WYTlohilo_h[ix] << endl
                 << "            "
                 << WYThilolo_h[ix] << "  " << WYTlololo_h[ix] << endl;
            ix = ix + 1;
         }
   }
}

void GPU_cmplx8_medium_VB_to_W
 ( int nrows, int ncols, int szt, int idx,
   double *Vrehihihi_h, double *Vrelohihi_h,
   double *Vrehilohi_h, double *Vrelolohi_h,
   double *Vrehihilo_h, double *Vrelohilo_h,
   double *Vrehilolo_h, double *Vrelololo_h,
   double *Vimhihihi_h, double *Vimlohihi_h,
   double *Vimhilohi_h, double *Vimlolohi_h,
   double *Vimhihilo_h, double *Vimlohilo_h,
   double *Vimhilolo_h, double *Vimlololo_h,
   double *Vrehihihi_d, double *Vrelohihi_d,
   double *Vrehilohi_d, double *Vrelolohi_d,
   double *Vrehihilo_d, double *Vrelohilo_d, 
   double *Vrehilolo_d, double *Vrelololo_d,
   double *Vimhihihi_d, double *Vimlohihi_d,
   double *Vimhilohi_d, double *Vimlolohi_d,
   double *Vimhihilo_d, double *Vimlohilo_d,
   double *Vimhilolo_d, double *Vimlololo_d,
   double *Wrehihihi_h, double *Wrelohihi_h,
   double *Wrehilohi_h, double *Wrelolohi_h,
   double *Wrehihilo_h, double *Wrelohilo_h,
   double *Wrehilolo_h, double *Wrelololo_h,
   double *Wimhihihi_h, double *Wimlohihi_h,
   double *Wimhilohi_h, double *Wimlolohi_h,
   double *Wimhihilo_h, double *Wimlohilo_h,
   double *Wimhilolo_h, double *Wimlololo_h,
   double *Wrehihihi_d, double *Wrelohihi_d,
   double *Wrehilohi_d, double *Wrelolohi_d,
   double *Wrehihilo_d, double *Wrelohilo_d,
   double *Wrehilolo_d, double *Wrelololo_d,
   double *Wimhihihi_d, double *Wimlohihi_d,
   double *Wimhilohi_d, double *Wimlolohi_d,
   double *Wimhihilo_d, double *Wimlohilo_d,
   double *Wimhilolo_d, double *Wimlololo_d,
   double *WYHrehihihi_h, double *WYHrelohihi_h,
   double *WYHrehilohi_h, double *WYHrelolohi_h,
   double *WYHrehihilo_h, double *WYHrelohilo_h,
   double *WYHrehilolo_h, double *WYHrelololo_h,
   double *WYHimhihihi_h, double *WYHimlohihi_h,
   double *WYHimhilohi_h, double *WYHimlolohi_h,
   double *WYHimhihilo_h, double *WYHimlohilo_h,
   double *WYHimhilolo_h, double *WYHimlololo_h,
   double *WYHrehihihi_d, double *WYHrelohihi_d,
   double *WYHrehilohi_d, double *WYHrelolohi_d,
   double *WYHrehihilo_d, double *WYHrelohilo_d,
   double *WYHrehilolo_d, double *WYHrelololo_d,
   double *WYHimhihihi_d, double *WYHimlohihi_d,
   double *WYHimhilohi_d, double *WYHimlolohi_d,
   double *WYHimhihilo_d, double *WYHimlohilo_d,
   double *WYHimhilolo_d, double *WYHimlololo_d,
   double *betahihihi_h, double *betalohihi_h,
   double *betahilohi_h, double *betalolohi_h,
   double *betahihilo_h, double *betalohilo_h,
   double *betahilolo_h, double *betalololo_h,
   double *betahihihi_d, double *betalohihi_d,
   double *betahilohi_d, double *betalolohi_d,
   double *betahihilo_d, double *betalohilo_d,
   double *betahilolo_d, double *betalololo_d,
   double *lapms, long long int *add, long long int *mul, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int rowdim = nrows - idx*szt;
   const int nbrblocks1 = (int) ceil(rowdim/((double) szt));

   cudaEventRecord(start);
   cmplx8_beta_times_V<<<nbrblocks1,szt>>>
      (rowdim,szt,betahihihi_d,betalohihi_d,betahilohi_d,betalolohi_d,
                  betahihilo_d,betalohilo_d,betahilolo_d,betalololo_d,
                   Vrehihihi_d, Vrelohihi_d, Vrehilohi_d, Vrelolohi_d,
                   Vrehihilo_d, Vrelohilo_d, Vrehilolo_d, Vrelololo_d,
                   Vimhihihi_d, Vimlohihi_d, Vimhilohi_d, Vimlolohi_d,
                   Vimhihilo_d, Vimlohilo_d, Vimhilolo_d, Vimlololo_d,
                   Wrehihihi_d, Wrelohihi_d, Wrehilohi_d, Wrelolohi_d,
                   Wrehihilo_d, Wrelohilo_d, Wrehilolo_d, Wrelololo_d,
                   Wimhihihi_d, Wimlohihi_d, Wimhilohi_d, Wimlolohi_d,
                   Wimhihilo_d, Wimlohilo_d, Wimhilolo_d, Wimlololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_cmplx_beta_times_V(rowdim,mul);

   const int nbrblocks2 = (int) ceil(rowdim*rowdim/((double) szt));

   cudaEventRecord(start);
   cmplx8_initialize_WYH<<<nbrblocks2,szt>>>
      (rowdim,szt,  Vrehihihi_d,  Vrelohihi_d,  Vrehilohi_d,  Vrelolohi_d,
                    Vrehihilo_d,  Vrelohilo_d,  Vrehilolo_d,  Vrelololo_d,
                    Vimhihihi_d,  Vimlohihi_d,  Vimhilohi_d,  Vimlolohi_d,
                    Vimhihilo_d,  Vimlohilo_d,  Vimhilolo_d,  Vimlololo_d,
                    Wrehihihi_d,  Wrelohihi_d,  Wrehilohi_d,  Wrelolohi_d,
                    Wrehihilo_d,  Wrelohilo_d,  Wrehilolo_d,  Wrelololo_d,
                    Wimhihihi_d,  Wimlohihi_d,  Wimhilohi_d,  Wimlolohi_d,
                    Wimhihilo_d,  Wimlohilo_d,  Wimhilolo_d,  Wimlololo_d,
                  WYHrehihihi_d,WYHrelohihi_d,WYHrehilohi_d,WYHrelolohi_d,
                  WYHrehihilo_d,WYHrelohilo_d,WYHrehilolo_d,WYHrelololo_d,
                  WYHimhihihi_d,WYHimlohihi_d,WYHimhilohi_d,WYHimlolohi_d,
                  WYHimhihilo_d,WYHimlohilo_d,WYHimhilolo_d,WYHimlololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_cmplx_initialize_WYH(rowdim,add,mul);

   for(int j=1; j<szt; j++)
   {
      cudaEventRecord(start);
      cmplx8_beta_next_W<<<nbrblocks1,szt>>>
         (rowdim,szt,
          &betahihihi_d[j],&betalohihi_d[j],&betahilohi_d[j],&betalolohi_d[j],
          &betahihilo_d[j],&betalohilo_d[j],&betahilolo_d[j],&betalololo_d[j],
          &Vrehihihi_d[j*rowdim],&Vrelohihi_d[j*rowdim],
          &Vrehilohi_d[j*rowdim],&Vrelolohi_d[j*rowdim],
          &Vrehihilo_d[j*rowdim],&Vrelohilo_d[j*rowdim],
          &Vrehilolo_d[j*rowdim],&Vrelololo_d[j*rowdim],
          &Vimhihihi_d[j*rowdim],&Vimlohihi_d[j*rowdim],
          &Vimhilohi_d[j*rowdim],&Vimlolohi_d[j*rowdim],
          &Vimhihilo_d[j*rowdim],&Vimlohilo_d[j*rowdim],
          &Vimhilolo_d[j*rowdim],&Vimlololo_d[j*rowdim],
          &Wrehihihi_d[j*rowdim],&Wrelohihi_d[j*rowdim],
          &Wrehilohi_d[j*rowdim],&Wrelolohi_d[j*rowdim],
          &Wrehihilo_d[j*rowdim],&Wrelohilo_d[j*rowdim],
          &Wrehilolo_d[j*rowdim],&Wrelololo_d[j*rowdim],
          &Wimhihihi_d[j*rowdim],&Wimlohihi_d[j*rowdim],
          &Wimhilohi_d[j*rowdim],&Wimlolohi_d[j*rowdim],
          &Wimhihilo_d[j*rowdim],&Wimlohilo_d[j*rowdim],
          &Wimhilolo_d[j*rowdim],&Wimlololo_d[j*rowdim],
          WYHrehihihi_d,WYHrelohihi_d,WYHrehilohi_d,WYHrelolohi_d,
          WYHrehihilo_d,WYHrelohilo_d,WYHrehilolo_d,WYHrelololo_d,
          WYHimhihihi_d,WYHimlohihi_d,WYHimhilohi_d,WYHimlolohi_d,
          WYHimhihilo_d,WYHimlohilo_d,WYHimhilolo_d,WYHimlololo_d);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds,start,stop);
      *lapms += milliseconds;
      flopcount_cmplx_beta_next_W(rowdim,add,mul);

      cudaEventRecord(start);
      cmplx8_update_WYH<<<nbrblocks2,szt>>>
         (rowdim,szt,&Vrehihihi_d[j*rowdim],&Vrelohihi_d[j*rowdim],
                     &Vrehilohi_d[j*rowdim],&Vrelolohi_d[j*rowdim],
                     &Vrehihilo_d[j*rowdim],&Vrelohilo_d[j*rowdim],
                     &Vrehilolo_d[j*rowdim],&Vrelololo_d[j*rowdim],
                     &Vimhihihi_d[j*rowdim],&Vimlohihi_d[j*rowdim],
                     &Vimhilohi_d[j*rowdim],&Vimlolohi_d[j*rowdim],
                     &Vimhihilo_d[j*rowdim],&Vimlohilo_d[j*rowdim],
                     &Vimhilolo_d[j*rowdim],&Vimlololo_d[j*rowdim],
                     &Wrehihihi_d[j*rowdim],&Wrelohihi_d[j*rowdim],
                     &Wrehilohi_d[j*rowdim],&Wrelolohi_d[j*rowdim],
                     &Wrehihilo_d[j*rowdim],&Wrelohilo_d[j*rowdim],
                     &Wrehilolo_d[j*rowdim],&Wrelololo_d[j*rowdim],
                     &Wimhihihi_d[j*rowdim],&Wimlohihi_d[j*rowdim],
                     &Wimhilohi_d[j*rowdim],&Wimlolohi_d[j*rowdim],
                     &Wimhihilo_d[j*rowdim],&Wimlohilo_d[j*rowdim],
                     &Wimhilolo_d[j*rowdim],&Wimlololo_d[j*rowdim],
                     WYHrehihihi_d,WYHrelohihi_d,WYHrehilohi_d,WYHrelolohi_d,
                     WYHrehihilo_d,WYHrelohilo_d,WYHrehilolo_d,WYHrelololo_d,
                     WYHimhihihi_d,WYHimlohihi_d,WYHimhilohi_d,WYHimlolohi_d,
                     WYHimhihilo_d,WYHimlohilo_d,WYHimhilolo_d,WYHimlololo_d);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds,start,stop);
      *lapms += milliseconds;
      flopcount_cmplx_update_WYH(rowdim,add,mul);
   }
   if(verbose)
   {
      const size_t szbeta = szt*sizeof(double);
      const size_t szhouse = rowdim*sizeof(double);
      const size_t szVandW = szt*szhouse;
      const size_t szmat = rowdim*rowdim*sizeof(double);

      cudaMemcpy(betahihihi_h,betahihihi_d,szbeta,cudaMemcpyDeviceToHost);
      cudaMemcpy(betalohihi_h,betalohihi_d,szbeta,cudaMemcpyDeviceToHost);
      cudaMemcpy(betahilohi_h,betahilohi_d,szbeta,cudaMemcpyDeviceToHost);
      cudaMemcpy(betalolohi_h,betalolohi_d,szbeta,cudaMemcpyDeviceToHost);
      cudaMemcpy(betahihilo_h,betahihilo_d,szbeta,cudaMemcpyDeviceToHost);
      cudaMemcpy(betalohilo_h,betalohilo_d,szbeta,cudaMemcpyDeviceToHost);
      cudaMemcpy(betahilolo_h,betahilolo_d,szbeta,cudaMemcpyDeviceToHost);
      cudaMemcpy(betalololo_h,betalololo_d,szbeta,cudaMemcpyDeviceToHost);

      cout << "the betas :" << endl;
      for(int j=0; j<szt; j++)
         cout << "beta[" << j << "] : "
              << betahihihi_h[j] << "  " << betalohihi_h[j] << endl
              << "          "
              << betahilohi_h[j] << "  " << betalolohi_h[j] << endl
              << "          "
              << betahihilo_h[j] << "  " << betalohilo_h[j] << endl
              << "          "
              << betahilolo_h[j] << "  " << betalololo_h[j] << endl;

      cudaMemcpy(Vrehihihi_h,Vrehihihi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Vrelohihi_h,Vrelohihi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Vrehilohi_h,Vrehilohi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Vrelolohi_h,Vrelolohi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Vrehihilo_h,Vrehihilo_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Vrelohilo_h,Vrelohilo_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Vrehilolo_h,Vrehilolo_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Vrelololo_h,Vrelololo_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Vimhihihi_h,Vimhihihi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Vimlohihi_h,Vimlohihi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Vimhilohi_h,Vimhilohi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Vimlolohi_h,Vimlolohi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Vimhihilo_h,Vimhihilo_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Vimlohilo_h,Vimlohilo_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Vimhilolo_h,Vimhilolo_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Vimlololo_h,Vimlololo_d,szVandW,cudaMemcpyDeviceToHost);

      cout << "the columns of the V matrix :" << endl;
      int ix = 0;
      for(int j=0; j<szt; j++) 
         for(int i=0; i<rowdim; i++) 
         {
            cout << "V[" << i << "][" << j << "]re : "
                 << Vrehihihi_h[ix] << "  " << Vrelohihi_h[ix] << endl
                 << "            "
                 << Vrehilohi_h[ix] << "  " << Vrelolohi_h[ix] << endl
                 << "            "
                 << Vrehihilo_h[ix] << "  " << Vrelohilo_h[ix] << endl
                 << "            "
                 << Vrehilolo_h[ix] << "  " << Vrelololo_h[ix] << endl;
            cout << "V[" << i << "][" << j << "]im : "
                 << Vimhihihi_h[ix] << "  " << Vimlohihi_h[ix] << endl
                 << "            "
                 << Vimhilohi_h[ix] << "  " << Vimlolohi_h[ix] << endl
                 << "            "
                 << Vimhihilo_h[ix] << "  " << Vimlohilo_h[ix] << endl
                 << "            "
                 << Vimhilolo_h[ix] << "  " << Vimlololo_h[ix] << endl;
            ix = ix + 1;
         }

      cudaMemcpy(Wrehihihi_h,Wrehihihi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Wrelohihi_h,Wrelohihi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Wrehilohi_h,Wrehilohi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Wrelolohi_h,Wrelolohi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Wrehihilo_h,Wrehihilo_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Wrelohilo_h,Wrelohilo_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Wrehilolo_h,Wrehilolo_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Wrelololo_h,Wrelololo_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Wimhihihi_h,Wimhihihi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Wimlohihi_h,Wimlohihi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Wimhilohi_h,Wimhilohi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Wimlolohi_h,Wimlolohi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Wimhihilo_h,Wimhihilo_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Wimlohilo_h,Wimlohilo_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Wimhilolo_h,Wimhilolo_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Wimlololo_h,Wimlololo_d,szVandW,cudaMemcpyDeviceToHost);

      cout << "the columns of the W matrix :" << endl;
      ix = 0;
      for(int j=0; j<szt; j++) 
         for(int i=0; i<rowdim; i++) 
         {
            cout << "W[" << i << "][" << j << "]re : "
                 << Wrehihihi_h[ix] << "  " << Wrelohihi_h[ix] << endl
                 << "            "
                 << Wrehilohi_h[ix] << "  " << Wrelolohi_h[ix] << endl
                 << "            "
                 << Wrehihilo_h[ix] << "  " << Wrelohilo_h[ix] << endl
                 << "            "
                 << Wrehilolo_h[ix] << "  " << Wrelololo_h[ix] << endl;
            cout << "W[" << i << "][" << j << "]im : "
                 << Wimhihihi_h[ix] << "  " << Wimlohihi_h[ix] << endl
                 << "            "
                 << Wimhilohi_h[ix] << "  " << Wimlolohi_h[ix] << endl
                 << "            "
                 << Wimhihilo_h[ix] << "  " << Wimlohilo_h[ix] << endl
                 << "            "
                 << Wimhilolo_h[ix] << "  " << Wimlololo_h[ix] << endl;
            ix = ix + 1;
         }

      cudaMemcpy(WYHrehihihi_h,WYHrehihihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYHrelohihi_h,WYHrelohihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYHrehilohi_h,WYHrehilohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYHrelolohi_h,WYHrelolohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYHrehihilo_h,WYHrehihilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYHrelohilo_h,WYHrelohilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYHrehilolo_h,WYHrehilolo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYHrelololo_h,WYHrelololo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYHimhihihi_h,WYHimhihihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYHimlohihi_h,WYHimlohihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYHimhilohi_h,WYHimhilohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYHimlolohi_h,WYHimlolohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYHimhihilo_h,WYHimhihilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYHimlohilo_h,WYHimlohilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYHimhilolo_h,WYHimhilolo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYHimlololo_h,WYHimlololo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the WYT matrix :" << endl;
      ix = 0;
      for(int i=0; i<rowdim; i++) 
         for(int j=0; j<rowdim; j++) 
         {
            cout << "WYH[" << i << "][" << j << "]re : "
                 << WYHrehihihi_h[ix] << "  " << WYHrelohihi_h[ix] << endl
                 << "              "
                 << WYHrehilohi_h[ix] << "  " << WYHrelolohi_h[ix] << endl
                 << "              "
                 << WYHrehihilo_h[ix] << "  " << WYHrelohilo_h[ix] << endl
                 << "              "
                 << WYHrehilolo_h[ix] << "  " << WYHrelololo_h[ix] << endl;
            cout << "WYH[" << i << "][" << j << "]im : "
                 << WYHimhihihi_h[ix] << "  " << WYHimlohihi_h[ix] << endl
                 << "              "
                 << WYHimhilohi_h[ix] << "  " << WYHimlolohi_h[ix] << endl
                 << "              "
                 << WYHimhihilo_h[ix] << "  " << WYHimlohilo_h[ix] << endl
                 << "              "
                 << WYHimhilolo_h[ix] << "  " << WYHimlololo_h[ix] << endl;
            ix = ix + 1;
         }
   }
}

void GPU_dbl8_small_WYT
 ( int nrows, int szt,
   double *Whihihi_d, double *Wlohihi_d,
   double *Whilohi_d, double *Wlolohi_d,
   double *Whihilo_d, double *Wlohilo_d,
   double *Whilolo_d, double *Wlololo_d,
   double *Yhihihi_d, double *Ylohihi_d,
   double *Yhilohi_d, double *Ylolohi_d,
   double *Yhihilo_d, double *Ylohilo_d,
   double *Yhilolo_d, double *Ylololo_d,
   double *WYThihihi_d, double *WYTlohihi_d,
   double *WYThilohi_d, double *WYTlolohi_d,
   double *WYThihilo_d, double *WYTlohilo_d,
   double *WYThilolo_d, double *WYTlololo_d,
   double *WYThihihi_h, double *WYTlohihi_h,
   double *WYThilohi_h, double *WYTlolohi_h,
   double *WYThihilo_h, double *WYTlohilo_h,
   double *WYThilolo_h, double *WYTlololo_h,
   double *lapms, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int nbrblocks = (int) ceil(nrows*nrows/((double) szt));

   cudaEventRecord(start);
   dbl8_small_WYT<<<nbrblocks,szt>>>
      (nrows,szt,  Whihihi_d,  Wlohihi_d,  Whilohi_d,  Wlolohi_d,
                   Whihilo_d,  Wlohilo_d,  Whilolo_d,  Wlololo_d,
                   Yhihihi_d,  Ylohihi_d,  Yhilohi_d,  Ylolohi_d,
                   Yhihilo_d,  Ylohilo_d,  Yhilolo_d,  Ylololo_d,
                 WYThihihi_d,WYTlohihi_d,WYThilohi_d,WYTlolohi_d,
                 WYThihilo_d,WYTlohilo_d,WYThilolo_d,WYTlololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   // flopcount_dbl_small_WYT(nrows,szt,add,mul);

   if(verbose)
   {
      const size_t szmat = nrows*nrows*sizeof(double);

      cudaMemcpy(WYThihihi_h,WYThihihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYTlohihi_h,WYTlohihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYThilohi_h,WYThilohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYTlolohi_h,WYTlolohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYThihilo_h,WYThihilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYTlohilo_h,WYTlohilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYThilolo_h,WYThilolo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYTlololo_h,WYTlololo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the WYT matrix :" << endl;
      int ix = 0;
      for(int i=0; i<nrows; i++) 
         for(int j=0; j<nrows; j++) 
         {
            cout << "WYT[" << i << "][" << j << "] : "
                 << WYThihihi_h[ix] << "  " << WYTlohihi_h[ix] << endl
                 << "            "
                 << WYThilohi_h[ix] << "  " << WYTlolohi_h[ix] << endl
                 << "            "
                 << WYThihilo_h[ix] << "  " << WYTlohilo_h[ix] << endl
                 << "            "
                 << WYThilolo_h[ix] << "  " << WYTlololo_h[ix] << endl;
            ix = ix + 1;
         }
   }
}

void GPU_cmplx8_small_WYH
 ( int nrows, int szt,
   double *Wrehihihi_d, double *Wrelohihi_d,
   double *Wrehilohi_d, double *Wrelolohi_d,
   double *Wrehihilo_d, double *Wrelohilo_d,
   double *Wrehilolo_d, double *Wrelololo_d,
   double *Wimhihihi_d, double *Wimlohihi_d,
   double *Wimhilohi_d, double *Wimlolohi_d,
   double *Wimhihilo_d, double *Wimlohilo_d,
   double *Wimhilolo_d, double *Wimlololo_d,
   double *Yrehihihi_d, double *Yrelohihi_d,
   double *Yrehilohi_d, double *Yrelolohi_d,
   double *Yrehihilo_d, double *Yrelohilo_d,
   double *Yrehilolo_d, double *Yrelololo_d,
   double *Yimhihihi_d, double *Yimlohihi_d,
   double *Yimhilohi_d, double *Yimlolohi_d,
   double *Yimhihilo_d, double *Yimlohilo_d,
   double *Yimhilolo_d, double *Yimlololo_d,
   double *WYTrehihihi_d, double *WYTrelohihi_d,
   double *WYTrehilohi_d, double *WYTrelolohi_d,
   double *WYTrehihilo_d, double *WYTrelohilo_d,
   double *WYTrehilolo_d, double *WYTrelololo_d,
   double *WYTimhihihi_d, double *WYTimlohihi_d,
   double *WYTimhilohi_d, double *WYTimlolohi_d,
   double *WYTimhihilo_d, double *WYTimlohilo_d,
   double *WYTimhilolo_d, double *WYTimlololo_d,
   double *WYTrehihihi_h, double *WYTrelohihi_h,
   double *WYTrehilohi_h, double *WYTrelolohi_h,
   double *WYTrehihilo_h, double *WYTrelohilo_h,
   double *WYTrehilolo_h, double *WYTrelololo_h,
   double *WYTimhihihi_h, double *WYTimlohihi_h,
   double *WYTimhilohi_h, double *WYTimlolohi_h,
   double *WYTimhihilo_h, double *WYTimlohilo_h,
   double *WYTimhilolo_h, double *WYTimlololo_h,
   double *lapms, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int nbrblocks = (int) ceil(nrows*nrows/((double) szt));

   cudaEventRecord(start);
   cmplx8_small_WYH<<<nbrblocks,szt>>>
      (nrows,szt, Wrehihihi_d,  Wrelohihi_d,  Wrehilohi_d,  Wrelolohi_d,
                  Wrehihilo_d,  Wrelohilo_d,  Wrehilolo_d,  Wrelololo_d,
                  Wimhihihi_d,  Wimlohihi_d,  Wimhilohi_d,  Wimlolohi_d,
                  Wimhihilo_d,  Wimlohilo_d,  Wimhilolo_d,  Wimlololo_d,
                  Yrehihihi_d,  Yrelohihi_d,  Yrehilohi_d,  Yrelolohi_d,
                  Yrehihilo_d,  Yrelohilo_d,  Yrehilolo_d,  Yrelololo_d,
                  Yimhihihi_d,  Yimlohihi_d,  Yimhilohi_d,  Yimlolohi_d,
                  Yimhihilo_d,  Yimlohilo_d,  Yimhilolo_d,  Yimlololo_d,
                WYTrehihihi_d,WYTrelohihi_d,WYTrehilohi_d,WYTrelolohi_d,
                WYTrehihilo_d,WYTrelohilo_d,WYTrehilolo_d,WYTrelololo_d,
                WYTimhihihi_d,WYTimlohihi_d,WYTimhilohi_d,WYTimlolohi_d,
                WYTimhihilo_d,WYTimlohilo_d,WYTimhilolo_d,WYTimlololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   // flopcount_cmplx_small_WYH(nrows,szt,add,mul);

   if(verbose)
   {
      const size_t szmat = nrows*nrows*sizeof(double);

      cudaMemcpy(WYTrehihihi_h,WYTrehihihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYTrelohihi_h,WYTrelohihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYTrehilohi_h,WYTrehilohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYTrelolohi_h,WYTrelolohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYTrehihilo_h,WYTrehihilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYTrelohilo_h,WYTrelohilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYTrehilolo_h,WYTrehilolo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYTrelololo_h,WYTrelololo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYTimhihihi_h,WYTimhihihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYTimlohihi_h,WYTimlohihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYTimhilohi_h,WYTimhilohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYTimlolohi_h,WYTimlolohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYTimhihilo_h,WYTimhihilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYTimlohilo_h,WYTimlohilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYTimhilolo_h,WYTimhilolo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYTimlololo_h,WYTimlololo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the WYT matrix :" << endl;
      int ix = 0;
      for(int i=0; i<nrows; i++) 
         for(int j=0; j<nrows; j++) 
         {
            cout << "WYT[" << i << "][" << j << "]re : "
                 << WYTrehihihi_h[ix] << "  " << WYTrelohihi_h[ix] << endl
                 << "              "
                 << WYTrehilohi_h[ix] << "  " << WYTrelolohi_h[ix] << endl
                 << "              "
                 << WYTrehihilo_h[ix] << "  " << WYTrelohilo_h[ix] << endl
                 << "              "
                 << WYTrehilolo_h[ix] << "  " << WYTrelololo_h[ix] << endl;
            cout << "WYT[" << i << "][" << j << "]im : "
                 << WYTimhihihi_h[ix] << "  " << WYTimlohihi_h[ix] << endl
                 << "              "
                 << WYTimhilohi_h[ix] << "  " << WYTimlolohi_h[ix] << endl
                 << "              "
                 << WYTimhihilo_h[ix] << "  " << WYTimlohilo_h[ix] << endl
                 << "              "
                 << WYTimhilolo_h[ix] << "  " << WYTimlololo_h[ix] << endl;
            ix = ix + 1;
         }
   }
}

void GPU_dbl8_small_YWT
 ( int nrows, int szt, int idx,
   double *Yhihihi_d, double *Ylohihi_d, double *Yhilohi_d, double *Ylolohi_d,
   double *Yhihilo_d, double *Ylohilo_d, double *Yhilolo_d, double *Ylololo_d,
   double *Whihihi_d, double *Wlohihi_d, double *Whilohi_d, double *Wlolohi_d,
   double *Whihilo_d, double *Wlohilo_d, double *Whilolo_d, double *Wlololo_d,
   double *YWThihihi_d, double *YWTlohihi_d,
   double *YWThilohi_d, double *YWTlolohi_d,
   double *YWThihilo_d, double *YWTlohilo_d,
   double *YWThilolo_d, double *YWTlololo_d,
   double *YWThihihi_h, double *YWTlohihi_h,
   double *YWThilohi_h, double *YWTlolohi_h,
   double *YWThihilo_h, double *YWTlohilo_h,
   double *YWThilolo_h, double *YWTlololo_h,
   double *lapms, long long int *add, long long int *mul, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int rowdim = nrows - idx*szt;
   int nbrblocks = (int) ceil(rowdim*rowdim/((double) szt));

   cudaEventRecord(start);
   dbl8_small_WYT<<<nbrblocks,szt>>>
      (rowdim,szt,  Yhihihi_d,  Ylohihi_d,  Yhilohi_d,  Ylolohi_d,
                    Yhihilo_d,  Ylohilo_d,  Yhilolo_d,  Ylololo_d,
                    Whihihi_d,  Wlohihi_d,  Whilohi_d,  Wlolohi_d,
                    Whihilo_d,  Wlohilo_d,  Whilolo_d,  Wlololo_d,
                  YWThihihi_d,YWTlohihi_d,YWThilohi_d,YWTlolohi_d,
                  YWThihilo_d,YWTlohilo_d,YWThilolo_d,YWTlololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_dbl_small_WYT(rowdim,szt,add,mul);

   if(verbose)
   {
      const size_t szmat = rowdim*rowdim*sizeof(double);

      cudaMemcpy(YWThihihi_h,YWThihihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTlohihi_h,YWTlohihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWThilohi_h,YWThilohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTlolohi_h,YWTlolohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWThihilo_h,YWThihilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTlohilo_h,YWTlohilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWThilolo_h,YWThilolo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTlololo_h,YWTlololo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the YWT matrix :" << endl;
      int ix = 0;
      for(int i=0; i<rowdim; i++) 
         for(int j=0; j<rowdim; j++) 
         {
            cout << "YWT[" << i << "][" << j << "] : "
                 << YWThihihi_h[ix] << "  " << YWTlohihi_h[ix] << endl
                 << "            "
                 << YWThilohi_h[ix] << "  " << YWTlolohi_h[ix] << endl
                 << "            "
                 << YWThihilo_h[ix] << "  " << YWTlohilo_h[ix] << endl
                 << "            "
                 << YWThilolo_h[ix] << "  " << YWTlololo_h[ix] << endl;
            ix = ix + 1;
         }
   }
}

void GPU_cmplx8_small_YWH
 ( int nrows, int szt, int idx,
   double *Yrehihihi_d, double *Yrelohihi_d,
   double *Yrehilohi_d, double *Yrelolohi_d,
   double *Yrehihilo_d, double *Yrelohilo_d,
   double *Yrehilolo_d, double *Yrelololo_d,
   double *Yimhihihi_d, double *Yimlohihi_d,
   double *Yimhilohi_d, double *Yimlolohi_d,
   double *Yimhihilo_d, double *Yimlohilo_d,
   double *Yimhilolo_d, double *Yimlololo_d,
   double *Wrehihihi_d, double *Wrelohihi_d,
   double *Wrehilohi_d, double *Wrelolohi_d,
   double *Wrehihilo_d, double *Wrelohilo_d,
   double *Wrehilolo_d, double *Wrelololo_d,
   double *Wimhihihi_d, double *Wimlohihi_d,
   double *Wimhilohi_d, double *Wimlolohi_d,
   double *Wimhihilo_d, double *Wimlohilo_d,
   double *Wimhilolo_d, double *Wimlololo_d,
   double *YWTrehihihi_d, double *YWTrelohihi_d,
   double *YWTrehilohi_d, double *YWTrelolohi_d,
   double *YWTrehihilo_d, double *YWTrelohilo_d,
   double *YWTrehilolo_d, double *YWTrelololo_d,
   double *YWTimhihihi_d, double *YWTimlohihi_d,
   double *YWTimhilohi_d, double *YWTimlolohi_d,
   double *YWTimhihilo_d, double *YWTimlohilo_d,
   double *YWTimhilolo_d, double *YWTimlololo_d,
   double *YWTrehihihi_h, double *YWTrelohihi_h,
   double *YWTrehilohi_h, double *YWTrelolohi_h,
   double *YWTrehihilo_h, double *YWTrelohilo_h,
   double *YWTrehilolo_h, double *YWTrelololo_h,
   double *YWTimhihihi_h, double *YWTimlohihi_h,
   double *YWTimhilohi_h, double *YWTimlolohi_h,
   double *YWTimhihilo_h, double *YWTimlohilo_h,
   double *YWTimhilolo_h, double *YWTimlololo_h,
   double *lapms, long long int *add, long long int *mul, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int rowdim = nrows - idx*szt;
   int nbrblocks = (int) ceil(rowdim*rowdim/((double) szt));

   cudaEventRecord(start);
   cmplx8_small_WYH<<<nbrblocks,szt>>>
      (rowdim,szt,  Yrehihihi_d,  Yrelohihi_d,  Yrehilohi_d,  Yrelolohi_d,
                    Yrehihilo_d,  Yrelohilo_d,  Yrehilolo_d,  Yrelololo_d,
                    Yimhihihi_d,  Yimlohihi_d,  Yimhilohi_d,  Yimlolohi_d,
                    Yimhihilo_d,  Yimlohilo_d,  Yimhilolo_d,  Yimlololo_d,
                    Wrehihihi_d,  Wrelohihi_d,  Wrehilohi_d,  Wrelolohi_d,
                    Wrehihilo_d,  Wrelohilo_d,  Wrehilolo_d,  Wrelololo_d,
                    Wimhihihi_d,  Wimlohihi_d,  Wimhilohi_d,  Wimlolohi_d,
                    Wimhihilo_d,  Wimlohilo_d,  Wimhilolo_d,  Wimlololo_d,
                  YWTrehihihi_d,YWTrelohihi_d,YWTrehilohi_d,YWTrelolohi_d,
                  YWTrehihilo_d,YWTrelohilo_d,YWTrehilolo_d,YWTrelololo_d,
                  YWTimhihihi_d,YWTimlohihi_d,YWTimhilohi_d,YWTimlolohi_d,
                  YWTimhihilo_d,YWTimlohilo_d,YWTimhilolo_d,YWTimlololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_cmplx_small_WYH(rowdim,szt,add,mul);

   if(verbose)
   {
      const size_t szmat = rowdim*rowdim*sizeof(double);

      cudaMemcpy(YWTrehihihi_h,YWTrehihihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTrelohihi_h,YWTrelohihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTrehilohi_h,YWTrehilohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTrelolohi_h,YWTrelolohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTrehihilo_h,YWTrehihilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTrelohilo_h,YWTrelohilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTrehilolo_h,YWTrehilolo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTrelololo_h,YWTrelololo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTimhihihi_h,YWTimhihihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTimlohihi_h,YWTimlohihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTimhilohi_h,YWTimhilohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTimlolohi_h,YWTimlolohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTimhihilo_h,YWTimhihilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTimlohilo_h,YWTimlohilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTimhilolo_h,YWTimhilolo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTimlololo_h,YWTimlololo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the YWT matrix :" << endl;
      int ix = 0;
      for(int i=0; i<rowdim; i++) 
         for(int j=0; j<rowdim; j++) 
         {
            cout << "YWT[" << i << "][" << j << "]re : "
                 << YWTrehihihi_h[ix] << "  " << YWTrelohihi_h[ix] << endl
                 << "              "
                 << YWTrehilohi_h[ix] << "  " << YWTrelolohi_h[ix] << endl
                 << "              "
                 << YWTrehihilo_h[ix] << "  " << YWTrelohilo_h[ix] << endl
                 << "              "
                 << YWTrehilolo_h[ix] << "  " << YWTrelololo_h[ix] << endl;
            cout << "YWT[" << i << "][" << j << "]im : "
                 << YWTimhihihi_h[ix] << "  " << YWTimlohihi_h[ix] << endl
                 << "              "
                 << YWTimhilohi_h[ix] << "  " << YWTimlolohi_h[ix] << endl
                 << "              "
                 << YWTimhihilo_h[ix] << "  " << YWTimlohilo_h[ix] << endl
                 << "              "
                 << YWTimhilolo_h[ix] << "  " << YWTimlololo_h[ix] << endl;
            ix = ix + 1;
         }
   }
}

void GPU_dbl8_small_QWYT
 ( int dim, int szt, int idx,
   double *Qhihihi_d, double *Qlohihi_d, double *Qhilohi_d, double *Qlolohi_d,
   double *Qhihilo_d, double *Qlohilo_d, double *Qhilolo_d, double *Qlololo_d,
   double *WYThihihi_d, double *WYTlohihi_d,
   double *WYThilohi_d, double *WYTlolohi_d,
   double *WYThihilo_d, double *WYTlohilo_d,
   double *WYThilolo_d, double *WYTlololo_d,
   double *QWYThihihi_d, double *QWYTlohihi_d,
   double *QWYThilohi_d, double *QWYTlolohi_d,
   double *QWYThihilo_d, double *QWYTlohilo_d,
   double *QWYThilolo_d, double *QWYTlololo_d,
   double *QWYThihihi_h, double *QWYTlohihi_h,
   double *QWYThilohi_h, double *QWYTlolohi_h,
   double *QWYThihilo_h, double *QWYTlohilo_h,
   double *QWYThilolo_h, double *QWYTlololo_h,
   double *Qhihihi_h, double *Qlohihi_h, double *Qhilohi_h, double *Qlolohi_h,
   double *Qhihilo_h, double *Qlohilo_h, double *Qhilolo_h, double *Qlololo_h,
   double *lapms, long long int *add, long long int *mul, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int coloff = idx*szt;
   const int rowdim = dim - coloff;
   const int nbrblocks = (int) ceil(dim*rowdim/((double) szt));

   if(verbose)
   {
      const size_t szmat = dim*dim*sizeof(double);

      cudaMemcpy(Qhihihi_h,Qhihihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qlohihi_h,Qlohihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qhilohi_h,Qhilohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qlolohi_h,Qlolohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qhihilo_h,Qhihilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qlohilo_h,Qlohilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qhilolo_h,Qhilolo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qlololo_h,Qlololo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the Q matrix :" << endl;
      int ix = 0;
      for(int i=0; i<dim; i++) 
         for(int j=0; j<dim; j++) 
         {
            cout << "Q[" << i << "][" << j << "] : "
                 << Qhihihi_h[ix] << "  " << Qlohihi_h[ix] << endl
                 << "          "
                 << Qhilohi_h[ix] << "  " << Qlolohi_h[ix] << endl
                 << "          "
                 << Qhihilo_h[ix] << "  " << Qlohilo_h[ix] << endl
                 << "          "
                 << Qhilolo_h[ix] << "  " << Qlololo_h[ix] << endl;
            ix = ix + 1;
         }
   }

   cudaEventRecord(start);
   dbl8_small_QWYT<<<nbrblocks,szt>>>
      (dim,rowdim,szt,coloff,
          Qhihihi_d,   Qlohihi_d,   Qhilohi_d,   Qlolohi_d,
          Qhihilo_d,   Qlohilo_d,   Qhilolo_d,   Qlololo_d,
        WYThihihi_d, WYTlohihi_d, WYThilohi_d, WYTlolohi_d,
        WYThihilo_d, WYTlohilo_d, WYThilolo_d, WYTlololo_d,
       QWYThihihi_d,QWYTlohihi_d,QWYThilohi_d,QWYTlolohi_d,
       QWYThihilo_d,QWYTlohilo_d,QWYThilolo_d,QWYTlololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_dbl_small_QWYT(dim,rowdim,szt,coloff,add,mul);

   if(verbose)
   {
      const size_t szmat = dim*rowdim*sizeof(double);

      cudaMemcpy(QWYThihihi_h,QWYThihihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(QWYTlohihi_h,QWYTlohihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(QWYThilohi_h,QWYThilohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(QWYTlolohi_h,QWYTlolohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(QWYThihilo_h,QWYThihilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(QWYTlohilo_h,QWYTlohilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(QWYThilolo_h,QWYThilolo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(QWYTlololo_h,QWYTlololo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the QWYT matrix :" << endl;
      int ix = 0;
      for(int i=0; i<dim; i++) 
         for(int j=0; j<rowdim; j++) 
         {
            cout << "QWYT[" << i << "][" << j << "] : "
                 << QWYThihihi_h[ix] << "  " << QWYTlohihi_h[ix] << endl
                 << "             "
                 << QWYThilohi_h[ix] << "  " << QWYTlolohi_h[ix] << endl
                 << "             "
                 << QWYThihilo_h[ix] << "  " << QWYTlohilo_h[ix] << endl
                 << "             "
                 << QWYThilolo_h[ix] << "  " << QWYTlololo_h[ix] << endl;
            ix = ix + 1;
         }
   }
}

void GPU_cmplx8_small_QWYH
 ( int dim, int szt, int idx,
   double *Qrehihihi_d, double *Qrelohihi_d,
   double *Qrehilohi_d, double *Qrelolohi_d,
   double *Qrehihilo_d, double *Qrelohilo_d,
   double *Qrehilolo_d, double *Qrelololo_d,
   double *Qimhihihi_d, double *Qimlohihi_d,
   double *Qimhilohi_d, double *Qimlolohi_d,
   double *Qimhihilo_d, double *Qimlohilo_d,
   double *Qimhilolo_d, double *Qimlololo_d,
   double *WYTrehihihi_d, double *WYTrelohihi_d,
   double *WYTrehilohi_d, double *WYTrelolohi_d,
   double *WYTrehihilo_d, double *WYTrelohilo_d,
   double *WYTrehilolo_d, double *WYTrelololo_d,
   double *WYTimhihihi_d, double *WYTimlohihi_d,
   double *WYTimhilohi_d, double *WYTimlolohi_d,
   double *WYTimhihilo_d, double *WYTimlohilo_d,
   double *WYTimhilolo_d, double *WYTimlololo_d,
   double *QWYTrehihihi_d, double *QWYTrelohihi_d,
   double *QWYTrehilohi_d, double *QWYTrelolohi_d,
   double *QWYTrehihilo_d, double *QWYTrelohilo_d,
   double *QWYTrehilolo_d, double *QWYTrelololo_d,
   double *QWYTimhihihi_d, double *QWYTimlohihi_d,
   double *QWYTimhilohi_d, double *QWYTimlolohi_d,
   double *QWYTimhihilo_d, double *QWYTimlohilo_d,
   double *QWYTimhilolo_d, double *QWYTimlololo_d,
   double *QWYTrehihihi_h, double *QWYTrelohihi_h,
   double *QWYTrehilohi_h, double *QWYTrelolohi_h,
   double *QWYTrehihilo_h, double *QWYTrelohilo_h,
   double *QWYTrehilolo_h, double *QWYTrelololo_h,
   double *QWYTimhihihi_h, double *QWYTimlohihi_h,
   double *QWYTimhilohi_h, double *QWYTimlolohi_h,
   double *QWYTimhihilo_h, double *QWYTimlohilo_h,
   double *QWYTimhilolo_h, double *QWYTimlololo_h,
   double *Qrehihihi_h, double *Qrelohihi_h,
   double *Qrehilohi_h, double *Qrelolohi_h,
   double *Qrehihilo_h, double *Qrelohilo_h,
   double *Qrehilolo_h, double *Qrelololo_h,
   double *Qimhihihi_h, double *Qimlohihi_h,
   double *Qimhilohi_h, double *Qimlolohi_h,
   double *Qimhihilo_h, double *Qimlohilo_h,
   double *Qimhilolo_h, double *Qimlololo_h,
   double *lapms, long long int *add, long long int *mul, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int coloff = idx*szt;
   const int rowdim = dim - coloff;
   const int nbrblocks = (int) ceil(dim*rowdim/((double) szt));

   if(verbose)
   {
      const size_t szmat = dim*dim*sizeof(double);

      cudaMemcpy(Qrehihihi_h,Qrehihihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qrelohihi_h,Qrelohihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qrehilohi_h,Qrehilohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qrelolohi_h,Qrelolohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qrehihilo_h,Qrehihilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qrelohilo_h,Qrelohilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qrehilolo_h,Qrehilolo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qrelololo_h,Qrelololo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qimhihihi_h,Qimhihihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qimlohihi_h,Qimlohihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qimhilohi_h,Qimhilohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qimlolohi_h,Qimlolohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qimhihilo_h,Qimhihilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qimlohilo_h,Qimlohilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qimhilolo_h,Qimhilolo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qimlololo_h,Qimlololo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the Q matrix :" << endl;
      int ix = 0;
      for(int i=0; i<dim; i++) 
         for(int j=0; j<dim; j++) 
         {
            cout << "Q[" << i << "][" << j << "]re : "
                 << Qrehihihi_h[ix] << "  " << Qrelohihi_h[ix] << endl
                 << "          "
                 << Qrehilohi_h[ix] << "  " << Qrelolohi_h[ix] << endl
                 << "          "
                 << Qrehihilo_h[ix] << "  " << Qrelohilo_h[ix] << endl
                 << "          "
                 << Qrehilolo_h[ix] << "  " << Qrelololo_h[ix] << endl;
            cout << "Q[" << i << "][" << j << "] : "
                 << Qimhihihi_h[ix] << "  " << Qimlohihi_h[ix] << endl
                 << "          "
                 << Qimhilohi_h[ix] << "  " << Qimlolohi_h[ix] << endl
                 << "          "
                 << Qimhihilo_h[ix] << "  " << Qimlohilo_h[ix] << endl
                 << "          "
                 << Qimhilolo_h[ix] << "  " << Qimlololo_h[ix] << endl;
            ix = ix + 1;
         }
   }

   cudaEventRecord(start);
   cmplx8_small_QWYH<<<nbrblocks,szt>>>
      (dim,rowdim,szt,coloff,
          Qrehihihi_d,   Qrelohihi_d,   Qrehilohi_d,   Qrelolohi_d,
          Qrehihilo_d,   Qrelohilo_d,   Qrehilolo_d,   Qrelololo_d,
          Qimhihihi_d,   Qimlohihi_d,   Qimhilohi_d,   Qimlolohi_d,
          Qimhihilo_d,   Qimlohilo_d,   Qimhilolo_d,   Qimlololo_d,
        WYTrehihihi_d, WYTrelohihi_d, WYTrehilohi_d, WYTrelolohi_d,
        WYTrehihilo_d, WYTrelohilo_d, WYTrehilolo_d, WYTrelololo_d,
        WYTimhihihi_d, WYTimlohihi_d, WYTimhilohi_d, WYTimlolohi_d,
        WYTimhihilo_d, WYTimlohilo_d, WYTimhilolo_d, WYTimlololo_d,
       QWYTrehihihi_d,QWYTrelohihi_d,QWYTrehilohi_d,QWYTrelolohi_d,
       QWYTrehihilo_d,QWYTrelohilo_d,QWYTrehilolo_d,QWYTrelololo_d,
       QWYTimhihihi_d,QWYTimlohihi_d,QWYTimhilohi_d,QWYTimlolohi_d,
       QWYTimhihilo_d,QWYTimlohilo_d,QWYTimhilolo_d,QWYTimlololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_cmplx_small_QWYH(dim,rowdim,szt,coloff,add,mul);

   if(verbose)
   {
      const size_t szmat = dim*rowdim*sizeof(double);

      cudaMemcpy(QWYTrehihihi_h,QWYTrehihihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(QWYTrelohihi_h,QWYTrelohihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(QWYTrehilohi_h,QWYTrehilohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(QWYTrelolohi_h,QWYTrelolohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(QWYTrehihilo_h,QWYTrehihilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(QWYTrelohilo_h,QWYTrelohilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(QWYTrehilolo_h,QWYTrehilolo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(QWYTrelololo_h,QWYTrelololo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(QWYTimhihihi_h,QWYTimhihihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(QWYTimlohihi_h,QWYTimlohihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(QWYTimhilohi_h,QWYTimhilohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(QWYTimlolohi_h,QWYTimlolohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(QWYTimhihilo_h,QWYTimhihilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(QWYTimlohilo_h,QWYTimlohilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(QWYTimhilolo_h,QWYTimhilolo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(QWYTimlololo_h,QWYTimlololo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the QWYT matrix :" << endl;
      int ix = 0;
      for(int i=0; i<dim; i++) 
         for(int j=0; j<rowdim; j++) 
         {
            cout << "QWYT[" << i << "][" << j << "]re : "
                 << QWYTrehihihi_h[ix] << "  " << QWYTrelohihi_h[ix] << endl
                 << "               "
                 << QWYTrehilohi_h[ix] << "  " << QWYTrelolohi_h[ix] << endl
                 << "               "
                 << QWYTrehihilo_h[ix] << "  " << QWYTrelohilo_h[ix] << endl
                 << "               "
                 << QWYTrehilolo_h[ix] << "  " << QWYTrelololo_h[ix] << endl;
            cout << "QWYT[" << i << "][" << j << "]im : "
                 << QWYTimhihihi_h[ix] << "  " << QWYTimlohihi_h[ix] << endl
                 << "               "
                 << QWYTimhilohi_h[ix] << "  " << QWYTimlolohi_h[ix] << endl
                 << "               "
                 << QWYTimhihilo_h[ix] << "  " << QWYTimlohilo_h[ix] << endl
                 << "               "
                 << QWYTimhilolo_h[ix] << "  " << QWYTimlololo_h[ix] << endl;
            ix = ix + 1;
         }
   }
}

void GPU_dbl8_small_YWTC
 ( int nrows, int ncols, int szt, int idx,
   double *YWThihihi_d, double *YWTlohihi_d,
   double *YWThilohi_d, double *YWTlolohi_d,
   double *YWThihilo_d, double *YWTlohilo_d,
   double *YWThilolo_d, double *YWTlololo_d,
   double *Chihihi_d, double *Clohihi_d, double *Chilohi_d, double *Clolohi_d,
   double *Chihilo_d, double *Clohilo_d, double *Chilolo_d, double *Clololo_d,
   double *YWTChihihi_d, double *YWTClohihi_d,
   double *YWTChilohi_d, double *YWTClolohi_d,
   double *YWTChihilo_d, double *YWTClohilo_d,
   double *YWTChilolo_d, double *YWTClololo_d,
   double *YWTChihihi_h, double *YWTClohihi_h,
   double *YWTChilohi_h, double *YWTClolohi_h,
   double *YWTChihilo_h, double *YWTClohilo_h,
   double *YWTChilolo_h, double *YWTClololo_h,
   double *lapms, long long int *add, long long int *mul, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int rowoff = idx*szt;
   const int rowdim = nrows - rowoff;
   const int coloff = (idx+1)*szt;
   const int coldim = ncols - coloff;
   const int nbrblocks = (int) ceil(rowdim*coldim/((double) szt));

   if(verbose)
   {
      cout << "in GPU_dbl2_small_YWTC ..." << endl;
      cout << "-> nrows : " << nrows
           << "  ncols : " << ncols
           << "  szt : " << szt
           << "  idx : " << idx << endl;
      cout << "   rowdim : " << rowdim
           << "  coldim : " << coldim
           << "  rowoff : " << rowoff
           << "  coloff : " << coloff
           << "  nbrblocks : " << nbrblocks << endl;

      double *Chihihi_h = new double[nrows*ncols];
      double *Clohihi_h = new double[nrows*ncols];
      double *Chilohi_h = new double[nrows*ncols];
      double *Clolohi_h = new double[nrows*ncols];
      double *Chihilo_h = new double[nrows*ncols];
      double *Clohilo_h = new double[nrows*ncols];
      double *Chilolo_h = new double[nrows*ncols];
      double *Clololo_h = new double[nrows*ncols];
      const size_t szmat = nrows*ncols*sizeof(double);

      cudaMemcpy(Chihihi_h,Chihihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Clohihi_h,Clohihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Chilohi_h,Chilohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Clolohi_h,Clolohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Chihilo_h,Chihilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Clohilo_h,Clohilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Chilolo_h,Chilolo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Clololo_h,Clololo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the matrix C : " << endl;
      for(int i=rowoff; i<nrows; i++)
         for(int j=coloff; j<ncols; j++)
            cout << "C_h[" << i << "][" << j << "] : "
                 << Chihihi_h[j*nrows+i] << "  "
                 << Clohihi_h[j*nrows+i] << endl
                 << "            "
                 << Chilohi_h[j*nrows+i] << "  "
                 << Clolohi_h[j*nrows+i] << endl
                 << "            "
                 << Chihilo_h[j*nrows+i] << "  "
                 << Clohilo_h[j*nrows+i] << endl
                 << "            "
                 << Chilolo_h[j*nrows+i] << "  "
                 << Clololo_h[j*nrows+i] << endl;

      free(Chihihi_h); free(Clohihi_h); free(Chilohi_h); free(Clolohi_h);
      free(Chihilo_h); free(Clohilo_h); free(Chilolo_h); free(Clololo_h);
   }

   cudaEventRecord(start);
   dbl8_small_YWTC<<<nbrblocks,szt>>>
      (nrows,ncols,rowdim,coldim,szt,rowoff,coloff,
        YWThihihi_d, YWTlohihi_d, YWThilohi_d, YWTlolohi_d,
        YWThihilo_d, YWTlohilo_d, YWThilolo_d, YWTlololo_d,
          Chihihi_d,   Clohihi_d,   Chilohi_d,   Clolohi_d,
          Chihilo_d,   Clohilo_d,   Chilolo_d,   Clololo_d,
       YWTChihihi_d,YWTClohihi_d,YWTChilohi_d,YWTClolohi_d,
       YWTChihilo_d,YWTClohilo_d,YWTChilolo_d,YWTClololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_dbl_small_YWTC(rowdim,coldim,add,mul);

   if(verbose)
   {
      const size_t szmat = nrows*ncols*sizeof(double);

      cudaMemcpy(YWTChihihi_h,YWTChihihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTClohihi_h,YWTClohihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTChilohi_h,YWTChilohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTClolohi_h,YWTClolohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTChihilo_h,YWTChihilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTClohilo_h,YWTClohilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTChilolo_h,YWTChilolo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTClololo_h,YWTClololo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the YWTC matrix :" << endl;
      for(int i=rowoff; i<nrows; i++) 
         for(int j=coloff; j<ncols; j++)
            cout << "YWTC[" << i << "][" << j << "] : "
                 << YWTChihihi_h[j*nrows + i] << "  "
                 << YWTClohihi_h[j*nrows + i] << endl
                 << "             "
                 << YWTChilohi_h[j*nrows + i] << "  "
                 << YWTClolohi_h[j*nrows + i] << endl
                 << "             "
                 << YWTChihilo_h[j*nrows + i] << "  "
                 << YWTClohilo_h[j*nrows + i] << endl
                 << "             "
                 << YWTChilolo_h[j*nrows + i] << "  "
                 << YWTClololo_h[j*nrows + i] << endl;
   }
}

void GPU_cmplx8_small_YWHC
 ( int nrows, int ncols, int szt, int idx,
   double *YWTrehihihi_d, double *YWTrelohihi_d,
   double *YWTrehilohi_d, double *YWTrelolohi_d,
   double *YWTrehihilo_d, double *YWTrelohilo_d,
   double *YWTrehilolo_d, double *YWTrelololo_d,
   double *YWTimhihihi_d, double *YWTimlohihi_d,
   double *YWTimhilohi_d, double *YWTimlolohi_d,
   double *YWTimhihilo_d, double *YWTimlohilo_d,
   double *YWTimhilolo_d, double *YWTimlololo_d,
   double *Crehihihi_d, double *Crelohihi_d,
   double *Crehilohi_d, double *Crelolohi_d,
   double *Crehihilo_d, double *Crelohilo_d,
   double *Crehilolo_d, double *Crelololo_d,
   double *Cimhihihi_d, double *Cimlohihi_d,
   double *Cimhilohi_d, double *Cimlolohi_d,
   double *Cimhihilo_d, double *Cimlohilo_d,
   double *Cimhilolo_d, double *Cimlololo_d,
   double *YWTCrehihihi_d, double *YWTCrelohihi_d,
   double *YWTCrehilohi_d, double *YWTCrelolohi_d,
   double *YWTCrehihilo_d, double *YWTCrelohilo_d,
   double *YWTCrehilolo_d, double *YWTCrelololo_d,
   double *YWTCimhihihi_d, double *YWTCimlohihi_d,
   double *YWTCimhilohi_d, double *YWTCimlolohi_d,
   double *YWTCimhihilo_d, double *YWTCimlohilo_d,
   double *YWTCimhilolo_d, double *YWTCimlololo_d,
   double *YWTCrehihihi_h, double *YWTCrelohihi_h,
   double *YWTCrehilohi_h, double *YWTCrelolohi_h,
   double *YWTCrehihilo_h, double *YWTCrelohilo_h,
   double *YWTCrehilolo_h, double *YWTCrelololo_h,
   double *YWTCimhihihi_h, double *YWTCimlohihi_h,
   double *YWTCimhilohi_h, double *YWTCimlolohi_h,
   double *YWTCimhihilo_h, double *YWTCimlohilo_h,
   double *YWTCimhilolo_h, double *YWTCimlololo_h,
   double *lapms, long long int *add, long long int *mul, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int rowoff = idx*szt;
   const int rowdim = nrows - rowoff;
   const int coloff = (idx+1)*szt;
   const int coldim = ncols - coloff;
   const int nbrblocks = (int) ceil(rowdim*coldim/((double) szt));

   if(verbose)
   {
      cout << "in GPU_dbl_small_YWTC ..." << endl;
      cout << "-> nrows : " << nrows
           << "  ncols : " << ncols
           << "  szt : " << szt
           << "  idx : " << idx << endl;
      cout << "   rowdim : " << rowdim
           << "  coldim : " << coldim
           << "  rowoff : " << rowoff
           << "  coloff : " << coloff
           << "  nbrblocks : " << nbrblocks << endl;

      double *Crehihihi_h = new double[nrows*ncols];
      double *Crelohihi_h = new double[nrows*ncols];
      double *Crehilohi_h = new double[nrows*ncols];
      double *Crelolohi_h = new double[nrows*ncols];
      double *Crehihilo_h = new double[nrows*ncols];
      double *Crelohilo_h = new double[nrows*ncols];
      double *Crehilolo_h = new double[nrows*ncols];
      double *Crelololo_h = new double[nrows*ncols];
      double *Cimhihihi_h = new double[nrows*ncols];
      double *Cimlohihi_h = new double[nrows*ncols];
      double *Cimhilohi_h = new double[nrows*ncols];
      double *Cimlolohi_h = new double[nrows*ncols];
      double *Cimhihilo_h = new double[nrows*ncols];
      double *Cimlohilo_h = new double[nrows*ncols];
      double *Cimhilolo_h = new double[nrows*ncols];
      double *Cimlololo_h = new double[nrows*ncols];
      const size_t szmat = nrows*ncols*sizeof(double);

      cudaMemcpy(Crehihihi_h,Crehihihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Crelohihi_h,Crelohihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Crehilohi_h,Crehilohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Crelolohi_h,Crelolohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Crehihilo_h,Crehihilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Crelohilo_h,Crelohilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Crehilolo_h,Crehilolo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Crelololo_h,Crelololo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Cimhihihi_h,Cimhihihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Cimlohihi_h,Cimlohihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Cimhilohi_h,Cimhilohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Cimlolohi_h,Cimlolohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Cimhihilo_h,Cimhihilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Cimlohilo_h,Cimlohilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Cimhilolo_h,Cimhilolo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Cimlololo_h,Cimlololo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the matrix C : " << endl;
      for(int i=rowoff; i<nrows; i++)
         for(int j=coloff; j<ncols; j++)
         {
            cout << "C_h[" << i << "][" << j << "]re : "
                 << Crehihihi_h[j*nrows+i] << "  "
                 << Crelohihi_h[j*nrows+i] << endl
                 << "              "
                 << Crehilohi_h[j*nrows+i] << "  "
                 << Crelolohi_h[j*nrows+i] << endl
                 << "              "
                 << Crehihilo_h[j*nrows+i] << "  "
                 << Crelohilo_h[j*nrows+i] << endl
                 << "              "
                 << Crehilolo_h[j*nrows+i] << "  "
                 << Crelololo_h[j*nrows+i] << endl;
            cout << "C_h[" << i << "][" << j << "]im : "
                 << Cimhihihi_h[j*nrows+i] << "  "
                 << Cimlohihi_h[j*nrows+i] << endl
                 << "              "
                 << Cimhilohi_h[j*nrows+i] << "  "
                 << Cimlolohi_h[j*nrows+i] << endl
                 << "              "
                 << Cimhihilo_h[j*nrows+i] << "  "
                 << Cimlohilo_h[j*nrows+i] << endl
                 << "              "
                 << Cimhilolo_h[j*nrows+i] << "  "
                 << Cimlololo_h[j*nrows+i] << endl;
         }

      free(Crehihihi_h); free(Crelohihi_h);
      free(Crehilohi_h); free(Crelolohi_h);
      free(Crehihilo_h); free(Crelohilo_h);
      free(Crehilolo_h); free(Crelololo_h);
      free(Cimhihihi_h); free(Cimlohihi_h);
      free(Cimhilohi_h); free(Cimlolohi_h);
      free(Cimhihilo_h); free(Cimlohilo_h);
      free(Cimhilolo_h); free(Cimlololo_h);
   }
   cudaEventRecord(start);
   cmplx8_small_YWHC<<<nbrblocks,szt>>>
      (nrows,ncols,rowdim,coldim,szt,rowoff,coloff,
        YWTrehihihi_d, YWTrelohihi_d,  YWTrehilohi_d, YWTrelolohi_d,
        YWTrehihilo_d, YWTrelohilo_d,  YWTrehilolo_d, YWTrelololo_d,
        YWTimhihihi_d, YWTimlohihi_d,  YWTimhilohi_d, YWTimlolohi_d,
        YWTimhihilo_d, YWTimlohilo_d,  YWTimhilolo_d, YWTimlololo_d,
          Crehihihi_d,   Crelohihi_d,    Crehilohi_d,   Crelolohi_d,
          Crehihilo_d,   Crelohilo_d,    Crehilolo_d,   Crelololo_d,
          Cimhihihi_d,   Cimlohihi_d,    Cimhilohi_d,   Cimlolohi_d,
          Cimhihilo_d,   Cimlohilo_d,    Cimhilolo_d,   Cimlololo_d,
       YWTCrehihihi_d,YWTCrelohihi_d, YWTCrehilohi_d,YWTCrelolohi_d,
       YWTCrehihilo_d,YWTCrelohilo_d, YWTCrehilolo_d,YWTCrelololo_d,
       YWTCimhihihi_d,YWTCimlohihi_d, YWTCimhilohi_d,YWTCimlolohi_d,
       YWTCimhihilo_d,YWTCimlohilo_d, YWTCimhilolo_d,YWTCimlololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_cmplx_small_YWHC(rowdim,coldim,add,mul);

   if(verbose)
   {
      const size_t szmat = nrows*ncols*sizeof(double);

      cudaMemcpy(YWTCrehihihi_h,YWTCrehihihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTCrelohihi_h,YWTCrelohihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTCrehilohi_h,YWTCrehilohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTCrelolohi_h,YWTCrelolohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTCrehihilo_h,YWTCrehihilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTCrelohilo_h,YWTCrelohilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTCrehilolo_h,YWTCrehilolo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTCrelololo_h,YWTCrelololo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTCimhihihi_h,YWTCimhihihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTCimlohihi_h,YWTCimlohihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTCimhilohi_h,YWTCimhilohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTCimlolohi_h,YWTCimlolohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTCimhihilo_h,YWTCimhihilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTCimlohilo_h,YWTCimlohilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTCimhilolo_h,YWTCimhilolo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTCimlololo_h,YWTCimlololo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the YWTC matrix :" << endl;
      for(int i=rowoff; i<nrows; i++) 
         for(int j=coloff; j<ncols; j++)
         {
            cout << "YWTC[" << i << "][" << j << "]re : "
                 << YWTCrehihihi_h[j*nrows + i] << "  "
                 << YWTCrelohihi_h[j*nrows + i] << endl
                 << "               "
                 << YWTCrehilohi_h[j*nrows + i] << "  "
                 << YWTCrelolohi_h[j*nrows + i] << endl
                 << "               "
                 << YWTCrehihilo_h[j*nrows + i] << "  "
                 << YWTCrelohilo_h[j*nrows + i] << endl
                 << "               "
                 << YWTCrehilolo_h[j*nrows + i] << "  "
                 << YWTCrelololo_h[j*nrows + i] << endl;
            cout << "YWTC[" << i << "][" << j << "]im : "
                 << YWTCimhihihi_h[j*nrows + i] << "  "
                 << YWTCimlohihi_h[j*nrows + i] << endl
                 << "               "
                 << YWTCimhilohi_h[j*nrows + i] << "  "
                 << YWTCimlolohi_h[j*nrows + i] << endl
                 << "               "
                 << YWTCimhihilo_h[j*nrows + i] << "  "
                 << YWTCimlohilo_h[j*nrows + i] << endl
                 << "               "
                 << YWTCimhilolo_h[j*nrows + i] << "  "
                 << YWTCimlololo_h[j*nrows + i] << endl;
         }
   }
}

void GPU_dbl8_small_Qupdate
 ( int dim, int szt, int idx,
   double *Qhihihi_d, double *Qlohihi_d, double *Qhilohi_d, double *Qlolohi_d,
   double *Qhihilo_d, double *Qlohilo_d, double *Qhilolo_d, double *Qlololo_d,
   double *QWYThihihi_d, double *QWYTlohihi_d,
   double *QWYThilohi_d, double *QWYTlolohi_d,
   double *QWYThihilo_d, double *QWYTlohilo_d,
   double *QWYThilolo_d, double *QWYTlololo_d,
   double *Qhihihi_h, double *Qlohihi_h, double *Qhilohi_h, double *Qlolohi_h,
   double *Qhihilo_h, double *Qlohilo_h, double *Qhilolo_h, double *Qlololo_h,
   double *lapms, long long int *add, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int coloff = idx*szt;
   const int rowdim = dim - coloff;
   const int nbrblocks = (int) ceil(dim*rowdim/((double) szt));

   cudaEventRecord(start);
   dbl8_small_Qupdate<<<nbrblocks,szt>>>
      (dim,rowdim,szt,coloff,
       Qhihihi_d,Qlohihi_d,Qhilohi_d,Qlolohi_d,
       Qhihilo_d,Qlohilo_d,Qhilolo_d,Qlololo_d,
       QWYThihihi_d,QWYTlohihi_d,QWYThilohi_d,QWYTlolohi_d,
       QWYThihilo_d,QWYTlohilo_d,QWYThilolo_d,QWYTlololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_dbl_small_Qupdate(dim,rowdim,add);

   if(verbose)
   {
      const size_t szmat = dim*dim*sizeof(double);

      cudaMemcpy(Qhihihi_h,Qhihihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qlohihi_h,Qlohihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qhilohi_h,Qhilohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qlolohi_h,Qlolohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qhihilo_h,Qhihilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qlohilo_h,Qlohilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qhilolo_h,Qhilolo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qlololo_h,Qlololo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the Q matrix :" << endl;
      int ix = 0;
      for(int i=0; i<dim; i++) 
         for(int j=0; j<dim; j++) 
         {
            cout << "Q[" << i << "][" << j << "] : "
                 << Qhihihi_h[ix] << "  " << Qlohihi_h[ix] << endl
                 << "          "
                 << Qhilohi_h[ix] << "  " << Qlolohi_h[ix] << endl
                 << "          "
                 << Qhihilo_h[ix] << "  " << Qlohilo_h[ix] << endl
                 << "          "
                 << Qhilolo_h[ix] << "  " << Qlololo_h[ix] << endl;
            ix = ix + 1;
         }
   }
}

void GPU_cmplx8_small_Qupdate
 ( int dim, int szt, int idx,
   double *Qrehihihi_d, double *Qrelohihi_d,
   double *Qrehilohi_d, double *Qrelolohi_d,
   double *Qrehihilo_d, double *Qrelohilo_d,
   double *Qrehilolo_d, double *Qrelololo_d,
   double *Qimhihihi_d, double *Qimlohihi_d,
   double *Qimhilohi_d, double *Qimlolohi_d,
   double *Qimhihilo_d, double *Qimlohilo_d,
   double *Qimhilolo_d, double *Qimlololo_d,
   double *QWYTrehihihi_d, double *QWYTrelohihi_d,
   double *QWYTrehilohi_d, double *QWYTrelolohi_d,
   double *QWYTrehihilo_d, double *QWYTrelohilo_d,
   double *QWYTrehilolo_d, double *QWYTrelololo_d,
   double *QWYTimhihihi_d, double *QWYTimlohihi_d,
   double *QWYTimhilohi_d, double *QWYTimlolohi_d,
   double *QWYTimhihilo_d, double *QWYTimlohilo_d,
   double *QWYTimhilolo_d, double *QWYTimlololo_d,
   double *Qrehihihi_h, double *Qrelohihi_h,
   double *Qrehilohi_h, double *Qrelolohi_h,
   double *Qrehihilo_h, double *Qrelohilo_h,
   double *Qrehilolo_h, double *Qrelololo_h,
   double *Qimhihihi_h, double *Qimlohihi_h,
   double *Qimhilohi_h, double *Qimlolohi_h,
   double *Qimhihilo_h, double *Qimlohilo_h,
   double *Qimhilolo_h, double *Qimlololo_h,
   double *lapms, long long int *add, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int coloff = idx*szt;
   const int rowdim = dim - coloff;
   const int nbrblocks = (int) ceil(dim*rowdim/((double) szt));

   cudaEventRecord(start);
   cmplx8_small_Qupdate<<<nbrblocks,szt>>>
      (dim,rowdim,szt,coloff,
          Qrehihihi_d,   Qrelohihi_d,   Qrehilohi_d,   Qrelolohi_d,
          Qrehihilo_d,   Qrelohilo_d,   Qrehilolo_d,   Qrelololo_d,
          Qimhihihi_d,   Qimlohihi_d,   Qimhilohi_d,   Qimlolohi_d,
          Qimhihilo_d,   Qimlohilo_d,   Qimhilolo_d,   Qimlololo_d,
       QWYTrehihihi_d,QWYTrelohihi_d,QWYTrehilohi_d,QWYTrelolohi_d,
       QWYTrehihilo_d,QWYTrelohilo_d,QWYTrehilolo_d,QWYTrelololo_d,
       QWYTimhihihi_d,QWYTimlohihi_d,QWYTimhilohi_d,QWYTimlolohi_d,
       QWYTimhihilo_d,QWYTimlohilo_d,QWYTimhilolo_d,QWYTimlololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_cmplx_small_Qupdate(dim,rowdim,add);

   if(verbose)
   {
      const size_t szmat = dim*dim*sizeof(double);

      cudaMemcpy(Qrehihihi_h,Qrehihihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qrelohihi_h,Qrelohihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qrehilohi_h,Qrehilohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qrelolohi_h,Qrelolohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qrehihilo_h,Qrehihilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qrelohilo_h,Qrelohilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qrehilolo_h,Qrehilolo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qrelololo_h,Qrelololo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qimhihihi_h,Qimhihihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qimlohihi_h,Qimlohihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qimhilohi_h,Qimhilohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qimlolohi_h,Qimlolohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qimhihilo_h,Qimhihilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qimlohilo_h,Qimlohilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qimhilolo_h,Qimhilolo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qimlololo_h,Qimlololo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the Q matrix :" << endl;
      int ix = 0;
      for(int i=0; i<dim; i++) 
         for(int j=0; j<dim; j++) 
         {
            cout << "Q[" << i << "][" << j << "]re : "
                 << Qrehihihi_h[ix] << "  " << Qrelohihi_h[ix] << endl
                 << "            "
                 << Qrehilohi_h[ix] << "  " << Qrelolohi_h[ix] << endl
                 << "            "
                 << Qrehihilo_h[ix] << "  " << Qrelohilo_h[ix] << endl
                 << "            "
                 << Qrehilolo_h[ix] << "  " << Qrelololo_h[ix] << endl;
            cout << "Q[" << i << "][" << j << "]im : "
                 << Qimhihihi_h[ix] << "  " << Qimlohihi_h[ix] << endl
                 << "            "
                 << Qimhilohi_h[ix] << "  " << Qimlolohi_h[ix] << endl
                 << "            "
                 << Qimhihilo_h[ix] << "  " << Qimlohilo_h[ix] << endl
                 << "            "
                 << Qimhilolo_h[ix] << "  " << Qimlololo_h[ix] << endl;
            ix = ix + 1;
         }
   }
}

void GPU_dbl8_small_R_add_YWTC
 ( int nrows, int ncols, int szt, int idx,
   double *Rhihihi_d, double *Rlohihi_d, double *Rhilohi_d, double *Rlolohi_d,
   double *Rhihilo_d, double *Rlohilo_d, double *Rhilolo_d, double *Rlololo_d,
   double *YWTChihihi_d, double *YWTClohihi_d,
   double *YWTChilohi_d, double *YWTClolohi_d,
   double *YWTChihilo_d, double *YWTClohilo_d,
   double *YWTChilolo_d, double *YWTClololo_d,
   double *Rhihihi_h, double *Rlohihi_h, double *Rhilohi_h, double *Rlolohi_h,
   double *Rhihilo_h, double *Rlohilo_h, double *Rhilolo_h, double *Rlololo_h,
   double *lapms, long long int *add, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int rowoff = idx*szt;
   const int rowdim = nrows - rowoff;
   const int coloff = (idx+1)*szt;
   const int coldim = ncols - coloff;
   const int nbrblocks = (int) ceil(rowdim*coldim/((double) szt));

   cudaEventRecord(start);
   dbl8_small_R_add_YWTC<<<nbrblocks,szt>>>
      (nrows,coldim,szt,rowoff,coloff,
       Rhihihi_d,Rlohihi_d,Rhilohi_d,Rlolohi_d,
       Rhihilo_d,Rlohilo_d,Rhilolo_d,Rlololo_d,
       YWTChihihi_d,YWTClohihi_d,YWTChilohi_d,YWTClolohi_d,
       YWTChihilo_d,YWTClohilo_d,YWTChilolo_d,YWTClololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_dbl_small_R_add_YWTC(nrows,coldim,szt,rowoff,coloff,add);

   if(verbose)
   {
      const size_t szmat = nrows*ncols*sizeof(double);

      cudaMemcpy(Rhihihi_h,Rhihihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Rlohihi_h,Rlohihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Rhilohi_h,Rhilohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Rlolohi_h,Rlolohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Rhihilo_h,Rhihilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Rlohilo_h,Rlohilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Rhilolo_h,Rhilolo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Rlololo_h,Rlololo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the R matrix :" << endl;
      for(int i=rowoff; i<nrows; i++) 
         for(int j=coloff; j<ncols; j++)
            cout << "R[" << i << "][" << j << "] : "
                 << Rhihihi_h[j*nrows + i] << "  "
                 << Rlohihi_h[j*nrows + i] << endl
                 << "          "
                 << Rhilohi_h[j*nrows + i] << "  "
                 << Rlolohi_h[j*nrows + i] << endl
                 << "          "
                 << Rhihilo_h[j*nrows + i] << "  "
                 << Rlohilo_h[j*nrows + i] << endl
                 << "          "
                 << Rhilolo_h[j*nrows + i] << "  "
                 << Rlololo_h[j*nrows + i] << endl;
   }
}

void GPU_cmplx8_small_R_add_YWHC
 ( int nrows, int ncols, int szt, int idx,
   double *Rrehihihi_d, double *Rrelohihi_d,
   double *Rrehilohi_d, double *Rrelolohi_d,
   double *Rrehihilo_d, double *Rrelohilo_d,
   double *Rrehilolo_d, double *Rrelololo_d,
   double *Rimhihihi_d, double *Rimlohihi_d,
   double *Rimhilohi_d, double *Rimlolohi_d,
   double *Rimhihilo_d, double *Rimlohilo_d,
   double *Rimhilolo_d, double *Rimlololo_d,
   double *YWTCrehihihi_d, double *YWTCrelohihi_d,
   double *YWTCrehilohi_d, double *YWTCrelolohi_d,
   double *YWTCrehihilo_d, double *YWTCrelohilo_d,
   double *YWTCrehilolo_d, double *YWTCrelololo_d,
   double *YWTCimhihihi_d, double *YWTCimlohihi_d,
   double *YWTCimhilohi_d, double *YWTCimlolohi_d,
   double *YWTCimhihilo_d, double *YWTCimlohilo_d,
   double *YWTCimhilolo_d, double *YWTCimlololo_d,
   double *Rrehihihi_h, double *Rrelohihi_h,
   double *Rrehilohi_h, double *Rrelolohi_h,
   double *Rrehihilo_h, double *Rrelohilo_h,
   double *Rrehilolo_h, double *Rrelololo_h,
   double *Rimhihihi_h, double *Rimlohihi_h,
   double *Rimhilohi_h, double *Rimlolohi_h,
   double *Rimhihilo_h, double *Rimlohilo_h,
   double *Rimhilolo_h, double *Rimlololo_h,
   double *lapms, long long int *add, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int rowoff = idx*szt;
   const int rowdim = nrows - rowoff;
   const int coloff = (idx+1)*szt;
   const int coldim = ncols - coloff;
   const int nbrblocks = (int) ceil(rowdim*coldim/((double) szt));

   cudaEventRecord(start);
   cmplx8_small_R_add_YWHC<<<nbrblocks,szt>>>
      (nrows,coldim,szt,rowoff,coloff,
          Rrehihihi_d,   Rrelohihi_d,   Rrehilohi_d,   Rrelolohi_d,
          Rrehihilo_d,   Rrelohilo_d,   Rrehilolo_d,   Rrelololo_d,
          Rimhihihi_d,   Rimlohihi_d,   Rimhilohi_d,   Rimlolohi_d,
          Rimhihilo_d,   Rimlohilo_d,   Rimhilolo_d,   Rimlololo_d,
       YWTCrehihihi_d,YWTCrelohihi_d,YWTCrehilohi_d,YWTCrelolohi_d,
       YWTCrehihilo_d,YWTCrelohilo_d,YWTCrehilolo_d,YWTCrelololo_d,
       YWTCimhihihi_d,YWTCimlohihi_d,YWTCimhilohi_d,YWTCimlolohi_d,
       YWTCimhihilo_d,YWTCimlohilo_d,YWTCimhilolo_d,YWTCimlololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_cmplx_small_R_add_YWHC(nrows,coldim,szt,rowoff,coloff,add);

   if(verbose)
   {
      const size_t szmat = nrows*ncols*sizeof(double);

      cudaMemcpy(Rrehihihi_h,Rrehihihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Rrelohihi_h,Rrelohihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Rrehilohi_h,Rrehilohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Rrelolohi_h,Rrelolohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Rrehihilo_h,Rrehihilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Rrelohilo_h,Rrelohilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Rrehilolo_h,Rrehilolo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Rrelololo_h,Rrelololo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Rimhihihi_h,Rimhihihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Rimlohihi_h,Rimlohihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Rimhilohi_h,Rimhilohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Rimlolohi_h,Rimlolohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Rimhihilo_h,Rimhihilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Rimlohilo_h,Rimlohilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Rimhilolo_h,Rimhilolo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Rimlololo_h,Rimlololo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the R matrix :" << endl;
      for(int i=rowoff; i<nrows; i++) 
         for(int j=coloff; j<ncols; j++)
         {
            cout << "R[" << i << "][" << j << "]re : "
                 << Rrehihihi_h[j*nrows + i] << "  "
                 << Rrelohihi_h[j*nrows + i] << endl
                 << "            "
                 << Rrehilohi_h[j*nrows + i] << "  "
                 << Rrelolohi_h[j*nrows + i] << endl
                 << "            "
                 << Rrehihilo_h[j*nrows + i] << "  "
                 << Rrelohilo_h[j*nrows + i] << endl
                 << "            "
                 << Rrehilolo_h[j*nrows + i] << "  "
                 << Rrelololo_h[j*nrows + i] << endl;
            cout << "R[" << i << "][" << j << "]im : "
                 << Rimhihihi_h[j*nrows + i] << "  "
                 << Rimlohihi_h[j*nrows + i] << endl
                 << "            "
                 << Rimhilohi_h[j*nrows + i] << "  "
                 << Rimlolohi_h[j*nrows + i] << endl
                 << "            "
                 << Rimhihilo_h[j*nrows + i] << "  "
                 << Rimlohilo_h[j*nrows + i] << endl
                 << "            "
                 << Rimhilolo_h[j*nrows + i] << "  "
                 << Rimlololo_h[j*nrows + i] << endl;
         }
   }
}

void GPU_dbl8_blocked_houseqr
 ( int nrows, int ncols, int szt, int nbt,
   double **Ahihihi, double **Alohihi, double **Ahilohi, double **Alolohi,
   double **Ahihilo, double **Alohilo, double **Ahilolo, double **Alololo,
   double **Qhihihi, double **Qlohihi, double **Qhilohi, double **Qlolohi,
   double **Qhihilo, double **Qlohilo, double **Qhilolo, double **Qlololo,
   double **Rhihihi, double **Rlohihi, double **Rhilohi, double **Rlolohi,
   double **Rhihilo, double **Rlohilo, double **Rhilolo, double **Rlololo,
   double *houselapms, double *RTvlapms, double *tileRlapms,
   double *vb2Wlapms, double *WYTlapms, double *QWYTlapms, double *Qaddlapms,
   double *YWTlapms, double *YWTClapms, double *Raddlapms,
   double *walltimesec, long long int *addcnt, long long int *mulcnt,
   long long int *divcnt, long long int *sqrtcnt, bool verbose )
{
   const int dim = nrows*ncols;          // total number of doubles
   const int nrows2 = nrows*nrows;
   double *Ahihihi_h = new double[dim];    // A on the host
   double *Alohihi_h = new double[dim]; 
   double *Ahilohi_h = new double[dim];
   double *Alolohi_h = new double[dim]; 
   double *Ahihilo_h = new double[dim]; 
   double *Alohilo_h = new double[dim]; 
   double *Ahilolo_h = new double[dim];
   double *Alololo_h = new double[dim]; 
   double *Ahihihi_d;                      // A on the device
   double *Alohihi_d; 
   double *Ahilohi_d; 
   double *Alolohi_d; 
   double *Ahihilo_d;
   double *Alohilo_d; 
   double *Ahilolo_d; 
   double *Alololo_d; 
   double *Qhihihi_h = new double[nrows2]; // Q on the host
   double *Qlohihi_h = new double[nrows2]; 
   double *Qhilohi_h = new double[nrows2]; 
   double *Qlolohi_h = new double[nrows2]; 
   double *Qhihilo_h = new double[nrows2];
   double *Qlohilo_h = new double[nrows2]; 
   double *Qhilolo_h = new double[nrows2]; 
   double *Qlololo_h = new double[nrows2]; 
   double *Qhihihi_d;                      // Q on the device
   double *Qlohihi_d;
   double *Qhilohi_d;
   double *Qlolohi_d;
   double *Qhihilo_d;                      // Q on the device
   double *Qlohilo_d;
   double *Qhilolo_d;
   double *Qlololo_d;
   double *vhihihi_h = new double[nrows];  // Householder vector
   double *vlohihi_h = new double[nrows];
   double *vhilohi_h = new double[nrows];
   double *vlolohi_h = new double[nrows];
   double *vhihilo_h = new double[nrows]; 
   double *vlohilo_h = new double[nrows];
   double *vhilolo_h = new double[nrows];
   double *vlololo_h = new double[nrows];
   double *betahihihi_h = new double[szt]; //  beta on the host
   double *betalohihi_h = new double[szt]; 
   double *betahilohi_h = new double[szt]; 
   double *betalolohi_h = new double[szt]; 
   double *betahihilo_h = new double[szt]; 
   double *betalohilo_h = new double[szt]; 
   double *betahilolo_h = new double[szt]; 
   double *betalololo_h = new double[szt]; 
   double *betahihihi_d;                   // beta on the device
   double *betalohihi_d;
   double *betahilohi_d;
   double *betalolohi_d;
   double *betahihilo_d;
   double *betalohilo_d;
   double *betahilolo_d;
   double *betalololo_d;
   double *Vhihihi_h = new double[nrows*szt]; // V matrix
   double *Vlohihi_h = new double[nrows*szt];
   double *Vhilohi_h = new double[nrows*szt];
   double *Vlolohi_h = new double[nrows*szt];
   double *Vhihilo_h = new double[nrows*szt];
   double *Vlohilo_h = new double[nrows*szt];
   double *Vhilolo_h = new double[nrows*szt];
   double *Vlololo_h = new double[nrows*szt];
   double *Vhihihi_d;                         // V on the device
   double *Vlohihi_d;
   double *Vhilohi_d;
   double *Vlolohi_d;
   double *Vhihilo_d;
   double *Vlohilo_d;
   double *Vhilolo_d;
   double *Vlololo_d;
   double *Whihihi_h = new double[nrows*szt]; // W on the host
   double *Wlohihi_h = new double[nrows*szt];
   double *Whilohi_h = new double[nrows*szt];
   double *Wlolohi_h = new double[nrows*szt];
   double *Whihilo_h = new double[nrows*szt];
   double *Wlohilo_h = new double[nrows*szt];
   double *Whilolo_h = new double[nrows*szt];
   double *Wlololo_h = new double[nrows*szt];
   double *Whihihi_d;                         // W on the device
   double *Wlohihi_d;
   double *Whilohi_d;
   double *Wlolohi_d;
   double *Whihilo_d;
   double *Wlohilo_d;
   double *Whilolo_d;
   double *Wlololo_d;
   double *WYThihihi_h = new double[nrows2];  // W*Y^T 
   double *WYTlohihi_h = new double[nrows2];
   double *WYThilohi_h = new double[nrows2];
   double *WYTlolohi_h = new double[nrows2];
   double *WYThihilo_h = new double[nrows2];
   double *WYTlohilo_h = new double[nrows2];
   double *WYThilolo_h = new double[nrows2];
   double *WYTlololo_h = new double[nrows2];
   double *WYThihihi_d;                       // WYT on the device
   double *WYTlohihi_d;
   double *WYThilohi_d;
   double *WYTlolohi_d;
   double *WYThihilo_d;
   double *WYTlohilo_d;
   double *WYThilolo_d;
   double *WYTlololo_d;
   double *YWThihihi_h = new double[nrows2];  // Y*W^T
   double *YWTlohihi_h = new double[nrows2];
   double *YWThilohi_h = new double[nrows2];
   double *YWTlolohi_h = new double[nrows2];
   double *YWThihilo_h = new double[nrows2];
   double *YWTlohilo_h = new double[nrows2];
   double *YWThilolo_h = new double[nrows2];
   double *YWTlololo_h = new double[nrows2];
   double *YWThihihi_d;                       // YWT on the device
   double *YWTlohihi_d;
   double *YWThilohi_d;
   double *YWTlolohi_d;
   double *YWThihilo_d;
   double *YWTlohilo_d;
   double *YWThilolo_d;
   double *YWTlololo_d;
   double *QWYThihihi_h = new double[nrows2]; // Q*WY^T
   double *QWYTlohihi_h = new double[nrows2];
   double *QWYThilohi_h = new double[nrows2];
   double *QWYTlolohi_h = new double[nrows2];
   double *QWYThihilo_h = new double[nrows2];
   double *QWYTlohilo_h = new double[nrows2];
   double *QWYThilolo_h = new double[nrows2];
   double *QWYTlololo_h = new double[nrows2];
   double *QWYThihihi_d;                      // QWYT on the device
   double *QWYTlohihi_d;
   double *QWYThilohi_d;
   double *QWYTlolohi_d;
   double *QWYThihilo_d;
   double *QWYTlohilo_d;
   double *QWYThilolo_d;
   double *QWYTlololo_d;
   double *YWTChihihi_h = new double[dim];    // YWT*C on the host
   double *YWTClohihi_h = new double[dim];
   double *YWTChilohi_h = new double[dim];
   double *YWTClolohi_h = new double[dim];
   double *YWTChihilo_h = new double[dim];
   double *YWTClohilo_h = new double[dim];
   double *YWTChilolo_h = new double[dim];
   double *YWTClololo_h = new double[dim];
   double *YWTChihihi_d;                      // YWTC on the device
   double *YWTClohihi_d;
   double *YWTChilohi_d;
   double *YWTClolohi_d;
   double *YWTChihilo_d;
   double *YWTClohilo_d;
   double *YWTChilolo_d;
   double *YWTClololo_d;
   double *RTdotvhihihi_h = new double[nrows2]; // R^T dotted with v
   double *RTdotvlohihi_h = new double[nrows2];
   double *RTdotvhilohi_h = new double[nrows2];
   double *RTdotvlolohi_h = new double[nrows2];
   double *RTdotvhihilo_h = new double[nrows2];
   double *RTdotvlohilo_h = new double[nrows2];
   double *RTdotvhilolo_h = new double[nrows2];
   double *RTdotvlololo_h = new double[nrows2];
   double *RTdotvhihihi_d;                      // RTdotv on the device
   double *RTdotvlohihi_d;
   double *RTdotvhilohi_d;
   double *RTdotvlolohi_d;
   double *RTdotvhihilo_d;                      // RTdotv on the device
   double *RTdotvlohilo_d;
   double *RTdotvhilolo_d;
   double *RTdotvlololo_d;
   double *bRTvhihihi_h = new double[nrows];  // beta*R^T*v
   double *bRTvlohihi_h = new double[nrows];
   double *bRTvhilohi_h = new double[nrows];
   double *bRTvlolohi_h = new double[nrows];
   double *bRTvhihilo_h = new double[nrows];
   double *bRTvlohilo_h = new double[nrows];
   double *bRTvhilolo_h = new double[nrows];
   double *bRTvlololo_h = new double[nrows];
   double *bRTvhihihi_d;                      // bRTv on the device
   double *bRTvlohihi_d;
   double *bRTvhilohi_d;
   double *bRTvlolohi_d;
   double *bRTvhihilo_d; 
   double *bRTvlohilo_d;
   double *bRTvhilolo_d;
   double *bRTvlololo_d;
   double *sumshihihi_h = new double[nrows];  // subsums for large house
   double *sumslohihi_h = new double[nrows];
   double *sumshilohi_h = new double[nrows];
   double *sumslolohi_h = new double[nrows];
   double *sumshihilo_h = new double[nrows];
   double *sumslohilo_h = new double[nrows];
   double *sumshilolo_h = new double[nrows];
   double *sumslololo_h = new double[nrows];
   double *sumshihihi_d;                      // sums on the device
   double *sumslohihi_d;
   double *sumshilohi_d;
   double *sumslolohi_d;
   double *sumshihilo_d;
   double *sumslohilo_d;
   double *sumshilolo_d;
   double *sumslololo_d;
   double sigmahihihi_h,sigmalohihi_h,sigmahilohi_h,sigmalolohi_h;
   double sigmahihilo_h,sigmalohilo_h,sigmahilolo_h,sigmalololo_h;
   double *sigmahihihi_d;                     // sigma on the device
   double *sigmalohihi_d;
   double *sigmahilohi_d;
   double *sigmalolohi_d;
   double *sigmahihilo_d;
   double *sigmalohilo_d;
   double *sigmahilolo_d;
   double *sigmalololo_d;

   int ix = 0;                          // copy the columns of A to A_h
   for(int j=0; j<ncols; j++)   
      for(int i=0; i<nrows; i++)
      {
         Ahihihi_h[ix]   = Ahihihi[i][j];
         Alohihi_h[ix]   = Alohihi[i][j];
         Ahilohi_h[ix]   = Ahilohi[i][j];
         Alolohi_h[ix]   = Alolohi[i][j];
         Ahihilo_h[ix]   = Ahihilo[i][j];
         Alohilo_h[ix]   = Alohilo[i][j];
         Ahilolo_h[ix]   = Ahilolo[i][j];
         Alololo_h[ix++] = Alololo[i][j];
      }

   ix = 0;                              // initialize Q with identity
   for(int i=0; i<nrows; i++)
   {
      for(int j=0; j<nrows; j++)
      {
         if(i == j)
         {
            Qhihihi_h[ix]   = 1.0;
            Qlohihi_h[ix]   = 0.0;
            Qhilohi_h[ix]   = 0.0;
            Qlolohi_h[ix]   = 0.0;
            Qhihilo_h[ix]   = 0.0;
            Qlohilo_h[ix]   = 0.0;
            Qhilolo_h[ix]   = 0.0;
            Qlololo_h[ix++] = 0.0;
         }
         else
         {
            Qhihihi_h[ix]   = 0.0;
            Qlohihi_h[ix]   = 0.0;
            Qhilohi_h[ix]   = 0.0;
            Qlolohi_h[ix]   = 0.0;
            Qhihilo_h[ix]   = 0.0;
            Qlohilo_h[ix]   = 0.0;
            Qhilolo_h[ix]   = 0.0;
            Qlololo_h[ix++] = 0.0;
         }
      }
   }
   const size_t sznum = dim*sizeof(double);
   cudaMalloc((void**)&Ahihihi_d,sznum);
   cudaMalloc((void**)&Alohihi_d,sznum);
   cudaMalloc((void**)&Ahilohi_d,sznum);
   cudaMalloc((void**)&Alolohi_d,sznum);
   cudaMalloc((void**)&Ahihilo_d,sznum);
   cudaMalloc((void**)&Alohilo_d,sznum);
   cudaMalloc((void**)&Ahilolo_d,sznum);
   cudaMalloc((void**)&Alololo_d,sznum);
   cudaMemcpy(Ahihihi_d,Ahihihi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Alohihi_d,Alohihi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Ahilohi_d,Ahilohi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Alolohi_d,Alolohi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Ahihilo_d,Ahihilo_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Alohilo_d,Alohilo_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Ahilolo_d,Ahilolo_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Alololo_d,Alololo_h,sznum,cudaMemcpyHostToDevice);

   const size_t szbeta = szt*sizeof(double);
   cudaMalloc((void**)&betahihihi_d,szbeta);
   cudaMalloc((void**)&betalohihi_d,szbeta);
   cudaMalloc((void**)&betahilohi_d,szbeta);
   cudaMalloc((void**)&betalolohi_d,szbeta);
   cudaMalloc((void**)&betahihilo_d,szbeta);
   cudaMalloc((void**)&betalohilo_d,szbeta);
   cudaMalloc((void**)&betahilolo_d,szbeta);
   cudaMalloc((void**)&betalololo_d,szbeta);

   for(int i=0; i<szt; i++)
   {
      betahihihi_h[i] = 0.0;
      betalohihi_h[i] = 0.0;
      betahilohi_h[i] = 0.0;
      betalolohi_h[i] = 0.0;
      betahihilo_h[i] = 0.0;
      betalohilo_h[i] = 0.0;
      betahilolo_h[i] = 0.0;
      betalololo_h[i] = 0.0;
   }
   cudaMemcpy(betahihihi_d,betahihihi_h,szbeta,cudaMemcpyHostToDevice);
   cudaMemcpy(betalohihi_d,betalohihi_h,szbeta,cudaMemcpyHostToDevice);
   cudaMemcpy(betahilohi_d,betahilohi_h,szbeta,cudaMemcpyHostToDevice);
   cudaMemcpy(betalolohi_d,betalolohi_h,szbeta,cudaMemcpyHostToDevice);
   cudaMemcpy(betahihilo_d,betahihilo_h,szbeta,cudaMemcpyHostToDevice);
   cudaMemcpy(betalohilo_d,betalohilo_h,szbeta,cudaMemcpyHostToDevice);
   cudaMemcpy(betahilolo_d,betahilolo_h,szbeta,cudaMemcpyHostToDevice);
   cudaMemcpy(betalololo_d,betalololo_h,szbeta,cudaMemcpyHostToDevice);

   const size_t szhouse = nrows*sizeof(double);
   const size_t szpad = szt*sizeof(double);  // padding for nonsquare tiles
   const size_t szVandW = szt*szhouse;
   cudaMalloc((void**)&Vhihihi_d,szVandW + szpad); // pad only in allocation
   cudaMalloc((void**)&Vlohihi_d,szVandW + szpad);
   cudaMalloc((void**)&Vhilohi_d,szVandW + szpad);
   cudaMalloc((void**)&Vlolohi_d,szVandW + szpad);
   cudaMalloc((void**)&Vhihilo_d,szVandW + szpad);
   cudaMalloc((void**)&Vlohilo_d,szVandW + szpad);
   cudaMalloc((void**)&Vhilolo_d,szVandW + szpad);
   cudaMalloc((void**)&Vlololo_d,szVandW + szpad);

   ix = 0;
   for(int i=0; i<nrows*szt; i++)
   {
      Vhihihi_h[ix] = 0.0;
      Vlohihi_h[ix] = 0.0; 
      Vhilohi_h[ix] = 0.0; 
      Vlolohi_h[ix] = 0.0; 
      Vhihilo_h[ix] = 0.0;
      Vlohilo_h[ix] = 0.0; 
      Vhilolo_h[ix] = 0.0; 
      Vlololo_h[ix++] = 0.0; 
   }
   Vhihihi_h[--ix] = 1.0; // initialize last vector for square tiles

   cudaMemcpy(Vhihihi_d,Vhihihi_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vlohihi_d,Vlohihi_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vhilohi_d,Vhilohi_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vlolohi_d,Vlolohi_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vhihilo_d,Vhihilo_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vlohilo_d,Vlohilo_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vhilolo_d,Vhilolo_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vlololo_d,Vlololo_h,szVandW,cudaMemcpyHostToDevice);
   cudaMalloc((void**)&Whihihi_d,szVandW + szpad); // pad only in allocation
   cudaMalloc((void**)&Wlohihi_d,szVandW + szpad); 
   cudaMalloc((void**)&Whilohi_d,szVandW + szpad); 
   cudaMalloc((void**)&Wlolohi_d,szVandW + szpad); 
   cudaMalloc((void**)&Whihilo_d,szVandW + szpad); 
   cudaMalloc((void**)&Wlohilo_d,szVandW + szpad); 
   cudaMalloc((void**)&Whilolo_d,szVandW + szpad); 
   cudaMalloc((void**)&Wlololo_d,szVandW + szpad); 

   cudaMalloc((void**)&RTdotvhihihi_d,szVandW + szpad);
   cudaMalloc((void**)&RTdotvlohihi_d,szVandW + szpad);
   cudaMalloc((void**)&RTdotvhilohi_d,szVandW + szpad);
   cudaMalloc((void**)&RTdotvlolohi_d,szVandW + szpad);
   cudaMalloc((void**)&RTdotvhihilo_d,szVandW + szpad);
   cudaMalloc((void**)&RTdotvlohilo_d,szVandW + szpad);
   cudaMalloc((void**)&RTdotvhilolo_d,szVandW + szpad);
   cudaMalloc((void**)&RTdotvlololo_d,szVandW + szpad);
   cudaMalloc((void**)&bRTvhihihi_d,szhouse + szpad);
   cudaMalloc((void**)&bRTvlohihi_d,szhouse + szpad);
   cudaMalloc((void**)&bRTvhilohi_d,szhouse + szpad);
   cudaMalloc((void**)&bRTvlolohi_d,szhouse + szpad);
   cudaMalloc((void**)&bRTvhihilo_d,szhouse + szpad);
   cudaMalloc((void**)&bRTvlohilo_d,szhouse + szpad);
   cudaMalloc((void**)&bRTvhilolo_d,szhouse + szpad);
   cudaMalloc((void**)&bRTvlololo_d,szhouse + szpad);

   cudaMalloc((void**)&sumshihihi_d,szhouse);
   cudaMalloc((void**)&sumslohihi_d,szhouse);
   cudaMalloc((void**)&sumshilohi_d,szhouse);
   cudaMalloc((void**)&sumslolohi_d,szhouse);
   cudaMalloc((void**)&sumshihilo_d,szhouse);
   cudaMalloc((void**)&sumslohilo_d,szhouse);
   cudaMalloc((void**)&sumshilolo_d,szhouse);
   cudaMalloc((void**)&sumslololo_d,szhouse);
   cudaMalloc((void**)&sigmahihihi_d,sizeof(double));
   cudaMalloc((void**)&sigmalohihi_d,sizeof(double));
   cudaMalloc((void**)&sigmahilohi_d,sizeof(double));
   cudaMalloc((void**)&sigmalolohi_d,sizeof(double));
   cudaMalloc((void**)&sigmahihilo_d,sizeof(double));
   cudaMalloc((void**)&sigmalohilo_d,sizeof(double));
   cudaMalloc((void**)&sigmahilolo_d,sizeof(double));
   cudaMalloc((void**)&sigmalololo_d,sizeof(double));

   const size_t szWYT = nrows2*sizeof(double);
   cudaMalloc((void**)&WYThihihi_d,szWYT + szpad); // pad for W*Y^T product
   cudaMalloc((void**)&WYTlohihi_d,szWYT + szpad); 
   cudaMalloc((void**)&WYThilohi_d,szWYT + szpad); 
   cudaMalloc((void**)&WYTlolohi_d,szWYT + szpad); 
   cudaMalloc((void**)&WYThihilo_d,szWYT + szpad);
   cudaMalloc((void**)&WYTlohilo_d,szWYT + szpad); 
   cudaMalloc((void**)&WYThilolo_d,szWYT + szpad); 
   cudaMalloc((void**)&WYTlololo_d,szWYT + szpad); 
   cudaMalloc((void**)&Qhihihi_d,szWYT + szpad);
   cudaMalloc((void**)&Qlohihi_d,szWYT + szpad);
   cudaMalloc((void**)&Qhilohi_d,szWYT + szpad);
   cudaMalloc((void**)&Qlolohi_d,szWYT + szpad);
   cudaMalloc((void**)&Qhihilo_d,szWYT + szpad);
   cudaMalloc((void**)&Qlohilo_d,szWYT + szpad);
   cudaMalloc((void**)&Qhilolo_d,szWYT + szpad);
   cudaMalloc((void**)&Qlololo_d,szWYT + szpad);
   cudaMemcpy(Qhihihi_d,Qhihihi_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qlohihi_d,Qlohihi_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qhilohi_d,Qhilohi_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qlolohi_d,Qlolohi_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qhihilo_d,Qhihilo_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qlohilo_d,Qlohilo_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qhilolo_d,Qhilolo_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qlololo_d,Qlololo_h,szWYT,cudaMemcpyHostToDevice);
   cudaMalloc((void**)&QWYThihihi_d,szWYT + szpad);
   cudaMalloc((void**)&QWYTlohihi_d,szWYT + szpad);
   cudaMalloc((void**)&QWYThilohi_d,szWYT + szpad);
   cudaMalloc((void**)&QWYTlolohi_d,szWYT + szpad);
   cudaMalloc((void**)&QWYThihilo_d,szWYT + szpad);
   cudaMalloc((void**)&QWYTlohilo_d,szWYT + szpad);
   cudaMalloc((void**)&QWYThilolo_d,szWYT + szpad);
   cudaMalloc((void**)&QWYTlololo_d,szWYT + szpad);

   const size_t szYWT = nrows2*sizeof(double);
   cudaMalloc((void**)&YWThihihi_d,szYWT + szpad); // pad for Y*W^T product
   cudaMalloc((void**)&YWTlohihi_d,szYWT + szpad);
   cudaMalloc((void**)&YWThilohi_d,szYWT + szpad);
   cudaMalloc((void**)&YWTlolohi_d,szYWT + szpad);
   cudaMalloc((void**)&YWThihilo_d,szYWT + szpad);
   cudaMalloc((void**)&YWTlohilo_d,szYWT + szpad);
   cudaMalloc((void**)&YWThilolo_d,szYWT + szpad);
   cudaMalloc((void**)&YWTlololo_d,szYWT + szpad);
   cudaMalloc((void**)&YWTChihihi_d,sznum + szpad);
   cudaMalloc((void**)&YWTClohihi_d,sznum + szpad);
   cudaMalloc((void**)&YWTChilohi_d,sznum + szpad);
   cudaMalloc((void**)&YWTClolohi_d,sznum + szpad);
   cudaMalloc((void**)&YWTChihilo_d,sznum + szpad);
   cudaMalloc((void**)&YWTClohilo_d,sznum + szpad);
   cudaMalloc((void**)&YWTChilolo_d,sznum + szpad);
   cudaMalloc((void**)&YWTClololo_d,sznum + szpad);

   *houselapms = 0.0; *RTvlapms = 0.0; *tileRlapms = 0.0; *vb2Wlapms = 0.0;
   *WYTlapms = 0.0; *QWYTlapms = 0.0; *Qaddlapms = 0.0;
   *YWTlapms = 0.0; *YWTClapms = 0.0; *Raddlapms = 0.0;
   *addcnt = 0; *mulcnt = 0; *divcnt = 0; *sqrtcnt = 0;
   struct timeval begintime,endtime; // wall clock time of computations

   gettimeofday(&begintime,0);

   for(int k=0; k<nbt; k++)       // k runs over the number of blocks
   {
      if(verbose)
         cout << "Tile k = " << k << " out of " << nbt << " ..." << endl;

      int colidx,nrows1;

      for(int L=0; L<szt; L++)  // L runs over the columns in one block
      {
         colidx = k*szt + L;              // index of the current column
         nrows1 = nrows - colidx - 1;     // #rows in Householder vector - 1

         if(verbose)
            cout << "-> current column : " << colidx << endl
                 << "-> #nrows in Householder vector - 1 : "
                 << nrows1 << endl;

         if(nrows1 <= szt)
         {
            GPU_dbl8_small_house
               (nrows,ncols,szt,nbt,colidx,nrows1,k,L,
                   Ahihihi_h,   Alohihi_h,   Ahilohi_h,   Alolohi_h,
                   Ahihilo_h,   Alohilo_h,   Ahilolo_h,   Alololo_h,
                   Ahihihi_d,   Alohihi_d,   Ahilohi_d,   Alolohi_d,
                   Ahihilo_d,   Alohilo_d,   Ahilolo_d,   Alololo_d,
                   vhihihi_h,   vlohihi_h,   vhilohi_h,   vlolohi_h,
                   vhihilo_h,   vlohilo_h,   vhilolo_h,   vlololo_h,
                   Vhihihi_d,   Vlohihi_d,   Vhilohi_d,   Vlolohi_d,
                   Vhihilo_d,   Vlohilo_d,   Vhilolo_d,   Vlololo_d,
                betahihihi_h,betalohihi_h,betahilohi_h,betalolohi_h,
                betahihilo_h,betalohilo_h,betahilolo_h,betalololo_h,
                betahihihi_d,betalohihi_d,betahilohi_d,betalolohi_d,
                betahihilo_d,betalohilo_d,betahilolo_d,betalololo_d,
                houselapms,addcnt,mulcnt,divcnt,sqrtcnt,verbose);

            if((betahihihi_h[L] == 0.0) && (betalohihi_h[L] == 0.0) &&
               (betahilohi_h[L] == 0.0) && (betalolohi_h[L] == 0.0) &&
               (betahihilo_h[L] == 0.0) && (betalohilo_h[L] == 0.0) &&
               (betahilolo_h[L] == 0.0) && (betalololo_h[L] == 0.0))
            {
               if(verbose) cout << "Zero beta detected." << endl;
            }
            else
            {
               GPU_dbl8_small_leftRupdate
                  (nrows,ncols,szt,colidx,k,L,
                      Ahihihi_h,   Alohihi_h,   Ahilohi_h,   Alolohi_h,
                      Ahihilo_h,   Alohilo_h,   Ahilolo_h,   Alololo_h,
                      Ahihihi_d,   Alohihi_d,   Ahilohi_d,   Alolohi_d,
                      Ahihilo_d,   Alohilo_d,   Ahilolo_d,   Alololo_d,
                      Vhihihi_d,   Vlohihi_d,   Vhilohi_d,   Vlolohi_d,
                      Vhihilo_d,   Vlohilo_d,   Vhilolo_d,   Vlololo_d,
                   betahihihi_h,betalohihi_h,betahilohi_h,betalolohi_h,
                   betahihilo_h,betalohilo_h,betahilolo_h,betalololo_h,
                   betahihihi_d,betalohihi_d,betahilohi_d,betalolohi_d,
                   betahihilo_d,betalohilo_d,betahilolo_d,betalololo_d,
                   tileRlapms,addcnt,mulcnt,verbose);
            }
         }
         else
         {
            GPU_dbl8_large_house
               (nrows,ncols,szt,nbt,colidx,nrows1,k,L,
                     Ahihihi_h,     Alohihi_h,     Ahilohi_h,     Alolohi_h,
                     Ahihilo_h,     Alohilo_h,     Ahilolo_h,     Alololo_h,
                     Ahihihi_d,     Alohihi_d,     Ahilohi_d,     Alolohi_d,
                     Ahihilo_d,     Alohilo_d,     Ahilolo_d,     Alololo_d,
                     vhihihi_h,     vlohihi_h,     vhilohi_h,     vlolohi_h,
                     vhihilo_h,     vlohilo_h,     vhilolo_h,     vlololo_h,
                     Vhihihi_d,     Vlohihi_d,     Vhilohi_d,     Vlolohi_d,
                     Vhihilo_d,     Vlohilo_d,     Vhilolo_d,     Vlololo_d,
                  betahihihi_h,  betalohihi_h,  betahilohi_h,  betalolohi_h,
                  betahihilo_h,  betalohilo_h,  betahilolo_h,  betalololo_h,
                  betahihihi_d,  betalohihi_d,  betahilohi_d,  betalolohi_d,
                  betahihilo_d,  betalohilo_d,  betahilolo_d,  betalololo_d,
                  sumshihihi_h,  sumslohihi_h,  sumshilohi_h,  sumslolohi_h,
                  sumshihilo_h,  sumslohilo_h,  sumshilolo_h,  sumslololo_h,
                  sumshihihi_d,  sumslohihi_d,  sumshilohi_d,  sumslolohi_d,
                  sumshihilo_d,  sumslohilo_d,  sumshilolo_d,  sumslololo_d,
                &sigmahihihi_h,&sigmalohihi_h,&sigmahilohi_h,&sigmalolohi_h,
                &sigmahihilo_h,&sigmalohilo_h,&sigmahilolo_h,&sigmalololo_h,
                 sigmahihihi_d, sigmalohihi_d, sigmahilohi_d, sigmalolohi_d,
                 sigmahihilo_d, sigmalohilo_d, sigmahilolo_d, sigmalololo_d,
                houselapms,addcnt,mulcnt,divcnt,sqrtcnt,verbose);

            if((betahihihi_h[L] == 0.0) && (betalohihi_h[L] == 0.0) &&
               (betahilohi_h[L] == 0.0) && (betalolohi_h[L] == 0.0) &&
               (betahihilo_h[L] == 0.0) && (betalohilo_h[L] == 0.0) &&
               (betahilolo_h[L] == 0.0) && (betalololo_h[L] == 0.0))
            {
               if(verbose) cout << "Zero beta detected." << endl;
            }
            else
            {
               GPU_dbl8_medium_leftRupdate
                  (nrows,ncols,szt,colidx,k,L,
                      Ahihihi_h,     Alohihi_h,     Ahilohi_h,     Alolohi_h,
                      Ahihilo_h,     Alohilo_h,     Ahilolo_h,     Alololo_h,
                      Ahihihi_d,     Alohihi_d,     Ahilohi_d,     Alolohi_d,
                      Ahihilo_d,     Alohilo_d,     Ahilolo_d,     Alololo_d,
                      Vhihihi_d,     Vlohihi_d,     Vhilohi_d,     Vlolohi_d,
                      Vhihilo_d,     Vlohilo_d,     Vhilolo_d,     Vlololo_d,
                   betahihihi_h,  betalohihi_h,  betahilohi_h,  betalolohi_h,
                   betahihilo_h,  betalohilo_h,  betahilolo_h,  betalololo_h,
                   betahihihi_d,  betalohihi_d,  betahilohi_d,  betalolohi_d,
                   betahihilo_d,  betalohilo_d,  betahilolo_d,  betalololo_d,
                 RTdotvhihihi_h,RTdotvlohihi_h,RTdotvhilohi_h,RTdotvlolohi_h,
                 RTdotvhihilo_h,RTdotvlohilo_h,RTdotvhilolo_h,RTdotvlololo_h,
                 RTdotvhihihi_d,RTdotvlohihi_d,RTdotvhilohi_d,RTdotvlolohi_d,
                 RTdotvhihilo_d,RTdotvlohilo_d,RTdotvhilolo_d,RTdotvlololo_d,
                   bRTvhihihi_h,  bRTvlohihi_h,  bRTvhilohi_h,  bRTvlolohi_h,
                   bRTvhihilo_h,  bRTvlohilo_h,  bRTvhilolo_h,  bRTvlololo_h,
                   bRTvhihihi_d,  bRTvlohihi_d,  bRTvhilohi_d,  bRTvlolohi_d,
                   bRTvhihilo_d,  bRTvlohilo_d,  bRTvhilolo_d,  bRTvlololo_d,
                 RTvlapms,tileRlapms,addcnt,mulcnt,verbose);
            }
         }
      }
      GPU_dbl8_medium_VB_to_W
         (nrows,szt,szt,k,
             Vhihihi_h,   Vlohihi_h,   Vhilohi_h,   Vlolohi_h,
             Vhihilo_h,   Vlohilo_h,   Vhilolo_h,   Vlololo_h,
             Vhihihi_d,   Vlohihi_d,   Vhilohi_d,   Vlolohi_d,
             Vhihilo_d,   Vlohilo_d,   Vhilolo_d,   Vlololo_d,
             Whihihi_h,   Wlohihi_h,   Whilohi_h,   Wlolohi_h,
             Whihilo_h,   Wlohilo_h,   Whilolo_h,   Wlololo_h,
             Whihihi_d,   Wlohihi_d,   Whilohi_d,   Wlolohi_d,
             Whihilo_d,   Wlohilo_d,   Whilolo_d,   Wlololo_d,
           WYThihihi_h, WYTlohihi_h, WYThilohi_h, WYTlolohi_h,
           WYThihilo_h, WYTlohilo_h, WYThilolo_h, WYTlololo_h,
           WYThihihi_d, WYTlohihi_d, WYThilohi_d, WYTlolohi_d,
           WYThihilo_d, WYTlohilo_d, WYThilolo_d, WYTlololo_d,
          betahihihi_h,betalohihi_h,betahilohi_h,betalolohi_h,
          betahihilo_h,betalohilo_h,betahilolo_h,betalololo_h,
          betahihihi_d,betalohihi_d,betahilohi_d,betalolohi_d,
          betahihilo_d,betalohilo_d,betahilolo_d,betalololo_d,
          vb2Wlapms,addcnt,mulcnt,verbose);
/*
      GPU_dbl2_small_WYT
         (nrows-k*szt,szt,Whi_d,Wlo_d,Vhi_d,Vlo_d,WYThi_d,WYTlo_d,
          WYThi_h,WYTlo_h,WYTlapms,verbose);
 */
      GPU_dbl8_small_QWYT
         (nrows,szt,k,
             Qhihihi_d,   Qlohihi_d,   Qhilohi_d,   Qlolohi_d,
             Qhihilo_d,   Qlohilo_d,   Qhilolo_d,   Qlololo_d,
           WYThihihi_d, WYTlohihi_d, WYThilohi_d, WYTlolohi_d,
           WYThihilo_d, WYTlohilo_d, WYThilolo_d, WYTlololo_d,
          QWYThihihi_d,QWYTlohihi_d,QWYThilohi_d,QWYTlolohi_d,
          QWYThihilo_d,QWYTlohilo_d,QWYThilolo_d,QWYTlololo_d,
          QWYThihihi_h,QWYTlohihi_h,QWYThilohi_h,QWYTlolohi_h,
          QWYThihilo_h,QWYTlohilo_h,QWYThilolo_h,QWYTlololo_h,
             Qhihihi_h,   Qlohihi_h,   Qhilohi_h,   Qlolohi_h,
             Qhihilo_h,   Qlohilo_h,   Qhilolo_h,   Qlololo_h,
          QWYTlapms,addcnt,mulcnt,verbose);

      GPU_dbl8_small_Qupdate
         (nrows,szt,k,
             Qhihihi_d,   Qlohihi_d,   Qhilohi_d,   Qlolohi_d,
             Qhihilo_d,   Qlohilo_d,   Qhilolo_d,   Qlololo_d,
          QWYThihihi_d,QWYTlohihi_d,QWYThilohi_d,QWYTlolohi_d,
          QWYThihilo_d,QWYTlohilo_d,QWYThilolo_d,QWYTlololo_d,
             Qhihihi_h,   Qlohihi_h,   Qhilohi_h,   Qlolohi_h,
             Qhihilo_h,   Qlohilo_h,   Qhilolo_h,   Qlololo_h,
          Qaddlapms,addcnt,verbose);

      if(k < nbt-1)                                           // update R
      {
         GPU_dbl8_small_YWT
            (nrows,szt,k,
               Vhihihi_d,  Vlohihi_d,  Vhilohi_d,  Vlolohi_d,
               Vhihilo_d,  Vlohilo_d,  Vhilolo_d,  Vlololo_d,
               Whihihi_d,  Wlohihi_d,  Whilohi_d,  Wlolohi_d,
               Whihilo_d,  Wlohilo_d,  Whilolo_d,  Wlololo_d,
             YWThihihi_d,YWTlohihi_d,YWThilohi_d,YWTlolohi_d,
             YWThihilo_d,YWTlohilo_d,YWThilolo_d,YWTlololo_d,
             YWThihihi_h,YWTlohihi_h,YWThilohi_h,YWTlolohi_h,
             YWThihilo_h,YWTlohilo_h,YWThilolo_h,YWTlololo_h,
             YWTlapms,addcnt,mulcnt,verbose);

         GPU_dbl8_small_YWTC
            (nrows,ncols,szt,k,
              YWThihihi_d, YWTlohihi_d, YWThilohi_d, YWTlolohi_d,
              YWThihilo_d, YWTlohilo_d, YWThilolo_d, YWTlololo_d,
                Ahihihi_d,   Alohihi_d,   Ahilohi_d,   Alolohi_d,
                Ahihilo_d,   Alohilo_d,   Ahilolo_d,   Alololo_d,
             YWTChihihi_d,YWTClohihi_d,YWTChilohi_d,YWTClolohi_d,
             YWTChihilo_d,YWTClohilo_d,YWTChilolo_d,YWTClololo_d,
             YWTChihihi_h,YWTClohihi_h,YWTChilohi_h,YWTClolohi_h,
             YWTChihilo_h,YWTClohilo_h,YWTChilolo_h,YWTClololo_h,
             YWTClapms,addcnt,mulcnt,verbose);

         GPU_dbl8_small_R_add_YWTC
            (nrows,ncols,szt,k,
                Ahihihi_d,   Alohihi_d,   Ahilohi_d,   Alolohi_d,
                Ahihilo_d,   Alohilo_d,   Ahilolo_d,   Alololo_d,
             YWTChihihi_d,YWTClohihi_d,YWTChilohi_d,YWTClolohi_d,
             YWTChihilo_d,YWTClohilo_d,YWTChilolo_d,YWTClololo_d,
                Ahihihi_h,   Alohihi_h,   Ahilohi_h,   Alolohi_h,
                Ahihilo_h,   Alohilo_h,   Ahilolo_h,   Alololo_h,
             Raddlapms,addcnt,verbose);
      }
   }
   gettimeofday(&endtime,0);
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   *walltimesec = seconds + microseconds*1.0e-6;

   cudaMemcpy(Qhihihi_h,Qhihihi_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qlohihi_h,Qlohihi_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qhilohi_h,Qhilohi_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qlolohi_h,Qlolohi_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qhihilo_h,Qhihilo_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qlohilo_h,Qlohilo_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qhilolo_h,Qhilolo_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qlololo_h,Qlololo_d,szWYT,cudaMemcpyDeviceToHost);
   ix = 0;                                           // copy rows of Q
   for(int i=0; i<nrows; i++)
      for(int j=0; j<nrows; j++)
      {
         Qhihihi[i][j] = Qhihihi_h[ix];
         Qlohihi[i][j] = Qlohihi_h[ix];
         Qhilohi[i][j] = Qhilohi_h[ix];
         Qlolohi[i][j] = Qlolohi_h[ix];
         Qhihilo[i][j] = Qhihilo_h[ix];
         Qlohilo[i][j] = Qlohilo_h[ix];
         Qhilolo[i][j] = Qhilolo_h[ix];
         Qlololo[i][j] = Qlololo_h[ix++];
      }

   cudaMemcpy(Ahihihi_h,Ahihihi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Alohihi_h,Alohihi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Ahilohi_h,Ahilohi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Alolohi_h,Alolohi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Ahihilo_h,Ahihilo_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Alohilo_h,Alohilo_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Ahilolo_h,Ahilolo_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Alololo_h,Alololo_d,sznum,cudaMemcpyDeviceToHost);

   for(int i=0; i<nrows; i++)                       // copy columns of R
      for(int j=0; j<ncols; j++)
      {
         Rhihihi[i][j] = Ahihihi_h[j*nrows+i];
         Rlohihi[i][j] = Alohihi_h[j*nrows+i];
         Rhilohi[i][j] = Ahilohi_h[j*nrows+i];
         Rlolohi[i][j] = Alolohi_h[j*nrows+i];
         Rhihilo[i][j] = Ahihilo_h[j*nrows+i];
         Rlohilo[i][j] = Alohilo_h[j*nrows+i];
         Rhilolo[i][j] = Ahilolo_h[j*nrows+i];
         Rlololo[i][j] = Alololo_h[j*nrows+i];
      }

   free(Ahihihi_h); free(Alohihi_h); free(Ahilohi_h); free(Alolohi_h);
   free(Ahihilo_h); free(Alohilo_h); free(Ahilolo_h); free(Alololo_h);
   free(Qhihihi_h); free(Qlohihi_h); free(Qhilohi_h); free(Qlolohi_h); 
   free(Qhihilo_h); free(Qlohilo_h); free(Qhilolo_h); free(Qlololo_h); 
   free(vhihihi_h); free(vlohihi_h); free(vhilohi_h); free(vlolohi_h);
   free(vhihilo_h); free(vlohilo_h); free(vhilolo_h); free(vlololo_h);
   free(Vhihihi_h); free(Vlohihi_h); free(Vhilohi_h); free(Vlolohi_h);
   free(Vhihilo_h); free(Vlohilo_h); free(Vhilolo_h); free(Vlololo_h);
   free(Whihihi_h); free(Wlohihi_h); free(Whilohi_h); free(Wlolohi_h);
   free(Whihilo_h); free(Wlohilo_h); free(Whilolo_h); free(Wlololo_h);
   free(sumshihihi_h); free(sumslohihi_h);
   free(sumshilohi_h); free(sumslolohi_h);
   free(sumshihilo_h); free(sumslohilo_h);
   free(sumshilolo_h); free(sumslololo_h);

   free(RTdotvhihihi_h); free(RTdotvlohihi_h);
   free(RTdotvhilohi_h); free(RTdotvlolohi_h);
   free(RTdotvhihilo_h); free(RTdotvlohilo_h);
   free(RTdotvhilolo_h); free(RTdotvlololo_h);
   free(bRTvhihihi_h); free(bRTvlohihi_h);
   free(bRTvhilohi_h); free(bRTvlolohi_h);
   free(bRTvhihilo_h); free(bRTvlohilo_h);
   free(bRTvhilolo_h); free(bRTvlololo_h);
   free(WYThihihi_h); free(QWYThihihi_h);
   free(WYThilohi_h); free(QWYThilohi_h);
   free(WYThihilo_h); free(QWYThihilo_h);
   free(WYThilolo_h); free(QWYThilolo_h);
   free(YWThihihi_h); free(YWTChihihi_h);
   free(YWThilohi_h); free(YWTChilohi_h);
   free(YWThihilo_h); free(YWTChihilo_h);
   free(YWThilolo_h); free(YWTChilolo_h);
   free(WYTlohihi_h); free(QWYTlohihi_h);
   free(WYTlolohi_h); free(QWYTlolohi_h);
   free(WYTlohilo_h); free(QWYTlohilo_h);
   free(WYTlololo_h); free(QWYTlololo_h);
   free(YWTlohihi_h); free(YWTClohihi_h);
   free(YWTlolohi_h); free(YWTClolohi_h);
   free(YWTlohilo_h); free(YWTClohilo_h);
   free(YWTlololo_h); free(YWTClololo_h);
}

void GPU_cmplx8_blocked_houseqr
 ( int nrows, int ncols, int szt, int nbt,
   double **Arehihihi, double **Arelohihi,
   double **Arehilohi, double **Arelolohi,
   double **Arehihilo, double **Arelohilo,
   double **Arehilolo, double **Arelololo,
   double **Aimhihihi, double **Aimlohihi,
   double **Aimhilohi, double **Aimlolohi,
   double **Aimhihilo, double **Aimlohilo,
   double **Aimhilolo, double **Aimlololo,
   double **Qrehihihi, double **Qrelohihi,
   double **Qrehilohi, double **Qrelolohi,
   double **Qrehihilo, double **Qrelohilo,
   double **Qrehilolo, double **Qrelololo,
   double **Qimhihihi, double **Qimlohihi,
   double **Qimhilohi, double **Qimlolohi,
   double **Qimhihilo, double **Qimlohilo,
   double **Qimhilolo, double **Qimlololo,
   double **Rrehihihi, double **Rrelohihi,
   double **Rrehilohi, double **Rrelolohi,
   double **Rrehihilo, double **Rrelohilo,
   double **Rrehilolo, double **Rrelololo,
   double **Rimhihihi, double **Rimlohihi,
   double **Rimhilohi, double **Rimlolohi,
   double **Rimhihilo, double **Rimlohilo,
   double **Rimhilolo, double **Rimlololo,
   double *houselapms, double *RHvlapms, double *tileRlapms,
   double *vb2Wlapms, double *WYHlapms, double *QWYHlapms, double *Qaddlapms,
   double *YWHlapms, double *YWHClapms, double *Raddlapms,
   double *walltimesec, long long int *addcnt, long long int *mulcnt,
   long long int *divcnt, long long int *sqrtcnt, bool verbose )
{
   const int dim = nrows*ncols;        // total number of doubles
   const int nrows2 = nrows*nrows;
   double *Arehihihi_h = new double[dim];  // the real parts of A
   double *Arelohihi_h = new double[dim];
   double *Arehilohi_h = new double[dim];
   double *Arelolohi_h = new double[dim];
   double *Arehihilo_h = new double[dim];
   double *Arelohilo_h = new double[dim];
   double *Arehilolo_h = new double[dim];
   double *Arelololo_h = new double[dim];
   double *Aimhihihi_h = new double[dim];  // the imaginary parts of A
   double *Aimlohihi_h = new double[dim]; 
   double *Aimhilohi_h = new double[dim]; 
   double *Aimlolohi_h = new double[dim]; 
   double *Aimhihilo_h = new double[dim];
   double *Aimlohilo_h = new double[dim]; 
   double *Aimhilolo_h = new double[dim]; 
   double *Aimlololo_h = new double[dim]; 
   double *Arehihihi_d;                    // Are on the device
   double *Arelohihi_d;
   double *Arehilohi_d;
   double *Arelolohi_d;
   double *Arehihilo_d;
   double *Arelohilo_d;
   double *Arehilolo_d;
   double *Arelololo_d;
   double *Aimhihihi_d;                    // Aim on the device
   double *Aimlohihi_d;
   double *Aimhilohi_d;
   double *Aimlolohi_d;
   double *Aimhihilo_d;
   double *Aimlohilo_d;
   double *Aimhilolo_d;
   double *Aimlololo_d;
   double *Qrehihihi_h = new double[nrows2]; // real parts of Q 
   double *Qrelohihi_h = new double[nrows2];
   double *Qrehilohi_h = new double[nrows2];
   double *Qrelolohi_h = new double[nrows2];
   double *Qrehihilo_h = new double[nrows2];
   double *Qrelohilo_h = new double[nrows2];
   double *Qrehilolo_h = new double[nrows2];
   double *Qrelololo_h = new double[nrows2];
   double *Qimhihihi_h = new double[nrows2]; // imaginary parts of Q
   double *Qimlohihi_h = new double[nrows2];
   double *Qimhilohi_h = new double[nrows2];
   double *Qimlolohi_h = new double[nrows2];
   double *Qimhihilo_h = new double[nrows2];
   double *Qimlohilo_h = new double[nrows2];
   double *Qimhilolo_h = new double[nrows2];
   double *Qimlololo_h = new double[nrows2];
   double *Qrehihihi_d;                      // Qre on the device
   double *Qrelohihi_d;
   double *Qrehilohi_d;
   double *Qrelolohi_d;
   double *Qrehihilo_d;
   double *Qrelohilo_d;
   double *Qrehilolo_d;
   double *Qrelololo_d;
   double *Qimhihihi_d;                      // Qim on the device
   double *Qimlohihi_d;
   double *Qimhilohi_d;
   double *Qimlolohi_d;
   double *Qimhihilo_d;
   double *Qimlohilo_d;
   double *Qimhilolo_d;
   double *Qimlololo_d;
   double *vrehihihi_h = new double[nrows];  // real parts of Householder v
   double *vrelohihi_h = new double[nrows];
   double *vrehilohi_h = new double[nrows];
   double *vrelolohi_h = new double[nrows];
   double *vrehihilo_h = new double[nrows];
   double *vrelohilo_h = new double[nrows];
   double *vrehilolo_h = new double[nrows];
   double *vrelololo_h = new double[nrows];
   double *vimhihihi_h = new double[nrows];  // imag parts of Householder v
   double *vimlohihi_h = new double[nrows]; 
   double *vimhilohi_h = new double[nrows]; 
   double *vimlolohi_h = new double[nrows]; 
   double *vimhihilo_h = new double[nrows];
   double *vimlohilo_h = new double[nrows]; 
   double *vimhilolo_h = new double[nrows]; 
   double *vimlololo_h = new double[nrows]; 
   double *betahihihi_h = new double[szt+1]; // beta
   double *betalohihi_h = new double[szt+1]; 
   double *betahilohi_h = new double[szt+1]; 
   double *betalolohi_h = new double[szt+1]; 
   double *betahihilo_h = new double[szt+1];
   double *betalohilo_h = new double[szt+1]; 
   double *betahilolo_h = new double[szt+1]; 
   double *betalololo_h = new double[szt+1]; 
   double *betahihihi_d;                     // beta on the device
   double *betalohihi_d;
   double *betahilohi_d;
   double *betalolohi_d;
   double *betahihilo_d;
   double *betalohilo_d;
   double *betahilolo_d;
   double *betalololo_d;
   double *Vrehihihi_h = new double[nrows*szt]; // real parts of V
   double *Vrelohihi_h = new double[nrows*szt]; 
   double *Vrehilohi_h = new double[nrows*szt]; 
   double *Vrelolohi_h = new double[nrows*szt]; 
   double *Vrehihilo_h = new double[nrows*szt];
   double *Vrelohilo_h = new double[nrows*szt]; 
   double *Vrehilolo_h = new double[nrows*szt]; 
   double *Vrelololo_h = new double[nrows*szt]; 
   double *Vimhihihi_h = new double[nrows*szt]; // imaginary parts of V
   double *Vimlohihi_h = new double[nrows*szt];
   double *Vimhilohi_h = new double[nrows*szt];
   double *Vimlolohi_h = new double[nrows*szt];
   double *Vimhihilo_h = new double[nrows*szt];
   double *Vimlohilo_h = new double[nrows*szt];
   double *Vimhilolo_h = new double[nrows*szt];
   double *Vimlololo_h = new double[nrows*szt];
   double *Vrehihihi_d;                         // Vre on device
   double *Vrelohihi_d;
   double *Vrehilohi_d;
   double *Vrelolohi_d;
   double *Vrehihilo_d;
   double *Vrelohilo_d;
   double *Vrehilolo_d;
   double *Vrelololo_d;
   double *Vimhihihi_d;                         // Vim on device
   double *Vimlohihi_d;
   double *Vimhilohi_d;
   double *Vimlolohi_d;
   double *Vimhihilo_d;
   double *Vimlohilo_d;
   double *Vimhilolo_d;
   double *Vimlololo_d;
   double *Wrehihihi_h = new double[nrows*szt]; // real parts of W
   double *Wrelohihi_h = new double[nrows*szt];
   double *Wrehilohi_h = new double[nrows*szt];
   double *Wrelolohi_h = new double[nrows*szt];
   double *Wrehihilo_h = new double[nrows*szt];
   double *Wrelohilo_h = new double[nrows*szt];
   double *Wrehilolo_h = new double[nrows*szt];
   double *Wrelololo_h = new double[nrows*szt];
   double *Wimhihihi_h = new double[nrows*szt]; // imaginary parts of W
   double *Wimlohihi_h = new double[nrows*szt];
   double *Wimhilohi_h = new double[nrows*szt];
   double *Wimlolohi_h = new double[nrows*szt];
   double *Wimhihilo_h = new double[nrows*szt];
   double *Wimlohilo_h = new double[nrows*szt];
   double *Wimhilolo_h = new double[nrows*szt];
   double *Wimlololo_h = new double[nrows*szt];
   double *Wrehihihi_d;                         // Wre on the device
   double *Wrelohihi_d;
   double *Wrehilohi_d;
   double *Wrelolohi_d;
   double *Wrehihilo_d;
   double *Wrelohilo_d;
   double *Wrehilolo_d;
   double *Wrelololo_d;
   double *Wimhihihi_d;                         // Wim on the device
   double *Wimlohihi_d;
   double *Wimhilohi_d;
   double *Wimlolohi_d;
   double *Wimhihilo_d;
   double *Wimlohilo_d;
   double *Wimhilolo_d;
   double *Wimlololo_d;
   double *WYTrehihihi_h = new double[nrows2];  // real parts of W*Y^H
   double *WYTrelohihi_h = new double[nrows2];
   double *WYTrehilohi_h = new double[nrows2];
   double *WYTrelolohi_h = new double[nrows2];
   double *WYTrehihilo_h = new double[nrows2];
   double *WYTrelohilo_h = new double[nrows2];
   double *WYTrehilolo_h = new double[nrows2];
   double *WYTrelololo_h = new double[nrows2];
   double *WYTimhihihi_h = new double[nrows2];  // imaginary parts of W*Y^H
   double *WYTimlohihi_h = new double[nrows2];
   double *WYTimhilohi_h = new double[nrows2];
   double *WYTimlolohi_h = new double[nrows2];
   double *WYTimhihilo_h = new double[nrows2];
   double *WYTimlohilo_h = new double[nrows2];
   double *WYTimhilolo_h = new double[nrows2];
   double *WYTimlololo_h = new double[nrows2];
   double *WYTrehihihi_d;                       // WYTre on the device 
   double *WYTrelohihi_d;
   double *WYTrehilohi_d;
   double *WYTrelolohi_d;
   double *WYTrehihilo_d;
   double *WYTrelohilo_d;
   double *WYTrehilolo_d;
   double *WYTrelololo_d;
   double *WYTimhihihi_d;                       // WYTim on the device
   double *WYTimlohihi_d;
   double *WYTimhilohi_d;
   double *WYTimlolohi_d;
   double *WYTimhihilo_d;
   double *WYTimlohilo_d;
   double *WYTimhilolo_d;
   double *WYTimlololo_d;
   double *YWTrehihihi_h = new double[nrows2];  // real parts of Y*W^H
   double *YWTrelohihi_h = new double[nrows2];
   double *YWTrehilohi_h = new double[nrows2];
   double *YWTrelolohi_h = new double[nrows2];
   double *YWTrehihilo_h = new double[nrows2];
   double *YWTrelohilo_h = new double[nrows2];
   double *YWTrehilolo_h = new double[nrows2];
   double *YWTrelololo_h = new double[nrows2];
   double *YWTimhihihi_h = new double[nrows2];  // imaginary parts of Y*W^H
   double *YWTimlohihi_h = new double[nrows2];
   double *YWTimhilohi_h = new double[nrows2];
   double *YWTimlolohi_h = new double[nrows2];
   double *YWTimhihilo_h = new double[nrows2];
   double *YWTimlohilo_h = new double[nrows2];
   double *YWTimhilolo_h = new double[nrows2];
   double *YWTimlololo_h = new double[nrows2];
   double *YWTrehihihi_d;                       // YWTre on the device
   double *YWTrelohihi_d;
   double *YWTrehilohi_d;
   double *YWTrelolohi_d;
   double *YWTrehihilo_d;
   double *YWTrelohilo_d;
   double *YWTrehilolo_d;
   double *YWTrelololo_d;
   double *YWTimhihihi_d;                       // YWTim on the device
   double *YWTimlohihi_d; 
   double *YWTimhilohi_d; 
   double *YWTimlolohi_d; 
   double *YWTimhihilo_d;
   double *YWTimlohilo_d; 
   double *YWTimhilolo_d; 
   double *YWTimlololo_d; 
   double *QWYTrehihihi_h = new double[nrows2]; // real parts of Q*WY^H
   double *QWYTrelohihi_h = new double[nrows2];
   double *QWYTrehilohi_h = new double[nrows2];
   double *QWYTrelolohi_h = new double[nrows2];
   double *QWYTrehihilo_h = new double[nrows2];
   double *QWYTrelohilo_h = new double[nrows2];
   double *QWYTrehilolo_h = new double[nrows2];
   double *QWYTrelololo_h = new double[nrows2];
   double *QWYTimhihihi_h = new double[nrows2]; // imaginary parts of Q*WY^H
   double *QWYTimlohihi_h = new double[nrows2];
   double *QWYTimhilohi_h = new double[nrows2];
   double *QWYTimlolohi_h = new double[nrows2];
   double *QWYTimhihilo_h = new double[nrows2];
   double *QWYTimlohilo_h = new double[nrows2];
   double *QWYTimhilolo_h = new double[nrows2];
   double *QWYTimlololo_h = new double[nrows2];
   double *QWYTrehihihi_d;                      // QWYTre on the device
   double *QWYTrelohihi_d;
   double *QWYTrehilohi_d;
   double *QWYTrelolohi_d;
   double *QWYTrehihilo_d;
   double *QWYTrelohilo_d;
   double *QWYTrehilolo_d;
   double *QWYTrelololo_d;
   double *QWYTimhihihi_d;                      // QWYTim on the device
   double *QWYTimlohihi_d;
   double *QWYTimhilohi_d;
   double *QWYTimlolohi_d;
   double *QWYTimhihilo_d;
   double *QWYTimlohilo_d;
   double *QWYTimhilolo_d;
   double *QWYTimlololo_d;
   double *YWTCrehihihi_h = new double[dim];    // real parts of YWT*C
   double *YWTCrelohihi_h = new double[dim]; 
   double *YWTCrehilohi_h = new double[dim]; 
   double *YWTCrelolohi_h = new double[dim]; 
   double *YWTCrehihilo_h = new double[dim];
   double *YWTCrelohilo_h = new double[dim]; 
   double *YWTCrehilolo_h = new double[dim]; 
   double *YWTCrelololo_h = new double[dim]; 
   double *YWTCimhihihi_h = new double[dim];    // imaginary parts of YWT*C
   double *YWTCimlohihi_h = new double[dim];
   double *YWTCimhilohi_h = new double[dim];
   double *YWTCimlolohi_h = new double[dim];
   double *YWTCimhihilo_h = new double[dim];
   double *YWTCimlohilo_h = new double[dim];
   double *YWTCimhilolo_h = new double[dim];
   double *YWTCimlololo_h = new double[dim];
   double *YWTCrehihihi_d;                      // YWTCre on the device
   double *YWTCrelohihi_d;
   double *YWTCrehilohi_d;
   double *YWTCrelolohi_d;
   double *YWTCrehihilo_d;
   double *YWTCrelohilo_d;
   double *YWTCrehilolo_d;
   double *YWTCrelololo_d;
   double *YWTCimhihihi_d;                      // YWTCim on the device
   double *YWTCimlohihi_d;
   double *YWTCimhilohi_d;
   double *YWTCimlolohi_d;
   double *YWTCimhihilo_d;
   double *YWTCimlohilo_d;
   double *YWTCimhilolo_d;
   double *YWTCimlololo_d;
   double *RHdotvrehihihi_h = new double[nrows2]; // real R^H dotted with v
   double *RHdotvrelohihi_h = new double[nrows2]; 
   double *RHdotvrehilohi_h = new double[nrows2]; 
   double *RHdotvrelolohi_h = new double[nrows2]; 
   double *RHdotvrehihilo_h = new double[nrows2];
   double *RHdotvrelohilo_h = new double[nrows2]; 
   double *RHdotvrehilolo_h = new double[nrows2]; 
   double *RHdotvrelololo_h = new double[nrows2]; 
   double *RHdotvimhihihi_h = new double[nrows2]; // imag R^H dotted with v
   double *RHdotvimlohihi_h = new double[nrows2]; 
   double *RHdotvimhilohi_h = new double[nrows2]; 
   double *RHdotvimlolohi_h = new double[nrows2]; 
   double *RHdotvimhihilo_h = new double[nrows2];
   double *RHdotvimlohilo_h = new double[nrows2]; 
   double *RHdotvimhilolo_h = new double[nrows2]; 
   double *RHdotvimlololo_h = new double[nrows2]; 
   double *RHdotvrehihihi_d;                      // RHdotvre on the device
   double *RHdotvrelohihi_d;
   double *RHdotvrehilohi_d;
   double *RHdotvrelolohi_d;
   double *RHdotvrehihilo_d;
   double *RHdotvrelohilo_d;
   double *RHdotvrehilolo_d;
   double *RHdotvrelololo_d;
   double *RHdotvimhihihi_d;                      // RHdotvim on the device
   double *RHdotvimlohihi_d;
   double *RHdotvimhilohi_d;
   double *RHdotvimlolohi_d;
   double *RHdotvimhihilo_d;
   double *RHdotvimlohilo_d;
   double *RHdotvimhilolo_d;
   double *RHdotvimlololo_d;
   double *bRHvrehihihi_h = new double[nrows];  // real parts of beta*R^H*v
   double *bRHvrelohihi_h = new double[nrows];
   double *bRHvrehilohi_h = new double[nrows];
   double *bRHvrelolohi_h = new double[nrows];
   double *bRHvrehihilo_h = new double[nrows];
   double *bRHvrelohilo_h = new double[nrows];
   double *bRHvrehilolo_h = new double[nrows];
   double *bRHvrelololo_h = new double[nrows];
   double *bRHvimhihihi_h = new double[nrows];  // imag parts of beta*R^H*v
   double *bRHvimlohihi_h = new double[nrows];
   double *bRHvimhilohi_h = new double[nrows];
   double *bRHvimlolohi_h = new double[nrows];
   double *bRHvimhihilo_h = new double[nrows];
   double *bRHvimlohilo_h = new double[nrows];
   double *bRHvimhilolo_h = new double[nrows];
   double *bRHvimlololo_h = new double[nrows];
   double *bRHvrehihihi_d;                      // bRHvre on the device
   double *bRHvrelohihi_d;
   double *bRHvrehilohi_d;
   double *bRHvrelolohi_d;
   double *bRHvrehihilo_d;
   double *bRHvrelohilo_d;
   double *bRHvrehilolo_d;
   double *bRHvrelololo_d;
   double *bRHvimhihihi_d;                      // bRHvim on the device
   double *bRHvimlohihi_d; 
   double *bRHvimhilohi_d; 
   double *bRHvimlolohi_d; 
   double *bRHvimhihilo_d;
   double *bRHvimlohilo_d; 
   double *bRHvimhilolo_d; 
   double *bRHvimlololo_d; 
   double *sumshihihi_h = new double[nrows];  // subsums for large house
   double *sumslohihi_h = new double[nrows]; 
   double *sumshilohi_h = new double[nrows]; 
   double *sumslolohi_h = new double[nrows]; 
   double *sumshihilo_h = new double[nrows];
   double *sumslohilo_h = new double[nrows]; 
   double *sumshilolo_h = new double[nrows]; 
   double *sumslololo_h = new double[nrows]; 
   double *sumshihihi_d;                      // sums on the device
   double *sumslohihi_d;
   double *sumshilohi_d;
   double *sumslolohi_d;
   double *sumshihilo_d;
   double *sumslohilo_d;
   double *sumshilolo_d;
   double *sumslololo_d;
   double sigmahihihi_h,sigmalohihi_h,sigmahilohi_h,sigmalolohi_h;
   double sigmahihilo_h,sigmalohilo_h,sigmahilolo_h,sigmalololo_h;
   double *sigmahihihi_d;                     // sigma on the device
   double *sigmalohihi_d;
   double *sigmahilohi_d;
   double *sigmalolohi_d;
   double *sigmahihilo_d;
   double *sigmalohilo_d;
   double *sigmahilolo_d;
   double *sigmalololo_d;

   int ix = 0;                            // copy the columns of A to A_h
   for(int j=0; j<ncols; j++)   
      for(int i=0; i<nrows; i++)
      {
         Arehihihi_h[ix]   = Arehihihi[i][j];
         Arelohihi_h[ix]   = Arelohihi[i][j];
         Arehilohi_h[ix]   = Arehilohi[i][j];
         Arelolohi_h[ix]   = Arelolohi[i][j];
         Arehihilo_h[ix]   = Arehihilo[i][j];
         Arelohilo_h[ix]   = Arelohilo[i][j];
         Arehilolo_h[ix]   = Arehilolo[i][j];
         Arelololo_h[ix]   = Arelololo[i][j];
         Aimhihihi_h[ix]   = Aimhihihi[i][j];
         Aimlohihi_h[ix]   = Aimlohihi[i][j];
         Aimhilohi_h[ix]   = Aimhilohi[i][j];
         Aimlolohi_h[ix]   = Aimlolohi[i][j];
         Aimhihilo_h[ix]   = Aimhihilo[i][j];
         Aimlohilo_h[ix]   = Aimlohilo[i][j];
         Aimhilolo_h[ix]   = Aimhilolo[i][j];
         Aimlololo_h[ix++] = Aimlololo[i][j];
      }

   ix = 0;                                // initialize Q with identity
   for(int i=0; i<nrows; i++)
   {
      for(int j=0; j<nrows; j++)
      {
         if(i == j)
         {
            Qrehihihi_h[ix]   = 1.0;
            Qrelohihi_h[ix]   = 0.0;
            Qrehilohi_h[ix]   = 0.0;
            Qrelolohi_h[ix]   = 0.0;
            Qrehihilo_h[ix]   = 0.0;
            Qrelohilo_h[ix]   = 0.0;
            Qrehilolo_h[ix]   = 0.0;
            Qrelololo_h[ix]   = 0.0;
            Qimhihihi_h[ix]   = 0.0;
            Qimlohihi_h[ix]   = 0.0;
            Qimhilohi_h[ix]   = 0.0;
            Qimlolohi_h[ix]   = 0.0;
            Qimhihilo_h[ix]   = 0.0;
            Qimlohilo_h[ix]   = 0.0;
            Qimhilolo_h[ix]   = 0.0;
            Qimlololo_h[ix++] = 0.0;
         }
         else
         {
            Qrehihihi_h[ix]   = 0.0;
            Qrelohihi_h[ix]   = 0.0;
            Qrehilohi_h[ix]   = 0.0;
            Qrelolohi_h[ix]   = 0.0;
            Qrehihilo_h[ix]   = 0.0;
            Qrelohilo_h[ix]   = 0.0;
            Qrehilolo_h[ix]   = 0.0;
            Qrelololo_h[ix]   = 0.0;
            Qimhihihi_h[ix]   = 0.0;
            Qimlohihi_h[ix]   = 0.0;
            Qimhilohi_h[ix]   = 0.0;
            Qimlolohi_h[ix]   = 0.0;
            Qimhihilo_h[ix]   = 0.0;
            Qimlohilo_h[ix]   = 0.0;
            Qimhilolo_h[ix]   = 0.0;
            Qimlololo_h[ix++] = 0.0;
         }
      }
   }
   const size_t sznum = dim*sizeof(double);
   cudaMalloc((void**)&Arehihihi_d,sznum);
   cudaMalloc((void**)&Arelohihi_d,sznum);
   cudaMalloc((void**)&Arehilohi_d,sznum);
   cudaMalloc((void**)&Arelolohi_d,sznum);
   cudaMalloc((void**)&Arehihilo_d,sznum);
   cudaMalloc((void**)&Arelohilo_d,sznum);
   cudaMalloc((void**)&Arehilolo_d,sznum);
   cudaMalloc((void**)&Arelololo_d,sznum);
   cudaMalloc((void**)&Aimhihihi_d,sznum);
   cudaMalloc((void**)&Aimlohihi_d,sznum);
   cudaMalloc((void**)&Aimhilohi_d,sznum);
   cudaMalloc((void**)&Aimlolohi_d,sznum);
   cudaMalloc((void**)&Aimhihilo_d,sznum);
   cudaMalloc((void**)&Aimlohilo_d,sznum);
   cudaMalloc((void**)&Aimhilolo_d,sznum);
   cudaMalloc((void**)&Aimlololo_d,sznum);
   cudaMemcpy(Arehihihi_d,Arehihihi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Arelohihi_d,Arelohihi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Arehilohi_d,Arehilohi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Arelolohi_d,Arelolohi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Arehihilo_d,Arehihilo_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Arelohilo_d,Arelohilo_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Arehilolo_d,Arehilolo_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Arelololo_d,Arelololo_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Aimhihihi_d,Aimhihihi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Aimlohihi_d,Aimlohihi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Aimhilohi_d,Aimhilohi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Aimlolohi_d,Aimlolohi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Aimhihilo_d,Aimhihilo_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Aimlohilo_d,Aimlohilo_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Aimhilolo_d,Aimhilolo_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Aimlololo_d,Aimlololo_h,sznum,cudaMemcpyHostToDevice);
   // allocate one extra beta for use in the cmplx8_normalize
   const size_t szbeta = (szt+1)*sizeof(double);
   cudaMalloc((void**)&betahihihi_d,szbeta);
   cudaMalloc((void**)&betalohihi_d,szbeta);
   cudaMalloc((void**)&betahilohi_d,szbeta);
   cudaMalloc((void**)&betalolohi_d,szbeta);
   cudaMalloc((void**)&betahihilo_d,szbeta);
   cudaMalloc((void**)&betalohilo_d,szbeta);
   cudaMalloc((void**)&betahilolo_d,szbeta);
   cudaMalloc((void**)&betalololo_d,szbeta);

   for(int i=0; i<=szt; i++)
   {
      betahihihi_h[i] = 0.0;
      betalohihi_h[i] = 0.0;
      betahilohi_h[i] = 0.0;
      betalolohi_h[i] = 0.0;
      betahihilo_h[i] = 0.0;
      betalohilo_h[i] = 0.0;
      betahilolo_h[i] = 0.0;
      betalololo_h[i] = 0.0;
   }
   cudaMemcpy(betahihihi_d,betahihihi_h,szbeta,cudaMemcpyHostToDevice);
   cudaMemcpy(betalohihi_d,betalohihi_h,szbeta,cudaMemcpyHostToDevice);
   cudaMemcpy(betahilohi_d,betahilohi_h,szbeta,cudaMemcpyHostToDevice);
   cudaMemcpy(betalolohi_d,betalolohi_h,szbeta,cudaMemcpyHostToDevice);
   cudaMemcpy(betahihilo_d,betahihilo_h,szbeta,cudaMemcpyHostToDevice);
   cudaMemcpy(betalohilo_d,betalohilo_h,szbeta,cudaMemcpyHostToDevice);
   cudaMemcpy(betahilolo_d,betahilolo_h,szbeta,cudaMemcpyHostToDevice);
   cudaMemcpy(betalololo_d,betalololo_h,szbeta,cudaMemcpyHostToDevice);

   const size_t szhouse = nrows*sizeof(double);
   const size_t szpad = szt*sizeof(double);    // padding for nonsquare tiles
   const size_t szVandW = szt*szhouse;
   cudaMalloc((void**)&Vrehihihi_d,szVandW + szpad); // padding added
   cudaMalloc((void**)&Vrelohihi_d,szVandW + szpad);
   cudaMalloc((void**)&Vrehilohi_d,szVandW + szpad);
   cudaMalloc((void**)&Vrelolohi_d,szVandW + szpad);
   cudaMalloc((void**)&Vrehihilo_d,szVandW + szpad);
   cudaMalloc((void**)&Vrelohilo_d,szVandW + szpad);
   cudaMalloc((void**)&Vrehilolo_d,szVandW + szpad);
   cudaMalloc((void**)&Vrelololo_d,szVandW + szpad);
   cudaMalloc((void**)&Vimhihihi_d,szVandW + szpad);
   cudaMalloc((void**)&Vimlohihi_d,szVandW + szpad);
   cudaMalloc((void**)&Vimhilohi_d,szVandW + szpad);
   cudaMalloc((void**)&Vimlolohi_d,szVandW + szpad);
   cudaMalloc((void**)&Vimhihilo_d,szVandW + szpad);
   cudaMalloc((void**)&Vimlohilo_d,szVandW + szpad);
   cudaMalloc((void**)&Vimhilolo_d,szVandW + szpad);
   cudaMalloc((void**)&Vimlololo_d,szVandW + szpad);

   ix = 0;
   for(int i=0; i<nrows*szt; i++)
   {
      Vrehihihi_h[ix]   = 0.0; 
      Vrelohihi_h[ix]   = 0.0; 
      Vrehilohi_h[ix]   = 0.0; 
      Vrelolohi_h[ix]   = 0.0; 
      Vrehihilo_h[ix]   = 0.0; 
      Vrelohilo_h[ix]   = 0.0; 
      Vrehilolo_h[ix]   = 0.0; 
      Vrelololo_h[ix]   = 0.0; 
      Vimhihihi_h[ix]   = 0.0; 
      Vimlohihi_h[ix]   = 0.0; 
      Vimhilohi_h[ix]   = 0.0; 
      Vimlolohi_h[ix]   = 0.0; 
      Vimhihilo_h[ix]   = 0.0; 
      Vimlohilo_h[ix]   = 0.0; 
      Vimhilolo_h[ix]   = 0.0; 
      Vimlololo_h[ix++] = 0.0; 
   }
   Vrehihihi_h[--ix] = 1.0; // initialize last vector for square tiles

   cudaMemcpy(Vrehihihi_d,Vrehihihi_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vrelohihi_d,Vrelohihi_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vrehilohi_d,Vrehilohi_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vrelolohi_d,Vrelolohi_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vrehihilo_d,Vrehihilo_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vrelohilo_d,Vrelohilo_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vrehilolo_d,Vrehilolo_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vrelololo_d,Vrelololo_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vimhihihi_d,Vimhihihi_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vimlohihi_d,Vimlohihi_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vimhilohi_d,Vimhilohi_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vimlolohi_d,Vimlolohi_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vimhihilo_d,Vimhihilo_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vimlohilo_d,Vimlohilo_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vimhilolo_d,Vimhilolo_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vimlololo_d,Vimlololo_h,szVandW,cudaMemcpyHostToDevice);

   cudaMalloc((void**)&Wrehihihi_d,szVandW + szpad); // padding added
   cudaMalloc((void**)&Wrelohihi_d,szVandW + szpad);
   cudaMalloc((void**)&Wrehilohi_d,szVandW + szpad);
   cudaMalloc((void**)&Wrelolohi_d,szVandW + szpad);
   cudaMalloc((void**)&Wrehihilo_d,szVandW + szpad);
   cudaMalloc((void**)&Wrelohilo_d,szVandW + szpad);
   cudaMalloc((void**)&Wrehilolo_d,szVandW + szpad);
   cudaMalloc((void**)&Wrelololo_d,szVandW + szpad);
   cudaMalloc((void**)&Wimhihihi_d,szVandW + szpad);
   cudaMalloc((void**)&Wimlohihi_d,szVandW + szpad);
   cudaMalloc((void**)&Wimhilohi_d,szVandW + szpad);
   cudaMalloc((void**)&Wimlolohi_d,szVandW + szpad);
   cudaMalloc((void**)&Wimhihilo_d,szVandW + szpad);
   cudaMalloc((void**)&Wimlohilo_d,szVandW + szpad);
   cudaMalloc((void**)&Wimhilolo_d,szVandW + szpad);
   cudaMalloc((void**)&Wimlololo_d,szVandW + szpad);

   cudaMalloc((void**)&RHdotvrehihihi_d,szVandW + szpad);
   cudaMalloc((void**)&RHdotvrelohihi_d,szVandW + szpad);
   cudaMalloc((void**)&RHdotvrehilohi_d,szVandW + szpad);
   cudaMalloc((void**)&RHdotvrelolohi_d,szVandW + szpad);
   cudaMalloc((void**)&RHdotvrehihilo_d,szVandW + szpad);
   cudaMalloc((void**)&RHdotvrelohilo_d,szVandW + szpad);
   cudaMalloc((void**)&RHdotvrehilolo_d,szVandW + szpad);
   cudaMalloc((void**)&RHdotvrelololo_d,szVandW + szpad);
   cudaMalloc((void**)&RHdotvimhihihi_d,szVandW + szpad);
   cudaMalloc((void**)&RHdotvimlohihi_d,szVandW + szpad);
   cudaMalloc((void**)&RHdotvimhilohi_d,szVandW + szpad);
   cudaMalloc((void**)&RHdotvimlolohi_d,szVandW + szpad);
   cudaMalloc((void**)&RHdotvimhihilo_d,szVandW + szpad);
   cudaMalloc((void**)&RHdotvimlohilo_d,szVandW + szpad);
   cudaMalloc((void**)&RHdotvimhilolo_d,szVandW + szpad);
   cudaMalloc((void**)&RHdotvimlololo_d,szVandW + szpad);
   cudaMalloc((void**)&bRHvrehihihi_d,szhouse + szpad);
   cudaMalloc((void**)&bRHvrelohihi_d,szhouse + szpad);
   cudaMalloc((void**)&bRHvrehilohi_d,szhouse + szpad);
   cudaMalloc((void**)&bRHvrelolohi_d,szhouse + szpad);
   cudaMalloc((void**)&bRHvrehihilo_d,szhouse + szpad);
   cudaMalloc((void**)&bRHvrelohilo_d,szhouse + szpad);
   cudaMalloc((void**)&bRHvrehilolo_d,szhouse + szpad);
   cudaMalloc((void**)&bRHvrelololo_d,szhouse + szpad);
   cudaMalloc((void**)&bRHvimhihihi_d,szhouse + szpad);
   cudaMalloc((void**)&bRHvimlohihi_d,szhouse + szpad);
   cudaMalloc((void**)&bRHvimhilohi_d,szhouse + szpad);
   cudaMalloc((void**)&bRHvimlolohi_d,szhouse + szpad);
   cudaMalloc((void**)&bRHvimhihilo_d,szhouse + szpad);
   cudaMalloc((void**)&bRHvimlohilo_d,szhouse + szpad);
   cudaMalloc((void**)&bRHvimhilolo_d,szhouse + szpad);
   cudaMalloc((void**)&bRHvimlololo_d,szhouse + szpad);

   cudaMalloc((void**)&sumshihihi_d,szhouse);
   cudaMalloc((void**)&sumslohihi_d,szhouse);
   cudaMalloc((void**)&sumshilohi_d,szhouse);
   cudaMalloc((void**)&sumslolohi_d,szhouse);
   cudaMalloc((void**)&sumshihilo_d,szhouse);
   cudaMalloc((void**)&sumslohilo_d,szhouse);
   cudaMalloc((void**)&sumshilolo_d,szhouse);
   cudaMalloc((void**)&sumslololo_d,szhouse);
   cudaMalloc((void**)&sigmahihihi_d,sizeof(double));
   cudaMalloc((void**)&sigmalohihi_d,sizeof(double));
   cudaMalloc((void**)&sigmahilohi_d,sizeof(double));
   cudaMalloc((void**)&sigmalolohi_d,sizeof(double));
   cudaMalloc((void**)&sigmahihilo_d,sizeof(double));
   cudaMalloc((void**)&sigmalohilo_d,sizeof(double));
   cudaMalloc((void**)&sigmahilolo_d,sizeof(double));
   cudaMalloc((void**)&sigmalololo_d,sizeof(double));

   const size_t szWYT = nrows2*sizeof(double);
   cudaMalloc((void**)&WYTrehihihi_d,szWYT + szpad); // padding for W*Y^H
   cudaMalloc((void**)&WYTrelohihi_d,szWYT + szpad);
   cudaMalloc((void**)&WYTrehilohi_d,szWYT + szpad);
   cudaMalloc((void**)&WYTrelolohi_d,szWYT + szpad);
   cudaMalloc((void**)&WYTrehihilo_d,szWYT + szpad);
   cudaMalloc((void**)&WYTrelohilo_d,szWYT + szpad);
   cudaMalloc((void**)&WYTrehilolo_d,szWYT + szpad);
   cudaMalloc((void**)&WYTrelololo_d,szWYT + szpad);
   cudaMalloc((void**)&WYTimhihihi_d,szWYT + szpad);
   cudaMalloc((void**)&WYTimlohihi_d,szWYT + szpad);
   cudaMalloc((void**)&WYTimhilohi_d,szWYT + szpad);
   cudaMalloc((void**)&WYTimlolohi_d,szWYT + szpad);
   cudaMalloc((void**)&WYTimhihilo_d,szWYT + szpad);
   cudaMalloc((void**)&WYTimlohilo_d,szWYT + szpad);
   cudaMalloc((void**)&WYTimhilolo_d,szWYT + szpad);
   cudaMalloc((void**)&WYTimlololo_d,szWYT + szpad);
   cudaMalloc((void**)&Qrehihihi_d,szWYT + szpad);
   cudaMalloc((void**)&Qrelohihi_d,szWYT + szpad);
   cudaMalloc((void**)&Qrehilohi_d,szWYT + szpad);
   cudaMalloc((void**)&Qrelolohi_d,szWYT + szpad);
   cudaMalloc((void**)&Qrehihilo_d,szWYT + szpad);
   cudaMalloc((void**)&Qrelohilo_d,szWYT + szpad);
   cudaMalloc((void**)&Qrehilolo_d,szWYT + szpad);
   cudaMalloc((void**)&Qrelololo_d,szWYT + szpad);
   cudaMalloc((void**)&Qimhihihi_d,szWYT + szpad);
   cudaMalloc((void**)&Qimlohihi_d,szWYT + szpad);
   cudaMalloc((void**)&Qimhilohi_d,szWYT + szpad);
   cudaMalloc((void**)&Qimlolohi_d,szWYT + szpad);
   cudaMalloc((void**)&Qimhihilo_d,szWYT + szpad);
   cudaMalloc((void**)&Qimlohilo_d,szWYT + szpad);
   cudaMalloc((void**)&Qimhilolo_d,szWYT + szpad);
   cudaMalloc((void**)&Qimlololo_d,szWYT + szpad);
   cudaMemcpy(Qrehihihi_d,Qrehihihi_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qrelohihi_d,Qrelohihi_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qrehilohi_d,Qrehilohi_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qrelolohi_d,Qrelolohi_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qrehihilo_d,Qrehihilo_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qrelohilo_d,Qrelohilo_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qrehilolo_d,Qrehilolo_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qrelololo_d,Qrelololo_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qimhihihi_d,Qimhihihi_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qimlohihi_d,Qimlohihi_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qimhilohi_d,Qimhilohi_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qimlolohi_d,Qimlolohi_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qimhihilo_d,Qimhihilo_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qimlohilo_d,Qimlohilo_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qimhilolo_d,Qimhilolo_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qimlololo_d,Qimlololo_h,szWYT,cudaMemcpyHostToDevice);
   cudaMalloc((void**)&QWYTrehihihi_d,szWYT + szpad);
   cudaMalloc((void**)&QWYTrelohihi_d,szWYT + szpad);
   cudaMalloc((void**)&QWYTrehilohi_d,szWYT + szpad);
   cudaMalloc((void**)&QWYTrelolohi_d,szWYT + szpad);
   cudaMalloc((void**)&QWYTrehihilo_d,szWYT + szpad);
   cudaMalloc((void**)&QWYTrelohilo_d,szWYT + szpad);
   cudaMalloc((void**)&QWYTrehilolo_d,szWYT + szpad);
   cudaMalloc((void**)&QWYTrelololo_d,szWYT + szpad);
   cudaMalloc((void**)&QWYTimhihihi_d,szWYT + szpad);
   cudaMalloc((void**)&QWYTimlohihi_d,szWYT + szpad);
   cudaMalloc((void**)&QWYTimhilohi_d,szWYT + szpad);
   cudaMalloc((void**)&QWYTimlolohi_d,szWYT + szpad);
   cudaMalloc((void**)&QWYTimhihilo_d,szWYT + szpad);
   cudaMalloc((void**)&QWYTimlohilo_d,szWYT + szpad);
   cudaMalloc((void**)&QWYTimhilolo_d,szWYT + szpad);
   cudaMalloc((void**)&QWYTimlololo_d,szWYT + szpad);

   const size_t szYWT = nrows2*sizeof(double);
   cudaMalloc((void**)&YWTrehihihi_d,szYWT + szpad); // padding for Y*W^H
   cudaMalloc((void**)&YWTrelohihi_d,szYWT + szpad);
   cudaMalloc((void**)&YWTrehilohi_d,szYWT + szpad);
   cudaMalloc((void**)&YWTrelolohi_d,szYWT + szpad);
   cudaMalloc((void**)&YWTrehihilo_d,szYWT + szpad);
   cudaMalloc((void**)&YWTrelohilo_d,szYWT + szpad);
   cudaMalloc((void**)&YWTrehilolo_d,szYWT + szpad);
   cudaMalloc((void**)&YWTrelololo_d,szYWT + szpad);
   cudaMalloc((void**)&YWTimhihihi_d,szYWT + szpad);
   cudaMalloc((void**)&YWTimlohihi_d,szYWT + szpad);
   cudaMalloc((void**)&YWTimhilohi_d,szYWT + szpad);
   cudaMalloc((void**)&YWTimlolohi_d,szYWT + szpad);
   cudaMalloc((void**)&YWTimhihilo_d,szYWT + szpad);
   cudaMalloc((void**)&YWTimlohilo_d,szYWT + szpad);
   cudaMalloc((void**)&YWTimhilolo_d,szYWT + szpad);
   cudaMalloc((void**)&YWTimlololo_d,szYWT + szpad);
   cudaMalloc((void**)&YWTCrehihihi_d,sznum + szpad);
   cudaMalloc((void**)&YWTCrelohihi_d,sznum + szpad);
   cudaMalloc((void**)&YWTCrehilohi_d,sznum + szpad);
   cudaMalloc((void**)&YWTCrelolohi_d,sznum + szpad);
   cudaMalloc((void**)&YWTCrehihilo_d,sznum + szpad);
   cudaMalloc((void**)&YWTCrelohilo_d,sznum + szpad);
   cudaMalloc((void**)&YWTCrehilolo_d,sznum + szpad);
   cudaMalloc((void**)&YWTCrelololo_d,sznum + szpad);
   cudaMalloc((void**)&YWTCimhihihi_d,sznum + szpad);
   cudaMalloc((void**)&YWTCimlohihi_d,sznum + szpad);
   cudaMalloc((void**)&YWTCimhilohi_d,sznum + szpad);
   cudaMalloc((void**)&YWTCimlolohi_d,sznum + szpad);
   cudaMalloc((void**)&YWTCimhihilo_d,sznum + szpad);
   cudaMalloc((void**)&YWTCimlohilo_d,sznum + szpad);
   cudaMalloc((void**)&YWTCimhilolo_d,sznum + szpad);
   cudaMalloc((void**)&YWTCimlololo_d,sznum + szpad);

   *houselapms = 0.0; *RHvlapms = 0.0; *tileRlapms = 0.0; *vb2Wlapms = 0.0;
   *WYHlapms = 0.0; *QWYHlapms = 0.0; *Qaddlapms = 0.0;
   *YWHlapms = 0.0; *YWHClapms = 0.0; *Raddlapms = 0.0;
   *addcnt = 0; *mulcnt = 0; *divcnt = 0; *sqrtcnt = 0;
   struct timeval begintime,endtime; // wall clock time of computations

   gettimeofday(&begintime,0);

   for(int k=0; k<nbt; k++)       // k runs over the number of blocks
   {
      if(verbose)
         cout << "Tile k = " << k << " out of " << nbt << " ..." << endl;

      int colidx,nrows1;

      for(int L=0; L<szt; L++)  // L runs over the columns in one block
      {
         colidx = k*szt + L;              // index of the current column
         nrows1 = nrows - colidx - 1;     // #rows in Householder vector - 1
         if(nrows1 <= szt)
         {
            GPU_cmplx8_small_house
               (nrows,ncols,szt,nbt,colidx,nrows1,k,L,
                 Arehihihi_h, Arelohihi_h, Arehilohi_h, Arelolohi_h,
                 Arehihilo_h, Arelohilo_h, Arehilolo_h, Arelololo_h,
                 Aimhihihi_h, Aimlohihi_h, Aimhilohi_h, Aimlolohi_h,
                 Aimhihilo_h, Aimlohilo_h, Aimhilolo_h, Aimlololo_h,
                 Arehihihi_d, Arelohihi_d, Arehilohi_d, Arelolohi_d,
                 Arehihilo_d, Arelohilo_d, Arehilolo_d, Arelololo_d,
                 Aimhihihi_d, Aimlohihi_d, Aimhilohi_d, Aimlolohi_d,
                 Aimhihilo_d, Aimlohilo_d, Aimhilolo_d, Aimlololo_d,
                 vrehihihi_h, vrelohihi_h, vrehilohi_h, vrelolohi_h,
                 vrehihilo_h, vrelohilo_h, vrehilolo_h, vrelololo_h,
                 vimhihihi_h, vimlohihi_h, vimhilohi_h, vimlolohi_h,
                 vimhihilo_h, vimlohilo_h, vimhilolo_h, vimlololo_h,
                 Vrehihihi_d, Vrelohihi_d, Vrehilohi_d, Vrelolohi_d,
                 Vrehihilo_d, Vrelohilo_d, Vrehilolo_d, Vrelololo_d,
                 Vimhihihi_d, Vimlohihi_d, Vimhilohi_d, Vimlolohi_d,
                 Vimhihilo_d, Vimlohilo_d, Vimhilolo_d, Vimlololo_d,
                betahihihi_h,betalohihi_h,betahilohi_h,betalolohi_h,
                betahihilo_h,betalohilo_h,betahilolo_h,betalololo_h,
                betahihihi_d,betalohihi_d,betahilohi_d,betalolohi_d,
                betahihilo_d,betalohilo_d,betahilolo_d,betalololo_d,
                houselapms,addcnt,mulcnt,divcnt,sqrtcnt,verbose);

            if((betahihihi_h[L] == 0.0) && (betalohihi_h[L] == 0.0) &&
               (betahilohi_h[L] == 0.0) && (betalolohi_h[L] == 0.0) &&
               (betahihilo_h[L] == 0.0) && (betalohilo_h[L] == 0.0) &&
               (betahilolo_h[L] == 0.0) && (betalololo_h[L] == 0.0))
            {
               if(verbose) cout << "Zero beta detected." << endl;
            }
            else
            {
               GPU_cmplx8_small_leftRupdate
                  (nrows,ncols,szt,colidx,k,L,
                    Arehihihi_h, Arelohihi_h, Arehilohi_h, Arelolohi_h,
                    Arehihilo_h, Arelohilo_h, Arehilolo_h, Arelololo_h,
                    Aimhihihi_h, Aimlohihi_h, Aimhilohi_h, Aimlolohi_h,
                    Aimhihilo_h, Aimlohilo_h, Aimhilolo_h, Aimlololo_h,
                    Arehihihi_d, Arelohihi_d, Arehilohi_d, Arelolohi_d,
                    Arehihilo_d, Arelohilo_d, Arehilolo_d, Arelololo_d,
                    Aimhihihi_d, Aimlohihi_d, Aimhilohi_d, Aimlolohi_d,
                    Aimhihilo_d, Aimlohilo_d, Aimhilolo_d, Aimlololo_d,
                    Vrehihihi_d, Vrelohihi_d, Vrehilohi_d, Vrelolohi_d,
                    Vrehihilo_d, Vrelohilo_d, Vrehilolo_d, Vrelololo_d,
                    Vimhihihi_d, Vimlohihi_d, Vimhilohi_d, Vimlolohi_d,
                    Vimhihilo_d, Vimlohilo_d, Vimhilolo_d, Vimlololo_d,
                   betahihihi_h,betalohihi_h,betahilohi_h,betalolohi_h,
                   betahihilo_h,betalohilo_h,betahilolo_h,betalololo_h,
                   betahihihi_d,betalohihi_d,betahilohi_d,betalolohi_d,
                   betahihilo_d,betalohilo_d,betahilolo_d,betalololo_d,
                   tileRlapms,addcnt,mulcnt,verbose);
            }
         }
         else
         {
            GPU_cmplx8_large_house
               (nrows,ncols,szt,nbt,colidx,nrows1,k,L,
                   Arehihihi_h,   Arelohihi_h,   Arehilohi_h,   Arelolohi_h,
                   Arehihilo_h,   Arelohilo_h,   Arehilolo_h,   Arelololo_h,
                   Aimhihihi_h,   Aimlohihi_h,   Aimhilohi_h,   Aimlolohi_h,
                   Aimhihilo_h,   Aimlohilo_h,   Aimhilolo_h,   Aimlololo_h,
                   Arehihihi_d,   Arelohihi_d,   Arehilohi_d,   Arelolohi_d,
                   Arehihilo_d,   Arelohilo_d,   Arehilolo_d,   Arelololo_d,
                   Aimhihihi_d,   Aimlohihi_d,   Aimhilohi_d,   Aimlolohi_d,
                   Aimhihilo_d,   Aimlohilo_d,   Aimhilolo_d,   Aimlololo_d,
                   vrehihihi_h,   vrelohihi_h,   vrehilohi_h,   vrelolohi_h,
                   vrehihilo_h,   vrelohilo_h,   vrehilolo_h,   vrelololo_h,
                   vimhihihi_h,   vimlohihi_h,   vimhilohi_h,   vimlolohi_h,
                   vimhihilo_h,   vimlohilo_h,   vimhilolo_h,   vimlololo_h,
                   Vrehihihi_d,   Vrelohihi_d,   Vrehilohi_d,   Vrelolohi_d,
                   Vrehihilo_d,   Vrelohilo_d,   Vrehilolo_d,   Vrelololo_d,
                   Vimhihihi_d,   Vimlohihi_d,   Vimhilohi_d,   Vimlolohi_d,
                   Vimhihilo_d,   Vimlohilo_d,   Vimhilolo_d,   Vimlololo_d,
                  betahihihi_h,  betalohihi_h,  betahilohi_h,  betalolohi_h,
                  betahihilo_h,  betalohilo_h,  betahilolo_h,  betalololo_h,
                  betahihihi_d,  betalohihi_d,  betahilohi_d,  betalolohi_d,
                  betahihilo_d,  betalohilo_d,  betahilolo_d,  betalololo_d,
                  sumshihihi_h,  sumslohihi_h,  sumshilohi_h,  sumslolohi_h,
                  sumshihilo_h,  sumslohilo_h,  sumshilolo_h,  sumslololo_h,
                  sumshihihi_d,  sumslohihi_d,  sumshilohi_d,  sumslolohi_d,
                  sumshihilo_d,  sumslohilo_d,  sumshilolo_d,  sumslololo_d,
                &sigmahihihi_h,&sigmalohihi_h,&sigmahilohi_h,&sigmalolohi_h,
                &sigmahihilo_h,&sigmalohilo_h,&sigmahilolo_h,&sigmalololo_h,
                 sigmahihihi_d, sigmalohihi_d, sigmahilohi_d, sigmalolohi_d,
                 sigmahihilo_d, sigmalohilo_d, sigmahilolo_d, sigmalololo_d,
                houselapms,addcnt,mulcnt,divcnt,sqrtcnt,verbose);

            if((betahihihi_h[L] == 0.0) && (betalohihi_h[L] == 0.0) &&
               (betahilohi_h[L] == 0.0) && (betalolohi_h[L] == 0.0) &&
               (betahihilo_h[L] == 0.0) && (betalohilo_h[L] == 0.0) &&
               (betahilolo_h[L] == 0.0) && (betalololo_h[L] == 0.0))
            {
               if(verbose) cout << "Zero beta detected." << endl;
            }
            else
            {
               GPU_cmplx8_medium_leftRupdate
                  (nrows,ncols,szt,colidx,k,L,
                     Arehihihi_h, Arelohihi_h, Arehilohi_h, Arelolohi_h,
                     Arehihilo_h, Arelohilo_h, Arehilolo_h, Arelololo_h,
                     Aimhihihi_h, Aimlohihi_h, Aimhilohi_h, Aimlolohi_h,
                     Aimhihilo_h, Aimlohilo_h, Aimhilolo_h, Aimlololo_h,
                     Arehihihi_d, Arelohihi_d, Arehilohi_d, Arelolohi_d,
                     Arehihilo_d, Arelohilo_d, Arehilolo_d, Arelololo_d,
                     Aimhihihi_d, Aimlohihi_d, Aimhilohi_d, Aimlolohi_d,
                     Aimhihilo_d, Aimlohilo_d, Aimhilolo_d, Aimlololo_d,
                     Vrehihihi_d, Vrelohihi_d, Vrehilohi_d, Vrelolohi_d,
                     Vrehihilo_d, Vrelohilo_d, Vrehilolo_d, Vrelololo_d,
                     Vimhihihi_d, Vimlohihi_d, Vimhilohi_d, Vimlolohi_d,
                     Vimhihilo_d, Vimlohilo_d, Vimhilolo_d, Vimlololo_d,
                    betahihihi_h,betalohihi_h,betahilohi_h,betalolohi_h,
                    betahihilo_h,betalohilo_h,betahilolo_h,betalololo_h,
                    betahihihi_d,betalohihi_d,betahilohi_d,betalolohi_d,
                    betahihilo_d,betalohilo_d,betahilolo_d,betalololo_d,
                RHdotvrehihihi_h,RHdotvrelohihi_h,
                RHdotvrehilohi_h,RHdotvrelolohi_h,
                RHdotvrehihilo_h,RHdotvrelohilo_h,
                RHdotvrehilolo_h,RHdotvrelololo_h,
                RHdotvimhihihi_h,RHdotvimlohihi_h,
                RHdotvimhilohi_h,RHdotvimlolohi_h,
                RHdotvimhihilo_h,RHdotvimlohilo_h,
                RHdotvimhilolo_h,RHdotvimlololo_h,
                RHdotvrehihihi_d,RHdotvrelohihi_d,
                RHdotvrehilohi_d,RHdotvrelolohi_d,
                RHdotvrehihilo_d,RHdotvrelohilo_d,
                RHdotvrehilolo_d,RHdotvrelololo_d,
                RHdotvimhihihi_d,RHdotvimlohihi_d,
                RHdotvimhilohi_d,RHdotvimlolohi_d,
                RHdotvimhihilo_d,RHdotvimlohilo_d,
                RHdotvimhilolo_d,RHdotvimlololo_d,
                  bRHvrehihihi_h,bRHvrelohihi_h,bRHvrehilohi_h,bRHvrelolohi_h,
                  bRHvrehihilo_h,bRHvrelohilo_h,bRHvrehilolo_h,bRHvrelololo_h,
                  bRHvimhihihi_h,bRHvimlohihi_h,bRHvimhilohi_h,bRHvimlolohi_h,
                  bRHvimhihilo_h,bRHvimlohilo_h,bRHvimhilolo_h,bRHvimlololo_h,
                  bRHvrehihihi_d,bRHvrelohihi_d,bRHvrehilohi_d,bRHvrelolohi_d,
                  bRHvrehihilo_d,bRHvrelohilo_d,bRHvrehilolo_d,bRHvrelololo_d,
                  bRHvimhihihi_d,bRHvimlohihi_d,bRHvimhilohi_d,bRHvimlolohi_d,
                  bRHvimhihilo_d,bRHvimlohilo_d,bRHvimhilolo_d,bRHvimlololo_d,
                RHvlapms,tileRlapms,addcnt,mulcnt,verbose);
            }
         }
      }
      GPU_cmplx8_medium_VB_to_W
         (nrows,szt,szt,k,
            Vrehihihi_h,  Vrelohihi_h,  Vrehilohi_h,  Vrelolohi_h,
            Vrehihilo_h,  Vrelohilo_h,  Vrehilolo_h,  Vrelololo_h,
            Vimhihihi_h,  Vimlohihi_h,  Vimhilohi_h,  Vimlolohi_h,
            Vimhihilo_h,  Vimlohilo_h,  Vimhilolo_h,  Vimlololo_h,
            Vrehihihi_d,  Vrelohihi_d,  Vrehilohi_d,  Vrelolohi_d,
            Vrehihilo_d,  Vrelohilo_d,  Vrehilolo_d,  Vrelololo_d,
            Vimhihihi_d,  Vimlohihi_d,  Vimhilohi_d,  Vimlolohi_d,
            Vimhihilo_d,  Vimlohilo_d,  Vimhilolo_d,  Vimlololo_d,
            Wrehihihi_h,  Wrelohihi_h,  Wrehilohi_h,  Wrelolohi_h,
            Wrehihilo_h,  Wrelohilo_h,  Wrehilolo_h,  Wrelololo_h,
            Wimhihihi_h,  Wimlohihi_h,  Wimhilohi_h,  Wimlolohi_h,
            Wimhihilo_h,  Wimlohilo_h,  Wimhilolo_h,  Wimlololo_h,
            Wrehihihi_d,  Wrelohihi_d,  Wrehilohi_d,  Wrelolohi_d,
            Wrehihilo_d,  Wrelohilo_d,  Wrehilolo_d,  Wrelololo_d,
            Wimhihihi_d,  Wimlohihi_d,  Wimhilohi_d,  Wimlolohi_d,
            Wimhihilo_d,  Wimlohilo_d,  Wimhilolo_d,  Wimlololo_d,
          WYTrehihihi_h,WYTrelohihi_h,WYTrehilohi_h,WYTrelolohi_h,
          WYTrehihilo_h,WYTrelohilo_h,WYTrehilolo_h,WYTrelololo_h,
          WYTimhihihi_h,WYTimlohihi_h,WYTimhilohi_h,WYTimlolohi_h,
          WYTimhihilo_h,WYTimlohilo_h,WYTimhilolo_h,WYTimlololo_h,
          WYTrehihihi_d,WYTrelohihi_d,WYTrehilohi_d,WYTrelolohi_d,
          WYTrehihilo_d,WYTrelohilo_d,WYTrehilolo_d,WYTrelololo_d,
          WYTimhihihi_d,WYTimlohihi_d,WYTimhilohi_d,WYTimlolohi_d,
          WYTimhihilo_d,WYTimlohilo_d,WYTimhilolo_d,WYTimlololo_d,
           betahihihi_h, betalohihi_h, betahilohi_h, betalolohi_h,
           betahihilo_h, betalohilo_h, betahilolo_h, betalololo_h,
           betahihihi_d, betalohihi_d, betahilohi_d, betalolohi_d,
           betahihilo_d, betalohilo_d, betahilolo_d, betalololo_d,
          vb2Wlapms,addcnt,mulcnt,verbose);

/*
      GPU_cmplx2_small_WYH(nrows-k*szt,szt,
            Wrehi_d,  Wrelo_d,  Wimhi_d,  Wimlo_d,
            Vrehi_d,  Vrelo_d,  Vimhi_d,  Vimlo_d,
          WYTrehi_d,WYTrelo_d,WYTimhi_d,WYTimlo_d,
          WYTrehi_h,WYTrelo_h,WYTimhi_h,WYTimlo_h,WYHlapms,verbose);
 */

      GPU_cmplx8_small_QWYH(nrows,szt,k,
             Qrehihihi_d,   Qrelohihi_d,   Qrehilohi_d,   Qrelolohi_d,
             Qrehihilo_d,   Qrelohilo_d,   Qrehilolo_d,   Qrelololo_d,
             Qimhihihi_d,   Qimlohihi_d,   Qimhilohi_d,   Qimlolohi_d,
             Qimhihilo_d,   Qimlohilo_d,   Qimhilolo_d,   Qimlololo_d,
           WYTrehihihi_d, WYTrelohihi_d, WYTrehilohi_d, WYTrelolohi_d,
           WYTrehihilo_d, WYTrelohilo_d, WYTrehilolo_d, WYTrelololo_d,
           WYTimhihihi_d, WYTimlohihi_d, WYTimhilohi_d, WYTimlolohi_d,
           WYTimhihilo_d, WYTimlohilo_d, WYTimhilolo_d, WYTimlololo_d,
          QWYTrehihihi_d,QWYTrelohihi_d,QWYTrehilohi_d,QWYTrelolohi_d,
          QWYTrehihilo_d,QWYTrelohilo_d,QWYTrehilolo_d,QWYTrelololo_d,
          QWYTimhihihi_d,QWYTimlohihi_d,QWYTimhilohi_d,QWYTimlolohi_d,
          QWYTimhihilo_d,QWYTimlohilo_d,QWYTimhilolo_d,QWYTimlololo_d,
          QWYTrehihihi_h,QWYTrelohihi_h,QWYTrehilohi_h,QWYTrelolohi_h,
          QWYTrehihilo_h,QWYTrelohilo_h,QWYTrehilolo_h,QWYTrelololo_h,
          QWYTimhihihi_h,QWYTimlohihi_h,QWYTimhilohi_h,QWYTimlolohi_h,
          QWYTimhihilo_h,QWYTimlohilo_h,QWYTimhilolo_h,QWYTimlololo_h,
             Qrehihihi_h,   Qrelohihi_h,   Qrehilohi_h,   Qrelolohi_h,
             Qrehihilo_h,   Qrelohilo_h,   Qrehilolo_h,   Qrelololo_h,
             Qimhihihi_h,   Qimlohihi_h,   Qimhilohi_h,   Qimlolohi_h,
             Qimhihilo_h,   Qimlohilo_h,   Qimhilolo_h,   Qimlololo_h,
          QWYHlapms,addcnt,mulcnt,verbose);

      GPU_cmplx8_small_Qupdate
         (nrows,szt,k,
             Qrehihihi_d,   Qrelohihi_d,   Qrehilohi_d,   Qrelolohi_d,
             Qrehihilo_d,   Qrelohilo_d,   Qrehilolo_d,   Qrelololo_d,
             Qimhihihi_d,   Qimlohihi_d,   Qimhilohi_d,   Qimlolohi_d,
             Qimhihilo_d,   Qimlohilo_d,   Qimhilolo_d,   Qimlololo_d,
          QWYTrehihihi_d,QWYTrelohihi_d,QWYTrehilohi_d,QWYTrelolohi_d,
          QWYTrehihilo_d,QWYTrelohilo_d,QWYTrehilolo_d,QWYTrelololo_d,
          QWYTimhihihi_d,QWYTimlohihi_d,QWYTimhilohi_d,QWYTimlolohi_d,
          QWYTimhihilo_d,QWYTimlohilo_d,QWYTimhilolo_d,QWYTimlololo_d,
             Qrehihihi_h,   Qrelohihi_h,   Qrehilohi_h,   Qrelolohi_h,
             Qrehihilo_h,   Qrelohilo_h,   Qrehilolo_h,   Qrelololo_h,
             Qimhihihi_h,   Qimlohihi_h,   Qimhilohi_h,   Qimlolohi_h,
             Qimhihilo_h,   Qimlohilo_h,   Qimhilolo_h,   Qimlololo_h,
          Qaddlapms,addcnt,verbose);

      if(k < nbt-1)                              // update R
      {
         GPU_cmplx8_small_YWH
            (nrows,szt,k,
               Vrehihihi_d,  Vrelohihi_d,  Vrehilohi_d,  Vrelolohi_d,
               Vrehihilo_d,  Vrelohilo_d,  Vrehilolo_d,  Vrelololo_d,
               Vimhihihi_d,  Vimlohihi_d,  Vimhilohi_d,  Vimlolohi_d,
               Vimhihilo_d,  Vimlohilo_d,  Vimhilolo_d,  Vimlololo_d,
               Wrehihihi_d,  Wrelohihi_d,  Wrehilohi_d,  Wrelolohi_d,
               Wrehihilo_d,  Wrelohilo_d,  Wrehilolo_d,  Wrelololo_d,
               Wimhihihi_d,  Wimlohihi_d,  Wimhilohi_d,  Wimlolohi_d,
               Wimhihilo_d,  Wimlohilo_d,  Wimhilolo_d,  Wimlololo_d,
             YWTrehihihi_d,YWTrelohihi_d,YWTrehilohi_d,YWTrelolohi_d,
             YWTrehihilo_d,YWTrelohilo_d,YWTrehilolo_d,YWTrelololo_d,
             YWTimhihihi_d,YWTimlohihi_d,YWTimhilohi_d,YWTimlolohi_d,
             YWTimhihilo_d,YWTimlohilo_d,YWTimhilolo_d,YWTimlololo_d,
             YWTrehihihi_h,YWTrelohihi_h,YWTrehilohi_h,YWTrelolohi_h,
             YWTrehihilo_h,YWTrelohilo_h,YWTrehilolo_h,YWTrelololo_h,
             YWTimhihihi_h,YWTimlohihi_h,YWTimhilohi_h,YWTimlolohi_h,
             YWTimhihilo_h,YWTimlohilo_h,YWTimhilolo_h,YWTimlololo_h,
             YWHlapms,addcnt,mulcnt,verbose);

         GPU_cmplx8_small_YWHC
            (nrows,ncols,szt,k,
              YWTrehihihi_d, YWTrelohihi_d, YWTrehilohi_d, YWTrelolohi_d,
              YWTrehihilo_d, YWTrelohilo_d, YWTrehilolo_d, YWTrelololo_d,
              YWTimhihihi_d, YWTimlohihi_d, YWTimhilohi_d, YWTimlolohi_d,
              YWTimhihilo_d, YWTimlohilo_d, YWTimhilolo_d, YWTimlololo_d,
                Arehihihi_d,   Arelohihi_d,   Arehilohi_d,   Arelolohi_d,
                Arehihilo_d,   Arelohilo_d,   Arehilolo_d,   Arelololo_d,
                Aimhihihi_d,   Aimlohihi_d,   Aimhilohi_d,   Aimlolohi_d,
                Aimhihilo_d,   Aimlohilo_d,   Aimhilolo_d,   Aimlololo_d,
             YWTCrehihihi_d,YWTCrelohihi_d,YWTCrehilohi_d,YWTCrelolohi_d,
             YWTCrehihilo_d,YWTCrelohilo_d,YWTCrehilolo_d,YWTCrelololo_d,
             YWTCimhihihi_d,YWTCimlohihi_d,YWTCimhilohi_d,YWTCimlolohi_d,
             YWTCimhihilo_d,YWTCimlohilo_d,YWTCimhilolo_d,YWTCimlololo_d,
             YWTCrehihihi_h,YWTCrelohihi_h,YWTCrehilohi_h,YWTCrelolohi_h,
             YWTCrehihilo_h,YWTCrelohilo_h,YWTCrehilolo_h,YWTCrelololo_h,
             YWTCimhihihi_h,YWTCimlohihi_h,YWTCimhilohi_h,YWTCimlolohi_h,
             YWTCimhihilo_h,YWTCimlohilo_h,YWTCimhilolo_h,YWTCimlololo_h,
             YWHClapms,addcnt,mulcnt,verbose);

         GPU_cmplx8_small_R_add_YWHC
            (nrows,ncols,szt,k,
                Arehihihi_d,   Arelohihi_d,   Arehilohi_d,   Arelolohi_d,
                Arehihilo_d,   Arelohilo_d,   Arehilolo_d,   Arelololo_d,
                Aimhihihi_d,   Aimlohihi_d,   Aimhilohi_d,   Aimlolohi_d,
                Aimhihilo_d,   Aimlohilo_d,   Aimhilolo_d,   Aimlololo_d,
             YWTCrehihihi_d,YWTCrelohihi_d,YWTCrehilohi_d,YWTCrelolohi_d,
             YWTCrehihilo_d,YWTCrelohilo_d,YWTCrehilolo_d,YWTCrelololo_d,
             YWTCimhihihi_d,YWTCimlohihi_d,YWTCimhilohi_d,YWTCimlolohi_d,
             YWTCimhihilo_d,YWTCimlohilo_d,YWTCimhilolo_d,YWTCimlololo_d,
                Arehihihi_h,   Arelohihi_h,   Arehilohi_h,   Arelolohi_h,
                Arehihilo_h,   Arelohilo_h,   Arehilolo_h,   Arelololo_h,
                Aimhihihi_h,   Aimlohihi_h,   Aimhilohi_h,   Aimlolohi_h,
                Aimhihilo_h,   Aimlohilo_h,   Aimhilolo_h,   Aimlololo_h,
             Raddlapms,addcnt,verbose);
      }
   }
   gettimeofday(&endtime,0);
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   *walltimesec = seconds + microseconds*1.0e-6;

   cudaMemcpy(Qrehihihi_h,Qrehihihi_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qrelohihi_h,Qrelohihi_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qrehilohi_h,Qrehilohi_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qrelolohi_h,Qrelolohi_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qrehihilo_h,Qrehihilo_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qrelohilo_h,Qrelohilo_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qrehilolo_h,Qrehilolo_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qrelololo_h,Qrelololo_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qimhihihi_h,Qimhihihi_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qimlohihi_h,Qimlohihi_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qimhilohi_h,Qimhilohi_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qimlolohi_h,Qimlolohi_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qimhihilo_h,Qimhihilo_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qimlohilo_h,Qimlohilo_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qimhilolo_h,Qimhilolo_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qimlololo_h,Qimlololo_d,szWYT,cudaMemcpyDeviceToHost);
   ix = 0;                                           // copy rows of Q
   for(int i=0; i<nrows; i++)
      for(int j=0; j<nrows; j++)
      {
         Qrehihihi[i][j] = Qrehihihi_h[ix];
         Qrelohihi[i][j] = Qrelohihi_h[ix];
         Qrehilohi[i][j] = Qrehilohi_h[ix];
         Qrelolohi[i][j] = Qrelolohi_h[ix];
         Qrehihilo[i][j] = Qrehihilo_h[ix];
         Qrelohilo[i][j] = Qrelohilo_h[ix];
         Qrehilolo[i][j] = Qrehilolo_h[ix];
         Qrelololo[i][j] = Qrelololo_h[ix];
         Qimhihihi[i][j] = Qimhihihi_h[ix];
         Qimlohihi[i][j] = Qimlohihi_h[ix];
         Qimhilohi[i][j] = Qimhilohi_h[ix];
         Qimlolohi[i][j] = Qimlolohi_h[ix];
         Qimhihilo[i][j] = Qimhihilo_h[ix];
         Qimlohilo[i][j] = Qimlohilo_h[ix];
         Qimhilolo[i][j] = Qimhilolo_h[ix];
         Qimlololo[i][j] = Qimlololo_h[ix++];
      }

   cudaMemcpy(Arehihihi_h,Arehihihi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Arelohihi_h,Arelohihi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Arehilohi_h,Arehilohi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Arelolohi_h,Arelolohi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Arehihilo_h,Arehihilo_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Arelohilo_h,Arelohilo_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Arehilolo_h,Arehilolo_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Arelololo_h,Arelololo_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Aimhihihi_h,Aimhihihi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Aimlohihi_h,Aimlohihi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Aimhilohi_h,Aimhilohi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Aimlolohi_h,Aimlolohi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Aimhihilo_h,Aimhihilo_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Aimlohilo_h,Aimlohilo_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Aimhilolo_h,Aimhilolo_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Aimlololo_h,Aimlololo_d,sznum,cudaMemcpyDeviceToHost);

   for(int i=0; i<nrows; i++)                       // copy columns of R
      for(int j=0; j<ncols; j++)
      {
         Rrehihihi[i][j] = Arehihihi_h[j*nrows+i];
         Rrelohihi[i][j] = Arelohihi_h[j*nrows+i];
         Rrehilohi[i][j] = Arehilohi_h[j*nrows+i];
         Rrelolohi[i][j] = Arelolohi_h[j*nrows+i];
         Rrehihilo[i][j] = Arehihilo_h[j*nrows+i];
         Rrelohilo[i][j] = Arelohilo_h[j*nrows+i];
         Rrehilolo[i][j] = Arehilolo_h[j*nrows+i];
         Rrelololo[i][j] = Arelololo_h[j*nrows+i];
         Rimhihihi[i][j] = Aimhihihi_h[j*nrows+i];
         Rimlohihi[i][j] = Aimlohihi_h[j*nrows+i];
         Rimhilohi[i][j] = Aimhilohi_h[j*nrows+i];
         Rimlolohi[i][j] = Aimlolohi_h[j*nrows+i];
         Rimhihilo[i][j] = Aimhihilo_h[j*nrows+i];
         Rimlohilo[i][j] = Aimlohilo_h[j*nrows+i];
         Rimhilolo[i][j] = Aimhilolo_h[j*nrows+i];
         Rimlololo[i][j] = Aimlololo_h[j*nrows+i];
      }

   free(Arehihihi_h); free(Arelohihi_h); free(Arehilohi_h); free(Arelolohi_h);
   free(Arehihilo_h); free(Arelohilo_h); free(Arehilolo_h); free(Arelololo_h);
   free(Aimhihihi_h); free(Aimlohihi_h); free(Aimhilohi_h); free(Aimlolohi_h);
   free(Aimhihilo_h); free(Aimlohilo_h); free(Aimhilolo_h); free(Aimlololo_h);
   free(Qrehihihi_h); free(Qrelohihi_h); free(Qrehilohi_h); free(Qrelolohi_h);
   free(Qrehihilo_h); free(Qrelohilo_h); free(Qrehilolo_h); free(Qrelololo_h);
   free(Qimhihihi_h); free(Qimlohihi_h); free(Qimhilohi_h); free(Qimlolohi_h);
   free(Qimhihilo_h); free(Qimlohilo_h); free(Qimhilolo_h); free(Qimlololo_h);
   free(vrehihihi_h); free(vrelohihi_h); free(vrehilohi_h); free(vrelolohi_h); 
   free(vrehihilo_h); free(vrelohilo_h); free(vrehilolo_h); free(vrelololo_h); 
   free(vimhihihi_h); free(vimlohihi_h); free(vimhilohi_h); free(vimlolohi_h);
   free(vimhihilo_h); free(vimlohilo_h); free(vimhilolo_h); free(vimlololo_h);
   free(Vrehihihi_h); free(Vrelohihi_h); free(Vrehilohi_h); free(Vrelolohi_h);
   free(Vrehihilo_h); free(Vrelohilo_h); free(Vrehilolo_h); free(Vrelololo_h);
   free(Vimhihihi_h); free(Vimlohihi_h); free(Vimhilohi_h); free(Vimlolohi_h);
   free(Vimhihilo_h); free(Vimlohilo_h); free(Vimhilolo_h); free(Vimlololo_h);
   free(RHdotvrehihihi_h); free(RHdotvrelohihi_h);
   free(RHdotvrehilohi_h); free(RHdotvrelolohi_h);
   free(RHdotvrehihilo_h); free(RHdotvrelohilo_h);
   free(RHdotvrehilolo_h); free(RHdotvrelololo_h);
   free(RHdotvimhihihi_h); free(RHdotvimlohihi_h);
   free(RHdotvimhilohi_h); free(RHdotvimlolohi_h);
   free(RHdotvimhihilo_h); free(RHdotvimlohilo_h);
   free(RHdotvimhilolo_h); free(RHdotvimlololo_h);
   free(bRHvrehihihi_h); free(bRHvrelohihi_h);
   free(bRHvrehilohi_h); free(bRHvrelolohi_h);
   free(bRHvrehihilo_h); free(bRHvrelohilo_h);
   free(bRHvrehilolo_h); free(bRHvrelololo_h);
   free(bRHvimhihihi_h); free(bRHvimlohihi_h);
   free(bRHvimhilohi_h); free(bRHvimlolohi_h);
   free(bRHvimhihilo_h); free(bRHvimlohilo_h);
   free(bRHvimhilolo_h); free(bRHvimlololo_h);
   free(sumshihihi_h); free(sumslohihi_h);
   free(sumshilohi_h); free(sumslolohi_h);
   free(sumshihilo_h); free(sumslohilo_h);
   free(sumshilolo_h); free(sumslololo_h);
   free(Wrehihihi_h); free(Wrelohihi_h); free(Wrehilohi_h); free(Wrelolohi_h);
   free(Wrehihilo_h); free(Wrelohilo_h); free(Wrehilolo_h); free(Wrelololo_h);
   free(Wimhihihi_h); free(Wimlohihi_h); free(Wimhilohi_h); free(Wimlolohi_h);
   free(Wimhihilo_h); free(Wimlohilo_h); free(Wimhilolo_h); free(Wimlololo_h);
   free(WYTrehihihi_h); free(WYTrelohihi_h);
   free(WYTrehilohi_h); free(WYTrelolohi_h);
   free(WYTrehihilo_h); free(WYTrelohilo_h);
   free(WYTrehilolo_h); free(WYTrelololo_h);
   free(WYTimhihihi_h); free(WYTimlohihi_h);
   free(WYTimhilohi_h); free(WYTimlolohi_h);
   free(WYTimhihilo_h); free(WYTimlohilo_h);
   free(WYTimhilolo_h); free(WYTimlololo_h);
   free(QWYTrehihihi_h); free(QWYTrelohihi_h);
   free(QWYTrehilohi_h); free(QWYTrelolohi_h);
   free(QWYTrehihilo_h); free(QWYTrelohilo_h);
   free(QWYTrehilolo_h); free(QWYTrelololo_h);
   free(QWYTimhihihi_h); free(QWYTimlohihi_h);
   free(QWYTimhilohi_h); free(QWYTimlolohi_h);
   free(QWYTimhihilo_h); free(QWYTimlohilo_h);
   free(QWYTimhilolo_h); free(QWYTimlololo_h);
   free(YWTrehihihi_h); free(YWTrelohihi_h);
   free(YWTrehilohi_h); free(YWTrelolohi_h);
   free(YWTrehihilo_h); free(YWTrelohilo_h);
   free(YWTrehilolo_h); free(YWTrelololo_h);
   free(YWTimhihihi_h); free(YWTimlohihi_h);
   free(YWTimhilohi_h); free(YWTimlolohi_h);
   free(YWTimhihilo_h); free(YWTimlohilo_h);
   free(YWTimhilolo_h); free(YWTimlololo_h);
   free(YWTCrehihihi_h); free(YWTCrelohihi_h);
   free(YWTCrehilohi_h); free(YWTCrelolohi_h);
   free(YWTCrehihilo_h); free(YWTCrelohilo_h);
   free(YWTCrehilolo_h); free(YWTCrelololo_h);
   free(YWTCimhihihi_h); free(YWTCimlohihi_h);
   free(YWTCimhilohi_h); free(YWTCimlolohi_h);
   free(YWTCimhihilo_h); free(YWTCimlohilo_h);
   free(YWTCimhilolo_h); free(YWTCimlololo_h);
}
