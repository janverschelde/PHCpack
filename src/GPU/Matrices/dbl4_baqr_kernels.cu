/* The file dbl4_baqr_kernels.cu defines the functions with prototypes in
 * the file dbl4_baqr_kernels.h. */

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
#endif
#include "dbl4_baqr_kernels.h"
#include "quad_double_functions.h"
#include "dbl_baqr_flopcounts.h"

using namespace std;

__global__ void dbl4_small_house
 ( double *x0hihi, double *x0lohi, double *x0hilo, double *x0lolo,
   double *x1hihi, double *x1lohi, double *x1hilo, double *x1lolo,
   int dim, int dimLog2,
   double *vhihi, double *vlohi, double *vhilo, double *vlolo,
   double *betahihi, double *betalohi, double *betahilo, double *betalolo )
{
   const int j = threadIdx.x;

   __shared__ double shvhihi[qd_shmemsize];
   __shared__ double shvlohi[qd_shmemsize];
   __shared__ double shvhilo[qd_shmemsize];
   __shared__ double shvlolo[qd_shmemsize];
   __shared__ double prdhihi[qd_shmemsize];
   __shared__ double prdlohi[qd_shmemsize];
   __shared__ double prdhilo[qd_shmemsize];
   __shared__ double prdlolo[qd_shmemsize];

   bool stopflag = false;
   double acchihi,acclohi,acchilo,acclolo;
   double muhihi,mulohi,muhilo,mulolo;
   double v0hihi,v0lohi,v0hilo,v0lolo;
   double v0p2hihi,v0p2lohi,v0p2hilo,v0p2lolo;

   shvhihi[j] = x1hihi[j];          // reading of vector into shared memory
   shvlohi[j] = x1lohi[j];
   shvhilo[j] = x1hilo[j];
   shvlolo[j] = x1lolo[j];
   // prd[j] = shv[j]*shv[j];   // for the 2-norm computation
   __syncthreads();
   qdg_sqr( shvhihi[j], shvlohi[j], shvhilo[j], shvlolo[j],
           &prdhihi[j],&prdlohi[j],&prdhilo[j],&prdlolo[j]);

   __syncthreads();
   vhihi[j+1] = shvhihi[j];         // copies x to v, in case beta is zero
   vlohi[j+1] = shvlohi[j];
   vhilo[j+1] = shvhilo[j];
   vlolo[j+1] = shvlolo[j];
   __syncthreads();
   if(j == 0)
   {
      vhihi[0] = 1.0;
      vlohi[0] = 0.0;
      vhilo[0] = 0.0;
      vlolo[0] = 0.0;
   }
   __syncthreads();
   int powTwo = 1;                          // sum reduction
   for(int k=0; k < dimLog2; k++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < dim) // prd[j] = prd[j] + prd[j+powTwo];
            qdg_inc(&prdhihi[j],&prdlohi[j],&prdhilo[j],&prdlolo[j],
                     prdhihi[j+powTwo],prdlohi[j+powTwo],
                     prdhilo[j+powTwo],prdlolo[j+powTwo]);
      powTwo = powTwo*2;
      __syncthreads();
   }
   // thread 0 computes the sqrt of the inner product, others wait
   if(j == 0)
   {                                       // prd[0] is sigma of house
      if((prdhihi[0] == 0.0) && (prdlohi[0] == 0.0) &&
         (prdhilo[0] == 0.0) && (prdlolo[0] == 0.0)) 
      {
         *betahihi = 0.0; *betalohi = 0.0;
         *betahilo = 0.0; *betalolo = 0.0; stopflag = true;
      }
   }
   __syncthreads();
   if(stopflag) return;                    // case when sigma is zero
   __syncthreads();
   if(j == 0)                              // thread zero sets beta
   {
      // mu = sqrt((*x0)*(*x0) + prd[0]);
      qdg_sqr( *x0hihi, *x0lohi, *x0hilo, *x0lolo,
              &acchihi,&acclohi,&acchilo,&acclolo);
      qdg_inc(&acchihi,  &acclohi,  &acchilo,  &acclolo,
               prdhihi[0],prdlohi[0],prdhilo[0],prdlolo[0]);
      qdg_sqrt(acchihi,acclohi,acchilo,acclolo,
               &muhihi,&mulohi,&muhilo,&mulolo);

      if(*x0hihi <= 0.0)
      {
         // v0 = *x0 - mu;
         qdg_sub(*x0hihi,*x0lohi,*x0hilo,*x0lolo,
                  muhihi, mulohi, muhilo, mulolo,
                 &v0hihi,&v0lohi,&v0hilo,&v0lolo);
      }
      else
      {
         // v0 = -prd[0]/(*x0 + mu);
         qdg_add( *x0hihi, *x0lohi, *x0hilo, *x0lolo,
                   muhihi,  mulohi,  muhilo,  mulolo,
                 &acchihi,&acclohi,&acchilo,&acclolo);
         qdg_div(prdhihi[0],prdlohi[0],prdhilo[0],prdlolo[0],
                 acchihi,   acclohi,   acchilo,   acclolo,
                 &v0hihi,   &v0lohi,   &v0hilo,   &v0lolo);
         qdg_minus(&v0hihi,&v0lohi,&v0hilo,&v0lolo);
      }
      // v0p2 = v0*v0;
      qdg_sqr(   v0hihi,   v0lohi,   v0hilo,   v0lolo,
              &v0p2hihi,&v0p2lohi,&v0p2hilo,&v0p2lolo);
      // *beta = 2.0*v0p2/(prd[0] + v0p2);
      qdg_add( prdhihi[0],prdlohi[0],prdhilo[0],prdlolo[0],
              v0p2hihi,  v0p2lohi,  v0p2hilo,  v0p2lolo,
              &acchihi,  &acclohi,  &acchilo,  &acclolo);
      qdg_div(v0p2hihi,v0p2lohi,v0p2hilo,v0p2lolo,
               acchihi, acclohi, acchilo, acclolo,
              betahihi,betalohi,betahilo,betalolo);
      qdg_mlt_d(betahihi,betalohi,betahilo,betalolo,2.0);
      prdhihi[0] = v0hihi;
      prdlohi[0] = v0lohi;
      prdhilo[0] = v0hilo;
      prdlolo[0] = v0lolo;                  // v0 needed for normalization
   }
   __syncthreads();
   // shv[j] = shv[j]/prd[0];
   vhihi[j+1] = 0.0;
   vlohi[j+1] = 0.0;  // initialize in case of zero beta ...
   vhilo[j+1] = 0.0;
   vlolo[j+1] = 0.0;
   if(1.0 + *betahihi + *betalohi + *betahilo + *betalolo != 1.0)
   {
      qdg_div(shvhihi[j],shvlohi[j],shvhilo[j],shvlolo[j],
              prdhihi[0],prdlohi[0],prdhilo[0],prdlolo[0],
             &acchihi,  &acclohi,  &acchilo  ,&acclolo);
      __syncthreads();
      vhihi[j+1] = acchihi;
      vlohi[j+1] = acclohi;
      vhilo[j+1] = acchilo;
      vlolo[j+1] = acclolo;
   }
   __syncthreads();
   if(j == 0)
   {
      vhihi[0] = 1.0;
      vlohi[0] = 0.0;
      vhilo[0] = 0.0;
      vlolo[0] = 0.0;
   }
}

__global__ void cmplx4_small_house
 ( double *x0rehihi, double *x0relohi, double *x0rehilo, double *x0relolo,
   double *x0imhihi, double *x0imlohi, double *x0imhilo, double *x0imlolo,
   double *x1rehihi, double *x1relohi, double *x1rehilo, double *x1relolo,
   double *x1imhihi, double *x1imlohi, double *x1imhilo, double *x1imlolo,
   int dim, int dimLog2,
   double *vrehihi, double *vrelohi, double *vrehilo, double *vrelolo,
   double *vimhihi, double *vimlohi, double *vimhilo, double *vimlolo,
   double *betahihi, double *betalohi, double *betahilo, double *betalolo )
{
   const int j = threadIdx.x;

   __shared__ double shvrehihi[cqd_shmemsize];
   __shared__ double shvrelohi[cqd_shmemsize];
   __shared__ double shvrehilo[cqd_shmemsize];
   __shared__ double shvrelolo[cqd_shmemsize];
   __shared__ double shvimhihi[cqd_shmemsize];
   __shared__ double shvimlohi[cqd_shmemsize];
   __shared__ double shvimhilo[cqd_shmemsize];
   __shared__ double shvimlolo[cqd_shmemsize];
   __shared__ double prdhihi[cqd_shmemsize];
   __shared__ double prdlohi[cqd_shmemsize];
   __shared__ double prdhilo[cqd_shmemsize];
   __shared__ double prdlolo[cqd_shmemsize];
   __shared__ double v0parts[8];

   bool stopflag = false;
   double muhihi,mulohi,muhilo,mulolo;
   double v0rehihi,v0relohi,v0rehilo,v0relolo;
   double v0imhihi,v0imlohi,v0imhilo,v0imlolo;
   double x0radhihi,x0radlohi,x0radhilo,x0radlolo;
   double sqrx0hihi,sqrx0lohi,sqrx0hilo,sqrx0lolo;
   double sqrv0hihi,sqrv0lohi,sqrv0hilo,sqrv0lolo;
   double inv0rehihi,inv0relohi,inv0rehilo,inv0relolo;
   double inv0imhihi,inv0imlohi,inv0imhilo,inv0imlolo;
   double zrehihi,zrelohi,zrehilo,zrelolo;
   double zimhihi,zimlohi,zimhilo,zimlolo;
   double acchihi,acclohi,acchilo,acclolo;

   shvrehihi[j] = x1rehihi[j];    // reading of vector into shared memory
   shvrelohi[j] = x1relohi[j];
   shvrehilo[j] = x1rehilo[j];
   shvrelolo[j] = x1relolo[j];
   shvimhihi[j] = x1imhihi[j];
   shvimlohi[j] = x1imlohi[j];
   shvimhilo[j] = x1imhilo[j];
   shvimlolo[j] = x1imlolo[j];
   // prd[j] = shv[j]*shv[j];   // for the 2-norm computation
   // prd[j] = shvre[j]*shvre[j] + shvim[j]*shvim[j];
   __syncthreads();
   qdg_sqr(shvrehihi[j],shvrelohi[j],shvrehilo[j],shvrelolo[j],
            &prdhihi[j], &prdlohi[j], &prdhilo[j], &prdlolo[j]);
   qdg_sqr(shvimhihi[j],shvimlohi[j],shvimhilo[j],shvimlolo[j],
            &acchihi,    &acclohi,    &acchilo,    &acclolo);
   qdg_inc(&prdhihi[j],&prdlohi[j],&prdhilo[j],&prdlolo[j],
            acchihi,    acclohi,    acchilo,    acclolo);
   __syncthreads();
   vrehihi[j+1] = shvrehihi[j];    // copies x to v, in case beta is zero
   vrelohi[j+1] = shvrelohi[j];
   vrehilo[j+1] = shvrehilo[j];
   vrelolo[j+1] = shvrelolo[j];
   vimhihi[j+1] = shvimhihi[j];
   vimlohi[j+1] = shvimlohi[j];
   vimhilo[j+1] = shvimhilo[j];
   vimlolo[j+1] = shvimlolo[j];
   __syncthreads();
   if(j == 0)
   {
      vrehihi[0] = 1.0; vrelohi[0] = 0.0;
      vrehilo[0] = 0.0; vrelolo[0] = 0.0;
      vimhihi[0] = 0.0; vimlohi[0] = 0.0;
      vimhilo[0] = 0.0; vimlolo[0] = 0.0;
   }
   __syncthreads();
   int powTwo = 1;                          // sum reduction
   for(int k=0; k < dimLog2; k++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < dim)     // prd[j] = prd[j] + prd[j+powTwo];
            qdg_inc(&prdhihi[j],&prdlohi[j],&prdhilo[j],&prdlolo[j],
                     prdhihi[j+powTwo],prdlohi[j+powTwo],
                     prdhilo[j+powTwo],prdlolo[j+powTwo]);
      powTwo = powTwo*2;
      __syncthreads();
   }
   // thread 0 computes the sqrt of the inner product, others wait
   if(j == 0)
   {                                        // prd[0] is sigma of house
      if((prdhihi[0] == 0.0) && (prdlohi[0] == 0.0) &&
         (prdhilo[0] == 0.0) && (prdlolo[0] == 0.0))
      {
         *betahihi = 0.0; *betalohi = 0.0;
         *betahilo = 0.0; *betalolo = 0.0; stopflag = true;
      }
   }
   __syncthreads();
   if(stopflag) return;                    // case when sigma is zero
   __syncthreads();
   if(j == 0)                              // thread zero sets beta
   {
      // sqrx0 = (*x0re)*(*x0re) + (*x0im)*(*x0im);
      qdg_sqr( *x0rehihi, *x0relohi, *x0rehilo, *x0relolo,
              &sqrx0hihi,&sqrx0lohi,&sqrx0hilo,&sqrx0lolo);
      qdg_sqr(*x0imhihi,*x0imlohi,*x0imhilo,*x0imlolo,
               &acchihi, &acclohi, &acchilo, &acclolo);
      qdg_inc(&sqrx0hihi,&sqrx0lohi,&sqrx0hilo,&sqrx0lolo,
                 acchihi,   acclohi,   acchilo,   acclolo);
      // x0rad = sqrt(sqrx0);
      qdg_sqrt( sqrx0hihi, sqrx0lohi, sqrx0hilo, sqrx0lolo,
               &x0radhihi,&x0radlohi,&x0radhilo,&x0radlolo);
      // mu = sqrt(sqrx0 + prd[0]);
      qdg_inc(&sqrx0hihi,&sqrx0lohi,&sqrx0hilo,&sqrx0lolo,
                 prdhihi[0],prdlohi[0],prdhilo[0],prdlolo[0]);
      qdg_sqrt(sqrx0hihi,sqrx0lohi,sqrx0hilo,sqrx0lolo,
                 &muhihi,  &mulohi,  &muhilo,  &mulolo);

      if((x0radhihi == 0.0) && (x0radlohi == 0.0) &&
         (x0radhilo == 0.0) && (x0radlolo == 0.0))
      {
         v0rehihi = muhihi; v0relohi = mulohi;
         v0rehilo = muhilo; v0relolo = mulolo;

         qdg_minus(&v0rehihi,&v0relohi,&v0rehilo,&v0relolo);

         v0imhihi = 0.0; v0imlohi = 0.0;
         v0imhilo = 0.0; v0imlolo = 0.0;
      }
      else
      {
         // mu = mu/x0rad;
         qdg_div(   muhihi,   mulohi,   muhilo,   mulolo,
                 x0radhihi,x0radlohi,x0radhilo,x0radlolo,
                  &acchihi, &acclohi, &acchilo, &acclolo);
         muhihi = acchihi; mulohi = acclohi;
         muhilo = acchilo; mulolo = acclolo;
         // v0re = (*x0re) - mu*(*x0re);
         qdg_mul(  muhihi,     mulohi,     muhilo,     mulolo,
                 x0rehihi[0],x0relohi[0],x0rehilo[0],x0relolo[0],
                 &acchihi,   &acclohi,   &acchilo,   &acclolo);
         qdg_sub(x0rehihi[0],x0relohi[0],x0rehilo[0],x0relolo[0],
                  acchihi,    acclohi,    acchilo,    acclolo,
                &v0rehihi,  &v0relohi,  &v0rehilo,  &v0relolo);
         // v0im = (*x0im) - mu*(*x0im);
         qdg_mul(  muhihi,     mulohi,     muhilo,     mulolo,
                 x0imhihi[0],x0imlohi[0],x0imhilo[0],x0imlolo[0],
                 &acchihi,   &acclohi,   &acchilo,   &acclolo);
         qdg_sub(x0imhihi[0],x0imlohi[0],x0imhilo[0],x0imlolo[0],
                  acchihi,    acclohi,    acchilo,    acclolo,
                &v0imhihi,  &v0imlohi,  &v0imhilo,  &v0imlolo);
      }
      // sqrv0 = v0re*v0re + v0im*v0im;
      qdg_sqr(  v0rehihi,  v0relohi,  v0rehilo,  v0relolo,
              &sqrv0hihi,&sqrv0lohi,&sqrv0hilo,&sqrv0lolo);
      qdg_sqr(v0imhihi,v0imlohi,v0imhilo,v0imlolo,
              &acchihi,&acclohi,&acchilo,&acclolo);
      qdg_inc(&sqrv0hihi,&sqrv0lohi,&sqrv0hilo,&sqrv0lolo,
                 acchihi,   acclohi,   acchilo,   acclolo);
      // *beta = 2.0*sqrv0/(prd[0] + sqrv0);
      qdg_add(  prdhihi[0],prdlohi[0],prdhilo[0],prdlolo[0],
              sqrv0hihi, sqrv0lohi, sqrv0hilo, sqrv0lolo,
               &acchihi,  &acclohi,  &acchilo,  &acclolo);
      qdg_div(sqrv0hihi,sqrv0lohi,sqrv0hilo,sqrv0lolo,
                acchihi,  acclohi,  acchilo,  acclolo,
               betahihi, betalohi, betahilo, betalolo);
      qdg_mlt_d(betahihi,betalohi,betahilo,betalolo,2.0);

      prdhihi[0] = sqrv0hihi;             // sqrv0 needed for normalization
      prdlohi[0] = sqrv0lohi;
      prdhilo[0] = sqrv0hilo;
      prdlolo[0] = sqrv0lolo;
      v0parts[0] = v0rehihi;              // share v0re with all threads
      v0parts[1] = v0relohi; 
      v0parts[2] = v0rehilo; 
      v0parts[3] = v0relolo; 
      v0parts[4] = v0imhihi;              // share v0im with all threads
      v0parts[5] = v0imlohi;
      v0parts[6] = v0imhilo;
      v0parts[7] = v0imlolo;
   }
   __syncthreads(); // important synchronization!
   // inv0re = v0parts[0]/prd[0];               // real part of 1/v[0]
   vrehihi[j+1] = 0.0;
   vrelohi[j+1] = 0.0;
   vrehilo[j+1] = 0.0;
   vrelolo[j+1] = 0.0;
   vimhihi[j+1] = 0.0;
   vimlohi[j+1] = 0.0;
   vimhilo[j+1] = 0.0;
   vimlolo[j+1] = 0.0;
   if(1.0 + prdhihi[0] + prdlohi[0] + prdhilo[0] + prdlolo[0] != 1.0)
   {
      qdg_div(v0parts[0],v0parts[1],v0parts[2],v0parts[3],
                  prdhihi[0], prdlohi[0], prdhilo[0], prdlolo[0],
              &inv0rehihi,&inv0relohi,&inv0rehilo,&inv0relolo);
      // inv0im = -v0parts[1]/prd[0];              // imag part of 1/v[0]
      qdg_div(v0parts[4],v0parts[5],v0parts[6],v0parts[7],
                  prdhihi[0], prdlohi[0], prdhilo[0], prdlolo[0],
              &inv0imhihi,&inv0imlohi,&inv0imhilo,&inv0imlolo);
      qdg_minus(&inv0imhihi,&inv0imlohi,&inv0imhilo,&inv0imlolo);
      // zre = shvre[j]*inv0re - shvim[j]*inv0im;  // real part of v[j]/v[0]
      __syncthreads();
      qdg_mul( shvrehihi[j],shvrelohi[j],shvrehilo[j],shvrelolo[j],
              inv0rehihi,  inv0relohi,  inv0rehilo,  inv0relolo,
                &zrehihi,    &zrelohi,    &zrehilo,    &zrelolo);
      qdg_mul( shvimhihi[j],shvimlohi[j],shvimhilo[j],shvimlolo[j],
              inv0imhihi,  inv0imlohi,  inv0imhilo,  inv0imlolo,
                &acchihi,    &acclohi,    &acchilo,    &acclolo);
      qdg_dec(&zrehihi,&zrelohi,&zrehilo,&zrelolo,
               acchihi, acclohi, acchilo, acclolo);
      // zim = shvim[j]*inv0re + shvre[j]*inv0im;  // imag part of v[j]/v[0]
      qdg_mul( shvimhihi[j],shvimlohi[j],shvimhilo[j],shvimlolo[j],
              inv0rehihi,  inv0relohi,  inv0rehilo,  inv0relolo,
                &zimhihi,    &zimlohi,    &zimhilo,    &zimlolo);
      qdg_mul( shvrehihi[j],shvrelohi[j],shvrehilo[j],shvrelolo[j],
              inv0imhihi,  inv0imlohi,  inv0imhilo,  inv0imlolo,
                &acchihi,    &acclohi,    &acchilo,    &acclolo);
      qdg_inc(&zimhihi,&zimlohi,&zimhilo,&zimlolo,
               acchihi, acclohi, acchilo, acclolo);
      __syncthreads();
      vrehihi[j+1] = zrehihi;
      vrelohi[j+1] = zrelohi;
      vrehilo[j+1] = zrehilo;
      vrelolo[j+1] = zrelolo;
      vimhihi[j+1] = zimhihi;
      vimlohi[j+1] = zimlohi;
      vimhilo[j+1] = zimhilo;
      vimlolo[j+1] = zimlolo;
   }
   __syncthreads();
   if(j == 0)
   {
      vrehihi[0] = 1.0; vrelohi[0] = 0.0;
      vrehilo[0] = 0.0; vrelolo[0] = 0.0;
      vimhihi[0] = 0.0; vimlohi[0] = 0.0;
      vimhilo[0] = 0.0; vimlolo[0] = 0.0;
   }
}

__global__ void dbl4_large_sum_of_squares
 ( double *vhihi, double *vlohi, double *vhilo, double *vlolo,
   double *sumshihi, double *sumslohi, double *sumshilo, double *sumslolo,
   int dim, int BS, int BSLog2 )
{
   const int i = blockIdx.x;
   const int j = threadIdx.x;
   const int k = i*BS + j;

   __shared__ double shvhihi[inner_qd_shmemsize];
   __shared__ double shvlohi[inner_qd_shmemsize];
   __shared__ double shvhilo[inner_qd_shmemsize];
   __shared__ double shvlolo[inner_qd_shmemsize];
   __shared__ double prdhihi[inner_qd_shmemsize];
   __shared__ double prdlohi[inner_qd_shmemsize];
   __shared__ double prdhilo[inner_qd_shmemsize];
   __shared__ double prdlolo[inner_qd_shmemsize];

   shvhihi[j] = vhihi[k];
   shvlohi[j] = vlohi[k];
   shvhilo[j] = vhilo[k];
   shvlolo[j] = vlolo[k];
   __syncthreads();
   if(k >= dim)
   {
      shvhihi[j] = 0.0;
      shvlohi[j] = 0.0;
      shvhilo[j] = 0.0;
      shvlolo[j] = 0.0;
   }
   __syncthreads();
   qdg_sqr( shvhihi[j], shvlohi[j], shvhilo[j], shvlolo[j],
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

__global__ void cmplx4_large_sum_of_squares
 ( double *vrehihi, double *vrelohi, double *vrehilo, double *vrelolo,
   double *vimhihi, double *vimlohi, double *vimhilo, double *vimlolo,
   double *sumshihi, double *sumslohi, double *sumshilo, double *sumslolo,
   int dim, int BS, int BSLog2 )
{
   const int i = blockIdx.x;
   const int j = threadIdx.x;
   const int k = i*BS + j;

   __shared__ double shvrehihi[inner_qd_shmemsize];
   __shared__ double shvrelohi[inner_qd_shmemsize];
   __shared__ double shvrehilo[inner_qd_shmemsize];
   __shared__ double shvrelolo[inner_qd_shmemsize];
   __shared__ double shvimhihi[inner_qd_shmemsize];
   __shared__ double shvimlohi[inner_qd_shmemsize];
   __shared__ double shvimhilo[inner_qd_shmemsize];
   __shared__ double shvimlolo[inner_qd_shmemsize];
   __shared__ double prdhihi[inner_qd_shmemsize];
   __shared__ double prdlohi[inner_qd_shmemsize];
   __shared__ double prdhilo[inner_qd_shmemsize];
   __shared__ double prdlolo[inner_qd_shmemsize];

   shvrehihi[j] = vrehihi[k];
   shvrelohi[j] = vrelohi[k];
   shvrehilo[j] = vrehilo[k];
   shvrelolo[j] = vrelolo[k];
   shvimhihi[j] = vimhihi[k];
   shvimlohi[j] = vimlohi[k];
   shvimhilo[j] = vimhilo[k];
   shvimlolo[j] = vimlolo[k];
   __syncthreads();
   if(k >= dim)
   {
      shvrehihi[j] = 0.0;
      shvrelohi[j] = 0.0;
      shvrehilo[j] = 0.0;
      shvrelolo[j] = 0.0;
      shvimhihi[j] = 0.0;
      shvimlohi[j] = 0.0;
      shvimhilo[j] = 0.0;
      shvimlolo[j] = 0.0;
   }
   double acchihi,acclohi,acchilo,acclolo;
   
   __syncthreads();
   qdg_sqr(shvrehihi[j],shvrelohi[j],shvrehilo[j],shvrelolo[j],
            &prdhihi[j], &prdlohi[j], &prdhilo[j], &prdlolo[j]);
   qdg_sqr(shvimhihi[j],shvimlohi[j],shvimhilo[j],shvimlolo[j],
            &acchihi,    &acclohi,    &acchilo,    &acclolo);
   qdg_inc(&prdhihi[j],&prdlohi[j],&prdhilo[j],&prdlolo[j],
            acchihi,    acclohi,    acchilo,    acclolo);

   __syncthreads();

   int powTwo = 1;                          // sum reduction
   for(int k=0; k < BSLog2; k++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < BS)     // prd[j] = prd[j] + prd[j+powTwo];
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

__global__ void dbl4_sum_accumulator
 ( double *sumshihi, double *sumslohi, double *sumshilo, double *sumslolo,
   int nbsums, int nbsumsLog2,
   double *acchihi, double *acclohi, double *acchilo, double *acclolo )
{
   const int j = threadIdx.x;

   __shared__ double shvhihi[outer_qd_shmemsize];
   __shared__ double shvlohi[outer_qd_shmemsize];
   __shared__ double shvhilo[outer_qd_shmemsize];
   __shared__ double shvlolo[outer_qd_shmemsize];

   shvhihi[j] = sumshihi[j];
   shvlohi[j] = sumslohi[j];
   shvhilo[j] = sumshilo[j];
   shvlolo[j] = sumslolo[j];

   __syncthreads();

   if(j >= nbsums)
   {
      shvhihi[j] = 0.0;
      shvlohi[j] = 0.0;
      shvhilo[j] = 0.0;
      shvlolo[j] = 0.0;
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
   __syncthreads();
   if(j == 0)
   {
      *acchihi = shvhihi[0];
      *acclohi = shvlohi[0];
      *acchilo = shvhilo[0];
      *acclolo = shvlolo[0];
   }
}

__global__ void dbl4_normalize
 ( int dim, int szt,
   double *xhihi, double *xlohi, double *xhilo, double *xlolo,
   double *v0hihi, double *v0lohi, double *v0hilo, double *v0lolo,
   double *vhihi, double *vlohi, double *vhilo, double *vlolo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int idx = bdx*szt + tdx;  // thread tdx scales idx

   __shared__ double shvhihi[inner_qd_shmemsize];
   __shared__ double shvlohi[inner_qd_shmemsize];
   __shared__ double shvhilo[inner_qd_shmemsize];
   __shared__ double shvlolo[inner_qd_shmemsize];

   shvhihi[tdx] = xhihi[idx];
   shvlohi[tdx] = xlohi[idx];
   shvhilo[tdx] = xhilo[idx];
   shvlolo[tdx] = xlolo[idx];
   __syncthreads();

   double resulthihi,resultlohi,resulthilo,resultlolo;

   // shv[j] = shv[j]/v0;
   qdg_div(    shvhihi[tdx],shvlohi[tdx],shvhilo[tdx],shvlolo[tdx],
                v0hihi[0],   v0lohi[0],   v0hilo[0],   v0lolo[0],
           &resulthihi, &resultlohi, &resulthilo, &resultlolo);

   __syncthreads();
   if(idx < dim)
   {
      vhihi[idx] = resulthihi;
      vlohi[idx] = resultlohi;
      vhilo[idx] = resulthilo;
      vlolo[idx] = resultlolo;
   }
}

__global__ void cmplx4_normalize
 ( int dim, int szt,
   double *xrehihi, double *xrelohi, double *xrehilo, double *xrelolo,
   double *ximhihi, double *ximlohi, double *ximhilo, double *ximlolo,
   double *inv0rehihi, double *inv0relohi,
   double *inv0rehilo, double *inv0relolo,
   double *inv0imhihi, double *inv0imlohi,
   double *inv0imhilo, double *inv0imlolo,
   double *vrehihi, double *vrelohi, double *vrehilo, double *vrelolo,
   double *vimhihi, double *vimlohi, double *vimhilo, double *vimlolo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int idx = bdx*szt + tdx;  // thread tdx scales idx

   const double invrehihi = *inv0rehihi;
   const double invrelohi = *inv0relohi;
   const double invrehilo = *inv0rehilo;
   const double invrelolo = *inv0relolo;
   const double invimhihi = *inv0imhihi;
   const double invimlohi = *inv0imlohi;
   const double invimhilo = *inv0imhilo;
   const double invimlolo = *inv0imlolo;

   __shared__ double shvrehihi[inner_qd_shmemsize];
   __shared__ double shvrelohi[inner_qd_shmemsize];
   __shared__ double shvrehilo[inner_qd_shmemsize];
   __shared__ double shvrelolo[inner_qd_shmemsize];
   __shared__ double shvimhihi[inner_qd_shmemsize];
   __shared__ double shvimlohi[inner_qd_shmemsize];
   __shared__ double shvimhilo[inner_qd_shmemsize];
   __shared__ double shvimlolo[inner_qd_shmemsize];

   shvrehihi[tdx] = xrehihi[idx];
   shvrelohi[tdx] = xrelohi[idx];
   shvrehilo[tdx] = xrehilo[idx];
   shvrelolo[tdx] = xrelolo[idx];
   shvimhihi[tdx] = ximhihi[idx];
   shvimlohi[tdx] = ximlohi[idx];
   shvimhilo[tdx] = ximhilo[idx];
   shvimlolo[tdx] = ximlolo[idx];
   __syncthreads();

   double resulthihi,resultlohi,resulthilo,resultlolo;
   double acchihi,acclohi,acchilo,acclolo;

   // shv[j] = shv[j]/v0;

   // resultre = vre[i]*inv0re - vim[i]*inv0im;
   qdg_mul(  shvrehihi[tdx],shvrelohi[tdx],shvrehilo[tdx],shvrelolo[tdx],
             invrehihi,     invrelohi,     invrehilo,     invrelolo,
           &resulthihi,   &resultlohi,   &resulthilo,   &resultlolo);
   qdg_mul(  shvimhihi[tdx],shvimlohi[tdx],shvimhilo[tdx],shvimlolo[tdx],
             invimhihi,     invimlohi,     invimhilo,     invimlolo,
              &acchihi,      &acclohi,      &acchilo,      &acclolo);
   qdg_dec(&resulthihi,&resultlohi,&resulthilo,&resultlolo,
               acchihi,    acclohi,    acchilo,    acclolo);
   if(idx < dim)
   {
      vrehihi[idx] = resulthihi;
      vrelohi[idx] = resultlohi;    
      vrehilo[idx] = resulthilo;
      vrelolo[idx] = resultlolo;    
   }
   // zim = vim[i]*inv0re + vre[i]*inv0im;

   qdg_mul(  shvimhihi[tdx],shvimlohi[tdx],shvimhilo[tdx],shvimlolo[tdx],
             invrehihi,     invrelohi,     invrehilo,     invrelolo,
           &resulthihi,   &resultlohi,   &resulthilo,   &resultlolo);
   qdg_mul(  shvrehihi[tdx],shvrelohi[tdx],shvrehilo[tdx],shvrelolo[tdx],
             invimhihi,     invimlohi,     invimhilo,     invimlolo,
              &acchihi,      &acclohi,      &acchilo,      &acclolo);
   qdg_inc(&resulthihi,&resultlohi,&resulthilo,&resultlolo,
               acchihi,    acclohi,    acchilo,    acclolo);

   __syncthreads();
   if(idx < dim)
   {
      vimhihi[idx] = resulthihi;
      vimlohi[idx] = resultlohi;    
      vimhilo[idx] = resulthilo;    
      vimlolo[idx] = resultlolo;
   }
}

__global__ void dbl4_small_leftRupdate
 ( int nrows, int ncols, int szt, int k,
   double *Rhihi, double *Rlohi, double *Rhilo, double *Rlolo,
   double *vhihi, double *vlohi, double *vhilo, double *vlolo,
   double *betahihi, double *betalohi,
   double *betahilo, double *betalolo )
{
   const int tdx = threadIdx.x;          // index of thread in block
   const int Roffset = k*nrows + k;
   int Rcolidx;
   double whihi,wlohi,whilo,wlolo;
   double Rtdxhihi,Rtdxlohi,Rtdxhilo,Rtdxlolo;
   double acchihi,acclohi,acchilo,acclolo;

   __shared__ double shvhihi[qd_shmemsize]; // slice of v
   __shared__ double shvlohi[qd_shmemsize]; 
   __shared__ double shvhilo[qd_shmemsize]; 
   __shared__ double shvlolo[qd_shmemsize]; 

   shvhihi[tdx] = vhihi[tdx];
   shvlohi[tdx] = vlohi[tdx];
   shvhilo[tdx] = vhilo[tdx];
   shvlolo[tdx] = vlolo[tdx];
   __syncthreads();
   whihi = 0.0;
   wlohi = 0.0;
   whilo = 0.0;
   wlolo = 0.0;

   for(int i=0; i<nrows-k; i++)   // loop through rows of R
   {
      Rcolidx = Roffset + i + tdx*nrows;
      __syncthreads();
      Rtdxhihi = Rhihi[Rcolidx];
      Rtdxlohi = Rlohi[Rcolidx];
      Rtdxhilo = Rhilo[Rcolidx];
      Rtdxlolo = Rlolo[Rcolidx];
      // w = w + Rtdx*shv[i];
      __syncthreads();
      qdg_mul(Rtdxhihi,  Rtdxlohi,  Rtdxhilo,  Rtdxlolo,
               shvhihi[i],shvlohi[i],shvhilo[i],shvlolo[i],
              &acchihi,  &acclohi,  &acchilo,  &acclolo);
      qdg_inc( &whihi, &wlohi, &whilo, &wlolo,
              acchihi,acclohi,acchilo,acclolo);
   }
   // w = (*beta)*w;
   // qdg_mlt(&whi,&wlo,*betahi,*betalo); <-- this does not work!
   __syncthreads();
   qdg_mul(*betahihi,*betalohi,*betahilo,*betalolo,
               whihi,    wlohi,    whilo,    wlolo,
            &acchihi, &acclohi, &acchilo, &acclolo);
   whihi = acchihi;
   wlohi = acclohi;
   whilo = acchilo;
   wlolo = acclolo;
   __syncthreads();
   for(int i=0; i<nrows-k; i++)   // update i-th row of R
   {
      Rcolidx = Roffset + i + tdx*nrows;
      __syncthreads();
      Rtdxhihi = Rhihi[Rcolidx];
      Rtdxlohi = Rlohi[Rcolidx];
      Rtdxhilo = Rhilo[Rcolidx];
      Rtdxlolo = Rlolo[Rcolidx];
      // Rtdx = Rtdx - shv[i]*w;
      __syncthreads();
      qdg_mul(shvhihi[i],shvlohi[i],shvhilo[i],shvlolo[i],
                whihi,     wlohi,     whilo,     wlolo,
             &acchihi,  &acclohi,  &acchilo,  &acclolo);
      qdg_dec(&Rtdxhihi,&Rtdxlohi,&Rtdxhilo,&Rtdxlolo,
                acchihi,  acclohi,  acchilo,  acclolo);
      __syncthreads();
      // changed nrows-k into ncols-k, where ncols = szt
      if(tdx < ncols-k)
      {
         Rhihi[Rcolidx] = Rtdxhihi;
         Rlohi[Rcolidx] = Rtdxlohi;
         Rhilo[Rcolidx] = Rtdxhilo;
         Rlolo[Rcolidx] = Rtdxlolo;
      }
      __syncthreads();
   }
}

__global__ void cmplx4_small_leftRupdate
 ( int nrows, int ncols, int szt, int k,
   double *Rrehihi, double *Rrelohi, double *Rrehilo, double *Rrelolo,
   double *Rimhihi, double *Rimlohi, double *Rimhilo, double *Rimlolo,
   double *vrehihi, double *vrelohi, double *vrehilo, double *vrelolo,
   double *vimhihi, double *vimlohi, double *vimhilo, double *vimlolo,
   double *betahihi, double *betalohi, double *betahilo, double *betalolo )
{
   const int tdx = threadIdx.x;          // index of thread in block
   const int Roffset = k*nrows + k;
   const double bthihi = *betahihi;
   const double btlohi = *betalohi;
   const double bthilo = *betahilo;
   const double btlolo = *betalolo;

   int Rcolidx;
   double w_rehihi,w_relohi,w_rehilo,w_relolo;
   double w_imhihi,w_imlohi,w_imhilo,w_imlolo;
   double Rtdx_hihi,Rtdx_lohi,Rtdx_hilo,Rtdx_lolo;
   double acchihi,acclohi,acchilo,acclolo;

   __shared__ double shvrehihi[cqd_shmemsize]; // slice of v
   __shared__ double shvrelohi[cqd_shmemsize];
   __shared__ double shvrehilo[cqd_shmemsize];
   __shared__ double shvrelolo[cqd_shmemsize];
   __shared__ double shvimhihi[cqd_shmemsize];
   __shared__ double shvimlohi[cqd_shmemsize];
   __shared__ double shvimhilo[cqd_shmemsize];
   __shared__ double shvimlolo[cqd_shmemsize];

   shvrehihi[tdx] = vrehihi[tdx];
   shvrelohi[tdx] = vrelohi[tdx];
   shvrehilo[tdx] = vrehilo[tdx];
   shvrelolo[tdx] = vrelolo[tdx];
   shvimhihi[tdx] = vimhihi[tdx];
   shvimlohi[tdx] = vimlohi[tdx];
   shvimhilo[tdx] = vimhilo[tdx];
   shvimlolo[tdx] = vimlolo[tdx];
   __syncthreads();
   w_rehihi = 0.0;
   w_relohi = 0.0;
   w_rehilo = 0.0;
   w_relolo = 0.0;
   w_imhihi = 0.0;
   w_imlohi = 0.0;
   w_imhilo = 0.0;
   w_imlolo = 0.0;

   for(int i=0; i<nrows-k; i++)   // loop through rows of R
   {
      Rcolidx = Roffset + i + tdx*nrows;
      __syncthreads();
      Rtdx_hihi = Rrehihi[Rcolidx];
      Rtdx_lohi = Rrelohi[Rcolidx];
      Rtdx_hilo = Rrehilo[Rcolidx];
      Rtdx_lolo = Rrelolo[Rcolidx];
      // w = w + Rtdx*shv[i]; beware of the Hermitian transpose!
      // w_re = w_re + Rtdx_re*shvre[i] + Rtdx_im*shvim[i];
      __syncthreads();
      qdg_mul(Rtdx_hihi,   Rtdx_lohi,   Rtdx_hilo,   Rtdx_lolo,
              shvrehihi[i],shvrelohi[i],shvrehilo[i],shvrelolo[i],
               &acchihi,    &acclohi,    &acchilo,    &acclolo);
      qdg_inc(&w_rehihi,&w_relohi,&w_rehilo,&w_relolo,
                acchihi,  acclohi,  acchilo,  acclolo);
      Rtdx_hihi = Rimhihi[Rcolidx];
      Rtdx_lohi = Rimlohi[Rcolidx];
      Rtdx_hilo = Rimhilo[Rcolidx];
      Rtdx_lolo = Rimlolo[Rcolidx];
      __syncthreads();
      qdg_mul(Rtdx_hihi,   Rtdx_lohi,   Rtdx_hilo,   Rtdx_lolo,
              shvimhihi[i],shvimlohi[i],shvimhilo[i],shvimlolo[i],
               &acchihi,    &acclohi,    &acchilo,    &acclolo);
      qdg_inc(&w_rehihi,&w_relohi,&w_rehilo,&w_relolo,
                acchihi,  acclohi,  acchilo,  acclolo);
      // w_im = w_im - Rtdx_im*shvre[i] + Rtdx_re*shvim[i];
      qdg_mul(Rtdx_hihi,   Rtdx_lohi,   Rtdx_hilo,   Rtdx_lolo,
              shvrehihi[i],shvrelohi[i],shvrehilo[i],shvrelolo[i],
               &acchihi,    &acclohi,    &acchilo,    &acclolo);
      qdg_dec(&w_imhihi,&w_imlohi,&w_imhilo,&w_imlolo,
                acchihi,  acclohi,  acchilo,  acclolo);

      Rtdx_hihi = Rrehihi[Rcolidx];
      Rtdx_lohi = Rrelohi[Rcolidx];
      Rtdx_hilo = Rrehilo[Rcolidx];
      Rtdx_lolo = Rrelolo[Rcolidx];
      qdg_mul(Rtdx_hihi,   Rtdx_lohi,   Rtdx_hilo,   Rtdx_lolo,
              shvimhihi[i],shvimlohi[i],shvimhilo[i],shvimlolo[i],
               &acchihi,    &acclohi,    &acchilo,    &acclolo);
      qdg_inc(&w_imhihi,&w_imlohi,&w_imhilo,&w_imlolo,
                acchihi,  acclohi,  acchilo,  acclolo);
   }
   // w_re = acc*w_re;
   __syncthreads();
   qdg_mul(w_rehihi,w_relohi,w_rehilo,w_relolo,
             bthihi,  btlohi,  bthilo,  btlolo,
           &acchihi,&acclohi,&acchilo,&acclolo);
   w_rehihi = acchihi;
   w_relohi = acclohi;
   w_rehilo = acchilo;
   w_relolo = acclolo;
   // w_im = acc*w_im;
   __syncthreads();
   qdg_mul(w_imhihi,w_imlohi,w_imhilo,w_imlolo,
             bthihi,  btlohi,  bthilo,  btlolo,
           &acchihi,&acclohi,&acchilo,&acclolo);
   w_imhihi = acchihi;
   w_imlohi = acclohi;
   w_imhilo = acchilo;
   w_imlolo = acclolo;
   __syncthreads();
   for(int i=0; i<nrows-k; i++)   // update i-th row of Rre
   {
      Rcolidx = Roffset + i + tdx*nrows;
      Rtdx_hihi = Rrehihi[Rcolidx];
      Rtdx_lohi = Rrelohi[Rcolidx];
      Rtdx_hilo = Rrehilo[Rcolidx];
      Rtdx_lolo = Rrelolo[Rcolidx];
      // Rtdx = Rtdx - shv[i]*w; beware of the Hermitian transpose!
      // Rtdx_re = Rtdx_re - (shvre[i]*w_re + shvim[i]*w_im);
      __syncthreads();
      qdg_mul(shvrehihi[i],shvrelohi[i],shvrehilo[i],shvrelolo[i],
               w_rehihi,    w_relohi,    w_rehilo,    w_relolo,
               &acchihi,    &acclohi,    &acchilo,    &acclolo);
      qdg_dec(&Rtdx_hihi,&Rtdx_lohi,&Rtdx_hilo,&Rtdx_lolo,
                 acchihi,   acclohi,   acchilo,   acclolo);
      qdg_mul(shvimhihi[i],shvimlohi[i],shvimhilo[i],shvimlolo[i],
               w_imhihi,    w_imlohi,    w_imhilo,    w_imlolo,
               &acchihi,    &acclohi,    &acchilo,    &acclolo);
      qdg_dec(&Rtdx_hihi,&Rtdx_lohi,&Rtdx_hilo,&Rtdx_lolo,
                 acchihi,   acclohi,   acchilo,   acclolo);
      __syncthreads();
      // changed nrows-k into ncols-k, where ncols = szt
      if(tdx < ncols-k)
      {
         Rrehihi[Rcolidx] = Rtdx_hihi;
         Rrelohi[Rcolidx] = Rtdx_lohi;
         Rrehilo[Rcolidx] = Rtdx_hilo;
         Rrelolo[Rcolidx] = Rtdx_lolo;
      }
   }
   __syncthreads();
   for(int i=0; i<nrows-k; i++)   // update i-th row of Rim
   {
      Rcolidx = Roffset + i + tdx*nrows;
      Rtdx_hihi = Rimhihi[Rcolidx];
      Rtdx_lohi = Rimlohi[Rcolidx];
      Rtdx_hilo = Rimhilo[Rcolidx];
      Rtdx_lolo = Rimlolo[Rcolidx];
      // Rtdx_im = Rtdx_im - (shvim[i]*w_re - shvre[i]*w_im);
      __syncthreads();
      qdg_mul(shvimhihi[i],shvimlohi[i],shvimhilo[i],shvimlolo[i],
               w_rehihi,    w_relohi,    w_rehilo,    w_relolo,
               &acchihi,    &acclohi,    &acchilo,    &acclolo);
      qdg_dec(&Rtdx_hihi,&Rtdx_lohi,&Rtdx_hilo,&Rtdx_lolo,
                 acchihi,   acclohi,   acchilo,   acclolo);
      qdg_mul(shvrehihi[i],shvrelohi[i],shvrehilo[i],shvrelolo[i],
               w_imhihi,    w_imlohi,    w_imhilo,    w_imlolo,
               &acchihi,    &acclohi,    &acchilo,    &acclolo);
      qdg_inc(&Rtdx_hihi,&Rtdx_lohi,&Rtdx_hilo,&Rtdx_lolo,
                 acchihi,   acclohi,   acchilo,   acclolo);
      __syncthreads();
      // changed nrows-k into ncols-k, where ncols = szt
      if(tdx < ncols-k)
      {
         Rimhihi[Rcolidx] = Rtdx_hihi;
         Rimlohi[Rcolidx] = Rtdx_lohi;
         Rimhilo[Rcolidx] = Rtdx_hilo;
         Rimlolo[Rcolidx] = Rtdx_lolo;
      }
   }
}

__global__ void dbl4_RTdotv
 ( int nrows, int szt, int colidx, int Roffset, int dim,
   double *Rhihi, double *Rlohi, double *Rhilo, double *Rlolo,
   double *vhihi, double *vlohi, double *vhilo, double *vlolo,
   double *RTdotvhihi, double *RTdotvlohi,
   double *RTdotvhilo, double *RTdotvlolo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int idx = bdx*szt + tdx;        // thread tdx computes RTv[idx]

   const int vdx = idx % nrows;          // index in v is column in R^T
   const int row = idx / nrows;          // R is stored column-by-column

   const int Rdx = Roffset + idx + (row+1)*colidx;

   const double Vvalhihi = vhihi[vdx];
   const double Vvallohi = vlohi[vdx];
   const double Vvalhilo = vhilo[vdx];
   const double Vvallolo = vlolo[vdx];
   const double Rvalhihi = Rhihi[Rdx];
   const double Rvallohi = Rlohi[Rdx];
   const double Rvalhilo = Rhilo[Rdx];
   const double Rvallolo = Rlolo[Rdx];
   // double result = Rval*Vval;
   double resulthihi,resultlohi,resulthilo,resultlolo;

   __syncthreads();
   qdg_mul(   Rvalhihi,   Rvallohi,   Rvalhilo,   Rvallolo,
              Vvalhihi,   Vvallohi,   Vvalhilo,   Vvallolo,
           &resulthihi,&resultlohi,&resulthilo,&resultlolo);
   __syncthreads();
   RTdotvhihi[idx] = resulthihi;
   RTdotvlohi[idx] = resultlohi;
   RTdotvhilo[idx] = resulthilo;
   RTdotvlolo[idx] = resultlolo;
}

__global__ void cmplx4_RHdotv
 ( int nrows, int szt, int colidx, int Roffset, int dim,
   double *Rrehihi, double *Rrelohi, double *Rrehilo, double *Rrelolo,
   double *Rimhihi, double *Rimlohi, double *Rimhilo, double *Rimlolo,
   double *vrehihi, double *vrelohi, double *vrehilo, double *vrelolo,
   double *vimhihi, double *vimlohi, double *vimhilo, double *vimlolo,
   double *RHdotvrehihi, double *RHdotvrelohi,
   double *RHdotvrehilo, double *RHdotvrelolo,
   double *RHdotvimhihi, double *RHdotvimlohi,
   double *RHdotvimhilo, double *RHdotvimlolo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int idx = bdx*szt + tdx;        // thread tdx computes RTv[idx]

   const int vdx = idx % nrows;          // index in v is column in R^T
   const int row = idx / nrows;          // R is stored column-by-column

   const int Rdx = Roffset + idx + (row+1)*colidx;

   double Vvalhihi,Vvallohi,Vvalhilo,Vvallolo;
   double Rvalhihi,Rvallohi,Rvalhilo,Rvallolo;
   // double result = Rval*Vval;
   double resulthihi,resultlohi,resulthilo,resultlolo;
   double acchihi,acclohi,acchilo,acclolo;

   Vvalhihi = vrehihi[vdx];
   Vvallohi = vrelohi[vdx];
   Vvalhilo = vrehilo[vdx];
   Vvallolo = vrelolo[vdx];
   Rvalhihi = Rrehihi[Rdx];
   Rvallohi = Rrelohi[Rdx];
   Rvalhilo = Rrehilo[Rdx];
   Rvallolo = Rrelolo[Rdx];
   __syncthreads();
   qdg_mul(   Rvalhihi,   Rvallohi,   Rvalhilo,   Rvallolo,
              Vvalhihi,   Vvallohi,   Vvalhilo,   Vvallolo,
           &resulthihi,&resultlohi,&resulthilo,&resultlolo);

   Vvalhihi = vimhihi[vdx];
   Vvallohi = vimlohi[vdx];
   Vvalhilo = vimhilo[vdx];
   Vvallolo = vimlolo[vdx];
   Rvalhihi = Rimhihi[Rdx];
   Rvallohi = Rimlohi[Rdx];
   Rvalhilo = Rimhilo[Rdx];
   Rvallolo = Rimlolo[Rdx];
   __syncthreads();
   qdg_mul(Rvalhihi,Rvallohi,Rvalhilo,Rvallolo,
           Vvalhihi,Vvallohi,Vvalhilo,Vvallolo,
           &acchihi,&acclohi,&acchilo,&acclolo);
   qdg_inc(&resulthihi,&resultlohi,&resulthilo,&resultlolo,
               acchihi,    acclohi,    acchilo,    acclolo);

   __syncthreads();
   RHdotvrehihi[idx] = resulthihi;
   RHdotvrelohi[idx] = resultlohi;
   RHdotvrehilo[idx] = resulthilo;
   RHdotvrelolo[idx] = resultlolo;

   // Vvalim is already available
   Rvalhihi = Rrehihi[Rdx];
   Rvallohi = Rrelohi[Rdx];
   Rvalhilo = Rrehilo[Rdx];
   Rvallolo = Rrelolo[Rdx];
   __syncthreads();
   qdg_mul(   Rvalhihi,   Rvallohi,   Rvalhilo,   Rvallolo,
              Vvalhihi,   Vvallohi,   Vvalhilo,   Vvallolo,
           &resulthihi,&resultlohi,&resulthilo,&resultlolo);

   Vvalhihi = vrehihi[vdx];
   Vvallohi = vrelohi[vdx];
   Vvalhilo = vrehilo[vdx];
   Vvallolo = vrelolo[vdx];
   Rvalhihi = Rimhihi[Rdx];
   Rvallohi = Rimlohi[Rdx];
   Rvalhilo = Rimhilo[Rdx];
   Rvallolo = Rimlolo[Rdx];
   __syncthreads();
   qdg_mul(Rvalhihi,Rvallohi,Rvalhilo,Rvallolo,
           Vvalhihi,Vvallohi,Vvalhilo,Vvallolo,
           &acchihi,&acclohi,&acchilo,&acclolo);
   qdg_dec(&resulthihi,&resultlohi,&resulthilo,&resultlolo,
               acchihi,    acclohi,    acchilo,    acclolo);

   __syncthreads();
   RHdotvimhihi[idx] = resulthihi;
   RHdotvimlohi[idx] = resultlohi;
   RHdotvimhilo[idx] = resulthilo;
   RHdotvimlolo[idx] = resultlolo;
}

__global__ void cmplx4_RHdotvre
 ( int nrows, int szt, int colidx, int Roffset, int dim,
   double *Rrehihi, double *Rrelohi, double *Rrehilo, double *Rrelolo,
   double *Rimhihi, double *Rimlohi, double *Rimhilo, double *Rimlolo,
   double *vrehihi, double *vrelohi, double *vrehilo, double *vrelolo,
   double *vimhihi, double *vimlohi, double *vimhilo, double *vimlolo,
   double *RHdotvrehihi, double *RHdotvrelohi,
   double *RHdotvrehilo, double *RHdotvrelolo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int idx = bdx*szt + tdx;        // thread tdx computes RTv[idx]

   const int vdx = idx % nrows;          // index in v is column in R^T
   const int row = idx / nrows;          // R is stored column-by-column

   const int Rdx = Roffset + idx + (row+1)*colidx;

   double Vvalhihi,Vvallohi,Vvalhilo,Vvallolo;
   double Rvalhihi,Rvallohi,Rvalhilo,Rvallolo;
   // double result = Rval*Vval;
   double resulthihi,resultlohi,resulthilo,resultlolo;
   double acchihi,acclohi,acchilo,acclolo;

   Vvalhihi = vrehihi[vdx];
   Vvallohi = vrelohi[vdx];
   Vvalhilo = vrehilo[vdx];
   Vvallolo = vrelolo[vdx];
   Rvalhihi = Rrehihi[Rdx];
   Rvallohi = Rrelohi[Rdx];
   Rvalhilo = Rrehilo[Rdx];
   Rvallolo = Rrelolo[Rdx];
   __syncthreads();
   qdg_mul(   Rvalhihi,   Rvallohi,   Rvalhilo,   Rvallolo,
              Vvalhihi,   Vvallohi,   Vvalhilo,   Vvallolo,
           &resulthihi,&resultlohi,&resulthilo,&resultlolo);
   __syncthreads();
   Vvalhihi = vimhihi[vdx];
   Vvallohi = vimlohi[vdx];
   Vvalhilo = vimhilo[vdx];
   Vvallolo = vimlolo[vdx];
   Rvalhihi = Rimhihi[Rdx];
   Rvallohi = Rimlohi[Rdx];
   Rvalhilo = Rimhilo[Rdx];
   Rvallolo = Rimlolo[Rdx];
   __syncthreads();
   qdg_mul(Rvalhihi,Rvallohi,Rvalhilo,Rvallolo,
           Vvalhihi,Vvallohi,Vvalhilo,Vvallolo,
           &acchihi,&acclohi,&acchilo,&acclolo);
   qdg_inc(&resulthihi,&resultlohi,&resulthilo,&resultlolo,
               acchihi,    acclohi,    acchilo,    acclolo);

   __syncthreads();
   RHdotvrehihi[idx] = resulthihi;
   RHdotvrelohi[idx] = resultlohi;
   RHdotvrehilo[idx] = resulthilo;
   RHdotvrelolo[idx] = resultlolo;
}

__global__ void cmplx4_RHdotvim
 ( int nrows, int szt, int colidx, int Roffset, int dim,
   double *Rrehihi, double *Rrelohi, double *Rrehilo, double *Rrelolo,
   double *Rimhihi, double *Rimlohi, double *Rimhilo, double *Rimlolo,
   double *vrehihi, double *vrelohi, double *vrehilo, double *vrelolo,
   double *vimhihi, double *vimlohi, double *vimhilo, double *vimlolo,
   double *RHdotvimhihi, double *RHdotvimlohi,
   double *RHdotvimhilo, double *RHdotvimlolo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int idx = bdx*szt + tdx;        // thread tdx computes RTv[idx]

   const int vdx = idx % nrows;          // index in v is column in R^T
   const int row = idx / nrows;          // R is stored column-by-column

   const int Rdx = Roffset + idx + (row+1)*colidx;

   double Vvalhihi,Vvallohi,Vvalhilo,Vvallolo;
   double Rvalhihi,Rvallohi,Rvalhilo,Rvallolo;
   // double result = Rval*Vval;
   double resulthihi,resultlohi,resulthilo,resultlolo;
   double acchihi,acclohi,acchilo,acclolo;

   Vvalhihi = vimhihi[vdx];
   Vvallohi = vimlohi[vdx];
   Vvalhilo = vimhilo[vdx];
   Vvallolo = vimlolo[vdx];
   Rvalhihi = Rrehihi[Rdx];
   Rvallohi = Rrelohi[Rdx];
   Rvalhilo = Rrehilo[Rdx];
   Rvallolo = Rrelolo[Rdx];
   __syncthreads();
   qdg_mul(   Rvalhihi,   Rvallohi,   Rvalhilo,   Rvallolo,
              Vvalhihi,   Vvallohi,   Vvalhilo,   Vvallolo,
           &resulthihi,&resultlohi,&resulthilo,&resultlolo);
   __syncthreads();
   Vvalhihi = vrehihi[vdx];
   Vvallohi = vrelohi[vdx];
   Vvalhilo = vrehilo[vdx];
   Vvallolo = vrelolo[vdx];
   Rvalhihi = Rimhihi[Rdx];
   Rvallohi = Rimlohi[Rdx];
   Rvalhilo = Rimhilo[Rdx];
   Rvallolo = Rimlolo[Rdx];
   __syncthreads();
   qdg_mul(Rvalhihi,Rvallohi,Rvalhilo,Rvallolo,
           Vvalhihi,Vvallohi,Vvalhilo,Vvallolo,
           &acchihi,&acclohi,&acchilo,&acclolo);
   qdg_dec(&resulthihi,&resultlohi,&resulthilo,&resultlolo,
               acchihi,    acclohi,    acchilo,    acclolo);

   __syncthreads();
   RHdotvimhihi[idx] = resulthihi;
   RHdotvimlohi[idx] = resultlohi;
   RHdotvimhilo[idx] = resulthilo;
   RHdotvimlolo[idx] = resultlolo;
}

__global__ void dbl4_sum_betaRTdotv
 ( int nrows,
   double *betahihi, double *betalohi, double *betahilo, double *betalolo,
   double *RTdotvhihi, double *RTdotvlohi,
   double *RTdotvhilo, double *RTdotvlolo,
   double *whihi, double *wlohi, double *whilo, double *wlolo )
{
   const int tdx = threadIdx.x;  // tdx sums elements on row tdx
   const int offset = tdx*nrows; // number of rows before current row
   int idx;

   double resulthihi = 0.0;
   double resultlohi = 0.0;
   double resulthilo = 0.0;
   double resultlolo = 0.0;
   double Rvalhihi,Rvallohi,Rvalhilo,Rvallolo;

   for(int i=0; i<nrows; i++)
   {
      idx = offset + i;
      __syncthreads();
      Rvalhihi = RTdotvhihi[idx];
      Rvallohi = RTdotvlohi[idx];
      Rvalhilo = RTdotvhilo[idx];
      Rvallolo = RTdotvlolo[idx];
      // result = result + Rval;
      __syncthreads();
      qdg_inc(&resulthihi,&resultlohi,&resulthilo,&resultlolo,
                 Rvalhihi,   Rvallohi,   Rvalhilo,   Rvallolo);
   }
   __syncthreads();
   Rvalhihi = *betahihi;
   Rvallohi = *betalohi;
   Rvalhilo = *betahilo;
   Rvallolo = *betalolo;
   // w[tdx] = Rval*result;
   __syncthreads();
   qdg_mul(  Rvalhihi,   Rvallohi,   Rvalhilo,   Rvallolo,
           resulthihi, resultlohi, resulthilo, resultlolo,
               &whihi[tdx],&wlohi[tdx],&whilo[tdx],&wlolo[tdx]);
}

__global__ void cmplx4_sum_betaRHdotv
 ( int nrows,
   double *betahihi, double *betalohi, double *betahilo, double *betalolo,
   double *RTdotvrehihi, double *RTdotvrelohi,
   double *RTdotvrehilo, double *RTdotvrelolo,
   double *RTdotvimhihi, double *RTdotvimlohi,
   double *RTdotvimhilo, double *RTdotvimlolo,
   double *wrehihi, double *wrelohi, double *wrehilo, double *wrelolo,
   double *wimhihi, double *wimlohi, double *wimhilo, double *wimlolo )
{
   const int tdx = threadIdx.x;  // tdx sums elements on row tdx
   const int offset = tdx*nrows; // number of rows before current row
   int idx;

   double resulthihi = 0.0;
   double resultlohi = 0.0;
   double resulthilo = 0.0;
   double resultlolo = 0.0;
   double Rvalhihi,Rvallohi,Rvalhilo,Rvallolo;
   double acchihi,acclohi,acchilo,acclolo;

   for(int i=0; i<nrows; i++)
   {
      idx = offset + i;
      Rvalhihi = RTdotvrehihi[idx];
      Rvallohi = RTdotvrelohi[idx];
      Rvalhilo = RTdotvrehilo[idx];
      Rvallolo = RTdotvrelolo[idx];
      // result = result + Rval;
      __syncthreads();
      qdg_inc(&resulthihi,&resultlohi,&resulthilo,&resultlolo,
                 Rvalhihi,   Rvallohi,   Rvalhilo,   Rvallolo);
   }
   __syncthreads();
   Rvalhihi = *betahihi;
   Rvallohi = *betalohi;
   Rvalhilo = *betahilo;
   Rvallolo = *betalolo;
   // w[tdx] = Rval*result;
   __syncthreads();
   qdg_mul(  Rvalhihi,  Rvallohi,  Rvalhilo,  Rvallolo,
           resulthihi,resultlohi,resulthilo,resultlolo,
             &acchihi,  &acclohi,  &acchilo,  &acclolo);

   wrehihi[tdx] = acchihi;
   wrelohi[tdx] = acclohi;
   wrehilo[tdx] = acchilo;
   wrelolo[tdx] = acclolo;

   resulthihi = 0.0;
   resultlohi = 0.0;
   resulthilo = 0.0;
   resultlolo = 0.0;

   for(int i=0; i<nrows; i++)
   {
      idx = offset + i;
      Rvalhihi = RTdotvimhihi[idx];
      Rvallohi = RTdotvimlohi[idx];
      Rvalhilo = RTdotvimhilo[idx];
      Rvallolo = RTdotvimlolo[idx];
      qdg_inc(&resulthihi,&resultlohi,&resulthilo,&resultlolo,
                 Rvalhihi,   Rvallohi,   Rvalhilo,   Rvallolo);
   }
   __syncthreads();
   Rvalhihi = *betahihi;
   Rvallohi = *betalohi;
   Rvalhilo = *betahilo;
   Rvallolo = *betalolo;
   // w[tdx] = Rval*result;
   __syncthreads();
   qdg_mul(  Rvalhihi,  Rvallohi,  Rvalhilo,  Rvallolo,
           resulthihi,resultlohi,resulthilo,resultlolo,
             &acchihi,  &acclohi,  &acchilo,  &acclolo);

   wimhihi[tdx] = acchihi;
   wimlohi[tdx] = acclohi;
   wimhilo[tdx] = acchilo;
   wimlolo[tdx] = acclolo;
}

__global__ void dbl4_medium_subvbetaRTv
 ( int nrows, int ncols, int szt, int k,
   double *Rhihi, double *Rlohi, double *Rhilo, double *Rlolo,
   double *vhihi, double *vlohi, double *vhilo, double *vlolo,
   double *betahihi, double *betalohi, double *betahilo, double *betalolo,
   double *whihi, double *wlohi, double *whilo, double *wlolo )
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

   __shared__ double shwhihi[qd_shmemsize];  // values in beta*R^T*v
   __shared__ double shwlohi[qd_shmemsize];  // are less in number than szt
   __shared__ double shwhilo[qd_shmemsize];
   __shared__ double shwlolo[qd_shmemsize];
   shwhihi[tdx] = whihi[tdx];
   shwlohi[tdx] = wlohi[tdx];
   shwhilo[tdx] = whilo[tdx];
   shwlolo[tdx] = wlolo[tdx];
   __syncthreads();

   double Rwidxhihi = Rhihi[Ridx];     // number that tdx updates
   double Rwidxlohi = Rlohi[Ridx];
   double Rwidxhilo = Rhilo[Ridx];
   double Rwidxlolo = Rlolo[Ridx];
   double vValhihi = vhihi[rowidx];    // value in Householder vector
   double vVallohi = vlohi[rowidx];
   double vValhilo = vhilo[rowidx];
   double vVallolo = vlolo[rowidx];
   double wValhihi = shwhihi[colidx];  // value in beta*R^T*v
   double wVallohi = shwlohi[colidx];
   double wValhilo = shwhilo[colidx];
   double wVallolo = shwlolo[colidx];
   double acchihi,acclohi,acchilo,acclolo;
 
   // Rwidx = Rwidx - vValue*wValue;   // update R[rowidx,colidx]
   __syncthreads();
   qdg_mul(vValhihi,vVallohi,vValhilo,vVallolo,
           wValhihi,wVallohi,wValhilo,wVallolo,
           &acchihi,&acclohi,&acchilo,&acclolo);
   qdg_dec(&Rwidxhihi,&Rwidxlohi,&Rwidxhilo,&Rwidxlolo,
              acchihi,   acclohi,   acchilo,   acclolo);
   __syncthreads();
   if(widx < bound)                    // if() takes care of padding
   {
      Rhihi[Ridx] = Rwidxhihi;
      Rlohi[Ridx] = Rwidxlohi;
      Rhilo[Ridx] = Rwidxhilo;
      Rlolo[Ridx] = Rwidxlolo;
   }
}

__global__ void cmplx4_medium_subvbetaRHv
 ( int nrows, int ncols, int szt, int k,
   double *Rrehihi, double *Rrelohi, double *Rrehilo, double *Rrelolo,
   double *Rimhihi, double *Rimlohi, double *Rimhilo, double *Rimlolo,
   double *vrehihi, double *vrelohi, double *vrehilo, double *vrelolo,
   double *vimhihi, double *vimlohi, double *vimhilo, double *vimlolo,
   double *betahihi, double *betalohi, double *betahilo, double *betalolo,
   double *wrehihi, double *wrelohi, double *wrehilo, double *wrelolo,
   double *wimhihi, double *wimlohi, double *wimhilo, double *wimlolo )
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

   __shared__ double shwrehihi[cqd_shmemsize];  // values in beta*R^H*v
   __shared__ double shwrelohi[cqd_shmemsize];  // are less in number than szt
   __shared__ double shwrehilo[cqd_shmemsize]; 
   __shared__ double shwrelolo[cqd_shmemsize]; 
   __shared__ double shwimhihi[cqd_shmemsize]; 
   __shared__ double shwimlohi[cqd_shmemsize];
   __shared__ double shwimhilo[cqd_shmemsize]; 
   __shared__ double shwimlolo[cqd_shmemsize];
   shwrehihi[tdx] = wrehihi[tdx];
   shwrelohi[tdx] = wrelohi[tdx];
   shwrehilo[tdx] = wrehilo[tdx];
   shwrelolo[tdx] = wrelolo[tdx];
   shwimhihi[tdx] = wimhihi[tdx];
   shwimlohi[tdx] = wimlohi[tdx];
   shwimhilo[tdx] = wimhilo[tdx];
   shwimlolo[tdx] = wimlolo[tdx];
   __syncthreads();

   double Rwidxhihi,Rwidxlohi,Rwidxhilo,Rwidxlolo;
   double vValhihi,vVallohi,vValhilo,vVallolo; // value in Householder vector
   double wValhihi,wVallohi,wValhilo,wVallolo; // value in beta*R^H*v
   double acchihi,acclohi,acchilo,acclolo;
 
   // Rwidx = Rwidx - vValue*wValue;   // update R[rowidx,colidx]
   // take the Hermitian transpose of w

   Rwidxhihi = Rrehihi[Ridx];     // number that tdx updates
   Rwidxlohi = Rrelohi[Ridx];
   Rwidxhilo = Rrehilo[Ridx];
   Rwidxlolo = Rrelolo[Ridx];

   __syncthreads();
   vValhihi = vrehihi[rowidx];    // value in Householder vector
   vVallohi = vrelohi[rowidx];
   vValhilo = vrehilo[rowidx];
   vVallolo = vrelolo[rowidx];
   wValhihi = shwrehihi[colidx];  // value in beta*R^H*v
   wVallohi = shwrelohi[colidx];
   wValhilo = shwrehilo[colidx];
   wVallolo = shwrelolo[colidx];
   __syncthreads();
   qdg_mul(vValhihi,vVallohi,vValhilo,vVallolo,
           wValhihi,wVallohi,wValhilo,wVallolo,
           &acchihi,&acclohi,&acchilo,&acclolo);
   qdg_dec(&Rwidxhihi,&Rwidxlohi,&Rwidxhilo,&Rwidxlolo,
              acchihi,   acclohi,   acchilo,   acclolo);

   vValhihi = vimhihi[rowidx];
   vVallohi = vimlohi[rowidx];
   vValhilo = vimhilo[rowidx];
   vVallolo = vimlolo[rowidx];
   wValhihi = shwimhihi[colidx];
   wVallohi = shwimlohi[colidx];
   wValhilo = shwimhilo[colidx];
   wVallolo = shwimlolo[colidx];
   __syncthreads();
   qdg_mul(vValhihi,vVallohi,vValhilo,vVallolo,
           wValhihi,wVallohi,wValhilo,wVallolo,
           &acchihi,&acclohi,&acchilo,&acclolo);
   qdg_dec(&Rwidxhihi,&Rwidxlohi,&Rwidxhilo,&Rwidxlolo,
              acchihi,   acclohi,   acchilo,   acclolo);

   __syncthreads();
   if(widx < bound)                    // if() takes care of padding
   {
      Rrehihi[Ridx] = Rwidxhihi;
      Rrelohi[Ridx] = Rwidxlohi;
      Rrehilo[Ridx] = Rwidxhilo;
      Rrelolo[Ridx] = Rwidxlolo;
   }

   __syncthreads();
   Rwidxhihi = Rimhihi[Ridx];
   Rwidxlohi = Rimlohi[Ridx];
   Rwidxhilo = Rimhilo[Ridx];
   Rwidxlolo = Rimlolo[Ridx];

   __syncthreads();
   // vValim is already loaded
   wValhihi = shwrehihi[colidx];  // value in beta*R^H*v
   wVallohi = shwrelohi[colidx];
   wValhilo = shwrehilo[colidx];
   wVallolo = shwrelolo[colidx];
   __syncthreads();
   qdg_mul(vValhihi,vVallohi,vValhilo,vVallolo,
           wValhihi,wVallohi,wValhilo,wVallolo,
           &acchihi,&acclohi,&acchilo,&acclolo);
   qdg_dec(&Rwidxhihi,&Rwidxlohi,&Rwidxhilo,&Rwidxlolo,
              acchihi,   acclohi,   acchilo,   acclolo);

   vValhihi = vrehihi[rowidx];    // value in Householder vector
   vVallohi = vrelohi[rowidx];
   vValhilo = vrehilo[rowidx];
   vVallolo = vrelolo[rowidx];
   wValhihi = shwimhihi[colidx];
   wVallohi = shwimlohi[colidx];
   wValhilo = shwimhilo[colidx];
   wVallolo = shwimlolo[colidx];
   __syncthreads();
   qdg_mul(vValhihi,vVallohi,vValhilo,vVallolo,
           wValhihi,wVallohi,wValhilo,wVallolo,
           &acchihi,&acclohi,&acchilo,&acclolo);
   qdg_inc(&Rwidxhihi,&Rwidxlohi,&Rwidxhilo,&Rwidxlolo,
              acchihi,   acclohi,   acchilo,   acclolo);

   __syncthreads();
   if(widx < bound)                    // if() takes care of padding
   {
      Rimhihi[Ridx] = Rwidxhihi;
      Rimlohi[Ridx] = Rwidxlohi;
      Rimhilo[Ridx] = Rwidxhilo;
      Rimlolo[Ridx] = Rwidxlolo;
   }
}

__global__ void cmplx4_medium_subvbetaRHvRe
 ( int nrows, int ncols, int szt, int k,
   double *Rrehihi, double *Rrelohi, double *Rrehilo, double *Rrelolo,
   double *vrehihi, double *vrelohi, double *vrehilo, double *vrelolo,
   double *vimhihi, double *vimlohi, double *vimhilo, double *vimlolo,
   double *betahihi, double *betalohi, double *betahilo, double *betalolo,
   double *wrehihi, double *wrelohi, double *wrehilo, double *wrelolo,
   double *wimhihi, double *wimlohi, double *wimhilo, double *wimlolo )
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

   __shared__ double shwrehihi[cqd_shmemsize];  // values in beta*R^H*v
   __shared__ double shwrelohi[cqd_shmemsize];  // are less in number than szt
   __shared__ double shwrehilo[cqd_shmemsize]; 
   __shared__ double shwrelolo[cqd_shmemsize]; 
   __shared__ double shwimhihi[cqd_shmemsize]; 
   __shared__ double shwimlohi[cqd_shmemsize];
   __shared__ double shwimhilo[cqd_shmemsize]; 
   __shared__ double shwimlolo[cqd_shmemsize];
   shwrehihi[tdx] = wrehihi[tdx];
   shwrelohi[tdx] = wrelohi[tdx];
   shwrehilo[tdx] = wrehilo[tdx];
   shwrelolo[tdx] = wrelolo[tdx];
   shwimhihi[tdx] = wimhihi[tdx];
   shwimlohi[tdx] = wimlohi[tdx];
   shwimhilo[tdx] = wimhilo[tdx];
   shwimlolo[tdx] = wimlolo[tdx];
   __syncthreads();

   double Rwidxhihi,Rwidxlohi,Rwidxhilo,Rwidxlolo;
   double vValhihi,vVallohi,vValhilo,vVallolo; // value in Householder vector
   double wValhihi,wVallohi,wValhilo,wVallolo; // value in beta*R^H*v
   double acchihi,acclohi,acchilo,acclolo;
 
   // Rwidx = Rwidx - vValue*wValue;   // update R[rowidx,colidx]
   // take the Hermitian transpose of w

   Rwidxhihi = Rrehihi[Ridx];     // number that tdx updates
   Rwidxlohi = Rrelohi[Ridx];
   Rwidxhilo = Rrehilo[Ridx];
   Rwidxlolo = Rrelolo[Ridx];

   __syncthreads();
   vValhihi = vrehihi[rowidx];    // value in Householder vector
   vVallohi = vrelohi[rowidx];
   vValhilo = vrehilo[rowidx];
   vVallolo = vrelolo[rowidx];
   wValhihi = shwrehihi[colidx];  // value in beta*R^H*v
   wVallohi = shwrelohi[colidx];
   wValhilo = shwrehilo[colidx];
   wVallolo = shwrelolo[colidx];
   __syncthreads();
   qdg_mul(vValhihi,vVallohi,vValhilo,vVallolo,
           wValhihi,wVallohi,wValhilo,wVallolo,
           &acchihi,&acclohi,&acchilo,&acclolo);
   qdg_dec(&Rwidxhihi,&Rwidxlohi,&Rwidxhilo,&Rwidxlolo,
              acchihi,   acclohi,   acchilo,   acclolo);

   vValhihi = vimhihi[rowidx];
   vVallohi = vimlohi[rowidx];
   vValhilo = vimhilo[rowidx];
   vVallolo = vimlolo[rowidx];
   wValhihi = shwimhihi[colidx];
   wVallohi = shwimlohi[colidx];
   wValhilo = shwimhilo[colidx];
   wVallolo = shwimlolo[colidx];
   __syncthreads();
   qdg_mul(vValhihi,vVallohi,vValhilo,vVallolo,
           wValhihi,wVallohi,wValhilo,wVallolo,
           &acchihi,&acclohi,&acchilo,&acclolo);
   qdg_dec(&Rwidxhihi,&Rwidxlohi,&Rwidxhilo,&Rwidxlolo,
              acchihi,   acclohi,   acchilo,   acclolo);

   __syncthreads();
   if(widx < bound)                    // if() takes care of padding
   {
      Rrehihi[Ridx] = Rwidxhihi;
      Rrelohi[Ridx] = Rwidxlohi;
      Rrehilo[Ridx] = Rwidxhilo;
      Rrelolo[Ridx] = Rwidxlolo;
   }
}

__global__ void cmplx4_medium_subvbetaRHvIm
 ( int nrows, int ncols, int szt, int k,
   double *Rimhihi, double *Rimlohi, double *Rimhilo, double *Rimlolo,
   double *vrehihi, double *vrelohi, double *vrehilo, double *vrelolo,
   double *vimhihi, double *vimlohi, double *vimhilo, double *vimlolo,
   double *betahihi, double *betalohi, double *betahilo, double *betalolo,
   double *wrehihi, double *wrelohi, double *wrehilo, double *wrelolo,
   double *wimhihi, double *wimlohi, double *wimhilo, double *wimlolo )
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

   __shared__ double shwrehihi[cqd_shmemsize];  // values in beta*R^H*v
   __shared__ double shwrelohi[cqd_shmemsize];  // are less in number than szt
   __shared__ double shwrehilo[cqd_shmemsize]; 
   __shared__ double shwrelolo[cqd_shmemsize]; 
   __shared__ double shwimhihi[cqd_shmemsize]; 
   __shared__ double shwimlohi[cqd_shmemsize];
   __shared__ double shwimhilo[cqd_shmemsize]; 
   __shared__ double shwimlolo[cqd_shmemsize];
   shwrehihi[tdx] = wrehihi[tdx];
   shwrelohi[tdx] = wrelohi[tdx];
   shwrehilo[tdx] = wrehilo[tdx];
   shwrelolo[tdx] = wrelolo[tdx];
   shwimhihi[tdx] = wimhihi[tdx];
   shwimlohi[tdx] = wimlohi[tdx];
   shwimhilo[tdx] = wimhilo[tdx];
   shwimlolo[tdx] = wimlolo[tdx];
   __syncthreads();

   double Rwidxhihi,Rwidxlohi,Rwidxhilo,Rwidxlolo;
   double vValhihi,vVallohi,vValhilo,vVallolo; // value in Householder vector
   double wValhihi,wVallohi,wValhilo,wVallolo; // value in beta*R^H*v
   double acchihi,acclohi,acchilo,acclolo;
 
   __syncthreads();
   Rwidxhihi = Rimhihi[Ridx];
   Rwidxlohi = Rimlohi[Ridx];
   Rwidxhilo = Rimhilo[Ridx];
   Rwidxlolo = Rimlolo[Ridx];

   __syncthreads();
   vValhihi = vimhihi[rowidx];
   vVallohi = vimlohi[rowidx];
   vValhilo = vimhilo[rowidx];
   vVallolo = vimlolo[rowidx];
   wValhihi = shwrehihi[colidx];  // value in beta*R^H*v
   wVallohi = shwrelohi[colidx];
   wValhilo = shwrehilo[colidx];
   wVallolo = shwrelolo[colidx];
   __syncthreads();
   qdg_mul(vValhihi,vVallohi,vValhilo,vVallolo,
           wValhihi,wVallohi,wValhilo,wVallolo,
           &acchihi,&acclohi,&acchilo,&acclolo);
   qdg_dec(&Rwidxhihi,&Rwidxlohi,&Rwidxhilo,&Rwidxlolo,
              acchihi,   acclohi,   acchilo,   acclolo);

   vValhihi = vrehihi[rowidx];    // value in Householder vector
   vVallohi = vrelohi[rowidx];
   vValhilo = vrehilo[rowidx];
   vVallolo = vrelolo[rowidx];
   wValhihi = shwimhihi[colidx];
   wVallohi = shwimlohi[colidx];
   wValhilo = shwimhilo[colidx];
   wVallolo = shwimlolo[colidx];
   __syncthreads();
   qdg_mul(vValhihi,vVallohi,vValhilo,vVallolo,
           wValhihi,wVallohi,wValhilo,wVallolo,
           &acchihi,&acclohi,&acchilo,&acclolo);
   qdg_inc(&Rwidxhihi,&Rwidxlohi,&Rwidxhilo,&Rwidxlolo,
              acchihi,   acclohi,   acchilo,   acclolo);

   __syncthreads();
   if(widx < bound)                    // if() takes care of padding
   {
      Rimhihi[Ridx] = Rwidxhihi;
      Rimlohi[Ridx] = Rwidxlohi;
      Rimhilo[Ridx] = Rwidxhilo;
      Rimlolo[Ridx] = Rwidxlolo;
   }
}

__global__ void dbl4_beta_times_V
 ( int nrows, int szt,
   double *Bhihi, double *Blohi, double *Bhilo, double *Blolo,
   double *Vhihi, double *Vlohi, double *Vhilo, double *Vlolo,
   double *Whihi, double *Wlohi, double *Whilo, double *Wlolo )
{
   const int bdx = blockIdx.x;        // index of block
   const int tdx = threadIdx.x;       // index of thread in block
   const int idx = bdx*szt + tdx;     // thread tdx computes W[idx]
   double resulthihi,resultlohi,resulthilo,resultlolo;

   __shared__ double shvhihi[qd_shmemsize]; // to store a slice of V
   __shared__ double shvlohi[qd_shmemsize];
   __shared__ double shvhilo[qd_shmemsize];
   __shared__ double shvlolo[qd_shmemsize];

   shvhihi[tdx] = Vhihi[idx]; // thread tdx loads the data at the global index
   shvlohi[tdx] = Vlohi[idx];
   shvhilo[tdx] = Vhilo[idx];
   shvlolo[tdx] = Vlolo[idx];

   // result = -B[0]*shv[tdx];
   __syncthreads();
   qdg_mul(     -Bhihi[0],   -Blohi[0],   -Bhilo[0],   -Blolo[0],
               shvhihi[tdx],shvlohi[tdx],shvhilo[tdx],shvlolo[tdx],
           &resulthihi, &resultlohi, &resulthilo, &resultlolo);

   __syncthreads();
   if(idx < nrows)
   {
      Whihi[idx] = resulthihi;
      Wlohi[idx] = resultlohi;
      Whilo[idx] = resulthilo;
      Wlolo[idx] = resultlolo;
   }
}

__global__ void cmplx4_beta_times_V
 ( int nrows, int szt,
   double *Bhihi, double *Blohi, double *Bhilo, double *Blolo,
   double *Vrehihi, double *Vrelohi, double *Vrehilo, double *Vrelolo,
   double *Vimhihi, double *Vimlohi, double *Vimhilo, double *Vimlolo,
   double *Wrehihi, double *Wrelohi, double *Wrehilo, double *Wrelolo,
   double *Wimhihi, double *Wimlohi, double *Wimhilo, double *Wimlolo )
{
   const int bdx = blockIdx.x;        // index of block
   const int tdx = threadIdx.x;       // index of thread in block
   const int idx = bdx*szt + tdx;     // thread tdx computes W[idx]
   double resulthihi,resultlohi,resulthilo,resultlolo;
   const double minBhihi = -Bhihi[0];
   const double minBlohi = -Blohi[0];
   const double minBhilo = -Bhilo[0];
   const double minBlolo = -Blolo[0];

   __shared__ double shvrehihi[cqd_shmemsize]; // to store a slice of V
   __shared__ double shvrelohi[cqd_shmemsize];
   __shared__ double shvrehilo[cqd_shmemsize];
   __shared__ double shvrelolo[cqd_shmemsize];
   __shared__ double shvimhihi[cqd_shmemsize];
   __shared__ double shvimlohi[cqd_shmemsize];
   __shared__ double shvimhilo[cqd_shmemsize];
   __shared__ double shvimlolo[cqd_shmemsize];

   shvrehihi[tdx] = Vrehihi[idx]; // thread tdx loads data at the global index
   shvrelohi[tdx] = Vrelohi[idx];
   shvrehilo[tdx] = Vrehilo[idx];
   shvrelolo[tdx] = Vrelolo[idx];

   // resultre = -B[0]*shvre[tdx];
   __syncthreads();
   qdg_mul(   minBhihi,      minBlohi,      minBhilo,      minBlolo,
             shvrehihi[tdx],shvrelohi[tdx],shvrehilo[tdx],shvrelolo[tdx],
           &resulthihi,   &resultlohi,   &resulthilo,   &resultlolo);

   __syncthreads();
   if(idx < nrows)
   {
      Wrehihi[idx] = resulthihi;
      Wrelohi[idx] = resultlohi;
      Wrehilo[idx] = resulthilo;
      Wrelolo[idx] = resultlolo;
   }
   // resultim = -B[0]*shvim[tdx];
   __syncthreads();
   shvimhihi[tdx] = Vimhihi[idx];
   shvimlohi[tdx] = Vimlohi[idx];
   shvimhilo[tdx] = Vimhilo[idx];
   shvimlolo[tdx] = Vimlolo[idx];
   __syncthreads();
   qdg_mul(   minBhihi,      minBlohi,      minBhilo,      minBlolo,
             shvimhihi[tdx],shvimlohi[tdx],shvimhilo[tdx],shvimlolo[tdx],
           &resulthihi,   &resultlohi,   &resulthilo,   &resultlolo);

   __syncthreads();
   if(idx < nrows)
   {
      Wimhihi[idx] = resulthihi;
      Wimlohi[idx] = resultlohi;
      Wimhilo[idx] = resulthilo;
      Wimlolo[idx] = resultlolo;
   }
}

__global__ void dbl4_initialize_WYT
 ( int dim, int szt,
   double *Vhihi, double *Vlohi, double *Vhilo, double *Vlolo,
   double *Whihi, double *Wlohi, double *Whilo, double *Wlolo,
   double *WYThihi, double *WYTlohi, double *WYThilo, double *WYTlolo )
{
   const int bdx = blockIdx.x;        // index of block
   const int tdx = threadIdx.x;       // index of thread in block
   const int idx = bdx*szt + tdx;     // global index of the thread
   const int row = idx / dim;         // row index in WYT
   const int col = idx % dim;         // column index in WYT

   const double Vvalhihi = Vhihi[col];
   const double Vvallohi = Vlohi[col];
   const double Vvalhilo = Vhilo[col];
   const double Vvallolo = Vlolo[col];
   const double Wvalhihi = Whihi[row];
   const double Wvallohi = Wlohi[row];
   const double Wvalhilo = Whilo[row];
   const double Wvallolo = Wlolo[row];
   // const double result = Vval*Wval;
   double resulthihi,resultlohi,resulthilo,resultlolo; 

   __syncthreads();
   qdg_mul(   Vvalhihi,   Vvallohi,   Vvalhilo,   Vvallolo,
              Wvalhihi,   Wvallohi,   Wvalhilo,   Wvallolo,
           &resulthihi,&resultlohi,&resulthilo,&resultlolo);

   __syncthreads();
   if(idx < dim*dim)
   {
      WYThihi[idx] = resulthihi;
      WYTlohi[idx] = resultlohi;
      WYThilo[idx] = resulthilo;
      WYTlolo[idx] = resultlolo;
   }
}

__global__ void cmplx4_initialize_WYH
 ( int dim, int szt,
   double *Vrehihi, double *Vrelohi, double *Vrehilo, double *Vrelolo,
   double *Vimhihi, double *Vimlohi, double *Vimhilo, double *Vimlolo,
   double *Wrehihi, double *Wrelohi, double *Wrehilo, double *Wrelolo,
   double *Wimhihi, double *Wimlohi, double *Wimhilo, double *Wimlolo,
   double *WYHrehihi, double *WYHrelohi,
   double *WYHrehilo, double *WYHrelolo,
   double *WYHimhihi, double *WYHimlohi,
   double *WYHimhilo, double *WYHimlolo )
{
   const int bdx = blockIdx.x;        // index of block
   const int tdx = threadIdx.x;       // index of thread in block
   const int idx = bdx*szt + tdx;     // global index of the thread
   const int row = idx / dim;         // row index in WYH
   const int col = idx % dim;         // column index in WYH

   const double Vvalrehihi = Vrehihi[col];
   const double Vvalrelohi = Vrelohi[col];
   const double Vvalrehilo = Vrehilo[col];
   const double Vvalrelolo = Vrelolo[col];
   const double Vvalimhihi = Vimhihi[col];
   const double Vvalimlohi = Vimlohi[col];
   const double Vvalimhilo = Vimhilo[col];
   const double Vvalimlolo = Vimlolo[col];
   const double Wvalrehihi = Wrehihi[row];
   const double Wvalrelohi = Wrelohi[row];
   const double Wvalrehilo = Wrehilo[row];
   const double Wvalrelolo = Wrelolo[row];
   const double Wvalimhihi = Wimhihi[row];
   const double Wvalimlohi = Wimlohi[row];
   const double Wvalimhilo = Wimhilo[row];
   const double Wvalimlolo = Wimlolo[row];
   // const double result = Vval*Wval;
   double resulthihi,resultlohi,resulthilo,resultlolo;
   double acchihi,acclohi,acchilo,acclolo;

   // take the Hermitian transpose of V
   __syncthreads();
   qdg_mul( Vvalrehihi, Vvalrelohi, Vvalrehilo, Vvalrelolo,
            Wvalrehihi, Wvalrelohi, Wvalrehilo, Wvalrelolo,
           &resulthihi,&resultlohi,&resulthilo,&resultlolo);
   qdg_mul(Vvalimhihi,Vvalimlohi,Vvalimhilo,Vvalimlolo,
           Wvalimhihi,Wvalimlohi,Wvalimhilo,Wvalimlolo,
             &acchihi,  &acclohi,  &acchilo,  &acclolo);
   qdg_inc(&resulthihi,&resultlohi,&resulthilo,&resultlolo,
               acchihi,    acclohi,    acchilo,    acclolo);
   __syncthreads();
   if(idx < dim*dim)
   {
      WYHrehihi[idx] = resulthihi;
      WYHrelohi[idx] = resultlohi;
      WYHrehilo[idx] = resulthilo;
      WYHrelolo[idx] = resultlolo;
   }
   __syncthreads();
   qdg_mul( Vvalrehihi, Vvalrelohi, Vvalrehilo, Vvalrelolo,
            Wvalimhihi, Wvalimlohi, Wvalimhilo, Wvalimlolo,
           &resulthihi,&resultlohi,&resulthilo,&resultlolo);
   qdg_mul(Vvalimhihi,Vvalimlohi,Vvalimhilo,Vvalimlolo,
           Wvalrehihi,Wvalrelohi,Wvalrehilo,Wvalrelolo,
             &acchihi,  &acclohi,  &acchilo,  &acclolo);
   qdg_dec(&resulthihi,&resultlohi,&resulthilo,&resultlolo,
               acchihi,    acclohi,    acchilo,    acclolo);

   __syncthreads();
   if(idx < dim*dim)
   {
      WYHimhihi[idx] = resulthihi;
      WYHimlohi[idx] = resultlohi;
      WYHimhilo[idx] = resulthilo;
      WYHimlolo[idx] = resultlolo;
   }
}

__global__ void dbl4_update_WYT
 ( int dim, int szt,
   double *Vhihi, double *Vlohi, double *Vhilo, double *Vlolo,
   double *Whihi, double *Wlohi, double *Whilo, double *Wlolo,
   double *WYThihi, double *WYTlohi, double *WYThilo, double *WYTlolo )
{
   const int bdx = blockIdx.x;        // index of block
   const int tdx = threadIdx.x;       // index of thread in block
   const int idx = bdx*szt + tdx;     // global index of the thread
   const int row = idx / dim;         // row index in WYT
   const int col = idx % dim;         // column index in WYT

   const double Vvalhihi = Vhihi[col];
   const double Vvallohi = Vlohi[col];
   const double Vvalhilo = Vhilo[col];
   const double Vvallolo = Vlolo[col];
   __syncthreads();
   const double Wvalhihi = Whihi[row];
   const double Wvallohi = Wlohi[row];
   const double Wvalhilo = Whilo[row];
   const double Wvallolo = Wlolo[row];
   __syncthreads();
   double resulthihi = WYThihi[idx];
   double resultlohi = WYTlohi[idx];
   double resulthilo = WYThilo[idx];
   double resultlolo = WYTlolo[idx];
   double acchihi,acclohi,acchilo,acclolo;

   // result = result + Vval*Wval;

   __syncthreads();
   qdg_mul(Vvalhihi,Vvallohi,Vvalhilo,Vvallolo,
           Wvalhihi,Wvallohi,Wvalhilo,Wvallolo,
           &acchihi,&acclohi,&acchilo,&acclolo);
   qdg_inc(&resulthihi,&resultlohi,&resulthilo,&resultlolo,
               acchihi,    acclohi,    acchilo,    acclolo);
   
   __syncthreads();
   if(idx < dim*dim)
   {
      WYThihi[idx] = resulthihi;
      WYTlohi[idx] = resultlohi;
      WYThilo[idx] = resulthilo;
      WYTlolo[idx] = resultlolo;
   }
}

__global__ void cmplx4_update_WYH
 ( int dim, int szt,
   double *Vrehihi, double *Vrelohi, double *Vrehilo, double *Vrelolo,
   double *Vimhihi, double *Vimlohi, double *Vimhilo, double *Vimlolo,
   double *Wrehihi, double *Wrelohi, double *Wrehilo, double *Wrelolo,
   double *Wimhihi, double *Wimlohi, double *Wimhilo, double *Wimlolo,
   double *WYHrehihi, double *WYHrelohi,
   double *WYHrehilo, double *WYHrelolo,
   double *WYHimhihi, double *WYHimlohi,
   double *WYHimhilo, double *WYHimlolo )
{
   const int bdx = blockIdx.x;        // index of block
   const int tdx = threadIdx.x;       // index of thread in block
   const int idx = bdx*szt + tdx;     // global index of the thread
   const int row = idx / dim;         // row index in WYH
   const int col = idx % dim;         // column index in WYH

   const double Vvalrehihi = Vrehihi[col];
   const double Vvalrelohi = Vrelohi[col];
   const double Vvalrehilo = Vrehilo[col];
   const double Vvalrelolo = Vrelolo[col];
   const double Vvalimhihi = Vimhihi[col];
   const double Vvalimlohi = Vimlohi[col];
   const double Vvalimhilo = Vimhilo[col];
   const double Vvalimlolo = Vimlolo[col];
   __syncthreads();
   const double Wvalrehihi = Wrehihi[row];
   const double Wvalrelohi = Wrelohi[row];
   const double Wvalrehilo = Wrehilo[row];
   const double Wvalrelolo = Wrelolo[row];
   const double Wvalimhihi = Wimhihi[row];
   const double Wvalimlohi = Wimlohi[row];
   const double Wvalimhilo = Wimhilo[row];
   const double Wvalimlolo = Wimlolo[row];
   // const double result = Vval*Wval;
   __syncthreads();
   double resulthihi = WYHrehihi[idx];
   double resultlohi = WYHrelohi[idx];
   double resulthilo = WYHrehilo[idx];
   double resultlolo = WYHrelolo[idx];
   double acchihi,acclohi,acchilo,acclolo;

   // take the Hermitian transpose of V
   __syncthreads();
   qdg_mul(Vvalrehihi,Vvalrelohi,Vvalrehilo,Vvalrelolo,
           Wvalrehihi,Wvalrelohi,Wvalrehilo,Wvalrelolo,
             &acchihi , &acclohi,  &acchilo,  &acclolo);
   qdg_inc(&resulthihi,&resultlohi,&resulthilo,&resultlolo,
               acchihi,    acclohi,    acchilo,    acclolo);
   qdg_mul(Vvalimhihi,Vvalimlohi,Vvalimhilo,Vvalimlolo,
           Wvalimhihi,Wvalimlohi,Wvalimhilo,Wvalimlolo,
             &acchihi,  &acclohi,  &acchilo,  &acclolo);
   qdg_inc(&resulthihi,&resultlohi,&resulthilo,&resultlolo,
               acchihi,    acclohi,    acchilo,    acclolo);
   __syncthreads();
   if(idx < dim*dim)
   {
      WYHrehihi[idx] = resulthihi;
      WYHrelohi[idx] = resultlohi;
      WYHrehilo[idx] = resulthilo;
      WYHrelolo[idx] = resultlolo;
   }
   __syncthreads();
   resulthihi = WYHimhihi[idx];
   resultlohi = WYHimlohi[idx];
   resulthilo = WYHimhilo[idx];
   resultlolo = WYHimlolo[idx];

   qdg_mul(Vvalrehihi,Vvalrelohi,Vvalrehilo,Vvalrelolo,
           Wvalimhihi,Wvalimlohi,Wvalimhilo,Wvalimlolo,
             &acchihi,  &acclohi,  &acchilo,  &acclolo);
   qdg_inc(&resulthihi,&resultlohi,&resulthilo,&resultlolo,
               acchihi,    acclohi,    acchilo,    acclolo);
   qdg_mul(Vvalimhihi,Vvalimlohi,Vvalimhilo,Vvalimlolo,
           Wvalrehihi,Wvalrelohi,Wvalrehilo,Wvalrelolo,
             &acchihi,  &acclohi,  &acchilo,  &acclolo);
   qdg_dec(&resulthihi,&resultlohi,&resulthilo,&resultlolo,
               acchihi,    acclohi,    acchilo,    acclolo);

   __syncthreads();
   if(idx < dim*dim)
   {
      WYHimhihi[idx] = resulthihi;
      WYHimlohi[idx] = resultlohi;
      WYHimhilo[idx] = resulthilo;
      WYHimlolo[idx] = resultlolo;
   }
}

__global__ void dbl4_beta_next_W
 ( int nrows, int szt,
   double *Bhihi, double *Blohi, double *Bhilo, double *Blolo,
   double *Vhihi, double *Vlohi, double *Vhilo, double *Vlolo,
   double *Whihi, double *Wlohi, double *Whilo, double *Wlolo,
   double *WYThihi, double *WYTlohi, double *WYThilo, double *WYTlolo )
{
   const int bdx = blockIdx.x;        // index of block
   const int tdx = threadIdx.x;       // index of thread in block
   const int idx = bdx*szt + tdx;     // global index of the thread
   const int WYToff = idx*nrows;      // start of idx row in YWT
   const double mybetahihi = Bhihi[0];
   const double mybetalohi = Blohi[0];
   const double mybetahilo = Bhilo[0];
   const double mybetalolo = Blolo[0];
   int vdx,ydx;
   double resulthihi,resultlohi,resulthilo,resultlolo;
   double WYTvalhihi,WYTvallohi,WYTvalhilo,WYTvallolo;
   double Vvalhihi,Vvallohi,Vvalhilo,Vvallolo;
   double acchihi,acclohi,acchilo,acclolo;

   __shared__ double shVhihi[qd_shmemsize];   // to store a slice of V
   __shared__ double shVlohi[qd_shmemsize];
   __shared__ double shVhilo[qd_shmemsize];
   __shared__ double shVlolo[qd_shmemsize];

   shVhihi[tdx] = Vhihi[idx]; // thread tdx loads the data at the global index
   shVlohi[tdx] = Vlohi[idx];
   shVhilo[tdx] = Vhilo[idx];
   shVlolo[tdx] = Vlolo[idx];

   __syncthreads();
   resulthihi = shVhihi[tdx]; // thread tdx computes the value at index idx
   resultlohi = shVlohi[tdx];
   resulthilo = shVhilo[tdx];
   resultlolo = shVlolo[tdx];

   for(int i=0; i<nrows/szt; i++)
   {
      vdx = i*szt + tdx;                 // index in V and in YWT
      __syncthreads();
      shVhihi[tdx] = Vhihi[vdx];         // threads load next szt values
      shVlohi[tdx] = Vlohi[vdx];
      shVhilo[tdx] = Vhilo[vdx];
      shVlolo[tdx] = Vlolo[vdx];

      __syncthreads();
      for(int j=0; j<szt; j++)           // multiply szt values with YWT
      {
         ydx = WYToff + i*szt + j;       // WYT is stored row by row
         __syncthreads();
         WYTvalhihi = WYThihi[ydx];
         WYTvallohi = WYTlohi[ydx];
         WYTvalhilo = WYThilo[ydx];
         WYTvallolo = WYTlolo[ydx];
         __syncthreads();
         Vvalhihi = shVhihi[j];
         Vvallohi = shVlohi[j];
         Vvalhilo = shVhilo[j];
         Vvallolo = shVlolo[j];
         // result = result + WYTval*Vvalue;
         __syncthreads();
         qdg_mul(  Vvalhihi,  Vvallohi,  Vvalhilo,  Vvallolo,
                 WYTvalhihi,WYTvallohi,WYTvalhilo,WYTvallolo,
                   &acchihi,  &acclohi,  &acchilo,  &acclolo);
         qdg_inc(&resulthihi,&resultlohi,&resulthilo,&resultlolo,
                     acchihi,    acclohi,    acchilo,    acclolo);
      }
      __syncthreads();
   }
   int quot = nrows/szt;
   int rest = nrows - quot*szt;          // remainder to compute

   vdx = quot*szt + tdx;                 // next index to compute
   __syncthreads();
   shVhihi[tdx] = Vhihi[vdx];
   shVlohi[tdx] = Vlohi[vdx];
   shVhilo[tdx] = Vhilo[vdx];
   shVlolo[tdx] = Vlolo[vdx];

   for(int j=0; j<rest; j++)            // rest < szt prevents overflow
   {
      __syncthreads();
      ydx = WYToff + quot*szt + j;
      WYTvalhihi = WYThihi[ydx];
      WYTvallohi = WYTlohi[ydx];
      WYTvalhilo = WYThilo[ydx];
      WYTvallolo = WYTlolo[ydx];
      __syncthreads();
      Vvalhihi = shVhihi[j];
      Vvallohi = shVlohi[j];
      Vvalhilo = shVhilo[j];
      Vvallolo = shVlolo[j];
      // result = result + WYTval*Vvalue;
      __syncthreads();
      qdg_mul(  Vvalhihi,  Vvallohi,  Vvalhilo,  Vvallolo,
              WYTvalhihi,WYTvallohi,WYTvalhilo,WYTvallolo,
                &acchihi,  &acclohi,  &acchilo,  &acclolo);
      qdg_inc(&resulthihi,&resultlohi,&resulthilo,&resultlolo,
                  acchihi,    acclohi,    acchilo,    acclolo);
   }
   // result = -mybeta*result;
   __syncthreads();
   qdg_mul(-mybetahihi,-mybetalohi,-mybetahilo,-mybetalolo,
            resulthihi, resultlohi, resulthilo, resultlolo,
              &acchihi,   &acclohi,   &acchilo,   &acclolo);

   __syncthreads();
   if(idx < nrows) 
   {
      Whihi[idx] = acchihi;
      Wlohi[idx] = acclohi;
      Whilo[idx] = acchilo;
      Wlolo[idx] = acclolo;
   }
}

__global__ void cmplx4_beta_next_W
 ( int nrows, int szt,
   double *Bhihi, double *Blohi, double *Bhilo, double *Blolo,
   double *Vrehihi, double *Vrelohi, double *Vrehilo, double *Vrelolo,
   double *Vimhihi, double *Vimlohi, double *Vimhilo, double *Vimlolo,
   double *Wrehihi, double *Wrelohi, double *Wrehilo, double *Wrelolo,
   double *Wimhihi, double *Wimlohi, double *Wimhilo, double *Wimlolo,
   double *WYHrehihi, double *WYHrelohi,
   double *WYHrehilo, double *WYHrelolo,
   double *WYHimhihi, double *WYHimlohi,
   double *WYHimhilo, double *WYHimlolo )
{
   const int bdx = blockIdx.x;        // index of block
   const int tdx = threadIdx.x;       // index of thread in block
   const int idx = bdx*szt + tdx;     // global index of the thread
   const int WYHoff = idx*nrows;      // start of idx row in YWH
   const double mybetahihi = Bhihi[0];
   const double mybetalohi = Blohi[0];
   const double mybetahilo = Bhilo[0];
   const double mybetalolo = Blolo[0];
   // __syncthreads();
   int vdx,ydx;
   double resultrehihi,resultrelohi,resultrehilo,resultrelolo;
   double resultimhihi,resultimlohi,resultimhilo,resultimlolo;
   double WYHvalhihi,WYHvallohi,WYHvalhilo,WYHvallolo;
   double Vvalhihi,Vvallohi,Vvalhilo,Vvallolo;
   double acchihi,acclohi,acchilo,acclolo;

   __shared__ double shVrehihi[cqd_shmemsize];   // to store a slice of V
   __shared__ double shVrelohi[cqd_shmemsize];
   __shared__ double shVrehilo[cqd_shmemsize];
   __shared__ double shVrelolo[cqd_shmemsize];
   __shared__ double shVimhihi[cqd_shmemsize];
   __shared__ double shVimlohi[cqd_shmemsize];
   __shared__ double shVimhilo[cqd_shmemsize];
   __shared__ double shVimlolo[cqd_shmemsize];

   shVrehihi[tdx] = Vrehihi[idx]; // thread tdx loads data at the global index
   shVrelohi[tdx] = Vrelohi[idx];
   shVrehilo[tdx] = Vrehilo[idx];
   shVrelolo[tdx] = Vrelolo[idx];
   shVimhihi[tdx] = Vimhihi[idx];
   shVimlohi[tdx] = Vimlohi[idx];
   shVimhilo[tdx] = Vimhilo[idx];
   shVimlolo[tdx] = Vimlolo[idx];

   __syncthreads();
   resultrehihi = shVrehihi[tdx]; // thread tdx computes the value at idx
   resultrelohi = shVrelohi[tdx];
   resultrehilo = shVrehilo[tdx];
   resultrelolo = shVrelolo[tdx];
   resultimhihi = shVimhihi[tdx];
   resultimlohi = shVimlohi[tdx];
   resultimhilo = shVimhilo[tdx];
   resultimlolo = shVimlolo[tdx];

   for(int i=0; i<nrows/szt; i++)
   {
      vdx = i*szt + tdx;                 // index in V and in YWT
      __syncthreads();
      shVrehihi[tdx] = Vrehihi[vdx];     // threads load next szt values
      shVrelohi[tdx] = Vrelohi[vdx];
      shVrehilo[tdx] = Vrehilo[vdx];
      shVrelolo[tdx] = Vrelolo[vdx];
      shVimhihi[tdx] = Vimhihi[vdx];
      shVimlohi[tdx] = Vimlohi[vdx];
      shVimhilo[tdx] = Vimhilo[vdx];
      shVimlolo[tdx] = Vimlolo[vdx];

      __syncthreads();
      for(int j=0; j<szt; j++)           // multiply szt values with YWT
      {
         ydx = WYHoff + i*szt + j;       // WYH is stored row by row
         // take the Hermitian transpose of V

         Vvalhihi = shVrehihi[j];
         Vvallohi = shVrelohi[j];
         Vvalhilo = shVrehilo[j];
         Vvallolo = shVrelolo[j];
         WYHvalhihi = WYHrehihi[ydx];
         WYHvallohi = WYHrelohi[ydx];
         WYHvalhilo = WYHrehilo[ydx];
         WYHvallolo = WYHrelolo[ydx];
         qdg_mul(  Vvalhihi,  Vvallohi,  Vvalhilo,  Vvallolo,
                 WYHvalhihi,WYHvallohi,WYHvalhilo,WYHvallolo,
                   &acchihi,  &acclohi,  &acchilo,  &acclolo);
         qdg_inc(&resultrehihi,&resultrelohi,&resultrehilo,&resultrelolo,
                       acchihi,      acclohi,      acchilo,      acclolo);

         Vvalhihi = shVimhihi[j];
         Vvallohi = shVimlohi[j];
         Vvalhilo = shVimhilo[j];
         Vvallolo = shVimlolo[j];
         WYHvalhihi = WYHimhihi[ydx];
         WYHvallohi = WYHimlohi[ydx];
         WYHvalhilo = WYHimhilo[ydx];
         WYHvallolo = WYHimlolo[ydx];
         qdg_mul  (Vvalhihi,  Vvallohi,  Vvalhilo,  Vvallolo,
                 WYHvalhihi,WYHvallohi,WYHvalhilo,WYHvallolo,
                   &acchihi,  &acclohi,  &acchilo,  &acclolo);
         qdg_dec(&resultrehihi,&resultrelohi,&resultrehilo,&resultrelolo,
                       acchihi,      acclohi,      acchilo,      acclolo);

         // have already imaginary values for V
         WYHvalhihi = WYHrehihi[ydx];
         WYHvallohi = WYHrelohi[ydx];
         WYHvalhilo = WYHrehilo[ydx];
         WYHvallolo = WYHrelolo[ydx];
         qdg_mul(  Vvalhihi,  Vvallohi,  Vvalhilo,  Vvallolo,
                 WYHvalhihi,WYHvallohi,WYHvalhilo,WYHvallolo,
                   &acchihi,  &acclohi,  &acchilo,  &acclolo);
         qdg_inc(&resultimhihi,&resultimlohi,&resultimhilo,&resultimlolo,
                       acchihi,      acclohi,      acchilo,      acclolo);

         Vvalhihi = shVrehihi[j];
         Vvallohi = shVrelohi[j];
         Vvalhilo = shVrehilo[j];
         Vvallolo = shVrelolo[j];
         WYHvalhihi = WYHimhihi[ydx];
         WYHvallohi = WYHimlohi[ydx];
         WYHvalhilo = WYHimhilo[ydx];
         WYHvallolo = WYHimlolo[ydx];
         qdg_mul(  Vvalhihi,  Vvallohi,  Vvalhilo,  Vvallolo,
                 WYHvalhihi,WYHvallohi,WYHvalhilo,WYHvallolo,
                   &acchihi,  &acclohi,  &acchilo,  &acclolo);
         qdg_inc(&resultimhihi,&resultimlohi,&resultimhilo,&resultimlolo,
                       acchihi,      acclohi,      acchilo,      acclolo);
      }
      __syncthreads();
   }
   int quot = nrows/szt;
   int rest = nrows - quot*szt;          // remainder to compute

   vdx = quot*szt + tdx;                 // next index to compute
   __syncthreads();
   shVrehihi[tdx] = Vrehihi[vdx];
   shVrelohi[tdx] = Vrelohi[vdx];
   shVrehilo[tdx] = Vrehilo[vdx];
   shVrelolo[tdx] = Vrelolo[vdx];
   shVimhihi[tdx] = Vimhihi[vdx];
   shVimlohi[tdx] = Vimlohi[vdx];
   shVimhilo[tdx] = Vimhilo[vdx];
   shVimlolo[tdx] = Vimlolo[vdx];

   for(int j=0; j<rest; j++)            // rest < szt prevents overflow
   {
      // result = result + WYTval*Vvalue;
      // take the Hermitian transpose of V

      Vvalhihi = shVrehihi[j];
      Vvallohi = shVrelohi[j];
      Vvalhilo = shVrehilo[j];
      Vvallolo = shVrelolo[j];
      WYHvalhihi = WYHrehihi[ydx];
      WYHvallohi = WYHrelohi[ydx];
      WYHvalhilo = WYHrehilo[ydx];
      WYHvallolo = WYHrelolo[ydx];
      qdg_mul(  Vvalhihi,  Vvallohi,  Vvalhilo,  Vvallolo,
              WYHvalhihi,WYHvallohi,WYHvalhilo,WYHvallolo,
                &acchihi,  &acclohi,  &acchilo,  &acclolo);
      qdg_inc(&resultrehihi,&resultrelohi,&resultrehilo,&resultrelolo,
                    acchihi,      acclohi,      acchilo,      acclolo);

      Vvalhihi = shVimhihi[j];
      Vvallohi = shVimlohi[j];
      Vvalhilo = shVimhilo[j];
      Vvallolo = shVimlolo[j];
      WYHvalhihi = WYHimhihi[ydx];
      WYHvallohi = WYHimlohi[ydx];
      WYHvalhilo = WYHimhilo[ydx];
      WYHvallolo = WYHimlolo[ydx];
      qdg_mul(  Vvalhihi,  Vvallohi,  Vvalhilo,  Vvallolo,
              WYHvalhihi,WYHvallohi,WYHvalhilo,WYHvallolo,
                &acchihi,  &acclohi,  &acchilo,  &acclolo);
      qdg_dec(&resultrehihi,&resultrelohi,&resultrehilo,&resultrelolo,
                    acchihi,      acclohi,      acchilo,      acclolo);

      // already have the imaginary values for V
      WYHvalhihi = WYHrehihi[ydx];
      WYHvallohi = WYHrelohi[ydx];
      WYHvalhilo = WYHrehilo[ydx];
      WYHvallolo = WYHrelolo[ydx];
      qdg_mul(  Vvalhihi,  Vvallohi,  Vvalhilo,  Vvallolo,
              WYHvalhihi,WYHvallohi,WYHvalhilo,WYHvallolo,
                &acchihi,  &acclohi,  &acchilo,  &acclolo);
      qdg_inc(&resultimhihi,&resultimlohi,&resultimhilo,&resultimlolo,
                    acchihi,      acclohi,      acchilo,      acclolo);

      Vvalhihi = shVrehihi[j];
      Vvallohi = shVrelohi[j];
      Vvalhilo = shVrehilo[j];
      Vvallolo = shVrelolo[j];
      WYHvalhihi = WYHimhihi[ydx];
      WYHvallohi = WYHimlohi[ydx];
      WYHvalhilo = WYHimhilo[ydx];
      WYHvallolo = WYHimlolo[ydx];
      qdg_mul(  Vvalhihi,  Vvallohi,  Vvalhilo,  Vvallolo,
              WYHvalhihi,WYHvallohi,WYHvalhilo,WYHvallolo,
                &acchihi,  &acclohi,  &acchilo,  &acclolo);
      qdg_inc(&resultimhihi,&resultimlohi,&resultimhilo,&resultimlolo,
                    acchihi,      acclohi,      acchilo,      acclolo);
   }
   // result = -mybeta*result;
   __syncthreads();
   qdg_mul( -mybetahihi, -mybetalohi, -mybetahilo, -mybetalolo,
           resultrehihi,resultrelohi,resultrehilo,resultrelolo,
               &acchihi,    &acclohi,    &acchilo,    &acclolo);

   // __syncthreads();
   if(idx < nrows) 
   {
      Wrehihi[idx] = acchihi;
      Wrelohi[idx] = acclohi;
      Wrehilo[idx] = acchilo;
      Wrelolo[idx] = acclolo;
   }
   __syncthreads();
   qdg_mul( -mybetahihi, -mybetalohi, -mybetahilo, -mybetalolo,
           resultimhihi,resultimlohi,resultimhilo,resultimlolo,
               &acchihi,    &acclohi,    &acchilo,    &acclolo);

   __syncthreads();
   if(idx < nrows) 
   {
      Wimhihi[idx] = acchihi;
      Wimlohi[idx] = acclohi;
      Wimhilo[idx] = acchilo;
      Wimlolo[idx] = acclolo;
   }
}

__global__ void dbl4_small_WYT
 ( int nrows, int szt,
   double *Whihi, double *Wlohi, double *Whilo, double *Wlolo,
   double *Vhihi, double *Vlohi, double *Vhilo, double *Vlolo,
   double *WYThihi, double *WYTlohi, double *WYThilo, double *WYTlolo )
{
   const int bdx = blockIdx.x;           // index of block
   const int tdx = threadIdx.x;          // index of thread in block
   const int offset = bdx*szt + tdx;     // offset in result
   const int row = offset / nrows;
   const int col = offset % nrows;       // thread 0 computes WYT[row][col]

   double resulthihi = 0.0;
   double resultlohi = 0.0;
   double resulthilo = 0.0;
   double resultlolo = 0.0;
   double ahihi,alohi,ahilo,alolo;
   double bhihi,blohi,bhilo,blolo;
   double chihi,clohi,chilo,clolo;

   for(int k=0; k<szt; k++)
   {
      __syncthreads();
      ahihi = Whihi[k*nrows + row];   // if(nrows == szt) then row = bdx
      alohi = Wlohi[k*nrows + row];
      ahilo = Whilo[k*nrows + row];
      alolo = Wlolo[k*nrows + row];
      __syncthreads();
      bhihi = Vhihi[k*nrows + col];   // if(nrows == szt) then col = tdx
      blohi = Vlohi[k*nrows + col]; 
      bhilo = Vhilo[k*nrows + col]; 
      blolo = Vlolo[k*nrows + col]; 
      // result = result + a*b;
      __syncthreads();
      qdg_mul( ahihi, alohi, ahilo, alolo,
               bhihi, blohi, bhilo, blolo,
              &chihi,&clohi,&chilo,&clolo);
      qdg_inc(&resulthihi,&resultlohi,&resulthilo,&resultlolo,
                    chihi,      clohi,      chilo,      clolo);
   }
   __syncthreads();
   WYThihi[offset] = resulthihi;
   WYTlohi[offset] = resultlohi;
   WYThilo[offset] = resulthilo;
   WYTlolo[offset] = resultlolo;
}

__global__ void cmplx4_small_WYH
 ( int nrows, int szt,
   double *Wrehihi, double *Wrelohi, double *Wrehilo, double *Wrelolo,
   double *Wimhihi, double *Wimlohi, double *Wimhilo, double *Wimlolo,
   double *Yrehihi, double *Yrelohi, double *Yrehilo, double *Yrelolo,
   double *Yimhihi, double *Yimlohi, double *Yimhilo, double *Yimlolo,
   double *WYTrehihi, double *WYTrelohi,
   double *WYTrehilo, double *WYTrelolo,
   double *WYTimhihi, double *WYTimlohi,
   double *WYTimhilo, double *WYTimlolo )
{
   const int bdx = blockIdx.x;           // index of block
   const int tdx = threadIdx.x;          // index of thread in block
   const int offset = bdx*szt + tdx;     // offset in result
   const int row = offset / nrows;
   const int col = offset % nrows;       // thread 0 computes WYT[row][col]

   double resultrehihi = 0.0;
   double resultrelohi = 0.0;
   double resultrehilo = 0.0;
   double resultrelolo = 0.0;
   double resultimhihi = 0.0;
   double resultimlohi = 0.0;
   double resultimhilo = 0.0;
   double resultimlolo = 0.0;
   double a_hihi,a_lohi,a_hilo,a_lolo;
   double b_hihi,b_lohi,b_hilo,b_lolo;
   double acchihi,acclohi,acchilo,acclolo;
   int Widx,Yidx;

   for(int k=0; k<szt; k++)
   {
      Widx = k*nrows + row;             // if(nrows == szt) then row = bdx
      Yidx = k*nrows + col;            // if(nrows == szt) then col = tdx
      // result = result + a*b; with Hermitian transpose of Y
      // resultre = resultre + a_re*b_re + a_im*b_im;
      __syncthreads();
      a_hihi = Wrehihi[Widx];
      a_lohi = Wrelohi[Widx];
      a_hilo = Wrehilo[Widx];
      a_lolo = Wrelolo[Widx];
      b_hihi = Yrehihi[Yidx];       
      b_lohi = Yrelohi[Yidx];
      b_hilo = Yrehilo[Yidx];
      b_lolo = Yrelolo[Yidx];
      qdg_mul(  a_hihi,  a_lohi,  a_hilo,  a_lolo,
                b_hihi,  b_lohi,  b_hilo,  b_lolo,
              &acchihi,&acclohi,&acchilo,&acclolo);
      qdg_inc(&resultrehihi,&resultrelohi,&resultrehilo,&resultrelolo,
                    acchihi,      acclohi,      acchilo,      acclolo);

      a_hihi = Wimhihi[Widx]; 
      a_lohi = Wimlohi[Widx]; 
      a_hilo = Wimhilo[Widx]; 
      a_lolo = Wimlolo[Widx]; 
      b_hihi = Yimhihi[Yidx];
      b_lohi = Yimlohi[Yidx];
      b_hilo = Yimhilo[Yidx];
      b_lolo = Yimlolo[Yidx];
      qdg_mul(  a_hihi,  a_lohi,  a_hilo,  a_lolo,
                b_hihi,  b_lohi,  b_hilo,  b_lolo,
              &acchihi,&acclohi,&acchilo,&acclolo);
      qdg_inc(&resultrehihi,&resultrelohi,&resultrehilo,&resultrelolo,
                    acchihi,      acclohi,      acchilo,      acclolo);
      // resultim = resultim + a_im*b_re - a_re*b_im;
      // a_im is already available
      b_hihi = Yrehihi[Yidx];
      b_lohi = Yrelohi[Yidx];
      b_hilo = Yrehilo[Yidx];
      b_lolo = Yrelolo[Yidx];
      qdg_mul(  a_hihi,  a_lohi,  a_hilo,  a_lolo,
                b_hihi,  b_lohi,  b_hilo,  b_lolo,
              &acchihi,&acclohi,&acchilo,&acclolo);
      qdg_inc(&resultimhihi,&resultimlohi,&resultimhilo,&resultimlolo,
                    acchihi,      acclohi,      acchilo,      acclolo);

      a_hihi = Wrehihi[Widx];
      a_lohi = Wrelohi[Widx];
      a_hilo = Wrehilo[Widx];
      a_lolo = Wrelolo[Widx];
      b_hihi = Yimhihi[Yidx];
      b_lohi = Yimlohi[Yidx];
      b_hilo = Yimhilo[Yidx];
      b_lolo = Yimlolo[Yidx];
      qdg_mul(  a_hihi,  a_lohi,  a_hilo,  a_lolo,
                b_hihi,  b_lohi,  b_hilo,  b_lolo,
              &acchihi,&acclohi,&acchilo,&acclolo);
      qdg_dec(&resultimhihi,&resultimlohi,&resultimhilo,&resultimlolo,
                    acchihi,      acclohi,      acchilo,      acclolo);
   }
   __syncthreads();
   WYTrehihi[offset] = resultrehihi;
   WYTrelohi[offset] = resultrelohi;
   WYTrehilo[offset] = resultrehilo;
   WYTrelolo[offset] = resultrelolo;
   WYTimhihi[offset] = resultimhihi;
   WYTimlohi[offset] = resultimlohi;
   WYTimhilo[offset] = resultimhilo;
   WYTimlolo[offset] = resultimlolo;
}

__global__ void dbl4_small_QWYT
 ( int dim, int rowdim, int szt, int coloff,
   double *Qhihi, double *Qlohi, double *Qhilo, double *Qlolo,
   double *WYThihi, double *WYTlohi, double *WYThilo, double *WYTlolo,
   double *QWYThihi, double *QWYTlohi, double *QWYThilo, double *QWYTlolo )
{
   const int bdx = blockIdx.x;         // index of block
   const int tdx = threadIdx.x;        // index of thread in block
   const int offset = bdx*szt + tdx;   // offset in result
   const int row = offset / rowdim;
   const int col = offset % rowdim;    // thread 0 computes QWYT[row][col]

   double resulthihi = 0.0;
   double resultlohi = 0.0;
   double resulthilo = 0.0;
   double resultlolo = 0.0;
   double ahihi,alohi,ahilo,alolo;
   double bhihi,blohi,bhilo,blolo;
   double chihi,clohi,chilo,clolo;
   int idx;

   for(int k=0; k<rowdim; k++)       // run over rowdim, not just szt
   {                                 // coloff shifts by col*row elements
      idx = row*dim + coloff + k;
      __syncthreads();
      ahihi = Qhihi[idx];            // row = bdx,
      alohi = Qlohi[idx];
      ahilo = Qhilo[idx];            // if dim == szt, coloff == 0
      alolo = Qlolo[idx];
      idx = k*rowdim + col;
      __syncthreads();
      bhihi = WYThihi[idx];          // if(dim == szt) then col = tdx
      blohi = WYTlohi[idx];
      bhilo = WYThilo[idx];
      blolo = WYTlolo[idx]; 
      // result = result + a*b;
      __syncthreads();
      qdg_mul( ahihi, alohi, ahilo, alolo,
               bhihi, blohi, bhilo, blolo,
              &chihi,&clohi,&chilo,&clolo);
      qdg_inc(&resulthihi,&resultlohi,&resulthilo,&resultlolo,
                    chihi,      clohi,      chilo,      clolo);
   }
   __syncthreads();
   QWYThihi[offset] = resulthihi;    // no column offset in saving QWYT
   QWYTlohi[offset] = resultlohi;
   QWYThilo[offset] = resulthilo;
   QWYTlolo[offset] = resultlolo;
}

__global__ void cmplx4_small_QWYH
 ( int dim, int rowdim, int szt, int coloff,
   double *Qrehihi, double *Qrelohi, double *Qrehilo, double *Qrelolo,
   double *Qimhihi, double *Qimlohi, double *Qimhilo, double *Qimlolo,
   double *WYTrehihi, double *WYTrelohi, double *WYTrehilo, double *WYTrelolo,
   double *WYTimhihi, double *WYTimlohi, double *WYTimhilo, double *WYTimlolo,
   double *QWYTrehihi, double *QWYTrelohi,
   double *QWYTrehilo, double *QWYTrelolo,
   double *QWYTimhihi, double *QWYTimlohi,
   double *QWYTimhilo, double *QWYTimlolo )
{
   const int bdx = blockIdx.x;         // index of block
   const int tdx = threadIdx.x;        // index of thread in block
   const int offset = bdx*szt + tdx;   // offset in result
   const int row = offset / rowdim;
   const int col = offset % rowdim;    // thread 0 computes QWYT[row][col]

   double resultrehihi = 0.0;
   double resultrelohi = 0.0;
   double resultrehilo = 0.0;
   double resultrelolo = 0.0;
   double resultimhihi = 0.0;
   double resultimlohi = 0.0;
   double resultimhilo = 0.0;
   double resultimlolo = 0.0;
   double a_rehihi,a_relohi,a_rehilo,a_relolo;
   double a_imhihi,a_imlohi,a_imhilo,a_imlolo;
   double b_rehihi,b_relohi,b_rehilo,b_relolo;
   double b_imhihi,b_imlohi,b_imhilo,b_imlolo;
   double acchihi,acclohi,acchilo,acclolo;
   int Qidx,WYTidx;

   for(int k=0; k<rowdim; k++)          // run over rowdim, not just szt
   {                                    // coloff shifts by col*row elements
      Qidx = row*dim + coloff + k;
      __syncthreads();
      a_rehihi = Qrehihi[Qidx];         // row = bdx,
      a_relohi = Qrelohi[Qidx];
      a_rehilo = Qrehilo[Qidx];
      a_relolo = Qrelolo[Qidx];
      a_imhihi = Qimhihi[Qidx];         // if dim == szt, coloff == 0
      a_imlohi = Qimlohi[Qidx];
      a_imhilo = Qimhilo[Qidx];
      a_imlolo = Qimlolo[Qidx];
      WYTidx = k*rowdim + col;
      __syncthreads();
      b_rehihi = WYTrehihi[WYTidx];     // if(dim == szt) then col = tdx
      b_relohi = WYTrelohi[WYTidx];
      b_rehilo = WYTrehilo[WYTidx];
      b_relolo = WYTrelolo[WYTidx];
      b_imhihi = WYTimhihi[WYTidx];
      b_imlohi = WYTimlohi[WYTidx];
      b_imhilo = WYTimhilo[WYTidx];
      b_imlolo = WYTimlolo[WYTidx];
      // result = result + a*b;
      // resultre = resultre + a_re*b_re - a_im*b_im;
      __syncthreads();
      qdg_mul(a_rehihi,a_relohi,a_rehilo,a_relolo,
              b_rehihi,b_relohi,b_rehilo,b_relolo,
              &acchihi,&acclohi,&acchilo,&acclolo);
      qdg_inc(&resultrehihi,&resultrelohi,&resultrehilo,&resultrelolo,
                    acchihi,      acclohi,      acchilo,      acclolo);
      qdg_mul(a_imhihi,a_imlohi,a_imhilo,a_imlolo,
              b_imhihi,b_imlohi,b_imhilo,b_imlolo,
              &acchihi,&acclohi,&acchilo,&acclolo);
      qdg_dec(&resultrehihi,&resultrelohi,&resultrehilo,&resultrelolo,
                    acchihi,      acclohi,      acchilo,      acclolo);
      // resultim = resultim + a_im*b_re + a_re*b_im;
      qdg_mul(a_imhihi,a_imlohi,a_imhilo,a_imlolo,
              b_rehihi,b_relohi,b_rehilo,b_relolo,
              &acchihi,&acclohi,&acchilo,&acclolo);
      qdg_inc(&resultimhihi,&resultimlohi,&resultimhilo,&resultimlolo,
                    acchihi,      acclohi,      acchilo,      acclolo);
      qdg_mul(a_rehihi,a_relohi,a_rehilo,a_relolo,
              b_imhihi,b_imlohi,b_imhilo,b_imlolo,
              &acchihi,&acclohi,&acchilo,&acclolo);
      qdg_inc(&resultimhihi,&resultimlohi,&resultimhilo,&resultimlolo,
                    acchihi,      acclohi,      acchilo,      acclolo);
   }
   __syncthreads();
   QWYTrehihi[offset] = resultrehihi;    // no column offset in saving QWYT
   QWYTrelohi[offset] = resultrelohi;
   QWYTrehilo[offset] = resultrehilo;
   QWYTrelolo[offset] = resultrelolo;
   QWYTimhihi[offset] = resultimhihi;
   QWYTimlohi[offset] = resultimlohi;
   QWYTimhilo[offset] = resultimhilo;
   QWYTimlolo[offset] = resultimlolo;
}

__global__ void dbl4_small_YWTC
 ( int nrows, int ncols, int rowdim, int coldim, int szt,
   int rowoff, int coloff,
   double *YWThihi, double *YWTlohi, double *YWThilo, double *YWTlolo,
   double *Chihi, double *Clohi, double *Chilo, double *Clolo,
   double *YWTChihi, double *YWTClohi, double *YWTChilo, double *YWTClolo )
{
   const int bdx = blockIdx.x;         // bdx*szt done by previous blocks
   const int tdx = threadIdx.x;        // index of thread in block
   const int offset = bdx*szt + tdx;   // offset in result
   const int row = offset / coldim;    // 1st thread does YWTC[row][col]
   const int col = offset % coldim;
   const int colCoff0 = (coloff+col)*nrows + rowoff; // 1st element in C

   double resulthihi = 0.0;
   double resultlohi = 0.0;
   double resulthilo = 0.0;
   double resultlolo = 0.0;
   double ahihi,alohi,ahilo,alolo;
   double bhihi,blohi,bhilo,blolo;
   double chihi,clohi,chilo,clolo;
   int idx;

   for(int k=0; k<rowdim; k++)         // innermost loop runs over rowdim
   {
      idx = row*rowdim + k;
      __syncthreads();
      ahihi = YWThihi[idx];            // YWT is stored row by row
      alohi = YWTlohi[idx];
      ahilo = YWThilo[idx];
      alolo = YWTlolo[idx];
      idx = colCoff0 + k;
      __syncthreads();
      bhihi = Chihi[idx];              // but C is stored column by column
      blohi = Clohi[idx];
      bhilo = Chilo[idx];
      blolo = Clolo[idx];
      // result = result + a*b;
      __syncthreads();
      qdg_mul( ahihi, alohi, ahilo, alolo,
               bhihi, blohi, bhilo, blolo,
              &chihi,&clohi,&chilo,&clolo);
      qdg_inc(&resulthihi,&resultlohi,&resulthilo,&resultlolo,
                    chihi,      clohi,      chilo,      clolo);
   }
   idx = (coloff + col)*nrows + (rowoff + row);
   __syncthreads();
   YWTChihi[idx] = resulthihi;
   YWTClohi[idx] = resultlohi;
   YWTChilo[idx] = resulthilo;
   YWTClolo[idx] = resultlolo;
}

__global__ void cmplx4_small_YWHC
 ( int nrows, int ncols, int rowdim, int coldim, int szt,
   int rowoff, int coloff,
   double *YWTrehihi, double *YWTrelohi, double *YWTrehilo, double *YWTrelolo,
   double *YWTimhihi, double *YWTimlohi, double *YWTimhilo, double *YWTimlolo,
   double *Crehihi, double *Crelohi, double *Crehilo, double *Crelolo,
   double *Cimhihi, double *Cimlohi, double *Cimhilo, double *Cimlolo,
   double *YWTCrehihi, double *YWTCrelohi,
   double *YWTCrehilo, double *YWTCrelolo,
   double *YWTCimhihi, double *YWTCimlohi,
   double *YWTCimhilo, double *YWTCimlolo )
{
   const int bdx = blockIdx.x;         // bdx*szt done by previous blocks
   const int tdx = threadIdx.x;        // index of thread in block
   const int offset = bdx*szt + tdx;   // offset in result
   const int row = offset / coldim;    // 1st thread does YWTC[row][col]
   const int col = offset % coldim;
   const int colCoff0 = (coloff+col)*nrows + rowoff; // 1st element in C
   int idx;

   double resultrehihi = 0.0;
   double resultrelohi = 0.0;
   double resultrehilo = 0.0;
   double resultrelolo = 0.0;
   double resultimhihi = 0.0;
   double resultimlohi = 0.0;
   double resultimhilo = 0.0;
   double resultimlolo = 0.0;
   double a_hihi,a_lohi,a_hilo,a_lolo;
   double b_hihi,b_lohi,b_hilo,b_lolo;
   double acchihi,acclohi,acchilo,acclolo;
   int YWTidx,Cidx;

   for(int k=0; k<rowdim; k++)         // innermost loop runs over rowdim
   {
      YWTidx = row*rowdim + k;         // YWT is stored row by row
      Cidx = colCoff0 + k;             // but C is stored column by column
      // result = result + a*b;
      // resultre = resultre + a_re*b_re - a_im*b_im;
      __syncthreads();
      a_hihi = YWTrehihi[YWTidx];    // YWT is stored row by row
      a_lohi = YWTrelohi[YWTidx];
      a_hilo = YWTrehilo[YWTidx];
      a_lolo = YWTrelolo[YWTidx];
      b_hihi = Crehihi[Cidx];        // but C is stored column by column
      b_lohi = Crelohi[Cidx];
      b_hilo = Crehilo[Cidx];
      b_lolo = Crelolo[Cidx];
      qdg_mul(  a_hihi,  a_lohi,  a_hilo,  a_lolo,
                b_hihi,  b_lohi,  b_hilo,  b_lolo,
              &acchihi,&acclohi,&acchilo,&acclolo);
      qdg_inc(&resultrehihi,&resultrelohi,&resultrehilo,&resultrelolo,
                    acchihi,      acclohi,      acchilo,      acclolo);

      a_hihi = YWTimhihi[YWTidx];
      a_lohi = YWTimlohi[YWTidx];
      a_hilo = YWTimhilo[YWTidx];
      a_lolo = YWTimlolo[YWTidx];
      b_hihi = Cimhihi[Cidx];
      b_lohi = Cimlohi[Cidx];
      b_hilo = Cimhilo[Cidx];
      b_lolo = Cimlolo[Cidx];
      qdg_mul(  a_hihi,  a_lohi,  a_hilo,  a_lolo,
                b_hihi,  b_lohi,  b_hilo,  b_lolo,
              &acchihi,&acclohi,&acchilo,&acclolo);
      qdg_dec(&resultrehihi,&resultrelohi,&resultrehilo,&resultrelolo,
                    acchihi,      acclohi,      acchilo,      acclolo);
      // resultim = resultim + a_im*b_re + a_re*b_im;
      // a_im is already available
      b_hihi = Crehihi[Cidx];
      b_lohi = Crelohi[Cidx];
      b_hilo = Crehilo[Cidx];
      b_lolo = Crelolo[Cidx];
      qdg_mul(  a_hihi,  a_lohi,  a_hilo,  a_lolo,
                b_hihi,  b_lohi,  b_hilo,  b_lolo,
              &acchihi,&acclohi,&acchilo,&acclolo);
      qdg_inc(&resultimhihi,&resultimlohi,&resultimhilo,&resultimlolo,
                    acchihi,      acclohi,      acchilo,      acclolo);

      a_hihi = YWTrehihi[YWTidx];
      a_lohi = YWTrelohi[YWTidx];
      a_hilo = YWTrehilo[YWTidx];
      a_lolo = YWTrelolo[YWTidx];
      b_hihi = Cimhihi[Cidx];
      b_lohi = Cimlohi[Cidx];
      b_hilo = Cimhilo[Cidx];
      b_lolo = Cimlolo[Cidx];
      qdg_mul(  a_hihi,  a_lohi,  a_hilo,  a_lolo,
                b_hihi,  b_lohi,  b_hilo,  b_lolo,
              &acchihi,&acclohi,&acchilo,&acclolo);
      qdg_inc(&resultimhihi,&resultimlohi,&resultimhilo,&resultimlolo,
                    acchihi,      acclohi,      acchilo,      acclolo);
   }
   __syncthreads();
   idx = (coloff + col)*nrows + (rowoff + row);
   YWTCrehihi[idx] = resultrehihi;
   YWTCrelohi[idx] = resultrelohi;
   YWTCrehilo[idx] = resultrehilo;
   YWTCrelolo[idx] = resultrelolo;
   YWTCimhihi[idx] = resultimhihi;
   YWTCimlohi[idx] = resultimlohi;
   YWTCimhilo[idx] = resultimhilo;
   YWTCimlolo[idx] = resultimlolo;
}

__global__ void dbl4_small_Qupdate
 ( int dim, int rowdim, int szt, int coloff,
   double *Qhihi, double *Qlohi, double *Qhilo, double *Qlolo,
   double *QWYThihi, double *QWYTlohi, double *QWYThilo, double *QWYTlolo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int offset = bdx*szt + tdx;   // offset in result
   const int row = offset / rowdim;
   const int col = offset % rowdim;
   const int idx1 = row*dim + coloff + col;

   double ahihi,alohi,ahilo,alolo;
   double bhihi,blohi,bhilo,blolo;

   ahihi = Qhihi[idx1];       // row = bdx, if dim == szt, coloff == 0
   alohi = Qlohi[idx1];
   ahilo = Qhilo[idx1];
   alolo = Qlolo[idx1];
   __syncthreads();
   bhihi = QWYThihi[offset];  // if(dim == szt) then col = tdx
   blohi = QWYTlohi[offset];
   bhilo = QWYThilo[offset];
   blolo = QWYTlolo[offset];
   // a = a + b;
   __syncthreads();
   qdg_inc(&ahihi,&alohi,&ahilo,&alolo,
            bhihi, blohi, bhilo, blolo);
   __syncthreads();
   Qhihi[idx1] = ahihi;
   Qlohi[idx1] = alohi;
   Qhilo[idx1] = ahilo;
   Qlolo[idx1] = alolo;
}

__global__ void cmplx4_small_Qupdate
 ( int dim, int rowdim, int szt, int coloff,
   double *Qrehihi, double *Qrelohi, double *Qrehilo, double *Qrelolo,
   double *Qimhihi, double *Qimlohi, double *Qimhilo, double *Qimlolo,
   double *QWYTrehihi, double *QWYTrelohi,
   double *QWYTrehilo, double *QWYTrelolo,
   double *QWYTimhihi, double *QWYTimlohi,
   double *QWYTimhilo, double *QWYTimlolo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int offset = bdx*szt + tdx;   // offset in result
   const int row = offset / rowdim;
   const int col = offset % rowdim;
   const int idx1 = row*dim + coloff + col;

   double a_rehihi,a_relohi,a_rehilo,a_relolo;
   double a_imhihi,a_imlohi,a_imhilo,a_imlolo;
   double b_rehihi,b_relohi,b_rehilo,b_relolo;
   double b_imhihi,b_imlohi,b_imhilo,b_imlolo;

   a_rehihi = Qrehihi[idx1];       // row = bdx, if dim == szt, coloff == 0
   a_relohi = Qrelohi[idx1];
   a_rehilo = Qrehilo[idx1];
   a_relolo = Qrelolo[idx1];
   a_imhihi = Qimhihi[idx1];
   a_imlohi = Qimlohi[idx1];
   a_imhilo = Qimhilo[idx1];
   a_imlolo = Qimlolo[idx1];
   __syncthreads();
   b_rehihi = QWYTrehihi[offset];  // if(dim == szt) then col = tdx
   b_relohi = QWYTrelohi[offset];
   b_rehilo = QWYTrehilo[offset];
   b_relolo = QWYTrelolo[offset];
   b_imhihi = QWYTimhihi[offset];
   b_imlohi = QWYTimlohi[offset];
   b_imhilo = QWYTimhilo[offset];
   b_imlolo = QWYTimlolo[offset];
   // a_re = a_re + b_re;
   __syncthreads();
   qdg_inc(&a_rehihi,&a_relohi,&a_rehilo,&a_relolo,
            b_rehihi, b_relohi, b_rehilo, b_relolo);
   // a_im = a_im + b_im;
   qdg_inc(&a_imhihi,&a_imlohi,&a_imhilo,&a_imlolo,
            b_imhihi, b_imlohi, b_imhilo, b_imlolo);

   __syncthreads();
   Qrehihi[idx1] = a_rehihi;
   Qrelohi[idx1] = a_relohi;
   Qrehilo[idx1] = a_rehilo;
   Qrelolo[idx1] = a_relolo;
   Qimhihi[idx1] = a_imhihi;
   Qimlohi[idx1] = a_imlohi;
   Qimhilo[idx1] = a_imhilo;
   Qimlolo[idx1] = a_imlolo;
}

__global__ void dbl4_small_R_add_YWTC
 ( int nrows, int coldim, int szt, int rowoff, int coloff,
   double *Rhihi, double *Rlohi, double *Rhilo, double *Rlolo,
   double *YWTChihi, double *YWTClohi, double *YWTChilo, double *YWTClolo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int offset = bdx*szt + tdx;   // offset in result
   const int row = offset / coldim;    // thread updates R[row][col]
   const int col = offset % coldim;
   const int idx = (coloff + col)*nrows + (rowoff + row);
 
   double ahihi,alohi,ahilo,alolo;
   double bhihi,blohi,bhilo,blolo;
   
   ahihi = Rhihi[idx];
   alohi = Rlohi[idx];
   ahilo = Rhilo[idx];
   alolo = Rlolo[idx];
   __syncthreads();
   bhihi = YWTChihi[idx];
   blohi = YWTClohi[idx];
   bhilo = YWTChilo[idx];
   blolo = YWTClolo[idx];
   // a = a + b;
   __syncthreads();
   qdg_inc(&ahihi,&alohi,&ahilo,&alolo,
            bhihi, blohi, bhilo, blolo);
  
   __syncthreads();
   Rhihi[idx] = ahihi;
   Rlohi[idx] = alohi;
   Rhilo[idx] = ahilo;
   Rlolo[idx] = alolo;
}

__global__ void cmplx4_small_R_add_YWHC
 ( int nrows, int coldim, int szt, int rowoff, int coloff,
   double *Rrehihi, double *Rrelohi, double *Rrehilo, double *Rrelolo,
   double *Rimhihi, double *Rimlohi, double *Rimhilo, double *Rimlolo,
   double *YWTCrehihi, double *YWTCrelohi,
   double *YWTCrehilo, double *YWTCrelolo,
   double *YWTCimhihi, double *YWTCimlohi,
   double *YWTCimhilo, double *YWTCimlolo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int offset = bdx*szt + tdx;   // offset in result
   const int row = offset / coldim;    // thread updates R[row][col]
   const int col = offset % coldim;
   const int idx = (coloff + col)*nrows + (rowoff + row);
 
   double a_hihi,a_lohi,a_hilo,a_lolo;
   double b_hihi,b_lohi,b_hilo,b_lolo;
   
   a_hihi = Rrehihi[idx];
   a_lohi = Rrelohi[idx];
   a_hilo = Rrehilo[idx];
   a_lolo = Rrelolo[idx];
   b_hihi = YWTCrehihi[idx];
   b_lohi = YWTCrelohi[idx];
   b_hilo = YWTCrehilo[idx];
   b_lolo = YWTCrelolo[idx];
   // a_re = a_re + b_re;
   __syncthreads();
   qdg_inc(&a_hihi,&a_lohi,&a_hilo,&a_lolo,
            b_hihi, b_lohi, b_hilo, b_lolo);
   __syncthreads();
   Rrehihi[idx] = a_hihi;
   Rrelohi[idx] = a_lohi;
   Rrehilo[idx] = a_hilo;
   Rrelolo[idx] = a_lolo;

   // a_im = a_im + b_im;
   a_hihi = Rimhihi[idx];
   a_lohi = Rimlohi[idx];
   a_hilo = Rimhilo[idx];
   a_lolo = Rimlolo[idx];
   b_hihi = YWTCimhihi[idx];
   b_lohi = YWTCimlohi[idx];
   b_hilo = YWTCimhilo[idx];
   b_lolo = YWTCimlolo[idx];
   __syncthreads();
   qdg_inc(&a_hihi,&a_lohi,&a_hilo,&a_lolo,
            b_hihi, b_lohi, b_hilo, b_lolo);
   __syncthreads();
   Rimhihi[idx] = a_hihi;
   Rimlohi[idx] = a_lohi;
   Rimhilo[idx] = a_hilo;
   Rimlolo[idx] = a_lolo;
}

void GPU_dbl4_small_house
 ( int nrows, int ncols, int szt, int nbt,
   int colidx, int nrows1, int k, int L,
   double *Ahihi_h, double *Alohi_h, double *Ahilo_h, double *Alolo_h,
   double *Ahihi_d, double *Alohi_d, double *Ahilo_d, double *Alolo_d,
   double *vhihi_h, double *vlohi_h, double *vhilo_h, double *vlolo_h,
   double *Vhihi_d, double *Vlohi_d, double *Vhilo_d, double *Vlolo_d,
   double *betahihi_h, double *betalohi_h,
   double *betahilo_h, double *betalolo_h,
   double *betahihi_d, double *betalohi_d,
   double *betahilo_d, double *betalolo_d,
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
         vhihi_h[i] = 0.0;
         vlohi_h[i] = 0.0;
         vhilo_h[i] = 0.0;
         vlolo_h[i] = 0.0;
      }
      cudaMemcpy(&Vhihi_d[L*nVrows],vhihi_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vlohi_d[L*nVrows],vlohi_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vhilo_d[L*nVrows],vhilo_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vlolo_d[L*nVrows],vlolo_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
   }
   if(nrows1 == 0)
   {
      betahihi_h[L] = 0.0;
      betalohi_h[L] = 0.0;
      betahilo_h[L] = 0.0;
      betalolo_h[L] = 0.0;
      vhihi_h[0] = 1.0;
      vlohi_h[0] = 0.0;
      vhilo_h[0] = 0.0;
      vlolo_h[0] = 0.0;
      cudaMemcpy(&betahihi_d[L],&betahihi_h[L],sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&betalohi_d[L],&betalohi_h[L],sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&betahilo_d[L],&betahilo_h[L],sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&betalolo_d[L],&betalolo_h[L],sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vhihi_d[L*nVrows+L],vhihi_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vlohi_d[L*nVrows+L],vlohi_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vhilo_d[L*nVrows+L],vhilo_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vlolo_d[L*nVrows+L],vlolo_h,sizeof(double),
                 cudaMemcpyHostToDevice);
   }
   else
   {
      cudaEvent_t start,stop;           // to measure time spent by kernels 
      cudaEventCreate(&start);
      cudaEventCreate(&stop);
      float milliseconds;

      cudaEventRecord(start);
      dbl4_small_house<<<1,nrows1>>>
         (&Ahihi_d[rowidx],&Alohi_d[rowidx],&Ahilo_d[rowidx],&Alolo_d[rowidx],
          &Ahihi_d[rowidx+1],&Alohi_d[rowidx+1],
          &Ahilo_d[rowidx+1],&Alolo_d[rowidx+1],
          nrows1,nrLog2,&Vhihi_d[L*nVrows+L],&Vlohi_d[L*nVrows+L],
                        &Vhilo_d[L*nVrows+L],&Vlolo_d[L*nVrows+L],
          &betahihi_d[L],&betalohi_d[L],&betahilo_d[L],&betalolo_d[L]);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds,start,stop);
      *lapms += milliseconds;
      flopcount_dbl_small_house(nrows1,nrLog2,add,mul,div,sqrtfun);
   }
   cudaMemcpy(&betahihi_h[L],&betahihi_d[L],sizeof(double),
              cudaMemcpyDeviceToHost);
   cudaMemcpy(&betalohi_h[L],&betalohi_d[L],sizeof(double),
              cudaMemcpyDeviceToHost);
   cudaMemcpy(&betahilo_h[L],&betahilo_d[L],sizeof(double),
              cudaMemcpyDeviceToHost);
   cudaMemcpy(&betalolo_h[L],&betalolo_d[L],sizeof(double),
              cudaMemcpyDeviceToHost);
   if(verbose)
   {
      const size_t szhouse = nVrows*sizeof(double);

      cudaMemcpy(vhihi_h,&Vhihi_d[L*nVrows],szhouse,cudaMemcpyDeviceToHost);
      cudaMemcpy(vlohi_h,&Vlohi_d[L*nVrows],szhouse,cudaMemcpyDeviceToHost);
      cudaMemcpy(vhilo_h,&Vhilo_d[L*nVrows],szhouse,cudaMemcpyDeviceToHost);
      cudaMemcpy(vlolo_h,&Vlolo_d[L*nVrows],szhouse,cudaMemcpyDeviceToHost);
      cout << scientific << setprecision(16)
           << "beta[" << colidx << "] : "
           << betahihi_h[L] << "  " << betalohi_h[L] << endl
           << "          "
           << betahilo_h[L] << "  " << betalolo_h[L] << endl;
      for(int i=0; i<nVrows; i++)
         cout << "v[" << i << "] : "
              << vhihi_h[i] << "  " << vlohi_h[i] << endl
              << "       "
              << vhilo_h[i] << "  " << vlolo_h[i] << endl;
   }
}

void GPU_cmplx4_small_house
 ( int nrows, int ncols, int szt, int nbt,
   int colidx, int nrows1, int k, int L,
   double *Arehihi_h, double *Arelohi_h, double *Arehilo_h, double *Arelolo_h,
   double *Aimhihi_h, double *Aimlohi_h, double *Aimhilo_h, double *Aimlolo_h,
   double *Arehihi_d, double *Arelohi_d, double *Arehilo_d, double *Arelolo_d,
   double *Aimhihi_d, double *Aimlohi_d, double *Aimhilo_d, double *Aimlolo_d,
   double *vrehihi_h, double *vrelohi_h, double *vrehilo_h, double *vrelolo_h,
   double *vimhihi_h, double *vimlohi_h, double *vimhilo_h, double *vimlolo_h,
   double *Vrehihi_d, double *Vrelohi_d, double *Vrehilo_d, double *Vrelolo_d,
   double *Vimhihi_d, double *Vimlohi_d, double *Vimhilo_d, double *Vimlolo_d,
   double *betahihi_h, double *betalohi_h,
   double *betahilo_h, double *betalolo_h,
   double *betahihi_d, double *betalohi_d,
   double *betahilo_d, double *betalolo_d,
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
         vrehihi_h[i] = 0.0; vrelohi_h[i] = 0.0;
         vrehilo_h[i] = 0.0; vrelolo_h[i] = 0.0;
         vimhihi_h[i] = 0.0; vimlohi_h[i] = 0.0;
         vimhilo_h[i] = 0.0; vimlolo_h[i] = 0.0;
      }
      cudaMemcpy(&Vrehihi_d[L*nVrows],vrehihi_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vrelohi_d[L*nVrows],vrelohi_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vrehilo_d[L*nVrows],vrehilo_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vrelolo_d[L*nVrows],vrelolo_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vimhihi_d[L*nVrows],vimhihi_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vimlohi_d[L*nVrows],vimlohi_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vimhilo_d[L*nVrows],vimhilo_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vimlolo_d[L*nVrows],vimlolo_h,L*sizeof(double),
                 cudaMemcpyHostToDevice);
   }
   if(nrows1 == 0)
   {
      betahihi_h[L] = 0.0; betalohi_h[L] = 0.0;
      betahilo_h[L] = 0.0; betalolo_h[L] = 0.0;
      vrehihi_h[0] = 1.0; vrelohi_h[0] = 0.0;
      vrehilo_h[0] = 0.0; vrelolo_h[0] = 0.0;
      vimhihi_h[0] = 0.0; vimlohi_h[0] = 0.0;
      vimhilo_h[0] = 0.0; vimlolo_h[0] = 0.0;
      cudaMemcpy(&betahihi_d[L],&betahihi_h[L],sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&betalohi_d[L],&betalohi_h[L],sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&betahilo_d[L],&betahilo_h[L],sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&betalolo_d[L],&betalolo_h[L],sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vrehihi_d[L*nVrows+L],vrehihi_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vrelohi_d[L*nVrows+L],vrelohi_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vrehilo_d[L*nVrows+L],vrehilo_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vrelolo_d[L*nVrows+L],vrelolo_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vimhihi_d[L*nVrows+L],vimhihi_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vimlohi_d[L*nVrows+L],vimlohi_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vimhilo_d[L*nVrows+L],vimhilo_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&Vimlolo_d[L*nVrows+L],vimlolo_h,sizeof(double),
                 cudaMemcpyHostToDevice);
   }
   else
   {
      if(verbose)  // to verify the input column is correct ...
      {
         cout << "The column of A : " << endl;
         cudaMemcpy(&Arehihi_h[rowidx],&Arehihi_d[rowidx],
                   (nrows1+1)*sizeof(double),cudaMemcpyDeviceToHost);
         cudaMemcpy(&Arelohi_h[rowidx],&Arelohi_d[rowidx],
                   (nrows1+1)*sizeof(double),cudaMemcpyDeviceToHost);
         cudaMemcpy(&Arehilo_h[rowidx],&Arehilo_d[rowidx],
                   (nrows1+1)*sizeof(double),cudaMemcpyDeviceToHost);
         cudaMemcpy(&Arelolo_h[rowidx],&Arelolo_d[rowidx],
                   (nrows1+1)*sizeof(double),cudaMemcpyDeviceToHost);
         cudaMemcpy(&Aimhihi_h[rowidx],&Aimhihi_d[rowidx],
                   (nrows1+1)*sizeof(double),cudaMemcpyDeviceToHost);
         cudaMemcpy(&Aimlohi_h[rowidx],&Aimlohi_d[rowidx],
                   (nrows1+1)*sizeof(double),cudaMemcpyDeviceToHost);
         cudaMemcpy(&Aimhilo_h[rowidx],&Aimhilo_d[rowidx],
                   (nrows1+1)*sizeof(double),cudaMemcpyDeviceToHost);
         cudaMemcpy(&Aimlolo_h[rowidx],&Aimlolo_d[rowidx],
                   (nrows1+1)*sizeof(double),cudaMemcpyDeviceToHost);
         for(int i=0; i<=nrows1; i++)
         {
            cout << "A[" << i << "]re : " 
                 << Arehihi_h[i] << "  " << Arelohi_h[i] << endl
                 << "         "
                 << Arehilo_h[i] << "  " << Arelolo_h[i] << endl;
            cout << "A[" << i << "]im : " 
                 << Aimhihi_h[i] << "  " << Aimlohi_h[i] << endl
                 << "         "
                 << Aimhilo_h[i] << "  " << Aimlolo_h[i] << endl;
         }
      }
      cudaEvent_t start,stop;           // to measure time spent by kernels 
      cudaEventCreate(&start);
      cudaEventCreate(&stop);
      float milliseconds;

      cudaEventRecord(start);
      cmplx4_small_house<<<1,nrows1>>>
         (&Arehihi_d[rowidx],&Arelohi_d[rowidx],
          &Arehilo_d[rowidx],&Arelolo_d[rowidx],
          &Aimhihi_d[rowidx],&Aimlohi_d[rowidx],
          &Aimhilo_d[rowidx],&Aimlolo_d[rowidx],
          &Arehihi_d[rowidx+1],&Arelohi_d[rowidx+1],
          &Arehilo_d[rowidx+1],&Arelolo_d[rowidx+1],
          &Aimhihi_d[rowidx+1],&Aimlohi_d[rowidx+1],
          &Aimhilo_d[rowidx+1],&Aimlolo_d[rowidx+1],nrows1,nrLog2,
          &Vrehihi_d[L*nVrows+L],&Vrelohi_d[L*nVrows+L],
          &Vrehilo_d[L*nVrows+L],&Vrelolo_d[L*nVrows+L],
          &Vimhihi_d[L*nVrows+L],&Vimlohi_d[L*nVrows+L],
          &Vimhilo_d[L*nVrows+L],&Vimlolo_d[L*nVrows+L],
          &betahihi_d[L],&betalohi_d[L],&betahilo_d[L],&betalolo_d[L]);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds,start,stop);
      *lapms += milliseconds;
      flopcount_cmplx_small_house(nrows1,nrLog2,add,mul,div,sqrtfun);
   }
   cudaMemcpy(&betahihi_h[L],&betahihi_d[L],sizeof(double),
              cudaMemcpyDeviceToHost);
   cudaMemcpy(&betalohi_h[L],&betalohi_d[L],sizeof(double),
              cudaMemcpyDeviceToHost);
   cudaMemcpy(&betahilo_h[L],&betahilo_d[L],sizeof(double),
              cudaMemcpyDeviceToHost);
   cudaMemcpy(&betalolo_h[L],&betalolo_d[L],sizeof(double),
              cudaMemcpyDeviceToHost);
   if(verbose)
   {
      const size_t szhouse = nVrows*sizeof(double);

      cudaMemcpy(vrehihi_h,&Vrehihi_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vrelohi_h,&Vrelohi_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vrehilo_h,&Vrehilo_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vrelolo_h,&Vrelolo_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vimhihi_h,&Vimhihi_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vimlohi_h,&Vimlohi_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vimhilo_h,&Vimhilo_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vimlolo_h,&Vimlolo_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cout << scientific << setprecision(16)
           << "beta[" << colidx << "] : "
           << betahihi_h[L] << "  " << betalohi_h[L] << endl
           << "          "
           << betahilo_h[L] << "  " << betalolo_h[L] << endl;
      for(int i=0; i<nVrows; i++)
      {
         cout << "v[" << i << "]re : "
              << vrehihi_h[i] << "  " << vrelohi_h[i] << endl
           << "          "
              << vrehilo_h[i] << "  " << vrelolo_h[i] << endl;
         cout << "v[" << i << "]im : "
              << vimhihi_h[i] << "  " << vimlohi_h[i] << endl
           << "          "
              << vimhilo_h[i] << "  " << vimlolo_h[i] << endl;
      }
   }
}

void GPU_dbl4_large_house
 ( int nrows, int ncols, int szt, int nbt,
   int colidx, int nrows1, int k, int L,
   double *Ahihi_h, double *Alohi_h, double *Ahilo_h, double *Alolo_h,
   double *Ahihi_d, double *Alohi_d, double *Ahilo_d, double *Alolo_d,
   double *vhihi_h, double *vlohi_h, double *vhilo_h, double *vlolo_h,
   double *Vhihi_d, double *Vlohi_d, double *Vhilo_d, double *Vlolo_d,
   double *betahihi_h, double *betalohi_h,
   double *betahilo_h, double *betalolo_h,
   double *betahihi_d, double *betalohi_d,
   double *betahilo_d, double *betalolo_d,
   double *sumshihi_h, double *sumslohi_h,
   double *sumshilo_h, double *sumslolo_h,
   double *sumshihi_d, double *sumslohi_d,
   double *sumshilo_d, double *sumslolo_d,
   double *sigmahihi_h, double *sigmalohi_h, 
   double *sigmahilo_h, double *sigmalolo_h, 
   double *sigmahihi_d, double *sigmalohi_d,
   double *sigmahilo_d, double *sigmalolo_d,
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
         vhihi_h[i] = 0.0;
         vlohi_h[i] = 0.0;
         vhilo_h[i] = 0.0;
         vlolo_h[i] = 0.0;
      }
   }
   vhihi_h[L] = 1.0;                    // set one on the diagonal
   vlohi_h[L] = 0.0;
   vhilo_h[L] = 0.0;
   vlolo_h[L] = 0.0;

   cudaMemcpy(&Vhihi_d[L*nVrows],vhihi_h,(L+1)*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&Vlohi_d[L*nVrows],vlohi_h,(L+1)*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&Vhilo_d[L*nVrows],vhilo_h,(L+1)*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&Vlolo_d[L*nVrows],vlolo_h,(L+1)*sizeof(double),
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
      sumshihi_h[i] = 0.0;
      sumslohi_h[i] = 0.0;
      sumshilo_h[i] = 0.0;
      sumslolo_h[i] = 0.0;
   }
   cudaMemcpy(sumshihi_d,sumshihi_h,nblocks*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(sumslohi_d,sumslohi_h,nblocks*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(sumshilo_d,sumshilo_h,nblocks*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(sumslolo_d,sumslolo_h,nblocks*sizeof(double),
              cudaMemcpyHostToDevice);

   cudaEventRecord(start);
   dbl4_large_sum_of_squares<<<nblocks,szt>>>
      (&Ahihi_d[rowidx+1],&Alohi_d[rowidx+1],
       &Ahilo_d[rowidx+1],&Alolo_d[rowidx+1],
       sumshihi_d,sumslohi_d,sumshilo_d,sumslolo_d,nrows1,szt,sztLog2);
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
   dbl4_sum_accumulator<<<1,nblocks>>>
      (sumshihi_d,sumslohi_d,sumshilo_d,sumslolo_d,
       nblocks,nblLog2,sigmahihi_d,sigmalohi_d,sigmahilo_d,sigmalolo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_dbl_sum_accumulator(nblocks,nblLog2,add);

   cudaMemcpy(sigmahihi_h,sigmahihi_d,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(sigmalohi_h,sigmalohi_d,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(sigmahilo_h,sigmahilo_d,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(sigmalolo_h,sigmalolo_d,sizeof(double),cudaMemcpyDeviceToHost);

   bool done = false;

   if((sigmahihi_h[0] == 0.0) && (sigmalohi_h[0] == 0.0) &&
      (sigmahilo_h[0] == 0.0) && (sigmalolo_h[0] == 0.0))
   {
      betahihi_h[L] = 0.0; betalohi_h[L] = 0.0;
      betahilo_h[L] = 0.0; betalolo_h[L] = 0.0; done = true;

      if(verbose)
         cout << "Zero sigma value encountered." << endl;
   }
   else // beta is computed on the host instead of by one GPU thread
   {
      // const double x0hi = Ahi_h[rowidx];
      // const double x0lo = Alo_h[rowidx];
      double acchihi,acclohi,acchilo,acclolo;
      double muhihi,mulohi,muhilo,mulolo;
      double v0hihi,v0lohi,v0hilo,v0lolo;
      double v0p2hihi,v0p2lohi,v0p2hilo,v0p2lolo;
      double x0hihi,x0lohi,x0hilo,x0lolo;

      cudaMemcpy(&x0hihi,&Ahihi_d[rowidx],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&x0lohi,&Alohi_d[rowidx],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&x0hilo,&Ahilo_d[rowidx],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&x0lolo,&Alolo_d[rowidx],sizeof(double),
                 cudaMemcpyDeviceToHost);

      // mu = sqrt((*x0)*(*x0) + sigma[0]);
      qdf_sqr(x0hihi,x0lohi,x0hilo,x0lolo,
              &acchihi,&acclohi,&acchilo,&acclolo);
      qdf_inc(&acchihi,&acclohi,&acchilo,&acclolo,
              sigmahihi_h[0],sigmalohi_h[0],sigmahilo_h[0],sigmalolo_h[0]);
      qdf_sqrt(acchihi,acclohi,acchilo,acclolo,
               &muhihi,&mulohi,&muhilo,&mulolo);
      if(x0hihi <= 0.0)
      {
         // v0 = *x0 - mu;
         qdf_sub( x0hihi, x0lohi, x0hilo, x0lolo,
                  muhihi, mulohi, muhilo, mulolo,
                 &v0hihi,&v0lohi,&v0hilo,&v0lolo);
      }
      else
      {
         // v0 = -sigma[0]/(*x0 + mu);
         qdf_add(  x0hihi,  x0lohi,  x0hilo,  x0lolo,
                   muhihi,  mulohi,  muhilo,  mulolo,
                 &acchihi,&acclohi,&acchilo,&acclolo);
         qdf_div(sigmahihi_h[0],sigmalohi_h[0],sigmahilo_h[0],sigmalolo_h[0],
                   acchihi,       acclohi,       acchilo,       acclolo,
                   &v0hihi,       &v0lohi,       &v0hilo,       &v0lolo);
         qdf_minus(&v0hihi,&v0lohi,&v0hilo,&v0lolo);
      }
      // v0p2 = v0*v0;
      qdf_sqr(   v0hihi,   v0lohi,   v0hilo,   v0lolo,
              &v0p2hihi,&v0p2lohi,&v0p2hilo,&v0p2lolo);
      // *beta = 2.0*v0p2/(sigma[0] + v0p2);
      qdf_add(sigmahihi_h[0],sigmalohi_h[0],sigmahilo_h[0],sigmalolo_h[0],
               v0p2hihi,      v0p2lohi,      v0p2hilo,      v0p2lolo,
               &acchihi,      &acclohi,      &acchilo,      &acclolo);
      qdf_div( v0p2hihi,      v0p2lohi,      v0p2hilo,      v0p2lolo,
                acchihi,       acclohi,       acchilo,       acclolo,
              &betahihi_h[L],&betalohi_h[L],&betahilo_h[L],&betalolo_h[L]);
      qdf_mlt_d(&betahihi_h[L],&betalohi_h[L],
                &betahilo_h[L],&betalolo_h[L],2.0);
      sigmahihi_h[0] = v0hihi;
      sigmalohi_h[0] = v0lohi;
      sigmahilo_h[0] = v0hilo;
      sigmalolo_h[0] = v0lolo;                // v0 needed for normalization
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
           << betahihi_h[L] << "  " << betalohi_h[L] << endl
           << "          "
           << betahilo_h[L] << "  " << betalolo_h[L] << endl;
   }
   cudaMemcpy(&betahihi_d[L],&betahihi_h[L],sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&betalohi_d[L],&betalohi_h[L],sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&betahilo_d[L],&betahilo_h[L],sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&betalolo_d[L],&betalolo_h[L],sizeof(double),
              cudaMemcpyHostToDevice);

   if(!done)  // normalization needed
   {
      // (sigmahi_h, sigmalo_h) has the values for (v0hi, v0lo).
      cudaMemcpy(sigmahihi_d,sigmahihi_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(sigmalohi_d,sigmalohi_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(sigmahilo_d,sigmahilo_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(sigmalolo_d,sigmalolo_h,sizeof(double),
                 cudaMemcpyHostToDevice);

      if(verbose)
      {
         cout << "-> launching " << nblocks << " blocks of "
              << szt << " threads to normalize ..." << endl;
         cout << "   nrows1 : " << nrows1
              << "  rowidx : " << rowidx << "  nVrows : " << nVrows << endl;
      }
      cudaEventRecord(start);
      dbl4_normalize<<<nblocks,szt>>>
         (nrows1,szt,
          &Ahihi_d[rowidx+1],&Alohi_d[rowidx+1],
          &Ahilo_d[rowidx+1],&Alolo_d[rowidx+1],
          sigmahihi_d,sigmalohi_d,sigmahilo_d,sigmalolo_d,
          &Vhihi_d[L*nVrows+L+1],&Vlohi_d[L*nVrows+L+1],
          &Vhilo_d[L*nVrows+L+1],&Vlolo_d[L*nVrows+L+1]);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds,start,stop);
      *lapms += milliseconds;
      flopcount_dbl_normalize(nblocks,szt,div);
   }
   if(verbose)
   {
      const size_t szhouse = nVrows*sizeof(double);

      cudaMemcpy(&betahihi_h[L],&betahihi_d[L],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&betalohi_h[L],&betalohi_d[L],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&betahilo_h[L],&betahilo_d[L],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&betalolo_h[L],&betalolo_d[L],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vhihi_h,&Vhihi_d[L*nVrows],szhouse,cudaMemcpyDeviceToHost);
      cudaMemcpy(vlohi_h,&Vlohi_d[L*nVrows],szhouse,cudaMemcpyDeviceToHost);
      cudaMemcpy(vhilo_h,&Vhilo_d[L*nVrows],szhouse,cudaMemcpyDeviceToHost);
      cudaMemcpy(vlolo_h,&Vlolo_d[L*nVrows],szhouse,cudaMemcpyDeviceToHost);
      cout << scientific << setprecision(16)
           << "beta[" << colidx << "] : "
           << betahihi_h[L] << "  " << betalohi_h[L] << endl
           << "           "
           << betahilo_h[L] << "  " << betalolo_h[L] << endl;
      for(int i=0; i<nVrows; i++)
         cout << "v[" << i << "] : "
              << vhihi_h[i] << "  " << vlohi_h[i] << endl
              << "       "
              << vhilo_h[i] << "  " << vlolo_h[i] << endl;
   }
}

void GPU_cmplx4_large_house
 ( int nrows, int ncols, int szt, int nbt,
   int colidx, int nrows1, int k, int L,
   double *Arehihi_h, double *Arelohi_h, double *Arehilo_h, double *Arelolo_h,
   double *Aimhihi_h, double *Aimlohi_h, double *Aimhilo_h, double *Aimlolo_h,
   double *Arehihi_d, double *Arelohi_d, double *Arehilo_d, double *Arelolo_d,
   double *Aimhihi_d, double *Aimlohi_d, double *Aimhilo_d, double *Aimlolo_d,
   double *vrehihi_h, double *vrelohi_h, double *vrehilo_h, double *vrelolo_h,
   double *vimhihi_h, double *vimlohi_h, double *vimhilo_h, double *vimlolo_h,
   double *Vrehihi_d, double *Vrelohi_d, double *Vrehilo_d, double *Vrelolo_d,
   double *Vimhihi_d, double *Vimlohi_d, double *Vimhilo_d, double *Vimlolo_d,
   double *betahihi_h, double *betalohi_h,
   double *betahilo_h, double *betalolo_h,
   double *betahihi_d, double *betalohi_d,
   double *betahilo_d, double *betalolo_d,
   double *sumshihi_h, double *sumslohi_h,
   double *sumshilo_h, double *sumslolo_h,
   double *sumshihi_d, double *sumslohi_d,
   double *sumshilo_d, double *sumslolo_d,
   double *sigmahihi_h, double *sigmalohi_h,
   double *sigmahilo_h, double *sigmalolo_h,
   double *sigmahihi_d, double *sigmalohi_d,
   double *sigmahilo_d, double *sigmalolo_d,
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
         vrehihi_h[i] = 0.0;
         vrelohi_h[i] = 0.0;
         vrehilo_h[i] = 0.0;
         vrelolo_h[i] = 0.0;
         vimhihi_h[i] = 0.0;
         vimlohi_h[i] = 0.0;
         vimhilo_h[i] = 0.0;
         vimlolo_h[i] = 0.0;
      }
   }
   vrehihi_h[L] = 1.0;                    // set one on the diagonal
   vrelohi_h[L] = 0.0;
   vrehilo_h[L] = 0.0;
   vrelolo_h[L] = 0.0;
   vimhihi_h[L] = 0.0;
   vimlohi_h[L] = 0.0;
   vimhilo_h[L] = 0.0;
   vimlolo_h[L] = 0.0;

   cudaMemcpy(&Vrehihi_d[L*nVrows],vrehihi_h,(L+1)*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&Vrelohi_d[L*nVrows],vrelohi_h,(L+1)*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&Vrehilo_d[L*nVrows],vrehilo_h,(L+1)*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&Vrelolo_d[L*nVrows],vrelolo_h,(L+1)*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&Vimhihi_d[L*nVrows],vimhihi_h,(L+1)*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&Vimlohi_d[L*nVrows],vimlohi_h,(L+1)*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&Vimhilo_d[L*nVrows],vimhilo_h,(L+1)*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&Vimlolo_d[L*nVrows],vimlolo_h,(L+1)*sizeof(double),
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
      sumshihi_h[i] = 0.0;
      sumslohi_h[i] = 0.0;
      sumshilo_h[i] = 0.0;
      sumslolo_h[i] = 0.0;
   }
   cudaMemcpy(sumshihi_d,sumshihi_h,nblocks*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(sumslohi_d,sumslohi_h,nblocks*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(sumshilo_d,sumshilo_h,nblocks*sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(sumslolo_d,sumslolo_h,nblocks*sizeof(double),
              cudaMemcpyHostToDevice);

   cudaEventRecord(start);
   cmplx4_large_sum_of_squares<<<nblocks,szt>>>
      (&Arehihi_d[rowidx+1],&Arelohi_d[rowidx+1],
       &Arehilo_d[rowidx+1],&Arelolo_d[rowidx+1],
       &Aimhihi_d[rowidx+1],&Aimlohi_d[rowidx+1],
       &Aimhilo_d[rowidx+1],&Aimlolo_d[rowidx+1],
       sumshihi_d,sumslohi_d,sumshilo_d,sumslolo_d,nrows1,szt,sztLog2);
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
   dbl4_sum_accumulator<<<1,nblocks>>>
      (sumshihi_d,sumslohi_d,sumshilo_d,sumslolo_d,
       nblocks,nblLog2,sigmahihi_d,sigmalohi_d,sigmahilo_d,sigmalolo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_dbl_sum_accumulator(nblocks,nblLog2,add);

   cudaMemcpy(sigmahihi_h,sigmahihi_d,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(sigmalohi_h,sigmalohi_d,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(sigmahilo_h,sigmahilo_d,sizeof(double),cudaMemcpyDeviceToHost);
   cudaMemcpy(sigmalolo_h,sigmalolo_d,sizeof(double),cudaMemcpyDeviceToHost);

   bool done = false;

   if((sigmahihi_h[0] == 0.0) && (sigmalohi_h[0] == 0.0) &&
      (sigmahilo_h[0] == 0.0) && (sigmalolo_h[0] == 0.0))
   {
      betahihi_h[L] = 0.0; betalohi_h[L] = 0.0;
      betahilo_h[L] = 0.0; betalolo_h[L] = 0.0; done = true;
      if(verbose)
         cout << "Zero sigma value encountered." << endl;
   }
   else // beta is computed on the host instead of by one GPU thread
   {
      double acchihi,acclohi,acchilo,acclolo;
      double muhihi,mulohi,muhilo,mulolo;
      double sqrv0hihi,sqrv0lohi,sqrv0hilo,sqrv0lolo;
      double x0rehihi,x0relohi,x0rehilo,x0relolo;
      double x0imhihi,x0imlohi,x0imhilo,x0imlolo;
      double sqrx0hihi,sqrx0lohi,sqrx0hilo,sqrx0lolo;
      double v0rehihi,v0relohi,v0rehilo,v0relolo;
      double v0imhihi,v0imlohi,v0imhilo,v0imlolo;
      double x0radhihi,x0radlohi,x0radhilo,x0radlolo;
      double inv0rehihi,inv0relohi,inv0rehilo,inv0relolo;
      double inv0imhihi,inv0imlohi,inv0imhilo,inv0imlolo;

      cudaMemcpy(&x0rehihi,&Arehihi_d[rowidx],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&x0relohi,&Arelohi_d[rowidx],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&x0rehilo,&Arehilo_d[rowidx],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&x0relolo,&Arelolo_d[rowidx],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&x0imhihi,&Aimhihi_d[rowidx],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&x0imlohi,&Aimlohi_d[rowidx],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&x0imhilo,&Aimhilo_d[rowidx],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&x0imlolo,&Aimlolo_d[rowidx],sizeof(double),
                 cudaMemcpyDeviceToHost);

      // sqrx0 = xre[0]*xre[0] + xim[0]*xim[0];
      qdf_sqr(  x0rehihi,  x0relohi,  x0rehilo,  x0relolo,
              &sqrx0hihi,&sqrx0lohi,&sqrx0hilo,&sqrx0lolo);
      qdf_sqr(x0imhihi,x0imlohi,x0imhilo,x0imlolo,
              &acchihi,&acclohi,&acchilo,&acclolo);
      qdf_inc(&sqrx0hihi,&sqrx0lohi,&sqrx0hilo,&sqrx0lolo,
                 acchihi,   acclohi,   acchilo,   acclolo);
      // x0rad = sqrt(sqrx0);
      qdf_sqrt( sqrx0hihi, sqrx0lohi, sqrx0hilo, sqrx0lolo,
               &x0radhihi,&x0radlohi,&x0radhilo,&x0radlolo);
      // mu = sqrt(sqrx0 + sigma); // norm of the vector x
      qdf_inc(&sqrx0hihi,  &sqrx0lohi,  &sqrx0hilo,  &sqrx0lolo,
              *sigmahihi_h,*sigmalohi_h,*sigmahilo_h,*sigmalolo_h);
      qdf_sqrt(sqrx0hihi,sqrx0lohi,sqrx0hilo,sqrx0lolo,
                 &muhihi,  &mulohi,  &muhilo,  &mulolo);

      if((x0radhihi == 0.0) && (x0radlohi == 0.0) &&
         (x0radhilo == 0.0) && (x0radlolo == 0.0))
      {
         v0rehihi = -muhihi; 
         v0relohi = -mulohi; 
         v0rehilo = -muhilo; 
         v0relolo = -mulolo; 
         v0imhihi = 0.0;
         v0imlohi = 0.0;
         v0imhilo = 0.0;
         v0imlolo = 0.0;
      }
      else // if(x0rad /= 0.0)   // xre[0]/xrad = cos(angle)
      {                          // xim[0]/xrad = sin(angle)
         // mu = mu/x0rad;
         qdf_div(   muhihi,   mulohi,   muhilo,   mulolo,
                 x0radhihi,x0radlohi,x0radhilo,x0radlolo,
                  &acchihi, &acclohi, &acchilo, &acclolo);
         muhihi = acchihi; mulohi = acclohi;
         muhilo = acchilo; mulolo = acclolo;
         // vre[0] = xre[0] - mu*xre[0];
         qdf_mul(  muhihi,  mulohi,  muhilo,  mulolo,
                 x0rehihi,x0relohi,x0rehilo,x0relolo,
                 &acchihi,&acclohi,&acchilo,&acclolo);
         qdf_sub( x0rehihi, x0relohi, x0rehilo, x0relolo,
                   acchihi,  acclohi,  acchilo,  acclolo,
                 &v0rehihi,&v0relohi,&v0rehilo,&v0relolo);
         // vim[0] = xim[0] - mu*xim[0];
         qdf_mul(  muhihi,  mulohi,  muhilo,  mulolo,
                 x0imhihi,x0imlohi,x0imhilo,x0imlolo,
                 &acchihi,&acclohi,&acchilo,&acclolo);
         qdf_sub( x0imhihi, x0imlohi, x0imhilo, x0imlolo,
                   acchihi,  acclohi,  acchilo,  acclolo,
                 &v0imhihi,&v0imlohi,&v0imhilo,&v0imlolo);
      }
      // sqrv0 = vre[0]*vre[0] + vim[0]*vim[0];
      qdf_sqr(  v0rehihi,  v0relohi,  v0rehilo,  v0relolo,
              &sqrv0hihi,&sqrv0lohi,&sqrv0hilo,&sqrv0lolo);
      qdf_sqr(v0imhihi,v0imlohi,v0imhilo,v0imlolo,
              &acchihi,&acclohi,&acchilo,&acclolo);
      qdf_inc(&sqrv0hihi,&sqrv0lohi,&sqrv0hilo,&sqrv0lolo,
                 acchihi,   acclohi,   acchilo,   acclolo);
      // *beta = 2.0*sqrv0/(sigma + sqrv0);
      qdf_inc(sigmahihi_h,sigmalohi_h,sigmahilo_h,sigmalolo_h,
              sqrv0hihi,  sqrv0lohi,  sqrv0hilo,  sqrv0lolo);
      qdf_div( sqrv0hihi,     sqrv0lohi,     sqrv0hilo,     sqrv0lolo,
              *sigmahihi_h,  *sigmalohi_h,  *sigmahilo_h,  *sigmalolo_h,
               &betahihi_h[L],&betalohi_h[L],&betahilo_h[L],&betalolo_h[L]);
      qdf_mlt_d(&betahihi_h[L],&betalohi_h[L],
                &betahilo_h[L],&betalolo_h[L],2.0);
      // inv0re = vre[0]/sqrv0;  // real part of 1/v[0]
      qdf_div(   v0rehihi,   v0relohi,   v0rehilo,   v0relolo,
                sqrv0hihi,  sqrv0lohi,  sqrv0hilo,  sqrv0lolo,
              &inv0rehihi,&inv0relohi,&inv0rehilo,&inv0relolo);
      // inv0im = -vim[0]/sqrv0; // imaginary part of 1/v[0]
      qdf_div(   v0imhihi,   v0imlohi,   v0imhilo,   v0imlolo,
                sqrv0hihi,  sqrv0lohi,  sqrv0hilo,  sqrv0lolo,
              &inv0imhihi,&inv0imlohi,&inv0imhilo,&inv0imlolo);
      qdf_minus(&inv0imhihi,&inv0imlohi,&inv0imhilo,&inv0imlolo);
      *sigmahihi_h = inv0rehihi;
      *sigmalohi_h = inv0relohi;
      *sigmahilo_h = inv0rehilo;
      *sigmalolo_h = inv0relolo;
      betahihi_h[szt] = inv0imhihi;
      betalohi_h[szt] = inv0imlohi;
      betahilo_h[szt] = inv0imhilo;
      betalolo_h[szt] = inv0imlolo;
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
           << betahihi_h[L] << "  " << betalohi_h[L] << endl
           << "          "
           << betahilo_h[L] << "  " << betalolo_h[L] << endl;
   }
   cudaMemcpy(&betahihi_d[L],&betahihi_h[L],sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&betalohi_d[L],&betalohi_h[L],sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&betahilo_d[L],&betahilo_h[L],sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(&betalolo_d[L],&betalolo_h[L],sizeof(double),
              cudaMemcpyHostToDevice);

   if(!done)  // normalization needed
   {
      // (sigmahi_h, sigmalo_h) has the values for (v0rehi, v0relo)
      cudaMemcpy(sigmahihi_d,sigmahihi_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(sigmalohi_d,sigmalohi_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(sigmahilo_d,sigmahilo_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(sigmalolo_d,sigmalolo_h,sizeof(double),
                 cudaMemcpyHostToDevice);
      // (betahi_h[szt], betalo_h[szt]) has the values for (v0imhi, v0rimlo).
      cudaMemcpy(&betahihi_d[szt],&betahihi_h[szt],sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&betalohi_d[szt],&betalohi_h[szt],sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&betahilo_d[szt],&betahilo_h[szt],sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(&betalolo_d[szt],&betalolo_h[szt],sizeof(double),
                 cudaMemcpyHostToDevice);

      if(verbose)
      {
         cout << "-> launching " << nblocks << " blocks of "
              << szt << " threads to normalize ..." << endl;
         cout << "   nrows1 : " << nrows1
              << "  rowidx : " << rowidx << "  nVrows : " << nVrows << endl;
      }
      cudaEventRecord(start);
      cmplx4_normalize<<<nblocks,szt>>>
         (nrows1,szt,&Arehihi_d[rowidx+1],&Arelohi_d[rowidx+1],
                     &Arehilo_d[rowidx+1],&Arelolo_d[rowidx+1],
                     &Aimhihi_d[rowidx+1],&Aimlohi_d[rowidx+1],
                     &Aimhilo_d[rowidx+1],&Aimlolo_d[rowidx+1],
          sigmahihi_d,sigmalohi_d,sigmahilo_d,sigmalolo_d,
          &betahihi_d[szt],&betalohi_d[szt],&betahilo_d[szt],&betalolo_d[szt],
          &Vrehihi_d[L*nVrows+L+1],&Vrelohi_d[L*nVrows+L+1],
          &Vrehilo_d[L*nVrows+L+1],&Vrelolo_d[L*nVrows+L+1],
          &Vimhihi_d[L*nVrows+L+1],&Vimlohi_d[L*nVrows+L+1],
          &Vimhilo_d[L*nVrows+L+1],&Vimlolo_d[L*nVrows+L+1]);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds,start,stop);
      *lapms += milliseconds;
      flopcount_cmplx_normalize(nblocks,szt,add,mul);
   }
   if(verbose)
   {
      const size_t szhouse = nVrows*sizeof(double);

      cudaMemcpy(&betahihi_h[L],&betahihi_d[L],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&betalohi_h[L],&betalohi_d[L],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&betahilo_h[L],&betahilo_d[L],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&betalolo_h[L],&betalolo_d[L],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vrehihi_h,&Vrehihi_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vrelohi_h,&Vrelohi_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vrehilo_h,&Vrehilo_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vrelolo_h,&Vrelolo_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vimhihi_h,&Vimhihi_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vimlohi_h,&Vimlohi_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vimhilo_h,&Vimhilo_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(vimlolo_h,&Vimlolo_d[L*nVrows],szhouse,
                 cudaMemcpyDeviceToHost);
      cout << scientific << setprecision(16)
           << "beta[" << colidx << "] : "
           << betahihi_h[L] << "  " << betalohi_h[L] << endl
           << "          "
           << betahilo_h[L] << "  " << betalolo_h[L] << endl;

      for(int i=0; i<nVrows; i++)
      {
         cout << "v[" << i << "]re : "
              << vrehihi_h[i] << "  " << vrelohi_h[i] << endl
              << "         "
              << vrehilo_h[i] << "  " << vrelolo_h[i] << endl;
         cout << "v[" << i << "]im : "
              << vimhihi_h[i] << "  " << vimlohi_h[i] << endl
              << "         "
              << vimhilo_h[i] << "  " << vimlolo_h[i] << endl;
      }
   }
}

void GPU_dbl4_small_leftRupdate
 ( int nrows, int ncols, int szt, int colidx, int k, int L,
   double *Ahihi_h, double *Alohi_h, double *Ahilo_h, double *Alolo_h,
   double *Ahihi_d, double *Alohi_d, double *Ahilo_d, double *Alolo_d,
   double *Vhihi_d, double *Vlohi_d, double *Vhilo_d, double *Vlolo_d,
   double *betahihi_h, double *betalohi_h,
   double *betahilo_h, double *betalolo_h,
   double *betahihi_d, double *betalohi_d,
   double *betahilo_d, double *betalolo_d,
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
   dbl4_small_leftRupdate<<<1,nrows-colidx>>>
      (nrows,endcol,szt,colidx,
       Ahihi_d,Alohi_d,Ahilo_d,Alolo_d,
       &Vhihi_d[L*nVrows+L],&Vlohi_d[L*nVrows+L],
       &Vhilo_d[L*nVrows+L],&Vlolo_d[L*nVrows+L],
       &betahihi_d[L],&betalohi_d[L],&betahilo_d[L],&betalolo_d[L]);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_dbl_small_leftRupdate(nrows,ncols,szt,colidx,add,mul);

   if(verbose)
   {
      const int dim = nrows*ncols;
      const size_t sznum = dim*sizeof(double);

      cudaMemcpy(Ahihi_h,Ahihi_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Alohi_h,Alohi_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Ahilo_h,Ahilo_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Alolo_h,Alolo_d,sznum,cudaMemcpyDeviceToHost);
      cout << "the matrix after the update :" << endl;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
            cout << "A_d[" << i << "][" << j << "] : "
                 << Ahihi_h[j*nrows+i] << "  "
                 << Alohi_h[j*nrows+i] << endl
                 << "            "
                 << Ahilo_h[j*nrows+i] << "  "
                 << Alolo_h[j*nrows+i] << endl;
   }
}

void GPU_cmplx4_small_leftRupdate
 ( int nrows, int ncols, int szt, int colidx, int k, int L,
   double *Arehihi_h, double *Arelohi_h, double *Arehilo_h, double *Arelolo_h,
   double *Aimhihi_h, double *Aimlohi_h, double *Aimhilo_h, double *Aimlolo_h,
   double *Arehihi_d, double *Arelohi_d, double *Arehilo_d, double *Arelolo_d,
   double *Aimhihi_d, double *Aimlohi_d, double *Aimhilo_d, double *Aimlolo_d,
   double *Vrehihi_d, double *Vrelohi_d, double *Vrehilo_d, double *Vrelolo_d,
   double *Vimhihi_d, double *Vimlohi_d, double *Vimhilo_d, double *Vimlolo_d,
   double *betahihi_h, double *betalohi_h,
   double *betahilo_h, double *betalolo_h,
   double *betahihi_d, double *betalohi_d,
   double *betahilo_d, double *betalolo_d,
   double *lapms, long long int *add, long long int *mul, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int endcol = (k+1)*szt;     // 1 + last column index in tile
   const int nVrows = nrows - k*szt;          // dimension of V matrix

   cudaEventRecord(start);
   cmplx4_small_leftRupdate<<<1,nrows-colidx>>>
      (nrows,endcol,szt,colidx,
       Arehihi_d,Arelohi_d,Arehilo_d,Arelolo_d,
       Aimhihi_d,Aimlohi_d,Aimhilo_d,Aimlolo_d,
       &Vrehihi_d[L*nVrows+L],&Vrelohi_d[L*nVrows+L],
       &Vrehilo_d[L*nVrows+L],&Vrelolo_d[L*nVrows+L],
       &Vimhihi_d[L*nVrows+L],&Vimlohi_d[L*nVrows+L],
       &Vimhilo_d[L*nVrows+L],&Vimlolo_d[L*nVrows+L],
       &betahihi_d[L],&betalohi_d[L],&betahilo_d[L],&betalolo_d[L]);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_cmplx_small_leftRupdate(nrows,ncols,szt,colidx,add,mul);

   if(verbose)
   {
      const int dim = nrows*ncols;
      const size_t sznum = dim*sizeof(double);

      cudaMemcpy(Arehihi_h,Arehihi_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Arelohi_h,Arelohi_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Arehilo_h,Arehilo_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Arelolo_h,Arelolo_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Aimhihi_h,Aimhihi_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Aimlohi_h,Aimlohi_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Aimhilo_h,Aimhilo_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Aimlolo_h,Aimlolo_d,sznum,cudaMemcpyDeviceToHost);
      cout << "the matrix after the update :" << endl;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
         {
            cout << "A_d[" << i << "][" << j << "]re : "
                 << Arehihi_h[j*nrows+i] << "  "
                 << Arelohi_h[j*nrows+i] << endl
                 << "              "
                 << Arehilo_h[j*nrows+i] << "  "
                 << Arelolo_h[j*nrows+i] << endl;
            cout << "A_d[" << i << "][" << j << "]im : "
                 << Aimhihi_h[j*nrows+i] << "  "
                 << Aimlohi_h[j*nrows+i] << endl
                 << "              "
                 << Aimhilo_h[j*nrows+i] << "  "
                 << Aimlolo_h[j*nrows+i] << endl;
         }
   }
}

void GPU_dbl4_medium_leftRupdate
 ( int nrows, int ncols, int szt, int colidx, int k, int L,
   double *Ahihi_h, double *Alohi_h, double *Ahilo_h, double *Alolo_h,
   double *Ahihi_d, double *Alohi_d, double *Ahilo_d, double *Alolo_d,
   double *Vhihi_d, double *Vlohi_d, double *Vhilo_d, double *Vlolo_d,
   double *betahihi_h, double *betalohi_h,
   double *betahilo_h, double *betalolo_h,
   double *betahihi_d, double *betalohi_d,
   double *betahilo_d, double *betalolo_d,
   double *RTdotvhihi_h, double *RTdotvlohi_h,
   double *RTdotvhilo_h, double *RTdotvlolo_h,
   double *RTdotvhihi_d, double *RTdotvlohi_d,
   double *RTdotvhilo_d, double *RTdotvlolo_d,
   double *whihi_h, double *wlohi_h, double *whilo_h, double *wlolo_h,
   double *whihi_d, double *wlohi_d, double *whilo_d, double *wlolo_d,
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
   dbl4_RTdotv<<<nbrblocks,szt>>>
      (nhouse,szt,colidx,RToffset,dimRTdotv,
       Ahihi_d,Alohi_d,Ahilo_d,Alolo_d,
       &Vhihi_d[L*nVrows+L],&Vlohi_d[L*nVrows+L],
       &Vhilo_d[L*nVrows+L],&Vlolo_d[L*nVrows+L],
       RTdotvhihi_d,RTdotvlohi_d,RTdotvhilo_d,RTdotvlolo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *RTvlapms += milliseconds;
   cudaEventRecord(start);
   dbl4_sum_betaRTdotv<<<1,dimRTdotv>>>
      (nhouse,&betahihi_d[L],&betalohi_d[L],&betahilo_d[L],&betalolo_d[L],
       RTdotvhihi_d,RTdotvlohi_d,RTdotvhilo_d,RTdotvlolo_d,
       whihi_d,wlohi_d,whilo_d,wlolo_d);
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
   dbl4_medium_subvbetaRTv<<<nbrblocks,szt>>>
      (nrows,endcol,szt,colidx,
       Ahihi_d,Alohi_d,Ahilo_d,Alolo_d,
       &Vhihi_d[L*nVrows+L],&Vlohi_d[L*nVrows+L],
       &Vhilo_d[L*nVrows+L],&Vlolo_d[L*nVrows+L],
       &betahihi_d[L],&betalohi_d[L],&betahilo_d[L],&betalolo_d[L],
       whihi_d,wlohi_d,whilo_d,wlolo_d);
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

      cudaMemcpy(RTdotvhihi_h,RTdotvhihi_d,szRTdotv,cudaMemcpyDeviceToHost);
      cudaMemcpy(RTdotvlohi_h,RTdotvlohi_d,szRTdotv,cudaMemcpyDeviceToHost);
      cudaMemcpy(RTdotvhilo_h,RTdotvhilo_d,szRTdotv,cudaMemcpyDeviceToHost);
      cudaMemcpy(RTdotvlolo_h,RTdotvlolo_d,szRTdotv,cudaMemcpyDeviceToHost);
      cout << "the matrix R^T dot v : " << endl;
      int ix = 0;
      for(int i=0; i<endcol-colidx; i++)
      {
         for(int j=0; j<nhouse; j++)      // must use nhouse
         {
            cout << "RTdotv[" << i << "][" << j << "] : "
                 << RTdotvhihi_h[ix] << "  "
                 << RTdotvlohi_h[ix] << endl
                 << "               "
                 << RTdotvhilo_h[ix] << "  "
                 << RTdotvlolo_h[ix] << endl;
            ix = ix + 1;
         }
      }
      cudaMemcpy(whihi_h,whihi_d,szbRTv,cudaMemcpyDeviceToHost);
      cudaMemcpy(wlohi_h,wlohi_d,szbRTv,cudaMemcpyDeviceToHost);
      cudaMemcpy(whilo_h,whilo_d,szbRTv,cudaMemcpyDeviceToHost);
      cudaMemcpy(wlolo_h,wlolo_d,szbRTv,cudaMemcpyDeviceToHost);
      cout << "the vector w = beta*R^T*v : " << endl;
      for(int i=0; i<endcol-colidx; i++)
         cout << "w[" << i << "] : "
              << whihi_h[i] << "  " << wlohi_h[i] << endl
              << "       "
              << whilo_h[i] << "  " << wlolo_h[i] << endl;

      cudaMemcpy(Ahihi_h,Ahihi_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Alohi_h,Alohi_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Ahilo_h,Ahilo_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Alolo_h,Alolo_d,sznum,cudaMemcpyDeviceToHost);
      cout << "the matrix after the update :" << endl;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
            cout << "A_d[" << i << "][" << j << "] : "
                 << Ahihi_h[j*nrows+i] << "  "
                 << Alohi_h[j*nrows+i] << endl
                 << "            "
                 << Ahilo_h[j*nrows+i] << "  "
                 << Alolo_h[j*nrows+i] << endl;
   }
}

void GPU_cmplx4_medium_leftRupdate
 ( int nrows, int ncols, int szt, int colidx, int k, int L,
   double *Arehihi_h, double *Arelohi_h, double *Arehilo_h, double *Arelolo_h,
   double *Aimhihi_h, double *Aimlohi_h, double *Aimhilo_h, double *Aimlolo_h,
   double *Arehihi_d, double *Arelohi_d, double *Arehilo_d, double *Arelolo_d,
   double *Aimhihi_d, double *Aimlohi_d, double *Aimhilo_d, double *Aimlolo_d,
   double *Vrehihi_d, double *Vrelohi_d, double *Vrehilo_d, double *Vrelolo_d,
   double *Vimhihi_d, double *Vimlohi_d, double *Vimhilo_d, double *Vimlolo_d,
   double *betahihi_h, double *betalohi_h,
   double *betahilo_h, double *betalolo_h,
   double *betahihi_d, double *betalohi_d,
   double *betahilo_d, double *betalolo_d,
   double *RHdotvrehihi_h, double *RHdotvrelohi_h,
   double *RHdotvrehilo_h, double *RHdotvrelolo_h,
   double *RHdotvimhihi_h, double *RHdotvimlohi_h,
   double *RHdotvimhilo_h, double *RHdotvimlolo_h,
   double *RHdotvrehihi_d, double *RHdotvrelohi_d,
   double *RHdotvrehilo_d, double *RHdotvrelolo_d,
   double *RHdotvimhihi_d, double *RHdotvimlohi_d,
   double *RHdotvimhilo_d, double *RHdotvimlolo_d,
   double *wrehihi_h, double *wrelohi_h, double *wrehilo_h, double *wrelolo_h,
   double *wimhihi_h, double *wimlohi_h, double *wimhilo_h, double *wimlolo_h,
   double *wrehihi_d, double *wrelohi_d, double *wrehilo_d, double *wrelolo_d,
   double *wimhihi_d, double *wimlohi_d, double *wimhilo_d, double *wimlolo_d,
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
   cmplx4_RHdotv<<<nbrblocks,szt>>>
      (nhouse,szt,colidx,RHoffset,dimRHdotv,
       Arehihi_d,Arelohi_d,Arehilo_d,Arelolo_d,
       Aimhihi_d,Aimlohi_d,Aimhilo_d,Aimlolo_d,
       &Vrehihi_d[L*nVrows+L],&Vrelohi_d[L*nVrows+L],
       &Vrehilo_d[L*nVrows+L],&Vrelolo_d[L*nVrows+L],
       &Vimhihi_d[L*nVrows+L],&Vimlohi_d[L*nVrows+L],
       &Vimhilo_d[L*nVrows+L],&Vimlolo_d[L*nVrows+L],
       RHdotvrehihi_d,RHdotvrelohi_d,RHdotvrehilo_d,RHdotvrelolo_d,
       RHdotvimhihi_d,RHdotvimlohi_d,RHdotvimhilo_d,RHdotvimlolo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *RHvlapms += milliseconds;
 */
   cudaEventRecord(start);
   cmplx4_RHdotvre<<<nbrblocks,szt>>>
      (nhouse,szt,colidx,RHoffset,dimRHdotv,
       Arehihi_d,Arelohi_d,Arehilo_d,Arelolo_d,
       Aimhihi_d,Aimlohi_d,Aimhilo_d,Aimlolo_d,
       &Vrehihi_d[L*nVrows+L],&Vrelohi_d[L*nVrows+L],
       &Vrehilo_d[L*nVrows+L],&Vrelolo_d[L*nVrows+L],
       &Vimhihi_d[L*nVrows+L],&Vimlohi_d[L*nVrows+L],
       &Vimhilo_d[L*nVrows+L],&Vimlolo_d[L*nVrows+L],
       RHdotvrehihi_d,RHdotvrelohi_d,RHdotvrehilo_d,RHdotvrelolo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *RHvlapms += milliseconds;

   cudaEventRecord(start);
   cmplx4_RHdotvim<<<nbrblocks,szt>>>
      (nhouse,szt,colidx,RHoffset,dimRHdotv,
       Arehihi_d,Arelohi_d,Arehilo_d,Arelolo_d,
       Aimhihi_d,Aimlohi_d,Aimhilo_d,Aimlolo_d,
       &Vrehihi_d[L*nVrows+L],&Vrelohi_d[L*nVrows+L],
       &Vrehilo_d[L*nVrows+L],&Vrelolo_d[L*nVrows+L],
       &Vimhihi_d[L*nVrows+L],&Vimlohi_d[L*nVrows+L],
       &Vimhilo_d[L*nVrows+L],&Vimlolo_d[L*nVrows+L],
       RHdotvimhihi_d,RHdotvimlohi_d,RHdotvimhilo_d,RHdotvimlolo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *RHvlapms += milliseconds;

   cudaEventRecord(start);
   cmplx4_sum_betaRHdotv<<<1,dimRHdotv>>>
      (nhouse,
          &betahihi_d[L],&betalohi_d[L],&betahilo_d[L],&betalolo_d[L],
       RHdotvrehihi_d,RHdotvrelohi_d,RHdotvrehilo_d,RHdotvrelolo_d,
       RHdotvimhihi_d,RHdotvimlohi_d,RHdotvimhilo_d,RHdotvimlolo_d,
            wrehihi_d,     wrelohi_d,     wrehilo_d,     wrelolo_d,
            wimhihi_d,     wimlohi_d,     wimhilo_d,     wimlolo_d);
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
   cmplx4_medium_subvbetaRHv<<<nbrblocks,szt>>>
      (nrows,endcol,szt,colidx,
       Arehihi_d,Arelohi_d,Arehilo_d,Arelolo_d,
       Aimhihi_d,Aimlohi_d,Aimhilo_d,Aimlolo_d,
       &Vrehihi_d[L*nVrows+L],&Vrelohi_d[L*nVrows+L],
       &Vrehilo_d[L*nVrows+L],&Vrelolo_d[L*nVrows+L],
       &Vimhihi_d[L*nVrows+L],&Vimlohi_d[L*nVrows+L],
       &Vimhilo_d[L*nVrows+L],&Vimlolo_d[L*nVrows+L],
       &betahihi_d[L],&betalohi_d[L],&betahilo_d[L],&betalolo_d[L],
       wrehihi_d,wrelohi_d,wrehilo_d,wrelolo_d,
       wimhihi_d,wimlohi_d,wimhilo_d,wimlolo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *redlapms += milliseconds;
 */
   cudaEventRecord(start);
   cmplx4_medium_subvbetaRHvRe<<<nbrblocks,szt>>>
      (nrows,endcol,szt,colidx,
       Arehihi_d,Arelohi_d,Arehilo_d,Arelolo_d,
       &Vrehihi_d[L*nVrows+L],&Vrelohi_d[L*nVrows+L],
       &Vrehilo_d[L*nVrows+L],&Vrelolo_d[L*nVrows+L],
       &Vimhihi_d[L*nVrows+L],&Vimlohi_d[L*nVrows+L],
       &Vimhilo_d[L*nVrows+L],&Vimlolo_d[L*nVrows+L],
       &betahihi_d[L],&betalohi_d[L],&betahilo_d[L],&betalolo_d[L],
       wrehihi_d,wrelohi_d,wrehilo_d,wrelolo_d,
       wimhihi_d,wimlohi_d,wimhilo_d,wimlolo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *redlapms += milliseconds;

   cudaEventRecord(start);
   cmplx4_medium_subvbetaRHvIm<<<nbrblocks,szt>>>
      (nrows,endcol,szt,colidx,
       Aimhihi_d,Aimlohi_d,Aimhilo_d,Aimlolo_d,
       &Vrehihi_d[L*nVrows+L],&Vrelohi_d[L*nVrows+L],
       &Vrehilo_d[L*nVrows+L],&Vrelolo_d[L*nVrows+L],
       &Vimhihi_d[L*nVrows+L],&Vimlohi_d[L*nVrows+L],
       &Vimhilo_d[L*nVrows+L],&Vimlolo_d[L*nVrows+L],
       &betahihi_d[L],&betalohi_d[L],&betahilo_d[L],&betalolo_d[L],
       wrehihi_d,wrelohi_d,wrehilo_d,wrelolo_d,
       wimhihi_d,wimlohi_d,wimhilo_d,wimlolo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *redlapms += milliseconds;

   flopcount_cmplx_medium_subvbetaRHv(nrows,endcol,szt,colidx,add,mul);

   if(verbose)
   {
      const int dim = nrows*ncols;
      const size_t sznum = dim*sizeof(double);
      const size_t szbRHv = dimRHdotv*sizeof(double);
      const size_t szRHdotv = nVrows*szbRHv;

      cudaMemcpy(RHdotvrehihi_h,RHdotvrehihi_d,szRHdotv,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(RHdotvrelohi_h,RHdotvrelohi_d,szRHdotv,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(RHdotvrehilo_h,RHdotvrehilo_d,szRHdotv,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(RHdotvrelolo_h,RHdotvrelolo_d,szRHdotv,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(RHdotvimhihi_h,RHdotvimhihi_d,szRHdotv,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(RHdotvimlohi_h,RHdotvimlohi_d,szRHdotv,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(RHdotvimhilo_h,RHdotvimhilo_d,szRHdotv,
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(RHdotvimlolo_h,RHdotvimlolo_d,szRHdotv,
                 cudaMemcpyDeviceToHost);
      cout << "the matrix R^H dot v : " << endl;
      int ix = 0;
      for(int i=0; i<endcol-colidx; i++)
      {
         for(int j=0; j<nhouse; j++)      // must use nhouse
         {
            cout << "RHdotv[" << i << "][" << j << "]re : "
                 << RHdotvrehihi_h[ix] << "  "
                 << RHdotvrelohi_h[ix] << endl
                 << "                 "
                 << RHdotvrehilo_h[ix] << "  "
                 << RHdotvrelolo_h[ix] << endl;
            cout << "RHdotv[" << i << "][" << j << "]im : "
                 << RHdotvimhihi_h[ix] << "  "
                 << RHdotvimlohi_h[ix] << endl
                 << "                 "
                 << RHdotvimhilo_h[ix] << "  "
                 << RHdotvimlolo_h[ix] << endl;
            ix = ix + 1;
         }
      }
      cudaMemcpy(wrehihi_h,wrehihi_d,szbRHv,cudaMemcpyDeviceToHost);
      cudaMemcpy(wrelohi_h,wrelohi_d,szbRHv,cudaMemcpyDeviceToHost);
      cudaMemcpy(wrehilo_h,wrehilo_d,szbRHv,cudaMemcpyDeviceToHost);
      cudaMemcpy(wrelolo_h,wrelolo_d,szbRHv,cudaMemcpyDeviceToHost);
      cudaMemcpy(wimhihi_h,wimhihi_d,szbRHv,cudaMemcpyDeviceToHost);
      cudaMemcpy(wimlohi_h,wimlohi_d,szbRHv,cudaMemcpyDeviceToHost);
      cudaMemcpy(wimhilo_h,wimhilo_d,szbRHv,cudaMemcpyDeviceToHost);
      cudaMemcpy(wimlolo_h,wimlolo_d,szbRHv,cudaMemcpyDeviceToHost);
      cout << "the vector w = beta*R^H*v : " << endl;
      for(int i=0; i<endcol-colidx; i++)
      {
         cout << "w[" << i << "]re : "
              << wrehihi_h[i] << "  " << wimlohi_h[i] << endl
              << "         "
              << wrehilo_h[i] << "  " << wimlolo_h[i] << endl;
         cout << "w[" << i << "]im : "
              << wimhihi_h[i] << "  " << wimlohi_h[i] << endl
              << "         "
              << wimhilo_h[i] << "  " << wimlolo_h[i] << endl;
      }
      cudaMemcpy(Arehihi_h,Arehihi_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Arelohi_h,Arelohi_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Arehilo_h,Arehilo_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Arelolo_h,Arelolo_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Aimhihi_h,Aimhihi_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Aimlohi_h,Aimlohi_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Aimhilo_h,Aimhilo_d,sznum,cudaMemcpyDeviceToHost);
      cudaMemcpy(Aimlolo_h,Aimlolo_d,sznum,cudaMemcpyDeviceToHost);
      cout << "the matrix after the update :" << endl;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
         {
            cout << "A_d[" << i << "][" << j << "]re : "
                 << Arehihi_h[j*nrows+i] << "  "
                 << Arelohi_h[j*nrows+i] << endl
                 << "              "
                 << Arehilo_h[j*nrows+i] << "  "
                 << Arelolo_h[j*nrows+i] << endl;
            cout << "A_d[" << i << "][" << j << "]im : "
                 << Aimhihi_h[j*nrows+i] << "  "
                 << Aimlohi_h[j*nrows+i] << endl
                 << "              "
                 << Aimhilo_h[j*nrows+i] << "  "
                 << Aimlolo_h[j*nrows+i] << endl;
         }
   }
}

void GPU_dbl4_medium_VB_to_W
 ( int nrows, int ncols, int szt, int idx,
   double *Vhihi_h, double *Vlohi_h, double *Vhilo_h, double *Vlolo_h,
   double *Vhihi_d, double *Vlohi_d, double *Vhilo_d, double *Vlolo_d,
   double *Whihi_h, double *Wlohi_h, double *Whilo_h, double *Wlolo_h,
   double *Whihi_d, double *Wlohi_d, double *Whilo_d, double *Wlolo_d,
   double *WYThihi_h, double *WYTlohi_h, double *WYThilo_h, double *WYTlolo_h,
   double *WYThihi_d, double *WYTlohi_d, double *WYThilo_d, double *WYTlolo_d,
   double *betahihi_h, double *betalohi_h,
   double *betahilo_h, double *betalolo_h,
   double *betahihi_d, double *betalohi_d,
   double *betahilo_d, double *betalolo_d,
   double *lapms, long long int *add, long long int *mul, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int rowdim = nrows - idx*szt;
   const int nbrblocks1 = (int) ceil(rowdim/((double) szt));

   cudaEventRecord(start);
   dbl4_beta_times_V<<<nbrblocks1,szt>>>
      (rowdim,szt,betahihi_d,betalohi_d,betahilo_d,betalolo_d,
                     Vhihi_d,   Vlohi_d,   Vhilo_d,   Vlolo_d,
                     Whihi_d,   Wlohi_d,   Whilo_d,   Wlolo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_dbl_beta_times_V(rowdim,mul);

   const int nbrblocks2 = (int) ceil(rowdim*rowdim/((double) szt));

   cudaEventRecord(start);
   dbl4_initialize_WYT<<<nbrblocks2,szt>>>
      (rowdim,szt,  Vhihi_d,  Vlohi_d,  Vhilo_d,  Vlolo_d,
                    Whihi_d,  Wlohi_d,  Whilo_d,  Wlolo_d,
                  WYThihi_d,WYTlohi_d,WYThilo_d,WYTlolo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_dbl_initialize_WYT(rowdim,mul);

   for(int j=1; j<szt; j++)
   {
      cudaEventRecord(start);
      dbl4_beta_next_W<<<nbrblocks1,szt>>>
         (rowdim,szt,
          &betahihi_d[j],&betalohi_d[j],&betahilo_d[j],&betalolo_d[j],
          &Vhihi_d[j*rowdim],&Vlohi_d[j*rowdim],
          &Vhilo_d[j*rowdim],&Vlolo_d[j*rowdim],
          &Whihi_d[j*rowdim],&Wlohi_d[j*rowdim],
          &Whilo_d[j*rowdim],&Wlolo_d[j*rowdim],
          WYThihi_d,WYTlohi_d,WYThilo_d,WYTlolo_d);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds,start,stop);
      *lapms += milliseconds;
      flopcount_dbl_beta_next_W(rowdim,add,mul);

      cudaEventRecord(start);
      dbl4_update_WYT<<<nbrblocks2,szt>>>
         (rowdim,szt,
          &Vhihi_d[j*rowdim],&Vlohi_d[j*rowdim],
          &Vhilo_d[j*rowdim],&Vlolo_d[j*rowdim],
          &Whihi_d[j*rowdim],&Wlohi_d[j*rowdim],
          &Whilo_d[j*rowdim],&Wlolo_d[j*rowdim],
          WYThihi_d,WYTlohi_d,WYThilo_d,WYTlolo_d);
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

      cudaMemcpy(betahihi_h,betahihi_d,szbeta,cudaMemcpyDeviceToHost);
      cudaMemcpy(betalohi_h,betalohi_d,szbeta,cudaMemcpyDeviceToHost);
      cudaMemcpy(betahilo_h,betahilo_d,szbeta,cudaMemcpyDeviceToHost);
      cudaMemcpy(betalolo_h,betalolo_d,szbeta,cudaMemcpyDeviceToHost);
      cout << "the betas :" << endl;
      for(int j=0; j<szt; j++)
         cout << "beta[" << j << "] : "
              << betahihi_h[j] << "  " << betalohi_h[j] << endl
              << "          "
              << betahilo_h[j] << "  " << betalolo_h[j] << endl;

      cudaMemcpy(Vhihi_h,Vhihi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Vlohi_h,Vlohi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Vhilo_h,Vhilo_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Vlolo_h,Vlolo_d,szVandW,cudaMemcpyDeviceToHost);
      cout << "the columns of the V matrix :" << endl;
      int ix = 0;
      for(int j=0; j<szt; j++) 
         for(int i=0; i<rowdim; i++) 
         {
            cout << "V[" << i << "][" << j << "] : "
                 << Vhihi_h[ix] << "  " << Vlohi_h[ix] << endl
                 << "          "
                 << Vhilo_h[ix] << "  " << Vlolo_h[ix] << endl;
            ix = ix + 1;
         }

      cudaMemcpy(Whihi_h,Whihi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Wlohi_h,Wlohi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Whilo_h,Whilo_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Wlolo_h,Wlolo_d,szVandW,cudaMemcpyDeviceToHost);
      cout << "the columns of the W matrix :" << endl;
      ix = 0;
      for(int j=0; j<szt; j++) 
         for(int i=0; i<rowdim; i++) 
         {
            cout << "W[" << i << "][" << j << "] : "
                 << Whihi_h[ix] << "  " << Wlohi_h[ix] << endl
                 << "          "
                 << Whilo_h[ix] << "  " << Wlolo_h[ix] << endl;
            ix = ix + 1;
         }

      cudaMemcpy(WYThihi_h,WYThihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYTlohi_h,WYTlohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYThilo_h,WYThilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYTlolo_h,WYTlolo_d,szmat,cudaMemcpyDeviceToHost);
      cout << "the WYT matrix :" << endl;
      ix = 0;
      for(int i=0; i<rowdim; i++) 
         for(int j=0; j<rowdim; j++) 
         {
            cout << "WYT[" << i << "][" << j << "] : "
                 << WYThihi_h[ix] << "  " << WYTlohi_h[ix] << endl
                 << "            "
                 << WYThilo_h[ix] << "  " << WYTlolo_h[ix] << endl;
            ix = ix + 1;
         }
   }
}

void GPU_cmplx4_medium_VB_to_W
 ( int nrows, int ncols, int szt, int idx,
   double *Vrehihi_h, double *Vrelohi_h, double *Vrehilo_h, double *Vrelolo_h,
   double *Vimhihi_h, double *Vimlohi_h, double *Vimhilo_h, double *Vimlolo_h,
   double *Vrehihi_d, double *Vrelohi_d, double *Vrehilo_d, double *Vrelolo_d,
   double *Vimhihi_d, double *Vimlohi_d, double *Vimhilo_d, double *Vimlolo_d,
   double *Wrehihi_h, double *Wrelohi_h, double *Wrehilo_h, double *Wrelolo_h,
   double *Wimhihi_h, double *Wimlohi_h, double *Wimhilo_h, double *Wimlolo_h,
   double *Wrehihi_d, double *Wrelohi_d, double *Wrehilo_d, double *Wrelolo_d,
   double *Wimhihi_d, double *Wimlohi_d, double *Wimhilo_d, double *Wimlolo_d,
   double *WYHrehihi_h, double *WYHrelohi_h,
   double *WYHrehilo_h, double *WYHrelolo_h,
   double *WYHimhihi_h, double *WYHimlohi_h,
   double *WYHimhilo_h, double *WYHimlolo_h,
   double *WYHrehihi_d, double *WYHrelohi_d,
   double *WYHrehilo_d, double *WYHrelolo_d,
   double *WYHimhihi_d, double *WYHimlohi_d,
   double *WYHimhilo_d, double *WYHimlolo_d,
   double *betahihi_h, double *betalohi_h,
   double *betahilo_h, double *betalolo_h,
   double *betahihi_d, double *betalohi_d,
   double *betahilo_d, double *betalolo_d,
   double *lapms, long long int *add, long long int *mul, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int rowdim = nrows - idx*szt;
   const int nbrblocks1 = (int) ceil(rowdim/((double) szt));

   cudaEventRecord(start);
   cmplx4_beta_times_V<<<nbrblocks1,szt>>>
      (rowdim,szt,betahihi_d,betalohi_d,betahilo_d,betalolo_d,
                   Vrehihi_d, Vrelohi_d, Vrehilo_d, Vrelolo_d,
                   Vimhihi_d, Vimlohi_d, Vimhilo_d, Vimlolo_d,
                   Wrehihi_d, Wrelohi_d, Wrehilo_d, Wrelolo_d,
                   Wimhihi_d, Wimlohi_d, Wimhilo_d, Wimlolo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_cmplx_beta_times_V(rowdim,mul);

   const int nbrblocks2 = (int) ceil(rowdim*rowdim/((double) szt));

   cudaEventRecord(start);
   cmplx4_initialize_WYH<<<nbrblocks2,szt>>>
      (rowdim,szt,  Vrehihi_d,  Vrelohi_d,  Vrehilo_d,  Vrelolo_d,
                    Vimhihi_d,  Vimlohi_d,  Vimhilo_d,  Vimlolo_d,
                    Wrehihi_d,  Wrelohi_d,  Wrehilo_d,  Wrelolo_d,
                    Wimhihi_d,  Wimlohi_d,  Wimhilo_d,  Wimlolo_d,
                  WYHrehihi_d,WYHrelohi_d,WYHrehilo_d,WYHrelolo_d,
                  WYHimhihi_d,WYHimlohi_d,WYHimhilo_d,WYHimlolo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_cmplx_initialize_WYH(rowdim,add,mul);

   for(int j=1; j<szt; j++)
   {
      cudaEventRecord(start);
      cmplx4_beta_next_W<<<nbrblocks1,szt>>>
         (rowdim,szt,
          &betahihi_d[j],&betalohi_d[j],&betahilo_d[j],&betalolo_d[j],
          &Vrehihi_d[j*rowdim],&Vrelohi_d[j*rowdim],
          &Vrehilo_d[j*rowdim],&Vrelolo_d[j*rowdim],
          &Vimhihi_d[j*rowdim],&Vimlohi_d[j*rowdim],
          &Vimhilo_d[j*rowdim],&Vimlolo_d[j*rowdim],
          &Wrehihi_d[j*rowdim],&Wrelohi_d[j*rowdim],
          &Wrehilo_d[j*rowdim],&Wrelolo_d[j*rowdim],
          &Wimhihi_d[j*rowdim],&Wimlohi_d[j*rowdim],
          &Wimhilo_d[j*rowdim],&Wimlolo_d[j*rowdim],
          WYHrehihi_d,WYHrelohi_d,WYHrehilo_d,WYHrelolo_d,
          WYHimhihi_d,WYHimlohi_d,WYHimhilo_d,WYHimlolo_d);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds,start,stop);
      *lapms += milliseconds;
      flopcount_cmplx_beta_next_W(rowdim,add,mul);

      cudaEventRecord(start);
      cmplx4_update_WYH<<<nbrblocks2,szt>>>
         (rowdim,szt,&Vrehihi_d[j*rowdim],&Vrelohi_d[j*rowdim],
                     &Vrehilo_d[j*rowdim],&Vrelolo_d[j*rowdim],
                     &Vimhihi_d[j*rowdim],&Vimlohi_d[j*rowdim],
                     &Vimhilo_d[j*rowdim],&Vimlolo_d[j*rowdim],
                     &Wrehihi_d[j*rowdim],&Wrelohi_d[j*rowdim],
                     &Wrehilo_d[j*rowdim],&Wrelolo_d[j*rowdim],
                     &Wimhihi_d[j*rowdim],&Wimlohi_d[j*rowdim],
                     &Wimhilo_d[j*rowdim],&Wimlolo_d[j*rowdim],
                     WYHrehihi_d,WYHrelohi_d,WYHrehilo_d,WYHrelolo_d,
                     WYHimhihi_d,WYHimlohi_d,WYHimhilo_d,WYHimlolo_d);
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

      cudaMemcpy(betahihi_h,betahihi_d,szbeta,cudaMemcpyDeviceToHost);
      cudaMemcpy(betalohi_h,betalohi_d,szbeta,cudaMemcpyDeviceToHost);
      cudaMemcpy(betahilo_h,betahilo_d,szbeta,cudaMemcpyDeviceToHost);
      cudaMemcpy(betalolo_h,betalolo_d,szbeta,cudaMemcpyDeviceToHost);
      cout << "the betas :" << endl;
      for(int j=0; j<szt; j++)
         cout << "beta[" << j << "] : "
              << betahihi_h[j] << "  " << betalohi_h[j] << endl
              << "          "
              << betahilo_h[j] << "  " << betalolo_h[j] << endl;

      cudaMemcpy(Vrehihi_h,Vrehihi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Vrelohi_h,Vrelohi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Vrehilo_h,Vrehilo_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Vrelolo_h,Vrelolo_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Vimhihi_h,Vimhihi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Vimlohi_h,Vimlohi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Vimhilo_h,Vimhilo_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Vimlolo_h,Vimlolo_d,szVandW,cudaMemcpyDeviceToHost);
      cout << "the columns of the V matrix :" << endl;
      int ix = 0;
      for(int j=0; j<szt; j++) 
         for(int i=0; i<rowdim; i++) 
         {
            cout << "V[" << i << "][" << j << "]re : "
                 << Vrehihi_h[ix] << "  " << Vrelohi_h[ix] << endl
                 << "            "
                 << Vrehilo_h[ix] << "  " << Vrelolo_h[ix] << endl;
            cout << "V[" << i << "][" << j << "]im : "
                 << Vimhihi_h[ix] << "  " << Vimlohi_h[ix] << endl
                 << "            "
                 << Vimhilo_h[ix] << "  " << Vimlolo_h[ix] << endl;
            ix = ix + 1;
         }

      cudaMemcpy(Wrehihi_h,Wrehihi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Wrelohi_h,Wrelohi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Wrehilo_h,Wrehilo_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Wrelolo_h,Wrelolo_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Wimhihi_h,Wimhihi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Wimlohi_h,Wimlohi_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Wimhilo_h,Wimhilo_d,szVandW,cudaMemcpyDeviceToHost);
      cudaMemcpy(Wimlolo_h,Wimlolo_d,szVandW,cudaMemcpyDeviceToHost);
      cout << "the columns of the W matrix :" << endl;
      ix = 0;
      for(int j=0; j<szt; j++) 
         for(int i=0; i<rowdim; i++) 
         {
            cout << "W[" << i << "][" << j << "]re : "
                 << Wrehihi_h[ix] << "  " << Wrelohi_h[ix] << endl
                 << "            "
                 << Wrehilo_h[ix] << "  " << Wrelolo_h[ix] << endl;
            cout << "W[" << i << "][" << j << "]im : "
                 << Wimhihi_h[ix] << "  " << Wimlohi_h[ix] << endl
                 << "            "
                 << Wimhilo_h[ix] << "  " << Wimlolo_h[ix] << endl;
            ix = ix + 1;
         }

      cudaMemcpy(WYHrehihi_h,WYHrehihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYHrelohi_h,WYHrelohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYHrehilo_h,WYHrehilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYHrelolo_h,WYHrelolo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYHimhihi_h,WYHimhihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYHimlohi_h,WYHimlohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYHimhilo_h,WYHimhilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYHimlolo_h,WYHimlolo_d,szmat,cudaMemcpyDeviceToHost);
      cout << "the WYT matrix :" << endl;
      ix = 0;
      for(int i=0; i<rowdim; i++) 
         for(int j=0; j<rowdim; j++) 
         {
            cout << "WYH[" << i << "][" << j << "]re : "
                 << WYHrehihi_h[ix] << "  " << WYHrelohi_h[ix] << endl
                 << "              "
                 << WYHrehilo_h[ix] << "  " << WYHrelolo_h[ix] << endl;
            cout << "WYH[" << i << "][" << j << "]im : "
                 << WYHimhihi_h[ix] << "  " << WYHimlohi_h[ix] << endl
                 << "              "
                 << WYHimhilo_h[ix] << "  " << WYHimlolo_h[ix] << endl;
            ix = ix + 1;
         }
   }
}

void GPU_dbl4_small_WYT
 ( int nrows, int szt,
   double *Whihi_d, double *Wlohi_d, double *Whilo_d, double *Wlolo_d,
   double *Yhihi_d, double *Ylohi_d, double *Yhilo_d, double *Ylolo_d,
   double *WYThihi_d, double *WYTlohi_d, double *WYThilo_d, double *WYTlolo_d,
   double *WYThihi_h, double *WYTlohi_h, double *WYThilo_h, double *WYTlolo_h,
   double *lapms, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int nbrblocks = (int) ceil(nrows*nrows/((double) szt));

   cudaEventRecord(start);
   dbl4_small_WYT<<<nbrblocks,szt>>>
      (nrows,szt,  Whihi_d,  Wlohi_d,  Whilo_d,  Wlolo_d,
                   Yhihi_d,  Ylohi_d,  Yhilo_d,  Ylolo_d,
                 WYThihi_d,WYTlohi_d,WYThilo_d,WYTlolo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   // flopcount_dbl_small_WYT(nrows,szt,add,mul);

   if(verbose)
   {
      const size_t szmat = nrows*nrows*sizeof(double);

      cudaMemcpy(WYThihi_h,WYThihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYTlohi_h,WYTlohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYThilo_h,WYThilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYTlolo_h,WYTlolo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the WYT matrix :" << endl;
      int ix = 0;
      for(int i=0; i<nrows; i++) 
         for(int j=0; j<nrows; j++) 
         {
            cout << "WYT[" << i << "][" << j << "] : "
                 << WYThihi_h[ix] << "  " << WYTlohi_h[ix] << endl
                 << "            "
                 << WYThilo_h[ix] << "  " << WYTlolo_h[ix] << endl;
            ix = ix + 1;
         }
   }
}

void GPU_cmplx4_small_WYH
 ( int nrows, int szt,
   double *Wrehihi_d, double *Wrelohi_d, double *Wrehilo_d, double *Wrelolo_d,
   double *Wimhihi_d, double *Wimlohi_d, double *Wimhilo_d, double *Wimlolo_d,
   double *Yrehihi_d, double *Yrelohi_d, double *Yrehilo_d, double *Yrelolo_d,
   double *Yimhihi_d, double *Yimlohi_d, double *Yimhilo_d, double *Yimlolo_d,
   double *WYTrehihi_d, double *WYTrelohi_d,
   double *WYTrehilo_d, double *WYTrelolo_d,
   double *WYTimhihi_d, double *WYTimlohi_d,
   double *WYTimhilo_d, double *WYTimlolo_d,
   double *WYTrehihi_h, double *WYTrelohi_h,
   double *WYTrehilo_h, double *WYTrelolo_h,
   double *WYTimhihi_h, double *WYTimlohi_h,
   double *WYTimhilo_h, double *WYTimlolo_h,
   double *lapms, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int nbrblocks = (int) ceil(nrows*nrows/((double) szt));

   cudaEventRecord(start);
   cmplx4_small_WYH<<<nbrblocks,szt>>>
      (nrows,szt, Wrehihi_d,  Wrelohi_d,  Wrehilo_d,  Wrelolo_d,
                  Wimhihi_d,  Wimlohi_d,  Wimhilo_d,  Wimlolo_d,
                  Yrehihi_d,  Yrelohi_d,  Yrehilo_d,  Yrelolo_d,
                  Yimhihi_d,  Yimlohi_d,  Yimhilo_d,  Yimlolo_d,
                WYTrehihi_d,WYTrelohi_d,WYTrehilo_d,WYTrelolo_d,
                WYTimhihi_d,WYTimlohi_d,WYTimhilo_d,WYTimlolo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   // flopcount_cmplx_small_WYH(nrows,szt,add,mul);

   if(verbose)
   {
      const size_t szmat = nrows*nrows*sizeof(double);

      cudaMemcpy(WYTrehihi_h,WYTrehihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYTrelohi_h,WYTrelohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYTrehilo_h,WYTrehilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYTrelolo_h,WYTrelolo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYTimhihi_h,WYTimhihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYTimlohi_h,WYTimlohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYTimhilo_h,WYTimhilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(WYTimlolo_h,WYTimlolo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the WYT matrix :" << endl;
      int ix = 0;
      for(int i=0; i<nrows; i++) 
         for(int j=0; j<nrows; j++) 
         {
            cout << "WYT[" << i << "][" << j << "]re : "
                 << WYTrehihi_h[ix] << "  " << WYTrelohi_h[ix] << endl
                 << "              "
                 << WYTrehilo_h[ix] << "  " << WYTrelolo_h[ix] << endl;
            cout << "WYT[" << i << "][" << j << "]im : "
                 << WYTimhihi_h[ix] << "  " << WYTimlohi_h[ix] << endl
                 << "              "
                 << WYTimhilo_h[ix] << "  " << WYTimlolo_h[ix] << endl;
            ix = ix + 1;
         }
   }
}

void GPU_dbl4_small_YWT
 ( int nrows, int szt, int idx,
   double *Yhihi_d, double *Ylohi_d, double *Yhilo_d, double *Ylolo_d,
   double *Whihi_d, double *Wlohi_d, double *Whilo_d, double *Wlolo_d,
   double *YWThihi_d, double *YWTlohi_d, double *YWThilo_d, double *YWTlolo_d,
   double *YWThihi_h, double *YWTlohi_h, double *YWThilo_h, double *YWTlolo_h,
   double *lapms, long long int *add, long long int *mul, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int rowdim = nrows - idx*szt;
   int nbrblocks = (int) ceil(rowdim*rowdim/((double) szt));

   cudaEventRecord(start);
   dbl4_small_WYT<<<nbrblocks,szt>>>
      (rowdim,szt,  Yhihi_d,  Ylohi_d,  Yhilo_d,  Ylolo_d,
                    Whihi_d,  Wlohi_d,  Whilo_d,  Wlolo_d,
                  YWThihi_d,YWTlohi_d,YWThilo_d,YWTlolo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_dbl_small_WYT(rowdim,szt,add,mul);

   if(verbose)
   {
      const size_t szmat = rowdim*rowdim*sizeof(double);

      cudaMemcpy(YWThihi_h,YWThihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTlohi_h,YWTlohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWThilo_h,YWThilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTlolo_h,YWTlolo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the YWT matrix :" << endl;
      int ix = 0;
      for(int i=0; i<rowdim; i++) 
         for(int j=0; j<rowdim; j++) 
         {
            cout << "YWT[" << i << "][" << j << "] : "
                 << YWThihi_h[ix] << "  " << YWTlohi_h[ix] << endl
                 << "            "
                 << YWThilo_h[ix] << "  " << YWTlolo_h[ix] << endl;
            ix = ix + 1;
         }
   }
}

void GPU_cmplx4_small_YWH
 ( int nrows, int szt, int idx,
   double *Yrehihi_d, double *Yrelohi_d, double *Yrehilo_d, double *Yrelolo_d,
   double *Yimhihi_d, double *Yimlohi_d, double *Yimhilo_d, double *Yimlolo_d,
   double *Wrehihi_d, double *Wrelohi_d, double *Wrehilo_d, double *Wrelolo_d,
   double *Wimhihi_d, double *Wimlohi_d, double *Wimhilo_d, double *Wimlolo_d,
   double *YWTrehihi_d, double *YWTrelohi_d,
   double *YWTrehilo_d, double *YWTrelolo_d,
   double *YWTimhihi_d, double *YWTimlohi_d,
   double *YWTimhilo_d, double *YWTimlolo_d,
   double *YWTrehihi_h, double *YWTrelohi_h,
   double *YWTrehilo_h, double *YWTrelolo_h,
   double *YWTimhihi_h, double *YWTimlohi_h,
   double *YWTimhilo_h, double *YWTimlolo_h,
   double *lapms, long long int *add, long long int *mul, bool verbose )
{
   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;
   const int rowdim = nrows - idx*szt;
   int nbrblocks = (int) ceil(rowdim*rowdim/((double) szt));

   cudaEventRecord(start);
   cmplx4_small_WYH<<<nbrblocks,szt>>>
      (rowdim,szt,  Yrehihi_d,  Yrelohi_d,  Yrehilo_d,  Yrelolo_d,
                    Yimhihi_d,  Yimlohi_d,  Yimhilo_d,  Yimlolo_d,
                    Wrehihi_d,  Wrelohi_d,  Wrehilo_d,  Wrelolo_d,
                    Wimhihi_d,  Wimlohi_d,  Wimhilo_d,  Wimlolo_d,
                  YWTrehihi_d,YWTrelohi_d,YWTrehilo_d,YWTrelolo_d,
                  YWTimhihi_d,YWTimlohi_d,YWTimhilo_d,YWTimlolo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_cmplx_small_WYH(rowdim,szt,add,mul);

   if(verbose)
   {
      const size_t szmat = rowdim*rowdim*sizeof(double);

      cudaMemcpy(YWTrehihi_h,YWTrehihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTrelohi_h,YWTrelohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTrehilo_h,YWTrehilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTrelolo_h,YWTrelolo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTimhihi_h,YWTimhihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTimlohi_h,YWTimlohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTimhilo_h,YWTimhilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTimlolo_h,YWTimlolo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the YWT matrix :" << endl;
      int ix = 0;
      for(int i=0; i<rowdim; i++) 
         for(int j=0; j<rowdim; j++) 
         {
            cout << "YWT[" << i << "][" << j << "]re : "
                 << YWTrehihi_h[ix] << "  " << YWTrelohi_h[ix] << endl
                 << "              "
                 << YWTrehilo_h[ix] << "  " << YWTrelolo_h[ix] << endl;
            cout << "YWT[" << i << "][" << j << "]im : "
                 << YWTimhihi_h[ix] << "  " << YWTimlohi_h[ix] << endl
                 << "              "
                 << YWTimhilo_h[ix] << "  " << YWTimlolo_h[ix] << endl;
            ix = ix + 1;
         }
   }
}

void GPU_dbl4_small_QWYT
 ( int dim, int szt, int idx,
   double *Qhihi_d, double *Qlohi_d, double *Qhilo_d, double *Qlolo_d,
   double *WYThihi_d, double *WYTlohi_d, double *WYThilo_d, double *WYTlolo_d,
   double *QWYThihi_d, double *QWYTlohi_d,
   double *QWYThilo_d, double *QWYTlolo_d,
   double *QWYThihi_h, double *QWYTlohi_h,
   double *QWYThilo_h, double *QWYTlolo_h,
   double *Qhihi_h, double *Qlohi_h, double *Qhilo_h, double *Qlolo_h,
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

      cudaMemcpy(Qhihi_h,Qhihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qlohi_h,Qlohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qhilo_h,Qhilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qlolo_h,Qlolo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the Q matrix :" << endl;
      int ix = 0;
      for(int i=0; i<dim; i++) 
         for(int j=0; j<dim; j++) 
         {
            cout << "Q[" << i << "][" << j << "] : "
                 << Qhihi_h[ix] << "  " << Qlohi_h[ix] << endl
                 << "          "
                 << Qhilo_h[ix] << "  " << Qlolo_h[ix] << endl;
            ix = ix + 1;
         }
   }

   cudaEventRecord(start);
   dbl4_small_QWYT<<<nbrblocks,szt>>>
      (dim,rowdim,szt,coloff,
          Qhihi_d,   Qlohi_d,   Qhilo_d,   Qlolo_d,
        WYThihi_d, WYTlohi_d, WYThilo_d, WYTlolo_d,
       QWYThihi_d,QWYTlohi_d,QWYThilo_d,QWYTlolo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_dbl_small_QWYT(dim,rowdim,szt,coloff,add,mul);

   if(verbose)
   {
      const size_t szmat = dim*rowdim*sizeof(double);

      cudaMemcpy(QWYThihi_h,QWYThihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(QWYTlohi_h,QWYTlohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(QWYThilo_h,QWYThilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(QWYTlolo_h,QWYTlolo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the QWYT matrix :" << endl;
      int ix = 0;
      for(int i=0; i<dim; i++) 
         for(int j=0; j<rowdim; j++) 
         {
            cout << "QWYT[" << i << "][" << j << "] : "
                 << QWYThihi_h[ix] << "  " << QWYTlohi_h[ix] << endl
                 << "             "
                 << QWYThilo_h[ix] << "  " << QWYTlolo_h[ix] << endl;
            ix = ix + 1;
         }
   }
}

void GPU_cmplx4_small_QWYH
 ( int dim, int szt, int idx,
   double *Qrehihi_d, double *Qrelohi_d, double *Qrehilo_d, double *Qrelolo_d,
   double *Qimhihi_d, double *Qimlohi_d, double *Qimhilo_d, double *Qimlolo_d,
   double *WYTrehihi_d, double *WYTrelohi_d,
   double *WYTrehilo_d, double *WYTrelolo_d,
   double *WYTimhihi_d, double *WYTimlohi_d,
   double *WYTimhilo_d, double *WYTimlolo_d,
   double *QWYTrehihi_d, double *QWYTrelohi_d,
   double *QWYTrehilo_d, double *QWYTrelolo_d,
   double *QWYTimhihi_d, double *QWYTimlohi_d,
   double *QWYTimhilo_d, double *QWYTimlolo_d,
   double *QWYTrehihi_h, double *QWYTrelohi_h,
   double *QWYTrehilo_h, double *QWYTrelolo_h,
   double *QWYTimhihi_h, double *QWYTimlohi_h,
   double *QWYTimhilo_h, double *QWYTimlolo_h,
   double *Qrehihi_h, double *Qrelohi_h, double *Qrehilo_h, double *Qrelolo_h,
   double *Qimhihi_h, double *Qimlohi_h, double *Qimhilo_h, double *Qimlolo_h,
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

      cudaMemcpy(Qrehihi_h,Qrehihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qrelohi_h,Qrelohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qrehilo_h,Qrehilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qrelolo_h,Qrelolo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qimhihi_h,Qimhihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qimlohi_h,Qimlohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qimhilo_h,Qimhilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qimlolo_h,Qimlolo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the Q matrix :" << endl;
      int ix = 0;
      for(int i=0; i<dim; i++) 
         for(int j=0; j<dim; j++) 
         {
            cout << "Q[" << i << "][" << j << "]re : "
                 << Qrehihi_h[ix] << "  " << Qrelohi_h[ix] << endl
                 << "          "
                 << Qrehilo_h[ix] << "  " << Qrelolo_h[ix] << endl;
            cout << "Q[" << i << "][" << j << "]im : "
                 << Qimhihi_h[ix] << "  " << Qimlohi_h[ix] << endl
                 << "          "
                 << Qimhilo_h[ix] << "  " << Qimlolo_h[ix] << endl;
            ix = ix + 1;
         }
   }

   cudaEventRecord(start);
   cmplx4_small_QWYH<<<nbrblocks,szt>>>
      (dim,rowdim,szt,coloff,
          Qrehihi_d,   Qrelohi_d,   Qrehilo_d,   Qrelolo_d,
          Qimhihi_d,   Qimlohi_d,   Qimhilo_d,   Qimlolo_d,
        WYTrehihi_d, WYTrelohi_d, WYTrehilo_d, WYTrelolo_d,
        WYTimhihi_d, WYTimlohi_d, WYTimhilo_d, WYTimlolo_d,
       QWYTrehihi_d,QWYTrelohi_d,QWYTrehilo_d,QWYTrelolo_d,
       QWYTimhihi_d,QWYTimlohi_d,QWYTimhilo_d,QWYTimlolo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_cmplx_small_QWYH(dim,rowdim,szt,coloff,add,mul);

   if(verbose)
   {
      const size_t szmat = dim*rowdim*sizeof(double);

      cudaMemcpy(QWYTrehihi_h,QWYTrehihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(QWYTrelohi_h,QWYTrelohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(QWYTrehilo_h,QWYTrehilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(QWYTrelolo_h,QWYTrelolo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(QWYTimhihi_h,QWYTimhihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(QWYTimlohi_h,QWYTimlohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(QWYTimhilo_h,QWYTimhilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(QWYTimlolo_h,QWYTimlolo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the QWYT matrix :" << endl;
      int ix = 0;
      for(int i=0; i<dim; i++) 
         for(int j=0; j<rowdim; j++) 
         {
            cout << "QWYT[" << i << "][" << j << "]re : "
                 << QWYTrehihi_h[ix] << "  " << QWYTrelohi_h[ix] << endl
                 << "               "
                 << QWYTrehilo_h[ix] << "  " << QWYTrelolo_h[ix] << endl;
            cout << "QWYT[" << i << "][" << j << "]im : "
                 << QWYTimhihi_h[ix] << "  " << QWYTimlohi_h[ix] << endl
                 << "               "
                 << QWYTimhilo_h[ix] << "  " << QWYTimlolo_h[ix] << endl;
            ix = ix + 1;
         }
   }
}

void GPU_dbl4_small_YWTC
 ( int nrows, int ncols, int szt, int idx,
   double *YWThihi_d, double *YWTlohi_d, double *YWThilo_d, double *YWTlolo_d,
   double *Chihi_d, double *Clohi_d, double *Chilo_d, double *Clolo_d,
   double *YWTChihi_d, double *YWTClohi_d,
   double *YWTChilo_d, double *YWTClolo_d,
   double *YWTChihi_h, double *YWTClohi_h,
   double *YWTChilo_h, double *YWTClolo_h,
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

      double *Chihi_h = new double[nrows*ncols];
      double *Clohi_h = new double[nrows*ncols];
      double *Chilo_h = new double[nrows*ncols];
      double *Clolo_h = new double[nrows*ncols];
      const size_t szmat = nrows*ncols*sizeof(double);

      cudaMemcpy(Chihi_h,Chihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Clohi_h,Clohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Chilo_h,Chilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Clolo_h,Clolo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the matrix C : " << endl;
      for(int i=rowoff; i<nrows; i++)
         for(int j=coloff; j<ncols; j++)
            cout << "C_h[" << i << "][" << j << "] : "
                 << Chihi_h[j*nrows+i] << "  "
                 << Clohi_h[j*nrows+i] << endl
                 << "            "
                 << Chilo_h[j*nrows+i] << "  "
                 << Clolo_h[j*nrows+i] << endl;

      free(Chihi_h); free(Clohi_h); free(Chilo_h); free(Clolo_h);
   }

   cudaEventRecord(start);
   dbl4_small_YWTC<<<nbrblocks,szt>>>
      (nrows,ncols,rowdim,coldim,szt,rowoff,coloff,
        YWThihi_d, YWTlohi_d, YWThilo_d, YWTlolo_d,
          Chihi_d,   Clohi_d,   Chilo_d,   Clolo_d,
       YWTChihi_d,YWTClohi_d,YWTChilo_d,YWTClolo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_dbl_small_YWTC(rowdim,coldim,add,mul);

   if(verbose)
   {
      const size_t szmat = nrows*ncols*sizeof(double);

      cudaMemcpy(YWTChihi_h,YWTChihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTClohi_h,YWTClohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTChilo_h,YWTChilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTClolo_h,YWTClolo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the YWTC matrix :" << endl;
      for(int i=rowoff; i<nrows; i++) 
         for(int j=coloff; j<ncols; j++)
            cout << "YWTC[" << i << "][" << j << "] : "
                 << YWTChihi_h[j*nrows + i] << "  "
                 << YWTClohi_h[j*nrows + i] << endl
                 << "             "
                 << YWTChilo_h[j*nrows + i] << "  "
                 << YWTClolo_h[j*nrows + i] << endl;
   }
}

void GPU_cmplx4_small_YWHC
 ( int nrows, int ncols, int szt, int idx,
   double *YWTrehihi_d, double *YWTrelohi_d,
   double *YWTrehilo_d, double *YWTrelolo_d,
   double *YWTimhihi_d, double *YWTimlohi_d,
   double *YWTimhilo_d, double *YWTimlolo_d,
   double *Crehihi_d, double *Crelohi_d, double *Crehilo_d, double *Crelolo_d,
   double *Cimhihi_d, double *Cimlohi_d, double *Cimhilo_d, double *Cimlolo_d,
   double *YWTCrehihi_d, double *YWTCrelohi_d,
   double *YWTCrehilo_d, double *YWTCrelolo_d,
   double *YWTCimhihi_d, double *YWTCimlohi_d,
   double *YWTCimhilo_d, double *YWTCimlolo_d,
   double *YWTCrehihi_h, double *YWTCrelohi_h,
   double *YWTCrehilo_h, double *YWTCrelolo_h,
   double *YWTCimhihi_h, double *YWTCimlohi_h,
   double *YWTCimhilo_h, double *YWTCimlolo_h,
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

      double *Crehihi_h = new double[nrows*ncols];
      double *Crelohi_h = new double[nrows*ncols];
      double *Crehilo_h = new double[nrows*ncols];
      double *Crelolo_h = new double[nrows*ncols];
      double *Cimhihi_h = new double[nrows*ncols];
      double *Cimlohi_h = new double[nrows*ncols];
      double *Cimhilo_h = new double[nrows*ncols];
      double *Cimlolo_h = new double[nrows*ncols];
      const size_t szmat = nrows*ncols*sizeof(double);

      cudaMemcpy(Crehihi_h,Crehihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Crelohi_h,Crelohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Crehilo_h,Crehilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Crelolo_h,Crelolo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Cimhihi_h,Cimhihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Cimlohi_h,Cimlohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Cimhilo_h,Cimhilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Cimlolo_h,Cimlolo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the matrix C : " << endl;
      for(int i=rowoff; i<nrows; i++)
         for(int j=coloff; j<ncols; j++)
         {
            cout << "C_h[" << i << "][" << j << "]re : "
                 << Crehihi_h[j*nrows+i] << "  "
                 << Crelohi_h[j*nrows+i] << endl
                 << "              "
                 << Crehilo_h[j*nrows+i] << "  "
                 << Crelolo_h[j*nrows+i] << endl;
            cout << "C_h[" << i << "][" << j << "]im : "
                 << Cimhihi_h[j*nrows+i] << "  "
                 << Cimlohi_h[j*nrows+i] << endl
                 << "              "
                 << Cimhilo_h[j*nrows+i] << "  "
                 << Cimlolo_h[j*nrows+i] << endl;
         }

      free(Crehihi_h); free(Crelohi_h); free(Crehilo_h); free(Crelolo_h);
      free(Cimhihi_h); free(Cimlohi_h); free(Cimhilo_h); free(Cimlolo_h);
   }
   cudaEventRecord(start);
   cmplx4_small_YWHC<<<nbrblocks,szt>>>
      (nrows,ncols,rowdim,coldim,szt,rowoff,coloff,
        YWTrehihi_d, YWTrelohi_d,  YWTrehilo_d, YWTrelolo_d,
        YWTimhihi_d, YWTimlohi_d,  YWTimhilo_d, YWTimlolo_d,
          Crehihi_d,   Crelohi_d,    Crehilo_d,   Crelolo_d,
          Cimhihi_d,   Cimlohi_d,    Cimhilo_d,   Cimlolo_d,
       YWTCrehihi_d,YWTCrelohi_d, YWTCrehilo_d,YWTCrelolo_d,
       YWTCimhihi_d,YWTCimlohi_d, YWTCimhilo_d,YWTCimlolo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_cmplx_small_YWHC(rowdim,coldim,add,mul);

   if(verbose)
   {
      const size_t szmat = nrows*ncols*sizeof(double);

      cudaMemcpy(YWTCrehihi_h,YWTCrehihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTCrelohi_h,YWTCrelohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTCrehilo_h,YWTCrehilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTCrelolo_h,YWTCrelolo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTCimhihi_h,YWTCimhihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTCimlohi_h,YWTCimlohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTCimhilo_h,YWTCimhilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(YWTCimlolo_h,YWTCimlolo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the YWTC matrix :" << endl;
      for(int i=rowoff; i<nrows; i++) 
         for(int j=coloff; j<ncols; j++)
         {
            cout << "YWTC[" << i << "][" << j << "]re : "
                 << YWTCrehihi_h[j*nrows + i] << "  "
                 << YWTCrelohi_h[j*nrows + i] << endl
                 << "               "
                 << YWTCrehilo_h[j*nrows + i] << "  "
                 << YWTCrelolo_h[j*nrows + i] << endl;
            cout << "YWTC[" << i << "][" << j << "]im : "
                 << YWTCimhihi_h[j*nrows + i] << "  "
                 << YWTCimlohi_h[j*nrows + i] << endl
                 << "               "
                 << YWTCimhilo_h[j*nrows + i] << "  "
                 << YWTCimlolo_h[j*nrows + i] << endl;
         }
   }
}

void GPU_dbl4_small_Qupdate
 ( int dim, int szt, int idx,
   double *Qhihi_d, double *Qlohi_d, double *Qhilo_d, double *Qlolo_d,
   double *QWYThihi_d, double *QWYTlohi_d,
   double *QWYThilo_d, double *QWYTlolo_d,
   double *Qhihi_h, double *Qlohi_h, double *Qhilo_h, double *Qlolo_h,
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
   dbl4_small_Qupdate<<<nbrblocks,szt>>>
      (dim,rowdim,szt,coloff,Qhihi_d,Qlohi_d,Qhilo_d,Qlolo_d,
       QWYThihi_d,QWYTlohi_d,QWYThilo_d,QWYTlolo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_dbl_small_Qupdate(dim,rowdim,add);

   if(verbose)
   {
      const size_t szmat = dim*dim*sizeof(double);

      cudaMemcpy(Qhihi_h,Qhihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qlohi_h,Qlohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qhilo_h,Qhilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qlolo_h,Qlolo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the Q matrix :" << endl;
      int ix = 0;
      for(int i=0; i<dim; i++) 
         for(int j=0; j<dim; j++) 
         {
            cout << "Q[" << i << "][" << j << "] : "
                 << Qhihi_h[ix] << "  " << Qlohi_h[ix] << endl
                 << "          "
                 << Qhilo_h[ix] << "  " << Qlolo_h[ix] << endl;
            ix = ix + 1;
         }
   }
}

void GPU_cmplx4_small_Qupdate
 ( int dim, int szt, int idx,
   double *Qrehihi_d, double *Qrelohi_d, double *Qrehilo_d, double *Qrelolo_d,
   double *Qimhihi_d, double *Qimlohi_d, double *Qimhilo_d, double *Qimlolo_d,
   double *QWYTrehihi_d, double *QWYTrelohi_d,
   double *QWYTrehilo_d, double *QWYTrelolo_d,
   double *QWYTimhihi_d, double *QWYTimlohi_d,
   double *QWYTimhilo_d, double *QWYTimlolo_d,
   double *Qrehihi_h, double *Qrelohi_h, double *Qrehilo_h, double *Qrelolo_h,
   double *Qimhihi_h, double *Qimlohi_h, double *Qimhilo_h, double *Qimlolo_h,
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
   cmplx4_small_Qupdate<<<nbrblocks,szt>>>
      (dim,rowdim,szt,coloff,
          Qrehihi_d,   Qrelohi_d,   Qrehilo_d,   Qrelolo_d,
          Qimhihi_d,   Qimlohi_d,   Qimhilo_d,   Qimlolo_d,
       QWYTrehihi_d,QWYTrelohi_d,QWYTrehilo_d,QWYTrelolo_d,
       QWYTimhihi_d,QWYTimlohi_d,QWYTimhilo_d,QWYTimlolo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_cmplx_small_Qupdate(dim,rowdim,add);

   if(verbose)
   {
      const size_t szmat = dim*dim*sizeof(double);

      cudaMemcpy(Qrehihi_h,Qrehihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qrelohi_h,Qrelohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qrehilo_h,Qrehilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qrelolo_h,Qrelolo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qimhihi_h,Qimhihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qimlohi_h,Qimlohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qimhilo_h,Qimhilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Qimlolo_h,Qimlolo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the Q matrix :" << endl;
      int ix = 0;
      for(int i=0; i<dim; i++) 
         for(int j=0; j<dim; j++) 
         {
            cout << "Q[" << i << "][" << j << "]re : "
                 << Qrehihi_h[ix] << "  " << Qrelohi_h[ix] << endl
                 << "            "
                 << Qrehilo_h[ix] << "  " << Qrelolo_h[ix] << endl;
            cout << "Q[" << i << "][" << j << "]im : "
                 << Qimhihi_h[ix] << "  " << Qimlohi_h[ix] << endl
                 << "            "
                 << Qimhilo_h[ix] << "  " << Qimlolo_h[ix] << endl;
            ix = ix + 1;
         }
   }
}

void GPU_dbl4_small_R_add_YWTC
 ( int nrows, int ncols, int szt, int idx,
   double *Rhihi_d, double *Rlohi_d, double *Rhilo_d, double *Rlolo_d,
   double *YWTChihi_d, double *YWTClohi_d,
   double *YWTChilo_d, double *YWTClolo_d,
   double *Rhihi_h, double *Rlohi_h, double *Rhilo_h, double *Rlolo_h,
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
   dbl4_small_R_add_YWTC<<<nbrblocks,szt>>>
      (nrows,coldim,szt,rowoff,coloff,Rhihi_d,Rlohi_d,Rhilo_d,Rlolo_d,
       YWTChihi_d,YWTClohi_d,YWTChilo_d,YWTClolo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_dbl_small_R_add_YWTC(nrows,coldim,szt,rowoff,coloff,add);

   if(verbose)
   {
      const size_t szmat = nrows*ncols*sizeof(double);

      cudaMemcpy(Rhihi_h,Rhihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Rlohi_h,Rlohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Rhilo_h,Rhilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Rlolo_h,Rlolo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the R matrix :" << endl;
      for(int i=rowoff; i<nrows; i++) 
         for(int j=coloff; j<ncols; j++)
            cout << "R[" << i << "][" << j << "] : "
                 << Rhihi_h[j*nrows + i] << "  "
                 << Rlohi_h[j*nrows + i] << endl
                 << "          "
                 << Rhilo_h[j*nrows + i] << "  "
                 << Rlolo_h[j*nrows + i] << endl;
   }
}

void GPU_cmplx4_small_R_add_YWHC
 ( int nrows, int ncols, int szt, int idx,
   double *Rrehihi_d, double *Rrelohi_d, double *Rrehilo_d, double *Rrelolo_d,
   double *Rimhihi_d, double *Rimlohi_d, double *Rimhilo_d, double *Rimlolo_d,
   double *YWTCrehihi_d, double *YWTCrelohi_d,
   double *YWTCrehilo_d, double *YWTCrelolo_d,
   double *YWTCimhihi_d, double *YWTCimlohi_d,
   double *YWTCimhilo_d, double *YWTCimlolo_d,
   double *Rrehihi_h, double *Rrelohi_h, double *Rrehilo_h, double *Rrelolo_h,
   double *Rimhihi_h, double *Rimlohi_h, double *Rimhilo_h, double *Rimlolo_h,
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
   cmplx4_small_R_add_YWHC<<<nbrblocks,szt>>>
      (nrows,coldim,szt,rowoff,coloff,
          Rrehihi_d,   Rrelohi_d,   Rrehilo_d,   Rrelolo_d,
          Rimhihi_d,   Rimlohi_d,   Rimhilo_d,   Rimlolo_d,
       YWTCrehihi_d,YWTCrelohi_d,YWTCrehilo_d,YWTCrelolo_d,
       YWTCimhihi_d,YWTCimlohi_d,YWTCimhilo_d,YWTCimlolo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;
   flopcount_cmplx_small_R_add_YWHC(nrows,coldim,szt,rowoff,coloff,add);

   if(verbose)
   {
      const size_t szmat = nrows*ncols*sizeof(double);

      cudaMemcpy(Rrehihi_h,Rrehihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Rrelohi_h,Rrelohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Rrehilo_h,Rrehilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Rrelolo_h,Rrelolo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Rimhihi_h,Rimhihi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Rimlohi_h,Rimlohi_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Rimhilo_h,Rimhilo_d,szmat,cudaMemcpyDeviceToHost);
      cudaMemcpy(Rimlolo_h,Rimlolo_d,szmat,cudaMemcpyDeviceToHost);

      cout << "the R matrix :" << endl;
      for(int i=rowoff; i<nrows; i++) 
         for(int j=coloff; j<ncols; j++)
         {
            cout << "R[" << i << "][" << j << "]re : "
                 << Rrehihi_h[j*nrows + i] << "  "
                 << Rrelohi_h[j*nrows + i] << endl
                 << "            "
                 << Rrehilo_h[j*nrows + i] << "  "
                 << Rrelolo_h[j*nrows + i] << endl;
            cout << "R[" << i << "][" << j << "]im : "
                 << Rimhihi_h[j*nrows + i] << "  "
                 << Rimlohi_h[j*nrows + i] << endl
                 << "            "
                 << Rimhilo_h[j*nrows + i] << "  "
                 << Rimlolo_h[j*nrows + i] << endl;
         }
   }
}

void GPU_dbl4_blocked_houseqr
 ( int nrows, int ncols, int szt, int nbt,
   double **Ahihi, double **Alohi, double **Ahilo, double **Alolo,
   double **Qhihi, double **Qlohi, double **Qhilo, double **Qlolo,
   double **Rhihi, double **Rlohi, double **Rhilo, double **Rlolo,
   double *houselapms, double *RTvlapms, double *tileRlapms,
   double *vb2Wlapms, double *WYTlapms, double *QWYTlapms, double *Qaddlapms,
   double *YWTlapms, double *YWTClapms, double *Raddlapms,
   double *walltimesec, long long int *addcnt, long long int *mulcnt,
   long long int *divcnt, long long int *sqrtcnt, bool verbose )
{
   const int dim = nrows*ncols;          // total number of doubles
   const int nrows2 = nrows*nrows;
   double *Ahihi_h = new double[dim];    // A on the host
   double *Alohi_h = new double[dim]; 
   double *Ahilo_h = new double[dim];
   double *Alolo_h = new double[dim]; 
   double *Ahihi_d;                      // A on the device
   double *Alohi_d; 
   double *Ahilo_d; 
   double *Alolo_d; 
   double *Qhihi_h = new double[nrows2]; // Q on the host
   double *Qlohi_h = new double[nrows2]; 
   double *Qhilo_h = new double[nrows2]; 
   double *Qlolo_h = new double[nrows2]; 
   double *Qhihi_d;                      // Q on the device
   double *Qlohi_d;
   double *Qhilo_d;
   double *Qlolo_d;
   double *vhihi_h = new double[nrows];  // Householder vector
   double *vlohi_h = new double[nrows];
   double *vhilo_h = new double[nrows];
   double *vlolo_h = new double[nrows];
   double *betahihi_h = new double[szt]; //  beta on the host
   double *betalohi_h = new double[szt]; 
   double *betahilo_h = new double[szt]; 
   double *betalolo_h = new double[szt]; 
   double *betahihi_d;                   // beta on the device
   double *betalohi_d;
   double *betahilo_d;
   double *betalolo_d;
   double *Vhihi_h = new double[nrows*szt]; // V matrix
   double *Vlohi_h = new double[nrows*szt];
   double *Vhilo_h = new double[nrows*szt];
   double *Vlolo_h = new double[nrows*szt];
   double *Vhihi_d;                         // V on the device
   double *Vlohi_d;
   double *Vhilo_d;
   double *Vlolo_d;
   double *Whihi_h = new double[nrows*szt]; // W on the host
   double *Wlohi_h = new double[nrows*szt];
   double *Whilo_h = new double[nrows*szt];
   double *Wlolo_h = new double[nrows*szt];
   double *Whihi_d;                         // W on the device
   double *Wlohi_d;
   double *Whilo_d;
   double *Wlolo_d;
   double *WYThihi_h = new double[nrows2];  // W*Y^T 
   double *WYTlohi_h = new double[nrows2];
   double *WYThilo_h = new double[nrows2];
   double *WYTlolo_h = new double[nrows2];
   double *WYThihi_d;                       // WYT on the device
   double *WYTlohi_d;
   double *WYThilo_d;
   double *WYTlolo_d;
   double *YWThihi_h = new double[nrows2];  // Y*W^T
   double *YWTlohi_h = new double[nrows2];
   double *YWThilo_h = new double[nrows2];
   double *YWTlolo_h = new double[nrows2];
   double *YWThihi_d;                       // YWT on the device
   double *YWTlohi_d;
   double *YWThilo_d;
   double *YWTlolo_d;
   double *QWYThihi_h = new double[nrows2]; // Q*WY^T
   double *QWYTlohi_h = new double[nrows2];
   double *QWYThilo_h = new double[nrows2];
   double *QWYTlolo_h = new double[nrows2];
   double *QWYThihi_d;                      // QWYT on the device
   double *QWYTlohi_d;
   double *QWYThilo_d;
   double *QWYTlolo_d;
   double *YWTChihi_h = new double[dim];    // YWT*C on the host
   double *YWTClohi_h = new double[dim];
   double *YWTChilo_h = new double[dim];
   double *YWTClolo_h = new double[dim];
   double *YWTChihi_d;                      // YWTC on the device
   double *YWTClohi_d;
   double *YWTChilo_d;
   double *YWTClolo_d;
   double *RTdotvhihi_h = new double[nrows2]; // R^T dotted with v
   double *RTdotvlohi_h = new double[nrows2];
   double *RTdotvhilo_h = new double[nrows2];
   double *RTdotvlolo_h = new double[nrows2];
   double *RTdotvhihi_d;                      // RTdotv on the device
   double *RTdotvlohi_d;
   double *RTdotvhilo_d;
   double *RTdotvlolo_d;
   double *bRTvhihi_h = new double[nrows];  // beta*R^T*v
   double *bRTvlohi_h = new double[nrows];
   double *bRTvhilo_h = new double[nrows];
   double *bRTvlolo_h = new double[nrows];
   double *bRTvhihi_d;                      // bRTv on the device
   double *bRTvlohi_d;
   double *bRTvhilo_d;
   double *bRTvlolo_d;
   double *sumshihi_h = new double[nrows];  // subsums for large house
   double *sumslohi_h = new double[nrows];
   double *sumshilo_h = new double[nrows];
   double *sumslolo_h = new double[nrows];
   double *sumshihi_d;                      // sums on the device
   double *sumslohi_d;
   double *sumshilo_d;
   double *sumslolo_d;
   double sigmahihi_h,sigmalohi_h,sigmahilo_h,sigmalolo_h;
   double *sigmahihi_d;                     // sigma on the device
   double *sigmalohi_d;
   double *sigmahilo_d;
   double *sigmalolo_d;

   int ix = 0;                          // copy the columns of A to A_h
   for(int j=0; j<ncols; j++)   
      for(int i=0; i<nrows; i++)
      {
         Ahihi_h[ix]   = Ahihi[i][j];
         Alohi_h[ix]   = Alohi[i][j];
         Ahilo_h[ix]   = Ahilo[i][j];
         Alolo_h[ix++] = Alolo[i][j];
      }

   ix = 0;                              // initialize Q with identity
   for(int i=0; i<nrows; i++)
   {
      for(int j=0; j<nrows; j++)
      {
         if(i == j)
         {
            Qhihi_h[ix]   = 1.0;
            Qlohi_h[ix]   = 0.0;
            Qhilo_h[ix]   = 0.0;
            Qlolo_h[ix++] = 0.0;
         }
         else
         {
            Qhihi_h[ix]   = 0.0;
            Qlohi_h[ix]   = 0.0;
            Qhilo_h[ix]   = 0.0;
            Qlolo_h[ix++] = 0.0;
         }
      }
   }
   const size_t sznum = dim*sizeof(double);
   cudaMalloc((void**)&Ahihi_d,sznum);
   cudaMalloc((void**)&Alohi_d,sznum);
   cudaMalloc((void**)&Ahilo_d,sznum);
   cudaMalloc((void**)&Alolo_d,sznum);
   cudaMemcpy(Ahihi_d,Ahihi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Alohi_d,Alohi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Ahilo_d,Ahilo_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Alolo_d,Alolo_h,sznum,cudaMemcpyHostToDevice);

   const size_t szbeta = szt*sizeof(double);
   cudaMalloc((void**)&betahihi_d,szbeta);
   cudaMalloc((void**)&betalohi_d,szbeta);
   cudaMalloc((void**)&betahilo_d,szbeta);
   cudaMalloc((void**)&betalolo_d,szbeta);

   for(int i=0; i<szt; i++)
   {
      betahihi_h[i] = 0.0;
      betalohi_h[i] = 0.0;
      betahilo_h[i] = 0.0;
      betalolo_h[i] = 0.0;
   }
   cudaMemcpy(betahihi_d,betahihi_h,szbeta,cudaMemcpyHostToDevice);
   cudaMemcpy(betalohi_d,betalohi_h,szbeta,cudaMemcpyHostToDevice);
   cudaMemcpy(betahilo_d,betahilo_h,szbeta,cudaMemcpyHostToDevice);
   cudaMemcpy(betalolo_d,betalolo_h,szbeta,cudaMemcpyHostToDevice);

   const size_t szhouse = nrows*sizeof(double);
   const size_t szpad = szt*sizeof(double);  // padding for nonsquare tiles
   const size_t szVandW = szt*szhouse;
   cudaMalloc((void**)&Vhihi_d,szVandW + szpad); // padding only in allocation
   cudaMalloc((void**)&Vlohi_d,szVandW + szpad);
   cudaMalloc((void**)&Vhilo_d,szVandW + szpad);
   cudaMalloc((void**)&Vlolo_d,szVandW + szpad);

   ix = 0;
   for(int i=0; i<nrows*szt; i++)
   {
      Vhihi_h[ix] = 0.0; 
      Vlohi_h[ix] = 0.0; 
      Vhilo_h[ix] = 0.0; 
      Vlolo_h[ix++] = 0.0; 
   }
   Vhihi_h[--ix] = 1.0; // initialize last vector for square tiles

   cudaMemcpy(Vhihi_d,Vhihi_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vlohi_d,Vlohi_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vhilo_d,Vhilo_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vlolo_d,Vlolo_h,szVandW,cudaMemcpyHostToDevice);
   cudaMalloc((void**)&Whihi_d,szVandW + szpad); // padding only in allocation
   cudaMalloc((void**)&Wlohi_d,szVandW + szpad); 
   cudaMalloc((void**)&Whilo_d,szVandW + szpad); 
   cudaMalloc((void**)&Wlolo_d,szVandW + szpad); 

   cudaMalloc((void**)&RTdotvhihi_d,szVandW + szpad);
   cudaMalloc((void**)&RTdotvlohi_d,szVandW + szpad);
   cudaMalloc((void**)&RTdotvhilo_d,szVandW + szpad);
   cudaMalloc((void**)&RTdotvlolo_d,szVandW + szpad);
   cudaMalloc((void**)&bRTvhihi_d,szhouse + szpad);
   cudaMalloc((void**)&bRTvlohi_d,szhouse + szpad);
   cudaMalloc((void**)&bRTvhilo_d,szhouse + szpad);
   cudaMalloc((void**)&bRTvlolo_d,szhouse + szpad);

   cudaMalloc((void**)&sumshihi_d,szhouse);
   cudaMalloc((void**)&sumslohi_d,szhouse);
   cudaMalloc((void**)&sumshilo_d,szhouse);
   cudaMalloc((void**)&sumslolo_d,szhouse);
   cudaMalloc((void**)&sigmahihi_d,sizeof(double));
   cudaMalloc((void**)&sigmalohi_d,sizeof(double));
   cudaMalloc((void**)&sigmahilo_d,sizeof(double));
   cudaMalloc((void**)&sigmalolo_d,sizeof(double));

   const size_t szWYT = nrows2*sizeof(double);
   cudaMalloc((void**)&WYThihi_d,szWYT + szpad); // padding for W*Y^T product
   cudaMalloc((void**)&WYTlohi_d,szWYT + szpad); 
   cudaMalloc((void**)&WYThilo_d,szWYT + szpad); 
   cudaMalloc((void**)&WYTlolo_d,szWYT + szpad); 
   cudaMalloc((void**)&Qhihi_d,szWYT + szpad);
   cudaMalloc((void**)&Qlohi_d,szWYT + szpad);
   cudaMalloc((void**)&Qhilo_d,szWYT + szpad);
   cudaMalloc((void**)&Qlolo_d,szWYT + szpad);
   cudaMemcpy(Qhihi_d,Qhihi_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qlohi_d,Qlohi_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qhilo_d,Qhilo_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qlolo_d,Qlolo_h,szWYT,cudaMemcpyHostToDevice);
   cudaMalloc((void**)&QWYThihi_d,szWYT + szpad);
   cudaMalloc((void**)&QWYTlohi_d,szWYT + szpad);
   cudaMalloc((void**)&QWYThilo_d,szWYT + szpad);
   cudaMalloc((void**)&QWYTlolo_d,szWYT + szpad);

   const size_t szYWT = nrows2*sizeof(double);
   cudaMalloc((void**)&YWThihi_d,szYWT + szpad); // padding for Y*W^T product
   cudaMalloc((void**)&YWTlohi_d,szYWT + szpad);
   cudaMalloc((void**)&YWThilo_d,szYWT + szpad);
   cudaMalloc((void**)&YWTlolo_d,szYWT + szpad);
   cudaMalloc((void**)&YWTChihi_d,sznum + szpad);
   cudaMalloc((void**)&YWTClohi_d,sznum + szpad);
   cudaMalloc((void**)&YWTChilo_d,sznum + szpad);
   cudaMalloc((void**)&YWTClolo_d,sznum + szpad);

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
            GPU_dbl4_small_house
               (nrows,ncols,szt,nbt,colidx,nrows1,k,L,
                   Ahihi_h,   Alohi_h,   Ahilo_h,   Alolo_h,
                   Ahihi_d,   Alohi_d,   Ahilo_d,   Alolo_d,
                   vhihi_h,   vlohi_h,   vhilo_h,   vlolo_h,
                   Vhihi_d,   Vlohi_d,   Vhilo_d,   Vlolo_d,
                betahihi_h,betalohi_h,betahilo_h,betalolo_h,
                betahihi_d,betalohi_d,betahilo_d,betalolo_d,
                houselapms,addcnt,mulcnt,divcnt,sqrtcnt,verbose);

            if((betahihi_h[L] == 0.0) && (betalohi_h[L] == 0.0) &&
               (betahilo_h[L] == 0.0) && (betalolo_h[L] == 0.0))
            {
               if(verbose) cout << "Zero beta detected." << endl;
            }
            else
            {
               GPU_dbl4_small_leftRupdate
                  (nrows,ncols,szt,colidx,k,L,
                      Ahihi_h,   Alohi_h,   Ahilo_h,   Alolo_h,
                      Ahihi_d,   Alohi_d,   Ahilo_d,   Alolo_d,
                      Vhihi_d,   Vlohi_d,   Vhilo_d,   Vlolo_d,
                   betahihi_h,betalohi_h,betahilo_h,betalolo_h,
                   betahihi_d,betalohi_d,betahilo_d,betalolo_d,
                   tileRlapms,addcnt,mulcnt,verbose);
            }
         }
         else
         {
            GPU_dbl4_large_house
               (nrows,ncols,szt,nbt,colidx,nrows1,k,L,
                     Ahihi_h,     Alohi_h,     Ahilo_h,     Alolo_h,
                     Ahihi_d,     Alohi_d,     Ahilo_d,     Alolo_d,
                     vhihi_h,     vlohi_h,     vhilo_h,     vlolo_h,
                     Vhihi_d,     Vlohi_d,     Vhilo_d,     Vlolo_d,
                  betahihi_h,  betalohi_h,  betahilo_h,  betalolo_h,
                  betahihi_d,  betalohi_d,  betahilo_d,  betalolo_d,
                  sumshihi_h,  sumslohi_h,  sumshilo_h,  sumslolo_h,
                  sumshihi_d,  sumslohi_d,  sumshilo_d,  sumslolo_d,
                &sigmahihi_h,&sigmalohi_h,&sigmahilo_h,&sigmalolo_h,
                 sigmahihi_d, sigmalohi_d, sigmahilo_d, sigmalolo_d,
                houselapms,addcnt,mulcnt,divcnt,sqrtcnt,verbose);

            if((betahihi_h[L] == 0.0) && (betalohi_h[L] == 0.0) &&
               (betahilo_h[L] == 0.0) && (betalolo_h[L] == 0.0))
            {
               if(verbose) cout << "Zero beta detected." << endl;
            }
            else
            {
               GPU_dbl4_medium_leftRupdate
                  (nrows,ncols,szt,colidx,k,L,
                        Ahihi_h,     Alohi_h,     Ahilo_h,     Alolo_h,
                        Ahihi_d,     Alohi_d,     Ahilo_d,     Alolo_d,
                        Vhihi_d,     Vlohi_d,     Vhilo_d,     Vlolo_d,
                     betahihi_h,  betalohi_h,  betahilo_h,  betalolo_h,
                     betahihi_d,  betalohi_d,  betahilo_d,  betalolo_d,
                   RTdotvhihi_h,RTdotvlohi_h,RTdotvhilo_h,RTdotvlolo_h,
                   RTdotvhihi_d,RTdotvlohi_d,RTdotvhilo_d,RTdotvlolo_d,
                     bRTvhihi_h,  bRTvlohi_h,  bRTvhilo_h,  bRTvlolo_h,
                     bRTvhihi_d,  bRTvlohi_d,  bRTvhilo_d,  bRTvlolo_d,
                   RTvlapms,tileRlapms,addcnt,mulcnt,verbose);
            }
         }
      }
      GPU_dbl4_medium_VB_to_W
         (nrows,szt,szt,k,
             Vhihi_h,   Vlohi_h,   Vhilo_h,   Vlolo_h,
             Vhihi_d,   Vlohi_d,   Vhilo_d,   Vlolo_d,
             Whihi_h,   Wlohi_h,   Whilo_h,   Wlolo_h,
             Whihi_d,   Wlohi_d,   Whilo_d,   Wlolo_d,
           WYThihi_h, WYTlohi_h, WYThilo_h, WYTlolo_h,
           WYThihi_d, WYTlohi_d, WYThilo_d, WYTlolo_d,
          betahihi_h,betalohi_h,betahilo_h,betalolo_h,
          betahihi_d,betalohi_d,betahilo_d,betalolo_d,
          vb2Wlapms,addcnt,mulcnt,verbose);
/*
      GPU_dbl2_small_WYT
         (nrows-k*szt,szt,Whi_d,Wlo_d,Vhi_d,Vlo_d,WYThi_d,WYTlo_d,
          WYThi_h,WYTlo_h,WYTlapms,verbose);
 */
      GPU_dbl4_small_QWYT
         (nrows,szt,k,
             Qhihi_d,   Qlohi_d,   Qhilo_d,   Qlolo_d,
           WYThihi_d, WYTlohi_d, WYThilo_d, WYTlolo_d,
          QWYThihi_d,QWYTlohi_d,QWYThilo_d,QWYTlolo_d,
          QWYThihi_h,QWYTlohi_h,QWYThilo_h,QWYTlolo_h,
             Qhihi_h,   Qlohi_h,   Qhilo_h,   Qlolo_h,
          QWYTlapms,addcnt,mulcnt,verbose);

      GPU_dbl4_small_Qupdate
         (nrows,szt,k,
             Qhihi_d,   Qlohi_d,   Qhilo_d,   Qlolo_d,
          QWYThihi_d,QWYTlohi_d,QWYThilo_d,QWYTlolo_d,
             Qhihi_h,   Qlohi_h,   Qhilo_h,   Qlolo_h,
          Qaddlapms,addcnt,verbose);

      if(k < nbt-1)                                           // update R
      {
         GPU_dbl4_small_YWT
            (nrows,szt,k,
               Vhihi_d,  Vlohi_d,  Vhilo_d,  Vlolo_d,
               Whihi_d,  Wlohi_d,  Whilo_d,  Wlolo_d,
             YWThihi_d,YWTlohi_d,YWThilo_d,YWTlolo_d,
             YWThihi_h,YWTlohi_h,YWThilo_h,YWTlolo_h,
             YWTlapms,addcnt,mulcnt,verbose);

         GPU_dbl4_small_YWTC
            (nrows,ncols,szt,k,
              YWThihi_d, YWTlohi_d, YWThilo_d, YWTlolo_d,
                Ahihi_d,   Alohi_d,   Ahilo_d,   Alolo_d,
             YWTChihi_d,YWTClohi_d,YWTChilo_d,YWTClolo_d,
             YWTChihi_h,YWTClohi_h,YWTChilo_h,YWTClolo_h,
             YWTClapms,addcnt,mulcnt,verbose);

         GPU_dbl4_small_R_add_YWTC
            (nrows,ncols,szt,k,
                Ahihi_d,   Alohi_d,   Ahilo_d,   Alolo_d,
             YWTChihi_d,YWTClohi_d,YWTChilo_d,YWTClolo_d,
                Ahihi_h,   Alohi_h,   Ahilo_h,   Alolo_h,
             Raddlapms,addcnt,verbose);
      }
   }
   gettimeofday(&endtime,0);
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   *walltimesec = seconds + microseconds*1.0e-6;

   cudaMemcpy(Qhihi_h,Qhihi_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qlohi_h,Qlohi_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qhilo_h,Qhilo_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qlolo_h,Qlolo_d,szWYT,cudaMemcpyDeviceToHost);
   ix = 0;                                           // copy rows of Q
   for(int i=0; i<nrows; i++)
      for(int j=0; j<nrows; j++)
      {
         Qhihi[i][j] = Qhihi_h[ix];
         Qlohi[i][j] = Qlohi_h[ix];
         Qhilo[i][j] = Qhilo_h[ix];
         Qlolo[i][j] = Qlolo_h[ix++];
      }

   cudaMemcpy(Ahihi_h,Ahihi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Alohi_h,Alohi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Ahilo_h,Ahilo_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Alolo_h,Alolo_d,sznum,cudaMemcpyDeviceToHost);

   for(int i=0; i<nrows; i++)                       // copy columns of R
      for(int j=0; j<ncols; j++)
      {
         Rhihi[i][j] = Ahihi_h[j*nrows+i];
         Rlohi[i][j] = Alohi_h[j*nrows+i];
         Rhilo[i][j] = Ahilo_h[j*nrows+i];
         Rlolo[i][j] = Alolo_h[j*nrows+i];
      }

   free(Ahihi_h); free(Alohi_h); free(Ahilo_h); free(Alolo_h);
   free(Qhihi_h); free(Qlohi_h); free(Qhilo_h); free(Qlolo_h); 
   free(vhihi_h); free(vlohi_h); free(vhilo_h); free(vlolo_h);
   free(Vhihi_h); free(Vlohi_h); free(Vhilo_h); free(Vlolo_h);
   free(Whihi_h); free(Wlohi_h); free(Whilo_h); free(Wlolo_h);
   free(sumshihi_h); free(sumslohi_h); free(sumshilo_h); free(sumslolo_h);

   free(RTdotvhihi_h); free(RTdotvlohi_h);
   free(RTdotvhilo_h); free(RTdotvlolo_h);
   free(bRTvhihi_h); free(bRTvlohi_h); free(bRTvhilo_h); free(bRTvlolo_h);
   free(WYThihi_h); free(QWYThihi_h); free(WYThilo_h); free(QWYThilo_h);
   free(YWThihi_h); free(YWTChihi_h); free(YWThilo_h); free(YWTChilo_h);
   free(WYTlohi_h); free(QWYTlohi_h); free(WYTlolo_h); free(QWYTlolo_h);
   free(YWTlohi_h); free(YWTClohi_h); free(YWTlolo_h); free(YWTClolo_h);
}


void GPU_cmplx4_blocked_houseqr
 ( int nrows, int ncols, int szt, int nbt,
   double **Arehihi, double **Arelohi, double **Arehilo, double **Arelolo,
   double **Aimhihi, double **Aimlohi, double **Aimhilo, double **Aimlolo,
   double **Qrehihi, double **Qrelohi, double **Qrehilo, double **Qrelolo,
   double **Qimhihi, double **Qimlohi, double **Qimhilo, double **Qimlolo,
   double **Rrehihi, double **Rrelohi, double **Rrehilo, double **Rrelolo,
   double **Rimhihi, double **Rimlohi, double **Rimhilo, double **Rimlolo,
   double *houselapms, double *RHvlapms, double *tileRlapms,
   double *vb2Wlapms, double *WYHlapms, double *QWYHlapms, double *Qaddlapms,
   double *YWHlapms, double *YWHClapms, double *Raddlapms,
   double *walltimesec, long long int *addcnt, long long int *mulcnt,
   long long int *divcnt, long long int *sqrtcnt, bool verbose )
{
   const int dim = nrows*ncols;        // total number of doubles
   const int nrows2 = nrows*nrows;
   double *Arehihi_h = new double[dim];  // the real parts of A
   double *Arelohi_h = new double[dim];
   double *Arehilo_h = new double[dim];
   double *Arelolo_h = new double[dim];
   double *Aimhihi_h = new double[dim];  // the imaginary parts of A
   double *Aimlohi_h = new double[dim]; 
   double *Aimhilo_h = new double[dim]; 
   double *Aimlolo_h = new double[dim]; 
   double *Arehihi_d;                    // Are on the device
   double *Arelohi_d;
   double *Arehilo_d;
   double *Arelolo_d;
   double *Aimhihi_d;                    // Aim on the device
   double *Aimlohi_d;
   double *Aimhilo_d;
   double *Aimlolo_d;
   double *Qrehihi_h = new double[nrows2]; // real parts of Q 
   double *Qrelohi_h = new double[nrows2];
   double *Qrehilo_h = new double[nrows2];
   double *Qrelolo_h = new double[nrows2];
   double *Qimhihi_h = new double[nrows2]; // imaginary parts of Q
   double *Qimlohi_h = new double[nrows2];
   double *Qimhilo_h = new double[nrows2];
   double *Qimlolo_h = new double[nrows2];
   double *Qrehihi_d;                      // Qre on the device
   double *Qrelohi_d;
   double *Qrehilo_d;
   double *Qrelolo_d;
   double *Qimhihi_d;                      // Qim on the device
   double *Qimlohi_d;
   double *Qimhilo_d;
   double *Qimlolo_d;
   double *vrehihi_h = new double[nrows];  // real parts of Householder v
   double *vrelohi_h = new double[nrows];
   double *vrehilo_h = new double[nrows];
   double *vrelolo_h = new double[nrows];
   double *vimhihi_h = new double[nrows];  // imaginary parts of Householder v
   double *vimlohi_h = new double[nrows]; 
   double *vimhilo_h = new double[nrows]; 
   double *vimlolo_h = new double[nrows]; 
   double *betahihi_h = new double[szt+1]; // beta
   double *betalohi_h = new double[szt+1]; 
   double *betahilo_h = new double[szt+1]; 
   double *betalolo_h = new double[szt+1]; 
   double *betahihi_d;                     // beta on the device
   double *betalohi_d;
   double *betahilo_d;
   double *betalolo_d;
   double *Vrehihi_h = new double[nrows*szt]; // real parts of V
   double *Vrelohi_h = new double[nrows*szt]; 
   double *Vrehilo_h = new double[nrows*szt]; 
   double *Vrelolo_h = new double[nrows*szt]; 
   double *Vimhihi_h = new double[nrows*szt]; // imaginary parts of V
   double *Vimlohi_h = new double[nrows*szt];
   double *Vimhilo_h = new double[nrows*szt];
   double *Vimlolo_h = new double[nrows*szt];
   double *Vrehihi_d;                         // Vre on device
   double *Vrelohi_d;
   double *Vrehilo_d;
   double *Vrelolo_d;
   double *Vimhihi_d;                         // Vim on device
   double *Vimlohi_d;
   double *Vimhilo_d;
   double *Vimlolo_d;
   double *Wrehihi_h = new double[nrows*szt]; // real parts of W
   double *Wrelohi_h = new double[nrows*szt];
   double *Wrehilo_h = new double[nrows*szt];
   double *Wrelolo_h = new double[nrows*szt];
   double *Wimhihi_h = new double[nrows*szt]; // imaginary parts of W
   double *Wimlohi_h = new double[nrows*szt];
   double *Wimhilo_h = new double[nrows*szt];
   double *Wimlolo_h = new double[nrows*szt];
   double *Wrehihi_d;                         // Wre on the device
   double *Wrelohi_d;
   double *Wrehilo_d;
   double *Wrelolo_d;
   double *Wimhihi_d;                         // Wim on the device
   double *Wimlohi_d;
   double *Wimhilo_d;
   double *Wimlolo_d;
   double *WYTrehihi_h = new double[nrows2];  // real parts of W*Y^T
   double *WYTrelohi_h = new double[nrows2];
   double *WYTrehilo_h = new double[nrows2];
   double *WYTrelolo_h = new double[nrows2];
   double *WYTimhihi_h = new double[nrows2];  // imaginary parts of W*Y^T
   double *WYTimlohi_h = new double[nrows2];
   double *WYTimhilo_h = new double[nrows2];
   double *WYTimlolo_h = new double[nrows2];
   double *WYTrehihi_d;                       // WYTre on the device 
   double *WYTrelohi_d;
   double *WYTrehilo_d;
   double *WYTrelolo_d;
   double *WYTimhihi_d;                       // WYTim on the device
   double *WYTimlohi_d;
   double *WYTimhilo_d;
   double *WYTimlolo_d;
   double *YWTrehihi_h = new double[nrows2];  // real parts of Y*W^T
   double *YWTrelohi_h = new double[nrows2];
   double *YWTrehilo_h = new double[nrows2];
   double *YWTrelolo_h = new double[nrows2];
   double *YWTimhihi_h = new double[nrows2];  // imaginary parts of Y*W^T
   double *YWTimlohi_h = new double[nrows2];
   double *YWTimhilo_h = new double[nrows2];
   double *YWTimlolo_h = new double[nrows2];
   double *YWTrehihi_d;                       // YWTre on the device
   double *YWTrelohi_d;
   double *YWTrehilo_d;
   double *YWTrelolo_d;
   double *YWTimhihi_d;                       // YWTim on the device
   double *YWTimlohi_d; 
   double *YWTimhilo_d; 
   double *YWTimlolo_d; 
   double *QWYTrehihi_h = new double[nrows2]; // real parts of Q*WY^T
   double *QWYTrelohi_h = new double[nrows2];
   double *QWYTrehilo_h = new double[nrows2];
   double *QWYTrelolo_h = new double[nrows2];
   double *QWYTimhihi_h = new double[nrows2]; // imaginary parts of Q*WY^T
   double *QWYTimlohi_h = new double[nrows2];
   double *QWYTimhilo_h = new double[nrows2];
   double *QWYTimlolo_h = new double[nrows2];
   double *QWYTrehihi_d;                      // QWYTre on the device
   double *QWYTrelohi_d;
   double *QWYTrehilo_d;
   double *QWYTrelolo_d;
   double *QWYTimhihi_d;                      // QWYTim on the device
   double *QWYTimlohi_d;
   double *QWYTimhilo_d;
   double *QWYTimlolo_d;
   double *YWTCrehihi_h = new double[dim];    // real parts of YWT*C
   double *YWTCrelohi_h = new double[dim]; 
   double *YWTCrehilo_h = new double[dim]; 
   double *YWTCrelolo_h = new double[dim]; 
   double *YWTCimhihi_h = new double[dim];    // imaginary parts of YWT*C
   double *YWTCimlohi_h = new double[dim];
   double *YWTCimhilo_h = new double[dim];
   double *YWTCimlolo_h = new double[dim];
   double *YWTCrehihi_d;                      // YWTCre on the device
   double *YWTCrelohi_d;
   double *YWTCrehilo_d;
   double *YWTCrelolo_d;
   double *YWTCimhihi_d;                      // YWTCim on the device
   double *YWTCimlohi_d;
   double *YWTCimhilo_d;
   double *YWTCimlolo_d;
   double *RHdotvrehihi_h = new double[nrows2]; // real R^H dotted with v
   double *RHdotvrelohi_h = new double[nrows2]; 
   double *RHdotvrehilo_h = new double[nrows2]; 
   double *RHdotvrelolo_h = new double[nrows2]; 
   double *RHdotvimhihi_h = new double[nrows2]; // imaginary R^H dotted with v
   double *RHdotvimlohi_h = new double[nrows2]; 
   double *RHdotvimhilo_h = new double[nrows2]; 
   double *RHdotvimlolo_h = new double[nrows2]; 
   double *RHdotvrehihi_d;                      // RHdotvre on the device
   double *RHdotvrelohi_d;
   double *RHdotvrehilo_d;
   double *RHdotvrelolo_d;
   double *RHdotvimhihi_d;                      // RHdotvim on the device
   double *RHdotvimlohi_d;
   double *RHdotvimhilo_d;
   double *RHdotvimlolo_d;
   double *bRHvrehihi_h = new double[nrows];  // real parts of beta*R^H*v
   double *bRHvrelohi_h = new double[nrows];
   double *bRHvrehilo_h = new double[nrows];
   double *bRHvrelolo_h = new double[nrows];
   double *bRHvimhihi_h = new double[nrows];  // imaginary parts of beta*R^H*v
   double *bRHvimlohi_h = new double[nrows];
   double *bRHvimhilo_h = new double[nrows];
   double *bRHvimlolo_h = new double[nrows];
   double *bRHvrehihi_d;                      // bRHvre on the device
   double *bRHvrelohi_d;
   double *bRHvrehilo_d;
   double *bRHvrelolo_d;
   double *bRHvimhihi_d;                      // bRHvim on the device
   double *bRHvimlohi_d; 
   double *bRHvimhilo_d; 
   double *bRHvimlolo_d; 
   double *sumshihi_h = new double[nrows];  // subsums for large house
   double *sumslohi_h = new double[nrows]; 
   double *sumshilo_h = new double[nrows]; 
   double *sumslolo_h = new double[nrows]; 
   double *sumshihi_d;                      // sums on the device
   double *sumslohi_d;
   double *sumshilo_d;
   double *sumslolo_d;
   double sigmahihi_h,sigmalohi_h,sigmahilo_h,sigmalolo_h;
   double *sigmahihi_d;                     // sigma on the device
   double *sigmalohi_d;
   double *sigmahilo_d;
   double *sigmalolo_d;

   int ix = 0;                            // copy the columns of A to A_h
   for(int j=0; j<ncols; j++)   
      for(int i=0; i<nrows; i++)
      {
         Arehihi_h[ix]   = Arehihi[i][j];
         Arelohi_h[ix]   = Arelohi[i][j];
         Arehilo_h[ix]   = Arehilo[i][j];
         Arelolo_h[ix]   = Arelolo[i][j];
         Aimhihi_h[ix]   = Aimhihi[i][j];
         Aimlohi_h[ix]   = Aimlohi[i][j];
         Aimhilo_h[ix]   = Aimhilo[i][j];
         Aimlolo_h[ix++] = Aimlolo[i][j];
      }

   ix = 0;                                // initialize Q with identity
   for(int i=0; i<nrows; i++)
   {
      for(int j=0; j<nrows; j++)
      {
         if(i == j)
         {
            Qrehihi_h[ix]   = 1.0;
            Qrelohi_h[ix]   = 0.0;
            Qrehilo_h[ix]   = 0.0;
            Qrelolo_h[ix]   = 0.0;
            Qimhihi_h[ix]   = 0.0;
            Qimlohi_h[ix]   = 0.0;
            Qimhilo_h[ix]   = 0.0;
            Qimlolo_h[ix++] = 0.0;
         }
         else
         {
            Qrehihi_h[ix]   = 0.0;
            Qrelohi_h[ix]   = 0.0;
            Qrehilo_h[ix]   = 0.0;
            Qrelolo_h[ix]   = 0.0;
            Qimhihi_h[ix]   = 0.0;
            Qimlohi_h[ix]   = 0.0;
            Qimhilo_h[ix]   = 0.0;
            Qimlolo_h[ix++] = 0.0;
         }
         // cout << "Q[" << ix-1 << "] : "
         //      << Qre_h[ix-1] << "  " << Qim_h[ix-1] << endl;
      }
   }
   const size_t sznum = dim*sizeof(double);
   cudaMalloc((void**)&Arehihi_d,sznum);
   cudaMalloc((void**)&Arelohi_d,sznum);
   cudaMalloc((void**)&Arehilo_d,sznum);
   cudaMalloc((void**)&Arelolo_d,sznum);
   cudaMalloc((void**)&Aimhihi_d,sznum);
   cudaMalloc((void**)&Aimlohi_d,sznum);
   cudaMalloc((void**)&Aimhilo_d,sznum);
   cudaMalloc((void**)&Aimlolo_d,sznum);
   cudaMemcpy(Arehihi_d,Arehihi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Arelohi_d,Arelohi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Arehilo_d,Arehilo_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Arelolo_d,Arelolo_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Aimhihi_d,Aimhihi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Aimlohi_d,Aimlohi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Aimhilo_d,Aimhilo_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Aimlolo_d,Aimlolo_h,sznum,cudaMemcpyHostToDevice);
   // allocate one extra beta for use in the cmplx2_normalize
   const size_t szbeta = (szt+1)*sizeof(double);
   cudaMalloc((void**)&betahihi_d,szbeta);
   cudaMalloc((void**)&betalohi_d,szbeta);
   cudaMalloc((void**)&betahilo_d,szbeta);
   cudaMalloc((void**)&betalolo_d,szbeta);

   for(int i=0; i<=szt; i++)
   {
      betahihi_h[i] = 0.0;
      betalohi_h[i] = 0.0;
      betahilo_h[i] = 0.0;
      betalolo_h[i] = 0.0;
   }
   cudaMemcpy(betahihi_d,betahihi_h,szbeta,cudaMemcpyHostToDevice);
   cudaMemcpy(betalohi_d,betalohi_h,szbeta,cudaMemcpyHostToDevice);
   cudaMemcpy(betahilo_d,betahilo_h,szbeta,cudaMemcpyHostToDevice);
   cudaMemcpy(betalolo_d,betalolo_h,szbeta,cudaMemcpyHostToDevice);

   const size_t szhouse = nrows*sizeof(double);
   const size_t szpad = szt*sizeof(double);    // padding for nonsquare tiles
   // const size_t szpad = 0;
   const size_t szVandW = szt*szhouse;
   cudaMalloc((void**)&Vrehihi_d,szVandW + szpad); // padding added
   cudaMalloc((void**)&Vrelohi_d,szVandW + szpad);
   cudaMalloc((void**)&Vrehilo_d,szVandW + szpad);
   cudaMalloc((void**)&Vrelolo_d,szVandW + szpad);
   cudaMalloc((void**)&Vimhihi_d,szVandW + szpad);
   cudaMalloc((void**)&Vimlohi_d,szVandW + szpad);
   cudaMalloc((void**)&Vimhilo_d,szVandW + szpad);
   cudaMalloc((void**)&Vimlolo_d,szVandW + szpad);

   ix = 0;
   for(int i=0; i<nrows*szt; i++)
   {
      Vrehihi_h[ix]   = 0.0; 
      Vrelohi_h[ix]   = 0.0; 
      Vrehilo_h[ix]   = 0.0; 
      Vrelolo_h[ix]   = 0.0; 
      Vimhihi_h[ix]   = 0.0; 
      Vimlohi_h[ix]   = 0.0; 
      Vimhilo_h[ix]   = 0.0; 
      Vimlolo_h[ix++] = 0.0; 
   }
   Vrehihi_h[--ix] = 1.0; // initialize last vector for square tiles

   cudaMemcpy(Vrehihi_d,Vrehihi_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vrelohi_d,Vrelohi_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vrehilo_d,Vrehilo_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vrelolo_d,Vrelolo_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vimhihi_d,Vimhihi_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vimlohi_d,Vimlohi_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vimhilo_d,Vimhilo_h,szVandW,cudaMemcpyHostToDevice);
   cudaMemcpy(Vimlolo_d,Vimlolo_h,szVandW,cudaMemcpyHostToDevice);
   cudaMalloc((void**)&Wrehihi_d,szVandW + szpad); // padding added
   cudaMalloc((void**)&Wrelohi_d,szVandW + szpad);
   cudaMalloc((void**)&Wrehilo_d,szVandW + szpad);
   cudaMalloc((void**)&Wrelolo_d,szVandW + szpad);
   cudaMalloc((void**)&Wimhihi_d,szVandW + szpad);
   cudaMalloc((void**)&Wimlohi_d,szVandW + szpad);
   cudaMalloc((void**)&Wimhilo_d,szVandW + szpad);
   cudaMalloc((void**)&Wimlolo_d,szVandW + szpad);

   cudaMalloc((void**)&RHdotvrehihi_d,szVandW + szpad);
   cudaMalloc((void**)&RHdotvrelohi_d,szVandW + szpad);
   cudaMalloc((void**)&RHdotvrehilo_d,szVandW + szpad);
   cudaMalloc((void**)&RHdotvrelolo_d,szVandW + szpad);
   cudaMalloc((void**)&RHdotvimhihi_d,szVandW + szpad);
   cudaMalloc((void**)&RHdotvimlohi_d,szVandW + szpad);
   cudaMalloc((void**)&RHdotvimhilo_d,szVandW + szpad);
   cudaMalloc((void**)&RHdotvimlolo_d,szVandW + szpad);
   cudaMalloc((void**)&bRHvrehihi_d,szhouse + szpad);
   cudaMalloc((void**)&bRHvrelohi_d,szhouse + szpad);
   cudaMalloc((void**)&bRHvrehilo_d,szhouse + szpad);
   cudaMalloc((void**)&bRHvrelolo_d,szhouse + szpad);
   cudaMalloc((void**)&bRHvimhihi_d,szhouse + szpad);
   cudaMalloc((void**)&bRHvimlohi_d,szhouse + szpad);
   cudaMalloc((void**)&bRHvimhilo_d,szhouse + szpad);
   cudaMalloc((void**)&bRHvimlolo_d,szhouse + szpad);

   cudaMalloc((void**)&sumshihi_d,szhouse);
   cudaMalloc((void**)&sumslohi_d,szhouse);
   cudaMalloc((void**)&sumshilo_d,szhouse);
   cudaMalloc((void**)&sumslolo_d,szhouse);
   cudaMalloc((void**)&sigmahihi_d,sizeof(double));
   cudaMalloc((void**)&sigmalohi_d,sizeof(double));
   cudaMalloc((void**)&sigmahilo_d,sizeof(double));
   cudaMalloc((void**)&sigmalolo_d,sizeof(double));

   const size_t szWYT = nrows2*sizeof(double);
   cudaMalloc((void**)&WYTrehihi_d,szWYT + szpad); // padding for W*Y^T 
   cudaMalloc((void**)&WYTrelohi_d,szWYT + szpad);
   cudaMalloc((void**)&WYTrehilo_d,szWYT + szpad);
   cudaMalloc((void**)&WYTrelolo_d,szWYT + szpad);
   cudaMalloc((void**)&WYTimhihi_d,szWYT + szpad);
   cudaMalloc((void**)&WYTimlohi_d,szWYT + szpad);
   cudaMalloc((void**)&WYTimhilo_d,szWYT + szpad);
   cudaMalloc((void**)&WYTimlolo_d,szWYT + szpad);
   cudaMalloc((void**)&Qrehihi_d,szWYT + szpad);
   cudaMalloc((void**)&Qrelohi_d,szWYT + szpad);
   cudaMalloc((void**)&Qrehilo_d,szWYT + szpad);
   cudaMalloc((void**)&Qrelolo_d,szWYT + szpad);
   cudaMalloc((void**)&Qimhihi_d,szWYT + szpad);
   cudaMalloc((void**)&Qimlohi_d,szWYT + szpad);
   cudaMalloc((void**)&Qimhilo_d,szWYT + szpad);
   cudaMalloc((void**)&Qimlolo_d,szWYT + szpad);
   cudaMemcpy(Qrehihi_d,Qrehihi_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qrelohi_d,Qrelohi_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qrehilo_d,Qrehilo_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qrelolo_d,Qrelolo_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qimhihi_d,Qimhihi_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qimlohi_d,Qimlohi_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qimhilo_d,Qimhilo_h,szWYT,cudaMemcpyHostToDevice);
   cudaMemcpy(Qimlolo_d,Qimlolo_h,szWYT,cudaMemcpyHostToDevice);
   cudaMalloc((void**)&QWYTrehihi_d,szWYT + szpad);
   cudaMalloc((void**)&QWYTrelohi_d,szWYT + szpad);
   cudaMalloc((void**)&QWYTrehilo_d,szWYT + szpad);
   cudaMalloc((void**)&QWYTrelolo_d,szWYT + szpad);
   cudaMalloc((void**)&QWYTimhihi_d,szWYT + szpad);
   cudaMalloc((void**)&QWYTimlohi_d,szWYT + szpad);
   cudaMalloc((void**)&QWYTimhilo_d,szWYT + szpad);
   cudaMalloc((void**)&QWYTimlolo_d,szWYT + szpad);

   const size_t szYWT = nrows2*sizeof(double);
   cudaMalloc((void**)&YWTrehihi_d,szYWT + szpad); // padding for Y*W^T
   cudaMalloc((void**)&YWTrelohi_d,szYWT + szpad);
   cudaMalloc((void**)&YWTrehilo_d,szYWT + szpad);
   cudaMalloc((void**)&YWTrelolo_d,szYWT + szpad);
   cudaMalloc((void**)&YWTimhihi_d,szYWT + szpad);
   cudaMalloc((void**)&YWTimlohi_d,szYWT + szpad);
   cudaMalloc((void**)&YWTimhilo_d,szYWT + szpad);
   cudaMalloc((void**)&YWTimlolo_d,szYWT + szpad);
   cudaMalloc((void**)&YWTCrehihi_d,sznum + szpad);
   cudaMalloc((void**)&YWTCrelohi_d,sznum + szpad);
   cudaMalloc((void**)&YWTCrehilo_d,sznum + szpad);
   cudaMalloc((void**)&YWTCrelolo_d,sznum + szpad);
   cudaMalloc((void**)&YWTCimhihi_d,sznum + szpad);
   cudaMalloc((void**)&YWTCimlohi_d,sznum + szpad);
   cudaMalloc((void**)&YWTCimhilo_d,sznum + szpad);
   cudaMalloc((void**)&YWTCimlolo_d,sznum + szpad);

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
            GPU_cmplx4_small_house
               (nrows,ncols,szt,nbt,colidx,nrows1,k,L,
                 Arehihi_h, Arelohi_h, Arehilo_h, Arelolo_h,
                 Aimhihi_h, Aimlohi_h, Aimhilo_h, Aimlolo_h,
                 Arehihi_d, Arelohi_d, Arehilo_d, Arelolo_d,
                 Aimhihi_d, Aimlohi_d, Aimhilo_d, Aimlolo_d,
                 vrehihi_h, vrelohi_h, vrehilo_h, vrelolo_h,
                 vimhihi_h, vimlohi_h, vimhilo_h, vimlolo_h,
                 Vrehihi_d, Vrelohi_d, Vrehilo_d, Vrelolo_d,
                 Vimhihi_d, Vimlohi_d, Vimhilo_d, Vimlolo_d,
                betahihi_h,betalohi_h,betahilo_h,betalolo_h,
                betahihi_d,betalohi_d,betahilo_d,betalolo_d,
                houselapms,addcnt,mulcnt,divcnt,sqrtcnt,verbose);

            if((betahihi_h[L] == 0.0) && (betalohi_h[L] == 0.0) &&
               (betahilo_h[L] == 0.0) && (betalolo_h[L] == 0.0))
            {
               if(verbose) cout << "Zero beta detected." << endl;
            }
            else
            {
               GPU_cmplx4_small_leftRupdate
                  (nrows,ncols,szt,colidx,k,L,
                    Arehihi_h, Arelohi_h, Arehilo_h, Arelolo_h,
                    Aimhihi_h, Aimlohi_h, Aimhilo_h, Aimlolo_h,
                    Arehihi_d, Arelohi_d, Arehilo_d, Arelolo_d,
                    Aimhihi_d, Aimlohi_d, Aimhilo_d, Aimlolo_d,
                    Vrehihi_d, Vrelohi_d, Vrehilo_d, Vrelolo_d,
                    Vimhihi_d, Vimlohi_d, Vimhilo_d, Vimlolo_d,
                   betahihi_h,betalohi_h,betahilo_h,betalolo_h,
                   betahihi_d,betalohi_d,betahilo_d,betalolo_d,
                   tileRlapms,addcnt,mulcnt,verbose);
            }
         }
         else
         {
            GPU_cmplx4_large_house
               (nrows,ncols,szt,nbt,colidx,nrows1,k,L,
                   Arehihi_h,   Arelohi_h,   Arehilo_h,   Arelolo_h,
                   Aimhihi_h,   Aimlohi_h,   Aimhilo_h,   Aimlolo_h,
                   Arehihi_d,   Arelohi_d,   Arehilo_d,   Arelolo_d,
                   Aimhihi_d,   Aimlohi_d,   Aimhilo_d,   Aimlolo_d,
                   vrehihi_h,   vrelohi_h,   vrehilo_h,   vrelolo_h,
                   vimhihi_h,   vimlohi_h,   vimhilo_h,   vimlolo_h,
                   Vrehihi_d,   Vrelohi_d,   Vrehilo_d,   Vrelolo_d,
                   Vimhihi_d,   Vimlohi_d,   Vimhilo_d,   Vimlolo_d,
                  betahihi_h,  betalohi_h,  betahilo_h,  betalolo_h,
                  betahihi_d,  betalohi_d,  betahilo_d,  betalolo_d,
                  sumshihi_h,  sumslohi_h,  sumshilo_h,  sumslolo_h,
                  sumshihi_d,  sumslohi_d,  sumshilo_d,  sumslolo_d,
                &sigmahihi_h,&sigmalohi_h,&sigmahilo_h,&sigmalolo_h,
                 sigmahihi_d, sigmalohi_d, sigmahilo_d, sigmalolo_d,
                houselapms,addcnt,mulcnt,divcnt,sqrtcnt,verbose);

            if((betahihi_h[L] == 0.0) && (betalohi_h[L] == 0.0) &&
               (betahilo_h[L] == 0.0) && (betalolo_h[L] == 0.0))
            {
               if(verbose) cout << "Zero beta detected." << endl;
            }
            else
            {
               GPU_cmplx4_medium_leftRupdate
                  (nrows,ncols,szt,colidx,k,L,
                      Arehihi_h,     Arelohi_h,     Arehilo_h,     Arelolo_h,
                      Aimhihi_h,     Aimlohi_h,     Aimhilo_h,     Aimlolo_h,
                      Arehihi_d,     Arelohi_d,     Arehilo_d,     Arelolo_d,
                      Aimhihi_d,     Aimlohi_d,     Aimhilo_d,     Aimlolo_d,
                      Vrehihi_d,     Vrelohi_d,     Vrehilo_d,     Vrelolo_d,
                      Vimhihi_d,     Vimlohi_d,     Vimhilo_d,     Vimlolo_d,
                     betahihi_h,    betalohi_h,    betahilo_h,    betalolo_h,
                     betahihi_d,    betalohi_d,    betahilo_d,    betalolo_d,
                 RHdotvrehihi_h,RHdotvrelohi_h,RHdotvrehilo_h,RHdotvrelolo_h,
                 RHdotvimhihi_h,RHdotvimlohi_h,RHdotvimhilo_h,RHdotvimlolo_h,
                 RHdotvrehihi_d,RHdotvrelohi_d,RHdotvrehilo_d,RHdotvrelolo_d,
                 RHdotvimhihi_d,RHdotvimlohi_d,RHdotvimhilo_d,RHdotvimlolo_d,
                   bRHvrehihi_h,  bRHvrelohi_h,  bRHvrehilo_h,  bRHvrelolo_h,
                   bRHvimhihi_h,  bRHvimlohi_h,  bRHvimhilo_h,  bRHvimlolo_h,
                   bRHvrehihi_d,  bRHvrelohi_d,  bRHvrehilo_d,  bRHvrelolo_d,
                   bRHvimhihi_d,  bRHvimlohi_d,  bRHvimhilo_d,  bRHvimlolo_d,
                 RHvlapms,tileRlapms,addcnt,mulcnt,verbose);
            }
         }
      }
      GPU_cmplx4_medium_VB_to_W
         (nrows,szt,szt,k,
            Vrehihi_h,  Vrelohi_h,  Vrehilo_h,  Vrelolo_h,
            Vimhihi_h,  Vimlohi_h,  Vimhilo_h,  Vimlolo_h,
            Vrehihi_d,  Vrelohi_d,  Vrehilo_d,  Vrelolo_d,
            Vimhihi_d,  Vimlohi_d,  Vimhilo_d,  Vimlolo_d,
            Wrehihi_h,  Wrelohi_h,  Wrehilo_h,  Wrelolo_h,
            Wimhihi_h,  Wimlohi_h,  Wimhilo_h,  Wimlolo_h,
            Wrehihi_d,  Wrelohi_d,  Wrehilo_d,  Wrelolo_d,
            Wimhihi_d,  Wimlohi_d,  Wimhilo_d,  Wimlolo_d,
          WYTrehihi_h,WYTrelohi_h,WYTrehilo_h,WYTrelolo_h,
          WYTimhihi_h,WYTimlohi_h,WYTimhilo_h,WYTimlolo_h,
          WYTrehihi_d,WYTrelohi_d,WYTrehilo_d,WYTrelolo_d,
          WYTimhihi_d,WYTimlohi_d,WYTimhilo_d,WYTimlolo_d,
           betahihi_h, betalohi_h, betahilo_h, betalolo_h,
           betahihi_d, betalohi_d, betahilo_d, betalolo_d,
          vb2Wlapms,addcnt,mulcnt,verbose);

/*
      GPU_cmplx2_small_WYH(nrows-k*szt,szt,
            Wrehi_d,  Wrelo_d,  Wimhi_d,  Wimlo_d,
            Vrehi_d,  Vrelo_d,  Vimhi_d,  Vimlo_d,
          WYTrehi_d,WYTrelo_d,WYTimhi_d,WYTimlo_d,
          WYTrehi_h,WYTrelo_h,WYTimhi_h,WYTimlo_h,WYHlapms,verbose);
 */

      GPU_cmplx4_small_QWYH(nrows,szt,k,
             Qrehihi_d,   Qrelohi_d,   Qrehilo_d,   Qrelolo_d,
             Qimhihi_d,   Qimlohi_d,   Qimhilo_d,   Qimlolo_d,
           WYTrehihi_d, WYTrelohi_d, WYTrehilo_d, WYTrelolo_d,
           WYTimhihi_d, WYTimlohi_d, WYTimhilo_d, WYTimlolo_d,
          QWYTrehihi_d,QWYTrelohi_d,QWYTrehilo_d,QWYTrelolo_d,
          QWYTimhihi_d,QWYTimlohi_d,QWYTimhilo_d,QWYTimlolo_d,
          QWYTrehihi_h,QWYTrelohi_h,QWYTrehilo_h,QWYTrelolo_h,
          QWYTimhihi_h,QWYTimlohi_h,QWYTimhilo_h,QWYTimlolo_h,
             Qrehihi_h,   Qrelohi_h,   Qrehilo_h,   Qrelolo_h,
             Qimhihi_h,   Qimlohi_h,   Qimhilo_h,   Qimlolo_h,
          QWYHlapms,addcnt,mulcnt,verbose);

      GPU_cmplx4_small_Qupdate
         (nrows,szt,k,
             Qrehihi_d,   Qrelohi_d,   Qrehilo_d,   Qrelolo_d,
             Qimhihi_d,   Qimlohi_d,   Qimhilo_d,   Qimlolo_d,
          QWYTrehihi_d,QWYTrelohi_d,QWYTrehilo_d,QWYTrelolo_d,
          QWYTimhihi_d,QWYTimlohi_d,QWYTimhilo_d,QWYTimlolo_d,
             Qrehihi_h,   Qrelohi_h,   Qrehilo_h,   Qrelolo_h,
             Qimhihi_h,   Qimlohi_h,   Qimhilo_h,   Qimlolo_h,
          Qaddlapms,addcnt,verbose);

      if(k < nbt-1)                              // update R
      {
         GPU_cmplx4_small_YWH
            (nrows,szt,k,
               Vrehihi_d,  Vrelohi_d,  Vrehilo_d,  Vrelolo_d,
               Vimhihi_d,  Vimlohi_d,  Vimhilo_d,  Vimlolo_d,
               Wrehihi_d,  Wrelohi_d,  Wrehilo_d,  Wrelolo_d,
               Wimhihi_d,  Wimlohi_d,  Wimhilo_d,  Wimlolo_d,
             YWTrehihi_d,YWTrelohi_d,YWTrehilo_d,YWTrelolo_d,
             YWTimhihi_d,YWTimlohi_d,YWTimhilo_d,YWTimlolo_d,
             YWTrehihi_h,YWTrelohi_h,YWTrehilo_h,YWTrelolo_h,
             YWTimhihi_h,YWTimlohi_h,YWTimhilo_h,YWTimlolo_h,
             YWHlapms,addcnt,mulcnt,verbose);

         GPU_cmplx4_small_YWHC
            (nrows,ncols,szt,k,
              YWTrehihi_d, YWTrelohi_d, YWTrehilo_d, YWTrelolo_d,
              YWTimhihi_d, YWTimlohi_d, YWTimhilo_d, YWTimlolo_d,
                Arehihi_d,   Arelohi_d,   Arehilo_d,   Arelolo_d,
                Aimhihi_d,   Aimlohi_d,   Aimhilo_d,   Aimlolo_d,
             YWTCrehihi_d,YWTCrelohi_d,YWTCrehilo_d,YWTCrelolo_d,
             YWTCimhihi_d,YWTCimlohi_d,YWTCimhilo_d,YWTCimlolo_d,
             YWTCrehihi_h,YWTCrelohi_h,YWTCrehilo_h,YWTCrelolo_h,
             YWTCimhihi_h,YWTCimlohi_h,YWTCimhilo_h,YWTCimlolo_h,
             YWHClapms,addcnt,mulcnt,verbose);

         GPU_cmplx4_small_R_add_YWHC
            (nrows,ncols,szt,k,
                Arehihi_d,   Arelohi_d,   Arehilo_d,   Arelolo_d,
                Aimhihi_d,   Aimlohi_d,   Aimhilo_d,   Aimlolo_d,
             YWTCrehihi_d,YWTCrelohi_d,YWTCrehilo_d,YWTCrelolo_d,
             YWTCimhihi_d,YWTCimlohi_d,YWTCimhilo_d,YWTCimlolo_d,
                Arehihi_h,   Arelohi_h,   Arehilo_h,   Arelolo_h,
                Aimhihi_h,   Aimlohi_h,   Aimhilo_h,   Aimlolo_h,
             Raddlapms,addcnt,verbose);
      }
   }
   gettimeofday(&endtime,0);
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   *walltimesec = seconds + microseconds*1.0e-6;

   cudaMemcpy(Qrehihi_h,Qrehihi_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qrelohi_h,Qrelohi_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qrehilo_h,Qrehilo_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qrelolo_h,Qrelolo_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qimhihi_h,Qimhihi_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qimlohi_h,Qimlohi_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qimhilo_h,Qimhilo_d,szWYT,cudaMemcpyDeviceToHost);
   cudaMemcpy(Qimlolo_h,Qimlolo_d,szWYT,cudaMemcpyDeviceToHost);
   ix = 0;                                           // copy rows of Q
   for(int i=0; i<nrows; i++)
      for(int j=0; j<nrows; j++)
      {
         Qrehihi[i][j] = Qrehihi_h[ix];
         Qrelohi[i][j] = Qrelohi_h[ix];
         Qrehilo[i][j] = Qrehilo_h[ix];
         Qrelolo[i][j] = Qrelolo_h[ix];
         Qimhihi[i][j] = Qimhihi_h[ix];
         Qimlohi[i][j] = Qimlohi_h[ix];
         Qimhilo[i][j] = Qimhilo_h[ix];
         Qimlolo[i][j] = Qimlolo_h[ix++];
      }

   cudaMemcpy(Arehihi_h,Arehihi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Arelohi_h,Arelohi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Arehilo_h,Arehilo_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Arelolo_h,Arelolo_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Aimhihi_h,Aimhihi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Aimlohi_h,Aimlohi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Aimhilo_h,Aimhilo_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(Aimlolo_h,Aimlolo_d,sznum,cudaMemcpyDeviceToHost);

   for(int i=0; i<nrows; i++)                       // copy columns of R
      for(int j=0; j<ncols; j++)
      {
         Rrehihi[i][j] = Arehihi_h[j*nrows+i];
         Rrelohi[i][j] = Arelohi_h[j*nrows+i];
         Rrehilo[i][j] = Arehilo_h[j*nrows+i];
         Rrelolo[i][j] = Arelolo_h[j*nrows+i];
         Rimhihi[i][j] = Aimhihi_h[j*nrows+i];
         Rimlohi[i][j] = Aimlohi_h[j*nrows+i];
         Rimhilo[i][j] = Aimhilo_h[j*nrows+i];
         Rimlolo[i][j] = Aimlolo_h[j*nrows+i];
      }

   free(Arehihi_h); free(Arelohi_h); free(Arehilo_h); free(Arelolo_h);
   free(Aimhihi_h); free(Aimlohi_h); free(Aimhilo_h); free(Aimlolo_h);
   free(Qrehihi_h); free(Qrelohi_h); free(Qrehilo_h); free(Qrelolo_h);
   free(Qimhihi_h); free(Qimlohi_h); free(Qimhilo_h); free(Qimlolo_h);
   free(vrehihi_h); free(vrelohi_h); free(vrehilo_h); free(vrelolo_h); 
   free(vimhihi_h); free(vimlohi_h); free(vimhilo_h); free(vimlolo_h);
   free(Vrehihi_h); free(Vrelohi_h); free(Vrehilo_h); free(Vrelolo_h);
   free(Vimhihi_h); free(Vimlohi_h); free(Vimhilo_h); free(Vimlolo_h);
   free(RHdotvrehihi_h); free(RHdotvrelohi_h);
   free(RHdotvrehilo_h); free(RHdotvrelolo_h);
   free(RHdotvimhihi_h); free(RHdotvimlohi_h);
   free(RHdotvimhilo_h); free(RHdotvimlolo_h);
   free(bRHvrehihi_h); free(bRHvrelohi_h);
   free(bRHvrehilo_h); free(bRHvrelolo_h);
   free(bRHvimhihi_h); free(bRHvimlohi_h);
   free(bRHvimhilo_h); free(bRHvimlolo_h);
   free(sumshihi_h); free(sumslohi_h);
   free(sumshilo_h); free(sumslolo_h);
   free(Wrehihi_h); free(Wrelohi_h); free(Wrehilo_h); free(Wrelolo_h);
   free(Wimhihi_h); free(Wimlohi_h); free(Wimhilo_h); free(Wimlolo_h);
   free(WYTrehihi_h); free(WYTrelohi_h); free(WYTrehilo_h); free(WYTrelolo_h);
   free(WYTimhihi_h); free(WYTimlohi_h); free(WYTimhilo_h); free(WYTimlolo_h);
   free(QWYTrehihi_h); free(QWYTrelohi_h);
   free(QWYTrehilo_h); free(QWYTrelolo_h);
   free(QWYTimhihi_h); free(QWYTimlohi_h);
   free(QWYTimhilo_h); free(QWYTimlolo_h);
   free(YWTrehihi_h); free(YWTrelohi_h); free(YWTrehilo_h); free(YWTrelolo_h);
   free(YWTimhihi_h); free(YWTimlohi_h); free(YWTimhilo_h); free(YWTimlolo_h);
   free(YWTCrehihi_h); free(YWTCrelohi_h);
   free(YWTCrehilo_h); free(YWTCrelolo_h);
   free(YWTCimhihi_h); free(YWTCimlohi_h);
   free(YWTCimhilo_h); free(YWTCimlolo_h);
}
