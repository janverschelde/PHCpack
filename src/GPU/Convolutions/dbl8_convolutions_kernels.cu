// The file dbl8_convolutions_kernels.cu defines kernels with prototypes
// in dbl8_convolution_kernels.h.

#ifdef gpufun
#include "double_double_gpufun.cu"
#include "quad_double_gpufun.cu"
#include "octo_double_gpufun.cu"
#endif
#include "dbl8_convolutions_kernels.h"

__global__ void dbl8_increment
 ( double *xhihihi, double *xlohihi, double *xhilohi, double *xlolohi,
   double *xhihilo, double *xlohilo, double *xhilolo, double *xlololo,
   double *yhihihi, double *ylohihi, double *yhilohi, double *ylolohi,
   double *yhihilo, double *ylohilo, double *yhilolo, double *ylololo,
   double *zhihihi, double *zlohihi, double *zhilohi, double *zlolohi,
   double *zhihilo, double *zlohilo, double *zhilolo, double *zlololo,
   int dim )
{
   int k = threadIdx.x;                 // thread k computes z[k]

   __shared__ double xvhihihi[od_shmemsize];
   __shared__ double xvlohihi[od_shmemsize];
   __shared__ double xvhilohi[od_shmemsize];
   __shared__ double xvlolohi[od_shmemsize];
   __shared__ double xvhihilo[od_shmemsize];
   __shared__ double xvlohilo[od_shmemsize];
   __shared__ double xvhilolo[od_shmemsize];
   __shared__ double xvlololo[od_shmemsize];
   __shared__ double yvhihihi[od_shmemsize];
   __shared__ double yvlohihi[od_shmemsize];
   __shared__ double yvhilohi[od_shmemsize];
   __shared__ double yvlolohi[od_shmemsize];
   __shared__ double yvhihilo[od_shmemsize];
   __shared__ double yvlohilo[od_shmemsize];
   __shared__ double yvhilolo[od_shmemsize];
   __shared__ double yvlololo[od_shmemsize];
   __shared__ double zvhihihi[od_shmemsize];
   __shared__ double zvlohihi[od_shmemsize];
   __shared__ double zvhilohi[od_shmemsize];
   __shared__ double zvlolohi[od_shmemsize];
   __shared__ double zvhihilo[od_shmemsize];
   __shared__ double zvlohilo[od_shmemsize];
   __shared__ double zvhilolo[od_shmemsize];
   __shared__ double zvlololo[od_shmemsize];

   xvhihihi[k] = xhihihi[k]; xvlohihi[k] = xlohihi[k];
   xvhilohi[k] = xhilohi[k]; xvlolohi[k] = xlolohi[k];
   xvhihilo[k] = xhihilo[k]; xvlohilo[k] = xlohilo[k];
   xvhilolo[k] = xhilolo[k]; xvlololo[k] = xlololo[k];
   yvhihihi[k] = yhihihi[k]; yvlohihi[k] = ylohihi[k];
   yvhilohi[k] = yhilohi[k]; yvlolohi[k] = ylolohi[k];
   yvhihilo[k] = yhihilo[k]; yvlohilo[k] = ylohilo[k];
   yvhilolo[k] = yhilolo[k]; yvlololo[k] = ylololo[k];

   __syncthreads();

   odg_add( xvhihihi[k], xvlohihi[k], xvhilohi[k], xvlolohi[k],
            xvhihilo[k], xvlohilo[k], xvhilolo[k], xvlololo[k],
            yvhihihi[k], yvlohihi[k], yvhilohi[k], yvlolohi[k],
            yvhihilo[k], yvlohilo[k], yvhilolo[k], yvlololo[k],
           &zvhihihi[k],&zvlohihi[k],&zvhilohi[k],&zvlolohi[k],
           &zvhihilo[k],&zvlohilo[k],&zvhilolo[k],&zvlololo[k]);

   __syncthreads();

   zhihihi[k] = zvhihihi[k]; zlohihi[k] = zvlohihi[k];
   zhilohi[k] = zvhilohi[k]; zlolohi[k] = zvlolohi[k];
   zhihilo[k] = zvhihilo[k]; zlohilo[k] = zvlohilo[k];
   zhilolo[k] = zvhilolo[k]; zlololo[k] = zvlololo[k];
}

__global__ void dbl8_decrement
 ( double *xhihihi, double *xlohihi, double *xhilohi, double *xlolohi,
   double *xhihilo, double *xlohilo, double *xhilolo, double *xlololo,
   double *yhihihi, double *ylohihi, double *yhilohi, double *ylolohi,
   double *yhihilo, double *ylohilo, double *yhilolo, double *ylololo,
   double *zhihihi, double *zlohihi, double *zhilohi, double *zlolohi,
   double *zhihilo, double *zlohilo, double *zhilolo, double *zlololo,
   int dim )
{
   int k = threadIdx.x;                 // thread k computes z[k]

   __shared__ double xvhihihi[od_shmemsize];
   __shared__ double xvlohihi[od_shmemsize];
   __shared__ double xvhilohi[od_shmemsize];
   __shared__ double xvlolohi[od_shmemsize];
   __shared__ double xvhihilo[od_shmemsize];
   __shared__ double xvlohilo[od_shmemsize];
   __shared__ double xvhilolo[od_shmemsize];
   __shared__ double xvlololo[od_shmemsize];
   __shared__ double yvhihihi[od_shmemsize];
   __shared__ double yvlohihi[od_shmemsize];
   __shared__ double yvhilohi[od_shmemsize];
   __shared__ double yvlolohi[od_shmemsize];
   __shared__ double yvhihilo[od_shmemsize];
   __shared__ double yvlohilo[od_shmemsize];
   __shared__ double yvhilolo[od_shmemsize];
   __shared__ double yvlololo[od_shmemsize];
   __shared__ double zvhihihi[od_shmemsize];
   __shared__ double zvlohihi[od_shmemsize];
   __shared__ double zvhilohi[od_shmemsize];
   __shared__ double zvlolohi[od_shmemsize];
   __shared__ double zvhihilo[od_shmemsize];
   __shared__ double zvlohilo[od_shmemsize];
   __shared__ double zvhilolo[od_shmemsize];
   __shared__ double zvlololo[od_shmemsize];

   xvhihihi[k] = xhihihi[k]; xvlohihi[k] = xlohihi[k];
   xvhilohi[k] = xhilohi[k]; xvlolohi[k] = xlolohi[k];
   xvhihilo[k] = xhihilo[k]; xvlohilo[k] = xlohilo[k];
   xvhilolo[k] = xhilolo[k]; xvlololo[k] = xlololo[k];
   yvhihihi[k] = yhihihi[k]; yvlohihi[k] = ylohihi[k];
   yvhilohi[k] = yhilohi[k]; yvlolohi[k] = ylolohi[k];
   yvhihilo[k] = yhihilo[k]; yvlohilo[k] = ylohilo[k];
   yvhilolo[k] = yhilolo[k]; yvlololo[k] = ylololo[k];

   __syncthreads();

   odg_sub( xvhihihi[k], xvlohihi[k], xvhilohi[k], xvlolohi[k],
            xvhihilo[k], xvlohilo[k], xvhilolo[k], xvlololo[k],
            yvhihihi[k], yvlohihi[k], yvhilohi[k], yvlolohi[k],
            yvhihilo[k], yvlohilo[k], yvhilolo[k], yvlololo[k],
           &zvhihihi[k],&zvlohihi[k],&zvhilohi[k],&zvlolohi[k],
           &zvhihilo[k],&zvlohilo[k],&zvhilolo[k],&zvlololo[k]);

   __syncthreads();

   zhihihi[k] = zvhihihi[k]; zlohihi[k] = zvlohihi[k];
   zhilohi[k] = zvhilohi[k]; zlolohi[k] = zvlolohi[k];
   zhihilo[k] = zvhihilo[k]; zlohilo[k] = zvlohilo[k];
   zhilolo[k] = zvhilolo[k]; zlololo[k] = zvlololo[k];
}

__global__ void dbl8_padded_convolute
 ( double *xhihihi, double *xlohihi, double *xhilohi, double *xlolohi,
   double *xhihilo, double *xlohilo, double *xhilolo, double *xlololo,
   double *yhihihi, double *ylohihi, double *yhilohi, double *ylolohi,
   double *yhihilo, double *ylohilo, double *yhilolo, double *ylololo,
   double *zhihihi, double *zlohihi, double *zhilohi, double *zlolohi,
   double *zhihilo, double *zlohilo, double *zhilolo, double *zlololo,
   int dim )
{
   int k = threadIdx.x;                 // thread k computes z[k]

   __shared__ double xvhihihi[od_shmemsize];
   __shared__ double xvlohihi[od_shmemsize];
   __shared__ double xvhilohi[od_shmemsize];
   __shared__ double xvlolohi[od_shmemsize];
   __shared__ double xvhihilo[od_shmemsize];
   __shared__ double xvlohilo[od_shmemsize];
   __shared__ double xvhilolo[od_shmemsize];
   __shared__ double xvlololo[od_shmemsize];
   __shared__ double yvhihihi[od_shmemsize];
   __shared__ double yvlohihi[od_shmemsize];
   __shared__ double yvhilohi[od_shmemsize];
   __shared__ double yvlolohi[od_shmemsize];
   __shared__ double yvhihilo[od_shmemsize];
   __shared__ double yvlohilo[od_shmemsize];
   __shared__ double yvhilolo[od_shmemsize];
   __shared__ double yvlololo[od_shmemsize];
   __shared__ double zvhihihi[od_shmemsize];
   __shared__ double zvlohihi[od_shmemsize];
   __shared__ double zvhilohi[od_shmemsize];
   __shared__ double zvlolohi[od_shmemsize];
   __shared__ double zvhihilo[od_shmemsize];
   __shared__ double zvlohilo[od_shmemsize];
   __shared__ double zvhilolo[od_shmemsize];
   __shared__ double zvlololo[od_shmemsize];

   double prdhihihi,prdlohihi,prdhilohi,prdlolohi;
   double prdhihilo,prdlohilo,prdhilolo,prdlololo;
   int idx = dim + k;

   xvhihihi[k] = xhihihi[k]; xvlohihi[k] = xlohihi[k];
   xvhilohi[k] = xhilohi[k]; xvlolohi[k] = xlolohi[k];
   xvhihilo[k] = xhihilo[k]; xvlohilo[k] = xlohilo[k];
   xvhilolo[k] = xhilolo[k]; xvlololo[k] = xlololo[k];
   yvhihihi[k] = 0.0; yvlohihi[k] = 0.0;
   yvhilohi[k] = 0.0; yvlolohi[k] = 0.0;
   yvhihilo[k] = 0.0; yvlohilo[k] = 0.0;
   yvhilolo[k] = 0.0; yvlololo[k] = 0.0;
   yvhihihi[idx] = yhihihi[k]; yvlohihi[idx] = ylohihi[k];
   yvhilohi[idx] = yhilohi[k]; yvlolohi[idx] = ylolohi[k];
   yvhihilo[idx] = yhihilo[k]; yvlohilo[idx] = ylohilo[k];
   yvhilolo[idx] = yhilolo[k]; yvlololo[idx] = ylololo[k];

   __syncthreads();

   // zv[k] = xv[0]*yv[k];
   odg_mul( xvhihihi[0], xvlohihi[0], xvhilohi[0], xvlolohi[0],
            xvhihilo[0], xvlohilo[0], xvhilolo[0], xvlololo[0],
            yvhihihi[idx], yvlohihi[idx], yvhilohi[idx], yvlolohi[idx],
            yvhihilo[idx], yvlohilo[idx], yvhilolo[idx], yvlololo[idx],
           &zvhihihi[k],&zvlohihi[k],&zvhilohi[k],&zvlolohi[k],
           &zvhihilo[k],&zvlohilo[k],&zvhilolo[k],&zvlololo[k]);

   for(int i=1; i<=k; i++) // zv[k] = zv[k] + xv[i]*yv[k-i];
   {
      idx = dim + k - i;
      odg_mul( xvhihihi[i],  xvlohihi[i],  xvhilohi[i],  xvlolohi[i],
               xvhihilo[i],  xvlohilo[i],  xvhilolo[i],  xvlololo[i],
               yvhihihi[idx],yvlohihi[idx],yvhilohi[idx],yvlolohi[idx],
               yvhihilo[idx],yvlohilo[idx],yvhilolo[idx],yvlololo[idx],
             &prdhihihi,   &prdlohihi,   &prdhilohi,   &prdlolohi,
             &prdhihilo,   &prdlohilo,   &prdhilolo,   &prdlololo);
      odg_inc(&zvhihihi[k], &zvlohihi[k], &zvhilohi[k], &zvlolohi[k],
              &zvhihilo[k] ,&zvlohilo[k], &zvhilolo[k], &zvlololo[k],
              prdhihihi,    prdlohihi,    prdhilohi,    prdlolohi,
              prdhihilo,    prdlohilo,    prdhilolo,    prdlololo);
   }

   __syncthreads();

   zhihihi[k] = zvhihihi[k]; zlohihi[k] = zvlohihi[k];
   zhilohi[k] = zvhilohi[k]; zlolohi[k] = zvlolohi[k];
   zhihilo[k] = zvhihilo[k]; zlohilo[k] = zvlohilo[k];
   zhilolo[k] = zvhilolo[k]; zlololo[k] = zvlololo[k];
}

__global__ void dbl8_convolute
 ( double *xhihihi, double *xlohihi, double *xhilohi, double *xlolohi,
   double *xhihilo, double *xlohilo, double *xhilolo, double *xlololo,
   double *yhihihi, double *ylohihi, double *yhilohi, double *ylolohi,
   double *yhihilo, double *ylohilo, double *yhilolo, double *ylololo,
   double *zhihihi, double *zlohihi, double *zhilohi, double *zlolohi,
   double *zhihilo, double *zlohilo, double *zhilolo, double *zlololo,
   int dim )
{
   int k = threadIdx.x;                 // thread k computes z[k]

   __shared__ double xvhihihi[od_shmemsize];
   __shared__ double xvlohihi[od_shmemsize];
   __shared__ double xvhilohi[od_shmemsize];
   __shared__ double xvlolohi[od_shmemsize];
   __shared__ double xvhihilo[od_shmemsize];
   __shared__ double xvlohilo[od_shmemsize];
   __shared__ double xvhilolo[od_shmemsize];
   __shared__ double xvlololo[od_shmemsize];
   __shared__ double yvhihihi[od_shmemsize];
   __shared__ double yvlohihi[od_shmemsize];
   __shared__ double yvhilohi[od_shmemsize];
   __shared__ double yvlolohi[od_shmemsize];
   __shared__ double yvhihilo[od_shmemsize];
   __shared__ double yvlohilo[od_shmemsize];
   __shared__ double yvhilolo[od_shmemsize];
   __shared__ double yvlololo[od_shmemsize];
   __shared__ double zvhihihi[od_shmemsize];
   __shared__ double zvlohihi[od_shmemsize];
   __shared__ double zvhilohi[od_shmemsize];
   __shared__ double zvlolohi[od_shmemsize];
   __shared__ double zvhihilo[od_shmemsize];
   __shared__ double zvlohilo[od_shmemsize];
   __shared__ double zvhilolo[od_shmemsize];
   __shared__ double zvlololo[od_shmemsize];

   double prdhihihi,prdlohihi,prdhilohi,prdlolohi;
   double prdhihilo,prdlohilo,prdhilolo,prdlololo;

   xvhihihi[k] = xhihihi[k]; xvlohihi[k] = xlohihi[k];
   xvhilohi[k] = xhilohi[k]; xvlolohi[k] = xlolohi[k];
   xvhihilo[k] = xhihilo[k]; xvlohilo[k] = xlohilo[k];
   xvhilolo[k] = xhilolo[k]; xvlololo[k] = xlololo[k];
   yvhihihi[k] = yhihihi[k]; yvlohihi[k] = ylohihi[k];
   yvhilohi[k] = yhilohi[k]; yvlolohi[k] = ylolohi[k];
   yvhihilo[k] = yhihilo[k]; yvlohilo[k] = ylohilo[k];
   yvhilolo[k] = yhilolo[k]; yvlololo[k] = ylololo[k];

   __syncthreads();

   // zv[k] = xv[0]*yv[k];
   odg_mul( xvhihihi[0], xvlohihi[0], xvhilohi[0], xvlolohi[0],
            xvhihilo[0], xvlohilo[0], xvhilolo[0], xvlololo[0],
            yvhihihi[k], yvlohihi[k], yvhilohi[k], yvlolohi[k],
            yvhihilo[k], yvlohilo[k], yvhilolo[k], yvlololo[k],
           &zvhihihi[k],&zvlohihi[k],&zvhilohi[k],&zvlolohi[k],
           &zvhihilo[k],&zvlohilo[k],&zvhilolo[k],&zvlololo[k]);

   for(int i=1; i<=k; i++) // zv[k] = zv[k] + xv[i]*yv[k-i];
   {
      odg_mul( xvhihihi[i],  xvlohihi[i],  xvhilohi[i],  xvlolohi[i],
               xvhihilo[i],  xvlohilo[i],  xvhilolo[i],  xvlololo[i],
               yvhihihi[k-i],yvlohihi[k-i],yvhilohi[k-i],yvlolohi[k-i],
               yvhihilo[k-i],yvlohilo[k-i],yvhilolo[k-i],yvlololo[k-i],
             &prdhihihi,   &prdlohihi,   &prdhilohi,   &prdlolohi,
             &prdhihilo,   &prdlohilo,   &prdhilolo,   &prdlololo);
      odg_inc(&zvhihihi[k], &zvlohihi[k], &zvhilohi[k], &zvlolohi[k],
              &zvhihilo[k] ,&zvlohilo[k], &zvhilolo[k], &zvlololo[k],
              prdhihihi,    prdlohihi,    prdhilohi,    prdlolohi,
              prdhihilo,    prdlohilo,    prdhilolo,    prdlololo);
   }

   __syncthreads();

   zhihihi[k] = zvhihihi[k]; zlohihi[k] = zvlohihi[k];
   zhilohi[k] = zvhilohi[k]; zlolohi[k] = zvlolohi[k];
   zhihilo[k] = zvhihilo[k]; zlohilo[k] = zvlohilo[k];
   zhilolo[k] = zvhilolo[k]; zlololo[k] = zvlololo[k];
}

__global__ void cmplx8_convolute
 ( double *xrehihihi, double *xrelohihi, double *xrehilohi, double *xrelolohi,
   double *xrehihilo, double *xrelohilo, double *xrehilolo, double *xrelololo,
   double *ximhihihi, double *ximlohihi, double *ximhilohi, double *ximlolohi,
   double *ximhihilo, double *ximlohilo, double *ximhilolo, double *ximlololo,
   double *yrehihihi, double *yrelohihi, double *yrehilohi, double *yrelolohi,
   double *yrehihilo, double *yrelohilo, double *yrehilolo, double *yrelololo,
   double *yimhihihi, double *yimlohihi, double *yimhilohi, double *yimlolohi,
   double *yimhihilo, double *yimlohilo, double *yimhilolo, double *yimlololo,
   double *zrehihihi, double *zrelohihi, double *zrehilohi, double *zrelolohi,
   double *zrehihilo, double *zrelohilo, double *zrehilolo, double *zrelololo,
   double *zimhihihi, double *zimlohihi, double *zimhilohi, double *zimlolohi,
   double *zimhihilo, double *zimlohilo, double *zimhilolo, double *zimlololo,
   int dim )
{
   int k = threadIdx.x;       // thread k computes zre[k] and zim[k]

   __shared__ double xvrehihihi[od_shmemsize];
   __shared__ double xvrelohihi[od_shmemsize];
   __shared__ double xvrehilohi[od_shmemsize];
   __shared__ double xvrelolohi[od_shmemsize];
   __shared__ double xvrehihilo[od_shmemsize];
   __shared__ double xvrelohilo[od_shmemsize];
   __shared__ double xvrehilolo[od_shmemsize];
   __shared__ double xvrelololo[od_shmemsize];
   __shared__ double xvimhihihi[od_shmemsize];
   __shared__ double xvimlohihi[od_shmemsize];
   __shared__ double xvimhilohi[od_shmemsize];
   __shared__ double xvimlolohi[od_shmemsize];
   __shared__ double xvimhihilo[od_shmemsize];
   __shared__ double xvimlohilo[od_shmemsize];
   __shared__ double xvimhilolo[od_shmemsize];
   __shared__ double xvimlololo[od_shmemsize];
   __shared__ double yvrehihihi[od_shmemsize];
   __shared__ double yvrelohihi[od_shmemsize];
   __shared__ double yvrehilohi[od_shmemsize];
   __shared__ double yvrelolohi[od_shmemsize];
   __shared__ double yvrehihilo[od_shmemsize];
   __shared__ double yvrelohilo[od_shmemsize];
   __shared__ double yvrehilolo[od_shmemsize];
   __shared__ double yvrelololo[od_shmemsize];
   __shared__ double yvimhihihi[od_shmemsize];
   __shared__ double yvimlohihi[od_shmemsize];
   __shared__ double yvimhilohi[od_shmemsize];
   __shared__ double yvimlolohi[od_shmemsize];
   __shared__ double yvimhihilo[od_shmemsize];
   __shared__ double yvimlohilo[od_shmemsize];
   __shared__ double yvimhilolo[od_shmemsize];
   __shared__ double yvimlololo[od_shmemsize];
   __shared__ double zvrehihihi[od_shmemsize];
   __shared__ double zvrelohihi[od_shmemsize];
   __shared__ double zvrehilohi[od_shmemsize];
   __shared__ double zvrelolohi[od_shmemsize];
   __shared__ double zvrehihilo[od_shmemsize];
   __shared__ double zvrelohilo[od_shmemsize];
   __shared__ double zvrehilolo[od_shmemsize];
   __shared__ double zvrelololo[od_shmemsize];
   __shared__ double zvimhihihi[od_shmemsize];
   __shared__ double zvimlohihi[od_shmemsize];
   __shared__ double zvimhilohi[od_shmemsize];
   __shared__ double zvimlolohi[od_shmemsize];
   __shared__ double zvimhihilo[od_shmemsize];
   __shared__ double zvimlohilo[od_shmemsize];
   __shared__ double zvimhilolo[od_shmemsize];
   __shared__ double zvimlololo[od_shmemsize];

   double xrhihihi,xihihihi,yrhihihi,yihihihi,zrhihihi,zihihihi,acchihihi;
   double xrhihilo,xihihilo,yrhihilo,yihihilo,zrhihilo,zihihilo,acchihilo;
   double xrlohihi,xilohihi,yrlohihi,yilohihi,zrlohihi,zilohihi,acclohihi;
   double xrlohilo,xilohilo,yrlohilo,yilohilo,zrlohilo,zilohilo,acclohilo;
   double xrhilohi,xihilohi,yrhilohi,yihilohi,zrhilohi,zihilohi,acchilohi;
   double xrhilolo,xihilolo,yrhilolo,yihilolo,zrhilolo,zihilolo,acchilolo;
   double xrlolohi,xilolohi,yrlolohi,yilolohi,zrlolohi,zilolohi,acclolohi;
   double xrlololo,xilololo,yrlololo,yilololo,zrlololo,zilololo,acclololo;

   xvrehihihi[k] = xrehihihi[k]; xvimhihihi[k] = ximhihihi[k];
   xvrelohihi[k] = xrelohihi[k]; xvimlohihi[k] = ximlohihi[k];
   xvrehilohi[k] = xrehilohi[k]; xvimhilohi[k] = ximhilohi[k];
   xvrelolohi[k] = xrelolohi[k]; xvimlolohi[k] = ximlolohi[k];
   xvrehihilo[k] = xrehihilo[k]; xvimhihilo[k] = ximhihilo[k];
   xvrelohilo[k] = xrelohilo[k]; xvimlohilo[k] = ximlohilo[k];
   xvrehilolo[k] = xrehilolo[k]; xvimhilolo[k] = ximhilolo[k];
   xvrelololo[k] = xrelololo[k]; xvimlololo[k] = ximlololo[k];
   yvrehihihi[k] = yrehihihi[k]; yvimhihihi[k] = yimhihihi[k];
   yvrelohihi[k] = yrelohihi[k]; yvimlohihi[k] = yimlohihi[k];
   yvrehilohi[k] = yrehilohi[k]; yvimhilohi[k] = yimhilohi[k];
   yvrelolohi[k] = yrelolohi[k]; yvimlolohi[k] = yimlolohi[k];
   yvrehihilo[k] = yrehihilo[k]; yvimhihilo[k] = yimhihilo[k];
   yvrelohilo[k] = yrelohilo[k]; yvimlohilo[k] = yimlohilo[k];
   yvrehilolo[k] = yrehilolo[k]; yvimhilolo[k] = yimhilolo[k];
   yvrelololo[k] = yrelololo[k]; yvimlololo[k] = yimlololo[k];

   __syncthreads();

   // z[k] = x[0]*y[k]
   xrhihihi = xvrehihihi[0]; xrlohihi = xvrelohihi[0];
   xrhilohi = xvrehilohi[0]; xrlolohi = xvrelolohi[0];
   xrhihilo = xvrehihilo[0]; xrlohilo = xvrelohilo[0];
   xrhilolo = xvrehilolo[0]; xrlololo = xvrelololo[0];
   xihihihi = xvimhihihi[0]; xilohihi = xvimlohihi[0];
   xihilohi = xvimhilohi[0]; xilolohi = xvimlolohi[0];
   xihihilo = xvimhihilo[0]; xilohilo = xvimlohilo[0];
   xihilolo = xvimhilolo[0]; xilololo = xvimlololo[0];
   yrhihihi = yvrehihihi[k]; yrlohihi = yvrelohihi[k];
   yrhilohi = yvrehilohi[k]; yrlolohi = yvrelolohi[k];
   yrhihilo = yvrehihilo[k]; yrlohilo = yvrelohilo[k];
   yrhilolo = yvrehilolo[k]; yrlololo = yvrelololo[k];
   yihihihi = yvimhihihi[k]; yilohihi = yvimlohihi[k];
   yihilohi = yvimhilohi[k]; yilolohi = yvimlolohi[k];
   yihihilo = yvimhihilo[k]; yilohilo = yvimlohilo[k];
   yihilolo = yvimhilolo[k]; yilololo = yvimlololo[k];

   odg_mul(xrhihihi,xrlohihi,xrhilohi,xrlolohi,
           xrhihilo,xrlohilo,xrhilolo,xrlololo,
           yrhihihi,yrlohihi,yrhilohi,yrlolohi,
           yrhihilo,yrlohilo,yrhilolo,yrlololo,
           &zrhihihi,&zrlohihi,&zrhilohi,&zrlolohi,
           &zrhihilo,&zrlohilo,&zrhilolo,&zrlololo);       // zr = xr*yr
   odg_mul(xihihihi,xilohihi,xihilohi,xilolohi,
           xihihilo,xilohilo,xihilolo,xilololo,
           yihihihi,yilohihi,yihilohi,yilolohi,
           yihihilo,yilohilo,yihilolo,yilololo,
           &acchihihi,&acclohihi,&acchilohi,&acclolohi,
           &acchihilo,&acclohilo,&acchilolo,&acclololo);   // acc = xi*yi
   odg_minus(&acchihihi,&acclohihi,&acchilohi,&acclolohi,
             &acchihilo,&acclohilo,&acchilolo,&acclololo);
   odg_inc(&zrhihihi,&zrlohihi,&zrhilohi,&zrlolohi,
           &zrhihilo,&zrlohilo,&zrhilolo,&zrlololo,
           acchihihi,acclohihi,acchilohi,acclolohi,
           acchihilo,acclohilo,acchilolo,acclololo);  // zr = xr*yr - xi*yi
   odg_mul(xrhihihi,xrlohihi,xrhilohi,xrlolohi,
           xrhihilo,xrlohilo,xrhilolo,xrlololo,
           yihihihi,yilohihi,yihilohi,yilolohi,
           yihihilo,yilohilo,yihilolo,yilololo,
           &zihihihi,&zilohihi,&zihilohi,&zilolohi,
           &zihihilo,&zilohilo,&zihilolo,&zilololo);       // zi = xr*yi
   odg_mul(xihihihi,xilohihi,xihilohi,xilolohi,
           xihihilo,xilohilo,xihilolo,xilololo,
           yrhihihi,yrlohihi,yrhilohi,yrlolohi,
           yrhihilo,yrlohilo,yrhilolo,yrlololo,
           &acchihihi,&acclohihi,&acchilohi,&acclolohi,
           &acchihilo,&acclohilo,&acchilolo,&acclololo);   // acc = xi*yr
   odg_inc(&zihihihi,&zilohihi,&zihilohi,&zilolohi,
           &zihihilo,&zilohilo,&zihilolo,&zilololo,
           acchihihi,acclohihi,acchilohi,acclolohi,
           acchihilo,acclohilo,acchilolo,acclololo); // zr = xr*yr + xi*yi

   zvrehihihi[k] = zrhihihi; zvrelohihi[k] = zrlohihi;
   zvrehilohi[k] = zrhilohi; zvrelolohi[k] = zrlolohi;
   zvrehihilo[k] = zrhihilo; zvrelohilo[k] = zrlohilo;
   zvrehilolo[k] = zrhilolo; zvrelololo[k] = zrlololo;
   zvimhihihi[k] = zihihihi; zvimlohihi[k] = zilohihi;
   zvimhilohi[k] = zihilohi; zvimlolohi[k] = zilolohi;
   zvimhihilo[k] = zihihilo; zvimlohilo[k] = zilohilo;
   zvimhilolo[k] = zihilolo; zvimlololo[k] = zilololo;

   for(int i=1; i<=k; i++) // z[k] = z[k] + x[i]*y[k-i]
   {
      xrhihihi = xvrehihihi[i]; xrlohihi = xvrelohihi[i];
      xrhilohi = xvrehilohi[i]; xrlolohi = xvrelolohi[i];
      xrhihilo = xvrehihilo[i]; xrlohilo = xvrelohilo[i];
      xrhilolo = xvrehilolo[i]; xrlololo = xvrelololo[i];
      xihihihi = xvimhihihi[i]; xilohihi = xvimlohihi[i];
      xihilohi = xvimhilohi[i]; xilolohi = xvimlolohi[i];
      xihihilo = xvimhihilo[i]; xilohilo = xvimlohilo[i];
      xihilolo = xvimhilolo[i]; xilololo = xvimlololo[i];
      yrhihihi = yvrehihihi[k-i]; yrlohihi = yvrelohihi[k-i];
      yrhilohi = yvrehilohi[k-i]; yrlolohi = yvrelolohi[k-i];
      yrhihilo = yvrehihilo[k-i]; yrlohilo = yvrelohilo[k-i];
      yrhilolo = yvrehilolo[k-i]; yrlololo = yvrelololo[k-i];
      yihihihi = yvimhihihi[k-i]; yilohihi = yvimlohihi[k-i];
      yihilohi = yvimhilohi[k-i]; yilolohi = yvimlolohi[k-i];
      yihihilo = yvimhihilo[k-i]; yilohilo = yvimlohilo[k-i];
      yihilolo = yvimhilolo[k-i]; yilololo = yvimlololo[k-i];

      odg_mul(xrhihihi,xrlohihi,xrhilohi,xrlolohi,
              xrhihilo,xrlohilo,xrhilolo,xrlololo,
              yrhihihi,yrlohihi,yrhilohi,yrlolohi,
              yrhihilo,yrlohilo,yrhilolo,yrlololo,
              &zrhihihi,&zrlohihi,&zrhilohi,&zrlolohi,
              &zrhihilo,&zrlohilo,&zrhilolo,&zrlololo);       // zr = xr*yr
      odg_mul(xihihihi,xilohihi,xihilohi,xilolohi,
              xihihilo,xilohilo,xihilolo,xilololo,
              yihihihi,yilohihi,yihilohi,yilolohi,
              yihihilo,yilohilo,yihilolo,yilololo,
              &acchihihi,&acclohihi,&acchilohi,&acclolohi,
              &acchihilo,&acclohilo,&acchilolo,&acclololo);   // acc = xi*yi
      odg_minus(&acchihihi,&acclohihi,&acchilohi,&acclolohi,
                &acchihilo,&acclohilo,&acchilolo,&acclololo);
      odg_inc(&zrhihihi,&zrlohihi,&zrhilohi,&zrlolohi,
              &zrhihilo,&zrlohilo,&zrhilolo,&zrlololo,
              acchihihi,acclohihi,acchilohi,acclolohi,
              acchihilo,acclohilo,acchilolo,acclololo); // zr = xr*yr - xi*yi
      odg_mul(xrhihihi,xrlohihi,xrhilohi,xrlolohi,
              xrhihilo,xrlohilo,xrhilolo,xrlololo,
              yihihihi,yilohihi,yihilohi,yilolohi,
              yihihilo,yilohilo,yihilolo,yilololo,
              &zihihihi,&zilohihi,&zihilohi,&zilolohi,
              &zihihilo,&zilohilo,&zihilolo,&zilololo);       // zi = xr*yi
      odg_mul(xihihihi,xilohihi,xihilohi,xilolohi,
              xihihilo,xilohilo,xihilolo,xilololo,
              yrhihihi,yrlohihi,yrhilohi,yrlolohi,
              yrhihilo,yrlohilo,yrhilolo,yrlololo,
              &acchihihi,&acclohihi,&acchilohi,&acclolohi,
              &acchihilo,&acclohilo,&acchilolo,&acclololo);   // acc = xi*yr
      odg_inc(&zihihihi,&zilohihi,&zihilohi,&zilolohi,
              &zihihilo,&zilohilo,&zihilolo,&zilololo,
              acchihihi,acclohihi,acchilohi,acclolohi,
              acchihilo,acclohilo,acchilolo,acclololo); // zr = xr*yi + xi*yr

      odg_inc(&zvrehihihi[k],&zvrelohihi[k],&zvrehilohi[k],&zvrelolohi[k],
              &zvrehihilo[k],&zvrelohilo[k],&zvrehilolo[k],&zvrelololo[k],
              zrhihihi,zrlohihi,zrhilohi,zrlolohi,
              zrhihilo,zrlohilo,zrhilolo,zrlololo);      // zvre[k] += zr
      odg_inc(&zvimhihihi[k],&zvimlohihi[k],&zvimhilohi[k],&zvimlolohi[k],
              &zvimhihilo[k],&zvimlohilo[k],&zvimhilolo[k],&zvimlololo[k],
              zihihihi,zilohihi,zihilohi,zilolohi,
              zihihilo,zilohilo,zihilolo,zilololo);      // zvim[k] += zi
   }

   __syncthreads();

   zrehihihi[k] = zvrehihihi[k]; zrelohihi[k] = zvrelohihi[k];
   zrehilohi[k] = zvrehilohi[k]; zrelolohi[k] = zvrelolohi[k];
   zrehihilo[k] = zvrehihilo[k]; zrelohilo[k] = zvrelohilo[k];
   zrehilolo[k] = zvrehilolo[k]; zrelololo[k] = zvrelololo[k];
   zimhihihi[k] = zvimhihihi[k]; zimlohihi[k] = zvimlohihi[k];
   zimhilohi[k] = zvimhilohi[k]; zimlolohi[k] = zvimlolohi[k];
   zimhihilo[k] = zvimhihilo[k]; zimlohilo[k] = zvimlohilo[k];
   zimhilolo[k] = zvimhilolo[k]; zimlololo[k] = zvimlololo[k];
}

__global__ void cmplx8_padded_convolute
 ( double *xrehihihi, double *xrelohihi, double *xrehilohi, double *xrelolohi,
   double *xrehihilo, double *xrelohilo, double *xrehilolo, double *xrelololo,
   double *ximhihihi, double *ximlohihi, double *ximhilohi, double *ximlolohi,
   double *ximhihilo, double *ximlohilo, double *ximhilolo, double *ximlololo,
   double *yrehihihi, double *yrelohihi, double *yrehilohi, double *yrelolohi,
   double *yrehihilo, double *yrelohilo, double *yrehilolo, double *yrelololo,
   double *yimhihihi, double *yimlohihi, double *yimhilohi, double *yimlolohi,
   double *yimhihilo, double *yimlohilo, double *yimhilolo, double *yimlololo,
   double *zrehihihi, double *zrelohihi, double *zrehilohi, double *zrelolohi,
   double *zrehihilo, double *zrelohilo, double *zrehilolo, double *zrelololo,
   double *zimhihihi, double *zimlohihi, double *zimhilohi, double *zimlolohi,
   double *zimhihilo, double *zimlohilo, double *zimhilolo, double *zimlololo,
   int dim )
{
   int k = threadIdx.x;       // thread k computes zre[k] and zim[k]

   __shared__ double xvrehihihi[od_shmemsize];
   __shared__ double xvrelohihi[od_shmemsize];
   __shared__ double xvrehilohi[od_shmemsize];
   __shared__ double xvrelolohi[od_shmemsize];
   __shared__ double xvrehihilo[od_shmemsize];
   __shared__ double xvrelohilo[od_shmemsize];
   __shared__ double xvrehilolo[od_shmemsize];
   __shared__ double xvrelololo[od_shmemsize];
   __shared__ double xvimhihihi[od_shmemsize];
   __shared__ double xvimlohihi[od_shmemsize];
   __shared__ double xvimhilohi[od_shmemsize];
   __shared__ double xvimlolohi[od_shmemsize];
   __shared__ double xvimhihilo[od_shmemsize];
   __shared__ double xvimlohilo[od_shmemsize];
   __shared__ double xvimhilolo[od_shmemsize];
   __shared__ double xvimlololo[od_shmemsize];
   __shared__ double yvrehihihi[od_shmemsize];
   __shared__ double yvrelohihi[od_shmemsize];
   __shared__ double yvrehilohi[od_shmemsize];
   __shared__ double yvrelolohi[od_shmemsize];
   __shared__ double yvrehihilo[od_shmemsize];
   __shared__ double yvrelohilo[od_shmemsize];
   __shared__ double yvrehilolo[od_shmemsize];
   __shared__ double yvrelololo[od_shmemsize];
   __shared__ double yvimhihihi[od_shmemsize];
   __shared__ double yvimlohihi[od_shmemsize];
   __shared__ double yvimhilohi[od_shmemsize];
   __shared__ double yvimlolohi[od_shmemsize];
   __shared__ double yvimhihilo[od_shmemsize];
   __shared__ double yvimlohilo[od_shmemsize];
   __shared__ double yvimhilolo[od_shmemsize];
   __shared__ double yvimlololo[od_shmemsize];
   __shared__ double zvrehihihi[od_shmemsize];
   __shared__ double zvrelohihi[od_shmemsize];
   __shared__ double zvrehilohi[od_shmemsize];
   __shared__ double zvrelolohi[od_shmemsize];
   __shared__ double zvrehihilo[od_shmemsize];
   __shared__ double zvrelohilo[od_shmemsize];
   __shared__ double zvrehilolo[od_shmemsize];
   __shared__ double zvrelololo[od_shmemsize];
   __shared__ double zvimhihihi[od_shmemsize];
   __shared__ double zvimlohihi[od_shmemsize];
   __shared__ double zvimhilohi[od_shmemsize];
   __shared__ double zvimlolohi[od_shmemsize];
   __shared__ double zvimhihilo[od_shmemsize];
   __shared__ double zvimlohilo[od_shmemsize];
   __shared__ double zvimhilolo[od_shmemsize];
   __shared__ double zvimlololo[od_shmemsize];

   double xrhihihi,xihihihi,yrhihihi,yihihihi,zrhihihi,zihihihi,acchihihi;
   double xrhihilo,xihihilo,yrhihilo,yihihilo,zrhihilo,zihihilo,acchihilo;
   double xrlohihi,xilohihi,yrlohihi,yilohihi,zrlohihi,zilohihi,acclohihi;
   double xrlohilo,xilohilo,yrlohilo,yilohilo,zrlohilo,zilohilo,acclohilo;
   double xrhilohi,xihilohi,yrhilohi,yihilohi,zrhilohi,zihilohi,acchilohi;
   double xrhilolo,xihilolo,yrhilolo,yihilolo,zrhilolo,zihilolo,acchilolo;
   double xrlolohi,xilolohi,yrlolohi,yilolohi,zrlolohi,zilolohi,acclolohi;
   double xrlololo,xilololo,yrlololo,yilololo,zrlololo,zilololo,acclololo;
   int idx = dim+k;

   xvrehihihi[k] = xrehihihi[k]; xvimhihihi[k] = ximhihihi[k];
   xvrelohihi[k] = xrelohihi[k]; xvimlohihi[k] = ximlohihi[k];
   xvrehilohi[k] = xrehilohi[k]; xvimhilohi[k] = ximhilohi[k];
   xvrelolohi[k] = xrelolohi[k]; xvimlolohi[k] = ximlolohi[k];
   xvrehihilo[k] = xrehihilo[k]; xvimhihilo[k] = ximhihilo[k];
   xvrelohilo[k] = xrelohilo[k]; xvimlohilo[k] = ximlohilo[k];
   xvrehilolo[k] = xrehilolo[k]; xvimhilolo[k] = ximhilolo[k];
   xvrelololo[k] = xrelololo[k]; xvimlololo[k] = ximlololo[k];
   yvrehihihi[k] = 0.0; yvimhihihi[k] = 0.0;
   yvrelohihi[k] = 0.0; yvimlohihi[k] = 0.0;
   yvrehilohi[k] = 0.0; yvimhilohi[k] = 0.0;
   yvrelolohi[k] = 0.0; yvimlolohi[k] = 0.0;
   yvrehihilo[k] = 0.0; yvimhihilo[k] = 0.0;
   yvrelohilo[k] = 0.0; yvimlohilo[k] = 0.0;
   yvrehilolo[k] = 0.0; yvimhilolo[k] = 0.0;
   yvrelololo[k] = 0.0; yvimlololo[k] = 0.0;
   yvrehihihi[idx] = yrehihihi[k]; yvimhihihi[idx] = yimhihihi[k];
   yvrelohihi[idx] = yrelohihi[k]; yvimlohihi[idx] = yimlohihi[k];
   yvrehilohi[idx] = yrehilohi[k]; yvimhilohi[idx] = yimhilohi[k];
   yvrelolohi[idx] = yrelolohi[k]; yvimlolohi[idx] = yimlolohi[k];
   yvrehihilo[idx] = yrehihilo[k]; yvimhihilo[idx] = yimhihilo[k];
   yvrelohilo[idx] = yrelohilo[k]; yvimlohilo[idx] = yimlohilo[k];
   yvrehilolo[idx] = yrehilolo[k]; yvimhilolo[idx] = yimhilolo[k];
   yvrelololo[idx] = yrelololo[k]; yvimlololo[idx] = yimlololo[k];

   __syncthreads();

   // z[k] = x[0]*y[k]
   xrhihihi = xvrehihihi[0]; xrlohihi = xvrelohihi[0];
   xrhilohi = xvrehilohi[0]; xrlolohi = xvrelolohi[0];
   xrhihilo = xvrehihilo[0]; xrlohilo = xvrelohilo[0];
   xrhilolo = xvrehilolo[0]; xrlololo = xvrelololo[0];
   xihihihi = xvimhihihi[0]; xilohihi = xvimlohihi[0];
   xihilohi = xvimhilohi[0]; xilolohi = xvimlolohi[0];
   xihihilo = xvimhihilo[0]; xilohilo = xvimlohilo[0];
   xihilolo = xvimhilolo[0]; xilololo = xvimlololo[0];
   yrhihihi = yvrehihihi[idx]; yrlohihi = yvrelohihi[idx];
   yrhilohi = yvrehilohi[idx]; yrlolohi = yvrelolohi[idx];
   yrhihilo = yvrehihilo[idx]; yrlohilo = yvrelohilo[idx];
   yrhilolo = yvrehilolo[idx]; yrlololo = yvrelololo[idx];
   yihihihi = yvimhihihi[idx]; yilohihi = yvimlohihi[idx];
   yihilohi = yvimhilohi[idx]; yilolohi = yvimlolohi[idx];
   yihihilo = yvimhihilo[idx]; yilohilo = yvimlohilo[idx];
   yihilolo = yvimhilolo[idx]; yilololo = yvimlololo[idx];

   odg_mul(xrhihihi,xrlohihi,xrhilohi,xrlolohi,
           xrhihilo,xrlohilo,xrhilolo,xrlololo,
           yrhihihi,yrlohihi,yrhilohi,yrlolohi,
           yrhihilo,yrlohilo,yrhilolo,yrlololo,
           &zrhihihi,&zrlohihi,&zrhilohi,&zrlolohi,
           &zrhihilo,&zrlohilo,&zrhilolo,&zrlololo);       // zr = xr*yr
   __syncthreads();
   odg_mul(xihihihi,xilohihi,xihilohi,xilolohi,
           xihihilo,xilohilo,xihilolo,xilololo,
           yihihihi,yilohihi,yihilohi,yilolohi,
           yihihilo,yilohilo,yihilolo,yilololo,
           &acchihihi,&acclohihi,&acchilohi,&acclolohi,
           &acchihilo,&acclohilo,&acchilolo,&acclololo);   // acc = xi*yi
   __syncthreads();
   odg_minus(&acchihihi,&acclohihi,&acchilohi,&acclolohi,
             &acchihilo,&acclohilo,&acchilolo,&acclololo);
   __syncthreads();
   odg_inc(&zrhihihi,&zrlohihi,&zrhilohi,&zrlolohi,
           &zrhihilo,&zrlohilo,&zrhilolo,&zrlololo,
           acchihihi,acclohihi,acchilohi,acclolohi,
           acchihilo,acclohilo,acchilolo,acclololo);  // zr = xr*yr - xi*yi
   __syncthreads();
   odg_mul(xrhihihi,xrlohihi,xrhilohi,xrlolohi,
           xrhihilo,xrlohilo,xrhilolo,xrlololo,
           yihihihi,yilohihi,yihilohi,yilolohi,
           yihihilo,yilohilo,yihilolo,yilololo,
           &zihihihi,&zilohihi,&zihilohi,&zilolohi,
           &zihihilo,&zilohilo,&zihilolo,&zilololo);       // zi = xr*yi
   __syncthreads();
   odg_mul(xihihihi,xilohihi,xihilohi,xilolohi,
           xihihilo,xilohilo,xihilolo,xilololo,
           yrhihihi,yrlohihi,yrhilohi,yrlolohi,
           yrhihilo,yrlohilo,yrhilolo,yrlololo,
           &acchihihi,&acclohihi,&acchilohi,&acclolohi,
           &acchihilo,&acclohilo,&acchilolo,&acclololo);   // acc = xi*yr
   __syncthreads();
   odg_inc(&zihihihi,&zilohihi,&zihilohi,&zilolohi,
           &zihihilo,&zilohilo,&zihilolo,&zilololo,
           acchihihi,acclohihi,acchilohi,acclolohi,
           acchihilo,acclohilo,acchilolo,acclololo); // zr = xr*yr + xi*yi
   __syncthreads();

   zvrehihihi[k] = zrhihihi; zvrelohihi[k] = zrlohihi;
   zvrehilohi[k] = zrhilohi; zvrelolohi[k] = zrlolohi;
   zvrehihilo[k] = zrhihilo; zvrelohilo[k] = zrlohilo;
   zvrehilolo[k] = zrhilolo; zvrelololo[k] = zrlololo;
   zvimhihihi[k] = zihihihi; zvimlohihi[k] = zilohihi;
   zvimhilohi[k] = zihilohi; zvimlolohi[k] = zilolohi;
   zvimhihilo[k] = zihihilo; zvimlohilo[k] = zilohilo;
   zvimhilolo[k] = zihilolo; zvimlololo[k] = zilololo;

   for(int i=1; i<dim; i++) // z[k] = z[k] + x[i]*y[k-i]
   {
      idx = dim + k - i;
      xrhihihi = xvrehihihi[i]; xrlohihi = xvrelohihi[i];
      xrhilohi = xvrehilohi[i]; xrlolohi = xvrelolohi[i];
      xrhihilo = xvrehihilo[i]; xrlohilo = xvrelohilo[i];
      xrhilolo = xvrehilolo[i]; xrlololo = xvrelololo[i];
      xihihihi = xvimhihihi[i]; xilohihi = xvimlohihi[i];
      xihilohi = xvimhilohi[i]; xilolohi = xvimlolohi[i];
      xihihilo = xvimhihilo[i]; xilohilo = xvimlohilo[i];
      xihilolo = xvimhilolo[i]; xilololo = xvimlololo[i];
      yrhihihi = yvrehihihi[idx]; yrlohihi = yvrelohihi[idx];
      yrhilohi = yvrehilohi[idx]; yrlolohi = yvrelolohi[idx];
      yrhihilo = yvrehihilo[idx]; yrlohilo = yvrelohilo[idx];
      yrhilolo = yvrehilolo[idx]; yrlololo = yvrelololo[idx];
      yihihihi = yvimhihihi[idx]; yilohihi = yvimlohihi[idx];
      yihilohi = yvimhilohi[idx]; yilolohi = yvimlolohi[idx];
      yihihilo = yvimhihilo[idx]; yilohilo = yvimlohilo[idx];
      yihilolo = yvimhilolo[idx]; yilololo = yvimlololo[idx];

      odg_mul(xrhihihi,xrlohihi,xrhilohi,xrlolohi,
              xrhihilo,xrlohilo,xrhilolo,xrlololo,
              yrhihihi,yrlohihi,yrhilohi,yrlolohi,
              yrhihilo,yrlohilo,yrhilolo,yrlololo,
              &zrhihihi,&zrlohihi,&zrhilohi,&zrlolohi,
              &zrhihilo,&zrlohilo,&zrhilolo,&zrlololo);       // zr = xr*yr
      __syncthreads();
      odg_mul(xihihihi,xilohihi,xihilohi,xilolohi,
              xihihilo,xilohilo,xihilolo,xilololo,
              yihihihi,yilohihi,yihilohi,yilolohi,
              yihihilo,yilohilo,yihilolo,yilololo,
              &acchihihi,&acclohihi,&acchilohi,&acclolohi,
              &acchihilo,&acclohilo,&acchilolo,&acclololo);   // acc = xi*yi
      __syncthreads();
      odg_minus(&acchihihi,&acclohihi,&acchilohi,&acclolohi,
                &acchihilo,&acclohilo,&acchilolo,&acclololo);
      __syncthreads();
      odg_inc(&zrhihihi,&zrlohihi,&zrhilohi,&zrlolohi,
              &zrhihilo,&zrlohilo,&zrhilolo,&zrlololo,
              acchihihi,acclohihi,acchilohi,acclolohi,
              acchihilo,acclohilo,acchilolo,acclololo); // zr = xr*yr - xi*yi
      __syncthreads();
      odg_mul(xrhihihi,xrlohihi,xrhilohi,xrlolohi,
              xrhihilo,xrlohilo,xrhilolo,xrlololo,
              yihihihi,yilohihi,yihilohi,yilolohi,
              yihihilo,yilohilo,yihilolo,yilololo,
              &zihihihi,&zilohihi,&zihilohi,&zilolohi,
              &zihihilo,&zilohilo,&zihilolo,&zilololo);       // zi = xr*yi
      __syncthreads();
      odg_mul(xihihihi,xilohihi,xihilohi,xilolohi,
              xihihilo,xilohilo,xihilolo,xilololo,
              yrhihihi,yrlohihi,yrhilohi,yrlolohi,
              yrhihilo,yrlohilo,yrhilolo,yrlololo,
              &acchihihi,&acclohihi,&acchilohi,&acclolohi,
              &acchihilo,&acclohilo,&acchilolo,&acclololo);   // acc = xi*yr
      __syncthreads();
      odg_inc(&zihihihi,&zilohihi,&zihilohi,&zilolohi,
              &zihihilo,&zilohilo,&zihilolo,&zilololo,
              acchihihi,acclohihi,acchilohi,acclolohi,
              acchihilo,acclohilo,acchilolo,acclololo); // zr = xr*yi + xi*yr
      __syncthreads();
      odg_inc(&zvrehihihi[k],&zvrelohihi[k],&zvrehilohi[k],&zvrelolohi[k],
              &zvrehihilo[k],&zvrelohilo[k],&zvrehilolo[k],&zvrelololo[k],
              zrhihihi,zrlohihi,zrhilohi,zrlolohi,
              zrhihilo,zrlohilo,zrhilolo,zrlololo);      // zvre[k] += zr
      __syncthreads();
      odg_inc(&zvimhihihi[k],&zvimlohihi[k],&zvimhilohi[k],&zvimlolohi[k],
              &zvimhihilo[k],&zvimlohilo[k],&zvimhilolo[k],&zvimlololo[k],
              zihihihi,zilohihi,zihilohi,zilolohi,
              zihihilo,zilohilo,zihilolo,zilololo);      // zvim[k] += zi
      __syncthreads();
   }
   zrehihihi[k] = zvrehihihi[k]; zrelohihi[k] = zvrelohihi[k];
   zrehilohi[k] = zvrehilohi[k]; zrelolohi[k] = zvrelolohi[k];
   zrehihilo[k] = zvrehihilo[k]; zrelohilo[k] = zvrelohilo[k];
   zrehilolo[k] = zvrehilolo[k]; zrelololo[k] = zvrelololo[k];
   zimhihihi[k] = zvimhihihi[k]; zimlohihi[k] = zvimlohihi[k];
   zimhilohi[k] = zvimhilohi[k]; zimlolohi[k] = zvimlolohi[k];
   zimhihilo[k] = zvimhihilo[k]; zimlohilo[k] = zvimlohilo[k];
   zimhilolo[k] = zvimhilolo[k]; zimlololo[k] = zvimlololo[k];
}

void GPU_dbl8_product
 ( double *xhihihi_h, double *xlohihi_h, double *xhilohi_h, double *xlolohi_h,
   double *xhihilo_h, double *xlohilo_h, double *xhilolo_h, double *xlololo_h,
   double *yhihihi_h, double *ylohihi_h, double *yhilohi_h, double *ylolohi_h,
   double *yhihilo_h, double *ylohilo_h, double *yhilolo_h, double *ylololo_h,
   double *zhihihi_h, double *zlohihi_h, double *zhilohi_h, double *zlolohi_h,
   double *zhihilo_h, double *zlohilo_h, double *zhilolo_h, double *zlololo_h,
   int deg, int freq, int BS, int padded )
{
   const int dim = deg+1;         // length of all vectors
   double* xhihihi_d;             // xhihihi_d is xhihihi_h on the device
   double* xlohihi_d;             // xlohihi_d is xlohihi_h on the device
   double* xhilohi_d;             // xhilohi_d is xhilohi_h on the device
   double* xlolohi_d;             // xlolohi_d is xlolohi_h on the device
   double* xhihilo_d;             // xhihilo_d is xhihilo_h on the device
   double* xlohilo_d;             // xlohilo_d is xlohilo_h on the device
   double* xhilolo_d;             // xhilolo_d is xhilolo_h on the device
   double* xlololo_d;             // xlololo_d is xlololo_h on the device
   double* yhihihi_d;             // yhihihi_d is yhihihi_h on the device
   double* ylohihi_d;             // ylohihi_d is ylohihi_h on the device
   double* yhilohi_d;             // yhilohi_d is yhilohi_h on the device
   double* ylolohi_d;             // ylolohi_d is ylolohi_h on the device
   double* yhihilo_d;             // yhihilo_d is yhihilo_h on the device
   double* ylohilo_d;             // ylohilo_d is ylohilo_h on the device
   double* yhilolo_d;             // yhilolo_d is yhilolo_h on the device
   double* ylololo_d;             // ylololo_d is ylololo_h on the device
   double* zhihihi_d;             // zhihihi_d is zhihihi_h on the device
   double* zlohihi_d;             // zlohihi_d is zlohihi_h on the device
   double* zhilohi_d;             // zhilohi_d is zhilohi_h on the device
   double* zlolohi_d;             // zlolohi_d is zlolohi_h on the device
   double* zhihilo_d;             // zhihilo_d is zhihilo_h on the device
   double* zlohilo_d;             // zlohilo_d is zlohilo_h on the device
   double* zhilolo_d;             // zhilolo_d is zhilolo_h on the device
   double* zlololo_d;             // zlololo_d is zlololo_h on the device
   size_t size = dim*sizeof(double); // number of bytes for each vector

   cudaMalloc((void**)&xhihihi_d,size);
   cudaMalloc((void**)&xlohihi_d,size);
   cudaMalloc((void**)&xhilohi_d,size);
   cudaMalloc((void**)&xlolohi_d,size);
   cudaMalloc((void**)&xhihilo_d,size);
   cudaMalloc((void**)&xlohilo_d,size);
   cudaMalloc((void**)&xhilolo_d,size);
   cudaMalloc((void**)&xlololo_d,size);
   cudaMalloc((void**)&yhihihi_d,size);
   cudaMalloc((void**)&ylohihi_d,size);
   cudaMalloc((void**)&yhilohi_d,size);
   cudaMalloc((void**)&ylolohi_d,size);
   cudaMalloc((void**)&yhihilo_d,size);
   cudaMalloc((void**)&ylohilo_d,size);
   cudaMalloc((void**)&yhilolo_d,size);
   cudaMalloc((void**)&ylololo_d,size);
   cudaMalloc((void**)&zhihihi_d,size);
   cudaMalloc((void**)&zlohihi_d,size);
   cudaMalloc((void**)&zhilohi_d,size);
   cudaMalloc((void**)&zlolohi_d,size);
   cudaMalloc((void**)&zhihilo_d,size);
   cudaMalloc((void**)&zlohilo_d,size);
   cudaMalloc((void**)&zhilolo_d,size);
   cudaMalloc((void**)&zlololo_d,size);
   cudaMemcpy(xhihihi_d,xhihihi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xlohihi_d,xlohihi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xhilohi_d,xhilohi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xlolohi_d,xlolohi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xhihilo_d,xhihilo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xlohilo_d,xlohilo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xhilolo_d,xhilolo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xlololo_d,xlololo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yhihihi_d,yhihihi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(ylohihi_d,ylohihi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yhilohi_d,yhilohi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(ylolohi_d,ylolohi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yhihilo_d,yhihilo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(ylohilo_d,ylohilo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yhilolo_d,yhilolo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(ylololo_d,ylololo_h,size,cudaMemcpyHostToDevice);

   if(dim == BS)
   {
      if(padded == 1)
      {
         for(int i=0; i<freq; i++)
            dbl8_padded_convolute<<<1,BS>>>
               (xhihihi_d,xlohihi_d,xhilohi_d,xlolohi_d,
                xhihilo_d,xlohilo_d,xhilolo_d,xlololo_d,
                yhihihi_d,ylohihi_d,yhilohi_d,ylolohi_d,
                yhihilo_d,ylohilo_d,yhilolo_d,ylololo_d,
                zhihihi_d,zlohihi_d,zhilohi_d,zlolohi_d,
                zhihilo_d,zlohilo_d,zhilolo_d,zlololo_d,dim);
      }
      else
      {
         for(int i=0; i<freq; i++)
            dbl8_convolute<<<1,BS>>>
               (xhihihi_d,xlohihi_d,xhilohi_d,xlolohi_d,
                xhihilo_d,xlohilo_d,xhilolo_d,xlololo_d,
                yhihihi_d,ylohihi_d,yhilohi_d,ylolohi_d,
                yhihilo_d,ylohilo_d,yhilolo_d,ylololo_d,
                zhihihi_d,zlohihi_d,zhilohi_d,zlolohi_d,
                zhihilo_d,zlohilo_d,zhilolo_d,zlololo_d,dim);
      }
   }
   cudaMemcpy(zhihihi_h,zhihihi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zlohihi_h,zlohihi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zhilohi_h,zhilohi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zlolohi_h,zlolohi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zhihilo_h,zhihilo_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zlohilo_h,zlohilo_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zhilolo_h,zhilolo_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zlololo_h,zlololo_d,size,cudaMemcpyDeviceToHost);
}

void GPU_cmplx8_product
 ( double *xrehihihi_h, double *xrelohihi_h,
   double *xrehilohi_h, double *xrelolohi_h,
   double *xrehihilo_h, double *xrelohilo_h,
   double *xrehilolo_h, double *xrelololo_h,
   double *ximhihihi_h, double *ximlohihi_h,
   double *ximhilohi_h, double *ximlolohi_h,
   double *ximhihilo_h, double *ximlohilo_h,
   double *ximhilolo_h, double *ximlololo_h,
   double *yrehihihi_h, double *yrelohihi_h,
   double *yrehilohi_h, double *yrelolohi_h,
   double *yrehihilo_h, double *yrelohilo_h,
   double *yrehilolo_h, double *yrelololo_h,
   double *yimhihihi_h, double *yimlohihi_h,
   double *yimhilohi_h, double *yimlolohi_h,
   double *yimhihilo_h, double *yimlohilo_h,
   double *yimhilolo_h, double *yimlololo_h,
   double *zrehihihi_h, double *zrelohihi_h,
   double *zrehilohi_h, double *zrelolohi_h,
   double *zrehihilo_h, double *zrelohilo_h,
   double *zrehilolo_h, double *zrelololo_h,
   double *zimhihihi_h, double *zimlohihi_h,
   double *zimhilohi_h, double *zimlolohi_h,
   double *zimhihilo_h, double *zimlohilo_h,
   double *zimhilolo_h, double *zimlololo_h,
   int deg, int freq, int BS, int mode )
{
   const int dim = deg+1;       // length of all vectors
   double* xrehihihi_d;         // xrehihihi_d is xrehihihi_h on the device
   double* xrelohihi_d;         // xrelohihi_d is xrelohihi_h on the device
   double* xrehilohi_d;         // xrehilohi_d is xrehilohi_h on the device
   double* xrelolohi_d;         // xrelolohi_d is xrelolohi_h on the device
   double* xrehihilo_d;         // xrehihilo_d is xrehihilo_h on the device
   double* xrelohilo_d;         // xrelohilo_d is xrelohilo_h on the device
   double* xrehilolo_d;         // xrehilolo_d is xrehilolo_h on the device
   double* xrelololo_d;         // xrelololo_d is xrelololo_h on the device
   double* ximhihihi_d;         // ximhihihi_d is ximhihihi_h on the device
   double* ximlohihi_d;         // ximlohihi_d is ximlohihi_h on the device
   double* ximhilohi_d;         // ximhilohi_d is ximhilohi_h on the device
   double* ximlolohi_d;         // ximlolohi_d is ximlolohi_h on the device
   double* ximhihilo_d;         // ximhihilo_d is ximhihilo_h on the device
   double* ximlohilo_d;         // ximlohilo_d is ximlohilo_h on the device
   double* ximhilolo_d;         // ximhilolo_d is ximhilolo_h on the device
   double* ximlololo_d;         // ximlololo_d is ximlololo_h on the device
   double* yrehihihi_d;         // yrehihihi_d is yrehihihi_h on the device
   double* yrelohihi_d;         // yrelohihi_d is yrelohihi_h on the device
   double* yrehilohi_d;         // yrehilohi_d is yrehilohi_h on the device
   double* yrelolohi_d;         // yrelolohi_d is yrelolohi_h on the device
   double* yrehihilo_d;         // yrehihilo_d is yrehihilo_h on the device
   double* yrelohilo_d;         // yrelohilo_d is yrelohilo_h on the device
   double* yrehilolo_d;         // yrehilolo_d is yrehilolo_h on the device
   double* yrelololo_d;         // yrelololo_d is yrelololo_h on the device
   double* yimhihihi_d;         // yimhihihi_d is yimhihihi_h on the device
   double* yimlohihi_d;         // yimlohihi_d is yimlohihi_h on the device
   double* yimhilohi_d;         // yimhilohi_d is yimhilohi_h on the device
   double* yimlolohi_d;         // yimlolohi_d is yimlolohi_h on the device
   double* yimhihilo_d;         // yimhihilo_d is yimhihilo_h on the device
   double* yimlohilo_d;         // yimlohilo_d is yimlohilo_h on the device
   double* yimhilolo_d;         // yimhilolo_d is yimhilolo_h on the device
   double* yimlololo_d;         // yimlololo_d is yimlololo_h on the device
   double* zrehihihi_d;         // zrehihihi_d is zrehihihi_h on the device
   double* zrelohihi_d;         // zrelohihi_d is zrelohihi_h on the device
   double* zrehilohi_d;         // zrehilohi_d is zrehilohi_h on the device
   double* zrelolohi_d;         // zrelolohi_d is zrelolohi_h on the device
   double* zrehihilo_d;         // zrehihilo_d is zrehihilo_h on the device
   double* zrelohilo_d;         // zrelohilo_d is zrelohilo_h on the device
   double* zrehilolo_d;         // zrehilolo_d is zrehilolo_h on the device
   double* zrelololo_d;         // zrelololo_d is zrelololo_h on the device
   double* zimhihihi_d;         // zimhihihi_d is zimhihihi_h on the device
   double* zimlohihi_d;         // zimlohihi_d is zimlohihi_h on the device
   double* zimhilohi_d;         // zimhilohi_d is zimhilohi_h on the device
   double* zimlolohi_d;         // zimlolohi_d is zimlolohi_h on the device
   double* zimhihilo_d;         // zimhihilo_d is zimhihilo_h on the device
   double* zimlohilo_d;         // zimlohilo_d is zimlohilo_h on the device
   double* zimhilolo_d;         // zimhilolo_d is zimhilolo_h on the device
   double* zimlololo_d;         // zimlololo_d is zimlololo_h on the device
   double* acchihihi_d;         // accumulates highest doubles
   double* acclohihi_d;         // accumulates second highest doubles
   double* acchilohi_d;         // accumulates third highest doubles
   double* acclolohi_d;         // accumulates fourth highest doubles
   double* acchihilo_d;         // accumulates fourth lowest doubles
   double* acclohilo_d;         // accumulates third lowest doubles
   double* acchilolo_d;         // accumulates second lowest doubles
   double* acclololo_d;         // accumulates lowest doubles
   size_t size = dim*sizeof(double); // number of bytes for each vector

   cudaMalloc((void**)&xrehihihi_d,size);
   cudaMalloc((void**)&xrelohihi_d,size);
   cudaMalloc((void**)&xrehilohi_d,size);
   cudaMalloc((void**)&xrelolohi_d,size);
   cudaMalloc((void**)&xrehihilo_d,size);
   cudaMalloc((void**)&xrelohilo_d,size);
   cudaMalloc((void**)&xrehilolo_d,size);
   cudaMalloc((void**)&xrelololo_d,size);
   cudaMalloc((void**)&ximhihihi_d,size);
   cudaMalloc((void**)&ximlohihi_d,size);
   cudaMalloc((void**)&ximhilohi_d,size);
   cudaMalloc((void**)&ximlolohi_d,size);
   cudaMalloc((void**)&ximhihilo_d,size);
   cudaMalloc((void**)&ximlohilo_d,size);
   cudaMalloc((void**)&ximhilolo_d,size);
   cudaMalloc((void**)&ximlololo_d,size);
   cudaMalloc((void**)&yrehihihi_d,size);
   cudaMalloc((void**)&yrelohihi_d,size);
   cudaMalloc((void**)&yrehilohi_d,size);
   cudaMalloc((void**)&yrelolohi_d,size);
   cudaMalloc((void**)&yrehihilo_d,size);
   cudaMalloc((void**)&yrelohilo_d,size);
   cudaMalloc((void**)&yrehilolo_d,size);
   cudaMalloc((void**)&yrelololo_d,size);
   cudaMalloc((void**)&yimhihihi_d,size);
   cudaMalloc((void**)&yimlohihi_d,size);
   cudaMalloc((void**)&yimhilohi_d,size);
   cudaMalloc((void**)&yimlolohi_d,size);
   cudaMalloc((void**)&yimhihilo_d,size);
   cudaMalloc((void**)&yimlohilo_d,size);
   cudaMalloc((void**)&yimhilolo_d,size);
   cudaMalloc((void**)&yimlololo_d,size);
   cudaMalloc((void**)&zrehihihi_d,size);
   cudaMalloc((void**)&zrelohihi_d,size);
   cudaMalloc((void**)&zrehilohi_d,size);
   cudaMalloc((void**)&zrelolohi_d,size);
   cudaMalloc((void**)&zrehihilo_d,size);
   cudaMalloc((void**)&zrelohilo_d,size);
   cudaMalloc((void**)&zrehilolo_d,size);
   cudaMalloc((void**)&zrelololo_d,size);
   cudaMalloc((void**)&zimhihihi_d,size);
   cudaMalloc((void**)&zimlohihi_d,size);
   cudaMalloc((void**)&zimhilohi_d,size);
   cudaMalloc((void**)&zimlolohi_d,size);
   cudaMalloc((void**)&zimhihilo_d,size);
   cudaMalloc((void**)&zimlohilo_d,size);
   cudaMalloc((void**)&zimhilolo_d,size);
   cudaMalloc((void**)&zimlololo_d,size);
   cudaMalloc((void**)&acchihihi_d,size);
   cudaMalloc((void**)&acclohihi_d,size);
   cudaMalloc((void**)&acchilohi_d,size);
   cudaMalloc((void**)&acclolohi_d,size);
   cudaMalloc((void**)&acchihilo_d,size);
   cudaMalloc((void**)&acclohilo_d,size);
   cudaMalloc((void**)&acchilolo_d,size);
   cudaMalloc((void**)&acclololo_d,size);
   cudaMemcpy(xrehihihi_d,xrehihihi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xrelohihi_d,xrelohihi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xrehilohi_d,xrehilohi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xrelolohi_d,xrelolohi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xrehihilo_d,xrehihilo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xrelohilo_d,xrelohilo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xrehilolo_d,xrehilolo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xrelololo_d,xrelololo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(ximhihihi_d,ximhihihi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(ximlohihi_d,ximlohihi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(ximhilohi_d,ximhilohi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(ximlolohi_d,ximlolohi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(ximhihilo_d,ximhihilo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(ximlohilo_d,ximlohilo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(ximhilolo_d,ximhilolo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(ximlololo_d,ximlololo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yrehihihi_d,yrehihihi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yrelohihi_d,yrelohihi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yrehilohi_d,yrehilohi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yrelolohi_d,yrelolohi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yrehihilo_d,yrehihilo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yrelohilo_d,yrelohilo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yrehilolo_d,yrehilolo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yrelololo_d,yrelololo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yimhihihi_d,yimhihihi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yimlohihi_d,yimlohihi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yimhilohi_d,yimhilohi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yimlolohi_d,yimlolohi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yimhihilo_d,yimhihilo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yimlohilo_d,yimlohilo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yimhilolo_d,yimhilolo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yimlololo_d,yimlololo_h,size,cudaMemcpyHostToDevice);

   if(dim == BS)
   {
      if(mode == 2)
      {
         for(int i=0; i<freq; i++)
            cmplx8_padded_convolute<<<1,BS>>>
               (xrehihihi_d,xrelohihi_d,xrehilohi_d,xrelolohi_d,
                xrehihilo_d,xrelohilo_d,xrehilolo_d,xrelololo_d,
                ximhihihi_d,ximlohihi_d,ximhilohi_d,ximlolohi_d,
                ximhihilo_d,ximlohilo_d,ximhilolo_d,ximlololo_d,
                yrehihihi_d,yrelohihi_d,yrehilohi_d,yrelolohi_d,
                yrehihilo_d,yrelohilo_d,yrehilolo_d,yrelololo_d,
                yimhihihi_d,yimlohihi_d,yimhilohi_d,yimlolohi_d,
                yimhihilo_d,yimlohilo_d,yimhilolo_d,yimlololo_d,
                zrehihihi_d,zrelohihi_d,zrehilohi_d,zrelolohi_d,
                zrehihilo_d,zrelohilo_d,zrehilolo_d,zrelololo_d,
                zimhihihi_d,zimlohihi_d,zimhilohi_d,zimlolohi_d,
                zimhihilo_d,zimlohilo_d,zimhilolo_d,zimlololo_d,dim);
      }
      else if(mode == 1)
      {
         dbl8_padded_convolute<<<1,BS>>>
            (xrehihihi_d,xrelohihi_d,xrehilohi_d,xrelolohi_d,
             xrehihilo_d,xrelohilo_d,xrehilolo_d,xrelololo_d,
             yrehihihi_d,yrelohihi_d,yrehilohi_d,yrelolohi_d,
             yrehihilo_d,yrelohilo_d,yrehilolo_d,yrelololo_d,
             zrehihihi_d,zrelohihi_d,zrehilohi_d,zrelolohi_d,
             zrehihilo_d,zrelohilo_d,zrehilolo_d,zrelololo_d,dim);
         dbl8_padded_convolute<<<1,BS>>>
            (ximhihihi_d,ximlohihi_d,ximhilohi_d,ximlolohi_d,
             ximhihilo_d,ximlohilo_d,ximhilolo_d,ximlololo_d,
             yimhihihi_d,yimlohihi_d,yimhilohi_d,yimlolohi_d,
             yimhihilo_d,yimlohilo_d,yimhilolo_d,yimlololo_d,
             acchihihi_d,acclohihi_d,acchilohi_d,acclolohi_d,
             acchihilo_d,acclohilo_d,acchilolo_d,acclololo_d,dim);
         dbl8_decrement<<<1,BS>>>
            (zrehihihi_d,zrelohihi_d,zrehilohi_d,zrelolohi_d,
             zrehihilo_d,zrelohilo_d,zrehilolo_d,zrelololo_d,
             acchihihi_d,acclohihi_d,acchilohi_d,acclolohi_d,
             acchihilo_d,acclohilo_d,acchilolo_d,acclololo_d,
             zrehihihi_d,zrelohihi_d,zrehilohi_d,zrelolohi_d,
             zrehihilo_d,zrelohilo_d,zrehilolo_d,zrelololo_d,dim);
         dbl8_padded_convolute<<<1,BS>>>
            (xrehihihi_d,xrelohihi_d,xrehilohi_d,xrelolohi_d,
             xrehihilo_d,xrelohilo_d,xrehilolo_d,xrelololo_d,
             yimhihihi_d,yimlohihi_d,yimhilohi_d,yimlolohi_d,
             yimhihilo_d,yimlohilo_d,yimhilolo_d,yimlololo_d,
             zimhihihi_d,zimlohihi_d,zimhilohi_d,zimlolohi_d,
             zimhihilo_d,zimlohilo_d,zimhilolo_d,zimlololo_d,dim);
         dbl8_padded_convolute<<<1,BS>>>
            (ximhihihi_d,ximlohihi_d,ximhilohi_d,ximlolohi_d,
             ximhihilo_d,ximlohilo_d,ximhilolo_d,ximlololo_d,
             yrehihihi_d,yrelohihi_d,yrehilohi_d,yrelolohi_d,
             yrehihilo_d,yrelohilo_d,yrehilolo_d,yrelololo_d,
             acchihihi_d,acclohihi_d,acchilohi_d,acclolohi_d,
             acchihilo_d,acclohilo_d,acchilolo_d,acclololo_d,dim);
         dbl8_increment<<<1,BS>>>
            (zimhihihi_d,zimlohihi_d,zimhilohi_d,zimlolohi_d,
             zimhihilo_d,zimlohilo_d,zimhilolo_d,zimlololo_d,
             acchihihi_d,acclohihi_d,acchilohi_d,acclolohi_d,
             acchihilo_d,acclohilo_d,acchilolo_d,acclololo_d,
             zimhihihi_d,zimlohihi_d,zimhilohi_d,zimlolohi_d,
             zimhihilo_d,zimlohilo_d,zimhilolo_d,zimlololo_d,dim);
      }
      else
      {
         for(int i=0; i<freq; i++)
            cmplx8_convolute<<<1,BS>>>
               (xrehihihi_d,xrelohihi_d,xrehilohi_d,xrelolohi_d,
                xrehihilo_d,xrelohilo_d,xrehilolo_d,xrelololo_d,
                ximhihihi_d,ximlohihi_d,ximhilohi_d,ximlolohi_d,
                ximhihilo_d,ximlohilo_d,ximhilolo_d,ximlololo_d,
                yrehihihi_d,yrelohihi_d,yrehilohi_d,yrelolohi_d,
                yrehihilo_d,yrelohilo_d,yrehilolo_d,yrelololo_d,
                yimhihihi_d,yimlohihi_d,yimhilohi_d,yimlolohi_d,
                yimhihilo_d,yimlohilo_d,yimhilolo_d,yimlololo_d,
                zrehihihi_d,zrelohihi_d,zrehilohi_d,zrelolohi_d,
                zrehihilo_d,zrelohilo_d,zrehilolo_d,zrelololo_d,
                zimhihihi_d,zimlohihi_d,zimhilohi_d,zimlolohi_d,
                zimhihilo_d,zimlohilo_d,zimhilolo_d,zimlololo_d,dim);
      }
   }
   cudaMemcpy(zrehihihi_h,zrehihihi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zrelohihi_h,zrelohihi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zrehilohi_h,zrehilohi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zrelolohi_h,zrelolohi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zrehihilo_h,zrehihilo_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zrelohilo_h,zrelohilo_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zrehilolo_h,zrehilolo_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zrelololo_h,zrelololo_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zimhihihi_h,zimhihihi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zimlohihi_h,zimlohihi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zimhilohi_h,zimhilohi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zimlolohi_h,zimlolohi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zimhihilo_h,zimhihilo_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zimlohilo_h,zimlohilo_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zimhilolo_h,zimhilolo_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zimlololo_h,zimlololo_d,size,cudaMemcpyDeviceToHost);
}
