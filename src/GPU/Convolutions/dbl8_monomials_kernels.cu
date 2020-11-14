// The file dbl8_monomials_kernels.cu defines kernels specified
// in dbl8_monomials_kernels.h.

/* The algorithm to compute forward, backward, and cross products
 * (denoted respectively by arrays f, b, and c)
 * for a monomial cff*x[0]*x[1]* .. *x[n-1] goes as follows:
 *
 * f[0] := cff*x[0]
 * for i from 1 to n-1 do f[i] := f[i-1]*x[i]
 * if n > 2 then
 *    b[0] := x[n-1]*x[n-2]
 *    for i from 1 to n-3 do b[i] := b[i-1]*x[n-2-i]
 *    b[n-3] := b[n-3]*cff
 *    if n = 3 then
 *       c[0] = f[0]*x[2]
 *    else
 *       for i from 0 to n-4 do c[i] := f[i]*b[n-4-i]
 *       c[n-3] := f[n-3]*x[n-1]
 *
 * Compared to the evaluation and differentiation of a product of variables,
 * (without coefficient cff), two extra multiplications must be done,
 * but this is better than n+1 multiplications with cff afterwards. */

#include "double_double_gpufun.cu"
#include "quad_double_gpufun.cu"
#include "octo_double_gpufun.cu"
#include "dbl8_convolutions_kernels.h"
#include "dbl8_monomials_kernels.h"

__device__ void dbl8_convolute
 ( double *xhihihi, double *xlohihi, double *xhilohi, double *xlolohi,
   double *xhihilo, double *xlohilo, double *xhilolo, double *xlololo,
   double *yhihihi, double *ylohihi, double *yhilohi, double *ylolohi,
   double *yhihilo, double *ylohilo, double *yhilolo, double *ylololo,
   double *zhihihi, double *zlohihi, double *zhilohi, double *zlolohi,
   double *zhihilo, double *zlohilo, double *zhilolo, double *zlololo,
   int dim, int k )
{
   double prdhihihi,prdlohihi,prdhilohi,prdlolohi;
   double prdhihilo,prdlohilo,prdhilolo,prdlololo;
   int idx;

   // z[k] = x[0]*y[k];
   odg_mul(xhihihi[0],xlohihi[0],xhilohi[0],xlolohi[0],
           xhihilo[0],xlohilo[0],xhilolo[0],xlololo[0],
           yhihihi[k],ylohihi[k],yhilohi[k],ylolohi[k],
           yhihilo[k],ylohilo[k],yhilolo[k],ylololo[k],
           &zhihihi[k],&zlohihi[k],&zhilohi[k],&zlolohi[k],
           &zhihilo[k],&zlohilo[k],&zhilolo[k],&zlololo[k]);

   for(int i=1; i<=k; i++) // z[k] = z[k] + x[i]*y[k-i];
   {
      idx = k-i;
      odg_mul(xhihihi[i],xlohihi[i],xhilohi[i],xlolohi[i],
              xhihilo[i],xlohilo[i],xhilolo[i],xlololo[i],
              yhihihi[idx],ylohihi[idx],yhilohi[idx],ylolohi[idx],
              yhihilo[idx],ylohilo[idx],yhilolo[idx],ylololo[idx],
              &prdhihihi,&prdlohihi,&prdhilohi,&prdlolohi,
              &prdhihilo,&prdlohilo,&prdhilolo,&prdlololo);
      odg_inc(&zhihihi[k],&zlohihi[k],&zhilohi[k],&zlolohi[k],
              &zhihilo[k],&zlohilo[k],&zhilolo[k],&zlololo[k],
              prdhihihi,prdlohihi,prdhilohi,prdlolohi,
              prdhihilo,prdlohilo,prdhilolo,prdlololo);
   }
}

__device__ void cmplx8_convolute
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
   int dim, int k )
{
   double xrhihihi,xihihihi,yrhihihi,yihihihi,zrhihihi,zihihihi,acchihihi;
   double xrhihilo,xihihilo,yrhihilo,yihihilo,zrhihilo,zihihilo,acchihilo;
   double xrlohihi,xilohihi,yrlohihi,yilohihi,zrlohihi,zilohihi,acclohihi;
   double xrlohilo,xilohilo,yrlohilo,yilohilo,zrlohilo,zilohilo,acclohilo;
   double xrhilohi,xihilohi,yrhilohi,yihilohi,zrhilohi,zihilohi,acchilohi;
   double xrhilolo,xihilolo,yrhilolo,yihilolo,zrhilolo,zihilolo,acchilolo;
   double xrlolohi,xilolohi,yrlolohi,yilolohi,zrlolohi,zilolohi,acclolohi;
   double xrlololo,xilololo,yrlololo,yilololo,zrlololo,zilololo,acclololo;
   int idx;

   // z[k] = x[0]*y[k]
   xrhihihi = xrehihihi[0]; xrlohihi = xrelohihi[0];
   xrhilohi = xrehilohi[0]; xrlolohi = xrelolohi[0];
   xrhihilo = xrehihilo[0]; xrlohilo = xrelohilo[0];
   xrhilolo = xrehilolo[0]; xrlololo = xrelololo[0];
   xihihihi = ximhihihi[0]; xilohihi = ximlohihi[0];
   xihilohi = ximhilohi[0]; xilolohi = ximlolohi[0];
   xihihilo = ximhihilo[0]; xilohilo = ximlohilo[0];
   xihilolo = ximhilolo[0]; xilololo = ximlololo[0];
   yrhihihi = yrehihihi[k]; yrlohihi = yrelohihi[k];
   yrhilohi = yrehilohi[k]; yrlolohi = yrelolohi[k];
   yrhihilo = yrehihilo[k]; yrlohilo = yrelohilo[k];
   yrhilolo = yrehilolo[k]; yrlololo = yrelololo[k];
   yihihihi = yimhihihi[k]; yilohihi = yimlohihi[k];
   yihilohi = yimhilohi[k]; yilolohi = yimlolohi[k];
   yihihilo = yimhihilo[k]; yilohilo = yimlohilo[k];
   yihilolo = yimhilolo[k]; yilololo = yimlololo[k];

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

   zrehihihi[k] = zrhihihi; zrelohihi[k] = zrlohihi;
   zrehilohi[k] = zrhilohi; zrelolohi[k] = zrlolohi;
   zrehihilo[k] = zrhihilo; zrelohilo[k] = zrlohilo;
   zrehilolo[k] = zrhilolo; zrelololo[k] = zrlololo;
   zimhihihi[k] = zihihihi; zimlohihi[k] = zilohihi;
   zimhilohi[k] = zihilohi; zimlolohi[k] = zilolohi;
   zimhihilo[k] = zihihilo; zimlohilo[k] = zilohilo;
   zimhilolo[k] = zihilolo; zimlololo[k] = zilololo;

   for(int i=1; i<=k; i++) // z[k] = z[k] + x[i]*y[k-i]
   {
      idx = k-i;
      xrhihihi = xrehihihi[i]; xrlohihi = xrelohihi[i];
      xrhilohi = xrehilohi[i]; xrlolohi = xrelolohi[i];
      xrhihilo = xrehihilo[i]; xrlohilo = xrelohilo[i];
      xrhilolo = xrehilolo[i]; xrlololo = xrelololo[i];
      xihihihi = ximhihihi[i]; xilohihi = ximlohihi[i];
      xihilohi = ximhilohi[i]; xilolohi = ximlolohi[i];
      xihihilo = ximhihilo[i]; xilohilo = ximlohilo[i];
      xihilolo = ximhilolo[i]; xilololo = ximlololo[i];
      yrhihihi = yrehihihi[idx]; yrlohihi = yrelohihi[idx];
      yrhilohi = yrehilohi[idx]; yrlolohi = yrelolohi[idx];
      yrhihilo = yrehihilo[idx]; yrlohilo = yrelohilo[idx];
      yrhilolo = yrehilolo[idx]; yrlololo = yrelololo[idx];
      yihihihi = yimhihihi[idx]; yilohihi = yimlohihi[idx];
      yihilohi = yimhilohi[idx]; yilolohi = yimlolohi[idx];
      yihihilo = yimhihilo[idx]; yilohilo = yimlohilo[idx];
      yihilolo = yimhilolo[idx]; yilololo = yimlololo[idx];

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

      odg_inc(&zrehihihi[k],&zrelohihi[k],&zrehilohi[k],&zrelolohi[k],
              &zrehihilo[k],&zrelohilo[k],&zrehilolo[k],&zrelololo[k],
              zrhihihi,zrlohihi,zrhilohi,zrlolohi,
              zrhihilo,zrlohilo,zrhilolo,zrlololo);      // zre[k] += zr
      odg_inc(&zimhihihi[k],&zimlohihi[k],&zimhilohi[k],&zimlolohi[k],
              &zimhihilo[k],&zimlohilo[k],&zimhilolo[k],&zimlololo[k],
              zihihihi,zilohihi,zihilohi,zilolohi,
              zihihilo,zilohilo,zihilolo,zilololo);      // zim[k] += zi
   }
}

__global__ void GPU_dbl8_speel
 ( int nvr, int deg, int *idx,
   double *cffhihihi, double *cfflohihi,
   double *cffhilohi, double *cfflolohi,
   double *cffhihilo, double *cfflohilo,
   double *cffhilolo, double *cfflololo,
   double *inputhihihi, double *inputlohihi,
   double *inputhilohi, double *inputlolohi,
   double *inputhihilo, double *inputlohilo,
   double *inputhilolo, double *inputlololo,
   double *forwardhihihi, double *forwardlohihi,
   double *forwardhilohi, double *forwardlolohi,
   double *forwardhihilo, double *forwardlohilo,
   double *forwardhilolo, double *forwardlololo,
   double *backwardhihihi, double *backwardlohihi,
   double *backwardhilohi, double *backwardlolohi,
   double *backwardhihilo, double *backwardlohilo,
   double *backwardhilolo, double *backwardlololo,
   double *crosshihihi, double *crosslohihi,
   double *crosshilohi, double *crosslolohi,
   double *crosshihilo, double *crosslohilo,
   double *crosshilolo, double *crosslololo )
{
   const int k = threadIdx.x;
   const int deg1 = deg+1;
   int ix1,ix2;

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
  
   xvhihihi[k] = cffhihihi[k]; xvlohihi[k] = cfflohihi[k];
   xvhilohi[k] = cffhilohi[k]; xvlolohi[k] = cfflolohi[k];
   xvhihilo[k] = cffhihilo[k]; xvlohilo[k] = cfflohilo[k];
   xvhilolo[k] = cffhilolo[k]; xvlololo[k] = cfflololo[k];
   ix1 = idx[0]*deg1+k;
   yvhihihi[k] = inputhihihi[ix1]; yvlohihi[k] = inputlohihi[ix1];
   yvhilohi[k] = inputhilohi[ix1]; yvlolohi[k] = inputlolohi[ix1]; 
   yvhihilo[k] = inputhihilo[ix1]; yvlohilo[k] = inputlohilo[ix1];
   yvhilolo[k] = inputhilolo[ix1]; yvlololo[k] = inputlololo[ix1]; 
   __syncthreads();                                       // f[0] = cff*x[0]
   dbl8_convolute(xvhihihi,xvlohihi,xvhilohi,xvlolohi,
                  xvhihilo,xvlohilo,xvhilolo,xvlololo,
                  yvhihihi,yvlohihi,yvhilohi,yvlolohi,
                  yvhihilo,yvlohilo,yvhilolo,yvlololo,
                  zvhihihi,zvlohihi,zvhilohi,zvlolohi,
                  zvhihilo,zvlohilo,zvhilolo,zvlololo,deg1,k);
   __syncthreads();
   forwardhihihi[k] = zvhihihi[k]; forwardlohihi[k] = zvlohihi[k];
   forwardhilohi[k] = zvhilohi[k]; forwardlolohi[k] = zvlolohi[k];
   forwardhihilo[k] = zvhihilo[k]; forwardlohilo[k] = zvlohilo[k];
   forwardhilolo[k] = zvhilolo[k]; forwardlololo[k] = zvlololo[k];

   for(int i=1; i<nvr; i++)
   {
      xvhihihi[k] = zvhihihi[k]; xvlohihi[k] = zvlohihi[k];
      xvhilohi[k] = zvhilohi[k]; xvlolohi[k] = zvlolohi[k];
      xvhihilo[k] = zvhihilo[k]; xvlohilo[k] = zvlohilo[k];
      xvhilolo[k] = zvhilolo[k]; xvlololo[k] = zvlololo[k];
      ix2 = idx[i]*deg1+k;
      yvhihihi[k] = inputhihihi[ix2]; yvlohihi[k] = inputlohihi[ix2];
      yvhilohi[k] = inputhilohi[ix2]; yvlolohi[k] = inputlolohi[ix2];
      yvhihilo[k] = inputhihilo[ix2]; yvlohilo[k] = inputlohilo[ix2];
      yvhilolo[k] = inputhilolo[ix2]; yvlololo[k] = inputlololo[ix2];
      __syncthreads();                                 // f[i] = f[i-1]*x[i]
      dbl8_convolute(xvhihihi,xvlohihi,xvhilohi,xvlolohi,
                     xvhihilo,xvlohilo,xvhilolo,xvlololo,
                     yvhihihi,yvlohihi,yvhilohi,yvlolohi,
                     yvhihilo,yvlohilo,yvhilolo,yvlololo,
                     zvhihihi,zvlohihi,zvhilohi,zvlolohi,
                     zvhihilo,zvlohilo,zvhilolo,zvlololo,deg1,k);
      __syncthreads();
      ix1 = i*deg1+k;
      forwardhihihi[ix1] = zvhihihi[k]; forwardlohihi[ix1] = zvlohihi[k];
      forwardhilohi[ix1] = zvhilohi[k]; forwardlolohi[ix1] = zvlolohi[k];
      forwardhihilo[ix1] = zvhihilo[k]; forwardlohilo[ix1] = zvlohilo[k];
      forwardhilolo[ix1] = zvhilolo[k]; forwardlololo[ix1] = zvlololo[k];
   }
   if(nvr > 2)
   {
      ix1 = idx[nvr-1]*deg1+k;
      xvhihihi[k] = inputhihihi[ix1]; xvlohihi[k] = inputlohihi[ix1]; 
      xvhilohi[k] = inputhilohi[ix1]; xvlolohi[k] = inputlolohi[ix1];
      xvhihilo[k] = inputhihilo[ix1]; xvlohilo[k] = inputlohilo[ix1]; 
      xvhilolo[k] = inputhilolo[ix1]; xvlololo[k] = inputlololo[ix1];
      ix2 = idx[nvr-2]*deg1+k;
      yvhihihi[k] = inputhihihi[ix2]; yvlohihi[k] = inputlohihi[ix2];
      yvhilohi[k] = inputhilohi[ix2]; yvlolohi[k] = inputlolohi[ix2];
      yvhihilo[k] = inputhihilo[ix2]; yvlohilo[k] = inputlohilo[ix2];
      yvhilolo[k] = inputhilolo[ix2]; yvlololo[k] = inputlololo[ix2];
      __syncthreads();                               // b[0] = x[n-1]*x[n-2]
      dbl8_convolute(xvhihihi,xvlohihi,xvhilohi,xvlolohi,
                     xvhihilo,xvlohilo,xvhilolo,xvlololo,
                     yvhihihi,yvlohihi,yvhilohi,yvlolohi,
                     yvhihilo,yvlohilo,yvhilolo,yvlololo,
                     zvhihihi,zvlohihi,zvhilohi,zvlolohi,
                     zvhihilo,zvlohilo,zvhilolo,zvlololo,deg1,k);
      __syncthreads();
      backwardhihihi[k] = zvhihihi[k]; backwardlohihi[k] = zvlohihi[k];
      backwardhilohi[k] = zvhilohi[k]; backwardlolohi[k] = zvlolohi[k];
      backwardhihilo[k] = zvhihilo[k]; backwardlohilo[k] = zvlohilo[k];
      backwardhilolo[k] = zvhilolo[k]; backwardlololo[k] = zvlololo[k];
      for(int i=1; i<nvr-2; i++)
      {
         xvhihihi[k] = zvhihihi[k]; xvlohihi[k] = zvlohihi[k];
         xvhilohi[k] = zvhilohi[k]; xvlolohi[k] = zvlolohi[k];
         xvhihilo[k] = zvhihilo[k]; xvlohilo[k] = zvlohilo[k];
         xvhilolo[k] = zvhilolo[k]; xvlololo[k] = zvlololo[k];
         ix2 = idx[nvr-2-i]*deg1+k;
         yvhihihi[k] = inputhihihi[ix2]; yvlohihi[k] = inputlohihi[ix2];
         yvhilohi[k] = inputhilohi[ix2]; yvlolohi[k] = inputlolohi[ix2];
         yvhihilo[k] = inputhihilo[ix2]; yvlohilo[k] = inputlohilo[ix2];
         yvhilolo[k] = inputhilolo[ix2]; yvlololo[k] = inputlololo[ix2];
         __syncthreads();                          // b[i] = b[i-1]*x[n-2-i]
         dbl8_convolute(xvhihihi,xvlohihi,xvhilohi,xvlolohi,
                        xvhihilo,xvlohilo,xvhilolo,xvlololo,
                        yvhihihi,yvlohihi,yvhilohi,yvlolohi,
                        yvhihilo,yvlohilo,yvhilolo,yvlololo,
                        zvhihihi,zvlohihi,zvhilohi,zvlolohi,
                        zvhihilo,zvlohilo,zvhilolo,zvlololo,deg1,k);
         __syncthreads();
         ix1 = i*deg1+k;
         backwardhihihi[ix1] = zvhihihi[k]; backwardlohihi[ix1] = zvlohihi[k];
         backwardhilohi[ix1] = zvhilohi[k]; backwardlolohi[ix1] = zvlolohi[k];
         backwardhihilo[ix1] = zvhihilo[k]; backwardlohilo[ix1] = zvlohilo[k];
         backwardhilolo[ix1] = zvhilolo[k]; backwardlololo[ix1] = zvlololo[k];
      }
      xvhihihi[k] = zvhihihi[k]; xvlohihi[k] = zvlohihi[k];
      xvhilohi[k] = zvhilohi[k]; xvlolohi[k] = zvlolohi[k];
      xvhihilo[k] = zvhihilo[k]; xvlohilo[k] = zvlohilo[k];
      xvhilolo[k] = zvhilolo[k]; xvlololo[k] = zvlololo[k];
      yvhihihi[k] = cffhihihi[k]; yvlohihi[k] = cfflohihi[k];
      yvhilohi[k] = cffhilohi[k]; yvlolohi[k] = cfflolohi[k];
      yvhihilo[k] = cffhihilo[k]; yvlohilo[k] = cfflohilo[k];
      yvhilolo[k] = cffhilolo[k]; yvlololo[k] = cfflololo[k];
      __syncthreads();                                // b[n-3] = b[n-3]*cff
      dbl8_convolute(xvhihihi,xvlohihi,xvhilohi,xvlolohi,
                     xvhihilo,xvlohilo,xvhilolo,xvlololo,
                     yvhihihi,yvlohihi,yvhilohi,yvlolohi,
                     yvhihilo,yvlohilo,yvhilolo,yvlololo,
                     zvhihihi,zvlohihi,zvhilohi,zvlolohi,
                     zvhihilo,zvlohilo,zvhilolo,zvlololo,deg1,k);
      __syncthreads();
      ix2 = (nvr-3)*deg1+k;
      backwardhihihi[ix2] = zvhihihi[k]; backwardlohihi[ix2] = zvlohihi[k];
      backwardhilohi[ix2] = zvhilohi[k]; backwardlolohi[ix2] = zvlolohi[k];
      backwardhihilo[ix2] = zvhihilo[k]; backwardlohilo[ix2] = zvlohilo[k];
      backwardhilolo[ix2] = zvhilolo[k]; backwardlololo[ix2] = zvlololo[k];

      if(nvr == 3)
      {
         xvhihihi[k] = forwardhihihi[k]; xvlohihi[k] = forwardlohihi[k];
         xvhilohi[k] = forwardhilohi[k]; xvlolohi[k] = forwardlolohi[k];
         xvhihilo[k] = forwardhihilo[k]; xvlohilo[k] = forwardlohilo[k];
         xvhilolo[k] = forwardhilolo[k]; xvlololo[k] = forwardlololo[k];
         ix2 = idx[2]*deg1+k;
         yvhihihi[k] = inputhihihi[ix2]; yvlohihi[k] = inputlohihi[ix2];
         yvhilohi[k] = inputhilohi[ix2]; yvlolohi[k] = inputlolohi[ix2];
         yvhihilo[k] = inputhihilo[ix2]; yvlohilo[k] = inputlohilo[ix2];
         yvhilolo[k] = inputhilolo[ix2]; yvlololo[k] = inputlololo[ix2];
         __syncthreads();                                // c[0] = f[0]*x[2]
         dbl8_convolute(xvhihihi,xvlohihi,xvhilohi,xvlolohi,
                        xvhihilo,xvlohilo,xvhilolo,xvlololo,
                        yvhihihi,yvlohihi,yvhilohi,yvlolohi,
                        yvhihilo,yvlohilo,yvhilolo,yvlololo,
                        zvhihihi,zvlohihi,zvhilohi,zvlolohi,
                        zvhihilo,zvlohilo,zvhilolo,zvlololo,deg1,k);
         __syncthreads();
         crosshihihi[k] = zvhihihi[k]; crosslohihi[k] = zvlohihi[k];
         crosshilohi[k] = zvhilohi[k]; crosslolohi[k] = zvlolohi[k];
         crosshihilo[k] = zvhihilo[k]; crosslohilo[k] = zvlohilo[k];
         crosshilolo[k] = zvhilolo[k]; crosslololo[k] = zvlololo[k];
      }
      else
      {
         for(int i=0; i<nvr-3; i++)
         {
            ix1 = i*deg1+k; 
            xvhihihi[k] = forwardhihihi[ix1];
            xvlohihi[k] = forwardlohihi[ix1];
            xvhilohi[k] = forwardhilohi[ix1];
            xvlolohi[k] = forwardlolohi[ix1];
            xvhihilo[k] = forwardhihilo[ix1];
            xvlohilo[k] = forwardlohilo[ix1];
            xvhilolo[k] = forwardhilolo[ix1];
            xvlololo[k] = forwardlololo[ix1];
            ix2 = (nvr-4-i)*deg1+k;
            yvhihihi[k] = backwardhihihi[ix2];
            yvlohihi[k] = backwardlohihi[ix2];
            yvhilohi[k] = backwardhilohi[ix2];
            yvlolohi[k] = backwardlolohi[ix2];
            yvhihilo[k] = backwardhihilo[ix2];
            yvlohilo[k] = backwardlohilo[ix2];
            yvhilolo[k] = backwardhilolo[ix2];
            yvlololo[k] = backwardlololo[ix2];
            __syncthreads();                         // c[i] = f[i]*b[n-4-i]
            dbl8_convolute(xvhihihi,xvlohihi,xvhilohi,xvlolohi,
                           xvhihilo,xvlohilo,xvhilolo,xvlololo,
                           yvhihihi,yvlohihi,yvhilohi,yvlolohi,
                           yvhihilo,yvlohilo,yvhilolo,yvlololo,
                           zvhihihi,zvlohihi,zvhilohi,zvlolohi,
                           zvhihilo,zvlohilo,zvhilolo,zvlololo,deg1,k);
            __syncthreads();
            crosshihihi[ix1] = zvhihihi[k]; crosslohihi[ix1] = zvlohihi[k];
            crosshilohi[ix1] = zvhilohi[k]; crosslolohi[ix1] = zvlolohi[k]; 
            crosshihilo[ix1] = zvhihilo[k]; crosslohilo[ix1] = zvlohilo[k];
            crosshilolo[ix1] = zvhilolo[k]; crosslololo[ix1] = zvlololo[k]; 
         }
         ix1 = (nvr-3)*deg1+k;
         xvhihihi[k] = forwardhihihi[ix1]; xvlohihi[k] = forwardlohihi[ix1];
         xvhilohi[k] = forwardhilohi[ix1]; xvlolohi[k] = forwardlolohi[ix1];
         xvhihilo[k] = forwardhihilo[ix1]; xvlohilo[k] = forwardlohilo[ix1];
         xvhilolo[k] = forwardhilolo[ix1]; xvlololo[k] = forwardlololo[ix1];
         ix2 = idx[nvr-1]*deg1+k;
         yvhihihi[k] = inputhihihi[ix2]; yvlohihi[k] = inputlohihi[ix2];
         yvhilohi[k] = inputhilohi[ix2]; yvlolohi[k] = inputlolohi[ix2];
         yvhihilo[k] = inputhihilo[ix2]; yvlohilo[k] = inputlohilo[ix2];
         yvhilolo[k] = inputhilolo[ix2]; yvlololo[k] = inputlololo[ix2];
         __syncthreads();                          // c[n-3] = f[n-3]*x[n-1]
         dbl8_convolute(xvhihihi,xvlohihi,xvhilohi,xvlolohi,
                        xvhihilo,xvlohilo,xvhilolo,xvlololo,
                        yvhihihi,yvlohihi,yvhilohi,yvlolohi,
                        yvhihilo,yvlohilo,yvhilolo,yvlololo,
                        zvhihihi,zvlohihi,zvhilohi,zvlolohi,
                        zvhihilo,zvlohilo,zvhilolo,zvlololo,deg1,k);
         __syncthreads();
         crosshihihi[ix1] = zvhihihi[k]; crosslohihi[ix1] = zvlohihi[k];
         crosshilohi[ix1] = zvhilohi[k]; crosslolohi[ix1] = zvlolohi[k];
         crosshihilo[ix1] = zvhihilo[k]; crosslohilo[ix1] = zvlohilo[k];
         crosshilolo[ix1] = zvhilolo[k]; crosslololo[ix1] = zvlololo[k];
      }
   }
}

__global__ void GPU_cmplx8_speel
 ( int nvr, int deg, int *idx,
   double *cffrehihihi, double *cffrelohihi,
   double *cffrehilohi, double *cffrelolohi,
   double *cffrehihilo, double *cffrelohilo,
   double *cffrehilolo, double *cffrelololo,
   double *cffimhihihi, double *cffimlohihi,
   double *cffimhilohi, double *cffimlolohi,
   double *cffimhihilo, double *cffimlohilo,
   double *cffimhilolo, double *cffimlololo,
   double *inputrehihihi, double *inputrelohihi,
   double *inputrehilohi, double *inputrelolohi,
   double *inputrehihilo, double *inputrelohilo,
   double *inputrehilolo, double *inputrelololo,
   double *inputimhihihi, double *inputimlohihi,
   double *inputimhilohi, double *inputimlolohi,
   double *inputimhihilo, double *inputimlohilo,
   double *inputimhilolo, double *inputimlololo,
   double *forwardrehihihi, double *forwardrelohihi,
   double *forwardrehilohi, double *forwardrelolohi,
   double *forwardrehihilo, double *forwardrelohilo,
   double *forwardrehilolo, double *forwardrelololo,
   double *forwardimhihihi, double *forwardimlohihi,
   double *forwardimhilohi, double *forwardimlolohi,
   double *forwardimhihilo, double *forwardimlohilo,
   double *forwardimhilolo, double *forwardimlololo,
   double *backwardrehihihi, double *backwardrelohihi,
   double *backwardrehilohi, double *backwardrelolohi,
   double *backwardrehihilo, double *backwardrelohilo,
   double *backwardrehilolo, double *backwardrelololo,
   double *backwardimhihihi, double *backwardimlohihi,
   double *backwardimhilohi, double *backwardimlolohi,
   double *backwardimhihilo, double *backwardimlohilo,
   double *backwardimhilolo, double *backwardimlololo,
   double *crossrehihihi, double *crossrelohihi,
   double *crossrehilohi, double *crossrelolohi,
   double *crossrehihilo, double *crossrelohilo,
   double *crossrehilolo, double *crossrelololo,
   double *crossimhihihi, double *crossimlohihi,
   double *crossimhilohi, double *crossimlolohi,
   double *crossimhihilo, double *crossimlohilo,
   double *crossimhilolo, double *crossimlololo )
{
}

void GPU_dbl8_evaldiff
 ( int BS, int dim, int nvr, int deg, int *idx,
   double *cffhihihi, double *cfflohihi, double *cffhilohi, double *cfflolohi,
   double *cffhihilo, double *cfflohilo, double *cffhilolo, double *cfflololo,
   double **inputhihihi, double **inputlohihi,
   double **inputhilohi, double **inputlolohi,
   double **inputhihilo, double **inputlohilo,
   double **inputhilolo, double **inputlololo,
   double **outputhihihi, double **outputlohihi,
   double **outputhilohi, double **outputlolohi,
   double **outputhihilo, double **outputlohilo,
   double **outputhilolo, double **outputlololo )
{
   const int deg1 = deg+1;   // length of all vectors
   double *inputhihihi_d;    // inputhihihi_d is inputhihihi on the device
   double *inputlohihi_d;    // inputlohihi_d is inputlohihi on the device
   double *inputhilohi_d;    // inputhilohi_d is inputhilohi on the device
   double *inputlolohi_d;    // inputlolohi_d is inputlolohi on the device
   double *inputhihilo_d;    // inputhihilo_d is inputhihilo on the device
   double *inputlohilo_d;    // inputlohilo_d is inputlohilo on the device
   double *inputhilolo_d;    // inputhilolo_d is inputhilolo on the device
   double *inputlololo_d;    // inputlololo_d is inputlololo on the device
   double *forwardhihihi_d;  // highest forward products on the device
   double *forwardlohihi_d;  // second highest forward products on the device
   double *forwardhilohi_d;  // third highest forward products on the device
   double *forwardlolohi_d;  // fourth highest forward products on the device
   double *forwardhihilo_d;  // fourth lowest forward products on the device
   double *forwardlohilo_d;  // third lowest forward products on the device
   double *forwardhilolo_d;  // second lowest forward products on the device
   double *forwardlololo_d;  // lowest forward products on the device
   double *backwardhihihi_d; // highest backward products on the device
   double *backwardlohihi_d; // second highest backward products on the device
   double *backwardhilohi_d; // third highest backward products on the device
   double *backwardlolohi_d; // fourth highest backward products on the device
   double *backwardhihilo_d; // fourth lowest backward products on the device
   double *backwardlohilo_d; // third lowest backward products on the device
   double *backwardhilolo_d; // second lowest backward products on the device
   double *backwardlololo_d; // lowest backward products on the device
   double *crosshihihi_d;    // highest cross products on the device
   double *crosslohihi_d;    // second highest cross products on the device
   double *crosshilohi_d;    // third highest cross products on the device
   double *crosslolohi_d;    // fourth highest cross products on the device
   double *crosshihilo_d;    // fourth lowest cross products on the device
   double *crosslohilo_d;    // third lowest cross products on the device
   double *crosshilolo_d;    // second lowest cross products on the device
   double *crosslololo_d;    // lowest cross products on the device
   double *cffhihihi_d;      // cffhihihi_d is cffhihihi on device
   double *cfflohihi_d;      // cfflohihi_d is cfflohihi on device
   double *cffhilohi_d;      // cffhilohi_d is cffhilohi on device
   double *cfflolohi_d;      // cfflolohi_d is cfflolohi on device
   double *cffhihilo_d;      // cffhihilo_d is cffhihilo on device
   double *cfflohilo_d;      // cfflohilo_d is cfflohilo on device
   double *cffhilolo_d;      // cffhilolo_d is cffhilolo on device
   double *cfflololo_d;      // cfflololo_d is cfflololo on device
   int *idx_d;               // idx_d is idx on device

   size_t szcff = deg1*sizeof(double);
   size_t szdim = dim*(deg1)*sizeof(double);
   size_t sznvr = nvr*(deg1)*sizeof(double);
   size_t sznvr2 = (nvr-2)*(deg1)*sizeof(double);
   size_t szidx = nvr*sizeof(int);

   cudaMalloc((void**)&idx_d,szidx);
   cudaMalloc((void**)&cffhihihi_d,szcff);
   cudaMalloc((void**)&cfflohihi_d,szcff);
   cudaMalloc((void**)&cffhilohi_d,szcff);
   cudaMalloc((void**)&cfflolohi_d,szcff);
   cudaMalloc((void**)&cffhihilo_d,szcff);
   cudaMalloc((void**)&cfflohilo_d,szcff);
   cudaMalloc((void**)&cffhilolo_d,szcff);
   cudaMalloc((void**)&cfflololo_d,szcff);
   cudaMalloc((void**)&inputhihihi_d,szdim);
   cudaMalloc((void**)&inputlohihi_d,szdim);
   cudaMalloc((void**)&inputhilohi_d,szdim);
   cudaMalloc((void**)&inputlolohi_d,szdim);
   cudaMalloc((void**)&inputhihilo_d,szdim);
   cudaMalloc((void**)&inputlohilo_d,szdim);
   cudaMalloc((void**)&inputhilolo_d,szdim);
   cudaMalloc((void**)&inputlololo_d,szdim);
   cudaMalloc((void**)&forwardhihihi_d,sznvr);
   cudaMalloc((void**)&forwardlohihi_d,sznvr);
   cudaMalloc((void**)&forwardhilohi_d,sznvr);
   cudaMalloc((void**)&forwardlolohi_d,sznvr);
   cudaMalloc((void**)&forwardhihilo_d,sznvr);
   cudaMalloc((void**)&forwardlohilo_d,sznvr);
   cudaMalloc((void**)&forwardhilolo_d,sznvr);
   cudaMalloc((void**)&forwardlololo_d,sznvr);
   cudaMalloc((void**)&backwardhihihi_d,sznvr2);
   cudaMalloc((void**)&backwardlohihi_d,sznvr2);
   cudaMalloc((void**)&backwardhilohi_d,sznvr2);
   cudaMalloc((void**)&backwardlolohi_d,sznvr2);
   cudaMalloc((void**)&backwardhihilo_d,sznvr2);
   cudaMalloc((void**)&backwardlohilo_d,sznvr2);
   cudaMalloc((void**)&backwardhilolo_d,sznvr2);
   cudaMalloc((void**)&backwardlololo_d,sznvr2);
   cudaMalloc((void**)&crosshihihi_d,sznvr2);
   cudaMalloc((void**)&crosslohihi_d,sznvr2);
   cudaMalloc((void**)&crosshilohi_d,sznvr2);
   cudaMalloc((void**)&crosslolohi_d,sznvr2);
   cudaMalloc((void**)&crosshihilo_d,sznvr2);
   cudaMalloc((void**)&crosslohilo_d,sznvr2);
   cudaMalloc((void**)&crosshilolo_d,sznvr2);
   cudaMalloc((void**)&crosslololo_d,sznvr2);

   double *inputhihihi_h = new double[dim*(deg1)];
   double *inputlohihi_h = new double[dim*(deg1)];
   double *inputhilohi_h = new double[dim*(deg1)];
   double *inputlolohi_h = new double[dim*(deg1)];
   double *inputhihilo_h = new double[dim*(deg1)];
   double *inputlohilo_h = new double[dim*(deg1)];
   double *inputhilolo_h = new double[dim*(deg1)];
   double *inputlololo_h = new double[dim*(deg1)];
   int ix = 0;
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++)
      {
         inputhihihi_h[ix] = inputhihihi[i][j];
         inputlohihi_h[ix] = inputlohihi[i][j];
         inputhilohi_h[ix] = inputhilohi[i][j];
         inputlolohi_h[ix] = inputlolohi[i][j];
         inputhihilo_h[ix] = inputhihilo[i][j];
         inputlohilo_h[ix] = inputlohilo[i][j];
         inputhilolo_h[ix] = inputhilolo[i][j];
         inputlololo_h[ix++] = inputlololo[i][j];
      }

   cudaMemcpy(idx_d,idx,szidx,cudaMemcpyHostToDevice);
   cudaMemcpy(cffhihihi_d,cffhihihi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cfflohihi_d,cfflohihi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffhilohi_d,cffhilohi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cfflolohi_d,cfflolohi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffhihilo_d,cffhihilo,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cfflohilo_d,cfflohilo,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffhilolo_d,cffhilolo,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cfflololo_d,cfflololo,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(inputhihihi_d,inputhihihi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputlohihi_d,inputlohihi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputhilohi_d,inputhilohi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputlolohi_d,inputlolohi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputhihilo_d,inputhihilo_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputlohilo_d,inputlohilo_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputhilolo_d,inputhilolo_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputlololo_d,inputlololo_h,szdim,cudaMemcpyHostToDevice);

   if(BS == deg1)
   {
      GPU_dbl8_speel<<<1,BS>>>
         (nvr,deg,idx_d,
          cffhihihi_d,cfflohihi_d,cffhilohi_d,cfflolohi_d,
          cffhihilo_d,cfflohilo_d,cffhilolo_d,cfflololo_d,
          inputhihihi_d,inputlohihi_d,inputhilohi_d,inputlolohi_d,
          inputhihilo_d,inputlohilo_d,inputhilolo_d,inputlololo_d,
          forwardhihihi_d,forwardlohihi_d,forwardhilohi_d,forwardlolohi_d,
          forwardhihilo_d,forwardlohilo_d,forwardhilolo_d,forwardlololo_d,
          backwardhihihi_d,backwardlohihi_d,backwardhilohi_d,backwardlolohi_d,
          backwardhihilo_d,backwardlohilo_d,backwardhilolo_d,backwardlololo_d,
          crosshihihi_d,crosslohihi_d,crosshilohi_d,crosslolohi_d,
          crosshihilo_d,crosslohilo_d,crosshilolo_d,crosslololo_d);
   }
   double *forwardhihihi_h = new double[(deg1)*nvr];
   double *forwardlohihi_h = new double[(deg1)*nvr];
   double *forwardhilohi_h = new double[(deg1)*nvr];
   double *forwardlolohi_h = new double[(deg1)*nvr];
   double *forwardhihilo_h = new double[(deg1)*nvr];
   double *forwardlohilo_h = new double[(deg1)*nvr];
   double *forwardhilolo_h = new double[(deg1)*nvr];
   double *forwardlololo_h = new double[(deg1)*nvr];
   double *backwardhihihi_h = new double[(deg1)*(nvr-2)];
   double *backwardlohihi_h = new double[(deg1)*(nvr-2)];
   double *backwardhilohi_h = new double[(deg1)*(nvr-2)];
   double *backwardlolohi_h = new double[(deg1)*(nvr-2)];
   double *backwardhihilo_h = new double[(deg1)*(nvr-2)];
   double *backwardlohilo_h = new double[(deg1)*(nvr-2)];
   double *backwardhilolo_h = new double[(deg1)*(nvr-2)];
   double *backwardlololo_h = new double[(deg1)*(nvr-2)];
   double *crosshihihi_h = new double[(deg1)*(nvr-2)];
   double *crosslohihi_h = new double[(deg1)*(nvr-2)];
   double *crosshilohi_h = new double[(deg1)*(nvr-2)];
   double *crosslolohi_h = new double[(deg1)*(nvr-2)];
   double *crosshihilo_h = new double[(deg1)*(nvr-2)];
   double *crosslohilo_h = new double[(deg1)*(nvr-2)];
   double *crosshilolo_h = new double[(deg1)*(nvr-2)];
   double *crosslololo_h = new double[(deg1)*(nvr-2)];
  
   cudaMemcpy(forwardhihihi_h,forwardhihihi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardlohihi_h,forwardlohihi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardhilohi_h,forwardhilohi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardlolohi_h,forwardlolohi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardhihilo_h,forwardhihilo_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardlohilo_h,forwardlohilo_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardhilolo_h,forwardhilolo_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardlololo_h,forwardlololo_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardhihihi_h,backwardhihihi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardlohihi_h,backwardlohihi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardhilohi_h,backwardhilohi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardlolohi_h,backwardlolohi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardhihilo_h,backwardhihilo_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardlohilo_h,backwardlohilo_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardhilolo_h,backwardhilolo_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardlololo_h,backwardlololo_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crosshihihi_h,crosshihihi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crosslohihi_h,crosslohihi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crosshilohi_h,crosshilohi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crosslolohi_h,crosslolohi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crosshihilo_h,crosshihilo_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crosslohilo_h,crosslohilo_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crosshilolo_h,crosshilolo_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crosslololo_h,crosslololo_d,sznvr2,cudaMemcpyDeviceToHost);

   int offset = (nvr-1)*deg1;            // assign value of the monomial
   for(int i=0; i<deg1; i++)
   {
      outputhihihi[dim][i] = forwardhihihi_h[offset+i];
      outputlohihi[dim][i] = forwardlohihi_h[offset+i];
      outputhilohi[dim][i] = forwardhilohi_h[offset+i];
      outputlolohi[dim][i] = forwardlolohi_h[offset+i];
      outputhihilo[dim][i] = forwardhihilo_h[offset+i];
      outputlohilo[dim][i] = forwardlohilo_h[offset+i];
      outputhilolo[dim][i] = forwardhilolo_h[offset+i];
      outputlololo[dim][i] = forwardlololo_h[offset+i];
   }
   ix = idx[nvr-1];                      // derivative with respect to x[n-1]
   offset = (nvr-2)*deg1;
   for(int i=0; i<deg1; i++)
   {
      outputhihihi[ix][i] = forwardhihihi_h[offset+i];
      outputlohihi[ix][i] = forwardlohihi_h[offset+i];
      outputhilohi[ix][i] = forwardhilohi_h[offset+i];
      outputlolohi[ix][i] = forwardlolohi_h[offset+i];
      outputhihilo[ix][i] = forwardhihilo_h[offset+i];
      outputlohilo[ix][i] = forwardlohilo_h[offset+i];
      outputhilolo[ix][i] = forwardhilolo_h[offset+i];
      outputlololo[ix][i] = forwardlololo_h[offset+i];
   }
   ix = idx[0];                          // derivative with respect to x[0]
   offset = (nvr-3)*deg1;
   for(int i=0; i<deg1; i++)
   {
      outputhihihi[ix][i] = backwardhihihi_h[offset+i];
      outputlohihi[ix][i] = backwardlohihi_h[offset+i];
      outputhilohi[ix][i] = backwardhilohi_h[offset+i];
      outputlolohi[ix][i] = backwardlolohi_h[offset+i];
      outputhihilo[ix][i] = backwardhihilo_h[offset+i];
      outputlohilo[ix][i] = backwardlohilo_h[offset+i];
      outputhilolo[ix][i] = backwardhilolo_h[offset+i];
      outputlololo[ix][i] = backwardlololo_h[offset+i];
   }
   for(int k=1; k<nvr-1; k++)            // derivative with respect to x[k]
   {
      ix = idx[k]; offset = (k-1)*deg1;
      for(int i=0; i<deg1; i++)
      {
         outputhihihi[ix][i] = crosshihihi_h[offset+i];
         outputlohihi[ix][i] = crosslohihi_h[offset+i];
         outputhilohi[ix][i] = crosshilohi_h[offset+i];
         outputlolohi[ix][i] = crosslolohi_h[offset+i];
         outputhihilo[ix][i] = crosshihilo_h[offset+i];
         outputlohilo[ix][i] = crosslohilo_h[offset+i];
         outputhilolo[ix][i] = crosshilolo_h[offset+i];
         outputlololo[ix][i] = crosslololo_h[offset+i];
      }
   }
}

void GPU_cmplx8_evaldiff
 ( int BS, int dim, int nvr, int deg, int *idx,
   double *cffrehihihi, double *cffrelohihi,
   double *cffrehilohi, double *cffrelolohi,
   double *cffrehihilo, double *cffrelohilo,
   double *cffrehilolo, double *cffrelololo,
   double *cffimhihihi, double *cffimlohihi,
   double *cffimhilohi, double *cffimlolohi,
   double *cffimhihilo, double *cffimlohilo,
   double *cffimhilolo, double *cffimlololo,
   double **inputrehihihi, double **inputrelohihi,
   double **inputrehilohi, double **inputrelolohi,
   double **inputrehihilo, double **inputrelohilo,
   double **inputrehilolo, double **inputrelololo,
   double **inputimhihihi, double **inputimlohihi,
   double **inputimhilohi, double **inputimlolohi,
   double **inputimhihilo, double **inputimlohilo,
   double **inputimhilolo, double **inputimlololo,
   double **outputrehihihi, double **outputrelohihi,
   double **outputrehilohi, double **outputrelolohi,
   double **outputrehihilo, double **outputrelohilo,
   double **outputrehilolo, double **outputrelololo,
   double **outputimhihihi, double **outputimlohihi,
   double **outputimhilohi, double **outputimlolohi,
   double **outputimhihilo, double **outputimlohilo,
   double **outputimhilolo, double **outputimlololo )
{
}
