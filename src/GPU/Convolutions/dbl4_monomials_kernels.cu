// The file dbl4_monomials_kernels.cu defines kernels specified
// in dbl4_monomials_kernels.h.

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
#include "dbl4_convolutions_kernels.h"
#include "dbl4_monomials_kernels.h"

__device__ void dbl4_convolute
 ( double *xhihi, double *xlohi, double *xhilo, double *xlolo,
   double *yhihi, double *ylohi, double *yhilo, double *ylolo,
   double *zhihi, double *zlohi, double *zhilo, double *zlolo,
   int dim, int k )
{
   double prdhihi,prdlohi,prdhilo,prdlolo;

   // z[k] = x[0]*y[k];
   qdg_mul(xhihi[0],xlohi[0],xhilo[0],xlolo[0],
           yhihi[k],ylohi[k],yhilo[k],ylolo[k],
           &zhihi[k],&zlohi[k],&zhilo[k],&zlolo[k]);

   for(int i=1; i<=k; i++) // z[k] = z[k] + x[i]*y[k-i];
   {
      qdg_mul(xhihi[i],xlohi[i],xhilo[i],xlolo[i],
              yhihi[k-i],ylohi[k-i],yhilo[k-i],ylolo[k-i],
              &prdhihi,&prdlohi,&prdhilo,&prdlolo);
      qdg_inc(&zhihi[k],&zlohi[k],&zhilo[k],&zlolo[k],
              prdhihi,prdlohi,prdhilo,prdlolo);
   }
}

__device__ void cmplx4_convolute
 ( double *xrehihi, double *xrelohi, double *xrehilo, double *xrelolo,
   double *ximhihi, double *ximlohi, double *ximhilo, double *ximlolo,
   double *yrehihi, double *yrelohi, double *yrehilo, double *yrelolo,
   double *yimhihi, double *yimlohi, double *yimhilo, double *yimlolo,
   double *zrehihi, double *zrelohi, double *zrehilo, double *zrelolo,
   double *zimhihi, double *zimlohi, double *zimhilo, double *zimlolo,
   int dim, int k )
{
   double xrhihi,xihihi,yrhihi,yihihi,zrhihi,zihihi,acchihi;
   double xrlohi,xilohi,yrlohi,yilohi,zrlohi,zilohi,acclohi;
   double xrhilo,xihilo,yrhilo,yihilo,zrhilo,zihilo,acchilo;
   double xrlolo,xilolo,yrlolo,yilolo,zrlolo,zilolo,acclolo;

   // z[k] = x[0]*y[k]
   xrhihi = xrehihi[0]; xrlohi = xrelohi[0];
   xrhilo = xrehilo[0]; xrlolo = xrelolo[0];
   xihihi = ximhihi[0]; xilohi = ximlohi[0];
   xihilo = ximhilo[0]; xilolo = ximlolo[0];
   yrhihi = yrehihi[k]; yrlohi = yrelohi[k];
   yrhilo = yrehilo[k]; yrlolo = yrelolo[k];
   yihihi = yimhihi[k]; yilohi = yimlohi[k];
   yihilo = yimhilo[k]; yilolo = yimlolo[k];

   qdg_mul(xrhihi,xrlohi,xrhilo,xrlolo,
           yrhihi,yrlohi,yrhilo,yrlolo,
           &zrhihi,&zrlohi,&zrhilo,&zrlolo);       // zr = xr*yr
   qdg_mul(xihihi,xilohi,xihilo,xilolo,
           yihihi,yilohi,yihilo,yilolo,
           &acchihi,&acclohi,&acchilo,&acclolo);   // acc = xi*yi
   qdg_minus(&acchihi,&acclohi,&acchilo,&acclolo);
   qdg_inc(&zrhihi,&zrlohi,&zrhilo,&zrlolo,
           acchihi,acclohi,acchilo,acclolo);       // zr = xr*yr - xi*yi
   qdg_mul(xrhihi,xrlohi,xrhilo,xrlolo,
           yihihi,yilohi,yihilo,yilolo,
           &zihihi,&zilohi,&zihilo,&zilolo);       // zi = xr*yi
   qdg_mul(xihihi,xilohi,xihilo,xilolo,
           yrhihi,yrlohi,yrhilo,yrlolo,
           &acchihi,&acclohi,&acchilo,&acclolo);   // acc = xi*yr
   qdg_inc(&zihihi,&zilohi,&zihilo,&zilolo,
           acchihi,acclohi,acchilo,acclolo);       // zr = xr*yr + xi*yi

   zrehihi[k] = zrhihi; zrelohi[k] = zrlohi;
   zrehilo[k] = zrhilo; zrelolo[k] = zrlolo;
   zimhihi[k] = zihihi; zimlohi[k] = zilohi;
   zimhilo[k] = zihilo; zimlolo[k] = zilolo;

   for(int i=1; i<=k; i++) // z[k] = z[k] + x[i]*y[k-i]
   {
      xrhihi = xrehihi[i]; xrlohi = xrelohi[i];
      xrhilo = xrehilo[i]; xrlolo = xrelolo[i];
      xihihi = ximhihi[i]; xilohi = ximlohi[i];
      xihilo = ximhilo[i]; xilolo = ximlolo[i];
      yrhihi = yrehihi[k-i]; yrlohi = yrelohi[k-i];
      yrhilo = yrehilo[k-i]; yrlolo = yrelolo[k-i];
      yihihi = yimhihi[k-i]; yilohi = yimlohi[k-i];
      yihilo = yimhilo[k-i]; yilolo = yimlolo[k-i];

      qdg_mul(xrhihi,xrlohi,xrhilo,xrlolo,
              yrhihi,yrlohi,yrhilo,yrlolo,
              &zrhihi,&zrlohi,&zrhilo,&zrlolo);       // zr = xr*yr
      qdg_mul(xihihi,xilohi,xihilo,xilolo,
              yihihi,yilohi,yihilo,yilolo,
              &acchihi,&acclohi,&acchilo,&acclolo);   // acc = xi*yi
      qdg_minus(&acchihi,&acclohi,&acchilo,&acclolo);
      qdg_inc(&zrhihi,&zrlohi,&zrhilo,&zrlolo,
              acchihi,acclohi,acchilo,acclolo);       // zr = xr*yr - xi*yi
      qdg_mul(xrhihi,xrlohi,xrhilo,xrlolo,
              yihihi,yilohi,yihilo,yilolo,
              &zihihi,&zilohi,&zihilo,&zilolo);       // zi = xr*yi
      qdg_mul(xihihi,xilohi,xihilo,xilolo,
              yrhihi,yrlohi,yrhilo,yrlolo,
              &acchihi,&acclohi,&acchilo,&acclolo);   // acc = xi*yr
      qdg_inc(&zihihi,&zilohi,&zihilo,&zilolo,
              acchihi,acclohi,acchilo,acclolo);       // zr = xr*yr + xi*yi

      qdg_inc(&zrehihi[k],&zrelohi[k],&zrehilo[k],&zrelolo[k],
              zrhihi,zrlohi,zrhilo,zrlolo);           // zre[k] += zr;
      qdg_inc(&zimhihi[k],&zimlohi[k],&zimhilo[k],&zimlolo[k],
              zihihi,zilohi,zihilo,zilolo);           // zim[k] += zi;
   }
}

__global__ void GPU_dbl4_speel
 ( int nvr, int deg, int *idx,
   double *cffhihi, double *cfflohi, double *cffhilo, double *cfflolo,
   double *inputhihi, double *inputlohi, double *inputhilo, double *inputlolo,
   double *forwardhihi, double *forwardlohi,
   double *forwardhilo, double *forwardlolo,
   double *backwardhihi, double *backwardlohi,
   double *backwardhilo, double *backwardlolo,
   double *crosshihi, double *crosslohi,
   double *crosshilo, double *crosslolo )
{
   const int k = threadIdx.x;
   const int deg1 = deg+1;
   int ix1,ix2;

   __shared__ double xvhihi[qd_shmemsize];
   __shared__ double xvlohi[qd_shmemsize];
   __shared__ double xvhilo[qd_shmemsize];
   __shared__ double xvlolo[qd_shmemsize];
   __shared__ double yvhihi[qd_shmemsize];
   __shared__ double yvlohi[qd_shmemsize];
   __shared__ double yvhilo[qd_shmemsize];
   __shared__ double yvlolo[qd_shmemsize];
   __shared__ double zvhihi[qd_shmemsize];
   __shared__ double zvlohi[qd_shmemsize];
   __shared__ double zvhilo[qd_shmemsize];
   __shared__ double zvlolo[qd_shmemsize];
  
   xvhihi[k] = cffhihi[k]; xvlohi[k] = cfflohi[k];
   xvhilo[k] = cffhilo[k]; xvlolo[k] = cfflolo[k];
   ix1 = idx[0]*deg1+k;
   yvhihi[k] = inputhihi[ix1]; yvlohi[k] = inputlohi[ix1];
   yvhilo[k] = inputhilo[ix1]; yvlolo[k] = inputlolo[ix1]; 
   __syncthreads();                                       // f[0] = cff*x[0]
   dbl4_convolute(xvhihi,xvlohi,xvhilo,xvlolo,
                  yvhihi,yvlohi,yvhilo,yvlolo,
                  zvhihi,zvlohi,zvhilo,zvlolo,deg1,k);
   forwardhihi[k] = zvhihi[k]; forwardlohi[k] = zvlohi[k];
   forwardhilo[k] = zvhilo[k]; forwardlolo[k] = zvlolo[k];

   for(int i=1; i<nvr; i++)
   {
      xvhihi[k] = zvhihi[k]; xvlohi[k] = zvlohi[k];
      xvhilo[k] = zvhilo[k]; xvlolo[k] = zvlolo[k];
      ix2 = idx[i]*deg1+k;
      yvhihi[k] = inputhihi[ix2]; yvlohi[k] = inputlohi[ix2];
      yvhilo[k] = inputhilo[ix2]; yvlolo[k] = inputlolo[ix2];
      __syncthreads();                                 // f[i] = f[i-1]*x[i]
      dbl4_convolute(xvhihi,xvlohi,xvhilo,xvlolo,
                     yvhihi,yvlohi,yvhilo,yvlolo,
                     zvhihi,zvlohi,zvhilo,zvlolo,deg1,k);
      ix1 = i*deg1+k;
      forwardhihi[ix1] = zvhihi[k]; forwardlohi[ix1] = zvlohi[k]; 
      forwardhilo[ix1] = zvhilo[k]; forwardlolo[ix1] = zvlolo[k];           
   }
   if(nvr > 2)
   {
      ix1 = idx[nvr-1]*deg1+k;
      xvhihi[k] = inputhihi[ix1]; xvlohi[k] = inputlohi[ix1]; 
      xvhilo[k] = inputhilo[ix1]; xvlolo[k] = inputlolo[ix1];
      ix2 = idx[nvr-2]*deg1+k;
      yvhihi[k] = inputhihi[ix2]; yvlohi[k] = inputlohi[ix2];
      yvhilo[k] = inputhilo[ix2]; yvlolo[k] = inputlolo[ix2];
      __syncthreads();                               // b[0] = x[n-1]*x[n-2]
      dbl4_convolute(xvhihi,xvlohi,xvhilo,xvlolo,
                     yvhihi,yvlohi,yvhilo,yvlolo,
                     zvhihi,zvlohi,zvhilo,zvlolo,deg1,k);
      backwardhihi[k] = zvhihi[k]; backwardlohi[k] = zvlohi[k];
      backwardhilo[k] = zvhilo[k]; backwardlolo[k] = zvlolo[k];                  
      for(int i=1; i<nvr-2; i++)
      {
         xvhihi[k] = zvhihi[k]; xvlohi[k] = zvlohi[k];
         xvhilo[k] = zvhilo[k]; xvlolo[k] = zvlolo[k];
         ix2 = idx[nvr-2-i]*deg1+k;
         yvhihi[k] = inputhihi[ix2]; yvlohi[k] = inputlohi[ix2];
         yvhilo[k] = inputhilo[ix2]; yvlolo[k] = inputlolo[ix2];
         __syncthreads();                          // b[i] = b[i-1]*x[n-2-i]
         dbl4_convolute(xvhihi,xvlohi,xvhilo,xvlolo,
                        yvhihi,yvlohi,yvhilo,yvlolo,
                        zvhihi,zvlohi,zvhilo,zvlolo,deg1,k);
         ix1 = i*deg1+k;
         backwardhihi[ix1] = zvhihi[k]; backwardlohi[ix1] = zvlohi[k];
         backwardhilo[ix1] = zvhilo[k]; backwardlolo[ix1] = zvlolo[k];        
      }
      xvhihi[k] = zvhihi[k]; xvlohi[k] = zvlohi[k];
      xvhilo[k] = zvhilo[k]; xvlolo[k] = zvlolo[k];
      yvhihi[k] = cffhihi[k]; yvlohi[k] = cfflohi[k];
      yvhilo[k] = cffhilo[k]; yvlolo[k] = cfflolo[k];
      __syncthreads();                                // b[n-3] = b[n-3]*cff
      dbl4_convolute(xvhihi,xvlohi,xvhilo,xvlolo,
                     yvhihi,yvlohi,yvhilo,yvlolo,
                     zvhihi,zvlohi,zvhilo,zvlolo,deg1,k);
      ix2 = (nvr-3)*deg1+k;
      backwardhihi[ix2] = zvhihi[k]; backwardlohi[ix2] = zvlohi[k];
      backwardhilo[ix2] = zvhilo[k]; backwardlolo[ix2] = zvlolo[k];

      if(nvr == 3)
      {
         xvhihi[k] = forwardhihi[k]; xvlohi[k] = forwardlohi[k];
         xvhilo[k] = forwardhilo[k]; xvlolo[k] = forwardlolo[k];
         ix2 = idx[2]*deg1+k;
         yvhihi[k] = inputhihi[ix2]; yvlohi[k] = inputlohi[ix2];
         yvhilo[k] = inputhilo[ix2]; yvlolo[k] = inputlolo[ix2];
         __syncthreads();                                // c[0] = f[0]*x[2]
         dbl4_convolute(xvhihi,xvlohi,xvhilo,xvlolo,
                        yvhihi,yvlohi,yvhilo,yvlolo,
                        zvhihi,zvlohi,zvhilo,zvlolo,deg1,k);
         crosshihi[k] = zvhihi[k]; crosslohi[k] = zvlohi[k];
         crosshilo[k] = zvhilo[k]; crosslolo[k] = zvlolo[k];
      }
      else
      {
         for(int i=0; i<nvr-3; i++)
         {
            ix1 = i*deg1+k; 
            xvhihi[k] = forwardhihi[ix1]; xvlohi[k] = forwardlohi[ix1];
            xvhilo[k] = forwardhilo[ix1]; xvlolo[k] = forwardlolo[ix1];
            ix2 = (nvr-4-i)*deg1+k;
            yvhihi[k] = backwardhihi[ix2]; yvlohi[k] = backwardlohi[ix2];
            yvhilo[k] = backwardhilo[ix2]; yvlolo[k] = backwardlolo[ix2];
            __syncthreads();                         // c[i] = f[i]*b[n-4-i]
            dbl4_convolute(xvhihi,xvlohi,xvhilo,xvlolo,
                           yvhihi,yvlohi,yvhilo,yvlolo,
                           zvhihi,zvlohi,zvhilo,zvlolo,deg1,k);
            crosshihi[ix1] = zvhihi[k]; crosslohi[ix1] = zvlohi[k];
            crosshilo[ix1] = zvhilo[k]; crosslolo[ix1] = zvlolo[k]; 
         }
         ix1 = (nvr-3)*deg1+k;
         xvhihi[k] = forwardhihi[ix1]; xvlohi[k] = forwardlohi[ix1];
         xvhilo[k] = forwardhilo[ix1]; xvlolo[k] = forwardlolo[ix1];
         ix2 = idx[nvr-1]*deg1+k;
         yvhihi[k] = inputhihi[ix2]; yvlohi[k] = inputlohi[ix2];
         yvhilo[k] = inputhilo[ix2]; yvlolo[k] = inputlolo[ix2];
         __syncthreads();                          // c[n-3] = f[n-3]*x[n-1]
         dbl4_convolute(xvhihi,xvlohi,xvhilo,xvlolo,
                        yvhihi,yvlohi,yvhilo,yvlolo,
                        zvhihi,zvlohi,zvhilo,zvlolo,deg1,k);
         crosshihi[ix1] = zvhihi[k]; crosslohi[ix1] = zvlohi[k];
         crosshilo[ix1] = zvhilo[k]; crosslolo[ix1] = zvlolo[k];
      }
   }
}

__global__ void GPU_cmplx4_speel
 ( int nvr, int deg, int *idx,
   double *cffrehihi, double *cffrelohi, double *cffrehilo, double *cffrelolo,
   double *cffimhihi, double *cffimlohi, double *cffimhilo, double *cffimlolo,
   double *inputrehihi, double *inputrelohi,
   double *inputrehilo, double *inputrelolo,
   double *inputimhihi, double *inputimlohi,
   double *inputimhilo, double *inputimlolo,
   double *forwardrehihi, double *forwardrelohi,
   double *forwardrehilo, double *forwardrelolo,
   double *forwardimhihi, double *forwardimlohi,
   double *forwardimhilo, double *forwardimlolo,
   double *backwardrehihi, double *backwardrelohi,
   double *backwardrehilo, double *backwardrelolo,
   double *backwardimhihi, double *backwardimlohi,
   double *backwardimhilo, double *backwardimlolo,
   double *crossrehihi, double *crossrelohi,
   double *crossrehilo, double *crossrelolo,
   double *crossimhihi, double *crossimlohi,
   double *crossimhilo, double *crossimlolo )
{
   const int k = threadIdx.x;
   const int deg1 = deg+1;
   int ix1,ix2;

   __shared__ double xvrehihi[qd_shmemsize];
   __shared__ double xvrelohi[qd_shmemsize];
   __shared__ double xvrehilo[qd_shmemsize];
   __shared__ double xvrelolo[qd_shmemsize];
   __shared__ double xvimhihi[qd_shmemsize];
   __shared__ double xvimlohi[qd_shmemsize];
   __shared__ double xvimhilo[qd_shmemsize];
   __shared__ double xvimlolo[qd_shmemsize];
   __shared__ double yvrehihi[qd_shmemsize];
   __shared__ double yvrelohi[qd_shmemsize];
   __shared__ double yvrehilo[qd_shmemsize];
   __shared__ double yvrelolo[qd_shmemsize];
   __shared__ double yvimhihi[qd_shmemsize];
   __shared__ double yvimlohi[qd_shmemsize];
   __shared__ double yvimhilo[qd_shmemsize];
   __shared__ double yvimlolo[qd_shmemsize];
   __shared__ double zvrehihi[qd_shmemsize];
   __shared__ double zvrelohi[qd_shmemsize];
   __shared__ double zvrehilo[qd_shmemsize];
   __shared__ double zvrelolo[qd_shmemsize];
   __shared__ double zvimhihi[qd_shmemsize];
   __shared__ double zvimlohi[qd_shmemsize];
   __shared__ double zvimhilo[qd_shmemsize];
   __shared__ double zvimlolo[qd_shmemsize];

   xvrehihi[k] = cffrehihi[k]; xvrelohi[k] = cffrelohi[k];
   xvrehilo[k] = cffrehilo[k]; xvrelolo[k] = cffrelolo[k];
   xvimhihi[k] = cffimhihi[k]; xvimlohi[k] = cffimlohi[k];
   xvimhilo[k] = cffimhilo[k]; xvimlolo[k] = cffimlolo[k];
   ix1 = idx[0]*deg1+k;
   yvrehihi[k] = inputrehihi[ix1]; yvrelohi[k] = inputrelohi[ix1];
   yvrehilo[k] = inputrehilo[ix1]; yvrelolo[k] = inputrelolo[ix1];
   yvimhihi[k] = inputimhihi[ix1]; yvimlohi[k] = inputimlohi[ix1];
   yvimhilo[k] = inputimhilo[ix1]; yvimlolo[k] = inputimlolo[ix1];
   __syncthreads();                                       // f[0] = cff*x[0] 
   cmplx4_convolute(xvrehihi,xvrelohi,xvrehilo,xvrelolo,
                    xvimhihi,xvimlohi,xvimhilo,xvimlolo,
                    yvrehihi,yvrelohi,yvrehilo,yvrelolo,
                    yvimhihi,yvimlohi,yvimhilo,yvimlolo,
                    zvrehihi,zvrelohi,zvrehilo,zvrelolo,
                    zvimhihi,zvimlohi,zvimhilo,zvimlolo,deg1,k);
   forwardrehihi[k] = zvrehihi[k]; forwardrelohi[k] = zvrelohi[k];
   forwardrehilo[k] = zvrehilo[k]; forwardrelolo[k] = zvrelolo[k];
   forwardimhihi[k] = zvimhihi[k]; forwardimlohi[k] = zvimlohi[k];
   forwardimhilo[k] = zvimhilo[k]; forwardimlolo[k] = zvimlolo[k];

   for(int i=1; i<nvr; i++)
   {
      xvrehihi[k] = zvrehihi[k]; xvrelohi[k] = zvrelohi[k];
      xvrehilo[k] = zvrehilo[k]; xvrelolo[k] = zvrelolo[k];
      xvimhihi[k] = zvimhihi[k]; xvimlohi[k] = zvimlohi[k];
      xvimhilo[k] = zvimhilo[k]; xvimlolo[k] = zvimlolo[k];
      ix2 = idx[i]*deg1+k;
      yvrehihi[k] = inputrehihi[ix2]; yvrelohi[k] = inputrelohi[ix2];
      yvrehilo[k] = inputrehilo[ix2]; yvrelolo[k] = inputrelolo[ix2];
      yvimhihi[k] = inputimhihi[ix2]; yvimlohi[k] = inputimlohi[ix2];
      yvimhilo[k] = inputimhilo[ix2]; yvimlolo[k] = inputimlolo[ix2];
      __syncthreads();                                 // f[i] = f[i-i]*x[i]
      cmplx4_convolute(xvrehihi,xvrelohi,xvrehilo,xvrelolo,
                       xvimhihi,xvimlohi,xvimhilo,xvimlolo,
                       yvrehihi,yvrelohi,yvrehilo,yvrelolo,
                       yvimhihi,yvimlohi,yvimhilo,yvimlolo,
                       zvrehihi,zvrelohi,zvrehilo,zvrelolo,
                       zvimhihi,zvimlohi,zvimhilo,zvimlolo,deg1,k);
      ix1 = i*deg1+k;                                   
      forwardrehihi[ix1] = zvrehihi[k]; forwardrelohi[ix1] = zvrelohi[k];
      forwardrehilo[ix1] = zvrehilo[k]; forwardrelolo[ix1] = zvrelolo[k];
      forwardimhihi[ix1] = zvimhihi[k]; forwardimlohi[ix1] = zvimlohi[k];
      forwardimhilo[ix1] = zvimhilo[k]; forwardimlolo[ix1] = zvimlolo[k]; 
   }
   if(nvr > 2)
   {
      ix1 = idx[nvr-1]*deg1+k;
      xvrehihi[k] = inputrehihi[ix1]; xvrelohi[k] = inputrelohi[ix1];
      xvrehilo[k] = inputrehilo[ix1]; xvrelolo[k] = inputrelolo[ix1];
      xvimhihi[k] = inputimhihi[ix1]; xvimlohi[k] = inputimlohi[ix1];
      xvimhilo[k] = inputimhilo[ix1]; xvimlolo[k] = inputimlolo[ix1];
      ix2 = idx[nvr-2]*deg1+k;
      yvrehihi[k] = inputrehihi[ix2]; yvrelohi[k] = inputrelohi[ix2];
      yvrehilo[k] = inputrehilo[ix2]; yvrelolo[k] = inputrelolo[ix2];
      yvimhihi[k] = inputimhihi[ix2]; yvimlohi[k] = inputimlohi[ix2];
      yvimhilo[k] = inputimhilo[ix2]; yvimlolo[k] = inputimlolo[ix2];
      __syncthreads();                               // b[0] = x[n-1]*x[n-2]
      cmplx4_convolute(xvrehihi,xvrelohi,xvrehilo,xvrelolo,
                       xvimhihi,xvimlohi,xvimhilo,xvimlolo,
                       yvrehihi,yvrelohi,yvrehilo,yvrelolo,
                       yvimhihi,yvimlohi,yvimhilo,yvimlolo,
                       zvrehihi,zvrelohi,zvrehilo,zvrelolo,
                       zvimhihi,zvimlohi,zvimhilo,zvimlolo,deg1,k);
      backwardrehihi[k] = zvrehihi[k]; backwardrelohi[k] = zvrelohi[k];
      backwardrehilo[k] = zvrehilo[k]; backwardrelolo[k] = zvrelolo[k];
      backwardimhihi[k] = zvimhihi[k]; backwardimlohi[k] = zvimlohi[k];
      backwardimhilo[k] = zvimhilo[k]; backwardimlolo[k] = zvimlolo[k];

      for(int i=1; i<nvr-2; i++)
      {
         xvrehihi[k] = zvrehihi[k]; xvrelohi[k] = zvrelohi[k];
         xvrehilo[k] = zvrehilo[k]; xvrelolo[k] = zvrelolo[k];
         xvimhihi[k] = zvimhihi[k]; xvimlohi[k] = zvimlohi[k];
         xvimhilo[k] = zvimhilo[k]; xvimlolo[k] = zvimlolo[k];
         ix2 = idx[nvr-2-i]*deg1+k;
         yvrehihi[k] = inputrehihi[ix2]; yvrelohi[k] = inputrelohi[ix2];
         yvrehilo[k] = inputrehilo[ix2]; yvrelolo[k] = inputrelolo[ix2];
         yvimhihi[k] = inputimhihi[ix2]; yvimlohi[k] = inputimlohi[ix2];
         yvimhilo[k] = inputimhilo[ix2]; yvimlolo[k] = inputimlolo[ix2];
         __syncthreads();                           // b[i] = b[i]*x[n-2-i]
         cmplx4_convolute(xvrehihi,xvrelohi,xvrehilo,xvrelolo,
                          xvimhihi,xvimlohi,xvimhilo,xvimlolo,
                          yvrehihi,yvrelohi,yvrehilo,yvrelolo,
                          yvimhihi,yvimlohi,yvimhilo,yvimlolo,
                          zvrehihi,zvrelohi,zvrehilo,zvrelolo,
                          zvimhihi,zvimlohi,zvimhilo,zvimlolo,deg1,k);
         ix1 = i*deg1+k;
         backwardrehihi[ix1] = zvrehihi[k]; backwardrelohi[ix1] = zvrelohi[k];
         backwardrehilo[ix1] = zvrehilo[k]; backwardrelolo[ix1] = zvrelolo[k];
         backwardimhihi[ix1] = zvimhihi[k]; backwardimlohi[ix1] = zvimlohi[k];
         backwardimhilo[ix1] = zvimhilo[k]; backwardimlolo[ix1] = zvimlolo[k];
      }
      xvrehihi[k] = zvrehihi[k]; xvrelohi[k] = zvrelohi[k];
      xvrehilo[k] = zvrehilo[k]; xvrelolo[k] = zvrelolo[k];
      xvimhihi[k] = zvimhihi[k]; xvimlohi[k] = zvimlohi[k];
      xvimhilo[k] = zvimhilo[k]; xvimlolo[k] = zvimlolo[k];
      yvrehihi[k] = cffrehihi[k]; yvrelohi[k] = cffrelohi[k];
      yvrehilo[k] = cffrehilo[k]; yvrelolo[k] = cffrelolo[k];
      yvimhihi[k] = cffimhihi[k]; yvimlohi[k] = cffimlohi[k];
      yvimhilo[k] = cffimhilo[k]; yvimlolo[k] = cffimlolo[k];
      __syncthreads();                               // b[n-3] = b[n-3]*cff
      cmplx4_convolute(xvrehihi,xvrelohi,xvrehilo,xvrelolo,
                       xvimhihi,xvimlohi,xvimhilo,xvimlolo,
                       yvrehihi,yvrelohi,yvrehilo,yvrelolo,
                       yvimhihi,yvimlohi,yvimhilo,yvimlolo,
                       zvrehihi,zvrelohi,zvrehilo,zvrelolo,
                       zvimhihi,zvimlohi,zvimhilo,zvimlolo,deg1,k);
      ix1 = (nvr-3)*deg1+k;
      backwardrehihi[ix1] = zvrehihi[k]; backwardrelohi[ix1] = zvrelohi[k];
      backwardrehilo[ix1] = zvrehilo[k]; backwardrelolo[ix1] = zvrelolo[k];
      backwardimhihi[ix1] = zvimhihi[k]; backwardimlohi[ix1] = zvimlohi[k];
      backwardimhilo[ix1] = zvimhilo[k]; backwardimlolo[ix1] = zvimlolo[k];

      if(nvr == 3)
      {
         xvrehihi[k] = forwardrehihi[k]; xvrelohi[k] = forwardrelohi[k];
         xvrehilo[k] = forwardrehilo[k]; xvrelolo[k] = forwardrelolo[k];
         xvimhihi[k] = forwardimhihi[k]; xvimlohi[k] = forwardimlohi[k];
         xvimhilo[k] = forwardimhilo[k]; xvimlolo[k] = forwardimlolo[k];
         ix2 = idx[2]*deg1+k;
         yvrehihi[k] = inputrehihi[ix2]; yvrelohi[k] = inputrelohi[ix2];
         yvrehilo[k] = inputrehilo[ix2]; yvrelolo[k] = inputrelolo[ix2];
         yvimhihi[k] = inputimhihi[ix2]; yvimlohi[k] = inputimlohi[ix2];
         yvimhilo[k] = inputimhilo[ix2]; yvimlolo[k] = inputimlolo[ix2];
         __syncthreads();                               // c[0] = f[0]*x[2]
         cmplx4_convolute(xvrehihi,xvrelohi,xvrehilo,xvrelolo,
                          xvimhihi,xvimlohi,xvimhilo,xvimlolo,
                          yvrehihi,yvrelohi,yvrehilo,yvrelolo,
                          yvimhihi,yvimlohi,yvimhilo,yvimlolo,
                          zvrehihi,zvrelohi,zvrehilo,zvrelolo,
                          zvimhihi,zvimlohi,zvimhilo,zvimlolo,deg1,k);
         crossrehihi[k] = zvrehihi[k]; crossrelohi[k] = zvrelohi[k];
         crossrehilo[k] = zvrehilo[k]; crossrelolo[k] = zvrelolo[k];
         crossimhihi[k] = zvimhihi[k]; crossimlohi[k] = zvimlohi[k];
         crossimhilo[k] = zvimhilo[k]; crossimlolo[k] = zvimlolo[k];
      }
      else
      {
         for(int i=0; i<nvr-3; i++)
         {
            ix1 = i*deg1+k;   
            xvrehihi[k] = forwardrehihi[ix1]; xvrelohi[k] = forwardrelohi[ix1];
            xvrehilo[k] = forwardrehilo[ix1]; xvrelolo[k] = forwardrelolo[ix1];
            xvimhihi[k] = forwardimhihi[ix1]; xvimlohi[k] = forwardimlohi[ix1];
            xvimhilo[k] = forwardimhilo[ix1]; xvimlolo[k] = forwardimlolo[ix1];
            ix2 = (nvr-4-i)*deg1+k;
            yvrehihi[k] = backwardrehihi[ix2];
            yvrelohi[k] = backwardrelohi[ix2];
            yvrehilo[k] = backwardrehilo[ix2];
            yvrelolo[k] = backwardrelolo[ix2];
            yvimhihi[k] = backwardimhihi[ix2];
            yvimlohi[k] = backwardimlohi[ix2];
            yvimhilo[k] = backwardimhilo[ix2];
            yvimlolo[k] = backwardimlolo[ix2];
            __syncthreads();                        // c[i] = f[i]*b[n-4-i]
            cmplx4_convolute(xvrehihi,xvrelohi,xvrehilo,xvrelolo,
                             xvimhihi,xvimlohi,xvimhilo,xvimlolo,
                             yvrehihi,yvrelohi,yvrehilo,yvrelolo,
                             yvimhihi,yvimlohi,yvimhilo,yvimlolo,
                             zvrehihi,zvrelohi,zvrehilo,zvrelolo,
                             zvimhihi,zvimlohi,zvimhilo,zvimlolo,deg1,k);
            ix1 = i*deg1+k;
            crossrehihi[ix1] = zvrehihi[k]; crossrelohi[ix1] = zvrelohi[k];
            crossrehilo[ix1] = zvrehilo[k]; crossrelolo[ix1] = zvrelolo[k];
            crossimhihi[ix1] = zvimhihi[k]; crossimlohi[ix1] = zvimlohi[k];
            crossimhilo[ix1] = zvimhilo[k]; crossimlolo[ix1] = zvimlolo[k];
         }
         ix1 = (nvr-3)*deg1+k;
         xvrehihi[k] = forwardrehihi[ix1]; xvrelohi[k] = forwardrelohi[ix1];
         xvrehilo[k] = forwardrehilo[ix1]; xvrelolo[k] = forwardrelolo[ix1];
         xvimhihi[k] = forwardimhihi[ix1]; xvimlohi[k] = forwardimlohi[ix1];
         xvimhilo[k] = forwardimhilo[ix1]; xvimlolo[k] = forwardimlolo[ix1];
         ix2 = idx[nvr-1]*deg1+k;
         yvrehihi[k] = inputrehihi[ix2]; yvrelohi[k] = inputrelohi[ix2];
         yvrehilo[k] = inputrehilo[ix2]; yvrelolo[k] = inputrelolo[ix2];
         yvimhihi[k] = inputimhihi[ix2]; yvimlohi[k] = inputimlohi[ix2];
         yvimhilo[k] = inputimhilo[ix2]; yvimlolo[k] = inputimlolo[ix2];
         __syncthreads();                         // c[n-3] = f[n-3]*x[n-1]
         cmplx4_convolute(xvrehihi,xvrelohi,xvrehilo,xvrelolo,
                          xvimhihi,xvimlohi,xvimhilo,xvimlolo,
                          yvrehihi,yvrelohi,yvrehilo,yvrelolo,
                          yvimhihi,yvimlohi,yvimhilo,yvimlolo,
                          zvrehihi,zvrelohi,zvrehilo,zvrelolo,
                          zvimhihi,zvimlohi,zvimhilo,zvimlolo,deg1,k);
         ix1 = (nvr-3)*deg1+k;
         crossrehihi[ix1] = zvrehihi[k]; crossrelohi[ix1] = zvrelohi[k];
         crossrehilo[ix1] = zvrehilo[k]; crossrelolo[ix1] = zvrelolo[k];
         crossimhihi[ix1] = zvimhihi[k]; crossimlohi[ix1] = zvimlohi[k];
         crossimhilo[ix1] = zvimhilo[k]; crossimlolo[ix1] = zvimlolo[k];
      }
   }
}

void GPU_dbl4_evaldiff
 ( int BS, int dim, int nvr, int deg, int *idx,
   double *cffhihi, double *cfflohi, double *cffhilo, double *cfflolo,
   double **inputhihi, double **inputlohi,
   double **inputhilo, double **inputlolo,
   double **outputhihi, double **outputlohi,
   double **outputhilo, double **outputlolo )
{
   const int deg1 = deg+1;   // length of all vectors
   double *inputhihi_d;      // inputhihi_d is input on the device
   double *inputlohi_d;      // inputlohi_d is input on the device
   double *inputhilo_d;      // inputhilo_d is input on the device
   double *inputlolo_d;      // inputlolo_d is input on the device
   double *forwardhihi_d;    // highest forward products on the device
   double *forwardlohi_d;    // second highest forward products on the device
   double *forwardhilo_d;    // second lowest forward products on the device
   double *forwardlolo_d;    // lowest forward products on the device
   double *backwardhihi_d;   // highest backward products on the device
   double *backwardlohi_d;   // second highest backward products on the device
   double *backwardhilo_d;   // second lowest backward products on the device
   double *backwardlolo_d;   // lowest backward products on the device
   double *crosshihi_d;      // highest cross products on the device
   double *crosslohi_d;      // second highest cross products on the device
   double *crosshilo_d;      // second lowest cross products on the device
   double *crosslolo_d;      // lowest cross products on the device
   double *cffhihi_d;        // cffhihi_d is cffhihi on device
   double *cfflohi_d;        // cfflohi_d is cfflohi on device
   double *cffhilo_d;        // cffhilo_d is cffhilo on device
   double *cfflolo_d;        // cfflolo_d is cfflolo on device
   int *idx_d;               // idx_d is idx on device

   size_t szcff = deg1*sizeof(double);
   size_t szdim = dim*(deg1)*sizeof(double);
   size_t sznvr = nvr*(deg1)*sizeof(double);
   size_t sznvr2 = (nvr-2)*(deg1)*sizeof(double);
   size_t szidx = nvr*sizeof(int);

   cudaMalloc((void**)&idx_d,szidx);
   cudaMalloc((void**)&cffhihi_d,szcff);
   cudaMalloc((void**)&cfflohi_d,szcff);
   cudaMalloc((void**)&cffhilo_d,szcff);
   cudaMalloc((void**)&cfflolo_d,szcff);
   cudaMalloc((void**)&inputhihi_d,szdim);
   cudaMalloc((void**)&inputlohi_d,szdim);
   cudaMalloc((void**)&inputhilo_d,szdim);
   cudaMalloc((void**)&inputlolo_d,szdim);
   cudaMalloc((void**)&forwardhihi_d,sznvr);
   cudaMalloc((void**)&forwardlohi_d,sznvr);
   cudaMalloc((void**)&forwardhilo_d,sznvr);
   cudaMalloc((void**)&forwardlolo_d,sznvr);
   cudaMalloc((void**)&backwardhihi_d,sznvr2);
   cudaMalloc((void**)&backwardlohi_d,sznvr2);
   cudaMalloc((void**)&backwardhilo_d,sznvr2);
   cudaMalloc((void**)&backwardlolo_d,sznvr2);
   cudaMalloc((void**)&crosshihi_d,sznvr2);
   cudaMalloc((void**)&crosslohi_d,sznvr2);
   cudaMalloc((void**)&crosshilo_d,sznvr2);
   cudaMalloc((void**)&crosslolo_d,sznvr2);

   double *inputhihi_h = new double[dim*(deg1)];
   double *inputlohi_h = new double[dim*(deg1)];
   double *inputhilo_h = new double[dim*(deg1)];
   double *inputlolo_h = new double[dim*(deg1)];
   int ix = 0;
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++)
      {
         inputhihi_h[ix] = inputhihi[i][j];
         inputlohi_h[ix] = inputlohi[i][j];
         inputhilo_h[ix] = inputhilo[i][j];
         inputlolo_h[ix++] = inputlolo[i][j];
      }

   cudaMemcpy(idx_d,idx,szidx,cudaMemcpyHostToDevice);
   cudaMemcpy(cffhihi_d,cffhihi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cfflohi_d,cfflohi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffhilo_d,cffhilo,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cfflolo_d,cfflolo,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(inputhihi_d,inputhihi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputlohi_d,inputlohi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputhilo_d,inputhilo_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputlolo_d,inputlolo_h,szdim,cudaMemcpyHostToDevice);

   if(BS == deg1)
   {
      GPU_dbl4_speel<<<1,BS>>>
         (nvr,deg,idx_d,cffhihi_d,cfflohi_d,cffhilo_d,cfflolo_d,
          inputhihi_d,inputlohi_d,inputhilo_d,inputlolo_d,
          forwardhihi_d,forwardlohi_d,forwardhilo_d,forwardlolo_d,
          backwardhihi_d,backwardlohi_d,backwardhilo_d,backwardlolo_d,
          crosshihi_d,crosslohi_d,crosshilo_d,crosslolo_d);
   }
   double *forwardhihi_h = new double[(deg1)*nvr];
   double *forwardlohi_h = new double[(deg1)*nvr];
   double *forwardhilo_h = new double[(deg1)*nvr];
   double *forwardlolo_h = new double[(deg1)*nvr];
   double *backwardhihi_h = new double[(deg1)*(nvr-2)];
   double *backwardlohi_h = new double[(deg1)*(nvr-2)];
   double *backwardhilo_h = new double[(deg1)*(nvr-2)];
   double *backwardlolo_h = new double[(deg1)*(nvr-2)];
   double *crosshihi_h = new double[(deg1)*(nvr-2)];
   double *crosslohi_h = new double[(deg1)*(nvr-2)];
   double *crosshilo_h = new double[(deg1)*(nvr-2)];
   double *crosslolo_h = new double[(deg1)*(nvr-2)];
  
   cudaMemcpy(forwardhihi_h,forwardhihi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardlohi_h,forwardlohi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardhilo_h,forwardhilo_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardlolo_h,forwardlolo_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardhihi_h,backwardhihi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardlohi_h,backwardlohi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardhilo_h,backwardhilo_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardlolo_h,backwardlolo_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crosshihi_h,crosshihi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crosslohi_h,crosslohi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crosshilo_h,crosshilo_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crosslolo_h,crosslolo_d,sznvr2,cudaMemcpyDeviceToHost);

   int offset = (nvr-1)*deg1;            // assign value of the monomial
   for(int i=0; i<deg1; i++)
   {
      outputhihi[dim][i] = forwardhihi_h[offset+i];
      outputlohi[dim][i] = forwardlohi_h[offset+i];
      outputhilo[dim][i] = forwardhilo_h[offset+i];
      outputlolo[dim][i] = forwardlolo_h[offset+i];
   }
   ix = idx[nvr-1];                      // derivative with respect to x[n-1]
   offset = (nvr-2)*deg1;
   for(int i=0; i<deg1; i++)
   {
      outputhihi[ix][i] = forwardhihi_h[offset+i];
      outputlohi[ix][i] = forwardlohi_h[offset+i];
      outputhilo[ix][i] = forwardhilo_h[offset+i];
      outputlolo[ix][i] = forwardlolo_h[offset+i];
   }
   ix = idx[0];                          // derivative with respect to x[0]
   offset = (nvr-3)*deg1;
   for(int i=0; i<deg1; i++)
   {
      outputhihi[ix][i] = backwardhihi_h[offset+i];
      outputlohi[ix][i] = backwardlohi_h[offset+i];
      outputhilo[ix][i] = backwardhilo_h[offset+i];
      outputlolo[ix][i] = backwardlolo_h[offset+i];
   }
   for(int k=1; k<nvr-1; k++)            // derivative with respect to x[k]
   {
      ix = idx[k]; offset = (k-1)*deg1;
      for(int i=0; i<deg1; i++)
      {
         outputhihi[ix][i] = crosshihi_h[offset+i];
         outputlohi[ix][i] = crosslohi_h[offset+i];
         outputhilo[ix][i] = crosshilo_h[offset+i];
         outputlolo[ix][i] = crosslolo_h[offset+i];
      }
   }
}

void GPU_cmplx4_evaldiff
 ( int BS, int dim, int nvr, int deg, int *idx,
   double *cffrehihi, double *cffrelohi, double *cffrehilo, double *cffrelolo,
   double *cffimhihi, double *cffimlohi, double *cffimhilo, double *cffimlolo,
   double **inputrehihi, double **inputrelohi,
   double **inputrehilo, double **inputrelolo,
   double **inputimhihi, double **inputimlohi,
   double **inputimhilo, double **inputimlolo,
   double **outputrehihi, double **outputrelohi,
   double **outputrehilo, double **outputrelolo,
   double **outputimhihi, double **outputimlohi,
   double **outputimhilo, double **outputimlolo )
{
   const int deg1 = deg+1;      // length of all vectors
   double *inputrehihi_d;       // inputrehihi_d is inputrehihi on the device
   double *inputrelohi_d;       // inputrelohi_d is inputrelohi on the device
   double *inputrehilo_d;       // inputrehilo_d is inputrehilo on the device
   double *inputrelolo_d;       // inputrelolo_d is inputrelolo on the device
   double *inputimhihi_d;       // inputimhihi_d is inputrehihi on the device
   double *inputimlohi_d;       // inputimlohi_d is inputrelohi on the device
   double *inputimhilo_d;       // inputimhilo_d is inputrehilo on the device
   double *inputimlolo_d;       // inputimlolo_d is inputrelolo on the device
   double *forwardrehihi_d;
   double *forwardrelohi_d;
   double *forwardrehilo_d;
   double *forwardrelolo_d;
   double *forwardimhihi_d;
   double *forwardimlohi_d;
   double *forwardimhilo_d;
   double *forwardimlolo_d;
   double *backwardrehihi_d;
   double *backwardrelohi_d;
   double *backwardrehilo_d;
   double *backwardrelolo_d;
   double *backwardimhihi_d;
   double *backwardimlohi_d;
   double *backwardimhilo_d;
   double *backwardimlolo_d;
   double *crossrehihi_d;
   double *crossrelohi_d;
   double *crossrehilo_d;
   double *crossrelolo_d;
   double *crossimhihi_d;
   double *crossimlohi_d;
   double *crossimhilo_d;
   double *crossimlolo_d;
   double *cffrehihi_d;             // cffrehihi_d is cffrehihi on the device
   double *cffrelohi_d;             // cffrelohi_d is cffrelohi on the device
   double *cffrehilo_d;             // cffrehilo_d is cffrehilo on the device
   double *cffrelolo_d;             // cffrelolo_d is cffrelolo on the device
   double *cffimhihi_d;             // cffimhihi_d is cffimhihi on the device
   double *cffimlohi_d;             // cffimlohi_d is cffimlohi on the device
   double *cffimhilo_d;             // cffimhilo_d is cffimhilo on the device
   double *cffimlolo_d;             // cffimlolo_d is cffimlolo on the device
   int *idx_d;                      // idx_d is idx on the device

   size_t szdim = dim*(deg1)*sizeof(double);
   size_t sznvr = nvr*(deg1)*sizeof(double);
   size_t sznvr2 = (nvr-2)*(deg1)*sizeof(double);
   size_t szidx = nvr*sizeof(int);
   size_t szcff = deg1*sizeof(double);

   cudaMalloc((void**)&idx_d,szidx);
   cudaMalloc((void**)&cffrehihi_d,szcff);
   cudaMalloc((void**)&cffrelohi_d,szcff);
   cudaMalloc((void**)&cffrehilo_d,szcff);
   cudaMalloc((void**)&cffrelolo_d,szcff);
   cudaMalloc((void**)&cffimhihi_d,szcff);
   cudaMalloc((void**)&cffimlohi_d,szcff);
   cudaMalloc((void**)&cffimhilo_d,szcff);
   cudaMalloc((void**)&cffimlolo_d,szcff);
   cudaMalloc((void**)&inputrehihi_d,szdim);
   cudaMalloc((void**)&inputrelohi_d,szdim);
   cudaMalloc((void**)&inputrehilo_d,szdim);
   cudaMalloc((void**)&inputrelolo_d,szdim);
   cudaMalloc((void**)&inputimhihi_d,szdim);
   cudaMalloc((void**)&inputimlohi_d,szdim);
   cudaMalloc((void**)&inputimhilo_d,szdim);
   cudaMalloc((void**)&inputimlolo_d,szdim);
   cudaMalloc((void**)&forwardrehihi_d,sznvr);
   cudaMalloc((void**)&forwardrelohi_d,sznvr);
   cudaMalloc((void**)&forwardrehilo_d,sznvr);
   cudaMalloc((void**)&forwardrelolo_d,sznvr);
   cudaMalloc((void**)&forwardimhihi_d,sznvr);
   cudaMalloc((void**)&forwardimlohi_d,sznvr);
   cudaMalloc((void**)&forwardimhilo_d,sznvr);
   cudaMalloc((void**)&forwardimlolo_d,sznvr);
   cudaMalloc((void**)&backwardrehihi_d,sznvr2);
   cudaMalloc((void**)&backwardrelohi_d,sznvr2);
   cudaMalloc((void**)&backwardrehilo_d,sznvr2);
   cudaMalloc((void**)&backwardrelolo_d,sznvr2);
   cudaMalloc((void**)&backwardimhihi_d,sznvr2);
   cudaMalloc((void**)&backwardimlohi_d,sznvr2);
   cudaMalloc((void**)&backwardimhilo_d,sznvr2);
   cudaMalloc((void**)&backwardimlolo_d,sznvr2);
   cudaMalloc((void**)&crossrehihi_d,sznvr2);
   cudaMalloc((void**)&crossrelohi_d,sznvr2);
   cudaMalloc((void**)&crossrehilo_d,sznvr2);
   cudaMalloc((void**)&crossrelolo_d,sznvr2);
   cudaMalloc((void**)&crossimhihi_d,sznvr2);
   cudaMalloc((void**)&crossimlohi_d,sznvr2);
   cudaMalloc((void**)&crossimhilo_d,sznvr2);
   cudaMalloc((void**)&crossimlolo_d,sznvr2);

   double *inputrehihi_h = new double[dim*(deg1)];
   double *inputrelohi_h = new double[dim*(deg1)];
   double *inputrehilo_h = new double[dim*(deg1)];
   double *inputrelolo_h = new double[dim*(deg1)];
   double *inputimhihi_h = new double[dim*(deg1)];
   double *inputimlohi_h = new double[dim*(deg1)];
   double *inputimhilo_h = new double[dim*(deg1)];
   double *inputimlolo_h = new double[dim*(deg1)];
   int ix = 0;
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++)
      {
         inputrehihi_h[ix] = inputrehihi[i][j];
         inputrelohi_h[ix] = inputrelohi[i][j];
         inputrehilo_h[ix] = inputrehilo[i][j];
         inputrelolo_h[ix] = inputrelolo[i][j];
         inputimhihi_h[ix] = inputimhihi[i][j];
         inputimlohi_h[ix] = inputimlohi[i][j];
         inputimhilo_h[ix] = inputimhilo[i][j];
         inputimlolo_h[ix++] = inputimlolo[i][j];
      }

   cudaMemcpy(idx_d,idx,szidx,cudaMemcpyHostToDevice);
   cudaMemcpy(cffrehihi_d,cffrehihi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffrelohi_d,cffrelohi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffrehilo_d,cffrehilo,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffrelolo_d,cffrelolo,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffimhihi_d,cffimhihi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffimlohi_d,cffimlohi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffimhilo_d,cffimhilo,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffimlolo_d,cffimlolo,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(inputrehihi_d,inputrehihi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputrelohi_d,inputrelohi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputrehilo_d,inputrehilo_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputrelolo_d,inputrelolo_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputimhihi_d,inputimhihi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputimlohi_d,inputimlohi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputimhilo_d,inputimhilo_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputimlolo_d,inputimlolo_h,szdim,cudaMemcpyHostToDevice);

   if(BS == deg1)
   {
      GPU_cmplx4_speel<<<1,BS>>>
         (nvr,deg,idx_d,
          cffrehihi_d,cffrelohi_d,cffrehilo_d,cffrelolo_d,
          cffimhihi_d,cffimlohi_d,cffimhilo_d,cffimlolo_d,
          inputrehihi_d,inputrelohi_d,inputrehilo_d,inputrelolo_d,
          inputimhihi_d,inputimlohi_d,inputimhilo_d,inputimlolo_d,
          forwardrehihi_d,forwardrelohi_d,forwardrehilo_d,forwardrelolo_d,
          forwardimhihi_d,forwardimlohi_d,forwardimhilo_d,forwardimlolo_d,
          backwardrehihi_d,backwardrelohi_d,backwardrehilo_d,backwardrelolo_d,
          backwardimhihi_d,backwardimlohi_d,backwardimhilo_d,backwardimlolo_d,
          crossrehihi_d,crossrelohi_d,crossrehilo_d,crossrelolo_d,
          crossimhihi_d,crossimlohi_d,crossimhilo_d,crossimlolo_d);
   }
   double *forwardrehihi_h = new double[(deg1)*nvr];
   double *forwardrelohi_h = new double[(deg1)*nvr];
   double *forwardrehilo_h = new double[(deg1)*nvr];
   double *forwardrelolo_h = new double[(deg1)*nvr];
   double *forwardimhihi_h = new double[(deg1)*nvr];
   double *forwardimlohi_h = new double[(deg1)*nvr];
   double *forwardimhilo_h = new double[(deg1)*nvr];
   double *forwardimlolo_h = new double[(deg1)*nvr];
   double *backwardrehihi_h = new double[(deg1)*(nvr-2)];
   double *backwardrelohi_h = new double[(deg1)*(nvr-2)];
   double *backwardrehilo_h = new double[(deg1)*(nvr-2)];
   double *backwardrelolo_h = new double[(deg1)*(nvr-2)];
   double *backwardimhihi_h = new double[(deg1)*(nvr-2)];
   double *backwardimlohi_h = new double[(deg1)*(nvr-2)];
   double *backwardimhilo_h = new double[(deg1)*(nvr-2)];
   double *backwardimlolo_h = new double[(deg1)*(nvr-2)];
   double *crossrehihi_h = new double[(deg1)*(nvr-2)];
   double *crossrelohi_h = new double[(deg1)*(nvr-2)];
   double *crossrehilo_h = new double[(deg1)*(nvr-2)];
   double *crossrelolo_h = new double[(deg1)*(nvr-2)];
   double *crossimhihi_h = new double[(deg1)*(nvr-2)];
   double *crossimlohi_h = new double[(deg1)*(nvr-2)];
   double *crossimhilo_h = new double[(deg1)*(nvr-2)];
   double *crossimlolo_h = new double[(deg1)*(nvr-2)];
  
   cudaMemcpy(forwardrehihi_h,forwardrehihi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardrelohi_h,forwardrelohi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardrehilo_h,forwardrehilo_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardrelolo_h,forwardrelolo_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardimhihi_h,forwardimhihi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardimlohi_h,forwardimlohi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardimhilo_h,forwardimhilo_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardimlolo_h,forwardimlolo_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardrehihi_h,backwardrehihi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardrelohi_h,backwardrelohi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardrehilo_h,backwardrehilo_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardrelolo_h,backwardrelolo_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardimhihi_h,backwardimhihi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardimlohi_h,backwardimlohi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardimhilo_h,backwardimhilo_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardimlolo_h,backwardimlolo_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossrehihi_h,crossrehihi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossrelohi_h,crossrelohi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossrehilo_h,crossrehilo_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossrelolo_h,crossrelolo_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossimhihi_h,crossimhihi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossimlohi_h,crossimlohi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossimhilo_h,crossimhilo_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossimlolo_h,crossimlolo_d,sznvr2,cudaMemcpyDeviceToHost);

   int offset = (nvr-1)*deg1;
   for(int i=0; i<deg1; i++)   // assign value of the monomial
   {
      outputrehihi[dim][i] = forwardrehihi_h[offset+i];
      outputrelohi[dim][i] = forwardrelohi_h[offset+i];
      outputrehilo[dim][i] = forwardrehilo_h[offset+i];
      outputrelolo[dim][i] = forwardrelolo_h[offset+i];
      outputimhihi[dim][i] = forwardimhihi_h[offset+i];
      outputimlohi[dim][i] = forwardimlohi_h[offset+i];
      outputimhilo[dim][i] = forwardimhilo_h[offset+i];
      outputimlolo[dim][i] = forwardimlolo_h[offset+i];
   }
   ix = idx[nvr-1];
   offset = (nvr-2)*deg1;
   for(int i=0; i<deg1; i++)  // derivative with respect to x[n-1]
   {
      outputrehihi[ix][i] = forwardrehihi_h[offset+i];
      outputrelohi[ix][i] = forwardrelohi_h[offset+i];
      outputrehilo[ix][i] = forwardrehilo_h[offset+i];
      outputrelolo[ix][i] = forwardrelolo_h[offset+i];
      outputimhihi[ix][i] = forwardimhihi_h[offset+i];
      outputimlohi[ix][i] = forwardimlohi_h[offset+i];
      outputimhilo[ix][i] = forwardimhilo_h[offset+i];
      outputimlolo[ix][i] = forwardimlolo_h[offset+i];
   }
   ix = idx[0]; 
   offset = (nvr-3)*deg1;
   for(int i=0; i<deg1; i++)   // derivative with respect to x[0]
   {
      outputrehihi[ix][i] = backwardrehihi_h[offset+i];
      outputrelohi[ix][i] = backwardrelohi_h[offset+i];
      outputrehilo[ix][i] = backwardrehilo_h[offset+i];
      outputrelolo[ix][i] = backwardrelolo_h[offset+i];
      outputimhihi[ix][i] = backwardimhihi_h[offset+i];
      outputimlohi[ix][i] = backwardimlohi_h[offset+i];
      outputimhilo[ix][i] = backwardimhilo_h[offset+i];
      outputimlolo[ix][i] = backwardimlolo_h[offset+i];
   }
   for(int k=1; k<nvr-1; k++)  // derivative with respect to x[k]
   {
      ix = idx[k]; offset = (k-1)*deg1;
      for(int i=0; i<deg1; i++)
      {
         outputrehihi[ix][i] = crossrehihi_h[offset+i];
         outputrelohi[ix][i] = crossrelohi_h[offset+i];
         outputrehilo[ix][i] = crossrehilo_h[offset+i];
         outputrelolo[ix][i] = crossrelolo_h[offset+i];
         outputimhihi[ix][i] = crossimhihi_h[offset+i];
         outputimlohi[ix][i] = crossimlohi_h[offset+i];
         outputimhilo[ix][i] = crossimhilo_h[offset+i];
         outputimlolo[ix][i] = crossimlolo_h[offset+i];
      }
   }
}
