// The file dbl4_convolutions_kernels.cu defines kernels with prototypes
// in dbl4_convolution_kernels.h.

#ifdef gpufun
#include "double_double_gpufun.cu"
#include "quad_double_gpufun.cu"
#endif
#include "dbl4_convolutions_kernels.h"

__global__ void dbl4_increment
 ( double *xhihi, double *xlohi, double *xhilo, double *xlolo,
   double *yhihi, double *ylohi, double *yhilo, double *ylolo,
   double *zhihi, double *zlohi, double *zhilo, double *zlolo, int dim )
{
   int k = threadIdx.x;                 // thread k computes z[k]

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

   xvhihi[k] = xhihi[k]; xvlohi[k] = xlohi[k]; 
   xvhilo[k] = xhilo[k]; xvlolo[k] = xlolo[k];
   yvhihi[k] = yhihi[k]; yvlohi[k] = ylohi[k]; 
   yvhilo[k] = yhilo[k]; yvlolo[k] = ylolo[k];

   qdg_add(xvhihi[k],xvlohi[k],xvhilo[k],xvlolo[k],
           yvhihi[k],yvlohi[k],yvhilo[k],yvlolo[k],
           &zvhihi[k],&zvlohi[k],&zvhilo[k],&zvlolo[k]);

   __syncthreads();

   zhihi[k] = zvhihi[k]; zlohi[k] = zvlohi[k];
   zhilo[k] = zvhilo[k]; zlolo[k] = zvlolo[k];
}

__global__ void dbl4_decrement
 ( double *xhihi, double *xlohi, double *xhilo, double *xlolo,
   double *yhihi, double *ylohi, double *yhilo, double *ylolo,
   double *zhihi, double *zlohi, double *zhilo, double *zlolo, int dim )
{
   int k = threadIdx.x;                 // thread k computes z[k]

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

   xvhihi[k] = xhihi[k]; xvlohi[k] = xlohi[k]; 
   xvhilo[k] = xhilo[k]; xvlolo[k] = xlolo[k];
   yvhihi[k] = yhihi[k]; yvlohi[k] = ylohi[k]; 
   yvhilo[k] = yhilo[k]; yvlolo[k] = ylolo[k];

   qdg_sub(xvhihi[k],xvlohi[k],xvhilo[k],xvlolo[k],
           yvhihi[k],yvlohi[k],yvhilo[k],yvlolo[k],
           &zvhihi[k],&zvlohi[k],&zvhilo[k],&zvlolo[k]);

   __syncthreads();

   zhihi[k] = zvhihi[k]; zlohi[k] = zvlohi[k];
   zhilo[k] = zvhilo[k]; zlolo[k] = zvlolo[k];
}

__global__ void dbl4_convolute
 ( double *xhihi, double *xlohi, double *xhilo, double *xlolo,
   double *yhihi, double *ylohi, double *yhilo, double *ylolo,
   double *zhihi, double *zlohi, double *zhilo, double *zlolo, int dim )
{
   int k = threadIdx.x;                 // thread k computes z[k]

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

   double prdhihi,prdlohi,prdhilo,prdlolo;

   xvhihi[k] = xhihi[k]; xvlohi[k] = xlohi[k]; 
   xvhilo[k] = xhilo[k]; xvlolo[k] = xlolo[k];
   yvhihi[k] = yhihi[k]; yvlohi[k] = ylohi[k]; 
   yvhilo[k] = yhilo[k]; yvlolo[k] = ylolo[k];

   __syncthreads();

   // zv[k] = xv[0]*yv[k];
   qdg_mul(xvhihi[0],xvlohi[0],xvhilo[0],xvlolo[0],
           yvhihi[k],yvlohi[k],yvhilo[k],yvlolo[k],
           &zvhihi[k],&zvlohi[k],&zvhilo[k],&zvlolo[k]);
   __syncthreads();

   for(int i=1; i<=k; i++) // zv[k] = zv[k] + xv[i]*yv[k-i];
   {
      qdg_mul(xvhihi[i],xvlohi[i],xvhilo[i],xvlolo[i],
              yvhihi[k-i],yvlohi[k-i],yvhilo[k-i],yvlolo[k-i],
              &prdhihi,&prdlohi,&prdhilo,&prdlolo);
      __syncthreads();
      qdg_inc(&zvhihi[k],&zvlohi[k],&zvhilo[k],&zvlolo[k],
              prdhihi,prdlohi,prdhilo,prdlolo);
      __syncthreads();
   }
   __syncthreads();

   zhihi[k] = zvhihi[k]; zlohi[k] = zvlohi[k];
   zhilo[k] = zvhilo[k]; zlolo[k] = zvlolo[k];
}

__global__ void dbl4_padded_convolute
 ( double *xhihi, double *xlohi, double *xhilo, double *xlolo,
   double *yhihi, double *ylohi, double *yhilo, double *ylolo,
   double *zhihi, double *zlohi, double *zhilo, double *zlolo, int dim )
{
   int k = threadIdx.x;                 // thread k computes z[k]

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

   double prdhihi,prdlohi,prdhilo,prdlolo;
   int idx = dim+k;

   xvhihi[k] = xhihi[k]; xvlohi[k] = xlohi[k]; 
   xvhilo[k] = xhilo[k]; xvlolo[k] = xlolo[k];
   yvhihi[k] = 0.0; yvlohi[k] = 0.0; 
   yvhilo[k] = 0.0; yvlolo[k] = 0.0;
   yvhihi[idx] = yhihi[k]; yvlohi[idx] = ylohi[k]; 
   yvhilo[idx] = yhilo[k]; yvlolo[idx] = ylolo[k];

   __syncthreads();

   // zv[k] = xv[0]*yv[k];
   qdg_mul(xvhihi[0],xvlohi[0],xvhilo[0],xvlolo[0],
           yvhihi[idx],yvlohi[idx],yvhilo[idx],yvlolo[idx],
           &zvhihi[k],&zvlohi[k],&zvhilo[k],&zvlolo[k]);
   __syncthreads();

   for(int i=1; i<dim; i++) // zv[k] = zv[k] + xv[i]*yv[k-i];
   {
      int idx = dim + k - i;
      qdg_mul(xvhihi[i],xvlohi[i],xvhilo[i],xvlolo[i],
              yvhihi[idx],yvlohi[idx],yvhilo[idx],yvlolo[idx],
              &prdhihi,&prdlohi,&prdhilo,&prdlolo);
      __syncthreads();
      qdg_inc(&zvhihi[k],&zvlohi[k],&zvhilo[k],&zvlolo[k],
              prdhihi,prdlohi,prdhilo,prdlolo);
      __syncthreads();
   }
   zhihi[k] = zvhihi[k]; zlohi[k] = zvlohi[k];
   zhilo[k] = zvhilo[k]; zlolo[k] = zvlolo[k];
}

__global__ void cmplx4_convolute
 ( double *xrehihi, double *xrelohi, double *xrehilo, double *xrelolo,
   double *ximhihi, double *ximlohi, double *ximhilo, double *ximlolo,
   double *yrehihi, double *yrelohi, double *yrehilo, double *yrelolo,
   double *yimhihi, double *yimlohi, double *yimhilo, double *yimlolo,
   double *zrehihi, double *zrelohi, double *zrehilo, double *zrelolo,
   double *zimhihi, double *zimlohi, double *zimhilo, double *zimlolo,
   int dim )
{
   int k = threadIdx.x;       // thread k computes zre[k] and zim[k]

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

   double xrhihi,xihihi,yrhihi,yihihi,zrhihi,zihihi,acchihi;
   double xrlohi,xilohi,yrlohi,yilohi,zrlohi,zilohi,acclohi;
   double xrhilo,xihilo,yrhilo,yihilo,zrhilo,zihilo,acchilo;
   double xrlolo,xilolo,yrlolo,yilolo,zrlolo,zilolo,acclolo;

   xvrehihi[k] = xrehihi[k]; xvimhihi[k] = ximhihi[k];
   xvrelohi[k] = xrelohi[k]; xvimlohi[k] = ximlohi[k];
   xvrehilo[k] = xrehilo[k]; xvimhilo[k] = ximhilo[k];
   xvrelolo[k] = xrelolo[k]; xvimlolo[k] = ximlolo[k];
   yvrehihi[k] = yrehihi[k]; yvimhihi[k] = yimhihi[k];
   yvrelohi[k] = yrelohi[k]; yvimlohi[k] = yimlohi[k];
   yvrehilo[k] = yrehilo[k]; yvimhilo[k] = yimhilo[k];
   yvrelolo[k] = yrelolo[k]; yvimlolo[k] = yimlolo[k];

   __syncthreads();

   // z[k] = x[0]*y[k]
   xrhihi = xvrehihi[0]; xrlohi = xvrelohi[0];
   xrhilo = xvrehilo[0]; xrlolo = xvrelolo[0];
   xihihi = xvimhihi[0]; xilohi = xvimlohi[0];
   xihilo = xvimhilo[0]; xilolo = xvimlolo[0];
   yrhihi = yvrehihi[k]; yrlohi = yvrelohi[k];
   yrhilo = yvrehilo[k]; yrlolo = yvrelolo[k];
   yihihi = yvimhihi[k]; yilohi = yvimlohi[k];
   yihilo = yvimhilo[k]; yilolo = yvimlolo[k];

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

   zvrehihi[k] = zrhihi; zvrelohi[k] = zrlohi;
   zvrehilo[k] = zrhilo; zvrelolo[k] = zrlolo;
   zvimhihi[k] = zihihi; zvimlohi[k] = zilohi;
   zvimhilo[k] = zihilo; zvimlolo[k] = zilolo;

   for(int i=1; i<=k; i++) // z[k] = z[k] + x[i]*y[k-i]
   {
      xrhihi = xvrehihi[i]; xrlohi = xvrelohi[i];
      xrhilo = xvrehilo[i]; xrlolo = xvrelolo[i];
      xihihi = xvimhihi[i]; xilohi = xvimlohi[i];
      xihilo = xvimhilo[i]; xilolo = xvimlolo[i];
      yrhihi = yvrehihi[k-i]; yrlohi = yvrelohi[k-i];
      yrhilo = yvrehilo[k-i]; yrlolo = yvrelolo[k-i];
      yihihi = yvimhihi[k-i]; yilohi = yvimlohi[k-i];
      yihilo = yvimhilo[k-i]; yilolo = yvimlolo[k-i];

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

      qdg_inc(&zvrehihi[k],&zvrelohi[k],&zvrehilo[k],&zvrelolo[k],
              zrhihi,zrlohi,zrhilo,zrlolo);           // zvre[k] += zr;
      qdg_inc(&zvimhihi[k],&zvimlohi[k],&zvimhilo[k],&zvimlolo[k],
              zihihi,zilohi,zihilo,zilolo);           // zvim[k] += zi;
   }

   __syncthreads();

   zrehihi[k] = zvrehihi[k]; zrelohi[k] = zvrelohi[k];
   zrehilo[k] = zvrehilo[k]; zrelolo[k] = zvrelolo[k];
   zimhihi[k] = zvimhihi[k]; zimlohi[k] = zvimlohi[k];
   zimhilo[k] = zvimhilo[k]; zimlolo[k] = zvimlolo[k];
}

__global__ void cmplx4_padded_convolute
 ( double *xrehihi, double *xrelohi, double *xrehilo, double *xrelolo,
   double *ximhihi, double *ximlohi, double *ximhilo, double *ximlolo,
   double *yrehihi, double *yrelohi, double *yrehilo, double *yrelolo,
   double *yimhihi, double *yimlohi, double *yimhilo, double *yimlolo,
   double *zrehihi, double *zrelohi, double *zrehilo, double *zrelolo,
   double *zimhihi, double *zimlohi, double *zimhilo, double *zimlolo,
   int dim )
{
   int k = threadIdx.x;       // thread k computes zre[k] and zim[k]

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

   double xrhihi,xihihi,yrhihi,yihihi,zrhihi,zihihi,acchihi;
   double xrlohi,xilohi,yrlohi,yilohi,zrlohi,zilohi,acclohi;
   double xrhilo,xihilo,yrhilo,yihilo,zrhilo,zihilo,acchilo;
   double xrlolo,xilolo,yrlolo,yilolo,zrlolo,zilolo,acclolo;
   int idx = dim+k;

   xvrehihi[k] = xrehihi[k]; xvimhihi[k] = ximhihi[k];
   xvrelohi[k] = xrelohi[k]; xvimlohi[k] = ximlohi[k];
   xvrehilo[k] = xrehilo[k]; xvimhilo[k] = ximhilo[k];
   xvrelolo[k] = xrelolo[k]; xvimlolo[k] = ximlolo[k];
   yvrehihi[k] = 0.0; yvimhihi[k] = 0.0;
   yvrelohi[k] = 0.0; yvimlohi[k] = 0.0;
   yvrehilo[k] = 0.0; yvimhilo[k] = 0.0;
   yvrelolo[k] = 0.0; yvimlolo[k] = 0.0;
   yvrehihi[idx] = yrehihi[k]; yvimhihi[idx] = yimhihi[k];
   yvrelohi[idx] = yrelohi[k]; yvimlohi[idx] = yimlohi[k];
   yvrehilo[idx] = yrehilo[k]; yvimhilo[idx] = yimhilo[k];
   yvrelolo[idx] = yrelolo[k]; yvimlolo[idx] = yimlolo[k];

   __syncthreads();

   // z[k] = x[0]*y[k]
   xrhihi = xvrehihi[0]; xrlohi = xvrelohi[0];
   xrhilo = xvrehilo[0]; xrlolo = xvrelolo[0];
   xihihi = xvimhihi[0]; xilohi = xvimlohi[0];
   xihilo = xvimhilo[0]; xilolo = xvimlolo[0];
   yrhihi = yvrehihi[idx]; yrlohi = yvrelohi[idx];
   yrhilo = yvrehilo[idx]; yrlolo = yvrelolo[idx];
   yihihi = yvimhihi[idx]; yilohi = yvimlohi[idx];
   yihilo = yvimhilo[idx]; yilolo = yvimlolo[idx];

   __syncthreads();

   qdg_mul(xrhihi,xrlohi,xrhilo,xrlolo,
           yrhihi,yrlohi,yrhilo,yrlolo,
           &zrhihi,&zrlohi,&zrhilo,&zrlolo);       // zr = xr*yr
   __syncthreads();
   qdg_mul(xihihi,xilohi,xihilo,xilolo,
           yihihi,yilohi,yihilo,yilolo,
           &acchihi,&acclohi,&acchilo,&acclolo);   // acc = xi*yi
   __syncthreads();
   qdg_minus(&acchihi,&acclohi,&acchilo,&acclolo);
   qdg_inc(&zrhihi,&zrlohi,&zrhilo,&zrlolo,
           acchihi,acclohi,acchilo,acclolo);       // zr = xr*yr - xi*yi
   __syncthreads();
   qdg_mul(xrhihi,xrlohi,xrhilo,xrlolo,
           yihihi,yilohi,yihilo,yilolo,
           &zihihi,&zilohi,&zihilo,&zilolo);       // zi = xr*yi
   __syncthreads();
   qdg_mul(xihihi,xilohi,xihilo,xilolo,
           yrhihi,yrlohi,yrhilo,yrlolo,
           &acchihi,&acclohi,&acchilo,&acclolo);   // acc = xi*yr
   __syncthreads();
   qdg_inc(&zihihi,&zilohi,&zihilo,&zilolo,
           acchihi,acclohi,acchilo,acclolo);       // zr = xr*yr + xi*yi
   __syncthreads();

   zvrehihi[k] = zrhihi; zvrelohi[k] = zrlohi;
   zvrehilo[k] = zrhilo; zvrelolo[k] = zrlolo;
   zvimhihi[k] = zihihi; zvimlohi[k] = zilohi;
   zvimhilo[k] = zihilo; zvimlolo[k] = zilolo;

   __syncthreads();

   for(int i=1; i<dim; i++) // z[k] = z[k] + x[i]*y[k-i]
   {
      idx = dim + k - i;
      xrhihi = xvrehihi[i]; xrlohi = xvrelohi[i];
      xrhilo = xvrehilo[i]; xrlolo = xvrelolo[i];
      xihihi = xvimhihi[i]; xilohi = xvimlohi[i];
      xihilo = xvimhilo[i]; xilolo = xvimlolo[i];
      yrhihi = yvrehihi[idx]; yrlohi = yvrelohi[idx];
      yrhilo = yvrehilo[idx]; yrlolo = yvrelolo[idx];
      yihihi = yvimhihi[idx]; yilohi = yvimlohi[idx];
      yihilo = yvimhilo[idx]; yilolo = yvimlolo[idx];

      __syncthreads();

      qdg_mul(xrhihi,xrlohi,xrhilo,xrlolo,
              yrhihi,yrlohi,yrhilo,yrlolo,
              &zrhihi,&zrlohi,&zrhilo,&zrlolo);       // zr = xr*yr
      __syncthreads();
      qdg_mul(xihihi,xilohi,xihilo,xilolo,
              yihihi,yilohi,yihilo,yilolo,
              &acchihi,&acclohi,&acchilo,&acclolo);   // acc = xi*yi
      __syncthreads();
      qdg_minus(&acchihi,&acclohi,&acchilo,&acclolo);
      qdg_inc(&zrhihi,&zrlohi,&zrhilo,&zrlolo,
              acchihi,acclohi,acchilo,acclolo);       // zr = xr*yr - xi*yi
      __syncthreads();
      qdg_mul(xrhihi,xrlohi,xrhilo,xrlolo,
              yihihi,yilohi,yihilo,yilolo,
              &zihihi,&zilohi,&zihilo,&zilolo);       // zi = xr*yi
      __syncthreads();
      qdg_mul(xihihi,xilohi,xihilo,xilolo,
              yrhihi,yrlohi,yrhilo,yrlolo,
              &acchihi,&acclohi,&acchilo,&acclolo);   // acc = xi*yr
      __syncthreads();
      qdg_inc(&zihihi,&zilohi,&zihilo,&zilolo,
              acchihi,acclohi,acchilo,acclolo);       // zr = xr*yr + xi*yi
      __syncthreads();
      qdg_inc(&zvrehihi[k],&zvrelohi[k],&zvrehilo[k],&zvrelolo[k],
              zrhihi,zrlohi,zrhilo,zrlolo);           // zvre[k] += zr;
      __syncthreads();
      qdg_inc(&zvimhihi[k],&zvimlohi[k],&zvimhilo[k],&zvimlolo[k],
              zihihi,zilohi,zihilo,zilolo);           // zvim[k] += zi;
      __syncthreads();
   }
   zrehihi[k] = zvrehihi[k]; zrelohi[k] = zvrelohi[k];
   zrehilo[k] = zvrehilo[k]; zrelolo[k] = zvrelolo[k];
   zimhihi[k] = zvimhihi[k]; zimlohi[k] = zvimlohi[k];
   zimhilo[k] = zvimhilo[k]; zimlolo[k] = zvimlolo[k];
   __syncthreads();
}

void GPU_dbl4_product
 ( double *xhihi_h, double *xlohi_h, double *xhilo_h, double *xlolo_h,
   double *yhihi_h, double *ylohi_h, double *yhilo_h, double *ylolo_h,
   double *zhihi_h, double *zlohi_h, double *zhilo_h, double *zlolo_h,
   int deg, int freq, int BS, int padded )
{
   const int dim = deg+1;            // length of all vectors
   double* xhihi_d;                  // xhihi_d is xhihi_h on the device
   double* xlohi_d;                  // xlohi_d is xlohi_h on the device
   double* xhilo_d;                  // xhilo_d is xhilo_h on the device
   double* xlolo_d;                  // xlolo_d is xlolo_h on the device
   double* yhihi_d;                  // yhihi_d is yhihi_h on the device
   double* ylohi_d;                  // ylohi_d is ylohi_h on the device
   double* yhilo_d;                  // yhilo_d is yhilo_h on the device
   double* ylolo_d;                  // ylolo_d is ylolo_h on the device
   double* zhihi_d;                  // zhihi_d is zhihi_h on the device
   double* zlohi_d;                  // zlohi_d is zlohi_h on the device
   double* zhilo_d;                  // zhilo_d is zhilo_h on the device
   double* zlolo_d;                  // zlolo_d is zlolo_h on the device
   size_t size = dim*sizeof(double); // number of bytes for each vector

   cudaMalloc((void**)&xhihi_d,size);
   cudaMalloc((void**)&xlohi_d,size);
   cudaMalloc((void**)&xhilo_d,size);
   cudaMalloc((void**)&xlolo_d,size);
   cudaMalloc((void**)&yhihi_d,size);
   cudaMalloc((void**)&ylohi_d,size);
   cudaMalloc((void**)&yhilo_d,size);
   cudaMalloc((void**)&ylolo_d,size);
   cudaMalloc((void**)&zhihi_d,size);
   cudaMalloc((void**)&zlohi_d,size);
   cudaMalloc((void**)&zhilo_d,size);
   cudaMalloc((void**)&zlolo_d,size);
   cudaMemcpy(xhihi_d,xhihi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xlohi_d,xlohi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xhilo_d,xhilo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xlolo_d,xlolo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yhihi_d,yhihi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(ylohi_d,ylohi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yhilo_d,yhilo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(ylolo_d,ylolo_h,size,cudaMemcpyHostToDevice);

   if(dim == BS)
   {
      if(padded == 1)
      {
         for(int i=0; i<freq; i++)
            dbl4_padded_convolute<<<1,BS>>>
               (xhihi_d,xlohi_d,xhilo_d,xlolo_d,
                yhihi_d,ylohi_d,yhilo_d,ylolo_d,
                zhihi_d,zlohi_d,zhilo_d,zlolo_d,dim);
      }
      else
      {
         for(int i=0; i<freq; i++)
            dbl4_convolute<<<1,BS>>>
               (xhihi_d,xlohi_d,xhilo_d,xlolo_d,
                yhihi_d,ylohi_d,yhilo_d,ylolo_d,
                zhihi_d,zlohi_d,zhilo_d,zlolo_d,dim);
      }
   }

   cudaMemcpy(zhihi_h,zhihi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zlohi_h,zlohi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zhilo_h,zhilo_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zlolo_h,zlolo_d,size,cudaMemcpyDeviceToHost);
}

void GPU_cmplx4_product
 ( double *xrehihi_h, double *xrelohi_h, double *xrehilo_h, double *xrelolo_h,
   double *ximhihi_h, double *ximlohi_h, double *ximhilo_h, double *ximlolo_h,
   double *yrehihi_h, double *yrelohi_h, double *yrehilo_h, double *yrelolo_h,
   double *yimhihi_h, double *yimlohi_h, double *yimhilo_h, double *yimlolo_h,
   double *zrehihi_h, double *zrelohi_h, double *zrehilo_h, double *zrelolo_h,
   double *zimhihi_h, double *zimlohi_h, double *zimhilo_h, double *zimlolo_h,
   int deg, int freq, int BS, int mode )
{
   const int dim = deg+1;            // length of all vectors
   double* xrehihi_d;                // xrehihi_d is xrehihi_h on the device
   double* xrelohi_d;                // xrelohi_d is xrelohi_h on the device
   double* xrehilo_d;                // xrehilo_d is xrehilo_h on the device
   double* xrelolo_d;                // xrelolo_d is xrelolo_h on the device
   double* ximhihi_d;                // ximhihi_d is ximhihi_h on the device
   double* ximlohi_d;                // ximlohi_d is ximlohi_h on the device
   double* ximhilo_d;                // ximhilo_d is ximhilo_h on the device
   double* ximlolo_d;                // ximlolo_d is ximlolo_h on the device
   double* yrehihi_d;                // yrehihi_d is yrehihi_h on the device
   double* yrelohi_d;                // yrelohi_d is yrelohi_h on the device
   double* yrehilo_d;                // yrehilo_d is yrehilo_h on the device
   double* yrelolo_d;                // yrelolo_d is yrelolo_h on the device
   double* yimhihi_d;                // yimhihi_d is yimhihi_h on the device
   double* yimlohi_d;                // yimlohi_d is yimlohi_h on the device
   double* yimhilo_d;                // yimhilo_d is yimhilo_h on the device
   double* yimlolo_d;                // yimlolo_d is yimlolo_h on the device
   double* zrehihi_d;                // zrehihi_d is zrehihi_h on the device
   double* zrelohi_d;                // zrelohi_d is zrelohi_h on the device
   double* zrehilo_d;                // zrehilo_d is zrehilo_h on the device
   double* zrelolo_d;                // zrelolo_d is zrelolo_h on the device
   double* zimhihi_d;                // zimhihi_d is zimhihi_h on the device
   double* zimlohi_d;                // zimlohi_d is zimlohi_h on the device
   double* zimhilo_d;                // zimhilo_d is zimhilo_h on the device
   double* zimlolo_d;                // zimlolo_d is zimlolo_h on the device
   double* acchihi_d;                // accumulates highest doubles
   double* acclohi_d;                // accumulates second highest doubles
   double* acchilo_d;                // accumulates second lowest doubles
   double* acclolo_d;                // accumulates lowest doubles
   size_t size = dim*sizeof(double); // number of bytes for each vector

   cudaMalloc((void**)&xrehihi_d,size);
   cudaMalloc((void**)&xrelohi_d,size);
   cudaMalloc((void**)&xrehilo_d,size);
   cudaMalloc((void**)&xrelolo_d,size);
   cudaMalloc((void**)&ximhihi_d,size);
   cudaMalloc((void**)&ximlohi_d,size);
   cudaMalloc((void**)&ximhilo_d,size);
   cudaMalloc((void**)&ximlolo_d,size);
   cudaMalloc((void**)&yrehihi_d,size);
   cudaMalloc((void**)&yrelohi_d,size);
   cudaMalloc((void**)&yrehilo_d,size);
   cudaMalloc((void**)&yrelolo_d,size);
   cudaMalloc((void**)&yimhihi_d,size);
   cudaMalloc((void**)&yimlohi_d,size);
   cudaMalloc((void**)&yimhilo_d,size);
   cudaMalloc((void**)&yimlolo_d,size);
   cudaMalloc((void**)&zrehihi_d,size);
   cudaMalloc((void**)&zrelohi_d,size);
   cudaMalloc((void**)&zrehilo_d,size);
   cudaMalloc((void**)&zrelolo_d,size);
   cudaMalloc((void**)&zimhihi_d,size);
   cudaMalloc((void**)&zimlohi_d,size);
   cudaMalloc((void**)&zimhilo_d,size);
   cudaMalloc((void**)&zimlolo_d,size);
   cudaMalloc((void**)&acchihi_d,size);
   cudaMalloc((void**)&acclohi_d,size);
   cudaMalloc((void**)&acchilo_d,size);
   cudaMalloc((void**)&acclolo_d,size);
   cudaMemcpy(xrehihi_d,xrehihi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xrelohi_d,xrelohi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xrehilo_d,xrehilo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xrelolo_d,xrelolo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(ximhihi_d,ximhihi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(ximlohi_d,ximlohi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(ximhilo_d,ximhilo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(ximlolo_d,ximlolo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yrehihi_d,yrehihi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yrelohi_d,yrelohi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yrehilo_d,yrehilo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yrelolo_d,yrelolo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yimhihi_d,yimhihi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yimlohi_d,yimlohi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yimhilo_d,yimhilo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yimlolo_d,yimlolo_h,size,cudaMemcpyHostToDevice);

   if(dim == BS)
   {
      if(mode == 2)
      { 
         for(int i=0; i<freq; i++)
            cmplx4_padded_convolute<<<1,BS>>>
               (xrehihi_d,xrelohi_d,xrehilo_d,xrelolo_d,
                ximhihi_d,ximlohi_d,ximhilo_d,ximlolo_d,
                yrehihi_d,yrelohi_d,yrehilo_d,yrelolo_d,
                yimhihi_d,yimlohi_d,yimhilo_d,yimlolo_d,
                zrehihi_d,zrelohi_d,zrehilo_d,zrelolo_d,
                zimhihi_d,zimlohi_d,zimhilo_d,zimlolo_d,dim);
      }
      else if(mode == 1)
      {
         for(int i=0; i<freq; i++)
         {
            dbl4_padded_convolute<<<1,BS>>>
               (xrehihi_d,xrelohi_d,xrehilo_d,xrelolo_d,
                yrehihi_d,yrelohi_d,yrehilo_d,yrelolo_d,
                zrehihi_d,zrelohi_d,zrehilo_d,zrelolo_d,dim);
            dbl4_padded_convolute<<<1,BS>>>
               (ximhihi_d,ximlohi_d,ximhilo_d,ximlolo_d,
                yimhihi_d,yimlohi_d,yimhilo_d,yimlolo_d,
                acchihi_d,acclohi_d,acchilo_d,acclolo_d,dim);
            dbl4_decrement<<<1,BS>>>
               (zrehihi_d,zrelohi_d,zrehilo_d,zrelolo_d,
                acchihi_d,acclohi_d,acchilo_d,acclolo_d,
                zrehihi_d,zrelohi_d,zrehilo_d,zrelolo_d,dim);
            dbl4_padded_convolute<<<1,BS>>>
               (xrehihi_d,xrelohi_d,xrehilo_d,xrelolo_d,
                yimhihi_d,yimlohi_d,yimhilo_d,yimlolo_d,
                zimhihi_d,zimlohi_d,zimhilo_d,zimlolo_d,dim);
            dbl4_padded_convolute<<<1,BS>>>
               (ximhihi_d,ximlohi_d,ximhilo_d,ximlolo_d,
                yrehihi_d,yrelohi_d,yrehilo_d,yrelolo_d,
                acchihi_d,acclohi_d,acchilo_d,acclolo_d,dim);
            dbl4_increment<<<1,BS>>>
               (zimhihi_d,zimlohi_d,zimhilo_d,zimlolo_d,
                acchihi_d,acclohi_d,acchilo_d,acclolo_d,
                zimhihi_d,zimlohi_d,zimhilo_d,zimlolo_d,dim);
         }
      }
      else
      { 
         for(int i=0; i<freq; i++)
            cmplx4_convolute<<<1,BS>>>
               (xrehihi_d,xrelohi_d,xrehilo_d,xrelolo_d,
                ximhihi_d,ximlohi_d,ximhilo_d,ximlolo_d,
                yrehihi_d,yrelohi_d,yrehilo_d,yrelolo_d,
                yimhihi_d,yimlohi_d,yimhilo_d,yimlolo_d,
                zrehihi_d,zrelohi_d,zrehilo_d,zrelolo_d,
                zimhihi_d,zimlohi_d,zimhilo_d,zimlolo_d,dim);
      }
   }

   cudaMemcpy(zrehihi_h,zrehihi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zrelohi_h,zrelohi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zrehilo_h,zrehilo_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zrelolo_h,zrelolo_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zimhihi_h,zimhihi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zimlohi_h,zimlohi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zimhilo_h,zimhilo_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zimlolo_h,zimlolo_d,size,cudaMemcpyDeviceToHost);
}
