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
}
