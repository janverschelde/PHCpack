// The file dbl4_polynomials_kernels.cu defines the kernels with prototypes
// in dbl4_polynomials_kernels.h.

#include <iostream>
#include <iomanip>
#ifdef winwalltime
#include "gettimeofday4win.h"
#else
#include <sys/time.h>
#endif
#include "job_coordinates.h"
#include "quad_double_functions.h"
#ifdef gpufun
#include "double_double_gpufun.cu"
#include "quad_double_gpufun.cu"
#endif
#include "dbl4_polynomials_kernels.h"
#include "write_gpu_timings.h"

// The constant qd_shmemsize is the bound on the shared memory size.

#define qd_shmemsize 192

using namespace std;

__global__ void dbl4_padded_convjobs
 ( double *datahihi, double *datalohi, double *datahilo, double *datalolo,
   int *in1idx, int *in2idx, int *outidx, int dim )
{
   const int bdx = blockIdx.x;           // index to the convolution job
   const int tdx = threadIdx.x;          // index to the output of the job
   const int idx1 = in1idx[bdx] + tdx;
   const int idx2 = in2idx[bdx] + tdx;
   const int idx3 = outidx[bdx] + tdx;

   __shared__ double xvhihi[qd_shmemsize];
   __shared__ double xvlohi[qd_shmemsize];
   __shared__ double xvhilo[qd_shmemsize];
   __shared__ double xvlolo[qd_shmemsize];
   __shared__ double yvhihi[2*qd_shmemsize];
   __shared__ double yvlohi[2*qd_shmemsize];
   __shared__ double yvhilo[2*qd_shmemsize];
   __shared__ double yvlolo[2*qd_shmemsize];
   __shared__ double zvhihi[qd_shmemsize];
   __shared__ double zvlohi[qd_shmemsize];
   __shared__ double zvhilo[qd_shmemsize];
   __shared__ double zvlolo[qd_shmemsize];

   double prdhihi,prdlohi,prdhilo,prdlolo;
   int ydx = dim + tdx;

   xvhihi[tdx] = datahihi[idx1];  // loading first input
   xvlohi[tdx] = datalohi[idx1]; 
   xvhilo[tdx] = datahilo[idx1]; 
   xvlolo[tdx] = datalolo[idx1]; 
   yvhihi[tdx] = 0.0;             // padded with zeros
   yvlohi[tdx] = 0.0;
   yvhilo[tdx] = 0.0;
   yvlolo[tdx] = 0.0;
   yvhihi[ydx] = datahihi[idx2];  // loading second input
   yvlohi[ydx] = datalohi[idx2];
   yvhilo[ydx] = datahilo[idx2];
   yvlolo[ydx] = datalolo[idx2];

   __syncthreads();

   // zv[tdx] = xv[0]*yv[tdx];
   qdg_mul( xvhihi[0],   xvlohi[0],   xvhilo[0],   xvlolo[0],
            yvhihi[ydx], yvlohi[ydx], yvhilo[ydx], yvlolo[ydx],
           &zvhihi[tdx],&zvlohi[tdx],&zvhilo[tdx],&zvlolo[tdx]);
   __syncthreads();

   for(int i=1; i<dim; i++) // zv[tdx] = zv[tdx] + xv[i]*yv[dim+tdx-i];
   {
      ydx = dim + tdx - i;

      qdg_mul( xvhihi[i],  xvlohi[i],  xvhilo[i],  xvlolo[i],
               yvhihi[ydx],yvlohi[ydx],yvhilo[ydx],yvlolo[ydx],
              &prdhihi,  &prdlohi,   &prdhilo,   &prdlolo);
      __syncthreads();

      qdg_inc(&zvhihi[tdx],&zvlohi[tdx],&zvhilo[tdx],&zvlolo[tdx],
              prdhihi,     prdlohi,     prdhilo,     prdlolo);
      __syncthreads();
   }
   __syncthreads();

   datahihi[idx3] = zvhihi[tdx]; // storing the output
   datalohi[idx3] = zvlohi[tdx];
   datahilo[idx3] = zvhilo[tdx];
   datalolo[idx3] = zvlolo[tdx];
}

__global__ void dbl4_increment_jobs
 ( double *datahihi, double *datalohi, double *datahilo, double *datalolo,
   int *in1idx, int *in2idx, int *outidx, int dim )
{
   const int bdx = blockIdx.x;           // index to the increment job
   const int tdx = threadIdx.x;          // index to the output of the job
   const int idx1 = in1idx[bdx] + tdx;
   const int idx2 = in2idx[bdx] + tdx;
   const int idx3 = outidx[bdx] + tdx;

   __shared__ double yvhihi[qd_shmemsize];
   __shared__ double yvlohi[qd_shmemsize];
   __shared__ double yvhilo[qd_shmemsize];
   __shared__ double yvlolo[qd_shmemsize];
   __shared__ double zvhihi[qd_shmemsize];
   __shared__ double zvlohi[qd_shmemsize];
   __shared__ double zvhilo[qd_shmemsize];
   __shared__ double zvlolo[qd_shmemsize];

   zvhihi[tdx] = datahihi[idx1];  // loading first input
   zvlohi[tdx] = datalohi[idx1]; 
   zvhilo[tdx] = datahilo[idx1]; 
   zvlolo[tdx] = datalolo[idx1]; 
   yvhihi[tdx] = datahihi[idx2];  // loading second input
   yvlohi[tdx] = datalohi[idx2];
   yvhilo[tdx] = datahilo[idx2];
   yvlolo[tdx] = datalolo[idx2];

   __syncthreads();

   qdg_inc(&zvhihi[tdx],&zvlohi[tdx],&zvhilo[tdx],&zvlolo[tdx],
            yvhihi[tdx], yvlohi[tdx], yvhilo[tdx], yvlolo[tdx]);

   __syncthreads();

   datahihi[idx3] = zvhihi[tdx]; // storing the output
   datalohi[idx3] = zvlohi[tdx];
   datahilo[idx3] = zvhilo[tdx];
   datalolo[idx3] = zvlolo[tdx];
}

/*
__global__ void cmplx4_padded_convjobs
 ( double *datarehihi, double *datarelohi,
   double *datarehilo, double *datarelolo,
   double *dataimhihi, double *dataimlohi,
   double *dataimhilo, double *dataimlolo,
   int *in1idx, int *in2idx, int *outidx, int dim )
{
   const int bdx = blockIdx.x;           // index to the convolution job
   const int tdx = threadIdx.x;          // index to the output of the job
   const int idx1 = in1idx[bdx] + tdx;
   const int idx2 = in2idx[bdx] + tdx;
   const int idx3 = outidx[bdx] + tdx;

   __shared__ double xvrehihi[qd_shmemsize];
   __shared__ double xvrelohi[qd_shmemsize];
   __shared__ double xvrehilo[qd_shmemsize];
   __shared__ double xvrelolo[qd_shmemsize];
   __shared__ double xvimhihi[qd_shmemsize];
   __shared__ double xvimlohi[qd_shmemsize];
   __shared__ double xvimhilo[qd_shmemsize];
   __shared__ double xvimlolo[qd_shmemsize];
   __shared__ double yvrehihi[2*qd_shmemsize];
   __shared__ double yvrelohi[2*qd_shmemsize];
   __shared__ double yvrehilo[2*qd_shmemsize];
   __shared__ double yvrelolo[2*qd_shmemsize];
   __shared__ double yvimhihi[2*qd_shmemsize];
   __shared__ double yvimlohi[2*qd_shmemsize];
   __shared__ double yvimhilo[2*qd_shmemsize];
   __shared__ double yvimlolo[2*qd_shmemsize];
   __shared__ double zvrehihi[qd_shmemsize];
   __shared__ double zvrelohi[qd_shmemsize];
   __shared__ double zvrehilo[qd_shmemsize];
   __shared__ double zvrelolo[qd_shmemsize];
   __shared__ double zvimhihi[qd_shmemsize];
   __shared__ double zvimlohi[qd_shmemsize];
   __shared__ double zvimhilo[qd_shmemsize];
   __shared__ double zvimlolo[qd_shmemsize];

   double prodhihi,prodlohi,prodhilo,prodlolo;
   int ydx = dim + tdx;

   xvrehihi[tdx] = datarehihi[idx1];  // loading first input
   xvrelohi[tdx] = datarelohi[idx1]; 
   xvrehilo[tdx] = datarehilo[idx1]; 
   xvrelolo[tdx] = datarelolo[idx1]; 
   xvimhihi[tdx] = dataimhihi[idx1];
   xvimlohi[tdx] = dataimlohi[idx1];
   xvimhilo[tdx] = dataimhilo[idx1]; 
   xvimlolo[tdx] = dataimlolo[idx1]; 
   yvrehihi[tdx] = 0.0;           // padded with zeros
   yvrelohi[tdx] = 0.0;
   yvrehilo[tdx] = 0.0;
   yvrelolo[tdx] = 0.0;
   yvimhihi[tdx] = 0.0;
   yvimlohi[tdx] = 0.0;
   yvimhilo[tdx] = 0.0;
   yvimlolo[tdx] = 0.0;
   yvrehihi[ydx] = datarehihi[idx2];  // loading second input
   yvrelohi[ydx] = datarelohi[idx2];
   yvrehilo[ydx] = datarehilo[idx2];
   yvrelolo[ydx] = datarelolo[idx2];
   yvimhihi[ydx] = dataimhihi[idx2];
   yvimlohi[ydx] = dataimlohi[idx2];
   yvimhilo[ydx] = dataimhilo[idx2];
   yvimlolo[ydx] = dataimlolo[idx2];

   __syncthreads();

   // zv[tdx] = xv[0]*yv[tdx];
   qdg_mul( xvrehihi[0],   xvrelohi[0],   xvrehilo[0],   xvrelolo[0],
            yvrehihi[ydx], yvrelohi[ydx], yvrehilo[ydx], yvrelolo[ydx],
           &zvrehihi[tdx],&zvrelohi[tdx],&zvrehilo[tdx],&zvrelolo[tdx]);
   __syncthreads();
   qdg_mul( xvimhihi[0],  xvimlohi[0],  xvimhilo[0],  xvimlolo[0],
            yvimhihi[ydx],yvimlohi[ydx],yvimhilo[ydx],yvimlolo[ydx],
           &prodhihi,    &prodlohi,    &prodhilo,    &prodlolo);
   __syncthreads();
   qdg_minus(&prodhihi,&prodlohi,&prodhilo,&prodlolo);
   qdg_inc(&zvrehihi[tdx],&zvrelohi[tdx],&zvrehilo[tdx],&zvrelolo[tdx],
             prodhihi,     prodlohi,      prodhilo,      prodlolo);
   __syncthreads();

   qdg_mul( xvrehihi[0],   xvrelohi[0],   xvrehilo[0],   xvrelolo[0],
            yvimhihi[ydx], yvimlohi[ydx], yvimhilo[ydx], yvimlolo[ydx],
           &zvimhihi[tdx],&zvimlohi[tdx],&zvimhilo[tdx],&zvimlolo[tdx]);
   __syncthreads();
   qdg_mul( xvimhihi[0],  xvimlohi[0],  xvimhilo[0],  xvimlolo[0],
            yvrehihi[ydx],yvrelohi[ydx],yvrehilo[ydx],yvrelolo[ydx],
           &prodhihi,    &prodlohi,    &prodhilo,    &prodlolo);
   __syncthreads();
   qdg_inc(&zvimhihi[tdx],&zvimlohi[tdx],&zvimhilo[tdx],&zvimlolo[tdx],
            prodhihi,      prodlohi,      prodhilo,      prodlolo);
   __syncthreads();

   for(int i=1; i<dim; i++) // zv[tdx] = zv[tdx] + xv[i]*yv[dim+tdx-i];
   {
      ydx = dim + tdx - i;

      qdg_mul( xvrehihi[i],  xvrelohi[i],  xvrehilo[i],  xvrelolo[i],
               yvrehihi[ydx],yvrelohi[ydx],yvrehilo[ydx],yvrelolo[ydx],
              &prodhihi,    &prodlohi,    &prodhilo,    &prodlolo);
      __syncthreads();
      qdg_inc(&zvrehihi[tdx],&zvrelohi[tdx],&zvrehilo[tdx],&zvrelolo[tdx],
               prodhihi,      prodlohi,      prodhilo,      prodlolo);
      __syncthreads();
      qdg_mul( xvimhihi[i],  xvimlohi[i],  xvimhilo[i],  xvimlolo[i],
               yvimhihi[ydx],yvimlohi[ydx],yvimhilo[ydx],yvimlolo[ydx],
              &prodhihi,    &prodlohi,    &prodhilo,    &prodlolo);
      __syncthreads();
      qdg_minus(&prodhihi,&prodlohi,&prodhilo,&prodlolo);
      qdg_inc(&zvrehihi[tdx],&zvrelohi[tdx],&zvrehilo[tdx],&zvrelolo[tdx],
               prodhihi,      prodlohi,      prodhilo,      prodlolo);
      __syncthreads();

      qdg_mul( xvrehihi[i],  xvrelohi[i],  xvrehilo[i],  xvrelolo[i],
               yvimhihi[ydx],yvimlohi[ydx],yvimhilo[ydx],yvimlolo[ydx],
              &prodhihi,    &prodlohi,    &prodhilo,    &prodlolo);
      __syncthreads();
      qdg_inc(&zvimhihi[tdx],&zvimlohi[tdx],&zvimhilo[tdx],&zvimlolo[tdx],
                prodhihi,     prodlohi,      prodhilo,      prodlolo);
      __syncthreads();
      qdg_mul( xvimhihi[i],  xvimlohi[i],  xvimhilo[i],  xvimlolo[i],
               yvrehihi[ydx],yvrelohi[ydx],yvrehilo[ydx],yvrelolo[ydx],
              &prodhihi,    &prodlohi,    &prodhilo,    &prodlolo);
      __syncthreads();
      qdg_inc(&zvimhihi[tdx],&zvimlohi[tdx],&zvimhilo[tdx],&zvimlolo[tdx],
               prodhihi,      prodlohi,      prodhilo,      prodlolo);
      __syncthreads();
   }
   __syncthreads();

   datarehihi[idx3] = zvrehihi[tdx]; // storing the output
   datarelohi[idx3] = zvrelohi[tdx];
   datarehilo[idx3] = zvrehilo[tdx];
   datarelolo[idx3] = zvrelolo[tdx];
   dataimhihi[idx3] = zvimhihi[tdx]; 
   dataimlohi[idx3] = zvimlohi[tdx]; 
   dataimhilo[idx3] = zvimhilo[tdx];
   dataimlolo[idx3] = zvimlolo[tdx];
   __syncthreads();
}
*/

__global__ void cmplx4_padded_convjobs
 ( double *datarehihi, double *datarelohi,
   double *datarehilo, double *datarelolo,
   double *dataimhihi, double *dataimlohi,
   double *dataimhilo, double *dataimlolo,
   int *in1idx, int *in2idx, int *outidx, int dim )
{
   const int bdx = blockIdx.x;           // index to the convolution job
   const int tdx = threadIdx.x;          // index to the output of the job
   const int idx1 = in1idx[bdx] + tdx;
   const int idx2 = in2idx[bdx] + tdx;
   const int idx3 = outidx[bdx] + tdx;

   __shared__ double xvrehihi[qd_shmemsize];
   __shared__ double xvrelohi[qd_shmemsize];
   __shared__ double xvrehilo[qd_shmemsize];
   __shared__ double xvrelolo[qd_shmemsize];
   __shared__ double xvimhihi[qd_shmemsize];
   __shared__ double xvimlohi[qd_shmemsize];
   __shared__ double xvimhilo[qd_shmemsize];
   __shared__ double xvimlolo[qd_shmemsize];
   __shared__ double yvrehihi[2*qd_shmemsize];
   __shared__ double yvrelohi[2*qd_shmemsize];
   __shared__ double yvrehilo[2*qd_shmemsize];
   __shared__ double yvrelolo[2*qd_shmemsize];
   __shared__ double yvimhihi[2*qd_shmemsize];
   __shared__ double yvimlohi[2*qd_shmemsize];
   __shared__ double yvimhilo[2*qd_shmemsize];
   __shared__ double yvimlolo[2*qd_shmemsize];
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
   int ydx = dim + tdx;

   xvrehihi[tdx] = datarehihi[idx1];  // loading first input
   xvrelohi[tdx] = datarelohi[idx1]; 
   xvrehilo[tdx] = datarehilo[idx1]; 
   xvrelolo[tdx] = datarelolo[idx1]; 
   xvimhihi[tdx] = dataimhihi[idx1];
   xvimlohi[tdx] = dataimlohi[idx1];
   xvimhilo[tdx] = dataimhilo[idx1]; 
   xvimlolo[tdx] = dataimlolo[idx1]; 
   yvrehihi[tdx] = 0.0;           // padded with zeros
   yvrelohi[tdx] = 0.0;
   yvrehilo[tdx] = 0.0;
   yvrelolo[tdx] = 0.0;
   yvimhihi[tdx] = 0.0;
   yvimlohi[tdx] = 0.0;
   yvimhilo[tdx] = 0.0;
   yvimlolo[tdx] = 0.0;
   yvrehihi[ydx] = datarehihi[idx2];  // loading second input
   yvrelohi[ydx] = datarelohi[idx2];
   yvrehilo[ydx] = datarehilo[idx2];
   yvrelolo[ydx] = datarelolo[idx2];
   yvimhihi[ydx] = dataimhihi[idx2];
   yvimlohi[ydx] = dataimlohi[idx2];
   yvimhilo[ydx] = dataimhilo[idx2];
   yvimlolo[ydx] = dataimlolo[idx2];

   __syncthreads();

   // z[tdx] = x[0]*y[tdx]
   xrhihi = xvrehihi[0]; xrlohi = xvrelohi[0];
   xrhilo = xvrehilo[0]; xrlolo = xvrelolo[0];
   xihihi = xvimhihi[0]; xilohi = xvimlohi[0];
   xihilo = xvimhilo[0]; xilolo = xvimlolo[0];
   yrhihi = yvrehihi[ydx]; yrlohi = yvrelohi[ydx];
   yrhilo = yvrehilo[ydx]; yrlolo = yvrelolo[ydx];
   yihihi = yvimhihi[ydx]; yilohi = yvimlohi[ydx];
   yihilo = yvimhilo[ydx]; yilolo = yvimlolo[ydx];

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

   zvrehihi[tdx] = zrhihi; zvrelohi[tdx] = zrlohi;
   zvrehilo[tdx] = zrhilo; zvrelolo[tdx] = zrlolo;
   zvimhihi[tdx] = zihihi; zvimlohi[tdx] = zilohi;
   zvimhilo[tdx] = zihilo; zvimlolo[tdx] = zilolo;

   __syncthreads();

   for(int i=1; i<dim; i++) // z[tdx] = z[tdx] + x[i]*y[tdx-i]
   {
      ydx = dim + tdx - i;
      xrhihi = xvrehihi[i]; xrlohi = xvrelohi[i];
      xrhilo = xvrehilo[i]; xrlolo = xvrelolo[i];
      xihihi = xvimhihi[i]; xilohi = xvimlohi[i];
      xihilo = xvimhilo[i]; xilolo = xvimlolo[i];
      yrhihi = yvrehihi[ydx]; yrlohi = yvrelohi[ydx];
      yrhilo = yvrehilo[ydx]; yrlolo = yvrelolo[ydx];
      yihihi = yvimhihi[ydx]; yilohi = yvimlohi[ydx];
      yihilo = yvimhilo[ydx]; yilolo = yvimlolo[ydx];

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
      qdg_inc(&zvrehihi[tdx],&zvrelohi[tdx],&zvrehilo[tdx],&zvrelolo[tdx],
              zrhihi,zrlohi,zrhilo,zrlolo);           // zvre[k] += zr;
      __syncthreads();
      qdg_inc(&zvimhihi[tdx],&zvimlohi[tdx],&zvimhilo[tdx],&zvimlolo[tdx],
              zihihi,zilohi,zihilo,zilolo);           // zvim[k] += zi;
      __syncthreads();
   }
   __syncthreads();

   datarehihi[idx3] = zvrehihi[tdx]; // storing the output
   datarelohi[idx3] = zvrelohi[tdx];
   datarehilo[idx3] = zvrehilo[tdx];
   datarelolo[idx3] = zvrelolo[tdx];
   dataimhihi[idx3] = zvimhihi[tdx]; 
   dataimlohi[idx3] = zvimlohi[tdx]; 
   dataimhilo[idx3] = zvimhilo[tdx];
   dataimlolo[idx3] = zvimlolo[tdx];
   __syncthreads();
}

__global__ void cmplx4vectorized_flipsigns
 ( double *datarihihi, double *datarilohi,
   double *datarihilo, double *datarilolo, int *flpidx, int dim )
{
   const int bdx = blockIdx.x;    // index to the series to flip
   const int tdx = threadIdx.x;
   const int idx = flpidx[bdx] + tdx; // which number to flip

   double x; // register to load data from global memory

   x = datarihihi[idx]; x = -x; datarihihi[idx] = x;
   x = datarilohi[idx]; x = -x; datarilohi[idx] = x;
   x = datarihilo[idx]; x = -x; datarihilo[idx] = x;
   x = datarilolo[idx]; x = -x; datarilolo[idx] = x;
}

void GPU_cmplx4vectorized_flipsigns
 ( int deg, int nbrflips, int *flipidx,
   double *datarihihi, double *datarilohi,
   double *datarihilo, double *datarilolo,
   double *elapsedms, int vrblvl )
{
   const int deg1 = deg+1;

   int *flipidx_d;                   // flip indices on device
   const size_t szjobidx = nbrflips*sizeof(int);
   cudaMalloc((void**)&flipidx_d,szjobidx);
   cudaMemcpy(flipidx_d,flipidx,szjobidx,cudaMemcpyHostToDevice);

   cudaEvent_t start,stop;           // to measure time spent by kernels
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;

   if(vrblvl > 1)
   {
      cout << "flip indices :";
      for(int i=0; i<nbrflips; i++) cout << " " << flipidx[i];
      cout << endl;
      cout << "launching " << nbrflips << " flip signs blocks of "
                           << deg1 << " threads ..." << endl;
   }
   cudaEventRecord(start);
   cmplx4vectorized_flipsigns<<<nbrflips,deg1>>>
      (datarihihi,datarilohi,datarihilo,datarilolo,flipidx_d,deg1);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);

   *elapsedms = (double) milliseconds;

   if(vrblvl > 1)
   {
       cout << fixed << setprecision(2);
       cout << "Time spent by flip sign kernels : ";
       cout << *elapsedms << " milliseconds." << endl;
       cout << scientific << setprecision(16);
   }
   cudaFree(flipidx_d);
}

__global__ void dbl4_update_addjobs
 ( double *datahihi, double *datalohi, double *datahilo, double *datalolo,
   int *in1idx, int *in2idx, int *outidx, int dim )
{
   const int bdx = blockIdx.x;           // index to the convolution job
   const int tdx = threadIdx.x;          // index to the output of the job
   const int idx1 = in1idx[bdx] + tdx;
   const int idx2 = in2idx[bdx] + tdx;
   const int idx3 = outidx[bdx] + tdx;

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

   xvhihi[tdx] = datahihi[idx1];  // loading first input
   xvlohi[tdx] = datalohi[idx1];
   xvhilo[tdx] = datahilo[idx1];
   xvlolo[tdx] = datalolo[idx1];
   yvhihi[tdx] = datahihi[idx2];  // loading second input
   yvlohi[tdx] = datalohi[idx2];
   yvhilo[tdx] = datahilo[idx2];
   yvlolo[tdx] = datalolo[idx2];

   // zv[tdx] = xv[tdx] + yv[tdx];

   __syncthreads();

   qdg_add( xvhihi[tdx], xvlohi[tdx], xvhilo[tdx], xvlolo[tdx],
            yvhihi[tdx], yvlohi[tdx], yvhilo[tdx], yvlolo[tdx],
           &zvhihi[tdx],&zvlohi[tdx],&zvhilo[tdx],&zvlolo[tdx]);

   __syncthreads();

   datahihi[idx3] = zvhihi[tdx]; // storing the output
   datalohi[idx3] = zvlohi[tdx];
   datahilo[idx3] = zvhilo[tdx];
   datalolo[idx3] = zvlolo[tdx];
}

__global__ void cmplx4_update_addjobs
 ( double *datarehihi, double *datarelohi,
   double *datarehilo, double *datarelolo,
   double *dataimhihi, double *dataimlohi,
   double *dataimhilo, double *dataimlolo,
   int *in1idx, int *in2idx, int *outidx, int dim )
{
   const int bdx = blockIdx.x;           // index to the addition job
   const int tdx = threadIdx.x;          // index to the output of the job
   const int idx1 = in1idx[bdx] + tdx;
   const int idx2 = in2idx[bdx] + tdx;
   const int idx3 = outidx[bdx] + tdx;

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

   xvrehihi[tdx] = datarehihi[idx1];  // loading first input
   xvrelohi[tdx] = datarelohi[idx1];
   xvrehilo[tdx] = datarehilo[idx1];
   xvrelolo[tdx] = datarelolo[idx1];
   xvimhihi[tdx] = dataimhihi[idx1];
   xvimlohi[tdx] = dataimlohi[idx1];
   xvimhilo[tdx] = dataimhilo[idx1];
   xvimlolo[tdx] = dataimlolo[idx1];
   yvrehihi[tdx] = datarehihi[idx2];  // loading second input
   yvrelohi[tdx] = datarelohi[idx2];
   yvrehilo[tdx] = datarehilo[idx2];
   yvrelolo[tdx] = datarelolo[idx2];
   yvimhihi[tdx] = dataimhihi[idx2];
   yvimlohi[tdx] = dataimlohi[idx2];
   yvimhilo[tdx] = dataimhilo[idx2];
   yvimlolo[tdx] = dataimlolo[idx2];

   // zv[tdx] = xv[tdx] + yv[tdx];

   qdg_add(xvrehihi[tdx],xvrelohi[tdx],xvrehilo[tdx],xvrelolo[tdx],
           yvrehihi[tdx],yvrelohi[tdx],yvrehilo[tdx],yvrelolo[tdx],
           &zvrehihi[tdx],&zvrelohi[tdx],&zvrehilo[tdx],&zvrelolo[tdx]);
   __syncthreads();

   qdg_add(xvimhihi[tdx],xvimlohi[tdx],xvimhilo[tdx],xvimlolo[tdx],
           yvimhihi[tdx],yvimlohi[tdx],yvimhilo[tdx],yvimlolo[tdx],
           &zvimhihi[tdx],&zvimlohi[tdx],&zvimhilo[tdx],&zvimlolo[tdx]);
   __syncthreads();

   datarehihi[idx3] = zvrehihi[tdx]; // storing the output
   datarelohi[idx3] = zvrelohi[tdx];
   datarehilo[idx3] = zvrehilo[tdx];
   datarelolo[idx3] = zvrelolo[tdx];
   dataimhihi[idx3] = zvimhihi[tdx];
   dataimlohi[idx3] = zvimlohi[tdx];
   dataimhilo[idx3] = zvimhilo[tdx];
   dataimlolo[idx3] = zvimlolo[tdx];
}

void dbl_convoluted_data4_to_output
 ( double *datahihi, double *datalohi, double *datahilo, double *datalolo,
   double **outputhihi, double **outputlohi,
   double **outputhilo, double **outputlolo,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart, int vrblvl )
{
   const int deg1 = deg+1;
   int ix0,ix1,ix2;

   for(int i=0; i<=deg; i++) // output[dim][i] = data[i];
   {
      outputhihi[dim][i] = datahihi[i];
      outputlohi[dim][i] = datalohi[i];
      outputhilo[dim][i] = datahilo[i];
      outputlolo[dim][i] = datalolo[i];
   }
   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++) // output[i][j] = 0.0;
      {
         outputhihi[i][j] = 0.0;
         outputlohi[i][j] = 0.0;
         outputhilo[i][j] = 0.0;
         outputlolo[i][j] = 0.0;
      }

   for(int k=0; k<nbr; k++)
   {
      ix1 = fstart[k] + (nvr[k]-1)*deg1;
      
      if(vrblvl > 1)
         cout << "monomial " << k << " update starts at " << ix1 << endl;

      for(int i=0; i<=deg; i++) // output[dim][i] += data[ix1++];
         qdf_inc(&outputhihi[dim][i],&outputlohi[dim][i],
                 &outputhilo[dim][i],&outputlolo[dim][i],
                    datahihi[ix1],      datalohi[ix1],
                    datahilo[ix1],      datalolo[ix1++]);
     
      ix0 = idx[k][0];
      if(nvr[k] == 1)
      {
         ix1 = (1 + k)*deg1;
            
         for(int i=0; i<=deg; i++) // output[ix0][i] += data[ix1++];
            qdf_inc(&outputhihi[ix0][i],&outputlohi[ix0][i],
                    &outputhilo[ix0][i],&outputlolo[ix0][i],
                       datahihi[ix1],      datalohi[ix1],
                       datahilo[ix1],      datalolo[ix1++]);
      }
      else
      {                               // update first and last derivative
         ix2 = nvr[k]-3;
         if(ix2 < 0) ix2 = 0;
         ix1 = bstart[k] + ix2*deg1;

         for(int i=0; i<=deg; i++) // output[ix0][i] += data[ix1++];
            qdf_inc(&outputhihi[ix0][i],&outputlohi[ix0][i],
                    &outputhilo[ix0][i],&outputlolo[ix0][i],
                       datahihi[ix1],      datalohi[ix1],
                       datahilo[ix1],      datalolo[ix1++]);

         ix2 = nvr[k]-2;
         ix1 = fstart[k] + ix2*deg1;
         ix0 = idx[k][ix2+1];

         for(int i=0; i<=deg; i++) // output[ix0][i] += data[ix1++];
            qdf_inc(&outputhihi[ix0][i],&outputlohi[ix0][i],
                    &outputhilo[ix0][i],&outputlolo[ix0][i],
                       datahihi[ix1],      datalohi[ix1],
                       datahilo[ix1],      datalolo[ix1++]);
 
         if(nvr[k] > 2)                   // update all other derivatives
         {
            for(int j=1; j<nvr[k]-1; j++)
            {
               ix0 = idx[k][j];            // j-th variable in monomial k
               ix1 = cstart[k] + (j-1)*deg1;

               if(vrblvl > 1)
                  cout << "monomial " << k << " derivative " << ix0
                       << " update starts at " << ix1 << endl;

               for(int i=0; i<=deg; i++) // output[ix0][i] += data[ix1++];
                  qdf_inc(&outputhihi[ix0][i],&outputlohi[ix0][i],
                          &outputhilo[ix0][i],&outputlolo[ix0][i],
                             datahihi[ix1],      datalohi[ix1],
                             datahilo[ix1],      datalolo[ix1++]);
            }
         }
      }
   }
}

void cmplx_convoluted_data4_to_output
 ( double *datarehihi, double *datarelohi,
   double *datarehilo, double *datarelolo,
   double *dataimhihi, double *dataimlohi,
   double *dataimhilo, double *dataimlolo,
   double **outputrehihi, double **outputrelohi,
   double **outputrehilo, double **outputrelolo,
   double **outputimhihi, double **outputimlohi,
   double **outputimhilo, double **outputimlolo,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart, int vrblvl )
{
   const int deg1 = deg+1;
   int ix0,ix1,ix2;

   for(int i=0; i<=deg; i++) // output[dim][i] = data[i];
   {
      outputrehihi[dim][i] = datarehihi[i];
      outputrelohi[dim][i] = datarelohi[i];
      outputrehilo[dim][i] = datarehilo[i];
      outputrelolo[dim][i] = datarelolo[i];
      outputimhihi[dim][i] = dataimhihi[i];
      outputimlohi[dim][i] = dataimlohi[i];
      outputimhilo[dim][i] = dataimhilo[i];
      outputimlolo[dim][i] = dataimlolo[i];
   }
   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++) // output[i][j] = 0.0;
      {
         outputrehihi[i][j] = 0.0;
         outputrelohi[i][j] = 0.0;
         outputrehilo[i][j] = 0.0;
         outputrelolo[i][j] = 0.0;
         outputimhihi[i][j] = 0.0;
         outputimlohi[i][j] = 0.0;
         outputimhilo[i][j] = 0.0;
         outputimlolo[i][j] = 0.0;
      }

   for(int k=0; k<nbr; k++)
   {
      ix1 = fstart[k] + (nvr[k]-1)*deg1;
      
      if(vrblvl > 1)
         cout << "monomial " << k << " update starts at " << ix1 << endl;

      for(int i=0; i<=deg; i++) // output[dim][i] += data[ix1++];
      {
         qdf_inc(&outputrehihi[dim][i],&outputrelohi[dim][i],
                 &outputrehilo[dim][i],&outputrelolo[dim][i],
                 datarehihi[ix1],datarelohi[ix1],
                 // datarehilo[ix1],datarelolo[ix1++]);
                 datarehilo[ix1],datarelolo[ix1]);
         qdf_inc(&outputimhihi[dim][i],&outputimlohi[dim][i],
                 &outputimhilo[dim][i],&outputimlolo[dim][i],
                 dataimhihi[ix1],dataimlohi[ix1],
                 dataimhilo[ix1],dataimlolo[ix1++]);
      }
      ix0 = idx[k][0];
      if(nvr[k] == 1)
      {
         ix1 = (1 + k)*deg1;
            
         for(int i=0; i<=deg; i++) // output[ix0][i] += data[ix1++];
         {
            qdf_inc(&outputrehihi[ix0][i],&outputrelohi[ix0][i],
                    &outputrehilo[ix0][i],&outputrelolo[ix0][i],
                    datarehihi[ix1],datarelohi[ix1],
                    // datarehilo[ix1],datarelolo[ix1++]);
                    datarehilo[ix1],datarelolo[ix1]);
            qdf_inc(&outputimhihi[ix0][i],&outputimlohi[ix0][i],
                    &outputimhilo[ix0][i],&outputimlolo[ix0][i],
                    dataimhihi[ix1],dataimlohi[ix1],
                    dataimhilo[ix1],dataimlolo[ix1++]);
         }
      }
      else
      {                               // update first and last derivative
         ix2 = nvr[k]-3;
         if(ix2 < 0) ix2 = 0;
         ix1 = bstart[k] + ix2*deg1;

         for(int i=0; i<=deg; i++) // output[ix0][i] += data[ix1++];
         {
            qdf_inc(&outputrehihi[ix0][i],&outputrelohi[ix0][i],
                    &outputrehilo[ix0][i],&outputrelolo[ix0][i],
                    datarehihi[ix1],datarelohi[ix1],
                    // datarehilo[ix1],datarelolo[ix1++]);
                    datarehilo[ix1],datarelolo[ix1]);
            qdf_inc(&outputimhihi[ix0][i],&outputimlohi[ix0][i],
                    &outputimhilo[ix0][i],&outputimlolo[ix0][i],
                    dataimhihi[ix1],dataimlohi[ix1],
                    dataimhilo[ix1],dataimlolo[ix1++]);
         }
         ix2 = nvr[k]-2;
         ix1 = fstart[k] + ix2*deg1;
         ix0 = idx[k][ix2+1];

         for(int i=0; i<=deg; i++) // output[ix0][i] += data[ix1++];
         {
            qdf_inc(&outputrehihi[ix0][i],&outputrelohi[ix0][i],
                    &outputrehilo[ix0][i],&outputrelolo[ix0][i],
                    datarehihi[ix1],datarelohi[ix1],
                    // datarehilo[ix1],datarelolo[ix1++]);
                    datarehilo[ix1],datarelolo[ix1]);
            qdf_inc(&outputimhihi[ix0][i],&outputimlohi[ix0][i],
                    &outputimhilo[ix0][i],&outputimlolo[ix0][i],
                    dataimhihi[ix1],dataimlohi[ix1],
                    dataimhilo[ix1],dataimlolo[ix1++]);
         }
         if(nvr[k] > 2)                   // update all other derivatives
         {
            for(int j=1; j<nvr[k]-1; j++)
            {
               ix0 = idx[k][j];            // j-th variable in monomial k
               ix1 = cstart[k] + (j-1)*deg1;

               if(vrblvl > 1)
                  cout << "monomial " << k << " derivative " << ix0
                       << " update starts at " << ix1 << endl;

               for(int i=0; i<=deg; i++) // output[ix0][i] += data[ix1++];
               {
                  qdf_inc(&outputrehihi[ix0][i],&outputrelohi[ix0][i],
                          &outputrehilo[ix0][i],&outputrelolo[ix0][i],
                          datarehihi[ix1],datarelohi[ix1],
                          // datarehilo[ix1],datarelolo[ix1++]);
                          datarehilo[ix1],datarelolo[ix1]);
                  qdf_inc(&outputimhihi[ix0][i],&outputimlohi[ix0][i],
                          &outputimhilo[ix0][i],&outputimlolo[ix0][i],
                          dataimhihi[ix1],dataimlohi[ix1],
                          dataimhilo[ix1],dataimlolo[ix1++]);
               }
            }
         }
      }
   }
}

void dbl_added_data4_to_output
 ( double *datahihi, double *datalohi, double *datahilo, double *datalolo,
   double **outputhihi, double **outputlohi,
   double **outputhilo, double **outputlolo,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart,
   AdditionJobs jobs, int vrblvl )
{
   const int deg1 = deg + 1;
   const int lastmon = nbr-1;
   const int lastidx = nvr[lastmon]-1;
   int ix;

   ix = fstart[lastmon] + lastidx*deg1;

   if(vrblvl > 1)
      cout << "Updating value starting at " << ix << " in data." << endl;

   for(int i=0; i<=deg; i++) // output[dim][i] = data[ix++];
   {
      outputhihi[dim][i] = datahihi[ix];
      outputlohi[dim][i] = datalohi[ix];
      outputhilo[dim][i] = datahilo[ix];
      outputlolo[dim][i] = datalolo[ix++];
   }
   int cnt = jobs.get_differential_count(0);

   if(vrblvl > 1)
      cout << "Differential count for variable 0 : " << cnt << endl;

   if(cnt == 0) // it could be there is no first variable anywhere ...
   {
      const int difidx = jobs.get_differential_index(0,0);

      if(vrblvl > 1)
         cout << "Differential index for variable 0 : " << difidx << endl;

      if(difidx < 0)
      {
         for(int i=0; i<=deg; i++) // output[0][i] = 0.0;
         {
            outputhihi[0][i] = 0.0; outputlohi[0][i] = 0.0;
            outputhilo[0][i] = 0.0; outputlolo[0][i] = 0.0;
         }
      }
      else
      {
         int ix1re = (1 + difidx)*deg1;

         if(vrblvl > 1)
            cout << "updating derivative with coefficient ..." << endl;

         for(int i=0; i<=deg; i++)
         {
            outputhihi[0][i] = datahihi[ix1re];
            outputlohi[0][i] = datalohi[ix1re];
            outputhilo[0][i] = datahilo[ix1re];
            outputlolo[0][i] = datalolo[ix1re++];
         }
      }
   }
   else
   {
      int ix0 = jobs.get_differential_index(0,cnt);
      int ix2 = nvr[ix0]-3;
      if(ix2 < 0) ix2 = 0; // on GPU, one backward item less

      ix = bstart[ix0] + ix2*deg1;
      
      if(vrblvl > 1)
         cout << "Updating derivative 0 at " << ix << " in data." << endl;

      for(int i=0; i<=deg; i++) // output[0][i] = data[ix++];
      {
         outputhihi[0][i] = datahihi[ix];
         outputlohi[0][i] = datalohi[ix];
         outputhilo[0][i] = datahilo[ix];
         outputlolo[0][i] = datalolo[ix++];
      }
   }
   for(int k=1; k<dim; k++) // updating all other derivatives
   {
      int cnt = jobs.get_differential_count(k);

      if(vrblvl > 1)
         cout << "Differential count for variable " << k
              << " : " << cnt << endl;

      if(cnt == 0) // it could be there is no variable k anywhere ...
      {
         const int difidx = jobs.get_differential_index(k,0);

         if(vrblvl > 1)
            cout << "Differential index for variable " << k 
                 << " : " << difidx << endl;

         if(difidx < 0)
         {
            for(int i=0; i<=deg; i++) // output[k][i] = 0.0;
            {
               outputhihi[k][i] = 0.0; outputlohi[k][i] = 0.0;
               outputhilo[k][i] = 0.0; outputlolo[k][i] = 0.0;
            }
         }
         else
         {
            int ix1re = (1 + difidx)*deg1;

            if(vrblvl > 1)
               cout << "updating derivative with coefficient ..." << endl;

            for(int i=0; i<=deg; i++)
            {
               outputhihi[k][i] = datahihi[ix1re];
               outputlohi[k][i] = datalohi[ix1re];
               outputhilo[k][i] = datahilo[ix1re];
               outputlolo[k][i] = datalolo[ix1re++];
            }
         }
      }
      else
      {
         int ix0 = jobs.get_differential_index(k,cnt);

         if(idx[ix0][0] == k) // k is first variable of monomial
         {
            int ix2 = nvr[ix0]-3;
            if(ix2 < 0) ix2 = 0;

            if(vrblvl > 1)
               cout << "Updating derivative " << k 
                    << " at " << ix << " in data." << endl;

            ix = bstart[ix0] + ix2*deg1;

            for(int i=0; i<=deg; i++) // output[k][i] = data[ix++];
            {
               outputhihi[k][i] = datahihi[ix];
               outputlohi[k][i] = datalohi[ix];
               outputhilo[k][i] = datahilo[ix];
               outputlolo[k][i] = datalolo[ix++];
            }
         }
         else if(idx[ix0][nvr[ix0]-1] == k) // k is last variable
         {
            int ix2 = nvr[ix0]-2;
 
            if(vrblvl > 1)
               cout << "Updating derivative " << k 
                    << " at " << ix << " in data." << endl;

            ix = fstart[ix0] + ix2*deg1;

            for(int i=0; i<=deg; i++) // output[k][i] = data[ix++];
            {
               outputhihi[k][i] = datahihi[ix];
               outputlohi[k][i] = datalohi[ix];
               outputhilo[k][i] = datahilo[ix];
               outputlolo[k][i] = datalolo[ix++];
            }
         }
         else // derivative is in some cross product
         {
            int ix2 = jobs.position(nvr[ix0],idx[ix0],k) - 1;
   
            if(vrblvl > 1)
               cout << "Updating derivative " << k 
                    << " at " << ix << " in data." << endl;

            ix = cstart[ix0] + ix2*deg1;

            for(int i=0; i<=deg; i++) // output[k][i] = data[ix++];
            {
               outputhihi[k][i] = datahihi[ix];
               outputlohi[k][i] = datalohi[ix];
               outputhilo[k][i] = datahilo[ix];
               outputlolo[k][i] = datalolo[ix++];
            }
         }
      }
   }
}

void cmplx_added_data4_to_output
 ( double *datarehihi, double *datarelohi,
   double *datarehilo, double *datarelolo,
   double *dataimhihi, double *dataimlohi,
   double *dataimhilo, double *dataimlolo,
   double **outputrehihi, double **outputrelohi,
   double **outputrehilo, double **outputrelolo,
   double **outputimhihi, double **outputimlohi,
   double **outputimhilo, double **outputimlolo,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart, AdditionJobs jobs,
   int vrblvl )
{
   const int deg1 = deg + 1;
   const int lastmon = nbr-1;
   const int lastidx = nvr[lastmon]-1;
   int ix;

   ix = fstart[lastmon] + lastidx*deg1;

   if(vrblvl > 1)
      cout << "Updating value starting at " << ix << " in data." << endl;

   for(int i=0; i<=deg; i++) // output[dim][i] = data[ix++];
   {
      outputrehihi[dim][i] = datarehihi[ix];
      outputrelohi[dim][i] = datarelohi[ix];
      outputrehilo[dim][i] = datarehilo[ix];
      outputrelolo[dim][i] = datarelolo[ix];
      outputimhihi[dim][i] = dataimhihi[ix];
      outputimlohi[dim][i] = dataimlohi[ix];
      outputimhilo[dim][i] = dataimhilo[ix];
      outputimlolo[dim][i] = dataimlolo[ix++];
   }
   int cnt = jobs.get_differential_count(0);

   if(vrblvl > 1)
      cout << "Differential count for variable 0 : " << cnt << endl;

   if(cnt == 0) // it could be there is no first variable anywhere ...
   {
      const int difidx = jobs.get_differential_index(0,0);

      if(vrblvl > 1)
         cout << "Differential index for variable 0 : " << difidx << endl;

      if(difidx < 0)
      {
         for(int i=0; i<=deg; i++) // output[0][i] = 0.0;
         {
            outputrehihi[0][i] = 0.0; outputrelohi[0][i] = 0.0;
            outputrehilo[0][i] = 0.0; outputrelolo[0][i] = 0.0;
            outputimhihi[0][i] = 0.0; outputimlohi[0][i] = 0.0; 
            outputimhilo[0][i] = 0.0; outputimlolo[0][i] = 0.0;
         }
      }
      else
      {
         int cffidx = (1 + difidx)*deg1;

         if(vrblvl > 1)
            cout << "updating derivative with coefficient ..." << endl;

         for(int i=0; i<=deg; i++)
         {
            outputrehihi[0][i] = datarehihi[cffidx];
            outputrelohi[0][i] = datarelohi[cffidx];
            outputrehilo[0][i] = datarehilo[cffidx];
            outputrelolo[0][i] = datarelolo[cffidx];
            outputimhihi[0][i] = dataimhihi[cffidx];
            outputimlohi[0][i] = dataimlohi[cffidx];
            outputimhilo[0][i] = dataimhilo[cffidx];
            outputimlolo[0][i] = dataimlolo[cffidx++];
         }
      }
   }
   else
   {
      int ix0 = jobs.get_differential_index(0,cnt);
      int ix2 = nvr[ix0]-3;
      if(ix2 < 0) ix2 = 0; // on GPU, one backward item less

      ix = bstart[ix0] + ix2*deg1;
      
      if(vrblvl > 1)
         cout << "Updating derivative 0 at " << ix << " in data." << endl;

      for(int i=0; i<=deg; i++) // output[0][i] = data[ix++];
      {
         outputrehihi[0][i] = datarehihi[ix];
         outputrelohi[0][i] = datarelohi[ix];
         outputrehilo[0][i] = datarehilo[ix];
         outputrelolo[0][i] = datarelolo[ix];
         outputimhihi[0][i] = dataimhihi[ix];
         outputimlohi[0][i] = dataimlohi[ix];
         outputimhilo[0][i] = dataimhilo[ix];
         outputimlolo[0][i] = dataimlolo[ix++];
      }
   }
   for(int k=1; k<dim; k++) // updating all other derivatives
   {
      int cnt = jobs.get_differential_count(k);

      if(vrblvl > 1)
         cout << "Differential count for variable " << k
              << " : " << cnt << endl;

      if(cnt == 0) // it could be there is no variable k anywhere ...
      {
         const int difidx = jobs.get_differential_index(k,0);

         if(vrblvl > 1)
            cout << "Differential index for variable " << k 
                 << " : " << difidx << endl;

         if(difidx < 0)
         {
            for(int i=0; i<=deg; i++) // output[k][i] = 0.0;
            {
               outputrehihi[k][i] = 0.0; outputrelohi[k][i] = 0.0;
               outputrehilo[k][i] = 0.0; outputrelolo[k][i] = 0.0;
               outputimhihi[k][i] = 0.0; outputimlohi[k][i] = 0.0;
               outputimhilo[k][i] = 0.0; outputimlolo[k][i] = 0.0;
            }
         }
         else
         {
            int cffidx = (1 + difidx)*deg1;

            if(vrblvl > 1)
               cout << "updating derivative with coefficient ..." << endl;

            for(int i=0; i<=deg; i++)
            {
               outputrehihi[k][i] = datarehihi[cffidx];
               outputrelohi[k][i] = datarelohi[cffidx];
               outputrehilo[k][i] = datarehilo[cffidx];
               outputrelolo[k][i] = datarelolo[cffidx];
               outputimhihi[k][i] = dataimhihi[cffidx];
               outputimlohi[k][i] = dataimlohi[cffidx];
               outputimhilo[k][i] = dataimhilo[cffidx];
               outputimlolo[k][i] = dataimlolo[cffidx++];
            }
         }
      }
      else
      {
         int ix0 = jobs.get_differential_index(k,cnt);
 
         if(idx[ix0][0] == k) // k is first variable of monomial
         {
            int ix2 = nvr[ix0]-3;
            if(ix2 < 0) ix2 = 0;

            if(vrblvl > 1)
               cout << "Updating derivative " << k 
                    << " at " << ix << " in data." << endl;

            ix = bstart[ix0] + ix2*deg1;

            for(int i=0; i<=deg; i++) // output[k][i] = data[ix++];
            {
               outputrehihi[k][i] = datarehihi[ix];
               outputrelohi[k][i] = datarelohi[ix];
               outputrehilo[k][i] = datarehilo[ix];
               outputrelolo[k][i] = datarelolo[ix];
               outputimhihi[k][i] = dataimhihi[ix];
               outputimlohi[k][i] = dataimlohi[ix];
               outputimhilo[k][i] = dataimhilo[ix];
               outputimlolo[k][i] = dataimlolo[ix++];
            }
         }
         else if(idx[ix0][nvr[ix0]-1] == k) // k is last variable
         {
            int ix2 = nvr[ix0]-2;

            if(vrblvl > 1)
               cout << "Updating derivative " << k 
                    << " at " << ix << " in data." << endl;

            ix = fstart[ix0] + ix2*deg1;

            for(int i=0; i<=deg; i++) // output[k][i] = data[ix++];
            {
               outputrehihi[k][i] = datarehihi[ix];
               outputrelohi[k][i] = datarelohi[ix];
               outputrehilo[k][i] = datarehilo[ix];
               outputrelolo[k][i] = datarelolo[ix];
               outputimhihi[k][i] = dataimhihi[ix];
               outputimlohi[k][i] = dataimlohi[ix];
               outputimhilo[k][i] = dataimhilo[ix];
               outputimlolo[k][i] = dataimlolo[ix++];
            }
         }
         else // derivative is in some cross product
         {
            int ix2 = jobs.position(nvr[ix0],idx[ix0],k) - 1;

            if(vrblvl > 1)
               cout << "Updating derivative " << k 
                    << " at " << ix << " in data." << endl;

            ix = cstart[ix0] + ix2*deg1;

            for(int i=0; i<=deg; i++) // output[k][i] = data[ix++];
            {
               outputrehihi[k][i] = datarehihi[ix];
               outputrelohi[k][i] = datarelohi[ix];
               outputrehilo[k][i] = datarehilo[ix];
               outputrelolo[k][i] = datarelolo[ix];
               outputimhihi[k][i] = dataimhihi[ix];
               outputimlohi[k][i] = dataimlohi[ix];
               outputimhilo[k][i] = dataimhilo[ix];
               outputimlolo[k][i] = dataimlolo[ix++];
            }
         }
      }
   }
}

void cmplx_added_data4vectorized_to_output
 ( double *datarihihi, double *datarilohi,
   double *datarihilo, double *datarilolo,
   double **outputrehihi, double **outputrelohi,
   double **outputrehilo, double **outputrelolo,
   double **outputimhihi, double **outputimlohi,
   double **outputimhilo, double **outputimlolo,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart,
   int totcff, int offsetri, ComplexAdditionJobs jobs,
   int vrblvl )
{
   const int deg1 = deg + 1;
   const int lastmon = nbr-1;
   const int lastidx = nvr[lastmon]-1;
   const int totcffoffset = totcff + offsetri;
   int ix1re,ix2im;

/*
   for(int i=0; i<=dim; i++)    // initialize the entire output
      for(int j=0; j<deg1; j++)
      {
         outputrehihi[i][j] = 0.0; outputrelohi[i][j] = 0.0;
         outputrehilo[i][j] = 0.0; outputrelolo[i][j] = 0.0;
         outputimhihi[i][j] = 0.0; outputimlohi[i][j] = 0.0; 
         outputimhilo[i][j] = 0.0; outputimlolo[i][j] = 0.0;
      }
 */
   ix1re = fstart[lastmon] + lastidx*deg1;
   ix2im = fstart[lastmon] + lastidx*deg1 + totcffoffset;

   if(vrblvl > 1)
      cout << "Updating value starting at " << ix1re << " in data." << endl;

   for(int i=0; i<=deg; i++) // output[dim][i] = data[ix++];
   {
      outputrehihi[dim][i] = datarihihi[ix1re];
      outputrelohi[dim][i] = datarilohi[ix1re];
      outputrehilo[dim][i] = datarihilo[ix1re];
      outputrelolo[dim][i] = datarilolo[ix1re++];
      outputimhihi[dim][i] = datarihihi[ix2im];
      outputimlohi[dim][i] = datarilohi[ix2im];
      outputimhilo[dim][i] = datarihilo[ix2im];
      outputimlolo[dim][i] = datarilolo[ix2im++];
   }
   int cnt = jobs.get_differential_count(0);

   if(vrblvl > 1)
      cout << "Differential count for variable 0 : " << cnt << endl;

   if(cnt == 0) // it could be there is no first variable anywhere ...
   {
      const int difidx = jobs.get_differential_index(0,0);

      if(vrblvl > 1)
         cout << "Differential index for variable 0 : " << difidx << endl;

      if(difidx < 0)
      {
         for(int i=0; i<=deg; i++) // output[0][i] = 0.0;
         {
            outputrehihi[0][i] = 0.0; outputrelohi[0][i] = 0.0;
            outputrehilo[0][i] = 0.0; outputrelolo[0][i] = 0.0;
            outputimhihi[0][i] = 0.0; outputimlohi[0][i] = 0.0; 
            outputimhilo[0][i] = 0.0; outputimlolo[0][i] = 0.0;
         }
      }
      else
      {
         ix1re = (1 + difidx)*deg1;
         ix2im = (1 + difidx)*deg1 + totcffoffset;

         if(vrblvl > 1)
            cout << "updating derivative with coefficient ..." << endl;

         for(int i=0; i<=deg; i++)
         {
            outputrehihi[0][i] = datarihihi[ix1re];
            outputrelohi[0][i] = datarilohi[ix1re];
            outputrehilo[0][i] = datarihilo[ix1re];
            outputrelolo[0][i] = datarilolo[ix1re++];
            outputimhihi[0][i] = datarihihi[ix2im];
            outputimlohi[0][i] = datarilohi[ix2im];
            outputimhilo[0][i] = datarihilo[ix2im];
            outputimlolo[0][i] = datarilolo[ix2im++];
         }
      }
   }
   else
   {
      int ix0 = jobs.get_differential_index(0,cnt);
      // int ix2 = nvr[ix0]-3;
      // if(ix2 < 0) ix2 = 0; // on GPU, one backward item less
      int ix2 = nvr[ix0]-2; // no longer the case for vectorized

      ix1re = bstart[ix0] + ix2*deg1;
      ix2im = bstart[ix0] + ix2*deg1 + totcffoffset;
      
      if(vrblvl > 1)
         cout << "Updating derivative 0 at " << ix1re << " in data." << endl;

      for(int i=0; i<=deg; i++) // output[0][i] = data[ix++];
      {
         outputrehihi[0][i] = datarihihi[ix1re];
         outputrelohi[0][i] = datarilohi[ix1re];
         outputrehilo[0][i] = datarihilo[ix1re];
         outputrelolo[0][i] = datarilolo[ix1re++];
         outputimhihi[0][i] = datarihihi[ix2im];
         outputimlohi[0][i] = datarilohi[ix2im];
         outputimhilo[0][i] = datarihilo[ix2im];
         outputimlolo[0][i] = datarilolo[ix2im++];
      }
   }
   for(int k=1; k<dim; k++) // updating all other derivatives
   {
      int cnt = jobs.get_differential_count(k);

      if(vrblvl > 1)
         cout << "Differential count for variable " << k
              << " : " << cnt << endl;

      if(cnt == 0) // it could be there is no variable k anywhere ...
      {
         const int difidx = jobs.get_differential_index(k,0);

         if(vrblvl > 1)
            cout << "Differential index for variable " << k 
                 << " : " << difidx << endl;

         if(difidx < 0)
         {
            for(int i=0; i<=deg; i++) // output[k][i] = 0.0;
            {
               outputrehihi[k][i] = 0.0; outputrelohi[k][i] = 0.0;
               outputrehilo[k][i] = 0.0; outputrelolo[k][i] = 0.0;
               outputimhihi[k][i] = 0.0; outputimlohi[k][i] = 0.0;
               outputimhilo[k][i] = 0.0; outputimlolo[k][i] = 0.0;
            }
         }
         else
         {
            ix1re = (1 + difidx)*deg1;
            ix2im = (1 + difidx)*deg1 + totcffoffset;

            if(vrblvl > 1)
               cout << "updating derivative with coefficient ..." << endl;

            for(int i=0; i<=deg; i++)
            {
               outputrehihi[k][i] = datarihihi[ix1re];
               outputrelohi[k][i] = datarilohi[ix1re];
               outputrehilo[k][i] = datarihilo[ix1re];
               outputrelolo[k][i] = datarilolo[ix1re++];
               outputimhihi[k][i] = datarihihi[ix2im];
               outputimlohi[k][i] = datarilohi[ix2im];
               outputimhilo[k][i] = datarihilo[ix2im];
               outputimlolo[k][i] = datarilolo[ix2im++];
            }
         }
      }
      else
      {
         int ix0 = jobs.get_differential_index(k,cnt);

         if(idx[ix0][0] == k) // k is first variable of monomial
         {
            // int ix2 = nvr[ix0]-3;
            // if(ix2 < 0) ix2 = 0;
            int ix2 = nvr[ix0] - 2;

            ix1re = bstart[ix0] + ix2*deg1;
            ix2im = bstart[ix0] + ix2*deg1 + totcffoffset;

            if(vrblvl > 1)
               cout << "Updating derivative " << k 
                    << " at " << ix1re << " in data." << endl;

            for(int i=0; i<=deg; i++) // output[k][i] = data[ix++];
            {
               outputrehihi[k][i] = datarihihi[ix1re];
               outputrelohi[k][i] = datarilohi[ix1re];
               outputrehilo[k][i] = datarihilo[ix1re];
               outputrelolo[k][i] = datarilolo[ix1re++];
               outputimhihi[k][i] = datarihihi[ix2im];
               outputimlohi[k][i] = datarilohi[ix2im];
               outputimhilo[k][i] = datarihilo[ix2im];
               outputimlolo[k][i] = datarilolo[ix2im++];
            }
         }
         else if(idx[ix0][nvr[ix0]-1] == k) // k is last variable
         {
            int ix2 = nvr[ix0]-2;

            ix1re = fstart[ix0] + ix2*deg1;
            ix2im = fstart[ix0] + ix2*deg1 + totcffoffset;
 
            if(vrblvl > 1)
               cout << "Updating derivative " << k 
                    << " at " << ix1re << " in data." << endl;

            for(int i=0; i<=deg; i++) // output[k][i] = data[ix++];
            {
               outputrehihi[k][i] = datarihihi[ix1re];
               outputrelohi[k][i] = datarilohi[ix1re];
               outputrehilo[k][i] = datarihilo[ix1re];
               outputrelolo[k][i] = datarilolo[ix1re++];
               outputimhihi[k][i] = datarihihi[ix2im];
               outputimlohi[k][i] = datarilohi[ix2im];
               outputimhilo[k][i] = datarihilo[ix2im];
               outputimlolo[k][i] = datarilolo[ix2im++];
            }
         }
         else // derivative is in some cross product
         {
            int ix2 = jobs.position(nvr[ix0],idx[ix0],k) - 1;

            if(vrblvl > 1)
               cout << "Updating derivative " << k 
                    << " at " << ix1re << " in data." << endl;

            ix1re = cstart[ix0] + ix2*deg1;
            ix2im = cstart[ix0] + ix2*deg1 + totcffoffset;

            for(int i=0; i<=deg; i++) // output[k][i] = data[ix++];
            {
               outputrehihi[k][i] = datarihihi[ix1re];
               outputrelohi[k][i] = datarilohi[ix1re];
               outputrehilo[k][i] = datarihilo[ix1re];
               outputrelolo[k][i] = datarilolo[ix1re++];
               outputimhihi[k][i] = datarihihi[ix2im];
               outputimlohi[k][i] = datarilohi[ix2im];
               outputimhilo[k][i] = datarihilo[ix2im];
               outputimlolo[k][i] = datarilolo[ix2im++];
            }
         }
      }
   }
}

void dbl4_data_setup
 ( int dim, int nbr, int deg, int totcff,
   double *datahihi, double *datalohi, double *datahilo, double *datalolo,
   double *csthihi, double *cstlohi, double *csthilo, double *cstlolo,
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo,
   double **inputhihi, double **inputlohi,
   double **inputhilo, double **inputlolo )
{
   const int deg1 = deg+1;
   int ix = 0;

   for(int i=0; i<deg1; i++)
   {
      datahihi[ix]   = csthihi[i];
      datalohi[ix]   = cstlohi[i];
      datahilo[ix]   = csthilo[i];
      datalolo[ix++] = cstlolo[i];
   }
   for(int i=0; i<nbr; i++)
      for(int j=0; j<deg1; j++)
      {
         datahihi[ix]   = cffhihi[i][j];
         datalohi[ix]   = cfflohi[i][j];
         datahilo[ix]   = cffhilo[i][j];
         datalolo[ix++] = cfflolo[i][j];
      }
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++)
      {
         datahihi[ix]   = inputhihi[i][j];
         datalohi[ix]   = inputlohi[i][j];
         datahilo[ix]   = inputhilo[i][j];
         datalolo[ix++] = inputlolo[i][j];
      }

   for(int i=ix; i<totcff; i++)
   {
      datahihi[i] = 0.0; datalohi[i] = 0.0;
      datahilo[i] = 0.0; datalolo[i] = 0.0;
   }
}

void cmplx4_data_setup
 ( int dim, int nbr, int deg, int totcff,
   double *datarehihi, double *datarelohi,
   double *datarehilo, double *datarelolo,
   double *dataimhihi, double *dataimlohi,
   double *dataimhilo, double *dataimlolo,
   double *cstrehihi, double *cstrelohi,
   double *cstrehilo, double *cstrelolo,
   double *cstimhihi, double *cstimlohi,
   double *cstimhilo, double *cstimlolo,
   double **cffrehihi, double **cffrelohi,
   double **cffrehilo, double **cffrelolo,
   double **cffimhihi, double **cffimlohi,
   double **cffimhilo, double **cffimlolo,
   double **inputrehihi, double **inputrelohi,
   double **inputrehilo, double **inputrelolo,
   double **inputimhihi, double **inputimlohi,
   double **inputimhilo, double **inputimlolo )
{
   const int deg1 = deg+1;
   int ix = 0;

   for(int i=0; i<deg1; i++)
   {
      datarehihi[ix]   = cstrehihi[i];
      datarelohi[ix]   = cstrelohi[i];
      datarehilo[ix]   = cstrehilo[i];
      datarelolo[ix]   = cstrelolo[i];
      dataimhihi[ix]   = cstimhihi[i];
      dataimlohi[ix]   = cstimlohi[i];
      dataimhilo[ix]   = cstimhilo[i];
      dataimlolo[ix++] = cstimlolo[i];
   }
   for(int i=0; i<nbr; i++)
      for(int j=0; j<deg1; j++)
      {
         datarehihi[ix]   = cffrehihi[i][j];
         datarelohi[ix]   = cffrelohi[i][j];
         datarehilo[ix]   = cffrehilo[i][j];
         datarelolo[ix]   = cffrelolo[i][j];
         dataimhihi[ix]   = cffimhihi[i][j];
         dataimlohi[ix]   = cffimlohi[i][j];
         dataimhilo[ix]   = cffimhilo[i][j];
         dataimlolo[ix++] = cffimlolo[i][j];
      }
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++)
      {
         datarehihi[ix]   = inputrehihi[i][j];
         datarelohi[ix]   = inputrelohi[i][j];
         datarehilo[ix]   = inputrehilo[i][j];
         datarelolo[ix]   = inputrelolo[i][j];
         dataimhihi[ix]   = inputimhihi[i][j];
         dataimlohi[ix]   = inputimlohi[i][j];
         dataimhilo[ix]   = inputimhilo[i][j];
         dataimlolo[ix++] = inputimlolo[i][j];
      }

   for(int i=ix; i<totcff; i++)
   {
      datarehihi[i] = 0.0; datarelohi[i] = 0.0;
      datarehilo[i] = 0.0; datarelolo[i] = 0.0;
      dataimhihi[i] = 0.0; dataimlohi[i] = 0.0;
      dataimhilo[i] = 0.0; dataimlolo[i] = 0.0;
   }
}

void cmplx4vectorized_data_setup
 ( int dim, int nbr, int deg, int totcff, int offsetri,
   double *datarihihi, double *datarilohi,
   double *datarihilo, double *datarilolo,
   double *cstrehihi, double *cstrelohi,
   double *cstrehilo, double *cstrelolo,
   double *cstimhihi, double *cstimlohi,
   double *cstimhilo, double *cstimlolo,
   double **cffrehihi, double **cffrelohi,
   double **cffrehilo, double **cffrelolo,
   double **cffimhihi, double **cffimlohi,
   double **cffimhilo, double **cffimlolo,
   double **inputrehihi, double **inputrelohi,
   double **inputrehilo, double **inputrelolo,
   double **inputimhihi, double **inputimlohi,
   double **inputimhilo, double **inputimlolo )
{
   const int deg1 = deg+1;

   int ix1 = 0;
   int ix2 = totcff + offsetri;

   for(int i=0; i<deg1; i++)
   {
      datarihihi[ix1]   = cstrehihi[i];
      datarilohi[ix1]   = cstrelohi[i];
      datarihilo[ix1]   = cstrehilo[i];
      datarilolo[ix1++] = cstrelolo[i];
      datarihihi[ix2]   = cstimhihi[i];
      datarilohi[ix2]   = cstimlohi[i];
      datarihilo[ix2]   = cstimhilo[i];
      datarilolo[ix2++] = cstimlolo[i];
   }
   for(int i=0; i<nbr; i++)
      for(int j=0; j<deg1; j++)
      {
         datarihihi[ix1]   = cffrehihi[i][j];
         datarilohi[ix1]   = cffrelohi[i][j];
         datarihilo[ix1]   = cffrehilo[i][j];
         datarilolo[ix1++] = cffrelolo[i][j];
         datarihihi[ix2]   = cffimhihi[i][j];
         datarilohi[ix2]   = cffimlohi[i][j];
         datarihilo[ix2]   = cffimhilo[i][j];
         datarilolo[ix2++] = cffimlolo[i][j];
      }
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++)
      {
         datarihihi[ix1]   = inputrehihi[i][j];
         datarilohi[ix1]   = inputrelohi[i][j];
         datarihilo[ix1]   = inputrehilo[i][j];
         datarilolo[ix1++] = inputrelolo[i][j];
         datarihihi[ix2]   = inputimhihi[i][j];
         datarilohi[ix2]   = inputimlohi[i][j];
         datarihilo[ix2]   = inputimhilo[i][j];
         datarilolo[ix2++] = inputimlolo[i][j];
      }

   for(int i=0; i<2*offsetri; i++)
   {
      datarihihi[ix1]   = 0.0;
      datarilohi[ix1]   = 0.0;
      datarihilo[ix1]   = 0.0;
      datarilolo[ix1++] = 0.0;
      datarihihi[ix2]   = 0.0;
      datarilohi[ix2]   = 0.0;
      datarihilo[ix2]   = 0.0;
      datarilolo[ix2++] = 0.0;
   }
}

void dbl4_convolution_jobs
 ( int dim, int nbr, int deg, int *nvr, ConvolutionJobs cnvjobs,
   int *fstart, int *bstart, int *cstart,
   double *datahihi, double *datalohi, double *datahilo, double *datalolo,
   double *cnvlapms, int vrblvl )
{
   const int deg1 = deg+1;

   cudaEvent_t start,stop;           // to measure time spent by kernels
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   *cnvlapms = 0.0;
   float milliseconds;
   const bool vrb = (vrblvl > 1);

   for(int k=0; k<cnvjobs.get_depth(); k++)
   {
      const int jobnbr = cnvjobs.get_layer_count(k);
      int *in1ix_h = new int[jobnbr];
      int *in2ix_h = new int[jobnbr];
      int *outix_h = new int[jobnbr];

      if(vrb) cout << "preparing convolution jobs at layer "
                   << k << " ..." << endl;

      convjobs_coordinates(cnvjobs,k,in1ix_h,in2ix_h,outix_h,dim,nbr,deg,nvr,
                           fstart,bstart,cstart,vrb);
      // if(deg1 == BS)
      {
         int *in1ix_d; // first input on device
         int *in2ix_d; // second input on device
         int *outix_d; // output indices on device
         const size_t szjobidx = jobnbr*sizeof(int);
         cudaMalloc((void**)&in1ix_d,szjobidx);
         cudaMalloc((void**)&in2ix_d,szjobidx);
         cudaMalloc((void**)&outix_d,szjobidx);
         cudaMemcpy(in1ix_d,in1ix_h,szjobidx,cudaMemcpyHostToDevice);
         cudaMemcpy(in2ix_d,in2ix_h,szjobidx,cudaMemcpyHostToDevice);
         cudaMemcpy(outix_d,outix_h,szjobidx,cudaMemcpyHostToDevice);

         if(vrb)
            cout << "launching " << jobnbr << " blocks of " << deg1
                 << " threads for convolutions ..." << endl;

         cudaEventRecord(start);
         dbl4_padded_convjobs<<<jobnbr,deg1>>>
            (datahihi,datalohi,datahilo,datalolo,
             in1ix_d,in2ix_d,outix_d,deg1);
         cudaEventRecord(stop);
         cudaEventSynchronize(stop);
         cudaEventElapsedTime(&milliseconds,start,stop);
         *cnvlapms += milliseconds;
         cudaFree(in1ix_d); cudaFree(in2ix_d); cudaFree(outix_d);
      }
      free(in1ix_h); free(in2ix_h); free(outix_h);
   }
}

void cmplx4_convolution_jobs
 ( int dim, int nbr, int deg, int *nvr, ConvolutionJobs cnvjobs,
   int *fstart, int *bstart, int *cstart,
   double *datarehihi, double *datarelohi,
   double *datarehilo, double *datarelolo,
   double *dataimhihi, double *dataimlohi,
   double *dataimhilo, double *dataimlolo,
   double *cnvlapms, int vrblvl )
{
   const int deg1 = deg+1;

   cudaEvent_t start,stop;           // to measure time spent by kernels
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   *cnvlapms = 0.0;
   float milliseconds;
   const bool vrb = (vrblvl > 1);

   for(int k=0; k<cnvjobs.get_depth(); k++)
   {
      const int jobnbr = cnvjobs.get_layer_count(k);
      int *in1ix_h = new int[jobnbr];
      int *in2ix_h = new int[jobnbr];
      int *outix_h = new int[jobnbr];

      if(vrb) cout << "preparing convolution jobs at layer "
                   << k << " ..." << endl;

      convjobs_coordinates(cnvjobs,k,in1ix_h,in2ix_h,outix_h,dim,nbr,deg,nvr,
                           fstart,bstart,cstart,vrb);
      // if(deg1 == BS)
      {
         int *in1ix_d; // first input on device
         int *in2ix_d; // second input on device
         int *outix_d; // output indices on device
         const size_t szjobidx = jobnbr*sizeof(int);
         cudaMalloc((void**)&in1ix_d,szjobidx);
         cudaMalloc((void**)&in2ix_d,szjobidx);
         cudaMalloc((void**)&outix_d,szjobidx);
         cudaMemcpy(in1ix_d,in1ix_h,szjobidx,cudaMemcpyHostToDevice);
         cudaMemcpy(in2ix_d,in2ix_h,szjobidx,cudaMemcpyHostToDevice);
         cudaMemcpy(outix_d,outix_h,szjobidx,cudaMemcpyHostToDevice);

         if(vrb)
            cout << "launching " << jobnbr << " blocks of " << deg1
                 << " threads for convolutions ..." << endl;

         cudaEventRecord(start);
         cmplx4_padded_convjobs<<<jobnbr,deg1>>>
            (datarehihi,datarelohi,datarehilo,datarelolo,
             dataimhihi,dataimlohi,dataimhilo,dataimlolo,
             in1ix_d,in2ix_d,outix_d,deg1);
         cudaEventRecord(stop);
         cudaEventSynchronize(stop);
         cudaEventElapsedTime(&milliseconds,start,stop);
         *cnvlapms += milliseconds;
         cudaFree(in1ix_d); cudaFree(in2ix_d); cudaFree(outix_d);
      }
      free(in1ix_h); free(in2ix_h); free(outix_h);
   }
}

void cmplx4vectorized_convolution_jobs
 ( int dim, int nbr, int deg, int *nvr, int totcff, int offsetri,
   ComplexConvolutionJobs cnvjobs, ComplexIncrementJobs incjobs,
   int *fstart, int *bstart, int *cstart,
   double *datarihihi, double *datarilohi,
   double *datarihilo, double *datarilolo,
   double *cnvlapms, int vrblvl )
{
   const int deg1 = deg+1;

   cudaEvent_t start,stop;           // to measure time spent by kernels
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   *cnvlapms = 0.0;
   float milliseconds;
   double fliplapms;
   const bool vrb = (vrblvl > 1);

   for(int k=0; k<cnvjobs.get_depth(); k++)
   {
      int jobnbr = cnvjobs.get_layer_count(k);
      int *in1ix_h = new int[jobnbr];
      int *in2ix_h = new int[jobnbr];
      int *outix_h = new int[jobnbr];

      if(vrb) cout << "preparing convolution jobs at layer "
                   << k << " ..." << endl;

      complex_convjobs_coordinates
         (cnvjobs,k,in1ix_h,in2ix_h,outix_h,dim,nbr,deg,nvr,totcff,offsetri,
          fstart,bstart,cstart,vrb);

      // if(deg1 == BS)
      {
         int *in1ix_d; // first input on device
         int *in2ix_d; // second input on device
         int *outix_d; // output indices on device
         const size_t szjobidx = jobnbr*sizeof(int);
         cudaMalloc((void**)&in1ix_d,szjobidx);
         cudaMalloc((void**)&in2ix_d,szjobidx);
         cudaMalloc((void**)&outix_d,szjobidx);
         cudaMemcpy(in1ix_d,in1ix_h,szjobidx,cudaMemcpyHostToDevice);
         cudaMemcpy(in2ix_d,in2ix_h,szjobidx,cudaMemcpyHostToDevice);
         cudaMemcpy(outix_d,outix_h,szjobidx,cudaMemcpyHostToDevice);

         if(vrb)
            cout << "launching " << jobnbr << " blocks of " << deg1
                 << " threads for convolutions ..." << endl;

         cudaEventRecord(start);
         dbl4_padded_convjobs<<<jobnbr,deg1>>>
            (datarihihi,datarilohi,datarihilo,datarilolo,
             in1ix_d,in2ix_d,outix_d,deg1);
         cudaEventRecord(stop);
         cudaEventSynchronize(stop);
         cudaEventElapsedTime(&milliseconds,start,stop);
         *cnvlapms += milliseconds;
         cudaFree(in1ix_d); cudaFree(in2ix_d); cudaFree(outix_d);
      }
      jobnbr = incjobs.get_layer_count(k);
      // note: only half the number of increment jobs

      if(vrb) cout << "preparing increment jobs at layer "
                   << k << " ..." << endl;

      complex_incjobs_coordinates
         (incjobs,k,in1ix_h,in2ix_h,outix_h,dim,nbr,deg,nvr,totcff,offsetri,
          fstart,bstart,cstart,vrb);

      const int nbrflips = jobnbr/2;
      int *rebidx = new int[nbrflips];
      for(int i=0, j=0; i<jobnbr; i=i+2, j++) rebidx[j] = in2ix_h[i];

      GPU_cmplx4vectorized_flipsigns
        (deg,nbrflips,rebidx,datarihihi,datarilohi,datarihilo,datarilolo,
         &fliplapms,vrblvl);
      *cnvlapms += fliplapms;

      // if(BS == deg1)
      {
         int *in1ix_d; // first input on device
         int *in2ix_d; // second input on device
         int *outix_d; // output indices on device
         const size_t szjobidx = jobnbr*sizeof(int);
         cudaMalloc((void**)&in1ix_d,szjobidx);
         cudaMalloc((void**)&in2ix_d,szjobidx);
         cudaMalloc((void**)&outix_d,szjobidx);
         cudaMemcpy(in1ix_d,in1ix_h,szjobidx,cudaMemcpyHostToDevice);
         cudaMemcpy(in2ix_d,in2ix_h,szjobidx,cudaMemcpyHostToDevice);
         cudaMemcpy(outix_d,outix_h,szjobidx,cudaMemcpyHostToDevice);

         if(vrb)
            cout << "launching " << jobnbr << " blocks of " << deg1
                 << " threads for increments ..." << endl;

         cudaEventRecord(start);
         dbl4_increment_jobs<<<jobnbr,deg1>>>
            (datarihihi,datarilohi,datarihilo,datarilolo,
             in1ix_d,in2ix_d,outix_d,deg1);
         cudaEventRecord(stop);
         cudaEventSynchronize(stop);
         cudaEventElapsedTime(&milliseconds,start,stop);
         *cnvlapms += milliseconds;
         cudaFree(in1ix_d); cudaFree(in2ix_d); cudaFree(outix_d);
      }
      free(in1ix_h); free(in2ix_h); free(outix_h);
      free(rebidx);
   }
}

void dbl4_addition_jobs
 ( int dim, int nbr, int deg, int *nvr, AdditionJobs addjobs,
   int *fstart, int *bstart, int *cstart,
   double *datahihi, double *datalohi, double *datahilo, double *datalolo,
   double *addlapms, int vrblvl )
{
   const int deg1 = deg+1;

   cudaEvent_t start,stop;           // to measure time spent by kernels
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   *addlapms = 0.0;
   float milliseconds;
   const bool vrb = (vrblvl > 1);

   for(int k=0; k<addjobs.get_depth(); k++)
   {
      const int jobnbr = addjobs.get_layer_count(k);
      int *in1ix_h = new int[jobnbr];
      int *in2ix_h = new int[jobnbr];
      int *outix_h = new int[jobnbr];

      if(vrb) cout << "preparing addition jobs at layer "
                   << k << " ..." << endl;

      addjobs_coordinates(addjobs,k,in1ix_h,in2ix_h,outix_h,dim,nbr,deg,nvr,
                          fstart,bstart,cstart,vrb);
      // if(deg1 == BS)
      {
         int *in1ix_d; // first input on device
         int *in2ix_d; // second input on device
         int *outix_d; // output indices on device
         const size_t szjobidx = jobnbr*sizeof(int);
         cudaMalloc((void**)&in1ix_d,szjobidx);
         cudaMalloc((void**)&in2ix_d,szjobidx);
         cudaMalloc((void**)&outix_d,szjobidx);
         cudaMemcpy(in1ix_d,in1ix_h,szjobidx,cudaMemcpyHostToDevice);
         cudaMemcpy(in2ix_d,in2ix_h,szjobidx,cudaMemcpyHostToDevice);
         cudaMemcpy(outix_d,outix_h,szjobidx,cudaMemcpyHostToDevice);

         if(vrb)
            cout << "launching " << jobnbr << " blocks of " << deg1
                 << " threads for additions ..." << endl;

         cudaEventRecord(start);
         dbl4_update_addjobs<<<jobnbr,deg1>>>
            (datahihi,datalohi,datahilo,datalolo,
             in1ix_d,in2ix_d,outix_d,deg1);
         cudaEventRecord(stop);
         cudaEventSynchronize(stop);
         cudaEventElapsedTime(&milliseconds,start,stop);
         *addlapms += milliseconds;
         cudaFree(in1ix_d); cudaFree(in2ix_d); cudaFree(outix_d);
      }
      free(in1ix_h); free(in2ix_h); free(outix_h);
   }
}

void cmplx4_addition_jobs
 ( int dim, int nbr, int deg, int *nvr, AdditionJobs addjobs,
   int *fstart, int *bstart, int *cstart,
   double *datarehihi, double *datarelohi,
   double *datarehilo, double *datarelolo,
   double *dataimhihi, double *dataimlohi,
   double *dataimhilo, double *dataimlolo,
   double *addlapms, int vrblvl )
{
   const int deg1 = deg+1;

   cudaEvent_t start,stop;           // to measure time spent by kernels
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   *addlapms = 0.0;
   float milliseconds;
   const bool vrb = (vrblvl > 1);

   for(int k=0; k<addjobs.get_depth(); k++)
   {
      const int jobnbr = addjobs.get_layer_count(k);
      int *in1ix_h = new int[jobnbr];
      int *in2ix_h = new int[jobnbr];
      int *outix_h = new int[jobnbr];

      if(vrb) cout << "preparing addition jobs at layer "
                   << k << " ..." << endl;

      addjobs_coordinates(addjobs,k,in1ix_h,in2ix_h,outix_h,dim,nbr,deg,nvr,
                          fstart,bstart,cstart,vrb);
      // if(deg1 == BS)
      {
         int *in1ix_d; // first input on device
         int *in2ix_d; // second input on device
         int *outix_d; // output indices on device
         const size_t szjobidx = jobnbr*sizeof(int);
         cudaMalloc((void**)&in1ix_d,szjobidx);
         cudaMalloc((void**)&in2ix_d,szjobidx);
         cudaMalloc((void**)&outix_d,szjobidx);
         cudaMemcpy(in1ix_d,in1ix_h,szjobidx,cudaMemcpyHostToDevice);
         cudaMemcpy(in2ix_d,in2ix_h,szjobidx,cudaMemcpyHostToDevice);
         cudaMemcpy(outix_d,outix_h,szjobidx,cudaMemcpyHostToDevice);

         if(vrb)
            cout << "launching " << jobnbr << " blocks of " << deg1
                 << " threads for additions ..." << endl;

         cudaEventRecord(start);
         cmplx4_update_addjobs<<<jobnbr,deg1>>>
            (datarehihi,datarelohi,datarehilo,datarelolo,
             dataimhihi,dataimlohi,dataimhilo,dataimlolo,
             in1ix_d,in2ix_d,outix_d,deg1);
         cudaEventRecord(stop);
         cudaEventSynchronize(stop);
         cudaEventElapsedTime(&milliseconds,start,stop);
         *addlapms += milliseconds;
         cudaFree(in1ix_d); cudaFree(in2ix_d); cudaFree(outix_d);
      }
      free(in1ix_h); free(in2ix_h); free(outix_h);
   }
}

void cmplx4vectorized_addition_jobs
 ( int dim, int nbr, int deg, int *nvr, int totcff, int offsetri,
   ComplexAdditionJobs addjobs,
   int *fstart, int *bstart, int *cstart,
   double *datarihihi, double *datarilohi,
   double *datarihilo, double *datarilolo,
   double *addlapms, int vrblvl )
{
   const int deg1 = deg+1;

   cudaEvent_t start,stop;           // to measure time spent by kernels
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   *addlapms = 0.0;
   float milliseconds;
   const bool vrb = (vrblvl > 1);

   for(int k=0; k<addjobs.get_depth(); k++)
   {
      const int jobnbr = addjobs.get_layer_count(k);
      int *in1ix_h = new int[jobnbr];
      int *in2ix_h = new int[jobnbr];
      int *outix_h = new int[jobnbr];

      if(vrb) cout << "preparing addition jobs at layer "
                   << k << " ..." << endl;

      complex_addjobs_coordinates
         (addjobs,k,in1ix_h,in2ix_h,outix_h,dim,nbr,deg,nvr,
          totcff,offsetri,fstart,bstart,cstart,vrb);

      // if(deg1 == BS)
      {
         int *in1ix_d; // first input on device
         int *in2ix_d; // second input on device
         int *outix_d; // output indices on device
         const size_t szjobidx = jobnbr*sizeof(int);
         cudaMalloc((void**)&in1ix_d,szjobidx);
         cudaMalloc((void**)&in2ix_d,szjobidx);
         cudaMalloc((void**)&outix_d,szjobidx);
         cudaMemcpy(in1ix_d,in1ix_h,szjobidx,cudaMemcpyHostToDevice);
         cudaMemcpy(in2ix_d,in2ix_h,szjobidx,cudaMemcpyHostToDevice);
         cudaMemcpy(outix_d,outix_h,szjobidx,cudaMemcpyHostToDevice);

         if(vrb)
            cout << "launching " << jobnbr << " blocks of " << deg1
                 << " threads for additions ..." << endl;

         cudaEventRecord(start);
         dbl4_update_addjobs<<<jobnbr,deg1>>>
            (datarihihi,datarilohi,datarihilo,datarilolo,
             in1ix_d,in2ix_d,outix_d,deg1);
         cudaEventRecord(stop);
         cudaEventSynchronize(stop);
         cudaEventElapsedTime(&milliseconds,start,stop);
         *addlapms += milliseconds;
         cudaFree(in1ix_d); cudaFree(in2ix_d); cudaFree(outix_d);
      }
      free(in1ix_h); free(in2ix_h); free(outix_h);
   }
}

void GPU_dbl4_poly_evaldiff
 ( int BS, int dim, int nbr, int deg, int *nvr, int **idx,
   double *csthihi, double *cstlohi, double *csthilo, double *cstlolo,
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo,
   double **inputhihi, double **inputlohi,
   double **inputhilo, double **inputlolo,
   double **outputhihi, double **outputlohi,
   double **outputhilo, double **outputlolo,
   ConvolutionJobs cnvjobs, AdditionJobs addjobs,
   double *cnvlapms, double *addlapms, double *elapsedms,
   double *walltimesec, int vrblvl )
{
   const int totalcff = coefficient_count(dim,nbr,deg,nvr);

   int *fstart = new int[nbr];
   int *bstart = new int[nbr];
   int *cstart = new int[nbr];
   int *fsums = new int[nbr];
   int *bsums = new int[nbr];
   int *csums = new int[nbr];

   coefficient_indices
      (dim,nbr,deg,nvr,fsums,bsums,csums,fstart,bstart,cstart);
   if(vrblvl > 1)
      write_coefficient_indices
         (totalcff,nbr,fsums,fstart,bsums,bstart,csums,cstart);

   double *datahihi_h = new double[totalcff];        // data on host
   double *datalohi_h = new double[totalcff];
   double *datahilo_h = new double[totalcff];
   double *datalolo_h = new double[totalcff];

   dbl4_data_setup
      (dim,nbr,deg,totalcff,datahihi_h,datalohi_h,datahilo_h,datalolo_h,
       csthihi,cstlohi,csthilo,cstlolo,cffhihi,cfflohi,cffhilo,cfflolo,
       inputhihi,inputlohi,inputhilo,inputlolo);

   double *datahihi_d;                               // device data
   double *datalohi_d;
   double *datahilo_d;
   double *datalolo_d;
   const size_t szdata = totalcff*sizeof(double);
   cudaMalloc((void**)&datahihi_d,szdata);
   cudaMalloc((void**)&datalohi_d,szdata);
   cudaMalloc((void**)&datahilo_d,szdata);
   cudaMalloc((void**)&datalolo_d,szdata);
   cudaMemcpy(datahihi_d,datahihi_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datalohi_d,datalohi_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datahilo_d,datahilo_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datalolo_d,datalolo_h,szdata,cudaMemcpyHostToDevice);

   struct timeval begintime,endtime; // wall clock time of computations
   gettimeofday(&begintime,0);

   dbl4_convolution_jobs
      (dim,nbr,deg,nvr,cnvjobs,fstart,bstart,cstart,
       datahihi_d,datalohi_d,datahilo_d,datalolo_d,cnvlapms,vrblvl);

   dbl4_addition_jobs
      (dim,nbr,deg,nvr,addjobs,fstart,bstart,cstart,
       datahihi_d,datalohi_d,datahilo_d,datalolo_d,addlapms,vrblvl);

   gettimeofday(&endtime,0);
   cudaMemcpy(datahihi_h,datahihi_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datalohi_h,datalohi_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datahilo_h,datahilo_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datalolo_h,datalolo_d,szdata,cudaMemcpyDeviceToHost);
   *elapsedms = *cnvlapms + *addlapms;
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   *walltimesec = seconds + microseconds*1.0e-6;

   // dbl_convoluted_data2_to_output
   //   (datahihi_h,datalohi_h,datahilo_h,datalolo_h,
   //    outputhihi,outputlohi,outputhilo,outputlolo,
   //    dim,nbr,deg,nvr,idx,fstart,bstart,cstart,vrblvl);
   dbl_added_data4_to_output
      (datahihi_h,datalohi_h,datahilo_h,datalolo_h,
       outputhihi,outputlohi,outputhilo,outputlolo,
       dim,nbr,deg,nvr,idx,fstart,bstart,cstart,addjobs,vrblvl);

   if(vrblvl > 0)
   {
      cout << fixed << setprecision(2);
      cout << "Time spent by convolution kernels : ";
      cout << *cnvlapms << " milliseconds." << endl;
      cout << "Time spent by addition kernels    : ";
      cout << *addlapms << " milliseconds." << endl;
      cout << "Time spent by all kernels         : ";
      cout << *elapsedms << " milliseconds." << endl;
      cout << "Total wall clock computation time : ";
      cout << fixed << setprecision(3) << *walltimesec
           << " seconds." << endl;
      cout << scientific << setprecision(16);
   }
   cudaFree(datahihi_d); cudaFree(datalohi_d);
   cudaFree(datahilo_d); cudaFree(datalolo_d);

   free(datahihi_h); free(datalohi_h);
   free(datahilo_h); free(datalolo_h);

   free(fstart); free(bstart); free(cstart);
   free(fsums); free(bsums); free(csums);
}

void GPU_cmplx4_poly_evaldiff
 ( int BS, int dim, int nbr, int deg, int *nvr, int **idx,
   double *cstrehihi, double *cstrelohi,
   double *cstrehilo, double *cstrelolo,
   double *cstimhihi, double *cstimlohi,
   double *cstimhilo, double *cstimlolo,
   double **cffrehihi, double **cffrelohi,
   double **cffrehilo, double **cffrelolo,
   double **cffimhihi, double **cffimlohi,
   double **cffimhilo, double **cffimlolo,
   double **inputrehihi, double **inputrelohi,
   double **inputrehilo, double **inputrelolo,
   double **inputimhihi, double **inputimlohi,
   double **inputimhilo, double **inputimlolo,
   double **outputrehihi, double **outputrelohi,
   double **outputrehilo, double **outputrelolo,
   double **outputimhihi, double **outputimlohi,
   double **outputimhilo, double **outputimlolo,
   ConvolutionJobs cnvjobs, AdditionJobs addjobs,
   double *cnvlapms, double *addlapms, double *elapsedms,
   double *walltimesec, int vrblvl )
{
   const int totalcff = coefficient_count(dim,nbr,deg,nvr);

   int *fstart = new int[nbr];
   int *bstart = new int[nbr];
   int *cstart = new int[nbr];
   int *fsums = new int[nbr];
   int *bsums = new int[nbr];
   int *csums = new int[nbr];

   coefficient_indices
      (dim,nbr,deg,nvr,fsums,bsums,csums,fstart,bstart,cstart);
   if(vrblvl > 1)
      write_coefficient_indices
         (totalcff,nbr,fsums,fstart,bsums,bstart,csums,cstart);

   double *datarehihi_h = new double[totalcff];      // data on host
   double *datarelohi_h = new double[totalcff];
   double *datarehilo_h = new double[totalcff];
   double *datarelolo_h = new double[totalcff];
   double *dataimhihi_h = new double[totalcff]; 
   double *dataimlohi_h = new double[totalcff]; 
   double *dataimhilo_h = new double[totalcff];
   double *dataimlolo_h = new double[totalcff];

   cmplx4_data_setup
      (dim,nbr,deg,totalcff,
       datarehihi_h,datarelohi_h,datarehilo_h,datarelolo_h,
       dataimhihi_h,dataimlohi_h,dataimhilo_h,dataimlolo_h,
       cstrehihi,cstrelohi,cstrehilo,cstrelolo,
       cstimhihi,cstimlohi,cstimhilo,cstimlolo,
       cffrehihi,cffrelohi,cffrehilo,cffrelolo,
       cffimhihi,cffimlohi,cffimhilo,cffimlolo,
       inputrehihi,inputrelohi,inputrehilo,inputrelolo,
       inputimhihi,inputimlohi,inputimhilo,inputimlolo);

   double *datarehihi_d;                               // device data
   double *datarelohi_d;
   double *datarehilo_d;
   double *datarelolo_d;
   double *dataimhihi_d;
   double *dataimlohi_d;
   double *dataimhilo_d;
   double *dataimlolo_d;
   const size_t szdata = totalcff*sizeof(double);
   cudaMalloc((void**)&datarehihi_d,szdata);
   cudaMalloc((void**)&datarelohi_d,szdata);
   cudaMalloc((void**)&datarehilo_d,szdata);
   cudaMalloc((void**)&datarelolo_d,szdata);
   cudaMalloc((void**)&dataimhihi_d,szdata);
   cudaMalloc((void**)&dataimlohi_d,szdata);
   cudaMalloc((void**)&dataimhilo_d,szdata);
   cudaMalloc((void**)&dataimlolo_d,szdata);
   cudaMemcpy(datarehihi_d,datarehihi_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datarelohi_d,datarelohi_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datarehilo_d,datarehilo_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datarelolo_d,datarelolo_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(dataimhihi_d,dataimhihi_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(dataimlohi_d,dataimlohi_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(dataimhilo_d,dataimhilo_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(dataimlolo_d,dataimlolo_h,szdata,cudaMemcpyHostToDevice);

   struct timeval begintime,endtime; // wall clock time of computations

   gettimeofday(&begintime,0);

   cmplx4_convolution_jobs
      (dim,nbr,deg,nvr,cnvjobs,fstart,bstart,cstart,
       datarehihi_d,datarelohi_d,datarehilo_d,datarelolo_d,
       dataimhihi_d,dataimlohi_d,dataimhilo_d,dataimlolo_d,
       cnvlapms,vrblvl);

   cmplx4_addition_jobs
      (dim,nbr,deg,nvr,addjobs,fstart,bstart,cstart,
       datarehihi_d,datarelohi_d,datarehilo_d,datarelolo_d,
       dataimhihi_d,dataimlohi_d,dataimhilo_d,dataimlolo_d,
       addlapms,vrblvl);

   gettimeofday(&endtime,0);
   cudaMemcpy(datarehihi_h,datarehihi_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datarelohi_h,datarelohi_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datarehilo_h,datarehilo_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datarelolo_h,datarelolo_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(dataimhihi_h,dataimhihi_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(dataimlohi_h,dataimlohi_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(dataimhilo_h,dataimhilo_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(dataimlolo_h,dataimlolo_d,szdata,cudaMemcpyDeviceToHost);
   *elapsedms = *cnvlapms + *addlapms;
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   *walltimesec = seconds + microseconds*1.0e-6;

   // cmplx_convoluted_data4_to_output
   //    (datarehihi_h,datarelohi_h,datarehilo_h,datarelolo_h,
   //     dataimhihi_h,dataimlohi_h,dataimhilo_h,dataimlolo_h,
   //     outputrehihi,outputrelohi,outputrehilo,outputrelolo,
   //     outputimhihi,outputimlohi,outputimhilo,outputimlolo,
   //     dim,nbr,deg,nvr,idx,fstart,bstart,cstart,vrblvl);
   cmplx_added_data4_to_output
      (datarehihi_h,datarelohi_h,datarehilo_h,datarelolo_h,
       dataimhihi_h,dataimlohi_h,dataimhilo_h,dataimlolo_h,
       outputrehihi,outputrelohi,outputrehilo,outputrelolo,
       outputimhihi,outputimlohi,outputimhilo,outputimlolo,
       dim,nbr,deg,nvr,idx,fstart,bstart,cstart,addjobs,vrblvl);

   if(vrblvl > 0)
      write_GPU_timings(*cnvlapms,*addlapms,*elapsedms,*walltimesec);

   cudaFree(datarehihi_d); cudaFree(datarelohi_d);
   cudaFree(datarehilo_d); cudaFree(datarelolo_d);
   cudaFree(dataimhihi_d); cudaFree(dataimlohi_d);
   cudaFree(dataimhilo_d); cudaFree(dataimlolo_d);

   free(datarehihi_h); free(datarelohi_h);
   free(datarehilo_h); free(datarelolo_h);
   free(dataimhihi_h); free(dataimlohi_h);
   free(dataimhilo_h); free(dataimlolo_h);

   free(fstart); free(bstart); free(cstart);
   free(fsums); free(bsums); free(csums);
}

void GPU_cmplx4vectorized_poly_evaldiff
 ( int BS, int dim, int nbr, int deg, int *nvr, int **idx,
   double *cstrehihi, double *cstrelohi,
   double *cstrehilo, double *cstrelolo,
   double *cstimhihi, double *cstimlohi,
   double *cstimhilo, double *cstimlolo,
   double **cffrehihi, double **cffrelohi,
   double **cffrehilo, double **cffrelolo,
   double **cffimhihi, double **cffimlohi,
   double **cffimhilo, double **cffimlolo,
   double **inputrehihi, double **inputrelohi,
   double **inputrehilo, double **inputrelolo,
   double **inputimhihi, double **inputimlohi,
   double **inputimhilo, double **inputimlolo,
   double **outputrehihi, double **outputrelohi,
   double **outputrehilo, double **outputrelolo,
   double **outputimhihi, double **outputimlohi,
   double **outputimhilo, double **outputimlolo,
   ComplexConvolutionJobs cnvjobs, ComplexIncrementJobs incjobs,
   ComplexAdditionJobs addjobs,
   double *cnvlapms, double *addlapms, double *elapsedms,
   double *walltimesec, int vrblvl )
{
   const int totalcff = complex_coefficient_count(dim,nbr,deg,nvr);
   const int diminput = (1 + nbr + dim)*(deg + 1); // dimension of input
   const int offsetri = totalcff - diminput; // offset for re/im operands
   const int cmplxtotcff = 2*(totalcff + offsetri);

   int *fstart = new int[nbr];
   int *bstart = new int[nbr];
   int *cstart = new int[nbr];
   int *fsums = new int[nbr];
   int *bsums = new int[nbr];
   int *csums = new int[nbr];

   complex_coefficient_indices
      (dim,nbr,deg,nvr,fsums,bsums,csums,fstart,bstart,cstart);

   if(vrblvl > 1)
   {
      cout << "        total count : " << totalcff << endl;
      cout << "offset for operands : " << offsetri << endl;
      cout << "complex total count : " << cmplxtotcff << endl;
      write_coefficient_indices
         (totalcff,nbr,fsums,fstart,bsums,bstart,csums,cstart);
   }
   double *datarihihi_h = new double[cmplxtotcff];      // data on host
   double *datarilohi_h = new double[cmplxtotcff];
   double *datarihilo_h = new double[cmplxtotcff];
   double *datarilolo_h = new double[cmplxtotcff];

   cmplx4vectorized_data_setup
      (dim,nbr,deg,totalcff,offsetri,
       datarihihi_h,datarilohi_h,datarihilo_h,datarilolo_h,
       cstrehihi,cstrelohi,cstrehilo,cstrelolo,
       cstimhihi,cstimlohi,cstimhilo,cstimlolo,
       cffrehihi,cffrelohi,cffrehilo,cffrelolo,
       cffimhihi,cffimlohi,cffimhilo,cffimlolo,
       inputrehihi,inputrelohi,inputrehilo,inputrelolo,
       inputimhihi,inputimlohi,inputimhilo,inputimlolo);

   double *datarihihi_d;                               // device data
   double *datarilohi_d;
   double *datarihilo_d;
   double *datarilolo_d;
   const size_t szdata = cmplxtotcff*sizeof(double);
   cudaMalloc((void**)&datarihihi_d,szdata);
   cudaMalloc((void**)&datarilohi_d,szdata);
   cudaMalloc((void**)&datarihilo_d,szdata);
   cudaMalloc((void**)&datarilolo_d,szdata);
   cudaMemcpy(datarihihi_d,datarihihi_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datarilohi_d,datarilohi_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datarihilo_d,datarihilo_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datarilolo_d,datarilolo_h,szdata,cudaMemcpyHostToDevice);

   struct timeval begintime,endtime; // wall clock time of computations

   gettimeofday(&begintime,0);

   cmplx4vectorized_convolution_jobs
      (dim,nbr,deg,nvr,totalcff,offsetri,cnvjobs,incjobs,fstart,bstart,cstart,
       datarihihi_d,datarilohi_d,datarihilo_d,datarilolo_d,
       cnvlapms,vrblvl);

   cmplx4vectorized_addition_jobs
      (dim,nbr,deg,nvr,totalcff,offsetri,addjobs,fstart,bstart,cstart,
       datarihihi_d,datarilohi_d,datarihilo_d,datarilolo_d,
       addlapms,vrblvl);

   gettimeofday(&endtime,0);
   cudaMemcpy(datarihihi_h,datarihihi_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datarilohi_h,datarilohi_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datarihilo_h,datarihilo_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datarilolo_h,datarilolo_d,szdata,cudaMemcpyDeviceToHost);
   *elapsedms = *cnvlapms + *addlapms;
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   *walltimesec = seconds + microseconds*1.0e-6;

   cmplx_added_data4vectorized_to_output
      (datarihihi_h,datarilohi_h,datarihilo_h,datarilolo_h,
       outputrehihi,outputrelohi,outputrehilo,outputrelolo,
       outputimhihi,outputimlohi,outputimhilo,outputimlolo,
       dim,nbr,deg,nvr,idx,fstart,bstart,cstart,
       totalcff,offsetri,addjobs,vrblvl);

   if(vrblvl > 0)
      write_GPU_timings(*cnvlapms,*addlapms,*elapsedms,*walltimesec);

   cudaFree(datarihihi_d); cudaFree(datarilohi_d);
   cudaFree(datarihilo_d); cudaFree(datarilolo_d);

   free(datarihihi_h); free(datarilohi_h);
   free(datarihilo_h); free(datarilolo_h);

   free(fstart); free(bstart); free(cstart);
   free(fsums); free(bsums); free(csums);
}
