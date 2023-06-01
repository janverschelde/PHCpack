// The file dbl8_polynomials_kernels.cu defines the kernels with prototypes
// in dbl8_polynomials_kernels.h.

#include <iostream>
#include <iomanip>
#ifdef winwalltime
#include "gettimeofday4win.h"
#else
#include <sys/time.h>
#endif
#include "job_coordinates.h"
#include "octo_double_functions.h"
#ifdef gpufun
#include "double_double_gpufun.cu"
#include "quad_double_gpufun.cu"
#include "octo_double_gpufun.cu"
#endif
#include "dbl8_polynomials_kernels.h"
#include "write_gpu_timings.h"

// The constant od_shmemsize is the bound on the shared memory size.

#define od_shmemsize 192
#define odh_shmemsize 96

// odh_shmesize is used on complex data, otherwise too much used ...

using namespace std;

__global__ void dbl8_padded_convjobs
 ( double *datahihihi, double *datalohihi,
   double *datahilohi, double *datalolohi,
   double *datahihilo, double *datalohilo,
   double *datahilolo, double *datalololo,
   int *in1idx, int *in2idx, int *outidx, int dim )
{
   const int bdx = blockIdx.x;           // index to the convolution job
   const int tdx = threadIdx.x;          // index to the output of the job
   const int idx1 = in1idx[bdx] + tdx;
   const int idx2 = in2idx[bdx] + tdx;
   const int idx3 = outidx[bdx] + tdx;

   __shared__ double xvhihihi[od_shmemsize];
   __shared__ double xvlohihi[od_shmemsize];
   __shared__ double xvhilohi[od_shmemsize];
   __shared__ double xvlolohi[od_shmemsize];
   __shared__ double xvhihilo[od_shmemsize];
   __shared__ double xvlohilo[od_shmemsize];
   __shared__ double xvhilolo[od_shmemsize];
   __shared__ double xvlololo[od_shmemsize];
   __shared__ double yvhihihi[2*od_shmemsize];
   __shared__ double yvlohihi[2*od_shmemsize];
   __shared__ double yvhilohi[2*od_shmemsize];
   __shared__ double yvlolohi[2*od_shmemsize];
   __shared__ double yvhihilo[2*od_shmemsize];
   __shared__ double yvlohilo[2*od_shmemsize];
   __shared__ double yvhilolo[2*od_shmemsize];
   __shared__ double yvlololo[2*od_shmemsize];
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
   int ydx = dim + tdx;

   xvhihihi[tdx] = datahihihi[idx1];  // loading first input
   xvlohihi[tdx] = datalohihi[idx1]; 
   xvhilohi[tdx] = datahilohi[idx1]; 
   xvlolohi[tdx] = datalolohi[idx1]; 
   xvhihilo[tdx] = datahihilo[idx1];
   xvlohilo[tdx] = datalohilo[idx1]; 
   xvhilolo[tdx] = datahilolo[idx1]; 
   xvlololo[tdx] = datalololo[idx1]; 
   yvhihihi[tdx] = 0.0;             // padded with zeros
   yvlohihi[tdx] = 0.0;
   yvhilohi[tdx] = 0.0;
   yvlolohi[tdx] = 0.0;
   yvhihilo[tdx] = 0.0;
   yvlohilo[tdx] = 0.0;
   yvhilolo[tdx] = 0.0;
   yvlololo[tdx] = 0.0;
   yvhihihi[ydx] = datahihihi[idx2];  // loading second input
   yvlohihi[ydx] = datalohihi[idx2];
   yvhilohi[ydx] = datahilohi[idx2];
   yvlolohi[ydx] = datalolohi[idx2];
   yvhihilo[ydx] = datahihilo[idx2];
   yvlohilo[ydx] = datalohilo[idx2];
   yvhilolo[ydx] = datahilolo[idx2];
   yvlololo[ydx] = datalololo[idx2];

   __syncthreads();

   // zv[tdx] = xv[0]*yv[tdx];
   odg_mul( xvhihihi[0],   xvlohihi[0],   xvhilohi[0],   xvlolohi[0],
            xvhihilo[0],   xvlohilo[0],   xvhilolo[0],   xvlololo[0],
            yvhihihi[ydx], yvlohihi[ydx], yvhilohi[ydx], yvlolohi[ydx],
            yvhihilo[ydx], yvlohilo[ydx], yvhilolo[ydx], yvlololo[ydx],
           &zvhihihi[tdx],&zvlohihi[tdx],&zvhilohi[tdx],&zvlolohi[tdx],
           &zvhihilo[tdx],&zvlohilo[tdx],&zvhilolo[tdx],&zvlololo[tdx]);
   __syncthreads();

   for(int i=1; i<dim; i++) // zv[tdx] = zv[tdx] + xv[i]*yv[dim+tdx-i];
   {
      ydx = dim + tdx - i;

      odg_mul( xvhihihi[i],  xvlohihi[i],  xvhilohi[i],  xvlolohi[i],
               xvhihilo[i],  xvlohilo[i],  xvhilolo[i],  xvlololo[i],
               yvhihihi[ydx],yvlohihi[ydx],yvhilohi[ydx],yvlolohi[ydx],
               yvhihilo[ydx],yvlohilo[ydx],yvhilolo[ydx],yvlololo[ydx],
             &prdhihihi,   &prdlohihi,   &prdhilohi,   &prdlolohi,
             &prdhihilo,   &prdlohilo,   &prdhilolo,   &prdlololo);
      __syncthreads();

      odg_inc(&zvhihihi[tdx],&zvlohihi[tdx],&zvhilohi[tdx],&zvlolohi[tdx],
              &zvhihilo[tdx],&zvlohilo[tdx],&zvhilolo[tdx],&zvlololo[tdx],
              prdhihihi,     prdlohihi,     prdhilohi,     prdlolohi,
              prdhihilo,     prdlohilo,     prdhilolo,     prdlololo);
      __syncthreads();
   }
   __syncthreads();

   datahihihi[idx3] = zvhihihi[tdx]; // storing the output
   datalohihi[idx3] = zvlohihi[tdx];
   datahilohi[idx3] = zvhilohi[tdx];
   datalolohi[idx3] = zvlolohi[tdx];
   datahihilo[idx3] = zvhihilo[tdx];
   datalohilo[idx3] = zvlohilo[tdx];
   datahilolo[idx3] = zvhilolo[tdx];
   datalololo[idx3] = zvlololo[tdx];
}

__global__ void dbl8_increment_jobs
 ( double *datahihihi, double *datalohihi,
   double *datahilohi, double *datalolohi,
   double *datahihilo, double *datalohilo,
   double *datahilolo, double *datalololo,
   int *in1idx, int *in2idx, int *outidx, int dim )
{
   const int bdx = blockIdx.x;           // index to the increment job
   const int tdx = threadIdx.x;          // index to the output of the job
   const int idx1 = in1idx[bdx] + tdx;
   const int idx2 = in2idx[bdx] + tdx;
   const int idx3 = outidx[bdx] + tdx;

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

   zvhihihi[tdx] = datahihihi[idx1];  // loading first input
   zvlohihi[tdx] = datalohihi[idx1]; 
   zvhilohi[tdx] = datahilohi[idx1]; 
   zvlolohi[tdx] = datalolohi[idx1]; 
   zvhihilo[tdx] = datahihilo[idx1];
   zvlohilo[tdx] = datalohilo[idx1]; 
   zvhilolo[tdx] = datahilolo[idx1]; 
   zvlololo[tdx] = datalololo[idx1]; 
   yvhihihi[tdx] = datahihihi[idx2];  // loading second input
   yvlohihi[tdx] = datalohihi[idx2];
   yvhilohi[tdx] = datahilohi[idx2];
   yvlolohi[tdx] = datalolohi[idx2];
   yvhihilo[tdx] = datahihilo[idx2];
   yvlohilo[tdx] = datalohilo[idx2];
   yvhilolo[tdx] = datahilolo[idx2];
   yvlololo[tdx] = datalololo[idx2];

   __syncthreads();

   odg_inc(&zvhihihi[tdx],&zvlohihi[tdx],&zvhilohi[tdx],&zvlolohi[tdx],
           &zvhihilo[tdx],&zvlohilo[tdx],&zvhilolo[tdx],&zvlololo[tdx],
            yvhihihi[tdx], yvlohihi[tdx], yvhilohi[tdx], yvlolohi[tdx],
            yvhihilo[tdx], yvlohilo[tdx], yvhilolo[tdx], yvlololo[tdx]);

   __syncthreads();

   datahihihi[idx3] = zvhihihi[tdx]; // storing the output
   datalohihi[idx3] = zvlohihi[tdx];
   datahilohi[idx3] = zvhilohi[tdx];
   datalolohi[idx3] = zvlolohi[tdx];
   datahihilo[idx3] = zvhihilo[tdx];
   datalohilo[idx3] = zvlohilo[tdx];
   datahilolo[idx3] = zvhilolo[tdx];
   datalololo[idx3] = zvlololo[tdx];
}

__global__ void cmplx8_padded_convjobs
 ( double *datarehihihi, double *datarelohihi,
   double *datarehilohi, double *datarelolohi,
   double *datarehihilo, double *datarelohilo,
   double *datarehilolo, double *datarelololo,
   double *dataimhihihi, double *dataimlohihi,
   double *dataimhilohi, double *dataimlolohi,
   double *dataimhihilo, double *dataimlohilo,
   double *dataimhilolo, double *dataimlololo,
   int *in1idx, int *in2idx, int *outidx, int dim )
{
   const int bdx = blockIdx.x;           // index to the convolution job
   const int tdx = threadIdx.x;          // index to the output of the job
   const int idx1 = in1idx[bdx] + tdx;
   const int idx2 = in2idx[bdx] + tdx;
   const int idx3 = outidx[bdx] + tdx;

   __shared__ double xvrehihihi[odh_shmemsize];
   __shared__ double xvrelohihi[odh_shmemsize];
   __shared__ double xvrehilohi[odh_shmemsize];
   __shared__ double xvrelolohi[odh_shmemsize];
   __shared__ double xvrehihilo[odh_shmemsize];
   __shared__ double xvrelohilo[odh_shmemsize];
   __shared__ double xvrehilolo[odh_shmemsize];
   __shared__ double xvrelololo[odh_shmemsize];
   __shared__ double xvimhihihi[odh_shmemsize];
   __shared__ double xvimlohihi[odh_shmemsize];
   __shared__ double xvimhilohi[odh_shmemsize];
   __shared__ double xvimlolohi[odh_shmemsize];
   __shared__ double xvimhihilo[odh_shmemsize];
   __shared__ double xvimlohilo[odh_shmemsize];
   __shared__ double xvimhilolo[odh_shmemsize];
   __shared__ double xvimlololo[odh_shmemsize];
   __shared__ double yvrehihihi[2*odh_shmemsize];
   __shared__ double yvrelohihi[2*odh_shmemsize];
   __shared__ double yvrehilohi[2*odh_shmemsize];
   __shared__ double yvrelolohi[2*odh_shmemsize];
   __shared__ double yvrehihilo[2*odh_shmemsize];
   __shared__ double yvrelohilo[2*odh_shmemsize];
   __shared__ double yvrehilolo[2*odh_shmemsize];
   __shared__ double yvrelololo[2*odh_shmemsize];
   __shared__ double yvimhihihi[2*odh_shmemsize];
   __shared__ double yvimlohihi[2*odh_shmemsize];
   __shared__ double yvimhilohi[2*odh_shmemsize];
   __shared__ double yvimlolohi[2*odh_shmemsize];
   __shared__ double yvimhihilo[2*odh_shmemsize];
   __shared__ double yvimlohilo[2*odh_shmemsize];
   __shared__ double yvimhilolo[2*odh_shmemsize];
   __shared__ double yvimlololo[2*odh_shmemsize];
   __shared__ double zvrehihihi[odh_shmemsize];
   __shared__ double zvrelohihi[odh_shmemsize];
   __shared__ double zvrehilohi[odh_shmemsize];
   __shared__ double zvrelolohi[odh_shmemsize];
   __shared__ double zvrehihilo[odh_shmemsize];
   __shared__ double zvrelohilo[odh_shmemsize];
   __shared__ double zvrehilolo[odh_shmemsize];
   __shared__ double zvrelololo[odh_shmemsize];
   __shared__ double zvimhihihi[odh_shmemsize];
   __shared__ double zvimlohihi[odh_shmemsize];
   __shared__ double zvimhilohi[odh_shmemsize];
   __shared__ double zvimlolohi[odh_shmemsize];
   __shared__ double zvimhihilo[odh_shmemsize];
   __shared__ double zvimlohilo[odh_shmemsize];
   __shared__ double zvimhilolo[odh_shmemsize];
   __shared__ double zvimlololo[odh_shmemsize];

   double prodhihihi,prodlohihi,prodhilohi,prodlolohi;
   double prodhihilo,prodlohilo,prodhilolo,prodlololo;
   int ydx = dim + tdx;

   xvrehihihi[tdx] = datarehihihi[idx1];  // loading first input
   xvrelohihi[tdx] = datarelohihi[idx1]; 
   xvrehilohi[tdx] = datarehilohi[idx1]; 
   xvrelolohi[tdx] = datarelolohi[idx1]; 
   xvrehihilo[tdx] = datarehihilo[idx1];
   xvrelohilo[tdx] = datarelohilo[idx1]; 
   xvrehilolo[tdx] = datarehilolo[idx1]; 
   xvrelololo[tdx] = datarelololo[idx1]; 
   xvimhihihi[tdx] = dataimhihihi[idx1];
   xvimlohihi[tdx] = dataimlohihi[idx1];
   xvimhilohi[tdx] = dataimhilohi[idx1]; 
   xvimlolohi[tdx] = dataimlolohi[idx1]; 
   xvimhihilo[tdx] = dataimhihilo[idx1];
   xvimlohilo[tdx] = dataimlohilo[idx1];
   xvimhilolo[tdx] = dataimhilolo[idx1]; 
   xvimlololo[tdx] = dataimlololo[idx1]; 
   yvrehihihi[tdx] = 0.0;                 // padded with zeros
   yvrelohihi[tdx] = 0.0;
   yvrehilohi[tdx] = 0.0;
   yvrelolohi[tdx] = 0.0;
   yvrehihilo[tdx] = 0.0;
   yvrelohilo[tdx] = 0.0;
   yvrehilolo[tdx] = 0.0;
   yvrelololo[tdx] = 0.0;
   yvimhihihi[tdx] = 0.0;
   yvimlohihi[tdx] = 0.0;
   yvimhilohi[tdx] = 0.0;
   yvimlolohi[tdx] = 0.0;
   yvimhihilo[tdx] = 0.0;
   yvimlohilo[tdx] = 0.0;
   yvimhilolo[tdx] = 0.0;
   yvimlololo[tdx] = 0.0;
   yvrehihihi[ydx] = datarehihihi[idx2]; // loading second input
   yvrelohihi[ydx] = datarelohihi[idx2];
   yvrehilohi[ydx] = datarehilohi[idx2];
   yvrelolohi[ydx] = datarelolohi[idx2];
   yvrehihilo[ydx] = datarehihilo[idx2];
   yvrelohilo[ydx] = datarelohilo[idx2];
   yvrehilolo[ydx] = datarehilolo[idx2];
   yvrelololo[ydx] = datarelololo[idx2];
   yvimhihihi[ydx] = dataimhihihi[idx2];
   yvimlohihi[ydx] = dataimlohihi[idx2];
   yvimhilohi[ydx] = dataimhilohi[idx2];
   yvimlolohi[ydx] = dataimlolohi[idx2];
   yvimhihilo[ydx] = dataimhihilo[idx2];
   yvimlohilo[ydx] = dataimlohilo[idx2];
   yvimhilolo[ydx] = dataimhilolo[idx2];
   yvimlololo[ydx] = dataimlololo[idx2];

   __syncthreads();

   // zv[tdx] = xv[0]*yv[tdx];
   odg_mul(xvrehihihi[0],xvrelohihi[0],xvrehilohi[0],xvrelolohi[0],
           xvrehihilo[0],xvrelohilo[0],xvrehilolo[0],xvrelololo[0],
           yvrehihihi[ydx],yvrelohihi[ydx],yvrehilohi[ydx],yvrelolohi[ydx],
           yvrehihilo[ydx],yvrelohilo[ydx],yvrehilolo[ydx],yvrelololo[ydx],
           &zvrehihihi[tdx],&zvrelohihi[tdx],
           &zvrehilohi[tdx],&zvrelolohi[tdx],
           &zvrehihilo[tdx],&zvrelohilo[tdx],
           &zvrehilolo[tdx],&zvrelololo[tdx]);
   __syncthreads();
   odg_mul(xvimhihihi[0],xvimlohihi[0],xvimhilohi[0],xvimlolohi[0],
           xvimhihilo[0],xvimlohilo[0],xvimhilolo[0],xvimlololo[0],
           yvimhihihi[ydx],yvimlohihi[ydx],yvimhilohi[ydx],yvimlolohi[ydx],
           yvimhihilo[ydx],yvimlohilo[ydx],yvimhilolo[ydx],yvimlololo[ydx],
           &prodhihihi,&prodlohihi,&prodhilohi,&prodlolohi,
           &prodhihilo,&prodlohilo,&prodhilolo,&prodlololo);
   __syncthreads();
   odg_minus(&prodhihihi,&prodlohihi,&prodhilohi,&prodlolohi,
             &prodhihilo,&prodlohilo,&prodhilolo,&prodlololo);
   odg_inc(&zvrehihihi[tdx],&zvrelohihi[tdx],&zvrehilohi[tdx],&zvrelolohi[tdx],
           &zvrehihilo[tdx],&zvrelohilo[tdx],&zvrehilolo[tdx],&zvrelololo[tdx],
           prodhihihi,prodlohihi,prodhilohi,prodlolohi,
           prodhihilo,prodlohilo,prodhilolo,prodlololo);
   __syncthreads();

   odg_mul(xvrehihihi[0],xvrelohihi[0],xvrehilohi[0],xvrelolohi[0],
           xvrehihilo[0],xvrelohilo[0],xvrehilolo[0],xvrelololo[0],
           yvimhihihi[ydx],yvimlohihi[ydx],yvimhilohi[ydx],yvimlolohi[ydx],
           yvimhihilo[ydx],yvimlohilo[ydx],yvimhilolo[ydx],yvimlololo[ydx],
           &zvimhihihi[tdx],&zvimlohihi[tdx],&zvimhilohi[tdx],&zvimlolohi[tdx],
           &zvimhihilo[tdx],&zvimlohilo[tdx],&zvimhilolo[tdx],&zvimlololo[tdx]);
   __syncthreads();
   odg_mul(xvimhihihi[0],xvimlohihi[0],xvimhilohi[0],xvimlolohi[0],
           xvimhihilo[0],xvimlohilo[0],xvimhilolo[0],xvimlololo[0],
           yvrehihihi[ydx],yvrelohihi[ydx],yvrehilohi[ydx],yvrelolohi[ydx],
           yvrehihilo[ydx],yvrelohilo[ydx],yvrehilolo[ydx],yvrelololo[ydx],
           &prodhihihi,&prodlohihi,&prodhilohi,&prodlolohi,
           &prodhihilo,&prodlohilo,&prodhilolo,&prodlololo);
   __syncthreads();
   odg_inc(&zvimhihihi[tdx],&zvimlohihi[tdx],&zvimhilohi[tdx],&zvimlolohi[tdx],
           &zvimhihilo[tdx],&zvimlohilo[tdx],&zvimhilolo[tdx],&zvimlololo[tdx],
           prodhihihi,prodlohihi,prodhilohi,prodlolohi,
           prodhihilo,prodlohilo,prodhilolo,prodlololo);
   __syncthreads();

   for(int i=1; i<dim; i++) // zv[tdx] = zv[tdx] + xv[i]*yv[dim+tdx-i];
   {
      ydx = dim + tdx - i;

      odg_mul(xvrehihihi[i],xvrelohihi[i],xvrehilohi[i],xvrelolohi[i],
              xvrehihilo[i],xvrelohilo[i],xvrehilolo[i],xvrelololo[i],
              yvrehihihi[ydx],yvrelohihi[ydx],yvrehilohi[ydx],yvrelolohi[ydx],
              yvrehihilo[ydx],yvrelohilo[ydx],yvrehilolo[ydx],yvrelololo[ydx],
              &prodhihihi,&prodlohihi,&prodhilohi,&prodlolohi,
              &prodhihilo,&prodlohilo,&prodhilolo,&prodlololo);
      __syncthreads();
      odg_inc(&zvrehihihi[tdx],&zvrelohihi[tdx],
              &zvrehilohi[tdx],&zvrelolohi[tdx],
              &zvrehihilo[tdx],&zvrelohilo[tdx],
              &zvrehilolo[tdx],&zvrelololo[tdx],
              prodhihihi,prodlohihi,prodhilohi,prodlolohi,
              prodhihilo,prodlohilo,prodhilolo,prodlololo);
      __syncthreads();
      odg_mul(xvimhihihi[i],xvimlohihi[i],xvimhilohi[i],xvimlolohi[i],
              xvimhihilo[i],xvimlohilo[i],xvimhilolo[i],xvimlololo[i],
              yvimhihihi[ydx],yvimlohihi[ydx],yvimhilohi[ydx],yvimlolohi[ydx],
              yvimhihilo[ydx],yvimlohilo[ydx],yvimhilolo[ydx],yvimlololo[ydx],
              &prodhihihi,&prodlohihi,&prodhilohi,&prodlolohi,
              &prodhihilo,&prodlohilo,&prodhilolo,&prodlololo);
      __syncthreads();
      odg_minus(&prodhihihi,&prodlohihi,&prodhilohi,&prodlolohi,
                &prodhihilo,&prodlohilo,&prodhilolo,&prodlololo);
      odg_inc(&zvrehihihi[tdx],&zvrelohihi[tdx],
              &zvrehilohi[tdx],&zvrelolohi[tdx],
              &zvrehihilo[tdx],&zvrelohilo[tdx],
              &zvrehilolo[tdx],&zvrelololo[tdx],
              prodhihihi,prodlohihi,prodhilohi,prodlolohi,
              prodhihilo,prodlohilo,prodhilolo,prodlololo);
      __syncthreads();

      odg_mul(xvrehihihi[i],xvrelohihi[i],xvrehilohi[i],xvrelolohi[i],
              xvrehihilo[i],xvrelohilo[i],xvrehilolo[i],xvrelololo[i],
              yvimhihihi[ydx],yvimlohihi[ydx],yvimhilohi[ydx],yvimlolohi[ydx],
              yvimhihilo[ydx],yvimlohilo[ydx],yvimhilolo[ydx],yvimlololo[ydx],
              &prodhihihi,&prodlohihi,&prodhilohi,&prodlolohi,
              &prodhihilo,&prodlohilo,&prodhilolo,&prodlololo);
      __syncthreads();
      odg_inc(&zvimhihihi[tdx],&zvimlohihi[tdx],
              &zvimhilohi[tdx],&zvimlolohi[tdx],
              &zvimhihilo[tdx],&zvimlohilo[tdx],
              &zvimhilolo[tdx],&zvimlololo[tdx],
              prodhihihi,prodlohihi,prodhilohi,prodlolohi,
              prodhihilo,prodlohilo,prodhilolo,prodlololo);
      __syncthreads();
      odg_mul(xvimhihihi[i],xvimlohihi[i],xvimhilohi[i],xvimlolohi[i],
              xvimhihilo[i],xvimlohilo[i],xvimhilolo[i],xvimlololo[i],
              yvrehihihi[ydx],yvrelohihi[ydx],yvrehilohi[ydx],yvrelolohi[ydx],
              yvrehihilo[ydx],yvrelohilo[ydx],yvrehilolo[ydx],yvrelololo[ydx],
              &prodhihihi,&prodlohihi,&prodhilohi,&prodlolohi,
              &prodhihilo,&prodlohilo,&prodhilolo,&prodlololo);
      __syncthreads();
      odg_inc(&zvimhihihi[tdx],&zvimlohihi[tdx],
              &zvimhilohi[tdx],&zvimlolohi[tdx],
              &zvimhihilo[tdx],&zvimlohilo[tdx],
              &zvimhilolo[tdx],&zvimlololo[tdx],
              prodhihihi,prodlohihi,prodhilohi,prodlolohi,
              prodhihilo,prodlohilo,prodhilolo,prodlololo);
      __syncthreads();
   }
   __syncthreads();

   datarehihihi[idx3] = zvrehihihi[tdx]; // storing the output
   datarelohihi[idx3] = zvrelohihi[tdx];
   datarehilohi[idx3] = zvrehilohi[tdx];
   datarelolohi[idx3] = zvrelolohi[tdx];
   datarehihilo[idx3] = zvrehihilo[tdx];
   datarelohilo[idx3] = zvrelohilo[tdx];
   datarehilolo[idx3] = zvrehilolo[tdx];
   datarelololo[idx3] = zvrelololo[tdx];
   dataimhihihi[idx3] = zvimhihihi[tdx]; 
   dataimlohihi[idx3] = zvimlohihi[tdx]; 
   dataimhilohi[idx3] = zvimhilohi[tdx];
   dataimlolohi[idx3] = zvimlolohi[tdx];
   dataimhihilo[idx3] = zvimhihilo[tdx]; 
   dataimlohilo[idx3] = zvimlohilo[tdx]; 
   dataimhilolo[idx3] = zvimhilolo[tdx];
   dataimlololo[idx3] = zvimlololo[tdx];
}

__global__ void cmplx8vectorized_flipsigns
 ( double *datarihihihi, double *datarilohihi,
   double *datarihilohi, double *datarilolohi,
   double *datarihihilo, double *datarilohilo,
   double *datarihilolo, double *datarilololo, int *flpidx, int dim )
{
   const int bdx = blockIdx.x;    // index to the series to flip
   const int tdx = threadIdx.x;
   const int idx = flpidx[bdx] + tdx; // which number to flip

   double x; // register to load data from global memory

   x = datarihihihi[idx]; x = -x; datarihihihi[idx] = x;
   x = datarilohihi[idx]; x = -x; datarilohihi[idx] = x;
   x = datarihilohi[idx]; x = -x; datarihilohi[idx] = x;
   x = datarilolohi[idx]; x = -x; datarilolohi[idx] = x;
   x = datarihihilo[idx]; x = -x; datarihihilo[idx] = x;
   x = datarilohilo[idx]; x = -x; datarilohilo[idx] = x;
   x = datarihilolo[idx]; x = -x; datarihilolo[idx] = x;
   x = datarilololo[idx]; x = -x; datarilololo[idx] = x;
}

void GPU_cmplx8vectorized_flipsigns
 ( int deg, int nbrflips, int *flipidx,
   double *datarihihihi, double *datarilohihi,
   double *datarihilohi, double *datarilolohi,
   double *datarihihilo, double *datarilohilo,
   double *datarihilolo, double *datarilololo,
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
/*
   for(int i=0; i<nbrflips; i++)
   {
      int in1idx = flipidx[i];

      cout << "Flip of job " << i
           << " at index " << in1idx << " :" << endl;

      double rihihihi,rilohihi,rihilohi,rilolohi;
      double rihihilo,rilohilo,rihilolo,rilololo;

      cudaMemcpy(&rihihihi,&datarihihihi[in1idx],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&rilohihi,&datarilohihi[in1idx],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&rihilohi,&datarihilohi[in1idx],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&rilolohi,&datarilolohi[in1idx],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&rihihilo,&datarihihilo[in1idx],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&rilohilo,&datarilohilo[in1idx],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&rihilolo,&datarihilolo[in1idx],sizeof(double),
                 cudaMemcpyDeviceToHost);
      cudaMemcpy(&rilololo,&datarilololo[in1idx],sizeof(double),
                 cudaMemcpyDeviceToHost);

      cout << scientific << setprecision(16);
      cout << rihihihi << "  " << rilohihi << endl;
      cout << rihilohi << "  " << rilolohi << endl;
      cout << rihihilo << "  " << rilohilo << endl;
      cout << rihilolo << "  " << rilololo << endl;
   }
 */
   cudaEventRecord(start);
   cmplx8vectorized_flipsigns<<<nbrflips,deg1>>>
      (datarihihihi,datarilohihi,datarihilohi,datarilolohi,
       datarihihilo,datarilohilo,datarihilolo,datarilololo,flipidx_d,deg1);
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

__global__ void dbl8_update_addjobs
 ( double *datahihihi, double *datalohihi,
   double *datahilohi, double *datalolohi,
   double *datahihilo, double *datalohilo,
   double *datahilolo, double *datalololo,
   int *in1idx, int *in2idx, int *outidx, int dim )
{
   const int bdx = blockIdx.x;           // index to the convolution job
   const int tdx = threadIdx.x;          // index to the output of the job
   const int idx1 = in1idx[bdx] + tdx;
   const int idx2 = in2idx[bdx] + tdx;
   const int idx3 = outidx[bdx] + tdx;

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

   xvhihihi[tdx] = datahihihi[idx1];  // loading first input
   xvlohihi[tdx] = datalohihi[idx1];
   xvhilohi[tdx] = datahilohi[idx1];
   xvlolohi[tdx] = datalolohi[idx1];
   xvhihilo[tdx] = datahihilo[idx1];
   xvlohilo[tdx] = datalohilo[idx1];
   xvhilolo[tdx] = datahilolo[idx1];
   xvlololo[tdx] = datalololo[idx1];
   yvhihihi[tdx] = datahihihi[idx2];  // loading second input
   yvlohihi[tdx] = datalohihi[idx2];
   yvhilohi[tdx] = datahilohi[idx2];
   yvlolohi[tdx] = datalolohi[idx2];
   yvhihilo[tdx] = datahihilo[idx2];
   yvlohilo[tdx] = datalohilo[idx2];
   yvhilolo[tdx] = datahilolo[idx2];
   yvlololo[tdx] = datalololo[idx2];

   // zv[tdx] = xv[tdx] + yv[tdx];

   __syncthreads();

   odg_add( xvhihihi[tdx], xvlohihi[tdx], xvhilohi[tdx], xvlolohi[tdx],
            xvhihilo[tdx], xvlohilo[tdx], xvhilolo[tdx], xvlololo[tdx],
            yvhihihi[tdx], yvlohihi[tdx], yvhilohi[tdx], yvlolohi[tdx],
            yvhihilo[tdx], yvlohilo[tdx], yvhilolo[tdx], yvlololo[tdx],
           &zvhihihi[tdx],&zvlohihi[tdx],&zvhilohi[tdx],&zvlolohi[tdx],
           &zvhihilo[tdx],&zvlohilo[tdx],&zvhilolo[tdx],&zvlololo[tdx]);

   __syncthreads();

   datahihihi[idx3] = zvhihihi[tdx]; // storing the output
   datalohihi[idx3] = zvlohihi[tdx];
   datahilohi[idx3] = zvhilohi[tdx];
   datalolohi[idx3] = zvlolohi[tdx];
   datahihilo[idx3] = zvhihilo[tdx];
   datalohilo[idx3] = zvlohilo[tdx];
   datahilolo[idx3] = zvhilolo[tdx];
   datalololo[idx3] = zvlololo[tdx];
}

void dbl_convoluted_data8_to_output
 ( double *datahihihi, double *datalohihi,
   double *datahilohi, double *datalolohi,
   double *datahihilo, double *datalohilo,
   double *datahilolo, double *datalololo,
   double **outputhihihi, double **outputlohihi,
   double **outputhilohi, double **outputlolohi,
   double **outputhihilo, double **outputlohilo,
   double **outputhilolo, double **outputlololo,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart, int vrblvl )
{
   const int deg1 = deg+1;
   int ix0,ix1,ix2;

   for(int i=0; i<=deg; i++) // output[dim][i] = data[i];
   {
      outputhihihi[dim][i] = datahihihi[i];
      outputlohihi[dim][i] = datalohihi[i];
      outputhilohi[dim][i] = datahilohi[i];
      outputlolohi[dim][i] = datalolohi[i];
      outputhihilo[dim][i] = datahihilo[i];
      outputlohilo[dim][i] = datalohilo[i];
      outputhilolo[dim][i] = datahilolo[i];
      outputlololo[dim][i] = datalololo[i];
   }
   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++) // output[i][j] = 0.0;
      {
         outputhihihi[i][j] = 0.0;
         outputlohihi[i][j] = 0.0;
         outputhilohi[i][j] = 0.0;
         outputlolohi[i][j] = 0.0;
         outputhihilo[i][j] = 0.0;
         outputlohilo[i][j] = 0.0;
         outputhilolo[i][j] = 0.0;
         outputlololo[i][j] = 0.0;
      }

   for(int k=0; k<nbr; k++)
   {
      ix1 = fstart[k] + (nvr[k]-1)*deg1;
      
      if(vrblvl > 1)
         cout << "monomial " << k << " update starts at " << ix1 << endl;

      for(int i=0; i<=deg; i++) // output[dim][i] += data[ix1++];
         odf_inc(&outputhihihi[dim][i],&outputlohihi[dim][i],
                 &outputhilohi[dim][i],&outputlolohi[dim][i],
                 &outputhihilo[dim][i],&outputlohilo[dim][i],
                 &outputhilolo[dim][i],&outputlololo[dim][i],
                    datahihihi[ix1],      datalohihi[ix1],
                    datahilohi[ix1],      datalolohi[ix1],
                    datahihilo[ix1],      datalohilo[ix1],
                    datahilolo[ix1],      datalololo[ix1++]);
     
      ix0 = idx[k][0];
      if(nvr[k] == 1)
      {
         ix1 = (1 + k)*deg1;
            
         for(int i=0; i<=deg; i++) // output[ix0][i] += data[ix1++];
            odf_inc(&outputhihihi[ix0][i],&outputlohihi[ix0][i],
                    &outputhilohi[ix0][i],&outputlolohi[ix0][i],
                    &outputhihilo[ix0][i],&outputlohilo[ix0][i],
                    &outputhilolo[ix0][i],&outputlololo[ix0][i],
                       datahihihi[ix1],      datalohihi[ix1],
                       datahilohi[ix1],      datalolohi[ix1],
                       datahihilo[ix1],      datalohilo[ix1],
                       datahilolo[ix1],      datalololo[ix1++]);
      }
      else
      {                               // update first and last derivative
         ix2 = nvr[k]-3;
         if(ix2 < 0) ix2 = 0;
         ix1 = bstart[k] + ix2*deg1;

         for(int i=0; i<=deg; i++) // output[ix0][i] += data[ix1++];
            odf_inc(&outputhihihi[ix0][i],&outputlohihi[ix0][i],
                    &outputhilohi[ix0][i],&outputlolohi[ix0][i],
                    &outputhihilo[ix0][i],&outputlohilo[ix0][i],
                    &outputhilolo[ix0][i],&outputlololo[ix0][i],
                       datahihihi[ix1],      datalohihi[ix1],
                       datahilohi[ix1],      datalolohi[ix1],
                       datahihilo[ix1],      datalohilo[ix1],
                       datahilolo[ix1],      datalololo[ix1++]);

         ix2 = nvr[k]-2;
         ix1 = fstart[k] + ix2*deg1;
         ix0 = idx[k][ix2+1];

         for(int i=0; i<=deg; i++) // output[ix0][i] += data[ix1++];
            odf_inc(&outputhihihi[ix0][i],&outputlohihi[ix0][i],
                    &outputhilohi[ix0][i],&outputlolohi[ix0][i],
                    &outputhihilo[ix0][i],&outputlohilo[ix0][i],
                    &outputhilolo[ix0][i],&outputlololo[ix0][i],
                       datahihihi[ix1],      datalohihi[ix1],
                       datahilolo[ix1],      datalolohi[ix1],
                       datahihilo[ix1],      datalohilo[ix1],
                       datahilolo[ix1],      datalololo[ix1++]);
 
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
                  odf_inc(&outputhihihi[ix0][i],&outputlohihi[ix0][i],
                          &outputhilohi[ix0][i],&outputlolohi[ix0][i],
                          &outputhihilo[ix0][i],&outputlohilo[ix0][i],
                          &outputhilolo[ix0][i],&outputlololo[ix0][i],
                             datahihihi[ix1],      datalohihi[ix1],
                             datahilohi[ix1],      datalolohi[ix1],
                             datahihilo[ix1],      datalohilo[ix1],
                             datahilolo[ix1],      datalololo[ix1++]);
            }
         }
      }
   }
}

void dbl_added_data8_to_output
 ( double *datahihihi, double *datalohihi,
   double *datahilohi, double *datalolohi,
   double *datahihilo, double *datalohilo,
   double *datahilolo, double *datalololo,
   double **outputhihihi, double **outputlohihi,
   double **outputhilohi, double **outputlolohi,
   double **outputhihilo, double **outputlohilo,
   double **outputhilolo, double **outputlololo,
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
      outputhihihi[dim][i] = datahihihi[ix];
      outputlohihi[dim][i] = datalohihi[ix];
      outputhilohi[dim][i] = datahilohi[ix];
      outputlolohi[dim][i] = datalolohi[ix];
      outputhihilo[dim][i] = datahihilo[ix];
      outputlohilo[dim][i] = datalohilo[ix];
      outputhilolo[dim][i] = datahilolo[ix];
      outputlololo[dim][i] = datalololo[ix++];
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
            outputhihihi[0][i] = 0.0;
            outputlohihi[0][i] = 0.0;
            outputhilohi[0][i] = 0.0;
            outputlolohi[0][i] = 0.0;
            outputhihilo[0][i] = 0.0;
            outputlohilo[0][i] = 0.0;
            outputhilolo[0][i] = 0.0;
            outputlololo[0][i] = 0.0;
         }
      }
      else
      {
         int cffidx = (1 + difidx)*deg1;

         if(vrblvl > 1)
            cout << "updating derivative with coefficient ..." << endl;

         for(int i=0; i<=deg; i++)
         {
            outputhihihi[0][i] = datahihihi[cffidx];
            outputlohihi[0][i] = datalohihi[cffidx];
            outputhilohi[0][i] = datahilohi[cffidx];
            outputlolohi[0][i] = datalolohi[cffidx];
            outputhihilo[0][i] = datahihilo[cffidx];
            outputlohilo[0][i] = datalohilo[cffidx];
            outputhilolo[0][i] = datahilolo[cffidx];
            outputlololo[0][i] = datalololo[cffidx++];
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
         outputhihihi[0][i] = datahihihi[ix];
         outputlohihi[0][i] = datalohihi[ix];
         outputhilohi[0][i] = datahilohi[ix];
         outputlolohi[0][i] = datalolohi[ix];
         outputhihilo[0][i] = datahihilo[ix];
         outputlohilo[0][i] = datalohilo[ix];
         outputhilolo[0][i] = datahilolo[ix];
         outputlololo[0][i] = datalololo[ix++];
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
               outputhihihi[k][i] = 0.0;
               outputlohihi[k][i] = 0.0;
               outputhilohi[k][i] = 0.0;
               outputlolohi[k][i] = 0.0;
               outputhihilo[k][i] = 0.0;
               outputlohilo[k][i] = 0.0;
               outputhilolo[k][i] = 0.0;
               outputlololo[k][i] = 0.0;
            }
         }
         else
         {
            int cffidx = (1 + difidx)*deg1;

            if(vrblvl > 1)
               cout << "updating derivative with coefficient ..." << endl;

            for(int i=0; i<=deg; i++)
            {
               outputhihihi[k][i] = datahihihi[cffidx];
               outputlohihi[k][i] = datalohihi[cffidx];
               outputhilohi[k][i] = datahilohi[cffidx];
               outputlolohi[k][i] = datalolohi[cffidx];
               outputhihilo[k][i] = datahihilo[cffidx];
               outputlohilo[k][i] = datalohilo[cffidx];
               outputhilolo[k][i] = datahilolo[cffidx];
               outputlololo[k][i] = datalololo[cffidx++];
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
               outputhihihi[k][i] = datahihihi[ix];
               outputlohihi[k][i] = datalohihi[ix];
               outputhilohi[k][i] = datahilohi[ix];
               outputlolohi[k][i] = datalolohi[ix];
               outputhihilo[k][i] = datahihilo[ix];
               outputlohilo[k][i] = datalohilo[ix];
               outputhilolo[k][i] = datahilolo[ix];
               outputlololo[k][i] = datalololo[ix++];
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
               outputhihihi[k][i] = datahihihi[ix];
               outputlohihi[k][i] = datalohihi[ix];
               outputhilohi[k][i] = datahilohi[ix];
               outputlolohi[k][i] = datalolohi[ix];
               outputhihilo[k][i] = datahihilo[ix];
               outputlohilo[k][i] = datalohilo[ix];
               outputhilolo[k][i] = datahilolo[ix];
               outputlololo[k][i] = datalololo[ix++];
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
               outputhihihi[k][i] = datahihihi[ix];
               outputlohihi[k][i] = datalohihi[ix];
               outputhilohi[k][i] = datahilohi[ix];
               outputlolohi[k][i] = datalolohi[ix];
               outputhihilo[k][i] = datahihilo[ix];
               outputlohilo[k][i] = datalohilo[ix];
               outputhilolo[k][i] = datahilolo[ix];
               outputlololo[k][i] = datalololo[ix++];
            }
         }
      }
   }
}

void cmplx_added_data8vectorized_to_output
 ( double *datarihihihi, double *datarilohihi,
   double *datarihilohi, double *datarilolohi,
   double *datarihihilo, double *datarilohilo,
   double *datarihilolo, double *datarilololo,
   double **outputrehihihi, double **outputrelohihi,
   double **outputrehilohi, double **outputrelolohi,
   double **outputrehihilo, double **outputrelohilo,
   double **outputrehilolo, double **outputrelololo,
   double **outputimhihihi, double **outputimlohihi,
   double **outputimhilohi, double **outputimlolohi,
   double **outputimhihilo, double **outputimlohilo,
   double **outputimhilolo, double **outputimlololo,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart,
   int totcff, int offsetri, ComplexAdditionJobs jobs, int vrblvl )
{
   const int deg1 = deg + 1;
   const int lastmon = nbr-1;
   const int lastidx = nvr[lastmon]-1;
   const int totcffoffset = totcff + offsetri;
   int ix1re,ix2im;

   ix1re = fstart[lastmon] + lastidx*deg1;
   ix2im = fstart[lastmon] + lastidx*deg1 + totcffoffset;

   if(vrblvl > 1)
      cout << "Updating value starting at "
           << ix1re << ", " << ix2im  << " in data." << endl;

   for(int i=0; i<=deg; i++) // output[dim][i] = data[ix++];
   {
      outputrehihihi[dim][i] = datarihihihi[ix1re];
      outputrelohihi[dim][i] = datarilohihi[ix1re];
      outputrehilohi[dim][i] = datarihilohi[ix1re];
      outputrelolohi[dim][i] = datarilolohi[ix1re];
      outputrehihilo[dim][i] = datarihihilo[ix1re];
      outputrelohilo[dim][i] = datarilohilo[ix1re];
      outputrehilolo[dim][i] = datarihilolo[ix1re];
      outputrelololo[dim][i] = datarilololo[ix1re++];
      outputimhihihi[dim][i] = datarihihihi[ix2im];
      outputimlohihi[dim][i] = datarilohihi[ix2im];
      outputimhilohi[dim][i] = datarihilohi[ix2im];
      outputimlolohi[dim][i] = datarilolohi[ix2im];
      outputimhihilo[dim][i] = datarihihilo[ix2im];
      outputimlohilo[dim][i] = datarilohilo[ix2im];
      outputimhilolo[dim][i] = datarihilolo[ix2im];
      outputimlololo[dim][i] = datarilololo[ix2im++];
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
            outputrehihihi[0][i] = 0.0; outputrelohihi[0][i] = 0.0;
            outputrehilohi[0][i] = 0.0; outputrelolohi[0][i] = 0.0;
            outputrehihilo[0][i] = 0.0; outputrelohilo[0][i] = 0.0;
            outputrehilolo[0][i] = 0.0; outputrelololo[0][i] = 0.0;
            outputimhihihi[0][i] = 0.0; outputimlohihi[0][i] = 0.0; 
            outputimhilohi[0][i] = 0.0; outputimlolohi[0][i] = 0.0;
            outputimhihilo[0][i] = 0.0; outputimlohilo[0][i] = 0.0; 
            outputimhilolo[0][i] = 0.0; outputimlololo[0][i] = 0.0;
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
            outputrehihihi[0][i] = datarihihihi[ix1re];
            outputrelohihi[0][i] = datarilohihi[ix1re];
            outputrehilohi[0][i] = datarihilohi[ix1re];
            outputrelolohi[0][i] = datarilolohi[ix1re];
            outputrehihilo[0][i] = datarihihilo[ix1re];
            outputrelohilo[0][i] = datarilohilo[ix1re];
            outputrehilolo[0][i] = datarihilolo[ix1re];
            outputrelololo[0][i] = datarilololo[ix1re++];
            outputimhihihi[0][i] = datarihihihi[ix2im];
            outputimlohihi[0][i] = datarilohihi[ix2im];
            outputimhilohi[0][i] = datarihilohi[ix2im];
            outputimlolohi[0][i] = datarilolohi[ix2im];
            outputimhihilo[0][i] = datarihihilo[ix2im];
            outputimlohilo[0][i] = datarilohilo[ix2im];
            outputimhilolo[0][i] = datarihilolo[ix2im];
            outputimlololo[0][i] = datarilololo[ix2im++];
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
         cout << "Updating derivative 0 at "
              << ix1re << ", " << ix2im << " in data." << endl;

      for(int i=0; i<=deg; i++) // output[0][i] = data[ix++];
      {
         outputrehihihi[0][i] = datarihihihi[ix1re];
         outputrelohihi[0][i] = datarilohihi[ix1re];
         outputrehilohi[0][i] = datarihilohi[ix1re];
         outputrelolohi[0][i] = datarilolohi[ix1re];
         outputrehihilo[0][i] = datarihihilo[ix1re];
         outputrelohilo[0][i] = datarilohilo[ix1re];
         outputrehilolo[0][i] = datarihilolo[ix1re];
         outputrelololo[0][i] = datarilololo[ix1re++];
         outputimhihihi[0][i] = datarihihihi[ix2im];
         outputimlohihi[0][i] = datarilohihi[ix2im];
         outputimhilohi[0][i] = datarihilohi[ix2im];
         outputimlolohi[0][i] = datarilolohi[ix2im];
         outputimhihilo[0][i] = datarihihilo[ix2im];
         outputimlohilo[0][i] = datarilohilo[ix2im];
         outputimhilolo[0][i] = datarihilolo[ix2im];
         outputimlololo[0][i] = datarilololo[ix2im++];
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
               outputrehihihi[k][i] = 0.0; outputrelohihi[k][i] = 0.0;
               outputrehilohi[k][i] = 0.0; outputrelolohi[k][i] = 0.0;
               outputrehihilo[k][i] = 0.0; outputrelohilo[k][i] = 0.0;
               outputrehilolo[k][i] = 0.0; outputrelololo[k][i] = 0.0;
               outputimhihihi[k][i] = 0.0; outputimlohihi[k][i] = 0.0;
               outputimhilohi[k][i] = 0.0; outputimlolohi[k][i] = 0.0;
               outputimhihilo[k][i] = 0.0; outputimlohilo[k][i] = 0.0;
               outputimhilolo[k][i] = 0.0; outputimlololo[k][i] = 0.0;
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
               outputrehihihi[k][i] = datarihihihi[ix1re];
               outputrelohihi[k][i] = datarilohihi[ix1re];
               outputrehilohi[k][i] = datarihilohi[ix1re];
               outputrelolohi[k][i] = datarilolohi[ix1re];
               outputrehihilo[k][i] = datarihihilo[ix1re];
               outputrelohilo[k][i] = datarilohilo[ix1re];
               outputrehilolo[k][i] = datarihilolo[ix1re];
               outputrelololo[k][i] = datarilololo[ix1re++];
               outputimhihihi[k][i] = datarihihihi[ix2im];
               outputimlohihi[k][i] = datarilohihi[ix2im];
               outputimhilohi[k][i] = datarihilohi[ix2im];
               outputimlolohi[k][i] = datarilolohi[ix2im];
               outputimhihilo[k][i] = datarihihilo[ix2im];
               outputimlohilo[k][i] = datarilohilo[ix2im];
               outputimhilolo[k][i] = datarihilolo[ix2im];
               outputimlololo[k][i] = datarilololo[ix2im++];
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
               cout << "Updating derivative " << k << " at "
                    << ix1re << ", " << ix2im << " in data." << endl;

            for(int i=0; i<=deg; i++) // output[k][i] = data[ix++];
            {
               outputrehihihi[k][i] = datarihihihi[ix1re];
               outputrelohihi[k][i] = datarilohihi[ix1re];
               outputrehilohi[k][i] = datarihilohi[ix1re];
               outputrelolohi[k][i] = datarilolohi[ix1re];
               outputrehihilo[k][i] = datarihihilo[ix1re];
               outputrelohilo[k][i] = datarilohilo[ix1re];
               outputrehilolo[k][i] = datarihilolo[ix1re];
               outputrelololo[k][i] = datarilololo[ix1re++];
               outputimhihihi[k][i] = datarihihihi[ix2im];
               outputimlohihi[k][i] = datarilohihi[ix2im];
               outputimhilohi[k][i] = datarihilohi[ix2im];
               outputimlolohi[k][i] = datarilolohi[ix2im];
               outputimhihilo[k][i] = datarihihilo[ix2im];
               outputimlohilo[k][i] = datarilohilo[ix2im];
               outputimhilolo[k][i] = datarihilolo[ix2im];
               outputimlololo[k][i] = datarilololo[ix2im++];
            }
         }
         else if(idx[ix0][nvr[ix0]-1] == k) // k is last variable
         {
            int ix2 = nvr[ix0]-2;

            ix1re = fstart[ix0] + ix2*deg1;
            ix2im = fstart[ix0] + ix2*deg1 + totcffoffset;
 
            if(vrblvl > 1)
               cout << "Updating derivative " << k << " at "
                    << ix1re << ", " << ix2im << " in data." << endl;

            for(int i=0; i<=deg; i++) // output[k][i] = data[ix++];
            {
               outputrehihihi[k][i] = datarihihihi[ix1re];
               outputrelohihi[k][i] = datarilohihi[ix1re];
               outputrehilohi[k][i] = datarihilohi[ix1re];
               outputrelolohi[k][i] = datarilolohi[ix1re];
               outputrehihilo[k][i] = datarihihilo[ix1re];
               outputrelohilo[k][i] = datarilohilo[ix1re];
               outputrehilolo[k][i] = datarihilolo[ix1re];
               outputrelololo[k][i] = datarilololo[ix1re++];
               outputimhihihi[k][i] = datarihihihi[ix2im];
               outputimlohihi[k][i] = datarilohihi[ix2im];
               outputimhilohi[k][i] = datarihilohi[ix2im];
               outputimlolohi[k][i] = datarilolohi[ix2im];
               outputimhihilo[k][i] = datarihihilo[ix2im];
               outputimlohilo[k][i] = datarilohilo[ix2im];
               outputimhilolo[k][i] = datarihilolo[ix2im];
               outputimlololo[k][i] = datarilololo[ix2im++];
            }
         }
         else // derivative is in some cross product
         {
            int ix2 = jobs.position(nvr[ix0],idx[ix0],k) - 1;

            ix1re = cstart[ix0] + ix2*deg1;
            ix2im = cstart[ix0] + ix2*deg1 + totcffoffset;

            if(vrblvl > 1)
               cout << "Updating derivative " << k << " at "
                    << ix1re << ", " << ix2im << " in data." << endl;

            for(int i=0; i<=deg; i++) // output[k][i] = data[ix++];
            {
               outputrehihihi[k][i] = datarihihihi[ix1re];
               outputrelohihi[k][i] = datarilohihi[ix1re];
               outputrehilohi[k][i] = datarihilohi[ix1re];
               outputrelolohi[k][i] = datarilolohi[ix1re];
               outputrehihilo[k][i] = datarihihilo[ix1re];
               outputrelohilo[k][i] = datarilohilo[ix1re];
               outputrehilolo[k][i] = datarihilolo[ix1re];
               outputrelololo[k][i] = datarilololo[ix1re++];
               outputimhihihi[k][i] = datarihihihi[ix2im];
               outputimlohihi[k][i] = datarilohihi[ix2im];
               outputimhilohi[k][i] = datarihilohi[ix2im];
               outputimlolohi[k][i] = datarilolohi[ix2im];
               outputimhihilo[k][i] = datarihihilo[ix2im];
               outputimlohilo[k][i] = datarilohilo[ix2im];
               outputimhilolo[k][i] = datarihilolo[ix2im];
               outputimlololo[k][i] = datarilololo[ix2im++];
            }
         }
      }
   }
}

void dbl8_data_setup
 ( int dim, int nbr, int deg, int totcff,
   double *datahihihi, double *datalohihi,
   double *datahilohi, double *datalolohi,
   double *datahihilo, double *datalohilo,
   double *datahilolo, double *datalololo,
   double *csthihihi, double *cstlohihi,
   double *csthilohi, double *cstlolohi,
   double *csthihilo, double *cstlohilo,
   double *csthilolo, double *cstlololo,
   double **cffhihihi, double **cfflohihi,
   double **cffhilohi, double **cfflolohi,
   double **cffhihilo, double **cfflohilo,
   double **cffhilolo, double **cfflololo,
   double **inputhihihi, double **inputlohihi,
   double **inputhilohi, double **inputlolohi,
   double **inputhihilo, double **inputlohilo,
   double **inputhilolo, double **inputlololo )
{
   const int deg1 = deg+1;
   int ix = 0;

   for(int i=0; i<deg1; i++)
   {
      datahihihi[ix]   = csthihihi[i];
      datalohihi[ix]   = cstlohihi[i];
      datahilohi[ix]   = csthilohi[i];
      datalolohi[ix]   = cstlolohi[i];
      datahihilo[ix]   = csthihilo[i];
      datalohilo[ix]   = cstlohilo[i];
      datahilolo[ix]   = csthilolo[i];
      datalololo[ix++] = cstlololo[i];
   }
   for(int i=0; i<nbr; i++)
      for(int j=0; j<deg1; j++)
      {
         datahihihi[ix]   = cffhihihi[i][j];
         datalohihi[ix]   = cfflohihi[i][j];
         datahilohi[ix]   = cffhilohi[i][j];
         datalolohi[ix]   = cfflolohi[i][j];
         datahihilo[ix]   = cffhihilo[i][j];
         datalohilo[ix]   = cfflohilo[i][j];
         datahilolo[ix]   = cffhilolo[i][j];
         datalololo[ix++] = cfflololo[i][j];
      }
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++)
      {
         datahihihi[ix]   = inputhihihi[i][j];
         datalohihi[ix]   = inputlohihi[i][j];
         datahilohi[ix]   = inputhilohi[i][j];
         datalolohi[ix]   = inputlolohi[i][j];
         datahihilo[ix]   = inputhihilo[i][j];
         datalohilo[ix]   = inputlohilo[i][j];
         datahilolo[ix]   = inputhilolo[i][j];
         datalololo[ix++] = inputlololo[i][j];
      }

   for(int i=ix; i<totcff; i++)
   {
      datahihihi[i] = 0.0; datalohihi[i] = 0.0;
      datahilohi[i] = 0.0; datalolohi[i] = 0.0;
      datahihilo[i] = 0.0; datalohilo[i] = 0.0;
      datahilolo[i] = 0.0; datalololo[i] = 0.0;
   }
}

void cmplx8vectorized_data_setup
 ( int dim, int nbr, int deg, int totcff, int offsetri,
   double *datarihihihi, double *datarilohihi,
   double *datarihilohi, double *datarilolohi,
   double *datarihihilo, double *datarilohilo,
   double *datarihilolo, double *datarilololo,
   double *cstrehihihi, double *cstrelohihi,
   double *cstrehilohi, double *cstrelolohi,
   double *cstrehihilo, double *cstrelohilo,
   double *cstrehilolo, double *cstrelololo,
   double *cstimhihihi, double *cstimlohihi,
   double *cstimhilohi, double *cstimlolohi,
   double *cstimhihilo, double *cstimlohilo,
   double *cstimhilolo, double *cstimlololo,
   double **cffrehihihi, double **cffrelohihi,
   double **cffrehilohi, double **cffrelolohi,
   double **cffrehihilo, double **cffrelohilo,
   double **cffrehilolo, double **cffrelololo,
   double **cffimhihihi, double **cffimlohihi,
   double **cffimhilohi, double **cffimlolohi,
   double **cffimhihilo, double **cffimlohilo,
   double **cffimhilolo, double **cffimlololo,
   double **inputrehihihi, double **inputrelohihi,
   double **inputrehilohi, double **inputrelolohi,
   double **inputrehihilo, double **inputrelohilo,
   double **inputrehilolo, double **inputrelololo,
   double **inputimhihihi, double **inputimlohihi,
   double **inputimhilohi, double **inputimlolohi,
   double **inputimhihilo, double **inputimlohilo,
   double **inputimhilolo, double **inputimlololo )
{
   const int deg1 = deg+1;

   int ix1 = 0;
   int ix2 = totcff + offsetri;
 
   // cout << "before copying the constant ";
   // cout << "ix1 : " << ix1 << ", ix2 : " << ix2 << endl;

   for(int i=0; i<deg1; i++)
   {
      datarihihihi[ix1]   = cstrehihihi[i];
      datarilohihi[ix1]   = cstrelohihi[i];
      datarihilohi[ix1]   = cstrehilohi[i];
      datarilolohi[ix1]   = cstrelolohi[i];
      datarihihilo[ix1]   = cstrehihilo[i];
      datarilohilo[ix1]   = cstrelohilo[i];
      datarihilolo[ix1]   = cstrehilolo[i];
      datarilololo[ix1++] = cstrelololo[i];
      datarihihihi[ix2]   = cstimhihihi[i];
      datarilohihi[ix2]   = cstimlohihi[i];
      datarihilohi[ix2]   = cstimhilohi[i];
      datarilolohi[ix2]   = cstimlolohi[i];
      datarihihilo[ix2]   = cstimhihilo[i];
      datarilohilo[ix2]   = cstimlohilo[i];
      datarihilolo[ix2]   = cstimhilolo[i];
      datarilololo[ix2++] = cstimlololo[i];
   }

   // cout << "before copying the coefficients ";
   // cout << "ix1 : " << ix1 << ", ix2 : " << ix2 << endl;

   for(int i=0; i<nbr; i++)
      for(int j=0; j<deg1; j++)
      {
         datarihihihi[ix1]   = cffrehihihi[i][j];
         datarilohihi[ix1]   = cffrelohihi[i][j];
         datarihilohi[ix1]   = cffrehilohi[i][j];
         datarilolohi[ix1]   = cffrelolohi[i][j];
         datarihihilo[ix1]   = cffrehihilo[i][j];
         datarilohilo[ix1]   = cffrelohilo[i][j];
         datarihilolo[ix1]   = cffrehilolo[i][j];
         datarilololo[ix1++] = cffrelololo[i][j];
         datarihihihi[ix2]   = cffimhihihi[i][j];
         datarilohihi[ix2]   = cffimlohihi[i][j];
         datarihilohi[ix2]   = cffimhilohi[i][j];
         datarilolohi[ix2]   = cffimlolohi[i][j];
         datarihihilo[ix2]   = cffimhihilo[i][j];
         datarilohilo[ix2]   = cffimlohilo[i][j];
         datarihilolo[ix2]   = cffimhilolo[i][j];
         datarilololo[ix2++] = cffimlololo[i][j];
      }

   // cout << "before copying the input series ";
   // cout << "ix1 : " << ix1 << ", ix2 : " << ix2 << endl;

   for(int i=0; i<dim; i++)
   {
      for(int j=0; j<deg1; j++)
      {
         datarihihihi[ix1]   = inputrehihihi[i][j];
         datarilohihi[ix1]   = inputrelohihi[i][j];
         datarihilohi[ix1]   = inputrehilohi[i][j];
         datarilolohi[ix1]   = inputrelolohi[i][j];
         datarihihilo[ix1]   = inputrehihilo[i][j];
         datarilohilo[ix1]   = inputrelohilo[i][j];
         datarihilolo[ix1]   = inputrehilolo[i][j];
         datarilololo[ix1++] = inputrelololo[i][j];
         datarihihihi[ix2]   = inputimhihihi[i][j];
         datarilohihi[ix2]   = inputimlohihi[i][j];
         datarihilohi[ix2]   = inputimhilohi[i][j];
         datarilolohi[ix2]   = inputimlolohi[i][j];
         datarihihilo[ix2]   = inputimhihilo[i][j];
         datarilohilo[ix2]   = inputimlohilo[i][j];
         datarihilolo[ix2]   = inputimhilolo[i][j];
         datarilololo[ix2++] = inputimlololo[i][j];
      }
   }
   // cout << "Initializing data with offsetri = " << offsetri
   //      << ", ix1 = " << ix1 << ", ix2 = " << ix2 << endl;

   for(int i=0; i<2*offsetri; i++)  // 2*offsetri because of 2nd operands
   {
      datarihihihi[ix1]   = 0.0;
      datarilohihi[ix1]   = 0.0;
      datarihilohi[ix1]   = 0.0;
      datarilolohi[ix1]   = 0.0;
      datarihihilo[ix1]   = 0.0;
      datarilohilo[ix1]   = 0.0;
      datarihilolo[ix1]   = 0.0;
      datarilololo[ix1++] = 0.0;
      datarihihihi[ix2]   = 0.0;
      datarilohihi[ix2]   = 0.0;
      datarihilohi[ix2]   = 0.0;
      datarilolohi[ix2]   = 0.0;
      datarihihilo[ix2]   = 0.0;
      datarilohilo[ix2]   = 0.0;
      datarihilolo[ix2]   = 0.0;
      datarilololo[ix2++] = 0.0;
   }
   // cout << "After initializing, ix1 = " << ix1 << ", ix2 = " << ix2 << endl;
}

void dbl8_convolution_jobs
 ( int dim, int nbr, int deg, int *nvr, ConvolutionJobs cnvjobs,
   int *fstart, int *bstart, int *cstart,
   double *datahihihi, double *datalohihi,
   double *datahilohi, double *datalolohi,
   double *datahihilo, double *datalohilo,
   double *datahilolo, double *datalololo,
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
         dbl8_padded_convjobs<<<jobnbr,deg1>>>
            (datahihihi,datalohihi,datahilohi,datalolohi,
             datahihilo,datalohilo,datahilolo,datalololo,
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

void cmplx8vectorized_convolution_jobs
 ( int dim, int nbr, int deg, int *nvr, int totcff, int offsetri,
   ComplexConvolutionJobs cnvjobs, ComplexIncrementJobs incjobs,
   int *fstart, int *bstart, int *cstart,
   double *datarihihihi, double *datarilohihi,
   double *datarihilohi, double *datarilolohi,
   double *datarihihilo, double *datarilohilo,
   double *datarihilolo, double *datarilololo,
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
/*
         for(int i=0; i<jobnbr; i++)
         {
            int in1idx = in1ix_h[i];

            cout << "First input of job " << i
                 << " at index " << in1idx << " :" << endl;

            double rihihihi,rilohihi,rihilohi,rilolohi;
            double rihihilo,rilohilo,rihilolo,rilololo;

            cudaMemcpy(&rihihihi,&datarihihihi[in1idx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rilohihi,&datarilohihi[in1idx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rihilohi,&datarihilohi[in1idx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rilolohi,&datarilolohi[in1idx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rihihilo,&datarihihilo[in1idx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rilohilo,&datarilohilo[in1idx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rihilolo,&datarihilolo[in1idx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rilololo,&datarilololo[in1idx],sizeof(double),
                       cudaMemcpyDeviceToHost);

            cout << scientific << setprecision(16);
            cout << rihihihi << "  " << rilohihi << endl;
            cout << rihilohi << "  " << rilolohi << endl;
            cout << rihihilo << "  " << rilohilo << endl;
            cout << rihilolo << "  " << rilololo << endl;

            int in2idx = in2ix_h[i];

            cout << "Second input of job " << i
                 << " at index " << in2idx << " :" << endl;

            cudaMemcpy(&rihihihi,&datarihihihi[in2idx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rilohihi,&datarilohihi[in2idx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rihilohi,&datarihilohi[in2idx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rilolohi,&datarilolohi[in2idx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rihihilo,&datarihihilo[in2idx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rilohilo,&datarilohilo[in2idx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rihilolo,&datarihilolo[in2idx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rilololo,&datarilololo[in2idx],sizeof(double),
                       cudaMemcpyDeviceToHost);

            cout << scientific << setprecision(16);
            cout << rihihihi << "  " << rilohihi << endl;
            cout << rihilohi << "  " << rilolohi << endl;
            cout << rihihilo << "  " << rilohilo << endl;
            cout << rihilolo << "  " << rilololo << endl;
         }
 */
         cudaEventRecord(start);
         dbl8_padded_convjobs<<<jobnbr,deg1>>>
            (datarihihihi,datarilohihi,datarihilohi,datarilolohi,
             datarihihilo,datarilohilo,datarihilolo,datarilololo,
             in1ix_d,in2ix_d,outix_d,deg1);
         cudaEventRecord(stop);
         cudaEventSynchronize(stop);
         cudaEventElapsedTime(&milliseconds,start,stop);
         *cnvlapms += milliseconds;
/*
         for(int i=0; i<jobnbr; i++)
         {
            int outidx = outix_h[i];

            cout << "Output of job " << i
                 << " at index " << outidx << " :" << endl;

            double rihihihi,rilohihi,rihilohi,rilolohi;
            double rihihilo,rilohilo,rihilolo,rilololo;

            cudaMemcpy(&rihihihi,&datarihihihi[outidx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rilohihi,&datarilohihi[outidx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rihilohi,&datarihilohi[outidx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rilolohi,&datarilolohi[outidx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rihihilo,&datarihihilo[outidx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rilohilo,&datarilohilo[outidx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rihilolo,&datarihilolo[outidx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rilololo,&datarilololo[outidx],sizeof(double),
                       cudaMemcpyDeviceToHost);

            cout << rihihihi << "  " << rilohihi << endl;
            cout << rihilohi << "  " << rilolohi << endl;
            cout << rihihilo << "  " << rilohilo << endl;
            cout << rihilolo << "  " << rilololo << endl;
         }
 */
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

      GPU_cmplx8vectorized_flipsigns
        (deg,nbrflips,rebidx,
         datarihihihi,datarilohihi,datarihilohi,datarilolohi,
         datarihihilo,datarilohilo,datarihilolo,datarilololo,
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
/*
         for(int i=0; i<jobnbr; i++)
         {
            int in1idx = in1ix_h[i];

            cout << "First input of job " << i
                 << " at index " << in1idx << " :" << endl;

            double rihihihi,rilohihi,rihilohi,rilolohi;
            double rihihilo,rilohilo,rihilolo,rilololo;

            cudaMemcpy(&rihihihi,&datarihihihi[in1idx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rilohihi,&datarilohihi[in1idx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rihilohi,&datarihilohi[in1idx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rilolohi,&datarilolohi[in1idx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rihihilo,&datarihihilo[in1idx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rilohilo,&datarilohilo[in1idx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rihilolo,&datarihilolo[in1idx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rilololo,&datarilololo[in1idx],sizeof(double),
                       cudaMemcpyDeviceToHost);

            cout << scientific << setprecision(16);
            cout << rihihihi << "  " << rilohihi << endl;
            cout << rihilohi << "  " << rilolohi << endl;
            cout << rihihilo << "  " << rilohilo << endl;
            cout << rihilolo << "  " << rilololo << endl;

            int in2idx = in2ix_h[i];

            cout << "Second input of job " << i
                 << " at index " << in2idx << " :" << endl;

            cudaMemcpy(&rihihihi,&datarihihihi[in2idx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rilohihi,&datarilohihi[in2idx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rihilohi,&datarihilohi[in2idx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rilolohi,&datarilolohi[in2idx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rihihilo,&datarihihilo[in2idx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rilohilo,&datarilohilo[in2idx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rihilolo,&datarihilolo[in2idx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rilololo,&datarilololo[in2idx],sizeof(double),
                       cudaMemcpyDeviceToHost);

            cout << scientific << setprecision(16);
            cout << rihihihi << "  " << rilohihi << endl;
            cout << rihilohi << "  " << rilolohi << endl;
            cout << rihihilo << "  " << rilohilo << endl;
            cout << rihilolo << "  " << rilololo << endl;
         }
 */
         cudaEventRecord(start);
         dbl8_increment_jobs<<<jobnbr,deg1>>>
            (datarihihihi,datarilohihi,datarihilohi,datarilolohi,
             datarihihilo,datarilohilo,datarihilolo,datarilololo,
             in1ix_d,in2ix_d,outix_d,deg1);
         cudaEventRecord(stop);
         cudaEventSynchronize(stop);
         cudaEventElapsedTime(&milliseconds,start,stop);
         *cnvlapms += milliseconds;
/*
         for(int i=0; i<jobnbr; i++)
         {
            int outidx = outix_h[i];

            cout << "Output of job " << i
                 << " at index " << outidx << " :" << endl;

            double rihihihi,rilohihi,rihilohi,rilolohi;
            double rihihilo,rilohilo,rihilolo,rilololo;

            cudaMemcpy(&rihihihi,&datarihihihi[outidx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rilohihi,&datarilohihi[outidx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rihilohi,&datarihilohi[outidx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rilolohi,&datarilolohi[outidx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rihihilo,&datarihihilo[outidx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rilohilo,&datarilohilo[outidx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rihilolo,&datarihilolo[outidx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rilololo,&datarilololo[outidx],sizeof(double),
                       cudaMemcpyDeviceToHost);

            cout << scientific << setprecision(16);
            cout << rihihihi << "  " << rilohihi << endl;
            cout << rihilohi << "  " << rilolohi << endl;
            cout << rihihilo << "  " << rilohilo << endl;
            cout << rihilolo << "  " << rilololo << endl;
         }
 */
         cudaFree(in1ix_d); cudaFree(in2ix_d); cudaFree(outix_d);
      }
      free(in1ix_h); free(in2ix_h); free(outix_h);
      free(rebidx);
   }
}

void dbl8_addition_jobs
 ( int dim, int nbr, int deg, int *nvr, AdditionJobs addjobs,
   int *fstart, int *bstart, int *cstart,
   double *datahihihi, double *datalohihi,
   double *datahilohi, double *datalolohi,
   double *datahihilo, double *datalohilo,
   double *datahilolo, double *datalololo, double *addlapms, int vrblvl )
{
   const int deg1 = deg+1;

   cudaEvent_t start,stop;           // to measure time spent by kernels
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
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
                 << " threads ..." << endl;

         cudaEventRecord(start);
         dbl8_update_addjobs<<<jobnbr,deg1>>>
            (datahihihi,datalohihi,datahilohi,datalolohi,
             datahihilo,datalohilo,datahilolo,datalololo,
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

void cmplx8vectorized_addition_jobs
 ( int dim, int nbr, int deg, int *nvr, int totcff, int offsetri,
   ComplexAdditionJobs addjobs,
   int *fstart, int *bstart, int *cstart,
   double *datarihihihi, double *datarilohihi,
   double *datarihilohi, double *datarilolohi,
   double *datarihihilo, double *datarilohilo,
   double *datarihilolo, double *datarilololo,
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
/*
         for(int i=0; i<jobnbr; i++)
         {
            int in1idx = in1ix_h[i];

            cout << "First input of job " << i
                 << " at index " << in1idx << " :" << endl;

            double rihihihi,rilohihi,rihilohi,rilolohi;
            double rihihilo,rilohilo,rihilolo,rilololo;

            cudaMemcpy(&rihihihi,&datarihihihi[in1idx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rilohihi,&datarilohihi[in1idx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rihilohi,&datarihilohi[in1idx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rilolohi,&datarilolohi[in1idx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rihihilo,&datarihihilo[in1idx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rilohilo,&datarilohilo[in1idx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rihilolo,&datarihilolo[in1idx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rilololo,&datarilololo[in1idx],sizeof(double),
                       cudaMemcpyDeviceToHost);

            cout << scientific << setprecision(16);
            cout << rihihihi << "  " << rilohihi << endl;
            cout << rihilohi << "  " << rilolohi << endl;
            cout << rihihilo << "  " << rilohilo << endl;
            cout << rihilolo << "  " << rilololo << endl;

            int in2idx = in2ix_h[i];

            cout << "Second input of job " << i
                 << " at index " << in2idx << " :" << endl;

            cudaMemcpy(&rihihihi,&datarihihihi[in2idx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rilohihi,&datarilohihi[in2idx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rihilohi,&datarihilohi[in2idx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rilolohi,&datarilolohi[in2idx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rihihilo,&datarihihilo[in2idx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rilohilo,&datarilohilo[in2idx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rihilolo,&datarihilolo[in2idx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rilololo,&datarilololo[in2idx],sizeof(double),
                       cudaMemcpyDeviceToHost);

            cout << scientific << setprecision(16);
            cout << rihihihi << "  " << rilohihi << endl;
            cout << rihilohi << "  " << rilolohi << endl;
            cout << rihihilo << "  " << rilohilo << endl;
            cout << rihilolo << "  " << rilololo << endl;
         }
 */
         cudaEventRecord(start);
         dbl8_update_addjobs<<<jobnbr,deg1>>>
            (datarihihihi,datarilohihi,datarihilohi,datarilolohi,
             datarihihilo,datarilohilo,datarihilolo,datarilololo,
             in1ix_d,in2ix_d,outix_d,deg1);
         cudaEventRecord(stop);
         cudaEventSynchronize(stop);
         cudaEventElapsedTime(&milliseconds,start,stop);
         *addlapms += milliseconds;
/*
         for(int i=0; i<jobnbr; i++)
         {
            int outidx = outix_h[i];

            cout << "Output of job " << i
                 << " at index " << outidx << " :" << endl;

            double rihihihi,rilohihi,rihilohi,rilolohi;
            double rihihilo,rilohilo,rihilolo,rilololo;

            cudaMemcpy(&rihihihi,&datarihihihi[outidx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rilohihi,&datarilohihi[outidx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rihilohi,&datarihilohi[outidx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rilolohi,&datarilolohi[outidx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rihihilo,&datarihihilo[outidx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rilohilo,&datarilohilo[outidx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rihilolo,&datarihilolo[outidx],sizeof(double),
                       cudaMemcpyDeviceToHost);
            cudaMemcpy(&rilololo,&datarilololo[outidx],sizeof(double),
                       cudaMemcpyDeviceToHost);

            cout << scientific << setprecision(16);
            cout << rihihihi << "  " << rilohihi << endl;
            cout << rihilohi << "  " << rilolohi << endl;
            cout << rihihilo << "  " << rilohilo << endl;
            cout << rihilolo << "  " << rilololo << endl;
         }
 */
         cudaFree(in1ix_d); cudaFree(in2ix_d); cudaFree(outix_d);
      }
      free(in1ix_h); free(in2ix_h); free(outix_h);
   }
}

void GPU_dbl8_poly_evaldiff
 ( int BS, int dim, int nbr, int deg, int *nvr, int **idx,
   double *csthihihi, double *cstlohihi,
   double *csthilohi, double *cstlolohi,
   double *csthihilo, double *cstlohilo,
   double *csthilolo, double *cstlololo,
   double **cffhihihi, double **cfflohihi,
   double **cffhilohi, double **cfflolohi,
   double **cffhihilo, double **cfflohilo,
   double **cffhilolo, double **cfflololo,
   double **inputhihihi, double **inputlohihi,
   double **inputhilohi, double **inputlolohi,
   double **inputhihilo, double **inputlohilo,
   double **inputhilolo, double **inputlololo,
   double **outputhihihi, double **outputlohihi,
   double **outputhilohi, double **outputlolohi,
   double **outputhihilo, double **outputlohilo,
   double **outputhilolo, double **outputlololo,
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
   {
      cout << "The output coefficient count : " << totalcff << endl;
      cout << "fsums :";
      for(int i=0; i<nbr; i++) cout << " " << fsums[i]; cout << endl;
      cout << "fstart :";
      for(int i=0; i<nbr; i++) cout << " " << fstart[i]; cout << endl;
      cout << "bsums :";
      for(int i=0; i<nbr; i++) cout << " " << bsums[i]; cout << endl;
      cout << "bstart :";
      for(int i=0; i<nbr; i++) cout << " " << bstart[i]; cout << endl;
      cout << "csums :";
      for(int i=0; i<nbr; i++) cout << " " << csums[i]; cout << endl;
      cout << "cstart :";
      for(int i=0; i<nbr; i++) cout << " " << cstart[i]; cout << endl;
   }
   double *datahihihi_h = new double[totalcff];        // data on host
   double *datalohihi_h = new double[totalcff];
   double *datahilohi_h = new double[totalcff];
   double *datalolohi_h = new double[totalcff];
   double *datahihilo_h = new double[totalcff];
   double *datalohilo_h = new double[totalcff];
   double *datahilolo_h = new double[totalcff];
   double *datalololo_h = new double[totalcff];

   dbl8_data_setup
      (dim,nbr,deg,totalcff,
       datahihihi_h,datalohihi_h,datahilohi_h,datalolohi_h,
       datahihilo_h,datalohilo_h,datahilolo_h,datalololo_h,
       csthihihi,cstlohihi,csthilohi,cstlolohi,
       csthihilo,cstlohilo,csthilolo,cstlololo,
       cffhihihi,cfflohihi,cffhilohi,cfflolohi,
       cffhihilo,cfflohilo,cffhilolo,cfflololo,
       inputhihihi,inputlohihi,inputhilohi,inputlolohi,
       inputhihilo,inputlohilo,inputhilolo,inputlololo);

   double *datahihihi_d;                               // device data
   double *datalohihi_d;
   double *datahilohi_d;
   double *datalolohi_d;
   double *datahihilo_d;
   double *datalohilo_d;
   double *datahilolo_d;
   double *datalololo_d;
   const size_t szdata = totalcff*sizeof(double);
   cudaMalloc((void**)&datahihihi_d,szdata);
   cudaMalloc((void**)&datalohihi_d,szdata);
   cudaMalloc((void**)&datahilohi_d,szdata);
   cudaMalloc((void**)&datalolohi_d,szdata);
   cudaMalloc((void**)&datahihilo_d,szdata);
   cudaMalloc((void**)&datalohilo_d,szdata);
   cudaMalloc((void**)&datahilolo_d,szdata);
   cudaMalloc((void**)&datalololo_d,szdata);
   cudaMemcpy(datahihihi_d,datahihihi_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datalohihi_d,datalohihi_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datahilohi_d,datahilohi_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datalolohi_d,datalolohi_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datahihilo_d,datahihilo_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datalohilo_d,datalohilo_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datahilolo_d,datahilolo_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datalololo_d,datalololo_h,szdata,cudaMemcpyHostToDevice);

   struct timeval begintime,endtime; // wall clock time of computations

   gettimeofday(&begintime,0);

   dbl8_convolution_jobs
      (dim,nbr,deg,nvr,cnvjobs,fstart,bstart,cstart,
       datahihihi_d,datalohihi_d,datahilohi_d,datalolohi_d,
       datahihilo_d,datalohilo_d,datahilolo_d,datalololo_d,
       cnvlapms,vrblvl);

   dbl8_addition_jobs
      (dim,nbr,deg,nvr,addjobs,fstart,bstart,cstart,
       datahihihi_d,datalohihi_d,datahilohi_d,datalolohi_d,
       datahihilo_d,datalohilo_d,datahilolo_d,datalololo_d,
       addlapms,vrblvl);

   gettimeofday(&endtime,0);
   cudaMemcpy(datahihihi_h,datahihihi_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datalohihi_h,datalohihi_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datahilohi_h,datahilohi_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datalolohi_h,datalolohi_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datahihilo_h,datahihilo_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datalohilo_h,datalohilo_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datahilolo_h,datahilolo_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datalololo_h,datalololo_d,szdata,cudaMemcpyDeviceToHost);
   *elapsedms = *cnvlapms + *addlapms;
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   *walltimesec = seconds + microseconds*1.0e-6;

   dbl_added_data8_to_output
      (datahihihi_h,datalohihi_h,datahilohi_h,datalolohi_h,
       datahihilo_h,datalohilo_h,datahilolo_h,datalololo_h,
       outputhihihi,outputlohihi,outputhilohi,outputlolohi,
       outputhihilo,outputlohilo,outputhilolo,outputlololo,
       dim,nbr,deg,nvr,idx,fstart,bstart,cstart,addjobs,vrblvl);

   if(vrblvl > 0)
      write_GPU_timings(*cnvlapms,*addlapms,*elapsedms,*walltimesec);

   cudaFree(datahihihi_d); cudaFree(datalohihi_d);
   cudaFree(datahilohi_d); cudaFree(datalolohi_d);
   cudaFree(datahihilo_d); cudaFree(datalohilo_d);
   cudaFree(datahilolo_d); cudaFree(datalololo_d);

   free(datahihihi_h); free(datalohihi_h);
   free(datahilohi_h); free(datalolohi_h);
   free(datahihilo_h); free(datalohilo_h);
   free(datahilolo_h); free(datalololo_h);

   free(fstart); free(bstart); free(cstart);
   free(fsums); free(bsums); free(csums);
}

void GPU_cmplx8vectorized_poly_evaldiff
 ( int BS, int dim, int nbr, int deg, int *nvr, int **idx,
   double *cstrehihihi, double *cstrelohihi,
   double *cstrehilohi, double *cstrelolohi,
   double *cstrehihilo, double *cstrelohilo,
   double *cstrehilolo, double *cstrelololo,
   double *cstimhihihi, double *cstimlohihi,
   double *cstimhilohi, double *cstimlolohi,
   double *cstimhihilo, double *cstimlohilo,
   double *cstimhilolo, double *cstimlololo,
   double **cffrehihihi, double **cffrelohihi,
   double **cffrehilohi, double **cffrelolohi,
   double **cffrehihilo, double **cffrelohilo,
   double **cffrehilolo, double **cffrelololo,
   double **cffimhihihi, double **cffimlohihi,
   double **cffimhilohi, double **cffimlolohi,
   double **cffimhihilo, double **cffimlohilo,
   double **cffimhilolo, double **cffimlololo,
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
   double **outputimhilolo, double **outputimlololo,
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
   double *datarihihihi_h = new double[cmplxtotcff];      // data on host
   double *datarilohihi_h = new double[cmplxtotcff];
   double *datarihilohi_h = new double[cmplxtotcff];
   double *datarilolohi_h = new double[cmplxtotcff];
   double *datarihihilo_h = new double[cmplxtotcff];
   double *datarilohilo_h = new double[cmplxtotcff];
   double *datarihilolo_h = new double[cmplxtotcff];
   double *datarilololo_h = new double[cmplxtotcff];

   cmplx8vectorized_data_setup
      (dim,nbr,deg,totalcff,offsetri,
       datarihihihi_h,datarilohihi_h,datarihilohi_h,datarilolohi_h,
       datarihihilo_h,datarilohilo_h,datarihilolo_h,datarilololo_h,
       cstrehihihi,cstrelohihi,cstrehilohi,cstrelolohi,
       cstrehihilo,cstrelohilo,cstrehilolo,cstrelololo,
       cstimhihihi,cstimlohihi,cstimhilohi,cstimlolohi,
       cstimhihilo,cstimlohilo,cstimhilolo,cstimlololo,
       cffrehihihi,cffrelohihi,cffrehilohi,cffrelolohi,
       cffrehihilo,cffrelohilo,cffrehilolo,cffrelololo,
       cffimhihihi,cffimlohihi,cffimhilohi,cffimlolohi,
       cffimhihilo,cffimlohilo,cffimhilolo,cffimlololo,
       inputrehihihi,inputrelohihi,inputrehilohi,inputrelolohi,
       inputrehihilo,inputrelohilo,inputrehilolo,inputrelololo,
       inputimhihihi,inputimlohihi,inputimhilohi,inputimlolohi,
       inputimhihilo,inputimlohilo,inputimhilolo,inputimlololo);

   double *datarihihihi_d;                               // device data
   double *datarilohihi_d;
   double *datarihilohi_d;
   double *datarilolohi_d;
   double *datarihihilo_d;
   double *datarilohilo_d;
   double *datarihilolo_d;
   double *datarilololo_d;
   const size_t szdata = cmplxtotcff*sizeof(double);
   cudaMalloc((void**)&datarihihihi_d,szdata);
   cudaMalloc((void**)&datarilohihi_d,szdata);
   cudaMalloc((void**)&datarihilohi_d,szdata);
   cudaMalloc((void**)&datarilolohi_d,szdata);
   cudaMalloc((void**)&datarihihilo_d,szdata);
   cudaMalloc((void**)&datarilohilo_d,szdata);
   cudaMalloc((void**)&datarihilolo_d,szdata);
   cudaMalloc((void**)&datarilololo_d,szdata);
   cudaMemcpy(datarihihihi_d,datarihihihi_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datarilohihi_d,datarilohihi_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datarihilohi_d,datarihilohi_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datarilolohi_d,datarilolohi_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datarihihilo_d,datarihihilo_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datarilohilo_d,datarilohilo_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datarihilolo_d,datarihilolo_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datarilololo_d,datarilololo_h,szdata,cudaMemcpyHostToDevice);

   struct timeval begintime,endtime; // wall clock time of computations

   gettimeofday(&begintime,0);

   cmplx8vectorized_convolution_jobs
      (dim,nbr,deg,nvr,totalcff,offsetri,cnvjobs,incjobs,fstart,bstart,cstart,
       datarihihihi_d,datarilohihi_d,datarihilohi_d,datarilolohi_d,
       datarihihilo_d,datarilohilo_d,datarihilolo_d,datarilololo_d,
       cnvlapms,vrblvl);

   cmplx8vectorized_addition_jobs
      (dim,nbr,deg,nvr,totalcff,offsetri,addjobs,fstart,bstart,cstart,
       datarihihihi_d,datarilohihi_d,datarihilohi_d,datarilolohi_d,
       datarihihilo_d,datarilohilo_d,datarihilolo_d,datarilololo_d,
       addlapms,vrblvl);

   gettimeofday(&endtime,0);
   cudaMemcpy(datarihihihi_h,datarihihihi_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datarilohihi_h,datarilohihi_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datarihilohi_h,datarihilohi_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datarilolohi_h,datarilolohi_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datarihihilo_h,datarihihilo_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datarilohilo_h,datarilohilo_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datarihilolo_h,datarihilolo_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datarilololo_h,datarilololo_d,szdata,cudaMemcpyDeviceToHost);
   *elapsedms = *cnvlapms + *addlapms;
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   *walltimesec = seconds + microseconds*1.0e-6;

   cmplx_added_data8vectorized_to_output
      (datarihihihi_h,datarilohihi_h,datarihilohi_h,datarilolohi_h,
       datarihihilo_h,datarilohilo_h,datarihilolo_h,datarilololo_h,
       outputrehihihi,outputrelohihi,outputrehilohi,outputrelolohi,
       outputrehihilo,outputrelohilo,outputrehilolo,outputrelololo,
       outputimhihihi,outputimlohihi,outputimhilohi,outputimlolohi,
       outputimhihilo,outputimlohilo,outputimhilolo,outputimlololo,
       dim,nbr,deg,nvr,idx,fstart,bstart,cstart,
       totalcff,offsetri,addjobs,vrblvl);

   if(vrblvl > 0)
      write_GPU_timings(*cnvlapms,*addlapms,*elapsedms,*walltimesec);

   cudaFree(datarihihihi_d); cudaFree(datarilohihi_d);
   cudaFree(datarihilohi_d); cudaFree(datarilolohi_d);
   cudaFree(datarihihilo_d); cudaFree(datarilohilo_d);
   cudaFree(datarihilolo_d); cudaFree(datarilololo_d);

   free(datarihihihi_h); free(datarilohihi_h);
   free(datarihilohi_h); free(datarilolohi_h);
   free(datarihihilo_h); free(datarilohilo_h);
   free(datarihilolo_h); free(datarilololo_h);

   free(fstart); free(bstart); free(cstart);
   free(fsums); free(bsums); free(csums);
}
