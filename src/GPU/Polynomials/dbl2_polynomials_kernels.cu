// The file dbl2_polynomials_kernels.cu defines the kernels with prototypes
// in dbl2_polynomials_kernels.h.

#include <iostream>
#include <iomanip>
#ifdef winwalltime
#include "gettimeofday4win.h"
#else
#include <sys/time.h>
#endif
#include "job_coordinates.h"
#include "double_double_functions.h"
#ifdef gpufun
#include "double_double_gpufun.cu"
#endif
#include "dbl2_polynomials_kernels.h"
#include "write_gpu_timings.h"

// The constant dd_shmemsize is the bound on the shared memory size.

#define dd_shmemsize 192

using namespace std;

__global__ void dbl2_padded_convjobs
 ( double *datahi, double *datalo,
   int *in1idx, int *in2idx, int *outidx, int dim )
{
   const int bdx = blockIdx.x;           // index to the convolution job
   const int tdx = threadIdx.x;          // index to the output of the job
   const int idx1 = in1idx[bdx] + tdx;
   const int idx2 = in2idx[bdx] + tdx;
   const int idx3 = outidx[bdx] + tdx;

   __shared__ double xvhi[dd_shmemsize];
   __shared__ double xvlo[dd_shmemsize];
   __shared__ double yvhi[2*dd_shmemsize];
   __shared__ double yvlo[2*dd_shmemsize];
   __shared__ double zvhi[dd_shmemsize];
   __shared__ double zvlo[dd_shmemsize];

   double prdhi,prdlo;
   int ydx = dim + tdx;

   xvhi[tdx] = datahi[idx1];  // loading first input
   xvlo[tdx] = datalo[idx1]; 
   yvhi[tdx] = 0.0;           // padded with zeros
   yvlo[tdx] = 0.0;
   yvhi[ydx] = datahi[idx2];  // loading second input
   yvlo[ydx] = datalo[idx2];

   __syncthreads();

   // zv[tdx] = xv[0]*yv[tdx];
   ddg_mul(xvhi[0],xvlo[0],yvhi[ydx],yvlo[ydx],&zvhi[tdx],&zvlo[tdx]);
   __syncthreads();

   for(int i=1; i<dim; i++) // zv[tdx] = zv[tdx] + xv[i]*yv[dim+tdx-i];
   {
      ydx = dim + tdx - i;
      ddg_mul(xvhi[i],xvlo[i],yvhi[ydx],yvlo[ydx],&prdhi,&prdlo);
      __syncthreads();

      ddg_inc(&zvhi[tdx],&zvlo[tdx],prdhi,prdlo);
      __syncthreads();
   }
   __syncthreads();

   datahi[idx3] = zvhi[tdx]; // storing the output
   datalo[idx3] = zvlo[tdx];
}

__global__ void dbl2_increment_jobs
 ( double *datahi, double *datalo,
   int *in1idx, int *in2idx, int *outidx, int dim )
{
   const int bdx = blockIdx.x;           // index to the increment job
   const int tdx = threadIdx.x;          // index to the output of the job
   const int idx1 = in1idx[bdx] + tdx;
   const int idx2 = in2idx[bdx] + tdx;
   const int idx3 = outidx[bdx] + tdx;

   __shared__ double yvhi[dd_shmemsize];
   __shared__ double yvlo[dd_shmemsize];
   __shared__ double zvhi[dd_shmemsize];
   __shared__ double zvlo[dd_shmemsize];

   zvhi[tdx] = datahi[idx1];  // loading first input
   zvlo[tdx] = datalo[idx1]; 
   yvhi[tdx] = datahi[idx2];  // loading second input
   yvlo[tdx] = datalo[idx2];

   __syncthreads();

   ddg_inc(&zvhi[tdx],&zvlo[tdx],yvhi[tdx],yvlo[tdx]);

   __syncthreads();

   datahi[idx3] = zvhi[tdx]; // storing the output
   datalo[idx3] = zvlo[tdx];
}

__global__ void cmplx2_padded_convjobs
 ( double *datarehi, double *datarelo,
   double *dataimhi, double *dataimlo,
   int *in1idx, int *in2idx, int *outidx, int dim )
{
   const int bdx = blockIdx.x;           // index to the convolution job
   const int tdx = threadIdx.x;          // index to the output of the job
   const int idx1 = in1idx[bdx] + tdx;
   const int idx2 = in2idx[bdx] + tdx;
   const int idx3 = outidx[bdx] + tdx;

   __shared__ double xvrehi[dd_shmemsize];
   __shared__ double xvrelo[dd_shmemsize];
   __shared__ double xvimhi[dd_shmemsize];
   __shared__ double xvimlo[dd_shmemsize];
   __shared__ double yvrehi[2*dd_shmemsize];
   __shared__ double yvrelo[2*dd_shmemsize];
   __shared__ double yvimhi[2*dd_shmemsize];
   __shared__ double yvimlo[2*dd_shmemsize];
   __shared__ double zvrehi[dd_shmemsize];
   __shared__ double zvrelo[dd_shmemsize];
   __shared__ double zvimhi[dd_shmemsize];
   __shared__ double zvimlo[dd_shmemsize];

   double prodhi,prodlo;
   int ydx = dim + tdx;

   xvrehi[tdx] = datarehi[idx1];  // loading first input
   xvrelo[tdx] = datarelo[idx1]; 
   xvimhi[tdx] = dataimhi[idx1];
   xvimlo[tdx] = dataimlo[idx1]; 
   yvrehi[tdx] = 0.0;           // padded with zeros
   yvrelo[tdx] = 0.0;
   yvimhi[tdx] = 0.0;
   yvimlo[tdx] = 0.0;
   yvrehi[ydx] = datarehi[idx2];  // loading second input
   yvrelo[ydx] = datarelo[idx2];
   yvimhi[ydx] = dataimhi[idx2];
   yvimlo[ydx] = dataimlo[idx2];

   __syncthreads();

   // zv[tdx] = xv[0]*yv[tdx];
   ddg_mul(xvrehi[0],xvrelo[0],yvrehi[ydx],yvrelo[ydx],
           &zvrehi[tdx],&zvrelo[tdx]);
   __syncthreads();
   ddg_mul(xvimhi[0],xvimlo[0],yvimhi[ydx],yvimlo[ydx],&prodhi,&prodlo);
   __syncthreads();
   ddg_dec(&zvrehi[tdx],&zvrelo[tdx],prodhi,prodlo);
   __syncthreads();

   ddg_mul(xvrehi[0],xvrelo[0],yvimhi[ydx],yvimlo[ydx],
           &zvimhi[tdx],&zvimlo[tdx]);
   __syncthreads();
   ddg_mul(xvimhi[0],xvimlo[0],yvrehi[ydx],yvrelo[ydx],&prodhi,&prodlo);
   __syncthreads();
   ddg_inc(&zvimhi[tdx],&zvimlo[tdx],prodhi,prodlo);
   __syncthreads();

   for(int i=1; i<dim; i++) // zv[tdx] = zv[tdx] + xv[i]*yv[dim+tdx-i];
   {
      ydx = dim + tdx - i;

      ddg_mul(xvrehi[i],xvrelo[i],yvrehi[ydx],yvrelo[ydx],&prodhi,&prodlo);
      __syncthreads();
      ddg_inc(&zvrehi[tdx],&zvrelo[tdx],prodhi,prodlo);
      __syncthreads();
      ddg_mul(xvimhi[i],xvimlo[i],yvimhi[ydx],yvimlo[ydx],&prodhi,&prodlo);
      __syncthreads();
      ddg_dec(&zvrehi[tdx],&zvrelo[tdx],prodhi,prodlo);
      __syncthreads();

      ddg_mul(xvrehi[i],xvrelo[i],yvimhi[ydx],yvimlo[ydx],&prodhi,&prodlo);
      __syncthreads();
      ddg_inc(&zvimhi[tdx],&zvimlo[tdx],prodhi,prodlo);
      __syncthreads();
      ddg_mul(xvimhi[i],xvimlo[i],yvrehi[ydx],yvrelo[ydx],&prodhi,&prodlo);
      __syncthreads();
      ddg_inc(&zvimhi[tdx],&zvimlo[tdx],prodhi,prodlo);
      __syncthreads();
   }
   __syncthreads();

   datarehi[idx3] = zvrehi[tdx]; // storing the output
   datarelo[idx3] = zvrelo[tdx];
   dataimhi[idx3] = zvimhi[tdx]; 
   dataimlo[idx3] = zvimlo[tdx];
}

__global__ void cmplx2vectorized_flipsigns
 ( double *datarihi, double *datarilo, int *flpidx, int dim )
{
   const int bdx = blockIdx.x;    // index to the series to flip
   const int tdx = threadIdx.x;
   const int idx = flpidx[bdx] + tdx; // which number to flip

   double x; // register to load data from global memory

   x = datarihi[idx]; x = -x; datarihi[idx] = x;
   x = datarilo[idx]; x = -x; datarilo[idx] = x;
}

void GPU_cmplx2vectorized_flipsigns
 ( int deg, int nbrflips, int *flipidx,
   double *datarihi, double *datarilo,
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
   cmplx2vectorized_flipsigns<<<nbrflips,deg1>>>
      (datarihi,datarilo,flipidx_d,deg1);
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

__global__ void dbl2_update_addjobs
 ( double *datahi, double *datalo,
   int *in1idx, int *in2idx, int *outidx, int dim )
{
   const int bdx = blockIdx.x;           // index to the addition job
   const int tdx = threadIdx.x;          // index to the output of the job
   const int idx1 = in1idx[bdx] + tdx;
   const int idx2 = in2idx[bdx] + tdx;
   const int idx3 = outidx[bdx] + tdx;

   __shared__ double xvhi[dd_shmemsize];
   __shared__ double xvlo[dd_shmemsize];
   __shared__ double yvhi[dd_shmemsize];
   __shared__ double yvlo[dd_shmemsize];
   __shared__ double zvhi[dd_shmemsize];
   __shared__ double zvlo[dd_shmemsize];

   xvhi[tdx] = datahi[idx1];  // loading first input
   xvlo[tdx] = datalo[idx1];
   yvhi[tdx] = datahi[idx2];  // loading second input
   yvlo[tdx] = datalo[idx2];

   // zv[tdx] = xv[tdx] + yv[tdx];

   ddg_add(xvhi[tdx],xvlo[tdx],yvhi[tdx],yvlo[tdx],&zvhi[tdx],&zvlo[tdx]);

   __syncthreads();

   datahi[idx3] = zvhi[tdx]; // storing the output
   datalo[idx3] = zvlo[tdx];
}

__global__ void cmplx2_update_addjobs
 ( double *datarehi, double *datarelo,
   double *dataimhi, double *dataimlo,
   int *in1idx, int *in2idx, int *outidx, int dim )
{
   const int bdx = blockIdx.x;           // index to the addition job
   const int tdx = threadIdx.x;          // index to the output of the job
   const int idx1 = in1idx[bdx] + tdx;
   const int idx2 = in2idx[bdx] + tdx;
   const int idx3 = outidx[bdx] + tdx;

   __shared__ double xvrehi[dd_shmemsize];
   __shared__ double xvrelo[dd_shmemsize];
   __shared__ double xvimhi[dd_shmemsize];
   __shared__ double xvimlo[dd_shmemsize];
   __shared__ double yvrehi[dd_shmemsize];
   __shared__ double yvrelo[dd_shmemsize];
   __shared__ double yvimhi[dd_shmemsize];
   __shared__ double yvimlo[dd_shmemsize];
   __shared__ double zvrehi[dd_shmemsize];
   __shared__ double zvrelo[dd_shmemsize];
   __shared__ double zvimhi[dd_shmemsize];
   __shared__ double zvimlo[dd_shmemsize];

   xvrehi[tdx] = datarehi[idx1];  // loading first input
   xvrelo[tdx] = datarelo[idx1];
   xvimhi[tdx] = dataimhi[idx1];
   xvimlo[tdx] = dataimlo[idx1];
   yvrehi[tdx] = datarehi[idx2];  // loading second input
   yvrelo[tdx] = datarelo[idx2];
   yvimhi[tdx] = dataimhi[idx2];
   yvimlo[tdx] = dataimlo[idx2];

   // zv[tdx] = xv[tdx] + yv[tdx];

   ddg_add(xvrehi[tdx],xvrelo[tdx],yvrehi[tdx],yvrelo[tdx],
           &zvrehi[tdx],&zvrelo[tdx]);
   __syncthreads();

   ddg_add(xvimhi[tdx],xvimlo[tdx],yvimhi[tdx],yvimlo[tdx],
           &zvimhi[tdx],&zvimlo[tdx]);
   __syncthreads();

   datarehi[idx3] = zvrehi[tdx]; // storing the output
   datarelo[idx3] = zvrelo[tdx];
   dataimhi[idx3] = zvimhi[tdx];
   dataimlo[idx3] = zvimlo[tdx];
}

void dbl_convoluted_data2_to_output
 ( double *datahi, double *datalo, double **outputhi, double **outputlo,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart, int vrblvl )
{
   const int deg1 = deg+1;
   int ix0,ix1,ix2;

   for(int i=0; i<=deg; i++) // output[dim][i] = data[i];
   {
      outputhi[dim][i] = datahi[i];
      outputlo[dim][i] = datalo[i];
   }
   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++) // output[i][j] = 0.0;
      {
         outputhi[i][j] = 0.0;
         outputlo[i][j] = 0.0;
      }

   for(int k=0; k<nbr; k++)
   {
      ix1 = fstart[k] + (nvr[k]-1)*deg1;
      
      if(vrblvl > 1)
         cout << "monomial " << k << " update starts at " << ix1 << endl;

      for(int i=0; i<=deg; i++) // output[dim][i] += data[ix1++];
         ddf_inc(&outputhi[dim][i],&outputlo[dim][i],
                 datahi[ix1],datalo[ix1++]);
     
      ix0 = idx[k][0];
      if(nvr[k] == 1)
      {
         ix1 = (1 + k)*deg1;
            
         for(int i=0; i<=deg; i++) // output[ix0][i] += data[ix1++];
            ddf_inc(&outputhi[ix0][i],&outputlo[ix0][i],
                    datahi[ix1],datalo[ix1++]);
      }
      else
      {                               // update first and last derivative
         ix2 = nvr[k]-3;
         if(ix2 < 0) ix2 = 0;
         ix1 = bstart[k] + ix2*deg1;

         for(int i=0; i<=deg; i++) // output[ix0][i] += data[ix1++];
            ddf_inc(&outputhi[ix0][i],&outputlo[ix0][i],
                    datahi[ix1],datalo[ix1++]);

         ix2 = nvr[k]-2;
         ix1 = fstart[k] + ix2*deg1;
         ix0 = idx[k][ix2+1];

         for(int i=0; i<=deg; i++) // output[ix0][i] += data[ix1++];
            ddf_inc(&outputhi[ix0][i],&outputlo[ix0][i],
                    datahi[ix1],datalo[ix1++]);
 
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
                  ddf_inc(&outputhi[ix0][i],&outputlo[ix0][i],
                          datahi[ix1],datalo[ix1++]);
            }
         }
      }
   }
}

void cmplx_convoluted_data2_to_output
 ( double *datarehi, double *datarelo,
   double *dataimhi, double *dataimlo,
   double **outputrehi, double **outputrelo,
   double **outputimhi, double **outputimlo,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart, int vrblvl )
{
   const int deg1 = deg+1;
   int ix0,ix1,ix2;

   for(int i=0; i<=deg; i++) // output[dim][i] = data[i];
   {
      outputrehi[dim][i] = datarehi[i];
      outputrelo[dim][i] = datarelo[i];
      outputimhi[dim][i] = dataimhi[i];
      outputimlo[dim][i] = dataimlo[i];
   }
   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++) // output[i][j] = 0.0;
      {
         outputrehi[i][j] = 0.0;
         outputrelo[i][j] = 0.0;
         outputimhi[i][j] = 0.0;
         outputimlo[i][j] = 0.0;
      }

   for(int k=0; k<nbr; k++)
   {
      ix1 = fstart[k] + (nvr[k]-1)*deg1;
      
      if(vrblvl > 1)
         cout << "monomial " << k << " update starts at " << ix1 << endl;

      for(int i=0; i<=deg; i++) // output[dim][i] += data[ix1++];
      {
         ddf_inc(&outputrehi[dim][i],&outputrelo[dim][i],
                 datarehi[ix1],datarelo[ix1++]);
         ddf_inc(&outputimhi[dim][i],&outputimlo[dim][i],
                 dataimhi[ix1],dataimlo[ix1++]);
      }
      ix0 = idx[k][0];
      if(nvr[k] == 1)
      {
         ix1 = (1 + k)*deg1;
            
         for(int i=0; i<=deg; i++) // output[ix0][i] += data[ix1++];
         {
            ddf_inc(&outputrehi[ix0][i],&outputrelo[ix0][i],
                    datarehi[ix1],datarelo[ix1++]);
            ddf_inc(&outputimhi[ix0][i],&outputimlo[ix0][i],
                    dataimhi[ix1],dataimlo[ix1++]);
         }
      }
      else
      {                               // update first and last derivative
         ix2 = nvr[k]-3;
         if(ix2 < 0) ix2 = 0;
         ix1 = bstart[k] + ix2*deg1;

         for(int i=0; i<=deg; i++) // output[ix0][i] += data[ix1++];
         {
            ddf_inc(&outputrehi[ix0][i],&outputrelo[ix0][i],
                    datarehi[ix1],datarelo[ix1++]);
            ddf_inc(&outputimhi[ix0][i],&outputimlo[ix0][i],
                    dataimhi[ix1],dataimlo[ix1++]);
         }
         ix2 = nvr[k]-2;
         ix1 = fstart[k] + ix2*deg1;
         ix0 = idx[k][ix2+1];

         for(int i=0; i<=deg; i++) // output[ix0][i] += data[ix1++];
         {
            ddf_inc(&outputrehi[ix0][i],&outputrelo[ix0][i],
                    datarehi[ix1],datarelo[ix1++]);
            ddf_inc(&outputimhi[ix0][i],&outputimlo[ix0][i],
                    dataimhi[ix1],dataimlo[ix1++]);
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
                  ddf_inc(&outputrehi[ix0][i],&outputrelo[ix0][i],
                          datarehi[ix1],datarelo[ix1++]);
                  ddf_inc(&outputimhi[ix0][i],&outputimlo[ix0][i],
                          dataimhi[ix1],dataimlo[ix1++]);
               }
            }
         }
      }
   }
}

void dbl_added_data2_to_output
 ( double *datahi, double *datalo, double **outputhi, double **outputlo,
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
      outputhi[dim][i] = datahi[ix];
      outputlo[dim][i] = datalo[ix++];
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
            outputhi[0][i] = 0.0;
            outputlo[0][i] = 0.0;
         }
      }
      else
      {
         int cffidx = (1 + difidx)*deg1;

         if(vrblvl > 1)
            cout << "updating derivative with coefficient ..." << endl;

         for(int i=0; i<=deg; i++)
         {
            outputhi[0][i] = datahi[cffidx];
            outputlo[0][i] = datalo[cffidx++];
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
         outputhi[0][i] = datahi[ix];
         outputlo[0][i] = datalo[ix++];
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
               outputhi[k][i] = 0.0;
               outputlo[k][i] = 0.0;
            }
         }
         else
         {
            int cffidx = (1 + difidx)*deg1;

            if(vrblvl > 1)
               cout << "updating derivative with coefficient ..." << endl;

            for(int i=0; i<=deg; i++)
            {
               outputhi[k][i] = datahi[cffidx];
               outputlo[k][i] = datalo[cffidx++];
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
               outputhi[k][i] = datahi[ix];
               outputlo[k][i] = datalo[ix++];
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
               outputhi[k][i] = datahi[ix];
               outputlo[k][i] = datalo[ix++];
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
               outputhi[k][i] = datahi[ix];
               outputlo[k][i] = datalo[ix++];
            }
         }
      }
   }
}

void cmplx_added_data2_to_output
 ( double *datarehi, double *datarelo,
   double *dataimhi, double *dataimlo,
   double **outputrehi, double **outputrelo,
   double **outputimhi, double **outputimlo,
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
      outputrehi[dim][i] = datarehi[ix];
      outputrelo[dim][i] = datarelo[ix];
      outputimhi[dim][i] = dataimhi[ix];
      outputimlo[dim][i] = dataimlo[ix++];
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
            outputrehi[0][i] = 0.0; outputrelo[0][i] = 0.0;
            outputimhi[0][i] = 0.0; outputimlo[0][i] = 0.0;
         }
      }
      else
      {
         int cffidx = (1 + difidx)*deg1;

         if(vrblvl > 1)
            cout << "updating derivative 0 with coefficient "
                 << difidx << endl;

         for(int i=0; i<=deg; i++)
         {
            outputrehi[0][i] = datarehi[cffidx];
            outputrelo[0][i] = datarelo[cffidx];
            outputimhi[0][i] = dataimhi[cffidx];
            outputimlo[0][i] = dataimlo[cffidx++];
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
         outputrehi[0][i] = datarehi[ix];
         outputrelo[0][i] = datarelo[ix];
         outputimhi[0][i] = dataimhi[ix];
         outputimlo[0][i] = dataimlo[ix++];
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
               outputrehi[k][i] = 0.0; outputrelo[k][i] = 0.0;
               outputimhi[k][i] = 0.0; outputimlo[k][i] = 0.0;
            }
         }
         else
         {
            int cffidx = (1 + difidx)*deg1;

            if(vrblvl > 1)
               cout << "updating derivative " << k
                    << " with coefficient " << difidx << endl;

            for(int i=0; i<=deg; i++)
            {
               outputrehi[k][i] = datarehi[cffidx];
               outputrelo[k][i] = datarelo[cffidx];
               outputimhi[k][i] = dataimhi[cffidx];
               outputimlo[k][i] = dataimlo[cffidx++];
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
               outputrehi[k][i] = datarehi[ix];
               outputrelo[k][i] = datarelo[ix];
               outputimhi[k][i] = dataimhi[ix];
               outputimlo[k][i] = dataimlo[ix++];
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
               outputrehi[k][i] = datarehi[ix];
               outputrelo[k][i] = datarelo[ix];
               outputimhi[k][i] = dataimhi[ix];
               outputimlo[k][i] = dataimlo[ix++];
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
               outputrehi[k][i] = datarehi[ix];
               outputrelo[k][i] = datarelo[ix];
               outputimhi[k][i] = dataimhi[ix];
               outputimlo[k][i] = dataimlo[ix++];
            }
         }
      }
   }
}

void cmplx_added_data2vectorized_to_output
 ( double *datarihi, double *datarilo,
   double **outputrehi, double **outputrelo,
   double **outputimhi, double **outputimlo,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart,
   int totcff, int offsetri, ComplexAdditionJobs jobs, int vrblvl )
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
         outputrehi[i][j] = 0.0; outputrelo[i][j] = 0.0;
         outputimhi[i][j] = 0.0; outputimlo[i][j] = 0.0; 
      }
 */
   ix1re = fstart[lastmon] + lastidx*deg1;
   ix2im = fstart[lastmon] + lastidx*deg1 + totcffoffset;

   if(vrblvl > 1)
      cout << "Updating value starting at " << ix1re << " in data." << endl;

   for(int i=0; i<=deg; i++) // output[dim][i] = data[ix++];
   {
      outputrehi[dim][i] = datarihi[ix1re];
      outputrelo[dim][i] = datarilo[ix1re++];
      outputimhi[dim][i] = datarihi[ix2im];
      outputimlo[dim][i] = datarilo[ix2im++];
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
            outputrehi[0][i] = 0.0; outputrelo[0][i] = 0.0;
            outputimhi[0][i] = 0.0; outputimlo[0][i] = 0.0;
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
            outputrehi[0][i] = datarihi[ix1re];
            outputrelo[0][i] = datarilo[ix1re++];
            outputimhi[0][i] = datarihi[ix2im];
            outputimlo[0][i] = datarilo[ix2im++];
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
         outputrehi[0][i] = datarihi[ix1re];
         outputrelo[0][i] = datarilo[ix1re++];
         outputimhi[0][i] = datarihi[ix2im];
         outputimlo[0][i] = datarilo[ix2im++];
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
               outputrehi[k][i] = 0.0; outputrelo[k][i] = 0.0;
               outputimhi[k][i] = 0.0; outputimlo[k][i] = 0.0;
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
               outputrehi[k][i] = datarihi[ix1re];
               outputrelo[k][i] = datarilo[ix1re++];
               outputimhi[k][i] = datarihi[ix2im];
               outputimlo[k][i] = datarilo[ix2im++];
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
               outputrehi[k][i] = datarihi[ix1re];
               outputrelo[k][i] = datarilo[ix1re++];
               outputimhi[k][i] = datarihi[ix2im];
               outputimlo[k][i] = datarilo[ix2im++];
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
               outputrehi[k][i] = datarihi[ix1re];
               outputrelo[k][i] = datarilo[ix1re++];
               outputimhi[k][i] = datarihi[ix2im];
               outputimlo[k][i] = datarilo[ix2im++];
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
               outputrehi[k][i] = datarihi[ix1re];
               outputrelo[k][i] = datarilo[ix1re++];
               outputimhi[k][i] = datarihi[ix2im];
               outputimlo[k][i] = datarilo[ix2im++];
            }
         }
      }
   }
}

void dbl2_data_setup
 ( int dim, int nbr, int deg, int totcff,
   double *datahi, double *datalo,
   double *csthi, double *cstlo, double **cffhi, double **cfflo,
   double **inputhi, double **inputlo )
{
   const int deg1 = deg+1;
   int ix = 0;

   for(int i=0; i<deg1; i++)
   {
      datahi[ix]   = csthi[i];
      datalo[ix++] = cstlo[i];
   }
   for(int i=0; i<nbr; i++)
      for(int j=0; j<deg1; j++)
      {
         datahi[ix]   = cffhi[i][j];
         datalo[ix++] = cfflo[i][j];
      }
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++)
      {
         datahi[ix]   = inputhi[i][j];
         datalo[ix++] = inputlo[i][j];
      }

   for(int i=ix; i<totcff; i++)
   {
      datahi[i] = 0.0; datalo[i] = 0.0;
   }
}

void cmplx2_data_setup
 ( int dim, int nbr, int deg, int totcff,
   double *datarehi, double *datarelo, double *dataimhi, double *dataimlo,
   double *cstrehi, double *cstrelo, double *cstimhi, double *cstimlo,
   double **cffrehi, double **cffrelo, double **cffimhi, double **cffimlo,
   double **inputrehi, double **inputrelo,
   double **inputimhi, double **inputimlo )
{
   const int deg1 = deg+1;
   int ix = 0;

   for(int i=0; i<deg1; i++)
   {
      datarehi[ix]   = cstrehi[i];
      datarelo[ix]   = cstrelo[i];
      dataimhi[ix]   = cstimhi[i];
      dataimlo[ix++] = cstimlo[i];
   }
   for(int i=0; i<nbr; i++)
      for(int j=0; j<deg1; j++)
      {
         datarehi[ix]   = cffrehi[i][j];
         datarelo[ix]   = cffrelo[i][j];
         dataimhi[ix]   = cffimhi[i][j];
         dataimlo[ix++] = cffimlo[i][j];
      }
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++)
      {
         datarehi[ix]   = inputrehi[i][j];
         datarelo[ix]   = inputrelo[i][j];
         dataimhi[ix]   = inputimhi[i][j];
         dataimlo[ix++] = inputimlo[i][j];
      }

   for(int i=ix; i<totcff; i++)
   {
      datarehi[i] = 0.0; datarelo[i] = 0.0;
      dataimhi[i] = 0.0; dataimlo[i] = 0.0;
   }
}

void cmplx2vectorized_data_setup
 ( int dim, int nbr, int deg, int totcff, int offsetri,
   double *datarihi, double *datarilo,
   double *cstrehi, double *cstrelo, double *cstimhi, double *cstimlo,
   double **cffrehi, double **cffrelo, double **cffimhi, double **cffimlo,
   double **inputrehi, double **inputrelo,
   double **inputimhi, double **inputimlo )
{
   const int deg1 = deg+1;

   int ix1 = 0;
   int ix2 = totcff + offsetri;

   for(int i=0; i<deg1; i++)
   {
      datarihi[ix1]   = cstrehi[i];
      datarilo[ix1++] = cstrelo[i];
      datarihi[ix2]   = cstimhi[i];
      datarilo[ix2++] = cstimlo[i];
   }
   for(int i=0; i<nbr; i++)
      for(int j=0; j<deg1; j++)
      {
         datarihi[ix1]   = cffrehi[i][j];
         datarilo[ix1++] = cffrelo[i][j];
         datarihi[ix2]   = cffimhi[i][j];
         datarilo[ix2++] = cffimlo[i][j];
      }
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++)
      {
         datarihi[ix1]   = inputrehi[i][j];
         datarilo[ix1++] = inputrelo[i][j];
         datarihi[ix2]   = inputimhi[i][j];
         datarilo[ix2++] = inputimlo[i][j];
      }

   for(int i=0; i<2*offsetri; i++)
   {
      datarihi[ix1]   = 0.0;
      datarilo[ix1++] = 0.0;
      datarihi[ix2]   = 0.0;
      datarilo[ix2++] = 0.0;
   }
}

void dbl2_convolution_jobs
 ( int dim, int nbr, int deg, int *nvr, ConvolutionJobs cnvjobs,
   int *fstart, int *bstart, int *cstart,
   double *datahi, double *datalo, double *cnvlapms, int vrblvl )
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
         dbl2_padded_convjobs<<<jobnbr,deg1>>>
            (datahi,datalo,in1ix_d,in2ix_d,outix_d,deg1);
         cudaEventRecord(stop);
         cudaEventSynchronize(stop);
         cudaEventElapsedTime(&milliseconds,start,stop);
         *cnvlapms += milliseconds;
         cudaFree(in1ix_d); cudaFree(in2ix_d); cudaFree(outix_d);
      }
      free(in1ix_h); free(in2ix_h); free(outix_h);
   }
}

void cmplx2_convolution_jobs
 ( int dim, int nbr, int deg, int *nvr, ConvolutionJobs cnvjobs,
   int *fstart, int *bstart, int *cstart,
   double *datarehi, double *datarelo, double *dataimhi, double *dataimlo,
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
         cmplx2_padded_convjobs<<<jobnbr,deg1>>>
            (datarehi,datarelo,dataimhi,dataimlo,
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

void cmplx2vectorized_convolution_jobs
 ( int dim, int nbr, int deg, int *nvr, int totcff, int offsetri,
   ComplexConvolutionJobs cnvjobs, ComplexIncrementJobs incjobs,
   int *fstart, int *bstart, int *cstart,
   double *datarihi, double *datarilo, double *cnvlapms, int vrblvl )
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
         dbl2_padded_convjobs<<<jobnbr,deg1>>>
            (datarihi,datarilo,in1ix_d,in2ix_d,outix_d,deg1);
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

      GPU_cmplx2vectorized_flipsigns
         (deg,nbrflips,rebidx,datarihi,datarilo,&fliplapms,vrblvl);
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
         dbl2_increment_jobs<<<jobnbr,deg1>>>
            (datarihi,datarilo,in1ix_d,in2ix_d,outix_d,deg1);
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

void dbl2_addition_jobs
 ( int dim, int nbr, int deg, int *nvr, AdditionJobs addjobs,
   int *fstart, int *bstart, int *cstart,
   double *datahi, double *datalo, double *addlapms, int vrblvl )
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
         dbl2_update_addjobs<<<jobnbr,deg1>>>
            (datahi,datalo,in1ix_d,in2ix_d,outix_d,deg1);
         cudaEventRecord(stop);
         cudaEventSynchronize(stop);
         cudaEventElapsedTime(&milliseconds,start,stop);
         *addlapms += milliseconds;
         cudaFree(in1ix_d); cudaFree(in2ix_d); cudaFree(outix_d);
      }
      free(in1ix_h); free(in2ix_h); free(outix_h);
   }
}

void cmplx2_addition_jobs
 ( int dim, int nbr, int deg, int *nvr, AdditionJobs addjobs,
   int *fstart, int *bstart, int *cstart,
   double *datarehi, double *datarelo, double *dataimhi, double *dataimlo,
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
         cmplx2_update_addjobs<<<jobnbr,deg1>>>
            (datarehi,datarelo,dataimhi,dataimlo,
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

void cmplx2vectorized_addition_jobs
 ( int dim, int nbr, int deg, int *nvr, int totcff, int offsetri,
   ComplexAdditionJobs addjobs,
   int *fstart, int *bstart, int *cstart,
   double *datarihi, double *datarilo, double *addlapms, int vrblvl )
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
         dbl2_update_addjobs<<<jobnbr,deg1>>>
            (datarihi,datarilo,in1ix_d,in2ix_d,outix_d,deg1);
         cudaEventRecord(stop);
         cudaEventSynchronize(stop);
         cudaEventElapsedTime(&milliseconds,start,stop);
         *addlapms += milliseconds;
         cudaFree(in1ix_d); cudaFree(in2ix_d); cudaFree(outix_d);
      }
      free(in1ix_h); free(in2ix_h); free(outix_h);
   }
}

void GPU_dbl2_poly_evaldiff
 ( int BS, int dim, int nbr, int deg, int *nvr, int **idx,
   double *csthi, double *cstlo, double **cffhi, double **cfflo,
   double **inputhi, double **inputlo,
   double **outputhi, double **outputlo,
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

   double *datahi_h = new double[totalcff];        // data on host
   double *datalo_h = new double[totalcff];

   dbl2_data_setup
      (dim,nbr,deg,totalcff,datahi_h,datalo_h,
       csthi,cstlo,cffhi,cfflo,inputhi,inputlo);

   double *datahi_d;                               // device data
   double *datalo_d;
   const size_t szdata = totalcff*sizeof(double);
   cudaMalloc((void**)&datahi_d,szdata);
   cudaMalloc((void**)&datalo_d,szdata);
   cudaMemcpy(datahi_d,datahi_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datalo_d,datalo_h,szdata,cudaMemcpyHostToDevice);

   struct timeval begintime,endtime; // wall clock time of computations

   gettimeofday(&begintime,0);

   dbl2_convolution_jobs
      (dim,nbr,deg,nvr,cnvjobs,fstart,bstart,cstart,
       datahi_d,datalo_d,cnvlapms,vrblvl);

   dbl2_addition_jobs
      (dim,nbr,deg,nvr,addjobs,fstart,bstart,cstart,
       datahi_d,datalo_d,addlapms,vrblvl);

   gettimeofday(&endtime,0);
   cudaMemcpy(datahi_h,datahi_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datalo_h,datalo_d,szdata,cudaMemcpyDeviceToHost);
   *elapsedms = *cnvlapms + *addlapms;
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   *walltimesec = seconds + microseconds*1.0e-6;

   // dbl_convoluted_data2_to_output
   //    (datahi_h,datalo_h,outputhi,outputlo,
   //     dim,nbr,deg,nvr,idx,fstart,bstart,cstart,vrblvl);
   dbl_added_data2_to_output
      (datahi_h,datalo_h,outputhi,outputlo,
       dim,nbr,deg,nvr,idx,fstart,bstart,cstart,addjobs,vrblvl);

   if(vrblvl > 0)
      write_GPU_timings(*cnvlapms,*addlapms,*elapsedms,*walltimesec);

   cudaFree(datahi_d); cudaFree(datalo_d);

   free(datahi_h); free(datalo_h);

   free(fstart); free(bstart); free(cstart);
   free(fsums); free(bsums); free(csums);
}

void GPU_cmplx2_poly_evaldiff
 ( int BS, int dim, int nbr, int deg, int *nvr, int **idx,
   double *cstrehi, double *cstrelo,
   double *cstimhi, double *cstimlo,
   double **cffrehi, double **cffrelo,
   double **cffimhi, double **cffimlo,
   double **inputrehi, double **inputrelo,
   double **inputimhi, double **inputimlo,
   double **outputrehi, double **outputrelo,
   double **outputimhi, double **outputimlo,
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

   double *datarehi_h = new double[totalcff];      // data on host
   double *datarelo_h = new double[totalcff];
   double *dataimhi_h = new double[totalcff]; 
   double *dataimlo_h = new double[totalcff];

   cmplx2_data_setup
      (dim,nbr,deg,totalcff,
       datarehi_h,datarelo_h,dataimhi_h,dataimlo_h,
       cstrehi,cstrelo,cstimhi,cstimlo,
       cffrehi,cffrelo,cffimhi,cffimlo,
       inputrehi,inputrelo,inputimhi,inputimlo);

   double *datarehi_d;                               // device data
   double *datarelo_d;
   double *dataimhi_d;
   double *dataimlo_d;
   const size_t szdata = totalcff*sizeof(double);
   cudaMalloc((void**)&datarehi_d,szdata);
   cudaMalloc((void**)&datarelo_d,szdata);
   cudaMalloc((void**)&dataimhi_d,szdata);
   cudaMalloc((void**)&dataimlo_d,szdata);
   cudaMemcpy(datarehi_d,datarehi_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datarelo_d,datarelo_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(dataimhi_d,dataimhi_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(dataimlo_d,dataimlo_h,szdata,cudaMemcpyHostToDevice);

   struct timeval begintime,endtime; // wall clock time of computations

   gettimeofday(&begintime,0);

   cmplx2_convolution_jobs
      (dim,nbr,deg,nvr,cnvjobs,fstart,bstart,cstart,
       datarehi_d,datarelo_d,dataimhi_d,dataimlo_d,cnvlapms,vrblvl);

   cmplx2_addition_jobs
      (dim,nbr,deg,nvr,addjobs,fstart,bstart,cstart,
       datarehi_d,datarelo_d,dataimhi_d,dataimlo_d,
       addlapms,vrblvl);

   gettimeofday(&endtime,0);
   cudaMemcpy(datarehi_h,datarehi_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datarelo_h,datarelo_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(dataimhi_h,dataimhi_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(dataimlo_h,dataimlo_d,szdata,cudaMemcpyDeviceToHost);
   *elapsedms = *cnvlapms + *addlapms;
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   *walltimesec = seconds + microseconds*1.0e-6;

   // cmplx_convoluted_data2_to_output
   //    (datarehi_h,datarelo_h,dataimhi_h,dataimlo_h,
   //     outputrehi,outputrelo,outputimhi,outputimlo,
   //     dim,nbr,deg,nvr,idx,fstart,bstart,cstart,vrblvl);
   cmplx_added_data2_to_output
      (datarehi_h,datarelo_h,dataimhi_h,dataimlo_h,
       outputrehi,outputrelo,outputimhi,outputimlo,
       dim,nbr,deg,nvr,idx,fstart,bstart,cstart,addjobs,vrblvl);

   if(vrblvl > 0)
      write_GPU_timings(*cnvlapms,*addlapms,*elapsedms,*walltimesec);

   cudaFree(datarehi_d); cudaFree(datarelo_d);
   cudaFree(dataimhi_d); cudaFree(dataimlo_d);

   free(datarehi_h); free(datarelo_h);
   free(dataimhi_h); free(dataimlo_h);

   free(fstart); free(bstart); free(cstart);
   free(fsums); free(bsums); free(csums);
}

void GPU_cmplx2vectorized_poly_evaldiff
 ( int BS, int dim, int nbr, int deg, int *nvr, int **idx,
   double *cstrehi, double *cstrelo, double *cstimhi, double *cstimlo,
   double **cffrehi, double **cffrelo, double **cffimhi, double **cffimlo,
   double **inputrehi, double **inputrelo,
   double **inputimhi, double **inputimlo,
   double **outputrehi, double **outputrelo,
   double **outputimhi, double **outputimlo,
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
   double *datarihi_h = new double[cmplxtotcff];      // data on host
   double *datarilo_h = new double[cmplxtotcff];

   cmplx2vectorized_data_setup
      (dim,nbr,deg,totalcff,offsetri,datarihi_h,datarilo_h,
       cstrehi,cstrelo,cstimhi,cstimlo, cffrehi,cffrelo,cffimhi,cffimlo,
       inputrehi,inputrelo,inputimhi,inputimlo);

   double *datarihi_d;                               // device data
   double *datarilo_d;
   const size_t szdata = cmplxtotcff*sizeof(double);
   cudaMalloc((void**)&datarihi_d,szdata);
   cudaMalloc((void**)&datarilo_d,szdata);
   cudaMemcpy(datarihi_d,datarihi_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datarilo_d,datarilo_h,szdata,cudaMemcpyHostToDevice);

   struct timeval begintime,endtime; // wall clock time of computations

   gettimeofday(&begintime,0);

   cmplx2vectorized_convolution_jobs
      (dim,nbr,deg,nvr,totalcff,offsetri,cnvjobs,incjobs,fstart,bstart,cstart,
       datarihi_d,datarilo_d,cnvlapms,vrblvl);

   cmplx2vectorized_addition_jobs
      (dim,nbr,deg,nvr,totalcff,offsetri,addjobs,fstart,bstart,cstart,
       datarihi_d,datarilo_d,addlapms,vrblvl);

   gettimeofday(&endtime,0);
   cudaMemcpy(datarihi_h,datarihi_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datarilo_h,datarilo_d,szdata,cudaMemcpyDeviceToHost);
   *elapsedms = *cnvlapms + *addlapms;
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   *walltimesec = seconds + microseconds*1.0e-6;

   cmplx_added_data2vectorized_to_output
      (datarihi_h,datarilo_h,outputrehi,outputrelo,outputimhi,outputimlo,
       dim,nbr,deg,nvr,idx,fstart,bstart,cstart,
       totalcff,offsetri,addjobs,vrblvl);

   if(vrblvl > 0)
      write_GPU_timings(*cnvlapms,*addlapms,*elapsedms,*walltimesec);

   cudaFree(datarihi_d); cudaFree(datarilo_d);

   free(datarihi_h); free(datarilo_h);

   free(fstart); free(bstart); free(cstart);
   free(fsums); free(bsums); free(csums);
}
