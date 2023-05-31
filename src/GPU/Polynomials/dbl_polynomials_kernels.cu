// The file dbl_polynomials_kernels.cu defines the kernels with prototypes
// in dbl_polynomials_kernels.h.

#include <iostream>
#include <iomanip>
#ifdef winwalltime
#include "gettimeofday4win.h"
#else
#include <sys/time.h>
#endif
#include "job_coordinates.h"
#include "write_gpu_timings.h"
#include "dbl_polynomials_kernels.h"

// The constant d_shmemsize is the bound on the shared memory size.

#define d_shmemsize 256

using namespace std;

__global__ void dbl_padded_convjobs
 ( double *data, int *in1idx, int *in2idx, int *outidx, int dim )
{
   const int bdx = blockIdx.x;           // index to the convolution job
   const int tdx = threadIdx.x;          // index to the output of the job
   const int idx1 = in1idx[bdx] + tdx;
   const int idx2 = in2idx[bdx] + tdx;
   const int idx3 = outidx[bdx] + tdx;

   __shared__ double xv[d_shmemsize];
   __shared__ double yv[2*d_shmemsize];
   __shared__ double zv[d_shmemsize];

   int ydx = dim + tdx;

   xv[tdx] = data[idx1];  // loading first input
   yv[tdx] = 0.0;         // padded with zeros
   yv[ydx] = data[idx2];  // loading second input

   zv[tdx] = xv[0]*yv[ydx];

   for(int i=1; i<dim; i++)
   {
      ydx = dim + tdx - i;
      zv[tdx] = zv[tdx] + xv[i]*yv[ydx];
   }
   data[idx3] = zv[tdx]; // storing the output
}

__global__ void dbl_increment_jobs
 ( double *data, int *in1idx, int *in2idx, int *outidx, int dim )
{
   const int bdx = blockIdx.x;           // index to the increment job
   const int tdx = threadIdx.x;          // index to the output of the job
   const int idx1 = in1idx[bdx] + tdx;
   const int idx2 = in2idx[bdx] + tdx;
   const int idx3 = outidx[bdx] + tdx;

   __shared__ double yv[d_shmemsize];
   __shared__ double zv[d_shmemsize];

   zv[tdx] = data[idx1];  // loading first input
   yv[tdx] = data[idx2];  // loading second input

   __syncthreads();

   zv[tdx] += yv[tdx];

   __syncthreads();

   data[idx3] = zv[tdx]; // storing the output
}

__global__ void cmplx_padded_convjobs
 ( double *datare, double *dataim,
   int *in1idx, int *in2idx, int *outidx, int dim )
{
   const int bdx = blockIdx.x;           // index to the convolution job
   const int tdx = threadIdx.x;          // index to the output of the job
   const int idx1 = in1idx[bdx] + tdx;
   const int idx2 = in2idx[bdx] + tdx;
   const int idx3 = outidx[bdx] + tdx;

   __shared__ double xvre[d_shmemsize];
   __shared__ double xvim[d_shmemsize];
   __shared__ double yvre[2*d_shmemsize];
   __shared__ double yvim[2*d_shmemsize];
   __shared__ double zvre[d_shmemsize];
   __shared__ double zvim[d_shmemsize];

   int ydx = dim + tdx;

   xvre[tdx] = datare[idx1];  // loading first input
   xvim[tdx] = dataim[idx1];
   yvre[tdx] = 0.0;           // padded with zeros
   yvim[tdx] = 0.0; 
   yvre[ydx] = datare[idx2];  // loading second input
   yvim[ydx] = dataim[idx2];

   zvre[tdx] = xvre[0]*yvre[ydx] - xvim[0]*yvim[ydx];
   zvim[tdx] = xvre[0]*yvim[ydx] + xvim[0]*yvre[ydx];

   for(int i=1; i<dim; i++)
   {
      ydx = dim + tdx - i;
      zvre[tdx] += xvre[i]*yvre[ydx] - xvim[i]*yvim[ydx];
      zvim[tdx] += xvre[i]*yvim[ydx] + xvim[i]*yvre[ydx];
   }
   datare[idx3] = zvre[tdx]; // storing the output
   dataim[idx3] = zvim[tdx];
}

__global__ void cmplxvectorized_flipsigns
 ( double *datari, int *flpidx, int dim )
{
   const int bdx = blockIdx.x;    // index to the series to flip
   const int tdx = threadIdx.x;
   const int idx = flpidx[bdx] + tdx; // which number to flip

   double x; // register to load data from global memory

   x = datari[idx]; x = -x; datari[idx] = x;
}

void GPU_cmplxvectorized_flipsigns
 ( int deg, int nbrflips, int *flipidx,
   double *datari, double *elapsedms, int vrblvl )
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
   cmplxvectorized_flipsigns<<<nbrflips,deg1>>>(datari,flipidx_d,deg1);
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

__global__ void dbl_update_addjobs
 ( double *data, int *in1idx, int *in2idx, int *outidx, int dim )
{
   const int bdx = blockIdx.x;           // index to the addition job
   const int tdx = threadIdx.x;          // index to the output of the job
   const int idx1 = in1idx[bdx] + tdx;
   const int idx2 = in2idx[bdx] + tdx;
   const int idx3 = outidx[bdx] + tdx;

   __shared__ double xv[d_shmemsize];
   __shared__ double yv[d_shmemsize];
   __shared__ double zv[d_shmemsize];

   xv[tdx] = data[idx1];  // loading first input
   yv[tdx] = data[idx2];  // loading second input

   zv[tdx] = xv[tdx] + yv[tdx];

   data[idx3] = zv[tdx]; // storing the output
}

__global__ void cmplx_update_addjobs
 ( double *datare, double *dataim,
   int *in1idx, int *in2idx, int *outidx, int dim )
{
   const int bdx = blockIdx.x;           // index to the addition job
   const int tdx = threadIdx.x;          // index to the output of the job
   const int idx1 = in1idx[bdx] + tdx;
   const int idx2 = in2idx[bdx] + tdx;
   const int idx3 = outidx[bdx] + tdx;

   __shared__ double xvre[d_shmemsize];
   __shared__ double xvim[d_shmemsize];
   __shared__ double yvre[d_shmemsize];
   __shared__ double yvim[d_shmemsize];
   __shared__ double zvre[d_shmemsize];
   __shared__ double zvim[d_shmemsize];

   xvre[tdx] = datare[idx1];  // loading first input
   xvim[tdx] = dataim[idx1];
   yvre[tdx] = datare[idx2];  // loading second input
   yvim[tdx] = dataim[idx2];

   zvre[tdx] = xvre[tdx] + yvre[tdx]; // adding real parts
   zvim[tdx] = xvim[tdx] + yvim[tdx]; // adding imaginary parts

   datare[idx3] = zvre[tdx]; // storing the output
   dataim[idx3] = zvim[tdx];
}

void dbl_convoluted_data_to_output
 ( double *data, double **output, int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart, int vrblvl )
{
   const int deg1 = deg+1;
   int ix0,ix1,ix2;

   for(int i=0; i<=deg; i++) output[dim][i] = data[i];
   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++) output[i][j] = 0.0;

   for(int k=0; k<nbr; k++)
   {
      ix1 = fstart[k] + (nvr[k]-1)*deg1;
      
      if(vrblvl > 1)
         cout << "monomial " << k << " update starts at " << ix1 << endl;

      for(int i=0; i<=deg; i++) output[dim][i] += data[ix1++];

      ix0 = idx[k][0];
      if(nvr[k] == 1)
      {
         ix1 = (1 + k)*deg1;
            
         for(int i=0; i<=deg; i++) output[ix0][i] += data[ix1++];
      }
      else if(nvr[k] > 1)
      {                               // update first and last derivative
         ix2 = nvr[k]-3;
         if(ix2 < 0) ix2 = 0;
         ix1 = bstart[k] + ix2*deg1;

         for(int i=0; i<=deg; i++) output[ix0][i] += data[ix1++];

         ix2 = nvr[k]-2;
         ix1 = fstart[k] + ix2*deg1;
         ix0 = idx[k][ix2+1];

         for(int i=0; i<=deg; i++) output[ix0][i] += data[ix1++];

         if(nvr[k] > 2)                   // update all other derivatives
         {
            for(int j=1; j<nvr[k]-1; j++)
            {
               ix0 = idx[k][j];            // j-th variable in monomial k
               ix1 = cstart[k] + (j-1)*deg1;

               if(vrblvl > 1)
                  cout << "monomial " << k << " derivative " << ix0
                       << " update starts at " << ix1 << endl;

               for(int i=0; i<=deg; i++) output[ix0][i] += data[ix1++];
            }
         }
      }
   }
}

void cmplx_convoluted_data_to_output
 ( double *datare, double *dataim, double **outputre, double **outputim,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart, int vrblvl )
{
   const int deg1 = deg+1;
   int ix0,ix1,ix2;

   for(int i=0; i<=deg; i++)
   {
      outputre[dim][i] = datare[i]; outputim[dim][i] = dataim[i];
   }
   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++)
      {
         outputre[i][j] = 0.0; outputim[i][j] = 0.0;
      }

   for(int k=0; k<nbr; k++)
   {
      ix1 = fstart[k] + (nvr[k]-1)*deg1;
      
      if(vrblvl > 1)
         cout << "monomial " << k << " update starts at " << ix1 << endl;

      for(int i=0; i<=deg; i++)
      {
         outputre[dim][i] += datare[ix1];
         outputim[dim][i] += dataim[ix1++];
      }
      ix0 = idx[k][0];
      if(nvr[k] == 1)
      {
         ix1 = (1 + k)*deg1;
            
         for(int i=0; i<=deg; i++)
         {
            outputre[ix0][i] += datare[ix1];
            outputim[ix0][i] += dataim[ix1++];
         }
      }
      else if(nvr[k] > 1)
      {                               // update first and last derivative
         ix2 = nvr[k]-3;
         if(ix2 < 0) ix2 = 0;
         ix1 = bstart[k] + ix2*deg1;

         for(int i=0; i<=deg; i++)
         {
            outputre[ix0][i] += datare[ix1];
            outputim[ix0][i] += dataim[ix1++];
         }
         ix2 = nvr[k]-2;
         ix1 = fstart[k] + ix2*deg1;
         ix0 = idx[k][ix2+1];

         for(int i=0; i<=deg; i++)
         {
            outputre[ix0][i] += datare[ix1];
            outputim[ix0][i] += dataim[ix1++];
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

               for(int i=0; i<=deg; i++)
               {
                  outputre[ix0][i] += datare[ix1];
                  outputim[ix0][i] += dataim[ix1++];
               }
            }
         }
      }
   }
}

void dbl_added_data_to_output
 ( double *data, double **output, int dim, int nbr, int deg, int *nvr,
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

   for(int i=0; i<=deg; i++) output[dim][i] = data[ix++];

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
         for(int i=0; i<=deg; i++) output[0][i] = 0.0;
      }
      else
      {
         int cffidx = (1 + difidx)*deg1;

         if(vrblvl > 1)
            cout << "updating derivative with coefficient ..." << endl;

         for(int i=0; i<=deg; i++) output[0][i] = data[cffidx++];
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

      for(int i=0; i<=deg; i++) output[0][i] = data[ix++];
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
            for(int i=0; i<=deg; i++) output[k][i] = 0.0;
         }
         else
         {
            int cffidx = (1 + difidx)*deg1;

            if(vrblvl > 1)
               cout << "updating derivative with coefficient ..." << endl;

            for(int i=0; i<=deg; i++) output[k][i] = data[cffidx++];
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

            for(int i=0; i<=deg; i++) output[k][i] = data[ix++];
         }
         else if(idx[ix0][nvr[ix0]-1] == k) // k is last variable
         {
            int ix2 = nvr[ix0]-2;
   
            if(vrblvl > 1)
               cout << "Updating derivative " << k 
                    << " at " << ix << " in data." << endl;

            ix = fstart[ix0] + ix2*deg1;

            for(int i=0; i<=deg; i++) output[k][i] = data[ix++];
         }
         else // derivative is in some cross product
         {
            int ix2 = jobs.position(nvr[ix0],idx[ix0],k) - 1;
   
            if(vrblvl > 1)
               cout << "Updating derivative " << k 
                    << " at " << ix << " in data." << endl;

            ix = cstart[ix0] + ix2*deg1;

            for(int i=0; i<=deg; i++) output[k][i] = data[ix++];
         }
      }
   }
}

void cmplx_added_data_to_output
 ( double *datare, double *dataim, double **outputre, double **outputim,
   int dim, int nbr, int deg, int *nvr, int **idx,
   int *fstart, int *bstart, int *cstart, AdditionJobs jobs, int vrblvl )
{
   const int deg1 = deg + 1;
   const int lastmon = nbr-1;
   const int lastidx = nvr[lastmon]-1;
   int ix;

   ix = fstart[lastmon] + lastidx*deg1;

   if(vrblvl > 1)
      cout << "Updating value starting at " << ix << " in data." << endl;

   for(int i=0; i<=deg; i++)
   {
      outputre[dim][i] = datare[ix]; outputim[dim][i] = dataim[ix++];
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
         for(int i=0; i<=deg; i++)
         {
            outputre[0][i] = 0.0; outputim[0][i] = 0.0;
         }
      }
      else
      {
         int cffidx = (1 + difidx)*deg1;

         if(vrblvl > 1)
            cout << "updating derivative with coefficient ..." << endl;

         for(int i=0; i<=deg; i++)
         {
            outputre[0][i] = datare[cffidx];
            outputim[0][i] = dataim[cffidx++];
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

      for(int i=0; i<=deg; i++)
      {
         outputre[0][i] = datare[ix]; outputim[0][i] = dataim[ix++];
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
            for(int i=0; i<=deg; i++)
            {
               outputre[k][i] = 0.0; outputim[k][i] = 0.0;
            }
         }
         else
         {
            int cffidx = (1 + difidx)*deg1;

            if(vrblvl > 1)
               cout << "updating derivative with coefficient ..." << endl;

            for(int i=0; i<=deg; i++)
            {
               outputre[k][i] = datare[cffidx];
               outputim[k][i] = dataim[cffidx++];
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

            for(int i=0; i<=deg; i++)
            {
               outputre[k][i] = datare[ix]; outputim[k][i] = dataim[ix++];
            }
         }
         else if(idx[ix0][nvr[ix0]-1] == k) // k is last variable
         {
            int ix2 = nvr[ix0]-2;
   
            if(vrblvl > 1)
               cout << "Updating derivative " << k 
                    << " at " << ix << " in data." << endl;

            ix = fstart[ix0] + ix2*deg1;

            for(int i=0; i<=deg; i++)
            {
               outputre[k][i] = datare[ix]; outputim[k][i] = dataim[ix++];
            }
         }
         else // derivative is in some cross product
         {
            int ix2 = jobs.position(nvr[ix0],idx[ix0],k) - 1;
   
            if(vrblvl > 1)
               cout << "Updating derivative " << k 
                    << " at " << ix << " in data." << endl;

            ix = cstart[ix0] + ix2*deg1;

            for(int i=0; i<=deg; i++)
            {
               outputre[k][i] = datare[ix]; outputim[k][i] = dataim[ix++];
            }
         }
      }
   }
}

void cmplx_added_datavectorized_to_output
 ( double *datari, double **outputre, double **outputim,
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
      cout << "Updating value starting at " << ix1re << " in data." << endl;

   for(int i=0; i<=deg; i++) // output[dim][i] = data[ix++];
   {
      outputre[dim][i] = datari[ix1re++];
      outputim[dim][i] = datari[ix2im++];
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
            outputre[0][i] = 0.0; outputim[0][i] = 0.0;
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
            outputre[0][i] = datari[ix1re++];
            outputim[0][i] = datari[ix2im++];
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
         outputre[0][i] = datari[ix1re++];
         outputim[0][i] = datari[ix2im++];
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
               outputre[k][i] = 0.0; outputim[k][i] = 0.0;
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
               outputre[k][i] = datari[ix1re++];
               outputim[k][i] = datari[ix2im++];
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
               outputre[k][i] = datari[ix1re++];
               outputim[k][i] = datari[ix2im++];
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
               outputre[k][i] = datari[ix1re++];
               outputim[k][i] = datari[ix2im++];
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
               outputre[k][i] = datari[ix1re++];
               outputim[k][i] = datari[ix2im++];
            }
         }
      }
   }
}

void dbl_data_setup
 ( int dim, int nbr, int deg, int totcff,
   double *data, double *cst, double **cff, double **input )
{
   const int deg1 = deg+1;
   int ix = 0;

   for(int i=0; i<deg1; i++) data[ix++] = cst[i];

   for(int i=0; i<nbr; i++)
      for(int j=0; j<deg1; j++) data[ix++] = cff[i][j];

   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++) data[ix++] = input[i][j];
}

void cmplx_data_setup
 ( int dim, int nbr, int deg, int totcff,
   double *datare, double *dataim, double *cstre, double *cstim,
   double **cffre, double **cffim, double **inputre, double **inputim )
{
   const int deg1 = deg+1;
   int ix = 0;

   for(int i=0; i<deg1; i++)
   {
      datare[ix] = cstre[i]; dataim[ix++] = cstim[i];
   }
   for(int i=0; i<nbr; i++)
      for(int j=0; j<deg1; j++)
      {
         datare[ix] = cffre[i][j]; dataim[ix++] = cffim[i][j];
      }
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++)
      {
         datare[ix] = inputre[i][j]; dataim[ix++] = inputim[i][j];
      }
}

void cmplxvectorized_data_setup
 ( int dim, int nbr, int deg, int totcff, int offsetri,
   double *datari, double *cstre, double *cstim,
   double **cffre, double **cffim, double **inputre, double **inputim )
{
   const int deg1 = deg+1;

   int ix1 = 0;
   int ix2 = totcff + offsetri;

   for(int i=0; i<deg1; i++)
   {
      datari[ix1++] = cstre[i];
      datari[ix2++] = cstim[i];
   }
   for(int i=0; i<nbr; i++)
      for(int j=0; j<deg1; j++)
      {
         datari[ix1++] = cffre[i][j];
         datari[ix2++] = cffim[i][j];
      }
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++)
      {
         datari[ix1++] = inputre[i][j];
         datari[ix2++] = inputim[i][j];
      }

   for(int i=0; i<2*offsetri; i++)
   {
      datari[ix1++] = 0.0;
      datari[ix2++] = 0.0;
   }
}

void dbl_convolution_jobs
 ( int dim, int nbr, int deg, int *nvr, ConvolutionJobs cnvjobs,
   int *fstart, int *bstart, int *cstart,
   double *data, double *cnvlapms, int vrblvl )
{
   const int deg1 = deg+1;

   cudaEvent_t start,stop;           // to measure time spent by kernels
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   *cnvlapms = 0.0;
   float milliseconds;
   const bool vrb  = (vrblvl > 1);

   for(int k=0; k<cnvjobs.get_depth(); k++)
   {
      const int jobnbr = cnvjobs.get_layer_count(k);
      int *in1ix_h = new int[jobnbr];
      int *in2ix_h = new int[jobnbr];
      int *outix_h = new int[jobnbr];

      if(vrblvl > 1) cout << "preparing convolution jobs at layer "
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

         if(vrblvl > 1)
            cout << "launching " << jobnbr << " blocks of " << deg1
                 << " threads for convolutions..." << endl;
         
         cudaEventRecord(start);
         dbl_padded_convjobs<<<jobnbr,deg1>>>
            (data,in1ix_d,in2ix_d,outix_d,deg1);
         cudaEventRecord(stop);
         cudaEventSynchronize(stop);
         cudaEventElapsedTime(&milliseconds,start,stop);
         *cnvlapms += milliseconds;
         cudaFree(in1ix_d); cudaFree(in2ix_d); cudaFree(outix_d);
      }
      free(in1ix_h); free(in2ix_h); free(outix_h);
   }
}

void cmplx_convolution_jobs
 ( int dim, int nbr, int deg, int *nvr, ConvolutionJobs cnvjobs,
   int *fstart, int *bstart, int *cstart,
   double *datare, double *dataim, double *cnvlapms, int vrblvl )
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

      if(vrblvl > 1) cout << "preparing convolution jobs at layer "
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

         if(vrblvl > 1)
            cout << "launching " << jobnbr << " blocks of " << deg1
                 << " threads for convolutions ..." << endl;
         
         cudaEventRecord(start);
         cmplx_padded_convjobs<<<jobnbr,deg1>>>
            (datare,dataim,in1ix_d,in2ix_d,outix_d,deg1);
         cudaEventRecord(stop);
         cudaEventSynchronize(stop);
         cudaEventElapsedTime(&milliseconds,start,stop);
         *cnvlapms += milliseconds;
         cudaFree(in1ix_d); cudaFree(in2ix_d); cudaFree(outix_d);
      }
      free(in1ix_h); free(in2ix_h); free(outix_h);
   }
}

void cmplxvectorized_convolution_jobs
 ( int dim, int nbr, int deg, int *nvr, int totcff, int offsetri,
   ComplexConvolutionJobs cnvjobs, ComplexIncrementJobs incjobs,
   int *fstart, int *bstart, int *cstart,
   double *datari, double *cnvlapms, int vrblvl )
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

      if(vrblvl > 1) cout << "preparing convolution jobs at layer "
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

         if(vrblvl > 1)
            cout << "launching " << jobnbr << " blocks of " << deg1
                 << " threads for convolutions ..." << endl;

         cudaEventRecord(start);
         dbl_padded_convjobs<<<jobnbr,deg1>>>
            (datari,in1ix_d,in2ix_d,outix_d,deg1);
         cudaEventRecord(stop);
         cudaEventSynchronize(stop);
         cudaEventElapsedTime(&milliseconds,start,stop);
         *cnvlapms += milliseconds;
         cudaFree(in1ix_d); cudaFree(in2ix_d); cudaFree(outix_d);
      }
      jobnbr = incjobs.get_layer_count(k);
      // note: only half the number of increment jobs

      if(vrblvl > 1) cout << "preparing increment jobs at layer "
                          << k << " ..." << endl;

      complex_incjobs_coordinates
         (incjobs,k,in1ix_h,in2ix_h,outix_h,dim,nbr,deg,nvr,totcff,offsetri,
          fstart,bstart,cstart,vrb);

      const int nbrflips = jobnbr/2;
      int *rebidx = new int[nbrflips];
      for(int i=0, j=0; i<jobnbr; i=i+2, j++) rebidx[j] = in2ix_h[i];

      GPU_cmplxvectorized_flipsigns
         (deg,nbrflips,rebidx,datari,&fliplapms,vrblvl);
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

         if(vrblvl > 1)
            cout << "launching " << jobnbr << " blocks of " << deg1
                 << " threads for increments ..." << endl;

         cudaEventRecord(start);
         dbl_increment_jobs<<<jobnbr,deg1>>>
            (datari,in1ix_d,in2ix_d,outix_d,deg1);
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

void dbl_addition_jobs
 ( int dim, int nbr, int deg, int *nvr, AdditionJobs addjobs,
   int *fstart, int *bstart, int *cstart,
   double *data, double *addlapms, int vrblvl )
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

      if(vrblvl > 1) cout << "preparing addition jobs at layer "
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

         if(vrblvl > 1)
            cout << "launching " << jobnbr << " blocks of " << deg1
                 << " threads for additions ..." << endl;

         cudaEventRecord(start);
         dbl_update_addjobs<<<jobnbr,deg1>>>
            (data,in1ix_d,in2ix_d,outix_d,deg1);
         cudaEventRecord(stop);
         cudaEventSynchronize(stop);
         cudaEventElapsedTime(&milliseconds,start,stop);
         *addlapms += milliseconds;
         cudaFree(in1ix_d); cudaFree(in2ix_d); cudaFree(outix_d);
      }
      free(in1ix_h); free(in2ix_h); free(outix_h);
   }
}

void cmplx_addition_jobs
 ( int dim, int nbr, int deg, int *nvr, AdditionJobs addjobs,
   int *fstart, int *bstart, int *cstart,
   double *datare, double *dataim, double *addlapms, int vrblvl )
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

      if(vrblvl > 1) cout << "preparing addition jobs at layer "
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

         if(vrblvl > 1)
            cout << "launching " << jobnbr << " blocks of " << deg1
                 << " threads for additions ..." << endl;

         cudaEventRecord(start);
         cmplx_update_addjobs<<<jobnbr,deg1>>>
            (datare,dataim,in1ix_d,in2ix_d,outix_d,deg1);
         cudaEventRecord(stop);
         cudaEventSynchronize(stop);
         cudaEventElapsedTime(&milliseconds,start,stop);
         *addlapms += milliseconds;
         cudaFree(in1ix_d); cudaFree(in2ix_d); cudaFree(outix_d);
      }
      free(in1ix_h); free(in2ix_h); free(outix_h);
   }
}

void cmplxvectorized_addition_jobs
 ( int dim, int nbr, int deg, int *nvr, int totcff, int offsetri,
   ComplexAdditionJobs addjobs,
   int *fstart, int *bstart, int *cstart,
   double *datari, double *addlapms, int vrblvl )
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

      if(vrblvl > 1) cout << "preparing addition jobs at layer "
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

         if(vrblvl > 1)
            cout << "launching " << jobnbr << " blocks of " << deg1
                 << " threads for additions ..." << endl;

         cudaEventRecord(start);
         dbl_update_addjobs<<<jobnbr,deg1>>>
            (datari,in1ix_d,in2ix_d,outix_d,deg1);
         cudaEventRecord(stop);
         cudaEventSynchronize(stop);
         cudaEventElapsedTime(&milliseconds,start,stop);
         *addlapms += milliseconds;
         cudaFree(in1ix_d); cudaFree(in2ix_d); cudaFree(outix_d);
      }
      free(in1ix_h); free(in2ix_h); free(outix_h);
   }
}

void GPU_dbl_poly_evaldiff
 ( int BS, int dim, int nbr, int deg, int *nvr, int **idx,
   double *cst, double **cff, double **input, double **output,
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

   double *data_h = new double[totalcff];        // data on host

   dbl_data_setup(dim,nbr,deg,totalcff,data_h,cst,cff,input);

   double *data_d;                               // device data
   const size_t szdata = totalcff*sizeof(double);
   cudaMalloc((void**)&data_d,szdata);
   cudaMemcpy(data_d,data_h,szdata,cudaMemcpyHostToDevice);

   struct timeval begintime,endtime; // wall clock time of computations

   gettimeofday(&begintime,0);

   dbl_convolution_jobs
      (dim,nbr,deg,nvr,cnvjobs,fstart,bstart,cstart,
       data_d,cnvlapms,vrblvl);

   dbl_addition_jobs
      (dim,nbr,deg,nvr,addjobs,fstart,bstart,cstart,
       data_d,addlapms,vrblvl);

   gettimeofday(&endtime,0);
   cudaMemcpy(data_h,data_d,szdata,cudaMemcpyDeviceToHost);
   *elapsedms = *cnvlapms + *addlapms;
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   *walltimesec = seconds + microseconds*1.0e-6;

   // dbl_convoluted_data_to_output
   //    (data_h,output,dim,nbr,deg,nvr,idx,fstart,bstart,cstart,vrblvl);

   dbl_added_data_to_output
      (data_h,output,dim,nbr,deg,nvr,idx,fstart,bstart,cstart,
       addjobs,vrblvl);

   if(vrblvl > 0)
       write_GPU_timings(*cnvlapms,*addlapms,*elapsedms,*walltimesec);

   cudaFree(data_d);

   free(data_h);

   free(fstart); free(bstart); free(cstart);
   free(fsums); free(bsums); free(csums);
}

void GPU_cmplx_poly_evaldiff
 ( int BS, int dim, int nbr, int deg, int *nvr, int **idx,
   double *cstre, double *cstim, double **cffre, double **cffim,
   double **inputre, double **inputim, double **outputre, double **outputim,
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

   double *datare_h = new double[totalcff];        // data on host
   double *dataim_h = new double[totalcff];

   cmplx_data_setup
      (dim,nbr,deg,totalcff,datare_h,dataim_h,
       cstre,cstim,cffre,cffim,inputre,inputim);

   double *datare_d;                               // device data
   double *dataim_d;
   const size_t szdata = totalcff*sizeof(double);
   cudaMalloc((void**)&datare_d,szdata);
   cudaMalloc((void**)&dataim_d,szdata);
   cudaMemcpy(datare_d,datare_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(dataim_d,dataim_h,szdata,cudaMemcpyHostToDevice);

   struct timeval begintime,endtime; // wall clock time of computations

   gettimeofday(&begintime,0);

   cmplx_convolution_jobs
      (dim,nbr,deg,nvr,cnvjobs,fstart,bstart,cstart,
       datare_d,dataim_d,cnvlapms,vrblvl);

   cmplx_addition_jobs
      (dim,nbr,deg,nvr,addjobs,fstart,bstart,cstart,
       datare_d,dataim_d,addlapms,vrblvl);

   gettimeofday(&endtime,0);
   cudaMemcpy(datare_h,datare_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(dataim_h,dataim_d,szdata,cudaMemcpyDeviceToHost);
   *elapsedms = *cnvlapms + *addlapms;
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   *walltimesec = seconds + microseconds*1.0e-6;

   // cmplx_convoluted_data_to_output
   //    (datare_h,dataim_h,outputre,outputim,
   //     dim,nbr,deg,nvr,idx,fstart,bstart,cstart,vrblvl);

   cmplx_added_data_to_output
      (datare_h,dataim_h,outputre,outputim,dim,nbr,deg,nvr,idx,
       fstart,bstart,cstart,addjobs,vrblvl);

   if(vrblvl > 0)
      write_GPU_timings(*cnvlapms,*addlapms,*elapsedms,*walltimesec);

   cudaFree(datare_d); cudaFree(dataim_d);

   free(datare_h); free(dataim_h);

   free(fstart); free(bstart); free(cstart);
   free(fsums); free(bsums); free(csums);
}

void GPU_cmplxvectorized_poly_evaldiff
 ( int BS, int dim, int nbr, int deg, int *nvr, int **idx,
   double *cstre, double *cstim, double **cffre, double **cffim,
   double **inputre, double **inputim, double **outputre, double **outputim,
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
   double *datari_h = new double[cmplxtotcff];      // data on host

   cmplxvectorized_data_setup
      (dim,nbr,deg,totalcff,offsetri,datari_h,
       cstre,cstim,cffre,cffim,inputre,inputim);

   double *datari_d;                               // device data
   const size_t szdata = cmplxtotcff*sizeof(double);
   cudaMalloc((void**)&datari_d,szdata);
   cudaMemcpy(datari_d,datari_h,szdata,cudaMemcpyHostToDevice);

   struct timeval begintime,endtime; // wall clock time of computations

   gettimeofday(&begintime,0);

   cmplxvectorized_convolution_jobs
      (dim,nbr,deg,nvr,totalcff,offsetri,cnvjobs,incjobs,fstart,bstart,cstart,
       datari_d,cnvlapms,vrblvl);

   cmplxvectorized_addition_jobs
      (dim,nbr,deg,nvr,totalcff,offsetri,addjobs,fstart,bstart,cstart,
       datari_d,addlapms,vrblvl);

   gettimeofday(&endtime,0);
   cudaMemcpy(datari_h,datari_d,szdata,cudaMemcpyDeviceToHost);
   *elapsedms = *cnvlapms + *addlapms;
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   *walltimesec = seconds + microseconds*1.0e-6;

   cmplx_added_datavectorized_to_output
      (datari_h,outputre,outputim,dim,nbr,deg,nvr,idx,fstart,bstart,cstart,
       totalcff,offsetri,addjobs,vrblvl);

   if(vrblvl > 0)
      write_GPU_timings(*cnvlapms,*addlapms,*elapsedms,*walltimesec);

   cudaFree(datari_d);

   free(datari_h);
 
   free(fstart); free(bstart); free(cstart);
   free(fsums); free(bsums); free(csums);
}
