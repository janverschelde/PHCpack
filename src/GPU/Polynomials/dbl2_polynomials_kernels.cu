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
   int **idx, int *fstart, int *bstart, int *cstart, bool verbose )
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
      
      if(verbose)
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

               if(verbose)
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
   int **idx, int *fstart, int *bstart, int *cstart, bool verbose )
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
      
      if(verbose)
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

               if(verbose)
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
   bool verbose )
{
   const int deg1 = deg + 1;
   const int lastmon = nbr-1;
   const int lastidx = nvr[lastmon]-1;
   int ix;

   ix = fstart[lastmon] + lastidx*deg1;

   if(verbose)
      cout << "Updating value starting at " << ix << " in data." << endl;

   for(int i=0; i<=deg; i++) // output[dim][i] = data[ix++];
   {
      outputhi[dim][i] = datahi[ix];
      outputlo[dim][i] = datalo[ix++];
   }
   int cnt = jobs.get_differential_count(0);
   if(cnt == 0) // it could be there is no first variable anywhere ...
   {
      for(int i=0; i<=deg; i++) // output[0][i] = 0.0;
      {
         outputhi[0][i] = 0.0;
         outputlo[0][i] = 0.0;
      }
   }
   else
   {
      int ix0 = jobs.get_differential_index(0,cnt);
      int ix2 = nvr[ix0]-3;
      if(ix2 < 0) ix2 = 0; // on GPU, one backward item less

      ix = bstart[ix0] + ix2*deg1;
      
      if(verbose)
         cout << "Updating derivative 0 at " << ix << " in data." << endl;

      for(int i=0; i<=deg; i++) // output[0][i] = data[ix++];
      {
         outputhi[0][i] = datahi[ix];
         outputlo[0][i] = datalo[ix++];
      }
      for(int k=1; k<dim; k++) // updating all other derivatives
      {
         int cnt = jobs.get_differential_count(k);
         if(cnt == 0) // it could be there is no variable k anywhere ...
         {
            for(int i=0; i<=deg; i++) // output[k][i] = 0.0;
            {
               outputhi[k][i] = 0.0;
               outputlo[k][i] = 0.0;
            }
         }
         else
         {
            int ix0 = jobs.get_differential_index(k,cnt);
   
            if(idx[ix0][0] == k) // k is first variable of monomial
            {
               int ix2 = nvr[ix0]-3;
               if(ix2 < 0) ix2 = 0;

               if(verbose)
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
   
               if(verbose)
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
   
               if(verbose)
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
}

void cmplx_added_data2_to_output
 ( double *datarehi, double *datarelo,
   double *dataimhi, double *dataimlo,
   double **outputrehi, double **outputrelo,
   double **outputimhi, double **outputimlo,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart, AdditionJobs jobs,
   bool verbose )
{
   const int deg1 = deg + 1;
   const int lastmon = nbr-1;
   const int lastidx = nvr[lastmon]-1;
   int ix;

   ix = fstart[lastmon] + lastidx*deg1;

   if(verbose)
      cout << "Updating value starting at " << ix << " in data." << endl;

   for(int i=0; i<=deg; i++) // output[dim][i] = data[ix++];
   {
      outputrehi[dim][i] = datarehi[ix];
      outputrelo[dim][i] = datarelo[ix];
      outputimhi[dim][i] = dataimhi[ix];
      outputimlo[dim][i] = dataimlo[ix++];
   }
   int cnt = jobs.get_differential_count(0);
   if(cnt == 0) // it could be there is no first variable anywhere ...
   {
      for(int i=0; i<=deg; i++) // output[0][i] = 0.0;
      {
         outputrehi[0][i] = 0.0; outputrelo[0][i] = 0.0;
         outputimhi[0][i] = 0.0; outputimlo[0][i] = 0.0;
      }
   }
   else
   {
      int ix0 = jobs.get_differential_index(0,cnt);
      int ix2 = nvr[ix0]-3;
      if(ix2 < 0) ix2 = 0; // on GPU, one backward item less

      ix = bstart[ix0] + ix2*deg1;
      
      if(verbose)
         cout << "Updating derivative 0 at " << ix << " in data." << endl;

      for(int i=0; i<=deg; i++) // output[0][i] = data[ix++];
      {
         outputrehi[0][i] = datarehi[ix];
         outputrelo[0][i] = datarelo[ix];
         outputimhi[0][i] = dataimhi[ix];
         outputimlo[0][i] = dataimlo[ix++];
      }
      for(int k=1; k<dim; k++) // updating all other derivatives
      {
         int cnt = jobs.get_differential_count(k);
         if(cnt == 0) // it could be there is no variable k anywhere ...
         {
            for(int i=0; i<=deg; i++) // output[k][i] = 0.0;
            {
               outputrehi[k][i] = 0.0; outputrelo[k][i] = 0.0;
               outputimhi[k][i] = 0.0; outputimlo[k][i] = 0.0;
            }
         }
         else
         {
            int ix0 = jobs.get_differential_index(k,cnt);
   
            if(idx[ix0][0] == k) // k is first variable of monomial
            {
               int ix2 = nvr[ix0]-3;
               if(ix2 < 0) ix2 = 0;

               if(verbose)
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
   
               if(verbose)
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
   
               if(verbose)
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
}

void GPU_dbl2_poly_evaldiff
 ( int BS, int dim, int nbr, int deg, int *nvr, int **idx,
   double *csthi, double *cstlo, double **cffhi, double **cfflo,
   double **inputhi, double **inputlo,
   double **outputhi, double **outputlo,
   ConvolutionJobs cnvjobs, AdditionJobs addjobs,
   double *cnvlapms, double *addlapms, double *elapsedms,
   double *walltimesec, bool verbose )
{
   const int deg1 = deg+1;
   const int totalcff = coefficient_count(dim,nbr,deg,nvr);

   int *fstart = new int[nbr];
   int *bstart = new int[nbr];
   int *cstart = new int[nbr];
   int *fsums = new int[nbr];
   int *bsums = new int[nbr];
   int *csums = new int[nbr];

   coefficient_indices
      (dim,nbr,deg,nvr,fsums,bsums,csums,fstart,bstart,cstart);

   if(verbose)
      write_coefficient_indices
         (totalcff,nbr,fsums,fstart,bsums,bstart,csums,cstart);

   double *datahi_h = new double[totalcff];        // data on host
   double *datalo_h = new double[totalcff];
   int ix = 0;
   for(int i=0; i<deg1; i++)
   {
      datahi_h[ix] = csthi[i];
      datalo_h[ix++] = cstlo[i];
   }
   for(int i=0; i<nbr; i++)
      for(int j=0; j<deg1; j++)
      {
         datahi_h[ix] = cffhi[i][j];
         datalo_h[ix++] = cfflo[i][j];
      }
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++)
      {
         datahi_h[ix] = inputhi[i][j];
         datalo_h[ix++] = inputlo[i][j];
      }

   double *datahi_d;                               // device data
   double *datalo_d;
   const size_t szdata = totalcff*sizeof(double);
   cudaMalloc((void**)&datahi_d,szdata);
   cudaMalloc((void**)&datalo_d,szdata);
   cudaMemcpy(datahi_d,datahi_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datalo_d,datalo_h,szdata,cudaMemcpyHostToDevice);

   cudaEvent_t start,stop;           // to measture time spent by kernels
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   *cnvlapms = 0.0;
   *addlapms = 0.0;
   float milliseconds;
   struct timeval begintime,endtime; // wall clock time of computations

   gettimeofday(&begintime,0);
   for(int k=0; k<cnvjobs.get_depth(); k++)
   {
      const int jobnbr = cnvjobs.get_layer_count(k);
      int *in1ix_h = new int[jobnbr];
      int *in2ix_h = new int[jobnbr];
      int *outix_h = new int[jobnbr];

      if(verbose) cout << "preparing convolution jobs at layer "
                       << k << " ..." << endl;

      convjobs_coordinates(cnvjobs,k,in1ix_h,in2ix_h,outix_h,dim,nbr,deg,nvr,
                           fstart,bstart,cstart,verbose);
      if(deg1 == BS)
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

         if(verbose)
            cout << "launching " << jobnbr << " blocks of " << BS
                 << " threads ..." << endl;

         cudaEventRecord(start);
         dbl2_padded_convjobs<<<jobnbr,BS>>>
            (datahi_d,datalo_d,in1ix_d,in2ix_d,outix_d,deg1);
         cudaEventRecord(stop);
         cudaEventSynchronize(stop);
         cudaEventElapsedTime(&milliseconds,start,stop);
         *cnvlapms += milliseconds;
      }
      free(in1ix_h); free(in2ix_h); free(outix_h);
   }
   for(int k=0; k<addjobs.get_depth(); k++)
   {
      const int jobnbr = addjobs.get_layer_count(k);
      int *in1ix_h = new int[jobnbr];
      int *in2ix_h = new int[jobnbr];
      int *outix_h = new int[jobnbr];

      if(verbose) cout << "preparing addition jobs at layer "
                       << k << " ..." << endl;

      addjobs_coordinates(addjobs,k,in1ix_h,in2ix_h,outix_h,dim,nbr,deg,nvr,
                          fstart,bstart,cstart,verbose);
      if(deg1 == BS)
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

         if(verbose)
            cout << "launching " << jobnbr << " blocks of " << BS
                 << " threads ..." << endl;

         cudaEventRecord(start);
         dbl2_update_addjobs<<<jobnbr,BS>>>
            (datahi_d,datalo_d,in1ix_d,in2ix_d,outix_d,deg1);
         cudaEventRecord(stop);
         cudaEventSynchronize(stop);
         cudaEventElapsedTime(&milliseconds,start,stop);
         *addlapms += milliseconds;
      }
      free(in1ix_h); free(in2ix_h); free(outix_h);
   }
   gettimeofday(&endtime,0);
   cudaMemcpy(datahi_h,datahi_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datalo_h,datalo_d,szdata,cudaMemcpyDeviceToHost);
   *elapsedms = *cnvlapms + *addlapms;
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   *walltimesec = seconds + microseconds*1.0e-6;

   // dbl_convoluted_data2_to_output
   //    (datahi_h,datalo_h,outputhi,outputlo,
   //     dim,nbr,deg,nvr,idx,fstart,bstart,cstart,verbose);
   dbl_added_data2_to_output
      (datahi_h,datalo_h,outputhi,outputlo,
       dim,nbr,deg,nvr,idx,fstart,bstart,cstart,addjobs,verbose);

   if(verbose) write_GPU_timings(*cnvlapms,*addlapms,*elapsedms,*walltimesec);
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
   double *walltimesec, bool verbose )
{
   const int deg1 = deg+1;
   const int totalcff = coefficient_count(dim,nbr,deg,nvr);

   int *fstart = new int[nbr];
   int *bstart = new int[nbr];
   int *cstart = new int[nbr];
   int *fsums = new int[nbr];
   int *bsums = new int[nbr];
   int *csums = new int[nbr];

   coefficient_indices
      (dim,nbr,deg,nvr,fsums,bsums,csums,fstart,bstart,cstart);

   if(verbose)
      write_coefficient_indices
         (totalcff,nbr,fsums,fstart,bsums,bstart,csums,cstart);

   double *datarehi_h = new double[totalcff];      // data on host
   double *datarelo_h = new double[totalcff];
   double *dataimhi_h = new double[totalcff]; 
   double *dataimlo_h = new double[totalcff];
   int ix = 0;
   for(int i=0; i<deg1; i++)
   {
      datarehi_h[ix] = cstrehi[i]; datarelo_h[ix] = cstrelo[i];
      dataimhi_h[ix] = cstimhi[i]; dataimlo_h[ix++] = cstimlo[i];
   }
   for(int i=0; i<nbr; i++)
      for(int j=0; j<deg1; j++)
      {
         datarehi_h[ix] = cffrehi[i][j]; datarelo_h[ix] = cffrelo[i][j];
         dataimhi_h[ix] = cffimhi[i][j]; dataimlo_h[ix++] = cffimlo[i][j];
      }
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++)
      {
         datarehi_h[ix] = inputrehi[i][j]; datarelo_h[ix] = inputrelo[i][j];
         dataimhi_h[ix] = inputimhi[i][j]; dataimlo_h[ix++] = inputimlo[i][j];
      }

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

   cudaEvent_t start,stop;           // to measture time spent by kernels
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   *cnvlapms = 0.0;
   *addlapms = 0.0;
   float milliseconds;
   struct timeval begintime,endtime; // wall clock time of computations

   gettimeofday(&begintime,0);
   for(int k=0; k<cnvjobs.get_depth(); k++)
   {
      const int jobnbr = cnvjobs.get_layer_count(k);
      int *in1ix_h = new int[jobnbr];
      int *in2ix_h = new int[jobnbr];
      int *outix_h = new int[jobnbr];

      if(verbose) cout << "preparing convolution jobs at layer "
                       << k << " ..." << endl;

      convjobs_coordinates(cnvjobs,k,in1ix_h,in2ix_h,outix_h,dim,nbr,deg,nvr,
                           fstart,bstart,cstart,verbose);
      if(deg1 == BS)
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

         if(verbose)
            cout << "launching " << jobnbr << " blocks of " << BS
                 << " threads ..." << endl;

         cudaEventRecord(start);
         cmplx2_padded_convjobs<<<jobnbr,BS>>>
            (datarehi_d,datarelo_d,dataimhi_d,dataimlo_d,
             in1ix_d,in2ix_d,outix_d,deg1);
         cudaEventRecord(stop);
         cudaEventSynchronize(stop);
         cudaEventElapsedTime(&milliseconds,start,stop);
         *cnvlapms += milliseconds;
      }
      free(in1ix_h); free(in2ix_h); free(outix_h);
   }
   for(int k=0; k<addjobs.get_depth(); k++)
   {
      const int jobnbr = addjobs.get_layer_count(k);
      int *in1ix_h = new int[jobnbr];
      int *in2ix_h = new int[jobnbr];
      int *outix_h = new int[jobnbr];

      if(verbose) cout << "preparing addition jobs at layer "
                       << k << " ..." << endl;

      addjobs_coordinates(addjobs,k,in1ix_h,in2ix_h,outix_h,dim,nbr,deg,nvr,
                          fstart,bstart,cstart,verbose);
      if(deg1 == BS)
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

         if(verbose)
            cout << "launching " << jobnbr << " blocks of " << BS
                 << " threads ..." << endl;

         cudaEventRecord(start);
         cmplx2_update_addjobs<<<jobnbr,BS>>>
            (datarehi_d,datarelo_d,dataimhi_d,dataimlo_d,
             in1ix_d,in2ix_d,outix_d,deg1);
         cudaEventRecord(stop);
         cudaEventSynchronize(stop);
         cudaEventElapsedTime(&milliseconds,start,stop);
         *addlapms += milliseconds;
      }
      free(in1ix_h); free(in2ix_h); free(outix_h);
   }
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
   //     dim,nbr,deg,nvr,idx,fstart,bstart,cstart,verbose);
   cmplx_added_data2_to_output
      (datarehi_h,datarelo_h,dataimhi_h,dataimlo_h,
       outputrehi,outputrelo,outputimhi,outputimlo,
       dim,nbr,deg,nvr,idx,fstart,bstart,cstart,addjobs,verbose);

   if(verbose) write_GPU_timings(*cnvlapms,*addlapms,*elapsedms,*walltimesec);
}
