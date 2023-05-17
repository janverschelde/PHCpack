// The file dbl3_polynomials_kernels.cu defines the kernels with prototypes
// in dbl3_polynomials_kernels.h.

#include <iostream>
#include <iomanip>
#ifdef winwalltime
#include "wingettimeofday.h"
#else
#include <sys/time.h>
#endif
#include "job_coordinates.h"
#include "triple_double_functions.h"
#ifdef gpufun
#include "double_double_gpufun.cu"
#include "triple_double_gpufun.cu"
#endif
#include "dbl3_polynomials_kernels.h"
#include "write_gpu_timings.h"

// The constant td_shmemsize is the bound on the shared memory size.

#define td_shmemsize 192

using namespace std;

__global__ void dbl3_padded_convjobs
 ( double *datahi, double *datami, double *datalo,
   int *in1idx, int *in2idx, int *outidx, int dim )
{
   const int bdx = blockIdx.x;           // index to the convolution job
   const int tdx = threadIdx.x;          // index to the output of the job
   const int idx1 = in1idx[bdx] + tdx;
   const int idx2 = in2idx[bdx] + tdx;
   const int idx3 = outidx[bdx] + tdx;

   __shared__ double xvhi[td_shmemsize];
   __shared__ double xvmi[td_shmemsize];
   __shared__ double xvlo[td_shmemsize];
   __shared__ double yvhi[2*td_shmemsize];
   __shared__ double yvmi[2*td_shmemsize];
   __shared__ double yvlo[2*td_shmemsize];
   __shared__ double zvhi[td_shmemsize];
   __shared__ double zvmi[td_shmemsize];
   __shared__ double zvlo[td_shmemsize];

   double prdhi,prdmi,prdlo;
   int ydx = dim + tdx;

   xvhi[tdx] = datahi[idx1];  // loading first input
   xvmi[tdx] = datami[idx1]; 
   xvlo[tdx] = datalo[idx1]; 
   yvhi[tdx] = 0.0;           // padded with zeros
   yvmi[tdx] = 0.0;
   yvlo[tdx] = 0.0;
   yvhi[ydx] = datahi[idx2];  // loading second input
   yvmi[ydx] = datami[idx2];
   yvlo[ydx] = datalo[idx2];

   __syncthreads();

   // zv[tdx] = xv[0]*yv[tdx];
   tdg_mul(xvhi[0],xvmi[0],xvlo[0],yvhi[ydx],yvmi[ydx],yvlo[ydx],
           &zvhi[tdx],&zvmi[tdx],&zvlo[tdx]);
   __syncthreads();

   for(int i=1; i<dim; i++) // zv[tdx] = zv[tdx] + xv[i]*yv[dim+tdx-i];
   {
      ydx = dim + tdx - i;
      tdg_mul(xvhi[i],xvmi[i],xvlo[i],yvhi[ydx],yvmi[ydx],yvlo[ydx],
              &prdhi,&prdmi,&prdlo);
      __syncthreads();
      tdg_inc(&zvhi[tdx],&zvmi[tdx],&zvlo[tdx],prdhi,prdmi,prdlo);
      __syncthreads();
   }

   __syncthreads();

   datahi[idx3] = zvhi[tdx]; // storing the output
   datami[idx3] = zvmi[tdx];
   datalo[idx3] = zvlo[tdx];
}

__global__ void cmplx3_padded_convjobs
 ( double *datarehi, double *dataremi, double *datarelo,
   double *dataimhi, double *dataimmi, double *dataimlo,
   int *in1idx, int *in2idx, int *outidx, int dim )
{
   const int bdx = blockIdx.x;           // index to the convolution job
   const int tdx = threadIdx.x;          // index to the output of the job
   const int idx1 = in1idx[bdx] + tdx;
   const int idx2 = in2idx[bdx] + tdx;
   const int idx3 = outidx[bdx] + tdx;

   __shared__ double xvrehi[td_shmemsize];
   __shared__ double xvremi[td_shmemsize];
   __shared__ double xvrelo[td_shmemsize];
   __shared__ double xvimhi[td_shmemsize];
   __shared__ double xvimmi[td_shmemsize];
   __shared__ double xvimlo[td_shmemsize];
   __shared__ double yvrehi[2*td_shmemsize];
   __shared__ double yvremi[2*td_shmemsize];
   __shared__ double yvrelo[2*td_shmemsize];
   __shared__ double yvimhi[2*td_shmemsize];
   __shared__ double yvimmi[2*td_shmemsize];
   __shared__ double yvimlo[2*td_shmemsize];
   __shared__ double zvrehi[td_shmemsize];
   __shared__ double zvremi[td_shmemsize];
   __shared__ double zvrelo[td_shmemsize];
   __shared__ double zvimhi[td_shmemsize];
   __shared__ double zvimmi[td_shmemsize];
   __shared__ double zvimlo[td_shmemsize];

   double prodhi,prodmi,prodlo;
   int ydx = dim + tdx;

   xvrehi[tdx] = datarehi[idx1];  // loading first input
   xvremi[tdx] = dataremi[idx1]; 
   xvrelo[tdx] = datarelo[idx1]; 
   xvimhi[tdx] = dataimhi[idx1];
   xvimmi[tdx] = dataimmi[idx1];
   xvimlo[tdx] = dataimlo[idx1]; 
   yvrehi[tdx] = 0.0;           // padded with zeros
   yvremi[tdx] = 0.0;
   yvrelo[tdx] = 0.0;
   yvimhi[tdx] = 0.0;
   yvimmi[tdx] = 0.0;
   yvimlo[tdx] = 0.0;
   yvrehi[ydx] = datarehi[idx2];  // loading second input
   yvremi[ydx] = dataremi[idx2];
   yvrelo[ydx] = datarelo[idx2];
   yvimhi[ydx] = dataimhi[idx2];
   yvimmi[ydx] = dataimmi[idx2];
   yvimlo[ydx] = dataimlo[idx2];

   __syncthreads();

   // zv[tdx] = xv[0]*yv[tdx];
   tdg_mul(xvrehi[0],xvremi[0],xvrelo[0],
           yvrehi[ydx],yvremi[ydx],yvrelo[ydx],
           &zvrehi[tdx],&zvremi[tdx],&zvrelo[tdx]);
   __syncthreads();
   tdg_mul(xvimhi[0],xvimmi[0],xvimlo[0],
           yvimhi[ydx],yvimmi[ydx],yvimlo[ydx],
           &prodhi,&prodmi,&prodlo);
   __syncthreads();
   tdg_minus(&prodhi,&prodmi,&prodlo);
   tdg_inc(&zvrehi[tdx],&zvremi[tdx],&zvrelo[tdx],
           prodhi,prodmi,prodlo);
   __syncthreads();

   tdg_mul(xvrehi[0],xvremi[0],xvrelo[0],
           yvimhi[ydx],yvimmi[ydx],yvimlo[ydx],
           &zvimhi[tdx],&zvimmi[tdx],&zvimlo[tdx]);
   __syncthreads();
   tdg_mul(xvimhi[0],xvimmi[0],xvimlo[0],
           yvrehi[ydx],yvremi[ydx],yvrelo[ydx],
           &prodhi,&prodmi,&prodlo);
   __syncthreads();
   tdg_inc(&zvimhi[tdx],&zvimmi[tdx],&zvimlo[tdx],
           prodhi,prodmi,prodlo);
   __syncthreads();

   for(int i=1; i<dim; i++) // zv[tdx] = zv[tdx] + xv[i]*yv[dim+tdx-i];
   {
      ydx = dim + tdx - i;

      tdg_mul(xvrehi[i],xvremi[i],xvrelo[i],
              yvrehi[ydx],yvremi[ydx],yvrelo[ydx],
              &prodhi,&prodmi,&prodlo);
      __syncthreads();
      tdg_inc(&zvrehi[tdx],&zvremi[tdx],&zvrelo[tdx],
              prodhi,prodmi,prodlo);
      __syncthreads();
      tdg_mul(xvimhi[i],xvimmi[i],xvimlo[i],
              yvimhi[ydx],yvimmi[ydx],yvimlo[ydx],
              &prodhi,&prodmi,&prodlo);
      __syncthreads();
      tdg_minus(&prodhi,&prodmi,&prodlo);
      tdg_inc(&zvrehi[tdx],&zvremi[tdx],&zvrelo[tdx],
              prodhi,prodmi,prodlo);
      __syncthreads();

      tdg_mul(xvrehi[i],xvremi[i],xvrelo[i],
              yvimhi[ydx],yvimmi[ydx],yvimlo[ydx],
              &prodhi,&prodmi,&prodlo);
      __syncthreads();
      tdg_inc(&zvimhi[tdx],&zvimmi[tdx],&zvimlo[tdx],
              prodhi,prodmi,prodlo);
      __syncthreads();
      tdg_mul(xvimhi[i],xvimmi[i],xvimlo[i],
              yvrehi[ydx],yvremi[ydx],yvrelo[ydx],
              &prodhi,&prodmi,&prodlo);
      __syncthreads();
      tdg_inc(&zvimhi[tdx],&zvimmi[tdx],&zvimlo[tdx],
              prodhi,prodmi,prodlo);
      __syncthreads();
   }
   __syncthreads();

   datarehi[idx3] = zvrehi[tdx]; // storing the output
   dataremi[idx3] = zvremi[tdx];
   datarelo[idx3] = zvrelo[tdx];
   dataimhi[idx3] = zvimhi[tdx]; 
   dataimmi[idx3] = zvimmi[tdx]; 
   dataimlo[idx3] = zvimlo[tdx];
}

__global__ void dbl3_update_addjobs
 ( double *datahi, double *datami, double *datalo,
   int *in1idx, int *in2idx, int *outidx, int dim )
{
   const int bdx = blockIdx.x;           // index to the convolution job
   const int tdx = threadIdx.x;          // index to the output of the job
   const int idx1 = in1idx[bdx] + tdx;
   const int idx2 = in2idx[bdx] + tdx;
   const int idx3 = outidx[bdx] + tdx;

   __shared__ double xvhi[td_shmemsize];
   __shared__ double xvmi[td_shmemsize];
   __shared__ double xvlo[td_shmemsize];
   __shared__ double yvhi[td_shmemsize];
   __shared__ double yvmi[td_shmemsize];
   __shared__ double yvlo[td_shmemsize];
   __shared__ double zvhi[td_shmemsize];
   __shared__ double zvmi[td_shmemsize];
   __shared__ double zvlo[td_shmemsize];

   xvhi[tdx] = datahi[idx1];  // loading first input
   xvmi[tdx] = datami[idx1];
   xvlo[tdx] = datalo[idx1];
   yvhi[tdx] = datahi[idx2];  // loading second input
   yvmi[tdx] = datami[idx2];
   yvlo[tdx] = datalo[idx2];

   // zv[tdx] = xv[tdx] + yv[tdx];

   __syncthreads();

   tdg_add(xvhi[tdx],xvmi[tdx],xvlo[tdx],yvhi[tdx],yvmi[tdx],yvlo[tdx],
           &zvhi[tdx],&zvmi[tdx],&zvlo[tdx]);

   __syncthreads();

   datahi[idx3] = zvhi[tdx]; // storing the output
   datami[idx3] = zvmi[tdx];
   datalo[idx3] = zvlo[tdx];
}

__global__ void cmplx3_update_addjobs
 ( double *datarehi, double *dataremi, double *datarelo,
   double *dataimhi, double *dataimmi, double *dataimlo,
   int *in1idx, int *in2idx, int *outidx, int dim )
{
   const int bdx = blockIdx.x;           // index to the addition job
   const int tdx = threadIdx.x;          // index to the output of the job
   const int idx1 = in1idx[bdx] + tdx;
   const int idx2 = in2idx[bdx] + tdx;
   const int idx3 = outidx[bdx] + tdx;

   __shared__ double xvrehi[td_shmemsize];
   __shared__ double xvremi[td_shmemsize];
   __shared__ double xvrelo[td_shmemsize];
   __shared__ double xvimhi[td_shmemsize];
   __shared__ double xvimmi[td_shmemsize];
   __shared__ double xvimlo[td_shmemsize];
   __shared__ double yvrehi[td_shmemsize];
   __shared__ double yvremi[td_shmemsize];
   __shared__ double yvrelo[td_shmemsize];
   __shared__ double yvimhi[td_shmemsize];
   __shared__ double yvimmi[td_shmemsize];
   __shared__ double yvimlo[td_shmemsize];
   __shared__ double zvrehi[td_shmemsize];
   __shared__ double zvremi[td_shmemsize];
   __shared__ double zvrelo[td_shmemsize];
   __shared__ double zvimhi[td_shmemsize];
   __shared__ double zvimmi[td_shmemsize];
   __shared__ double zvimlo[td_shmemsize];

   xvrehi[tdx] = datarehi[idx1];  // loading first input
   xvremi[tdx] = dataremi[idx1];
   xvrelo[tdx] = datarelo[idx1];
   xvimhi[tdx] = dataimhi[idx1];
   xvimmi[tdx] = dataimmi[idx1];
   xvimlo[tdx] = dataimlo[idx1];
   yvrehi[tdx] = datarehi[idx2];  // loading second input
   yvremi[tdx] = dataremi[idx2];
   yvrelo[tdx] = datarelo[idx2];
   yvimhi[tdx] = dataimhi[idx2];
   yvimmi[tdx] = dataimmi[idx2];
   yvimlo[tdx] = dataimlo[idx2];

   // zv[tdx] = xv[tdx] + yv[tdx];

   tdg_add(xvrehi[tdx],xvremi[tdx],xvrelo[tdx],
           yvrehi[tdx],yvremi[tdx],yvrelo[tdx],
           &zvrehi[tdx],&zvremi[tdx],&zvrelo[tdx]);
   __syncthreads();

   tdg_add(xvimhi[tdx],xvimmi[tdx],xvimlo[tdx],
           yvimhi[tdx],yvimmi[tdx],yvimlo[tdx],
           &zvimhi[tdx],&zvimmi[tdx],&zvimlo[tdx]);
   __syncthreads();

   datarehi[idx3] = zvrehi[tdx]; // storing the output
   dataremi[idx3] = zvremi[tdx];
   datarelo[idx3] = zvrelo[tdx];
   dataimhi[idx3] = zvimhi[tdx];
   dataimmi[idx3] = zvimmi[tdx];
   dataimlo[idx3] = zvimlo[tdx];
}

void dbl_convoluted_data3_to_output
 ( double *datahi, double *datami, double *datalo,
   double **outputhi, double **outputmi, double **outputlo,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart, bool verbose )
{
   const int deg1 = deg+1;
   int ix0,ix1,ix2;

   for(int i=0; i<=deg; i++) // output[dim][i] = data[i];
   {
      outputhi[dim][i] = datahi[i];
      outputmi[dim][i] = datami[i];
      outputlo[dim][i] = datalo[i];
   }
   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++) // output[i][j] = 0.0;
      {
         outputhi[i][j] = 0.0;
         outputmi[i][j] = 0.0;
         outputlo[i][j] = 0.0;
      }

   for(int k=0; k<nbr; k++)
   {
      ix1 = fstart[k] + (nvr[k]-1)*deg1;
      
      if(verbose)
         cout << "monomial " << k << " update starts at " << ix1 << endl;

      for(int i=0; i<=deg; i++) // output[dim][i] += data[ix1++];
         tdf_inc(&outputhi[dim][i],&outputmi[dim][i],&outputlo[dim][i],
                 datahi[ix1],datami[ix1],datalo[ix1++]);
     
      ix0 = idx[k][0];
      if(nvr[k] == 1)
      {
         ix1 = (1 + k)*deg1;
            
         for(int i=0; i<=deg; i++) // output[ix0][i] += data[ix1++];
            tdf_inc(&outputhi[ix0][i],&outputmi[ix0][i],&outputlo[ix0][i],
                    datahi[ix1],datami[ix1],datalo[ix1++]);
      }
      else
      {                               // update first and last derivative
         ix2 = nvr[k]-3;
         if(ix2 < 0) ix2 = 0;
         ix1 = bstart[k] + ix2*deg1;

         for(int i=0; i<=deg; i++) // output[ix0][i] += data[ix1++];
            tdf_inc(&outputhi[ix0][i],&outputmi[ix0][i],&outputlo[ix0][i],
                    datahi[ix1],datami[ix1],datalo[ix1++]);

         ix2 = nvr[k]-2;
         ix1 = fstart[k] + ix2*deg1;
         ix0 = idx[k][ix2+1];

         for(int i=0; i<=deg; i++) // output[ix0][i] += data[ix1++];
            tdf_inc(&outputhi[ix0][i],&outputmi[ix0][i],&outputlo[ix0][i],
                    datahi[ix1],datami[ix1],datalo[ix1++]);
 
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
                  tdf_inc(&outputhi[ix0][i],&outputmi[ix0][i],
                          &outputlo[ix0][i],
                          datahi[ix1],datami[ix1],datalo[ix1++]);
            }
         }
      }
   }
}

void cmplx_convoluted_data3_to_output
 ( double *datarehi, double *dataremi, double *datarelo,
   double *dataimhi, double *dataimmi, double *dataimlo,
   double **outputrehi, double **outputremi, double **outputrelo,
   double **outputimhi, double **outputimmi, double **outputimlo,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart, bool verbose )
{
   const int deg1 = deg+1;
   int ix0,ix1,ix2;

   for(int i=0; i<=deg; i++) // output[dim][i] = data[i];
   {
      outputrehi[dim][i] = datarehi[i];
      outputremi[dim][i] = dataremi[i];
      outputrelo[dim][i] = datarelo[i];
      outputimhi[dim][i] = dataimhi[i];
      outputimmi[dim][i] = dataimmi[i];
      outputimlo[dim][i] = dataimlo[i];
   }
   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++) // output[i][j] = 0.0;
      {
         outputrehi[i][j] = 0.0;
         outputremi[i][j] = 0.0;
         outputrelo[i][j] = 0.0;
         outputimhi[i][j] = 0.0;
         outputimmi[i][j] = 0.0;
         outputimlo[i][j] = 0.0;
      }

   for(int k=0; k<nbr; k++)
   {
      ix1 = fstart[k] + (nvr[k]-1)*deg1;
      
      if(verbose)
         cout << "monomial " << k << " update starts at " << ix1 << endl;

      for(int i=0; i<=deg; i++) // output[dim][i] += data[ix1++];
      {
         tdf_inc(&outputrehi[dim][i],&outputremi[dim][i],
                 &outputrelo[dim][i],
                 datarehi[ix1],dataremi[ix1],datarelo[ix1++]);
         tdf_inc(&outputimhi[dim][i],&outputimmi[dim][i],
                 &outputimlo[dim][i],
                 dataimhi[ix1],dataimmi[ix1],dataimlo[ix1++]);
      }
      ix0 = idx[k][0];
      if(nvr[k] == 1)
      {
         ix1 = (1 + k)*deg1;
            
         for(int i=0; i<=deg; i++) // output[ix0][i] += data[ix1++];
         {
            tdf_inc(&outputrehi[ix0][i],&outputremi[ix0][i],
                    &outputrelo[ix0][i],
                    datarehi[ix1],dataremi[ix1],datarelo[ix1++]);
            tdf_inc(&outputimhi[ix0][i],&outputimmi[ix0][i],
                    &outputimlo[ix0][i],
                    dataimhi[ix1],dataimmi[ix1],dataimlo[ix1++]);
         }
      }
      else
      {                               // update first and last derivative
         ix2 = nvr[k]-3;
         if(ix2 < 0) ix2 = 0;
         ix1 = bstart[k] + ix2*deg1;

         for(int i=0; i<=deg; i++) // output[ix0][i] += data[ix1++];
         {
            tdf_inc(&outputrehi[ix0][i],&outputremi[ix0][i],
                    &outputrelo[ix0][i],
                    datarehi[ix1],dataremi[ix1],datarelo[ix1++]);
            tdf_inc(&outputimhi[ix0][i],&outputimmi[ix0][i],
                    &outputimlo[ix0][i],
                    dataimhi[ix1],dataimmi[ix1],dataimlo[ix1++]);
         }
         ix2 = nvr[k]-2;
         ix1 = fstart[k] + ix2*deg1;
         ix0 = idx[k][ix2+1];

         for(int i=0; i<=deg; i++) // output[ix0][i] += data[ix1++];
         {
            tdf_inc(&outputrehi[ix0][i],&outputremi[ix0][i],
                    &outputrelo[ix0][i],
                    datarehi[ix1],dataremi[ix1],datarelo[ix1++]);
            tdf_inc(&outputimhi[ix0][i],&outputimmi[ix0][i],
                    &outputimlo[ix0][i],
                    dataimhi[ix1],dataimmi[ix1],dataimlo[ix1++]);
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
                  tdf_inc(&outputrehi[ix0][i],&outputremi[ix0][i],
                          &outputrelo[ix0][i],
                          datarehi[ix1],dataremi[ix1],datarelo[ix1++]);
                  tdf_inc(&outputimhi[ix0][i],&outputimmi[ix0][i],
                          &outputimlo[ix0][i],
                          dataimhi[ix1],dataimmi[ix1],dataimlo[ix1++]);
               }
            }
         }
      }
   }
}

void dbl_added_data3_to_output
 ( double *datahi, double *datami, double *datalo,
   double **outputhi, double **outputmi, double **outputlo,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart,
   AdditionJobs jobs, bool verbose )
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
      outputmi[dim][i] = datami[ix];
      outputlo[dim][i] = datalo[ix++];
   }
   int cnt = jobs.get_differential_count(0);

   if(verbose)
      cout << "Differential count for variable 0 : " << cnt << endl;

   if(cnt == 0) // it could be there is no first variable anywhere ...
   {
      const int difidx = jobs.get_differential_index(0,0);

      if(verbose)
         cout << "Differential index for variable 0 : " << difidx << endl;

      if(difidx < 0)
      {
         for(int i=0; i<=deg; i++) // output[0][i] = 0.0;
         {
            outputhi[0][i] = 0.0;
            outputmi[0][i] = 0.0;
            outputlo[0][i] = 0.0;
         }
      }
      else
      {
         int cffidx = (1 + difidx)*deg1;

         if(verbose)
            cout << "updating derivative with coefficient ..." << endl;

         for(int i=0; i<=deg; i++)
         {
            outputhi[0][i] = datahi[cffidx];
            outputmi[0][i] = datami[cffidx];
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
      
      if(verbose)
         cout << "Updating derivative 0 at " << ix << " in data." << endl;

      for(int i=0; i<=deg; i++) // output[0][i] = data[ix++];
      {
         outputhi[0][i] = datahi[ix];
         outputmi[0][i] = datami[ix];
         outputlo[0][i] = datalo[ix++];
      }
   }
   for(int k=1; k<dim; k++) // updating all other derivatives
   {
      int cnt = jobs.get_differential_count(k);

      if(verbose)
         cout << "Differential count for variable " << k
              << " : " << cnt << endl;

      if(cnt == 0) // it could be there is no variable k anywhere ...
      {
         const int difidx = jobs.get_differential_index(k,0);

         if(verbose)
            cout << "Differential index for variable " << k 
                 << " : " << difidx << endl;

         if(difidx < 0)
         {
            for(int i=0; i<=deg; i++) // output[k][i] = 0.0;
            {
               outputhi[k][i] = 0.0;
               outputmi[k][i] = 0.0;
               outputlo[k][i] = 0.0;
            }
         }
         else
         {
            int cffidx = (1 + difidx)*deg1;

            if(verbose)
               cout << "updating derivative with coefficient ..." << endl;

            for(int i=0; i<=deg; i++)
            {
               outputhi[k][i] = datahi[cffidx];
               outputmi[k][i] = datami[cffidx];
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

            if(verbose)
               cout << "Updating derivative " << k 
                    << " at " << ix << " in data." << endl;

            ix = bstart[ix0] + ix2*deg1;

            for(int i=0; i<=deg; i++) // output[k][i] = data[ix++];
            {
               outputhi[k][i] = datahi[ix];
               outputmi[k][i] = datami[ix];
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
               outputmi[k][i] = datami[ix];
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
               outputmi[k][i] = datami[ix];
               outputlo[k][i] = datalo[ix++];
            }
         }
      }
   }
}

void cmplx_added_data3_to_output
 ( double *datarehi, double *dataremi, double *datarelo,
   double *dataimhi, double *dataimmi, double *dataimlo,
   double **outputrehi, double **outputremi, double **outputrelo,
   double **outputimhi, double **outputimmi, double **outputimlo,
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
      outputremi[dim][i] = dataremi[ix];
      outputrelo[dim][i] = datarelo[ix];
      outputimhi[dim][i] = dataimhi[ix];
      outputimmi[dim][i] = dataimmi[ix];
      outputimlo[dim][i] = dataimlo[ix++];
   }
   int cnt = jobs.get_differential_count(0);

   if(verbose)
      cout << "Differential count for variable 0 : " << cnt << endl;

   if(cnt == 0) // it could be there is no first variable anywhere ...
   {
      const int difidx = jobs.get_differential_index(0,0);

      if(verbose)
         cout << "Differential index for variable 0 : " << difidx << endl;

      if(difidx < 0)
      {
         for(int i=0; i<=deg; i++) // output[0][i] = 0.0;
         {
            outputrehi[0][i] = 0.0; outputremi[0][i] = 0.0;
            outputrelo[0][i] = 0.0;
            outputimhi[0][i] = 0.0; outputimmi[0][i] = 0.0; 
            outputimlo[0][i] = 0.0;
         }
      }
      else
      {
         int cffidx = (1 + difidx)*deg1;

         if(verbose)
            cout << "updating derivative with coefficient ..." << endl;

         for(int i=0; i<=deg; i++)
         {
            outputrehi[0][i] = datarehi[cffidx];
            outputremi[0][i] = dataremi[cffidx];
            outputrelo[0][i] = datarelo[cffidx];
            outputimhi[0][i] = dataimhi[cffidx];
            outputimmi[0][i] = dataimmi[cffidx];
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
      
      if(verbose)
         cout << "Updating derivative 0 at " << ix << " in data." << endl;

      for(int i=0; i<=deg; i++) // output[0][i] = data[ix++];
      {
         outputrehi[0][i] = datarehi[ix]; outputremi[0][i] = dataremi[ix];
         outputrelo[0][i] = datarelo[ix];
         outputimhi[0][i] = dataimhi[ix]; outputimmi[0][i] = dataimmi[ix];
         outputimlo[0][i] = dataimlo[ix++];
      }
   }
   for(int k=1; k<dim; k++) // updating all other derivatives
   {
      int cnt = jobs.get_differential_count(k);

      if(verbose)
         cout << "Differential count for variable " << k
              << " : " << cnt << endl;

      if(cnt == 0) // it could be there is no variable k anywhere ...
      {
         const int difidx = jobs.get_differential_index(k,0);

         if(verbose)
            cout << "Differential index for variable " << k 
                 << " : " << difidx << endl;

         if(difidx < 0)
         {
            for(int i=0; i<=deg; i++) // output[k][i] = 0.0;
            {
               outputrehi[k][i] = 0.0; outputremi[k][i] = 0.0;
               outputrelo[k][i] = 0.0;
               outputimhi[k][i] = 0.0; outputimmi[k][i] = 0.0;
               outputimlo[k][i] = 0.0;
            }
         }
         else
         {
            int cffidx = (1 + difidx)*deg1;

            if(verbose)
               cout << "updating derivative with coefficient ..." << endl;

            for(int i=0; i<=deg; i++)
            {
               outputrehi[k][i] = datarehi[cffidx];
               outputremi[k][i] = dataremi[cffidx];
               outputrelo[k][i] = datarelo[cffidx];
               outputimhi[k][i] = dataimhi[cffidx];
               outputimmi[k][i] = dataimmi[cffidx];
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

            if(verbose)
               cout << "Updating derivative " << k 
                    << " at " << ix << " in data." << endl;

            ix = bstart[ix0] + ix2*deg1;

            for(int i=0; i<=deg; i++) // output[k][i] = data[ix++];
            {
               outputrehi[k][i] = datarehi[ix];
               outputremi[k][i] = dataremi[ix];
               outputrelo[k][i] = datarelo[ix];
               outputimhi[k][i] = dataimhi[ix];
               outputimmi[k][i] = dataimmi[ix];
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
               outputremi[k][i] = dataremi[ix];
               outputrelo[k][i] = datarelo[ix];
               outputimhi[k][i] = dataimhi[ix];
               outputimmi[k][i] = dataimmi[ix];
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
               outputremi[k][i] = dataremi[ix];
               outputrelo[k][i] = datarelo[ix];
               outputimhi[k][i] = dataimhi[ix];
               outputimmi[k][i] = dataimmi[ix];
               outputimlo[k][i] = dataimlo[ix++];
            }
         }
      }
   }
}

void GPU_dbl3_poly_evaldiff
 ( int BS, int dim, int nbr, int deg, int *nvr, int **idx,
   double *csthi, double *cstmi, double *cstlo,
   double **cffhi, double **cffmi, double **cfflo,
   double **inputhi, double **inputmi,  double **inputlo,
   double **outputhi, double **outputmi, double **outputlo,
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
   double *datami_h = new double[totalcff];
   double *datalo_h = new double[totalcff];
   int ix = 0;
   for(int i=0; i<deg1; i++)
   {
      datahi_h[ix] = csthi[i];
      datami_h[ix] = cstmi[i];
      datalo_h[ix++] = cstlo[i];
   }
   for(int i=0; i<nbr; i++)
      for(int j=0; j<deg1; j++)
      {
         datahi_h[ix] = cffhi[i][j];
         datami_h[ix] = cffmi[i][j];
         datalo_h[ix++] = cfflo[i][j];
      }
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++)
      {
         datahi_h[ix] = inputhi[i][j];
         datami_h[ix] = inputmi[i][j];
         datalo_h[ix++] = inputlo[i][j];
      }

   double *datahi_d;                               // device data
   double *datami_d;
   double *datalo_d;
   const size_t szdata = totalcff*sizeof(double);
   cudaMalloc((void**)&datahi_d,szdata);
   cudaMalloc((void**)&datami_d,szdata);
   cudaMalloc((void**)&datalo_d,szdata);
   cudaMemcpy(datahi_d,datahi_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datami_d,datami_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datalo_d,datalo_h,szdata,cudaMemcpyHostToDevice);

   cudaEvent_t start,stop;           // to measure time spent by kernels
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
         dbl3_padded_convjobs<<<jobnbr,BS>>>
            (datahi_d,datami_d,datalo_d,in1ix_d,in2ix_d,outix_d,deg1);
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
         dbl3_update_addjobs<<<jobnbr,BS>>>
            (datahi_d,datami_d,datalo_d,in1ix_d,in2ix_d,outix_d,deg1);
         cudaEventRecord(stop);
         cudaEventSynchronize(stop);
         cudaEventElapsedTime(&milliseconds,start,stop);
         *addlapms += milliseconds;
      }
      free(in1ix_h); free(in2ix_h); free(outix_h);
   }
   gettimeofday(&endtime,0);
   cudaMemcpy(datahi_h,datahi_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datami_h,datami_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datalo_h,datalo_d,szdata,cudaMemcpyDeviceToHost);
   *elapsedms = *cnvlapms + *addlapms;
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   *walltimesec = seconds + microseconds*1.0e-6;

   // dbl_convoluted_data2_to_output
   //    (data_h,output,dim,nbr,deg,nvr,idx,fstart,bstart,cstart,verbose);
   dbl_added_data3_to_output
      (datahi_h,datami_h,datalo_h,outputhi,outputmi,outputlo,
       dim,nbr,deg,nvr,idx,fstart,bstart,cstart,addjobs,verbose);

   if(verbose) write_GPU_timings(*cnvlapms,*addlapms,*elapsedms,*walltimesec);
}

void GPU_cmplx3_poly_evaldiff
 ( int BS, int dim, int nbr, int deg, int *nvr, int **idx,
   double *cstrehi, double *cstremi, double *cstrelo,
   double *cstimhi, double *cstimmi, double *cstimlo,
   double **cffrehi, double **cffremi, double **cffrelo,
   double **cffimhi, double **cffimmi, double **cffimlo,
   double **inputrehi, double **inputremi, double **inputrelo,
   double **inputimhi, double **inputimmi, double **inputimlo,
   double **outputrehi, double **outputremi, double **outputrelo,
   double **outputimhi, double **outputimmi, double **outputimlo,
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
   double *dataremi_h = new double[totalcff];
   double *datarelo_h = new double[totalcff];
   double *dataimhi_h = new double[totalcff]; 
   double *dataimmi_h = new double[totalcff]; 
   double *dataimlo_h = new double[totalcff];
   int ix = 0;
   for(int i=0; i<deg1; i++)
   {
      datarehi_h[ix] = cstrehi[i]; dataremi_h[ix] = cstremi[i];
      datarelo_h[ix] = cstrelo[i];
      dataimhi_h[ix] = cstimhi[i]; dataimmi_h[ix] = cstimmi[i];
      dataimlo_h[ix++] = cstimlo[i];
   }
   for(int i=0; i<nbr; i++)
      for(int j=0; j<deg1; j++)
      {
         datarehi_h[ix] = cffrehi[i][j]; dataremi_h[ix] = cffremi[i][j];
         datarelo_h[ix] = cffrelo[i][j];
         dataimhi_h[ix] = cffimhi[i][j]; dataimmi_h[ix] = cffimmi[i][j];
         dataimlo_h[ix++] = cffimlo[i][j];
      }
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++)
      {
         datarehi_h[ix] = inputrehi[i][j]; dataremi_h[ix] = inputremi[i][j];
         datarelo_h[ix] = inputrelo[i][j];
         dataimhi_h[ix] = inputimhi[i][j]; dataimmi_h[ix] = inputimmi[i][j];
         dataimlo_h[ix++] = inputimlo[i][j];
      }

   double *datarehi_d;                               // device data
   double *dataremi_d;
   double *datarelo_d;
   double *dataimhi_d;
   double *dataimmi_d;
   double *dataimlo_d;
   const size_t szdata = totalcff*sizeof(double);
   cudaMalloc((void**)&datarehi_d,szdata);
   cudaMalloc((void**)&dataremi_d,szdata);
   cudaMalloc((void**)&datarelo_d,szdata);
   cudaMalloc((void**)&dataimhi_d,szdata);
   cudaMalloc((void**)&dataimmi_d,szdata);
   cudaMalloc((void**)&dataimlo_d,szdata);
   cudaMemcpy(datarehi_d,datarehi_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(dataremi_d,dataremi_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datarelo_d,datarelo_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(dataimhi_d,dataimhi_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(dataimmi_d,dataimmi_h,szdata,cudaMemcpyHostToDevice);
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
         cmplx3_padded_convjobs<<<jobnbr,BS>>>
            (datarehi_d,dataremi_d,datarelo_d,
             dataimhi_d,dataimmi_d,dataimlo_d,
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
         cmplx3_update_addjobs<<<jobnbr,BS>>>
            (datarehi_d,dataremi_d,datarelo_d,
             dataimhi_d,dataimmi_d,dataimlo_d,
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
   cudaMemcpy(dataremi_h,dataremi_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datarelo_h,datarelo_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(dataimhi_h,dataimhi_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(dataimmi_h,dataimmi_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(dataimlo_h,dataimlo_d,szdata,cudaMemcpyDeviceToHost);
   *elapsedms = *cnvlapms + *addlapms;
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   *walltimesec = seconds + microseconds*1.0e-6;

   // cmplx_convoluted_data3_to_output
   //    (datarehi_h,dataremi_h,datarelo_h,dataimhi_h,dataimmi_h,dataimlo_h,
   //     outputrehi,outputremi,outputrelo,outputimhi,outputimmi,outputimlo,
   //     dim,nbr,deg,nvr,idx,fstart,bstart,cstart,verbose);
   cmplx_added_data3_to_output
      (datarehi_h,dataremi_h,datarelo_h,dataimhi_h,dataimmi_h,dataimlo_h,
       outputrehi,outputremi,outputrelo,outputimhi,outputimmi,outputimlo,
       dim,nbr,deg,nvr,idx,fstart,bstart,cstart,addjobs,verbose);

   if(verbose) write_GPU_timings(*cnvlapms,*addlapms,*elapsedms,*walltimesec);
}
