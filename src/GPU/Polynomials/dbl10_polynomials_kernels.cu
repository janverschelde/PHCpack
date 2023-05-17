// The file dbl10_polynomials_kernels.cu defines the kernels with prototypes
// in dbl10_polynomials_kernels.h.

#include <iostream>
#include <iomanip>
#ifdef winwalltime
#include "wingettimeofday.h"
#else
#include <sys/time.h>
#endif
#include "job_coordinates.h"
#include "deca_double_functions.h"
#ifdef gpufun
#include "double_double_gpufun.cu"
#include "quad_double_gpufun.cu"
#include "octo_double_gpufun.cu"
#include "deca_double_gpufun.cu"
#endif
#include "dbl10_polynomials_kernels.h"

// The constant da_shmemsize is the bound on the shared memory size.

#define da_shmemsize 153

using namespace std;

__global__ void dbl10_padded_convjobs
 ( double *datartb, double *datarix, double *datarmi, 
   double *datarrg, double *datarpk,
   double *dataltb, double *datalix, double *datalmi, 
   double *datalrg, double *datalpk,
   int *in1idx, int *in2idx, int *outidx, int dim )
{
   const int bdx = blockIdx.x;           // index to the convolution job
   const int tdx = threadIdx.x;          // index to the output of the job
   const int idx1 = in1idx[bdx] + tdx;
   const int idx2 = in2idx[bdx] + tdx;
   const int idx3 = outidx[bdx] + tdx;

   __shared__ double xvrtb[da_shmemsize];
   __shared__ double xvrix[da_shmemsize];
   __shared__ double xvrmi[da_shmemsize];
   __shared__ double xvrrg[da_shmemsize];
   __shared__ double xvrpk[da_shmemsize];
   __shared__ double xvltb[da_shmemsize];
   __shared__ double xvlix[da_shmemsize];
   __shared__ double xvlmi[da_shmemsize];
   __shared__ double xvlrg[da_shmemsize];
   __shared__ double xvlpk[da_shmemsize];
   __shared__ double yvrtb[2*da_shmemsize];
   __shared__ double yvrix[2*da_shmemsize];
   __shared__ double yvrmi[2*da_shmemsize];
   __shared__ double yvrrg[2*da_shmemsize];
   __shared__ double yvrpk[2*da_shmemsize];
   __shared__ double yvltb[2*da_shmemsize];
   __shared__ double yvlix[2*da_shmemsize];
   __shared__ double yvlmi[2*da_shmemsize];
   __shared__ double yvlrg[2*da_shmemsize];
   __shared__ double yvlpk[2*da_shmemsize];
   __shared__ double zvrtb[da_shmemsize];
   __shared__ double zvrix[da_shmemsize];
   __shared__ double zvrmi[da_shmemsize];
   __shared__ double zvrrg[da_shmemsize];
   __shared__ double zvrpk[da_shmemsize];
   __shared__ double zvltb[da_shmemsize];
   __shared__ double zvlix[da_shmemsize];
   __shared__ double zvlmi[da_shmemsize];
   __shared__ double zvlrg[da_shmemsize];
   __shared__ double zvlpk[da_shmemsize];

   double prdrtb,prdrix,prdrmi,prdrrg,prdrpk;
   double prdltb,prdlix,prdlmi,prdlrg,prdlpk;
   int ydx = dim + tdx;

   xvrtb[tdx] = datartb[idx1];  // loading first input
   xvrix[tdx] = datarix[idx1]; 
   xvrmi[tdx] = datarmi[idx1]; 
   xvrrg[tdx] = datarrg[idx1]; 
   xvrpk[tdx] = datarpk[idx1]; 
   xvltb[tdx] = dataltb[idx1];
   xvlix[tdx] = datalix[idx1]; 
   xvlmi[tdx] = datalmi[idx1]; 
   xvlrg[tdx] = datalrg[idx1]; 
   xvlpk[tdx] = datalpk[idx1]; 
   yvrtb[tdx] = 0.0;             // padded with zeros
   yvrix[tdx] = 0.0;
   yvrmi[tdx] = 0.0;
   yvrrg[tdx] = 0.0;
   yvrpk[tdx] = 0.0;
   yvltb[tdx] = 0.0;
   yvlix[tdx] = 0.0;
   yvlmi[tdx] = 0.0;
   yvlrg[tdx] = 0.0;
   yvlpk[tdx] = 0.0;
   yvrtb[ydx] = datartb[idx2];  // loading second input
   yvrix[ydx] = datarix[idx2];
   yvrmi[ydx] = datarmi[idx2];
   yvrrg[ydx] = datarrg[idx2];
   yvrpk[ydx] = datarpk[idx2];
   yvltb[ydx] = dataltb[idx2]; 
   yvlix[ydx] = datalix[idx2];
   yvlmi[ydx] = datalmi[idx2];
   yvlrg[ydx] = datalrg[idx2];
   yvlpk[ydx] = datalpk[idx2];

   __syncthreads();

   // zv[tdx] = xv[0]*yv[tdx];
   dag_mul( xvrtb[0],   xvrix[0],   xvrmi[0],   xvrrg[0],   xvrpk[0],
            xvltb[0],   xvlix[0],   xvlmi[0],   xvlrg[0],   xvlpk[0],
            yvrtb[ydx], yvrix[ydx], yvrmi[ydx], yvrrg[ydx], yvrpk[ydx],
            yvltb[ydx], yvlix[ydx], yvlmi[ydx], yvlrg[ydx], yvlpk[ydx],
           &zvrtb[tdx],&zvrix[tdx],&zvrmi[tdx],&zvrrg[tdx],&zvrpk[tdx],
           &zvltb[tdx],&zvlix[tdx],&zvlmi[tdx],&zvlrg[tdx],&zvlpk[tdx]);
   __syncthreads();

   for(int i=1; i<dim; i++) // zv[tdx] = zv[tdx] + xv[i]*yv[dim+tdx-i];
   {
      ydx = dim + tdx - i;

      dag_mul( xvrtb[i],  xvrix[i],  xvrmi[i],  xvrrg[i],  xvrpk[i],
               xvltb[i],  xvlix[i],  xvlmi[i],  xvlrg[i],  xvlpk[i],
               yvrtb[ydx],yvrix[ydx],yvrmi[ydx],yvrrg[ydx],yvrpk[ydx],
               yvltb[ydx],yvlix[ydx],yvlmi[ydx],yvlrg[ydx],yvlpk[ydx],
              &prdrtb,  &prdrix,   &prdrmi,   &prdrrg,   &prdrpk,
              &prdltb,  &prdlix,   &prdlmi,   &prdlrg,   &prdlpk);
      __syncthreads();

      dag_inc(&zvrtb[tdx],&zvrix[tdx],&zvrmi[tdx],&zvrrg[tdx],&zvrpk[tdx],
              &zvltb[tdx],&zvlix[tdx],&zvlmi[tdx],&zvlrg[tdx],&zvlpk[tdx],
              prdrtb,     prdrix,     prdrmi,     prdrrg,     prdrpk,
              prdltb,     prdlix,     prdlmi,     prdlrg,     prdlpk);
      __syncthreads();
   }
   __syncthreads();

   datartb[idx3] = zvrtb[tdx]; // storing the output
   datarix[idx3] = zvrix[tdx];
   datarmi[idx3] = zvrmi[tdx];
   datarrg[idx3] = zvrrg[tdx];
   datarpk[idx3] = zvrpk[tdx];
   dataltb[idx3] = zvltb[tdx];
   datalix[idx3] = zvlix[tdx];
   datalmi[idx3] = zvlmi[tdx];
   datalrg[idx3] = zvlrg[tdx];
   datalpk[idx3] = zvlpk[tdx];
}

__global__ void dbl10_update_addjobs
 ( double *datartb, double *datarix, double *datarmi,
   double *datarrg, double *datarpk,
   double *dataltb, double *datalix, double *datalmi,
   double *datalrg, double *datalpk,
   int *in1idx, int *in2idx, int *outidx, int dim )
{
   const int bdx = blockIdx.x;           // index to the convolution job
   const int tdx = threadIdx.x;          // index to the output of the job
   const int idx1 = in1idx[bdx] + tdx;
   const int idx2 = in2idx[bdx] + tdx;
   const int idx3 = outidx[bdx] + tdx;

   __shared__ double xvrtb[da_shmemsize];
   __shared__ double xvrix[da_shmemsize];
   __shared__ double xvrmi[da_shmemsize];
   __shared__ double xvrrg[da_shmemsize];
   __shared__ double xvrpk[da_shmemsize];
   __shared__ double xvltb[da_shmemsize];
   __shared__ double xvlix[da_shmemsize];
   __shared__ double xvlmi[da_shmemsize];
   __shared__ double xvlrg[da_shmemsize];
   __shared__ double xvlpk[da_shmemsize];
   __shared__ double yvrtb[da_shmemsize];
   __shared__ double yvrix[da_shmemsize];
   __shared__ double yvrmi[da_shmemsize];
   __shared__ double yvrrg[da_shmemsize];
   __shared__ double yvrpk[da_shmemsize];
   __shared__ double yvltb[da_shmemsize];
   __shared__ double yvlix[da_shmemsize];
   __shared__ double yvlmi[da_shmemsize];
   __shared__ double yvlrg[da_shmemsize];
   __shared__ double yvlpk[da_shmemsize];
   __shared__ double zvrtb[da_shmemsize];
   __shared__ double zvrix[da_shmemsize];
   __shared__ double zvrmi[da_shmemsize];
   __shared__ double zvrrg[da_shmemsize];
   __shared__ double zvrpk[da_shmemsize];
   __shared__ double zvltb[da_shmemsize];
   __shared__ double zvlix[da_shmemsize];
   __shared__ double zvlmi[da_shmemsize];
   __shared__ double zvlrg[da_shmemsize];
   __shared__ double zvlpk[da_shmemsize];

   xvrtb[tdx] = datartb[idx1];  // loading first input
   xvrix[tdx] = datarix[idx1];
   xvrmi[tdx] = datarmi[idx1];
   xvrrg[tdx] = datarrg[idx1];
   xvrpk[tdx] = datarpk[idx1];
   xvltb[tdx] = dataltb[idx1]; 
   xvlix[tdx] = datalix[idx1];
   xvlmi[tdx] = datalmi[idx1];
   xvlrg[tdx] = datalrg[idx1];
   xvlpk[tdx] = datalpk[idx1];
   yvrtb[tdx] = datartb[idx2];  // loading second input
   yvrix[tdx] = datarix[idx2];
   yvrmi[tdx] = datarmi[idx2];
   yvrrg[tdx] = datarrg[idx2];
   yvrpk[tdx] = datarpk[idx2];
   yvltb[tdx] = dataltb[idx2];
   yvlix[tdx] = datalix[idx2];
   yvlmi[tdx] = datalmi[idx2];
   yvlrg[tdx] = datalrg[idx2];
   yvlpk[tdx] = datalpk[idx2];

   // zv[tdx] = xv[tdx] + yv[tdx];

   __syncthreads();

   dag_add( xvrtb[tdx], xvrix[tdx], xvrmi[tdx], xvrrg[tdx], xvrpk[tdx],
            xvltb[tdx], xvlix[tdx], xvlmi[tdx], xvlrg[tdx], xvlpk[tdx],
            yvrtb[tdx], yvrix[tdx], yvrmi[tdx], yvrrg[tdx], yvrpk[tdx],
            yvltb[tdx], yvlix[tdx], yvlmi[tdx], yvlrg[tdx], yvlpk[tdx],
           &zvrtb[tdx],&zvrix[tdx],&zvrmi[tdx],&zvrrg[tdx],&zvrpk[tdx],
           &zvltb[tdx],&zvlix[tdx],&zvlmi[tdx],&zvlrg[tdx],&zvlpk[tdx]);

   __syncthreads();

   datartb[idx3] = zvrtb[tdx]; // storing the output
   datarix[idx3] = zvrix[tdx];
   datarmi[idx3] = zvrmi[tdx];
   datarrg[idx3] = zvrrg[tdx];
   datarpk[idx3] = zvrpk[tdx];
   dataltb[idx3] = zvltb[tdx];
   datalix[idx3] = zvlix[tdx];
   datalmi[idx3] = zvlmi[tdx];
   datalrg[idx3] = zvlrg[tdx];
   datalpk[idx3] = zvlpk[tdx];
}

void convoluted_data10_to_output
 ( double *datartb, double *datarix, double *datarmi,
   double *datarrg, double *datarpk,
   double *dataltb, double *datalix, double *datalmi,
   double *datalrg, double *datalpk,
   double **outputrtb, double **outputrix, double **outputrmi,
   double **outputrrg, double **outputrpk,
   double **outputltb, double **outputlix, double **outputlmi,
   double **outputlrg, double **outputlpk,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart, bool verbose )
{
   const int deg1 = deg+1;
   int ix0,ix1,ix2;

   for(int i=0; i<=deg; i++) // output[dim][i] = data[i];
   {
      outputrtb[dim][i] = datartb[i];
      outputrix[dim][i] = datarix[i];
      outputrmi[dim][i] = datarmi[i];
      outputrrg[dim][i] = datarrg[i];
      outputrpk[dim][i] = datarpk[i];
      outputltb[dim][i] = dataltb[i];
      outputlix[dim][i] = datalix[i];
      outputlmi[dim][i] = datalmi[i];
      outputlrg[dim][i] = datalrg[i];
      outputlpk[dim][i] = datalpk[i];
   }
   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++) // output[i][j] = 0.0;
      {
         outputrtb[i][j] = 0.0;
         outputrix[i][j] = 0.0;
         outputrmi[i][j] = 0.0;
         outputrrg[i][j] = 0.0;
         outputrpk[i][j] = 0.0;
         outputltb[i][j] = 0.0;
         outputlix[i][j] = 0.0;
         outputlmi[i][j] = 0.0;
         outputlrg[i][j] = 0.0;
         outputlpk[i][j] = 0.0;
      }

   for(int k=0; k<nbr; k++)
   {
      ix1 = fstart[k] + (nvr[k]-1)*deg1;
      
      if(verbose)
         cout << "monomial " << k << " update starts at " << ix1 << endl;

      for(int i=0; i<=deg; i++) // output[dim][i] += data[ix1++];
         daf_inc(&outputrtb[dim][i],&outputrix[dim][i],&outputrmi[dim][i],
                 &outputrrg[dim][i],&outputrpk[dim][i],
                 &outputltb[dim][i],&outputlix[dim][i],&outputlmi[dim][i],
                 &outputlrg[dim][i],&outputlpk[dim][i],
                    datartb[ix1],      datarix[ix1],      datarmi[ix1],
                    datarrg[ix1],      datarpk[ix1],
                    dataltb[ix1],      datalix[ix1],      datalmi[ix1],
                    datalrg[ix1],      datalpk[ix1++]);
     
      ix0 = idx[k][0];
      if(nvr[k] == 1)
      {
         ix1 = (1 + k)*deg1;
            
         for(int i=0; i<=deg; i++) // output[ix0][i] += data[ix1++];
            daf_inc(&outputrtb[ix0][i],&outputrix[ix0][i],&outputrmi[ix0][i],
                    &outputrrg[ix0][i],&outputrpk[ix0][i],
                    &outputltb[ix0][i],&outputlix[ix0][i],&outputlmi[ix0][i],
                    &outputlrg[ix0][i],&outputlpk[ix0][i],
                       datartb[ix1],      datarix[ix1],      datarmi[ix1],
                       datarrg[ix1],      datarpk[ix1],
                       dataltb[ix1],      datalix[ix1],      datalmi[ix1],
                       datalrg[ix1],      datalpk[ix1++]);
      }
      else
      {                               // update first and last derivative
         ix2 = nvr[k]-3;
         if(ix2 < 0) ix2 = 0;
         ix1 = bstart[k] + ix2*deg1;

         for(int i=0; i<=deg; i++) // output[ix0][i] += data[ix1++];
            daf_inc(&outputrtb[ix0][i],&outputrix[ix0][i],&outputrmi[ix0][i],
                    &outputrrg[ix0][i],&outputrpk[ix0][i],
                    &outputltb[ix0][i],&outputlix[ix0][i],&outputlmi[ix0][i],
                    &outputlrg[ix0][i],&outputlpk[ix0][i],
                       datartb[ix1],      datarix[ix1],      datarmi[ix1],
                       datarrg[ix1],      datarpk[ix1],
                       dataltb[ix1],      datalix[ix1],      datalmi[ix1],
                       datalrg[ix1],      datalpk[ix1++]);

         ix2 = nvr[k]-2;
         ix1 = fstart[k] + ix2*deg1;
         ix0 = idx[k][ix2+1];

         for(int i=0; i<=deg; i++) // output[ix0][i] += data[ix1++];
            daf_inc(&outputrtb[ix0][i],&outputrix[ix0][i],&outputrmi[ix0][i],
                    &outputrrg[ix0][i],&outputrpk[ix0][i],
                    &outputltb[ix0][i],&outputlix[ix0][i],&outputlmi[ix0][i],
                    &outputlrg[ix0][i],&outputlpk[ix0][i],
                       datartb[ix1],      datarix[ix1],      datarmi[ix1],
                       datarrg[ix1],      datarpk[ix1],
                       dataltb[ix1],      datalix[ix1],      datalmi[ix1],
                       datalrg[ix1],      datalpk[ix1++]);
 
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
                  daf_inc(&outputrtb[ix0][i],&outputrix[ix0][i],
                          &outputrmi[ix0][i],
                          &outputrrg[ix0][i],&outputrpk[ix0][i],
                          &outputltb[ix0][i],&outputlix[ix0][i],
                          &outputlmi[ix0][i],
                          &outputlrg[ix0][i],&outputlpk[ix0][i],
                             datartb[ix1],      datarix[ix1],
                             datarmi[ix1],
                             datarrg[ix1],      datarpk[ix1],
                             dataltb[ix1],      datalix[ix1],
                             datalmi[ix1],
                             datalrg[ix1],      datalpk[ix1++]);
            }
         }
      }
   }
}

void added_data10_to_output
 ( double *datartb, double *datarix, double *datarmi,
   double *datarrg, double *datarpk,
   double *dataltb, double *datalix, double *datalmi,
   double *datalrg, double *datalpk,
   double **outputrtb, double **outputrix, double **outputrmi,
   double **outputrrg, double **outputrpk,
   double **outputltb, double **outputlix, double **outputlmi,
   double **outputlrg, double **outputlpk,
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
      outputrtb[dim][i] = datartb[ix];
      outputrix[dim][i] = datarix[ix];
      outputrmi[dim][i] = datarmi[ix];
      outputrrg[dim][i] = datarrg[ix];
      outputrpk[dim][i] = datarpk[ix];
      outputltb[dim][i] = dataltb[ix];
      outputlix[dim][i] = datalix[ix];
      outputlmi[dim][i] = datalmi[ix];
      outputlrg[dim][i] = datalrg[ix];
      outputlpk[dim][i] = datalpk[ix++];
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
            outputrtb[0][i] = 0.0;
            outputrix[0][i] = 0.0;
            outputrmi[0][i] = 0.0;
            outputrrg[0][i] = 0.0;
            outputrpk[0][i] = 0.0;
            outputltb[0][i] = 0.0;
            outputlix[0][i] = 0.0;
            outputlmi[0][i] = 0.0;
            outputlrg[0][i] = 0.0;
            outputlpk[0][i] = 0.0;
         }
      }
      else
      {
         int cffidx = (1 + difidx)*deg1;

         if(verbose)
            cout << "updating derivative with coefficient ..." << endl;

         for(int i=0; i<=deg; i++)
         {
            outputrtb[0][i] = datartb[cffidx];
            outputrix[0][i] = datarix[cffidx];
            outputrmi[0][i] = datarmi[cffidx];
            outputrrg[0][i] = datarrg[cffidx];
            outputrpk[0][i] = datarpk[cffidx];
            outputltb[0][i] = dataltb[cffidx];
            outputlix[0][i] = datalix[cffidx];
            outputlmi[0][i] = datalmi[cffidx];
            outputlrg[0][i] = datalrg[cffidx];
            outputlpk[0][i] = datalpk[cffidx++];
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
         outputrtb[0][i] = datartb[ix];
         outputrix[0][i] = datarix[ix];
         outputrmi[0][i] = datarmi[ix];
         outputrrg[0][i] = datarrg[ix];
         outputrpk[0][i] = datarpk[ix];
         outputltb[0][i] = dataltb[ix];
         outputlix[0][i] = datalix[ix];
         outputlmi[0][i] = datalmi[ix];
         outputlrg[0][i] = datalrg[ix];
         outputlpk[0][i] = datalpk[ix++];
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
               outputrtb[k][i] = 0.0;
               outputrix[k][i] = 0.0;
               outputrmi[k][i] = 0.0;
               outputrrg[k][i] = 0.0;
               outputrpk[k][i] = 0.0;
               outputltb[k][i] = 0.0;
               outputlix[k][i] = 0.0;
               outputlmi[k][i] = 0.0;
               outputlrg[k][i] = 0.0;
               outputlpk[k][i] = 0.0;
            }
         }
         else
         {
            int cffidx = (1 + difidx)*deg1;

            if(verbose)
               cout << "updating derivative with coefficient ..." << endl;

            for(int i=0; i<=deg; i++)
            {
               outputrtb[k][i] = datartb[cffidx];
               outputrix[k][i] = datarix[cffidx];
               outputrmi[k][i] = datarmi[cffidx];
               outputrrg[k][i] = datarrg[cffidx];
               outputrpk[k][i] = datarpk[cffidx];
               outputltb[k][i] = dataltb[cffidx];
               outputlix[k][i] = datalix[cffidx];
               outputlmi[k][i] = datalmi[cffidx];
               outputlrg[k][i] = datalrg[cffidx];
               outputlpk[k][i] = datalpk[cffidx++];
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
               outputrtb[k][i] = datartb[ix];
               outputrix[k][i] = datarix[ix];
               outputrmi[k][i] = datarmi[ix];
               outputrrg[k][i] = datarrg[ix];
               outputrpk[k][i] = datarpk[ix];
               outputltb[k][i] = dataltb[ix];
               outputlix[k][i] = datalix[ix];
               outputlmi[k][i] = datalmi[ix];
               outputlrg[k][i] = datalrg[ix];
               outputlpk[k][i] = datalpk[ix++];
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
               outputrtb[k][i] = datartb[ix];
               outputrix[k][i] = datarix[ix];
               outputrmi[k][i] = datarmi[ix];
               outputrrg[k][i] = datarrg[ix];
               outputrpk[k][i] = datarpk[ix];
               outputltb[k][i] = dataltb[ix];
               outputlix[k][i] = datalix[ix];
               outputlmi[k][i] = datalmi[ix];
               outputlrg[k][i] = datalrg[ix];
               outputlpk[k][i] = datalpk[ix++];
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
               outputrtb[k][i] = datartb[ix];
               outputrix[k][i] = datarix[ix];
               outputrmi[k][i] = datarmi[ix];
               outputrrg[k][i] = datarrg[ix];
               outputrpk[k][i] = datarpk[ix];
               outputltb[k][i] = dataltb[ix];
               outputlix[k][i] = datalix[ix];
               outputlmi[k][i] = datalmi[ix];
               outputlrg[k][i] = datalrg[ix];
               outputlpk[k][i] = datalpk[ix++];
            }
         }
      }
   }
}

void GPU_dbl10_poly_evaldiff
 ( int BS, int dim, int nbr, int deg, int *nvr, int **idx,
   double *cstrtb, double *cstrix, double *cstrmi,
   double *cstrrg, double *cstrpk,
   double *cstltb, double *cstlix, double *cstlmi,
   double *cstlrg, double *cstlpk,
   double **cffrtb, double **cffrix, double **cffrmi,
   double **cffrrg, double **cffrpk,
   double **cffltb, double **cfflix, double **cfflmi,
   double **cfflrg, double **cfflpk,
   double **inputrtb, double **inputrix, double **inputrmi,
   double **inputrrg, double **inputrpk,
   double **inputltb, double **inputlix, double **inputlmi,
   double **inputlrg, double **inputlpk,
   double **outputrtb, double **outputrix, double **outputrmi,
   double **outputrrg, double **outputrpk,
   double **outputltb, double **outputlix, double **outputlmi,
   double **outputlrg, double **outputlpk,
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
   double *datartb_h = new double[totalcff];        // data on host
   double *datarix_h = new double[totalcff];
   double *datarmi_h = new double[totalcff];
   double *datarrg_h = new double[totalcff];
   double *datarpk_h = new double[totalcff];
   double *dataltb_h = new double[totalcff];
   double *datalix_h = new double[totalcff];
   double *datalmi_h = new double[totalcff];
   double *datalrg_h = new double[totalcff];
   double *datalpk_h = new double[totalcff];
   int ix = 0;
   for(int i=0; i<deg1; i++)
   {
      datartb_h[ix] = cstrtb[i];
      datarix_h[ix] = cstrix[i];
      datarmi_h[ix] = cstrmi[i];
      datarrg_h[ix] = cstrrg[i];
      datarpk_h[ix] = cstrpk[i];
      dataltb_h[ix] = cstltb[i];
      datalix_h[ix] = cstlix[i];
      datalmi_h[ix] = cstlmi[i];
      datalrg_h[ix] = cstlrg[i];
      datalpk_h[ix++] = cstlpk[i];
   }
   for(int i=0; i<nbr; i++)
      for(int j=0; j<deg1; j++)
      {
         datartb_h[ix] = cffrtb[i][j];
         datarix_h[ix] = cffrix[i][j];
         datarmi_h[ix] = cffrmi[i][j];
         datarrg_h[ix] = cffrrg[i][j];
         datarpk_h[ix] = cffrpk[i][j];
         dataltb_h[ix] = cffltb[i][j];
         datalix_h[ix] = cfflix[i][j];
         datalmi_h[ix] = cfflmi[i][j];
         datalrg_h[ix] = cfflrg[i][j];
         datalpk_h[ix++] = cfflpk[i][j];
      }
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++)
      {
         datartb_h[ix] = inputrtb[i][j];
         datarix_h[ix] = inputrix[i][j];
         datarmi_h[ix] = inputrmi[i][j];
         datarrg_h[ix] = inputrrg[i][j];
         datarpk_h[ix] = inputrpk[i][j];
         dataltb_h[ix] = inputltb[i][j];
         datalix_h[ix] = inputlix[i][j];
         datalmi_h[ix] = inputlmi[i][j];
         datalrg_h[ix] = inputlrg[i][j];
         datalpk_h[ix++] = inputlpk[i][j];
      }

   double *datartb_d;                               // device data
   double *datarix_d;
   double *datarmi_d;
   double *datarrg_d;
   double *datarpk_d;
   double *dataltb_d;
   double *datalix_d;
   double *datalmi_d;
   double *datalrg_d;
   double *datalpk_d;
   const size_t szdata = totalcff*sizeof(double);
   cudaMalloc((void**)&datartb_d,szdata);
   cudaMalloc((void**)&datarix_d,szdata);
   cudaMalloc((void**)&datarmi_d,szdata);
   cudaMalloc((void**)&datarrg_d,szdata);
   cudaMalloc((void**)&datarpk_d,szdata);
   cudaMalloc((void**)&dataltb_d,szdata);
   cudaMalloc((void**)&datalix_d,szdata);
   cudaMalloc((void**)&datalmi_d,szdata);
   cudaMalloc((void**)&datalrg_d,szdata);
   cudaMalloc((void**)&datalpk_d,szdata);
   cudaMemcpy(datartb_d,datartb_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datarix_d,datarix_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datarmi_d,datarmi_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datarrg_d,datarrg_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datarpk_d,datarpk_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(dataltb_d,dataltb_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datalix_d,datalix_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datalmi_d,datalmi_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datalrg_d,datalrg_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datalpk_d,datalpk_h,szdata,cudaMemcpyHostToDevice);

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
         dbl10_padded_convjobs<<<jobnbr,BS>>>
            (datartb_d,datarix_d,datarmi_d,datarrg_d,datarpk_d,
             dataltb_d,datalix_d,datalmi_d,datalrg_d,datalpk_d,
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
         dbl10_update_addjobs<<<jobnbr,BS>>>
            (datartb_d,datarix_d,datarmi_d,datarrg_d,datarpk_d,
             dataltb_d,datalix_d,datalmi_d,datalrg_d,datalpk_d,
             in1ix_d,in2ix_d,outix_d,deg1);
         cudaEventRecord(stop);
         cudaEventSynchronize(stop);
         cudaEventElapsedTime(&milliseconds,start,stop);
         *addlapms += milliseconds;
      }
      free(in1ix_h); free(in2ix_h); free(outix_h);
   }
   gettimeofday(&endtime,0);
   cudaMemcpy(datartb_h,datartb_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datarix_h,datarix_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datarmi_h,datarmi_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datarrg_h,datarrg_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datarpk_h,datarpk_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(dataltb_h,dataltb_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datalix_h,datalix_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datalmi_h,datalmi_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datalrg_h,datalrg_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datalpk_h,datalpk_d,szdata,cudaMemcpyDeviceToHost);
   *elapsedms = *cnvlapms + *addlapms;
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   *walltimesec = seconds + microseconds*1.0e-6;

   // convoluted_data2_to_output
   //    (data_h,output,dim,nbr,deg,nvr,idx,fstart,bstart,cstart,verbose);
   added_data10_to_output
      (datartb_h,datarix_h,datarmi_h,datarrg_h,datarpk_h,
       dataltb_h,datalix_h,datalmi_h,datalrg_h,datalpk_h,
       outputrtb,outputrix,outputrmi,outputrrg,outputrpk,
       outputltb,outputlix,outputlmi,outputlrg,outputlpk,
       dim,nbr,deg,nvr,idx,fstart,bstart,cstart,addjobs,verbose);
   if(verbose)
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
}
