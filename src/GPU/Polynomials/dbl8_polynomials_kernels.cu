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

// The constant od_shmemsize is the bound on the shared memory size.

#define od_shmemsize 192

using namespace std;

__global__ void dbl8_padded_convjobs
 ( double *datahihihi, double *datahilohi,
   double *datahihilo, double *datahilolo,
   double *datalohihi, double *datalolohi,
   double *datalohilo, double *datalololo,
   int *in1idx, int *in2idx, int *outidx, int dim )
{
   const int bdx = blockIdx.x;           // index to the convolution job
   const int tdx = threadIdx.x;          // index to the output of the job
   const int idx1 = in1idx[bdx] + tdx;
   const int idx2 = in2idx[bdx] + tdx;
   const int idx3 = outidx[bdx] + tdx;

   __shared__ double xvhihihi[od_shmemsize];
   __shared__ double xvhilohi[od_shmemsize];
   __shared__ double xvhihilo[od_shmemsize];
   __shared__ double xvhilolo[od_shmemsize];
   __shared__ double xvlohihi[od_shmemsize];
   __shared__ double xvlolohi[od_shmemsize];
   __shared__ double xvlohilo[od_shmemsize];
   __shared__ double xvlololo[od_shmemsize];
   __shared__ double yvhihihi[2*od_shmemsize];
   __shared__ double yvhilohi[2*od_shmemsize];
   __shared__ double yvhihilo[2*od_shmemsize];
   __shared__ double yvhilolo[2*od_shmemsize];
   __shared__ double yvlohihi[2*od_shmemsize];
   __shared__ double yvlolohi[2*od_shmemsize];
   __shared__ double yvlohilo[2*od_shmemsize];
   __shared__ double yvlololo[2*od_shmemsize];
   __shared__ double zvhihihi[od_shmemsize];
   __shared__ double zvhilohi[od_shmemsize];
   __shared__ double zvhihilo[od_shmemsize];
   __shared__ double zvhilolo[od_shmemsize];
   __shared__ double zvlohihi[od_shmemsize];
   __shared__ double zvlolohi[od_shmemsize];
   __shared__ double zvlohilo[od_shmemsize];
   __shared__ double zvlololo[od_shmemsize];

   double prdhihihi,prdhilohi,prdhihilo,prdhilolo;
   double prdlohihi,prdlolohi,prdlohilo,prdlololo;
   int ydx = dim + tdx;

   xvhihihi[tdx] = datahihihi[idx1];  // loading first input
   xvhilohi[tdx] = datahilohi[idx1]; 
   xvhihilo[tdx] = datahihilo[idx1]; 
   xvhilolo[tdx] = datahilolo[idx1]; 
   xvlohihi[tdx] = datalohihi[idx1];
   xvlolohi[tdx] = datalolohi[idx1]; 
   xvlohilo[tdx] = datalohilo[idx1]; 
   xvlololo[tdx] = datalololo[idx1]; 
   yvhihihi[tdx] = 0.0;             // padded with zeros
   yvhilohi[tdx] = 0.0;
   yvhihilo[tdx] = 0.0;
   yvhilolo[tdx] = 0.0;
   yvlohihi[tdx] = 0.0;
   yvlolohi[tdx] = 0.0;
   yvlohilo[tdx] = 0.0;
   yvlololo[tdx] = 0.0;
   yvhihihi[ydx] = datahihihi[idx2];  // loading second input
   yvhilohi[ydx] = datahilohi[idx2];
   yvhihilo[ydx] = datahihilo[idx2];
   yvhilolo[ydx] = datahilolo[idx2];
   yvlohihi[ydx] = datalohihi[idx2];
   yvlolohi[ydx] = datalolohi[idx2];
   yvlohilo[ydx] = datalohilo[idx2];
   yvlololo[ydx] = datalololo[idx2];

   __syncthreads();

   // zv[tdx] = xv[0]*yv[tdx];
   odg_mul( xvhihihi[0],   xvhilohi[0],   xvhihilo[0],   xvhilolo[0],
            xvlohihi[0],   xvlolohi[0],   xvlohilo[0],   xvlololo[0],
            yvhihihi[ydx], yvhilohi[ydx], yvhihilo[ydx], yvhilolo[ydx],
            yvlohihi[ydx], yvlolohi[ydx], yvlohilo[ydx], yvlololo[ydx],
           &zvhihihi[tdx],&zvhilohi[tdx],&zvhihilo[tdx],&zvhilolo[tdx],
           &zvlohihi[tdx],&zvlolohi[tdx],&zvlohilo[tdx],&zvlololo[tdx]);
   __syncthreads();

   for(int i=1; i<dim; i++) // zv[tdx] = zv[tdx] + xv[i]*yv[dim+tdx-i];
   {
      ydx = dim + tdx - i;

      odg_mul( xvhihihi[i],  xvhilohi[i],  xvhihilo[i],  xvhilolo[i],
               xvlohihi[i],  xvlolohi[i],  xvlohilo[i],  xvlololo[i],
               yvhihihi[ydx],yvhilohi[ydx],yvhihilo[ydx],yvhilolo[ydx],
               yvlohihi[ydx],yvlolohi[ydx],yvlohilo[ydx],yvlololo[ydx],
              &prdhihihi,  &prdhilohi,   &prdhihilo,   &prdhilolo,
              &prdlohihi,  &prdlolohi,   &prdlohilo,   &prdlololo);
      __syncthreads();

      odg_inc(&zvhihihi[tdx],&zvhilohi[tdx],&zvhihilo[tdx],&zvhilolo[tdx],
              &zvlohihi[tdx],&zvlolohi[tdx],&zvlohilo[tdx],&zvlololo[tdx],
              prdhihihi,     prdhilohi,     prdhihilo,     prdhilolo,
              prdlohihi,     prdlolohi,     prdlohilo,     prdlololo);
      __syncthreads();
   }
   __syncthreads();

   datahihihi[idx3] = zvhihihi[tdx]; // storing the output
   datahilohi[idx3] = zvhilohi[tdx];
   datahihilo[idx3] = zvhihilo[tdx];
   datahilolo[idx3] = zvhilolo[tdx];
   datalohihi[idx3] = zvlohihi[tdx];
   datalolohi[idx3] = zvlolohi[tdx];
   datalohilo[idx3] = zvlohilo[tdx];
   datalololo[idx3] = zvlololo[tdx];
}

__global__ void dbl8_update_addjobs
 ( double *datahihihi, double *datahilohi,
   double *datahihilo, double *datahilolo,
   double *datalohihi, double *datalolohi,
   double *datalohilo, double *datalololo,
   int *in1idx, int *in2idx, int *outidx, int dim )
{
   const int bdx = blockIdx.x;           // index to the convolution job
   const int tdx = threadIdx.x;          // index to the output of the job
   const int idx1 = in1idx[bdx] + tdx;
   const int idx2 = in2idx[bdx] + tdx;
   const int idx3 = outidx[bdx] + tdx;

   __shared__ double xvhihihi[od_shmemsize];
   __shared__ double xvhilohi[od_shmemsize];
   __shared__ double xvhihilo[od_shmemsize];
   __shared__ double xvhilolo[od_shmemsize];
   __shared__ double xvlohihi[od_shmemsize];
   __shared__ double xvlolohi[od_shmemsize];
   __shared__ double xvlohilo[od_shmemsize];
   __shared__ double xvlololo[od_shmemsize];
   __shared__ double yvhihihi[od_shmemsize];
   __shared__ double yvhilohi[od_shmemsize];
   __shared__ double yvhihilo[od_shmemsize];
   __shared__ double yvhilolo[od_shmemsize];
   __shared__ double yvlohihi[od_shmemsize];
   __shared__ double yvlolohi[od_shmemsize];
   __shared__ double yvlohilo[od_shmemsize];
   __shared__ double yvlololo[od_shmemsize];
   __shared__ double zvhihihi[od_shmemsize];
   __shared__ double zvhilohi[od_shmemsize];
   __shared__ double zvhihilo[od_shmemsize];
   __shared__ double zvhilolo[od_shmemsize];
   __shared__ double zvlohihi[od_shmemsize];
   __shared__ double zvlolohi[od_shmemsize];
   __shared__ double zvlohilo[od_shmemsize];
   __shared__ double zvlololo[od_shmemsize];

   xvhihihi[tdx] = datahihihi[idx1];  // loading first input
   xvhilohi[tdx] = datahilohi[idx1];
   xvhihilo[tdx] = datahihilo[idx1];
   xvhilolo[tdx] = datahilolo[idx1];
   xvlohihi[tdx] = datalohihi[idx1];
   xvlolohi[tdx] = datalolohi[idx1];
   xvlohilo[tdx] = datalohilo[idx1];
   xvlololo[tdx] = datalololo[idx1];
   yvhihihi[tdx] = datahihihi[idx2];  // loading second input
   yvhilohi[tdx] = datahilohi[idx2];
   yvhihilo[tdx] = datahihilo[idx2];
   yvhilolo[tdx] = datahilolo[idx2];
   yvlohihi[tdx] = datalohihi[idx2];
   yvlolohi[tdx] = datalolohi[idx2];
   yvlohilo[tdx] = datalohilo[idx2];
   yvlololo[tdx] = datalololo[idx2];

   // zv[tdx] = xv[tdx] + yv[tdx];

   __syncthreads();

   odg_add( xvhihihi[tdx], xvhilohi[tdx], xvhihilo[tdx], xvhilolo[tdx],
            xvlohihi[tdx], xvlolohi[tdx], xvlohilo[tdx], xvlololo[tdx],
            yvhihihi[tdx], yvhilohi[tdx], yvhihilo[tdx], yvhilolo[tdx],
            yvlohihi[tdx], yvlolohi[tdx], yvlohilo[tdx], yvlololo[tdx],
           &zvhihihi[tdx],&zvhilohi[tdx],&zvhihilo[tdx],&zvhilolo[tdx],
           &zvlohihi[tdx],&zvlolohi[tdx],&zvlohilo[tdx],&zvlololo[tdx]);

   __syncthreads();

   datahihihi[idx3] = zvhihihi[tdx]; // storing the output
   datahilohi[idx3] = zvhilohi[tdx];
   datahihilo[idx3] = zvhihilo[tdx];
   datahilolo[idx3] = zvhilolo[tdx];
   datalohihi[idx3] = zvlohihi[tdx];
   datalolohi[idx3] = zvlolohi[tdx];
   datalohilo[idx3] = zvlohilo[tdx];
   datalololo[idx3] = zvlololo[tdx];
}

void convoluted_data8_to_output
 ( double *datahihihi, double *datahilohi,
   double *datahihilo, double *datahilolo,
   double *datalohihi, double *datalolohi,
   double *datalohilo, double *datalololo,
   double **outputhihihi, double **outputhilohi,
   double **outputhihilo, double **outputhilolo,
   double **outputlohihi, double **outputlolohi,
   double **outputlohilo, double **outputlololo,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart, bool verbose )
{
   const int deg1 = deg+1;
   int ix0,ix1,ix2;

   for(int i=0; i<=deg; i++) // output[dim][i] = data[i];
   {
      outputhihihi[dim][i] = datahihihi[i];
      outputhilohi[dim][i] = datahilohi[i];
      outputhihilo[dim][i] = datahihilo[i];
      outputhilolo[dim][i] = datahilolo[i];
      outputlohihi[dim][i] = datalohihi[i];
      outputlolohi[dim][i] = datalolohi[i];
      outputlohilo[dim][i] = datalohilo[i];
      outputlololo[dim][i] = datalololo[i];
   }
   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++) // output[i][j] = 0.0;
      {
         outputhihihi[i][j] = 0.0;
         outputhilohi[i][j] = 0.0;
         outputhihilo[i][j] = 0.0;
         outputhilolo[i][j] = 0.0;
         outputlohihi[i][j] = 0.0;
         outputlolohi[i][j] = 0.0;
         outputlohilo[i][j] = 0.0;
         outputlololo[i][j] = 0.0;
      }

   for(int k=0; k<nbr; k++)
   {
      ix1 = fstart[k] + (nvr[k]-1)*deg1;
      
      if(verbose)
         cout << "monomial " << k << " update starts at " << ix1 << endl;

      for(int i=0; i<=deg; i++) // output[dim][i] += data[ix1++];
         odf_inc(&outputhihihi[dim][i],&outputhilohi[dim][i],
                 &outputhihilo[dim][i],&outputhilolo[dim][i],
                 &outputlohihi[dim][i],&outputlolohi[dim][i],
                 &outputlohilo[dim][i],&outputlololo[dim][i],
                    datahihihi[ix1],      datahilohi[ix1],
                    datahihilo[ix1],      datahilolo[ix1],
                    datalohihi[ix1],      datalolohi[ix1],
                    datalohilo[ix1],      datalololo[ix1++]);
     
      ix0 = idx[k][0];
      if(nvr[k] == 1)
      {
         ix1 = (1 + k)*deg1;
            
         for(int i=0; i<=deg; i++) // output[ix0][i] += data[ix1++];
            odf_inc(&outputhihihi[ix0][i],&outputhilohi[ix0][i],
                    &outputhihilo[ix0][i],&outputhilolo[ix0][i],
                    &outputlohihi[ix0][i],&outputlolohi[ix0][i],
                    &outputlohilo[ix0][i],&outputlololo[ix0][i],
                       datahihihi[ix1],      datahilohi[ix1],
                       datahihilo[ix1],      datahilolo[ix1],
                       datalohihi[ix1],      datalolohi[ix1],
                       datalohilo[ix1],      datalololo[ix1++]);
      }
      else
      {                               // update first and last derivative
         ix2 = nvr[k]-3;
         if(ix2 < 0) ix2 = 0;
         ix1 = bstart[k] + ix2*deg1;

         for(int i=0; i<=deg; i++) // output[ix0][i] += data[ix1++];
            odf_inc(&outputhihihi[ix0][i],&outputhilohi[ix0][i],
                    &outputhihilo[ix0][i],&outputhilolo[ix0][i],
                    &outputlohihi[ix0][i],&outputlolohi[ix0][i],
                    &outputlohilo[ix0][i],&outputlololo[ix0][i],
                       datahihihi[ix1],      datahilohi[ix1],
                       datahihilo[ix1],      datahilolo[ix1],
                       datalohihi[ix1],      datalolohi[ix1],
                       datalohilo[ix1],      datalololo[ix1++]);

         ix2 = nvr[k]-2;
         ix1 = fstart[k] + ix2*deg1;
         ix0 = idx[k][ix2+1];

         for(int i=0; i<=deg; i++) // output[ix0][i] += data[ix1++];
            odf_inc(&outputhihihi[ix0][i],&outputhilohi[ix0][i],
                    &outputhihilo[ix0][i],&outputhilolo[ix0][i],
                    &outputlohihi[ix0][i],&outputlolohi[ix0][i],
                    &outputlohilo[ix0][i],&outputlololo[ix0][i],
                       datahihihi[ix1],      datahilohi[ix1],
                       datahihilo[ix1],      datahilolo[ix1],
                       datalohihi[ix1],      datalolohi[ix1],
                       datalohilo[ix1],      datalololo[ix1++]);
 
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
                  odf_inc(&outputhihihi[ix0][i],&outputhilohi[ix0][i],
                          &outputhihilo[ix0][i],&outputhilolo[ix0][i],
                          &outputlohihi[ix0][i],&outputlolohi[ix0][i],
                          &outputlohilo[ix0][i],&outputlololo[ix0][i],
                             datahihihi[ix1],      datahilohi[ix1],
                             datahihilo[ix1],      datahilolo[ix1],
                             datalohihi[ix1],      datalolohi[ix1],
                             datalohilo[ix1],      datalololo[ix1++]);
            }
         }
      }
   }
}

void added_data8_to_output
 ( double *datahihihi, double *datahilohi,
   double *datahihilo, double *datahilolo,
   double *datalohihi, double *datalolohi,
   double *datalohilo, double *datalololo,
   double **outputhihihi, double **outputhilohi,
   double **outputhihilo, double **outputhilolo,
   double **outputlohihi, double **outputlolohi,
   double **outputlohilo, double **outputlololo,
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
      outputhihihi[dim][i] = datahihihi[ix];
      outputhilohi[dim][i] = datahilohi[ix];
      outputhihilo[dim][i] = datahihilo[ix];
      outputhilolo[dim][i] = datahilolo[ix];
      outputlohihi[dim][i] = datalohihi[ix];
      outputlolohi[dim][i] = datalolohi[ix];
      outputlohilo[dim][i] = datalohilo[ix];
      outputlololo[dim][i] = datalololo[ix++];
   }
   int cnt = jobs.get_differential_count(0);
   if(cnt == 0) // it could be there is no first variable anywhere ...
   {
      for(int i=0; i<=deg; i++) // output[0][i] = 0.0;
      {
         outputhihihi[0][i] = 0.0;
         outputhilohi[0][i] = 0.0;
         outputhihilo[0][i] = 0.0;
         outputhilolo[0][i] = 0.0;
         outputlohihi[0][i] = 0.0;
         outputlolohi[0][i] = 0.0;
         outputlohilo[0][i] = 0.0;
         outputlololo[0][i] = 0.0;
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
         outputhihihi[0][i] = datahihihi[ix];
         outputhilohi[0][i] = datahilohi[ix];
         outputhihilo[0][i] = datahihilo[ix];
         outputhilolo[0][i] = datahilolo[ix];
         outputlohihi[0][i] = datalohihi[ix];
         outputlolohi[0][i] = datalolohi[ix];
         outputlohilo[0][i] = datalohilo[ix];
         outputlololo[0][i] = datalololo[ix++];
      }
      for(int k=1; k<dim; k++) // updating all other derivatives
      {
         int cnt = jobs.get_differential_count(k);
         if(cnt == 0) // it could be there is no variable k anywhere ...
         {
            for(int i=0; i<=deg; i++) // output[k][i] = 0.0;
            {
               outputhihihi[k][i] = 0.0;
               outputhilohi[k][i] = 0.0;
               outputhihilo[k][i] = 0.0;
               outputhilolo[k][i] = 0.0;
               outputlohihi[k][i] = 0.0;
               outputlolohi[k][i] = 0.0;
               outputlohilo[k][i] = 0.0;
               outputlololo[k][i] = 0.0;
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
                  outputhihihi[k][i] = datahihihi[ix];
                  outputhilohi[k][i] = datahilohi[ix];
                  outputhihilo[k][i] = datahihilo[ix];
                  outputhilolo[k][i] = datahilolo[ix];
                  outputlohihi[k][i] = datalohihi[ix];
                  outputlolohi[k][i] = datalolohi[ix];
                  outputlohilo[k][i] = datalohilo[ix];
                  outputlololo[k][i] = datalololo[ix++];
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
                  outputhihihi[k][i] = datahihihi[ix];
                  outputhilohi[k][i] = datahilohi[ix];
                  outputhihilo[k][i] = datahihilo[ix];
                  outputhilolo[k][i] = datahilolo[ix];
                  outputlohihi[k][i] = datalohihi[ix];
                  outputlolohi[k][i] = datalolohi[ix];
                  outputlohilo[k][i] = datalohilo[ix];
                  outputlololo[k][i] = datalololo[ix++];
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
                  outputhihihi[k][i] = datahihihi[ix];
                  outputhilohi[k][i] = datahilohi[ix];
                  outputhihilo[k][i] = datahihilo[ix];
                  outputhilolo[k][i] = datahilolo[ix];
                  outputlohihi[k][i] = datalohihi[ix];
                  outputlolohi[k][i] = datalolohi[ix];
                  outputlohilo[k][i] = datalohilo[ix];
                  outputlololo[k][i] = datalololo[ix++];
               }
            }
         }
      }
   }
}

void GPU_dbl8_poly_evaldiff
 ( int BS, int dim, int nbr, int deg, int *nvr, int **idx,
   double *csthihihi, double *csthilohi,
   double *csthihilo, double *csthilolo,
   double *cstlohihi, double *cstlolohi,
   double *cstlohilo, double *cstlololo,
   double **cffhihihi, double **cffhilohi,
   double **cffhihilo, double **cffhilolo,
   double **cfflohihi, double **cfflolohi,
   double **cfflohilo, double **cfflololo,
   double **inputhihihi, double **inputhilohi,
   double **inputhihilo, double **inputhilolo,
   double **inputlohihi, double **inputlolohi,
   double **inputlohilo, double **inputlololo,
   double **outputhihihi, double **outputhilohi,
   double **outputhihilo, double **outputhilolo,
   double **outputlohihi, double **outputlolohi,
   double **outputlohilo, double **outputlololo,
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
   double *datahihihi_h = new double[totalcff];        // data on host
   double *datahilohi_h = new double[totalcff];
   double *datahihilo_h = new double[totalcff];
   double *datahilolo_h = new double[totalcff];
   double *datalohihi_h = new double[totalcff];
   double *datalolohi_h = new double[totalcff];
   double *datalohilo_h = new double[totalcff];
   double *datalololo_h = new double[totalcff];
   int ix = 0;
   for(int i=0; i<deg1; i++)
   {
      datahihihi_h[ix] = csthihihi[i];
      datahilohi_h[ix] = csthilohi[i];
      datahihilo_h[ix] = csthihilo[i];
      datahilolo_h[ix] = csthilolo[i];
      datalohihi_h[ix] = cstlohihi[i];
      datalolohi_h[ix] = cstlolohi[i];
      datalohilo_h[ix] = cstlohilo[i];
      datalololo_h[ix++] = cstlololo[i];
   }
   for(int i=0; i<nbr; i++)
      for(int j=0; j<deg1; j++)
      {
         datahihihi_h[ix] = cffhihihi[i][j];
         datahilohi_h[ix] = cffhilohi[i][j];
         datahihilo_h[ix] = cffhihilo[i][j];
         datahilolo_h[ix] = cffhilolo[i][j];
         datalohihi_h[ix] = cfflohihi[i][j];
         datalolohi_h[ix] = cfflolohi[i][j];
         datalohilo_h[ix] = cfflohilo[i][j];
         datalololo_h[ix++] = cfflololo[i][j];
      }
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++)
      {
         datahihihi_h[ix] = inputhihihi[i][j];
         datahilohi_h[ix] = inputhilohi[i][j];
         datahihilo_h[ix] = inputhihilo[i][j];
         datahilolo_h[ix] = inputhilolo[i][j];
         datalohihi_h[ix] = inputlohihi[i][j];
         datalolohi_h[ix] = inputlolohi[i][j];
         datalohilo_h[ix] = inputlohilo[i][j];
         datalololo_h[ix++] = inputlololo[i][j];
      }

   double *datahihihi_d;                               // device data
   double *datahilohi_d;
   double *datahihilo_d;
   double *datahilolo_d;
   double *datalohihi_d;
   double *datalolohi_d;
   double *datalohilo_d;
   double *datalololo_d;
   const size_t szdata = totalcff*sizeof(double);
   cudaMalloc((void**)&datahihihi_d,szdata);
   cudaMalloc((void**)&datahilohi_d,szdata);
   cudaMalloc((void**)&datahihilo_d,szdata);
   cudaMalloc((void**)&datahilolo_d,szdata);
   cudaMalloc((void**)&datalohihi_d,szdata);
   cudaMalloc((void**)&datalolohi_d,szdata);
   cudaMalloc((void**)&datalohilo_d,szdata);
   cudaMalloc((void**)&datalololo_d,szdata);
   cudaMemcpy(datahihihi_d,datahihihi_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datahilohi_d,datahilohi_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datahihilo_d,datahihilo_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datahilolo_d,datahilolo_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datalohihi_d,datalohihi_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datalolohi_d,datalolohi_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datalohilo_d,datalohilo_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datalololo_d,datalololo_h,szdata,cudaMemcpyHostToDevice);

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
         dbl8_padded_convjobs<<<jobnbr,BS>>>
            (datahihihi_d,datahilohi_d,datahihilo_d,datahilolo_d,
             datalohihi_d,datalolohi_d,datalohilo_d,datalololo_d,
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
         dbl8_update_addjobs<<<jobnbr,BS>>>
            (datahihihi_d,datahilohi_d,datahihilo_d,datahilolo_d,
             datalohihi_d,datalolohi_d,datalohilo_d,datalololo_d,
             in1ix_d,in2ix_d,outix_d,deg1);
         cudaEventRecord(stop);
         cudaEventSynchronize(stop);
         cudaEventElapsedTime(&milliseconds,start,stop);
         *addlapms += milliseconds;
      }
      free(in1ix_h); free(in2ix_h); free(outix_h);
   }
   gettimeofday(&endtime,0);
   cudaMemcpy(datahihihi_h,datahihihi_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datahilohi_h,datahilohi_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datahihilo_h,datahihilo_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datahilolo_h,datahilolo_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datalohihi_h,datalohihi_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datalolohi_h,datalolohi_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datalohilo_h,datalohilo_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datalololo_h,datalololo_d,szdata,cudaMemcpyDeviceToHost);
   *elapsedms = *cnvlapms + *addlapms;
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   *walltimesec = seconds + microseconds*1.0e-6;

   // convoluted_data2_to_output
   //    (data_h,output,dim,nbr,deg,nvr,idx,fstart,bstart,cstart,verbose);
   added_data8_to_output
      (datahihihi_h,datahilohi_h,datahihilo_h,datahilolo_h,
       datalohihi_h,datalolohi_h,datalohilo_h,datalololo_h,
       outputhihihi,outputhilohi,outputhihilo,outputhilolo,
       outputlohihi,outputlolohi,outputlohilo,outputlololo,
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
