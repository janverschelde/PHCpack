// The file dbl2_polynomials_kernels.cu defines the kernels with prototypes
// in dbl2_polynomials_kernels.h.

#include <iostream>
#include "job_coordinates.h"
#include "double_double_functions.h"
#ifdef gpufun
#include "double_double_gpufun.cu"
#endif
#include "dbl2_convolutions_kernels.h"
#include "dbl2_polynomials_kernels.h"

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

void convoluted_data2_to_output
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

void added_data2_to_output
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

void GPU_dbl2_poly_evaldiff
 ( int BS, int dim, int nbr, int deg, int *nvr, int **idx,
   double *csthi, double *cstlo, double **cffhi, double **cfflo,
   double **inputhi, double **inputlo,
   double **outputhi, double **outputlo,
   ConvolutionJobs cnvjobs, AdditionJobs addjobs, bool verbose )
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

         dbl2_padded_convjobs<<<jobnbr,BS>>>
            (datahi_d,datalo_d,in1ix_d,in2ix_d,outix_d,deg1);
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

         dbl2_update_addjobs<<<jobnbr,BS>>>
            (datahi_d,datalo_d,in1ix_d,in2ix_d,outix_d,deg1);
      }
      free(in1ix_h); free(in2ix_h); free(outix_h);
   }
   cudaMemcpy(datahi_h,datahi_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datalo_h,datalo_d,szdata,cudaMemcpyDeviceToHost);

   // convoluted_data2_to_output
   //    (data_h,output,dim,nbr,deg,nvr,idx,fstart,bstart,cstart,verbose);
   added_data2_to_output
      (datahi_h,datalo_h,outputhi,outputlo,
       dim,nbr,deg,nvr,idx,fstart,bstart,cstart,addjobs,verbose);
}
