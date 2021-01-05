// The file dbl_polynomials_kernels.cu defines the kernels with prototypes
// in dbl_polynomials_kernels.h.

#include <iostream>
#include "job_coordinates.h"
#include "dbl_convolutions_kernels.h"
#include "dbl_polynomials_kernels.h"

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

void convoluted_data_to_output
 ( double *data, double **output, int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart, bool verbose )
{
   const int deg1 = deg+1;
   int ix0,ix1,ix2;

   for(int i=0; i<=deg; i++) output[dim][i] = data[i];
   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++) output[i][j] = 0.0;

   for(int k=0; k<nbr; k++)
   {
      ix1 = fstart[k] + (nvr[k]-1)*deg1;
      
      if(verbose)
         cout << "monomial " << k << " update starts at " << ix1 << endl;

      for(int i=0; i<=deg; i++) output[dim][i] += data[ix1++];

      ix0 = idx[k][0];
      if(nvr[k] == 1)
      {
         ix1 = (1 + k)*deg1;
            
         for(int i=0; i<=deg; i++) output[ix0][i] += data[ix1++];
      }
      else
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

               if(verbose)
                  cout << "monomial " << k << " derivative " << ix0
                       << " update starts at " << ix1 << endl;

               for(int i=0; i<=deg; i++) output[ix0][i] += data[ix1++];
            }
         }
      }
   }
}

void added_data_to_output
 ( double *data, double **output, int dim, int nbr, int deg, int *nvr,
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

   for(int i=0; i<=deg; i++) output[dim][i] = data[ix++];

   int cnt = jobs.get_differential_count(0);
   if(cnt == 0) // it could be there is no first variable anywhere ...
   {
      for(int i=0; i<=deg; i++) output[0][i] = 0.0;
   }
   else
   {
      int ix0 = jobs.get_differential_index(0,cnt);
      int ix2 = nvr[ix0]-3;
      if(ix2 < 0) ix2 = 0; // on GPU, one backward item less

      ix = bstart[ix0] + ix2*deg1;
      
      if(verbose)
         cout << "Updating derivative 0 at " << ix << " in data." << endl;

      for(int i=0; i<=deg; i++) output[0][i] = data[ix++];

      for(int k=1; k<dim; k++) // updating all other derivatives
      {
         int cnt = jobs.get_differential_count(k);
         if(cnt == 0) // it could be there is no variable k anywhere ...
         {
            for(int i=0; i<=deg; i++) output[k][i] = 0.0;
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

               for(int i=0; i<=deg; i++) output[k][i] = data[ix++];
            }
            else if(idx[ix0][nvr[ix0]-1] == k) // k is last variable
            {
               int ix2 = nvr[ix0]-2;
   
               if(verbose)
                  cout << "Updating derivative " << k 
                       << " at " << ix << " in data." << endl;

               ix = fstart[ix0] + ix2*deg1;

               for(int i=0; i<=deg; i++) output[k][i] = data[ix++];
            }
            else // derivative is in some cross product
            {
               int ix2 = jobs.position(nvr[ix0],idx[ix0],k) - 1;
   
               if(verbose)
                  cout << "Updating derivative " << k 
                       << " at " << ix << " in data." << endl;

               ix = cstart[ix0] + ix2*deg1;

               for(int i=0; i<=deg; i++) output[k][i] = data[ix++];
            }
         }
      }
   }
}

void GPU_dbl_poly_evaldiff
 ( int BS, int dim, int nbr, int deg, int *nvr, int **idx,
   double *cst, double **cff, double **input, double **output,
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

   double *data_h = new double[totalcff];        // data on host
   int ix = 0;
   for(int i=0; i<deg1; i++) data_h[ix++] = cst[i];
   for(int i=0; i<nbr; i++)
      for(int j=0; j<deg1; j++) data_h[ix++] = cff[i][j];
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++) data_h[ix++] = input[i][j];

   double *data_d;                               // device data
   const size_t szdata = totalcff*sizeof(double);
   cudaMalloc((void**)&data_d,szdata);
   cudaMemcpy(data_d,data_h,szdata,cudaMemcpyHostToDevice);

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

         dbl_padded_convjobs<<<jobnbr,BS>>>
            (data_d,in1ix_d,in2ix_d,outix_d,deg1);
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

         dbl_update_addjobs<<<jobnbr,BS>>>
            (data_d,in1ix_d,in2ix_d,outix_d,deg1);
      }
      free(in1ix_h); free(in2ix_h); free(outix_h);
   }
   cudaMemcpy(data_h,data_d,szdata,cudaMemcpyDeviceToHost);

   // convoluted_data_to_output
   //    (data_h,output,dim,nbr,deg,nvr,idx,fstart,bstart,cstart,verbose);
   added_data_to_output
      (data_h,output,dim,nbr,deg,nvr,idx,fstart,bstart,cstart,addjobs,verbose);
}
