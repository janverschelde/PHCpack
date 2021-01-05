// The file dbl3_polynomials_kernels.cu defines the kernels with prototypes
// in dbl3_polynomials_kernels.h.

#include <iostream>
#include <iomanip>
#include "job_coordinates.h"
#include "triple_double_functions.h"
#ifdef gpufun
#include "double_double_gpufun.cu"
#include "triple_double_gpufun.cu"
#endif
#include "dbl3_convolutions_kernels.h"
#include "dbl3_polynomials_kernels.h"

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

void convoluted_data3_to_output
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

void added_data3_to_output
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
   if(cnt == 0) // it could be there is no first variable anywhere ...
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
      for(int k=1; k<dim; k++) // updating all other derivatives
      {
         int cnt = jobs.get_differential_count(k);
         if(cnt == 0) // it could be there is no variable k anywhere ...
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
}

void GPU_dbl3_poly_evaldiff
 ( int BS, int dim, int nbr, int deg, int *nvr, int **idx,
   double *csthi, double *cstmi, double *cstlo,
   double **cffhi, double **cffmi, double **cfflo,
   double **inputhi, double **inputmi,  double **inputlo,
   double **outputhi, double **outputmi, double **outputlo,
   ConvolutionJobs cnvjobs, AdditionJobs addjobs, double *elapsedms,
   bool verbose )
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

   cudaEvent_t start,stop;
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   *elapsedms = 0.0;
   float milliseconds;

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
         *elapsedms += milliseconds;
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
         *elapsedms += milliseconds;
      }
      free(in1ix_h); free(in2ix_h); free(outix_h);
   }
   cudaMemcpy(datahi_h,datahi_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datami_h,datami_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datalo_h,datalo_d,szdata,cudaMemcpyDeviceToHost);

   // convoluted_data2_to_output
   //    (data_h,output,dim,nbr,deg,nvr,idx,fstart,bstart,cstart,verbose);
   added_data3_to_output
      (datahi_h,datami_h,datalo_h,outputhi,outputmi,outputlo,
       dim,nbr,deg,nvr,idx,fstart,bstart,cstart,addjobs,verbose);
   if(verbose)
   {
      cout << "Time spent by all kernels in milliseconds : ";
      cout << fixed << setprecision(2) << *elapsedms << endl;
      cout << scientific << setprecision(16);
   }
}
