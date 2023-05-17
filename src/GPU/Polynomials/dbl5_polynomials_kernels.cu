// The file dbl5_polynomials_kernels.cu defines the kernels with prototypes
// in dbl5_polynomials_kernels.h.

#include <iostream>
#include <iomanip>
#ifdef winwalltime
#include "wingettimeofday.h"
#else
#include <sys/time.h>
#endif
#include "job_coordinates.h"
#include "penta_double_functions.h"
#ifdef gpufun
#include "double_double_gpufun.cu"
#include "penta_double_gpufun.cu"
#endif
#include "dbl5_polynomials_kernels.h"

// The constant pd_shmemsize is the bound on the shared memory size.

#define pd_shmemsize 192

using namespace std;

__global__ void dbl5_padded_convjobs
 ( double *datatb, double *dataix, double *datami, 
   double *datarg, double *datapk,
   int *in1idx, int *in2idx, int *outidx, int dim )
{
   const int bdx = blockIdx.x;           // index to the convolution job
   const int tdx = threadIdx.x;          // index to the output of the job
   const int idx1 = in1idx[bdx] + tdx;
   const int idx2 = in2idx[bdx] + tdx;
   const int idx3 = outidx[bdx] + tdx;

   __shared__ double xvtb[pd_shmemsize];
   __shared__ double xvix[pd_shmemsize];
   __shared__ double xvmi[pd_shmemsize];
   __shared__ double xvrg[pd_shmemsize];
   __shared__ double xvpk[pd_shmemsize];
   __shared__ double yvtb[2*pd_shmemsize];
   __shared__ double yvix[2*pd_shmemsize];
   __shared__ double yvmi[2*pd_shmemsize];
   __shared__ double yvrg[2*pd_shmemsize];
   __shared__ double yvpk[2*pd_shmemsize];
   __shared__ double zvtb[pd_shmemsize];
   __shared__ double zvix[pd_shmemsize];
   __shared__ double zvmi[pd_shmemsize];
   __shared__ double zvrg[pd_shmemsize];
   __shared__ double zvpk[pd_shmemsize];

   double prdtb,prdix,prdmi,prdrg,prdpk;
   int ydx = dim + tdx;

   xvtb[tdx] = datatb[idx1];  // loading first input
   xvix[tdx] = dataix[idx1]; 
   xvmi[tdx] = datami[idx1]; 
   xvrg[tdx] = datarg[idx1]; 
   xvpk[tdx] = datapk[idx1]; 
   yvtb[tdx] = 0.0;             // padded with zeros
   yvix[tdx] = 0.0;
   yvmi[tdx] = 0.0;
   yvrg[tdx] = 0.0;
   yvpk[tdx] = 0.0;
   yvtb[ydx] = datatb[idx2];  // loading second input
   yvix[ydx] = dataix[idx2];
   yvmi[ydx] = datami[idx2];
   yvrg[ydx] = datarg[idx2];
   yvpk[ydx] = datapk[idx2];

   __syncthreads();

   // zv[tdx] = xv[0]*yv[tdx];
   pdg_mul( xvtb[0],   xvix[0],   xvmi[0],   xvrg[0],   xvpk[0],
            yvtb[ydx], yvix[ydx], yvmi[ydx], yvrg[ydx], yvpk[ydx],
           &zvtb[tdx],&zvix[tdx],&zvmi[tdx],&zvrg[tdx],&zvpk[tdx]);
   __syncthreads();

   for(int i=1; i<dim; i++) // zv[tdx] = zv[tdx] + xv[i]*yv[dim+tdx-i];
   {
      ydx = dim + tdx - i;

      pdg_mul( xvtb[i],  xvix[i],  xvmi[i],  xvrg[i],  xvpk[i],
               yvtb[ydx],yvix[ydx],yvmi[ydx],yvrg[ydx],yvpk[ydx],
              &prdtb,  &prdix,   &prdmi,   &prdrg,   &prdpk);
      __syncthreads();

      pdg_inc(&zvtb[tdx],&zvix[tdx],&zvmi[tdx],&zvrg[tdx],&zvpk[tdx],
              prdtb,     prdix,     prdmi,     prdrg,     prdpk);
      __syncthreads();
   }
   __syncthreads();

   datatb[idx3] = zvtb[tdx]; // storing the output
   dataix[idx3] = zvix[tdx];
   datami[idx3] = zvmi[tdx];
   datarg[idx3] = zvrg[tdx];
   datapk[idx3] = zvpk[tdx];
}

__global__ void dbl5_update_addjobs
 ( double *datatb, double *dataix, double *datami,
   double *datarg, double *datapk,
   int *in1idx, int *in2idx, int *outidx, int dim )
{
   const int bdx = blockIdx.x;           // index to the convolution job
   const int tdx = threadIdx.x;          // index to the output of the job
   const int idx1 = in1idx[bdx] + tdx;
   const int idx2 = in2idx[bdx] + tdx;
   const int idx3 = outidx[bdx] + tdx;

   __shared__ double xvtb[pd_shmemsize];
   __shared__ double xvix[pd_shmemsize];
   __shared__ double xvmi[pd_shmemsize];
   __shared__ double xvrg[pd_shmemsize];
   __shared__ double xvpk[pd_shmemsize];
   __shared__ double yvtb[pd_shmemsize];
   __shared__ double yvix[pd_shmemsize];
   __shared__ double yvmi[pd_shmemsize];
   __shared__ double yvrg[pd_shmemsize];
   __shared__ double yvpk[pd_shmemsize];
   __shared__ double zvtb[pd_shmemsize];
   __shared__ double zvix[pd_shmemsize];
   __shared__ double zvmi[pd_shmemsize];
   __shared__ double zvrg[pd_shmemsize];
   __shared__ double zvpk[pd_shmemsize];

   xvtb[tdx] = datatb[idx1];  // loading first input
   xvix[tdx] = dataix[idx1];
   xvmi[tdx] = datami[idx1];
   xvrg[tdx] = datarg[idx1];
   xvpk[tdx] = datapk[idx1];
   yvtb[tdx] = datatb[idx2];  // loading second input
   yvix[tdx] = dataix[idx2];
   yvmi[tdx] = datami[idx2];
   yvrg[tdx] = datarg[idx2];
   yvpk[tdx] = datapk[idx2];

   // zv[tdx] = xv[tdx] + yv[tdx];

   __syncthreads();

   pdg_add( xvtb[tdx], xvix[tdx], xvmi[tdx], xvrg[tdx], xvpk[tdx],
            yvtb[tdx], yvix[tdx], yvmi[tdx], yvrg[tdx], yvpk[tdx],
           &zvtb[tdx],&zvix[tdx],&zvmi[tdx],&zvrg[tdx],&zvpk[tdx]);

   __syncthreads();

   datatb[idx3] = zvtb[tdx]; // storing the output
   dataix[idx3] = zvix[tdx];
   datami[idx3] = zvmi[tdx];
   datarg[idx3] = zvrg[tdx];
   datapk[idx3] = zvpk[tdx];
}

void convoluted_data5_to_output
 ( double *datatb, double *dataix, double *datami,
   double *datarg, double *datapk,
   double **outputtb, double **outputix, double **outputmi,
   double **outputrg, double **outputpk,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart, bool verbose )
{
   const int deg1 = deg+1;
   int ix0,ix1,ix2;

   for(int i=0; i<=deg; i++) // output[dim][i] = data[i];
   {
      outputtb[dim][i] = datatb[i];
      outputix[dim][i] = dataix[i];
      outputmi[dim][i] = datami[i];
      outputrg[dim][i] = datarg[i];
      outputpk[dim][i] = datapk[i];
   }
   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++) // output[i][j] = 0.0;
      {
         outputtb[i][j] = 0.0;
         outputix[i][j] = 0.0;
         outputmi[i][j] = 0.0;
         outputrg[i][j] = 0.0;
         outputpk[i][j] = 0.0;
      }

   for(int k=0; k<nbr; k++)
   {
      ix1 = fstart[k] + (nvr[k]-1)*deg1;
      
      if(verbose)
         cout << "monomial " << k << " update starts at " << ix1 << endl;

      for(int i=0; i<=deg; i++) // output[dim][i] += data[ix1++];
         pdf_inc(&outputtb[dim][i],&outputix[dim][i],&outputmi[dim][i],
                 &outputrg[dim][i],&outputpk[dim][i],
                    datatb[ix1],      dataix[ix1],      datami[ix1],
                    datarg[ix1],      datapk[ix1++]);
     
      ix0 = idx[k][0];
      if(nvr[k] == 1)
      {
         ix1 = (1 + k)*deg1;
            
         for(int i=0; i<=deg; i++) // output[ix0][i] += data[ix1++];
            pdf_inc(&outputtb[ix0][i],&outputix[ix0][i],&outputmi[ix0][i],
                    &outputrg[ix0][i],&outputpk[ix0][i],
                       datatb[ix1],      dataix[ix1],      datami[ix1],
                       datarg[ix1],      datapk[ix1++]);
      }
      else
      {                               // update first and last derivative
         ix2 = nvr[k]-3;
         if(ix2 < 0) ix2 = 0;
         ix1 = bstart[k] + ix2*deg1;

         for(int i=0; i<=deg; i++) // output[ix0][i] += data[ix1++];
            pdf_inc(&outputtb[ix0][i],&outputix[ix0][i],&outputmi[ix0][i],
                    &outputrg[ix0][i],&outputpk[ix0][i],
                       datatb[ix1],      dataix[ix1],      datami[ix1],
                       datarg[ix1],      datapk[ix1++]);

         ix2 = nvr[k]-2;
         ix1 = fstart[k] + ix2*deg1;
         ix0 = idx[k][ix2+1];

         for(int i=0; i<=deg; i++) // output[ix0][i] += data[ix1++];
            pdf_inc(&outputtb[ix0][i],&outputix[ix0][i],&outputmi[ix0][i],
                    &outputrg[ix0][i],&outputpk[ix0][i],
                       datatb[ix1],      dataix[ix1],      datami[ix1],
                       datarg[ix1],      datapk[ix1++]);
 
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
                  pdf_inc(&outputtb[ix0][i],&outputix[ix0][i],
                          &outputmi[ix0][i],
                          &outputrg[ix0][i],&outputpk[ix0][i],
                             datatb[ix1],      dataix[ix1],
                             datami[ix1],
                             datarg[ix1],      datapk[ix1++]);
            }
         }
      }
   }
}

void added_data5_to_output
 ( double *datatb, double *dataix, double *datami,
   double *datarg, double *datapk,
   double **outputtb, double **outputix, double **outputmi,
   double **outputrg, double **outputpk,
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
      outputtb[dim][i] = datatb[ix];
      outputix[dim][i] = dataix[ix];
      outputmi[dim][i] = datami[ix];
      outputrg[dim][i] = datarg[ix];
      outputpk[dim][i] = datapk[ix++];
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
            outputtb[0][i] = 0.0;
            outputix[0][i] = 0.0;
            outputmi[0][i] = 0.0;
            outputrg[0][i] = 0.0;
            outputpk[0][i] = 0.0;
         }
      }
      else
      {
         int cffidx = (1 + difidx)*deg1;

         if(verbose)
            cout << "updating derivative with coefficient ..." << endl;

         for(int i=0; i<=deg; i++)
         {
            outputtb[0][i] = datatb[cffidx];
            outputix[0][i] = dataix[cffidx];
            outputmi[0][i] = datami[cffidx];
            outputrg[0][i] = datarg[cffidx];
            outputpk[0][i] = datapk[cffidx++];
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
         outputtb[0][i] = datatb[ix];
         outputix[0][i] = dataix[ix];
         outputmi[0][i] = datami[ix];
         outputrg[0][i] = datarg[ix];
         outputpk[0][i] = datapk[ix++];
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
               outputtb[k][i] = 0.0;
               outputix[k][i] = 0.0;
               outputmi[k][i] = 0.0;
               outputrg[k][i] = 0.0;
               outputpk[k][i] = 0.0;
            }
         }
         else
         {
            int cffidx = (1 + difidx)*deg1;

            if(verbose)
               cout << "updating derivative with coefficient ..." << endl;

            for(int i=0; i<=deg; i++)
            {
               outputtb[k][i] = datatb[cffidx];
               outputix[k][i] = dataix[cffidx];
               outputmi[k][i] = datami[cffidx];
               outputrg[k][i] = datarg[cffidx];
               outputpk[k][i] = datapk[cffidx++];
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
               outputtb[k][i] = datatb[ix];
               outputix[k][i] = dataix[ix];
               outputmi[k][i] = datami[ix];
               outputrg[k][i] = datarg[ix];
               outputpk[k][i] = datapk[ix++];
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
               outputtb[k][i] = datatb[ix];
               outputix[k][i] = dataix[ix];
               outputmi[k][i] = datami[ix];
               outputrg[k][i] = datarg[ix];
               outputpk[k][i] = datapk[ix++];
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
               outputtb[k][i] = datatb[ix];
               outputix[k][i] = dataix[ix];
               outputmi[k][i] = datami[ix];
               outputrg[k][i] = datarg[ix];
               outputpk[k][i] = datapk[ix++];
            }
         }
      }
   }
}

void GPU_dbl5_poly_evaldiff
 ( int BS, int dim, int nbr, int deg, int *nvr, int **idx,
   double *csttb, double *cstix, double *cstmi,
   double *cstrg, double *cstpk,
   double **cfftb, double **cffix, double **cffmi,
   double **cffrg, double **cffpk,
   double **inputtb, double **inputix, double **inputmi,
   double **inputrg, double **inputpk,
   double **outputtb, double **outputix, double **outputmi,
   double **outputrg, double **outputpk,
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
   double *datatb_h = new double[totalcff];        // data on host
   double *dataix_h = new double[totalcff];
   double *datami_h = new double[totalcff];
   double *datarg_h = new double[totalcff];
   double *datapk_h = new double[totalcff];
   int ix = 0;
   for(int i=0; i<deg1; i++)
   {
      datatb_h[ix] = csttb[i];
      dataix_h[ix] = cstix[i];
      datami_h[ix] = cstmi[i];
      datarg_h[ix] = cstrg[i];
      datapk_h[ix++] = cstpk[i];
   }
   for(int i=0; i<nbr; i++)
      for(int j=0; j<deg1; j++)
      {
         datatb_h[ix] = cfftb[i][j];
         dataix_h[ix] = cffix[i][j];
         datami_h[ix] = cffmi[i][j];
         datarg_h[ix] = cffrg[i][j];
         datapk_h[ix++] = cffpk[i][j];
      }
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++)
      {
         datatb_h[ix] = inputtb[i][j];
         dataix_h[ix] = inputix[i][j];
         datami_h[ix] = inputmi[i][j];
         datarg_h[ix] = inputrg[i][j];
         datapk_h[ix++] = inputpk[i][j];
      }

   double *datatb_d;                               // device data
   double *dataix_d;
   double *datami_d;
   double *datarg_d;
   double *datapk_d;
   const size_t szdata = totalcff*sizeof(double);
   cudaMalloc((void**)&datatb_d,szdata);
   cudaMalloc((void**)&dataix_d,szdata);
   cudaMalloc((void**)&datami_d,szdata);
   cudaMalloc((void**)&datarg_d,szdata);
   cudaMalloc((void**)&datapk_d,szdata);
   cudaMemcpy(datatb_d,datatb_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(dataix_d,dataix_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datami_d,datami_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datarg_d,datarg_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datapk_d,datapk_h,szdata,cudaMemcpyHostToDevice);

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
         dbl5_padded_convjobs<<<jobnbr,BS>>>
            (datatb_d,dataix_d,datami_d,datarg_d,datapk_d,
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
         dbl5_update_addjobs<<<jobnbr,BS>>>
            (datatb_d,dataix_d,datami_d,datarg_d,datapk_d,
             in1ix_d,in2ix_d,outix_d,deg1);
         cudaEventRecord(stop);
         cudaEventSynchronize(stop);
         cudaEventElapsedTime(&milliseconds,start,stop);
         *addlapms += milliseconds;
      }
      free(in1ix_h); free(in2ix_h); free(outix_h);
   }
   gettimeofday(&endtime,0);
   cudaMemcpy(datatb_h,datatb_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(dataix_h,dataix_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datami_h,datami_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datarg_h,datarg_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datapk_h,datapk_d,szdata,cudaMemcpyDeviceToHost);
   *elapsedms = *cnvlapms + *addlapms;
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   *walltimesec = seconds + microseconds*1.0e-6;

   // convoluted_data2_to_output
   //    (data_h,output,dim,nbr,deg,nvr,idx,fstart,bstart,cstart,verbose);
   added_data5_to_output
      (datatb_h,dataix_h,datami_h,datarg_h,datapk_h,
       outputtb,outputix,outputmi,outputrg,outputpk,
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
