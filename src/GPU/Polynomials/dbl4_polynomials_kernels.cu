// The file dbl4_polynomials_kernels.cu defines the kernels with prototypes
// in dbl4_polynomials_kernels.h.

#include <iostream>
#include "job_coordinates.h"
#include "quad_double_functions.h"
#ifdef gpufun
#include "double_double_gpufun.cu"
#include "quad_double_gpufun.cu"
#endif
#include "dbl4_convolutions_kernels.h"
#include "dbl4_polynomials_kernels.h"

using namespace std;

__global__ void dbl4_padded_convjobs
 ( double *datahihi, double *datalohi, double *datahilo, double *datalolo,
   int *in1idx, int *in2idx, int *outidx, int dim )
{
   const int bdx = blockIdx.x;           // index to the convolution job
   const int tdx = threadIdx.x;          // index to the output of the job
   const int idx1 = in1idx[bdx] + tdx;
   const int idx2 = in2idx[bdx] + tdx;
   const int idx3 = outidx[bdx] + tdx;

   __shared__ double xvhihi[qd_shmemsize];
   __shared__ double xvlohi[qd_shmemsize];
   __shared__ double xvhilo[qd_shmemsize];
   __shared__ double xvlolo[qd_shmemsize];
   __shared__ double yvhihi[2*qd_shmemsize];
   __shared__ double yvlohi[2*qd_shmemsize];
   __shared__ double yvhilo[2*qd_shmemsize];
   __shared__ double yvlolo[2*qd_shmemsize];
   __shared__ double zvhihi[qd_shmemsize];
   __shared__ double zvlohi[qd_shmemsize];
   __shared__ double zvhilo[qd_shmemsize];
   __shared__ double zvlolo[qd_shmemsize];

   double prdhihi,prdlohi,prdhilo,prdlolo;
   int ydx = dim + tdx;

   xvhihi[tdx] = datahihi[idx1];  // loading first input
   xvlohi[tdx] = datalohi[idx1]; 
   xvhilo[tdx] = datahilo[idx1]; 
   xvlolo[tdx] = datalolo[idx1]; 
   yvhihi[tdx] = 0.0;             // padded with zeros
   yvlohi[tdx] = 0.0;
   yvhilo[tdx] = 0.0;
   yvlolo[tdx] = 0.0;
   yvhihi[ydx] = datahihi[idx2];  // loading second input
   yvlohi[ydx] = datalohi[idx2];
   yvhilo[ydx] = datahilo[idx2];
   yvlolo[ydx] = datalolo[idx2];

   __syncthreads();

   // zv[tdx] = xv[0]*yv[tdx];
   qdg_mul( xvhihi[0],   xvlohi[0],   xvhilo[0],   xvlolo[0],
            yvhihi[ydx], yvlohi[ydx], yvhilo[ydx], yvlolo[ydx],
           &zvhihi[tdx],&zvlohi[tdx],&zvhilo[tdx],&zvlolo[tdx]);
   __syncthreads();

   for(int i=1; i<dim; i++) // zv[tdx] = zv[tdx] + xv[i]*yv[dim+tdx-i];
   {
      ydx = dim + tdx - i;

      qdg_mul( xvhihi[i],  xvlohi[i],  xvhilo[i],  xvlolo[i],
               yvhihi[ydx],yvlohi[ydx],yvhilo[ydx],yvlolo[ydx],
              &prdhihi,  &prdlohi,   &prdhilo,   &prdlolo);
      __syncthreads();

      qdg_inc(&zvhihi[tdx],&zvlohi[tdx],&zvhilo[tdx],&zvlolo[tdx],
              prdhihi,     prdlohi,     prdhilo,     prdlolo);
      __syncthreads();
   }
   __syncthreads();

   datahihi[idx3] = zvhihi[tdx]; // storing the output
   datalohi[idx3] = zvlohi[tdx];
   datahilo[idx3] = zvhilo[tdx];
   datalolo[idx3] = zvlolo[tdx];
}

__global__ void dbl4_update_addjobs
 ( double *datahihi, double *datalohi, double *datahilo, double *datalolo,
   int *in1idx, int *in2idx, int *outidx, int dim )
{
   const int bdx = blockIdx.x;           // index to the convolution job
   const int tdx = threadIdx.x;          // index to the output of the job
   const int idx1 = in1idx[bdx] + tdx;
   const int idx2 = in2idx[bdx] + tdx;
   const int idx3 = outidx[bdx] + tdx;

   __shared__ double xvhihi[qd_shmemsize];
   __shared__ double xvlohi[qd_shmemsize];
   __shared__ double xvhilo[qd_shmemsize];
   __shared__ double xvlolo[qd_shmemsize];
   __shared__ double yvhihi[qd_shmemsize];
   __shared__ double yvlohi[qd_shmemsize];
   __shared__ double yvhilo[qd_shmemsize];
   __shared__ double yvlolo[qd_shmemsize];
   __shared__ double zvhihi[qd_shmemsize];
   __shared__ double zvlohi[qd_shmemsize];
   __shared__ double zvhilo[qd_shmemsize];
   __shared__ double zvlolo[qd_shmemsize];

   xvhihi[tdx] = datahihi[idx1];  // loading first input
   xvlohi[tdx] = datalohi[idx1];
   xvhilo[tdx] = datahilo[idx1];
   xvlolo[tdx] = datalolo[idx1];
   yvhihi[tdx] = datahihi[idx2];  // loading second input
   yvlohi[tdx] = datalohi[idx2];
   yvhilo[tdx] = datahilo[idx2];
   yvlolo[tdx] = datalolo[idx2];

   // zv[tdx] = xv[tdx] + yv[tdx];

   __syncthreads();

   qdg_add( xvhihi[tdx], xvlohi[tdx], xvhilo[tdx], xvlolo[tdx],
            yvhihi[tdx], yvlohi[tdx], yvhilo[tdx], yvlolo[tdx],
           &zvhihi[tdx],&zvlohi[tdx],&zvhilo[tdx],&zvlolo[tdx]);

   __syncthreads();

   datahihi[idx3] = zvhihi[tdx]; // storing the output
   datalohi[idx3] = zvlohi[tdx];
   datahilo[idx3] = zvhilo[tdx];
   datalolo[idx3] = zvlolo[tdx];
}

void convoluted_data4_to_output
 ( double *datahihi, double *datalohi, double *datahilo, double *datalolo,
   double **outputhihi, double **outputlohi,
   double **outputhilo, double **outputlolo,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart, bool verbose )
{
   const int deg1 = deg+1;
   int ix0,ix1,ix2;

   for(int i=0; i<=deg; i++) // output[dim][i] = data[i];
   {
      outputhihi[dim][i] = datahihi[i];
      outputlohi[dim][i] = datalohi[i];
      outputhilo[dim][i] = datahilo[i];
      outputlolo[dim][i] = datalolo[i];
   }
   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++) // output[i][j] = 0.0;
      {
         outputhihi[i][j] = 0.0;
         outputlohi[i][j] = 0.0;
         outputhilo[i][j] = 0.0;
         outputlolo[i][j] = 0.0;
      }

   for(int k=0; k<nbr; k++)
   {
      ix1 = fstart[k] + (nvr[k]-1)*deg1;
      
      if(verbose)
         cout << "monomial " << k << " update starts at " << ix1 << endl;

      for(int i=0; i<=deg; i++) // output[dim][i] += data[ix1++];
         qdf_inc(&outputhihi[dim][i],&outputlohi[dim][i],
                 &outputhilo[dim][i],&outputlolo[dim][i],
                    datahihi[ix1],      datalohi[ix1],
                    datahilo[ix1],      datalolo[ix1++]);
     
      ix0 = idx[k][0];
      if(nvr[k] == 1)
      {
         ix1 = (1 + k)*deg1;
            
         for(int i=0; i<=deg; i++) // output[ix0][i] += data[ix1++];
            qdf_inc(&outputhihi[ix0][i],&outputlohi[ix0][i],
                    &outputhilo[ix0][i],&outputlolo[ix0][i],
                       datahihi[ix1],      datalohi[ix1],
                       datahilo[ix1],      datalolo[ix1++]);
      }
      else
      {                               // update first and last derivative
         ix2 = nvr[k]-3;
         if(ix2 < 0) ix2 = 0;
         ix1 = bstart[k] + ix2*deg1;

         for(int i=0; i<=deg; i++) // output[ix0][i] += data[ix1++];
            qdf_inc(&outputhihi[ix0][i],&outputlohi[ix0][i],
                    &outputhilo[ix0][i],&outputlolo[ix0][i],
                       datahihi[ix1],      datalohi[ix1],
                       datahilo[ix1],      datalolo[ix1++]);

         ix2 = nvr[k]-2;
         ix1 = fstart[k] + ix2*deg1;
         ix0 = idx[k][ix2+1];

         for(int i=0; i<=deg; i++) // output[ix0][i] += data[ix1++];
            qdf_inc(&outputhihi[ix0][i],&outputlohi[ix0][i],
                    &outputhilo[ix0][i],&outputlolo[ix0][i],
                       datahihi[ix1],      datalohi[ix1],
                       datahilo[ix1],      datalolo[ix1++]);
 
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
                  qdf_inc(&outputhihi[ix0][i],&outputlohi[ix0][i],
                          &outputhilo[ix0][i],&outputlolo[ix0][i],
                             datahihi[ix1],      datalohi[ix1],
                             datahilo[ix1],      datalolo[ix1++]);
            }
         }
      }
   }
}

void added_data4_to_output
 ( double *datahihi, double *datalohi, double *datahilo, double *datalolo,
   double **outputhihi, double **outputlohi,
   double **outputhilo, double **outputlolo,
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
      outputhihi[dim][i] = datahihi[ix];
      outputlohi[dim][i] = datalohi[ix];
      outputhilo[dim][i] = datahilo[ix];
      outputlolo[dim][i] = datalolo[ix++];
   }
   int cnt = jobs.get_differential_count(0);
   if(cnt == 0) // it could be there is no first variable anywhere ...
   {
      for(int i=0; i<=deg; i++) // output[0][i] = 0.0;
      {
         outputhihi[0][i] = 0.0;
         outputlohi[0][i] = 0.0;
         outputhilo[0][i] = 0.0;
         outputlolo[0][i] = 0.0;
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
         outputhihi[0][i] = datahihi[ix];
         outputlohi[0][i] = datalohi[ix];
         outputhilo[0][i] = datahilo[ix];
         outputlolo[0][i] = datalolo[ix++];
      }
      for(int k=1; k<dim; k++) // updating all other derivatives
      {
         int cnt = jobs.get_differential_count(k);
         if(cnt == 0) // it could be there is no variable k anywhere ...
         {
            for(int i=0; i<=deg; i++) // output[k][i] = 0.0;
            {
               outputhihi[k][i] = 0.0;
               outputlohi[k][i] = 0.0;
               outputhilo[k][i] = 0.0;
               outputlolo[k][i] = 0.0;
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
                  outputhihi[k][i] = datahihi[ix];
                  outputlohi[k][i] = datalohi[ix];
                  outputhilo[k][i] = datahilo[ix];
                  outputlolo[k][i] = datalolo[ix++];
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
                  outputhihi[k][i] = datahihi[ix];
                  outputlohi[k][i] = datalohi[ix];
                  outputhilo[k][i] = datahilo[ix];
                  outputlolo[k][i] = datalolo[ix++];
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
                  outputhihi[k][i] = datahihi[ix];
                  outputlohi[k][i] = datalohi[ix];
                  outputhilo[k][i] = datahilo[ix];
                  outputlolo[k][i] = datalolo[ix++];
               }
            }
         }
      }
   }
}

void GPU_dbl4_poly_evaldiff
 ( int BS, int dim, int nbr, int deg, int *nvr, int **idx,
   double *csthihi, double *cstlohi, double *csthilo, double *cstlolo,
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo,
   double **inputhihi, double **inputlohi,
   double **inputhilo, double **inputlolo,
   double **outputhihi, double **outputlohi,
   double **outputhilo, double **outputlolo,
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

   double *datahihi_h = new double[totalcff];        // data on host
   double *datalohi_h = new double[totalcff];
   double *datahilo_h = new double[totalcff];
   double *datalolo_h = new double[totalcff];
   int ix = 0;
   for(int i=0; i<deg1; i++)
   {
      datahihi_h[ix] = csthihi[i];
      datalohi_h[ix] = cstlohi[i];
      datahilo_h[ix] = csthilo[i];
      datalolo_h[ix++] = cstlolo[i];
   }
   for(int i=0; i<nbr; i++)
      for(int j=0; j<deg1; j++)
      {
         datahihi_h[ix] = cffhihi[i][j];
         datalohi_h[ix] = cfflohi[i][j];
         datahilo_h[ix] = cffhilo[i][j];
         datalolo_h[ix++] = cfflolo[i][j];
      }
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++)
      {
         datahihi_h[ix] = inputhihi[i][j];
         datalohi_h[ix] = inputlohi[i][j];
         datahilo_h[ix] = inputhilo[i][j];
         datalolo_h[ix++] = inputlolo[i][j];
      }

   double *datahihi_d;                               // device data
   double *datalohi_d;
   double *datahilo_d;
   double *datalolo_d;
   const size_t szdata = totalcff*sizeof(double);
   cudaMalloc((void**)&datahihi_d,szdata);
   cudaMalloc((void**)&datalohi_d,szdata);
   cudaMalloc((void**)&datahilo_d,szdata);
   cudaMalloc((void**)&datalolo_d,szdata);
   cudaMemcpy(datahihi_d,datahihi_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datalohi_d,datalohi_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datahilo_d,datahilo_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datalolo_d,datalolo_h,szdata,cudaMemcpyHostToDevice);

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

         dbl4_padded_convjobs<<<jobnbr,BS>>>
            (datahihi_d,datalohi_d,datahilo_d,datalolo_d,
             in1ix_d,in2ix_d,outix_d,deg1);
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

         dbl4_update_addjobs<<<jobnbr,BS>>>
            (datahihi_d,datalohi_d,datahilo_d,datalolo_d,
             in1ix_d,in2ix_d,outix_d,deg1);
      }
      free(in1ix_h); free(in2ix_h); free(outix_h);
   }
   cudaMemcpy(datahihi_h,datahihi_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datalohi_h,datalohi_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datahilo_h,datahilo_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datalolo_h,datalolo_d,szdata,cudaMemcpyDeviceToHost);

   // convoluted_data2_to_output
   //    (data_h,output,dim,nbr,deg,nvr,idx,fstart,bstart,cstart,verbose);
   added_data4_to_output
      (datahihi_h,datalohi_h,datahilo_h,datalolo_h,
       outputhihi,outputlohi,outputhilo,outputlolo,
       dim,nbr,deg,nvr,idx,fstart,bstart,cstart,addjobs,verbose);
}
