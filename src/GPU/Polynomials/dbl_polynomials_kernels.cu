// The file dbl_polynomials_kernels.cu defines the kernels with prototypes
// in dbl_polynomials_kernels.h.

#include <iostream>
#include "dbl_convolutions_kernels.h"
#include "dbl_polynomials_kernels.h"

using namespace std;

int coefficient_count ( int dim, int nbr, int deg, int *nvr )
{
   int count = nbr + dim;

   for(int i=0; i<nbr; i++)
   {
      count = count + nvr[i];
      if(nvr[i] == 2)
         count = count + 1;
      else if(nvr[i] > 2)
         count = count + 2*(nvr[i]-2);
   }
   count = count*(deg+1);

   return count;
}

void coefficient_indices
 ( int dim, int nbr, int deg, int *nvr,
   int *fsums, int *bsums, int *csums,
   int *fstart, int *bstart, int *cstart )
{
   fsums[0] = nvr[0]; bsums[0] = 0; csums[0] = 0;

   if(nvr[0] == 2)
   {
      bsums[0] = 1;
   }
   else if(nvr[0] > 2)
   {
      bsums[0] = nvr[0] - 2;
      csums[0] = nvr[0] - 2;
   }
   for(int i=1; i<nbr; i++)
   {
      fsums[i] = fsums[i-1] + nvr[i];
      if(nvr[i] < 2)
      {
         bsums[i] = bsums[i-1];
         csums[i] = csums[i-1];
      }
      else if(nvr[i] == 2)
      {
         bsums[i] = bsums[i-1] + 1;
         csums[i] = csums[i-1];
      }
      else // nvr[i] > 2
      {
         bsums[i] = bsums[i-1] + nvr[i] - 2;
         csums[i] = csums[i-1] + nvr[i] - 2;
      }
   }
   fstart[0] = (nbr+dim)*(deg+1);
   for(int i=1; i<nbr; i++)
      fstart[i] = fstart[0] + fsums[i-1]*(deg+1);

   bstart[0] = fstart[0] + fsums[nbr-1]*(deg+1);
   for(int i=1; i<nbr; i++)
      bstart[i] = bstart[0] + bsums[i-1]*(deg+1);

   cstart[0] = bstart[0] + bsums[nbr-1]*(deg+1);
   for(int i=1; i<nbr; i++)
      cstart[i] = cstart[0] + csums[i-1]*(deg+1);
}

void job_indices
 ( ConvolutionJob job, int *inp1ix, int *inp2ix, int *outidx,
   int dim, int nbr, int deg,
   int *fstart, int *bstart, int *cstart, bool verbose )
{
   const int monidx = job.get_monomial_index();
   const int jobinp1tp = job.get_first_type();
   const int jobinp1ix = job.get_first_input();
   const int jobinp2tp = job.get_second_type();
   const int jobinp2ix = job.get_second_input();
   const int joboutptp = job.get_output_type();
   const int joboutidx = job.get_output_index();
   const int deg1 = deg+1;

   if(verbose)
   {
      cout << "  in1 type : " << jobinp1tp << ", idx : " << jobinp1ix;
      cout << "  in2 type : " << jobinp2tp << ", idx : " << jobinp2ix;
      cout << "  out type : " << joboutptp << ", idx : " << joboutidx;
      cout << endl;
   }
   if(jobinp1tp < 0)           // first input is coefficient
      *inp1ix = monidx*deg1;
   else if(jobinp1tp == 0)     // first input is input series
      *inp1ix = (nbr + jobinp1ix)*deg1;
   else if(jobinp1tp == 1)     // first input is forward product
      *inp1ix = fstart[monidx] + jobinp1ix*deg1;
   else if(jobinp1tp == 2)     // first input is backward product
      *inp1ix = bstart[monidx] + jobinp1ix*deg1;
   else if(jobinp1tp == 3)     // first input is cross product
      *inp1ix = cstart[monidx] + jobinp1ix*deg1;

   if(jobinp2tp < 0)           // second input is coefficient
      *inp2ix = monidx*deg1;
   else if(jobinp2tp == 0)     // second input is input series
      *inp2ix = (nbr + jobinp2ix)*deg1;
   else if(jobinp2tp == 1)     // second input is forward product
      *inp2ix = fstart[monidx] + jobinp2ix*deg1;
   else if(jobinp2tp == 2)     // second input is backward product
      *inp2ix = bstart[monidx] + jobinp2ix*deg1;
   else if(jobinp2tp == 3)     // second input is cross product
      *outidx = cstart[monidx] + jobinp2ix*deg1;

   if(joboutptp == 1)          // output is forward product
      *outidx = fstart[monidx] + joboutidx*deg1;
   else if(joboutptp == 2)    // output is backward product
      *outidx = bstart[monidx] + joboutidx*deg1;
   else if(joboutptp == 3)    // output is cross product
      *outidx = cstart[monidx] + joboutidx*deg1;

   if(verbose)
   {
      cout << "-> inp1ix : " << *inp1ix
           << ", inp2ix : " << *inp2ix
           << ", outidx : " << *outidx << endl;
   }
}

void jobs_coordinates
 ( ConvolutionJobs jobs, int layer,
   int *inp1ix, int *inp2ix, int *outidx,
   int dim, int nbr, int deg,
   int *fstart, int *bstart, int *cstart, bool verbose )
{ 
   for(int i=0; i<jobs.get_layer_count(layer); i++)
      job_indices(jobs.get_job(layer,i),&inp1ix[i],&inp2ix[i],&outidx[i],
                  dim,nbr,deg,fstart,bstart,cstart,verbose);
}

__global__ void dbl_padded_convjobs
 ( double *data, int *in1idx, int *in2idx, int *outidx, int dim )
{
   const int bdx = blockIdx.x;           // index to the convolution job
   const int tdx = threadIdx.x;          // index to the output of the job
   const int idx1 = in1idx[bdx] + tdx;
   const int idx2 = in2idx[bdx] + tdx;
   const int idx3 = outidx[bdx] + tdx;

   __shared__ double xv[d_shmemsize];
   __shared__ double yv[d_shmemsize];
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

void data_to_output
 ( double *data, double *cst, double **output,
   int dim, int nbr, int deg, int *nvr,
   int *fsums, int *bsums, int *csums,
   int *fstart, int *bstart, int *cstart, bool verbose )
{
   const int deg1 = deg+1;

   for(int i=0; i<=deg; i++) output[dim][i] = cst[i];
   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++) output[i][j] = 0.0;

   for(int k=0; k<nbr; k++)
   {
      int idx = fstart[k] + (nvr[k]-1)*deg1;
      
      if(verbose)
         cout << "monomial " << k << " update starts at " << idx << endl;

      for(int i=0; i<=deg; i++) output[dim][i] += data[idx++];
   }
}

void GPU_dbl_poly_evaldiff
 ( int BS, int dim, int nbr, int deg, int *nvr, int **idx,
   double *cst, double **cff, double **input, double **output,
   ConvolutionJobs jobs, bool verbose )
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
   for(int i=0; i<nbr; i++)
      for(int j=0; j<deg1; j++) data_h[ix++] = cff[i][j];
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++) data_h[ix++] = input[i][j];

   double *data_d;                               // device data
   const size_t szdata = totalcff*sizeof(double);
   cudaMalloc((void**)&data_d,szdata);
   cudaMemcpy(data_d,data_h,szdata,cudaMemcpyHostToDevice);

   for(int k=0; k<jobs.get_depth(); k++)
   {
      const int jobnbr = jobs.get_layer_count(k);
      int *in1ix_h = new int[jobnbr];
      int *in2ix_h = new int[jobnbr];
      int *outix_h = new int[jobnbr];

      if(verbose) cout << "preparing jobs at layer " << k << " ..." << endl;

      jobs_coordinates(jobs,k,in1ix_h,in2ix_h,outix_h,dim,nbr,deg,
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
   cudaMemcpy(data_h,data_d,szdata,cudaMemcpyDeviceToHost);

   data_to_output(data_h,cst,output,dim,nbr,deg,nvr,
                  fsums,bsums,csums,fstart,bstart,cstart,verbose);
}
