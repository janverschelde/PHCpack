// The file dbl_polynomials_kernels.cu defines the kernels with prototypes
// in dbl_polynomials_kernels.h.

#include <iostream>
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
   int dim, int nbr, int deg, int *fsums, int *bsums, int *csums,
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
   else if(jobinp1tp == 1) 
   {
      if(jobinp1ix == 0)       // first input is forward product
         *inp1ix = fstart[monidx];
      else
         *inp1ix = fstart[monidx] + fsums[jobinp1ix-1]*deg1;
   }
   else if(jobinp1tp == 2)     // first input is backward product
   {
      if(jobinp1ix == 0)
         *inp1ix = bstart[monidx];
      else
         *inp1ix = bstart[monidx] + bsums[jobinp1ix-1]*deg1;
   }
   else if(jobinp1tp == 3)     // first input is cross product
   {
      if(jobinp1ix == 0)
         *inp1ix = cstart[monidx];
      else
         *inp1ix = cstart[monidx] + csums[jobinp1ix-1]*deg1;
   }
   if(jobinp2tp < 0)           // second input is coefficient
      *inp2ix = monidx*deg1;
   else if(jobinp2tp == 0)     // second input is input series
      *inp2ix = (nbr + jobinp2ix)*deg1;
   else if(jobinp2tp == 1)     // second input is forward product
   {
      if(jobinp2ix == 0)
         *inp2ix = fstart[monidx];
      else
         *inp2ix = fstart[monidx] + fsums[jobinp2ix-1]*deg1;
   }
   else if(jobinp2tp == 2)     // second input is backward product
   {
      if(jobinp2ix == 0)
         *inp2ix = bstart[monidx];
      else
         *inp2ix = bstart[monidx] + bsums[jobinp2ix-1]*deg1;
   }
   else if(jobinp2tp == 3)     // second input is cross product
   {
      if(jobinp2ix == 0)
         *outidx = cstart[monidx];
      else
         *outidx = cstart[monidx] + csums[jobinp2ix-1]*deg1;
   }
   if(joboutptp == 1)          // output is forward product
   {
      if(joboutidx == 0)
         *outidx = fstart[monidx];
      else
         *outidx = fstart[monidx] + fsums[joboutidx-1]*deg1;
   }
   else if(joboutptp == 2)    // output is backward product
   {
      if(joboutidx == 0)
         *outidx = bstart[monidx];
      else
         *outidx = bstart[monidx] + bsums[joboutidx-1]*deg1;
   }
   else if(joboutptp == 3)    // output is cross product
   {
      if(joboutidx == 0)
         *outidx = cstart[monidx];
      else
         *outidx = cstart[monidx] + csums[joboutidx-1]*deg1;
   }
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
   int dim, int nbr, int deg, int *fsums, int *bsums, int *csums,
   int *fstart, int *bstart, int *cstart, bool verbose )
{ 
   for(int i=0; i<jobs.get_layer_count(layer); i++)
      job_indices(jobs.get_job(layer,i),&inp1ix[i],&inp2ix[i],&outidx[i],
                  dim,nbr,deg,fsums,bsums,csums,fstart,bstart,cstart,verbose);
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
      cout << "The total coefficient count : " << totalcff << endl;
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

   double *data_h = new double[totalcff];       // linearized data
   int ix = 0;
   for(int i=0; i<nbr; i++)
      for(int j=0; j<deg1; j++) data_h[ix++] = cff[i][j];
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++) data_h[ix++] = input[i][j];

   double *data_d;                             // device data
   size_t szdata = totalcff*sizeof(double);
   cudaMalloc((void**)&data_d,szdata);
   cudaMemcpy(data_d,data_h,szdata,cudaMemcpyHostToDevice);

   for(int k=0; k<jobs.get_depth(); k++)
   {
      const int jobnbr = jobs.get_layer_count(k);
      int *in1ix = new int[jobnbr];
      int *in2ix = new int[jobnbr];
      int *outix = new int[jobnbr];

      if(verbose) cout << "preparing jobs at layer " << k << " ..." << endl;

      jobs_coordinates(jobs,k,in1ix,in2ix,outix,dim,nbr,deg,
                       fsums,bsums,csums,fstart,bstart,cstart,verbose);

      free(in1ix); free(in2ix); free(outix);
   }
}
