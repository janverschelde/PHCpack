// The file dbl2_systems_kernels.cu defines the functions with prototypes in
// the file dbl2_systems_kernels.h.

#include <iostream>
#include <iomanip>
#ifdef winwalltime
#include "gettimeofday4win.h"
#else
#include <sys/time.h>
#endif
#include "write_gpu_timings.h"
#include "job_coordinates.h"
#include "dbl2_convolutions_host.h"
#include "dbl2_monomials_host.h"
#include "dbl2_polynomials_kernels.h"

using namespace std;

// The code in dbl2_evaldiffdata_to_output is an adaptation of the
// function dbl2_convoluted_data_to_output in dbl2_polynomials_kernels.cu.

void dbl2_evaldiffdata_to_output
 ( double *datahi, double *datalo, double ***outputhi, double ***outputlo,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart, bool verbose )
{
   const int deg1 = deg+1;
   int ix0,ix1,ix2;

   for(int k=0; k<nbr; k++)
   {
      ix1 = fstart[k] + (nvr[k]-1)*deg1;
      
      if(verbose)
         cout << "monomial " << k << " update starts at " << ix1 << endl;

      for(int i=0; i<=deg; i++)
      {
         outputhi[k][dim][i] = datahi[ix1];
         outputlo[k][dim][i] = datalo[ix1++];
      }
      ix0 = idx[k][0];
      if(nvr[k] == 1)
      {
         ix1 = (1 + k)*deg1;
            
         for(int i=0; i<=deg; i++)
         {
            outputhi[k][ix0][i] = datahi[ix1];
            outputlo[k][ix0][i] = datalo[ix1++];
         }
      }
      else
      {                               // first and last derivative
         ix2 = nvr[k]-3;
         if(ix2 < 0) ix2 = 0;
         ix1 = bstart[k] + ix2*deg1;

         for(int i=0; i<=deg; i++)
         {
            outputhi[k][ix0][i] = datahi[ix1];
            outputlo[k][ix0][i] = datalo[ix1++];
         }
         ix2 = nvr[k]-2;
         ix1 = fstart[k] + ix2*deg1;
         ix0 = idx[k][ix2+1];

         for(int i=0; i<=deg; i++)
         {
            outputhi[k][ix0][i] = datahi[ix1];
            outputlo[k][ix0][i] = datalo[ix1++];
         }
         if(nvr[k] > 2)                   // all other derivatives
         {
            for(int j=1; j<nvr[k]-1; j++)
            {
               ix0 = idx[k][j];            // j-th variable in monomial k
               ix1 = cstart[k] + (j-1)*deg1;

               if(verbose)
                  cout << "monomial " << k << " derivative " << ix0
                       << " update starts at " << ix1 << endl;

               for(int i=0; i<=deg; i++)
               {
                  outputhi[k][ix0][i] = datahi[ix1];
                  outputlo[k][ix0][i] = datalo[ix1++];
               }
            }
         }
      }
   }
}

// The code for GPU_dbl2_mon_evaldiff is an adaptation of the
// function GPU_dbl2_poly_evaldiff of dbl_polynomials_kernels.cu.

void GPU_dbl2_mon_evaldiff
 ( int szt, int dim, int nbr, int deg, int *nvr, int **idx,
   double **cffhi, double **cfflo, double **inputhi, double **inputlo,
   double ***outputhi, double ***outputlo, ConvolutionJobs cnvjobs,
   double *cnvlapms, double *elapsedms, double *walltimesec, bool verbose )
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
   double *datalo_h = new double[totalcff];

   int ix = 0;
   for(int i=0; i<deg1; i++)
   {
      datahi_h[ix]   = 0.0; // cst[i]; no constant
      datalo_h[ix++] = 0.0;
   }
   for(int i=0; i<nbr; i++)
      for(int j=0; j<deg1; j++)
      {
         datahi_h[ix]   = cffhi[i][j];
         datalo_h[ix++] = cfflo[i][j];
      }

   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++)
      {
         datahi_h[ix]   = inputhi[i][j];
         datalo_h[ix++] = inputlo[i][j];
      }

   double *datahi_d;                               // device data
   double *datalo_d;
   const size_t szdata = totalcff*sizeof(double);
   cudaMalloc((void**)&datahi_d,szdata);
   cudaMalloc((void**)&datalo_d,szdata);
   cudaMemcpy(datahi_d,datahi_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datalo_d,datalo_h,szdata,cudaMemcpyHostToDevice);

   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   *cnvlapms = 0.0;
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
      // if(deg1 == szt)
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
            cout << "launching " << jobnbr << " blocks of " << deg1
                 << " threads ..." << endl;
         
         cudaEventRecord(start);
         dbl2_padded_convjobs<<<jobnbr,deg1>>>
            (datahi_d,datalo_d,in1ix_d,in2ix_d,outix_d,deg1);
         cudaEventRecord(stop);
         cudaEventSynchronize(stop);
         cudaEventElapsedTime(&milliseconds,start,stop);
         *cnvlapms += milliseconds;
      }
      free(in1ix_h); free(in2ix_h); free(outix_h);
   }
   gettimeofday(&endtime,0);
   cudaMemcpy(datahi_h,datahi_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datalo_h,datalo_d,szdata,cudaMemcpyDeviceToHost);
   *elapsedms = *cnvlapms;
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   *walltimesec = seconds + microseconds*1.0e-6;

   dbl2_evaldiffdata_to_output
      (datahi_h,datalo_h,outputhi,outputlo,dim,nbr,deg,nvr,idx,
       fstart,bstart,cstart,verbose);

   if(verbose)
      write_GPU_timings(*cnvlapms,0.0,*elapsedms,*walltimesec);
}

void GPU_dbl2_evaluate_monomials
 ( int dim, int deg, int szt, int nbt,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffhi, double **cfflo, double *acchi, double *acclo,
   double **inputhi, double **inputlo, double ***outputhi, double ***outputlo,
   int vrblvl )
{
   for(int i=0; i<dim; i++) // common factors in the coefficients
   {
      if(nbrfac[i] > 0) // there are common factors in monomial i
      {
         for(int j=0; j<nvr[i]; j++) // run over all exponents
         {
            if(expfac[i][j] > 0) // the j-th exponent with variable idx[i][j]
            {
               int idxvar = idx[i][j];

               for(int k=0; k<expfac[i][j]; k++)
               {
                  CPU_dbl2_product(deg,inputhi[idxvar],inputlo[idxvar],
                                   cffhi[i],cfflo[i],acchi,acclo);
                  for(int L=0; L<=deg; L++)
                  {
                     cffhi[i][L] = acchi[L];
                     cfflo[i][L] = acclo[L];
                  }
               }
            }
         }
      }
   }
   if(vrblvl > 0)
   {
      for(int i=0; i<dim; i++)
      {
         cout << "coefficients for monomial " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << cffhi[i][j] << "  " << cfflo[i][j] << endl;
      }
      cout << "dim : " << dim << "  nvr :";
      for(int i=0; i<dim; i++) cout << " " << nvr[i];
      cout << endl;
      cout << "deg : " << deg;
      for(int i=0; i<dim; i++) 
      {
         cout << "  idx[" << i << "] :";
         for(int j=0; j<nvr[i]; j++) cout << " " << idx[i][j];
      }
      cout << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "input series for variable " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << inputhi[i][j] << "  " << inputlo[i][j] << endl;
      }
   }
   // for(int i=0; i<dim; i++)
   //    CPU_dbl_evaldiff(dim,nvr[i],deg,idx[i],cff[i],input,output[i]);

   bool verbose = (vrblvl > 0);
   double cnvlapms,elapsedms,walltimesec;

   ConvolutionJobs jobs(dim);

   jobs.make(dim,nvr,idx,verbose);

   if(verbose)
   {
      for(int k=0; k<jobs.get_depth(); k++)
      {
         cout << "jobs at layer " << k << " :" << endl;
         for(int i=0; i<jobs.get_layer_count(k); i++)
            cout << jobs.get_job(k,i) << endl;
      }
      cout << "dimension : " << dim << endl;
      cout << "number of monomials : " << dim << endl;
      cout << "number of convolution jobs : " << jobs.get_count() << endl;
      cout << "number of layers : " << jobs.get_depth() << endl;
      cout << "frequency of layer counts :" << endl;
      int checksum = 0;
      for(int i=0; i<jobs.get_depth(); i++)
      {
         cout << i << " : " << jobs.get_layer_count(i) << endl;
         checksum = checksum + jobs.get_layer_count(i); 
      }
      cout << "layer count sum : " << checksum << endl;
   }
   GPU_dbl2_mon_evaldiff
      (szt,dim,dim,deg,nvr,idx,
       cffhi,cfflo,inputhi,inputlo,outputhi,outputlo,jobs,
       &cnvlapms,&elapsedms,&walltimesec,verbose);

   if(vrblvl > 0)
   {
      for(int i=0; i<dim; i++)
      {
         cout << "output series for monomial " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << outputhi[i][dim][j] << "  "
                 << outputlo[i][dim][j] << endl;
      }
   }

   for(int i=0; i<dim; i++) // multiply derivatives with the powers
   {
      if(nbrfac[i] > 0) // there are common factors in monomial i
      {
         for(int j=0; j<nvr[i]; j++) // run over all exponents
         {
            if(expfac[i][j] > 0) // the j-th exponent with variable idx[i][j]
            {
               int idxvar = idx[i][j];
               double factor = (double) exp[i][j];

               // multiply derivative w.r.t. idxvar with factor
               for(int k=0; k<=deg; k++)
               {
                  outputhi[i][idxvar][k] = factor*outputhi[i][idxvar][k];
                  outputlo[i][idxvar][k] = factor*outputlo[i][idxvar][k];
               }
            }
         }
      }
   }
}
