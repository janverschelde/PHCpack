// The file dbl2_systems_kernels.cu defines the functions with prototypes in
// the file dbl2_systems_kernels.h.

#include <iostream>
#include <iomanip>
#ifdef winwalltime
#include "gettimeofday4win.h"
#else
#include <sys/time.h>
#endif
#include "job_coordinates.h"
#include "write_job_counts.h"
#include "write_gpu_timings.h"
#include "double_double_functions.h"
#include "dbl2_convolutions_host.h"
#include "dbl2_monomials_host.h"
#include "dbl2_polynomials_kernels.h"

using namespace std;

void write_dbl2_cnvflops
 ( int dim, int deg, int ctype,
   ConvolutionJobs cnvjobs, double kernms, double wallsec )
{
   long long int addcnt,mulcnt,flopcnt;

   convolution_operation_counts(deg,cnvjobs,&addcnt,&mulcnt,1);

   if(ctype == 0)
      flopcnt = 20*addcnt + 23*mulcnt;
   else
      flopcnt = 4*20*addcnt + 4*23*mulcnt;
   /*
      1 complex addition takes 2 floating-point (fp) additions
      1 complex multiplication takes 2 fp additions and 4 fp multiplications
   => quadruple the number of fp additions and multiplications */

   long long int bytecnt;

   if(ctype == 0)
      bytecnt = 2*dim*(deg+1);
   else
      bytecnt = 4*dim*(deg+1);

   cout << "    Total number of bytes : " << bytecnt << endl;

   const double intensity = ((double) flopcnt)/bytecnt;
   cout << "     Arithmetic intensity : "
        << scientific << setprecision(3) << intensity
        << " #flops/#bytes" << endl;

   const double kernflops = 1000.0*((double) flopcnt)/kernms;
   const double wallflops = ((double) flopcnt)/wallsec;
   const int gigacnt = pow(2.0,30);

   cout << "Kernel Time Flops : "
        << scientific << setprecision(3) << kernflops;
   cout << fixed << setprecision(3)
        << " = " << kernflops/gigacnt << " Gigaflops" << endl;
   cout << " Wall Clock Flops : "
        << scientific << setprecision(3) << wallflops;
   cout << fixed << setprecision(3)
        << " = " << wallflops/gigacnt << " Gigaflops" << endl;
}

void write_vectorized2_cnvincflops
 ( int dim, int deg,
   ComplexConvolutionJobs cnvjobs, ComplexIncrementJobs incjobs,
   double kernms, double wallsec )
{
   long long int addcnt1,addcnt2,mulcnt,flopcnt;

   complexconv_operation_counts(deg,cnvjobs,&addcnt1,&mulcnt,1);
   complexinc_operation_counts(deg,incjobs,&addcnt2,1);

   flopcnt = 20*(addcnt1 + addcnt2) + 23*mulcnt;
   /*
      1 complex addition takes 2 floating-point (fp) additions
      1 complex multiplication takes 2 fp additions and 4 fp multiplications
   => quadruple the number of fp additions and multiplications */

   long long int bytecnt = 4*dim*(deg+1);

   cout << "    Total number of bytes : " << bytecnt << endl;

   const double intensity = ((double) flopcnt)/bytecnt;
   cout << "     Arithmetic intensity : "
        << scientific << setprecision(3) << intensity
        << " #flops/#bytes" << endl;

   const double kernflops = 1000.0*((double) flopcnt)/kernms;
   const double wallflops = ((double) flopcnt)/wallsec;
   const int gigacnt = pow(2.0,30);

   cout << "Kernel Time Flops : "
        << scientific << setprecision(3) << kernflops;
   cout << fixed << setprecision(3)
        << " = " << kernflops/gigacnt << " Gigaflops" << endl;
   cout << " Wall Clock Flops : "
        << scientific << setprecision(3) << wallflops;
   cout << fixed << setprecision(3)
        << " = " << wallflops/gigacnt << " Gigaflops" << endl;
}

// The code in dbl2_evaldiffdata_to_output is an adaptation of the
// function dbl2_convoluted_data_to_output in dbl2_polynomials_kernels.cu.

void dbl2_evaldiffdata_to_output
 ( double *datahi, double *datalo, double ***outputhi, double ***outputlo,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart, int vrblvl )
{
   const int deg1 = deg+1;
   int ix0,ix1,ix2;

   for(int k=0; k<nbr; k++)
   {
      ix1 = fstart[k] + (nvr[k]-1)*deg1;
      
      if(vrblvl > 1)
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
      else if(nvr[k] > 1)
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

               if(vrblvl > 1)
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

void cmplx2_evaldiffdata_to_output
 ( double *datarehi, double *datarelo, double *dataimhi, double *dataimlo,
   double ***outputrehi, double ***outputrelo,
   double ***outputimhi, double ***outputimlo,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart, int vrblvl )
{
   const int deg1 = deg+1;
   int ix0,ix1,ix2;

   for(int k=0; k<nbr; k++)
   {
      ix1 = fstart[k] + (nvr[k]-1)*deg1;
      
      if(vrblvl > 1)
         cout << "monomial " << k << " update starts at " << ix1 << endl;

      for(int i=0; i<=deg; i++)
      {
         outputrehi[k][dim][i] = datarehi[ix1];
         outputrelo[k][dim][i] = datarelo[ix1];
         outputimhi[k][dim][i] = dataimhi[ix1];
         outputimlo[k][dim][i] = dataimlo[ix1++];
      }
      ix0 = idx[k][0];
      if(nvr[k] == 1)
      {
         ix1 = (1 + k)*deg1;
            
         for(int i=0; i<=deg; i++)
         {
            outputrehi[k][ix0][i] = datarehi[ix1];
            outputrelo[k][ix0][i] = datarelo[ix1];
            outputimhi[k][ix0][i] = dataimhi[ix1];
            outputimlo[k][ix0][i] = dataimlo[ix1++];
         }
      }
      else if(nvr[k] > 1)
      {                               // first and last derivative
         ix2 = nvr[k]-3;
         if(ix2 < 0) ix2 = 0;
         ix1 = bstart[k] + ix2*deg1;

         for(int i=0; i<=deg; i++)
         {
            outputrehi[k][ix0][i] = datarehi[ix1];
            outputrelo[k][ix0][i] = datarelo[ix1];
            outputimhi[k][ix0][i] = dataimhi[ix1];
            outputimlo[k][ix0][i] = dataimlo[ix1++];
         }
         ix2 = nvr[k]-2;
         ix1 = fstart[k] + ix2*deg1;
         ix0 = idx[k][ix2+1];

         for(int i=0; i<=deg; i++)
         {
            outputrehi[k][ix0][i] = datarehi[ix1];
            outputrelo[k][ix0][i] = datarelo[ix1];
            outputimhi[k][ix0][i] = dataimhi[ix1];
            outputimlo[k][ix0][i] = dataimlo[ix1++];
         }
         if(nvr[k] > 2)                   // all other derivatives
         {
            for(int j=1; j<nvr[k]-1; j++)
            {
               ix0 = idx[k][j];            // j-th variable in monomial k
               ix1 = cstart[k] + (j-1)*deg1;

               if(vrblvl > 1)
                  cout << "monomial " << k << " derivative " << ix0
                       << " update starts at " << ix1 << endl;

               for(int i=0; i<=deg; i++)
               {
                  outputrehi[k][ix0][i] = datarehi[ix1];
                  outputrelo[k][ix0][i] = datarelo[ix1];
                  outputimhi[k][ix0][i] = dataimhi[ix1];
                  outputimlo[k][ix0][i] = dataimlo[ix1++];
               }
            }
         }
      }
   }
}

void cmplx2vectorized_evaldiffdata_to_output
 ( double *datarihi, double *datarilo,
   double ***outputrehi, double ***outputrelo,
   double ***outputimhi, double ***outputimlo,
   int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart,
   int totalcff, int offsetri, int vrblvl )
{
   const int deg1 = deg+1;
   const int totcffoffset = totalcff + offsetri;
   int ix0,ix2;
   int ix1re,ix2im;

   for(int k=0; k<nbr; k++)
   {
      ix1re = fstart[k] + (nvr[k]-1)*deg1;
      ix2im = fstart[k] + (nvr[k]-1)*deg1 + totcffoffset;
      
      if(vrblvl > 1)
         cout << "monomial " << k << " update starts at "
              << ix1re << ", " << ix2im << endl;

      for(int i=0; i<=deg; i++)
      {
         outputrehi[k][dim][i] = datarihi[ix1re];
         outputrelo[k][dim][i] = datarilo[ix1re++];
         outputimhi[k][dim][i] = datarihi[ix2im];
         outputimlo[k][dim][i] = datarilo[ix2im++];
      }
      ix0 = idx[k][0];

      if(nvr[k] == 1)
      {
         ix1re = (1 + k)*deg1;
         ix2im = (1 + k)*deg1 + totcffoffset;
            
         for(int i=0; i<=deg; i++)
         {
            outputrehi[k][ix0][i] = datarihi[ix1re];
            outputrelo[k][ix0][i] = datarilo[ix1re++];
            outputimhi[k][ix0][i] = datarihi[ix2im];
            outputimlo[k][ix0][i] = datarilo[ix2im++];
         }
      }
      else if(nvr[k] > 1)
      {                               // first and last derivative
         // ix2 = nvr[k]-3;
         // if(ix2 < 0) ix2 = 0;
         ix2 = nvr[k] - 2;            // complex vectorized is different ...
         ix1re = bstart[k] + ix2*deg1;
         ix2im = bstart[k] + ix2*deg1 + totcffoffset;

         for(int i=0; i<=deg; i++)
         {
            outputrehi[k][ix0][i] = datarihi[ix1re];
            outputrelo[k][ix0][i] = datarilo[ix1re++];
            outputimhi[k][ix0][i] = datarihi[ix2im];
            outputimlo[k][ix0][i] = datarilo[ix2im++];
         }
         ix2 = nvr[k]-2;
         ix1re = fstart[k] + ix2*deg1;
         ix2im = fstart[k] + ix2*deg1 + totcffoffset;
         ix0 = idx[k][ix2+1];

         for(int i=0; i<=deg; i++)
         {
            outputrehi[k][ix0][i] = datarihi[ix1re];
            outputrelo[k][ix0][i] = datarilo[ix1re++];
            outputimhi[k][ix0][i] = datarihi[ix2im];
            outputimlo[k][ix0][i] = datarilo[ix2im++];
         }
         if(nvr[k] > 2)                   // all other derivatives
         {
            for(int j=1; j<nvr[k]-1; j++)
            {
               ix0 = idx[k][j];            // j-th variable in monomial k
               ix1re = cstart[k] + (j-1)*deg1;
               ix2im = cstart[k] + (j-1)*deg1 + totcffoffset;

               if(vrblvl > 1)
                  cout << "monomial " << k << " derivative " << ix0
                       << " update starts at "
                       << ix1re << ", " << ix2im << endl;

               for(int i=0; i<=deg; i++)
               {
                  outputrehi[k][ix0][i] = datarihi[ix1re];
                  outputrelo[k][ix0][i] = datarilo[ix1re++];
                  outputimhi[k][ix0][i] = datarihi[ix2im];
                  outputimlo[k][ix0][i] = datarilo[ix2im++];
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
   double *cnvlapms, double *elapsedms, double *walltimesec, int vrblvl )
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

   if(vrblvl > 1)
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
   bool verbose = (vrblvl > 1);

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

         cudaFree(in1ix_d); cudaFree(in2ix_d); cudaFree(outix_d);
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
       fstart,bstart,cstart,vrblvl);

   if(vrblvl > 0)
   {
      write_GPU_timings(*cnvlapms,0.0,*elapsedms,*walltimesec);
      write_dbl2_cnvflops(dim,deg,0,cnvjobs,*elapsedms,*walltimesec);
   }
   cudaFree(datahi_d); cudaFree(datalo_d);
   free(datahi_h); free(datalo_h);

   free(fstart); free(bstart); free(cstart);
   free(fsums); free(bsums); free(csums);
}

void GPU_cmplx2_mon_evaldiff
 ( int szt, int dim, int nbr, int deg, int *nvr, int **idx,
   double **cffrehi, double **cffrelo, double **cffimhi, double **cffimlo,
   double **inputrehi, double **inputrelo,
   double **inputimhi, double **inputimlo,
   double ***outputrehi, double ***outputrelo,
   double ***outputimhi, double ***outputimlo, ConvolutionJobs cnvjobs,
   double *cnvlapms, double *elapsedms, double *walltimesec, int vrblvl )
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

   if(vrblvl > 1)
      write_coefficient_indices
         (totalcff,nbr,fsums,fstart,bsums,bstart,csums,cstart);

   double *datarehi_h = new double[totalcff];   // real data on host
   double *datarelo_h = new double[totalcff];
   double *dataimhi_h = new double[totalcff];   // imaginary data on host
   double *dataimlo_h = new double[totalcff];
   int ix = 0;
   for(int i=0; i<deg1; i++)
   {
      datarehi_h[ix]   = 0.0; // cst[i]; no constant
      datarelo_h[ix]   = 0.0;
      dataimhi_h[ix]   = 0.0;
      dataimlo_h[ix++] = 0.0;
   }
   for(int i=0; i<nbr; i++)
      for(int j=0; j<deg1; j++)
      {
         datarehi_h[ix]   = cffrehi[i][j];
         datarelo_h[ix]   = cffrelo[i][j];
         dataimhi_h[ix]   = cffimhi[i][j];
         dataimlo_h[ix++] = cffimlo[i][j];
      }
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++)
      {
         datarehi_h[ix]   = inputrehi[i][j];
         datarelo_h[ix]   = inputrelo[i][j];
         dataimhi_h[ix]   = inputimhi[i][j];
         dataimlo_h[ix++] = inputimlo[i][j];
      }

   double *datarehi_d;                        // device real data
   double *datarelo_d;
   double *dataimhi_d;                        // device imag data
   double *dataimlo_d;
   const size_t szdata = totalcff*sizeof(double);
   cudaMalloc((void**)&datarehi_d,szdata);
   cudaMalloc((void**)&datarelo_d,szdata);
   cudaMalloc((void**)&dataimhi_d,szdata);
   cudaMalloc((void**)&dataimlo_d,szdata);
   cudaMemcpy(datarehi_d,datarehi_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datarelo_d,datarelo_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(dataimhi_d,dataimhi_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(dataimlo_d,dataimlo_h,szdata,cudaMemcpyHostToDevice);

   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   *cnvlapms = 0.0;
   float milliseconds;
   struct timeval begintime,endtime; // wall clock time of computations
   bool verbose = (vrblvl > 1);

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
         cmplx2_padded_convjobs<<<jobnbr,deg1>>>
            (datarehi_d,datarelo_d,dataimhi_d,dataimlo_d,
             in1ix_d,in2ix_d,outix_d,deg1);
         cudaEventRecord(stop);
         cudaEventSynchronize(stop);
         cudaEventElapsedTime(&milliseconds,start,stop);
         *cnvlapms += milliseconds;

         cudaFree(in1ix_d); cudaFree(in2ix_d); cudaFree(outix_d);
      }
      free(in1ix_h); free(in2ix_h); free(outix_h);
   }
   gettimeofday(&endtime,0);
   cudaMemcpy(datarehi_h,datarehi_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datarelo_h,datarelo_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(dataimhi_h,dataimhi_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(dataimlo_h,dataimlo_d,szdata,cudaMemcpyDeviceToHost);
   *elapsedms = *cnvlapms;
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   *walltimesec = seconds + microseconds*1.0e-6;

   cmplx2_evaldiffdata_to_output
      (datarehi_h,datarelo_h,dataimhi_h,dataimlo_h,
       outputrehi,outputrelo,outputimhi,outputimlo,
       dim,nbr,deg,nvr,idx,fstart,bstart,cstart,vrblvl);

   if(vrblvl > 0)
   {
      write_GPU_timings(*cnvlapms,0.0,*elapsedms,*walltimesec);
      write_dbl2_cnvflops(dim,deg,1,cnvjobs,*elapsedms,*walltimesec);
   }
   cudaFree(datarehi_d); cudaFree(datarelo_d);
   cudaFree(dataimhi_d); cudaFree(dataimlo_d);
   free(datarehi_h); free(datarelo_h);
   free(dataimhi_h); free(dataimlo_h);

   free(fstart); free(bstart); free(cstart);
   free(fsums); free(bsums); free(csums);
}

void GPU_cmplx2vectorized_mon_evaldiff
 ( int szt, int dim, int nbr, int deg, int *nvr, int **idx,
   double **cffrehi, double **cffrelo, double **cffimhi, double **cffimlo,
   double **inputrehi, double **inputrelo,
   double **inputimhi, double **inputimlo,
   double ***outputrehi, double ***outputrelo,
   double ***outputimhi, double ***outputimlo,
   ComplexConvolutionJobs cnvjobs, ComplexIncrementJobs incjobs,
   double *cnvlapms, double *elapsedms, double *walltimesec, int vrblvl )
{
   const int deg1 = deg+1;
   const int totalcff = complex_coefficient_count(dim,nbr,deg,nvr);
   const int diminput = (1 + nbr + dim)*deg1; // dimension of input
   const int offsetri = totalcff - diminput;  // offset for re/im operands
   const int cmplxtotcff = 2*(totalcff + offsetri);

   int *fstart = new int[nbr];
   int *bstart = new int[nbr];
   int *cstart = new int[nbr];
   int *fsums = new int[nbr];
   int *bsums = new int[nbr];
   int *csums = new int[nbr];

   // code is adjusted from the complex vectorized versions 
   // in dbl_polynomials_kernels.cu

   complex_coefficient_indices
      (dim,nbr,deg,nvr,fsums,bsums,csums,fstart,bstart,cstart);

   if(vrblvl > 1)
   {
      cout << "        total count : " << totalcff << endl;
      cout << "offset for operands : " << offsetri << endl;
      cout << "complex total count : " << cmplxtotcff << endl;
      write_coefficient_indices
         (totalcff,nbr,fsums,fstart,bsums,bstart,csums,cstart);
   }
   double *datarihi_h = new double[cmplxtotcff];   // data on host
   double *datarilo_h = new double[cmplxtotcff];

   int ix1 = 0;
   int ix2 = totalcff + offsetri;

   for(int i=0; i<deg1; i++)
   {
      datarihi_h[ix1]   = 0.0; // cst[i]; no constant
      datarilo_h[ix1++] = 0.0;
      datarihi_h[ix2]   = 0.0;
      datarilo_h[ix2++] = 0.0;
   }
   for(int i=0; i<nbr; i++)
      for(int j=0; j<deg1; j++)
      {
         datarihi_h[ix1]   = cffrehi[i][j];
         datarilo_h[ix1++] = cffrelo[i][j];
         datarihi_h[ix2]   = cffimhi[i][j];
         datarilo_h[ix2++] = cffimlo[i][j];
      }
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++)
      {
         datarihi_h[ix1]   = inputrehi[i][j];
         datarilo_h[ix1++] = inputrelo[i][j];
         datarihi_h[ix2]   = inputimhi[i][j];
         datarilo_h[ix2++] = inputimlo[i][j];
      }

   for(int i=0; i<2*offsetri; i++)
   {
      datarihi_h[ix1]   = 0.0;
      datarilo_h[ix1++] = 0.0;
      datarihi_h[ix2]   = 0.0;
      datarilo_h[ix2++] = 0.0;
   }

   double *datarihi_d;                        // device data
   double *datarilo_d;
   const size_t szdata = cmplxtotcff*sizeof(double);
   cudaMalloc((void**)&datarihi_d,szdata);
   cudaMalloc((void**)&datarilo_d,szdata);
   cudaMemcpy(datarihi_d,datarihi_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datarilo_d,datarilo_h,szdata,cudaMemcpyHostToDevice);

   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   *cnvlapms = 0.0;
   float milliseconds;
   double fliplapms;
   struct timeval begintime,endtime; // wall clock time of computations
   bool verbose = (vrblvl > 1);

   gettimeofday(&begintime,0);
   for(int k=0; k<cnvjobs.get_depth(); k++)
   {
      int jobnbr = cnvjobs.get_layer_count(k);
      int *in1ix_h = new int[jobnbr];
      int *in2ix_h = new int[jobnbr];
      int *outix_h = new int[jobnbr];

      if(verbose) cout << "preparing convolution jobs at layer "
                       << k << " ..." << endl;

      complex_convjobs_coordinates
         (cnvjobs,k,in1ix_h,in2ix_h,outix_h,dim,nbr,deg,nvr,totalcff,offsetri,
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
                 << " threads for convolutions ..." << endl;
         
         cudaEventRecord(start);
         dbl2_padded_convjobs<<<jobnbr,deg1>>>
            (datarihi_d,datarilo_d,in1ix_d,in2ix_d,outix_d,deg1);
         cudaEventRecord(stop);
         cudaEventSynchronize(stop);
         cudaEventElapsedTime(&milliseconds,start,stop);
         *cnvlapms += milliseconds;

         cudaFree(in1ix_d); cudaFree(in2ix_d); cudaFree(outix_d);
      }
      // code to flip signs and run the increment jobs
      jobnbr = incjobs.get_layer_count(k);
      // note: only half the number of increment jobs

      if(vrblvl > 1) cout << "preparing increment jobs at layer "
                          << k << " ..." << endl;

      complex_incjobs_coordinates
         (incjobs,k,in1ix_h,in2ix_h,outix_h,dim,nbr,deg,nvr,totalcff,offsetri,
          fstart,bstart,cstart,verbose);

      const int nbrflips = jobnbr/2;
      int *rebidx = new int[nbrflips];
      for(int i=0, j=0; i<jobnbr; i=i+2, j++) rebidx[j] = in2ix_h[i];

      GPU_cmplx2vectorized_flipsigns
         (deg,nbrflips,rebidx,datarihi_d,datarilo_d,&fliplapms,vrblvl);
      *cnvlapms += fliplapms;

      // if(BS == deg1)
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
                 << " threads for increments ..." << endl;

         cudaEventRecord(start);
         dbl2_increment_jobs<<<jobnbr,deg1>>>
            (datarihi_d,datarilo_d,in1ix_d,in2ix_d,outix_d,deg1);
         cudaEventRecord(stop);
         cudaEventSynchronize(stop);
         cudaEventElapsedTime(&milliseconds,start,stop);
         *cnvlapms += milliseconds;
         cudaFree(in1ix_d); cudaFree(in2ix_d); cudaFree(outix_d);
      }
      free(in1ix_h); free(in2ix_h); free(outix_h); free(rebidx);
   }
   gettimeofday(&endtime,0);
   cudaMemcpy(datarihi_h,datarihi_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datarilo_h,datarilo_d,szdata,cudaMemcpyDeviceToHost);
   *elapsedms = *cnvlapms;
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   *walltimesec = seconds + microseconds*1.0e-6;

   cmplx2vectorized_evaldiffdata_to_output
      (datarihi_h,datarilo_h,outputrehi,outputrelo,outputimhi,outputimlo,
       dim,nbr,deg,nvr,idx,fstart,bstart,cstart,totalcff,offsetri,vrblvl);

   if(vrblvl > 0)
   {
      write_GPU_timings(*cnvlapms,0.0,*elapsedms,*walltimesec);
      write_vectorized2_cnvincflops
         (dim,deg,cnvjobs,incjobs,*elapsedms,*walltimesec);
   }
   cudaFree(datarihi_d); cudaFree(datarilo_d);
   free(datarihi_h); free(datarilo_h);

   free(fstart); free(bstart); free(cstart);
   free(fsums); free(bsums); free(csums);
}

void GPU_dbl2_evaluate_monomials
 ( int dim, int deg, int szt, int nbt,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffhi, double **cfflo, double *acchi, double *acclo,
   double **inputhi, double **inputlo, double ***outputhi, double ***outputlo,
   double *totcnvlapsedms, int vrblvl )
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
   if(vrblvl > 1)
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

   bool verbose = (vrblvl > 1);
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
       &cnvlapms,&elapsedms,&walltimesec,vrblvl);

   *totcnvlapsedms += elapsedms;

   if(vrblvl > 1)
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
               double acchi,acclo;

               // multiply derivative w.r.t. idxvar with factor
               for(int k=0; k<=deg; k++)
               {
                  // output[i][idxvar][k] = factor*output[i][idxvar][k];
                  ddf_mul(outputhi[i][idxvar][k],outputlo[i][idxvar][k],
                          factor,0.0,&acchi,&acclo);
                  outputhi[i][idxvar][k] = acchi;
                  outputlo[i][idxvar][k] = acclo;
               }
            }
         }
      }
   }
}

void GPU_cmplx2_evaluate_monomials
 ( int dim, int deg, int szt, int nbt,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffrehi, double **cffrelo, double **cffimhi, double **cffimlo,
   double *accrehi, double *accrelo, double *accimhi, double *accimlo,
   double **inputrehi, double **inputrelo,
   double **inputimhi, double **inputimlo, 
   double ***outputrehi, double ***outputrelo, 
   double ***outputimhi, double ***outputimlo,
   double *totcnvlapsedms, int vrblvl )
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
                  CPU_cmplx2_product
                     (deg,inputrehi[idxvar],inputrelo[idxvar],
                          inputimhi[idxvar],inputimlo[idxvar],
                      cffrehi[i],cffrelo[i],cffimhi[i],cffimlo[i],
                      accrehi,accrelo,accimhi,accimlo);

                  for(int L=0; L<=deg; L++)
                  {
                     cffrehi[i][L] = accrehi[L]; cffrelo[i][L] = accrelo[L];
                     cffimhi[i][L] = accimhi[L]; cffimlo[i][L] = accimlo[L];
                  }
               }
            }
         }
      }
   }
   if(vrblvl > 1)
   {
      for(int i=0; i<dim; i++)
      {
         cout << "coefficients for monomial " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << cffrehi[i][j] << "  " << cffrelo[i][j] << endl << "  "
                 << cffimhi[i][j] << "  " << cffimlo[i][j] << endl;
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
            cout << inputrehi[i][j] << "  " << inputrelo[i][j] << endl << "  "
                 << inputimhi[i][j] << "  " << inputimlo[i][j] << endl;
      }
   }
   bool verbose = (vrblvl > 1);
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
   GPU_cmplx2_mon_evaldiff
      (szt,dim,dim,deg,nvr,idx,cffrehi,cffrelo,cffimhi,cffimlo,
       inputrehi,inputrelo,inputimhi,inputimlo,
       outputrehi,outputrelo,outputimhi,outputimlo,jobs,
       &cnvlapms,&elapsedms,&walltimesec,vrblvl);

   *totcnvlapsedms += elapsedms;

   if(vrblvl > 1)
   {
      for(int i=0; i<dim; i++)
      {
         cout << "output series for monomial " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << outputrehi[i][dim][j] << "  "
                 << outputrelo[i][dim][j] << endl << "  " 
                 << outputimhi[i][dim][j] << "  "
                 << outputimlo[i][dim][j] << endl;
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
               double fac = (double) exp[i][j];
               double acchi,acclo;

               // multiply derivative w.r.t. idxvar with factor
               for(int k=0; k<=deg; k++)
               {
                  // outputre[i][idxvar][k] = fac*outputre[i][idxvar][k];
                  ddf_mul(outputrehi[i][idxvar][k],
                          outputrelo[i][idxvar][k],fac,0.0,&acchi,&acclo);
                  outputrehi[i][idxvar][k] = acchi;
                  outputrelo[i][idxvar][k] = acclo;
                  // outputim[i][idxvar][k] = fac*outputim[i][idxvar][k];
                  ddf_mul(outputimhi[i][idxvar][k],
                          outputimlo[i][idxvar][k],fac,0.0,&acchi,&acclo);
                  outputimhi[i][idxvar][k] = acchi;
                  outputimlo[i][idxvar][k] = acclo;
               }
            }
         }
      }
   }
}

void GPU_cmplx2vectorized_evaluate_monomials
 ( int dim, int deg, int szt, int nbt,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffrehi, double **cffrelo, double **cffimhi, double **cffimlo,
   double *accrehi, double *accrelo, double *accimhi, double *accimlo,
   double **inputrehi, double **inputrelo,
   double **inputimhi, double **inputimlo, 
   double ***outputrehi, double ***outputrelo, 
   double ***outputimhi, double ***outputimlo,
   double *totcnvlapsedms, int vrblvl )
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
                  CPU_cmplx2_product
                     (deg,inputrehi[idxvar],inputrelo[idxvar],
                          inputimhi[idxvar],inputimlo[idxvar],
                      cffrehi[i],cffrelo[i],cffimhi[i],cffimlo[i],
                      accrehi,accrelo,accimhi,accimlo);

                  for(int L=0; L<=deg; L++)
                  {
                     cffrehi[i][L] = accrehi[L]; cffrelo[i][L] = accrelo[L];
                     cffimhi[i][L] = accimhi[L]; cffimlo[i][L] = accimlo[L];
                  }
               }
            }
         }
      }
   }
   if(vrblvl > 1)
   {
      for(int i=0; i<dim; i++)
      {
         cout << "coefficients for monomial " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << cffrehi[i][j] << "  " << cffrelo[i][j] << endl << "  "
                 << cffimhi[i][j] << "  " << cffimlo[i][j] << endl;
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
            cout << inputrehi[i][j] << "  " << inputrelo[i][j] << endl << "  "
                 << inputimhi[i][j] << "  " << inputimlo[i][j] << endl;
      }
   }
   bool verbose = (vrblvl > 1);
   double cnvlapms,elapsedms,walltimesec;

   // ConvolutionJobs jobs(dim);
   // jobs.make(dim,nvr,idx,verbose);
   ComplexConvolutionJobs cnvjobs(dim);
   cnvjobs.make(dim,nvr,idx,verbose);

   ComplexIncrementJobs incjobs(cnvjobs,verbose);

   if(verbose)
   {
      for(int k=0; k<cnvjobs.get_depth(); k++)
      {
         cout << "jobs at layer " << k << " :" << endl;
         for(int i=0; i<cnvjobs.get_layer_count(k); i++)
            cout << cnvjobs.get_job(k,i) << endl;
      }
      cout << "dimension : " << dim << endl;
      cout << "number of monomials : " << dim << endl;
      cout << "number of convolution jobs : "
           << cnvjobs.get_count() << endl;
      cout << "number of layers : " << cnvjobs.get_depth() << endl;
      cout << "frequency of layer counts :" << endl;
      int checksum = 0;
      for(int i=0; i<cnvjobs.get_depth(); i++)
      {
         cout << i << " : " << cnvjobs.get_layer_count(i) << endl;
         checksum = checksum + cnvjobs.get_layer_count(i); 
      }
      cout << "layer count sum : " << checksum << endl;
   }
   GPU_cmplx2vectorized_mon_evaldiff
      (szt,dim,dim,deg,nvr,idx,cffrehi,cffrelo,cffimhi,cffimlo,
       inputrehi,inputrelo,inputimhi,inputimlo,
       outputrehi,outputrelo,outputimhi,outputimlo,cnvjobs,incjobs,
       &cnvlapms,&elapsedms,&walltimesec,vrblvl);

   *totcnvlapsedms += elapsedms;

   if(vrblvl > 1)
   {
      for(int i=0; i<dim; i++)
      {
         cout << "output series for monomial " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << outputrehi[i][dim][j] << "  "
                 << outputrelo[i][dim][j] << endl << "  " 
                 << outputimhi[i][dim][j] << "  "
                 << outputimlo[i][dim][j] << endl;
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
               double fac = (double) exp[i][j];
               double acchi,acclo;

               // multiply derivative w.r.t. idxvar with factor
               for(int k=0; k<=deg; k++)
               {
                  // outputre[i][idxvar][k] = fac*outputre[i][idxvar][k];
                  ddf_mul(outputrehi[i][idxvar][k],
                          outputrelo[i][idxvar][k],fac,0.0,&acchi,&acclo);
                  outputrehi[i][idxvar][k] = acchi;
                  outputrelo[i][idxvar][k] = acclo;
                  // outputim[i][idxvar][k] = fac*outputim[i][idxvar][k];
                  ddf_mul(outputimhi[i][idxvar][k],
                          outputimlo[i][idxvar][k],fac,0.0,&acchi,&acclo);
                  outputimhi[i][idxvar][k] = acchi;
                  outputimlo[i][idxvar][k] = acclo;
               }
            }
         }
      }
   }
}

void GPU_dbl2_evaluate_columns
 ( int dim, int deg, int nbrcol, int szt, int nbt, int **nvr, int ***idx,
   double ***cffhi, double ***cfflo, double **inputhi, double **inputlo, 
   double ***outputhi, double ***outputlo,
   double **funvalhi, double **funvallo,
   double ***jacvalhi, double ***jacvallo,
   double *totcnvlapsedms, int vrblvl )
{
   const int degp1 = deg+1;

   for(int i=0; i<dim; i++)
      for(int j=0; j<degp1; j++)
      {
         funvalhi[i][j] = 0.0;
         funvallo[i][j] = 0.0;
      }
   funvalhi[dim-1][0] = -1.0; // constant of last eq in cyclic dim-roots

   for(int k=0; k<degp1; k++)  // the Jacobian is linearized
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
         {
            jacvalhi[k][i][j] = 0.0;
            jacvallo[k][i][j] = 0.0;
         }

   if(vrblvl > 1)
   {
      for(int i=0; i<nbrcol; i++)
      {
         cout << "coefficients for column " << i << " :" << endl;
         for(int j=0; j<dim; j++)
         {
            cout << "coefficients for monomial " << j << " :" << endl;
            for(int k=0; k<=deg; k++)
               cout << cffhi[i][j][k] << "  " << cfflo[i][j][k] << endl;
         }
      }
      cout << "dim : " << dim << endl;
      for(int i=0; i<nbrcol; i++)
      {
         for(int j=0; j<dim; j++)
         {
            cout << "nvr[" << i << "][" << j << "] : " << nvr[i][j] << " :";
            cout << "  idx[" << i << "][" << j << "] :";
            for(int k=0; k<nvr[i][j]; k++) cout << " " << idx[i][j][k];
            cout << endl;
         }
      }
      for(int i=0; i<dim; i++)
      {
         cout << "input series for variable " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << inputhi[i][j] << "  " << inputlo[i][j] << endl;
      }
   }
   bool verbose = (vrblvl > 1);

   double cnvlapms,elapsedms,walltimesec;

   for(int i=0; i<nbrcol; i++)
   {
      ConvolutionJobs jobs(dim);

      jobs.make(dim,nvr[i],idx[i],verbose);

      if(verbose)
      {
         for(int k=0; k<jobs.get_depth(); k++)
         {
            cout << "jobs at layer " << k << " :" << endl;
            for(int j=0; j<jobs.get_layer_count(k); j++)
               cout << jobs.get_job(k,j) << endl;
         }
         cout << "dimension : " << dim << endl;
         cout << "number of monomials : " << dim << endl;
         cout << "number of convolution jobs : " << jobs.get_count() << endl;
         cout << "number of layers : " << jobs.get_depth() << endl;
         cout << "frequency of layer counts :" << endl;
         int checksum = 0;
         for(int j=0; j<jobs.get_depth(); j++)
         {
            cout << j << " : " << jobs.get_layer_count(j) << endl;
            checksum = checksum + jobs.get_layer_count(j); 
         }
         cout << "layer count sum : " << checksum << endl;
      }
      GPU_dbl2_mon_evaldiff
         (szt,dim,dim,deg,nvr[i],idx[i],cffhi[i],cfflo[i],
          inputhi,inputlo,outputhi,outputlo,jobs,
          &cnvlapms,&elapsedms,&walltimesec,vrblvl);

      *totcnvlapsedms += elapsedms;

      for(int j=0; j<dim; j++)
         if(nvr[i][j] > 0)       // update values
         {
            for(int L=0; L<degp1; L++) // funval[j][L] += output[j][dim][L];
            {
               ddf_inc(&funvalhi[j][L],&funvallo[j][L],
                       outputhi[j][dim][L],outputlo[j][dim][L]);
            }

            int *indexes = idx[i][j];      // indices of the variables
            for(int k=0; k<nvr[i][j]; k++) // derivative w.r.t. idx[i][j][k]
            {                              // has j-th coefficient
               int idxval = indexes[k];
               for(int L=0; L<degp1; L++) 
               {
                  // jacval[L][j][idxval] += output[j][idxval][L];
                  ddf_inc(&jacvalhi[L][j][idxval],&jacvallo[L][j][idxval],
                          outputhi[j][idxval][L],outputlo[j][idxval][L]);
               }
            }
         }
   }
   if(vrblvl > 1)
   {
      for(int i=0; i<dim; i++)
      {
         cout << "output series for monomial " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << funvalhi[i][j] << "  " << funvallo[i][j] << endl;
      }
   }
}

void GPU_cmplx2_evaluate_columns
 ( int dim, int deg, int nbrcol, int szt, int nbt, int **nvr, int ***idx, 
   double ***cffrehi, double ***cffrelo, double ***cffimhi, double ***cffimlo, 
   double **inputrehi, double **inputrelo,
   double **inputimhi, double **inputimlo,
   double ***outputrehi, double ***outputrelo,
   double ***outputimhi, double ***outputimlo, 
   double **funvalrehi, double **funvalrelo,
   double **funvalimhi, double **funvalimlo,
   double ***jacvalrehi, double ***jacvalrelo,
   double ***jacvalimhi, double ***jacvalimlo,
   double *totcnvlapsedms, int vrblvl )
{
   const int degp1 = deg+1;

   for(int i=0; i<dim; i++)
      for(int j=0; j<degp1; j++)
      {
         funvalrehi[i][j] = 0.0; funvalrelo[i][j] = 0.0;
         funvalimhi[i][j] = 0.0; funvalimlo[i][j] = 0.0;
      }
   funvalrehi[dim-1][0] = -1.0; // constant of last eq in cyclic dim-roots

   for(int k=0; k<degp1; k++)  // the Jacobian is linearized
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
         {
            jacvalrehi[k][i][j] = 0.0; jacvalrelo[k][i][j] = 0.0;
            jacvalimhi[k][i][j] = 0.0; jacvalimlo[k][i][j] = 0.0;
         }

   if(vrblvl > 1)
   {
      for(int i=0; i<nbrcol; i++)
      {
         cout << "coefficients for column " << i << " :" << endl;
         for(int j=0; j<dim; j++)
         {
            cout << "coefficients for monomial " << j << " :" << endl;
            for(int k=0; k<=deg; k++)
               cout << cffrehi[i][j][k] << "  " << cffrelo[i][j][k] << endl
                    << "  "
                    << cffimhi[i][j][k] << "  " << cffimlo[i][j][k] << endl;
         }
      }
      cout << "dim : " << dim << endl;
      for(int i=0; i<nbrcol; i++)
      {
         for(int j=0; j<dim; j++)
         {
            cout << "nvr[" << i << "][" << j << "] : " << nvr[i][j] << " :";
            cout << "  idx[" << i << "][" << j << "] :";
            for(int k=0; k<nvr[i][j]; k++) cout << " " << idx[i][j][k];
            cout << endl;
         }
      }
      for(int i=0; i<dim; i++)
      {
         cout << "input series for variable " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << inputrehi[i][j] << "  " << inputrelo[i][j] << endl
                 << "  "
                 << inputimhi[i][j] << "  " << inputimlo[i][j] << endl;
      }
   }
   bool verbose = (vrblvl > 1);

   double cnvlapms,elapsedms,walltimesec;

   for(int i=0; i<nbrcol; i++)
   {
      ConvolutionJobs jobs(dim);

      jobs.make(dim,nvr[i],idx[i],verbose);

      if(verbose)
      {
         for(int k=0; k<jobs.get_depth(); k++)
         {
            cout << "jobs at layer " << k << " :" << endl;
            for(int j=0; j<jobs.get_layer_count(k); j++)
               cout << jobs.get_job(k,j) << endl;
         }
         cout << "dimension : " << dim << endl;
         cout << "number of monomials : " << dim << endl;
         cout << "number of convolution jobs : " << jobs.get_count() << endl;
         cout << "number of layers : " << jobs.get_depth() << endl;
         cout << "frequency of layer counts :" << endl;
         int checksum = 0;
         for(int j=0; j<jobs.get_depth(); j++)
         {
            cout << j << " : " << jobs.get_layer_count(j) << endl;
            checksum = checksum + jobs.get_layer_count(j); 
         }
         cout << "layer count sum : " << checksum << endl;
      }
      GPU_cmplx2_mon_evaldiff
         (szt,dim,dim,deg,nvr[i],idx[i],
          cffrehi[i],cffrelo[i],cffimhi[i],cffimlo[i],
          inputrehi,inputrelo,inputimhi,inputimlo,
          outputrehi,outputrelo,outputimhi,outputimlo,jobs,
          &cnvlapms,&elapsedms,&walltimesec,vrblvl);

      *totcnvlapsedms += elapsedms;

      for(int j=0; j<dim; j++)
         if(nvr[i][j] > 0)       // update values
         {
            for(int L=0; L<degp1; L++)
            {
               // funvalre[j][L] += outputre[j][dim][L];
               // funvalim[j][L] += outputim[j][dim][L];
               ddf_inc(&funvalrehi[j][L],&funvalrelo[j][L],
                       outputrehi[j][dim][L],outputrelo[j][dim][L]);
               ddf_inc(&funvalimhi[j][L],&funvalimlo[j][L],
                       outputimhi[j][dim][L],outputimlo[j][dim][L]);
            }

            int *indexes = idx[i][j];      // indices of the variables
            for(int k=0; k<nvr[i][j]; k++) // derivative w.r.t. idx[i][j][k]
            {                              // has j-th coefficient
               int idxval = indexes[k];
               for(int L=0; L<degp1; L++) 
               {
                  // jacvalre[L][j][idxval] += outputre[j][idxval][L];
                  // jacvalim[L][j][idxval] += outputim[j][idxval][L];
                  ddf_inc(&jacvalrehi[L][j][idxval],&jacvalrelo[L][j][idxval],
                          outputrehi[j][idxval][L],outputrelo[j][idxval][L]);
                  ddf_inc(&jacvalimhi[L][j][idxval],&jacvalimlo[L][j][idxval],
                          outputimhi[j][idxval][L],outputimlo[j][idxval][L]);
               }
            }
         }
   }
   if(vrblvl > 1)
   {
      for(int i=0; i<dim; i++)
      {
         cout << "output series for polynomial " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << funvalrehi[i][j] << "  " << funvalrelo[i][j] << endl
                 << "  "
                 << funvalimhi[i][j] << "  " << funvalimlo[i][j] << endl;
      }
   }
}
