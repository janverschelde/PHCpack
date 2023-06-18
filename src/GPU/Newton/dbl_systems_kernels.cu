// The file dbl_systems_kernels.cu defines the functions with prototypes in
// the file dbl_systems_kernels.h.

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
#include "dbl_convolutions_host.h"
#include "dbl_monomials_host.h"
#include "dbl_polynomials_kernels.h"

using namespace std;

void write_arithmetic_intensity
 ( long long int flopcnt, long long int bytecnt,
   double kernms, double wallsec )
{
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

void write_dbl_cnvflops
 ( int dim, int deg, int ctype,
   ConvolutionJobs cnvjobs, double kernms, double wallsec )
{
   long long int addcnt,mulcnt,flopcnt;

   convolution_operation_counts(deg,cnvjobs,&addcnt,&mulcnt,1);

   if(ctype == 0)
      flopcnt = addcnt + mulcnt;
   else
      flopcnt = 4*addcnt + 4*mulcnt;
   /*
      1 complex addition takes 2 floating-point (fp) additions
      1 complex multiplication takes 2 fp additions and 4 fp multiplications
   => quadruple the number of fp additions and multiplications */

   long long int bytecnt;

   if(ctype == 0)
      bytecnt = dim*(deg+1);
   else
      bytecnt = 2*dim*(deg+1);

   cout << "    Total number of bytes : " << bytecnt << endl;
 
   write_arithmetic_intensity(flopcnt,bytecnt,kernms,wallsec);
}

void write_vectorized_cnvincflops
 ( int dim, int deg,
   ComplexConvolutionJobs cnvjobs, ComplexIncrementJobs incjobs,
   double kernms, double wallsec )
{
   long long int addcnt1,addcnt2,mulcnt,flopcnt;

   complexconv_operation_counts(deg,cnvjobs,&addcnt1,&mulcnt,1);
   complexinc_operation_counts(deg,incjobs,&addcnt2,1);

   flopcnt = addcnt1 + addcnt2 + mulcnt;

   long long int bytecnt = 2*dim*(deg+1);

   cout << "    Total number of bytes : " << bytecnt << endl;

   write_arithmetic_intensity(flopcnt,bytecnt,kernms,wallsec);
}

// The code in dbl_evaldiffdata_to_output is an adaptation of the
// function dbl_convoluted_data_to_output in dbl_polynomials_kernels.cu.

void dbl_evaldiffdata_to_output
 ( double *data, double ***output, int dim, int nbr, int deg, int *nvr,
   int **idx, int *fstart, int *bstart, int *cstart, int vrblvl )
{
   const int deg1 = deg+1;
   int ix0,ix1,ix2;

   for(int k=0; k<nbr; k++)
   {
      ix1 = fstart[k] + (nvr[k]-1)*deg1;
      
      if(vrblvl > 1)
         cout << "monomial " << k << " update starts at " << ix1 << endl;

      for(int i=0; i<=deg; i++) output[k][dim][i] = data[ix1++];

      ix0 = idx[k][0];
      if(nvr[k] == 1)
      {
         ix1 = (1 + k)*deg1;
            
         for(int i=0; i<=deg; i++) output[k][ix0][i] = data[ix1++];
      }
      else if(nvr[k] > 1)
      {                               // first and last derivative
         ix2 = nvr[k]-3;
         if(ix2 < 0) ix2 = 0;
         ix1 = bstart[k] + ix2*deg1;

         for(int i=0; i<=deg; i++) output[k][ix0][i] = data[ix1++];

         ix2 = nvr[k]-2;
         ix1 = fstart[k] + ix2*deg1;
         ix0 = idx[k][ix2+1];

         for(int i=0; i<=deg; i++) output[k][ix0][i] = data[ix1++];

         if(nvr[k] > 2)                   // all other derivatives
         {
            for(int j=1; j<nvr[k]-1; j++)
            {
               ix0 = idx[k][j];            // j-th variable in monomial k
               ix1 = cstart[k] + (j-1)*deg1;

               if(vrblvl > 1)
                  cout << "monomial " << k << " derivative " << ix0
                       << " update starts at " << ix1 << endl;

               for(int i=0; i<=deg; i++) output[k][ix0][i] = data[ix1++];
            }
         }
      }
   }
}

void cmplx_evaldiffdata_to_output
 ( double *datare, double *dataim, double ***outputre, double ***outputim,
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
         outputre[k][dim][i] = datare[ix1];
         outputim[k][dim][i] = dataim[ix1++];
      }
      ix0 = idx[k][0];
      if(nvr[k] == 1)
      {
         ix1 = (1 + k)*deg1;
            
         for(int i=0; i<=deg; i++)
         {
            outputre[k][ix0][i] = datare[ix1];
            outputim[k][ix0][i] = dataim[ix1++];
         }
      }
      else if(nvr[k] > 1)
      {                               // first and last derivative
         ix2 = nvr[k]-3;
         if(ix2 < 0) ix2 = 0;
         ix1 = bstart[k] + ix2*deg1;

         for(int i=0; i<=deg; i++)
         {
            outputre[k][ix0][i] = datare[ix1];
            outputim[k][ix0][i] = dataim[ix1++];
         }
         ix2 = nvr[k]-2;
         ix1 = fstart[k] + ix2*deg1;
         ix0 = idx[k][ix2+1];

         for(int i=0; i<=deg; i++)
         {
            outputre[k][ix0][i] = datare[ix1];
            outputim[k][ix0][i] = dataim[ix1++];
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
                  outputre[k][ix0][i] = datare[ix1];
                  outputim[k][ix0][i] = dataim[ix1++];
               }
            }
         }
      }
   }
}

void cmplxvectorized_evaldiffdata_to_output
 ( double *datari, double ***outputre, double ***outputim,
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
         outputre[k][dim][i] = datari[ix1re++];
         outputim[k][dim][i] = datari[ix2im++];
      }
      ix0 = idx[k][0];

      if(nvr[k] == 1)
      {
         ix1re = (1 + k)*deg1;
         ix2im = (1 + k)*deg1 + totcffoffset;
            
         for(int i=0; i<=deg; i++)
         {
            outputre[k][ix0][i] = datari[ix1re++];
            outputim[k][ix0][i] = datari[ix2im++];
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
            outputre[k][ix0][i] = datari[ix1re++];
            outputim[k][ix0][i] = datari[ix2im++];
         }
         ix2 = nvr[k]-2;
         ix1re = fstart[k] + ix2*deg1;
         ix2im = fstart[k] + ix2*deg1 + totcffoffset;
         ix0 = idx[k][ix2+1];

         for(int i=0; i<=deg; i++)
         {
            outputre[k][ix0][i] = datari[ix1re++];
            outputim[k][ix0][i] = datari[ix2im++];
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
                  outputre[k][ix0][i] = datari[ix1re++];
                  outputim[k][ix0][i] = datari[ix2im++];
               }
            }
         }
      }
   }
}

// The code for GPU_dbl_mon_evaldiff is an adaptation of the
// function GPU_dbl_poly_evaldiff of dbl_polynomials_kernels.cu.

void GPU_dbl_mon_evaldiff
 ( int szt, int dim, int nbr, int deg, int *nvr, int **idx,
   double **cff, double **input, double ***output,
   ConvolutionJobs cnvjobs,
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

   double *data_h = new double[totalcff];        // data on host
   int ix = 0;
   for(int i=0; i<deg1; i++) data_h[ix++] = 0.0; // cst[i]; no constant
   for(int i=0; i<nbr; i++)
      for(int j=0; j<deg1; j++) data_h[ix++] = cff[i][j];
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++) data_h[ix++] = input[i][j];

   double *data_d;                               // device data
   const size_t szdata = totalcff*sizeof(double);
   cudaMalloc((void**)&data_d,szdata);
   cudaMemcpy(data_d,data_h,szdata,cudaMemcpyHostToDevice);

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
         // dbl_padded_convjobs<<<jobnbr,szt>>>
         dbl_padded_convjobs<<<jobnbr,deg1>>>
            (data_d,in1ix_d,in2ix_d,outix_d,deg1);
         cudaEventRecord(stop);
         cudaEventSynchronize(stop);
         cudaEventElapsedTime(&milliseconds,start,stop);
         *cnvlapms += milliseconds;

         cudaFree(in1ix_d); cudaFree(in2ix_d); cudaFree(outix_d);
      }
      free(in1ix_h); free(in2ix_h); free(outix_h);
   }
   gettimeofday(&endtime,0);
   cudaMemcpy(data_h,data_d,szdata,cudaMemcpyDeviceToHost);
   *elapsedms = *cnvlapms;
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   *walltimesec = seconds + microseconds*1.0e-6;

   dbl_evaldiffdata_to_output
      (data_h,output,dim,nbr,deg,nvr,idx,fstart,bstart,cstart,vrblvl);

   if(vrblvl > 0)
   {
      write_GPU_timings(*cnvlapms,0.0,*elapsedms,*walltimesec);
      write_dbl_cnvflops(dim,deg,0,cnvjobs,*elapsedms,*walltimesec);
   }
   cudaFree(data_d);
   free(data_h);

   free(fstart); free(bstart); free(cstart);
   free(fsums); free(bsums); free(csums);
}

void GPU_cmplx_mon_evaldiff
 ( int szt, int dim, int nbr, int deg, int *nvr, int **idx,
   double **cffre, double **cffim, double **inputre, double **inputim,
   double ***outputre, double ***outputim, ConvolutionJobs cnvjobs,
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

   double *datare_h = new double[totalcff];   // real data on host
   double *dataim_h = new double[totalcff];   // imaginary data on host
   int ix = 0;
   for(int i=0; i<deg1; i++)
   {
      datare_h[ix]   = 0.0; // cst[i]; no constant
      dataim_h[ix++] = 0.0;
   }
   for(int i=0; i<nbr; i++)
      for(int j=0; j<deg1; j++)
      {
         datare_h[ix]   = cffre[i][j];
         dataim_h[ix++] = cffim[i][j];
      }
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++)
      {
         datare_h[ix]   = inputre[i][j];
         dataim_h[ix++] = inputim[i][j];
      }

   double *datare_d;                        // device real data
   double *dataim_d;                        // device complex data
   const size_t szdata = totalcff*sizeof(double);
   cudaMalloc((void**)&datare_d,szdata);
   cudaMalloc((void**)&dataim_d,szdata);
   cudaMemcpy(datare_d,datare_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(dataim_d,dataim_h,szdata,cudaMemcpyHostToDevice);

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
         cmplx_padded_convjobs<<<jobnbr,deg1>>>
            (datare_d,dataim_d,in1ix_d,in2ix_d,outix_d,deg1);
         cudaEventRecord(stop);
         cudaEventSynchronize(stop);
         cudaEventElapsedTime(&milliseconds,start,stop);
         *cnvlapms += milliseconds;

         cudaFree(in1ix_d); cudaFree(in2ix_d); cudaFree(outix_d);
      }
      free(in1ix_h); free(in2ix_h); free(outix_h);
   }
   gettimeofday(&endtime,0);
   cudaMemcpy(datare_h,datare_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(dataim_h,dataim_d,szdata,cudaMemcpyDeviceToHost);
   *elapsedms = *cnvlapms;
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   *walltimesec = seconds + microseconds*1.0e-6;

   cmplx_evaldiffdata_to_output
      (datare_h,dataim_h,outputre,outputim,
       dim,nbr,deg,nvr,idx,fstart,bstart,cstart,vrblvl);

   if(vrblvl > 0)
   {
      write_GPU_timings(*cnvlapms,0.0,*elapsedms,*walltimesec);
      write_dbl_cnvflops(dim,deg,1,cnvjobs,*elapsedms,*walltimesec);
   }
   cudaFree(datare_d); cudaFree(dataim_d);
   free(datare_h); free(dataim_h);

   free(fstart); free(bstart); free(cstart);
   free(fsums); free(bsums); free(csums);
}

void GPU_cmplxvectorized_mon_evaldiff
 ( int szt, int dim, int nbr, int deg, int *nvr, int **idx,
   double **cffre, double **cffim, double **inputre, double **inputim,
   double ***outputre, double ***outputim,
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
   double *datari_h = new double[cmplxtotcff];   // data on host

   int ix1 = 0;
   int ix2 = totalcff + offsetri;

   for(int i=0; i<deg1; i++)
   {
      datari_h[ix1++] = 0.0; // cst[i]; no constant
      datari_h[ix2++] = 0.0;
   }
   for(int i=0; i<nbr; i++)
      for(int j=0; j<deg1; j++)
      {
         datari_h[ix1++] = cffre[i][j];
         datari_h[ix2++] = cffim[i][j];
      }
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++)
      {
         datari_h[ix1++] = inputre[i][j];
         datari_h[ix2++] = inputim[i][j];
      }

   for(int i=0; i<2*offsetri; i++)
   {
      datari_h[ix1++] = 0.0;
      datari_h[ix2++] = 0.0;
   }

   double *datari_d;                        // device data
   const size_t szdata = cmplxtotcff*sizeof(double);
   cudaMalloc((void**)&datari_d,szdata);
   cudaMemcpy(datari_d,datari_h,szdata,cudaMemcpyHostToDevice);

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
         dbl_padded_convjobs<<<jobnbr,deg1>>>
            (datari_d,in1ix_d,in2ix_d,outix_d,deg1);
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

      GPU_cmplxvectorized_flipsigns
         (deg,nbrflips,rebidx,datari_d,&fliplapms,vrblvl);
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

         if(vrblvl > 1)
            cout << "launching " << jobnbr << " blocks of " << deg1
                 << " threads for increments ..." << endl;

         cudaEventRecord(start);
         dbl_increment_jobs<<<jobnbr,deg1>>>
            (datari_d,in1ix_d,in2ix_d,outix_d,deg1);
         cudaEventRecord(stop);
         cudaEventSynchronize(stop);
         cudaEventElapsedTime(&milliseconds,start,stop);
         *cnvlapms += milliseconds;
         cudaFree(in1ix_d); cudaFree(in2ix_d); cudaFree(outix_d);
      }
      free(in1ix_h); free(in2ix_h); free(outix_h); free(rebidx);
   }
   gettimeofday(&endtime,0);
   cudaMemcpy(datari_h,datari_d,szdata,cudaMemcpyDeviceToHost);
   *elapsedms = *cnvlapms;
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   *walltimesec = seconds + microseconds*1.0e-6;

   cmplxvectorized_evaldiffdata_to_output
      (datari_h,outputre,outputim,dim,nbr,deg,nvr,idx,
       fstart,bstart,cstart,totalcff,offsetri,vrblvl);

   if(vrblvl > 0)
   {
      write_GPU_timings(*cnvlapms,0.0,*elapsedms,*walltimesec);
      write_vectorized_cnvincflops
          (dim,deg,cnvjobs,incjobs,*elapsedms,*walltimesec);
   }
   cudaFree(datari_d);
   free(datari_h);

   free(fstart); free(bstart); free(cstart);
   free(fsums); free(bsums); free(csums);
}

void GPU_dbl_evaluate_monomials
 ( int dim, int deg, int szt, int nbt,
   int *nvr, int **idx, int **exp, int *nbrfac,
   int **expfac, double **cff, double *acc, double **input,
   double ***output, double *totcnvlapsedms, int vrblvl )
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
                  CPU_dbl_product(deg,input[idxvar],cff[i],acc);
                  for(int L=0; L<=deg; L++) cff[i][L] = acc[L];
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
         for(int j=0; j<=deg; j++) cout << cff[i][j] << endl;
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
         for(int j=0; j<=deg; j++) cout << input[i][j] << endl;
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
   GPU_dbl_mon_evaldiff
      (szt,dim,dim,deg,nvr,idx,cff,input,output,jobs,
       &cnvlapms,&elapsedms,&walltimesec,vrblvl);

   *totcnvlapsedms += elapsedms;

   if(vrblvl > 1)
   {
      for(int i=0; i<dim; i++)
      {
         cout << "output series for monomial " << i << " :" << endl;
         for(int j=0; j<=deg; j++) cout << output[i][dim][j] << endl;
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
                  output[i][idxvar][k] = factor*output[i][idxvar][k];
            }
         }
      }
   }
}

void GPU_cmplx_evaluate_monomials
 ( int dim, int deg, int szt, int nbt,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffre, double **cffim, double *accre, double *accim,
   double **inputre, double **inputim, double ***outputre, double ***outputim,
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
                  CPU_cmplx_product
                     (deg,inputre[idxvar],inputim[idxvar],
                      cffre[i],cffim[i],accre,accim);

                  for(int L=0; L<=deg; L++)
                  {
                     cffre[i][L] = accre[L];
                     cffim[i][L] = accim[L];
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
            cout << cffre[i][j] << "  " << cffim[i][j] << endl;
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
            cout << inputre[i][j] << "  " << inputim[i][j] << endl;
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
   GPU_cmplx_mon_evaldiff
      (szt,dim,dim,deg,nvr,idx,cffre,cffim,inputre,inputim,
       outputre,outputim,jobs,&cnvlapms,&elapsedms,&walltimesec,vrblvl);

   *totcnvlapsedms += elapsedms;

   if(vrblvl > 1)
   {
      for(int i=0; i<dim; i++)
      {
         cout << "output series for monomial " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << outputre[i][dim][j] << "  "
                 << outputim[i][dim][j] << endl;
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

               // multiply derivative w.r.t. idxvar with factor
               for(int k=0; k<=deg; k++)
               {
                  outputre[i][idxvar][k] = fac*outputre[i][idxvar][k];
                  outputim[i][idxvar][k] = fac*outputim[i][idxvar][k];
               }
            }
         }
      }
   }
}

void GPU_cmplxvectorized_evaluate_monomials
 ( int dim, int deg, int szt, int nbt,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffre, double **cffim, double *accre, double *accim,
   double **inputre, double **inputim, double ***outputre, double ***outputim,
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
                  CPU_cmplx_product
                     (deg,inputre[idxvar],inputim[idxvar],
                      cffre[i],cffim[i],accre,accim);

                  for(int L=0; L<=deg; L++)
                  {
                     cffre[i][L] = accre[L];
                     cffim[i][L] = accim[L];
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
            cout << cffre[i][j] << "  " << cffim[i][j] << endl;
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
            cout << inputre[i][j] << "  " << inputim[i][j] << endl;
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
         cout << "convolution jobs at layer " << k << " :" << endl;
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
   GPU_cmplxvectorized_mon_evaldiff
      (szt,dim,dim,deg,nvr,idx,cffre,cffim,inputre,inputim,outputre,outputim,
       cnvjobs,incjobs,&cnvlapms,&elapsedms,&walltimesec,vrblvl);

   *totcnvlapsedms += elapsedms;

   if(vrblvl > 1)
   {
      for(int i=0; i<dim; i++)
      {
         cout << "output series for monomial " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << outputre[i][dim][j] << "  "
                 << outputim[i][dim][j] << endl;
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

               // multiply derivative w.r.t. idxvar with factor
               for(int k=0; k<=deg; k++)
               {
                  outputre[i][idxvar][k] = fac*outputre[i][idxvar][k];
                  outputim[i][idxvar][k] = fac*outputim[i][idxvar][k];
               }
            }
         }
      }
   }
}

void GPU_dbl_evaluate_columns
 ( int dim, int deg, int nbrcol, int szt, int nbt, int **nvr, int ***idx,
   double ***cff, double **input, double ***output,
   double **funval, double ***jacval, double *totcnvlapsedms, int vrblvl )
{
   const int degp1 = deg+1;

   for(int i=0; i<dim; i++)
      for(int j=0; j<degp1; j++) funval[i][j] = 0.0;

   funval[dim-1][0] = -1.0; // constant of last eq in cyclic dim-roots

   for(int k=0; k<degp1; k++)  // the Jacobian is linearized
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++) jacval[k][i][j] = 0.0;

   if(vrblvl > 1)
   {
      for(int i=0; i<nbrcol; i++)
      {
         cout << "coefficients for column " << i << " :" << endl;
         for(int j=0; j<dim; j++)
         {
            cout << "coefficients for monomial " << j << " :" << endl;
            for(int k=0; k<=deg; k++) cout << cff[i][j][k] << endl;
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
         for(int j=0; j<=deg; j++) cout << input[i][j] << endl;
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
      GPU_dbl_mon_evaldiff
         (szt,dim,dim,deg,nvr[i],idx[i],cff[i],input,output,jobs,
          &cnvlapms,&elapsedms,&walltimesec,vrblvl);

      *totcnvlapsedms += elapsedms;

      for(int j=0; j<dim; j++)
         if(nvr[i][j] > 0)       // update values
         {
            for(int L=0; L<degp1; L++) funval[j][L] += output[j][dim][L];

            int *indexes = idx[i][j];      // indices of the variables
            for(int k=0; k<nvr[i][j]; k++) // derivative w.r.t. idx[i][j][k]
            {                              // has j-th coefficient
               int idxval = indexes[k];
               for(int L=0; L<degp1; L++) 
                  jacval[L][j][idxval] += output[j][idxval][L];
            }
         }
   }
   if(vrblvl > 1)
   {
      for(int i=0; i<dim; i++)
      {
         cout << "output series for monomial " << i << " :" << endl;
         for(int j=0; j<=deg; j++) cout << funval[i][j] << endl;
      }
   }
}

void GPU_cmplx_evaluate_columns
 ( int dim, int deg, int nbrcol, int szt, int nbt, int **nvr, int ***idx, 
   double ***cffre, double ***cffim, double **inputre, double **inputim,
   double ***outputre, double ***outputim, 
   double **funvalre, double **funvalim,
   double ***jacvalre, double ***jacvalim, double *totcnvlapsedms,
   int vrblvl )
{
   const int degp1 = deg+1;

   for(int i=0; i<dim; i++)
      for(int j=0; j<degp1; j++)
      {
         funvalre[i][j] = 0.0;
         funvalim[i][j] = 0.0;
      }
   funvalre[dim-1][0] = -1.0; // constant of last eq in cyclic dim-roots

   for(int k=0; k<degp1; k++)  // the Jacobian is linearized
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
         {
            jacvalre[k][i][j] = 0.0;
            jacvalim[k][i][j] = 0.0;
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
               cout << cffre[i][j][k] << "  " << cffim[i][j][k] << endl;
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
            cout << inputre[i][j] << "  " << inputim[i][j] << endl;
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
      GPU_cmplx_mon_evaldiff
         (szt,dim,dim,deg,nvr[i],idx[i],cffre[i],cffim[i],inputre,inputim,
          outputre,outputim,jobs,&cnvlapms,&elapsedms,&walltimesec,vrblvl);

      *totcnvlapsedms += elapsedms;

      for(int j=0; j<dim; j++)
         if(nvr[i][j] > 0)       // update values
         {
            for(int L=0; L<degp1; L++)
            {
               funvalre[j][L] += outputre[j][dim][L];
               funvalim[j][L] += outputim[j][dim][L];
            }

            int *indexes = idx[i][j];      // indices of the variables
            for(int k=0; k<nvr[i][j]; k++) // derivative w.r.t. idx[i][j][k]
            {                              // has j-th coefficient
               int idxval = indexes[k];
               for(int L=0; L<degp1; L++) 
               {
                  jacvalre[L][j][idxval] += outputre[j][idxval][L];
                  jacvalim[L][j][idxval] += outputim[j][idxval][L];
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
            cout << funvalre[i][j] << "  "
                 << funvalim[i][j] << endl;
      }
   }
}
