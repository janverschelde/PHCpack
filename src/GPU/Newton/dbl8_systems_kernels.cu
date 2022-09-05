// The file dbl8_systems_kernels.cu defines the functions with prototypes in
// the file dbl8_systems_kernels.h.

#include <iostream>
#include <iomanip>
#ifdef winwalltime
#include "gettimeofday4win.h"
#else
#include <sys/time.h>
#endif
#include "write_gpu_timings.h"
#include "job_coordinates.h"
#include "dbl8_convolutions_host.h"
#include "dbl8_monomials_host.h"
#include "dbl8_polynomials_kernels.h"

using namespace std;

// The code in dbl8_evaldiffdata_to_output is an adaptation of the
// function dbl8_convoluted_data_to_output in dbl2_polynomials_kernels.cu.

void dbl8_evaldiffdata_to_output
 ( double *datahihihi, double *datalohihi,
   double *datahilohi, double *datalolohi,
   double *datahihilo, double *datalohilo,
   double *datahilolo, double *datalololo,
   double ***outputhihihi, double ***outputlohihi,
   double ***outputhilohi, double ***outputlolohi,
   double ***outputhihilo, double ***outputlohilo,
   double ***outputhilolo, double ***outputlololo,
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
         outputhihihi[k][dim][i] = datahihihi[ix1];
         outputlohihi[k][dim][i] = datalohihi[ix1];
         outputhilohi[k][dim][i] = datahilohi[ix1];
         outputlolohi[k][dim][i] = datalolohi[ix1];
         outputhihilo[k][dim][i] = datahihilo[ix1];
         outputlohilo[k][dim][i] = datalohilo[ix1];
         outputhilolo[k][dim][i] = datahilolo[ix1];
         outputlololo[k][dim][i] = datalololo[ix1++];
      }
      ix0 = idx[k][0];
      if(nvr[k] == 1)
      {
         ix1 = (1 + k)*deg1;
            
         for(int i=0; i<=deg; i++)
         {
            outputhihihi[k][ix0][i] = datahihihi[ix1];
            outputlohihi[k][ix0][i] = datalohihi[ix1];
            outputhilohi[k][ix0][i] = datahilohi[ix1];
            outputlolohi[k][ix0][i] = datalolohi[ix1];
            outputhihilo[k][ix0][i] = datahihilo[ix1];
            outputlohilo[k][ix0][i] = datalohilo[ix1];
            outputhilolo[k][ix0][i] = datahilolo[ix1];
            outputlololo[k][ix0][i] = datalololo[ix1++];
         }
      }
      else
      {                               // first and last derivative
         ix2 = nvr[k]-3;
         if(ix2 < 0) ix2 = 0;
         ix1 = bstart[k] + ix2*deg1;

         for(int i=0; i<=deg; i++)
         {
            outputhihihi[k][ix0][i] = datahihihi[ix1];
            outputlohihi[k][ix0][i] = datalohihi[ix1];
            outputhilohi[k][ix0][i] = datahilohi[ix1];
            outputlolohi[k][ix0][i] = datalolohi[ix1];
            outputhihilo[k][ix0][i] = datahihilo[ix1];
            outputlohilo[k][ix0][i] = datalohilo[ix1];
            outputhilolo[k][ix0][i] = datahilolo[ix1];
            outputlololo[k][ix0][i] = datalololo[ix1++];
         }
         ix2 = nvr[k]-2;
         ix1 = fstart[k] + ix2*deg1;
         ix0 = idx[k][ix2+1];

         for(int i=0; i<=deg; i++)
         {
            outputhihihi[k][ix0][i] = datahihihi[ix1];
            outputlohihi[k][ix0][i] = datalohihi[ix1];
            outputhilohi[k][ix0][i] = datahilohi[ix1];
            outputlolohi[k][ix0][i] = datalolohi[ix1];
            outputhihilo[k][ix0][i] = datahihilo[ix1];
            outputlohilo[k][ix0][i] = datalohilo[ix1];
            outputhilolo[k][ix0][i] = datahilolo[ix1];
            outputlololo[k][ix0][i] = datalololo[ix1++];
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
                  outputhihihi[k][ix0][i] = datahihihi[ix1];
                  outputlohihi[k][ix0][i] = datalohihi[ix1];
                  outputhilohi[k][ix0][i] = datahilohi[ix1];
                  outputlolohi[k][ix0][i] = datalolohi[ix1];
                  outputhihilo[k][ix0][i] = datahihilo[ix1];
                  outputlohilo[k][ix0][i] = datalohilo[ix1];
                  outputhilolo[k][ix0][i] = datahilolo[ix1];
                  outputlololo[k][ix0][i] = datalololo[ix1++];
               }
            }
         }
      }
   }
}

// The code for GPU_dbl8_mon_evaldiff is an adaptation of the
// function GPU_dbl8_poly_evaldiff of dbl_polynomials_kernels.cu.

void GPU_dbl8_mon_evaldiff
 ( int szt, int dim, int nbr, int deg, int *nvr, int **idx,
   double **cffhihihi, double **cfflohihi,
   double **cffhilohi, double **cfflolohi,
   double **cffhihilo, double **cfflohilo,
   double **cffhilolo, double **cfflololo,
   double **inputhihihi, double **inputlohihi,
   double **inputhilohi, double **inputlolohi,
   double **inputhihilo, double **inputlohilo,
   double **inputhilolo, double **inputlololo,
   double ***outputhihihi, double ***outputlohihi,
   double ***outputhilohi, double ***outputlolohi,
   double ***outputhihilo, double ***outputlohilo,
   double ***outputhilolo, double ***outputlololo, ConvolutionJobs cnvjobs,
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

   double *datahihihi_h = new double[totalcff];        // data on host
   double *datalohihi_h = new double[totalcff];
   double *datahilohi_h = new double[totalcff];
   double *datalolohi_h = new double[totalcff];
   double *datahihilo_h = new double[totalcff];
   double *datalohilo_h = new double[totalcff];
   double *datahilolo_h = new double[totalcff];
   double *datalololo_h = new double[totalcff];

   int ix = 0;
   for(int i=0; i<deg1; i++)
   {
      datahihihi_h[ix]   = 0.0; // cst[i]; no constant
      datalohihi_h[ix]   = 0.0; 
      datahilohi_h[ix]   = 0.0; 
      datalolohi_h[ix]   = 0.0;
      datahihilo_h[ix]   = 0.0;
      datalohilo_h[ix]   = 0.0; 
      datahilolo_h[ix]   = 0.0; 
      datalololo_h[ix++] = 0.0;
   }
   for(int i=0; i<nbr; i++)
      for(int j=0; j<deg1; j++)
      {
         datahihihi_h[ix]   = cffhihihi[i][j];
         datalohihi_h[ix]   = cfflohihi[i][j];
         datahilohi_h[ix]   = cffhilohi[i][j];
         datalolohi_h[ix]   = cfflolohi[i][j];
         datahihilo_h[ix]   = cffhihilo[i][j];
         datalohilo_h[ix]   = cfflohilo[i][j];
         datahilolo_h[ix]   = cffhilolo[i][j];
         datalololo_h[ix++] = cfflololo[i][j];
      }

   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++)
      {
         datahihihi_h[ix]   = inputhihihi[i][j];
         datalohihi_h[ix]   = inputlohihi[i][j];
         datahilohi_h[ix]   = inputhilohi[i][j];
         datalolohi_h[ix]   = inputlolohi[i][j];
         datahihilo_h[ix]   = inputhihilo[i][j];
         datalohilo_h[ix]   = inputlohilo[i][j];
         datahilolo_h[ix]   = inputhilolo[i][j];
         datalololo_h[ix++] = inputlololo[i][j];
      }

   double *datahihihi_d;                               // device data
   double *datalohihi_d;
   double *datahilohi_d;
   double *datalolohi_d;
   double *datahihilo_d;
   double *datalohilo_d;
   double *datahilolo_d;
   double *datalololo_d;
   const size_t szdata = totalcff*sizeof(double);
   cudaMalloc((void**)&datahihihi_d,szdata);
   cudaMalloc((void**)&datalohihi_d,szdata);
   cudaMalloc((void**)&datahilohi_d,szdata);
   cudaMalloc((void**)&datalolohi_d,szdata);
   cudaMalloc((void**)&datahihilo_d,szdata);
   cudaMalloc((void**)&datalohilo_d,szdata);
   cudaMalloc((void**)&datahilolo_d,szdata);
   cudaMalloc((void**)&datalololo_d,szdata);
   cudaMemcpy(datahihihi_d,datahihihi_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datalohihi_d,datalohihi_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datahilohi_d,datahilohi_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datalolohi_d,datalolohi_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datahihilo_d,datahihilo_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datalohilo_d,datalohilo_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datahilolo_d,datahilolo_h,szdata,cudaMemcpyHostToDevice);
   cudaMemcpy(datalololo_d,datalololo_h,szdata,cudaMemcpyHostToDevice);

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
         dbl8_padded_convjobs<<<jobnbr,deg1>>>
            (datahihihi_d,datalohihi_d,datahilohi_d,datalolohi_d,
             datahihilo_d,datalohilo_d,datahilolo_d,datalololo_d,
             in1ix_d,in2ix_d,outix_d,deg1);
         cudaEventRecord(stop);
         cudaEventSynchronize(stop);
         cudaEventElapsedTime(&milliseconds,start,stop);
         *cnvlapms += milliseconds;
      }
      free(in1ix_h); free(in2ix_h); free(outix_h);
   }
   gettimeofday(&endtime,0);
   cudaMemcpy(datahihihi_h,datahihihi_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datalohihi_h,datalohihi_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datahilohi_h,datahilohi_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datalolohi_h,datalolohi_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datahihilo_h,datahihilo_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datalohilo_h,datalohilo_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datahilolo_h,datahilolo_d,szdata,cudaMemcpyDeviceToHost);
   cudaMemcpy(datalololo_h,datalololo_d,szdata,cudaMemcpyDeviceToHost);
   *elapsedms = *cnvlapms;
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   *walltimesec = seconds + microseconds*1.0e-6;

   dbl8_evaldiffdata_to_output
      (datahihihi_h,datalohihi_h,datahilohi_h,datalolohi_h,
       datahihilo_h,datalohilo_h,datahilolo_h,datalololo_h,
       outputhihihi,outputlohihi,outputhilohi,outputlolohi,
       outputhihilo,outputlohilo,outputhilolo,outputlololo,
       dim,nbr,deg,nvr,idx,fstart,bstart,cstart,verbose);

   if(verbose)
      write_GPU_timings(*cnvlapms,0.0,*elapsedms,*walltimesec);
}

void GPU_dbl8_evaluate_monomials
 ( int dim, int deg, int szt, int nbt,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffhihihi, double **cfflohihi,
   double **cffhilohi, double **cfflolohi,
   double **cffhihilo, double **cfflohilo,
   double **cffhilolo, double **cfflololo,
   double *acchihihi, double *acclohihi, double *acchilohi, double *acclolohi,
   double *acchihilo, double *acclohilo, double *acchilolo, double *acclololo,
   double **inputhihihi, double **inputlohihi,
   double **inputhilohi, double **inputlolohi,
   double **inputhihilo, double **inputlohilo,
   double **inputhilolo, double **inputlololo,
   double ***outputhihihi, double ***outputlohihi,
   double ***outputhilohi, double ***outputlolohi,
   double ***outputhihilo, double ***outputlohilo,
   double ***outputhilolo, double ***outputlololo, int vrblvl )
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
                  CPU_dbl8_product
                     (deg,inputhihihi[idxvar],inputlohihi[idxvar],
                          inputhilohi[idxvar],inputlolohi[idxvar],
                          inputhihilo[idxvar],inputlohilo[idxvar],
                          inputhilolo[idxvar],inputlololo[idxvar],
                      cffhihihi[i],cfflohihi[i],cffhilohi[i],cfflolohi[i],
                      cffhihilo[i],cfflohilo[i],cffhilolo[i],cfflololo[i],
                      acchihihi,acclohihi,acchilohi,acclolohi,
                      acchihilo,acclohilo,acchilolo,acclololo);

                  for(int L=0; L<=deg; L++)
                  {
                     cffhihihi[i][L] = acchihihi[L];
                     cfflohihi[i][L] = acclohihi[L];
                     cffhilohi[i][L] = acchilohi[L];
                     cfflolohi[i][L] = acclolohi[L];
                     cffhihilo[i][L] = acchihilo[L];
                     cfflohilo[i][L] = acclohilo[L];
                     cffhilolo[i][L] = acchilolo[L];
                     cfflololo[i][L] = acclololo[L];
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
         {
            cout << cffhihihi[i][j] << "  " << cfflohihi[i][j] << endl;
            cout << cffhilohi[i][j] << "  " << cfflolohi[i][j] << endl;
            cout << cffhihilo[i][j] << "  " << cfflohilo[i][j] << endl;
            cout << cffhilolo[i][j] << "  " << cfflololo[i][j] << endl;
         }
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
         {
            cout << inputhihihi[i][j] << "  " << inputlohihi[i][j] << endl;
            cout << inputhilohi[i][j] << "  " << inputlolohi[i][j] << endl;
            cout << inputhihilo[i][j] << "  " << inputlohilo[i][j] << endl;
            cout << inputhilolo[i][j] << "  " << inputlololo[i][j] << endl;
         }
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
   GPU_dbl8_mon_evaldiff
      (szt,dim,dim,deg,nvr,idx,
       cffhihihi,cfflohihi,cffhilohi,cfflolohi,
       cffhihilo,cfflohilo,cffhilolo,cfflololo,
       inputhihihi,inputlohihi,inputhilohi,inputlolohi,
       inputhihilo,inputlohilo,inputhilolo,inputlololo,
       outputhihihi,outputlohihi,outputhilohi,outputlolohi,
       outputhihilo,outputlohilo,outputhilolo,outputlololo,jobs,
       &cnvlapms,&elapsedms,&walltimesec,verbose);

   if(vrblvl > 0)
   {
      for(int i=0; i<dim; i++)
      {
         cout << "output series for monomial " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
         {
            cout << outputhihihi[i][dim][j] << "  "
                 << outputlohihi[i][dim][j] << endl;
            cout << outputhilohi[i][dim][j] << "  "
                 << outputlolohi[i][dim][j] << endl;
            cout << outputhihilo[i][dim][j] << "  "
                 << outputlohilo[i][dim][j] << endl;
            cout << outputhilolo[i][dim][j] << "  "
                 << outputlololo[i][dim][j] << endl;
         }
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
                  outputhihihi[i][idxvar][k] = fac*outputhihihi[i][idxvar][k];
                  outputlohihi[i][idxvar][k] = fac*outputlohihi[i][idxvar][k];
                  outputhilohi[i][idxvar][k] = fac*outputhilohi[i][idxvar][k];
                  outputlolohi[i][idxvar][k] = fac*outputlolohi[i][idxvar][k];
                  outputhihilo[i][idxvar][k] = fac*outputhihilo[i][idxvar][k];
                  outputlohilo[i][idxvar][k] = fac*outputlohilo[i][idxvar][k];
                  outputhilolo[i][idxvar][k] = fac*outputhilolo[i][idxvar][k];
                  outputlololo[i][idxvar][k] = fac*outputlololo[i][idxvar][k];
               }
            }
         }
      }
   }
}
