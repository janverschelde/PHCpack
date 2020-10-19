/* This is the main program to run manual tests,
   with four input parameters at the command line.
   For example, typing at the command prompt
   ./run_dbl5_norm 32 64 1 2
   computes the norm once with 32 threads in a block,
   of a vector of dimension 64, on both GPU and CPU,
   and applies this norm to normalize a random vector.

Including the vector_types.h is needed for the include of
the headers for the dbl5_norm_kernels. */

#include <ctime>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "dbl5_norm_kernels.h"
#include "dbl5_norm_host.h"
#include "random5_vectors.h"
#include "parse_run_arguments.h"
#include "penta_double_functions.h"

using namespace std;

int main ( int argc, char *argv[] )
{
   // initialization of the execution parameters

   int BS,dim,freq,mode;
   if(parse_run_arguments(argc,argv,&BS,&dim,&freq,&mode) == 1) return 1;

   int timevalue;
   if(mode == 2)
      timevalue = time(NULL); // no fixed seed to verify correctness
   else
      timevalue = 1287178355; // fixed seed for timings
   srand(timevalue);

   double* vtb_host = new double[dim];   // highest parts on the host
   double* vix_host = new double[dim];   // 2nd highest parts on the host
   double* vmi_host = new double[dim];   // middle parts on the host
   double* vrg_host = new double[dim];   // 2nd lowest parts on the host
   double* vpk_host = new double[dim];   // lowest parts on the host
   double* vtb_device = new double[dim]; // highest parts on the device
   double* vix_device = new double[dim]; // 2nd highest parts on the device
   double* vmi_device = new double[dim]; // middle parts on the device
   double* vrg_device = new double[dim]; // 2nd lowest parts on the device
   double* vpk_device = new double[dim]; // lowest parts on the device
   double* wtb_host = new double[dim];   // highest parts copy
   double* wix_host = new double[dim];   // 2nd highest parts copy
   double* wmi_host = new double[dim];   // middle parts copy
   double* wrg_host = new double[dim];   // 2nd lowest parts copy
   double* wpk_host = new double[dim];   // lowest parts copy

   random_double5_vectors
      (dim,vtb_host,vix_host,vmi_host,vrg_host,vpk_host,
       vtb_device,vix_device,vmi_device,vrg_device,vpk_device);

   double vtbnorm_device,vixnorm_device,vrgnorm_device,vpknorm_device;
   double vminorm_device;
   double wtbnorm_device,wixnorm_device,wrgnorm_device,wpknorm_device;
   double wminorm_device;
   double vtbnorm_host,vixnorm_host,vminorm_host,vrgnorm_host,vpknorm_host;
   double wtbnorm_host,wixnorm_host,wminorm_host,wrgnorm_host,wpknorm_host;

   if(mode==0 || mode==2) // GPU computation of the norm
   {
      GPU_norm(vtb_device,vix_device,vmi_device,vrg_device,vpk_device,
               dim,1,BS,
               &vtbnorm_device,&vixnorm_device,&vminorm_device,
               &vrgnorm_device,&vpknorm_device,1);
      GPU_norm(vtb_device,vix_device,vmi_device,vrg_device,vpk_device,
               dim,freq,BS,
               &wtbnorm_device,&wixnorm_device,&wminorm_device,
               &wrgnorm_device,&wpknorm_device,1);
   }

   if(mode==1 || mode==2) // CPU computation of the norm
      for(int i=0; i<=freq; i++)
      {
         CPU_norm(vtb_host,vix_host,vmi_host,vrg_host,vpk_host,dim,
                  &vtbnorm_host,&vixnorm_host,&vminorm_host,
                  &vrgnorm_host,&vpknorm_host);
         make_copy(dim,vtb_host,vix_host,vmi_host,vrg_host,vpk_host,
                       wtb_host,wix_host,wmi_host,wrg_host,wpk_host);
         CPU_normalize(wtb_host,wix_host,wmi_host,wrg_host,wpk_host,dim,
                       vtbnorm_host,vixnorm_host,vminorm_host,
                       vrgnorm_host,vpknorm_host);
         CPU_norm(wtb_host,wix_host,wmi_host,wrg_host,wpk_host,dim,
                  &wtbnorm_host,&wixnorm_host,&wminorm_host,
                  &wrgnorm_host,&wpknorm_host);
      }

   if(mode == 2) // GPU vs CPU correctness verification
   {
      cout << scientific << setprecision(16);

      cout << "CPU norm : " << endl;
      pdf_write_doubles(vtbnorm_host,vixnorm_host,vminorm_host,
                        vrgnorm_host,vpknorm_host);
      cout << endl;
      cout << "GPU norm : " << endl;
      pdf_write_doubles(vtbnorm_device,vixnorm_device,vminorm_device,
                        vrgnorm_device,vpknorm_device);
      cout << endl;
      cout << "CPU norm after normalization : " << endl;
      pdf_write_doubles(wtbnorm_host,wixnorm_host,wminorm_host,
                        wrgnorm_host,wpknorm_host);
      cout << endl;
      cout << "GPU norm after normalization : " << endl;
      pdf_write_doubles(wtbnorm_device,wixnorm_device,wminorm_device,
                        wrgnorm_device,wpknorm_device);
      cout << endl;
   }

   return 0;
}
