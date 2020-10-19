/* This is the main program to run manual tests,
   with four input parameters at the command line.

For example, typing at the command prompt
  ./run_cmplx3_norm_d 64 256 1 2
computes the norm once with 64 threads in a block,
of a vector of dimension 128, on both GPU and CPU,
and applies this norm to normalize a random vector.
As all complex numbers are on the complex unit circle,
the 2-norm of a vector of dimension 64 equals 8.
If the dimension is 256, then the 2-norm is 16.

Including the vector_types.h is needed for the include of
the headers for the cmplx3_norm_kernels. */

#include <ctime>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "cmplx3_norm_kernels.h"
#include "cmplx3_norm_host.h"
#include "random3_vectors.h"
#include "parse_run_arguments.h"
#include "triple_double_functions.h"

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

   double* vrehi_host = new double[dim];   // high real parts on host
   double* vremi_host = new double[dim];   // middle real parts on host
   double* vrelo_host = new double[dim];   // low real parts on host
   double* vimhi_host = new double[dim];   // high imaginary parts on host
   double* vimmi_host = new double[dim];   // middle imaginary parts on host
   double* vimlo_host = new double[dim];   // low imaginary parts on host
   double* wrehi_host = new double[dim];   // a copy for normalization
   double* wremi_host = new double[dim];
   double* wrelo_host = new double[dim];
   double* wimhi_host = new double[dim];
   double* wimmi_host = new double[dim];
   double* wimlo_host = new double[dim];
   double* vrehi_device = new double[dim]; // high real parts on the device
   double* vremi_device = new double[dim]; // middle real parts on the device
   double* vrelo_device = new double[dim]; // low real parts on the device
   double* vimhi_device = new double[dim]; // high imaginary parts on device
   double* vimmi_device = new double[dim]; // middle imaginary parts on device
   double* vimlo_device = new double[dim]; // low imaginary parts on device

   random_complex3_vectors
      (dim,vrehi_host,vremi_host,vrelo_host,
           vimhi_host,vimmi_host,vimlo_host,
       vrehi_device,vremi_device,vrelo_device,
       vimhi_device,vimmi_device,vimlo_device);

   double vhinorm_device,vminorm_device,vlonorm_device;
   double whinorm_device,wminorm_device,wlonorm_device;
   double vhinorm_host,vminorm_host,vlonorm_host;
   double whinorm_host,wminorm_host,wlonorm_host;

   if(mode==0 || mode==2) // GPU computation of the norm
   {
      GPU_norm(vrehi_device,vremi_device,vrelo_device,
               vimhi_device,vimmi_device,vimlo_device,dim,1,BS,
               &vhinorm_device,&vminorm_device,&vlonorm_device,1);
      GPU_norm(vrehi_device,vremi_device,vrelo_device,
               vimhi_device,vimmi_device,vimlo_device,dim,freq,BS,
               &whinorm_device,&wminorm_device,&wlonorm_device,1);
   }

   if(mode==1 || mode==2) // CPU computation of the norm
      for(int i=0; i<=freq; i++)
      {
         CPU_norm(vrehi_host,vremi_host,vrelo_host,
                  vimhi_host,vimmi_host,vimlo_host,dim,
                  &vhinorm_host,&vminorm_host,&vlonorm_host);
         make_copy(dim,vrehi_host,vremi_host,vrelo_host,
                       vimhi_host,vimmi_host,vimlo_host,
                   wrehi_host,wremi_host,wrelo_host,
                   wimhi_host,wimmi_host,wimlo_host);
         CPU_normalize(wrehi_host,wremi_host,wrelo_host,
                       wimhi_host,wimmi_host,wimlo_host,dim,
                       vhinorm_host,vminorm_host,vlonorm_host);
         CPU_norm(wrehi_host,wremi_host,wrelo_host,
                  wimhi_host,wimmi_host,wimlo_host,dim,
                  &whinorm_host,&wminorm_host,&wlonorm_host);
      }

   if(mode == 2) // GPU vs CPU correctness verification
   {
      cout << scientific << setprecision(16);

      cout << "CPU norm : " << endl;
      tdf_write_doubles(vhinorm_host,vminorm_host,vlonorm_host);
      cout << endl;
      cout << "GPU norm : " << endl;
      tdf_write_doubles(vhinorm_device,vminorm_device,vlonorm_device);
      cout << endl;
      cout << "CPU norm after normalization : " << endl;
      tdf_write_doubles(whinorm_host,wminorm_host,wlonorm_host);
      cout << endl;
      cout << "GPU norm after normalization : " << endl;
      tdf_write_doubles(whinorm_device,wminorm_device,wlonorm_device);
      cout << endl;
   }

   return 0;
}
