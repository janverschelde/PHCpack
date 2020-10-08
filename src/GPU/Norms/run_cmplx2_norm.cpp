/* This is the main program to run manual tests,
   with four input parameters at the command line.

For example, typing at the command prompt
  ./run_cmplx2_norm_d 64 64 1 2
computes the norm once with 64 threads in a block,
of a vector of dimension 64, on both GPU and CPU,
and applies this norm to normalize a random vector.
As all complex numbers are on the complex unit circle,
the 2-norm of a vector of dimension 64 equals 8.
If the dimension is 256, then the 2-norm is 16.

Including the vector_types.h is needed for the include of
the headers for the cmplx2_norm_kernels. */

#include <ctime>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "double_double.h"
#include "cmplx2_norm_kernels.h"
#include "cmplx2_norm_host.h"
#include "random2_vectors.h"
#include "parse_run_arguments.h"

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
   double* vrelo_host = new double[dim];   // low real parts on host
   double* vimhi_host = new double[dim];   // high imaginary parts on host
   double* vimlo_host = new double[dim];   // low imaginary parts on host
   double* wrehi_host = new double[dim];   // a copy for normalization
   double* wrelo_host = new double[dim];
   double* wimhi_host = new double[dim];
   double* wimlo_host = new double[dim];
   double* vrehi_device = new double[dim]; // high real parts on the device
   double* vrelo_device = new double[dim]; // low real parts on the device
   double* vimhi_device = new double[dim]; // high imaginary parts on device
   double* vimlo_device = new double[dim]; // low imaginary parts on device

   random_complex2_vectors
      (dim,vrehi_host,vrelo_host,vimhi_host,vimlo_host,
       vrehi_device,vrelo_device,vimhi_device,vimlo_device);

   double vhinorm_device,vlonorm_device,whinorm_device,wlonorm_device;
   double vhinorm_host,vlonorm_host,whinorm_host,wlonorm_host;

   if(mode==0 || mode==2) // GPU computation of the norm
   {
      GPU_norm(vrehi_device,vrelo_device,vimhi_device,vimlo_device,
               dim,1,BS,&vhinorm_device,&vlonorm_device,1);
      GPU_norm(vrehi_device,vrelo_device,vimhi_device,vimlo_device,
               dim,freq,BS,&whinorm_device,&wlonorm_device,1);
   }

   if(mode==1 || mode==2) // CPU computation of the norm
      for(int i=0; i<=freq; i++)
      {
         CPU_norm(vrehi_host,vrelo_host,vimhi_host,vimlo_host,dim,
                  &vhinorm_host,&vlonorm_host);
         make_copy(dim,vrehi_host,vrelo_host,vimhi_host,vimlo_host,
                   wrehi_host,wrelo_host,wimhi_host,wimlo_host);
         CPU_normalize(wrehi_host,wrelo_host,wimhi_host,wimlo_host,dim,
                       vhinorm_host,vlonorm_host);
         CPU_norm(wrehi_host,wrelo_host,wimhi_host,wimlo_host,dim,
                  &whinorm_host,&wlonorm_host);
      }

   if(mode == 2) // GPU vs CPU correctness verification
   {
      cout << scientific << setprecision(16);
      cout << "GPU norm hi : " << vhinorm_device << endl;
      cout << "GPU norm lo : " << vlonorm_device << endl;
      cout << "CPU norm hi : " << vhinorm_host << endl;
      cout << "CPU norm lo : " << vlonorm_host << endl;
      cout << "GPU norm after normalization : " << endl;
      cout << "         hi : " << whinorm_device << endl;
      cout << "         lo : " << wlonorm_device << endl;
      cout << "CPU norm after normalization : " << endl;
      cout << "         hi : " << whinorm_host << endl;
      cout << "         lo : " << wlonorm_host << endl;

      double vnrm_h[2],vnrm_d[2],wnrm_h[2],wnrm_d[2];

      vnrm_h[0] = vhinorm_host; vnrm_h[1] = vlonorm_host;
      wnrm_h[0] = whinorm_host; wnrm_h[1] = wlonorm_host;
      vnrm_d[0] = vhinorm_device; vnrm_d[1] = vlonorm_device;
      wnrm_d[0] = whinorm_device; wnrm_d[1] = wlonorm_device;

      cout << "   CPU norm : ";
      dd_write(vnrm_h,32); cout << endl;
      cout << "   GPU norm : ";
      dd_write(vnrm_d,32); cout << endl;
      cout << "   CPU norm after normalization : ";
      dd_write(wnrm_h,32); cout << endl;
      cout << "   GPU norm after normalization : ";
      dd_write(wnrm_d,32); cout << endl;
   }

   return 0;
}
