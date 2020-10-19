/* This is the main program to run manual tests,
   with four input parameters at the command line.

For example, typing at the command prompt
  ./run_cmplx4_norm_d 64 256 1 2
computes the norm once with 64 threads in a block,
of a vector of dimension 128, on both GPU and CPU,
and applies this norm to normalize a random vector.
As all complex numbers are on the complex unit circle,
the 2-norm of a vector of dimension 64 equals 8.
If the dimension is 256, then the 2-norm is 16.

Including the vector_types.h is needed for the include of
the headers for the cmplx4_norm_kernels. */

#include <ctime>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "cmplx4_norm_kernels.h"
#include "cmplx4_norm_host.h"
#include "random4_vectors.h"
#include "parse_run_arguments.h"
#include "quad_double_functions.h"

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

   double* vrehihi_host = new double[dim];   // highest real parts on host
   double* vrelohi_host = new double[dim];   // 2nd highest real parts on host
   double* vrehilo_host = new double[dim];   // 2nd lowest real parts on host
   double* vrelolo_host = new double[dim];   // lowest real parts on host
   double* vimhihi_host = new double[dim];   // highest imaginary parts 
   double* vimlohi_host = new double[dim];   // 2nd highest imaginary parts
   double* vimhilo_host = new double[dim];   // 2nd lowest imaginary parts
   double* vimlolo_host = new double[dim];   // lowest imaginary parts
   double* wrehihi_host = new double[dim];   // a copy for normalization
   double* wrelohi_host = new double[dim];
   double* wrehilo_host = new double[dim];
   double* wrelolo_host = new double[dim];
   double* wimhihi_host = new double[dim];
   double* wimlohi_host = new double[dim];
   double* wimhilo_host = new double[dim];
   double* wimlolo_host = new double[dim];
   double* vrehihi_device = new double[dim]; // highest real parts on device
   double* vrelohi_device = new double[dim]; // 2nd highest real parts
   double* vrehilo_device = new double[dim]; // 2nd lowest real parts
   double* vrelolo_device = new double[dim]; // lowest real parts
   double* vimhihi_device = new double[dim]; // highest imaginary parts
   double* vimlohi_device = new double[dim]; // 2nd highest imaginary parts
   double* vimhilo_device = new double[dim]; // 2nd lowest imaginary parts
   double* vimlolo_device = new double[dim]; // lowest imaginary parts

   random_complex4_vectors
      (dim,vrehihi_host,vrelohi_host,vrehilo_host,vrelolo_host,
           vimhihi_host,vimlohi_host,vimhilo_host,vimlolo_host,
       vrehihi_device,vrelohi_device,vrehilo_device,vrelolo_device,
       vimhihi_device,vimlohi_device,vimhilo_device,vimlolo_device);

   double vhihinorm_device,vlohinorm_device,vhilonorm_device,vlolonorm_device;
   double whihinorm_device,wlohinorm_device,whilonorm_device,wlolonorm_device;
   double vhihinorm_host,vlohinorm_host,vhilonorm_host,vlolonorm_host;
   double whihinorm_host,wlohinorm_host,whilonorm_host,wlolonorm_host;

   if(mode==0 || mode==2) // GPU computation of the norm
   {
      GPU_norm(vrehihi_device,vrelohi_device,vrehilo_device,vrelolo_device,
               vimhihi_device,vimlohi_device,vimhilo_device,vimlolo_device,
               dim,1,BS,
               &vhihinorm_device,&vlohinorm_device,
               &vhilonorm_device,&vlolonorm_device,1);
      GPU_norm(vrehihi_device,vrelohi_device,vrehilo_device,vrelolo_device,
               vimhihi_device,vimlohi_device,vimhilo_device,vimlolo_device,
               dim,freq,BS,
               &whihinorm_device,&wlohinorm_device,
               &whilonorm_device,&wlolonorm_device,1);
   }

   if(mode==1 || mode==2) // CPU computation of the norm
      for(int i=0; i<=freq; i++)
      {
         CPU_norm(vrehihi_host,vrelohi_host,vrehilo_host,vrelolo_host,
                  vimhihi_host,vimlohi_host,vimhilo_host,vimlolo_host,dim,
                  &vhihinorm_host,&vlohinorm_host,
                  &vhilonorm_host,&vlolonorm_host);
         make_copy(dim,vrehihi_host,vrelohi_host,vrehilo_host,vrelolo_host,
                       vimhihi_host,vimlohi_host,vimhilo_host,vimlolo_host,
                   wrehihi_host,wrelohi_host,wrehilo_host,wrelolo_host,
                   wimhihi_host,wimlohi_host,wimhilo_host,wimlolo_host);
         CPU_normalize
            (wrehihi_host,wrelohi_host,wrehilo_host,wrelolo_host,
             wimhihi_host,wimlohi_host,wimhilo_host,wimlolo_host,dim,
             vhihinorm_host,vlohinorm_host,vhilonorm_host,vlolonorm_host);
         CPU_norm(wrehihi_host,wrelohi_host,wrehilo_host,wrelolo_host,
                  wimhihi_host,wimlohi_host,wimhilo_host,wimlolo_host,dim,
                  &whihinorm_host,&wlohinorm_host,
                  &whilonorm_host,&wlolonorm_host);
      }

   if(mode == 2) // GPU vs CPU correctness verification
   {
      cout << scientific << setprecision(16);

      cout << "CPU norm : " << endl;
      qdf_write_doubles(vhihinorm_host,vlohinorm_host,
                        vhilonorm_host,vlolonorm_host);
      cout << "GPU norm : " << endl;
      qdf_write_doubles(vhihinorm_device,vlohinorm_device,
                        vhilonorm_device,vlolonorm_device);
      cout << "CPU norm after normalization : " << endl;
      qdf_write_doubles(whihinorm_host,wlohinorm_host,
                        whilonorm_host,wlolonorm_host);
      cout << "GPU norm after normalization : " << endl;
      qdf_write_doubles(whihinorm_device,wlohinorm_device,
                        whilonorm_device,wlolonorm_device);
   }

   return 0;
}
