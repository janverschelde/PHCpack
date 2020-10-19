/* This is the main program to run manual tests,
   with four input parameters at the command line.

For example, typing at the command prompt
  ./run_cmplx8_norm_d 64 256 1 2
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
#include "cmplx8_norm_kernels.h"
#include "cmplx8_norm_host.h"
#include "random8_vectors.h"
#include "parse_run_arguments.h"
#include "octo_double_functions.h"

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

   double* vrehihihi_host = new double[dim]; // highest real parts on host
   double* vrelohihi_host = new double[dim]; // 2nd highest real parts on host
   double* vrehilohi_host = new double[dim]; // 3rd highest real parts on host
   double* vrelolohi_host = new double[dim]; // 4th highest real parts on host
   double* vrehihilo_host = new double[dim]; // 4th lowest real parts on host
   double* vrelohilo_host = new double[dim]; // 3rd lowest real parts on host
   double* vrehilolo_host = new double[dim]; // 2nd lowest real parts on host
   double* vrelololo_host = new double[dim]; // lowest real parts on host
   double* vimhihihi_host = new double[dim]; // highest imaginary parts 
   double* vimlohihi_host = new double[dim]; // 2nd highest imaginary parts
   double* vimhilohi_host = new double[dim]; // 3rd highest imaginary parts
   double* vimlolohi_host = new double[dim]; // 4th highest imaginary parts
   double* vimhihilo_host = new double[dim]; // 4th lowest imaginary parts 
   double* vimlohilo_host = new double[dim]; // 3rd lowest imaginary parts
   double* vimhilolo_host = new double[dim]; // 2nd lowest imaginary parts
   double* vimlololo_host = new double[dim]; // lowest imaginary parts
   double* wrehihihi_host = new double[dim];   // a copy for normalization
   double* wrelohihi_host = new double[dim];
   double* wrehilohi_host = new double[dim];
   double* wrelolohi_host = new double[dim];
   double* wrehihilo_host = new double[dim];
   double* wrelohilo_host = new double[dim];
   double* wrehilolo_host = new double[dim];
   double* wrelololo_host = new double[dim];
   double* wimhihihi_host = new double[dim];
   double* wimlohihi_host = new double[dim];
   double* wimhilohi_host = new double[dim];
   double* wimlolohi_host = new double[dim];
   double* wimhihilo_host = new double[dim];
   double* wimlohilo_host = new double[dim];
   double* wimhilolo_host = new double[dim];
   double* wimlololo_host = new double[dim];
   double* vrehihihi_device = new double[dim]; // highest real parts on device
   double* vrelohihi_device = new double[dim]; // 2nd highest real parts
   double* vrehilohi_device = new double[dim]; // 3rd highest real parts
   double* vrelolohi_device = new double[dim]; // 4th highest real parts
   double* vrehihilo_device = new double[dim]; // 4th lowest real parts
   double* vrelohilo_device = new double[dim]; // 3rd lowest real parts
   double* vrehilolo_device = new double[dim]; // 2nd lowest real parts
   double* vrelololo_device = new double[dim]; // lowest real parts
   double* vimhihihi_device = new double[dim]; // highest imaginary parts
   double* vimlohihi_device = new double[dim]; // 2nd highest imaginary parts
   double* vimhilohi_device = new double[dim]; // 3rd highest imaginary parts
   double* vimlolohi_device = new double[dim]; // 4th highest imaginary parts
   double* vimhihilo_device = new double[dim]; // 4th lowest imaginary parts
   double* vimlohilo_device = new double[dim]; // 3rd lowest imaginary parts
   double* vimhilolo_device = new double[dim]; // 2nd lowest imaginary parts
   double* vimlololo_device = new double[dim]; // lowest imaginary parts

   random_complex8_vectors
      (dim,vrehihihi_host,vrelohihi_host,vrehilohi_host,vrelolohi_host,
           vrehihilo_host,vrelohilo_host,vrehilolo_host,vrelololo_host,
           vimhihihi_host,vimlohihi_host,vimhilohi_host,vimlolohi_host,
           vimhihilo_host,vimlohilo_host,vimhilolo_host,vimlololo_host,
       vrehihihi_device,vrelohihi_device,vrehilohi_device,vrelolohi_device,
       vrehihilo_device,vrelohilo_device,vrehilolo_device,vrelololo_device,
       vimhihihi_device,vimlohihi_device,vimhilohi_device,vimlolohi_device,
       vimhihilo_device,vimlohilo_device,vimhilolo_device,vimlololo_device);

   double vhihihinorm_device,vlohihinorm_device;
   double vhilohinorm_device,vlolohinorm_device;
   double vhihilonorm_device,vlohilonorm_device;
   double vhilolonorm_device,vlololonorm_device;
   double whihihinorm_device,wlohihinorm_device;
   double whilohinorm_device,wlolohinorm_device;
   double whihilonorm_device,wlohilonorm_device;
   double whilolonorm_device,wlololonorm_device;
   double vhihihinorm_host,vlohihinorm_host,vhilohinorm_host,vlolohinorm_host;
   double vhihilonorm_host,vlohilonorm_host,vhilolonorm_host,vlololonorm_host;
   double whihihinorm_host,wlohihinorm_host,whilohinorm_host,wlolohinorm_host;
   double whihilonorm_host,wlohilonorm_host,whilolonorm_host,wlololonorm_host;

   if(mode==0 || mode==2) // GPU computation of the norm
   {
      GPU_norm
         (vrehihihi_device,vrelohihi_device,vrehilohi_device,vrelolohi_device,
          vrehihilo_device,vrelohilo_device,vrehilolo_device,vrelololo_device,
          vimhihihi_device,vimlohihi_device,vimhilohi_device,vimlolohi_device,
          vimhihilo_device,vimlohilo_device,vimhilolo_device,vimlololo_device,
          dim,1,BS,
          &vhihihinorm_device,&vlohihinorm_device,
          &vhilohinorm_device,&vlolohinorm_device,
          &vhihilonorm_device,&vlohilonorm_device,
          &vhilolonorm_device,&vlololonorm_device,1);
      GPU_norm
         (vrehihihi_device,vrelohihi_device,vrehilohi_device,vrelolohi_device,
          vrehihilo_device,vrelohilo_device,vrehilolo_device,vrelololo_device,
          vimhihihi_device,vimlohihi_device,vimhilohi_device,vimlolohi_device,
          vimhihilo_device,vimlohilo_device,vimhilolo_device,vimlololo_device,
          dim,freq,BS,
          &whihihinorm_device,&wlohihinorm_device,
          &whilohinorm_device,&wlolohinorm_device,
          &whihilonorm_device,&wlohilonorm_device,
          &whilolonorm_device,&wlololonorm_device,1);
   }

   if(mode==1 || mode==2) // CPU computation of the norm
      for(int i=0; i<=freq; i++)
      {
         CPU_norm
            (vrehihihi_host,vrelohihi_host,vrehilohi_host,vrelolohi_host,
             vrehihilo_host,vrelohilo_host,vrehilolo_host,vrelololo_host,
             vimhihihi_host,vimlohihi_host,vimhilohi_host,vimlolohi_host,
             vimhihilo_host,vimlohilo_host,vimhilolo_host,vimlololo_host,dim,
             &vhihihinorm_host,&vlohihinorm_host,
             &vhilohinorm_host,&vlolohinorm_host,
             &vhihilonorm_host,&vlohilonorm_host,
             &vhilolonorm_host,&vlololonorm_host);
         make_copy
            (dim,vrehihihi_host,vrelohihi_host,vrehilohi_host,vrelolohi_host,
                 vrehihilo_host,vrelohilo_host,vrehilolo_host,vrelololo_host,
                 vimhihihi_host,vimlohihi_host,vimhilohi_host,vimlolohi_host,
                 vimhihilo_host,vimlohilo_host,vimhilolo_host,vimlololo_host,
             wrehihihi_host,wrelohihi_host,wrehilohi_host,wrelolohi_host,
             wrehihilo_host,wrelohilo_host,wrehilolo_host,wrelololo_host,
             wimhihihi_host,wimlohihi_host,wimhilohi_host,wimlolohi_host,
             wimhihilo_host,wimlohilo_host,wimhilolo_host,wimlololo_host);
         CPU_normalize
            (wrehihihi_host,wrelohihi_host,wrehilohi_host,wrelolohi_host,
             wrehihilo_host,wrelohilo_host,wrehilolo_host,wrelololo_host,
             wimhihihi_host,wimlohihi_host,wimhilohi_host,wimlolohi_host,
             wimhihilo_host,wimlohilo_host,wimhilolo_host,wimlololo_host,dim,
             vhihihinorm_host,vlohihinorm_host,
             vhilohinorm_host,vlolohinorm_host,
             vhihilonorm_host,vlohilonorm_host,
             vhilolonorm_host,vlololonorm_host);
         CPU_norm
            (wrehihihi_host,wrelohihi_host,wrehilohi_host,wrelolohi_host,
             wrehihilo_host,wrelohilo_host,wrehilolo_host,wrelololo_host,
             wimhihihi_host,wimlohihi_host,wimhilohi_host,wimlolohi_host,
             wimhihilo_host,wimlohilo_host,wimhilolo_host,wimlololo_host,dim,
             &whihihinorm_host,&wlohihinorm_host,
             &whilohinorm_host,&wlolohinorm_host,
             &whihilonorm_host,&wlohilonorm_host,
             &whilolonorm_host,&wlololonorm_host);
      }

   if(mode == 2) // GPU vs CPU correctness verification
   {
      cout << scientific << setprecision(16);

      cout << "CPU norm : " << endl;
      odf_write_doubles(vhihihinorm_host,vlohihinorm_host,
                        vhilohinorm_host,vlolohinorm_host,
                        vhihilonorm_host,vlohilonorm_host,
                        vhilolonorm_host,vlololonorm_host);
      cout << endl;
      cout << "GPU norm : " << endl;
      odf_write_doubles(vhihihinorm_device,vlohihinorm_device,
                        vhilohinorm_device,vlolohinorm_device,
                        vhihilonorm_device,vlohilonorm_device,
                        vhilolonorm_device,vlololonorm_device);
      cout << endl;
      cout << "CPU norm after normalization : " << endl;
      odf_write_doubles(whihihinorm_host,wlohihinorm_host,
                        whilohinorm_host,wlolohinorm_host,
                        whihilonorm_host,wlohilonorm_host,
                        whilolonorm_host,wlololonorm_host);
      cout << endl;
      cout << "GPU norm after normalization : " << endl;
      odf_write_doubles(whihihinorm_device,wlohihinorm_device,
                        whilohinorm_device,wlolohinorm_device,
                        whihilonorm_device,wlohilonorm_device,
                        whilolonorm_device,wlololonorm_device);
      cout << endl;
   }

   return 0;
}
