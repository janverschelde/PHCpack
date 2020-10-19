/* This is the main program to run manual tests,
   with four input parameters at the command line.
   For example, typing at the command prompt
   ./run_dbl8_norm 32 64 1 2
   computes the norm once with 32 threads in a block,
   of a vector of dimension 64, on both GPU and CPU,
   and applies this norm to normalize a random vector.

Including the vector_types.h is needed for the include of
the headers for the dbl8_norm_kernels. */

#include <ctime>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "dbl8_norm_kernels.h"
#include "dbl8_norm_host.h"
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

   double* vhihihi_host = new double[dim];   // highest parts on the host
   double* vlohihi_host = new double[dim];   // 2nd highest parts on the host
   double* vhilohi_host = new double[dim];   // 3rd highest parts on the host
   double* vlolohi_host = new double[dim];   // 4th highest parts on the host
   double* vhihilo_host = new double[dim];   // 4th lowest parts on the host
   double* vlohilo_host = new double[dim];   // 3rd lowest parts on the host
   double* vhilolo_host = new double[dim];   // 2nd lowest parts on the host
   double* vlololo_host = new double[dim];   // lowest parts on the host
   double* vhihihi_device = new double[dim]; // highest parts on the device
   double* vlohihi_device = new double[dim]; // 2nd highest parts
   double* vhilohi_device = new double[dim]; // 3rd highest parts
   double* vlolohi_device = new double[dim]; // 4th highest parts
   double* vhihilo_device = new double[dim]; // 4th lowest parts on the device
   double* vlohilo_device = new double[dim]; // 3nd lowest parts
   double* vhilolo_device = new double[dim]; // 2nd lowest parts
   double* vlololo_device = new double[dim]; // lowest parts on the device
   double* whihihi_host = new double[dim];   // highest parts copy
   double* wlohihi_host = new double[dim];   // 2nd highest parts copy
   double* whilohi_host = new double[dim];   // 3rd highest parts copy
   double* wlolohi_host = new double[dim];   // 4th highest parts copy
   double* whihilo_host = new double[dim];   // 4th lowest parts copy
   double* wlohilo_host = new double[dim];   // 3rd lowest parts copy
   double* whilolo_host = new double[dim];   // 2nd lowest parts copy
   double* wlololo_host = new double[dim];   // lowest parts copy

   random_double8_vectors
      (dim,vhihihi_host,vlohihi_host,vhilohi_host,vlolohi_host,
           vhihilo_host,vlohilo_host,vhilolo_host,vlololo_host,
       vhihihi_device,vlohihi_device,vhilohi_device,vlolohi_device,
       vhihilo_device,vlohilo_device,vhilolo_device,vlololo_device);

   double vhihihinorm_device,vlohihinorm_device;
   double vhilohinorm_device,vlolohinorm_device;
   double vhihilonorm_device,vlohilonorm_device;
   double vhilolonorm_device,vlololonorm_device;
   double whihihinorm_device,wlohihinorm_device;
   double whilohinorm_device,wlolohinorm_device;
   double whihilonorm_device,wlohilonorm_device;
   double whilolonorm_device,wlololonorm_device;
   double vhihihinorm_host,vlohihinorm_host;
   double vhilohinorm_host,vlolohinorm_host;
   double vhihilonorm_host,vlohilonorm_host;
   double vhilolonorm_host,vlololonorm_host;
   double whihihinorm_host,wlohihinorm_host;
   double whilohinorm_host,wlolohinorm_host;
   double whihilonorm_host,wlohilonorm_host;
   double whilolonorm_host,wlololonorm_host;

   if(mode==0 || mode==2) // GPU computation of the norm
   {
      GPU_norm(vhihihi_device,vlohihi_device,vhilohi_device,vlolohi_device,
               vhihilo_device,vlohilo_device,vhilolo_device,vlololo_device,
               dim,1,BS,
               &vhihihinorm_device,&vlohihinorm_device,
               &vhilohinorm_device,&vlolohinorm_device,
               &vhihilonorm_device,&vlohilonorm_device,
               &vhilolonorm_device,&vlololonorm_device,1);
      GPU_norm(vhihihi_device,vlohihi_device,vhilohi_device,vlolohi_device,
               vhihilo_device,vlohilo_device,vhilolo_device,vlololo_device,
               dim,freq,BS,
               &whihihinorm_device,&wlohihinorm_device,
               &whilohinorm_device,&wlolohinorm_device,
               &whihilonorm_device,&wlohilonorm_device,
               &whilolonorm_device,&wlololonorm_device,1);
   }

   if(mode==1 || mode==2) // CPU computation of the norm
      for(int i=0; i<=freq; i++)
      {
         CPU_norm(vhihihi_host,vlohihi_host,vhilohi_host,vlolohi_host,
                  vhihilo_host,vlohilo_host,vhilolo_host,vlololo_host,dim,
                  &vhihihinorm_host,&vlohihinorm_host,
                  &vhilohinorm_host,&vlolohinorm_host,
                  &vhihilonorm_host,&vlohilonorm_host,
                  &vhilolonorm_host,&vlololonorm_host);
         make_copy
            (dim,vhihihi_host,vlohihi_host,vhilohi_host,vlolohi_host,
                 vhihilo_host,vlohilo_host,vhilolo_host,vlololo_host,
                 whihihi_host,wlohihi_host,whilohi_host,wlolohi_host,
                 whihilo_host,wlohilo_host,whilolo_host,wlololo_host);
         CPU_normalize
            (whihihi_host,wlohihi_host,whilohi_host,wlolohi_host,
             whihilo_host,wlohilo_host,whilolo_host,wlololo_host,dim,
             vhihihinorm_host,vlohihinorm_host,
             vhilohinorm_host,vlolohinorm_host,
             vhihilonorm_host,vlohilonorm_host,
             vhilolonorm_host,vlololonorm_host);
         CPU_norm(whihihi_host,wlohihi_host,whilohi_host,wlolohi_host,
                  whihilo_host,wlohilo_host,whilolo_host,wlololo_host,dim,
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
