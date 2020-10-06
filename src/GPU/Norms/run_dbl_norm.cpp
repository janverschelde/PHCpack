/* This is the main program to run manual tests,
   with four input parameters at the command line.
   For example, typing at the command prompt
   ./run_dbl_norm_d 32 32 1 2
   computes the norm once with 32 threads in a block,
   of a vector of dimension 32, on both GPU and CPU,
   and applies this norm to normalize a random vector.

Including the vector_types.h is needed for the include of
the headers for the dbl_norm_kernels. */

#include <ctime>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "dbl_norm_kernels.h"
#include "dbl_norm_host.h"
#include "random_vectors.h"
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

   double* v_host = new double[dim];   // vector on the host
   double* w_host = new double[dim];   // a copy for normalization
   double* v_device = new double[dim]; // vector on the device

   random_double_vectors(dim,v_host,v_device);

   double vnorm_device,wnorm_device;
   double vnorm_host,wnorm_host;

   if(mode==0 || mode==2) // GPU computation of the norm
   {
      GPU_norm(v_device,dim,1,BS,&vnorm_device,1);
      GPU_norm(v_device,dim,freq,BS,&wnorm_device,1);
   }

   if(mode==1 || mode==2) // CPU computation of the norm
      for(int i=0; i<=freq; i++)
      {
         CPU_norm(v_host,dim,&vnorm_host);
         make_copy(dim,v_host,w_host);
         CPU_normalize(w_host,dim,vnorm_host);
         CPU_norm(w_host,dim,&wnorm_host);
      }

   if(mode == 2) // GPU vs CPU correctness verification
   {
      cout << scientific << setprecision(16);
      cout << "GPU norm : " << vnorm_device << endl;
      cout << "GPU norm after normalization : " << wnorm_device << endl;
      cout << "CPU norm : " << vnorm_host << endl;
      cout << "CPU norm after normalization : " << wnorm_host << endl;
   }

   return 0;
}
