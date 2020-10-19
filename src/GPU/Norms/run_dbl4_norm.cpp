/* This is the main program to run manual tests,
   with four input parameters at the command line.
   For example, typing at the command prompt
   ./run_dbl4_norm 32 64 1 2
   computes the norm once with 32 threads in a block,
   of a vector of dimension 64, on both GPU and CPU,
   and applies this norm to normalize a random vector.

Including the vector_types.h is needed for the include of
the headers for the dbl4_norm_kernels. */

#include <ctime>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "dbl4_norm_kernels.h"
#include "dbl4_norm_host.h"
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

   double* vhihi_host = new double[dim];   // highest parts on the host
   double* vlohi_host = new double[dim];   // 2nd highest parts on the host
   double* vhilo_host = new double[dim];   // 2nd lowest parts on the host
   double* vlolo_host = new double[dim];   // lowest parts on the host
   double* vhihi_device = new double[dim]; // highest parts on the device
   double* vlohi_device = new double[dim]; // 2nd highest parts on the device
   double* vhilo_device = new double[dim]; // 2nd lowest parts on the device
   double* vlolo_device = new double[dim]; // lowest parts on the device
   double* whihi_host = new double[dim];   // highest parts copy
   double* wlohi_host = new double[dim];   // 2nd highest parts copy
   double* whilo_host = new double[dim];   // 2nd lowest parts copy
   double* wlolo_host = new double[dim];   // lowest parts copy

   random_double4_vectors
      (dim,vhihi_host,vlohi_host,vhilo_host,vlolo_host,
       vhihi_device,vlohi_device,vhilo_device,vlolo_device);

   double vhihinorm_device,vlohinorm_device,vhilonorm_device,vlolonorm_device;
   double whihinorm_device,wlohinorm_device,whilonorm_device,wlolonorm_device;
   double vhihinorm_host,vlohinorm_host,vhilonorm_host,vlolonorm_host;
   double whihinorm_host,wlohinorm_host,whilonorm_host,wlolonorm_host;

   if(mode==0 || mode==2) // GPU computation of the norm
   {
      GPU_norm(vhihi_device,vlohi_device,vhilo_device,vlolo_device,dim,1,BS,
               &vhihinorm_device,&vlohinorm_device,
               &vhilonorm_device,&vlolonorm_device,1);
      GPU_norm(vhihi_device,vlohi_device,vhilo_device,vlolo_device,
               dim,freq,BS,
               &whihinorm_device,&wlohinorm_device,
               &whilonorm_device,&vlolonorm_device,1);
   }

   if(mode==1 || mode==2) // CPU computation of the norm
      for(int i=0; i<=freq; i++)
      {
         CPU_norm(vhihi_host,vlohi_host,vhilo_host,vlolo_host,dim,
                  &vhihinorm_host,&vlohinorm_host,
                  &vhilonorm_host,&vlolonorm_host);
         make_copy(dim,vhihi_host,vlohi_host,vhilo_host,vlolo_host,
                       whihi_host,wlohi_host,whilo_host,wlolo_host);
         CPU_normalize(whihi_host,wlohi_host,whilo_host,wlolo_host,dim,
                       vhihinorm_host,vlohinorm_host,
                       vhilonorm_host,vlolonorm_host);
         CPU_norm(whihi_host,wlohi_host,whilo_host,wlolo_host,dim,
                  &whihinorm_host,&wlohinorm_host,
                  &whilonorm_host,&wlolonorm_host);
      }

   if(mode == 2) // GPU vs CPU correctness verification
   {
      cout << scientific << setprecision(16);

      cout << "CPU norm : " << endl;
      qdf_write_doubles(vhihinorm_host,vlohinorm_host,
                        vhilonorm_host,vlolonorm_host);
      cout << endl;
      cout << "GPU norm : " << endl;
      qdf_write_doubles(vhihinorm_device,vlohinorm_device,
                        vhilonorm_device,vlolonorm_device);
      cout << endl;
      cout << "CPU norm after normalization : " << endl;
      qdf_write_doubles(whihinorm_host,wlohinorm_host,
                        whilonorm_host,wlolonorm_host);
      cout << endl;
      cout << "GPU norm after normalization : " << endl;
      qdf_write_doubles(whihinorm_device,wlohinorm_device,
                        whilonorm_device,wlolonorm_device);
      cout << endl;
   }

   return 0;
}
