/* This is the main program to run manual tests,
   with four input parameters at the command line.
   For example, typing at the command prompt
   ./run_dbl2_norm 32 32 1 2
   computes the norm once with 32 threads in a block,
   of a vector of dimension 32, on both GPU and CPU,
   and applies this norm to normalize a random vector.

Including the vector_types.h is needed for the include of
the headers for the dbl2_norm_kernels. */

#include <ctime>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "dbl2_norm_kernels.h"
#include "dbl2_norm_host.h"
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

   double* vhi_host = new double[dim];   // high parts on the host
   double* vlo_host = new double[dim];   // low parts on the host
   double* vhi_device = new double[dim]; // high parts on the device
   double* vlo_device = new double[dim]; // low parts on the device
   double* whi_host = new double[dim];   // high parts copy for normalization
   double* wlo_host = new double[dim];   // low parts copy for normalization

   random_double2_vectors(dim,vhi_host,vlo_host,vhi_device,vlo_device);

   double vhinorm_device,vlonorm_device,whinorm_device,wlonorm_device;
   double vhinorm_host,vlonorm_host,whinorm_host,wlonorm_host;

   if(mode==0 || mode==2) // GPU computation of the norm
   {
      GPU_norm(vhi_device,vlo_device,dim,1,BS,
               &vhinorm_device,&vlonorm_device,1);
      GPU_norm(vhi_device,vlo_device,dim,freq,BS,
               &whinorm_device,&wlonorm_device,1);
   }

   if(mode==1 || mode==2) // CPU computation of the norm
      for(int i=0; i<=freq; i++)
      {
         CPU_norm(vhi_host,vlo_host,dim,&vhinorm_host,&vlonorm_host);
         make_copy(dim,vhi_host,vlo_host,whi_host,wlo_host);
         CPU_normalize(whi_host,wlo_host,dim,vhinorm_host,vlonorm_host);
         CPU_norm(whi_host,wlo_host,dim,&whinorm_host,&wlonorm_host);
      }

   if(mode == 2) // GPU vs CPU correctness verification
   {
      cout << scientific << setprecision(16);

      cout << "   CPU norm : " << endl;
      cout << "     hi : " << vhinorm_host;
      cout << "     lo : " << vlonorm_host << endl;
      cout << "   GPU norm : " << endl;
      cout << "     hi : " << vhinorm_device;
      cout << "     lo : " << vlonorm_device << endl;
      cout << "   CPU norm after normalization : " << endl;
      cout << "     hi : " << whinorm_host;
      cout << "     lo : " << wlonorm_host << endl;
      cout << "   GPU norm after normalization : " << endl;
      cout << "     hi : " << whinorm_device;
      cout << "     lo : " << wlonorm_device << endl;
   }

   return 0;
}
