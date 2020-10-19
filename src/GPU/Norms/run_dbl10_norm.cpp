/* This is the main program to run manual tests,
   with four input parameters at the command line.
   For example, typing at the command prompt
   ./run_dbl10_norm 32 64 1 2
   computes the norm once with 32 threads in a block,
   of a vector of dimension 64, on both GPU and CPU,
   and applies this norm to normalize a random vector.

Including the vector_types.h is needed for the include of
the headers for the dbl10_norm_kernels. */

#include <ctime>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "dbl10_norm_kernels.h"
#include "dbl10_norm_host.h"
#include "random10_vectors.h"
#include "parse_run_arguments.h"
#include "deca_double_functions.h"

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

   double* vrtb_host = new double[dim];   // highest parts on the host
   double* vrix_host = new double[dim];   // 2nd highest parts on the host
   double* vrmi_host = new double[dim];   // 3rd highest parts on the host
   double* vrrg_host = new double[dim];   // 4th highest parts on the host
   double* vrpk_host = new double[dim];   // 5th highest parts on the host
   double* vltb_host = new double[dim];   // 5th lowest parts on the host
   double* vlix_host = new double[dim];   // 4th lowest parts on the host
   double* vlmi_host = new double[dim];   // 3rd lowest parts on the host
   double* vlrg_host = new double[dim];   // 2nd lowest parts on the host
   double* vlpk_host = new double[dim];   // lowest parts on the host
   double* vrtb_device = new double[dim]; // highest parts on the device
   double* vrix_device = new double[dim]; // 2nd highest parts on the device
   double* vrmi_device = new double[dim]; // 3rd highest parts on the device
   double* vrrg_device = new double[dim]; // 4th highest parts on the device
   double* vrpk_device = new double[dim]; // 5th highest parts on the device
   double* vltb_device = new double[dim]; // 5th lowest parts on the device
   double* vlix_device = new double[dim]; // 4th lowest parts on the device
   double* vlmi_device = new double[dim]; // 3rd lowest parts on the device
   double* vlrg_device = new double[dim]; // 2nd lowest parts on the device
   double* vlpk_device = new double[dim]; // lowest parts on the device
   double* wrtb_host = new double[dim];   // highest parts copy
   double* wrix_host = new double[dim];   // 2nd highest parts copy
   double* wrmi_host = new double[dim];   // 3rd highest parts copy
   double* wrrg_host = new double[dim];   // 4th highest parts copy
   double* wrpk_host = new double[dim];   // 5th highest parts copy
   double* wltb_host = new double[dim];   // 5th lowest parts copy
   double* wlix_host = new double[dim];   // 4th lowest parts copy
   double* wlmi_host = new double[dim];   // 3rd lowest parts copy
   double* wlrg_host = new double[dim];   // 2nd lowest parts copy
   double* wlpk_host = new double[dim];   // lowest parts copy

   random_double10_vectors
      (dim,vrtb_host,vrix_host,vrmi_host,vrrg_host,vrpk_host,
           vltb_host,vlix_host,vlmi_host,vlrg_host,vlpk_host,
       vrtb_device,vrix_device,vrmi_device,vrrg_device,vrpk_device,
       vltb_device,vlix_device,vlmi_device,vlrg_device,vlpk_device);

   double vrtbnorm_device,vrixnorm_device,vrrgnorm_device,vrpknorm_device;
   double vrminorm_device;
   double vltbnorm_device,vlixnorm_device,vlrgnorm_device,vlpknorm_device;
   double vlminorm_device;
   double wrtbnorm_device,wrixnorm_device,wrrgnorm_device,wrpknorm_device;
   double wrminorm_device;
   double wltbnorm_device,wlixnorm_device,wlrgnorm_device,wlpknorm_device;
   double wlminorm_device;
   double vrtbnorm_host,vrixnorm_host,vrminorm_host,vrrgnorm_host;
   double vrpknorm_host;
   double vltbnorm_host,vlixnorm_host,vlminorm_host,vlrgnorm_host;
   double vlpknorm_host;
   double wrtbnorm_host,wrixnorm_host,wrminorm_host,wrrgnorm_host;
   double wrpknorm_host;
   double wltbnorm_host,wlixnorm_host,wlminorm_host,wlrgnorm_host;
   double wlpknorm_host;

   if(mode==0 || mode==2) // GPU computation of the norm
   {
      GPU_norm(vrtb_device,vrix_device,vrmi_device,vrrg_device,vrpk_device,
               vltb_device,vlix_device,vlmi_device,vlrg_device,vlpk_device,
               dim,1,BS,
               &vrtbnorm_device,&vrixnorm_device,&vrminorm_device,
               &vrrgnorm_device,&vrpknorm_device,
               &vltbnorm_device,&vlixnorm_device,&vlminorm_device,
               &vlrgnorm_device,&vlpknorm_device,1);
      GPU_norm(vrtb_device,vrix_device,vrmi_device,vrrg_device,vrpk_device,
               vltb_device,vlix_device,vlmi_device,vlrg_device,vlpk_device,
               dim,freq,BS,
               &wrtbnorm_device,&wrixnorm_device,&wrminorm_device,
               &wrrgnorm_device,&wrpknorm_device,
               &wltbnorm_device,&wlixnorm_device,&wlminorm_device,
               &wlrgnorm_device,&wlpknorm_device,1);
   }

   if(mode==1 || mode==2) // CPU computation of the norm
      for(int i=0; i<=freq; i++)
      {
         CPU_norm(vrtb_host,vrix_host,vrmi_host,vrrg_host,vrpk_host,
                  vltb_host,vlix_host,vlmi_host,vlrg_host,vlpk_host,dim,
                  &vrtbnorm_host,&vrixnorm_host,&vrminorm_host,
                  &vrrgnorm_host,&vrpknorm_host,
                  &vltbnorm_host,&vlixnorm_host,&vlminorm_host,
                  &vlrgnorm_host,&vlpknorm_host);
         make_copy(dim,vrtb_host,vrix_host,vrmi_host,vrrg_host,vrpk_host,
                       vltb_host,vlix_host,vlmi_host,vlrg_host,vlpk_host,
                       wrtb_host,wrix_host,wrmi_host,wrrg_host,wrpk_host,
                       wltb_host,wlix_host,wlmi_host,wlrg_host,wlpk_host);
         CPU_normalize(wrtb_host,wrix_host,wrmi_host,wrrg_host,wrpk_host,
                       wltb_host,wlix_host,wlmi_host,wlrg_host,wlpk_host,dim,
                       vrtbnorm_host,vrixnorm_host,vrminorm_host,
                       vrrgnorm_host,vrpknorm_host,
                       vltbnorm_host,vlixnorm_host,vlminorm_host,
                       vlrgnorm_host,vlpknorm_host);
         CPU_norm(wrtb_host,wrix_host,wrmi_host,wrrg_host,wrpk_host,
                  wltb_host,wlix_host,wlmi_host,wlrg_host,wlpk_host,dim,
                  &wrtbnorm_host,&wrixnorm_host,&wrminorm_host,
                  &wrrgnorm_host,&wrpknorm_host,
                  &wltbnorm_host,&wlixnorm_host,&wlminorm_host,
                  &wlrgnorm_host,&wlpknorm_host);
      }

   if(mode == 2) // GPU vs CPU correctness verification
   {
      cout << scientific << setprecision(16);

      cout << "CPU norm : " << endl;
      daf_write_doubles
         (vrtbnorm_host,vrixnorm_host,vrminorm_host,vrrgnorm_host,
          vrpknorm_host,vltbnorm_host,vlixnorm_host,vlminorm_host,
          vlrgnorm_host,vlpknorm_host);
      cout << endl;
      cout << "GPU norm : " << endl;
      daf_write_doubles
         (vrtbnorm_device,vrixnorm_device,vrminorm_device,vrrgnorm_device,
          vrpknorm_device,vltbnorm_device,vlixnorm_device,vlminorm_device,
          vlrgnorm_device,vlpknorm_device);
      cout << endl;
      cout << "CPU norm after normalization : " << endl;
      daf_write_doubles
         (wrtbnorm_host,wrixnorm_host,wrminorm_host,wrrgnorm_host,
          wrpknorm_host,wltbnorm_host,wlixnorm_host,wlminorm_host,
          wlrgnorm_host,wlpknorm_host);
      cout << endl;
      cout << "GPU norm after normalization : " << endl;
      daf_write_doubles
         (wrtbnorm_device,wrixnorm_device,wrminorm_device,wrrgnorm_device,
          wrpknorm_device,wltbnorm_device,wlixnorm_device,wlminorm_device,
          wlrgnorm_device,wlpknorm_device);
      cout << endl;
   }

   return 0;
}
