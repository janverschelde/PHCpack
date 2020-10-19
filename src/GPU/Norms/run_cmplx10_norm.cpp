/* This is the main program to run manual tests,
   with four input parameters at the command line.

For example, typing at the command prompt
  ./run_cmplx10_norm_d 64 256 1 2
computes the norm once with 64 threads in a block,
of a vector of dimension 128, on both GPU and CPU,
and applies this norm to normalize a random vector.
As all complex numbers are on the complex unit circle,
the 2-norm of a vector of dimension 64 equals 8.
If the dimension is 256, then the 2-norm is 16.

Including the vector_types.h is needed for the include of
the headers for the cmplx10_norm_kernels. */

#include <ctime>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "cmplx10_norm_kernels.h"
#include "cmplx10_norm_host.h"
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

   double* vrertb_host = new double[dim];   // highest real parts on host
   double* vrerix_host = new double[dim];   // 2nd highest real parts on host
   double* vrermi_host = new double[dim];   // 3rd highest real parts on host
   double* vrerrg_host = new double[dim];   // 4th highest real parts on host
   double* vrerpk_host = new double[dim];   // 5th highest real parts on host
   double* vreltb_host = new double[dim];   // 5th lowest real parts on host
   double* vrelix_host = new double[dim];   // 4th lowest real parts on host
   double* vrelmi_host = new double[dim];   // 3rd lowest real parts on host
   double* vrelrg_host = new double[dim];   // 2nd lowest real parts on host
   double* vrelpk_host = new double[dim];   // lowest real parts on host
   double* vimrtb_host = new double[dim];   // highest imaginary parts 
   double* vimrix_host = new double[dim];   // 2nd highest imaginary parts
   double* vimrmi_host = new double[dim];   // 3rd highest imaginary parts
   double* vimrrg_host = new double[dim];   // 4th highest imaginary parts
   double* vimrpk_host = new double[dim];   // 5th highest imaginary parts
   double* vimltb_host = new double[dim];   // 5th lowest imaginary parts 
   double* vimlix_host = new double[dim];   // 4th lowest imaginary parts
   double* vimlmi_host = new double[dim];   // 3rd lowest imaginary parts
   double* vimlrg_host = new double[dim];   // 2nd lowest imaginary parts
   double* vimlpk_host = new double[dim];   // lowest imaginary parts
   double* wrertb_host = new double[dim];   // a copy for normalization
   double* wrerix_host = new double[dim];
   double* wrermi_host = new double[dim];
   double* wrerrg_host = new double[dim];
   double* wrerpk_host = new double[dim];
   double* wreltb_host = new double[dim];
   double* wrelix_host = new double[dim];
   double* wrelmi_host = new double[dim];
   double* wrelrg_host = new double[dim];
   double* wrelpk_host = new double[dim];
   double* wimrtb_host = new double[dim];
   double* wimrix_host = new double[dim];
   double* wimrmi_host = new double[dim];
   double* wimrrg_host = new double[dim];
   double* wimrpk_host = new double[dim];
   double* wimltb_host = new double[dim];
   double* wimlix_host = new double[dim];
   double* wimlmi_host = new double[dim];
   double* wimlrg_host = new double[dim];
   double* wimlpk_host = new double[dim];
   double* vrertb_device = new double[dim]; // highest real parts on device
   double* vrerix_device = new double[dim]; // 2nd highest real parts
   double* vrermi_device = new double[dim]; // 3rd highest real parts
   double* vrerrg_device = new double[dim]; // 4th highest real parts
   double* vrerpk_device = new double[dim]; // 5th highest real parts
   double* vreltb_device = new double[dim]; // 5th lowest real parts on device
   double* vrelix_device = new double[dim]; // 4th lowest real parts
   double* vrelmi_device = new double[dim]; // 3rd lowest real parts
   double* vrelrg_device = new double[dim]; // 2nd lowest real parts
   double* vrelpk_device = new double[dim]; // lowest real parts
   double* vimrtb_device = new double[dim]; // highest imaginary parts
   double* vimrix_device = new double[dim]; // 2nd highest imaginary parts
   double* vimrmi_device = new double[dim]; // 3rd highest imaginary parts
   double* vimrrg_device = new double[dim]; // 4th highest imaginary parts
   double* vimrpk_device = new double[dim]; // 5th highest imaginary parts
   double* vimltb_device = new double[dim]; // 5th lowest imaginary parts
   double* vimlix_device = new double[dim]; // 4th lowest imaginary parts
   double* vimlmi_device = new double[dim]; // 3rd lowest imaginary parts
   double* vimlrg_device = new double[dim]; // 2nd lowest imaginary parts
   double* vimlpk_device = new double[dim]; // lowest imaginary parts

   random_complex10_vectors
      (dim,vrertb_host,vrerix_host,vrermi_host,vrerrg_host,vrerpk_host,
           vreltb_host,vrelix_host,vrelmi_host,vrelrg_host,vrelpk_host,
           vimrtb_host,vimrix_host,vimrmi_host,vimrrg_host,vimrpk_host,
           vimltb_host,vimlix_host,vimlmi_host,vimlrg_host,vimlpk_host,
       vrertb_device,vrerix_device,vrermi_device,vrerrg_device,vrerpk_device,
       vreltb_device,vrelix_device,vrelmi_device,vrelrg_device,vrelpk_device,
       vimrtb_device,vimrix_device,vimrmi_device,vimrrg_device,vimrpk_device,
       vimltb_device,vimlix_device,vimlmi_device,vimlrg_device,vimlpk_device);

   double vrtbnorm_device,vrixnorm_device,vrminorm_device;
   double vrrgnorm_device,vrpknorm_device;
   double vltbnorm_device,vlixnorm_device,vlminorm_device;
   double vlrgnorm_device,vlpknorm_device;
   double wrtbnorm_device,wrixnorm_device,wrminorm_device;
   double wrrgnorm_device,wrpknorm_device;
   double wltbnorm_device,wlixnorm_device,wlminorm_device;
   double wlrgnorm_device,wlpknorm_device;
   double vrtbnorm_host,vrixnorm_host,vrminorm_host;
   double vrrgnorm_host,vrpknorm_host;
   double vltbnorm_host,vlixnorm_host,vlminorm_host;
   double vlrgnorm_host,vlpknorm_host;
   double wrtbnorm_host,wrixnorm_host,wrminorm_host;
   double wrrgnorm_host,wrpknorm_host;
   double wltbnorm_host,wlixnorm_host,wlminorm_host;
   double wlrgnorm_host,wlpknorm_host;

   if(mode==0 || mode==2) // GPU computation of the norm
   {
      GPU_norm
         (vrertb_device,vrerix_device,vrermi_device,vrerrg_device,
          vrerpk_device,vreltb_device,vrelix_device,vrelmi_device,
          vrelrg_device,vrelpk_device,vimrtb_device,vimrix_device,
          vimrmi_device,vimrrg_device,vimrpk_device,vimltb_device,
          vimlix_device,vimlmi_device,vimlrg_device,vimlpk_device,
          dim,1,BS,
          &vrtbnorm_device,&vrixnorm_device,&vrminorm_device,
          &vrrgnorm_device,&vrpknorm_device,
          &vltbnorm_device,&vlixnorm_device,&vlminorm_device,
          &vlrgnorm_device,&vlpknorm_device,1);
      GPU_norm
         (vrertb_device,vrerix_device,vrermi_device,vrerrg_device,
          vrerpk_device,vreltb_device,vrelix_device,vrelmi_device,
          vrelrg_device,vrelpk_device,vimrtb_device,vimrix_device,
          vimrmi_device,vimrrg_device,vimrpk_device,vimltb_device,
          vimlix_device,vimlmi_device,vimlrg_device,vimlpk_device,
          dim,freq,BS,
          &wrtbnorm_device,&wrixnorm_device,&vrminorm_device,
          &wrrgnorm_device,&wrpknorm_device,
          &wltbnorm_device,&wlixnorm_device,&vlminorm_device,
          &wlrgnorm_device,&wlpknorm_device,1);
   }

   if(mode==1 || mode==2) // CPU computation of the norm
      for(int i=0; i<=freq; i++)
      {
         CPU_norm
            (vrertb_host,vrerix_host,vrermi_host,vrerrg_host,vrerpk_host,
             vreltb_host,vrelix_host,vrelmi_host,vrelrg_host,vrelpk_host,
             vimrtb_host,vimrix_host,vimrmi_host,vimrrg_host,vimrpk_host,
             vimltb_host,vimlix_host,vimlmi_host,vimlrg_host,vimlpk_host,dim,
             &vrtbnorm_host,&vrixnorm_host,&vrminorm_host,
             &vrrgnorm_host,&vrpknorm_host,
             &vltbnorm_host,&vlixnorm_host,&vlminorm_host,
             &vlrgnorm_host,&vlpknorm_host);
         make_copy
            (dim,vrertb_host,vrerix_host,vrermi_host,vrerrg_host,vrerpk_host,
                 vreltb_host,vrelix_host,vrelmi_host,vrelrg_host,vrelpk_host,
                 vimrtb_host,vimrix_host,vimrmi_host,vimrrg_host,vimrpk_host,
                 vimltb_host,vimlix_host,vimlmi_host,vimlrg_host,vimlpk_host,
             wrertb_host,wrerix_host,wrermi_host,wrerrg_host,wrerpk_host,
             wreltb_host,wrelix_host,wrelmi_host,wrelrg_host,wrelpk_host,
             wimrtb_host,wimrix_host,wimrmi_host,wimrrg_host,wimrpk_host,
             wimltb_host,wimlix_host,wimlmi_host,wimlrg_host,wimlpk_host);
         CPU_normalize
            (wrertb_host,wrerix_host,wrermi_host,wrerrg_host,wrerpk_host,
             wreltb_host,wrelix_host,wrelmi_host,wrelrg_host,wrelpk_host,
             wimrtb_host,wimrix_host,wimrmi_host,wimrrg_host,wimrpk_host,
             wimltb_host,wimlix_host,wimlmi_host,wimlrg_host,wimlpk_host,dim,
             vrtbnorm_host,vrixnorm_host,vrminorm_host,vrrgnorm_host,
             vrpknorm_host,vltbnorm_host,vlixnorm_host,vlminorm_host,
             vlrgnorm_host,vlpknorm_host);
         CPU_norm
            (wrertb_host,wrerix_host,wrermi_host,wrerrg_host,wrerpk_host,
             wreltb_host,wrelix_host,wrelmi_host,wrelrg_host,wrelpk_host,
             wimrtb_host,wimrix_host,wimrmi_host,wimrrg_host,wimrpk_host,
             wimltb_host,wimlix_host,wimlmi_host,wimlrg_host,wimlpk_host,dim,
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
