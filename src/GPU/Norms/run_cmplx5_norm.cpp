/* This is the main program to run manual tests,
   with four input parameters at the command line.

For example, typing at the command prompt
  ./run_cmplx5_norm_d 64 256 1 2
computes the norm once with 64 threads in a block,
of a vector of dimension 128, on both GPU and CPU,
and applies this norm to normalize a random vector.
As all complex numbers are on the complex unit circle,
the 2-norm of a vector of dimension 64 equals 8.
If the dimension is 256, then the 2-norm is 16.

Including the vector_types.h is needed for the include of
the headers for the cmplx5_norm_kernels. */

#include <ctime>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "cmplx5_norm_kernels.h"
#include "cmplx5_norm_host.h"
#include "random5_vectors.h"
#include "parse_run_arguments.h"
#include "penta_double_functions.h"

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

   double* vretb_host = new double[dim];   // highest real parts on host
   double* vreix_host = new double[dim];   // 2nd highest real parts on host
   double* vremi_host = new double[dim];   // middle real parts on host
   double* vrerg_host = new double[dim];   // 2nd lowest real parts on host
   double* vrepk_host = new double[dim];   // lowest real parts on host
   double* vimtb_host = new double[dim];   // highest imaginary parts 
   double* vimix_host = new double[dim];   // 2nd highest imaginary parts
   double* vimmi_host = new double[dim];   // middle imaginary parts
   double* vimrg_host = new double[dim];   // 2nd lowest imaginary parts
   double* vimpk_host = new double[dim];   // lowest imaginary parts
   double* wretb_host = new double[dim];   // a copy for normalization
   double* wreix_host = new double[dim];
   double* wremi_host = new double[dim];
   double* wrerg_host = new double[dim];
   double* wrepk_host = new double[dim];
   double* wimtb_host = new double[dim];
   double* wimix_host = new double[dim];
   double* wimmi_host = new double[dim];
   double* wimrg_host = new double[dim];
   double* wimpk_host = new double[dim];
   double* vretb_device = new double[dim]; // highest real parts on device
   double* vreix_device = new double[dim]; // 2nd highest real parts
   double* vremi_device = new double[dim]; // middel real parts
   double* vrerg_device = new double[dim]; // 2nd lowest real parts
   double* vrepk_device = new double[dim]; // lowest real parts
   double* vimtb_device = new double[dim]; // highest imaginary parts
   double* vimix_device = new double[dim]; // 2nd highest imaginary parts
   double* vimmi_device = new double[dim]; // middle imaginary parts
   double* vimrg_device = new double[dim]; // 2nd lowest imaginary parts
   double* vimpk_device = new double[dim]; // lowest imaginary parts

   random_complex5_vectors
      (dim,vretb_host,vreix_host,vremi_host,vrerg_host,vrepk_host,
           vimtb_host,vimix_host,vimmi_host,vimrg_host,vimpk_host,
       vretb_device,vreix_device,vremi_device,vrerg_device,vrepk_device,
       vimtb_device,vimix_device,vimmi_device,vimrg_device,vimpk_device);

   double vtbnorm_device,vixnorm_device,vminorm_device;
   double vrgnorm_device,vpknorm_device;
   double wtbnorm_device,wixnorm_device,wminorm_device;
   double wrgnorm_device,wpknorm_device;
   double vtbnorm_host,vixnorm_host,vminorm_host,vrgnorm_host,vpknorm_host;
   double wtbnorm_host,wixnorm_host,wminorm_host,wrgnorm_host,wpknorm_host;

   if(mode==0 || mode==2) // GPU computation of the norm
   {
      GPU_norm
         (vretb_device,vreix_device,vremi_device,vrerg_device,vrepk_device,
          vimtb_device,vimix_device,vimmi_device,vimrg_device,vimpk_device,
          dim,1,BS,
          &vtbnorm_device,&vixnorm_device,&vminorm_device,
          &vrgnorm_device,&vpknorm_device,1);
      GPU_norm
         (vretb_device,vreix_device,vremi_device,vrerg_device,vrepk_device,
          vimtb_device,vimix_device,vimmi_device,vimrg_device,vimpk_device,
          dim,freq,BS,
          &wtbnorm_device,&wixnorm_device,&wminorm_device,
          &wrgnorm_device,&wpknorm_device,1);
   }

   if(mode==1 || mode==2) // CPU computation of the norm
      for(int i=0; i<=freq; i++)
      {
         CPU_norm(vretb_host,vreix_host,vremi_host,vrerg_host,vrepk_host,
                  vimtb_host,vimix_host,vimmi_host,vimrg_host,vimpk_host,dim,
                  &vtbnorm_host,&vixnorm_host,&vminorm_host,
                  &vrgnorm_host,&vpknorm_host);
         make_copy(dim,vretb_host,vreix_host,vremi_host,vrerg_host,vrepk_host,
                       vimtb_host,vimix_host,vimmi_host,vimrg_host,vimpk_host,
                   wretb_host,wreix_host,wremi_host,wrerg_host,wrepk_host,
                   wimtb_host,wimix_host,wimmi_host,wimrg_host,wimpk_host);
         CPU_normalize
            (wretb_host,wreix_host,wremi_host,wrerg_host,wrepk_host,
             wimtb_host,wimix_host,wimmi_host,wimrg_host,wimpk_host,dim,
             vtbnorm_host,vixnorm_host,vminorm_host,vrgnorm_host,vpknorm_host);
         CPU_norm(wretb_host,wreix_host,wremi_host,wrerg_host,wrepk_host,
                  wimtb_host,wimix_host,wimmi_host,wimrg_host,wimpk_host,dim,
                  &wtbnorm_host,&wixnorm_host,&wminorm_host,
                  &wrgnorm_host,&wpknorm_host);
      }

   if(mode == 2) // GPU vs CPU correctness verification
   {
      cout << scientific << setprecision(16);

      cout << "CPU norm : " << endl;
      pdf_write_doubles(vtbnorm_host,vixnorm_host,vminorm_host,
                        vrgnorm_host,vpknorm_host);
      cout << endl;
      cout << "GPU norm : " << endl;
      pdf_write_doubles(vtbnorm_device,vixnorm_device,vminorm_device,
                        vrgnorm_device,vpknorm_device);
      cout << endl;
      cout << "CPU norm after normalization : " << endl;
      pdf_write_doubles(wtbnorm_host,wixnorm_host,wminorm_host,
                        wrgnorm_host,wpknorm_host);
      cout << endl;
      cout << "GPU norm after normalization : " << endl;
      pdf_write_doubles(wtbnorm_device,wixnorm_device,wminorm_device,
                        wrgnorm_device,wpknorm_device);
      cout << endl;
   }

   return 0;
}
