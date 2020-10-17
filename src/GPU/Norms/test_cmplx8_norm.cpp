/* This program runs automatic tests to compute norms of complex vectors
   in octo double precision, for preset values of the parameters.
   No input is required from the user. */

#include <ctime>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "cmplx8_norm_kernels.h"
#include "cmplx8_norm_host.h"
#include "random8_vectors.h"
#include "octo_double_functions.h"

using namespace std;

void run
 ( int dim, int BS, int freq, int mode, int blocked,
   double *vhihihinorm_device, double *vlohihinorm_device,
   double *vhilohinorm_device, double *vlolohinorm_device,
   double *vhihilonorm_device, double *vlohilonorm_device,
   double *vhilolonorm_device, double *vlololonorm_device,
   double *vhihihinorm_host,   double *vlohihinorm_host,
   double *vhilohinorm_host,   double *vlolohinorm_host,
   double *vhihilonorm_host,   double *vlohilonorm_host,
   double *vhilolonorm_host,   double *vlololonorm_host,
   double *whihihinorm_device, double *wlohihinorm_device,
   double *whilohinorm_device, double *wlolohinorm_device,
   double *whihilonorm_device, double *wlohilonorm_device,
   double *whilolonorm_device, double *wlololonorm_device,
   double *whihihinorm_host,   double *wlohihinorm_host,
   double *whilohinorm_host,   double *wlolohinorm_host,
   double *whihilonorm_host,   double *wlohilonorm_host,
   double *whilolonorm_host,   double *wlololonorm_host );
/*
 * DESCRIPTION :
 *   Computes norms for random complex vectors,
 *   in double doble precision.
 *
 * ON ENTRY :
 *   dim      dimension of the random vector;
 *   BS       block size;
 *   freq     frequency of runs for the tilohings;
 *   mode     if 0, then only GPU computes,
 *            if 1, then only CPU computes,
 *            if 2, then both GPU and CPU  compute;
 *   blocked  if 0, then the vector should be of medium size
 *            and only one block will compute,
 *            if 1, then as many as dim/BS blocks will compute.
 *
 * ON RETURN :
 *   vhihinorm_device  highest part of 2-norm compute by the device;
 *   vlohinorm_device  second highest part of 2-norm compute by the device;
 *   vhilonorm_device  second lowest part of 2-norm computed by the device;
 *   vlolonorm_device  lowest part of 2-norm computed by the device;
 *   vhihinorm_host    highest part of 2-norm compute by the host;
 *   vlohinorm_host    second highest part of 2-norm computed by the host;
 *   vhilonorm_host    second lowest part of 2-norm computed by the host;
 *   vlolonorm_host    lowest part of 2-norm computed by the host;
 *   whihinorm_device  highest part of 2-norm on device after normalization;
 *   wlohinorm_device  2nd highest part of 2-norm on device after normalization;
 *   whilonorm_device  2nd lowest part of 2-norm on device after normalization;
 *   wlolonorm_device  lowest part of 2-norm on device after normalization;
 *   whihinorm_host    highest part of 2-norm on host after normalization;
 *   wlohinorm_host    2nd highest part of 2-norm on host after normalization;
 *   whilonorm_host    2nd lowest part of 2-norm on host after normalization;
 *   wlolonorm_host    lowest part of 2-norm on host after normalization.  */

int verify_correctness ( int dim, int BS, int blocked );
/*
 * Computes norms of real vectors both on GPU and CPU
 * of dimension dim and verifies the correctness.
 * The number of threads in a block is given in the block size BS.
 * If blocked and dim is a multiple of BS, then many blocks compute. */

int main ( void )
{
   int fail;

   cout << "Verifying correctness for dimension and block size 64 ..."
        << endl;
   fail = verify_correctness(64,64,0);

   cout << "Nonblocked version form dimension 1024 and block size 128..."
        << endl;
   fail = verify_correctness(1024,128,0);

   cout << "Blocked version form dimension 4096 and block size 128..."
        << endl;
   fail = verify_correctness(4096,128,1);

   return 0;
}

void run
 ( int dim, int BS, int freq, int mode, int blocked,
   double *vhihihinorm_device, double *vlohihinorm_device,
   double *vhilohinorm_device, double *vlolohinorm_device,
   double *vhihilonorm_device, double *vlohilonorm_device,
   double *vhilolonorm_device, double *vlololonorm_device,
   double *vhihihinorm_host,   double *vlohihinorm_host,
   double *vhilohinorm_host,   double *vlolohinorm_host,
   double *vhihilonorm_host,   double *vlohilonorm_host,
   double *vhilolonorm_host,   double *vlololonorm_host,
   double *whihihinorm_device, double *wlohihinorm_device,
   double *whilohinorm_device, double *wlolohinorm_device,
   double *whihilonorm_device, double *wlohilonorm_device,
   double *whilolonorm_device, double *wlololonorm_device,
   double *whihihinorm_host,   double *wlohihinorm_host,
   double *whilohinorm_host,   double *wlolohinorm_host,
   double *whihilonorm_host,   double *wlohilonorm_host,
   double *whilolonorm_host,   double *wlololonorm_host )
{
   const int timevalue = time(NULL); // no fixed seed to verify correctness
   srand(timevalue);

   double* vrehihihi_host = new double[dim]; // highest real parts on host
   double* vrelohihi_host = new double[dim]; // 2nd highest real parts on host
   double* vrehilohi_host = new double[dim]; // 3rd highest real parts on host
   double* vrelolohi_host = new double[dim]; // 4th highest real parts on host
   double* vrehihilo_host = new double[dim]; // 4th lowest real parts on host
   double* vrelohilo_host = new double[dim]; // 3rd lowest real parts on host
   double* vrehilolo_host = new double[dim]; // 2nd lowest real parts on host
   double* vrelololo_host = new double[dim]; // lowest real parts on host
   double* vimhihihi_host = new double[dim]; // highest imaginary parts on host
   double* vimlohihi_host = new double[dim]; // 2nd highest imag parts on host
   double* vimhilohi_host = new double[dim]; // 3nd highest imag parts on host
   double* vimlolohi_host = new double[dim]; // 4th highest imag parts on host
   double* vimhihilo_host = new double[dim]; // 4th lowest imag parts on host
   double* vimlohilo_host = new double[dim]; // 3rd lowest imag parts on host
   double* vimhilolo_host = new double[dim]; // 2nd lowest imag parts on host
   double* vimlololo_host = new double[dim]; // lowest imaginary parts on host
   double* wrehihihi_host = new double[dim]; // a copy for normalization
   double* wrelohihi_host = new double[dim];
   double* wrehilohi_host = new double[dim];
   double* wrelolohi_host = new double[dim];
   double* wrehihilo_host = new double[dim]; // a copy for normalization
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
   double* vrelololo_device = new double[dim]; // lowest real parts on device
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

   if(mode == 0 || mode == 2)
   {
      GPU_norm
         (vrehihihi_device,vrelohihi_device,vrehilohi_device,vrelolohi_device,
          vrehihilo_device,vrelohilo_device,vrehilolo_device,vrelololo_device,
          vimhihihi_device,vimlohihi_device,vimhilohi_device,vimlolohi_device,
          vimhihilo_device,vimlohilo_device,vimhilolo_device,vimlololo_device,
          dim,1,BS,vhihihinorm_device,vlohihinorm_device,vhilohinorm_device,
          vlolohinorm_device,vhihilonorm_device,vlohilonorm_device,
          vhilolonorm_device,vlololonorm_device,blocked);
      GPU_norm
         (vrehihihi_device,vrelohihi_device,vrehilohi_device,vrelolohi_device,
          vrehihilo_device,vrelohilo_device,vrehilolo_device,vrelololo_device,
          vimhihihi_device,vimlohihi_device,vimhilohi_device,vimlolohi_device,
          vimhihilo_device,vimlohilo_device,vimhilolo_device,vimlololo_device,
          dim,freq,BS,whihihinorm_device,wlohihinorm_device,
          whilohinorm_device,wlolohinorm_device,whihilonorm_device,
          wlohilonorm_device,whilolonorm_device,wlololonorm_device,blocked);
   }

   if(mode == 1 || mode == 2)
   {
      for(int i=0; i<=freq; i++)
      {
         CPU_norm
            (vrehihihi_host,vrelohihi_host,vrehilohi_host,vrelolohi_host,
             vrehihilo_host,vrelohilo_host,vrehilolo_host,vrelololo_host,
             vimhihihi_host,vimlohihi_host,vimhilohi_host,vimlolohi_host,
             vimhihilo_host,vimlohilo_host,vimhilolo_host,vimlololo_host,
             dim,vhihihinorm_host,vlohihinorm_host,vhilohinorm_host,
                 vlolohinorm_host,vhihilonorm_host,vlohilonorm_host,
                 vhilolonorm_host,vlololonorm_host);
         make_copy(dim,
            vrehihihi_host,vrelohihi_host,vrehilohi_host,vrelolohi_host,
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
             *vhihihinorm_host,*vlohihinorm_host,*vhilohinorm_host,
             *vlolohinorm_host,*vhihilonorm_host,*vlohilonorm_host,
             *vhilolonorm_host,*vlololonorm_host);
         CPU_norm
            (wrehihihi_host,wrelohihi_host,wrehilohi_host,wrelolohi_host,
             wrehihilo_host,wrelohilo_host,wrehilolo_host,wrelololo_host,
             wimhihihi_host,wimlohihi_host,wimhilohi_host,wimlolohi_host,
             wimhihilo_host,wimlohilo_host,wimhilolo_host,wimlololo_host,
             dim,whihihinorm_host,wlohihinorm_host,whilohinorm_host,
             wlolohinorm_host,whihilonorm_host,wlohilonorm_host,
             whilolonorm_host,wlololonorm_host);
      }
   }
}

int verify_correctness ( int dim, int BS, int blocked )
{
   double vhihihinorm_device,vlohihinorm_device,vhilohinorm_device;
   double vlolohinorm_device,vhihilonorm_device,vlohilonorm_device;
   double vhilolonorm_device,vlololonorm_device;
   // norm before normalization on device
   double whihihinorm_device,wlohihinorm_device,whilohinorm_device;
   double wlolohinorm_device,whihilonorm_device,wlohilonorm_device;
   double whilolonorm_device,wlololonorm_device;
   // norm after normalization on device
   double vhihihinorm_host,vlohihinorm_host,vhilohinorm_host,vlolohinorm_host;
   double vhihilonorm_host,vlohilonorm_host,vhilolonorm_host,vlololonorm_host;
   // norm before normalization on host
   double whihihinorm_host,wlohihinorm_host,whilohinorm_host,wlolohinorm_host;
   double whihilonorm_host,wlohilonorm_host,whilolonorm_host,wlololonorm_host;
   // norm after normalization on host

   run(dim,BS,1,2,blocked,
       &vhihihinorm_device,&vlohihinorm_device,&vhilohinorm_device,
       &vlolohinorm_device,&vhihilonorm_device,&vlohilonorm_device,
       &vhilolonorm_device,&vlololonorm_device,
       &vhihihinorm_host,  &vlohihinorm_host,  &vhilohinorm_host,
       &vlolohinorm_host,  &vhihilonorm_host,  &vlohilonorm_host,
       &vhilolonorm_host,  &vlololonorm_host,
       &whihihinorm_device,&wlohihinorm_device,&whilohinorm_device,
       &wlolohinorm_device,&whihilonorm_device,&wlohilonorm_device,
       &whilolonorm_device,&wlololonorm_device,
       &whihihinorm_host,  &wlohihinorm_host,  &whilohinorm_host,
       &wlolohinorm_host,  &whihilonorm_host,  &wlohilonorm_host,
       &whilolonorm_host,  &wlololonorm_host);

   cout << scientific << setprecision(16);

   cout << "CPU norm : " << endl;
   odf_write_doubles
      (vhihihinorm_host,vlohihinorm_host,vhilohinorm_host,vlolohinorm_host,
       vhihilonorm_host,vlohilonorm_host,vhilolonorm_host,vlololonorm_host);
   cout << "GPU norm : " << endl;
   odf_write_doubles
      (vhihihinorm_device,vlohihinorm_device,vhilohinorm_device,
       vlolohinorm_device,vhihilonorm_device,vlohilonorm_device,
       vhilolonorm_device,vlololonorm_device);
   cout << "CPU norm after normalization : " << endl;
   odf_write_doubles
      (whihihinorm_host,wlohihinorm_host,whilohinorm_host,wlolohinorm_host,
       whihilonorm_host,wlohilonorm_host,whilolonorm_host,wlololonorm_host);
   cout << "GPU norm after normalization : " << endl;
   odf_write_doubles
      (whihihinorm_device,wlohihinorm_device,whilohinorm_device,
       wlolohinorm_device,whihilonorm_device,wlohilonorm_device,
       whilolonorm_device,wlololonorm_device);

   const double tol = 1.0e-100;
   double err = abs(vhihihinorm_device - vhihihinorm_host)
              + abs(whihihinorm_device - whihihinorm_host)
              + abs(vlohihinorm_device - vlohihinorm_host)
              + abs(vlohihinorm_device - vlohihinorm_host)
              + abs(vhilohinorm_device - vhilohinorm_host)
              + abs(whilohinorm_device - whilohinorm_host)
              + abs(vlolohinorm_device - vlolohinorm_host)
              + abs(wlolohinorm_device - wlolohinorm_host)
              + abs(vhihilonorm_device - vhihilonorm_host)
              + abs(whihilonorm_device - whihilonorm_host)
              + abs(vlohilonorm_device - vlohilonorm_host)
              + abs(vlohilonorm_device - vlohilonorm_host)
              + abs(vhilolonorm_device - vhilolonorm_host)
              + abs(whilolonorm_device - whilolonorm_host)
              + abs(vlololonorm_device - vlololonorm_host)
              + abs(wlololonorm_device - wlololonorm_host);

   cout << scientific << setprecision(4) << "error : " << err;
   if(err <= tol)
   {
      cout << " <= " << tol << "  okay" << endl;
      return 0;
   }
   else
   {
      cout << " > " << tol << "  failure" << endl;
      return 1;
   }
}
