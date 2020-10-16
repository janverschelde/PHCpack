/* This program runs automatic tests to compute norms of complex vectors
   in quad double precision, for preset values of the parameters.
   No input is required from the user. */

#include <ctime>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "cmplx4_norm_kernels.h"
#include "cmplx4_norm_host.h"
#include "random4_vectors.h"
#include "quad_double_functions.h"

using namespace std;

void run
 ( int dim, int BS, int freq, int mode, int blocked,
   double *vhihinorm_device, double *vlohinorm_device,
   double *vhilonorm_device, double *vlolonorm_device,
   double *vhihinorm_host,   double *vlohinorm_host,
   double *vhilonorm_host,   double *vlolonorm_host,
   double *whihinorm_device, double *wlohinorm_device,
   double *whilonorm_device, double *wlolonorm_device,
   double *whihinorm_host,   double *wlohinorm_host,
   double *whilonorm_host,   double *wlolonorm_host );
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

   cout << "Verifying correctness for dimension and block size 256 ..."
        << endl;
   fail = verify_correctness(256,256,0);

   cout << "Nonblocked version form dimension 1024 and block size 256..."
        << endl;
   fail = verify_correctness(1024,256,0);

   cout << "Blocked version form dimension 4096 and block size 256..."
        << endl;
   fail = verify_correctness(4096,256,1);

   return 0;
}

void run
 ( int dim, int BS, int freq, int mode, int blocked,
   double *vhihinorm_device, double *vlohinorm_device,
   double *vhilonorm_device, double *vlolonorm_device,
   double *vhihinorm_host,   double *vlohinorm_host,
   double *vhilonorm_host,   double *vlolonorm_host,
   double *whihinorm_device, double *wlohinorm_device,
   double *whilonorm_device, double *wlolonorm_device,
   double *whihinorm_host,   double *wlohinorm_host,
   double *whilonorm_host,   double *wlolonorm_host )
{
   const int timevalue = time(NULL); // no fixed seed to verify correctness
   srand(timevalue);

   double* vrehihi_host = new double[dim]; // highest real parts on host
   double* vrelohi_host = new double[dim]; // 2nd highest real parts on host
   double* vrehilo_host = new double[dim]; // 2nd lowest real parts on host
   double* vrelolo_host = new double[dim]; // lowest real parts on host
   double* vimhihi_host = new double[dim]; // highest imaginary parts on host
   double* vimlohi_host = new double[dim]; // 2nd highest imag parts on host
   double* vimhilo_host = new double[dim]; // 2nd lowest imag parts on host
   double* vimlolo_host = new double[dim]; // lowest imaginary parts on host
   double* wrehihi_host = new double[dim]; // a copy for normalization
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
   double* vrelolo_device = new double[dim]; // lowest real parts on device
   double* vimhihi_device = new double[dim]; // highest imaginary parts
   double* vimlohi_device = new double[dim]; // 2nd highest imaginary parts
   double* vimhilo_device = new double[dim]; // 2nd lowest imaginary parts
   double* vimlolo_device = new double[dim]; // lowest imaginary parts

   random_complex4_vectors
      (dim,vrehihi_host,vrelohi_host,vrehilo_host,vrelolo_host,
           vimhihi_host,vimlohi_host,vimhilo_host,vimlolo_host,
       vrehihi_device,vrelohi_device,vrehilo_device,vrelolo_device,
       vimhihi_device,vimlohi_device,vimhilo_device,vimlolo_device);

   if(mode == 0 || mode == 2)
   {
      GPU_norm(vrehihi_device,vrelohi_device,vrehilo_device,vrelolo_device,
               vimhihi_device,vimlohi_device,vimhilo_device,vimlolo_device,
         dim,1,BS,
         vhihinorm_device,vlohinorm_device,vhilonorm_device,vlolonorm_device,
         blocked);
      GPU_norm(vrehihi_device,vrelohi_device,vrehilo_device,vrelolo_device,
               vimhihi_device,vimlohi_device,vimhilo_device,vimlolo_device,
         dim,freq,BS,
         whihinorm_device,wlohinorm_device,whilonorm_device,wlolonorm_device,
         blocked);
   }

   if(mode == 1 || mode == 2)
   {
      for(int i=0; i<=freq; i++)
      {
         CPU_norm(vrehihi_host,vrelohi_host,vrehilo_host,vrelolo_host,
                  vimhihi_host,vimlohi_host,vimhilo_host,vimlolo_host,dim,
                  vhihinorm_host,vlohinorm_host,vhilonorm_host,vlolonorm_host);
         make_copy(dim,vrehihi_host,vrelohi_host,vrehilo_host,vrelolo_host,
                       vimhihi_host,vimlohi_host,vimhilo_host,vimlolo_host,
                   wrehihi_host,wrelohi_host,wrehilo_host,wrelolo_host,
                   wimhihi_host,wimlohi_host,wimhilo_host,wimlolo_host);
         CPU_normalize(wrehihi_host,wrelohi_host,wrehilo_host,wrelolo_host,
                       wimhihi_host,wimlohi_host,wimhilo_host,wimlolo_host,dim,
                       *vhihinorm_host,*vlohinorm_host,
                       *vhilonorm_host,*vlolonorm_host);
         CPU_norm(wrehihi_host,wrelohi_host,wrehilo_host,wrelolo_host,
                  wimhihi_host,wimlohi_host,wimhilo_host,wimlolo_host,dim,
                  whihinorm_host,wlohinorm_host,whilonorm_host,wlolonorm_host);
      }
   }
}

int verify_correctness ( int dim, int BS, int blocked )
{
   double vhihinorm_device,vlohinorm_device,vhilonorm_device,vlolonorm_device;
   // norm before normalization on device
   double whihinorm_device,wlohinorm_device,whilonorm_device,wlolonorm_device;
   // norm after normalization on device
   double vhihinorm_host,vlohinorm_host,vhilonorm_host,vlolonorm_host;
   // norm before normalization on host
   double whihinorm_host,wlohinorm_host,whilonorm_host,wlolonorm_host;
   // norm after normalization on host

   run(dim,BS,1,2,blocked,
       &vhihinorm_device,&vlohinorm_device,&vhilonorm_device,&vlolonorm_device,
       &vhihinorm_host,  &vlohinorm_host,  &vhilonorm_host,  &vlolonorm_host,
       &whihinorm_device,&wlohinorm_device,&whilonorm_device,&wlolonorm_device,
       &whihinorm_host,  &wlohinorm_host,  &whilonorm_host,  &wlolonorm_host);

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
   cout << endl;

   const double tol = 1.0e-36;
   double err = abs(vhihinorm_device - vhihinorm_host)
              + abs(whihinorm_device - whihinorm_host)
              + abs(vlohinorm_device - vlohinorm_host)
              + abs(vlohinorm_device - vlohinorm_host)
              + abs(vlohinorm_device - vlohinorm_host)
              + abs(wlohinorm_device - wlohinorm_host)
              + abs(vlolonorm_device - vlolonorm_host)
              + abs(wlolonorm_device - wlolonorm_host);

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
