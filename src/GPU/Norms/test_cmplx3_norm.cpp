/* This program runs automatic tests to compute norms of complex vectors
   in triple double precision, for preset values of the parameters.
   No input is required from the user. */

#include <ctime>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "cmplx3_norm_kernels.h"
#include "cmplx3_norm_host.h"
#include "random3_vectors.h"
#include "triple_double_functions.h"

using namespace std;

void run
 ( int dim, int BS, int freq, int mode, int blocked,
   double *vhinorm_device, double *vminorm_device, double *vlonorm_device,
   double *vhinorm_host,   double *vminorm_host,   double *vlonorm_host,
   double *whinorm_device, double *wminorm_device, double *wlonorm_device,
   double *whinorm_host,   double *wminorm_host,   double *wlonorm_host );
/*
 * DESCRIPTION :
 *   Computes norms for random complex vectors,
 *   in double doble precision.
 *
 * ON ENTRY :
 *   dim      dimension of the random vector;
 *   BS       block size;
 *   freq     frequency of runs for the timings;
 *   mode     if 0, then only GPU computes,
 *            if 1, then only CPU computes,
 *            if 2, then both GPU and CPU  compute;
 *   blocked  if 0, then the vector should be of medium size
 *            and only one block will compute,
 *            if 1, then as many as dim/BS blocks will compute.
 *
 * ON RETURN :
 *   vhinorm_device  high part of 2-norm compute by the device;
 *   vminorm_device  middle part of 2-norm compute by the device;
 *   vlonorm_device  low part of 2-norm computed by the device;
 *   vhinorm_host    high part of 2-norm compute by the host;
 *   vminorm_host    middle part of 2-norm computed by the host;
 *   vlonorm_host    low part of 2-norm computed by the host;
 *   whinorm_device  high part of 2-norm on device after normalization;
 *   wminorm_device  middle part of 2-norm on device after normalization;
 *   wlonorm_device  low part of 2-norm on device after normalization;
 *   whinorm_host    high part of 2-norm on host after normalization;
 *   wminorm_host    middle part of 2-norm on host after normalization;
 *   wlonorm_host    low part of 2-norm on host after normalization.  */

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

   cout << "nonblocked version form dimension 1024 and block size 256..."
        << endl;
   fail = verify_correctness(1024,256,0);

   cout << "blocked version ..." << endl;
   fail = verify_correctness(4096,256,1);

   return 0;
}

void run
 ( int dim, int BS, int freq, int mode, int blocked,
   double *vhinorm_device, double *vminorm_device, double *vlonorm_device,
   double *vhinorm_host,   double *vminorm_host,   double *vlonorm_host,
   double *whinorm_device, double *wminorm_device, double *wlonorm_device,
   double *whinorm_host,   double *wminorm_host,   double *wlonorm_host )
{
   const int timevalue = time(NULL); // no fixed seed to verify correctness
   srand(timevalue);

   double* vrehi_host = new double[dim]; // high real parts on the host
   double* vremi_host = new double[dim]; // middle real parts on the host
   double* vrelo_host = new double[dim]; // low real parts on the host
   double* vimhi_host = new double[dim]; // high imaginary parts on the host
   double* vimmi_host = new double[dim]; // middle imaginary parts on the host
   double* vimlo_host = new double[dim]; // low imaginary parts on the host
   double* wrehi_host = new double[dim]; // a copy for normalization
   double* wremi_host = new double[dim];
   double* wrelo_host = new double[dim];
   double* wimhi_host = new double[dim];
   double* wimmi_host = new double[dim];
   double* wimlo_host = new double[dim];
   double* vrehi_device = new double[dim]; // high real parts on the device
   double* vremi_device = new double[dim]; // middle real parts on the device
   double* vrelo_device = new double[dim]; // low real parts on the device
   double* vimhi_device = new double[dim]; // high imaginary parts
   double* vimmi_device = new double[dim]; // middle imaginary parts
   double* vimlo_device = new double[dim]; // low imaginary parts

   random_complex3_vectors
      (dim,vrehi_host,vremi_host,vrelo_host,
           vimhi_host,vimmi_host,vimlo_host,
       vrehi_device,vremi_device,vrelo_device,
       vimhi_device,vimmi_device,vimlo_device);

   if(mode == 0 || mode == 2)
   {
      GPU_norm(vrehi_device,vremi_device,vrelo_device,
               vimhi_device,vimmi_device,vimlo_device,dim,1,BS,
               vhinorm_device,vminorm_device,vlonorm_device,blocked);
      GPU_norm(vrehi_device,vremi_device,vrelo_device,
               vimhi_device,vimmi_device,vimlo_device,dim,freq,BS,
               whinorm_device,wminorm_device,wlonorm_device,blocked);
   }

   if(mode == 1 || mode == 2)
   {
      for(int i=0; i<=freq; i++)
      {
         CPU_norm(vrehi_host,vremi_host,vrelo_host,
                  vimhi_host,vimmi_host,vimlo_host,dim,
                  vhinorm_host,vminorm_host,vlonorm_host);
         make_copy(dim,vrehi_host,vremi_host,vrelo_host,
                       vimhi_host,vimmi_host,vimlo_host,
                   wrehi_host,wremi_host,wrelo_host,
                   wimhi_host,wimmi_host,wimlo_host);
         CPU_normalize(wrehi_host,wremi_host,wrelo_host,
                       wimhi_host,wimmi_host,wimlo_host,dim,
                       *vhinorm_host,*vminorm_host,*vlonorm_host);
         CPU_norm(wrehi_host,wremi_host,wrelo_host,
                  wimhi_host,wimmi_host,wimlo_host,dim,
                  whinorm_host,wminorm_host,wlonorm_host);
      }
   }
}

int verify_correctness ( int dim, int BS, int blocked )
{
   double vhinorm_device,vminorm_device,vlonorm_device;
   // norm before normalization on device
   double whinorm_device,wminorm_device,wlonorm_device;
   // norm after normalization on device
   double vhinorm_host,vminorm_host,vlonorm_host;
   // norm before normalization on host
   double whinorm_host,wminorm_host,wlonorm_host;
   // norm after normalization on host

   run(dim,BS,1,2,blocked,
       &vhinorm_device,&vminorm_device,&vlonorm_device,
       &vhinorm_host,  &vminorm_host,  &vlonorm_host,
       &whinorm_device,&wminorm_device,&wlonorm_device,
       &whinorm_host,  &wminorm_host,  &wlonorm_host);

   cout << scientific << setprecision(16);

   cout << "CPU norm : " << endl;
   tdf_write_doubles(vhinorm_host,vminorm_host,vlonorm_host);
   cout << endl;
   cout << "GPU norm : " << endl;
   tdf_write_doubles(vhinorm_device,vminorm_device,vlonorm_device);
   cout << endl;
   cout << "CPU norm after normalization : " << endl;
   tdf_write_doubles(whinorm_host,wminorm_host,wlonorm_host);
   cout << endl;
   cout << "GPU norm after normalization : " << endl;
   tdf_write_doubles(whinorm_device,wminorm_device,wlonorm_device);
   cout << endl;

   const double tol = 1.0e-24;
   double err = abs(vhinorm_device - vhinorm_host)
              + abs(whinorm_device - whinorm_host)
              + abs(vminorm_device - vminorm_host)
              + abs(wminorm_device - wminorm_host)
              + abs(vlonorm_device - vlonorm_host)
              + abs(wlonorm_device - wlonorm_host);

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
