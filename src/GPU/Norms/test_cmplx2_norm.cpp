/* This program runs automatic tests to compute norms of complex vectors
   in double double precision, for preset values of the parameters.
   No input is required from the user. */

#include <ctime>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
// #include "cmplx_norm_kernels.h"
#include "cmplx2_norm_host.h"
#include "random2_vectors.h"
#include "double_double.h"

using namespace std;

void run ( int dim, int BS, int freq, int mode,
           double *vhinorm_device, double *vlonorm_device,
           double *vhinorm_host, double *vlonorm_host,
           double *whinorm_device, double *wlonorm_device,
           double *whinorm_host, double *wlonorm_host );
/*
 * DESCRIPTION :
 *   Computes norms for random complex vectors,
 *   in double doble precision.
 *
 * ON ENTRY :
 *   dim   dimension of the random vector;
 *   BS    block size;
 *   freq  frequency of runs for the timings;
 *   mode  if 0, then only GPU computes,
 *         if 1, then only CPU computes,
 *         if 2, then both GPU and CPU  compute.
 *
 * ON RETURN :
 *   vhinorm_device  high part of 2-norm compute by the device;
 *   vlonorm_device  low part of 2-norm computed by the device;
 *   vhinorm_host    high part of 2-norm compute by the host;
 *   vlonorm_host    low part of 2-norm computed by the host;
 *   whinorm_device  high part of 2-norm on device after normalization;
 *   wlonorm_device  low norm of 2-norm on device after normalization;
 *   whinorm_host    high part of 2-norm on host after normalization;
 *   wlonorm_host    low norm of 2-norm on host after normalization.  */

int verify_correctness ( int dim, int BS );
/*
 * Computes norms of real vectors both on GPU and CPU
 * of dimension dim and verifies the correctness.
 * The number of threads in a block is given in the block size BS. */

int main ( void )
{
   int fail = verify_correctness(64,64);

   return 0;
}

void run ( int dim, int BS, int freq, int mode,
           double *vhinorm_device, double *vlonorm_device,
           double *vhinorm_host, double *vlonorm_host,
           double *whinorm_device, double *wlonorm_device,
           double *whinorm_host, double *wlonorm_host )
{
   const int timevalue = time(NULL); // no fixed seed to verify correctness
   srand(timevalue);

   double* vrehi_host = new double[dim]; // high real parts on the host
   double* vrelo_host = new double[dim]; // low real parts on the host
   double* vimhi_host = new double[dim]; // high imaginary parts on the host
   double* vimlo_host = new double[dim]; // low imaginary parts on the host
   double* wrehi_host = new double[dim]; // a copy for normalization
   double* wrelo_host = new double[dim];
   double* wimhi_host = new double[dim];
   double* wimlo_host = new double[dim];
   double* vrehi_device = new double[dim]; // high real parts on the device
   double* vrelo_device = new double[dim]; // low real parts on the device
   double* vimhi_device = new double[dim]; // high imaginary parts
   double* vimlo_device = new double[dim]; // low imaginary parts

   random_complex2_vectors
     (dim,vrehi_host,vrelo_host,vimhi_host,vimlo_host,
      vrehi_device,vrelo_device,vimhi_device,vimlo_device);
/*
   if(mode == 0 || mode == 2)
   {
      GPU_norm(vre_device,vim_device,dim,1,BS,vnorm_device,1);
      GPU_norm(vre_device,vim_device,dim,freq,BS,wnorm_device,1);
   }
   if(mode == 1 || mode == 2)
   {
      for(int i=0; i<=freq; i++)
      {
 */
         CPU_norm(vrehi_host,vrelo_host,vimhi_host,vimlo_host,dim,
                  vhinorm_host,vlonorm_host);
         make_copy(dim,vrehi_host,vrelo_host,vimhi_host,vimlo_host,
                   wrehi_host,wrelo_host,wimhi_host,wimlo_host);
         CPU_normalize(wrehi_host,wrelo_host,wimhi_host,wimlo_host,dim,
                       *vhinorm_host,*vlonorm_host);
         CPU_norm(wrehi_host,wrelo_host,wimhi_host,wimlo_host,dim,
                  whinorm_host,wlonorm_host);
/*
      }
   }
 */
}

int verify_correctness ( int dim, int BS )
{

   double vhinorm_device,vlonorm_device; // norm before normalization on device
   double whinorm_device,wlonorm_device; // norm after normalization on device
   double vhinorm_host,vlonorm_host;     // norm before normalization on host
   double whinorm_host,wlonorm_host;     // norm after normalization on host
   double vnrm_h[2],vnrm_d[2],wnrm_h[2],wnrm_d[2];

   run(dim,BS,1,2,
       &vhinorm_device,&vlonorm_device,&vhinorm_host,&vlonorm_host,
       &whinorm_device,&wlonorm_device,&whinorm_host,&wlonorm_host);

   cout << scientific << setprecision(16);

   cout << "   CPU norm : " << endl;
   cout << "     hi : " << vhinorm_host << endl;
   cout << "     lo : " << vlonorm_host << endl;
   cout << "   CPU norm after normalization : " << endl;
   cout << "     hi : " << whinorm_host << endl;
   cout << "     lo : " << wlonorm_host << endl;

   vnrm_h[0] = vhinorm_host; vnrm_h[1] = vlonorm_host;
   wnrm_h[0] = whinorm_host; wnrm_h[1] = wlonorm_host;

   cout << "   CPU norm : ";
   dd_write(vnrm_h,32); cout << endl;
   cout << "       after normalization : ";
   dd_write(wnrm_h,32); cout << endl;

   return 0;
}
