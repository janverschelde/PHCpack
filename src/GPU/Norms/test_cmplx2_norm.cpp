/* This program runs automatic tests to compute norms of complex vectors
   in double double precision, for preset values of the parameters.
   No input is required from the user. */

#include <ctime>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "cmplx2_norm_kernels.h"
#include "cmplx2_norm_host.h"
#include "random2_vectors.h"
#include "double_double.h"

using namespace std;

void run ( int dim, int BS, int freq, int mode, int blocked,
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
 *   vlonorm_device  low part of 2-norm computed by the device;
 *   vhinorm_host    high part of 2-norm compute by the host;
 *   vlonorm_host    low part of 2-norm computed by the host;
 *   whinorm_device  high part of 2-norm on device after normalization;
 *   wlonorm_device  low part of 2-norm on device after normalization;
 *   whinorm_host    high part of 2-norm on host after normalization;
 *   wlonorm_host    low part of 2-norm on host after normalization.  */

int verify_correctness ( int dim, int BS, int blocked );
/*
 * Computes norms of real vectors both on GPU and CPU
 * of dimension dim and verifies the correctness.
 * The number of threads in a block is given in the block size BS.
 * If blocked and dim is a multiple of BS, then many blocks compute. */

void time_host ( int dim, int freq );
/*
 * Computes norms on the host for random vectors of dimension dim
 * as many times as the frequency freq and shows the elapsed time. */

void time_device ( int dim, int freq, int bs );
/*
 * Computes norms on the device for random vectors of dimension dim
 * as many times as the frequency freq, for block size bs,
 * and shows the elapsed time. */

int main ( void )
{
   int fail;

   cout << "Verifying correctness for dimension and block size 64 ..."
        << endl;
   fail = verify_correctness(64,64,0);

   cout << "nonblocked version ..." << endl;
   fail = verify_correctness(256,64,0);

   cout << "blocked version ..." << endl;
   fail = verify_correctness(256,64,1);

   const int dim = 256*256; // largest dimension
   const int freq = 1024;
   const int bs = 256;

   fail = verify_correctness(dim,bs,1);

   cout << endl;
   cout << "Time on host for dimension " << dim
        << " and frequency " << freq << " ..." << endl;
   time_host(dim,freq);

   cout << "Time on device for dimension " << dim
        << ", block size " << bs
        << ", and frequency " << freq << " ..." << endl;
   time_device(dim,freq,bs);

   return 0;
}

void run ( int dim, int BS, int freq, int mode, int blocked,
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

   if(mode == 0 || mode == 2)
   {
      GPU_norm(vrehi_device,vrelo_device,vimhi_device,vimlo_device,
               dim,1,BS,vhinorm_device,vlonorm_device,blocked);
      GPU_norm(vrehi_device,vrelo_device,vimhi_device,vimlo_device,
               dim,freq,BS,whinorm_device,wlonorm_device,blocked);
   }
   if(mode == 1 || mode == 2)
   {
      for(int i=0; i<=freq; i++)
      {
         CPU_norm(vrehi_host,vrelo_host,vimhi_host,vimlo_host,dim,
                  vhinorm_host,vlonorm_host);
         make_copy(dim,vrehi_host,vrelo_host,vimhi_host,vimlo_host,
                   wrehi_host,wrelo_host,wimhi_host,wimlo_host);
         CPU_normalize(wrehi_host,wrelo_host,wimhi_host,wimlo_host,dim,
                       *vhinorm_host,*vlonorm_host);
         CPU_norm(wrehi_host,wrelo_host,wimhi_host,wimlo_host,dim,
                  whinorm_host,wlonorm_host);
      }
   }
}

int verify_correctness ( int dim, int BS, int blocked )
{
   double vhinorm_device,vlonorm_device; // norm before normalization on device
   double whinorm_device,wlonorm_device; // norm after normalization on device
   double vhinorm_host,vlonorm_host;     // norm before normalization on host
   double whinorm_host,wlonorm_host;     // norm after normalization on host
   double vnrm_h[2],vnrm_d[2],wnrm_h[2],wnrm_d[2];

   run(dim,BS,1,2,blocked,
       &vhinorm_device,&vlonorm_device,&vhinorm_host,&vlonorm_host,
       &whinorm_device,&wlonorm_device,&whinorm_host,&wlonorm_host);

   cout << scientific << setprecision(16);

   cout << "   CPU norm : " << endl;
   cout << "     hi : " << vhinorm_host << endl;
   cout << "     lo : " << vlonorm_host << endl;
   cout << "   GPU norm : " << endl;
   cout << "     hi : " << vhinorm_device << endl;
   cout << "     lo : " << vlonorm_device << endl;
   cout << "   CPU norm after normalization : " << endl;
   cout << "     hi : " << whinorm_host << endl;
   cout << "     lo : " << wlonorm_host << endl;
   cout << "   GPU norm after normalization : " << endl;
   cout << "     hi : " << whinorm_device << endl;
   cout << "     lo : " << wlonorm_device << endl;

   vnrm_h[0] = vhinorm_host; vnrm_h[1] = vlonorm_host;
   wnrm_h[0] = whinorm_host; wnrm_h[1] = wlonorm_host;
   vnrm_d[0] = vhinorm_device; vnrm_d[1] = vlonorm_device;
   wnrm_d[0] = whinorm_device; wnrm_d[1] = wlonorm_device;

   cout << "   CPU norm : ";
   dd_write(vnrm_h,32); cout << endl;
   cout << "   GPU norm : ";
   dd_write(vnrm_d,32); cout << endl;
   cout << "   CPU norm after normalization : ";
   dd_write(wnrm_h,32); cout << endl;
   cout << "   GPU norm after normalization : ";
   dd_write(wnrm_d,32); cout << endl;

   const double tol = 1.0e-24;
   double err = abs(vlonorm_device - vlonorm_host)
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

void time_host ( int dim, int freq )
{
   const int BS = 32;  // block size, not used anyway

   double vhinorm_device,vlonorm_device; // norm before normalization on device
   double whinorm_device,wlonorm_device; // norm after normalization on device
   double vhinorm_host,vlonorm_host;     // norm before normalization on host
   double whinorm_host,wlonorm_host;     // norm after normalization on host

   long int begin_time = clock();
   run(dim,BS,freq,1,1,
       &vhinorm_device,&vlonorm_device,&vhinorm_host,&vlonorm_host,
       &whinorm_device,&wlonorm_device,&whinorm_host,&wlonorm_host);
   long int end_time = clock();
   int elapsed = 1000*(end_time - begin_time)/CLOCKS_PER_SEC;
   cout << "elapsed time : " << elapsed << " milliseconds" << endl;
}

void time_device ( int dim, int freq, int BS )
{
   double vhinorm_device,vlonorm_device; // norm before normalization on device
   double whinorm_device,wlonorm_device; // norm after normalization on device
   double vhinorm_host,vlonorm_host;     // norm before normalization on host
   double whinorm_host,wlonorm_host;     // norm after normalization on host

   long int begin_time = clock();
   run(dim,BS,freq,0,1,
       &vhinorm_device,&vlonorm_device,&vhinorm_host,&vlonorm_host,
       &whinorm_device,&wlonorm_device,&whinorm_host,&wlonorm_host);
   long int end_time = clock();
   int elapsed = 1000*(end_time - begin_time)/CLOCKS_PER_SEC;
   cout << "elapsed time : " << elapsed << " milliseconds" << endl;
}
