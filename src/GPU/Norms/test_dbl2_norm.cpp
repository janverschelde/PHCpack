/* This program runs automatic tests to compute norms of real vectors
   in double double precision, for preset values of the parameters.
   No input is required from the user. */

#include <ctime>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "dbl2_norm_host.h"
#include "dbl2_norm_kernels.h"
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
 *   Computes norms for random real vectors.
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
 * Computes norms of real vectors of dimension dim, block size BS,
 * and verifies the correctness.  Returns 1 if the error was too large,
 * returns 0 if the 2-norm and normalization was within the tolerance. */

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
   const int maxbs = dd_shmemsize;

   for(int d=512; d<=4096; d=d+512) // for(int d=32; d<=1024; d=d+32)
   {
      for(int bs=256; bs<=d; bs=bs+32)
      {
         if(d % bs == 0) // dimension must be a multiple of block size
         {
            cout << "Verifying for dimension " << d
                 << " and for block size " << bs << " ..." << endl;
            fail = verify_correctness(d,bs);
            if(fail != 0)
            {
               cout << "Failed for dimension " << d
                    << " and block size " << bs << "!" << endl;
               break;
            }
         }
         if(bs >= maxbs) break;
      }
      if(fail != 0) break;
   }
   const int dimmaxbs = maxbs;      // dimension equals maxbs
   const int freqmaxbs = dimmaxbs*32;

   cout << endl;
   cout << "Time on host for dimension " << dimmaxbs
        << " and frequency " << freqmaxbs << " ..." << endl;
   time_host(dimmaxbs,freqmaxbs);

   cout << endl;
   cout << "Time on device for dimension " << dimmaxbs
        << ", block size " << maxbs
        << ", and frequency " << freqmaxbs << " ..." << endl;
   time_device(dimmaxbs,freqmaxbs,maxbs);

   const int dim = 512*512; // largest dimension
   const int freq = 1024;
   const int bs = 512;

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

void run ( int dim, int BS, int freq, int mode,
           double *vhinorm_device, double *vlonorm_device,
           double *vhinorm_host, double *vlonorm_host,
           double *whinorm_device, double *wlonorm_device,
           double *whinorm_host, double *wlonorm_host )
{
   const int timevalue = time(NULL); // no fixed seed to verify correctness
   srand(timevalue);

   double* vhi_host = new double[dim];   // high parts on the host
   double* vlo_host = new double[dim];   // low parts on the host
   double* vhi_device = new double[dim]; // high parts on the device
   double* vlo_device = new double[dim]; // low parts on the device
   double* whi_host = new double[dim];   // high parts copy for normalization
   double* wlo_host = new double[dim];   // low parts copy for normalization

   random_double2_vectors(dim,vhi_host,vlo_host,vhi_device,vlo_device);

   if(mode == 0 || mode == 2)
   {
      GPU_norm(vhi_device,vlo_device,dim,1,BS,
               vhinorm_device,vlonorm_device,1);
      GPU_norm(vhi_device,vlo_device,dim,freq,BS,
               whinorm_device,wlonorm_device,1);
   }
   if(mode == 1 || mode == 2)
   {
      for(int i=0; i<=freq; i++)
      {
         CPU_norm(vhi_host,vlo_host,dim,vhinorm_host,vlonorm_host);
         make_copy(dim,vhi_host,vlo_host,whi_host,wlo_host);
         CPU_normalize(whi_host,wlo_host,dim,*vhinorm_host,*vlonorm_host);
         CPU_norm(whi_host,wlo_host,dim,whinorm_host,wlonorm_host);
      }
   }
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
   cout << "   GPU norm : " << endl;
   cout << "     hi : " << vhinorm_device << endl;
   cout << "     lo : " << vlonorm_device << endl;
   cout << "   GPU norm after normalization : " << endl;
   cout << "     hi : " << whinorm_device << endl;
   cout << "     lo : " << wlonorm_device << endl;

   vnrm_d[0] = vhinorm_device; vnrm_d[1] = vlonorm_device;
   wnrm_d[0] = whinorm_device; wnrm_d[1] = wlonorm_device;
   vnrm_h[0] = vhinorm_host; vnrm_h[1] = vlonorm_host;
   wnrm_h[0] = whinorm_host; wnrm_h[1] = wlonorm_host;

   cout << "   GPU norm : ";
   dd_write(vnrm_d,32); cout << endl;
   cout << "       after normalization : ";
   dd_write(wnrm_d,32); cout << endl;
   cout << "   CPU norm : ";
   dd_write(vnrm_h,32); cout << endl;
   cout << "       after normalization : ";
   dd_write(wnrm_h,32); cout << endl;

   const double tol = 1.0e-24;
   double err = abs(vhinorm_device - vhinorm_host)
              + abs(whinorm_device - whinorm_host)
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

void time_host ( int dim, int freq )
{
   const int BS = 32;  // block size, not used anyway

   double vhinorm_device,vlonorm_device; // norm before normalization on device
   double whinorm_device,wlonorm_device; // norm after normalization on device
   double vhinorm_host,vlonorm_host;     // norm before normalization on host
   double whinorm_host,wlonorm_host;     // norm after normalization on host

   long int begin_time = clock();
   run(dim,BS,freq,1,
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
   run(dim,BS,freq,0,
       &vhinorm_device,&vlonorm_device,&vhinorm_host,&vlonorm_host,
       &whinorm_device,&wlonorm_device,&whinorm_host,&wlonorm_host);
   long int end_time = clock();
   int elapsed = 1000*(end_time - begin_time)/CLOCKS_PER_SEC;
   cout << "elapsed time : " << elapsed << " milliseconds" << endl;
}
