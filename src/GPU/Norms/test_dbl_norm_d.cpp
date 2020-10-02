/* This program runs automatic tests to compute norms of real vectors
   in double precision, for preset values of the parameters.
   No input is required from the user. */

#include <ctime>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "dbl_norm_kernels.h"
#include "dbl_norm_host.h"
#include "random_vectors.h"

using namespace std;

void run ( int dim, int BS, int freq, int mode,
           double *vnorm_device, double *vnorm_host,
           double *wnorm_device, double *wnorm_host );
/*
  DESCRIPTION :
    Computes norms for random real vectors.

  ON ENTRY :
    dim    dimension of the random vector;
    BS     block size, number of threads in a block;
    freq   frequency of the runs;
    mode   if 0, then only GPU computes,
           if 1, then only CPU computes,
           if 2, then both GPU and CPU  compute.

  ON RETURN :
    vnorm_device   norm computed by the device, if mode is 0 or 2;
    vnorm_host     norm computed by the host, if mode is 1 or 2;
    wnorm_device   norm on device after normalization, if mode is 0 or 2;
    wnorm_host     norm on host after normalization, if mode is 1 or 2.  */

int verify_correctness ( int dim, int BS );
/*
 * Computes norms of real vectors both on GPU and CPU
 * of dimension dim and verifies the correctness.
 * The number of threads in a block is given in the block size BS. */

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
   const int maxbs = d_shmemsize;

   for(int d=512; d<=4096; d=d+512) // for(int d=32; d<=1024; d=d+32)
   {
      for(int bs=512; bs<=d; bs=bs+32)
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

   const int dim = 1024*1024;
   const int freq = 4096;

   cout << endl;
   cout << "Time on host for dimension " << dim
        << " and frequency " << freq << " ..." << endl;
   time_host(dim,freq);

   cout << endl;
   for(int bs=128; bs<=maxbs; bs=bs+128)
   {
      if(dim % bs == 0) // dimension must be a multiple of BS
      {
         cout << "Time on device for dimension " << dim
               << ", block size " << bs
               << ", and frequency " << freq << " ..." << endl;
          time_device(dim,freq,bs);
      }
   }
   return 0;
}

void run ( int dim, int BS, int freq, int mode,
           double *vnorm_device, double *vnorm_host,
           double *wnorm_device, double *wnorm_host )
{
   const int timevalue = time(NULL); // no fixed seed to verify correctness
   srand(timevalue);

   double* v_host = new double[dim];   // vector on the host
   double* w_host = new double[dim];   // a copy for normalization
   double* v_device = new double[dim]; // vector on the device

   random_double_vectors(dim,v_host,v_device);

   if(mode == 0 || mode == 2)
   {
      GPU_norm(v_device,dim,1,BS,vnorm_device,1);
      GPU_norm(v_device,dim,freq,BS,wnorm_device,1);
   }
   if(mode == 1 || mode == 2)
   {
      for(int i=0; i<=freq; i++)
      {
         CPU_norm(v_host,dim,vnorm_host);
         make_copy(dim,v_host,w_host);
         CPU_normalize(w_host,dim,*vnorm_host);
         CPU_norm(w_host,dim,wnorm_host);
      }
   }
}

int verify_correctness ( int dim, int BS )
{
   const int freq = 1; // frequency, run only once
   const double tol = 1.0e-12; // tolerance on error

   double vnorm_device,vnorm_host; // norm before normalization
   double wnorm_device,wnorm_host; // norm after normalization

   run(dim,BS,freq,2,&vnorm_device,&vnorm_host,&wnorm_device,&wnorm_host);

   cout << scientific << setprecision(16);
   cout << "   GPU norm : " << vnorm_device << endl;
   cout << "       after normalization : " << wnorm_device << endl;
   cout << "   CPU norm : " << vnorm_host << endl;
   cout << "       after normalization : " << wnorm_host << endl;

   double err = abs(wnorm_device - wnorm_host);

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
   const int BS = 32;  // block size, number of threads in a a block

   double vnorm_device,vnorm_host; // norm before normalization
   double wnorm_device,wnorm_host; // norm after normalization

   long int begin_time = clock();
   run(dim,BS,freq,1,&vnorm_device,&vnorm_host,&wnorm_device,&wnorm_host);
   long int end_time = clock();
   int elapsed = 1000*(end_time - begin_time)/CLOCKS_PER_SEC;
   cout << "elapsed time : " << elapsed << " milliseconds" << endl;
}

void time_device ( int dim, int freq, int BS )
{
   // const int BS = 512;  // block size, number of threads in a a block

   double vnorm_device,vnorm_host; // norm before normalization
   double wnorm_device,wnorm_host; // norm after normalization

   long int begin_time = clock();
   run(dim,BS,freq,0,&vnorm_device,&vnorm_host,&wnorm_device,&wnorm_host);
   long int end_time = clock();
   int elapsed = 1000*(end_time - begin_time)/CLOCKS_PER_SEC;
   cout << "elapsed time : " << elapsed << " milliseconds" << endl;
}
