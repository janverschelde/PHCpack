/* This program runs automatic tests to compute norms of real vectors
   in octo double precision, for preset values of the parameters.
   No input is required from the user. */

#include <ctime>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "octo_double_functions.h"
#include "random8_vectors.h"
#include "dbl8_norm_host.h"
#include "dbl8_norm_kernels.h"

using namespace std;

void run
 ( int dim, int BS, int freq, int mode, int blocked,
   double *vhihihinorm_device, double *vlohihinorm_device,
   double *vhilohinorm_device, double *vlolohinorm_device,
   double *vhihilonorm_device, double *vlohilonorm_device,
   double *vhilolonorm_device, double *vlololonorm_device,
   double *vhihihinorm_host, double *vlohihinorm_host,
   double *vhilohinorm_host, double *vlolohinorm_host,
   double *vhihilonorm_host, double *vlohilonorm_host,
   double *vhilolonorm_host, double *vlololonorm_host,
   double *whihihinorm_device, double *wlohihinorm_device,
   double *whilohinorm_device, double *wlolohinorm_device,
   double *whihilonorm_device, double *wlohilonorm_device,
   double *whilolonorm_device, double *wlololonorm_device,
   double *whihihinorm_host, double *wlohihinorm_host,
   double *whilohinorm_host, double *wlolohinorm_host,
   double *whihilonorm_host, double *wlohilonorm_host,
   double *whilolonorm_host, double *wlololonorm_host );
/*
 * DESCRIPTION :
 *   Computes norms for random real vectors.
 *
 * ON ENTRY :
 *   dim      dimension of the random vector;
 *   BS       block size;
 *   freq     frequency of runs for the timings;
 *   mode     if 0, then only GPU computes,
 *            if 1, then only CPU computes,
 *            if 2, then both GPU and CPU compute;
 *   blocked  if 0, then the vector should be of medium size
 *            and only one block will compute,
 *            if 1, then as many as dim/BS blocks will compute.
 *
 * ON RETURN :
 *   vhihihinorm_device  highest part of 2-norm compute by the device;
 *   vlohihinorm_device  2nd highest part of 2-norm computed by the device;
 *   vhilohinorm_device  3rd highest part of 2-norm computed by the device;
 *   vlolohinorm_device  4th highest part of 2-norm computed by the device;
 *   vhihilonorm_device  4th lowest part of 2-norm computed by the device;
 *   vlohilonorm_device  3rd lowest part of 2-norm computed by the device;
 *   vhilolonorm_device  2nd lowest part of 2-norm computed by the device;
 *   vlololonorm_device  lowest part of 2-norm computed by the device;
 *   vhihihinorm_host    highest part of 2-norm compute by the host;
 *   vlohihinorm_host    2nd highest part of 2-norm computed by the host;
 *   vhilohinorm_host    3rd highest part of 2-norm computed by the host;
 *   vlolohinorm_host    4th highest part of 2-norm computed by the host;
 *   vhihilonorm_host    4th lowest part of 2-norm computed by the host;
 *   vlohilonorm_host    3rd lowest part of 2-norm computed by the host;
 *   vhilolonorm_host    2nd lowest part of 2-norm computed by the host;
 *   vlololonorm_host    lowest part of 2-norm computed by the host;
 *   whihihinorm_device  highest part of 2-norm on device after normalization;
 *   wlohihinorm_device  2nd highest of 2-norm on device after normalization;
 *   whilohinorm_device  3rd highest of 2-norm on device after normalization;
 *   wlolohinorm_device  4th highest of 2-norm on device after normalization;
 *   whihilonorm_device  4th lowest part of 2-norm on device;
 *   wlohilonorm_device  3rd lowest part of 2-norm on device;
 *   whilolonorm_device  2nd lowest part of 2-norm on device;
 *   wlololonorm_device  lowest part of 2-norm on device after normalization;
 *   whihihinorm_host    highest part of 2-norm on host after normalization;
 *   wlohihinorm_host    2nd highest part of 2-norm on host;
 *   whilohinorm_host    3rd highest part of 2-norm on host;
 *   wlolohinorm_host    4th highest part of 2-norm on host;
 *   whihilonorm_host    4th lowest part of 2-norm on host;
 *   wlohilonorm_host    3rd lowest part of 2-norm on host;
 *   whilolonorm_host    2nd lowest part of 2-norm on host;
 *   wlololonorm_host    lowest part of 2-norm on host after normalization. */

int verify_correctness ( int dim, int BS, int blocked );
/*
 * Computes norms of real vectors of dimension dim, block size BS,
 * and verifies the correctness.  Returns 1 if the error was too large,
 * returns 0 if the 2-norm and normalization was within the tolerance.
 * If blocked and dim is a multiple of BS, then many blocks compute. */

int main ( void )
{
   int fail;

   cout << "Verifying correctness for dimension and block size 128 ..."
        << endl;
   fail = verify_correctness(128,128,0);

   cout << "Nonblocked version for dimension 1024 and block size 128 ..."
        << endl;
   fail = verify_correctness(1024,128,0);

   cout << "Blocked version for dimension 4096 and block size 128 ..."
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
   double *vhihihinorm_host, double *vlohihinorm_host,
   double *vhilohinorm_host, double *vlolohinorm_host,
   double *vhihilonorm_host, double *vlohilonorm_host,
   double *vhilolonorm_host, double *vlololonorm_host,
   double *whihihinorm_device, double *wlohihinorm_device,
   double *whilohinorm_device, double *wlolohinorm_device,
   double *whihilonorm_device, double *wlohilonorm_device,
   double *whilolonorm_device, double *wlololonorm_device,
   double *whihihinorm_host, double *wlohihinorm_host,
   double *whilohinorm_host, double *wlolohinorm_host,
   double *whihilonorm_host, double *wlohilonorm_host,
   double *whilolonorm_host, double *wlololonorm_host )
{
   const int timevalue = time(NULL); // no fixed seed to verify correctness
   srand(timevalue);

   double* vhihihi_host = new double[dim];   // highest parts on the host
   double* vlohihi_host = new double[dim];   // 2nd highest parts on the host
   double* vhilohi_host = new double[dim];   // 3nd highest parts on the host
   double* vlolohi_host = new double[dim];   // 4th highest parts on the host
   double* vhihilo_host = new double[dim];   // 4th lowest parts on the host
   double* vlohilo_host = new double[dim];   // 3rd lowest parts on the host
   double* vhilolo_host = new double[dim];   // 2nd lowest parts on the host
   double* vlololo_host = new double[dim];   // lowest parts on the host
   double* vhihihi_device = new double[dim]; // highest parts on the device
   double* vlohihi_device = new double[dim]; // 2nd highest parts on the device
   double* vhilohi_device = new double[dim]; // 3rd highest parts on the device
   double* vlolohi_device = new double[dim]; // 4th highest parts on the device
   double* vhihilo_device = new double[dim]; // 4th lowest parts on the device
   double* vlohilo_device = new double[dim]; // 3nd lowest parts on the device
   double* vhilolo_device = new double[dim]; // 2nd lowest parts on the device
   double* vlololo_device = new double[dim]; // lowest parts on the device
   double* whihihi_host = new double[dim];   // highest parts copy
   double* wlohihi_host = new double[dim];   // second highest parts copy
   double* whilohi_host = new double[dim];   // third highest parts copy
   double* wlolohi_host = new double[dim];   // fourth highest parts copy
   double* whihilo_host = new double[dim];   // fourth lowest parts copy
   double* wlohilo_host = new double[dim];   // third lowest parts copy
   double* whilolo_host = new double[dim];   // second lowest parts copy
   double* wlololo_host = new double[dim];   // lowest parts copy

   random_double8_vectors
      (dim,vhihihi_host,vlohihi_host,vhilohi_host,vlolohi_host,
           vhihilo_host,vlohilo_host,vhilolo_host,vlololo_host,
       vhihihi_device,vlohihi_device,vhilohi_device,vlolohi_device,
       vhihilo_device,vlohilo_device,vhilolo_device,vlololo_device);

   if(mode == 0 || mode == 2)
   {
      GPU_norm(vhihihi_device,vlohihi_device,vhilohi_device,vlolohi_device,
               vhihilo_device,vlohilo_device,vhilolo_device,vlololo_device,
               dim,1,BS,vhihihinorm_device,vlohihinorm_device,
               vhilohinorm_device,vlolohinorm_device,
               vhihilonorm_device,vlohilonorm_device,
               vhilolonorm_device,vlololonorm_device,blocked);
      GPU_norm(vhihihi_device,vlohihi_device,vhilohi_device,vlolohi_device,
               vhihilo_device,vlohilo_device,vhilolo_device,vlololo_device,
               dim,freq,BS,whihihinorm_device,wlohihinorm_device,
               whilohinorm_device,wlolohinorm_device,
               whihilonorm_device,wlohilonorm_device,
               whilolonorm_device,wlololonorm_device,blocked);
   }
   if(mode == 1 || mode == 2)
   {
      for(int i=0; i<=freq; i++)
      {
         CPU_norm(vhihihi_host,vlohihi_host,vhilohi_host,vlolohi_host,
                  vhihilo_host,vlohilo_host,vhilolo_host,vlololo_host,dim,
                  vhihihinorm_host,vlohihinorm_host,
                  vhilohinorm_host,vlolohinorm_host,
                  vhihilonorm_host,vlohilonorm_host,
                  vhilolonorm_host,vlololonorm_host);
         make_copy(dim,vhihihi_host,vlohihi_host,vhilohi_host,vlolohi_host,
                       vhihilo_host,vlohilo_host,vhilolo_host,vlololo_host,
                       whihihi_host,wlohihi_host,whilohi_host,wlolohi_host,
                       whihilo_host,wlohilo_host,whilolo_host,wlololo_host);
         CPU_normalize(whihihi_host,wlohihi_host,whilohi_host,wlolohi_host,
                       whihilo_host,wlohilo_host,whilolo_host,wlololo_host,
                       dim,*vhihihinorm_host,*vlohihinorm_host,
                       *vhilohinorm_host,*vlolohinorm_host,
                       *vhihilonorm_host,*vlohilonorm_host,
                       *vhilolonorm_host,*vlololonorm_host);
         CPU_norm(whihihi_host,wlohihi_host,whilohi_host,wlolohi_host,
                  whihilo_host,wlohilo_host,whilolo_host,wlololo_host,dim,
                  whihihinorm_host,wlohihinorm_host,
                  whilohinorm_host,wlolohinorm_host,
                  whihilonorm_host,wlohilonorm_host,
                  whilolonorm_host,wlololonorm_host);
      }
   }
}

int verify_correctness ( int dim, int BS, int blocked )
{
   double vhihihinorm_device,vlohihinorm_device;
   double vhilohinorm_device,vlolohinorm_device;
   double vhihilonorm_device,vlohilonorm_device;
   double vhilolonorm_device,vlololonorm_device;
   // norm before normalization on device
   double whihihinorm_device,wlohihinorm_device;
   double whilohinorm_device,wlolohinorm_device;
   double whihilonorm_device,wlohilonorm_device;
   double whilolonorm_device,wlololonorm_device;
   // norm after normalization on device
   double vhihihinorm_host,vlohihinorm_host;
   double vhilohinorm_host,vlolohinorm_host;
   double vhihilonorm_host,vlohilonorm_host;
   double vhilolonorm_host,vlololonorm_host;
   // norm before normalization on host
   double whihihinorm_host,wlohihinorm_host;
   double whilohinorm_host,wlolohinorm_host;
   double whihilonorm_host,wlohilonorm_host;
   double whilolonorm_host,wlololonorm_host;
   // norm after normalization on host

   run(dim,BS,1,2,blocked,
       &vhihihinorm_device,&vlohihinorm_device,
       &vhilohinorm_device,&vlolohinorm_device,
       &vhihilonorm_device,&vlohilonorm_device,
       &vhilolonorm_device,&vlololonorm_device,
       &vhihihinorm_host,&vlohihinorm_host,
       &vhilohinorm_host,&vlolohinorm_host,
       &vhihilonorm_host,&vlohilonorm_host,
       &vhilolonorm_host,&vlololonorm_host,
       &whihihinorm_device,&wlohihinorm_device,
       &whilohinorm_device,&wlolohinorm_device,
       &whihilonorm_device,&wlohilonorm_device,
       &whilolonorm_device,&wlololonorm_device,
       &whihihinorm_host,&wlohihinorm_host,
       &whilohinorm_host,&wlolohinorm_host,
       &whihilonorm_host,&wlohihinorm_host,
       &whilolonorm_host,&wlolohinorm_host);

   cout << scientific << setprecision(16);

   cout << "CPU norm : " << endl;
   odf_write_doubles(vhihihinorm_host,vlohihinorm_host,
                     vhilohinorm_host,vlolohinorm_host,
                     vhihilonorm_host,vlohilonorm_host,
                     vhilolonorm_host,vlololonorm_host);
   cout << "CPU norm after normalization : " << endl;
   odf_write_doubles(whihihinorm_host,wlohihinorm_host,
                     whilohinorm_host,wlolohinorm_host,
                     whihilonorm_host,wlohilonorm_host,
                     whilolonorm_host,wlololonorm_host);

   cout << "GPU norm : " << endl;
   odf_write_doubles(vhihihinorm_device,vlohihinorm_device,
                     vhilohinorm_device,vlolohinorm_device,
                     vhihilonorm_device,vlohilonorm_device,
                     vhilolonorm_device,vlololonorm_device);
   cout << "GPU norm after normalization : " << endl;
   odf_write_doubles(whihihinorm_device,wlohihinorm_device,
                     whilohinorm_device,wlolohinorm_device,
                     whihilonorm_device,wlohilonorm_device,
                     whilolonorm_device,wlololonorm_device);

   const double tol = 1.0e-96;
   double err = fabs(vhihihinorm_device - vhihihinorm_host)
              + fabs(whihihinorm_device - whihihinorm_host)
              + fabs(vlohihinorm_device - vlohihinorm_host)
              + fabs(wlohihinorm_device - wlohihinorm_host)
              + fabs(vhilohinorm_device - vhilohinorm_host)
              + fabs(whilohinorm_device - whilohinorm_host)
              + fabs(vlolohinorm_device - vlolohinorm_host)
              + fabs(wlolohinorm_device - wlolohinorm_host)
              + fabs(vhihilonorm_device - vhihilonorm_host)
              + fabs(whihilonorm_device - whihilonorm_host)
              + fabs(vlohilonorm_device - vlohilonorm_host)
              + fabs(wlohilonorm_device - wlohilonorm_host)
              + fabs(vhilolonorm_device - vhilolonorm_host)
              + fabs(whilolonorm_device - whilolonorm_host)
              + fabs(vlololonorm_device - vlololonorm_host)
              + fabs(wlololonorm_device - wlololonorm_host);

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
