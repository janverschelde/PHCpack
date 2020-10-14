/* This program runs automatic tests to compute norms of real vectors
   in quad double precision, for preset values of the parameters.
   No input is required from the user. */

#include <ctime>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "quad_double_functions.h"
#include "random4_vectors.h"
#include "dbl4_norm_host.h"
#include "dbl4_norm_kernels.h"

using namespace std;

void run
 ( int dim, int BS, int freq, int mode, int blocked,
   double *vhihinorm_device, double *vlohinorm_device,
   double *vhilonorm_device, double *vlolonorm_device,
   double *vhihinorm_host, double *vlohinorm_host,
   double *vhilonorm_host, double *vlolonorm_host,
   double *whihinorm_device, double *wlohiinorm_device,
   double *whilonorm_device, double *wlolonorm_device,
   double *whihinorm_host, double *wlohinorm_host,
   double *whilonorm_host, double *wlolonorm_host );
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
 *   vhihinorm_device  highest part of 2-norm compute by the device;
 *   vlohinorm_device  2nd highest part of 2-norm computed by the device;
 *   vhilonorm_device  2nd lowest part of 2-norm computed by the device;
 *   vlolonorm_device  lowest part of 2-norm computed by the device;
 *   vhihinorm_host    highest part of 2-norm compute by the host;
 *   vlohinorm_host    2nd highest part of 2-norm computed by the host;
 *   vhilonorm_host    2nd lowest part of 2-norm computed by the host;
 *   vlolonorm_host    lowest part of 2-norm computed by the host;
 *   whihinorm_device  highest part of 2-norm on device after normalization;
 *   wlominorm_device  2nd highest of 2-norm on device after normalization;
 *   whilonorm_device  2nd lowest part of 2-norm on device after normalization;
 *   wlolonorm_device  lowest part of 2-norm on device after normalization;
 *   whihinorm_host    highest part of 2-norm on host after normalization;
 *   wlohinorm_host    2nd highest part of 2-norm on host after normalization;
 *   whilonorm_host    2nd lowest part of 2-norm on host after normalization;
 *   wlolonorm_host    lowest part of 2-norm on host after normalization. */

int verify_correctness ( int dim, int BS, int blocked );
/*
 * Computes norms of real vectors of dimension dim, block size BS,
 * and verifies the correctness.  Returns 1 if the error was too large,
 * returns 0 if the 2-norm and normalization was within the tolerance.
 * If blocked and dim is a multiple of BS, then many blocks compute. */

int main ( void )
{
   int fail;

   cout << "Verifying correctness for dimension and block size 256 ..."
        << endl;
   fail = verify_correctness(256,256,0);

   cout << "Nonblocked version for dimension 1024 and block size 256 ..."
        << endl;
   fail = verify_correctness(1024,256,0);

   cout << "Blocked version for dimension 4096 and block size 256 ..."
        << endl;
   fail = verify_correctness(4096,256,1);

   return 0;
}

void run 
 ( int dim, int BS, int freq, int mode, int blocked,
   double *vhihinorm_device, double *vlohinorm_device,
   double *vhilonorm_device, double *vlolonorm_device,
   double *vhihinorm_host, double *vlohinorm_host,
   double *vhilonorm_host, double *vlolonorm_host,
   double *whihinorm_device, double *wlohinorm_device,
   double *whilonorm_device, double *wlolonorm_device,
   double *whihinorm_host, double *wlohinorm_host,
   double *whilonorm_host, double *wlolonorm_host )
{
   const int timevalue = time(NULL); // no fixed seed to verify correctness
   srand(timevalue);

   double* vhihi_host = new double[dim];   // highest parts on the host
   double* vlohi_host = new double[dim];   // 2nd highest parts on the host
   double* vhilo_host = new double[dim];   // 2nd lowest parts on the host
   double* vlolo_host = new double[dim];   // lowest parts on the host
   double* vhihi_device = new double[dim]; // highest parts on the device
   double* vlohi_device = new double[dim]; // 2nd highest parts on the device
   double* vhilo_device = new double[dim]; // 2nd lowest parts on the device
   double* vlolo_device = new double[dim]; // lowest parts on the device
   double* whihi_host = new double[dim];   // highest parts copy
   double* wlohi_host = new double[dim];   // second highest parts copy
   double* whilo_host = new double[dim];   // second lowest parts copy
   double* wlolo_host = new double[dim];   // lowest parts copy

   random_double4_vectors
      (dim,vhihi_host,vlohi_host,vhilo_host,vlolo_host,
       vhihi_device,vlohi_device,vhilo_device,vlolo_device);

   if(mode == 0 || mode == 2)
   {
      GPU_norm(vhihi_device,vlohi_device,vhilo_device,vlolo_device,dim,1,BS,
               vhihinorm_device,vlohinorm_device,
               vhilonorm_device,vlolonorm_device,blocked);
      GPU_norm(vhihi_device,vlohi_device,vhilo_device,vlolo_device,
               dim,freq,BS,whihinorm_device,wlohinorm_device,
               whilonorm_device,wlolonorm_device,blocked);
   }
   if(mode == 1 || mode == 2)
   {
      for(int i=0; i<=freq; i++)
      {
         CPU_norm(vhihi_host,vlohi_host,vhilo_host,vlolo_host,dim,
                  vhihinorm_host,vlohinorm_host,
                  vhilonorm_host,vlolonorm_host);
         make_copy(dim,vhihi_host,vlohi_host,vhilo_host,vlolo_host,
                       whihi_host,wlohi_host,whilo_host,wlolo_host);
         CPU_normalize(whihi_host,wlohi_host,whilo_host,wlolo_host,dim,
                       *vhihinorm_host,*vlohinorm_host,
                       *vhilonorm_host,*vlolonorm_host);
         CPU_norm(whihi_host,wlohi_host,whilo_host,wlolo_host,dim,
                  whihinorm_host,wlohinorm_host,
                  whilonorm_host,wlolonorm_host);
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
       &vhihinorm_host,&vlohinorm_host,&vhilonorm_host,&vlolonorm_host,
       &whihinorm_device,&wlohinorm_device,&whilonorm_device,&wlolonorm_device,
       &whihinorm_host,&wlohinorm_host,&whilonorm_host,&wlolonorm_host);

   cout << scientific << setprecision(16);

   cout << "CPU norm : " << endl;
   qdf_write_doubles(vhihinorm_host,vlohinorm_host,
                     vhilonorm_host,vlolonorm_host);
   cout << "CPU norm after normalization : " << endl;
   qdf_write_doubles(whihinorm_host,wlohinorm_host,
                     whilonorm_host,wlolonorm_host);

   cout << "GPU norm : " << endl;
   qdf_write_doubles(vhihinorm_device,vlohinorm_device,
                     vhilonorm_device,vlolonorm_device);
   cout << "GPU norm after normalization : " << endl;
   qdf_write_doubles(whihinorm_device,wlohinorm_device,
                     whilonorm_device,wlolonorm_device);

   const double tol = 1.0e-48;
   double err = fabs(vhihinorm_device - vhihinorm_host)
              + fabs(whihinorm_device - whihinorm_host)
              + fabs(vlohinorm_device - vlohinorm_host)
              + fabs(wlohinorm_device - wlohinorm_host)
              + fabs(vhilonorm_device - vhilonorm_host)
              + fabs(whilonorm_device - whilonorm_host)
              + fabs(vlolonorm_device - vlolonorm_host)
              + fabs(wlolonorm_device - wlolonorm_host);

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
