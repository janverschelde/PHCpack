/* This program runs automatic tests to compute norms of real vectors
   in deca double precision, for preset values of the parameters.
   No input is required from the user. */

#include <ctime>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "deca_double_functions.h"
#include "random10_vectors.h"
#include "dbl10_norm_host.h"
#include "dbl10_norm_kernels.h"

using namespace std;

void run
 ( int dim, int BS, int freq, int mode, int blocked,
   double *vrtbnorm_device, double *vrixnorm_device, double *vrminorm_device,
   double *vrrgnorm_device, double *vrpknorm_device,
   double *vltbnorm_device, double *vlixnorm_device, double *vlminorm_device,
   double *vlrgnorm_device, double *vlpknorm_device,
   double *vrtbnorm_host, double *vrixnorm_host, double *vrminorm_host,
   double *vrrgnorm_host, double *vrpknorm_host,
   double *vltbnorm_host, double *vlixnorm_host, double *vlminorm_host,
   double *vlrgnorm_host, double *vlpknorm_host,
   double *wrtbnorm_device, double *wrixnorm_device, double *wrminorm_device,
   double *wrrgnorm_device, double *wrpknorm_device,
   double *wltbnorm_device, double *wlixnorm_device, double *wlminorm_device,
   double *wlrgnorm_device, double *wlpknorm_device,
   double *wrtbnorm_host, double *wrixnorm_host, double *wrminorm_host,
   double *wrrgnorm_host, double *wrpknorm_host,
   double *wltbnorm_host, double *wlixnorm_host, double *wlminorm_host,
   double *wlrgnorm_host, double *wlpknorm_host );
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
 *   vtbnorm_device  highest part of 2-norm compute by the device;
 *   vixnorm_device  second highest part of 2-norm compute by the device;
 *   vminorm_device  middle part of 2-norm computed by the device;
 *   vrgnorm_device  second lowest part of 2-norm computed by the device;
 *   vpknorm_device  lowest part of 2-norm computed by the device;
 *   vtbnorm_host    highest part of 2-norm compute by the host;
 *   vixnorm_host    second highest part of 2-norm compute by the host;
 *   vminorm_host    middle part of 2-norm computed by the host;
 *   vrgnorm_host    second lowest part of 2-norm computed by the host;
 *   vpknorm_host    lowest part of 2-norm computed by the host;
 *   wtbnorm_device  highest part of 2-norm on device after normalization;
 *   wixnorm_device  2nd highest part of 2-norm on device after normalization;
 *   wminorm_device  middle part of 2-norm on device after normalization;
 *   wrgnorm_device  2nd lowest part of 2-norm on device after normalization;
 *   wpknorm_device  lowest part of 2-norm on device after normalization;
 *   wtbnorm_host    highest part of 2-norm on host after normalization;
 *   wixnorm_host    2nd highest part of 2-norm on host after normalization;
 *   wminorm_host    middle part of 2-norm on host after normalization;
 *   wrgnorm_host    2nd lowest part of 2-norm on host after normalization;
 *   wpknorm_host    lowest part of 2-norm on host after normalization. */

int verify_correctness ( int dim, int BS, int blocked );
/*
 * Computes norms of real vectors of dimension dim, block size BS,
 * and verifies the correctness.  Returns 1 if the error was too large,
 * returns 0 if the 2-norm and normalization was within the tolerance.
 * If blocked and dim is a multiple of BS, then many blocks compute. */

int main ( void )
{
   int fail;

   cout << "Verifying correctness for dimension and block size 64 ..."
        << endl;
   fail = verify_correctness(64,64,0);


   cout << "Nonblocked version for dimension 1024 and block size 64 ..."
        << endl;
   fail = verify_correctness(1024,64,0);

   cout << "Blocked version for dimension 4096 and block size 64 ..."
        << endl;
   fail = verify_correctness(4096,64,1);

   return 0;
}

void run 
 ( int dim, int BS, int freq, int mode, int blocked,
   double *vrtbnorm_device, double *vrixnorm_device, double *vrminorm_device,
   double *vrrgnorm_device, double *vrpknorm_device,
   double *vltbnorm_device, double *vlixnorm_device, double *vlminorm_device,
   double *vlrgnorm_device, double *vlpknorm_device,
   double *vrtbnorm_host, double *vrixnorm_host, double *vrminorm_host,
   double *vrrgnorm_host, double *vrpknorm_host,
   double *vltbnorm_host, double *vlixnorm_host, double *vlminorm_host,
   double *vlrgnorm_host, double *vlpknorm_host,
   double *wrtbnorm_device, double *wrixnorm_device, double *wrminorm_device,
   double *wrrgnorm_device, double *wrpknorm_device,
   double *wltbnorm_device, double *wlixnorm_device, double *wlminorm_device,
   double *wlrgnorm_device, double *wlpknorm_device,
   double *wrtbnorm_host, double *wrixnorm_host, double *wrminorm_host,
   double *wrrgnorm_host, double *wrpknorm_host,
   double *wltbnorm_host, double *wlixnorm_host, double *wlminorm_host,
   double *wlrgnorm_host, double *wlpknorm_host )
{
   const int timevalue = time(NULL); // no fixed seed to verify correctness
   srand(timevalue);

   double* vrtb_host = new double[dim];   // highest parts on the host
   double* vrix_host = new double[dim];   // second highest parts on the host
   double* vrmi_host = new double[dim];   // third highest parts on the host
   double* vrrg_host = new double[dim];   // fourth highest parts on the host
   double* vrpk_host = new double[dim];   // fifth highest parts on the host
   double* vltb_host = new double[dim];   // fifth lowest parts the host
   double* vlix_host = new double[dim];   // fourth lowest parts on the host
   double* vlmi_host = new double[dim];   // third lowest parts on the host
   double* vlrg_host = new double[dim];   // second lowest parts on the host
   double* vlpk_host = new double[dim];   // lowest parts on the host
   double* vrtb_device = new double[dim]; // highest parts on the device
   double* vrix_device = new double[dim]; // second highest parts on the device
   double* vrmi_device = new double[dim]; // third highest parts on the device
   double* vrrg_device = new double[dim]; // fourth highest parts on the device
   double* vrpk_device = new double[dim]; // fifth highest parts on the device
   double* vltb_device = new double[dim]; // fifth lowest parts on the device
   double* vlix_device = new double[dim]; // fourth lowest parts on the device
   double* vlmi_device = new double[dim]; // third lowest parts on the device
   double* vlrg_device = new double[dim]; // second lowest parts on the device
   double* vlpk_device = new double[dim]; // lowest parts on the device
   double* wrtb_host = new double[dim];   // highest parts copy
   double* wrix_host = new double[dim];   // second highest parts copy
   double* wrmi_host = new double[dim];   // third highest parts copy
   double* wrrg_host = new double[dim];   // fourth highest parts copy
   double* wrpk_host = new double[dim];   // fifth highest parts copy
   double* wltb_host = new double[dim];   // fifth lowest parts copy
   double* wlix_host = new double[dim];   // fourth lowest parts copy
   double* wlmi_host = new double[dim];   // third lowest parts copy
   double* wlrg_host = new double[dim];   // second lowest parts copy
   double* wlpk_host = new double[dim];   // lowest parts copy

   random_double10_vectors
      (dim,vrtb_host,vrix_host,vrmi_host,vrrg_host,vrpk_host,
           vltb_host,vlix_host,vlmi_host,vlrg_host,vlpk_host,
           vrtb_device,vrix_device,vrmi_device,vrrg_device,vrpk_device,
           vltb_device,vlix_device,vlmi_device,vlrg_device,vlpk_device);
 
   if(mode == 0 || mode == 2)
   {
      GPU_norm(vrtb_device,vrix_device,vrmi_device,vrrg_device,vrpk_device,
               vltb_device,vlix_device,vlmi_device,vlrg_device,vlpk_device,
               dim,1,BS,vrtbnorm_device,vrixnorm_device,vrminorm_device,
               vrrgnorm_device,vrpknorm_device,vltbnorm_device,
               vlixnorm_device,vlminorm_device,vlrgnorm_device,
               vlpknorm_device,blocked);
      GPU_norm(vrtb_device,vrix_device,vrmi_device,vrrg_device,vrpk_device,
               vltb_device,vlix_device,vlmi_device,vlrg_device,vlpk_device,
               dim,freq,BS,wrtbnorm_device,wrixnorm_device,wrminorm_device,
               wrrgnorm_device,wrpknorm_device,wltbnorm_device,
               wlixnorm_device,wlminorm_device,wlrgnorm_device,
               wlpknorm_device,blocked);
   }

   if(mode == 1 || mode == 2)
   {
      for(int i=0; i<=freq; i++)
      {
         CPU_norm(vrtb_host,vrix_host,vrmi_host,vrrg_host,vrpk_host,
                  vltb_host,vlix_host,vlmi_host,vlrg_host,vlpk_host,dim,
                  vrtbnorm_host,vrixnorm_host,vrminorm_host,
                  vrrgnorm_host,vrpknorm_host,
                  vltbnorm_host,vlixnorm_host,vlminorm_host,
                  vlrgnorm_host,vlpknorm_host);
         make_copy(dim,vrtb_host,vrix_host,vrmi_host,vrrg_host,vrpk_host,
                       vltb_host,vlix_host,vlmi_host,vlrg_host,vlpk_host,
                       wrtb_host,wrix_host,wrmi_host,wrrg_host,wrpk_host,
                       wltb_host,wlix_host,wlmi_host,wlrg_host,wlpk_host);
         CPU_normalize(wrtb_host,wrix_host,wrmi_host,wrrg_host,wrpk_host,
                       wltb_host,wlix_host,wlmi_host,wlrg_host,wlpk_host,dim,
                       *vrtbnorm_host,*vrixnorm_host,*vrminorm_host,
                       *vrrgnorm_host,*vrpknorm_host,
                       *vltbnorm_host,*vlixnorm_host,*vlminorm_host,
                       *vlrgnorm_host,*vlpknorm_host);
         CPU_norm(wrtb_host,wrix_host,wrmi_host,wrrg_host,wrpk_host,
                  wltb_host,wlix_host,wlmi_host,wlrg_host,wlpk_host,dim,
                  wrtbnorm_host,wrixnorm_host,wrminorm_host,
                  wrrgnorm_host,wrpknorm_host,
                  wltbnorm_host,wlixnorm_host,wlminorm_host,
                  wlrgnorm_host,wlpknorm_host);
      }
   }
}

int verify_correctness ( int dim, int BS, int blocked )
{
   double vrtbnorm_device,vrixnorm_device,vrminorm_device;
   double vrrgnorm_device,vrpknorm_device;
   double vltbnorm_device,vlixnorm_device,vlminorm_device;
   double vlrgnorm_device,vlpknorm_device;
   // norm before normalization on device
   double wrtbnorm_device,wrixnorm_device,wrminorm_device;
   double wrrgnorm_device,wrpknorm_device;
   double wltbnorm_device,wlixnorm_device,wlminorm_device;
   double wlrgnorm_device,wlpknorm_device;
   // norm after normalization on device
   double vrtbnorm_host,vrixnorm_host,vrminorm_host;
   double vrrgnorm_host,vrpknorm_host;
   double vltbnorm_host,vlixnorm_host,vlminorm_host;
   double vlrgnorm_host,vlpknorm_host;
   // norm before normalization on host
   double wrtbnorm_host,wrixnorm_host,wrminorm_host;
   double wrrgnorm_host,wrpknorm_host;
   double wltbnorm_host,wlixnorm_host,wlminorm_host;
   double wlrgnorm_host,wlpknorm_host;
   // norm after normalization on host

   run(dim,BS,1,2,blocked,
       &vrtbnorm_device,&vrixnorm_device,&vrminorm_device,
       &vrrgnorm_device,&vrpknorm_device,
       &vltbnorm_device,&vlixnorm_device,&vlminorm_device,
       &vlrgnorm_device,&vlpknorm_device,
       &vrtbnorm_host,&vrixnorm_host,&vrminorm_host,
       &vrrgnorm_host,&vrpknorm_host,
       &vltbnorm_host,&vlixnorm_host,&vlminorm_host,
       &vlrgnorm_host,&vlpknorm_host,
       &wrtbnorm_device,&wrixnorm_device,&wrminorm_device,
       &wrrgnorm_device,&wrpknorm_device,
       &wltbnorm_device,&wlixnorm_device,&wlminorm_device,
       &wlrgnorm_device,&wlpknorm_device,
       &wrtbnorm_host,&wrixnorm_host,&wrminorm_host,
       &wrrgnorm_host,&wrpknorm_host,
       &wltbnorm_host,&wlixnorm_host,&wlminorm_host,
       &wlrgnorm_host,&wlpknorm_host);

   cout << scientific << setprecision(16);

   cout << "CPU norm : " << endl;
   daf_write_doubles(vrtbnorm_host,vrixnorm_host,vrminorm_host,
                     vrrgnorm_host,vrpknorm_host,
                     vltbnorm_host,vlixnorm_host,vlminorm_host,
                     vlrgnorm_host,vlpknorm_host);
   cout << "CPU norm after normalization : " << endl;
   daf_write_doubles(wrtbnorm_host,wrixnorm_host,wrminorm_host,
                     wrrgnorm_host,wrpknorm_host,
                     wltbnorm_host,wlixnorm_host,wlminorm_host,
                     wlrgnorm_host,wlpknorm_host);
   cout << "GPU norm : " << endl;
   daf_write_doubles(vrtbnorm_device,vrixnorm_device,vrminorm_device,
                     vrrgnorm_device,vrpknorm_device,
                     vltbnorm_device,vlixnorm_device,vlminorm_device,
                     vlrgnorm_device,vlpknorm_device);
   cout << "GPU norm after normalization : " << endl;
   daf_write_doubles(wrtbnorm_device,wrixnorm_device,wrminorm_device,
                     wrrgnorm_device,wrpknorm_device,
                     wltbnorm_device,wlixnorm_device,wlminorm_device,
                     wlrgnorm_device,wlpknorm_device);
   const double tol = 1.0e-150;
   double err = fabs(vrtbnorm_device - vrtbnorm_host)
              + fabs(wrtbnorm_device - wrtbnorm_host)
              + fabs(vrixnorm_device - vrixnorm_host)
              + fabs(wrixnorm_device - wrixnorm_host)
              + fabs(vrminorm_device - vrminorm_host)
              + fabs(wrminorm_device - wrminorm_host)
              + fabs(vrrgnorm_device - vrrgnorm_host)
              + fabs(wrrgnorm_device - wrrgnorm_host)
              + fabs(vrpknorm_device - vrpknorm_host)
              + fabs(wrpknorm_device - wrpknorm_host)
              + fabs(vltbnorm_device - vltbnorm_host)
              + fabs(wltbnorm_device - wltbnorm_host)
              + fabs(vlixnorm_device - vlixnorm_host)
              + fabs(wlixnorm_device - wlixnorm_host)
              + fabs(vlminorm_device - vlminorm_host)
              + fabs(wlminorm_device - wlminorm_host)
              + fabs(vlrgnorm_device - vlrgnorm_host)
              + fabs(wlrgnorm_device - wlrgnorm_host)
              + fabs(vlpknorm_device - vlpknorm_host)
              + fabs(wlpknorm_device - wlpknorm_host);

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
