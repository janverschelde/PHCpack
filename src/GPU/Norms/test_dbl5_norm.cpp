/* This program runs automatic tests to compute norms of real vectors
   in penta double precision, for preset values of the parameters.
   No input is required from the user. */

#include <ctime>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "penta_double_functions.h"
#include "random5_vectors.h"
#include "dbl5_norm_host.h"
#include "dbl5_norm_kernels.h"

using namespace std;

void run
 ( int dim, int BS, int freq, int mode, int blocked,
   double *vtbnorm_device, double *vixnorm_device, double *vminorm_device,
   double *vrgnorm_device, double *vpknorm_device,
   double *vtbnorm_host, double *vixnorm_host, double *vminorm_host,
   double *vrgnorm_host, double *vpknorm_host,
   double *wtbnorm_device, double *wixnorm_device, double *wminorm_device,
   double *wrgnorm_device, double *wpknorm_device,
   double *wtbnorm_host, double *wixnorm_host, double *wminorm_host,
   double *wrgnorm_host, double *wpknorm_host );
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
   double *vtbnorm_device, double *vixnorm_device, double *vminorm_device,
   double *vrgnorm_device, double *vpknorm_device,
   double *vtbnorm_host, double *vixnorm_host, double *vminorm_host,
   double *vrgnorm_host, double *vpknorm_host,
   double *wtbnorm_device, double *wixnorm_device, double *wminorm_device,
   double *wrgnorm_device, double *wpknorm_device,
   double *wtbnorm_host, double *wixnorm_host, double *wminorm_host,
   double *wrgnorm_host, double *wpknorm_host )
{
   const int timevalue = time(NULL); // no fixed seed to verify correctness
   srand(timevalue);

   double* vtb_host = new double[dim];   // highest parts on the host
   double* vix_host = new double[dim];   // second highest parts on the host
   double* vmi_host = new double[dim];   // middle parts on the host
   double* vrg_host = new double[dim];   // second lowest parts on the host
   double* vpk_host = new double[dim];   // lowest parts on the host
   double* vtb_device = new double[dim]; // highest parts on the device
   double* vix_device = new double[dim]; // second highest parts on the device
   double* vmi_device = new double[dim]; // middle parts on the device
   double* vrg_device = new double[dim]; // second lowest parts on the device
   double* vpk_device = new double[dim]; // lowest parts on the device
   double* wtb_host = new double[dim];   // highest parts copy
   double* wix_host = new double[dim];   // second highest parts copy
   double* wmi_host = new double[dim];   // middle parts copy
   double* wrg_host = new double[dim];   // second lowest parts copy
   double* wpk_host = new double[dim];   // lowest parts copy

   random_double5_vectors
      (dim,vtb_host,vix_host,vmi_host,vrg_host,vpk_host,
           vtb_device,vix_device,vmi_device,vrg_device,vpk_device);

   if(mode == 0 || mode == 2)
   {
      GPU_norm(vtb_device,vix_device,vmi_device,vrg_device,vpk_device,
               dim,1,BS, vtbnorm_device,vixnorm_device,vminorm_device,
               vrgnorm_device,vpknorm_device,blocked);
      GPU_norm(vtb_device,vix_device,vmi_device,vrg_device,vpk_device,
               dim,freq,BS,wtbnorm_device,wixnorm_device,wminorm_device,
               wrgnorm_device,wpknorm_device,blocked);
   }

   if(mode == 1 || mode == 2)
   {
      for(int i=0; i<=freq; i++)
      {
         CPU_norm(vtb_host,vix_host,vmi_host,vrg_host,vpk_host,dim,
                  vtbnorm_host,vixnorm_host,vminorm_host,
                  vrgnorm_host,vpknorm_host);
         make_copy(dim,vtb_host,vix_host,vmi_host,vrg_host,vpk_host,
                       wtb_host,wix_host,wmi_host,wrg_host,wpk_host);
         CPU_normalize(wtb_host,wix_host,wmi_host,wrg_host,wpk_host,dim,
                       *vtbnorm_host,*vixnorm_host,*vminorm_host,
                       *vrgnorm_host,*vpknorm_host);
         CPU_norm(wtb_host,wix_host,wmi_host,wrg_host,wpk_host,dim,
                  wtbnorm_host,wixnorm_host,wminorm_host,
                  wrgnorm_host,wpknorm_host);
      }
   }
}

int verify_correctness ( int dim, int BS, int blocked )
{
   double vtbnorm_device,vixnorm_device,vminorm_device;
   double vrgnorm_device,vpknorm_device;
   // norm before normalization on device
   double wtbnorm_device,wixnorm_device,wminorm_device;
   double wrgnorm_device,wpknorm_device;
   // norm after normalization on device
   double vtbnorm_host,vixnorm_host,vminorm_host;
   double vrgnorm_host,vpknorm_host;
   // norm before normalization on host
   double wtbnorm_host,wixnorm_host,wminorm_host;
   double wrgnorm_host,wpknorm_host;
   // norm after normalization on host

   run(dim,BS,1,2,blocked,
       &vtbnorm_device,&vixnorm_device,&vminorm_device,
       &vrgnorm_device,&vpknorm_device,
       &vtbnorm_host,&vixnorm_host,&vminorm_host,
       &vrgnorm_host,&vpknorm_host,
       &wtbnorm_device,&wixnorm_device,&wminorm_device,
       &wrgnorm_device,&wpknorm_device,
       &wtbnorm_host,&wixnorm_host,&wminorm_host,
       &wrgnorm_host,&wpknorm_host);

   cout << scientific << setprecision(16);

   cout << "CPU norm : " << endl;
   pdf_write_doubles(vtbnorm_host,vixnorm_host,vminorm_host,
                     vrgnorm_host,vpknorm_host);
   cout << endl;
   cout << "CPU norm after normalization : " << endl;
   pdf_write_doubles(wtbnorm_host,wixnorm_host,wminorm_host,
                     wrgnorm_host,wpknorm_host);
   cout << endl;

   cout << "GPU norm : " << endl;
   pdf_write_doubles(vtbnorm_device,vixnorm_device,vminorm_device,
                     vrgnorm_device,vpknorm_device);
   cout << endl;
   cout << "GPU norm after normalization : " << endl;
   pdf_write_doubles(wtbnorm_device,wixnorm_device,wminorm_device,
                     wrgnorm_device,wpknorm_device);
   cout << endl;

   const double tol = 1.0e-60;
   double err = fabs(vtbnorm_device - vtbnorm_host)
              + fabs(wtbnorm_device - wtbnorm_host)
              + fabs(vixnorm_device - vixnorm_host)
              + fabs(wixnorm_device - wixnorm_host)
              + fabs(vminorm_device - vminorm_host)
              + fabs(wminorm_device - wminorm_host)
              + fabs(vrgnorm_device - vrgnorm_host)
              + fabs(wrgnorm_device - wrgnorm_host)
              + fabs(vpknorm_device - vpknorm_host)
              + fabs(wpknorm_device - wpknorm_host);

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
