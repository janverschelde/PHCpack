/* This program runs automatic tests to compute norms of complex vectors
   in penta double precision, for preset values of the parameters.
   No input is required from the user. */

#include <ctime>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "cmplx5_norm_kernels.h"
#include "cmplx5_norm_host.h"
#include "random5_vectors.h"
#include "penta_double_functions.h"

using namespace std;

void run
 ( int dim, int BS, int freq, int mode, int blocked,
   double *vtbnorm_device, double *vixnorm_device, double *vminorm_device,
   double *vrgnorm_device, double *vpknorm_device,
   double *vtbnorm_host,   double *vixnorm_host,   double *vminorm_host,
   double *vrgnorm_host,   double *vpknorm_host,
   double *wtbnorm_device, double *wixnorm_device, double *wminorm_device,
   double *wrgnorm_device, double *wpknorm_device,
   double *wtbnorm_host,   double *wixnorm_host,   double *wminorm_host,
   double *wrgnorm_host,   double *wpknorm_host );
/*
 * DESCRIPTION :
 *   Computes norms for random complex vectors,
 *   in double doble precision.
 *
 * ON ENTRY :
 *   dim      dimension of the random vector;
 *   BS       block size;
 *   freq     frequency of runs for the tiixngs;
 *   mode     if 0, then only GPU computes,
 *            if 1, then only CPU computes,
 *            if 2, then both GPU and CPU  compute;
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
 *   vixnorm_host    second highest part of 2-norm computed by the host;
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
 *   wpknorm_host    lowest part of 2-norm on host after normalization.  */

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
   double *vtbnorm_device, double *vixnorm_device, double *vminorm_device,
   double *vrgnorm_device, double *vpknorm_device,
   double *vtbnorm_host,   double *vixnorm_host,   double *vminorm_host,
   double *vrgnorm_host,   double *vpknorm_host,
   double *wtbnorm_device, double *wixnorm_device, double *wminorm_device,
   double *wrgnorm_device, double *wpknorm_device,
   double *wtbnorm_host,   double *wixnorm_host,   double *wminorm_host,
   double *wrgnorm_host,   double *wpknorm_host )
{
   const int timevalue = time(NULL); // no fixed seed to verify correctness
   srand(timevalue);

   double* vretb_host = new double[dim]; // highest real parts on host
   double* vreix_host = new double[dim]; // 2nd highest real parts on host
   double* vremi_host = new double[dim]; // middle real parts on host
   double* vrerg_host = new double[dim]; // 2nd lowest real parts on host
   double* vrepk_host = new double[dim]; // lowest real parts on host
   double* vimtb_host = new double[dim]; // highest imaginary parts on host
   double* vimix_host = new double[dim]; // 2nd highest imag parts on host
   double* vimmi_host = new double[dim]; // middle imag parts on host
   double* vimrg_host = new double[dim]; // 2nd lowest imag parts on host
   double* vimpk_host = new double[dim]; // lowest imaginary parts on host
   double* wretb_host = new double[dim]; // a copy for normalization
   double* wreix_host = new double[dim];
   double* wremi_host = new double[dim];
   double* wrerg_host = new double[dim];
   double* wrepk_host = new double[dim];
   double* wimtb_host = new double[dim];
   double* wimix_host = new double[dim];
   double* wimmi_host = new double[dim];
   double* wimrg_host = new double[dim];
   double* wimpk_host = new double[dim];
   double* vretb_device = new double[dim]; // highest real parts on device
   double* vreix_device = new double[dim]; // 2nd highest real parts
   double* vremi_device = new double[dim]; // middle real parts
   double* vrerg_device = new double[dim]; // 2nd lowest real parts 
   double* vrepk_device = new double[dim]; // lowest real parts on device
   double* vimtb_device = new double[dim]; // highest imaginary parts
   double* vimix_device = new double[dim]; // 2nd highest imaginary parts
   double* vimmi_device = new double[dim]; // middle imaginary parts
   double* vimrg_device = new double[dim]; // 2nd lowest imaginary parts
   double* vimpk_device = new double[dim]; // lowest imaginary parts

   random_complex5_vectors
      (dim,vretb_host,vreix_host,vremi_host,vrerg_host,vrepk_host,
           vimtb_host,vimix_host,vimmi_host,vimrg_host,vimpk_host,
       vretb_device,vreix_device,vremi_device,vrerg_device,vrepk_device,
       vimtb_device,vimix_device,vimmi_device,vimrg_device,vimpk_device);

   if(mode == 0 || mode == 2)
   {
      GPU_norm
         (vretb_device,vreix_device,vremi_device,vrerg_device,vrepk_device,
          vimtb_device,vimix_device,vimmi_device,vimrg_device,vimpk_device,
          dim,1,BS,
          vtbnorm_device,vixnorm_device,vminorm_device,vrgnorm_device,
          vpknorm_device,blocked);
      GPU_norm
         (vretb_device,vreix_device,vremi_device,vrerg_device,vrepk_device,
          vimtb_device,vimix_device,vimmi_device,vimrg_device,vimpk_device,
          dim,freq,BS,
          wtbnorm_device,wixnorm_device,wminorm_device,wrgnorm_device,
          wpknorm_device,blocked);
   }

   if(mode == 1 || mode == 2)
   {
      for(int i=0; i<=freq; i++)
      {
         CPU_norm(vretb_host,vreix_host,vremi_host,vrerg_host,vrepk_host,
                  vimtb_host,vimix_host,vimmi_host,vimrg_host,vimpk_host,dim,
                  vtbnorm_host,vixnorm_host,vminorm_host,
                  vrgnorm_host,vpknorm_host);
         make_copy(dim,vretb_host,vreix_host,vremi_host,vrerg_host,vrepk_host,
                       vimtb_host,vimix_host,vimmi_host,vimrg_host,vimpk_host,
                   wretb_host,wreix_host,wremi_host,wrerg_host,wrepk_host,
                   wimtb_host,wimix_host,wimmi_host,wimrg_host,wimpk_host);
         CPU_normalize(wretb_host,wreix_host,wremi_host,wrerg_host,wrepk_host,
                       wimtb_host,wimix_host,wimmi_host,wimrg_host,wimpk_host,
                       dim,*vtbnorm_host,*vixnorm_host,*vminorm_host,
                           *vrgnorm_host,*vpknorm_host);
         CPU_norm(wretb_host,wreix_host,wremi_host,wrerg_host,wrepk_host,
                  wimtb_host,wimix_host,wimmi_host,wimrg_host,wimpk_host,dim,
                  wtbnorm_host,wixnorm_host,wminorm_host,wrgnorm_host,
                  wpknorm_host);
      }
   }
}

int verify_correctness ( int dim, int BS, int blocked )
{
   double vtbnorm_device,vixnorm_device,vminorm_device,vrgnorm_device;
   double vpknorm_device;
   // norm before normalization on device
   double wtbnorm_device,wixnorm_device,wminorm_device,wrgnorm_device;
   double wpknorm_device;
   // norm after normalization on device
   double vtbnorm_host,vixnorm_host,vminorm_host,vrgnorm_host,vpknorm_host;
   // norm before normalization on host
   double wtbnorm_host,wixnorm_host,wminorm_host,wrgnorm_host,wpknorm_host;
   // norm after normalization on host

   run(dim,BS,1,2,blocked,
       &vtbnorm_device,&vixnorm_device,&vminorm_device,&vrgnorm_device,
       &vpknorm_device,
       &vtbnorm_host,  &vixnorm_host,  &vminorm_host,  &vrgnorm_host,
       &vpknorm_host,
       &wtbnorm_device,&wixnorm_device,&wminorm_device,&wrgnorm_device,
       &wpknorm_device,
       &wtbnorm_host,  &wixnorm_host,  &wminorm_host,  &wrgnorm_host,
       &wpknorm_host);

   cout << scientific << setprecision(16);

   cout << "CPU norm : " << endl;
   pdf_write_doubles(vtbnorm_host,vixnorm_host,vminorm_host,
                     vrgnorm_host,vpknorm_host);
   cout << endl;
   cout << "GPU norm : " << endl;
   pdf_write_doubles(vtbnorm_device,vixnorm_device,vminorm_device,
                     vrgnorm_device,vpknorm_device);
   cout << endl;
   cout << "CPU norm after normalization : " << endl;
   pdf_write_doubles(wtbnorm_host,wixnorm_host,wminorm_host,
                     wrgnorm_host,wpknorm_host);
   cout << endl;
   cout << "GPU norm after normalization : " << endl;
   pdf_write_doubles(wtbnorm_device,wixnorm_device,wminorm_device,
                     wrgnorm_device,wpknorm_device);
   cout << endl;

   const double tol = 1.0e-70;
   double err = abs(vtbnorm_device - vtbnorm_host)
              + abs(wtbnorm_device - wtbnorm_host)
              + abs(vixnorm_device - vixnorm_host)
              + abs(wixnorm_device - wixnorm_host)
              + abs(vminorm_device - vminorm_host)
              + abs(wminorm_device - wminorm_host)
              + abs(vrgnorm_device - vrgnorm_host)
              + abs(wrgnorm_device - wrgnorm_host)
              + abs(vpknorm_device - vpknorm_host)
              + abs(wpknorm_device - wpknorm_host);

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
