/* This program runs automatic tests to compute norms of complex vectors
   in deca double precision, for preset values of the parameters.
   No input is required from the user. */

#include <ctime>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "cmplx10_norm_kernels.h"
#include "cmplx10_norm_host.h"
#include "random10_vectors.h"
#include "deca_double_functions.h"

using namespace std;

void run
 ( int dim, int BS, int freq, int mode, int blocked,
   double *vrtbnorm_device, double *vrixnorm_device, double *vrminorm_device,
   double *vrrgnorm_device, double *vrpknorm_device,
   double *vltbnorm_device, double *vlixnorm_device, double *vlminorm_device,
   double *vlrgnorm_device, double *vlpknorm_device,
   double *vrtbnorm_host,   double *vrixnorm_host,   double *vrminorm_host,
   double *vrrgnorm_host,   double *vrpknorm_host,
   double *vltbnorm_host,   double *vlixnorm_host,   double *vlminorm_host,
   double *vlrgnorm_host,   double *vlpknorm_host,
   double *wrtbnorm_device, double *wrixnorm_device, double *wrminorm_device,
   double *wrrgnorm_device, double *wrpknorm_device,
   double *wltbnorm_device, double *wlixnorm_device, double *wlminorm_device,
   double *wlrgnorm_device, double *wlpknorm_device,
   double *wrtbnorm_host,   double *wrixnorm_host,   double *wrminorm_host,
   double *wrrgnorm_host,   double *wrpknorm_host,
   double *wltbnorm_host,   double *wlixnorm_host,   double *wlminorm_host,
   double *wlrgnorm_host,   double *wlpknorm_host );
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

   cout << "Verifying correctness for dimension and block size 128 ..."
        << endl;
   fail = verify_correctness(128,128,0);

   cout << "Nonblocked version form dimension 1024 and block size 128 ..."
        << endl;
   fail = verify_correctness(1024,128,0);

   cout << "Blocked version form dimension 4096 and block size 128 ..."
        << endl;
   fail = verify_correctness(4096,128,1);

   return 0;
}

void run
 ( int dim, int BS, int freq, int mode, int blocked,
   double *vrtbnorm_device, double *vrixnorm_device, double *vrminorm_device,
   double *vrrgnorm_device, double *vrpknorm_device,
   double *vltbnorm_device, double *vlixnorm_device, double *vlminorm_device,
   double *vlrgnorm_device, double *vlpknorm_device,
   double *vrtbnorm_host,   double *vrixnorm_host,   double *vrminorm_host,
   double *vrrgnorm_host,   double *vrpknorm_host,
   double *vltbnorm_host,   double *vlixnorm_host,   double *vlminorm_host,
   double *vlrgnorm_host,   double *vlpknorm_host,
   double *wrtbnorm_device, double *wrixnorm_device, double *wrminorm_device,
   double *wrrgnorm_device, double *wrpknorm_device,
   double *wltbnorm_device, double *wlixnorm_device, double *wlminorm_device,
   double *wlrgnorm_device, double *wlpknorm_device,
   double *wrtbnorm_host,   double *wrixnorm_host,   double *wrminorm_host,
   double *wrrgnorm_host,   double *wrpknorm_host,
   double *wltbnorm_host,   double *wlixnorm_host,   double *wlminorm_host,
   double *wlrgnorm_host,   double *wlpknorm_host )
{
   const int timevalue = time(NULL); // no fixed seed to verify correctness
   srand(timevalue);

   double* vrertb_host = new double[dim]; // highest real parts on host
   double* vrerix_host = new double[dim]; // 2nd highest real parts on host
   double* vrermi_host = new double[dim]; // 3rd highest real parts on host
   double* vrerrg_host = new double[dim]; // 4th highest real parts on host
   double* vrerpk_host = new double[dim]; // 5th highest real parts on host
   double* vreltb_host = new double[dim]; // 5th lowest real parts on host
   double* vrelix_host = new double[dim]; // 4th lowest real parts on host
   double* vrelmi_host = new double[dim]; // 3rd lowest real parts on host
   double* vrelrg_host = new double[dim]; // 2nd lowest real parts on host
   double* vrelpk_host = new double[dim]; // lowest real parts on host
   double* vimrtb_host = new double[dim]; // highest imaginary parts on host
   double* vimrix_host = new double[dim]; // 2nd highest imag parts on host
   double* vimrmi_host = new double[dim]; // 3rd highest imag parts on host
   double* vimrrg_host = new double[dim]; // 4th highest imag parts on host
   double* vimrpk_host = new double[dim]; // 5th highest imag parts on host
   double* vimltb_host = new double[dim]; // 5th lowest imag parts on host
   double* vimlix_host = new double[dim]; // 4th lowest imag parts on host
   double* vimlmi_host = new double[dim]; // 3rd lowest imag parts on host
   double* vimlrg_host = new double[dim]; // 2nd lowest imag parts on host
   double* vimlpk_host = new double[dim]; // lowest imaginary parts on host
   double* wrertb_host = new double[dim]; // a copy for normalization
   double* wrerix_host = new double[dim];
   double* wrermi_host = new double[dim];
   double* wrerrg_host = new double[dim];
   double* wrerpk_host = new double[dim];
   double* wreltb_host = new double[dim]; 
   double* wrelix_host = new double[dim];
   double* wrelmi_host = new double[dim];
   double* wrelrg_host = new double[dim];
   double* wrelpk_host = new double[dim];
   double* wimrtb_host = new double[dim];
   double* wimrix_host = new double[dim];
   double* wimrmi_host = new double[dim];
   double* wimrrg_host = new double[dim];
   double* wimrpk_host = new double[dim];
   double* wimltb_host = new double[dim];
   double* wimlix_host = new double[dim];
   double* wimlmi_host = new double[dim];
   double* wimlrg_host = new double[dim];
   double* wimlpk_host = new double[dim];
   double* vrertb_device = new double[dim]; // highest real parts on device
   double* vrerix_device = new double[dim]; // 2nd highest real parts
   double* vrermi_device = new double[dim]; // 3rd highest real parts
   double* vrerrg_device = new double[dim]; // 4th highest real parts 
   double* vrerpk_device = new double[dim]; // 5th highest real parts 
   double* vreltb_device = new double[dim]; // 5th lowest real parts
   double* vrelix_device = new double[dim]; // 4th lowest real parts
   double* vrelmi_device = new double[dim]; // 3rd lowest real parts
   double* vrelrg_device = new double[dim]; // 2nd lowest real parts 
   double* vrelpk_device = new double[dim]; // lowest real parts on device
   double* vimrtb_device = new double[dim]; // highest imaginary parts
   double* vimrix_device = new double[dim]; // 2nd highest imaginary parts
   double* vimrmi_device = new double[dim]; // 3rd highest imaginary parts
   double* vimrrg_device = new double[dim]; // 4th highest imaginary parts
   double* vimrpk_device = new double[dim]; // 5th highest imaginary parts
   double* vimltb_device = new double[dim]; // 5th lowest imaginary parts
   double* vimlix_device = new double[dim]; // 4th lowest imaginary parts
   double* vimlmi_device = new double[dim]; // 3rd lowest imaginary parts
   double* vimlrg_device = new double[dim]; // 2nd lowest imaginary parts
   double* vimlpk_device = new double[dim]; // lowest imaginary parts

   random_complex10_vectors
      (dim,vrertb_host,vrerix_host,vrermi_host,vrerrg_host,vrerpk_host,
           vreltb_host,vrelix_host,vrelmi_host,vrelrg_host,vrelpk_host,
           vimrtb_host,vimrix_host,vimrmi_host,vimrrg_host,vimrpk_host,
           vimltb_host,vimlix_host,vimlmi_host,vimlrg_host,vimlpk_host,
       vrertb_device,vrerix_device,vrermi_device,vrerrg_device,vrerpk_device,
       vreltb_device,vrelix_device,vrelmi_device,vrelrg_device,vrelpk_device,
       vimrtb_device,vimrix_device,vimrmi_device,vimrrg_device,vimrpk_device,
       vimltb_device,vimlix_device,vimlmi_device,vimlrg_device,vimlpk_device);
/*
   for(int k=0; k<dim; k++)
   {
      vimrtb_host[k] = 0.0; vimrtb_device[k] = 0.0;
      vimrix_host[k] = 0.0; vimrix_device[k] = 0.0;
      vimrmi_host[k] = 0.0; vimrmi_device[k] = 0.0;
      vimrrg_host[k] = 0.0; vimrrg_device[k] = 0.0;
      vimrpk_host[k] = 0.0; vimrpk_device[k] = 0.0;
      vimltb_host[k] = 0.0; vimltb_device[k] = 0.0;
      vimlix_host[k] = 0.0; vimlix_device[k] = 0.0;
      vimlmi_host[k] = 0.0; vimlmi_device[k] = 0.0;
      vimlrg_host[k] = 0.0; vimlrg_device[k] = 0.0;
      vrertb_host[k] = 0.0; vrertb_device[k] = 0.0;
      vrerix_host[k] = 0.0; vrerix_device[k] = 0.0;
      vrermi_host[k] = 0.0; vrermi_device[k] = 0.0;
      vrerrg_host[k] = 0.0; vrerrg_device[k] = 0.0;
      vrerpk_host[k] = 0.0; vrerpk_device[k] = 0.0;
      vreltb_host[k] = 0.0; vreltb_device[k] = 0.0;
      vrelix_host[k] = 0.0; vrelix_device[k] = 0.0;
      vrelmi_host[k] = 0.0; vrelmi_device[k] = 0.0;
      vrelrg_host[k] = 0.0; vrelrg_device[k] = 0.0;
      vrelpk_host[k] = 0.0; vrelpk_device[k] = 0.0; 
   }
 */
   // vrertb_host[0] = 1.0; vrertb_device[0] = 1.0;

   double diffsum = 0.0;
   for(int k=0; k<dim; k++)
      diffsum += fabs(vrertb_host[k] - vrertb_device[k])
               + fabs(vrerix_host[k] - vrerix_device[k])
               + fabs(vrermi_host[k] - vrermi_device[k])
               + fabs(vrerrg_host[k] - vrerrg_device[k])
               + fabs(vrerpk_host[k] - vrerpk_device[k])
               + fabs(vreltb_host[k] - vreltb_device[k])
               + fabs(vrelix_host[k] - vrelix_device[k])
               + fabs(vrelmi_host[k] - vrelmi_device[k])
               + fabs(vrelrg_host[k] - vrelrg_device[k])
               + fabs(vrelpk_host[k] - vrelpk_device[k])
               + fabs(vimrtb_host[k] - vimrtb_device[k])
               + fabs(vimrix_host[k] - vimrix_device[k])
               + fabs(vimrmi_host[k] - vimrmi_device[k])
               + fabs(vimrrg_host[k] - vimrrg_device[k])
               + fabs(vimrpk_host[k] - vimrpk_device[k])
               + fabs(vimltb_host[k] - vimltb_device[k])
               + fabs(vimlix_host[k] - vimlix_device[k])
               + fabs(vimlmi_host[k] - vimlmi_device[k])
               + fabs(vimlrg_host[k] - vimlrg_device[k])
               + fabs(vimlpk_host[k] - vimlpk_device[k]);

   cout << "dim : " << dim << endl;
   cout << "BS : " << BS << endl;
   cout << "diffsum : " << diffsum << endl;

   if(mode == 0 || mode == 2)
   {
      GPU_norm
         (vrertb_device,vrerix_device,vrermi_device,vrerrg_device,vrerpk_device,
          vreltb_device,vrelix_device,vrelmi_device,vrelrg_device,vrelpk_device,
          vimrtb_device,vimrix_device,vimrmi_device,vimrrg_device,vimrpk_device,
          vimltb_device,vimlix_device,vimlmi_device,vimlrg_device,vimlpk_device,
          dim,1,BS,
          vrtbnorm_device,vrixnorm_device,vrminorm_device,vrrgnorm_device,
          vrpknorm_device,vltbnorm_device,vlixnorm_device,vlminorm_device,
          vlrgnorm_device,vlpknorm_device,blocked);

      cout << "vrixnorm_device : " << *vrixnorm_device << endl;

      GPU_norm
         (vrertb_device,vrerix_device,vrermi_device,vrerrg_device,vrerpk_device,
          vreltb_device,vrelix_device,vrelmi_device,vrelrg_device,vrelpk_device,
          vimrtb_device,vimrix_device,vimrmi_device,vimrrg_device,vimrpk_device,
          vimltb_device,vimlix_device,vimlmi_device,vimlrg_device,vimlpk_device,
          dim,freq,BS,
          wrtbnorm_device,wrixnorm_device,wrminorm_device,wrrgnorm_device,
          wrpknorm_device,wltbnorm_device,wlixnorm_device,wlminorm_device,
          wlrgnorm_device,wlpknorm_device,blocked);
   }

   if(mode == 1 || mode == 2)
   {
      for(int i=0; i<=freq; i++)
      {
         CPU_norm
            (vrertb_host,vrerix_host,vrermi_host,vrerrg_host,vrerpk_host,
             vreltb_host,vrelix_host,vrelmi_host,vrelrg_host,vrelpk_host,
             vimrtb_host,vimrix_host,vimrmi_host,vimrrg_host,vimrpk_host,
             vimltb_host,vimlix_host,vimlmi_host,vimlrg_host,vimlpk_host,dim,
             vrtbnorm_host,vrixnorm_host,vrminorm_host,vrrgnorm_host,
             vrpknorm_host,vltbnorm_host,vlixnorm_host,vlminorm_host,
             vlrgnorm_host,vlpknorm_host);
         make_copy
            (dim,vrertb_host,vrerix_host,vrermi_host,vrerrg_host,vrerpk_host,
                 vreltb_host,vrelix_host,vrelmi_host,vrelrg_host,vrelpk_host,
                 vimrtb_host,vimrix_host,vimrmi_host,vimrrg_host,vimrpk_host,
                 vimltb_host,vimlix_host,vimlmi_host,vimlrg_host,vimlpk_host,
             wrertb_host,wrerix_host,wrermi_host,wrerrg_host,wrerpk_host,
             wreltb_host,wrelix_host,wrelmi_host,wrelrg_host,wrelpk_host,
             wimrtb_host,wimrix_host,wimrmi_host,wimrrg_host,wimrpk_host,
             wimltb_host,wimlix_host,wimlmi_host,wimlrg_host,wimlpk_host);
         CPU_normalize
            (wrertb_host,wrerix_host,wrermi_host,wrerrg_host,wrerpk_host,
             wreltb_host,wrelix_host,wrelmi_host,wrelrg_host,wrelpk_host,
             wimrtb_host,wimrix_host,wimrmi_host,wimrrg_host,wimrpk_host,
             wimltb_host,wimlix_host,wimlmi_host,wimlrg_host,wimlpk_host,
             dim,*vrtbnorm_host,*vrixnorm_host,*vrminorm_host,
                 *vrrgnorm_host,*vrpknorm_host,*vltbnorm_host,
             *vlixnorm_host,*vlminorm_host,*vlrgnorm_host,*vlpknorm_host);
         CPU_norm
            (wrertb_host,wrerix_host,wrermi_host,wrerrg_host,wrerpk_host,
             wreltb_host,wrelix_host,wrelmi_host,wrelrg_host,wrelpk_host,
             wimrtb_host,wimrix_host,wimrmi_host,wimrrg_host,wimrpk_host,
             wimltb_host,wimlix_host,wimlmi_host,wimlrg_host,wimlpk_host,dim,
             wrtbnorm_host,wrixnorm_host,wrminorm_host,wrrgnorm_host,
             wrpknorm_host,wltbnorm_host,wlixnorm_host,wlminorm_host,
             wlrgnorm_host,wlpknorm_host);
      }
   }
}

int verify_correctness ( int dim, int BS, int blocked )
{
   double vrtbnorm_device,vrixnorm_device,vrminorm_device,vrrgnorm_device;
   double vrpknorm_device,vltbnorm_device,vlixnorm_device,vlminorm_device;
   double vlrgnorm_device,vlpknorm_device;
   // norm before normalization on device
   double wrtbnorm_device,wrixnorm_device,wrminorm_device,wrrgnorm_device;
   double wrpknorm_device,wltbnorm_device,wlixnorm_device,wlminorm_device;
   double wlrgnorm_device,wlpknorm_device;
   // norm after normalization on device
   double vrtbnorm_host,vrixnorm_host,vrminorm_host,vrrgnorm_host;
   double vrpknorm_host,vltbnorm_host,vlixnorm_host,vlminorm_host;
   double vlrgnorm_host,vlpknorm_host;
   // norm before normalization on host
   double wrtbnorm_host,wrixnorm_host,wrminorm_host,wrrgnorm_host;
   double wrpknorm_host,wltbnorm_host,wlixnorm_host,wlminorm_host;
   double wlrgnorm_host,wlpknorm_host;
   // norm after normalization on host

   run(dim,BS,1,2,blocked,
       &vrtbnorm_device,&vrixnorm_device,&vrminorm_device,&vrrgnorm_device,
       &vrpknorm_device,&vltbnorm_device,&vlixnorm_device,&vlminorm_device,
       &vlrgnorm_device,&vlpknorm_device,
       &vrtbnorm_host,  &vrixnorm_host,  &vrminorm_host,  &vrrgnorm_host,
       &vrpknorm_host,  &vltbnorm_host,  &vlixnorm_host,  &vlminorm_host,
       &vlrgnorm_host,  &vlpknorm_host,
       &wrtbnorm_device,&wrixnorm_device,&wrminorm_device,&wrrgnorm_device,
       &wrpknorm_device,&wltbnorm_device,&wlixnorm_device,&wlminorm_device,
       &wlrgnorm_device,&wlpknorm_device,
       &wrtbnorm_host,  &wrixnorm_host,  &wrminorm_host,  &wrrgnorm_host,
       &wrpknorm_host,  &wltbnorm_host,  &wlixnorm_host,  &wlminorm_host,
       &wlrgnorm_host,  &wlpknorm_host);

   cout << scientific << setprecision(16);

   cout << "CPU norm : " << endl;
   daf_write_doubles(vrtbnorm_host,vrixnorm_host,vrminorm_host,
                     vrrgnorm_host,vrpknorm_host,
                     vltbnorm_host,vlixnorm_host,vlminorm_host,
                     vlrgnorm_host,vlpknorm_host);
   cout << "GPU norm : " << endl;
   daf_write_doubles(vrtbnorm_device,vrixnorm_device,vrminorm_device,
                     vrrgnorm_device,vrpknorm_device,
                     vltbnorm_device,vlixnorm_device,vlminorm_device,
                     vlrgnorm_device,vlpknorm_device);
   cout << "CPU norm after normalization : " << endl;
   daf_write_doubles(wrtbnorm_host,wrixnorm_host,wrminorm_host,
                     wrrgnorm_host,wrpknorm_host,
                     wltbnorm_host,wlixnorm_host,wlminorm_host,
                     wlrgnorm_host,wlpknorm_host);
   cout << "GPU norm after normalization : " << endl;
   daf_write_doubles(wrtbnorm_device,wrixnorm_device,wrminorm_device,
                     wrrgnorm_device,wrpknorm_device,
                     wltbnorm_device,wlixnorm_device,wlminorm_device,
                     wlrgnorm_device,wlpknorm_device);

   const double tol = 1.0e-140;
   double err = abs(vrtbnorm_device - vrtbnorm_host)
              + abs(wrtbnorm_device - wrtbnorm_host)
              + abs(vrixnorm_device - vrixnorm_host)
              + abs(wrixnorm_device - wrixnorm_host)
              + abs(vrminorm_device - vrminorm_host)
              + abs(wrminorm_device - wrminorm_host)
              + abs(vrrgnorm_device - vrrgnorm_host)
              + abs(wrrgnorm_device - wrrgnorm_host)
              + abs(vrpknorm_device - vrpknorm_host)
              + abs(wrpknorm_device - wrpknorm_host)
              + abs(vltbnorm_device - vltbnorm_host)
              + abs(wltbnorm_device - wltbnorm_host)
              + abs(vlixnorm_device - vlixnorm_host)
              + abs(wlixnorm_device - wlixnorm_host)
              + abs(vlminorm_device - vlminorm_host)
              + abs(wlminorm_device - wlminorm_host)
              + abs(vlrgnorm_device - vlrgnorm_host)
              + abs(wlrgnorm_device - wlrgnorm_host)
              + abs(vlpknorm_device - vlpknorm_host)
              + abs(wlpknorm_device - wlpknorm_host);

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
