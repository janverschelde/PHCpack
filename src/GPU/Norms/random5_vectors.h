// The file random5_vectors.h specifies functions 
// to generate random vectors in penta double precision.

#ifndef __random5_vectors_h__
#define __random5_vectors_h__

void random_penta_double
 ( double *x_tb, double *x_ix, double *x_mi, double *x_rg, double *x_pk );
/*
 * DESCRIPTION :
 *   Returns a random penta double x in [-1, +1],
 *   with the generation of five random doubles
 *   so all five parts of the random penta double are filled.
 *
 * ON RETURN :
 *   x_tb     highest part of the random penta double x;
 *   x_ix     second highest part of the random penta double x;
 *   x_mi     middle part of the random penta double x;
 *   x_rg     second lowest part of the random penta double x;
 *   x_pk     lowest part of the random penta double x. */

void random_double5_vectors
 ( int dim, double *vtb_host, double *vix_host, double *vmi_host,
   double *vrg_host, double *vpk_host,
   double *vtb_device, double *vix_device, double *vmi_device,
   double *vrg_device, double *vpk_device );
/*
 * DESCRIPTION :
 *   Generates two instances of the same random penta double vector,
 *   one for the host and another for the device.
 *
 * REQUIRED :
 *   Space has been allocated for all double arrays of size at least
 *   equal to the dimension dim.
 *
 * ON ENTRY :
 *   dim      dimension of the random real vector.
 *
 * ON RETURN :
 *   vtb_host holds as many randomly generated doubles as the value of dim;
 *   vix_host holds as many randomly generated doubles as the value of dim;
 *   vmi_host holds as many randomly generated doubles as the value of dim;
 *   vrg_host holds as many randomly generated doubles as the value of dim;
 *   vpk_host holds as many randomly generated doubles as the value of dim;
 *   vtb_device is the same vector as vtb_host;
 *   vix_device is the same vector as vix_host;
 *   vmi_device is the same vector as vmi_host;
 *   vrg_device is the same vector as vrg_host;
 *   vpk_device is the same vector as vpk_host. */

void random_complex5_vectors
 ( int dim, double *vretb_host, double *vreix_host, double *vremi_host,
            double *vrerg_host, double *vrepk_host,
            double *vimtb_host, double *vimix_host, double *vimmi_host,
            double *vimrg_host, double *vimpk_host,
   double *vretb_device, double *vreix_device, double *vremi_device,
   double *vrerg_device, double *vrepk_device,
   double *vimtb_device, double *vimix_device, double *vimmi_device,
   double *vimrg_device, double *vimpk_device );
/*
 * DESCRIPTION :
 *   Generates two instances of the same random complex vector,
 *   one for the host and another for the device.
 *   The complex vector is represented by ten real arrays,
 *   five double arrays for the real parts and five double arrays
 *   for the imaginary parts of the complex numbers.
 *   The complex numbers are on the complex unit circle.
 *
 * REQUIRED :
 *   Space has been allocated for all arrays of size at least dim.
 *
 * ON ENTRY :
 *   dim          dimension of the random complex vector.
 *
 * ON RETURN :
 *   vretb_host   highest real parts of the complex vectors;
 *   vreix_host   second highest real parts of the complex vectors;
 *   vremi_host   middle real parts of the complex vectors;
 *   vrerg_host   second lowest real parts of the complex vectors;
 *   vrepk_host   lowest real parts of the complex vectors;
 *   vimtb_host   highest imaginary parts of the complex vectors;
 *   vimix_host   second highest imaginary parts of the complex vectors;
 *   vimmi_host   middle imaginary parts of the complex vectors;
 *   vimrg_host   second lowest imaginary parts of the complex vectors;
 *   vimpk_host   lowest imaginary parts of the complex vectors;
 *   vretb_device is the same vector as vretb_host;
 *   vreix_device is the same vector as vreix_host;
 *   vremi_device is the same vector as vremi_host;
 *   vrerg_device is the same vector as vrerg_host;
 *   vrepk_device is the same vector as vrepk_host;
 *   vimtb_device is the same vector as vimtb_host;
 *   vimix_device is the same vector as vimix_host;
 *   vimmi_device is the same vector as vimmi_host;
 *   vimrg_device is the same vector as vimrg_host;
 *   vimpk_device is the same vector as vimpk_host. */

#endif
