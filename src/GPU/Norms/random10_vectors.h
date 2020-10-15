// The file random10_vectors.h specifies functions 
// to generate random vectors in deca double precision.

#ifndef __random10_vectors_h__
#define __random10_vectors_h__

void random_deca_double
 ( double *x_rtb, double *x_rix, double *x_rmi, double *x_rrg, double *x_rpk,
   double *x_ltb, double *x_lix, double *x_lmi, double *x_lrg, double *x_lpk );
/*
 * DESCRIPTION :
 *   Returns a random deca double x in [-1, +1],
 *   with the generation of ten random doubles
 *   so all ten parts of the random deca double are filled.
 *
 * ON RETURN :
 *   x_rtb     highest part of the random deca double x;
 *   x_rix     second highest part of the random deca double x;
 *   x_rmi     third highest part of the random deca double x;
 *   x_rrg     fourth highest part of the random deca double x;
 *   x_rpk     fifth highest part of the random deca double x;
 *   x_ltb     fifth lowest part of the random deca double x;
 *   x_lix     fourth lowest part of the random deca double x;
 *   x_lmi     third lowest part of the random deca double x;
 *   x_lrg     second lowest part of the random deca double x;
 *   x_lpk     lowest part of the random deca double x. */

void random_double10_vectors
 ( int dim, double *vrtb_host, double *vrix_host, double *vrmi_host,
   double *vrrg_host, double *vrpk_host, double *vltb_host,
   double *vlix_host, double *vlmi_host, double *vlrg_host,
   double *vlpk_host, double *vrtb_device, double *vrix_device,
   double *vrmi_device, double *vrrg_device, double *vrpk_device,
   double *vltb_device, double *vlix_device, double *vlmi_device,
   double *vlrg_device, double *vlpk_device );
/*
 * DESCRIPTION :
 *   Generates two instances of the same random deca double vector,
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
 *   vrtb_host holds as many randomly generated doubles as the value of dim;
 *   vrix_host holds as many randomly generated doubles as the value of dim;
 *   vrmi_host holds as many randomly generated doubles as the value of dim;
 *   vrrg_host holds as many randomly generated doubles as the value of dim;
 *   vrpk_host holds as many randomly generated doubles as the value of dim;
 *   vltb_host holds as many randomly generated doubles as the value of dim;
 *   vlix_host holds as many randomly generated doubles as the value of dim;
 *   vlmi_host holds as many randomly generated doubles as the value of dim;
 *   vlrg_host holds as many randomly generated doubles as the value of dim;
 *   vlpk_host holds as many randomly generated doubles as the value of dim;
 *   vrtb_device is the same vector as vtb_host;
 *   vrix_device is the same vector as vix_host;
 *   vrmi_device is the same vector as vmi_host;
 *   vrrg_device is the same vector as vrg_host;
 *   vrpk_device is the same vector as vpk_host;
 *   vltb_device is the same vector as vtb_host;
 *   vlix_device is the same vector as vix_host;
 *   vlmi_device is the same vector as vmi_host;
 *   vlrg_device is the same vector as vrg_host;
 *   vlpk_device is the same vector as vpk_host. */

void random_complex10_vectors
 ( int dim,
   double *vrertb_host, double *vrerix_host, double *vrermi_host,
   double *vrerrg_host, double *vrerpk_host,
   double *vreltb_host, double *vrelix_host, double *vrelmi_host,
   double *vrelrg_host, double *vrelpk_host,
   double *vimrtb_host, double *vimrix_host, double *vimrmi_host,
   double *vimrrg_host, double *vimrpk_host,
   double *vimltb_host, double *vimlix_host, double *vimlmi_host,
   double *vimlrg_host, double *vimlpk_host,
   double *vrertb_device, double *vrerix_device, double *vrermi_device,
   double *vrerrg_device, double *vrerpk_device,
   double *vreltb_device, double *vrelix_device, double *vrelmi_device,
   double *vrelrg_device, double *vrelpk_device,
   double *vimrtb_device, double *vimrix_device, double *vimrmi_device,
   double *vimrrg_device, double *vimrpk_device,
   double *vimltb_device, double *vimlix_device, double *vimlmi_device,
   double *vimlrg_device, double *vimlpk_device );
/*
 * DESCRIPTION :
 *   Generates two instances of the same random complex vector,
 *   one for the host and another for the device.
 *   The complex vector is represented by twenty real arrays,
 *   ten double arrays for the real parts and ten double arrays
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
 *   vrertb_host  highest real parts of the complex vectors;
 *   vrerix_host  second highest real parts of the complex vectors;
 *   vrermi_host  third highest real parts of the complex vectors;
 *   vrerrg_host  fourth highest real parts of the complex vectors;
 *   vrerpk_host  fifth highest real parts of the complex vectors;
 *   vreltb_host  fifth lowest real parts of the complex vectors;
 *   vrelix_host  fourth lowest real parts of the complex vectors;
 *   vrelmi_host  third lowest real parts of the complex vectors;
 *   vrelrg_host  second lowest real parts of the complex vectors;
 *   vrelpk_host  lowest real parts of the complex vectors;
 *   vimrtb_host  highest imaginary parts of the complex vectors;
 *   vimrix_host  second highest imaginary parts of the complex vectors;
 *   vimrmi_host  third highest imaginary parts of the complex vectors;
 *   vimrrg_host  fourth highest imaginary parts of the complex vectors;
 *   vimrpk_host  fifth highest imaginary parts of the complex vectors;
 *   vimltb_host  fifth lowest imaginary parts of the complex vectors;
 *   vimlix_host  fourth lowest imaginary parts of the complex vectors;
 *   vimlmi_host  third lowest imaginary parts of the complex vectors;
 *   vimlrg_host  second lowest imaginary parts of the complex vectors;
 *   vimlpk_host  lowest imaginary parts of the complex vectors;
 *   vrertb_device is the same vector as vrertb_host;
 *   vrerix_device is the same vector as vrerix_host;
 *   vrermi_device is the same vector as vrermi_host;
 *   vrerrg_device is the same vector as vrerrg_host;
 *   vrerpk_device is the same vector as vrerpk_host;
 *   vreltb_device is the same vector as vreltb_host;
 *   vrelix_device is the same vector as vrelix_host;
 *   vrelmi_device is the same vector as vrelmi_host;
 *   vrelrg_device is the same vector as vrelrg_host;
 *   vrelpk_device is the same vector as vrelpk_host;
 *   vimrtb_device is the same vector as vimrtb_host;
 *   vimrix_device is the same vector as vimrix_host;
 *   vimrmi_device is the same vector as vimrmi_host;
 *   vimrrg_device is the same vector as vimrrg_host;
 *   vimrpk_device is the same vector as vimrpk_host;
 *   vimltb_device is the same vector as vimltb_host;
 *   vimlix_device is the same vector as vimlix_host;
 *   vimlmi_device is the same vector as vimlmi_host;
 *   vimlrg_device is the same vector as vimlrg_host;
 *   vimlpk_device is the same vector as vimlpk_host. */

#endif
