// The file random3_vectors.h specifies functions
// to generate random vectors in triple double precision.

#ifndef __random3_vectors_h__
#define __random3_vectors_h__

void random_triple_double
 ( double *x_hi, double *x_mi, double *x_lo );
/*
 * DESCRIPTION :
 *   Returns a random triple double x in [-1, +1],
 *   with the generation of three random doubles
 *   so all three parts of the random triple double are filled.
 *
 * ON RETURN :
 *   x_hi     highest part of the random triple double x;
 *   x_mi     middle part of the random triple double x;
 *   x_lo     lowest part of the random triple double x. */

void random_double3_vectors
 ( int dim, double *vhi_host, double *vmi_host, double *vlo_host,
   double *vhi_device, double *vmi_device, double *vlo_device );
/*
 * DESCRIPTION :
 *   Generates two instances of the same random triple double vector,
 *   one for the host and another for the device.
 *
 * REQUIRED :
 *   Space has been allocated for vhi_host, vmi_host, vlo_host, v_hidevice,
 *   vmi_device, and vlo_device, for double arrays of size at least equal 
 *   to the dimension dim.
 *
 * ON ENTRY :
 *   dim      dimension of the random real vector.
 *
 * ON RETURN :
 *   vhi_host holds as many randomly generated doubles as the value of dim;
 *   vmi_host holds as many randomly generated doubles as the value of dim;
 *   vlo_host holds as many randomly generated doubles as the value of dim;
 *   vhi_device is the same vector as vhi_host;
 *   vmi_device is the same vector as vmi_host;
 *   vlo_device is the same vector as vlo_host. */

void random_complex3_vectors
 ( int dim, double *vrehi_host, double *vremi_host, double *vrelo_host,
            double *vimhi_host, double *vimmi_host, double *vimlo_host,
   double *vrehi_device, double *vremi_device, double *vrelo_device,
   double *vimhi_device, double *vimmi_device, double *vimlo_device );
/*
 * DESCRIPTION :
 *   Generates two instances of the same random complex vector,
 *   one for the host and another for the device.
 *   The complex vector is represented by six real vectors,
 *   with the high, middle, and low real parts and with high, middle,
 *   and low imaginary parts.
 *   The complex numbers are on the complex unit circle.
 *
 * REQUIRED :
 *   Space has been allocated for all arrays of size at least dim.
 *
 * ON ENTRY :
 *   dim          dimension of the random complex vector.
 *
 * ON RETURN :
 *   vrehi_host   high real parts of the complex vectors;
 *   vremi_host   middle real parts of the complex vectors;
 *   vrelo_host   low real parts of the complex vectors;
 *   vimhi_host   high imaginary parts of the complex vectors;
 *   vimmi_host   middle imaginary parts of the complex vectors;
 *   vimlo_host   low imaginary parts of the complex vectors;
 *   vrehi_device is the same vector as vrehi_host;
 *   vremi_device is the same vector as vremi_host;
 *   vrelo_device is the same vector as vrelo_host;
 *   vimhi_device is the same vector as vimhi_host;
 *   vimmi_device is the same vector as vimmi_host;
 *   vimlo_device is the same vector as vimlo_host. */

#endif
