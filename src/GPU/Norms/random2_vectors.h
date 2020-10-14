// The file random2_vectors.h specifies functions
// to generate random vectors in double double precision.

#ifndef __random2_vectors_h__
#define __random2_vectors_h__

void random_double_double ( double *x_hi, double *x_lo );
/*
 * DESCRIPTION :
 *   Returns a random double x in [-1, +1],
 *   with the generation of two random doubles
 *   so all three parts of the random double double are filled.
 *
 * ON RETURN :
 *   x_hi     highest part of the random double double x;
 *   x_lo     lowest part of the random double double x. */

void random_double2_vectors
 ( int dim, double *vhi_host, double *vlo_host,
            double *vhi_device, double *vlo_device );
/*
 * DESCRIPTION :
 *   Generates two instances of the same random double double vector,
 *   one for the host and another for the device.
 *
 * REQUIRED :
 *   Space has been allocated for vhi_host, vlo_host, v_hidevice, and
 *   vlo_device, for double arrays of size at least equal to dim.
 *
 * ON ENTRY :
 *   dim      dimension of the random real vector.
 *
 * ON RETURN :
 *   vhi_host holds as many randomly generated doubles as the value of dim;
 *   vlo_host holds as many randomly generated doubles as the value of dim;
 *   vhi_device is the same vector as vhi_host;
 *   vlo_device is the same vector as vlo_host. */

void random_complex2_vectors
 ( int dim, double *vrehi_host, double *vrelo_host,
            double *vimhi_host, double *vimlo_host,
   double *vrehi_device, double *vrelo_device,
   double *vimhi_device, double *vimlo_device );
/*
 * DESCRIPTION :
 *   Generates two instances of the same random complex vector,
 *   one for the host and another for the device.
 *   The complex vector is represented by four real vectors,
 *   with the high and low real parts and with high and low imaginary parts.
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
 *   vrelo_host   low real parts of the complex vectors;
 *   vimhi_host   high imaginary parts of the complex vectors;
 *   vimlo_host   low imaginary parts of the complex vectors;
 *   vrehi_device is the same vector as vrehi_host;
 *   vrelo_device is the same vector as vrelo_host;
 *   vimhi_device is the same vector as vimhi_host;
 *   vimlo_device is the same vector as vimlo_host. */

#endif
