// The file random4_vectors.h specifies functions
// to generate random vectors in quad double precision.

#ifndef __random4_vectors_h__
#define __random4_vectors_h__

void random_quad_double
 ( double *x_hihi, double *x_lohi, double *x_hilo, double *x_lolo );
/*
 * DESCRIPTION :
 *   Returns a random quad double x in [-1, +1],
 *   with the generation of four random doubles
 *   so all four parts of the random quad double are filled.
 *
 * ON RETURN :
 *   x_hihi   highest part of the random quad double x;
 *   x_lohi   second highest part of the random quad double x;
 *   x_hilo   second lowest part of the random quad double x;
 *   x_lolo   lowest part of the random quad double x. */

void random_double4_vectors
 ( int dim, double *vhihi_host, double *vlohi_host,
            double *vhilo_host, double *vlolo_host,
   double *vhihi_device, double *vlohi_device,
   double *vhilo_device, double *vlolo_device );
/*
 * DESCRIPTION :
 *   Generates two instances of the same random quad double vector,
 *   one for the host and another for the device.
 *
 * REQUIRED :
 *   Space has been allocated for all double arrays of size
 *   at least equal to the dimension dim.
 *
 * ON ENTRY :
 *   dim      dimension of the random real vector.
 *
 * ON RETURN :
 *   vhihi_host holds as many randomly generated doubles as the value of dim;
 *   vlohi_host holds as many randomly generated doubles as the value of dim;
 *   vhilo_host holds as many randomly generated doubles as the value of dim;
 *   vlolo_host holds as many randomly generated doubles as the value of dim;
 *   vhihi_device is the same vector as vhihi_host;
 *   vlohi_device is the same vector as vlohi_host;
 *   vhilo_device is the same vector as vlohi_host;
 *   vlolo_device is the same vector as vlolo_host. */

void random_complex4_vectors
 ( int dim, double *vrehihi_host, double *vrelohi_host,
            double *vrehilo_host, double *vrelolo_host,
            double *vimhihi_host, double *vimlohi_host,
            double *vimhilo_host, double *vimlolo_host,
   double *vrehihi_device, double *vrelohi_device,
   double *vrehilo_device, double *vrelolo_device,
   double *vimhihi_device, double *vimlohi_device,
   double *vimhilo_device, double *vimlolo_device );
/*
 * DESCRIPTION :
 *   Generates two instances of the same random complex vector,
 *   one for the host and another for the device.
 *   The complex vector is represented by eight arrays,
 *   with the highest, second highest, second lowest, and lowest doubles
 *   both for real and imaginary parts of the complex numbers.
 *   The complex numbers are on the complex unit circle.
 *
 * REQUIRED :
 *   Space has been allocated for all arrays of size at least dim.
 *
 * ON ENTRY :
 *   dim            dimension of the random complex vector.
 *
 * ON RETURN :
 *   vrehihi_host   highest real parts of the complex vectors;
 *   vrelohi_host   second highest real parts of the complex vectors;
 *   vrehilo_host   second lowest real parts of the complex vectors;
 *   vrelolo_host   lowest real parts of the complex vectors;
 *   vimhihi_host   highest imaginary parts of the complex vectors;
 *   vimlohi_host   second highest imaginary parts of the complex vectors;
 *   vimhilo_host   second lowest imaginary parts of the complex vectors;
 *   vimlolo_host   lowest imaginary parts of the complex vectors;
 *   vrehihi_device is the same vector as vrehihi_host;
 *   vrelohi_device is the same vector as vrelohi_host;
 *   vrehilo_device is the same vector as vrehilo_host;
 *   vrelolo_device is the same vector as vrelolo_host;
 *   vimhihi_device is the same vector as vimhihi_host;
 *   vimlohi_device is the same vector as vimlohi_host;
 *   vimhilo_device is the same vector as vimhilo_host;
 *   vimlolo_device is the same vector as vimlolo_host. */

#endif
