// The file random8_vectors.h specifies functions to generate random vectors
// in octo double precision.

#ifndef __random8_vectors_h__
#define __random8_vectors_h__

void random_double8_vectors
 ( int dim, double *vhihihi_host, double *vlohihi_host,
            double *vhilohi_host, double *vlolohi_host,
            double *vhihilo_host, double *vlohilo_host,
            double *vhilolo_host, double *vlololo_host,
   double *vhihihi_device, double *vlohihi_device,
   double *vhilohi_device, double *vlolohi_device,
   double *vhihilo_device, double *vlohilo_device,
   double *vhilolo_device, double *vlololo_device );
/*
 * DESCRIPTION :
 *   Generates two instances of the same random octo double vector,
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
 *   vhihihi_host   as many randomly generated doubles as the value of dim;
 *   vlohihi_host   as many randomly generated doubles as the value of dim;
 *   vhilohi_host   as many randomly generated doubles as the value of dim;
 *   vlolohi_host   as many randomly generated doubles as the value of dim;
 *   vhihilo_host   as many randomly generated doubles as the value of dim;
 *   vlohilo_host   as many randomly generated doubles as the value of dim;
 *   vhilolo_host   as many randomly generated doubles as the value of dim;
 *   vlololo_host   as many randomly generated doubles as the value of dim;
 *   vhihihi_device is the same vector as vhihihi_host;
 *   vlohihi_device is the same vector as vlohihi_host;
 *   vhilohi_device is the same vector as vlohihi_host;
 *   vlolohi_device is the same vector as vlolohi_host;
 *   vhihilo_device is the same vector as vhihilo_host;
 *   vlohilo_device is the same vector as vlohilo_host;
 *   vhilolo_device is the same vector as vlohilo_host; 
 *   vlololo_device is the same vector as vlololo_host. */

void random_complex8_vectors
 ( int dim, double *vrehihihi_host, double *vrelohihi_host,
            double *vrehilohi_host, double *vrelolohi_host,
            double *vrehihilo_host, double *vrelohilo_host,
            double *vrehilolo_host, double *vrelololo_host,
            double *vimhihihi_host, double *vimlohihi_host,
            double *vimhilohi_host, double *vimlolohi_host,
            double *vimhihilo_host, double *vimlohilo_host,
            double *vimhilolo_host, double *vimlololo_host,
   double *vrehihihi_device, double *vrelohihi_device,
   double *vrehilohi_device, double *vrelolohi_device,
   double *vimhihihi_device, double *vimlohihi_device,
   double *vimhilohi_device, double *vimlolohi_device,
   double *vrehihilo_device, double *vrelohilo_device,
   double *vrehilolo_device, double *vrelololo_device,
   double *vimhihilo_device, double *vimlohilo_device,
   double *vimhilolo_device, double *vimlololo_device );
/*
 * DESCRIPTION :
 *   Generates two instances of the same random complex vector,
 *   one for the host and another for the device.
 *   The complex vector is represented by sixteen arrays,
 *   eight vectors for both real and imaginary parts of the complex numbers.
 *   both for real and imaginary parts of the complex numbers.
 *   The complex numbers are on the complex unit circle.
 *
 * REQUIRED :
 *   Space has been allocated for all arrays of size at least dim.
 *
 * ON ENTRY :
 *   dim        dimension of the random complex vector.
 *
 * ON RETURN :
 *   vrehihihi_host   highest real parts of the complex vectors;
 *   vrelohihi_host   second highest real parts of the complex vectors;
 *   vrehilohi_host   third highest real parts of the complex vectors;
 *   vrelolohi_host   fourth highest real parts of the complex vectors;
 *   vrehihilo_host   fourth lowest real parts of the complex vectors;
 *   vrehilolo_host   third lowest real parts of the complex vectors;
 *   vrelohilo_host   second lowest real parts of the complex vectors;
 *   vrelololo_host   lowest real parts of the complex vectors;
 *   vimhihihi_host   highest imaginary parts of the complex vectors;
 *   vimhilohi_host   second highest imaginary parts of the complex vectors;
 *   vimlohihi_host   third highest imaginary parts of the complex vectors;
 *   vimlolohi_host   fourth highest imaginary parts of the complex vectors;
 *   vimhihilo_host   fourth lowest imaginary parts of the complex vectors;
 *   vimlohilo_host   third lowest imaginary parts of the complex vectors;
 *   vimhilolo_host   second lowest imaginary parts of the complex vectors;
 *   vimlololo_host   lowest imaginary parts of the complex vectors;
 *   vrehihihi_device is the same vector as vrehihihi_host;
 *   vrelohihi_device is the same vector as vrelohihi_host;
 *   vrehilohi_device is the same vector as vrehilohi_host;
 *   vrelolohi_device is the same vector as vrelolohi_host;
 *   vrehihilo_device is the same vector as vrehihilo_host;
 *   vrelohilo_device is the same vector as vrelohilo_host;
 *   vrehilolo_device is the same vector as vrehilolo_host;
 *   vrelololo_device is the same vector as vrelololo_host;
 *   vimhihihi_device is the same vector as vimhihihi_host;
 *   vimlohihi_device is the same vector as vimlohihi_host;
 *   vimhilohi_device is the same vector as vimhilohi_host;
 *   vimlolohi_device is the same vector as vimlolohi_host;
 *   vimhihilo_device is the same vector as vimhihilo_host;
 *   vimlohilo_device is the same vector as vimlohilo_host;
 *   vimhilolo_device is the same vector as vimhilolo_host;
 *   vimlololo_device is the same vector as vimlolo_host. */

#endif
