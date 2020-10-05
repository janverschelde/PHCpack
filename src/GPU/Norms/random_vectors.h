// The file random_vectors.h specifies functions to generate random vectors.

double random_double ( void );
// Returns a random double in [-1,+1].

double random_angle ( void );
// Returns a random double in [0, 2*pi].

void random_double_vectors ( int dim, double* v_host, double* v_device );
/*
   DESCRIPTION :
     Generates two instances of the same random vector,
     one for the host and another for the device.

   REQUIRED :
     Space has been allocated for both v_host and v_device,
     for double arrays of size at least equal to dim.

   ON ENTRY :
     dim      dimension of the random real vector.

   ON RETURN :
     v_host   as many randomly generated doubles as the value of dim;
     v_device is the same vector as v_host.                            */

void random_double2_vectors
 ( int dim, double* vhi_host, double* vlo_host,
   double* vhi_device, double* vlo_device );
/*
   DESCRIPTION :
     Generates two instances of the same random double double vector,
     one for the host and another for the device.

   REQUIRED :
     Space has been allocated for vhi_host, vlo_host, v_hidevice, and
     vlo_device, for double arrays of size at least equal to dim.

   ON ENTRY :
     dim      dimension of the random real vector.

   ON RETURN :
     vhi_host as many randomly generated doubles as the value of dim;
     vlo_host as many randomly generated doubles as the value of dim;
     vhi_device is the same vector as vhi_host;
     vlo_device is the same vector as vlo_host.                            */

void random_complex_vectors
 ( int dim, double* vre_host, double* vim_host,
   double* vre_device, double* vim_device );
/*
   DESCRIPTION :
     Generates two instances of the same random complex vector,
     one for the host and another for the device.
     The complex vector is represented by two real vectors,
     one with the real parts and another with the imaginary parts.
     The complex numbers are on the complex unit circle.

   REQUIRED :
     Space has been allocated for both vre_host, vim_host, vre_device,
     and vim_device for double arrays of size at least equal to dim.

   ON ENTRY :
     dim        dimension of the random complex vector.

   ON RETURN :
     vre_host   real parts of the complex vectors;
     vim_host   imaginary parts of the complex vectors;
     vre_device is the same vector as vre_host;
     vim_device is the same vector as vim_host.                            */
