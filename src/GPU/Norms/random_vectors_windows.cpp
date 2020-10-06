// The file random_vectors.cpp defines the code for the functions
// specified in random_vectors.h.

#include <cstdlib>
#include <corecrt_math_defines.h>
#include "random_vectors.h"

double random_double ( void )
{
   double r = (double) rand();

   r = r/RAND_MAX;        // r is in [0,1]

   return (2.0*r - 1.0);  // result is in [-1, +1]
}

double random_angle ( void )
{
   double r = (double) rand();

   r = r/RAND_MAX;

   return (2.0*M_PI*r);
}

void random_double_vectors ( int dim, double* v_host, double* v_device )
{
   double r;

   for(int k=0; k<dim; k++)
   {
      r = random_double();
      v_host[k] = r;
      v_device[k] = r;
   }
}

void random_double2_vectors
 ( int dim, double* vhi_host, double* vlo_host,
   double* vhi_device, double* vlo_device )
{
   double r;

   for(int k=0; k<dim; k++)
   {
      r = random_double();
      vhi_host[k] = r;
      vhi_device[k] = r;
      vlo_host[k] = 0.0;
      vlo_device[k] = 0.0;
   }
}

void random_complex_vectors
 ( int dim, double* vre_host, double* vim_host,
   double* vre_device, double* vim_device )
{
   double r,cr,sr;

   for(int k=0; k<dim; k++)
   {
      r = random_angle();
      cr = cos(r); sr = sin(r);
      vre_host[k] = cr; vre_device[k] = cr;
      vim_host[k] = sr; vim_device[k] = sr;
   }
}
