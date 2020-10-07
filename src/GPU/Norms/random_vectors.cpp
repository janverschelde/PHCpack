// The file random_vectors.cpp defines the code for the functions
// specified in random_vectors.h.

#include <cstdlib>
#include <cmath>
#include "random_numbers.h"
#include "random_vectors.h"

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

void random_complex_vectors
 ( int dim, double* vre_host, double* vim_host,
   double* vre_device, double* vim_device )
{
   double r,cr,sr;

   for(int k=0; k<dim; k++)
   {
      r = random_angle(); cr = cos(r); sr = sin(r);
      vre_host[k] = cr; vre_device[k] = cr;
      vim_host[k] = sr; vim_device[k] = sr;
   }
}
