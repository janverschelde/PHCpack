// The file random_vectors.cpp defines the code for the functions
// specified in random_numbers.h.

#include <cstdlib>
#include <cmath>
#include "random_numbers.h"

#define M_PI 3.141592653589793

double random_double ( void )
{
   double r = (double) rand();

   r = r/RAND_MAX;        // r is in [0,1]

   return (2.0*r - 1.0);  // result is in [-1, +1]
}

double random_angle ( void )
{
   double r = random_double();

   return (M_PI*r);
}
