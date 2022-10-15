// The file dbl_monomial_systems.cpp defines the functions specified in
// the file dbl_monomial_systems.h.

#include <iostream>
#include <cmath>
#include "random_numbers.h"
#include "random_series.h"
#include "dbl_monomial_systems.h"

using namespace std;

void make_complex_exponentials
 ( int dim, int deg, double *angles, double **sre, double **sim )
{
   double xre,xim;

   for(int i=0; i<dim; i++)
   {
      angles[i] = random_angle();
      xre = cos(angles[i]);
      xim = sin(angles[i]);

      cmplx_exponential(deg,xre,xim,sre[i],sim[i]);
   }
}
