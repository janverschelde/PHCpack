// The file dbl_norm_host.cpp defines the code for the functions
// specified in dbl_norm_host.h.

#include "dbl_norm_host.h"
#include <cmath>

void make_copy ( int dim, double* org, double* dup )
{
   for(int i=0; i<dim; i++) dup[i] = org[i];
}

void CPU_norm ( double* v, int dim, double* twonorm )
{
   double sum = 0.0;

   for(int i=0; i<dim; i++) sum = sum + v[i]*v[i];

   *twonorm = sqrt(sum);
}

void CPU_normalize ( double* v, int dim, double norm )
{
   for(int i=0; i<dim; i++) v[i] = v[i]/norm;
}
