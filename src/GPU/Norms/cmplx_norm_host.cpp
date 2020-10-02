// The file cmplx_norm_host.cpp defines the code for the functions
// specified in cmplx_norm_host.h.

#include <cmath>
#include "cmplx_norm_host.h"

void make_copy
 ( int dim, double* orgre, double* orgim, double* dupre, double* dupim )
{
   for(int i=0; i<dim; i++)
   {
      dupre[i] = orgre[i];
      dupim[i] = orgim[i];
   }
}

void CPU_norm ( double* vre, double* vim, int dim, double* twonorm )
{
   double sum = 0.0;

   for(int k=0; k<dim; k++) sum = sum + vre[k]*vre[k] + vim[k]*vim[k];

   *twonorm = sqrt(sum);
}

void CPU_normalize ( double* vre, double* vim, int dim, double norm )
{
   for(int i=0; i<dim; i++)
   {
      vre[i] = vre[i]/norm;
      vim[i] = vim[i]/norm;
   }
}
