// The file dbl_onenorms_host.cpp defines functions with prototypes in
// the file dbl_onenorms_host.h.

#include <cmath>
#include "dbl_onenorms_host.h"

using namespace std;

void CPU_dbl_onenorm ( int dim, double *v, double *nrm )
{
   double result = abs(v[0]);

   for(int i=1; i<dim; i++)
      if(abs(v[i]) > result) result = abs(v[i]);

   *nrm = result;
}

void CPU_cmplx_onenorm ( int dim, double *vre, double *vim, double *nrm )
{
   double result = abs(vre[0]) + abs(vim[0]);

   for(int i=1; i<dim;i ++)
      if(abs(vre[i]) + abs(vim[i]) > result)
         result = abs(vre[i]) + abs(vim[i]);
 
   *nrm = result;
}
