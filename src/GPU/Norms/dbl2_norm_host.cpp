// The file dbl2_norm_host.cpp defines the code for the functions
// specified in dbl2_norm_host.h.

#include "dbl2_norm_host.h"
#include "double_double_functions.h"
#include <cmath>

void make_copy
 ( int dim, double* orghi, double* orglo, double* duphi, double* duplo )
{
   for(int i=0; i<dim; i++)
   {
      duphi[i] = orghi[i];
      duplo[i] = orglo[i];
   }
}

void CPU_norm
 ( double* vhi, double* vlo, int dim, double* normhi, double* normlo )
{
   double sumhi,sumlo = 0.0;
   double prodhi,prodlo;

   for(int i=0; i<dim; i++) // sum = sum + v[i]*v[i]
   {
      dd_sqr(vhi[i],vlo[i],&prodhi,&prodlo);
      dd_add(prodhi,prodlo,sumhi,sumlo,&sumhi,&sumlo);
   }
   dd_sqrt(sumhi,sumlo,normhi,normlo);
}

void CPU_normalize
 ( double* vhi, double* vlo, int dim, double normhi, double normlo )
{
   for(int i=0; i<dim; i++)
      dd_div(vhi[i],vlo[i],normhi,normlo,&vhi[i],&vlo[i]);
}
