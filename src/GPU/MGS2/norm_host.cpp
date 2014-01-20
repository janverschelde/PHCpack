// defines code of the norm computations for execution on the host

#include "norm_host.h"

void make_copy ( int dim, complex<T>* org, complex<T>* dup )
{
   for(int i=0; i<dim; i++) dup[i] = org[i];
}

void make_copyH ( int dim, complexH<T1>* org, complexH<T1>* dup )
{
   for(int i=0; i<dim; i++) dup[i] = org[i];
}

void CPU_norm ( complexH<T1>* v, int dim, T1* twonorm )
{
   T1 sum(0.0);

   for(int i=0; i<dim; i++)
      sum = sum + v[i].real*v[i].real + v[i].imag*v[i].imag;

   *twonorm = sqrt(sum);
}

void CPU_normalize ( complexH<T1>* v, int dim, T1 norm )
{
   for(int i=0; i<dim; i++) v[i] = v[i]/norm;
}
