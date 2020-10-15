// The file cmplx2_norm_host.cpp defines the code for the functions
// specified in cmplx2_norm_host.h.

#include "double_double_functions.h"
#include "cmplx2_norm_host.h"

void make_copy
 ( int dim,
   double *orgrehi, double *orgrelo, double *orgimhi, double *orgimlo,
   double *duprehi, double *duprelo, double *dupimhi, double *dupimlo )
{
   for(int i=0; i<dim; i++)
   {
      duprehi[i] = orgrehi[i];
      duprelo[i] = orgrelo[i];
      dupimhi[i] = orgimhi[i];
      dupimlo[i] = orgimlo[i];
   }
}

void CPU_norm
 ( double *vrehi, double *vrelo, double *vimhi, double *vimlo,
   int dim, double *normhi, double *normlo )
{
   double sumhi = 0.0;
   double sumlo = 0.0;
   double prodhi,prodlo;

   for(int k=0; k<dim; k++) // sum = sum + vre[k]*vre[k] + vim[k]*vim[k]
   {
      ddf_sqr(vrehi[k],vrelo[k],&prodhi,&prodlo);
      ddf_inc(&sumhi,&sumlo,prodhi,prodlo);
      ddf_sqr(vimhi[k],vimlo[k],&prodhi,&prodlo);
      ddf_inc(&sumhi,&sumlo,prodhi,prodlo);
   }
   ddf_sqrt(sumhi,sumlo,normhi,normlo);
}

void CPU_normalize
 ( double *vrehi, double *vrelo, double *vimhi, double *vimlo,
   int dim, double normhi, double normlo )
{
   for(int i=0; i<dim; i++)
   {
      ddf_div(vrehi[i],vrelo[i],normhi,normlo,&vrehi[i],&vrelo[i]);
      ddf_div(vimhi[i],vimlo[i],normhi,normlo,&vimhi[i],&vimlo[i]);
   }
}
