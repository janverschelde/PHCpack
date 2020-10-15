// The file cmplx3_norm_host.cpp defines the code for the functions
// specified in cmplx3_norm_host.h.

#include "triple_double_functions.h"
#include "cmplx3_norm_host.h"

void make_copy
 ( int dim,
   double *orgrehi, double *orgremi, double *orgrelo,
   double *orgimhi, double *orgimmi, double *orgimlo,
   double *duprehi, double *dupremi, double *duprelo,
   double *dupimhi, double *dupimmi, double *dupimlo )
{
   for(int i=0; i<dim; i++)
   {
      duprehi[i] = orgrehi[i];
      dupremi[i] = orgremi[i];
      duprelo[i] = orgrelo[i];
      dupimhi[i] = orgimhi[i];
      dupimmi[i] = orgimmi[i];
      dupimlo[i] = orgimlo[i];
   }
}

void CPU_norm
 ( double *vrehi, double *vremi, double *vrelo,
   double *vimhi, double *vimmi, double *vimlo, int dim,
   double *normhi, double *normmi, double *normlo )
{
   double sumhi = 0.0;
   double summi = 0.0;
   double sumlo = 0.0;
   double prodhi,prodmi,prodlo;

   for(int k=0; k<dim; k++) // sum = sum + vre[k]*vre[k] + vim[k]*vim[k]
   {
      tdf_sqr(vrehi[k],vremi[k],vrelo[k],&prodhi,&prodmi,&prodlo);
      tdf_inc(&sumhi,&summi,&sumlo,prodhi,prodmi,prodlo);
      tdf_sqr(vimhi[k],vimmi[k],vimlo[k],&prodhi,&prodmi,&prodlo);
      tdf_inc(&sumhi,&summi,&sumlo,prodhi,prodmi,prodlo);
   }
   tdf_sqrt(sumhi,summi,sumlo,normhi,normmi,normlo);
}

void CPU_normalize
 ( double *vrehi, double *vremi, double *vrelo,
   double *vimhi, double *vimmi, double *vimlo, int dim,
   double normhi, double normmi, double normlo )
{
   for(int i=0; i<dim; i++)
   {
      tdf_div(vrehi[i],vremi[i],vrelo[i],normhi,normmi,normlo,
              &vrehi[i],&vremi[i],&vrelo[i]);
      tdf_div(vimhi[i],vimmi[i],vimlo[i],normhi,normmi,normlo,
              &vimhi[i],&vimmi[i],&vimlo[i]);
   }
}
