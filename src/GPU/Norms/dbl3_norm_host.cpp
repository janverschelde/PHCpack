// The file dbl3_norm_host.cpp defines the code for the functions
// specified in dbl3_norm_host.h.

#include "dbl3_norm_host.h"
#include "triple_double_functions.h"

void make_copy
 ( int dim, double *orghi, double *orgmi, double *orglo,
            double *duphi, double *dupmi, double* duplo )
{
   for(int i=0; i<dim; i++)
   {
      duphi[i] = orghi[i];
      dupmi[i] = orgmi[i];
      duplo[i] = orglo[i];
   }
}

void CPU_norm
 ( double *vhi, double *vmi, double *vlo, int dim,
   double *normhi, double *normmi, double *normlo )
{
   double sumhi = 0.0;
   double summi = 0.0;
   double sumlo = 0.0;
   double prodhi,prodmi,prodlo;

   for(int i=0; i<dim; i++) // sum = sum + v[i]*v[i]
   {
      tdf_sqr(vhi[i],vmi[i],vlo[i],&prodhi,&prodmi,&prodlo);
      tdf_inc(&sumhi,&summi,&sumlo,prodhi,prodmi,prodlo);
   }
   tdf_sqrt(sumhi,summi,sumlo,normhi,normmi,normlo);
}

void CPU_normalize
 ( double *vhi, double *vmi, double *vlo, int dim,
   double normhi, double normmi, double normlo )
{
   for(int i=0; i<dim; i++)
      tdf_div(vhi[i],vmi[i],vlo[i],normhi,normmi,normlo,
              &vhi[i],&vmi[i],&vlo[i]);
}
