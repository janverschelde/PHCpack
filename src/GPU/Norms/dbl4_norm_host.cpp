// The file dbl4_norm_host.cpp defines the code for the functions
// specified in dbl4_norm_host.h.

#include "dbl4_norm_host.h"
#include "quad_double_functions.h"

void make_copy
 ( int dim,
   double *orghihi, double *orglohi, double *orghilo, double *orglolo,
   double *duphihi, double *duplohi, double *duphilo, double *duplolo)
{
   for(int i=0; i<dim; i++)
   {
      duphihi[i] = orghihi[i];
      duplohi[i] = orglohi[i];
      duphilo[i] = orghilo[i];
      duplolo[i] = orglolo[i];
   }
}

void CPU_norm
 ( double *vhihi, double *vlohi, double *vhilo, double *vlolo, int dim,
   double *normhihi, double *normlohi, double *normhilo, double *normlolo )
{
   double sumhihi = 0.0;
   double sumlohi = 0.0;
   double sumhilo = 0.0;
   double sumlolo = 0.0;
   double prodhihi,prodlohi,prodhilo,prodlolo;

   for(int i=0; i<dim; i++) // sum = sum + v[i]*v[i]
   {
      qdf_sqr(vhihi[i],vlohi[i],vhilo[i],vlolo[i],
              &prodhihi,&prodlohi,&prodhilo,&prodlolo);
      qdf_inc(&sumhihi,&sumlohi,&sumhilo,&sumlolo,
              prodhihi,prodlohi,prodhilo,prodlolo);
   }
   qdf_sqrt(sumhihi,sumlohi,sumhilo,sumlolo,
            normhihi,normlohi,normhilo,normlolo);
}

void CPU_normalize
 ( double *vhihi, double *vlohi, double *vhilo, double *vlolo, int dim,
   double normhihi, double normlohi, double normhilo, double normlolo )
{
   for(int i=0; i<dim; i++)
      qdf_div(vhihi[i],vlohi[i],vhilo[i],vlolo[i],
              normhihi,normlohi,normhilo,normlolo,
              &vhihi[i],&vlohi[i],&vhilo[i],&vlolo[i]);
}
