// The file cmplx4_norm_host.cpp defines the code for the functions
// specified in cmplx4_norm_host.h.

#include "quad_double_functions.h"
#include "cmplx4_norm_host.h"

void make_copy
 ( int dim,
   double *orgrehihi, double *orgrelohi, double *orgrehilo, double *orgrelolo,
   double *orgimhihi, double *orgimlohi, double *orgimhilo, double *orgimlolo,
   double *duprehihi, double *duprelohi, double *duprehilo, double *duprelolo,
   double *dupimhihi, double *dupimlohi, double *dupimhilo, double *dupimlolo )
{
   for(int i=0; i<dim; i++)
   {
      duprehihi[i] = orgrehihi[i];
      duprelohi[i] = orgrelohi[i];
      duprehilo[i] = orgrehilo[i];
      duprelolo[i] = orgrelolo[i];
      dupimhihi[i] = orgimhihi[i];
      dupimlohi[i] = orgimlohi[i];
      dupimhilo[i] = orgimhilo[i];
      dupimlolo[i] = orgimlolo[i];
   }
}

void CPU_norm
 ( double *vrehihi, double *vrelohi, double *vrehilo, double *vrelolo,
   double *vimhihi, double *vimlohi, double *vimhilo, double *vimlolo,
   int dim,
   double *normhihi, double *normlohi, double *normhilo, double *normlolo )
{
   double sumhihi = 0.0;
   double sumlohi = 0.0;
   double sumhilo = 0.0;
   double sumlolo = 0.0;
   double prodhihi,prodlohi,prodhilo,prodlolo;

   for(int k=0; k<dim; k++) // sum = sum + vre[k]*vre[k] + vim[k]*vim[k]
   {
      qdf_sqr(vrehihi[k],vrelohi[k],vrehilo[k],vrelolo[k],
              &prodhihi,&prodlohi,&prodhilo,&prodlolo);
      qdf_inc(&sumhihi,&sumlohi,&sumhilo,&sumlolo,
              prodhihi,prodlohi,prodhilo,prodlolo);
      qdf_sqr(vimhihi[k],vimlohi[k],vimhilo[k],vimlolo[k],
              &prodhihi,&prodlohi,&prodhilo,&prodlolo);
      qdf_inc(&sumhihi,&sumlohi,&sumhilo,&sumlolo,
              prodhihi,prodlohi,prodhilo,prodlolo);
   }
   qdf_sqrt(sumhihi,sumlohi,sumhilo,sumlolo,
            normhihi,normlohi,normhilo,normlolo);
}

void CPU_normalize
 ( double *vrehihi, double *vrelohi, double *vrehilo, double *vrelolo,
   double *vimhihi, double *vimlohi, double *vimhilo, double *vimlolo,
   int dim,
   double normhihi, double normlohi, double normhilo, double normlolo )
{
   for(int i=0; i<dim; i++)
   {
      qdf_div(vrehihi[i],vrelohi[i],vrehilo[i],vrelolo[i],
              normhihi,normlohi,normhilo,normlolo,
              &vrehihi[i],&vrelohi[i],&vrehilo[i],&vrelolo[i]);
      qdf_div(vimhihi[i],vimlohi[i],vimhilo[i],vimlolo[i],
              normhihi,normlohi,normhilo,normlolo,
              &vimhihi[i],&vimlohi[i],&vimhilo[i],&vimlolo[i]);
   }
}
