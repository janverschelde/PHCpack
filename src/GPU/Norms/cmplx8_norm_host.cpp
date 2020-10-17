// The file cmplx8_norm_host.cpp defines the code for the functions
// specified in cmplx8_norm_host.h.

#include "octo_double_functions.h"
#include "cmplx8_norm_host.h"

void make_copy
 ( int dim,
   double *orgrehihihi, double *orgrelohihi, double *orgrehilohi,
   double *orgrelolohi, double *orgrehihilo, double *orgrelohilo,
   double *orgrehilolo, double *orgrelololo,
   double *orgimhihihi, double *orgimlohihi, double *orgimhilohi,
   double *orgimlolohi, double *orgimhihilo, double *orgimlohilo,
   double *orgimhilolo, double *orgimlololo,
   double *duprehihihi, double *duprelohihi, double *duprehilohi,
   double *duprelolohi, double *duprehihilo, double *duprelohilo,
   double *duprehilolo, double *duprelololo,
   double *dupimhihihi, double *dupimlohihi, double *dupimhilohi,
   double *dupimlolohi, double *dupimhihilo, double *dupimlohilo,
   double *dupimhilolo, double *dupimlololo )
{
   for(int i=0; i<dim; i++)
   {
      duprehihihi[i] = orgrehihihi[i];
      duprelohihi[i] = orgrelohihi[i];
      duprehilohi[i] = orgrehilohi[i];
      duprelolohi[i] = orgrelolohi[i];
      duprehihilo[i] = orgrehihilo[i];
      duprelohilo[i] = orgrelohilo[i];
      duprehilolo[i] = orgrehilolo[i];
      duprelololo[i] = orgrelololo[i];
      dupimhihihi[i] = orgimhihihi[i];
      dupimlohihi[i] = orgimlohihi[i];
      dupimhilohi[i] = orgimhilohi[i];
      dupimlolohi[i] = orgimlolohi[i];
      dupimhihilo[i] = orgimhihilo[i];
      dupimlohilo[i] = orgimlohilo[i];
      dupimhilolo[i] = orgimhilolo[i];
      dupimlololo[i] = orgimlololo[i];
   }
}

void CPU_norm
 ( double *vrehihihi, double *vrelohihi, double *vrehilohi,
   double *vrelolohi, double *vrehihilo, double *vrelohilo,
   double *vrehilolo, double *vrelololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi,
   double *vimlolohi, double *vimhihilo, double *vimlohilo,
   double *vimhilolo, double *vimlololo,
   int dim,
   double *normhihihi, double *normlohihi, double *normhilohi,
   double *normlolohi, double *normhihilo, double *normlohilo,
   double *normhilolo, double *normlololo )
{
   double sumhihihi = 0.0;
   double sumlohihi = 0.0;
   double sumhilohi = 0.0;
   double sumlolohi = 0.0;
   double sumhihilo = 0.0;
   double sumlohilo = 0.0;
   double sumhilolo = 0.0;
   double sumlololo = 0.0;
   double prodhihihi,prodlohihi,prodhilohi,prodlolohi;
   double prodhihilo,prodlohilo,prodhilolo,prodlololo;

   for(int k=0; k<dim; k++) // sum = sum + vre[k]*vre[k] + vim[k]*vim[k]
   {
      odf_sqr(vrehihihi[k],vrelohihi[k],vrehilohi[k],vrelolohi[k],
              vrehihilo[k],vrelohilo[k],vrehilolo[k],vrelololo[k],
              &prodhihihi,&prodlohihi,&prodhilohi,&prodlolohi,
              &prodhihilo,&prodlohilo,&prodhilolo,&prodlololo);
      odf_inc(&sumhihihi,&sumlohihi,&sumhilohi,&sumlolohi,
              &sumhihilo,&sumlohilo,&sumhilolo,&sumlololo,
              prodhihihi,prodlohihi,prodhilohi,prodlolohi,
              prodhihilo,prodlohilo,prodhilolo,prodlololo);
      odf_sqr(vimhihihi[k],vimlohihi[k],vimhilohi[k],vimlolohi[k],
              vimhihilo[k],vimlohilo[k],vimhilolo[k],vimlololo[k],
              &prodhihihi,&prodlohihi,&prodhilohi,&prodlolohi,
              &prodhihilo,&prodlohilo,&prodhilolo,&prodlololo);
      odf_inc(&sumhihihi,&sumlohihi,&sumhilohi,&sumlolohi,
              &sumhihilo,&sumlohilo,&sumhilolo,&sumlololo,
              prodhihihi,prodlohihi,prodhilohi,prodlolohi,
              prodhihilo,prodlohilo,prodhilolo,prodlololo);
   }
   odf_sqrt(sumhihihi,sumlohihi,sumhilohi,sumlolohi,
            sumhihilo,sumlohilo,sumhilolo,sumlololo,
            normhihihi,normlohihi,normhilohi,normlolohi,
            normhihilo,normlohilo,normhilolo,normlololo);
}

void CPU_normalize
 ( double *vrehihihi, double *vrelohihi, double *vrehilohi,
   double *vrelolohi, double *vrehihilo, double *vrelohilo,
   double *vrehilolo, double *vrelololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi,
   double *vimlolohi, double *vimhihilo, double *vimlohilo,
   double *vimhilolo, double *vimlololo,
   int dim,
   double normhihihi, double normlohihi, double normhilohi,
   double normlolohi, double normhihilo, double normlohilo,
   double normhilolo, double normlololo )
{
   for(int i=0; i<dim; i++)
   {
      odf_div(vrehihihi[i],vrelohihi[i],vrehilohi[i],vrelolohi[i],
              vrehihilo[i],vrelohilo[i],vrehilolo[i],vrelololo[i],
              normhihihi,normlohihi,normhilohi,normlolohi,
              normhihilo,normlohilo,normhilolo,normlololo,
              &vrehihihi[i],&vrelohihi[i],&vrehilohi[i],&vrelolohi[i],
              &vrehihilo[i],&vrelohilo[i],&vrehilolo[i],&vrelololo[i]);
      odf_div(vimhihihi[i],vimlohihi[i],vimhilohi[i],vimlolohi[i],
              vimhihilo[i],vimlohilo[i],vimhilolo[i],vimlololo[i],
              normhihihi,normlohihi,normhilohi,normlolohi,
              normhihilo,normlohilo,normhilolo,normlololo,
              &vimhihihi[i],&vimlohihi[i],&vimhilohi[i],&vimlolohi[i],
              &vimhihilo[i],&vimlohilo[i],&vimhilolo[i],&vimlololo[i]);
   }
}
