// The file dbl8_norm_host.cpp defines the code for the functions
// specified in dbl8_norm_host.h.

#include "dbl8_norm_host.h"
#include "octo_double_functions.h"

void make_copy
 ( int dim,
   double *orghihihi, double *orglohihi, double *orghilohi, double *orglolohi,
   double *orghihilo, double *orglohilo, double *orghilolo, double *orglololo,
   double *duphihihi, double *duplohihi, double *duphilohi, double *duplolohi,
   double *duphihilo, double *duplohilo, double *duphilolo, double *duplololo)
{
   for(int i=0; i<dim; i++)
   {
      duphihihi[i] = orghihihi[i];
      duplohihi[i] = orglohihi[i];
      duphilohi[i] = orghilohi[i];
      duplolohi[i] = orglolohi[i];
      duphihilo[i] = orghihilo[i];
      duplohilo[i] = orglohilo[i];
      duphilolo[i] = orghilolo[i];
      duplololo[i] = orglololo[i];
   }
}

void CPU_norm
 ( double *vhihihi, double *vlohihi, double *vhilohi, double *vlolohi,
   double *vhihilo, double *vlohilo, double *vhilolo, double *vlololo,
   int dim, double *normhihihi, double *normlohihi, double *normhilohi,
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

   for(int i=0; i<dim; i++) // sum = sum + v[i]*v[i]
   {
      odf_sqr(vhihihi[i],vlohihi[i],vhilohi[i],vlolohi[i],
              vhihilo[i],vlohilo[i],vhilolo[i],vlololo[i],
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
 ( double *vhihihi, double *vlohihi, double *vhilohi, double *vlolohi,
   double *vhihilo, double *vlohilo, double *vhilolo, double *vlololo,
   int dim, double normhihihi, double normlohihi, double normhilohi,
   double normlolohi, double normhihilo, double normlohilo,
   double normhilolo, double normlololo )
{
   for(int i=0; i<dim; i++)
      odf_div(vhihihi[i],vlohihi[i],vhilohi[i],vlolohi[i],
              vhihilo[i],vlohilo[i],vhilolo[i],vlololo[i],
              normhihihi,normlohihi,normhilohi,normlolohi,
              normhihilo,normlohilo,normhilolo,normlololo,
              &vhihihi[i],&vlohihi[i],&vhilohi[i],&vlolohi[i],
              &vhihilo[i],&vlohilo[i],&vhilolo[i],&vlololo[i]);
}
