// The file random8_vectors.cpp defines the code for the functions
// specified in random8_vectors.h.

#include "random_numbers.h"
#include "random8_vectors.h"
#include "octo_double_functions.h"

void random_double8_vectors
 ( int dim, double *vhihihi_host, double *vlohihi_host,
            double *vhilohi_host, double *vlolohi_host,
            double *vhihilo_host, double *vlohilo_host,
            double *vhilolo_host, double *vlololo_host,
   double *vhihihi_device, double *vlohihi_device,
   double *vhilohi_device, double *vlolohi_device,
   double *vhihilo_device, double *vlohilo_device,
   double *vhilolo_device, double *vlololo_device )
{
   double r;

   for(int k=0; k<dim; k++)
   {
      r = random_double();
      vhihihi_host[k] = r;   vhihihi_device[k] = r;
      vlohihi_host[k] = 0.0; vlohihi_device[k] = 0.0;
      vhilohi_host[k] = 0.0; vhilohi_device[k] = 0.0;
      vlolohi_host[k] = 0.0; vlolohi_device[k] = 0.0;
      vhihilo_host[k] = 0.0; vhihilo_device[k] = 0.0;
      vlohilo_host[k] = 0.0; vlohilo_device[k] = 0.0;
      vhilolo_host[k] = 0.0; vhilolo_device[k] = 0.0;
      vlololo_host[k] = 0.0; vlololo_device[k] = 0.0;
   }
}

void random_complex8_vectors
 ( int dim, double *vrehihihi_host, double *vrelohihi_host,
            double *vrehilohi_host, double *vrelolohi_host,
            double *vrehihilo_host, double *vrelohilo_host,
            double *vrehilolo_host, double *vrelololo_host,
            double *vimhihihi_host, double *vimlohihi_host,
            double *vimhilohi_host, double *vimlolohi_host,
            double *vimhihilo_host, double *vimlohilo_host,
            double *vimhilolo_host, double *vimlololo_host,
   double *vrehihihi_device, double *vrelohihi_device,
   double *vrehilohi_device, double *vrelolohi_device,
   double *vimhihihi_device, double *vimlohihi_device,
   double *vimhilohi_device, double *vimlolohi_device,
   double *vrehihilo_device, double *vrelohilo_device,
   double *vrehilolo_device, double *vrelololo_device,
   double *vimhihilo_device, double *vimlohilo_device,
   double *vimhilolo_device, double *vimlololo_device )
{
   double rnd;
   double cosrnd_hihihi,cosrnd_lohihi,cosrnd_hilohi,cosrnd_lolohi;
   double cosrnd_hihilo,cosrnd_lohilo,cosrnd_hilolo,cosrnd_lololo;
   double sinrnd_hihihi,sinrnd_lohihi,sinrnd_hilohi,sinrnd_lolohi;
   double sinrnd_hihilo,sinrnd_lohilo,sinrnd_hilolo,sinrnd_lololo;

   for(int k=0; k<dim; k++)
   {
      rnd = random_double(); // rnd is in [-1, +1]
      sinrnd_hihihi = rnd; sinrnd_lohihi = 0.0;
      sinrnd_hilohi = 0.0; sinrnd_lolohi = 0.0;
      sinrnd_hihilo = 0.0; sinrnd_lohilo = 0.0;
      sinrnd_hilolo = 0.0; sinrnd_lololo = 0.0;

      double y_hihihi,y_lohihi,y_hilohi,y_lolohi;   // to compute cos
      double y_hihilo,y_lohilo,y_hilolo,y_lololo;

      odf_sqr(sinrnd_hihihi,sinrnd_lohihi,sinrnd_hilohi,sinrnd_lolohi,
              sinrnd_hihilo,sinrnd_lohilo,sinrnd_hilolo,sinrnd_lololo,
              &y_hihihi,&y_lohihi,&y_hilohi,&y_lolohi,
              &y_hihilo,&y_lohilo,&y_hilolo,&y_lololo);
      odf_minus(&y_hihihi,&y_lohihi,&y_hilohi,&y_lolohi,
                &y_hihilo,&y_lohilo,&y_hilolo,&y_lololo);     // y = -sin^2
      odf_inc_d(&y_hihihi,&y_lohihi,&y_hilohi,&y_lolohi,
                &y_hihilo,&y_lohilo,&y_hilolo,&y_lololo,1.0); // y = 1 - sin^2
      odf_sqrt(y_hihihi,y_lohihi,y_hilohi,y_lolohi,
               y_hihilo,y_lohilo,y_hilolo,y_lololo,
               &cosrnd_hihihi,&cosrnd_lohihi,&cosrnd_hilohi,&cosrnd_lolohi,
               &cosrnd_hihilo,&cosrnd_lohilo,&cosrnd_hilolo,&cosrnd_lololo);
      // cos is sqrt(1-sin^2)

      vrehihihi_host[k] = cosrnd_hihihi; vrehihihi_device[k] = cosrnd_hihihi;
      vrelohihi_host[k] = cosrnd_lohihi; vrelohihi_device[k] = cosrnd_lohihi;
      vrehilohi_host[k] = cosrnd_hilohi; vrehilohi_device[k] = cosrnd_hilohi;
      vrelolohi_host[k] = cosrnd_lolohi; vrelolohi_device[k] = cosrnd_lolohi;
      vrehihilo_host[k] = cosrnd_hihilo; vrehihilo_device[k] = cosrnd_hihilo;
      vrelohilo_host[k] = cosrnd_lohilo; vrelohilo_device[k] = cosrnd_lohilo;
      vrehilolo_host[k] = cosrnd_hilolo; vrehilolo_device[k] = cosrnd_hilolo;
      vrelololo_host[k] = cosrnd_lololo; vrelololo_device[k] = cosrnd_lololo;

      vimhihihi_host[k] = sinrnd_hihihi; vimhihihi_device[k] = sinrnd_hihihi;
      vimlohihi_host[k] = sinrnd_lohihi; vimlohihi_device[k] = sinrnd_lohihi;
      vimhilohi_host[k] = sinrnd_hilohi; vimhilohi_device[k] = sinrnd_hilohi;
      vimlolohi_host[k] = sinrnd_lolohi; vimlolohi_device[k] = sinrnd_lolohi;
      vimhihilo_host[k] = sinrnd_hihilo; vimhihilo_device[k] = sinrnd_hihilo;
      vimlohilo_host[k] = sinrnd_lohilo; vimlohilo_device[k] = sinrnd_lohilo;
      vimhilolo_host[k] = sinrnd_hilolo; vimhilolo_device[k] = sinrnd_hilolo;
      vimlololo_host[k] = sinrnd_lololo; vimlololo_device[k] = sinrnd_lololo;
   }
}
