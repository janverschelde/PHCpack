// The file random4_vectors.cpp defines the code for the functions
// specified in random4_vectors.h.

#include "random_numbers.h"
#include "random4_vectors.h"
#include "quad_double_functions.h"

void random_double4_vectors
 ( int dim, double *vhihi_host, double *vlohi_host,
            double *vhilo_host, double *vlolo_host,
   double *vhihi_device, double *vlohi_device,
   double *vhilo_device, double *vlolo_device )
{
   double r;

   for(int k=0; k<dim; k++)
   {
      r = random_double();
      vhihi_host[k] = r;   vhihi_device[k] = r;
      vlohi_host[k] = 0.0; vlohi_device[k] = 0.0;
      vhilo_host[k] = 0.0; vhilo_device[k] = 0.0;
      vlolo_host[k] = 0.0; vlolo_device[k] = 0.0;
   }
}

void random_complex4_vectors
 ( int dim, double *vrehihi_host, double *vrelohi_host,
            double *vrehilo_host, double *vrelolo_host,
            double *vimhihi_host, double *vimlohi_host,
            double *vimhilo_host, double *vimlolo_host,
   double *vrehihi_device, double *vrelohi_device,
   double *vrehilo_device, double *vrelolo_device,
   double *vimhihi_device, double *vimlohi_device,
   double *vimhilo_device, double *vimlolo_device )
{
   double rnd;
   double cosrnd_hihi,cosrnd_lohi,cosrnd_hilo,cosrnd_lolo;
   double sinrnd_hihi,sinrnd_lohi,sinrnd_hilo,sinrnd_lolo;

   for(int k=0; k<dim; k++)
   {
      rnd = random_double(); // rnd is in [-1, +1]
      sinrnd_hihi = rnd; sinrnd_lohi = 0.0;
      sinrnd_hilo = 0.0; sinrnd_lolo = 0.0;

      double y_hihi,y_lohi,y_hilo,y_lolo;   // work around to compute cos

      qdf_sqr(sinrnd_hihi,sinrnd_lohi,sinrnd_hilo,sinrnd_lolo,
              &y_hihi,&y_lohi,&y_hilo,&y_lolo);
      qdf_minus(&y_hihi,&y_lohi,&y_hilo,&y_lolo);      // y = -sin^2
      qdf_inc_d(&y_hihi,&y_lohi,&y_hilo,&y_lolo,1.0);  // y = 1 - sin^2
      qdf_sqrt(y_hihi,y_lohi,y_hilo,y_lolo,
               &cosrnd_hihi,&cosrnd_lohi,&cosrnd_hilo,&cosrnd_lolo);
      // cos is sqrt(1-sin^2)

      vrehihi_host[k] = cosrnd_hihi; vrehihi_device[k] = cosrnd_hihi;
      vrelohi_host[k] = cosrnd_lohi; vrelohi_device[k] = cosrnd_lohi;
      vrehilo_host[k] = cosrnd_hilo; vrehilo_device[k] = cosrnd_hilo;
      vrelolo_host[k] = cosrnd_lolo; vrelolo_device[k] = cosrnd_lolo;

      vimhihi_host[k] = sinrnd_hihi; vimhihi_device[k] = sinrnd_hihi;
      vimlohi_host[k] = sinrnd_lohi; vimlohi_device[k] = sinrnd_lohi;
      vimhilo_host[k] = sinrnd_hilo; vimhilo_device[k] = sinrnd_hilo;
      vimlolo_host[k] = sinrnd_lolo; vimlolo_device[k] = sinrnd_lolo;
   }
}
