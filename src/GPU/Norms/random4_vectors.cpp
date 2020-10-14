// The file random4_vectors.cpp defines the code for the functions
// specified in random4_vectors.h.

#include "random_numbers.h"
#include "random4_vectors.h"
#include "quad_double_functions.h"

void random_quad_double
 ( double *x_hihi, double *x_lohi, double *x_hilo, double *x_lolo )
{
   const double eps = 2.220446049250313e-16; // 2^(-52)
   const double eps2 = eps*eps;
   const double eps3 = eps*eps2;
   const double r0 = random_double();
   const double r1 = random_double();
   const double r2 = random_double();
   const double r3 = random_double();

   *x_hihi = r0; *x_lohi = 0.0; *x_hilo = 0.0; *x_lolo = 0.0;

   qdf_inc_d(x_hihi,x_lohi,x_hilo,x_lolo,r1*eps);
   qdf_inc_d(x_hihi,x_lohi,x_hilo,x_lolo,r2*eps2);
   qdf_inc_d(x_hihi,x_lohi,x_hilo,x_lolo,r3*eps3);
}

void random_double4_vectors
 ( int dim, double *vhihi_host, double *vlohi_host,
            double *vhilo_host, double *vlolo_host,
   double *vhihi_device, double *vlohi_device,
   double *vhilo_device, double *vlolo_device )
{
   double r_hihi,r_lohi,r_hilo,r_lolo;

   for(int k=0; k<dim; k++)
   {
      random_quad_double(&r_hihi,&r_lohi,&r_hilo,&r_lolo);
      vhihi_host[k] = r_hihi; vhihi_device[k] = r_hihi;
      vlohi_host[k] = r_lohi; vlohi_device[k] = r_lohi;
      vhilo_host[k] = r_hilo; vhilo_device[k] = r_hilo;
      vlolo_host[k] = r_lolo; vlolo_device[k] = r_lolo;
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
   double rnd_hihi,rnd_lohi,rnd_hilo,rnd_lolo;
   double cosrnd_hihi,cosrnd_lohi,cosrnd_hilo,cosrnd_lolo;
   double sinrnd_hihi,sinrnd_lohi,sinrnd_hilo,sinrnd_lolo;

   for(int k=0; k<dim; k++)
   {
      random_quad_double(&rnd_hihi,&rnd_lohi,&rnd_hilo,&rnd_lolo);
      sinrnd_hihi = rnd_hihi; sinrnd_lohi = rnd_lohi;
      sinrnd_hilo = rnd_hilo; sinrnd_lolo = rnd_lolo;

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
