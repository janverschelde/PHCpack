// The file random8_vectors.cpp defines the code for the functions
// specified in random8_vectors.h.

#include "random_numbers.h"
#include "random8_vectors.h"
#include "octo_double_functions.h"

void random_octo_double
 ( double *x_hihihi, double *x_lohihi, double *x_hilohi, double *x_lolohi,
   double *x_hihilo, double *x_lohilo, double *x_hilolo, double *x_lololo )
{
   const double eps = 2.220446049250313e-16; // 2^(-52)
   const double eps2 = eps*eps;
   const double eps3 = eps*eps2;
   const double eps4 = eps*eps3;
   const double eps5 = eps*eps4;
   const double eps6 = eps*eps5;
   const double eps7 = eps*eps6;
   const double r0 = random_double();
   const double r1 = random_double();
   const double r2 = random_double();
   const double r3 = random_double();
   const double r4 = random_double();
   const double r5 = random_double();
   const double r6 = random_double();
   const double r7 = random_double();

   *x_hihihi = r0;  *x_lohihi = 0.0; *x_hilohi = 0.0; *x_lolohi = 0.0;
   *x_hihilo = 0.0; *x_lohilo = 0.0; *x_hilolo = 0.0; *x_lololo = 0.0;

   odf_inc_d(x_hihihi,x_lohihi,x_hilohi,x_lolohi,
             x_hihilo,x_lohilo,x_hilolo,x_lololo,r1*eps);
   odf_inc_d(x_hihihi,x_lohihi,x_hilohi,x_lolohi,
             x_hihilo,x_lohilo,x_hilolo,x_lololo,r2*eps2);
   odf_inc_d(x_hihihi,x_lohihi,x_hilohi,x_lolohi,
             x_hihilo,x_lohilo,x_hilolo,x_lololo,r3*eps3);
   odf_inc_d(x_hihihi,x_lohihi,x_hilohi,x_lolohi,
             x_hihilo,x_lohilo,x_hilolo,x_lololo,r4*eps4);
   odf_inc_d(x_hihihi,x_lohihi,x_hilohi,x_lolohi,
             x_hihilo,x_lohilo,x_hilolo,x_lololo,r5*eps5);
   odf_inc_d(x_hihihi,x_lohihi,x_hilohi,x_lolohi,
             x_hihilo,x_lohilo,x_hilolo,x_lololo,r6*eps6);
   odf_inc_d(x_hihihi,x_lohihi,x_hilohi,x_lolohi,
             x_hihilo,x_lohilo,x_hilolo,x_lololo,r7*eps7);
}

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
   double r_hihihi,r_lohihi,r_hilohi,r_lolohi;
   double r_hihilo,r_lohilo,r_hilolo,r_lololo;

   for(int k=0; k<dim; k++)
   {
      random_octo_double(&r_hihihi,&r_lohihi,&r_hilohi,&r_lolohi,
                         &r_hihilo,&r_lohilo,&r_hilolo,&r_lololo);

      vhihihi_host[k] = r_hihihi; vhihihi_device[k] = r_hihihi;
      vlohihi_host[k] = r_lohihi; vlohihi_device[k] = r_lohihi;
      vhilohi_host[k] = r_hilohi; vhilohi_device[k] = r_hilohi;
      vlolohi_host[k] = r_lolohi; vlolohi_device[k] = r_lolohi;
      vhihilo_host[k] = r_hihilo; vhihilo_device[k] = r_hihilo;
      vlohilo_host[k] = r_lohilo; vlohilo_device[k] = r_lohilo;
      vhilolo_host[k] = r_hilolo; vhilolo_device[k] = r_hilolo;
      vlololo_host[k] = r_lololo; vlololo_device[k] = r_lololo;
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
   double *vrehihilo_device, double *vrelohilo_device,
   double *vrehilolo_device, double *vrelololo_device,
   double *vimhihihi_device, double *vimlohihi_device,
   double *vimhilohi_device, double *vimlolohi_device,
   double *vimhihilo_device, double *vimlohilo_device,
   double *vimhilolo_device, double *vimlololo_device )
{
   double rnd_hihihi,rnd_lohihi,rnd_hilohi,rnd_lolohi;
   double rnd_hihilo,rnd_lohilo,rnd_hilolo,rnd_lololo;
   double cosrnd_hihihi,cosrnd_lohihi,cosrnd_hilohi,cosrnd_lolohi;
   double cosrnd_hihilo,cosrnd_lohilo,cosrnd_hilolo,cosrnd_lololo;
   double sinrnd_hihihi,sinrnd_lohihi,sinrnd_hilohi,sinrnd_lolohi;
   double sinrnd_hihilo,sinrnd_lohilo,sinrnd_hilolo,sinrnd_lololo;

   for(int k=0; k<dim; k++)
   {
      random_octo_double(&rnd_hihihi,&rnd_lohihi,&rnd_hilohi,&rnd_lolohi,
                         &rnd_hihilo,&rnd_lohilo,&rnd_hilolo,&rnd_lololo);

      sinrnd_hihihi = rnd_hihihi; sinrnd_lohihi = rnd_lohihi;
      sinrnd_hilohi = rnd_hilohi; sinrnd_lolohi = rnd_lolohi;
      sinrnd_hihilo = rnd_hihilo; sinrnd_lohilo = rnd_lohilo;
      sinrnd_hilolo = rnd_hilolo; sinrnd_lololo = rnd_lololo;

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
