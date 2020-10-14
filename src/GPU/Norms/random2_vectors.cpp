// The file random2_vectors.cpp defines the code for the functions
// specified in random2_vectors.h.

#include "random_numbers.h"
#include "random2_vectors.h"
#include "double_double_functions.h"

void random_double_double ( double *x_hi, double *x_lo )
{
   const double eps = 2.220446049250313e-16; // 2^(-52)
   const double r0 = random_double();
   const double r1 = random_double();

   *x_hi = r0; *x_lo = 0.0;

   ddf_inc_d(x_hi,x_lo,r1*eps);
}

void random_double2_vectors
 ( int dim, double *vhi_host, double *vlo_host,
            double *vhi_device, double *vlo_device )
{
   double r_hi,r_lo;

   for(int k=0; k<dim; k++)
   {
      random_double_double(&r_hi,&r_lo);
      vhi_host[k] = r_hi;   vhi_device[k] = r_hi;
      vlo_host[k] = r_lo; vlo_device[k] = r_lo;
   }
}

void random_complex2_vectors
 ( int dim, double *vrehi_host, double *vrelo_host,
            double *vimhi_host, double *vimlo_host,
   double *vrehi_device, double *vrelo_device,
   double *vimhi_device, double *vimlo_device )
{
   double rnd_hi,rnd_lo;
   double cosrnd_hi,cosrnd_lo,sinrnd_hi,sinrnd_lo;

   for(int k=0; k<dim; k++)
   {
      random_double_double(&rnd_hi,&rnd_lo);
      sinrnd_hi = rnd_hi; sinrnd_lo = rnd_lo;

      double y_hi,y_lo;          // work around to compute cos

      ddf_sqr(sinrnd_hi,sinrnd_lo,&y_hi,&y_lo);
      ddf_minus(&y_hi,&y_lo);                    // y = -sin^2
      ddf_inc_d(&y_hi,&y_lo,1.0);                // y = 1 - sin^2
      ddf_sqrt(y_hi,y_lo,&cosrnd_hi,&cosrnd_lo); // cos is sqrt(1-sin^2)

      vrehi_host[k] = cosrnd_hi; vrehi_device[k] = cosrnd_hi;
      vrelo_host[k] = cosrnd_lo; vrelo_device[k] = cosrnd_lo;
      vimhi_host[k] = sinrnd_hi; vimhi_device[k] = sinrnd_hi;
      vimlo_host[k] = sinrnd_lo; vimlo_device[k] = sinrnd_lo;
   }
}
