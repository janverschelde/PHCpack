// The file random3_vectors.cpp defines the code for the functions
// specified in random3_vectors.h.

#include "random_numbers.h"
#include "random3_vectors.h"
#include "triple_double_functions.h"

void random_triple_double
 ( double *x_hi, double *x_mi, double *x_lo )
{
   const double eps = 2.220446049250313e-16; // 2^(-52)
   const double eps2 = eps*eps;
   const double r0 = random_double();
   const double r1 = random_double();
   const double r2 = random_double();

   *x_hi = r0; *x_mi = 0.0; *x_lo = 0.0;

   tdf_inc_d(x_hi,x_mi,x_lo,r1*eps);
   tdf_inc_d(x_hi,x_mi,x_lo,r2*eps2);
}

void random_double3_vectors
 ( int dim, double *vhi_host, double *vmi_host, double *vlo_host,
   double *vhi_device, double *vmi_device, double *vlo_device )
{
   double r_hi,r_mi,r_lo;

   for(int k=0; k<dim; k++)
   {
      random_triple_double(&r_hi,&r_mi,&r_lo);
      vhi_host[k] = r_hi; vhi_device[k] = r_hi;
      vmi_host[k] = r_mi; vmi_device[k] = r_mi;
      vlo_host[k] = r_lo; vlo_device[k] = r_lo;
   }
}

void random_complex3_vectors
 ( int dim, double *vrehi_host, double *vremi_host, double *vrelo_host,
            double *vimhi_host, double *vimmi_host, double *vimlo_host,
   double *vrehi_device, double *vremi_device, double *vrelo_device,
   double *vimhi_device, double *vimmi_device, double *vimlo_device )
{
   double rnd_hi,rnd_mi,rnd_lo;
   double cosrnd_hi,cosrnd_mi,cosrnd_lo,sinrnd_hi,sinrnd_mi,sinrnd_lo;

   for(int k=0; k<dim; k++)
   {
      random_triple_double(&rnd_hi,&rnd_mi,&rnd_lo);
      sinrnd_hi = rnd_hi; sinrnd_mi = rnd_mi; sinrnd_lo = rnd_lo;

      double y_hi,y_mi,y_lo;          // work around to compute cos

      tdf_sqr(sinrnd_hi,sinrnd_mi,sinrnd_lo,&y_hi,&y_mi,&y_lo);
      tdf_minus(&y_hi,&y_mi,&y_lo);              // y = -sin^2
      tdf_inc_d(&y_hi,&y_mi,&y_lo,1.0);          // y = 1 - sin^2
      tdf_sqrt(y_hi,y_mi,y_lo,&cosrnd_hi,&cosrnd_mi,&cosrnd_lo);
      // cos is sqrt(1-sin^2)

      vrehi_host[k] = cosrnd_hi; vrehi_device[k] = cosrnd_hi;
      vremi_host[k] = cosrnd_mi; vremi_device[k] = cosrnd_mi;
      vrelo_host[k] = cosrnd_lo; vrelo_device[k] = cosrnd_lo;
      vimhi_host[k] = sinrnd_hi; vimhi_device[k] = sinrnd_hi;
      vimmi_host[k] = sinrnd_mi; vimmi_device[k] = sinrnd_mi;
      vimlo_host[k] = sinrnd_lo; vimlo_device[k] = sinrnd_lo;
   }
}
