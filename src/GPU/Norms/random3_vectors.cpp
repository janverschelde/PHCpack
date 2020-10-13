// The file random3_vectors.cpp defines the code for the functions
// specified in random3_vectors.h.

#include "random_numbers.h"
#include "random3_vectors.h"
#include "triple_double_functions.h"

void random_double3_vectors
 ( int dim, double *vhi_host, double *vmi_host, double *vlo_host,
   double *vhi_device, double *vmi_device, double *vlo_device )
{
   double r;

   for(int k=0; k<dim; k++)
   {
      r = random_double();
      vhi_host[k] = r;   vhi_device[k] = r;
      vmi_host[k] = 0.0; vmi_device[k] = 0.0;
      vlo_host[k] = 0.0; vlo_device[k] = 0.0;
   }
}

void random_complex3_vectors
 ( int dim, double *vrehi_host, double *vremi_host, double *vrelo_host,
            double *vimhi_host, double *vimmi_host, double *vimlo_host,
   double *vrehi_device, double *vremi_device, double *vrelo_device,
   double *vimhi_device, double *vimmi_device, double *vimlo_device )
{
   double rnd;
   double cosrnd_hi,cosrnd_mi,cosrnd_lo,sinrnd_hi,sinrnd_mi,sinrnd_lo;

   for(int k=0; k<dim; k++)
   {
      rnd = random_double(); // rnd is in [-1, +1]
      sinrnd_hi = rnd; sinrnd_mi = 0.0; sinrnd_lo = 0.0;

      double y_hi,y_mi,y_lo;          // work around to compute cos

      tdf_sqr(sinrnd_hi,sinrnd_mi,sinrnd_lo,&y_hi,&y_mi,&y_lo);
      tdf_minus(&y_hi,&y_mi,&y_lo);              // y = -sin^2
      tdf_inc_d(&y_hi,&y_mi,&y_lo,1.0);          // y = 1 - sin^2
      tdf_sqrt(y_hi,y_mi,y_lo,&cosrnd_hi,&cosrnd_mi,&cosrnd_lo);
      // cos is sqrt(1-sin^2)

      vrehi_host[k] = cosrnd_hi; vrehi_device[k] = cosrnd_hi;
      vremi_host[k] = cosrnd_mi; vrelo_device[k] = cosrnd_mi;
      vrelo_host[k] = cosrnd_lo; vrelo_device[k] = cosrnd_lo;
      vimhi_host[k] = sinrnd_hi; vimhi_device[k] = sinrnd_hi;
      vimlo_host[k] = sinrnd_mi; vimlo_device[k] = sinrnd_mi;
      vimlo_host[k] = sinrnd_lo; vimlo_device[k] = sinrnd_lo;
   }
}
