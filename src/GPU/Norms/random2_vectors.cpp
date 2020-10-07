// The file random_vectors.cpp defines the code for the functions
// specified in random_vectors.h.

#include <cstdlib>
#include <cmath>
#include "random_numbers.h"
#include "random2_vectors.h"
#include "double_double.h"

void random_double2_vectors
 ( int dim, double* vhi_host, double* vlo_host,
   double* vhi_device, double* vlo_device )
{
   double r;

   for(int k=0; k<dim; k++)
   {
      r = random_double();
      vhi_host[k] = r;   vhi_device[k] = r;
      vlo_host[k] = 0.0; vlo_device[k] = 0.0;
   }
}

void random_complex2_vectors
 ( int dim, double* vrehi_host, double* vrelo_host,
            double* vimhi_host, double* vimlo_host,
   double* vrehi_device, double* vrelo_device,
   double* vimhi_device, double* vimlo_device )
{
   double rnd;
   double ddrnd[2],cosrnd[2],sinrnd[2];

   for(int k=0; k<dim; k++)
   {
      rnd = random_angle();

      ddrnd[0] = rnd; ddrnd[1] = 0.0;

      dd_sin(ddrnd,sinrnd); dd_cos(ddrnd,cosrnd);

      // sinrnd[0] = sin(rnd); sinrnd[1] = 0.0;
      // cosrnd[0] = cos(rnd); cosrnd[1] = 0.0;

      vrehi_host[k] = cosrnd[0]; vrehi_device[k] = cosrnd[0];
      vrelo_host[k] = cosrnd[1]; vrelo_device[k] = cosrnd[1];
      vimhi_host[k] = sinrnd[0]; vimhi_device[k] = sinrnd[0];
      vimlo_host[k] = sinrnd[1]; vimlo_device[k] = cosrnd[1];
   }
}
