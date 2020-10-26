/* The file dbl_convolutions_host.cpp defines the functions 
 * specified in dbl_convolutions_host.h. */

#include "dbl_convolutions_host.h"

void CPU_dbl_product ( int deg, double *x, double *y, double *z )
{
   z[0] = x[0]*y[0];
   for(int k=0; k<=deg; k++)
   {
      z[k] = x[0]*y[k];
      for(int i=1; i<=k; i++)
         z[k] = z[k] + x[i]*y[k-i];
   }
}

void CPU_cmplx_product
 ( int deg, double *xre, double *xim, double *yre, double *yim,
            double *zre, double *zim )
{
   double rpa,ipa; // accumulates real and imaginary parts
   double xr0,xi0; // to hold values in xr and xi
   double yr0,yi0; // to hold values in yr and yi
   int idx;

   xr0 = xre[0]; xi0 = xim[0];
   yr0 = yre[0]; yi0 = yim[0];
   zre[0] = xr0*yr0 - xi0*yi0;
   zim[0] = xi0*yr0 + xr0*yi0;

   for(int k=0; k<=deg; k++)
   {

      xr0 = xre[0]; xi0 = xim[0];
      yr0 = yre[k]; yi0 = yim[k];
      rpa = xr0*yr0 - xi0*yi0;
      ipa = xi0*yr0 + xr0*yi0;

      for(int i=1; i<=k; i++)
      {
         idx = k-i;
         xr0 = xre[i];   xi0 = xim[i];
         yr0 = yre[idx]; yi0 = yim[idx];
         rpa = rpa + xr0*yr0 - xi0*yi0;
         ipa = ipa + xi0*yr0 + xr0*yi0;
      }
      zre[k] = rpa; zim[k] = ipa;
   }
}
