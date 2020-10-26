/* The file dbl2_convolutions_host.cpp defines the functions
 * specified in dbl2_convolutions_host.h */

#include "dbl2_convolutions_host.h"
#include "double_double_functions.h"

void CPU_dbl2_product
 ( int deg, double *xhi, double *xlo, double *yhi, double *ylo,
            double *zhi, double *zlo )
{
   int idx;
   double phi,plo;

   ddf_mul(xhi[0],xlo[0],yhi[0],ylo[0],&zhi[0],&zlo[0]);    // z[0] = x[0]*y[0]

   for(int k=0; k<=deg; k++)
   {
      ddf_mul(xhi[0],xlo[0],yhi[k],ylo[k],&zhi[k],&zlo[k]); // z[k] = x[0]*y[k]
      for(int i=1; i<=k; i++)
      {
         idx = k-i;
         ddf_mul(xhi[i],xlo[i],yhi[idx],ylo[idx],&phi,&plo);     // x[i]*y[k-i]
         ddf_inc(&zhi[k],&zlo[k],phi,plo);         // z[k] = z[k] + x[i]*y[k-i]
      }
   }
}

void CPU_cmplx2_product
 ( int deg, double *xrehi, double *xrelo, double *ximhi, double *ximlo,
            double *yrehi, double *yrelo, double *yimhi, double *yimlo,
            double *zrehi, double *zrelo, double *zimhi, double *zimlo )
{
   double rpahi,rpalo,ipahi,ipalo; // accumulates real and imaginary parts
   double tmphi,tmplo;             // temporary high and low doubles 
   double xr0hi,xr0lo,xi0hi,xi0lo; // to hold values in xr and xi
   double yr0hi,yr0lo,yi0hi,yi0lo; // to hold values in yr and yi
   int idx;

   xr0hi = xrehi[0]; xr0lo = xrelo[0];
   xi0hi = ximhi[0]; xi0lo = ximlo[0];
   yr0hi = yrehi[0]; yr0lo = yrelo[0];
   yi0hi = yimhi[0]; yi0lo = yimlo[0];
   // zre[0] = xr0*yr0 - xi0*yi0;
   ddf_mul(xr0hi,xr0lo,yr0hi,yr0lo,&zrehi[0],&zrelo[0]);
   ddf_mul(xi0hi,xi0lo,yi0hi,yi0lo,&rpahi,&rpalo);
   ddf_dec(&zrehi[0],&zrelo[0],rpahi,rpalo);
   // zim[0] = xi0*yr0 + xr0*yi0;
   ddf_mul(xi0hi,xi0lo,yr0hi,yr0lo,&zimhi[0],&zimlo[0]);
   ddf_mul(xr0hi,xr0lo,yi0hi,yi0lo,&ipahi,&ipalo);
   ddf_inc(&zimhi[0],&zimlo[0],ipahi,ipalo);

   for(int k=0; k<=deg; k++)
   {
      xr0hi = xrehi[0]; xr0lo = xrelo[0];
      xi0hi = ximhi[0]; xi0lo = ximlo[0];
      yr0hi = yrehi[k]; yr0lo = yrelo[k];
      yi0hi = yimhi[k]; yi0lo = yimlo[k];
      // rpa = xr0*yr0 - xi0*yi0;
      ddf_mul(xr0hi,xr0lo,yr0hi,yr0lo,&rpahi,&rpalo);
      ddf_mul(xi0hi,xi0lo,yi0hi,yi0lo,&tmphi,&tmplo);
      ddf_dec(&rpahi,&rpalo,tmphi,tmplo);
      // ipa = xi0*yr0 + xr0*yi0;
      ddf_mul(xi0hi,xi0lo,yr0hi,yr0lo,&ipahi,&ipalo);
      ddf_mul(xr0hi,xr0lo,yi0hi,yi0lo,&tmphi,&tmplo);
      ddf_inc(&ipahi,&ipalo,tmphi,tmplo);

      for(int i=1; i<=k; i++)
      {
         idx = k-i;
         xr0hi = xrehi[i]; xr0lo = xrelo[i];
         xi0hi = ximhi[i]; xi0lo = ximlo[i];
         yr0hi = yrehi[idx]; yr0lo = yrelo[idx];
         yi0hi = yimhi[idx]; yi0lo = yimlo[idx];
         // rpa = rpa + xr0*yr0 - xi0*yi0;
         ddf_mul(xr0hi,xr0lo,yr0hi,yr0lo,&tmphi,&tmplo);
         ddf_inc(&rpahi,&rpalo,tmphi,tmplo);
         ddf_mul(xi0hi,xi0lo,yi0hi,yi0lo,&tmphi,&tmplo);
         ddf_dec(&rpahi,&rpalo,tmphi,tmplo);
         // ipa = ipa + xi0*yr0 + xr0*yi0;
         ddf_mul(xi0hi,xi0lo,yr0hi,yr0lo,&tmphi,&tmplo);
         ddf_inc(&ipahi,&ipalo,tmphi,tmplo);
         ddf_mul(xr0hi,xr0lo,yi0hi,yi0lo,&tmphi,&tmplo);
         ddf_inc(&ipahi,&ipalo,tmphi,tmplo);
      }
      zrehi[k] = rpahi; zrelo[k] = rpalo;
      zimhi[k] = ipahi; zimlo[k] = ipalo;
   }
}
