/* The file dbl3_convolutions_host.cpp defines functions
 * specified in dbl3_convolutions_host.h. */

#include "dbl3_convolutions_host.h"
#include "triple_double_functions.h"

void CPU_dbl3_product
 ( int deg, double *xhi, double *xmi, double *xlo,
            double *yhi, double *ymi, double *ylo,
            double *zhi, double *zmi, double *zlo )
{
   int idx;
   double phi,pmi,plo;

   tdf_mul(xhi[0],xmi[0],xlo[0],yhi[0],ymi[0],ylo[0],
           &zhi[0],&zmi[0],&zlo[0]);                    // z[0] = x[0]*y[0]

   for(int k=1; k<=deg; k++)
   {
      tdf_mul(xhi[0],xmi[0],xlo[0],yhi[k],ymi[k],ylo[k],
              &zhi[k],&zmi[k],&zlo[k]);                 // z[k] = x[0]*y[k]
      for(int i=1; i<=k; i++)
      {
         idx = k-i;
         tdf_mul(xhi[i],xmi[i],xlo[i],yhi[idx],ymi[idx],ylo[idx],
                 &phi,&pmi,&plo);                       // x[i]*y[k-i]
         // z[k] = z[k] + x[i]*y[k-i]
         tdf_inc(&zhi[k],&zmi[k],&zlo[k],phi,pmi,plo); 
      }
   }
}

void CPU_cmplx3_product
 ( int deg, double *xrehi, double *xremi, double *xrelo,
            double *ximhi, double *ximmi, double *ximlo,
            double *yrehi, double *yremi, double *yrelo,
            double *yimhi, double *yimmi, double *yimlo,
            double *zrehi, double *zremi, double *zrelo,
            double *zimhi, double *zimmi, double *zimlo )
{
   double rpahi,rpami,rpalo;    // accumulates real parts
   double ipahi,ipami,ipalo;    // accumulates imaginary parts
   double tmphi,tmpmi,tmplo;    // temporary high, middle, and low doubles 
   double xr0hi,xr0mi,xr0lo,xi0hi,xi0mi,xi0lo; // values in xr and xi
   double yr0hi,yr0mi,yr0lo,yi0hi,yi0mi,yi0lo; // values in yr and yi
   int idx;

   xr0hi = xrehi[0]; xr0mi = xremi[0]; xr0lo = xrelo[0];
   xi0hi = ximhi[0]; xi0mi = ximmi[0]; xi0lo = ximlo[0];
   yr0hi = yrehi[0]; yr0mi = yremi[0]; yr0lo = yrelo[0];
   yi0hi = yimhi[0]; yi0mi = yimmi[0]; yi0lo = yimlo[0];
   // zre[0] = xr0*yr0 - xi0*yi0;
   tdf_mul(xr0hi,xr0mi,xr0lo,yr0hi,yr0mi,yr0lo,
           &zrehi[0],&zremi[0],&zrelo[0]);
   tdf_mul(xi0hi,xi0mi,xi0lo,yi0hi,yi0mi,yi0lo,
           &rpahi,&rpami,&rpalo);
   tdf_minus(&rpahi,&rpami,&rpalo);
   tdf_inc(&zrehi[0],&zremi[0],&zrelo[0],rpahi,rpami,rpalo);
   // zim[0] = xi0*yr0 + xr0*yi0;
   tdf_mul(xi0hi,xi0mi,xi0lo,yr0hi,yr0mi,yr0lo,
           &zimhi[0],&zimmi[0],&zimlo[0]);
   tdf_mul(xr0hi,xr0mi,xr0lo,yi0hi,yi0mi,yi0lo,
           &ipahi,&ipami,&ipalo);
   tdf_inc(&zimhi[0],&zimmi[0],&zimlo[0],ipahi,ipami,ipalo);

   for(int k=1; k<=deg; k++)
   {
      xr0hi = xrehi[0]; xr0mi = xremi[0]; xr0lo = xrelo[0];
      xi0hi = ximhi[0]; xi0mi = ximmi[0]; xi0lo = ximlo[0];
      yr0hi = yrehi[k]; yr0mi = yremi[k]; yr0lo = yrelo[k];
      yi0hi = yimhi[k]; yi0mi = yimmi[k]; yi0lo = yimlo[k];
      // rpa = xr0*yr0 - xi0*yi0;
      tdf_mul(xr0hi,xr0mi,xr0lo,yr0hi,yr0mi,yr0lo,
              &rpahi,&rpami,&rpalo);
      tdf_mul(xi0hi,xi0mi,xi0lo,yi0hi,yi0mi,yi0lo,
              &tmphi,&tmpmi,&tmplo);
      tdf_minus(&tmphi,&tmpmi,&tmplo);
      tdf_inc(&rpahi,&rpami,&rpalo,tmphi,tmpmi,tmplo);
      // ipa = xi0*yr0 + xr0*yi0;
      tdf_mul(xi0hi,xi0mi,xi0lo,yr0hi,yr0mi,yr0lo,
              &ipahi,&ipami,&ipalo);
      tdf_mul(xr0hi,xr0mi,xr0lo,yi0hi,yi0mi,yi0lo,
              &tmphi,&tmpmi,&tmplo);
      tdf_inc(&ipahi,&ipami,&ipalo,tmphi,tmpmi,tmplo);

      for(int i=1; i<=k; i++)
      {
         idx = k-i;
         xr0hi = xrehi[i]; xr0mi = xremi[i]; xr0lo = xrelo[i];
         xi0hi = ximhi[i]; xi0mi = ximmi[i]; xi0lo = ximlo[i];
         yr0hi = yrehi[idx]; yr0mi = yremi[idx]; yr0lo = yrelo[idx];
         yi0hi = yimhi[idx]; yi0mi = yimmi[idx]; yi0lo = yimlo[idx];
         // rpa = rpa + xr0*yr0 - xi0*yi0;
         tdf_mul(xr0hi,xr0mi,xr0lo,yr0hi,yr0mi,yr0lo,
                 &tmphi,&tmpmi,&tmplo);
         tdf_inc(&rpahi,&rpami,&rpalo,tmphi,tmpmi,tmplo);
         tdf_mul(xi0hi,xi0mi,xi0lo,yi0hi,yi0mi,yi0lo,
                 &tmphi,&tmpmi,&tmplo);
         tdf_minus(&tmphi,&tmpmi,&tmplo);
         tdf_inc(&rpahi,&rpami,&rpalo,tmphi,tmpmi,tmplo);
         // ipa = ipa + xi0*yr0 + xr0*yi0;
         tdf_mul(xi0hi,xi0mi,xi0lo,yr0hi,yr0mi,yr0lo,
                 &tmphi,&tmpmi,&tmplo);
         tdf_inc(&ipahi,&ipami,&ipalo,tmphi,tmpmi,tmplo);
         tdf_mul(xr0hi,xr0mi,xr0lo,yi0hi,yi0mi,yi0lo,
                 &tmphi,&tmpmi,&tmplo);
         tdf_inc(&ipahi,&ipami,&ipalo,tmphi,tmpmi,tmplo);
      }
      zrehi[k] = rpahi; zremi[k] = rpami; zrelo[k] = rpalo;
      zimhi[k] = ipahi; zimmi[k] = ipami; zimlo[k] = ipalo;
   }
}
