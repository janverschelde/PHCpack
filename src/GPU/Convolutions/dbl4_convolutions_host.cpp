/* The file dbl4_convolutions_host.cpp defines the functions
 * specified in dbl4_convolutions_host.h. */

#include "dbl4_convolutions_host.h"
#include "quad_double_functions.h"

void CPU_dbl4_product
 ( int deg, double *xhihi, double *xlohi, double *xhilo, double *xlolo,
            double *yhihi, double *ylohi, double *yhilo, double *ylolo,
            double *zhihi, double *zlohi, double *zhilo, double *zlolo )
{
   int idx;
   double phihi,plohi,philo,plolo;

   qdf_mul(xhihi[0],xlohi[0],xhilo[0],xlolo[0],
           yhihi[0],ylohi[0],yhilo[0],ylolo[0],
           &zhihi[0],&zlohi[0],&zhilo[0],&zlolo[0]);     // z[0] = x[0]*y[0]

   for(int k=1; k<=deg; k++)
   {
      qdf_mul(xhihi[0],xlohi[0],xhilo[0],xlolo[0],
              yhihi[k],ylohi[k],yhilo[k],ylolo[k],
              &zhihi[k],&zlohi[k],&zhilo[k],&zlolo[k]);  // z[k] = x[0]*y[k]
      for(int i=1; i<=k; i++)
      {
         idx = k-i;
         qdf_mul(xhihi[i],xlohi[i],xhilo[i],xlolo[i],
                 yhihi[idx],ylohi[idx],yhilo[idx],ylolo[idx],
                 &phihi,&plohi,&philo,&plolo);                // x[i]*y[k-i]
         // z[k] = z[k] + x[i]*y[k-i]
         qdf_inc(&zhihi[k],&zlohi[k],&zhilo[k],&zlolo[k],
                 phihi,plohi,philo,plolo); 
      }
   }
}

void CPU_cmplx4_product
 ( int deg,
   double *xrehihi, double *xrelohi, double *xrehilo, double *xrelolo,
   double *ximhihi, double *ximlohi, double *ximhilo, double *ximlolo,
   double *yrehihi, double *yrelohi, double *yrehilo, double *yrelolo,
   double *yimhihi, double *yimlohi, double *yimhilo, double *yimlolo,
   double *zrehihi, double *zrelohi, double *zrehilo, double *zrelolo,
   double *zimhihi, double *zimlohi, double *zimhilo,  double *zimlolo )
{
   double rpahihi,rpalohi,rpahilo,rpalolo;    // accumulates real parts
   double ipahihi,ipalohi,ipahilo,ipalolo;    // accumulates imaginary parts
   double tmphihi,tmplohi,tmphilo,tmplolo;    // temporary quad double
   double xr0hihi,xr0lohi,xr0hilo,xr0lolo;    // real values of x
   double xi0hihi,xi0lohi,xi0hilo,xi0lolo;    // imaginary values of x
   double yr0hihi,yr0lohi,yr0hilo,yr0lolo;    // real values of y
   double yi0hihi,yi0lohi,yi0hilo,yi0lolo;    // imaginary values of y
   int idx;

   xr0hihi = xrehihi[0]; xr0lohi = xrelohi[0];
   xr0hilo = xrehilo[0]; xr0lolo = xrelolo[0];
   xi0hihi = ximhihi[0]; xi0lohi = ximlohi[0];
   xi0hilo = ximhilo[0]; xi0lolo = ximlolo[0];
   yr0hihi = yrehihi[0]; yr0lohi = yrelohi[0];
   yr0hilo = yrehilo[0]; yr0lolo = yrelolo[0];
   yi0hihi = yimhihi[0]; yi0lohi = yimlohi[0];
   yi0hilo = yimhilo[0]; yi0lolo = yimlolo[0];
   // zre[0] = xr0*yr0 - xi0*yi0;
   qdf_mul(xr0hihi,xr0lohi,xr0hilo,xr0lolo,
           yr0hihi,yr0lohi,yr0hilo,yr0lolo,
           &zrehihi[0],&zrelohi[0],&zrehilo[0],&zrelolo[0]);
   qdf_mul(xi0hihi,xi0lohi,xi0hilo,xi0lolo,
           yi0hihi,yi0lohi,yi0hilo,yi0lolo,
           &rpahihi,&rpalohi,&rpahilo,&rpalolo);
   qdf_minus(&rpahihi,&rpalohi,&rpahilo,&rpalolo);
   qdf_inc(&zrehihi[0],&zrelohi[0],&zrehilo[0],&zrelolo[0],
           rpahihi,rpalohi,rpahilo,rpalolo);
   // zim[0] = xi0*yr0 + xr0*yi0;
   qdf_mul(xi0hihi,xi0lohi,xi0hilo,xi0lolo,
           yr0hihi,yr0lohi,yr0hilo,yr0lolo,
           &zimhihi[0],&zimlohi[0],&zimhilo[0],&zimlolo[0]);
   qdf_mul(xr0hihi,xr0lohi,xr0hilo,xr0lolo,
           yi0hihi,yi0lohi,yi0hilo,yi0lolo,
           &ipahihi,&ipalohi,&ipahilo,&ipalolo);
   qdf_inc(&zimhihi[0],&zimlohi[0],&zimhilo[0],&zimlolo[0],
           ipahihi,ipalohi,ipahilo,ipalolo);

   for(int k=1; k<=deg; k++)
   {
      xr0hihi = xrehihi[0]; xr0lohi = xrelohi[0];
      xr0hilo = xrehilo[0]; xr0lolo = xrelolo[0];
      xi0hihi = ximhihi[0]; xi0lohi = ximlohi[0];
      xi0hilo = ximhilo[0]; xi0lolo = ximlolo[0];
      yr0hihi = yrehihi[k]; yr0lohi = yrelohi[k];
      yr0hilo = yrehilo[k]; yr0lolo = yrelolo[k];
      yi0hihi = yimhihi[k]; yi0lohi = yimlohi[k];
      yi0hilo = yimhilo[k]; yi0lolo = yimlolo[k];
      // rpa = xr0*yr0 - xi0*yi0;
      qdf_mul(xr0hihi,xr0lohi,xr0hilo,xr0lolo,
              yr0hihi,yr0lohi,yr0hilo,yr0lolo,
              &rpahihi,&rpalohi,&rpahilo,&rpalolo);
      qdf_mul(xi0hihi,xi0lohi,xi0hilo,xi0lolo,
              yi0hihi,yi0lohi,yi0hilo,yi0lolo,
              &tmphihi,&tmplohi,&tmphilo,&tmplolo);
      qdf_minus(&tmphihi,&tmplohi,&tmphilo,&tmplolo);
      qdf_inc(&rpahihi,&rpalohi,&rpahilo,&rpalolo,
              tmphihi,tmplohi,tmphilo,tmplolo);
      // ipa = xi0*yr0 + xr0*yi0;
      qdf_mul(xi0hihi,xi0lohi,xi0hilo,xi0lolo,
              yr0hihi,yr0lohi,yr0hilo,yr0lolo,
              &ipahihi,&ipalohi,&ipahilo,&ipalolo);
      qdf_mul(xr0hihi,xr0lohi,xr0hilo,xr0lolo,
              yi0hihi,yi0lohi,yi0hilo,yi0lolo,
              &tmphihi,&tmplohi,&tmphilo,&tmplolo);
      qdf_inc(&ipahihi,&ipalohi,&ipahilo,&ipalolo,
              tmphihi,tmplohi,tmphilo,tmplolo);

      for(int i=1; i<=k; i++)
      {
         idx = k-i;
         xr0hihi = xrehihi[i]; xr0lohi = xrelohi[i];
         xr0hilo = xrehilo[i]; xr0lolo = xrelolo[i];
         xi0hihi = ximhihi[i]; xi0lohi = ximlohi[i];
         xi0hilo = ximhilo[i]; xi0lolo = ximlolo[i];
         yr0hihi = yrehihi[idx]; yr0lohi = yrelohi[idx];
         yr0hilo = yrehilo[idx]; yr0lolo = yrelolo[idx];
         yi0hihi = yimhihi[idx]; yi0lohi = yimlohi[idx];
         yi0hilo = yimhilo[idx]; yi0lolo = yimlolo[idx];
         // rpa = rpa + xr0*yr0 - xi0*yi0;
         qdf_mul(xr0hihi,xr0lohi,xr0hilo,xr0lolo,
                 yr0hihi,yr0lohi,yr0hilo,yr0lolo,
                 &tmphihi,&tmplohi,&tmphilo,&tmplolo);
         qdf_inc(&rpahihi,&rpalohi,&rpahilo,&rpalolo,
                 tmphihi,tmplohi,tmphilo,tmplolo);
         qdf_mul(xi0hihi,xi0lohi,xi0hilo,xi0lolo,
                 yi0hihi,yi0lohi,yi0hilo,yi0lolo,
                 &tmphihi,&tmplohi,&tmphilo,&tmplolo);
         qdf_minus(&tmphihi,&tmplohi,&tmphilo,&tmplolo);
         qdf_inc(&rpahihi,&rpalohi,&rpahilo,&rpalolo,
                 tmphihi,tmplohi,tmphilo,tmplolo);
         // ipa = ipa + xi0*yr0 + xr0*yi0;
         qdf_mul(xi0hihi,xi0lohi,xi0hilo,xi0lolo,
                 yr0hihi,yr0lohi,yr0hilo,yr0lolo,
                 &tmphihi,&tmplohi,&tmphilo,&tmplolo);
         qdf_inc(&ipahihi,&ipalohi,&ipahilo,&ipalolo,
                 tmphihi,tmplohi,tmphilo,tmplolo);
         qdf_mul(xr0hihi,xr0lohi,xr0hilo,xr0lolo,
                 yi0hihi,yi0lohi,yi0hilo,yi0lolo,
                 &tmphihi,&tmplohi,&tmphilo,&tmplolo);
         qdf_inc(&ipahihi,&ipalohi,&ipahilo,&ipalolo,
                 tmphihi,tmplohi,tmphilo,tmplolo);
      }
      zrehihi[k] = rpahihi; zrelohi[k] = rpalohi;
      zrehilo[k] = rpahilo; zrelolo[k] = rpalolo;
      zimhihi[k] = ipahihi; zimlohi[k] = ipalohi;
      zimhilo[k] = ipahilo; zimlolo[k] = ipalolo;
   }
}
