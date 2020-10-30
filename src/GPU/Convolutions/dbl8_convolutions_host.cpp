/* The file dbl8_convolutions_host.cpp defines the functions
 * specified in dbl8_convolutions_host.h. */

#include "dbl8_convolutions_host.h"
#include "octo_double_functions.h"

void CPU_dbl8_product
 ( int deg,
   double *xhihihi, double *xlohihi, double *xhilohi, double *xlolohi,
   double *xhihilo, double *xlohilo, double *xhilolo, double *xlololo,
   double *yhihihi, double *ylohihi, double *yhilohi, double *ylolohi,
   double *yhihilo, double *ylohilo, double *yhilolo, double *ylololo,
   double *zhihihi, double *zlohihi, double *zhilohi, double *zlolohi,
   double *zhihilo, double *zlohilo, double *zhilolo, double *zlololo )
{
   int idx;
   double phihihi,plohihi,philohi,plolohi;
   double phihilo,plohilo,philolo,plololo;

   // z[0] = x[0]*y[0]
   odf_mul(xhihihi[0],xlohihi[0],xhilohi[0],xlolohi[0],
           xhihilo[0],xlohilo[0],xhilolo[0],xlololo[0],
           yhihihi[0],ylohihi[0],yhilohi[0],ylolohi[0],
           yhihilo[0],ylohilo[0],yhilolo[0],ylololo[0],
           &zhihihi[0],&zlohihi[0],&zhilohi[0],&zlolohi[0],
           &zhihilo[0],&zlohilo[0],&zhilolo[0],&zlololo[0]);

   for(int k=1; k<=deg; k++)
   {
      // z[k] = x[0]*y[k]
      odf_mul(xhihihi[0],xlohihi[0],xhilohi[0],xlolohi[0],
              xhihilo[0],xlohilo[0],xhilolo[0],xlololo[0],
              yhihihi[k],ylohihi[k],yhilohi[k],ylolohi[k],
              yhihilo[k],ylohilo[k],yhilolo[k],ylololo[k],
              &zhihihi[k],&zlohihi[k],&zhilohi[k],&zlolohi[k],
              &zhihilo[k],&zlohilo[k],&zhilolo[k],&zlololo[k]);
      for(int i=1; i<=k; i++)
      {
         idx = k-i;
         // x[i]*y[k-i]
         odf_mul(xhihihi[i],xlohihi[i],xhilohi[i],xlolohi[i],
                 xhihilo[i],xlohilo[i],xhilolo[i],xlololo[i],
                 yhihihi[idx],ylohihi[idx],yhilohi[idx],ylolohi[idx],
                 yhihilo[idx],ylohilo[idx],yhilolo[idx],ylololo[idx],
                 &phihihi,&plohihi,&philohi,&plolohi,
                 &phihilo,&plohilo,&philolo,&plololo);
         // z[k] = z[k] + x[i]*y[k-i]
         odf_inc(&zhihihi[k],&zlohihi[k],&zhilohi[k],&zlolohi[k],
                 &zhihilo[k],&zlohilo[k],&zhilolo[k],&zlololo[k],
                 phihihi,plohihi,philohi,plolohi,
                 phihilo,plohilo,philolo,plololo); 
      }
   }
}

void CPU_cmplx8_product
 ( int deg,
   double *xrehihihi, double *xrelohihi, double *xrehilohi, double *xrelolohi,
   double *xrehihilo, double *xrelohilo, double *xrehilolo, double *xrelololo,
   double *ximhihihi, double *ximlohihi, double *ximhilohi, double *ximlolohi,
   double *ximhihilo, double *ximlohilo, double *ximhilolo, double *ximlololo,
   double *yrehihihi, double *yrelohihi, double *yrehilohi, double *yrelolohi,
   double *yrehihilo, double *yrelohilo, double *yrehilolo, double *yrelololo,
   double *yimhihihi, double *yimlohihi, double *yimhilohi, double *yimlolohi,
   double *yimhihilo, double *yimlohilo, double *yimhilolo, double *yimlololo,
   double *zrehihihi, double *zrelohihi, double *zrehilohi, double *zrelolohi,
   double *zrehihilo, double *zrelohilo, double *zrehilolo, double *zrelololo,
   double *zimhihihi, double *zimlohihi, double *zimhilohi, double *zimlolohi,
   double *zimhihilo, double *zimlohilo, double *zimhilolo, double *zimlololo )
{
   double rpahihihi,rpalohihi,rpahilohi,rpalolohi;  // high real parts
   double rpahihilo,rpalohilo,rpahilolo,rpalololo;  // low real parts
   double ipahihihi,ipalohihi,ipahilohi,ipalolohi;  // high imaginary parts
   double ipahihilo,ipalohilo,ipahilolo,ipalololo;  // low imaginary parts
   double tmphihihi,tmplohihi,tmphilohi,tmplolohi;  // temporary high parts
   double tmphihilo,tmplohilo,tmphilolo,tmplololo;  // temporary low parts
   double xr0hihihi,xr0lohihi,xr0hilohi,xr0lolohi;  // high real values of x
   double xr0hihilo,xr0lohilo,xr0hilolo,xr0lololo;  // low real values of x
   double xi0hihihi,xi0lohihi,xi0hilohi,xi0lolohi;  // high imag values of x
   double xi0hihilo,xi0lohilo,xi0hilolo,xi0lololo;  // low imag values of x
   double yr0hihihi,yr0lohihi,yr0hilohi,yr0lolohi;  // high real values of y
   double yr0hihilo,yr0lohilo,yr0hilolo,yr0lololo;  // low real values of y
   double yi0hihihi,yi0lohihi,yi0hilohi,yi0lolohi;  // high imag values of y
   double yi0hihilo,yi0lohilo,yi0hilolo,yi0lololo;  // low imagi values of y
   int idx;

   xr0hihihi = xrehihihi[0]; xr0lohihi = xrelohihi[0];
   xr0hilohi = xrehilohi[0]; xr0lolohi = xrelolohi[0];
   xr0hihilo = xrehihilo[0]; xr0lohilo = xrelohilo[0];
   xr0hilolo = xrehilolo[0]; xr0lololo = xrelololo[0];
   xi0hihihi = ximhihihi[0]; xi0lohihi = ximlohihi[0];
   xi0hilohi = ximhilohi[0]; xi0lolohi = ximlolohi[0];
   xi0hihilo = ximhihilo[0]; xi0lohilo = ximlohilo[0];
   xi0hilolo = ximhilolo[0]; xi0lololo = ximlololo[0];
   yr0hihihi = yrehihihi[0]; yr0lohihi = yrelohihi[0];
   yr0hilohi = yrehilohi[0]; yr0lolohi = yrelolohi[0];
   yr0hihilo = yrehihilo[0]; yr0lohilo = yrelohilo[0];
   yr0hilolo = yrehilolo[0]; yr0lololo = yrelololo[0];
   yi0hihihi = yimhihihi[0]; yi0lohihi = yimlohihi[0];
   yi0hilohi = yimhilohi[0]; yi0lolohi = yimlolohi[0];
   yi0hihilo = yimhihilo[0]; yi0lohilo = yimlohilo[0];
   yi0hilolo = yimhilolo[0]; yi0lololo = yimlololo[0];
   // zre[0] = xr0*yr0 - xi0*yi0;
   odf_mul(xr0hihihi,xr0lohihi,xr0hilohi,xr0lolohi,
           xr0hihilo,xr0lohilo,xr0hilolo,xr0lololo,
           yr0hihihi,yr0lohihi,yr0hilohi,yr0lolohi,
           yr0hihilo,yr0lohilo,yr0hilolo,yr0lololo,
           &zrehihihi[0],&zrelohihi[0],&zrehilohi[0],&zrelolohi[0],
           &zrehihilo[0],&zrelohilo[0],&zrehilolo[0],&zrelololo[0]);
   odf_mul(xi0hihihi,xi0lohihi,xi0hilohi,xi0lolohi,
           xi0hihilo,xi0lohilo,xi0hilolo,xi0lololo,
           yi0hihihi,yi0lohihi,yi0hilohi,yi0lolohi,
           yi0hihilo,yi0lohilo,yi0hilolo,yi0lololo,
           &rpahihihi,&rpalohihi,&rpahilohi,&rpalolohi,
           &rpahihilo,&rpalohilo,&rpahilolo,&rpalololo);
   odf_minus(&rpahihihi,&rpalohihi,&rpahilohi,&rpalolohi,
             &rpahihilo,&rpalohilo,&rpahilolo,&rpalololo);
   odf_inc(&zrehihihi[0],&zrelohihi[0],&zrehilohi[0],&zrelolohi[0],
           &zrehihilo[0],&zrelohilo[0],&zrehilolo[0],&zrelololo[0],
           rpahihihi,rpalohihi,rpahilohi,rpalolohi,
           rpahihilo,rpalohilo,rpahilolo,rpalololo);
   // zim[0] = xi0*yr0 + xr0*yi0;
   odf_mul(xi0hihihi,xi0lohihi,xi0hilohi,xi0lolohi,
           xi0hihilo,xi0lohilo,xi0hilolo,xi0lololo,
           yr0hihihi,yr0lohihi,yr0hilohi,yr0lolohi,
           yr0hihilo,yr0lohilo,yr0hilolo,yr0lololo,
           &zimhihihi[0],&zimlohihi[0],&zimhilohi[0],&zimlolohi[0],
           &zimhihilo[0],&zimlohilo[0],&zimhilolo[0],&zimlololo[0]);
   odf_mul(xr0hihihi,xr0lohihi,xr0hilohi,xr0lolohi,
           xr0hihilo,xr0lohilo,xr0hilolo,xr0lololo,
           yi0hihihi,yi0lohihi,yi0hilohi,yi0lolohi,
           yi0hihilo,yi0lohilo,yi0hilolo,yi0lololo,
           &ipahihihi,&ipalohihi,&ipahilohi,&ipalolohi,
           &ipahihilo,&ipalohilo,&ipahilolo,&ipalololo);
   odf_inc(&zimhihihi[0],&zimlohihi[0],&zimhilohi[0],&zimlolohi[0],
           &zimhihilo[0],&zimlohilo[0],&zimhilolo[0],&zimlololo[0],
           ipahihihi,ipalohihi,ipahilohi,ipalolohi,
           ipahihilo,ipalohilo,ipahilolo,ipalololo);

   for(int k=1; k<=deg; k++)
   {
      xr0hihihi = xrehihihi[0]; xr0lohihi = xrelohihi[0];
      xr0hilohi = xrehilohi[0]; xr0lolohi = xrelolohi[0];
      xr0hihilo = xrehihilo[0]; xr0lohilo = xrelohilo[0];
      xr0hilolo = xrehilolo[0]; xr0lololo = xrelololo[0];
      xi0hihihi = ximhihihi[0]; xi0lohihi = ximlohihi[0];
      xi0hilohi = ximhilohi[0]; xi0lolohi = ximlolohi[0];
      xi0hihilo = ximhihilo[0]; xi0lohilo = ximlohilo[0];
      xi0hilolo = ximhilolo[0]; xi0lololo = ximlololo[0];
      yr0hihihi = yrehihihi[k]; yr0lohihi = yrelohihi[k];
      yr0hilohi = yrehilohi[k]; yr0lolohi = yrelolohi[k];
      yr0hihilo = yrehihilo[k]; yr0lohilo = yrelohilo[k];
      yr0hilolo = yrehilolo[k]; yr0lololo = yrelololo[k];
      yi0hihihi = yimhihihi[k]; yi0lohihi = yimlohihi[k];
      yi0hilohi = yimhilohi[k]; yi0lolohi = yimlolohi[k];
      yi0hihilo = yimhihilo[k]; yi0lohilo = yimlohilo[k];
      yi0hilolo = yimhilolo[k]; yi0lololo = yimlololo[k];
      // rpa = xr0*yr0 - xi0*yi0;
      odf_mul(xr0hihihi,xr0lohihi,xr0hilohi,xr0lolohi,
              xr0hihilo,xr0lohilo,xr0hilolo,xr0lololo,
              yr0hihihi,yr0lohihi,yr0hilohi,yr0lolohi,
              yr0hihilo,yr0lohilo,yr0hilolo,yr0lololo,
              &rpahihihi,&rpalohihi,&rpahilohi,&rpalolohi,
              &rpahihilo,&rpalohilo,&rpahilolo,&rpalololo);
      odf_mul(xi0hihihi,xi0lohihi,xi0hilohi,xi0lolohi,
              xi0hihilo,xi0lohilo,xi0hilolo,xi0lololo,
              yi0hihihi,yi0lohihi,yi0hilohi,yi0lolohi,
              yi0hihilo,yi0lohilo,yi0hilolo,yi0lololo,
              &tmphihihi,&tmplohihi,&tmphilohi,&tmplolohi,
              &tmphihilo,&tmplohilo,&tmphilolo,&tmplololo);
      odf_minus(&tmphihihi,&tmplohihi,&tmphilohi,&tmplolohi,
                &tmphihilo,&tmplohilo,&tmphilolo,&tmplololo);
      odf_inc(&rpahihihi,&rpalohihi,&rpahilohi,&rpalolohi,
              &rpahihilo,&rpalohilo,&rpahilolo,&rpalololo,
              tmphihihi,tmplohihi,tmphilohi,tmplolohi,
              tmphihilo,tmplohilo,tmphilolo,tmplololo);
      // ipa = xi0*yr0 + xr0*yi0;
      odf_mul(xi0hihihi,xi0lohihi,xi0hilohi,xi0lolohi,
              xi0hihilo,xi0lohilo,xi0hilolo,xi0lololo,
              yr0hihihi,yr0lohihi,yr0hilohi,yr0lolohi,
              yr0hihilo,yr0lohilo,yr0hilolo,yr0lololo,
              &ipahihihi,&ipalohihi,&ipahilohi,&ipalolohi,
              &ipahihilo,&ipalohilo,&ipahilolo,&ipalololo);
      odf_mul(xr0hihihi,xr0lohihi,xr0hilohi,xr0lolohi,
              xr0hihilo,xr0lohilo,xr0hilolo,xr0lololo,
              yi0hihihi,yi0lohihi,yi0hilohi,yi0lolohi,
              yi0hihilo,yi0lohilo,yi0hilolo,yi0lololo,
              &tmphihihi,&tmplohihi,&tmphilohi,&tmplolohi,
              &tmphihilo,&tmplohilo,&tmphilolo,&tmplololo);
      odf_inc(&ipahihihi,&ipalohihi,&ipahilohi,&ipalolohi,
              &ipahihilo,&ipalohilo,&ipahilolo,&ipalololo,
              tmphihihi,tmplohihi,tmphilohi,tmplolohi,
              tmphihilo,tmplohilo,tmphilolo,tmplololo);

      for(int i=1; i<=k; i++)
      {
         idx = k-i;
         xr0hihihi = xrehihihi[i]; xr0lohihi = xrelohihi[i];
         xr0hilohi = xrehilohi[i]; xr0lolohi = xrelolohi[i];
         xr0hihilo = xrehihilo[i]; xr0lohilo = xrelohilo[i];
         xr0hilolo = xrehilolo[i]; xr0lololo = xrelololo[i];
         xi0hihihi = ximhihihi[i]; xi0lohihi = ximlohihi[i];
         xi0hilohi = ximhilohi[i]; xi0lolohi = ximlolohi[i];
         xi0hihilo = ximhihilo[i]; xi0lohilo = ximlohilo[i];
         xi0hilolo = ximhilolo[i]; xi0lololo = ximlololo[i];
         yr0hihihi = yrehihihi[idx]; yr0lohihi = yrelohihi[idx];
         yr0hilohi = yrehilohi[idx]; yr0lolohi = yrelolohi[idx];
         yr0hihilo = yrehihilo[idx]; yr0lohilo = yrelohilo[idx];
         yr0hilolo = yrehilolo[idx]; yr0lololo = yrelololo[idx];
         yi0hihihi = yimhihihi[idx]; yi0lohihi = yimlohihi[idx];
         yi0hilohi = yimhilohi[idx]; yi0lolohi = yimlolohi[idx];
         yi0hihilo = yimhihilo[idx]; yi0lohilo = yimlohilo[idx];
         yi0hilolo = yimhilolo[idx]; yi0lololo = yimlololo[idx];
         // rpa = rpa + xr0*yr0 - xi0*yi0;
         odf_mul(xr0hihihi,xr0lohihi,xr0hilohi,xr0lolohi,
                 xr0hihilo,xr0lohilo,xr0hilolo,xr0lololo,
                 yr0hihihi,yr0lohihi,yr0hilohi,yr0lolohi,
                 yr0hihilo,yr0lohilo,yr0hilolo,yr0lololo,
                 &tmphihihi,&tmplohihi,&tmphilohi,&tmplolohi,
                 &tmphihilo,&tmplohilo,&tmphilolo,&tmplololo);
         odf_inc(&rpahihihi,&rpalohihi,&rpahilohi,&rpalolohi,
                 &rpahihilo,&rpalohilo,&rpahilolo,&rpalololo,
                 tmphihihi,tmplohihi,tmphilohi,tmplolohi,
                 tmphihilo,tmplohilo,tmphilolo,tmplololo);
         odf_mul(xi0hihihi,xi0lohihi,xi0hilohi,xi0lolohi,
                 xi0hihilo,xi0lohilo,xi0hilolo,xi0lololo,
                 yi0hihihi,yi0lohihi,yi0hilohi,yi0lolohi,
                 yi0hihilo,yi0lohilo,yi0hilolo,yi0lololo,
                 &tmphihihi,&tmplohihi,&tmphilohi,&tmplolohi,
                 &tmphihilo,&tmplohilo,&tmphilolo,&tmplololo);
         odf_minus(&tmphihihi,&tmplohihi,&tmphilohi,&tmplolohi,
                   &tmphihilo,&tmplohilo,&tmphilolo,&tmplololo);
         odf_inc(&rpahihihi,&rpalohihi,&rpahilohi,&rpalolohi,
                 &rpahihilo,&rpalohilo,&rpahilolo,&rpalololo,
                 tmphihihi,tmplohihi,tmphilohi,tmplolohi,
                 tmphihilo,tmplohilo,tmphilolo,tmplololo);
         // ipa = ipa + xi0*yr0 + xr0*yi0;
         odf_mul(xi0hihihi,xi0lohihi,xi0hilohi,xi0lolohi,
                 xi0hihilo,xi0lohilo,xi0hilolo,xi0lololo,
                 yr0hihihi,yr0lohihi,yr0hilohi,yr0lolohi,
                 yr0hihilo,yr0lohilo,yr0hilolo,yr0lololo,
                 &tmphihihi,&tmplohihi,&tmphilohi,&tmplolohi,
                 &tmphihilo,&tmplohilo,&tmphilolo,&tmplololo);
         odf_inc(&ipahihihi,&ipalohihi,&ipahilohi,&ipalolohi,
                 &ipahihilo,&ipalohilo,&ipahilolo,&ipalololo,
                 tmphihihi,tmplohihi,tmphilohi,tmplolohi,
                 tmphihilo,tmplohilo,tmphilolo,tmplololo);
         odf_mul(xr0hihihi,xr0lohihi,xr0hilohi,xr0lolohi,
                 xr0hihilo,xr0lohilo,xr0hilolo,xr0lololo,
                 yi0hihihi,yi0lohihi,yi0hilohi,yi0lolohi,
                 yi0hihilo,yi0lohilo,yi0hilolo,yi0lololo,
                 &tmphihihi,&tmplohihi,&tmphilohi,&tmplolohi,
                 &tmphihilo,&tmplohilo,&tmphilolo,&tmplololo);
         odf_inc(&ipahihihi,&ipalohihi,&ipahilohi,&ipalolohi,
                 &ipahihilo,&ipalohilo,&ipahilolo,&ipalololo,
                 tmphihihi,tmplohihi,tmphilohi,tmplolohi,
                 tmphihilo,tmplohilo,tmphilolo,tmplololo);
      }
      zrehihihi[k] = rpahihihi; zrelohihi[k] = rpalohihi;
      zrehilohi[k] = rpahilohi; zrelolohi[k] = rpalolohi;
      zrehihilo[k] = rpahihilo; zrelohilo[k] = rpalohilo;
      zrehilolo[k] = rpahilolo; zrelololo[k] = rpalololo;
      zimhihihi[k] = ipahihihi; zimlohihi[k] = ipalohihi;
      zimhilohi[k] = ipahilohi; zimlolohi[k] = ipalolohi;
      zimhihilo[k] = ipahihilo; zimlohilo[k] = ipalohilo;
      zimhilolo[k] = ipahilolo; zimlololo[k] = ipalololo;
   }
}
