// The file random8_series.cpp defines functions specified
// in random8_series.h.

#include "octo_double_functions.h"
#include "random8_vectors.h"
#include "random8_series.h"

void dbl8_exponential
 ( int deg, double xhihihi, double xlohihi, double xhilohi, double xlolohi,
            double xhihilo, double xlohilo, double xhilolo, double xlololo,
   double *shihihi, double *slohihi, double *shilohi, double *slolohi,
   double *shihilo, double *slohilo, double *shilolo, double *slololo )
{
   double fhihihi,flohihi,fhilohi,flolohi;
   double fhihilo,flohilo,fhilolo,flololo;

   shihihi[0] = 1.0; slohihi[0] = 0.0;
   shilohi[0] = 0.0; slolohi[0] = 0.0;
   shihilo[0] = 0.0; slohilo[0] = 0.0;
   shilolo[0] = 0.0; slololo[0] = 0.0;
   shihihi[1] = xhihihi; slohihi[1] = xlohihi;
   shilohi[1] = xhilohi; slolohi[1] = xlolohi;
   shihilo[1] = xhihilo; slohilo[1] = xlohilo;
   shilolo[1] = xhilolo; slololo[1] = xlololo;

   for(int k=2; k<=deg; k++)
   {
      odf_mul(shihihi[k-1],slohihi[k-1],shilohi[k-1],slolohi[k-1],
              shihilo[k-1],slohilo[k-1],shilolo[k-1],slololo[k-1],
              xhihihi,xlohihi,xhilohi,xlolohi,
              xhihilo,xlohilo,xhilolo,xlololo,
              &shihihi[k],&slohihi[k],&shilohi[k],&slolohi[k],
              &shihilo[k],&slohilo[k],&shilolo[k],&slololo[k]);
      // x[k] = x[k-1]*r
      fhihihi = (double) k; flohihi = 0.0; fhilohi = 0.0; flolohi = 0.0;
      fhihilo = 0.0; flohilo = 0.0; fhilolo = 0.0; flololo = 0.0;
      odf_div(shihihi[k],slohihi[k],shilohi[k],slolohi[k],
              shihilo[k],slohilo[k],shilolo[k],slololo[k],
              fhihihi,flohihi,fhilohi,flolohi,
              fhihilo,flohilo,fhilolo,flololo,
              &shihihi[k],&slohihi[k],&shilohi[k],&slolohi[k],
              &shihilo[k],&slohilo[k],&shilolo[k],&slololo[k]);
   }
}

void dbl8_exponentials
 ( int deg, double xhihihi, double xlohihi, double xhilohi, double xlolohi, 
            double xhihilo, double xlohilo, double xhilolo, double xlololo, 
   double *pluxhihihi, double *pluxlohihi, double *pluxhilohi,
   double *pluxlolohi, double *pluxhihilo, double *pluxlohilo,
   double *pluxhilolo, double *pluxlololo,
   double *minxhihihi, double *minxlohihi, double *minxhilohi,
   double *minxlolohi, double *minxhihilo, double *minxlohilo,
   double *minxhilolo, double *minxlololo )
{
   double fhihihi,flohihi,fhilohi,flolohi;
   double fhihilo,flohilo,fhilolo,flololo;

   pluxhihihi[0] = 1.0; pluxlohihi[0] = 0.0;
   pluxhilohi[0] = 0.0; pluxlolohi[0] = 0.0;
   pluxhihilo[0] = 0.0; pluxlohilo[0] = 0.0;
   pluxhilolo[0] = 0.0; pluxlololo[0] = 0.0;
   minxhihihi[0] = 1.0; minxlohihi[0] = 0.0;
   minxhilohi[0] = 0.0; minxlolohi[0] = 0.0;
   minxhihilo[0] = 0.0; minxlohilo[0] = 0.0;
   minxhilolo[0] = 0.0; minxlololo[0] = 0.0;
   pluxhihihi[1] = xhihihi; pluxlohihi[1] = xlohihi;
   pluxhilohi[1] = xhilohi; pluxlolohi[1] = xlolohi;
   pluxhihilo[1] = xhihilo; pluxlohilo[1] = xlohilo;
   pluxhilolo[1] = xhilolo; pluxlololo[1] = xlololo;
   minxhihihi[1] = -xhihihi; minxlohihi[1] = -xlohihi;
   minxhilohi[1] = -xhilohi; minxlolohi[1] = -xlolohi;
   minxhihilo[1] = -xhihilo; minxlohilo[1] = -xlohilo;
   minxhilolo[1] = -xhilolo; minxlololo[1] = -xlololo;

   for(int k=2; k<=deg; k++)
   {
      odf_mul(pluxhihihi[k-1],pluxlohihi[k-1],pluxhilohi[k-1],pluxlolohi[k-1],
              pluxhihilo[k-1],pluxlohilo[k-1],pluxhilolo[k-1],pluxlololo[k-1],
              xhihihi,xlohihi,xhilohi,xlolohi,
              xhihilo,xlohilo,xhilolo,xlololo,
              &pluxhihihi[k],&pluxlohihi[k],&pluxhilohi[k],&pluxlolohi[k],
              &pluxhihilo[k],&pluxlohilo[k],&pluxhilolo[k],&pluxlololo[k]);
     // x[k] = x[k-1]*r
      odf_mul(minxhihihi[k-1],minxlohihi[k-1],minxhilohi[k-1],minxlolohi[k-1],
              minxhihilo[k-1],minxlohilo[k-1],minxhilolo[k-1],minxlololo[k-1],
              -xhihihi,-xlohihi,-xhilohi,-xlolohi,
              -xhihilo,-xlohilo,-xhilolo,-xlololo,
              &minxhihihi[k],&minxlohihi[k],&minxhilohi[k],&minxlolohi[k],
              &minxhihilo[k],&minxlohilo[k],&minxhilolo[k],&minxlololo[k]); 
      // y[k] = y[k-1]*(-r)
      fhihihi = (double) k; flohihi = 0.0; fhilohi = 0.0; flolohi = 0.0;
      fhihilo = 0.0; flohilo = 0.0; fhilolo = 0.0; flololo = 0.0;
      odf_div(pluxhihihi[k],pluxlohihi[k],pluxhilohi[k],pluxlolohi[k],
              pluxhihilo[k],pluxlohilo[k],pluxhilolo[k],pluxlololo[k],
              fhihihi,flohihi,fhilohi,flolohi,
              fhihilo,flohilo,fhilolo,flololo,
              &pluxhihihi[k],&pluxlohihi[k],&pluxhilohi[k],&pluxlolohi[k],
              &pluxhihilo[k],&pluxlohilo[k],&pluxhilolo[k],&pluxlololo[k]);
      odf_div(minxhihihi[k],minxlohihi[k],minxhilohi[k],minxlolohi[k],
              minxhihilo[k],minxlohilo[k],minxhilolo[k],minxlololo[k],
              fhihihi,flohihi,fhilohi,flolohi,
              fhihilo,flohilo,fhilolo,flololo,
              &minxhihihi[k],&minxlohihi[k],&minxhilohi[k],&minxlolohi[k],
              &minxhihilo[k],&minxlohilo[k],&minxhilolo[k],&minxlololo[k]);
   }
}

void random_dbl8_exponential
 ( int deg,
   double *xhihihi, double *xlohihi, double *xhilohi, double *xlolohi,
   double *xhihilo, double *xlohilo, double *xhilolo, double *xlololo,
   double *shihihi, double *slohihi, double *shilohi, double *slolohi,
   double *shihilo, double *slohilo, double *shilolo, double *slololo )
{
   random_octo_double
      (xhihihi,xlohihi,xhilohi,xlolohi,xhihilo,xlohilo,xhilolo,xlololo);

   dbl8_exponential
      (deg,*xhihihi,*xlohihi,*xhilohi,*xlolohi,
           *xhihilo,*xlohilo,*xhilolo,*xlololo,
       shihihi,slohihi,shilohi,slolohi,shihilo,slohilo,shilolo,slololo);
}

void random_dbl8_exponentials
 ( int deg,
   double *xhihihi, double *xlohihi, double *xhilohi, double *xlolohi,
   double *xhihilo, double *xlohilo, double *xhilolo, double *xlololo,
   double *pluxhihihi, double *pluxlohihi, double *pluxhilohi,
   double *pluxlolohi, double *pluxhihilo, double *pluxlohilo,
   double *pluxhilolo, double *pluxlololo,
   double *minxhihihi, double *minxlohihi, double *minxhilohi,
   double *minxlolohi, double *minxhihilo, double *minxlohilo,
   double *minxhilolo, double *minxlololo )
{
   random_octo_double
      (xhihihi,xlohihi,xhilohi,xlolohi,xhihilo,xlohilo,xhilolo,xlololo);

   dbl8_exponentials
      (deg,*xhihihi,*xlohihi,*xhilohi,*xlolohi,
           *xhihilo,*xlohilo,*xhilolo,*xlololo,
           pluxhihihi,pluxlohihi,pluxhilohi,pluxlolohi,
           pluxhihilo,pluxlohilo,pluxhilolo,pluxlololo,
           minxhihihi,minxlohihi,minxhilohi,minxlolohi,
           minxhihilo,minxlohilo,minxhilolo,minxlololo);
}

void cmplx8_exponential
 ( int deg,
   double xrehihihi, double xrelohihi, double xrehilohi, double xrelolohi, 
   double xrehihilo, double xrelohilo, double xrehilolo, double xrelololo, 
   double ximhihihi, double ximlohihi, double ximhilohi, double ximlolohi,
   double ximhihilo, double ximlohilo, double ximhilolo, double ximlololo,
   double *srehihihi, double *srelohihi, double *srehilohi, double *srelolohi,
   double *srehihilo, double *srelohilo, double *srehilolo, double *srelololo,
   double *simhihihi, double *simlohihi, double *simhilohi, double *simlolohi,
   double *simhihilo, double *simlohilo, double *simhilolo, double *simlololo )
{
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   srehihihi[0] = 1.0; srelohihi[0] = 0.0;
   srehilohi[0] = 0.0; srelolohi[0] = 0.0;
   srehihilo[0] = 0.0; srelohilo[0] = 0.0;
   srehilolo[0] = 0.0; srelololo[0] = 0.0;
   simhihihi[0] = 0.0; simlohihi[0] = 0.0;
   simhilohi[0] = 0.0; simlolohi[0] = 0.0;
   simhihilo[0] = 0.0; simlohilo[0] = 0.0;
   simhilolo[0] = 0.0; simlololo[0] = 0.0;
   srehihihi[1] = xrehihihi; srelohihi[1] = xrelohihi;
   srehilohi[1] = xrehilohi; srelolohi[1] = xrelolohi;
   srehihilo[1] = xrehihilo; srelohilo[1] = xrelohilo;
   srehilolo[1] = xrehilolo; srelololo[1] = xrelololo;
   simhihihi[1] = ximhihihi; simlohihi[1] = ximlohihi;
   simhilohi[1] = ximhilohi; simlolohi[1] = ximlolohi;
   simhihilo[1] = ximhihilo; simlohilo[1] = ximlohilo;
   simhilolo[1] = ximhilolo; simlololo[1] = ximlololo;

   for(int k=2; k<=deg; k++)
   {
      // sre[k] = (sre[k-1]*xre - sim[k-1]*xim)/k;
      odf_mul(srehihihi[k-1],srelohihi[k-1],srehilohi[k-1],srelolohi[k-1],
              srehihilo[k-1],srelohilo[k-1],srehilolo[k-1],srelololo[k-1],
              xrehihihi,     xrelohihi,     xrehilohi,     xrelolohi,
              xrehihilo,     xrelohilo,     xrehilolo,     xrelololo,
             &srehihihi[k], &srelohihi[k], &srehilohi[k], &srelolohi[k],
             &srehihilo[k], &srelohilo[k], &srehilolo[k], &srelololo[k]);
      odf_mul(simhihihi[k-1],simlohihi[k-1],simhilohi[k-1],simlolohi[k-1],
              simhihilo[k-1],simlohilo[k-1],simhilolo[k-1],simlololo[k-1],
              ximhihihi,     ximlohihi,     ximhilohi,     ximlolohi,
              ximhihilo,     ximlohilo,     ximhilolo,     ximlololo,
             &acchihihi,    &acclohihi,    &acchilohi,    &acclolohi,
             &acchihilo,    &acclohilo,    &acchilolo,    &acclololo);
      odf_minus(&acchihihi,&acclohihi,&acchilohi,&acclolohi,
                &acchihilo,&acclohilo,&acchilolo,&acclololo);
      odf_inc(&srehihihi[k],&srelohihi[k],&srehilohi[k],&srelolohi[k],
              &srehihilo[k],&srelohilo[k],&srehilolo[k],&srelololo[k],
               acchihihi,    acclohihi,    acchilohi,    acclolohi,
               acchihilo,    acclohilo,    acchilolo,    acclololo);
      acchihihi = (double) k;
                       acclohihi = 0.0; acchilohi = 0.0; acclolohi = 0.0;
      acchihilo = 0.0; acclohilo = 0.0; acchilolo = 0.0; acclololo = 0.0;
      odf_div(srehihihi[k], srelohihi[k], srehilohi[k], srelolohi[k],
              srehihilo[k], srelohilo[k], srehilolo[k], srelololo[k],
              acchihihi,    acclohihi,    acchilohi,    acclolohi,
              acchihilo,    acclohilo,    acchilolo,    acclololo,
             &srehihihi[k],&srelohihi[k],&srehilohi[k],&srelolohi[k],
             &srehihilo[k],&srelohilo[k],&srehilolo[k],&srelololo[k]);
      // sim[k] = (sre[k-1]*xim + sim[k-1]*xre)/k;
      odf_mul(srehihihi[k-1],srelohihi[k-1],srehilohi[k-1],srelolohi[k-1],
              srehihilo[k-1],srelohilo[k-1],srehilolo[k-1],srelololo[k-1],
              ximhihihi,     ximlohihi,     ximhilohi,     ximlolohi,
              ximhihilo,     ximlohilo,     ximhilolo,     ximlololo,
             &simhihihi[k], &simlohihi[k], &simhilohi[k], &simlolohi[k],
             &simhihilo[k], &simlohilo[k], &simhilolo[k], &simlololo[k]);
      odf_mul(simhihihi[k-1],simlohihi[k-1],simhilohi[k-1],simlolohi[k-1],
              simhihilo[k-1],simlohilo[k-1],simhilolo[k-1],simlololo[k-1],
              xrehihihi,     xrelohihi,     xrehilohi,     xrelolohi,
              xrehihilo,     xrelohilo,     xrehilolo,     xrelololo,
             &acchihihi,    &acclohihi,    &acchilohi,    &acclolohi,
             &acchihilo,    &acclohilo,    &acchilolo,    &acclololo);
      odf_inc(&simhihihi[k],&simlohihi[k],&simhilohi[k],&simlolohi[k],
              &simhihilo[k],&simlohilo[k],&simhilolo[k],&simlololo[k],
               acchihihi,    acclohihi,    acchilohi,    acclolohi,
               acchihilo,    acclohilo,    acchilolo,    acclololo);
      acchihihi = (double) k;
                       acclohihi = 0.0; acchilohi = 0.0; acclolohi = 0.0;
      acchihilo = 0.0; acclohilo = 0.0; acchilolo = 0.0; acclololo = 0.0;
      odf_div(simhihihi[k], simlohihi[k], simhilohi[k], simlolohi[k],
              simhihilo[k], simlohilo[k], simhilolo[k], simlololo[k],
              acchihihi,    acclohihi,    acchilohi,    acclolohi,
              acchihilo,    acclohilo,    acchilolo,    acclololo,
             &simhihihi[k],&simlohihi[k],&simhilohi[k],&simlolohi[k],
             &simhihilo[k],&simlohilo[k],&simhilolo[k],&simlololo[k]);
   }
}

void cmplx8_exponentials
 ( int deg,
   double xrehihihi, double xrelohihi, double xrehilohi, double xrelolohi,
   double xrehihilo, double xrelohilo, double xrehilolo, double xrelololo,
   double ximhihihi, double ximlohihi, double ximhilohi, double ximlolohi,
   double ximhihilo, double ximlohilo, double ximhilolo, double ximlololo,
   double *pluxrehihihi, double *pluxrelohihi,
   double *pluxrehilohi, double *pluxrelolohi,
   double *pluxrehihilo, double *pluxrelohilo,
   double *pluxrehilolo, double *pluxrelololo,
   double *pluximhihihi, double *pluximlohihi,
   double *pluximhilohi, double *pluximlolohi,
   double *pluximhihilo, double *pluximlohilo,
   double *pluximhilolo, double *pluximlololo,
   double *minxrehihihi, double *minxrelohihi,
   double *minxrehilohi, double *minxrelolohi,
   double *minxrehihilo, double *minxrelohilo,
   double *minxrehilolo, double *minxrelololo,
   double *minximhihihi, double *minximlohihi,
   double *minximhilohi, double *minximlolohi,
   double *minximhihilo, double *minximlohilo,
   double *minximhilolo, double *minximlololo )
{
   double tmphihihi,tmplohihi,tmphilohi,tmplolohi;
   double tmphihilo,tmplohilo,tmphilolo,tmplololo;

   pluxrehihihi[0] = 1.0; pluxrelohihi[0] = 0.0;
   pluxrehilohi[0] = 0.0; pluxrelolohi[0] = 0.0;
   pluxrehihilo[0] = 0.0; pluxrelohilo[0] = 0.0;
   pluxrehilolo[0] = 0.0; pluxrelololo[0] = 0.0;
   minxrehihihi[0] = 1.0; minxrelohihi[0] = 0.0;
   minxrehilohi[0] = 0.0; minxrelolohi[0] = 0.0;
   minxrehihilo[0] = 0.0; minxrelohilo[0] = 0.0;
   minxrehilolo[0] = 0.0; minxrelololo[0] = 0.0;
   pluximhihihi[0] = 0.0; pluximlohihi[0] = 0.0;
   pluximhilohi[0] = 0.0; pluximlolohi[0] = 0.0;
   pluximhihilo[0] = 0.0; pluximlohilo[0] = 0.0;
   pluximhilolo[0] = 0.0; pluximlololo[0] = 0.0;
   minximhihihi[0] = 0.0; minximlohihi[0] = 0.0;
   minximhilohi[0] = 0.0; minximlolohi[0] = 0.0;
   minximhihilo[0] = 0.0; minximlohilo[0] = 0.0;
   minximhilolo[0] = 0.0; minximlololo[0] = 0.0;
   pluxrehihihi[1] = xrehihihi; pluxrelohihi[1] = xrelohihi;
   pluxrehilohi[1] = xrehilohi; pluxrelolohi[1] = xrelolohi;
   pluxrehihilo[1] = xrehihilo; pluxrelohilo[1] = xrelohilo;
   pluxrehilolo[1] = xrehilolo; pluxrelololo[1] = xrelololo;
   pluximhihihi[1] = ximhihihi; pluximlohihi[1] = ximlohihi;
   pluximhilohi[1] = ximhilohi; pluximlolohi[1] = ximlolohi;
   pluximhihilo[1] = ximhihilo; pluximlohilo[1] = ximlohilo;
   pluximhilolo[1] = ximhilolo; pluximlololo[1] = ximlololo;
   minxrehihihi[1] = -xrehihihi; minxrelohihi[1] = -xrelohihi;
   minxrehilohi[1] = -xrehilohi; minxrelolohi[1] = -xrelolohi;
   minxrehihilo[1] = -xrehihilo; minxrelohilo[1] = -xrelohilo;
   minxrehilolo[1] = -xrehilolo; minxrelololo[1] = -xrelololo;
   minximhihihi[1] = -ximhihihi; minximlohihi[1] = -ximlohihi;
   minximhilohi[1] = -ximhilohi; minximlolohi[1] = -ximlolohi;
   minximhihilo[1] = -ximhihilo; minximlohilo[1] = -ximlohilo;
   minximhilolo[1] = -ximhilolo; minximlololo[1] = -ximlololo;

   for(int k=2; k<=deg; k++)
   {
      // xre[k] = (xre[k-1]*cr - xim[k-1]*sr)/k;
      odf_mul(pluxrehihihi[k-1],pluxrelohihi[k-1],pluxrehilohi[k-1],
              pluxrelolohi[k-1],pluxrehihilo[k-1],pluxrelohilo[k-1],
              pluxrehilolo[k-1],pluxrelololo[k-1],
              xrehihihi,xrelohihi,xrehilohi,xrelolohi,
              xrehihilo,xrelohilo,xrehilolo,xrelololo,
              &pluxrehihihi[k],&pluxrelohihi[k],&pluxrehilohi[k],
              &pluxrelolohi[k],&pluxrehihilo[k],&pluxrelohilo[k],
              &pluxrehilolo[k],&pluxrelololo[k]);
      odf_mul(pluximhihihi[k-1],pluximlohihi[k-1],pluximhilohi[k-1],
              pluximlolohi[k-1],pluximhihilo[k-1],pluximlohilo[k-1],
              pluximhilolo[k-1],pluximlololo[k-1],
              ximhihihi,ximlohihi,ximhilohi,ximlolohi,
              ximhihilo,ximlohilo,ximhilolo,ximlololo,
              &tmphihihi,&tmplohihi,&tmphilohi,&tmplolohi,
              &tmphihilo,&tmplohilo,&tmphilolo,&tmplololo);
      odf_minus(&tmphihihi,&tmplohihi,&tmphilohi,&tmplolohi,
                &tmphihilo,&tmplohilo,&tmphilolo,&tmplololo);
      odf_inc(&pluxrehihihi[k],&pluxrelohihi[k],&pluxrehilohi[k],
              &pluxrelolohi[k],&pluxrehihilo[k],&pluxrelohilo[k],
              &pluxrehilolo[k],&pluxrelololo[k],
              tmphihihi,tmplohihi,tmphilohi,tmplolohi,
              tmphihilo,tmplohilo,tmphilolo,tmplololo);
      tmphihihi = (double) k; tmplohihi = 0.0; tmphilohi = 0.0;
      tmplolohi = 0.0; tmphihilo = 0.0; tmplohilo = 0.0;
      tmphilolo = 0.0; tmplololo = 0.0;
      odf_div(pluxrehihihi[k],pluxrelohihi[k],pluxrehilohi[k],pluxrelolohi[k],
              pluxrehihilo[k],pluxrelohilo[k],pluxrehilolo[k],pluxrelololo[k],
              tmphihihi,tmplohihi,tmphilohi,tmplolohi,
              tmphihilo,tmplohilo,tmphilolo,tmplololo,
              &pluxrehihihi[k],&pluxrelohihi[k],&pluxrehilohi[k],
              &pluxrelolohi[k],&pluxrehihilo[k],&pluxrelohilo[k],
              &pluxrehilolo[k],&pluxrelololo[k]);
      // xim[k] = (xre[k-1]*sr + xim[k-1]*cr)/k;
      odf_mul(pluxrehihihi[k-1],pluxrelohihi[k-1],pluxrehilohi[k-1],
              pluxrelolohi[k-1],pluxrehihilo[k-1],pluxrelohilo[k-1],
              pluxrehilolo[k-1],pluxrelololo[k-1],
              ximhihihi,ximlohihi,ximhilohi,ximlolohi,
              ximhihilo,ximlohilo,ximhilolo,ximlololo,
              &pluximhihihi[k],&pluximlohihi[k],&pluximhilohi[k],
              &pluximlolohi[k],&pluximhihilo[k],&pluximlohilo[k],
              &pluximhilolo[k],&pluximlololo[k]);
      odf_mul(pluximhihihi[k-1],pluximlohihi[k-1],pluximhilohi[k-1],
              pluximlolohi[k-1],pluximhihilo[k-1],pluximlohilo[k-1],
              pluximhilolo[k-1],pluximlololo[k-1],
              xrehihihi,xrelohihi,xrehilohi,xrelolohi,
              xrehihilo,xrelohilo,xrehilolo,xrelololo,
              &tmphihihi,&tmplohihi,&tmphilohi,&tmplolohi,
              &tmphihilo,&tmplohilo,&tmphilolo,&tmplololo);
      odf_inc(&pluximhihihi[k],&pluximlohihi[k],&pluximhilohi[k],
              &pluximlolohi[k],&pluximhihilo[k],&pluximlohilo[k],
              &pluximhilolo[k],&pluximlololo[k],
              tmphihihi,tmplohihi,tmphilohi,tmplolohi,
              tmphihilo,tmplohilo,tmphilolo,tmplololo);
      tmphihihi = (double) k; tmplohihi = 0.0; tmphilohi = 0.0;
      tmplolohi = 0.0; tmphihilo = 0.0; tmplohilo = 0.0;
      tmphilolo = 0.0; tmplololo = 0.0;
      odf_div(pluximhihihi[k],pluximlohihi[k],pluximhilohi[k],pluximlolohi[k],
              pluximhihilo[k],pluximlohilo[k],pluximhilolo[k],pluximlololo[k],
              tmphihihi,tmplohihi,tmphilohi,tmplolohi,
              tmphihilo,tmplohilo,tmphilolo,tmplololo,
              &pluximhihihi[k],&pluximlohihi[k],&pluximhilohi[k],
              &pluximlolohi[k],&pluximhihilo[k],&pluximlohilo[k],
              &pluximhilolo[k],&pluximlololo[k]);
      // yre[k] = (yre[k-1]*(-cr) - yim[k-1]*(-sr))/k;
      odf_mul(minxrehihihi[k-1],minxrelohihi[k-1],minxrehilohi[k-1],
              minxrelolohi[k-1],minxrehihilo[k-1],minxrelohilo[k-1],
              minxrehilolo[k-1],minxrelololo[k-1],
              -xrehihihi,-xrelohihi,-xrehilohi,-xrelolohi,
              -xrehihilo,-xrelohilo,-xrehilolo,-xrelololo,
              &minxrehihihi[k],&minxrelohihi[k],&minxrehilohi[k],
              &minxrelolohi[k],&minxrehihilo[k],&minxrelohilo[k],
              &minxrehilolo[k],&minxrelololo[k]);
      odf_mul(minximhihihi[k-1],minximlohihi[k-1],minximhilohi[k-1],
              minximlolohi[k-1],minximhihilo[k-1],minximlohilo[k-1],
              minximhilolo[k-1],minximlololo[k-1],
              -ximhihihi,-ximlohihi,-ximhilohi,-ximlolohi,
              -ximhihilo,-ximlohilo,-ximhilolo,-ximlololo,
              &tmphihihi,&tmplohihi,&tmphilohi,&tmplolohi,
              &tmphihilo,&tmplohilo,&tmphilolo,&tmplololo);
      odf_minus(&tmphihihi,&tmplohihi,&tmphilohi,&tmplolohi,
                &tmphihilo,&tmplohilo,&tmphilolo,&tmplololo);
      odf_inc(&minxrehihihi[k],&minxrelohihi[k],&minxrehilohi[k],
              &minxrelolohi[k],&minxrehihilo[k],&minxrelohilo[k],
              &minxrehilolo[k],&minxrelololo[k],
              tmphihihi,tmplohihi,tmphilohi,tmplolohi,
              tmphihilo,tmplohilo,tmphilolo,tmplololo);
      tmphihihi = (double) k; tmplohihi = 0.0; tmphilohi = 0.0;
      tmplolohi = 0.0; tmphihilo = 0.0; tmplohilo = 0.0;
      tmphilolo = 0.0; tmplololo = 0.0;
      odf_div(minxrehihihi[k],minxrelohihi[k],minxrehilohi[k],minxrelolohi[k],
              minxrehihilo[k],minxrelohilo[k],minxrehilolo[k],minxrelololo[k],
              tmphihihi,tmplohihi,tmphilohi,tmplolohi,
              tmphihilo,tmplohilo,tmphilolo,tmplololo,
              &minxrehihihi[k],&minxrelohihi[k],&minxrehilohi[k],
              &minxrelolohi[k],&minxrehihilo[k],&minxrelohilo[k],
              &minxrehilolo[k],&minxrelololo[k]);
      // yim[k] = (yre[k-1]*(-sr) + yim[k-1]*(-cr))/k;
      odf_mul(minxrehihihi[k-1],minxrelohihi[k-1],minxrehilohi[k-1],
              minxrelolohi[k-1],minxrehihilo[k-1],minxrelohilo[k-1],
              minxrehilolo[k-1],minxrelololo[k-1],
              -ximhihihi,-ximlohihi,-ximhilohi,-ximlolohi,
              -ximhihilo,-ximlohilo,-ximhilolo,-ximlololo,
              &minximhihihi[k],&minximlohihi[k],&minximhilohi[k],
              &minximlolohi[k],&minximhihilo[k],&minximlohilo[k],
              &minximhilolo[k],&minximlololo[k]);
      odf_mul(minximhihihi[k-1],minximlohihi[k-1],minximhilohi[k-1],
              minximlolohi[k-1],minximhihilo[k-1],minximlohilo[k-1],
              minximhilolo[k-1],minximlololo[k-1],
              -xrehihihi,-xrelohihi,-xrehilohi,-xrelolohi,
              -xrehihilo,-xrelohilo,-xrehilolo,-xrelololo,
              &tmphihihi,&tmplohihi,&tmphilohi,&tmplolohi,
              &tmphihilo,&tmplohilo,&tmphilolo,&tmplololo);
      odf_inc(&minximhihihi[k],&minximlohihi[k],&minximhilohi[k],
              &minximlolohi[k],&minximhihilo[k],&minximlohilo[k],
              &minximhilolo[k],&minximlololo[k],
              tmphihihi,tmplohihi,tmphilohi,tmplolohi,
              tmphihilo,tmplohilo,tmphilolo,tmplololo);
      tmphihihi = (double) k; tmplohihi = 0.0; tmphilohi = 0.0;
      tmplolohi = 0.0; tmphihilo = 0.0; tmplohilo = 0.0;
      tmphilolo = 0.0; tmplololo = 0.0;
      odf_div(minximhihihi[k],minximlohihi[k],minximhilohi[k],minximlolohi[k],
              minximhihilo[k],minximlohilo[k],minximhilolo[k],minximlololo[k],
              tmphihihi,tmplohihi,tmphilohi,tmplolohi,
              tmphihilo,tmplohilo,tmphilolo,tmplololo,
              &minximhihihi[k],&minximlohihi[k],&minximhilohi[k],
              &minximlolohi[k],&minximhihilo[k],&minximlohilo[k],
              &minximhilolo[k],&minximlololo[k]);
   }
}

void random_cmplx8_exponential
 ( int deg,
   double *xrehihihi, double *xrelohihi, double *xrehilohi, double *xrelolohi,
   double *xrehihilo, double *xrelohilo, double *xrehilolo, double *xrelololo,
   double *ximhihihi, double *ximlohihi, double *ximhilohi, double *ximlolohi,
   double *ximhihilo, double *ximlohilo, double *ximhilolo, double *ximlololo,
   double *srehihihi, double *srelohihi, double *srehilohi, double *srelolohi,
   double *srehihilo, double *srelohilo, double *srehilolo, double *srelololo,
   double *simhihihi, double *simlohihi, double *simhilohi, double *simlolohi,
   double *simhihilo, double *simlohilo, double *simhilolo, double *simlololo )
{
   double tmphihihi,tmplohihi,tmphilohi,tmplolohi;
   double tmphihilo,tmplohilo,tmphilolo,tmplololo;

   random_octo_double
      (xrehihihi,xrelohihi,xrehilohi,xrelolohi,
       xrehihilo,xrelohilo,xrehilolo,xrelololo);               // cos(a)

   odf_sqr(*xrehihihi,*xrelohihi,*xrehilohi,*xrelolohi,
           *xrehihilo,*xrelohilo,*xrehilolo,*xrelololo,
           &tmphihihi,&tmplohihi,&tmphilohi,&tmplolohi,
           &tmphihilo,&tmplohilo,&tmphilolo,&tmplololo);       // cos^2(a)
   odf_minus(&tmphihihi,&tmplohihi,&tmphilohi,&tmplolohi,
             &tmphihilo,&tmplohilo,&tmphilolo,&tmplololo);     // -cos^2(a)
   odf_inc_d(&tmphihihi,&tmplohihi,&tmphilohi,&tmplolohi,
             &tmphihilo,&tmplohilo,&tmphilolo,&tmplololo,1.0); // 1-cos^2(a)
   odf_sqrt(tmphihihi,tmplohihi,tmphilohi,tmplolohi,
            tmphihilo,tmplohilo,tmphilolo,tmplololo,
            ximhihihi,ximlohihi,ximhilohi,ximlolohi,
            ximhihilo,ximlohilo,ximhilolo,ximlololo);          // sin is sqrt

   cmplx8_exponential
      (deg,*xrehihihi,*xrelohihi,*xrehilohi,*xrelolohi,
           *xrehihilo,*xrelohilo,*xrehilolo,*xrelololo,
           *ximhihihi,*ximlohihi,*ximhilohi,*ximlolohi,
           *ximhihilo,*ximlohilo,*ximhilolo,*ximlololo,
       srehihihi,srelohihi,srehilohi,srelolohi,
       srehihilo,srelohilo,srehilolo,srelololo,
       simhihihi,simlohihi,simhilohi,simlolohi,
       simhihilo,simlohilo,simhilolo,simlololo);
}

void random_cmplx8_exponentials
 ( int deg, 
   double *xrehihihi, double *xrelohihi, double *xrehilohi, double *xrelolohi,
   double *xrehihilo, double *xrelohilo, double *xrehilolo, double *xrelololo,
   double *ximhihihi, double *ximlohihi, double *ximhilohi, double *ximlolohi,
   double *ximhihilo, double *ximlohilo, double *ximhilolo, double *ximlololo,
   double *pluxrehihihi, double *pluxrelohihi,
   double *pluxrehilohi, double *pluxrelolohi,
   double *pluxrehihilo, double *pluxrelohilo,
   double *pluxrehilolo, double *pluxrelololo,
   double *pluximhihihi, double *pluximlohihi,
   double *pluximhilohi, double *pluximlolohi,
   double *pluximhihilo, double *pluximlohilo,
   double *pluximhilolo, double *pluximlololo,
   double *minxrehihihi, double *minxrelohihi,
   double *minxrehilohi, double *minxrelolohi,
   double *minxrehihilo, double *minxrelohilo,
   double *minxrehilolo, double *minxrelololo,
   double *minximhihihi, double *minximlohihi,
   double *minximhilohi, double *minximlolohi,
   double *minximhihilo, double *minximlohilo,
   double *minximhilolo, double *minximlololo )
{
   double tmphihihi,tmplohihi,tmphilohi,tmplolohi;
   double tmphihilo,tmplohilo,tmphilolo,tmplololo;

   random_octo_double
      (xrehihihi,xrelohihi,xrehilohi,xrelolohi,
       xrehihilo,xrelohilo,xrehilolo,xrelololo);               // cos(a)

   odf_sqr(*xrehihihi,*xrelohihi,*xrehilohi,*xrelolohi,
           *xrehihilo,*xrelohilo,*xrehilolo,*xrelololo,
           &tmphihihi,&tmplohihi,&tmphilohi,&tmplolohi,
           &tmphihilo,&tmplohilo,&tmphilolo,&tmplololo);       // cos^2(a)
   odf_minus(&tmphihihi,&tmplohihi,&tmphilohi,&tmplolohi,
             &tmphihilo,&tmplohilo,&tmphilolo,&tmplololo);     // -cos^2(a)
   odf_inc_d(&tmphihihi,&tmplohihi,&tmphilohi,&tmplolohi,
             &tmphihilo,&tmplohilo,&tmphilolo,&tmplololo,1.0); // 1-cos^2(a)
   odf_sqrt(tmphihihi,tmplohihi,tmphilohi,tmplolohi,
            tmphihilo,tmplohilo,tmphilolo,tmplololo,
            ximhihihi,ximlohihi,ximhilohi,ximlolohi,
            ximhihilo,ximlohilo,ximhilolo,ximlololo);          // sin is sqrt

   cmplx8_exponentials
      (deg,*xrehihihi,*xrelohihi,*xrehilohi,*xrelolohi,
           *xrehihilo,*xrelohilo,*xrehilolo,*xrelololo,
           *ximhihihi,*ximlohihi,*ximhilohi,*ximlolohi,
           *ximhihilo,*ximlohilo,*ximhilolo,*ximlololo,
           pluxrehihihi,pluxrelohihi,pluxrehilohi,pluxrelolohi,
           pluxrehihilo,pluxrelohilo,pluxrehilolo,pluxrelololo,
           pluximhihihi,pluximlohihi,pluximhilohi,pluximlolohi,
           pluximhihilo,pluximlohilo,pluximhilolo,pluximlololo,
           minxrehihihi,minxrelohihi,minxrehilohi,minxrelolohi,
           minxrehihilo,minxrelohilo,minxrehilolo,minxrelololo,
           minximhihihi,minximlohihi,minximhilohi,minximlolohi,
           minximhihilo,minximlohilo,minximhilolo,minximlololo);
}
