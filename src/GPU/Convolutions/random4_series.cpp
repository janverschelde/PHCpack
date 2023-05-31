// The file random4_series.cpp defines functions specified
// in random4_series.h.

#include "random4_vectors.h"
#include "quad_double_functions.h"

void dbl4_exponential
 ( int deg, double xhihi, double xlohi, double xhilo, double xlolo,
   double *shihi, double *slohi, double *shilo, double *slolo )
{
   double fhihi,flohi,fhilo,flolo;

   shihi[0] = 1.0; slohi[0] = 0.0; shilo[0] = 0.0; slolo[0] = 0.0;
   shihi[1] = xhihi; slohi[1] = xlohi;
   shilo[1] = xhilo; slolo[1] = xlolo;

   for(int k=2; k<=deg; k++)
   {
      qdf_mul(shihi[k-1],slohi[k-1],shilo[k-1],slolo[k-1],
              xhihi,xlohi,xhilo,xlolo,
              &shihi[k],&slohi[k],&shilo[k],&slolo[k]); 
      // x[k] = x[k-1]*r;
      fhihi = (double) k; flohi = 0.0; fhilo = 0.0; flolo = 0.0;
      qdf_div(shihi[k],slohi[k],shilo[k],slolo[k],
              fhihi,flohi,fhilo,flolo,
              &shihi[k],&slohi[k],&shilo[k],&slolo[k]);
   }
}

void dbl4_exponentials
 ( int deg, double xhihi, double xlohi, double xhilo, double xlolo, 
   double *pluxhihi, double *pluxlohi, double *pluxhilo, double *pluxlolo,
   double *minxhihi, double *minxlohi, double *minxhilo, double *minxlolo )
{
   double fhihi,flohi,fhilo,flolo;

   pluxhihi[0] = 1.0; pluxlohi[0] = 0.0; pluxhilo[0] = 0.0; pluxlolo[0] = 0.0;
   minxhihi[0] = 1.0; minxlohi[0] = 0.0; minxhilo[0] = 0.0; minxlolo[0] = 0.0;
   pluxhihi[1] = xhihi; pluxlohi[1] = xlohi;
   pluxhilo[1] = xhilo; pluxlolo[1] = xlolo;
   minxhihi[1] = -xhihi; minxlohi[1] = -xlohi;
   minxhilo[1] = -xhilo; minxlolo[1] = -xlolo;

   for(int k=2; k<=deg; k++)
   {
      qdf_mul(pluxhihi[k-1],pluxlohi[k-1],pluxhilo[k-1],pluxlolo[k-1],
              xhihi,xlohi,xhilo,xlolo,
              &pluxhihi[k],&pluxlohi[k],&pluxhilo[k],&pluxlolo[k]); 
      // x[k] = x[k-1]*r;
      qdf_mul(minxhihi[k-1],minxlohi[k-1],minxhilo[k-1],minxlolo[k-1],
              -xhihi,-xlohi,-xhilo,-xlolo,
              &minxhihi[k],&minxlohi[k],&minxhilo[k],&minxlolo[k]); 
      // y[k] = y[k-1]*(-r);
      fhihi = (double) k; flohi = 0.0; fhilo = 0.0; flolo = 0.0;
      qdf_div(pluxhihi[k],pluxlohi[k],pluxhilo[k],pluxlolo[k],
              fhihi,flohi,fhilo,flolo,
              &pluxhihi[k],&pluxlohi[k],&pluxhilo[k],&pluxlolo[k]);
      qdf_div(minxhihi[k],minxlohi[k],minxhilo[k],minxlolo[k],
              fhihi,flohi,fhilo,flolo,
              &minxhihi[k],&minxlohi[k],&minxhilo[k],&minxlolo[k]);
   }
}

void random_dbl4_exponential
 ( int deg, double *xhihi, double *xlohi, double *xhilo, double *xlolo,
   double *shihi, double *slohi, double *shilo, double *slolo )
{
   random_quad_double(xhihi,xlohi,xhilo,xlolo);

   dbl4_exponential
      (deg,*xhihi,*xlohi,*xhilo,*xlolo,shihi,slohi,shilo,slolo);
}

void random_dbl4_exponentials
 ( int deg, double *xhihi, double *xlohi, double *xhilo, double *xlolo,
   double *pluxhihi, double *pluxlohi, double *pluxhilo, double *pluxlolo,
   double *minxhihi, double *minxlohi, double *minxhilo, double *minxlolo )
{
   random_quad_double(xhihi,xlohi,xhilo,xlolo);

   dbl4_exponentials
      (deg,*xhihi,*xlohi,*xhilo,*xlolo,
       pluxhihi,pluxlohi,pluxhilo,pluxlolo,
       minxhihi,minxlohi,minxhilo,minxlolo);
}

void cmplx4_exponential
 ( int deg, double xrehihi, double xrelohi, double xrehilo, double xrelolo, 
            double ximhihi, double ximlohi, double ximhilo, double ximlolo,
   double *srehihi, double *srelohi, double *srehilo, double *srelolo,
   double *simhihi, double *simlohi, double *simhilo, double *simlolo )
{
   double acchihi,acclohi,acchilo,acclolo;

   srehihi[0] = 1.0; srelohi[0] = 0.0;
   srehilo[0] = 0.0; srelolo[0] = 0.0;
   simhihi[0] = 0.0; simlohi[0] = 0.0;
   simhilo[0] = 0.0; simlolo[0] = 0.0;
   srehihi[1] = xrehihi; srelohi[1] = xrelohi;
   srehilo[1] = xrehilo; srelolo[1] = xrelolo;
   simhihi[1] = ximhihi; simlohi[1] = ximlohi;
   simhilo[1] = ximhilo; simlolo[1] = ximlolo;

   for(int k=2; k<=deg; k++)
   {
      // sre[k] = (sre[k-1]*xre - sim[k-1]*xim)/k;
      qdf_mul(srehihi[k-1],srelohi[k-1],srehilo[k-1],srelolo[k-1],
              xrehihi,     xrelohi,     xrehilo,     xrelolo,
             &srehihi[k], &srelohi[k], &srehilo[k], &srelolo[k]);
      qdf_mul(simhihi[k-1],simlohi[k-1],simhilo[k-1],simlolo[k-1],
              ximhihi,     ximlohi,     ximhilo,     ximlolo,
             &acchihi,    &acclohi,    &acchilo,    &acclolo);
      qdf_minus(&acchihi,&acclohi,&acchilo,&acclolo);
      qdf_inc(&srehihi[k],&srelohi[k],&srehilo[k],&srelolo[k],
               acchihi,    acclohi,    acchilo,    acclolo);
      acchihi = (double) k; acclohi = 0.0; acchilo = 0.0; acclolo = 0.0;
      qdf_div(srehihi[k], srelohi[k], srehilo[k], srelolo[k],
              acchihi,    acclohi,    acchilo,    acclolo,
             &srehihi[k],&srelohi[k],&srehilo[k],&srelolo[k]);
      // sim[k] = (sre[k-1]*xim + sim[k-1]*xre)/k;
      qdf_mul(srehihi[k-1],srelohi[k-1],srehilo[k-1],srelolo[k-1],
              ximhihi,     ximlohi,     ximhilo,     ximlolo,
             &simhihi[k], &simlohi[k], &simhilo[k], &simlolo[k]);
      qdf_mul(simhihi[k-1],simlohi[k-1],simhilo[k-1],simlolo[k-1],
              xrehihi,     xrelohi,     xrehilo,     xrelolo,
             &acchihi,    &acclohi,    &acchilo,    &acclolo);
      qdf_inc(&simhihi[k],&simlohi[k],&simhilo[k],&simlolo[k],
               acchihi,    acclohi,    acchilo,    acclolo);
      acchihi = (double) k; acclohi = 0.0; acchilo = 0.0; acclolo = 0.0;
      qdf_div(simhihi[k], simlohi[k], simhilo[k], simlolo[k],
              acchihi,    acclohi,    acchilo,    acclolo,
             &simhihi[k],&simlohi[k],&simhilo[k],&simlolo[k]);
   }
}

void cmplx4_exponentials
 ( int deg, double xrehihi, double xrelohi, double xrehilo, double xrelolo,
            double ximhihi, double ximlohi, double ximhilo, double ximlolo,
   double *pluxrehihi, double *pluxrelohi, double *pluxrehilo,
   double *pluxrelolo, double *pluximhihi, double *pluximlohi,
   double *pluximhilo, double *pluximlolo,
   double *minxrehihi, double *minxrelohi, double *minxrehilo,
   double *minxrelolo, double *minximhihi, double *minximlohi,
   double *minximhilo, double *minximlolo )
{
   double tmphihi,tmplohi,tmphilo,tmplolo;

   pluxrehihi[0] = 1.0; pluxrelohi[0] = 0.0;
   pluxrehilo[0] = 0.0; pluxrelolo[0] = 0.0;
   minxrehihi[0] = 1.0; minxrelohi[0] = 0.0;
   minxrehilo[0] = 0.0; minxrelolo[0] = 0.0;
   pluximhihi[0] = 0.0; pluximlohi[0] = 0.0;
   pluximhilo[0] = 0.0; pluximlolo[0] = 0.0;
   minximhihi[0] = 0.0; minximlohi[0] = 0.0;
   minximhilo[0] = 0.0; minximlolo[0] = 0.0;
   pluxrehihi[1] = xrehihi; pluxrelohi[1] = xrelohi;
   pluxrehilo[1] = xrehilo; pluxrelolo[1] = xrelolo;
   pluximhihi[1] = ximhihi; pluximlohi[1] = ximlohi;
   pluximhilo[1] = ximhilo; pluximlolo[1] = ximlolo;
   minxrehihi[1] = -xrehihi; minxrelohi[1] = -xrelohi;
   minxrehilo[1] = -xrehilo; minxrelolo[1] = -xrelolo;
   minximhihi[1] = -ximhihi; minximlohi[1] = -ximlohi;
   minximhilo[1] = -ximhilo; minximlolo[1] = -ximlolo;

   for(int k=2; k<=deg; k++)
   {
      // pluxre[k] = (pluxre[k-1]*cr - pluxim[k-1]*sr)/k;
      qdf_mul(pluxrehihi[k-1],pluxrelohi[k-1],pluxrehilo[k-1],pluxrelolo[k-1],
              xrehihi,xrelohi,xrehilo,xrelolo,
              &pluxrehihi[k],&pluxrelohi[k],&pluxrehilo[k],&pluxrelolo[k]);
      qdf_mul(pluximhihi[k-1],pluximlohi[k-1],pluximhilo[k-1],pluximlolo[k-1],
              ximhihi,ximlohi,ximhilo,ximlolo,
              &tmphihi,&tmplohi,&tmphilo,&tmplolo);
      qdf_minus(&tmphihi,&tmplohi,&tmphilo,&tmplolo);
      qdf_inc(&pluxrehihi[k],&pluxrelohi[k],&pluxrehilo[k],&pluxrelolo[k],
              tmphihi,tmplohi,tmphilo,tmplolo);
      tmphihi = (double) k; tmplohi = 0.0; tmphilo = 0.0; tmplolo = 0.0;
      qdf_div(pluxrehihi[k],pluxrelohi[k],pluxrehilo[k],pluxrelolo[k],
              tmphihi,tmplohi,tmphilo,tmplolo,
              &pluxrehihi[k],&pluxrelohi[k],&pluxrehilo[k],&pluxrelolo[k]);
      // pluxim[k] = (pluxre[k-1]*sr + pluxim[k-1]*cr)/k;
      qdf_mul(pluxrehihi[k-1],pluxrelohi[k-1],pluxrehilo[k-1],pluxrelolo[k-1],
              ximhihi,ximlohi,ximhilo,ximlolo,
              &pluximhihi[k],&pluximlohi[k],&pluximhilo[k],&pluximlolo[k]);
      qdf_mul(pluximhihi[k-1],pluximlohi[k-1],pluximhilo[k-1],pluximlolo[k-1],
              xrehihi,xrelohi,xrehilo,xrelolo,
              &tmphihi,&tmplohi,&tmphilo,&tmplolo);
      qdf_inc(&pluximhihi[k],&pluximlohi[k],&pluximhilo[k],&pluximlolo[k],
              tmphihi,tmplohi,tmphilo,tmplolo);
      tmphihi = (double) k; tmplohi = 0.0; tmphilo = 0.0; tmplolo = 0.0;
      qdf_div(pluximhihi[k],pluximlohi[k],pluximhilo[k],pluximlolo[k],
              tmphihi,tmplohi,tmphilo,tmplolo,
              &pluximhihi[k],&pluximlohi[k],&pluximhilo[k],&pluximlolo[k]);
      // minxre[k] = (minxre[k-1]*(-cr) - minxim[k-1]*(-sr))/k;
      qdf_mul(minxrehihi[k-1],minxrelohi[k-1],minxrehilo[k-1],minxrelolo[k-1],
              -xrehihi,-xrelohi,-xrehilo,-xrelolo,
              &minxrehihi[k],&minxrelohi[k],&minxrehilo[k],&minxrelolo[k]);
      qdf_mul(minximhihi[k-1],minximlohi[k-1],minximhilo[k-1],minximlolo[k-1],
              -ximhihi,-ximlohi,-ximhilo,-ximlolo,
              &tmphihi,&tmplohi,&tmphilo,&tmplolo);
      qdf_minus(&tmphihi,&tmplohi,&tmphilo,&tmplolo);
      qdf_inc(&minxrehihi[k],&minxrelohi[k],&minxrehilo[k],&minxrelolo[k],
              tmphihi,tmplohi,tmphilo,tmplolo);
      tmphihi = (double) k; tmplohi = 0.0; tmphilo = 0.0; tmplolo = 0.0;
      qdf_div(minxrehihi[k],minxrelohi[k],minxrehilo[k],minxrelolo[k],
              tmphihi,tmplohi,tmphilo,tmplolo,
              &minxrehihi[k],&minxrelohi[k],&minxrehilo[k],&minxrelolo[k]);
      // minxim[k] = (minxre[k-1]*(-sr) + minxim[k-1]*(-cr))/k;
      qdf_mul(minxrehihi[k-1],minxrelohi[k-1],minxrehilo[k-1],minxrelolo[k-1],
              -ximhihi,-ximlohi,-ximhilo,-ximlolo,
              &minximhihi[k],&minximlohi[k],&minximhilo[k],&minximlolo[k]);
      qdf_mul(minximhihi[k-1],minximlohi[k-1],minximhilo[k-1],minximlolo[k-1],
              -xrehihi,-xrelohi,-xrehilo,-xrelolo,
              &tmphihi,&tmplohi,&tmphilo,&tmplolo);
      qdf_inc(&minximhihi[k],&minximlohi[k],&minximhilo[k],&minximlolo[k],
              tmphihi,tmplohi,tmphilo,tmplolo);
      tmphihi = (double) k; tmplohi = 0.0; tmphilo = 0.0; tmplolo = 0.0;
      qdf_div(minximhihi[k],minximlohi[k],minximhilo[k],minximlolo[k],
              tmphihi,tmplohi,tmphilo,tmplolo,
              &minximhihi[k],&minximlohi[k],&minximhilo[k],&minximlolo[k]);
   }
}

void random_cmplx4_exponential
 ( int deg,
   double *xrehihi, double *xrelohi, double *xrehilo, double *xrelolo,
   double *ximhihi, double *ximlohi, double *ximhilo, double *ximlolo,
   double *srehihi, double *srelohi, double *srehilo, double *srelolo,
   double *simhihi, double *simlohi, double *simhilo, double *simlolo )
{
   double tmphihi,tmplohi,tmphilo,tmplolo;

   random_quad_double(xrehihi,xrelohi,xrehilo,xrelolo);    // cos(a)

   qdf_sqr(*xrehihi,*xrelohi,*xrehilo,*xrelolo,
           &tmphihi,&tmplohi,&tmphilo,&tmplolo);           // cos^2(a)
   qdf_minus(&tmphihi,&tmplohi,&tmphilo,&tmplolo);         // -cos^2(a)
   qdf_inc_d(&tmphihi,&tmplohi,&tmphilo,&tmplolo,1.0);     // 1-cos^2(a)
   qdf_sqrt(tmphihi,tmplohi,tmphilo,tmplolo,
            ximhihi,ximlohi,ximhilo,ximlolo);              // sin is sqrt

   cmplx4_exponential
      (deg,*xrehihi,*xrelohi,*xrehilo,*xrelolo,
           *ximhihi,*ximlohi,*ximhilo,*ximlolo,
       srehihi,srelohi,srehilo,srelolo,simhihi,simlohi,simhilo,simlolo);
}

void random_cmplx4_exponentials
 ( int deg, 
   double *xrehihi, double *xrelohi, double *xrehilo, double *xrelolo,
   double *ximhihi, double *ximlohi, double *ximhilo, double *ximlolo,
   double *pluxrehihi, double *pluxrelohi, double *pluxrehilo,
   double *pluxrelolo, double *pluximhihi, double *pluximlohi,
   double *pluximhilo, double *pluximlolo,
   double *minxrehihi, double *minxrelohi, double *minxrehilo,
   double *minxrelolo, double *minximhihi, double *minximlohi,
   double *minximhilo, double *minximlolo )
{
   double tmphihi,tmplohi,tmphilo,tmplolo;

   random_quad_double(xrehihi,xrelohi,xrehilo,xrelolo);    // cos(a)

   qdf_sqr(*xrehihi,*xrelohi,*xrehilo,*xrelolo,
           &tmphihi,&tmplohi,&tmphilo,&tmplolo);           // cos^2(a)
   qdf_minus(&tmphihi,&tmplohi,&tmphilo,&tmplolo);         // -cos^2(a)
   qdf_inc_d(&tmphihi,&tmplohi,&tmphilo,&tmplolo,1.0);     // 1-cos^2(a)
   qdf_sqrt(tmphihi,tmplohi,tmphilo,tmplolo,
            ximhihi,ximlohi,ximhilo,ximlolo);              // sin is sqrt

   cmplx4_exponentials
      (deg,*xrehihi,*xrelohi,*xrehilo,*xrelolo,
           *ximhihi,*ximlohi,*ximhilo,*ximlolo,
           pluxrehihi,pluxrelohi,pluxrehilo,pluxrelolo,
           pluximhihi,pluximlohi,pluximhilo,pluximlolo,
           minxrehihi,minxrelohi,minxrehilo,minxrelolo,
           minximhihi,minximlohi,minximhilo,minximlolo);
}
