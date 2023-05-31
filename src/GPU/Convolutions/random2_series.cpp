// The file random2_series.cpp defines functions specified in
// in random2_series.h.

#include "random2_vectors.h"
#include "double_double_functions.h"

void dbl2_exponential
 ( int deg, double xhi, double xlo, double *shi, double *slo )
{
   double fhi,flo;

   shi[0] = 1.0; slo[0] = 0.0;
   shi[1] = xhi; slo[1] = xlo;

   for(int k=2; k<=deg; k++)
   {
      // x[k] = x[k-1]*r
      ddf_mul(shi[k-1],slo[k-1],xhi,xlo,&shi[k],&slo[k]);
      // x[k] = x[k]/k
      fhi = (double) k; flo = 0.0;
      ddf_div(shi[k],slo[k],fhi,flo,&shi[k],&slo[k]);
   }
}

void dbl2_exponentials
 ( int deg, double xhi, double xlo, 
   double *pluxhi, double *pluxlo, double *minxhi, double *minxlo )
{
   double fhi,flo;

   pluxhi[0] = 1.0; pluxlo[0] = 0.0; minxhi[0] = 1.0; minxlo[0] = 0.0;
   pluxhi[1] = xhi; pluxlo[1] = xlo; minxhi[1] = -xhi; minxlo[1] = -xlo;

   for(int k=2; k<=deg; k++)
   {
      // x[k] = x[k-1]*r;
      ddf_mul(pluxhi[k-1],pluxlo[k-1],xhi,xlo,&pluxhi[k],&pluxlo[k]);
      ddf_mul(minxhi[k-1],minxlo[k-1],-xhi,-xlo,&minxhi[k],&minxlo[k]); 
      // y[k] = y[k-1]*(-r);
      fhi = (double) k; flo = 0.0;
      ddf_div(pluxhi[k],pluxlo[k],fhi,flo,&pluxhi[k],&pluxlo[k]);
      ddf_div(minxhi[k],minxlo[k],fhi,flo,&minxhi[k],&minxlo[k]);
   }
}

void random_dbl2_exponential
 ( int deg, double *xhi, double *xlo, double *shi, double *slo )
{
   random_double_double(xhi,xlo);
   dbl2_exponential(deg,*xhi,*xlo,shi,slo);
}

void random_dbl2_exponentials
 ( int deg, double *xhi, double *xlo,
   double *pluxhi, double *pluxlo, double *minxhi, double *minxlo )
{
   random_double_double(xhi,xlo);
   dbl2_exponentials(deg,*xhi,*xlo,pluxhi,pluxlo,minxhi,minxlo);
}

void cmplx2_exponential
 ( int deg, double xrehi, double xrelo, double ximhi, double ximlo,
   double *srehi, double *srelo, double *simhi, double *simlo )
{
   double acchi,acclo;

   srehi[0] = 1.0; srelo[0] = 0.0;
   simhi[0] = 0.0; simlo[0] = 0.0;
   srehi[1] = xrehi; srelo[1] = xrelo;
   simhi[1] = ximhi; simlo[1] = ximlo;

   for(int k=2; k<=deg; k++)
   {
      // sre[k] = (sre[k-1]*xre - sim[k-1]*xim)/k;
      ddf_mul(srehi[k-1],srelo[k-1],xrehi,xrelo,&srehi[k],&srelo[k]);
      ddf_mul(simhi[k-1],simlo[k-1],ximhi,ximlo,&acchi,&acclo);
      ddf_minus(&acchi,&acclo);
      ddf_inc(&srehi[k],&srelo[k],acchi,acclo);
      acchi = (double) k; acclo = 0.0;
      ddf_div(srehi[k],srelo[k],acchi,acclo,&srehi[k],&srelo[k]);
      // sim[k] = (sre[k-1]*xim + sim[k-1]*xre)/k;
      ddf_mul(srehi[k-1],srelo[k-1],ximhi,ximlo,&simhi[k],&simlo[k]);
      ddf_mul(simhi[k-1],simlo[k-1],xrehi,xrelo,&acchi,&acclo);
      ddf_inc(&simhi[k],&simlo[k],acchi,acclo);
      acchi = (double) k; acclo = 0.0;
      ddf_div(simhi[k],simlo[k],acchi,acclo,&simhi[k],&simlo[k]);
   }
}

void cmplx2_exponentials
 ( int deg, double xrehi, double xrelo, double ximhi, double ximlo,
   double *pluxrehi, double *pluxrelo, double *pluximhi, double *pluximlo,
   double *minxrehi, double *minxrelo, double *minximhi, double *minximlo )
{
   double tmphi,tmplo;

   pluxrehi[0] = 1.0; pluxrelo[0] = 0.0;
   minxrehi[0] = 1.0; minxrelo[0] = 0.0;
   pluximhi[0] = 0.0; pluximlo[0] = 0.0;
   minximhi[0] = 0.0; minximlo[0] = 0.0;
   pluxrehi[1] = xrehi; pluxrelo[1] = xrelo;
   pluximhi[1] = ximhi; pluximlo[1] = ximlo;
   minxrehi[1] = -(xrehi); minxrelo[1] = -(xrelo);
   minximhi[1] = -(ximhi); minximlo[1] = -(ximlo);

   for(int k=2; k<=deg; k++)
   {
      // pluxre[k] = (pluxre[k-1]*cr - pluxim[k-1]*sr)/k;
      ddf_mul(pluxrehi[k-1],pluxrelo[k-1],xrehi,xrelo,
              &pluxrehi[k],&pluxrelo[k]);
      ddf_mul(pluximhi[k-1],pluximlo[k-1],ximhi,ximlo,&tmphi,&tmplo);
      ddf_minus(&tmphi,&tmplo);
      ddf_inc(&pluxrehi[k],&pluxrelo[k],tmphi,tmplo);
      tmphi = (double) k; tmplo = 0.0;
      ddf_div(pluxrehi[k],pluxrelo[k],tmphi,tmplo,&pluxrehi[k],&pluxrelo[k]);
      // pluxim[k] = (pluxre[k-1]*sr + pluxim[k-1]*cr)/k;
      ddf_mul(pluxrehi[k-1],pluxrelo[k-1],ximhi,ximlo,
              &pluximhi[k],&pluximlo[k]);
      ddf_mul(pluximhi[k-1],pluximlo[k-1],xrehi,xrelo,&tmphi,&tmplo);
      ddf_inc(&pluximhi[k],&pluximlo[k],tmphi,tmplo);
      tmphi = (double) k; tmplo = 0.0;
      ddf_div(pluximhi[k],pluximlo[k],tmphi,tmplo,&pluximhi[k],&pluximlo[k]);
      // minxre[k] = (minxre[k-1]*(-cr) - minxim[k-1]*(-sr))/k;
      ddf_mul(minxrehi[k-1],minxrelo[k-1],-xrehi,-xrelo,
              &minxrehi[k],&minxrelo[k]);
      ddf_mul(minximhi[k-1],minximlo[k-1],-ximhi,-ximlo,&tmphi,&tmplo);
      ddf_minus(&tmphi,&tmplo);
      ddf_inc(&minxrehi[k],&minxrelo[k],tmphi,tmplo);
      tmphi = (double) k; tmplo = 0.0;
      ddf_div(minxrehi[k],minxrelo[k],tmphi,tmplo,&minxrehi[k],&minxrelo[k]);
      // minxim[k] = (minxre[k-1]*(-sr) + minxim[k-1]*(-cr))/k;
      ddf_mul(minxrehi[k-1],minxrelo[k-1],-ximhi,-ximlo,
              &minximhi[k],&minximlo[k]);
      ddf_mul(minximhi[k-1],minximlo[k-1],-xrehi,-xrelo,&tmphi,&tmplo);
      ddf_inc(&minximhi[k],&minximlo[k],tmphi,tmplo);
      tmphi = (double) k; tmplo = 0.0;
      ddf_div(minximhi[k],minximlo[k],tmphi,tmplo,&minximhi[k],&minximlo[k]);
   }
}

void random_cmplx2_exponential
 ( int deg, double *xrehi, double *xrelo, double *ximhi, double *ximlo,
   double *srehi, double *srelo, double *simhi, double *simlo )
{
   double tmphi,tmplo;

   random_double_double(xrehi,xrelo);           // cos(a)

   ddf_sqr(*xrehi,*xrelo,&tmphi,&tmplo);        // cos^2(a)
   ddf_minus(&tmphi,&tmplo);                    // -cos^2(a)
   ddf_inc_d(&tmphi,&tmplo,1.0);                // 1-cos^2(a)
   ddf_sqrt(tmphi,tmplo,ximhi,ximlo);           // sin is sqrt

   cmplx2_exponential
      (deg,*xrehi,*xrelo,*ximhi,*ximlo,srehi,srelo,simhi,simlo);
}

void random_cmplx2_exponentials
 ( int deg, double *xrehi, double *xrelo, double *ximhi, double *ximlo,
   double *pluxrehi, double *pluxrelo, double *pluximhi, double *pluximlo,
   double *minxrehi, double *minxrelo, double *minximhi, double *minximlo )
{
   double tmphi,tmplo;

   random_double_double(xrehi,xrelo);           // cos(a)

   ddf_sqr(*xrehi,*xrelo,&tmphi,&tmplo);        // cos^2(a)
   ddf_minus(&tmphi,&tmplo);                    // -cos^2(a)
   ddf_inc_d(&tmphi,&tmplo,1.0);                // 1-cos^2(a)
   ddf_sqrt(tmphi,tmplo,ximhi,ximlo);           // sin is sqrt

   cmplx2_exponentials(deg,*xrehi,*xrelo,*ximhi,*ximlo,
                       pluxrehi,pluxrelo,pluximhi,pluximlo,
                       minxrehi,minxrelo,minximhi,minximlo);
}
