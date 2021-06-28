// The file random_series.cpp defines the functions with prototypes
// in random_series.h.

#include <cmath>
#include "random_numbers.h"

void dbl_exponential ( int deg, double x, double *s )
{
   s[0] = 1.0;

   for(int k=1; k<=deg; k++)
      s[k] = s[k-1]*x/k;
}

void dbl_exponentials ( int deg, double x, double *plux, double *minx )
{
   plux[0] = 1.0; minx[0] = 1.0;

   for(int k=1; k<=deg; k++)
   {
      plux[k] = plux[k-1]*x/k;
      minx[k] = minx[k-1]*(-x)/k;
   }
}

void random_dbl_exponential ( int deg, double *x, double *s )
{
   *x = random_double();
   dbl_exponential(deg,*x,s);
}

void random_dbl_exponentials
 ( int deg, double *x, double *plux, double *minx )
{
   *x = random_double();
   dbl_exponentials(deg,*x,plux,minx);
}

void cmplx_exponential
 ( int deg, double xre, double xim, double *sre, double *sim )
{
   sre[0] = 1.0; sim[0] = 0.0;
   sre[1] = xre; sim[1] = xim;

   for(int k=2; k<=deg; k++)
   {
      sre[k] = (sre[k-1]*xre - sim[k-1]*xim)/k;
      sim[k] = (sre[k-1]*xim + sim[k-1]*xre)/k;
   }
}

void cmplx_exponentials
 ( int deg, double xre, double xim,
   double *pluxre, double *pluxim, double *minxre, double *minxim )
{
   pluxre[0] = 1.0; pluxim[0] = 0.0; minxre[0] = 1.0; minxim[0] = 0.0;
   pluxre[1] = xre; pluxim[1] = xim;  minxre[1] = -xre; minxim[1] = -xim;

   for(int k=2; k<=deg; k++)
   {
      pluxre[k] = (pluxre[k-1]*xre - pluxim[k-1]*xim)/k;
      pluxim[k] = (pluxre[k-1]*xim + pluxim[k-1]*xre)/k;
      minxre[k] = (minxre[k-1]*(-xre) - minxim[k-1]*(-xim))/k;
      minxim[k] = (minxre[k-1]*(-xim) + minxim[k-1]*(-xre))/k;
   }
}

void random_cmplx_exponential
 ( int deg, double *xre, double *xim, double *sre, double *sim )
{
   const double r = random_angle();

   *xre = cos(r);
   *xim = sin(r);

   cmplx_exponential(deg,*xre,*xim,sre,sim);
}

void random_cmplx_exponentials
 ( int deg, double *xre, double *xim,
   double *pluxre, double *pluxim, double *minxre, double *minxim )
{
   const double r = random_angle();

   *xre = cos(r);
   *xim = sin(r);

   cmplx_exponentials(deg,*xre,*xim,pluxre,pluxim,minxre,minxim);
}

void dbl_logarithm ( int deg, double x, double *s )
{
   double pwrx = x;
   s[0] = x;

   for(int k=1; k<=deg; k++)
   {
      pwrx = -x*pwrx; // alternating sign
      s[k] = pwrx/k;
   }
}

void cmplx_logarithm 
 ( int deg, double xre, double xim, double *sre, double *sim )
{
   double pwrxre = xre;
   double pwrxim = xim;
   sre[0] = xre;
   sim[0] = xim;

   for(int k=1; k<=deg; k++)
   {
      pwrxre = xre*pwrxre - xim*pwrxim;
      pwrxim = xre*pwrxim + xim*pwrxre;
      pwrxre = -pwrxre;                 // alternating sign
      pwrxim = -pwrxim;
      sre[k] = pwrxre/k;
      sim[k] = pwrxim/k;
   }
}

void random_dbl_logarithm ( int deg, double *x, double *s )
{
   *x = random_double();
   dbl_logarithm(deg,*x,s);
}

void random_cmplx_logarithm
 ( int deg, double *xre, double *xim, double *sre, double *sim )
{
   const double r = random_angle();

   *xre = cos(r);
   *xim = sin(r);

   cmplx_logarithm(deg,*xre,*xim,sre,sim);
}
