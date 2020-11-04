// The file random3_series.cpp defines functions specified 
// in random3_series.h.

#include "random3_vectors.h"
#include "triple_double_functions.h"

void dbl3_exponentials
 ( int deg, double xhi, double xmi, double xlo, 
   double *pluxhi, double *pluxmi, double *pluxlo,
   double *minxhi, double *minxmi, double *minxlo )
{
   double fhi,fmi,flo;

   pluxhi[0] = 1.0; pluxmi[0] = 0.0; pluxlo[0] = 0.0;
   minxhi[0] = 1.0; minxmi[0] = 0.0; minxlo[0] = 0.0;
   pluxhi[1] = xhi; pluxmi[1] = xmi; pluxlo[1] = xlo;
   minxhi[1] = -xhi; minxmi[1] = -xmi; minxlo[1] = -xlo;

   for(int k=2; k<=deg; k++)
   {
      tdf_mul(pluxhi[k-1],pluxmi[k-1],pluxlo[k-1],xhi,xmi,xlo,
              &pluxhi[k],&pluxmi[k],&pluxlo[k]);  // x[k] = x[k-1]*r;
      tdf_mul(minxhi[k-1],minxmi[k-1],minxlo[k-1],-xhi,-xmi,-xlo,
              &minxhi[k],&minxmi[k],&minxlo[k]);  // y[k] = y[k-1]*(-r);
      fhi = (double) k; fmi = 0.0; flo = 0.0;
      tdf_div(pluxhi[k],pluxmi[k],pluxlo[k],fhi,fmi,flo,
              &pluxhi[k],&pluxmi[k],&pluxlo[k]);
      tdf_div(minxhi[k],minxmi[k],minxlo[k],fhi,fmi,flo,
              &minxhi[k],&minxmi[k],&minxlo[k]);
   }
}

void random_dbl3_exponentials
 ( int deg, double *xhi, double *xmi, double *xlo,
   double *pluxhi, double *pluxmi, double *pluxlo,
   double *minxhi, double *minxmi, double *minxlo )
{
   random_triple_double(xhi,xmi,xlo);

   dbl3_exponentials(deg,*xhi,*xmi,*xlo,
                     pluxhi,pluxmi,pluxlo,minxhi,minxmi,minxlo);
}

void cmplx3_exponentials
 ( int deg, double xrehi, double xremi, double xrelo,
            double ximhi, double ximmi, double ximlo,
   double *pluxrehi, double *pluxremi, double *pluxrelo,
   double *pluximhi, double *pluximmi, double *pluximlo,
   double *minxrehi, double *minxremi, double *minxrelo,
   double *minximhi, double *minximmi, double *minximlo )
{
   double tmphi,tmpmi,tmplo;

   pluxrehi[0] = 1.0; pluxremi[0] = 0.0; pluxrelo[0] = 0.0;
   minxrehi[0] = 1.0; minxremi[0] = 0.0; minxrelo[0] = 0.0;
   pluximhi[0] = 0.0; pluximmi[0] = 0.0; pluximlo[0] = 0.0;
   minximhi[0] = 0.0; minximmi[0] = 0.0; minximlo[0] = 0.0;
   pluxrehi[1] = xrehi; pluxremi[1] = xremi; pluxrelo[1] = xrelo;
   pluximhi[1] = ximhi; pluximmi[1] = ximmi; pluximlo[1] = ximlo;
   minxrehi[1] = -xrehi; minxremi[1] = -xremi; minxrelo[1] = -xrelo;
   minximhi[1] = -ximhi; minximmi[1] = -ximmi; minximlo[1] = -ximlo;

   for(int k=2; k<=deg; k++)
   {
      // pluxre[k] = (pluxre[k-1]*cr - xim[k-1]*sr)/k;
      tdf_mul(pluxrehi[k-1],pluxremi[k-1],pluxrelo[k-1],xrehi,xremi,xrelo,
              &pluxrehi[k],&pluxremi[k],&pluxrelo[k]);
      tdf_mul(pluximhi[k-1],pluximmi[k-1],pluximlo[k-1],ximhi,ximmi,ximlo,
              &tmphi,&tmpmi,&tmplo);
      tdf_minus(&tmphi,&tmpmi,&tmplo);
      tdf_inc(&pluxrehi[k],&pluxremi[k],&pluxrelo[k],tmphi,tmpmi,tmplo);
      tmphi = (double) k; tmpmi = 0.0; tmplo = 0.0;
      tdf_div(pluxrehi[k],pluxremi[k],pluxrelo[k],tmphi,tmpmi,tmplo,
              &pluxrehi[k],&pluxremi[k],&pluxrelo[k]);
      // xim[k] = (pluxre[k-1]*sr + xim[k-1]*cr)/k;
      tdf_mul(pluxrehi[k-1],pluxremi[k-1],pluxrelo[k-1],ximhi,ximmi,ximlo,
              &pluximhi[k],&pluximmi[k],&pluximlo[k]);
      tdf_mul(pluximhi[k-1],pluximmi[k-1],pluximlo[k-1],xrehi,xremi,xrelo,
              &tmphi,&tmpmi,&tmplo);
      tdf_inc(&pluximhi[k],&pluximmi[k],&pluximlo[k],tmphi,tmpmi,tmplo);
      tmphi = (double) k; tmpmi = 0.0; tmplo = 0.0;
      tdf_div(pluximhi[k],pluximmi[k],pluximlo[k],tmphi,tmpmi,tmplo,
              &pluximhi[k],&pluximmi[k],&pluximlo[k]);
      // minxre[k] = (minxre[k-1]*(-cr) - minxim[k-1]*(-sr))/k;
      tdf_mul(minxrehi[k-1],minxremi[k-1],minxrelo[k-1],-xrehi,-xremi,-xrelo,
              &minxrehi[k],&minxremi[k],&minxrelo[k]);
      tdf_mul(minximhi[k-1],minximmi[k-1],minximlo[k-1],-ximhi,-ximmi,-ximlo,
              &tmphi,&tmpmi,&tmplo);
      tdf_minus(&tmphi,&tmpmi,&tmplo);
      tdf_inc(&minxrehi[k],&minxremi[k],&minxrelo[k],tmphi,tmpmi,tmplo);
      tmphi = (double) k; tmpmi = 0.0; tmplo = 0.0;
      tdf_div(minxrehi[k],minxremi[k],minxrelo[k],tmphi,tmpmi,tmplo,
              &minxrehi[k],&minxremi[k],&minxrelo[k]);
      // minxim[k] = (minxre[k-1]*(-sr) + minxim[k-1]*(-cr))/k;
      tdf_mul(minxrehi[k-1],minxremi[k-1],minxrelo[k-1],-ximhi,-ximmi,-ximlo,
              &minximhi[k],&minximmi[k],&minximlo[k]);
      tdf_mul(minximhi[k-1],minximmi[k-1],minximlo[k-1],-xrehi,-xremi,-xrelo,
              &tmphi,&tmpmi,&tmplo);
      tdf_inc(&minximhi[k],&minximmi[k],&minximlo[k],tmphi,tmpmi,tmplo);
      tmphi = (double) k; tmpmi = 0.0; tmplo = 0.0;
      tdf_div(minximhi[k],minximmi[k],minximlo[k],tmphi,tmpmi,tmplo,
              &minximhi[k],&minximmi[k],&minximlo[k]);
   }
}

void random_cmplx3_exponentials
 ( int deg, double *xrehi, double *xremi, double *xrelo,
            double *ximhi, double *ximmi, double *ximlo,
   double *pluxrehi, double *pluxremi, double *pluxrelo,
   double *pluximhi, double *pluximmi, double *pluximlo,
   double *minxrehi, double *minxremi, double *minxrelo,
   double *minximhi, double *minximmi, double *minximlo )
{
   double tmphi,tmpmi,tmplo;

   random_triple_double(xrehi,xremi,xrelo);             // cos(a)

   tdf_sqr(*xrehi,*xremi,*xrelo,&tmphi,&tmpmi,&tmplo);  // cos^2(a)
   tdf_minus(&tmphi,&tmpmi,&tmplo);                     // -cos^2(a)
   tdf_inc_d(&tmphi,&tmpmi,&tmplo,1.0);                 // 1-cos^2(a)
   tdf_sqrt(tmphi,tmpmi,tmplo,ximhi,ximmi,ximlo);       // sin is sqrt

   cmplx3_exponentials
      (deg,*xrehi,*xremi,*xrelo,*ximhi,*ximmi,*ximlo,
       pluxrehi,pluxremi,pluxrelo,pluximhi,pluximmi,pluximlo,
       minxrehi,minxremi,minxrelo,minximhi,minximmi,minximlo);
}
