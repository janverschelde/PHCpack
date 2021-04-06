/* The file dbl_convolutions_host.cpp defines the functions 
 * specified in dbl_convolutions_host.h. */

#include "dbl_convolutions_host.h"

void CPU_dbl_product ( int deg, double *x, double *y, double *z )
{
   z[0] = x[0]*y[0];
   for(int k=1; k<=deg; k++)
   {
      z[k] = x[0]*y[k];
      for(int i=1; i<=k; i++)
         z[k] = z[k] + x[i]*y[k-i];
   }
}

void CPU_dbl_Laurent_product
 ( int deg, int xe, int ye, int *ze, double *x, double *y, double *z )
{
   *ze = xe + ye;
   CPU_dbl_product(deg,x,y,z);
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

   for(int k=1; k<=deg; k++)
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

void CPU_cmplx_Laurent_product
 ( int deg, int xe, int ye, int *ze, double *xre, double *xim,
   double *yre, double *yim, double *zre, double *zim )
{
   *ze = xe + ye;
   CPU_cmplx_product(deg,xre,xim,yre,yim,zre,zim);
}

void CPU_dbl_inverse ( int deg, double *x, double *y )
{
   y[0] = 1.0/x[0];
   for(int i=1; i<=deg; i++)
   {
      y[i] = -x[1]*y[i-1];
      for(int j=2; j<=i; j++) y[i] = y[i] - x[j]*y[i-j];
      y[i] = y[i]/x[0];
   }
}

void CPU_dbl_Laurent_inverse
 ( int deg, int xe, int *ye, double *x, double *y )
{
   *ye = -xe;
   CPU_dbl_inverse(deg,x,y);
}

void CPU_cmplx_inverse
 ( int deg, double *xre, double *xim, double *yre, double *yim )
{
   double rpa,ipa;
   double xr0,xi0;
   double yr0,yi0;
   int idx;

   xr0 = xre[0]; xi0 = xim[0];
   rpa = xr0*xr0 + xi0*xi0;
   yre[0] = xr0/rpa; yim[0] = -xi0/rpa;

   for(int i=1; i<=deg; i++)
   {
      xr0 = xre[1];   xi0 = xim[1]; 
      yr0 = yre[i-1]; yi0 = yim[i-1]; 
      rpa = - xr0*yr0 + xi0*yi0;
      ipa = - xi0*yr0 - xr0*yi0;

      for(int j=2; j<=i; j++)
      {
         idx = i-j;
         xr0 = xre[j];   xi0 = xim[j];
         yr0 = yre[idx]; yi0 = yim[idx];
         rpa = rpa - (xr0*yr0 - xi0*yi0);
         ipa = ipa - (xi0*yr0 + xr0*yi0);
      }
      yre[i] = rpa; yim[i] = ipa;

      xr0 = yre[0]; xi0 = yim[0]; // contains 1/x[0]
      yr0 = yre[i]; yi0 = yim[i];
      rpa = xr0*yr0 - xi0*yi0;
      ipa = xi0*yr0 + xr0*yi0;
      yre[i] = rpa; yim[i] = ipa;
   }
}

void CPU_cmplx_Laurent_inverse
 ( int deg, int xe, int *ye, double *xre, double *xim,
   double *yre, double *yim )
{
   *ye = -xe;
   CPU_cmplx_inverse(deg,xre,xim,yre,yim);
}
