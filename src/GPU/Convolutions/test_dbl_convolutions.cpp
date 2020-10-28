/* Tests the product of two series in double precision. */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector_types.h>
#include "random_numbers.h"
#include "dbl_convolutions_host.h"
#include "dbl_convolutions_kernels.h"

using namespace std;

void test_real ( int deg );
/*
 * DESCRIPTION :
 *   Multiplies the power series of 1/(1-x) with 1+x,
 *   truncated to degree deg, for real coefficients. */

void test_complex ( int deg );
/*
 * DESCRIPTION :
 *   Multiplies the power series of 1/(1-x) with 1+x,
 *   truncated to degree deg, for complex coefficients. */

void test_real_exponential ( int deg );
/*
 * DESCRIPTION :
 *   Multiplies the power series for exp(x) with exp(-x)
 *   for some random x in [-1,+1], for real coefficients
 *   of a series of degree truncated to deg. */

void test_complex_exponential ( int deg );
/*
 * DESCRIPTION :
 *   Multiplies the power series for exp(x) with exp(-x)
 *   for some random complex number on the unit circle,
 *   for series of degree truncated to deg. */

int main ( void )
{
   int deg;

   cout << "Give a degree larger than one : "; cin >> deg;

   if(deg > 0) test_real(deg);
   if(deg > 0) test_complex(deg);
   if(deg > 0) test_real_exponential(deg);
   if(deg > 0) test_complex_exponential(deg);

   return 0;
}

void test_real ( int deg )
{
   double *x = new double[deg+1];
   double *y = new double[deg+1];
   double *z_h = new double[deg+1];
   double *z_d = new double[deg+1];

   for(int k=0; k<=deg; k++)
   {
      x[k] = 1.0; x[k] = 1.0;
      y[k] = 0.0; y[k] = 0.0;
   }
   y[0] = 1.0; y[1] = -1.0;

   CPU_dbl_product(deg,x,y,z_h);

   cout << "Series of 1/(1-x) multiplied with 1-x :" << endl;

   for(int k=0; k<=deg; k++)
   {
      cout << "  z[" << k << "] : " << z_h[k];
      if((k+1) % 4 == 0) cout << endl;
   }
   cout << endl;

   GPU_dbl_product(x,y,z_d,deg,1,deg+1);

   cout << "GPU computed product :" << endl;

   for(int k=0; k<=deg; k++)
   {
      cout << "  z[" << k << "] : " << z_d[k];
      if((k+1) % 4 == 0) cout << endl;
   }
   cout << endl;
}

void test_complex ( int deg )
{
   double *xre = new double[deg+1];
   double *xim = new double[deg+1];
   double *yre = new double[deg+1];
   double *yim = new double[deg+1];
   double *zre_h = new double[deg+1];
   double *zim_h = new double[deg+1];
   double *zre_d = new double[deg+1];
   double *zim_d = new double[deg+1];

   for(int k=0; k<=deg; k++)
   {
      xre[k] = 1.0; xim[k] = 0.0;
      yre[k] = 0.0; yim[k] = 0.0;
   }
   yre[0] = 1.0; yre[1] = -1.0;

   CPU_cmplx_product(deg,xre,xim,yre,yim,zre_h,zim_h);

   cout << "Series of 1/(1-x) multiplied with 1-x : " << endl;

   for(int k=0; k<=deg; k++)
   {
      cout << "  zre[" << k << "] : " << zre_h[k];
      cout << "  zim[" << k << "] : " << zim_h[k];
      if((k+1) % 2 == 0) cout << endl;
   }
   cout << endl;

   GPU_cmplx_product(xre,xim,yre,yim,zre_d,zim_d,deg,1,deg+1);

   cout << "GPU computed product :" << endl;

   for(int k=0; k<=deg; k++)
   {
      cout << "  zre[" << k << "] : " << zre_d[k];
      cout << "  zim[" << k << "] : " << zim_d[k];
      if((k+1) % 2 == 0) cout << endl;
   }
   cout << endl;
}

void test_real_exponential ( int deg )
{
   double *x = new double[deg+1];
   double *y = new double[deg+1];
   double *z_h = new double[deg+1];
   double *z_d = new double[deg+1];
   double r = random_double();

   x[0] = 1.0; y[0] = 1.0;

   for(int k=1; k<=deg; k++)
   {
      x[k] = x[k-1]*r/k;
      y[k] = y[k-1]*(-r)/k;
   }

   CPU_dbl_product(deg,x,y,z_h);

   cout << scientific << setprecision(16);
   cout << "Series of exp(x)*exp(-x), for x = " << r << endl;

   double sum = 0.0;
   for(int k=0; k<=deg; k++) sum = sum + z_h[k];
   cout << "Summation of all coefficients of the product ..." << endl;
   cout << "  sum : " << sum << endl;

   GPU_dbl_product(x,y,z_d,deg,1,deg+1);

   sum = 0.0;
   for(int k=0; k<=deg; k++) sum = sum + z_d[k];
   cout << "Summation of all coefficients of the GPU computed product ..."
        << endl;
   cout << "  sum : " << sum << endl;
}

void test_complex_exponential ( int deg )
{
   double *xre = new double[deg+1];
   double *xim = new double[deg+1];
   double *yre = new double[deg+1];
   double *yim = new double[deg+1];
   double *zre_h = new double[deg+1];
   double *zim_h = new double[deg+1];
   double *zre_d = new double[deg+1];
   double *zim_d = new double[deg+1];
   double r = random_angle();
   double cr = cos(r);
   double sr = sin(r);

   xre[0] = 1.0; xim[0] = 0.0; yre[0] = 1.0; yim[0] = 0.0;
   xre[1] = cr;  xim[1] = sr;  yre[1] = -cr; yim[1] = -sr;

   for(int k=2; k<=deg; k++)
   {
      xre[k] = (xre[k-1]*cr - xim[k-1]*sr)/k;
      xim[k] = (xre[k-1]*sr + xim[k-1]*cr)/k;
      yre[k] = (yre[k-1]*(-cr) - yim[k-1]*(-sr))/k;
      yim[k] = (yre[k-1]*(-sr) + yim[k-1]*(-cr))/k;
   }

   CPU_cmplx_product(deg,xre,xim,yre,yim,zre_h,zim_h);

   cout << scientific << setprecision(16);

   cout << "Series of exp(x)*exp(-x), for xre = " << cr << endl;
   cout << "  and xim = " << sr << endl;

   double sumre = 0.0;
   double sumim = 0.0;

   for(int k=0; k<=deg; k++) 
   {
      sumre = sumre + zre_h[k];
      sumim = sumim + zim_h[k];
   }
   cout << "Summation of all coefficients of the product ..." << endl;
   cout << "  sumre : " << sumre << endl;
   cout << "  sumim : " << sumim << endl;

   GPU_cmplx_product(xre,xim,yre,yim,zre_d,zim_d,deg,1,deg+1);

   sumre = 0.0; sumim = 0.0;
   for(int k=0; k<=deg; k++) 
   {
      sumre = sumre + zre_d[k];
      sumim = sumim + zim_d[k];
   }
   cout << "Summation of all coefficients of the GPU computed product ..."
        << endl;
   cout << "  sumre : " << sumre << endl;
   cout << "  sumim : " << sumim << endl;

}
