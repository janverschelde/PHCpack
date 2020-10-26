/* Tests the product of two series in double precision. */

#include <iostream>
#include "dbl_convolutions_host.h"

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

int main ( void )
{
   int deg;

   cout << "Give a degree larger than one : "; cin >> deg;

   if(deg > 0) test_real(deg);
   if(deg > 0) test_complex(deg);

   return 0;
}

void test_real ( int deg )
{
   double *x = new double[deg+1];
   double *y = new double[deg+1];
   double *z = new double[deg+1];

   for(int k=0; k<=deg; k++)
   {
      x[k] = 1.0;
      y[k] = 0.0;
   }
   y[0] = 1.0; y[1] = -1.0;

   CPU_dbl_product(deg,x,y,z);

   cout << "Series of 1/(1-x) multiplied with 1-x : " << endl;

   for(int k=0; k<=deg; k++)
      cout << "z[" << k << "] : " << z[k] << endl;
}

void test_complex ( int deg )
{
   double *xre = new double[deg+1];
   double *xim = new double[deg+1];
   double *yre = new double[deg+1];
   double *yim = new double[deg+1];
   double *zre = new double[deg+1];
   double *zim = new double[deg+1];

   for(int k=0; k<=deg; k++)
   {
      xre[k] = 1.0; xim[k] = 0.0;
      yre[k] = 0.0; yim[k] = 0.0;
   }
   yre[0] = 1.0; yre[1] = -1.0;

   CPU_cmplx_product(deg,xre,xim,yre,yim,zre,zim);

   cout << "Series of 1/(1-x) multiplied with 1-x : " << endl;

   for(int k=0; k<=deg; k++)
   {
      cout << "zre[" << k << "] : " << zre[k] << endl;
      cout << "zim[" << k << "] : " << zim[k] << endl;
   }
}
