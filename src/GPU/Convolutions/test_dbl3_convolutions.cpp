/* Tests the product of two series in triple double precision. */

#include <iostream>
#include "dbl3_convolutions_host.h"

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
   double* xhi = new double[deg+1];
   double* xmi = new double[deg+1];
   double* xlo = new double[deg+1];
   double* yhi = new double[deg+1];
   double* ymi = new double[deg+1];
   double* ylo = new double[deg+1];
   double* zhi = new double[deg+1];
   double* zmi = new double[deg+1];
   double* zlo = new double[deg+1];

   for(int k=0; k<=deg; k++)
   {
      xhi[k] = 1.0;
      xmi[k] = 0.0;
      xlo[k] = 0.0;
      yhi[k] = 0.0;
      ymi[k] = 0.0;
      ylo[k] = 0.0;
   }
   yhi[0] = 1.0; yhi[1] = -1.0;

   CPU_dbl3_product(deg,xhi,xmi,xlo,yhi,ymi,ylo,zhi,zmi,zlo);

   cout << "Series of 1/(1-x) multiplied with 1-x : " << endl;

   for(int k=0; k<=deg; k++)
   {
      cout << "zhi[" << k << "] : " << zhi[k] << endl;
      cout << "zmi[" << k << "] : " << zmi[k] << endl;
      cout << "zlo[" << k << "] : " << zlo[k] << endl;
   }
}

void test_complex ( int deg )
{
   double* xrehi = new double[deg+1];
   double* xremi = new double[deg+1];
   double* xrelo = new double[deg+1];
   double* ximhi = new double[deg+1];
   double* ximmi = new double[deg+1];
   double* ximlo = new double[deg+1];
   double* yrehi = new double[deg+1];
   double* yremi = new double[deg+1];
   double* yrelo = new double[deg+1];
   double* yimhi = new double[deg+1];
   double* yimmi = new double[deg+1];
   double* yimlo = new double[deg+1];
   double* zrehi = new double[deg+1];
   double* zremi = new double[deg+1];
   double* zrelo = new double[deg+1];
   double* zimhi = new double[deg+1];
   double* zimmi = new double[deg+1];
   double* zimlo = new double[deg+1];

   for(int k=0; k<=deg; k++)
   {
      xrehi[k] = 1.0; xremi[k] = 0.0; xrelo[k] = 0.0;
      ximhi[k] = 0.0; ximmi[k] = 0.0; ximlo[k] = 0.0;
      yrehi[k] = 0.0; yremi[k] = 0.0; yrelo[k] = 0.0;
      yimhi[k] = 0.0; yimmi[k] = 0.0; yimlo[k] = 0.0;
   }
   yrehi[0] = 1.0; yrehi[1] = -1.0;

   CPU_cmplx3_product(deg,xrehi,xremi,xrelo,ximhi,ximmi,ximlo,
                          yrehi,yremi,yrelo,yimhi,yimmi,yimlo,
                          zrehi,zremi,zrelo,zimhi,zimmi,zimlo);

   cout << "Series of 1/(1-x) multiplied with 1-x : " << endl;

   for(int k=0; k<=deg; k++)
   {
      cout << "zrehi[" << k << "] : " << zrehi[k] << endl;
      cout << "zremi[" << k << "] : " << zremi[k] << endl;
      cout << "zrelo[" << k << "] : " << zrelo[k] << endl;
      cout << "zimhi[" << k << "] : " << zimhi[k] << endl;
      cout << "zimmi[" << k << "] : " << zimmi[k] << endl;
      cout << "zimlo[" << k << "] : " << zimlo[k] << endl;
   }
}
