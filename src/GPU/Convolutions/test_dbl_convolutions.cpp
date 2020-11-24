/* Tests the product of two series in double precision. */

#include <ctime>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector_types.h>
#include "random_series.h"
#include "dbl_convolutions_host.h"
#include "dbl_convolutions_kernels.h"

using namespace std;

double test_real ( int deg );
/*
 * DESCRIPTION :
 *   Multiplies the power series of 1/(1-x) with 1+x,
 *   truncated to degree deg, for real coefficients.
 *   Returns the sum of all errors. */

double test_complex ( int deg );
/*
 * DESCRIPTION :
 *   Multiplies the power series of 1/(1-x) with 1+x,
 *   truncated to degree deg, for complex coefficients.
 *   Returns the sum of all errors. */

double test_real_exponential ( int deg );
/*
 * DESCRIPTION :
 *   Multiplies the power series for exp(x) with exp(-x)
 *   for some random x in [-1,+1], for real coefficients
 *   of a series of degree truncated to deg.
 *   Returns the sum of all errors. */

double test_complex_exponential ( int deg );
/*
 * DESCRIPTION :
 *   Multiplies the power series for exp(x) with exp(-x)
 *   for some random complex number on the unit circle,
 *   for series of degree truncated to deg.
 *   Returns the sum of all errors. */

int main ( void )
{
   const int timevalue = time(NULL); // for a random seed
   srand(timevalue);

   ios_base::fmtflags f(cout.flags()); // to restore format flags
 
   cout << "Give a degree larger than one : ";
   int deg; cin >> deg;

   if(deg > 0) 
   {
      double realerror1 = test_real(deg);

      cout.flags(f);
      double complexerror1 = test_complex(deg);

      double realerror2 = test_real_exponential(deg);
      double complexerror2 = test_complex_exponential(deg);

      const double tol = 1.0e-12;

      cout << scientific << setprecision(2);
      cout << "First test on real data, sum of all errors : ";
      cout << realerror1;
      if(realerror1 < tol)
         cout << "  pass." << endl;
      else
         cout << "  fail!" << endl;

      cout << "First test on complex data, sum of all errors : ";
      cout << complexerror1;
      if(complexerror1 < tol)
         cout << "  pass." << endl;
      else
         cout << "  fail!" << endl;

      cout << "Second test on real data, sum of all errors : ";
      cout << realerror2;
      if(realerror2 < tol)
         cout << "  pass." << endl;
      else
         cout << "  fail!" << endl;

      cout << "Second test on complex data, sum of all errors : ";
      cout << complexerror2;
      if(complexerror2 < tol)
         cout << "  pass." << endl;
      else
         cout << "  fail!" << endl;
   }
   return 0;
}

double test_real ( int deg )
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

   double err = 0.0;

   for(int k=0; k<=deg; k++)
   {
      cout << "  z[" << k << "] : " << z_d[k];
      if((k+1) % 4 == 0) cout << endl;
      err = err + abs(z_h[k] - z_d[k]);
   }
   cout << endl;

   cout << scientific << setprecision(16);
   cout << "the error : " << err << endl;

   return err;
}

double test_complex ( int deg )
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

   GPU_cmplx_product(xre,xim,yre,yim,zre_d,zim_d,deg,1,deg+1,2);

   cout << "GPU computed product :" << endl;

   double err = 0.0;

   for(int k=0; k<=deg; k++)
   {
      cout << "  zre[" << k << "] : " << zre_d[k];
      cout << "  zim[" << k << "] : " << zim_d[k];
      if((k+1) % 2 == 0) cout << endl;
      err = err + abs(zre_h[k] - zre_d[k]) + abs(zim_h[k] - zim_d[k]);
   }
   cout << endl;

   cout << scientific << setprecision(16);
   cout << "the error : " << err << endl;

   return err;
}

double test_real_exponential ( int deg )
{
   double *x = new double[deg+1];
   double *y = new double[deg+1];
   double *z_h = new double[deg+1];
   double *z_d = new double[deg+1];
   double r;

   random_dbl_exponentials(deg,&r,x,y);

   CPU_dbl_product(deg,x,y,z_h);

   cout << scientific << setprecision(16);
   cout << "Series of exp(x)*exp(-x), for x = " << r << endl;

   double sum = 0.0;
   for(int k=0; k<=deg; k++) sum = sum + z_h[k];
   cout << "Summation of all coefficients of the product ..." << endl;
   cout << "  sum : " << sum << endl;

   GPU_dbl_product(x,y,z_d,deg,1,deg+1);

   double err = 0.0;

   sum = 0.0;
   for(int k=0; k<=deg; k++)
   {
      sum = sum + z_d[k];
      err = err + abs(z_h[k] - z_d[k]);
   }
   cout << "Summation of all coefficients of the GPU computed product ..."
        << endl;
   cout << "  sum : " << sum << endl;

   cout << scientific << setprecision(16);
   cout << "the error : " << err << endl;

   return err;
}

double test_complex_exponential ( int deg )
{
   double *xre = new double[deg+1];
   double *xim = new double[deg+1];
   double *yre = new double[deg+1];
   double *yim = new double[deg+1];
   double *zre_h = new double[deg+1];
   double *zim_h = new double[deg+1];
   double *zre_d = new double[deg+1];
   double *zim_d = new double[deg+1];
   double cr,sr;

   random_cmplx_exponentials(deg,&cr,&sr,xre,xim,yre,yim);

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

   GPU_cmplx_product(xre,xim,yre,yim,zre_d,zim_d,deg,1,deg+1,2);

   double err = 0.0;

   sumre = 0.0; sumim = 0.0;
   for(int k=0; k<=deg; k++) 
   {
      sumre = sumre + zre_d[k];
      sumim = sumim + zim_d[k];
      err = err + abs(zre_h[k] - zre_d[k]) + abs(zim_h[k] - zim_d[k]);
   }
   cout << "Summation of all coefficients of the GPU computed product ..."
        << endl;
   cout << "  sumre : " << sumre << endl;
   cout << "  sumim : " << sumim << endl;

   cout << scientific << setprecision(16);
   cout << "the error : " << err << endl;

   return err;
}
