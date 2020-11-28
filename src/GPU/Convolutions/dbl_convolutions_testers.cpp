/* The file dbl_convolutions_testers.cpp contains the definitions of
 * functions to test the product of two series in double precision. */

#include <ctime>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector_types.h>
#include "random_numbers.h"
#include "random_series.h"
#include "dbl_convolutions_host.h"
#include "dbl_convolutions_kernels.h"
#include "dbl_convolutions_testers.h"

using namespace std;

int main_dbl_test ( int seed, int deg, int vrblvl )
{
   int fail,seedused;

   if(seed != 0)
   {
      srand(seed);
      seedused = seed;
   }
   else
   {
      const int timevalue = time(NULL); // for a random seed
      srand(timevalue);
      seedused = timevalue;
   }
   ios_base::fmtflags f(cout.flags()); // to restore format flags

   if(deg > 0) 
   {
      double realerror1 = test_dbl_real(deg,vrblvl-1);
      double realerror2 = test_dbl_real_random(deg,vrblvl-1);

      cout.flags(f);
      double complexerror1 = test_dbl_complex(deg,vrblvl-1);
      double complexerror2 = test_dbl_complex_random(deg,vrblvl-1);

      double realerror3 = test_dbl_real_exponential(deg,vrblvl-1);
      double complexerror3 = test_dbl_complex_exponential(deg,vrblvl-1);

      const double tol = 1.0e-12;

      fail = int(realerror1 > tol)
           + int(realerror2 > tol)
           + int(realerror3 > tol)
           + int(complexerror1 > tol)
           + int(complexerror2 > tol)
           + int(complexerror3 > tol);

      if(vrblvl > 0)
      {
         cout << scientific << setprecision(2);
         cout << "-> First test on real data, sum of all errors : ";
         cout << realerror1;
         if(realerror1 < tol)
            cout << "  pass." << endl;
         else
            cout << "  fail!" << endl;

         cout << scientific << setprecision(2);
         cout << "-> Second test on real data, sum of all errors : ";
         cout << realerror2;
         if(realerror2 < tol)
            cout << "  pass." << endl;
         else
            cout << "  fail!" << endl;

         cout << "-> Second test on complex data, sum of all errors : ";
         cout << complexerror2;
         if(complexerror2 < tol)
            cout << "  pass." << endl;
         else
            cout << "  fail!" << endl;

         cout << "-> First test on complex data, sum of all errors : ";
         cout << complexerror1;
         if(complexerror1 < tol)
            cout << "  pass." << endl;
         else
            cout << "  fail!" << endl;

         cout << "-> Third test on real data, sum of all errors : ";
         cout << realerror3;
         if(realerror3 < tol)
            cout << "  pass." << endl;
         else
            cout << "  fail!" << endl;

         cout << "-> Third test on complex data, sum of all errors : ";
         cout << complexerror3;
         if(complexerror3 < tol)
            cout << "  pass." << endl;
         else
            cout << "  fail!" << endl;

         cout << "   Seed used : " <<  seedused << endl;
      }
   }
   return fail;
}

double test_dbl_real ( int deg, int verbose )
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

   if(verbose > 0)
   {
      cout << "Series of 1/(1-x) multiplied with 1-x :" << endl;

      for(int k=0; k<=deg; k++)
      {
         cout << "  z[" << k << "] : " << z_h[k];
         if((k+1) % 4 == 0) cout << endl;
      }
      cout << endl;
   }
   GPU_dbl_product(x,y,z_d,deg,1,deg+1,1);

   if(verbose > 0) cout << "GPU computed product :" << endl;

   double err = 0.0;

   for(int k=0; k<=deg; k++)
   {
      if(verbose > 0)
      {
         cout << "  z[" << k << "] : " << z_d[k];
         if((k+1) % 4 == 0) cout << endl;
      }
      err = err + abs(z_h[k] - z_d[k]);
   }
   if(verbose > 0)
   {
      cout << endl;
      cout << scientific << setprecision(16);
      cout << "the error : " << err << endl;
   }
   return err;
}

double test_dbl_real_random ( int deg, int verbose )
{
   double* x = new double[deg+1];
   double* y = new double[deg+1];
   double* z_h = new double[deg+1];
   double* z_d = new double[deg+1];

   for(int k=0; k<=deg; k++)
   {
      x[k] = random_double();
      y[k] = random_double();
   }

   CPU_dbl_product(deg,x,y,z_h);

   if(verbose > 0)
   {
      cout << "Product of two random real series : " << endl;
      cout << scientific << setprecision(16);

      for(int k=0; k<=deg; k++)
      {
         cout << "z[" << k << "] : " << z_h[k] << endl;
      }
   }
   GPU_dbl_product(x,y,z_d,deg,1,deg+1,1);

   if(verbose > 0) cout << "GPU computed product :" << endl;

   double err = 0.0;

   for(int k=0; k<=deg; k++)
   {
      if(verbose > 0)
      {
         cout << "z[" << k << "] : " << z_d[k] << endl;
      }
      err = err + abs(z_h[k] - z_d[k]);
   }
   if(verbose > 0)
   {
      cout << scientific << setprecision(16);
      cout << "the error : " << err << endl;
   }
   return err;
}

double test_dbl_complex ( int deg, int verbose )
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

   if(verbose > 0)
   {
      cout << "Series of 1/(1-x) multiplied with 1-x : " << endl;

      for(int k=0; k<=deg; k++)
      {
         cout << "  zre[" << k << "] : " << zre_h[k];
         cout << "  zim[" << k << "] : " << zim_h[k];
         if((k+1) % 2 == 0) cout << endl;
      }
      cout << endl;
   }
   GPU_cmplx_product(xre,xim,yre,yim,zre_d,zim_d,deg,1,deg+1,3);

   if(verbose > 0) cout << "GPU computed product :" << endl;

   double err = 0.0;

   for(int k=0; k<=deg; k++)
   {
      if(verbose > 0)
      {
         cout << "  zre[" << k << "] : " << zre_d[k];
         cout << "  zim[" << k << "] : " << zim_d[k];
         if((k+1) % 2 == 0) cout << endl;
      }
      err = err + abs(zre_h[k] - zre_d[k]) + abs(zim_h[k] - zim_d[k]);
   }
   if(verbose > 0)
   {
      cout << endl;

      cout << scientific << setprecision(16);
      cout << "the error : " << err << endl;
   }
   return err;
}

double test_dbl_complex_random ( int deg, int verbose )
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
      xre[k] = random_double(); xim[k] = random_double();
      yre[k] = random_double(); yim[k] = random_double();
   }
   CPU_cmplx_product(deg,xre,xim,yre,yim,zre_h,zim_h);

   if(verbose > 0)
   {
      cout << "Product of two random complex series :" << endl;

      for(int k=0; k<=deg; k++)
      {
         cout << "  zre[" << k << "] : " << zre_h[k];
         cout << "  zim[" << k << "] : " << zim_h[k];
         if((k+1) % 2 == 0) cout << endl;
      }
      cout << endl;
   }
   GPU_cmplx_product(xre,xim,yre,yim,zre_d,zim_d,deg,1,deg+1,3);

   if(verbose > 0) cout << "GPU computed product :" << endl;

   double err = 0.0;

   for(int k=0; k<=deg; k++)
   {
      if(verbose > 0)
      {
         cout << "  zre[" << k << "] : " << zre_d[k];
         cout << "  zim[" << k << "] : " << zim_d[k];
         if((k+1) % 2 == 0) cout << endl;
      }
      err = err + abs(zre_h[k] - zre_d[k]) + abs(zim_h[k] - zim_d[k]);
   }
   if(verbose > 0)
   {
      cout << endl;

      cout << scientific << setprecision(16);
      cout << "the error : " << err << endl;
   }
   return err;
}

double test_dbl_real_exponential ( int deg, int verbose )
{
   double *x = new double[deg+1];
   double *y = new double[deg+1];
   double *z_h = new double[deg+1];
   double *z_d = new double[deg+1];
   double r,sum;

   random_dbl_exponentials(deg,&r,x,y);

   CPU_dbl_product(deg,x,y,z_h);

   if(verbose > 0)
   {
      cout << scientific << setprecision(16);
      cout << "Series of exp(x)*exp(-x), for x = " << r << endl;

      sum = 0.0;
      for(int k=0; k<=deg; k++) sum = sum + z_h[k];
      cout << "Summation of all coefficients of the product ..." << endl;
      cout << "  sum : " << sum << endl;
   }
   GPU_dbl_product(x,y,z_d,deg,1,deg+1,1);

   double err = 0.0;

   if(verbose > 0) sum = 0.0;
   for(int k=0; k<=deg; k++)
   {
      if(verbose > 0) sum = sum + z_d[k];
      err = err + abs(z_h[k] - z_d[k]);
   }
   if(verbose > 0)
   {
      cout << "Summation of all coefficients of the GPU computed product ..."
            << endl;
      cout << "  sum : " << sum << endl;
      cout << "the error : " << err << endl;
   }
   return err;
}

double test_dbl_complex_exponential ( int deg, int verbose )
{
   double *xre = new double[deg+1];
   double *xim = new double[deg+1];
   double *yre = new double[deg+1];
   double *yim = new double[deg+1];
   double *zre_h = new double[deg+1];
   double *zim_h = new double[deg+1];
   double *zre_d = new double[deg+1];
   double *zim_d = new double[deg+1];
   double cr,sr,sumre,sumim;

   random_cmplx_exponentials(deg,&cr,&sr,xre,xim,yre,yim);

   CPU_cmplx_product(deg,xre,xim,yre,yim,zre_h,zim_h);

   if(verbose > 0)
   {
      cout << scientific << setprecision(16);
      cout << "Series of exp(x)*exp(-x), for xre = " << cr << endl;
      cout << "  and xim = " << sr << endl;

      sumre = 0.0; sumim = 0.0;

      for(int k=0; k<=deg; k++) 
      {
         sumre = sumre + zre_h[k];
         sumim = sumim + zim_h[k];
      }
      cout << "Summation of all coefficients of the product ..." << endl;
      cout << "  sumre : " << sumre << endl;
      cout << "  sumim : " << sumim << endl;
   }
   GPU_cmplx_product(xre,xim,yre,yim,zre_d,zim_d,deg,1,deg+1,3);

   double err = 0.0;

   if(verbose > 0)
   {
      sumre = 0.0; sumim = 0.0;
   }
   for(int k=0; k<=deg; k++) 
   {
      if(verbose > 0)
      {
         sumre = sumre + zre_d[k];
         sumim = sumim + zim_d[k];
      }
      err = err + abs(zre_h[k] - zre_d[k]) + abs(zim_h[k] - zim_d[k]);
   }
   if(verbose > 0)
   {
      cout << "Summation of all coefficients of the GPU computed product ..."
           << endl;
      cout << "  sumre : " << sumre << endl;
      cout << "  sumim : " << sumim << endl;
      cout << "the error : " << err << endl;
   }
   return err;
}
