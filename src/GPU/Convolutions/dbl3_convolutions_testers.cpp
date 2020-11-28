/* The file dbl3_convolutions_testers.cpp contains the definitions of
 * functions to test the product of two series in triple double precision. */

#include <ctime>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector_types.h>
#include "triple_double_functions.h"
#include "random3_vectors.h"
#include "random3_series.h"
#include "dbl3_convolutions_host.h"
#include "dbl3_convolutions_kernels.h"
#include "dbl3_convolutions_testers.h"

using namespace std;

int main_dbl3_test ( int seed, int deg, int vrblvl )
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
      double realerror1 = test_dbl3_real(deg,vrblvl-1);
      double realerror2 = test_dbl3_real_random(deg,vrblvl-1);

      cout.flags(f);
      double complexerror1 = test_dbl3_complex(deg,vrblvl-1);
      double complexerror2 = test_dbl3_complex_random(deg,vrblvl-1);

      double realerror3 = test_dbl3_real_exponential(deg,vrblvl-1);
      double complexerror3 = test_dbl3_complex_exponential(deg,vrblvl-1);

      const double tol = 1.0e-42;

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

         cout << "-> First test on complex data, sum of all errors : ";
         cout << complexerror1;
         if(complexerror1 < tol)
            cout << "  pass." << endl;
         else
            cout << "  fail!" << endl;

         cout << "-> Second test on complex data, sum of all errors : ";
         cout << complexerror2;
         if(complexerror2 < tol)
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

double test_dbl3_real ( int deg, int verbose )
{
   double* xhi = new double[deg+1];
   double* xmi = new double[deg+1];
   double* xlo = new double[deg+1];
   double* yhi = new double[deg+1];
   double* ymi = new double[deg+1];
   double* ylo = new double[deg+1];
   double* zhi_h = new double[deg+1];
   double* zmi_h = new double[deg+1];
   double* zlo_h = new double[deg+1];
   double* zhi_d = new double[deg+1];
   double* zmi_d = new double[deg+1];
   double* zlo_d = new double[deg+1];

   for(int k=0; k<=deg; k++)
   {
      xhi[k] = 1.0; xmi[k] = 0.0; xlo[k] = 0.0;
      yhi[k] = 0.0; ymi[k] = 0.0; ylo[k] = 0.0;
   }
   yhi[0] = 1.0; yhi[1] = -1.0;

   CPU_dbl3_product(deg,xhi,xmi,xlo,yhi,ymi,ylo,zhi_h,zmi_h,zlo_h);

   if(verbose > 0)
   {
      cout << "Series of 1/(1-x) multiplied with 1-x : " << endl;

      for(int k=0; k<=deg; k++)
      {
         cout << "zhi[" << k << "] : " << zhi_h[k];
         cout << "  zmi[" << k << "] : " << zmi_h[k];
         cout << "  zlo[" << k << "] : " << zlo_h[k] << endl;
      }
   }
   GPU_dbl3_product(xhi,xmi,xlo,yhi,ymi,ylo,zhi_d,zmi_d,zlo_d,deg,1,deg+1,1);

   if(verbose > 0) cout << "GPU computed product :" << endl;

   double err = 0.0;

   for(int k=0; k<=deg; k++)
   {
      if(verbose > 0)
      {
         cout << "zhi[" << k << "] : " << zhi_d[k];
         cout << "  zmi[" << k << "] : " << zmi_d[k];
         cout << "  zlo[" << k << "] : " << zlo_d[k] << endl;
      }
      err = err
          + abs(zhi_h[k] - zhi_d[k]) + abs(zmi_h[k] - zmi_d[k])
          + abs(zlo_h[k] - zlo_d[k]);
   }
   if(verbose > 0)
   {
      cout << endl;

      cout << scientific << setprecision(16);
      cout << "the error : " << err << endl;
   }
   return err;
}

double test_dbl3_real_random ( int deg, int verbose )
{
   double* xhi = new double[deg+1];
   double* xmi = new double[deg+1];
   double* xlo = new double[deg+1];
   double* yhi = new double[deg+1];
   double* ymi = new double[deg+1];
   double* ylo = new double[deg+1];
   double* zhi_h = new double[deg+1];
   double* zmi_h = new double[deg+1];
   double* zlo_h = new double[deg+1];
   double* zhi_d = new double[deg+1];
   double* zmi_d = new double[deg+1];
   double* zlo_d = new double[deg+1];

   for(int k=0; k<=deg; k++)
   {
      random_triple_double(&xhi[k],&xmi[k],&xlo[k]);
      random_triple_double(&yhi[k],&ymi[k],&ylo[k]);
   }

   CPU_dbl3_product(deg,xhi,xmi,xlo,yhi,ymi,ylo,zhi_h,zmi_h,zlo_h);

   if(verbose > 0)
   {
      cout << "Product of two random real series : " << endl;
      cout << scientific << setprecision(16);

      for(int k=0; k<=deg; k++)
      {
         cout << "zhi[" << k << "] : " << zhi_h[k];
         cout << "  zmi[" << k << "] : " << zmi_h[k] << endl;
         cout << "  zlo[" << k << "] : " << zlo_h[k] << endl;
      }
   }
   GPU_dbl3_product(xhi,xmi,xlo,yhi,ymi,ylo,zhi_d,zmi_d,zlo_d,deg,1,deg+1,1);

   if(verbose > 0) cout << "GPU computed product :" << endl;

   double err = 0.0;

   for(int k=0; k<=deg; k++)
   {
      if(verbose > 0)
      {
         cout << "zhi[" << k << "] : " << zhi_d[k];
         cout << "  zmi[" << k << "] : " << zmi_d[k] << endl;
         cout << "  zlo[" << k << "] : " << zlo_d[k] << endl;
      }
      err = err
          + abs(zhi_h[k] - zhi_d[k]) + abs(zmi_h[k] - zmi_d[k])
          + abs(zlo_h[k] - zlo_d[k]);
   }
   if(verbose > 0)
   {
      cout << scientific << setprecision(16);
      cout << "the error : " << err << endl;
   }
   return err;
}

double test_dbl3_complex ( int deg, int verbose )
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
   double* zrehi_h = new double[deg+1];
   double* zremi_h = new double[deg+1];
   double* zrelo_h = new double[deg+1];
   double* zimhi_h = new double[deg+1];
   double* zimmi_h = new double[deg+1];
   double* zimlo_h = new double[deg+1];
   double* zrehi_d = new double[deg+1];
   double* zremi_d = new double[deg+1];
   double* zrelo_d = new double[deg+1];
   double* zimhi_d = new double[deg+1];
   double* zimmi_d = new double[deg+1];
   double* zimlo_d = new double[deg+1];

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
                          zrehi_h,zremi_h,zrelo_h,zimhi_h,zimmi_h,zimlo_h);

   if(verbose > 0)
   {
      cout << "Series of 1/(1-x) multiplied with 1-x :" << endl;

      for(int k=0; k<=deg; k++)
      {
         cout << "zrehi[" << k << "] : " << zrehi_h[k];
         cout << "  zremi[" << k << "] : " << zremi_h[k];
         cout << "  zrelo[" << k << "] : " << zrelo_h[k] << endl;
         cout << "zimhi[" << k << "] : " << zimhi_h[k];
         cout << "  zimmi[" << k << "] : " << zimmi_h[k];
         cout << "  zimlo[" << k << "] : " << zimlo_h[k] << endl;
      }
   }
   GPU_cmplx3_product
      (xrehi,xremi,xrelo,ximhi,ximmi,ximlo,
       yrehi,yremi,yrelo,yimhi,yimmi,yimlo,
       zrehi_d,zremi_d,zrelo_d,zimhi_d,zimmi_d,zimlo_d,deg,1,deg+1,3);

   if(verbose > 0) cout << "GPU computed product :" << endl;

   double err = 0.0;

   for(int k=0; k<=deg; k++)
   {
      if(verbose > 0)
      {
         cout << "zrehi[" << k << "] : " << zrehi_d[k];
         cout << "  zremi[" << k << "] : " << zremi_d[k];
         cout << "  zrelo[" << k << "] : " << zrelo_d[k] << endl;
         cout << "zimhi[" << k << "] : " << zimhi_d[k];
         cout << "  zimmi[" << k << "] : " << zimmi_d[k];
         cout << "  zimlo[" << k << "] : " << zimlo_d[k] << endl;
      }
      err = err
          + abs(zrehi_h[k] - zrehi_d[k]) + abs(zremi_h[k] - zremi_d[k])
          + abs(zrelo_h[k] - zrelo_d[k])
          + abs(zimhi_h[k] - zimhi_d[k]) + abs(zimmi_h[k] - zimmi_d[k])
          + abs(zimlo_h[k] - zimlo_d[k]);
   }
   if(verbose > 0)
   {
      cout << scientific << setprecision(16);
      cout << "the error : " << err << endl;
   }
   return err;
}

double test_dbl3_complex_random ( int deg, int verbose )
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
   double* zrehi_h = new double[deg+1];
   double* zremi_h = new double[deg+1];
   double* zrelo_h = new double[deg+1];
   double* zimhi_h = new double[deg+1];
   double* zimmi_h = new double[deg+1];
   double* zimlo_h = new double[deg+1];
   double* zrehi_d = new double[deg+1];
   double* zremi_d = new double[deg+1];
   double* zrelo_d = new double[deg+1];
   double* zimhi_d = new double[deg+1];
   double* zimmi_d = new double[deg+1];
   double* zimlo_d = new double[deg+1];

   for(int k=0; k<=deg; k++)
   {
      random_triple_double(&xrehi[k],&xremi[k],&xrelo[k]);
      random_triple_double(&ximhi[k],&ximmi[k],&ximlo[k]);
      random_triple_double(&yrehi[k],&yremi[k],&yrelo[k]);
      random_triple_double(&yimhi[k],&yimmi[k],&yimlo[k]);
   }
   CPU_cmplx3_product(deg,xrehi,xremi,xrelo,ximhi,ximmi,ximlo,
                          yrehi,yremi,yrelo,yimhi,yimmi,yimlo,
                          zrehi_h,zremi_h,zrelo_h,zimhi_h,zimmi_h,zimlo_h);

   if(verbose > 0)
   {
      cout << "Product of two random complex series :" << endl;

      for(int k=0; k<=deg; k++)
      {
         cout << "zrehi[" << k << "] : " << zrehi_h[k];
         cout << "  zremi[" << k << "] : " << zremi_h[k];
         cout << "  zrelo[" << k << "] : " << zrelo_h[k] << endl;
         cout << "zimhi[" << k << "] : " << zimhi_h[k];
         cout << "  zimmi[" << k << "] : " << zimmi_h[k];
         cout << "  zimlo[" << k << "] : " << zimlo_h[k] << endl;
      }
   }
   GPU_cmplx3_product
      (xrehi,xremi,xrelo,ximhi,ximmi,ximlo,
       yrehi,yremi,yrelo,yimhi,yimmi,yimlo,
       zrehi_d,zremi_d,zrelo_d,zimhi_d,zimmi_d,zimlo_d,deg,1,deg+1,3);

   if(verbose > 0) cout << "GPU computed product :" << endl;

   double err = 0.0;

   for(int k=0; k<=deg; k++)
   {
      if(verbose > 0)
      {
         cout << "zrehi[" << k << "] : " << zrehi_d[k];
         cout << "  zremi[" << k << "] : " << zremi_d[k];
         cout << "  zrelo[" << k << "] : " << zrelo_d[k] << endl;
         cout << "zimhi[" << k << "] : " << zimhi_d[k];
         cout << "  zimmi[" << k << "] : " << zimmi_d[k];
         cout << "  zimlo[" << k << "] : " << zimlo_d[k] << endl;
      }
      err = err
          + abs(zrehi_h[k] - zrehi_d[k]) + abs(zremi_h[k] - zremi_d[k])
          + abs(zrelo_h[k] - zrelo_d[k])
          + abs(zimhi_h[k] - zimhi_d[k]) + abs(zimmi_h[k] - zimmi_d[k])
          + abs(zimlo_h[k] - zimlo_d[k]);
   }
   if(verbose > 0)
   {
      cout << scientific << setprecision(16);
      cout << "the error : " << err << endl;
   }
   return err;
}

double test_dbl3_real_exponential ( int deg, int verbose )
{
   double *xhi = new double[deg+1];
   double *xmi = new double[deg+1];
   double *xlo = new double[deg+1];
   double *yhi = new double[deg+1];
   double *ymi = new double[deg+1];
   double *ylo = new double[deg+1];
   double *zhi_h = new double[deg+1];
   double *zmi_h = new double[deg+1];
   double *zlo_h = new double[deg+1];
   double *zhi_d = new double[deg+1];
   double *zmi_d = new double[deg+1];
   double *zlo_d = new double[deg+1];
   double rhi,rmi,rlo;
   double sumhi,summi,sumlo;

   random_dbl3_exponentials(deg,&rhi,&rmi,&rlo,xhi,xmi,xlo,yhi,ymi,ylo);

   CPU_dbl3_product(deg,xhi,xmi,xlo,yhi,ymi,ylo,zhi_h,zmi_h,zlo_h);

   if(verbose > 0)
   {
      cout << scientific << setprecision(16);

      cout << "Product of series of exp(x) with series of exp(-x)," << endl;
      cout << "  for xhi = " << rhi << endl;
      cout << "      xmi = " << rmi << endl;
      cout << "  and xlo = " << rlo << endl;

      sumhi = 0.0; summi = 0.0; sumlo = 0.0;

      for(int k=0; k<=deg; k++)
         tdf_inc(&sumhi,&summi,&sumlo,zhi_h[k],zmi_h[k],zlo_h[k]);

      cout << "Summation of all coefficients in the product ..." << endl;
      cout << "  high part of sum : " << sumhi << endl;
      cout << "middle part of sum : " << summi << endl;
      cout << "   low part of sum : " << sumlo << endl;
   }
   GPU_dbl3_product(xhi,xmi,xlo,yhi,ymi,ylo,zhi_d,zmi_d,zlo_d,deg,1,deg+1,1);

   if(verbose > 0)
   {
      sumhi = 0.0; summi = 0.0; sumlo = 0.0;
   }
   double err = 0.0;

   for(int k=0; k<=deg; k++)
   {
      if(verbose > 0)
         tdf_inc(&sumhi,&summi,&sumlo,zhi_d[k],zmi_d[k],zlo_d[k]);

      err = err
          + abs(zhi_h[k] - zhi_d[k]) + abs(zmi_h[k] - zmi_d[k])
          + abs(zlo_h[k] - zlo_d[k]);
   }
   if(verbose > 0)
   {
      cout << "Summation of all coefficients in the GPU computed product ..."
           << endl;
      cout << "  high part of sum : " << sumhi << endl;
      cout << "middle part of sum : " << summi << endl;
      cout << "   low part of sum : " << sumlo << endl;

      cout << "the error : " << err << endl;
   }
   return err;
}

double test_dbl3_complex_exponential ( int deg, int verbose )
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
   double* zrehi_h = new double[deg+1];
   double* zremi_h = new double[deg+1];
   double* zrelo_h = new double[deg+1];
   double* zimhi_h = new double[deg+1];
   double* zimmi_h = new double[deg+1];
   double* zimlo_h = new double[deg+1];
   double* zrehi_d = new double[deg+1];
   double* zremi_d = new double[deg+1];
   double* zrelo_d = new double[deg+1];
   double* zimhi_d = new double[deg+1];
   double* zimmi_d = new double[deg+1];
   double* zimlo_d = new double[deg+1];
   double rndrehi,rndremi,rndrelo;
   double rndimhi,rndimmi,rndimlo;
   double tmphi,tmpmi,tmplo;
   double sumrehi,sumremi,sumrelo,sumimhi,sumimmi,sumimlo;

   random_cmplx3_exponentials
      (deg,&rndrehi,&rndremi,&rndrelo,&rndimhi,&rndimmi,&rndimlo,
           xrehi,xremi,xrelo,ximhi,ximmi,ximlo,
           yrehi,yremi,yrelo,yimhi,yimmi,yimlo);

   CPU_cmplx3_product(deg,xrehi,xremi,xrelo,ximhi,ximmi,ximlo,
                          yrehi,yremi,yrelo,yimhi,yimmi,yimlo,
                          zrehi_h,zremi_h,zrelo_h,zimhi_h,zimmi_h,zimlo_h);

   if(verbose > 0)
   {
      cout << scientific << setprecision(16);

      cout << "Product of series of exp(x) with series of exp(-x)," << endl;
      cout << "  for xrehi = " << rndrehi;
      cout << "      xremi = " << rndremi << endl;
      cout << "  and xrelo = " << rndrelo << endl;
      cout << "  for ximhi = " << rndimhi;
      cout << "      ximmi = " << rndimmi << endl;
      cout << "  and ximlo = " << rndimlo << endl;

      sumrehi = 0.0; sumremi = 0.0; sumrelo = 0.0;
      sumimhi = 0.0; sumimmi = 0.0; sumimlo = 0.0;

      for(int k=0; k<=deg; k++) 
      {
         tdf_inc(&sumrehi,&sumremi,&sumrelo,zrehi_h[k],zremi_h[k],zrelo_h[k]);
         tdf_inc(&sumimhi,&sumimmi,&sumimlo,zimhi_h[k],zimmi_h[k],zimlo_h[k]);
      }
      cout << "Summation of all coefficients of the product ..." << endl;
      cout << "  sumrehi : " << sumrehi;
      cout << "  sumremi : " << sumremi << endl;
      cout << "  sumrelo : " << sumrelo << endl;
      cout << "  sumimhi : " << sumimhi;
      cout << "  sumimmi : " << sumimmi << endl;
      cout << "  sumimlo : " << sumimlo << endl;
   }
   GPU_cmplx3_product
      (xrehi,xremi,xrelo,ximhi,ximmi,ximlo,
       yrehi,yremi,yrelo,yimhi,yimmi,yimlo,
       zrehi_d,zremi_d,zrelo_d,zimhi_d,zimmi_d,zimlo_d,deg,1,deg+1,3);

   if(verbose > 0)
   {
      sumrehi = 0.0; sumremi = 0.0; sumrelo = 0.0;
      sumimhi = 0.0; sumimmi = 0.0; sumimlo = 0.0;
   }
   double err = 0.0;

   for(int k=0; k<=deg; k++) 
   {
      if(verbose > 0)
      {
         tdf_inc(&sumrehi,&sumremi,&sumrelo,zrehi_d[k],zremi_d[k],zrelo_d[k]);
         tdf_inc(&sumimhi,&sumimmi,&sumimlo,zimhi_d[k],zimmi_d[k],zimlo_d[k]);
      }
      err = err
          + abs(zrehi_h[k] - zrehi_d[k]) + abs(zremi_h[k] - zremi_d[k])
          + abs(zrelo_h[k] - zrelo_d[k])
          + abs(zimhi_h[k] - zimhi_d[k]) + abs(zimmi_h[k] - zimmi_d[k])
          + abs(zimlo_h[k] - zimlo_d[k]);
   }
   if(verbose > 0)
   {
      cout << "Summation of all coefficients of the GPU computed product ..."
           << endl;
      cout << "  sumrehi : " << sumrehi;
      cout << "  sumremi : " << sumremi << endl;
      cout << "  sumrelo : " << sumrelo << endl;
      cout << "  sumimhi : " << sumimhi;
      cout << "  sumimmi : " << sumimmi << endl;
      cout << "  sumimlo : " << sumimlo << endl;

      cout << "the error : " << err << endl;
   }
   return err;
}
