/* The file dbl4_convolutions_testers.cpp contains the definitions of
 * functions to test the product of two series in quad double precision. */

#include <ctime>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector_types.h>
#include "quad_double_functions.h"
#include "random4_vectors.h"
#include "random4_series.h"
#include "dbl4_convolutions_host.h"
#include "dbl4_convolutions_kernels.h"
#include "dbl4_convolutions_testers.h"

using namespace std;

int main_dbl4_test ( int seed, int deg, int vrblvl )
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
      double realerror1 = test_dbl4_real(deg,vrblvl-1);
      double realerror2 = test_dbl4_real_random(deg,vrblvl-1);

      cout.flags(f);
      double complexerror1 = test_dbl4_complex(deg,vrblvl-1);
      double complexerror2 = test_dbl4_complex_random(deg,vrblvl-1);

      double realerror3 = test_dbl4_real_exponential(deg,vrblvl-1);
      double complexerror3 = test_dbl4_complex_exponential(deg,vrblvl-1);

      const double tol = 1.0e-60;

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

double test_dbl4_real ( int deg, int verbose )
{
   double* xhihi = new double[deg+1];
   double* xlohi = new double[deg+1];
   double* xhilo = new double[deg+1];
   double* xlolo = new double[deg+1];
   double* yhihi = new double[deg+1];
   double* ylohi = new double[deg+1];
   double* yhilo = new double[deg+1];
   double* ylolo = new double[deg+1];
   double* zhihi_h = new double[deg+1];
   double* zlohi_h = new double[deg+1];
   double* zhilo_h = new double[deg+1];
   double* zlolo_h = new double[deg+1];
   double* zhihi_d = new double[deg+1];
   double* zlohi_d = new double[deg+1];
   double* zhilo_d = new double[deg+1];
   double* zlolo_d = new double[deg+1];

   for(int k=0; k<=deg; k++)
   {
      xhihi[k] = 1.0; xlohi[k] = 0.0;
      xhilo[k] = 0.0; xlolo[k] = 0.0;
      yhihi[k] = 0.0; ylohi[k] = 0.0;
      yhilo[k] = 0.0; ylolo[k] = 0.0;
   }
   yhihi[0] = 1.0; yhihi[1] = -1.0;

   CPU_dbl4_product(deg,xhihi,xlohi,xhilo,xlolo,yhihi,ylohi,yhilo,ylolo,
                        zhihi_h,zlohi_h,zhilo_h,zlolo_h);
   if(verbose > 0)
   {
      cout << "Series of 1/(1-x) multiplied with 1-x : " << endl;

      for(int k=0; k<=deg; k++)
      {
         cout << "zhihi[" << k << "] : " << zhihi_h[k];
         cout << "  zlohi[" << k << "] : " << zlohi_h[k];
         cout << "  zhilo[" << k << "] : " << zhilo_h[k];
         cout << "  zlolo[" << k << "] : " << zlolo_h[k] << endl;
      }
   }
   GPU_dbl4_product(xhihi,xlohi,xhilo,xlolo,
                    yhihi,ylohi,yhilo,ylolo,
                    zhihi_d,zlohi_d,zhilo_d,zlolo_d,deg,1,deg+1,1);

   if(verbose > 0) cout << "GPU computed product :" << endl;

   double err = 0.0;

   for(int k=0; k<=deg; k++)
   {
      if(verbose > 0)
      {
         cout << "zhihi[" << k << "] : " << zhihi_d[k];
         cout << "  zlohi[" << k << "] : " << zlohi_d[k];
         cout << "  zhilo[" << k << "] : " << zhilo_d[k];
         cout << "  zlolo[" << k << "] : " << zlolo_d[k] << endl;
      }
      err = err
          + abs(zhihi_h[k] - zhihi_d[k]) + abs(zlohi_h[k] - zlohi_d[k])
          + abs(zhilo_h[k] - zhilo_d[k]) + abs(zlolo_h[k] - zlolo_d[k]);
   }
   if(verbose > 0)
   {
      cout << scientific << setprecision(16);
      cout << "the error : " << err << endl;
   }
   return err;
}

double test_dbl4_real_random ( int deg, int verbose )
{
   double* xhihi = new double[deg+1];
   double* xlohi = new double[deg+1];
   double* xhilo = new double[deg+1];
   double* xlolo = new double[deg+1];
   double* yhihi = new double[deg+1];
   double* ylohi = new double[deg+1];
   double* yhilo = new double[deg+1];
   double* ylolo = new double[deg+1];
   double* zhihi_h = new double[deg+1];
   double* zlohi_h = new double[deg+1];
   double* zhilo_h = new double[deg+1];
   double* zlolo_h = new double[deg+1];
   double* zhihi_d = new double[deg+1];
   double* zlohi_d = new double[deg+1];
   double* zhilo_d = new double[deg+1];
   double* zlolo_d = new double[deg+1];

   for(int k=0; k<=deg; k++)
   {
      random_quad_double
         (&xhihi[k],&xlohi[k],&xhilo[k],&xlolo[k]);
      random_quad_double
         (&yhihi[k],&ylohi[k],&yhilo[k],&ylolo[k]);
   }
   CPU_dbl4_product
      (deg,xhihi,xlohi,xhilo,xlolo,yhihi,ylohi,yhilo,ylolo,
           zhihi_h,zlohi_h,zhilo_h,zlolo_h);

   if(verbose > 0)
   {
      cout << "Product of two random real series : " << endl;
      cout << scientific << setprecision(16);

      for(int k=0; k<=deg; k++)
      {
         cout << "zhihi[" << k << "] : " << zhihi_h[k];
         cout << "  zlohi[" << k << "] : " << zlohi_h[k] << endl;
         cout << "zhilo[" << k << "] : " << zhilo_h[k];
         cout << "  zlolo[" << k << "] : " << zlolo_h[k] << endl;
      }
   }
   GPU_dbl4_product
      (xhihi,xlohi,xhilo,xlolo,yhihi,ylohi,yhilo,ylolo,
       zhihi_d,zlohi_d,zhilo_d,zlolo_d,deg,1,deg+1,1);

   if(verbose > 0) cout << "GPU computed product :" << endl;

   double err = 0.0;

   for(int k=0; k<=deg; k++)
   {
      if(verbose > 0)
      {
         cout << "zhihi[" << k << "] : " << zhihi_d[k];
         cout << "  zlohi[" << k << "] : " << zlohi_d[k] << endl;
         cout << "zhilo[" << k << "] : " << zhilo_d[k];
         cout << "  zlolo[" << k << "] : " << zlolo_d[k] << endl;
      }
      err = err
          + abs(zhihi_h[k] - zhihi_d[k])
          + abs(zlohi_h[k] - zlohi_d[k])
          + abs(zhilo_h[k] - zhilo_d[k])
          + abs(zlolo_h[k] - zlolo_d[k]);
   }
   if(verbose > 0)
   {
      cout << scientific << setprecision(16);
      cout << "the error : " << err << endl;
   }
   return err;
}

double test_dbl4_complex ( int deg, int verbose )
{
   double* xrehihi = new double[deg+1];
   double* xrelohi = new double[deg+1];
   double* xrehilo = new double[deg+1];
   double* xrelolo = new double[deg+1];
   double* ximhihi = new double[deg+1];
   double* ximlohi = new double[deg+1];
   double* ximhilo = new double[deg+1];
   double* ximlolo = new double[deg+1];
   double* yrehihi = new double[deg+1];
   double* yrelohi = new double[deg+1];
   double* yrehilo = new double[deg+1];
   double* yrelolo = new double[deg+1];
   double* yimhihi = new double[deg+1];
   double* yimlohi = new double[deg+1];
   double* yimhilo = new double[deg+1];
   double* yimlolo = new double[deg+1];
   double* zrehihi_h = new double[deg+1];
   double* zrelohi_h = new double[deg+1];
   double* zrehilo_h = new double[deg+1];
   double* zrelolo_h = new double[deg+1];
   double* zimhihi_h = new double[deg+1];
   double* zimlohi_h = new double[deg+1];
   double* zimhilo_h = new double[deg+1];
   double* zimlolo_h = new double[deg+1];
   double* zrehihi_d = new double[deg+1];
   double* zrelohi_d = new double[deg+1];
   double* zrehilo_d = new double[deg+1];
   double* zrelolo_d = new double[deg+1];
   double* zimhihi_d = new double[deg+1];
   double* zimlohi_d = new double[deg+1];
   double* zimhilo_d = new double[deg+1];
   double* zimlolo_d = new double[deg+1];

   for(int k=0; k<=deg; k++)
   {
      xrehihi[k] = 1.0; xrelohi[k] = 0.0;
      xrehilo[k] = 0.0; xrelolo[k] = 0.0;
      ximhihi[k] = 0.0; ximlohi[k] = 0.0;
      ximhilo[k] = 0.0; ximlolo[k] = 0.0;
      yrehihi[k] = 0.0; yrelohi[k] = 0.0;
      yrehilo[k] = 0.0; yrelolo[k] = 0.0;
      yimhihi[k] = 0.0; yimlohi[k] = 0.0;
      yimhilo[k] = 0.0; yimlolo[k] = 0.0;
   }
   yrehihi[0] = 1.0; yrehihi[1] = -1.0;

   CPU_cmplx4_product
      (deg,xrehihi,xrelohi,xrehilo,xrelolo,ximhihi,ximlohi,ximhilo,ximlolo,
           yrehihi,yrelohi,yrehilo,yimlolo,yimhihi,yimlohi,yimhilo,yimlolo,
           zrehihi_h,zrelohi_h,zrehilo_h,zrelolo_h,
           zimhihi_h,zimlohi_h,zimhilo_h,zimlolo_h);

   if(verbose > 0)
   {
      cout << "Series of 1/(1-x) multiplied with 1-x : " << endl;

      for(int k=0; k<=deg; k++)
      {
         cout << "zrehihi[" << k << "] : " << zrehihi_h[k];
         cout << "  zrelohi[" << k << "] : " << zrelohi_h[k];
         cout << "  zrehilo[" << k << "] : " << zrehilo_h[k];
         cout << "  zrelolo[" << k << "] : " << zrelolo_h[k] << endl;
         cout << "zimhihi[" << k << "] : " << zimhihi_h[k];
         cout << "  zimlohi[" << k << "] : " << zimlohi_h[k];
         cout << "  zimhilo[" << k << "] : " << zimhilo_h[k];
         cout << "  zimlolo[" << k << "] : " << zimlolo_h[k] << endl;
      }
   }
   GPU_cmplx4_product
      (xrehihi,xrelohi,xrehilo,xrelolo,ximhihi,ximlohi,ximhilo,ximlolo,
       yrehihi,yrelohi,yrehilo,yrelolo,yimhihi,yimlohi,yimhilo,yimlolo,
       zrehihi_d,zrelohi_d,zrehilo_d,zrelolo_d,
       zimhihi_d,zimlohi_d,zimhilo_d,zimlolo_d,deg,1,deg+1,2);

   if(verbose > 0) cout << "GPU computed product :" << endl;

   double err = 0.0;

   for(int k=0; k<=deg; k++)
   {
      if(verbose > 0)
      {
         cout << "zrehihi[" << k << "] : " << zrehihi_d[k];
         cout << "  zrelohi[" << k << "] : " << zrelohi_d[k];
         cout << "  zrehilo[" << k << "] : " << zrehilo_d[k];
         cout << "  zrelolo[" << k << "] : " << zrelolo_d[k] << endl;
         cout << "zimhihi[" << k << "] : " << zimhihi_d[k];
         cout << "  zimlohi[" << k << "] : " << zimlohi_d[k];
         cout << "  zimhilo[" << k << "] : " << zimhilo_d[k];
         cout << "  zimlolo[" << k << "] : " << zimlolo_d[k] << endl;
      }
      err = err
          + abs(zrehihi_h[k] - zrehihi_d[k])
          + abs(zrelohi_h[k] - zrelohi_d[k])
          + abs(zrehilo_h[k] - zrehilo_d[k])
          + abs(zrelolo_h[k] - zrelolo_d[k])
          + abs(zimhihi_h[k] - zimhihi_d[k])
          + abs(zimlohi_h[k] - zimlohi_d[k])
          + abs(zimhilo_h[k] - zimhilo_d[k])
          + abs(zimlolo_h[k] - zimlolo_d[k]);
   }
   if(verbose > 0)
   {
      cout << scientific << setprecision(16);
      cout << "the error : " << err << endl;
   }
   return err;
}

double test_dbl4_complex_random ( int deg, int verbose )
{
   double* xrehihi = new double[deg+1];
   double* xrelohi = new double[deg+1];
   double* xrehilo = new double[deg+1];
   double* xrelolo = new double[deg+1];
   double* ximhihi = new double[deg+1];
   double* ximlohi = new double[deg+1];
   double* ximhilo = new double[deg+1];
   double* ximlolo = new double[deg+1];
   double* yrehihi = new double[deg+1];
   double* yrelohi = new double[deg+1];
   double* yrehilo = new double[deg+1];
   double* yrelolo = new double[deg+1];
   double* yimhihi = new double[deg+1];
   double* yimlohi = new double[deg+1];
   double* yimhilo = new double[deg+1];
   double* yimlolo = new double[deg+1];
   double* zrehihi_h = new double[deg+1];
   double* zrelohi_h = new double[deg+1];
   double* zrehilo_h = new double[deg+1];
   double* zrelolo_h = new double[deg+1];
   double* zimhihi_h = new double[deg+1];
   double* zimlohi_h = new double[deg+1];
   double* zimhilo_h = new double[deg+1];
   double* zimlolo_h = new double[deg+1];
   double* zrehihi_d = new double[deg+1];
   double* zrelohi_d = new double[deg+1];
   double* zrehilo_d = new double[deg+1];
   double* zrelolo_d = new double[deg+1];
   double* zimhihi_d = new double[deg+1];
   double* zimlohi_d = new double[deg+1];
   double* zimhilo_d = new double[deg+1];
   double* zimlolo_d = new double[deg+1];

   for(int k=0; k<=deg; k++)
   {
      random_quad_double
         (&xrehihi[k],&xrelohi[k],&xrehilo[k],&xrelolo[k]);
      random_quad_double
         (&ximhihi[k],&ximlohi[k],&ximhilo[k],&ximlolo[k]);
      random_quad_double
         (&yrehihi[k],&yrelohi[k],&yrehilo[k],&yrelolo[k]);
      random_quad_double
         (&yimhihi[k],&yimlohi[k],&yimhilo[k],&yimlolo[k]);
   }
   CPU_cmplx4_product
      (deg,xrehihi,xrelohi,xrehilo,xrelolo,ximhihi,ximlohi,ximhilo,ximlolo,
           yrehihi,yrelohi,yrehilo,yimlolo,yimhihi,yimlohi,yimhilo,yimlolo,
           zrehihi_h,zrelohi_h,zrehilo_h,zrelolo_h,
           zimhihi_h,zimlohi_h,zimhilo_h,zimlolo_h);

   if(verbose > 0)
   {
      cout << "Product of two random complex series : " << endl;

      for(int k=0; k<=deg; k++)
      {
         cout << "zrehihi[" << k << "] : " << zrehihi_h[k];
         cout << "  zrelohi[" << k << "] : " << zrelohi_h[k] << endl;
         cout << "zrehilo[" << k << "] : " << zrehilo_h[k];
         cout << "  zrelolo[" << k << "] : " << zrelolo_h[k] << endl;
         cout << "zimhihi[" << k << "] : " << zimhihi_h[k];
         cout << "  zimlohi[" << k << "] : " << zimlohi_h[k] << endl;
         cout << "zimhilo[" << k << "] : " << zimhilo_h[k];
         cout << "  zimlolo[" << k << "] : " << zimlolo_h[k] << endl;
      }
   }
   GPU_cmplx4_product
      (xrehihi,xrelohi,xrehilo,xrelolo,ximhihi,ximlohi,ximhilo,ximlolo,
       yrehihi,yrelohi,yrehilo,yrelolo,yimhihi,yimlohi,yimhilo,yimlolo,
       zrehihi_d,zrelohi_d,zrehilo_d,zrelolo_d,
       zimhihi_d,zimlohi_d,zimhilo_d,zimlolo_d,deg,1,deg+1,2);

   if(verbose > 0) cout << "GPU computed product :" << endl;

   double err = 0.0;

   for(int k=0; k<=deg; k++)
   {
      if(verbose > 0)
      {
         cout << "zrehihi[" << k << "] : " << zrehihi_d[k];
         cout << "  zrelohi[" << k << "] : " << zrelohi_d[k] << endl;
         cout << "zrehilo[" << k << "] : " << zrehilo_d[k];
         cout << "  zrelolo[" << k << "] : " << zrelolo_d[k] << endl;
         cout << "zimhihi[" << k << "] : " << zimhihi_d[k];
         cout << "  zimlohi[" << k << "] : " << zimlohi_d[k] << endl;
         cout << "zimhilo[" << k << "] : " << zimhilo_d[k];
         cout << "  zimlolo[" << k << "] : " << zimlolo_d[k] << endl;
      }
      err = err
          + abs(zrehihi_h[k] - zrehihi_d[k])
          + abs(zrelohi_h[k] - zrelohi_d[k])
          + abs(zrehilo_h[k] - zrehilo_d[k])
          + abs(zrelolo_h[k] - zrelolo_d[k])
          + abs(zimhihi_h[k] - zimhihi_d[k])
          + abs(zimlohi_h[k] - zimlohi_d[k])
          + abs(zimhilo_h[k] - zimhilo_d[k])
          + abs(zimlolo_h[k] - zimlolo_d[k]);
   }
   if(verbose > 0)
   {
      cout << scientific << setprecision(16);
      cout << "the error : " << err << endl;
   }
   return err;
}

double test_dbl4_real_exponential ( int deg, int verbose )
{
   double *xhihi = new double[deg+1];
   double *xlohi = new double[deg+1];
   double *xhilo = new double[deg+1];
   double *xlolo = new double[deg+1];
   double *yhihi = new double[deg+1];
   double *ylohi = new double[deg+1];
   double *yhilo = new double[deg+1];
   double *ylolo = new double[deg+1];
   double *zhihi_h = new double[deg+1];
   double *zlohi_h = new double[deg+1];
   double *zhilo_h = new double[deg+1];
   double *zlolo_h = new double[deg+1];
   double *zhihi_d = new double[deg+1];
   double *zlohi_d = new double[deg+1];
   double *zhilo_d = new double[deg+1];
   double *zlolo_d = new double[deg+1];
   double rhihi,rlohi,rhilo,rlolo;
   double sumhihi,sumlohi,sumhilo,sumlolo;

   random_dbl4_exponentials
      (deg,&rhihi,&rlohi,&rhilo,&rlolo,
           xhihi,xlohi,xhilo,xlolo,yhihi,ylohi,yhilo,ylolo);

   CPU_dbl4_product(deg,xhihi,xlohi,xhilo,xlolo,yhihi,ylohi,yhilo,ylolo,
                        zhihi_h,zlohi_h,zhilo_h,zlolo_h);
   if(verbose > 0)
   {
      cout << scientific << setprecision(16);

      cout << "Product of series of exp(x) with series of exp(-x)," << endl;
      cout << "  for xhihi = " << rhihi;
      cout << "      xlohi = " << rlohi << endl;
      cout << "      xhilo = " << rhilo;
      cout << "  and xlolo = " << rlolo << endl;

      sumhihi = 0.0; sumlohi = 0.0; sumhilo = 0.0; sumlolo = 0.0;

      for(int k=0; k<=deg; k++)
         qdf_inc(&sumhihi,&sumlohi,&sumhilo,&sumlolo,
                 zhihi_h[k],zlohi_h[k],zhilo_h[k],zlolo_h[k]);

      cout << "Summation of all coefficients in the product ..." << endl;
      cout << "    highest part of sum : " << sumhihi << endl;
      cout << "2nd highest part of sum : " << sumlohi << endl;
      cout << " 2nd lowest part of sum : " << sumhilo << endl;
      cout << "     lowest part of sum : " << sumlolo << endl;
   }
   GPU_dbl4_product(xhihi,xlohi,xhilo,xlolo,
                    yhihi,ylohi,yhilo,ylolo,
                    zhihi_d,zlohi_d,zhilo_d,zlolo_d,deg,1,deg+1,1);
   if(verbose > 0)
   {
      sumhihi = 0.0; sumlohi = 0.0; sumhilo = 0.0; sumlolo = 0.0;
   }
   double err = 0.0;

   for(int k=0; k<=deg; k++)
   {
      if(verbose > 0)
         qdf_inc(&sumhihi,&sumlohi,&sumhilo,&sumlolo,
                 zhihi_d[k],zlohi_d[k],zhilo_d[k],zlolo_d[k]);

      err = err
          + abs(zhihi_h[k] - zhihi_d[k]) + abs(zlohi_h[k] - zlohi_d[k])
          + abs(zhilo_h[k] - zhilo_d[k]) + abs(zlolo_h[k] - zlolo_d[k]);
   }
   if(verbose > 0)
   {
      cout << "Summation of all coefficients in the GPU computed product ..."
           << endl;
      cout << "    highest part of sum : " << sumhihi << endl;
      cout << "2nd highest part of sum : " << sumlohi << endl;
      cout << " 2nd lowest part of sum : " << sumhilo << endl;
      cout << "     lowest part of sum : " << sumlolo << endl;

      cout << "the error : " << err << endl;
   }
   return err;
}

double test_dbl4_complex_exponential ( int deg, int verbose )
{
   double* xrehihi = new double[deg+1];
   double* xrelohi = new double[deg+1];
   double* xrehilo = new double[deg+1];
   double* xrelolo = new double[deg+1];
   double* ximhihi = new double[deg+1];
   double* ximlohi = new double[deg+1];
   double* ximhilo = new double[deg+1];
   double* ximlolo = new double[deg+1];
   double* yrehihi = new double[deg+1];
   double* yrelohi = new double[deg+1];
   double* yrehilo = new double[deg+1];
   double* yrelolo = new double[deg+1];
   double* yimhihi = new double[deg+1];
   double* yimlohi = new double[deg+1];
   double* yimhilo = new double[deg+1];
   double* yimlolo = new double[deg+1];
   double* zrehihi_h = new double[deg+1];
   double* zrelohi_h = new double[deg+1];
   double* zrehilo_h = new double[deg+1];
   double* zrelolo_h = new double[deg+1];
   double* zimhihi_h = new double[deg+1];
   double* zimlohi_h = new double[deg+1];
   double* zimhilo_h = new double[deg+1];
   double* zimlolo_h = new double[deg+1];
   double* zrehihi_d = new double[deg+1];
   double* zrelohi_d = new double[deg+1];
   double* zrehilo_d = new double[deg+1];
   double* zrelolo_d = new double[deg+1];
   double* zimhihi_d = new double[deg+1];
   double* zimlohi_d = new double[deg+1];
   double* zimhilo_d = new double[deg+1];
   double* zimlolo_d = new double[deg+1];
   double rndrehihi,rndrelohi,rndrehilo,rndrelolo;
   double rndimhihi,rndimlohi,rndimhilo,rndimlolo;
   double sumrehihi,sumrelohi,sumrehilo,sumrelolo;
   double sumimhihi,sumimlohi,sumimhilo,sumimlolo;

   random_cmplx4_exponentials
      (deg,&rndrehihi,&rndrelohi,&rndrehilo,&rndrelolo,
           &rndimhihi,&rndimlohi,&rndimhilo,&rndimlolo,
           xrehihi,xrelohi,xrehilo,xrelolo,ximhihi,ximlohi,ximhilo,ximlolo,
           yrehihi,yrelohi,yrehilo,yrelolo,yimhihi,yimlohi,yimhilo,yimlolo);

   CPU_cmplx4_product
      (deg,xrehihi,xrelohi,xrehilo,xrelolo,ximhihi,ximlohi,ximhilo,ximlolo,
           yrehihi,yrelohi,yrehilo,yrelolo,yimhihi,yimlohi,yimhilo,yimlolo,
           zrehihi_h,zrelohi_h,zrehilo_h,zrelolo_h,
           zimhihi_h,zimlohi_h,zimhilo_h,zimlolo_h);

   if(verbose > 0)
   {
      cout << scientific << setprecision(16);

      cout << "Product of series of exp(x) with series of exp(-x)," << endl;
      cout << "  for xrehihi = " << rndrehihi;
      cout << "      xrelohi = " << rndrelohi << endl;
      cout << "      xrehilo = " << rndrehilo;
      cout << "  and xrelolo = " << rndrelolo << endl;
      cout << "  for ximhihi = " << rndimhihi;
      cout << "      ximlohi = " << rndimlohi << endl;
      cout << "      ximhilo = " << rndimhilo;
      cout << "  and ximlolo = " << rndimlolo << endl;

      sumrehihi = 0.0; sumrelohi = 0.0; sumrehilo = 0.0; sumrelolo = 0.0;
      sumimhihi = 0.0; sumimlohi = 0.0; sumimhilo = 0.0; sumimlolo = 0.0;

      for(int k=0; k<=deg; k++) 
      {
         qdf_inc(&sumrehihi,&sumrelohi,&sumrehilo,&sumrelolo,
                 zrehihi_h[k],zrelohi_h[k],zrehilo_h[k],zrelolo_h[k]);
         qdf_inc(&sumimhihi,&sumimlohi,&sumimhilo,&sumimlolo,
                 zimhihi_h[k],zimlohi_h[k],zimhilo_h[k],zimlolo_h[k]);
      }
      cout << "Summation of all coefficients of the product ..." << endl;
      cout << "  sumrehihi : " << sumrehihi;
      cout << "  sumrelohi : " << sumrelohi << endl;
      cout << "  sumrehilo : " << sumrehilo;
      cout << "  sumrelolo : " << sumrelolo << endl;
      cout << "  sumimhihi : " << sumimhihi;
      cout << "  sumimlohi : " << sumimlohi << endl;
      cout << "  sumimhilo : " << sumimhilo;
      cout << "  sumimlolo : " << sumimlolo << endl;
   }
   GPU_cmplx4_product
      (xrehihi,xrelohi,xrehilo,xrelolo,ximhihi,ximlohi,ximhilo,ximlolo,
       yrehihi,yrelohi,yrehilo,yrelolo,yimhihi,yimlohi,yimhilo,yimlolo,
       zrehihi_d,zrelohi_d,zrehilo_d,zrelolo_d,
       zimhihi_d,zimlohi_d,zimhilo_d,zimlolo_d,deg,1,deg+1,2);

   if(verbose > 0)
   {
      sumrehihi = 0.0; sumrelohi = 0.0; sumimhihi = 0.0; sumimlohi = 0.0;
      sumrehilo = 0.0; sumrelolo = 0.0; sumimhilo = 0.0; sumimlolo = 0.0;
   }
   double err = 0.0;

   for(int k=0; k<=deg; k++) 
   {
      if(verbose > 0)
      {
         qdf_inc(&sumrehihi,&sumrelohi,&sumrehilo,&sumrelolo,
                 zrehihi_d[k],zrelohi_d[k],zrehilo_d[k],zrelolo_d[k]);
         qdf_inc(&sumimhihi,&sumimlohi,&sumimhilo,&sumimlolo,
                 zimhihi_d[k],zimlohi_d[k],zimhilo_d[k],zimlolo_d[k]);
      }
      err = err
          + abs(zrehihi_h[k] - zrehihi_d[k])
          + abs(zrelohi_h[k] - zrelohi_d[k])
          + abs(zrehilo_h[k] - zrehilo_d[k])
          + abs(zrelolo_h[k] - zrelolo_d[k])
          + abs(zimhihi_h[k] - zimhihi_d[k])
          + abs(zimlohi_h[k] - zimlohi_d[k])
          + abs(zimhilo_h[k] - zimhilo_d[k])
          + abs(zimlolo_h[k] - zimlolo_d[k]);
   }
   if(verbose > 0)
   {
      cout << "Summation of all coefficients of the GPU computed product ..."
           << endl;
      cout << "  sumrehihi : " << sumrehihi;
      cout << "  sumrelohi : " << sumrelohi << endl;
      cout << "  sumrehilo : " << sumrehilo;
      cout << "  sumrelolo : " << sumrelolo << endl;
      cout << "  sumimhihi : " << sumimhihi;
      cout << "  sumimlohi : " << sumimlohi << endl;
      cout << "  sumimhilo : " << sumimhilo;
      cout << "  sumimlolo : " << sumimlolo << endl;

      cout << "the error : " << err << endl;
   }
   return err;
}
