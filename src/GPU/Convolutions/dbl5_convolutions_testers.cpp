/* The file dbl5_convolutions_testers.cpp contains the definitions of
 * functions to test the product of two series in penta double precision. */

#include <ctime>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector_types.h>
#include "penta_double_functions.h"
#include "random5_vectors.h"
#include "random5_series.h"
#include "dbl5_convolutions_host.h"
#include "dbl5_convolutions_kernels.h"
#include "dbl5_convolutions_testers.h"

using namespace std;

int main_dbl5_test ( int seed, int deg, int vrblvl )
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
      double realerror1 = test_dbl5_real(deg,vrblvl-1);
      double realerror2 = test_dbl5_real_random(deg,vrblvl-1);

      cout.flags(f);
      double complexerror1 = test_dbl5_complex(deg,vrblvl-1);
      double complexerror2 = test_dbl5_complex_random(deg,vrblvl-1);

      double realerror3 = test_dbl5_real_exponential(deg,vrblvl-1);
      double complexerror3 = test_dbl5_complex_exponential(deg,vrblvl-1);

      const double tol = 1.0e-75;

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

double test_dbl5_real ( int deg, int verbose )
{
   double* xtb = new double[deg+1];
   double* xix = new double[deg+1];
   double* xmi = new double[deg+1];
   double* xrg = new double[deg+1];
   double* xpk = new double[deg+1];
   double* ytb = new double[deg+1];
   double* yix = new double[deg+1];
   double* ymi = new double[deg+1];
   double* yrg = new double[deg+1];
   double* ypk = new double[deg+1];
   double* ztb_h = new double[deg+1];
   double* zix_h = new double[deg+1];
   double* zmi_h = new double[deg+1];
   double* zrg_h = new double[deg+1];
   double* zpk_h = new double[deg+1];
   double* ztb_d = new double[deg+1];
   double* zix_d = new double[deg+1];
   double* zmi_d = new double[deg+1];
   double* zrg_d = new double[deg+1];
   double* zpk_d = new double[deg+1];

   for(int k=0; k<=deg; k++)
   {
      xtb[k] = 1.0; xix[k] = 0.0; xmi[k] = 0.0; xrg[k] = 0.0; xpk[k] = 0.0;
      ytb[k] = 0.0; yix[k] = 0.0; ymi[k] = 0.0; yrg[k] = 0.0; ypk[k] = 0.0;
   }
   ytb[0] = 1.0; ytb[1] = -1.0;

   CPU_dbl5_product
      (deg,xtb,xix,xmi,xrg,xpk,ytb,yix,ymi,yrg,ypk,
       ztb_h,zix_h,zmi_h,zrg_h,zpk_h);

   if(verbose > 0)
   {
      cout << "Series of 1/(1-x) multiplied with 1-x : " << endl;

      for(int k=0; k<=deg; k++)
      {
         cout << "ztb[" << k << "] : " << ztb_h[k];
         cout << "  zix[" << k << "] : " << zix_h[k];
         cout << "  zmi[" << k << "] : " << zmi_h[k];
         cout << "  zrg[" << k << "] : " << zrg_h[k];
         cout << "  zpk[" << k << "] : " << zpk_h[k] << endl;
      }
      /*
         for(int k=0; k<=deg; k++)
         {
            cout << "x[" << k << "] : " << endl;
            pdf_write_doubles(xtb[k],xix[k],xmi[k],xrg[k],xpk[k]);
            cout << endl;
            cout << "y[" << k << "] : " << endl;
            pdf_write_doubles(ytb[k],yix[k],ymi[k],yrg[k],ypk[k]);
            cout << endl;
         }
       */
   }
   GPU_dbl5_product
      (xtb,xix,xmi,xrg,xpk,ytb,yix,ymi,yrg,ypk,
       ztb_d,zix_d,zmi_d,zrg_d,zpk_d,deg,1,deg+1,1);

   if(verbose > 0) cout << "GPU computed product :" << endl;

   double err = 0.0;

   for(int k=0; k<=deg; k++)
   {
      if(verbose > 0)
      {
         cout << "ztb[" << k << "] : " << ztb_d[k];
         cout << "  zix[" << k << "] : " << zix_d[k];
         cout << "  zmi[" << k << "] : " << zmi_d[k];
         cout << "  zrg[" << k << "] : " << zrg_d[k];
         cout << "  zpk[" << k << "] : " << zpk_d[k] << endl;
      }
      err = err
          + abs(ztb_h[k] - ztb_d[k]) + abs(zix_h[k] - zix_d[k])
          + abs(zmi_h[k] - zmi_d[k]) + abs(zrg_h[k] - zrg_d[k])
          + abs(zpk_h[k] - zpk_d[k]);
   }
   if(verbose > 0)
   {
      cout << scientific << setprecision(16);
      cout << "the error : " << err << endl;
   }
   return err;
}

double test_dbl5_real_random ( int deg, int verbose )
{
   double* xtb = new double[deg+1];
   double* xix = new double[deg+1];
   double* xmi = new double[deg+1];
   double* xrg = new double[deg+1];
   double* xpk = new double[deg+1];
   double* ytb = new double[deg+1];
   double* yix = new double[deg+1];
   double* ymi = new double[deg+1];
   double* yrg = new double[deg+1];
   double* ypk = new double[deg+1];
   double* ztb_h = new double[deg+1];
   double* zix_h = new double[deg+1];
   double* zmi_h = new double[deg+1];
   double* zrg_h = new double[deg+1];
   double* zpk_h = new double[deg+1];
   double* ztb_d = new double[deg+1];
   double* zix_d = new double[deg+1];
   double* zmi_d = new double[deg+1];
   double* zrg_d = new double[deg+1];
   double* zpk_d = new double[deg+1];

   for(int k=0; k<=deg; k++)
   {
      random_penta_double
         (&xtb[k],&xix[k],&xmi[k],&xrg[k],&xpk[k]);
      random_penta_double
         (&ytb[k],&yix[k],&ymi[k],&yrg[k],&ypk[k]);
   }
   CPU_dbl5_product
      (deg,xtb,xix,xmi,xrg,xpk,ytb,yix,ymi,yrg,ypk,
           ztb_h,zix_h,zmi_h,zrg_h,zpk_h);

   if(verbose > 0)
   {
      cout << "Product of two random real series : " << endl;
      cout << scientific << setprecision(16);

      for(int k=0; k<=deg; k++)
      {
         cout << "ztb[" << k << "] : " << ztb_h[k];
         cout << "  zix[" << k << "] : " << zix_h[k] << endl;
         cout << "zmi[" << k << "] : " << zmi_h[k];
         cout << "  zrg[" << k << "] : " << zrg_h[k] << endl;
         cout << "zpk[" << k << "] : " << zpk_h[k] << endl;
      }
   }
   GPU_dbl5_product
      (xtb,xix,xmi,xrg,xpk,ytb,yix,ymi,yrg,ypk,
       ztb_d,zix_d,zmi_d,zrg_d,zpk_d,deg,1,deg+1,1);

   if(verbose > 0) cout << "GPU computed product :" << endl;

   double err = 0.0;

   for(int k=0; k<=deg; k++)
   {
      if(verbose > 0)
      {
         cout << "ztb[" << k << "] : " << ztb_d[k];
         cout << "  zix[" << k << "] : " << zix_d[k] << endl;
         cout << "zmi[" << k << "] : " << zmi_d[k];
         cout << "  zrg[" << k << "] : " << zrg_d[k] << endl;
         cout << "zpk[" << k << "] : " << zpk_d[k] << endl;
      }
      err = err
          + abs(ztb_h[k] - ztb_d[k])
          + abs(zix_h[k] - zix_d[k])
          + abs(zmi_h[k] - zmi_d[k])
          + abs(zrg_h[k] - zrg_d[k])
          + abs(zpk_h[k] - zpk_d[k]);
   }
   if(verbose > 0)
   {
      cout << scientific << setprecision(16);
      cout << "the error : " << err << endl;
   }
   return err;
}

double test_dbl5_complex ( int deg, int verbose )
{
   double* xretb = new double[deg+1];
   double* xreix = new double[deg+1];
   double* xremi = new double[deg+1];
   double* xrerg = new double[deg+1];
   double* xrepk = new double[deg+1];
   double* ximtb = new double[deg+1];
   double* ximix = new double[deg+1];
   double* ximmi = new double[deg+1];
   double* ximrg = new double[deg+1];
   double* ximpk = new double[deg+1];
   double* yretb = new double[deg+1];
   double* yreix = new double[deg+1];
   double* yremi = new double[deg+1];
   double* yrerg = new double[deg+1];
   double* yrepk = new double[deg+1];
   double* yimtb = new double[deg+1];
   double* yimix = new double[deg+1];
   double* yimmi = new double[deg+1];
   double* yimrg = new double[deg+1];
   double* yimpk = new double[deg+1];
   double* zretb_h = new double[deg+1];
   double* zreix_h = new double[deg+1];
   double* zremi_h = new double[deg+1];
   double* zrerg_h = new double[deg+1];
   double* zrepk_h = new double[deg+1];
   double* zimtb_h = new double[deg+1];
   double* zimix_h = new double[deg+1];
   double* zimmi_h = new double[deg+1];
   double* zimrg_h = new double[deg+1];
   double* zimpk_h = new double[deg+1];
   double* zretb_d = new double[deg+1];
   double* zreix_d = new double[deg+1];
   double* zremi_d = new double[deg+1];
   double* zrerg_d = new double[deg+1];
   double* zrepk_d = new double[deg+1];
   double* zimtb_d = new double[deg+1];
   double* zimix_d = new double[deg+1];
   double* zimmi_d = new double[deg+1];
   double* zimrg_d = new double[deg+1];
   double* zimpk_d = new double[deg+1];

   for(int k=0; k<=deg; k++)
   {
      xretb[k] = 1.0; xreix[k] = 0.0; xremi[k] = 0.0;
      xrerg[k] = 0.0; xrepk[k] = 0.0;
      ximtb[k] = 0.0; ximix[k] = 0.0; ximmi[k] = 0.0;
      ximrg[k] = 0.0; ximpk[k] = 0.0;
      yretb[k] = 0.0; yreix[k] = 0.0; yremi[k] = 0.0;
      yrerg[k] = 0.0; yrepk[k] = 0.0;
      yimtb[k] = 0.0; yimix[k] = 0.0; yimmi[k] = 0.0;
      yimrg[k] = 0.0; yimpk[k] = 0.0;
   }
   yretb[0] = 1.0; yretb[1] = -1.0;

   CPU_cmplx5_product
      (deg,xretb,xreix,xremi,xrerg,xrepk,ximtb,ximix,ximmi,ximrg,ximpk,
           yretb,yreix,yremi,yrerg,yrepk,yimtb,yimix,yimmi,yimrg,yimpk,
           zretb_h,zreix_h,zremi_h,zrerg_h,zrepk_h,
           zimtb_h,zimix_h,zimmi_h,zimrg_h,zimpk_h);

   if(verbose > 0)
   {
      cout << "Series of 1/(1-x) multiplied with 1-x :" << endl;

      for(int k=0; k<=deg; k++)
      {
         cout << "zretb[" << k << "] : " << zretb_h[k];
         cout << "  zreix[" << k << "] : " << zreix_h[k];
         cout << "  zremi[" << k << "] : " << zremi_h[k];
         cout << "  zrerg[" << k << "] : " << zrerg_h[k];
         cout << "  zrepk[" << k << "] : " << zrepk_h[k] << endl;
         cout << "zimtb[" << k << "] : " << zimtb_h[k];
         cout << "  zimix[" << k << "] : " << zimix_h[k];
         cout << "  zimmi[" << k << "] : " << zimmi_h[k];
         cout << "  zimrg[" << k << "] : " << zimrg_h[k];
         cout << "  zimpk[" << k << "] : " << zimpk_h[k] << endl;
      }
   }
   GPU_cmplx5_product
      (xretb,xreix,xremi,xrerg,xrepk,ximtb,ximix,ximmi,ximrg,ximpk,
       yretb,yreix,yremi,yrerg,yrepk,yimtb,yimix,yimmi,yimrg,yimpk,
       zretb_d,zreix_d,zremi_d,zrerg_d,zrepk_d,
       zimtb_d,zimix_d,zimmi_d,zimrg_d,zimpk_d,deg,1,deg+1,2);

   if(verbose > 0) cout << "GPU computed product :" << endl;

   double err = 0.0;

   for(int k=0; k<=deg; k++)
   {
      if(verbose > 0)
      {
         cout << "zretb[" << k << "] : " << zretb_d[k];
         cout << "  zreix[" << k << "] : " << zreix_d[k];
         cout << "  zremi[" << k << "] : " << zremi_d[k];
         cout << "  zrerg[" << k << "] : " << zrerg_d[k];
         cout << "  zrepk[" << k << "] : " << zrepk_d[k] << endl;
         cout << "zimtb[" << k << "] : " << zimtb_d[k];
         cout << "  zimix[" << k << "] : " << zimix_d[k];
         cout << "  zimmi[" << k << "] : " << zimmi_d[k];
         cout << "  zimrg[" << k << "] : " << zimrg_d[k];
         cout << "  zimpk[" << k << "] : " << zimpk_d[k] << endl;
      }
      err = err
          + abs(zretb_h[k] - zretb_d[k]) + abs(zreix_h[k] - zreix_d[k])
          + abs(zremi_h[k] - zremi_d[k]) + abs(zrerg_h[k] - zrerg_d[k])
          + abs(zrepk_h[k] - zrepk_d[k])
          + abs(zimtb_h[k] - zimtb_d[k]) + abs(zimix_h[k] - zimix_d[k])
          + abs(zimmi_h[k] - zimmi_d[k]) + abs(zimrg_h[k] - zimrg_d[k])
          + abs(zimpk_h[k] - zimpk_d[k]);
   }
   if(verbose > 0)
   {
      cout << scientific << setprecision(16);
      cout << "the error : " << err << endl;
   }
   return err;
}

double test_dbl5_complex_random ( int deg, int verbose )
{
   double* xretb = new double[deg+1];
   double* xreix = new double[deg+1];
   double* xremi = new double[deg+1];
   double* xrerg = new double[deg+1];
   double* xrepk = new double[deg+1];
   double* ximtb = new double[deg+1];
   double* ximix = new double[deg+1];
   double* ximmi = new double[deg+1];
   double* ximrg = new double[deg+1];
   double* ximpk = new double[deg+1];
   double* yretb = new double[deg+1];
   double* yreix = new double[deg+1];
   double* yremi = new double[deg+1];
   double* yrerg = new double[deg+1];
   double* yrepk = new double[deg+1];
   double* yimtb = new double[deg+1];
   double* yimix = new double[deg+1];
   double* yimmi = new double[deg+1];
   double* yimrg = new double[deg+1];
   double* yimpk = new double[deg+1];
   double* zretb_h = new double[deg+1];
   double* zreix_h = new double[deg+1];
   double* zremi_h = new double[deg+1];
   double* zrerg_h = new double[deg+1];
   double* zrepk_h = new double[deg+1];
   double* zimtb_h = new double[deg+1];
   double* zimix_h = new double[deg+1];
   double* zimmi_h = new double[deg+1];
   double* zimrg_h = new double[deg+1];
   double* zimpk_h = new double[deg+1];
   double* zretb_d = new double[deg+1];
   double* zreix_d = new double[deg+1];
   double* zremi_d = new double[deg+1];
   double* zrerg_d = new double[deg+1];
   double* zrepk_d = new double[deg+1];
   double* zimtb_d = new double[deg+1];
   double* zimix_d = new double[deg+1];
   double* zimmi_d = new double[deg+1];
   double* zimrg_d = new double[deg+1];
   double* zimpk_d = new double[deg+1];

   for(int k=0; k<=deg; k++)
   {
      random_penta_double
         (&xretb[k],&xreix[k],&xremi[k],&xrerg[k],&xrepk[k]);
      random_penta_double
         (&ximtb[k],&ximix[k],&ximmi[k],&ximrg[k],&ximpk[k]);
      random_penta_double
         (&yretb[k],&yreix[k],&yremi[k],&yrerg[k],&yrepk[k]);
      random_penta_double
         (&yimtb[k],&yimix[k],&yimmi[k],&yimrg[k],&yimpk[k]);
   }
   CPU_cmplx5_product
      (deg,xretb,xreix,xremi,xrerg,xrepk,ximtb,ximix,ximmi,ximrg,ximpk,
           yretb,yreix,yremi,yrerg,yrepk,yimtb,yimix,yimmi,yimrg,yimpk,
           zretb_h,zreix_h,zremi_h,zrerg_h,zrepk_h,
           zimtb_h,zimix_h,zimmi_h,zimrg_h,zimpk_h);

   if(verbose > 0)
   {
      cout << "Product of two random series :" << endl;

      for(int k=0; k<=deg; k++)
      {
         cout << "zretb[" << k << "] : " << zretb_h[k];
         cout << "  zreix[" << k << "] : " << zreix_h[k] << endl;
         cout << "zremi[" << k << "] : " << zremi_h[k];
         cout << "  zrerg[" << k << "] : " << zrerg_h[k] << endl;
         cout << "zrepk[" << k << "] : " << zrepk_h[k] << endl;
         cout << "zimtb[" << k << "] : " << zimtb_h[k];
         cout << "  zimix[" << k << "] : " << zimix_h[k] << endl;
         cout << "zimmi[" << k << "] : " << zimmi_h[k];
         cout << "  zimrg[" << k << "] : " << zimrg_h[k] << endl;
         cout << "zimpk[" << k << "] : " << zimpk_h[k] << endl;
      }
   }
   GPU_cmplx5_product
      (xretb,xreix,xremi,xrerg,xrepk,ximtb,ximix,ximmi,ximrg,ximpk,
       yretb,yreix,yremi,yrerg,yrepk,yimtb,yimix,yimmi,yimrg,yimpk,
       zretb_d,zreix_d,zremi_d,zrerg_d,zrepk_d,
       zimtb_d,zimix_d,zimmi_d,zimrg_d,zimpk_d,deg,1,deg+1,2);

   if(verbose > 0) cout << "GPU computed product :" << endl;

   double err = 0.0;

   for(int k=0; k<=deg; k++)
   {
      if(verbose > 0)
      {
         cout << "zretb[" << k << "] : " << zretb_d[k];
         cout << "  zreix[" << k << "] : " << zreix_d[k] << endl;
         cout << "zremi[" << k << "] : " << zremi_d[k];
         cout << "  zrerg[" << k << "] : " << zrerg_d[k] << endl;
         cout << "zrepk[" << k << "] : " << zrepk_d[k] << endl;
         cout << "zimtb[" << k << "] : " << zimtb_d[k];
         cout << "  zimix[" << k << "] : " << zimix_d[k] << endl;
         cout << "zimmi[" << k << "] : " << zimmi_d[k];
         cout << "  zimrg[" << k << "] : " << zimrg_d[k] << endl;
         cout << "zimpk[" << k << "] : " << zimpk_d[k] << endl;
      }
      err = err
          + abs(zretb_h[k] - zretb_d[k]) + abs(zreix_h[k] - zreix_d[k])
          + abs(zremi_h[k] - zremi_d[k]) + abs(zrerg_h[k] - zrerg_d[k])
          + abs(zrepk_h[k] - zrepk_d[k])
          + abs(zimtb_h[k] - zimtb_d[k]) + abs(zimix_h[k] - zimix_d[k])
          + abs(zimmi_h[k] - zimmi_d[k]) + abs(zimrg_h[k] - zimrg_d[k])
          + abs(zimpk_h[k] - zimpk_d[k]);
   }
   if(verbose > 0)
   {
      cout << scientific << setprecision(16);
      cout << "the error : " << err << endl;
   }
   return err;
}

double test_dbl5_real_exponential ( int deg, int verbose )
{
   double *xtb = new double[deg+1];
   double *xix = new double[deg+1];
   double *xmi = new double[deg+1];
   double *xrg = new double[deg+1];
   double *xpk = new double[deg+1];
   double *ytb = new double[deg+1];
   double *yix = new double[deg+1];
   double *ymi = new double[deg+1];
   double *yrg = new double[deg+1];
   double *ypk = new double[deg+1];
   double *ztb_h = new double[deg+1];
   double *zix_h = new double[deg+1];
   double *zmi_h = new double[deg+1];
   double *zrg_h = new double[deg+1];
   double *zpk_h = new double[deg+1];
   double *ztb_d = new double[deg+1];
   double *zix_d = new double[deg+1];
   double *zmi_d = new double[deg+1];
   double *zrg_d = new double[deg+1];
   double *zpk_d = new double[deg+1];
   double rtb,rix,rmi,rrg,rpk;
   double sumtb,sumix,summi,sumrg,sumpk;

   random_dbl5_exponentials
      (deg,&rtb,&rix,&rmi,&rrg,&rpk,
           xtb,xix,xmi,xrg,xpk,ytb,yix,ymi,yrg,ypk);

   CPU_dbl5_product
      (deg,xtb,xix,xmi,xrg,xpk,ytb,yix,ymi,yrg,ypk,
       ztb_h,zix_h,zmi_h,zrg_h,zpk_h);

   if(verbose > 0)
   {
      cout << scientific << setprecision(16);

      cout << "Product of series of exp(x) with series of exp(-x)," << endl;
      cout << "  for xtb = " << rtb << endl;
      cout << "      xix = " << rix << endl;
      cout << "      xmi = " << rmi << endl;
      cout << "      xrg = " << rrg << endl;
      cout << "  and xpk = " << rpk << endl;

      sumtb = 0.0; sumix = 0.0; summi = 0.0; sumrg = 0.0; sumpk = 0.0;

      for(int k=0; k<=deg; k++)
         pdf_inc(&sumtb,&sumix,&summi,&sumrg,&sumpk,
                 ztb_h[k],zix_h[k],zmi_h[k],zrg_h[k],zpk_h[k]);

      cout << "Summation of all coefficients in the product ..." << endl;
      cout << "       highest part of the sum : " << sumtb << endl;
      cout << "second highest part of the sum : " << sumix << endl;
      cout << "        middle part of the sum : " << summi << endl;
      cout << " second lowest part of the sum : " << sumrg << endl;
      cout << "        lowest part of the sum : " << sumpk << endl;
   }
   GPU_dbl5_product
      (xtb,xix,xmi,xrg,xpk,ytb,yix,ymi,yrg,ypk,
       ztb_d,zix_d,zmi_d,zrg_d,zpk_d,deg,1,deg+1,1);

   if(verbose > 0)
   {
      sumtb = 0.0; sumix = 0.0; summi = 0.0; sumrg = 0.0; sumpk = 0.0;
   }
   double err = 0.0;

   for(int k=0; k<=deg; k++)
   {
      if(verbose > 0)
         pdf_inc(&sumtb,&sumix,&summi,&sumrg,&sumpk,
                 ztb_d[k],zix_d[k],zmi_d[k],zrg_d[k],zpk_d[k]);

      err = err
          + abs(ztb_h[k] - ztb_d[k]) + abs(zix_h[k] - zix_d[k])
          + abs(zmi_h[k] - zmi_d[k]) + abs(zrg_h[k] - zrg_d[k])
          + abs(zpk_h[k] - zpk_d[k]);
   }
   if(verbose > 0)
   {
      cout << "Summation of all coefficients in the GPU computed product ..."
           << endl;
      cout << "       highest part of the sum : " << sumtb << endl;
      cout << "second highest part of the sum : " << sumix << endl;
      cout << "        middle part of the sum : " << summi << endl;
      cout << " second lowest part of the sum : " << sumrg << endl;
      cout << "        lowest part of the sum : " << sumpk << endl;

      cout << "the error : " << err << endl;
   }
   return err;
}

double test_dbl5_complex_exponential ( int deg, int verbose )
{
   double* xretb = new double[deg+1];
   double* xreix = new double[deg+1];
   double* xremi = new double[deg+1];
   double* xrerg = new double[deg+1];
   double* xrepk = new double[deg+1];
   double* ximtb = new double[deg+1];
   double* ximix = new double[deg+1];
   double* ximmi = new double[deg+1];
   double* ximrg = new double[deg+1];
   double* ximpk = new double[deg+1];
   double* yretb = new double[deg+1];
   double* yreix = new double[deg+1];
   double* yremi = new double[deg+1];
   double* yrerg = new double[deg+1];
   double* yrepk = new double[deg+1];
   double* yimtb = new double[deg+1];
   double* yimix = new double[deg+1];
   double* yimmi = new double[deg+1];
   double* yimrg = new double[deg+1];
   double* yimpk = new double[deg+1];
   double* zretb_h = new double[deg+1];
   double* zreix_h = new double[deg+1];
   double* zremi_h = new double[deg+1];
   double* zrerg_h = new double[deg+1];
   double* zrepk_h = new double[deg+1];
   double* zimtb_h = new double[deg+1];
   double* zimix_h = new double[deg+1];
   double* zimmi_h = new double[deg+1];
   double* zimrg_h = new double[deg+1];
   double* zimpk_h = new double[deg+1];
   double* zretb_d = new double[deg+1];
   double* zreix_d = new double[deg+1];
   double* zremi_d = new double[deg+1];
   double* zrerg_d = new double[deg+1];
   double* zrepk_d = new double[deg+1];
   double* zimtb_d = new double[deg+1];
   double* zimix_d = new double[deg+1];
   double* zimmi_d = new double[deg+1];
   double* zimrg_d = new double[deg+1];
   double* zimpk_d = new double[deg+1];
   double rndretb,rndreix,rndremi,rndrerg,rndrepk;
   double rndimtb,rndimix,rndimmi,rndimrg,rndimpk;
   double sumretb,sumreix,sumremi,sumrerg,sumrepk;
   double sumimtb,sumimix,sumimmi,sumimrg,sumimpk;

   random_cmplx5_exponentials
      (deg,&rndretb,&rndreix,&rndremi,&rndrerg,&rndrepk,
           &rndimtb,&rndimix,&rndimmi,&rndimrg,&rndimpk,
           xretb,xreix,xremi,xrerg,xrepk,ximtb,ximix,ximmi,ximrg,ximpk,
           yretb,yreix,yremi,yrerg,yrepk,yimtb,yimix,yimmi,yimrg,yimpk);

   CPU_cmplx5_product
      (deg,xretb,xreix,xremi,xrerg,xrepk,ximtb,ximix,ximmi,ximrg,ximpk,
           yretb,yreix,yremi,yrerg,yrepk,yimtb,yimix,yimmi,yimrg,yimpk,
           zretb_h,zreix_h,zremi_h,zrerg_h,zrepk_h,
           zimtb_h,zimix_h,zimmi_h,zimrg_h,zimpk_h);

   if(verbose > 0)
   {
      /*
         for(int k=0; k<=deg; k++)
         {
            cout << "zre[" << k << "] : " << endl;
            pdf_write_doubles
               (zretb_h[k],zreix_h[k],zremi_h[k],zrerg_h[k],zrepk_h[k]);
            cout << endl;
            cout << "zim[" << k << "] : " << endl;
            pdf_write_doubles
               (zimtb_h[k],zimix_h[k],zimmi_h[k],zimrg_h[k],zimpk_h[k]);
            cout << endl;
         }
      */
      cout << scientific << setprecision(16);

      cout << "Product of series of exp(x) with series of exp(-x)," << endl;
      cout << "  for xretb = " << rndretb;
      cout << "      xreix = " << rndreix << endl;
      cout << "      xremi = " << rndremi;
      cout << "      xrerg = " << rndrerg << endl;
      cout << "  and xrepk = " << rndrepk << endl;
      cout << "  for ximtb = " << rndimtb;
      cout << "      ximix = " << rndimix << endl;
      cout << "      ximmi = " << rndimmi;
      cout << "      ximrg = " << rndimrg << endl;
      cout << "  and ximpk = " << rndimpk << endl;

      sumretb = 0.0; sumreix = 0.0; sumremi = 0.0;
      sumrerg = 0.0; sumrepk = 0.0;
      sumimtb = 0.0; sumimix = 0.0; sumimmi = 0.0;
      sumimrg = 0.0; sumimpk = 0.0;

      for(int k=0; k<=deg; k++) 
      {
         pdf_inc(&sumretb,&sumreix,&sumremi,&sumrerg,&sumrepk,
                 zretb_h[k],zreix_h[k],zremi_h[k],zrerg_h[k],zrepk_h[k]);
         pdf_inc(&sumimtb,&sumimix,&sumimmi,&sumimrg,&sumimpk,
                 zimtb_h[k],zimix_h[k],zimmi_h[k],zimrg_h[k],zimpk_h[k]);
      }
      cout << "Summation of all coefficients of the product ..." << endl;
      cout << "  sumretb : " << sumretb;
      cout << "  sumreix : " << sumreix << endl;
      cout << "  sumremi : " << sumremi;
      cout << "  sumrerg : " << sumrerg << endl;
      cout << "  sumrepk : " << sumrepk << endl;
      cout << "  sumimtb : " << sumimtb;
      cout << "  sumimix : " << sumimix << endl;
      cout << "  sumimmi : " << sumimmi;
      cout << "  sumimrg : " << sumimrg << endl;
      cout << "  sumimpk : " << sumimpk << endl;

      for(int k=0; k<=deg; k++)
      {
         cout << "xre[" << k << "] : " << endl;
         pdf_write_doubles(xretb[k],xreix[k],xremi[k],xrerg[k],xrepk[k]);
         cout << endl;
         cout << "xim[" << k << "] : " << endl;
         pdf_write_doubles(ximtb[k],ximix[k],ximmi[k],ximrg[k],ximpk[k]);
         cout << endl;
         cout << "yre[" << k << "] : " << endl;
         pdf_write_doubles(yretb[k],yreix[k],yremi[k],yrerg[k],yrepk[k]);
         cout << endl;
         cout << "yim[" << k << "] : " << endl;
         pdf_write_doubles(yimtb[k],yimix[k],yimmi[k],yimrg[k],yimpk[k]);
         cout << endl;
      }
   }
   GPU_cmplx5_product
      (xretb,xreix,xremi,xrerg,xrepk,ximtb,ximix,ximmi,ximrg,ximpk,
       yretb,yreix,yremi,yrerg,yrepk,yimtb,yimix,yimmi,yimrg,yimpk,
       zretb_d,zreix_d,zremi_d,zrerg_d,zrepk_d,
       zimtb_d,zimix_d,zimmi_d,zimrg_d,zimpk_d,deg,1,deg+1,2);

   /*
      for(int k=0; k<=deg; k++)
      {
         cout << "zre[" << k << "] : " << endl;
         pdf_write_doubles
            (zretb_d[k],zreix_d[k],zremi_d[k],zrerg_d[k],zrepk_d[k]);
         cout << endl;
         cout << "zim[" << k << "] : " << endl;
         pdf_write_doubles
            (zimtb_d[k],zimix_d[k],zimmi_d[k],zimrg_d[k],zimpk_d[k]);
         cout << endl;
   }
   */
   if(verbose > 0)
   {
      sumretb = 0.0; sumreix = 0.0; sumremi = 0.0;
      sumrerg = 0.0; sumrepk = 0.0;
      sumimtb = 0.0; sumimix = 0.0; sumimmi = 0.0;
      sumimrg = 0.0; sumimpk = 0.0;
   }
   double err = 0.0;

   for(int k=0; k<=deg; k++) 
   {
      if(verbose > 0)
      {
         pdf_inc(&sumretb,&sumreix,&sumremi,&sumrerg,&sumrepk,
                 zretb_d[k],zreix_d[k],zremi_d[k],zrerg_d[k],zrepk_d[k]);
         pdf_inc(&sumimtb,&sumimix,&sumimmi,&sumimrg,&sumimpk,
                 zimtb_d[k],zimix_d[k],zimmi_d[k],zimrg_d[k],zimpk_d[k]);
      }
      err = err
          + abs(zretb_h[k] - zretb_d[k]) + abs(zreix_h[k] - zreix_d[k])
          + abs(zremi_h[k] - zremi_d[k]) + abs(zrerg_h[k] - zrerg_d[k])
          + abs(zrepk_h[k] - zrepk_d[k])
          + abs(zimtb_h[k] - zimtb_d[k]) + abs(zimix_h[k] - zimix_d[k])
          + abs(zimmi_h[k] - zimmi_d[k]) + abs(zimrg_h[k] - zimrg_d[k])
          + abs(zimpk_h[k] - zimpk_d[k]);
   }
   if(verbose > 0)
   {
      cout << "Summation of all coefficients of the GPU computed product ..."
           << endl;
      cout << "  sumretb : " << sumretb;
      cout << "  sumreix : " << sumreix << endl;
      cout << "  sumremi : " << sumremi;
      cout << "  sumrerg : " << sumrerg << endl;
      cout << "  sumrepk : " << sumrepk << endl;
      cout << "  sumimtb : " << sumimtb;
      cout << "  sumimix : " << sumimix << endl;
      cout << "  sumimmi : " << sumimmi;
      cout << "  sumimrg : " << sumimrg << endl;
      cout << "  sumimpk : " << sumimpk << endl;

      cout << "the error : " << err << endl;
   }
   return err;
}
