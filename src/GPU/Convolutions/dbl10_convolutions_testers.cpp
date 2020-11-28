/* The file dbl10_convolutions_testers.cpp contains the definitions of
 * functions to test the product of two series in deca double precision. */

#include <ctime>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector_types.h>
#include "deca_double_functions.h"
#include "random10_vectors.h"
#include "random10_series.h"
#include "dbl10_convolutions_host.h"
#include "dbl10_convolutions_kernels.h"
#include "dbl10_convolutions_testers.h"

using namespace std;

int main_dbl10_test ( int seed, int deg, int vrblvl )
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
      double realerror1 = test_dbl10_real(deg,vrblvl-1);
      double realerror2 = test_dbl10_real_random(deg,vrblvl-1);

      cout.flags(f);
      double complexerror1 = test_dbl10_complex(deg,vrblvl-1);
      double complexerror2 = test_dbl10_complex_random(deg,vrblvl-1);

      double realerror3 = test_dbl10_real_exponential(deg,vrblvl-1);
      double complexerror3 = test_dbl10_complex_exponential(deg,vrblvl-1);

      const double tol = 1.0e-154;

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

double test_dbl10_real ( int deg, int verbose )
{
   double* xrtb = new double[deg+1];
   double* xrix = new double[deg+1];
   double* xrmi = new double[deg+1];
   double* xrrg = new double[deg+1];
   double* xrpk = new double[deg+1];
   double* xltb = new double[deg+1];
   double* xlix = new double[deg+1];
   double* xlmi = new double[deg+1];
   double* xlrg = new double[deg+1];
   double* xlpk = new double[deg+1];
   double* yrtb = new double[deg+1];
   double* yrix = new double[deg+1];
   double* yrmi = new double[deg+1];
   double* yrrg = new double[deg+1];
   double* yrpk = new double[deg+1];
   double* yltb = new double[deg+1];
   double* ylix = new double[deg+1];
   double* ylmi = new double[deg+1];
   double* ylrg = new double[deg+1];
   double* ylpk = new double[deg+1];
   double* zrtb_h = new double[deg+1];
   double* zrix_h = new double[deg+1];
   double* zrmi_h = new double[deg+1];
   double* zrrg_h = new double[deg+1];
   double* zrpk_h = new double[deg+1];
   double* zltb_h = new double[deg+1];
   double* zlix_h = new double[deg+1];
   double* zlmi_h = new double[deg+1];
   double* zlrg_h = new double[deg+1];
   double* zlpk_h = new double[deg+1];
   double* zrtb_d = new double[deg+1];
   double* zrix_d = new double[deg+1];
   double* zrmi_d = new double[deg+1];
   double* zrrg_d = new double[deg+1];
   double* zrpk_d = new double[deg+1];
   double* zltb_d = new double[deg+1];
   double* zlix_d = new double[deg+1];
   double* zlmi_d = new double[deg+1];
   double* zlrg_d = new double[deg+1];
   double* zlpk_d = new double[deg+1];

   for(int k=0; k<=deg; k++)
   {
      xrtb[k] = 1.0; xrix[k] = 0.0; xrmi[k] = 0.0; xrrg[k] = 0.0;
      xrpk[k] = 0.0; xltb[k] = 0.0; xlix[k] = 0.0; xlmi[k] = 0.0;
      xlrg[k] = 0.0; xlpk[k] = 0.0;
      yrtb[k] = 0.0; yrix[k] = 0.0; yrmi[k] = 0.0; yrrg[k] = 0.0;
      yrpk[k] = 0.0; yltb[k] = 0.0; ylix[k] = 0.0; ylmi[k] = 0.0;
      ylrg[k] = 0.0; ylpk[k] = 0.0;
   }
   yrtb[0] = 1.0; yrtb[1] = -1.0;

   CPU_dbl10_product
      (deg,xrtb,xrix,xrmi,xrrg,xrpk,xltb,xlix,xlmi,xlrg,xlpk,
           yrtb,yrix,yrmi,yrrg,yrpk,yltb,ylix,ylmi,ylrg,ylpk,
       zrtb_h,zrix_h,zrmi_h,zrrg_h,zrpk_h,
       zltb_h,zlix_h,zlmi_h,zlrg_h,zlpk_h);

   if(verbose > 0)
   {
      cout << "Series of 1/(1-x) multiplied with 1-x : " << endl;

      for(int k=0; k<=deg; k++)
      {
         cout << "zrtb[" << k << "] : " << zrtb_h[k];
         cout << "  zrix[" << k << "] : " << zrix_h[k];
         cout << "  zrmi[" << k << "] : " << zrmi_h[k];
         cout << "  zrrg[" << k << "] : " << zrrg_h[k];
         cout << "  zrpk[" << k << "] : " << zrpk_h[k] << endl;
         cout << "zltb[" << k << "] : " << zltb_h[k];
         cout << "  zlix[" << k << "] : " << zlix_h[k];
         cout << "  zlmi[" << k << "] : " << zlmi_h[k];
         cout << "  zlrg[" << k << "] : " << zlrg_h[k];
         cout << "  zlpk[" << k << "] : " << zlpk_h[k] << endl;
      }
   }
   GPU_dbl10_product
      (xrtb,xrix,xrmi,xrrg,xrpk,xltb,xlix,xlmi,xlrg,xlpk,
       yrtb,yrix,yrmi,yrrg,yrpk,yltb,ylix,ylmi,ylrg,ylpk,
       zrtb_d,zrix_d,zrmi_d,zrrg_d,zrpk_d,
       zltb_d,zlix_d,zlmi_d,zlrg_d,zlpk_d,deg,1,deg+1,1);

   if(verbose > 0) cout << "GPU computed product :" << endl;

   double err = 0.0;

   for(int k=0; k<=deg; k++)
   {
      if(verbose > 0)
      {
         cout << "zrtb[" << k << "] : " << zrtb_d[k];
         cout << "  zrix[" << k << "] : " << zrix_d[k];
         cout << "  zrmi[" << k << "] : " << zrmi_d[k];
         cout << "  zrrg[" << k << "] : " << zrrg_d[k];
         cout << "  zrpk[" << k << "] : " << zrpk_d[k] << endl;
         cout << "zltb[" << k << "] : " << zltb_d[k];
         cout << "  zlix[" << k << "] : " << zlix_d[k];
         cout << "  zlmi[" << k << "] : " << zlmi_d[k];
         cout << "  zlrg[" << k << "] : " << zlrg_d[k];
         cout << "  zlpk[" << k << "] : " << zlpk_d[k] << endl;
      }
      err = err
          + abs(zrtb_h[k] - zrtb_d[k]) + abs(zrix_h[k] - zrix_d[k])
          + abs(zrmi_h[k] - zrmi_d[k]) + abs(zrrg_h[k] - zrrg_d[k])
          + abs(zrpk_h[k] - zrpk_d[k])
          + abs(zltb_h[k] - zltb_d[k]) + abs(zlix_h[k] - zlix_d[k])
          + abs(zlmi_h[k] - zlmi_d[k]) + abs(zlrg_h[k] - zlrg_d[k])
          + abs(zlpk_h[k] - zlpk_d[k]);
   }
   if(verbose > 0)
   {
      cout << endl;

      cout << scientific << setprecision(16);
      cout << "the error : " << err << endl;
   }
   return err;
}

double test_dbl10_real_random ( int deg, int verbose )
{
   double* xrtb = new double[deg+1];
   double* xrix = new double[deg+1];
   double* xrmi = new double[deg+1];
   double* xrrg = new double[deg+1];
   double* xrpk = new double[deg+1];
   double* xltb = new double[deg+1];
   double* xlix = new double[deg+1];
   double* xlmi = new double[deg+1];
   double* xlrg = new double[deg+1];
   double* xlpk = new double[deg+1];
   double* yrtb = new double[deg+1];
   double* yrix = new double[deg+1];
   double* yrmi = new double[deg+1];
   double* yrrg = new double[deg+1];
   double* yrpk = new double[deg+1];
   double* yltb = new double[deg+1];
   double* ylix = new double[deg+1];
   double* ylmi = new double[deg+1];
   double* ylrg = new double[deg+1];
   double* ylpk = new double[deg+1];
   double* zrtb_h = new double[deg+1];
   double* zrix_h = new double[deg+1];
   double* zrmi_h = new double[deg+1];
   double* zrrg_h = new double[deg+1];
   double* zrpk_h = new double[deg+1];
   double* zltb_h = new double[deg+1];
   double* zlix_h = new double[deg+1];
   double* zlmi_h = new double[deg+1];
   double* zlrg_h = new double[deg+1];
   double* zlpk_h = new double[deg+1];
   double* zrtb_d = new double[deg+1];
   double* zrix_d = new double[deg+1];
   double* zrmi_d = new double[deg+1];
   double* zrrg_d = new double[deg+1];
   double* zrpk_d = new double[deg+1];
   double* zltb_d = new double[deg+1];
   double* zlix_d = new double[deg+1];
   double* zlmi_d = new double[deg+1];
   double* zlrg_d = new double[deg+1];
   double* zlpk_d = new double[deg+1];

   for(int k=0; k<=deg; k++)
   {
      random_deca_double
         (&xrtb[k],&xrix[k],&xrmi[k],&xrrg[k],&xrpk[k],
          &xltb[k],&xlix[k],&xlmi[k],&xlrg[k],&xlpk[k]);
      random_deca_double
         (&yrtb[k],&yrix[k],&yrmi[k],&yrrg[k],&yrpk[k],
          &yltb[k],&ylix[k],&ylmi[k],&ylrg[k],&ylpk[k]);
   }
   CPU_dbl10_product
      (deg,xrtb,xrix,xrmi,xrrg,xrpk,xltb,xlix,xlmi,xlrg,xlpk,
           yrtb,yrix,yrmi,yrrg,yrpk,yltb,ylix,ylmi,ylrg,ylpk,
       zrtb_h,zrix_h,zrmi_h,zrrg_h,zrpk_h,
       zltb_h,zlix_h,zlmi_h,zlrg_h,zlpk_h);

   if(verbose > 0)
   {
      cout << "The product of two random real series : " << endl;

      for(int k=0; k<=deg; k++)
      {
         cout << "zrtb[" << k << "] : " << zrtb_h[k];
         cout << "  zrix[" << k << "] : " << zrix_h[k] << endl;
         cout << "zrmi[" << k << "] : " << zrmi_h[k];
         cout << "  zrrg[" << k << "] : " << zrrg_h[k] << endl;
         cout << "zrpk[" << k << "] : " << zrpk_h[k] << endl;
         cout << "zltb[" << k << "] : " << zltb_h[k];
         cout << "  zlix[" << k << "] : " << zlix_h[k] << endl;
         cout << "zlmi[" << k << "] : " << zlmi_h[k];
         cout << "  zlrg[" << k << "] : " << zlrg_h[k] << endl;
         cout << "zlpk[" << k << "] : " << zlpk_h[k] << endl;
      }
   }
   GPU_dbl10_product
      (xrtb,xrix,xrmi,xrrg,xrpk,xltb,xlix,xlmi,xlrg,xlpk,
       yrtb,yrix,yrmi,yrrg,yrpk,yltb,ylix,ylmi,ylrg,ylpk,
       zrtb_d,zrix_d,zrmi_d,zrrg_d,zrpk_d,
       zltb_d,zlix_d,zlmi_d,zlrg_d,zlpk_d,deg,1,deg+1,1);

   if(verbose > 0) cout << "GPU computed product :" << endl;

   double err = 0.0;

   for(int k=0; k<=deg; k++)
   {
      if(verbose > 0)
      {
         cout << "zrtb[" << k << "] : " << zrtb_d[k];
         cout << "  zrix[" << k << "] : " << zrix_d[k];
         cout << "  zrmi[" << k << "] : " << zrmi_d[k];
         cout << "  zrrg[" << k << "] : " << zrrg_d[k];
         cout << "  zrpk[" << k << "] : " << zrpk_d[k] << endl;
         cout << "zltb[" << k << "] : " << zltb_d[k];
         cout << "  zlix[" << k << "] : " << zlix_d[k];
         cout << "  zlmi[" << k << "] : " << zlmi_d[k];
         cout << "  zlrg[" << k << "] : " << zlrg_d[k];
         cout << "  zlpk[" << k << "] : " << zlpk_d[k] << endl;
      }
      err = err
          + abs(zrtb_h[k] - zrtb_d[k]) + abs(zrix_h[k] - zrix_d[k])
          + abs(zrmi_h[k] - zrmi_d[k]) + abs(zrrg_h[k] - zrrg_d[k])
          + abs(zrpk_h[k] - zrpk_d[k])
          + abs(zltb_h[k] - zltb_d[k]) + abs(zlix_h[k] - zlix_d[k])
          + abs(zlmi_h[k] - zlmi_d[k]) + abs(zlrg_h[k] - zlrg_d[k])
          + abs(zlpk_h[k] - zlpk_d[k]);
   }
   if(verbose > 0)
   {
      cout << endl;

      cout << scientific << setprecision(16);
      cout << "the error : " << err << endl;
   }
   return err;
}

double test_dbl10_complex ( int deg, int verbose )
{
   double* xrertb = new double[deg+1];
   double* xrerix = new double[deg+1];
   double* xrermi = new double[deg+1];
   double* xrerrg = new double[deg+1];
   double* xrerpk = new double[deg+1];
   double* xreltb = new double[deg+1];
   double* xrelix = new double[deg+1];
   double* xrelmi = new double[deg+1];
   double* xrelrg = new double[deg+1];
   double* xrelpk = new double[deg+1];
   double* ximrtb = new double[deg+1];
   double* ximrix = new double[deg+1];
   double* ximrmi = new double[deg+1];
   double* ximrrg = new double[deg+1];
   double* ximrpk = new double[deg+1];
   double* ximltb = new double[deg+1];
   double* ximlix = new double[deg+1];
   double* ximlmi = new double[deg+1];
   double* ximlrg = new double[deg+1];
   double* ximlpk = new double[deg+1];
   double* yrertb = new double[deg+1];
   double* yrerix = new double[deg+1];
   double* yrermi = new double[deg+1];
   double* yrerrg = new double[deg+1];
   double* yrerpk = new double[deg+1];
   double* yreltb = new double[deg+1];
   double* yrelix = new double[deg+1];
   double* yrelmi = new double[deg+1];
   double* yrelrg = new double[deg+1];
   double* yrelpk = new double[deg+1];
   double* yimrtb = new double[deg+1];
   double* yimrix = new double[deg+1];
   double* yimrmi = new double[deg+1];
   double* yimrrg = new double[deg+1];
   double* yimrpk = new double[deg+1];
   double* yimltb = new double[deg+1];
   double* yimlix = new double[deg+1];
   double* yimlmi = new double[deg+1];
   double* yimlrg = new double[deg+1];
   double* yimlpk = new double[deg+1];
   double* zrertb_h = new double[deg+1];
   double* zrerix_h = new double[deg+1];
   double* zrermi_h = new double[deg+1];
   double* zrerrg_h = new double[deg+1];
   double* zrerpk_h = new double[deg+1];
   double* zreltb_h = new double[deg+1];
   double* zrelix_h = new double[deg+1];
   double* zrelmi_h = new double[deg+1];
   double* zrelrg_h = new double[deg+1];
   double* zrelpk_h = new double[deg+1];
   double* zimrtb_h = new double[deg+1];
   double* zimrix_h = new double[deg+1];
   double* zimrmi_h = new double[deg+1];
   double* zimrrg_h = new double[deg+1];
   double* zimrpk_h = new double[deg+1];
   double* zimltb_h = new double[deg+1];
   double* zimlix_h = new double[deg+1];
   double* zimlmi_h = new double[deg+1];
   double* zimlrg_h = new double[deg+1];
   double* zimlpk_h = new double[deg+1];
   double* zrertb_d = new double[deg+1];
   double* zrerix_d = new double[deg+1];
   double* zrermi_d = new double[deg+1];
   double* zrerrg_d = new double[deg+1];
   double* zrerpk_d = new double[deg+1];
   double* zreltb_d = new double[deg+1];
   double* zrelix_d = new double[deg+1];
   double* zrelmi_d = new double[deg+1];
   double* zrelrg_d = new double[deg+1];
   double* zrelpk_d = new double[deg+1];
   double* zimrtb_d = new double[deg+1];
   double* zimrix_d = new double[deg+1];
   double* zimrmi_d = new double[deg+1];
   double* zimrrg_d = new double[deg+1];
   double* zimrpk_d = new double[deg+1];
   double* zimltb_d = new double[deg+1];
   double* zimlix_d = new double[deg+1];
   double* zimlmi_d = new double[deg+1];
   double* zimlrg_d = new double[deg+1];
   double* zimlpk_d = new double[deg+1];

   for(int k=0; k<=deg; k++)
   {
      xrertb[k] = 1.0; xrerix[k] = 0.0; xrermi[k] = 0.0;
      xrerrg[k] = 0.0; xrerpk[k] = 0.0;
      xreltb[k] = 0.0; xrelix[k] = 0.0; xrelmi[k] = 0.0;
      xrelrg[k] = 0.0; xrelpk[k] = 0.0;
      ximrtb[k] = 0.0; ximrix[k] = 0.0; ximrmi[k] = 0.0;
      ximrrg[k] = 0.0; ximrpk[k] = 0.0;
      ximltb[k] = 0.0; ximlix[k] = 0.0; ximlmi[k] = 0.0;
      ximlrg[k] = 0.0; ximlpk[k] = 0.0;
      yrertb[k] = 0.0; yrerix[k] = 0.0; yrermi[k] = 0.0;
      yrerrg[k] = 0.0; yrerpk[k] = 0.0;
      yreltb[k] = 0.0; yrelix[k] = 0.0; yrelmi[k] = 0.0;
      yrelrg[k] = 0.0; yrelpk[k] = 0.0;
      yimrtb[k] = 0.0; yimrix[k] = 0.0; yimrmi[k] = 0.0;
      yimrrg[k] = 0.0; yimrpk[k] = 0.0;
      yimltb[k] = 0.0; yimlix[k] = 0.0; yimlmi[k] = 0.0;
      yimlrg[k] = 0.0; yimlpk[k] = 0.0;
   }
   yrertb[0] = 1.0; yrertb[1] = -1.0;

   CPU_cmplx10_product
      (deg,
       xrertb,xrerix,xrermi,xrerrg,xrerpk,xreltb,xrelix,xrelmi,xrelrg,xrelpk,
       ximrtb,ximrix,ximrmi,ximrrg,ximrpk,ximltb,ximlix,ximlmi,ximlrg,ximlpk,
       yrertb,yrerix,yrermi,yrerrg,yrerpk,yreltb,yrelix,yrelmi,yrelrg,yrelpk,
       yimrtb,yimrix,yimrmi,yimrrg,yimrpk,yimltb,yimlix,yimlmi,yimlrg,yimlpk,
       zrertb_h,zrerix_h,zrermi_h,zrerrg_h,zrerpk_h,
       zreltb_h,zrelix_h,zrelmi_h,zrelrg_h,zrelpk_h,
       zimrtb_h,zimrix_h,zimrmi_h,zimrrg_h,zimrpk_h,
       zimltb_h,zimlix_h,zimlmi_h,zimlrg_h,zimlpk_h);

   if(verbose > 0)
   {
      cout << "Series of 1/(1-x) multiplied with 1-x :" << endl;

      for(int k=0; k<=deg; k++)
      {
         cout << "zrertb[" << k << "] : " << zrertb_h[k];
         cout << "  zrerix[" << k << "] : " << zrerix_h[k];
         cout << "  zrermi[" << k << "] : " << zrermi_h[k];
         cout << "  zrerrg[" << k << "] : " << zrerrg_h[k];
         cout << "  zrerpk[" << k << "] : " << zrerpk_h[k] << endl;
         cout << "zreltb[" << k << "] : " << zreltb_h[k];
         cout << "  zrelix[" << k << "] : " << zrelix_h[k];
         cout << "  zrelmi[" << k << "] : " << zrelmi_h[k];
         cout << "  zrelrg[" << k << "] : " << zrelrg_h[k];
         cout << "  zrelpk[" << k << "] : " << zrelpk_h[k] << endl;
         cout << "zimrtb[" << k << "] : " << zimrtb_h[k];
         cout << "  zimrix[" << k << "] : " << zimrix_h[k];
         cout << "  zimrmi[" << k << "] : " << zimrmi_h[k];
         cout << "  zimrrg[" << k << "] : " << zimrrg_h[k];
         cout << "  zimrpk[" << k << "] : " << zimrpk_h[k] << endl;
         cout << "zimltb[" << k << "] : " << zimltb_h[k];
         cout << "  zimlix[" << k << "] : " << zimlix_h[k];
         cout << "  zimlmi[" << k << "] : " << zimlmi_h[k];
         cout << "  zimlrg[" << k << "] : " << zimlrg_h[k];
         cout << "  zimlpk[" << k << "] : " << zimlpk_h[k] << endl;
      }
   }
   GPU_cmplx10_product
      (xrertb,xrerix,xrermi,xrerrg,xrerpk,xreltb,xrelix,xrelmi,xrelrg,xrelpk,
       ximrtb,ximrix,ximrmi,ximrrg,ximrpk,ximltb,ximlix,ximlmi,ximlrg,ximlpk,
       yrertb,yrerix,yrermi,yrerrg,yrerpk,yreltb,yrelix,yrelmi,yrelrg,yrelpk,
       yimrtb,yimrix,yimrmi,yimrrg,yimrpk,yimltb,yimlix,yimlmi,yimlrg,yimlpk,
       zrertb_d,zrerix_d,zrermi_d,zrerrg_d,zrerpk_d,
       zreltb_d,zrelix_d,zrelmi_d,zrelrg_d,zrelpk_d,
       zimrtb_d,zimrix_d,zimrmi_d,zimrrg_d,zimrpk_d,
       zimltb_d,zimlix_d,zimlmi_d,zimlrg_d,zimlpk_d,deg,1,deg+1,2);

   if(verbose > 0) cout << "GPU computed product :" << endl;

   double err = 0.0;

   for(int k=0; k<=deg; k++)
   {
      if(verbose > 0)
      {
         cout << "zrertb[" << k << "] : " << zrertb_d[k];
         cout << "  zrerix[" << k << "] : " << zrerix_d[k];
         cout << "  zrermi[" << k << "] : " << zrermi_d[k];
         cout << "  zrerrg[" << k << "] : " << zrerrg_d[k];
         cout << "  zrerpk[" << k << "] : " << zrerpk_d[k] << endl;
         cout << "zreltb[" << k << "] : " << zreltb_d[k];
         cout << "  zrelix[" << k << "] : " << zrelix_d[k];
         cout << "  zrelmi[" << k << "] : " << zrelmi_d[k];
         cout << "  zrelrg[" << k << "] : " << zrelrg_d[k];
         cout << "  zrelpk[" << k << "] : " << zrelpk_d[k] << endl;
         cout << "zimrtb[" << k << "] : " << zimrtb_d[k];
         cout << "  zimrix[" << k << "] : " << zimrix_d[k];
         cout << "  zimrmi[" << k << "] : " << zimrmi_d[k];
         cout << "  zimrrg[" << k << "] : " << zimrrg_d[k];
         cout << "  zimrpk[" << k << "] : " << zimrpk_d[k] << endl;
         cout << "zimltb[" << k << "] : " << zimltb_d[k];
         cout << "  zimlix[" << k << "] : " << zimlix_d[k];
         cout << "  zimlmi[" << k << "] : " << zimlmi_d[k];
         cout << "  zimlrg[" << k << "] : " << zimlrg_d[k];
         cout << "  zimlpk[" << k << "] : " << zimlpk_d[k] << endl;
      }
      err = err
          + abs(zrertb_h[k] - zrertb_d[k]) + abs(zrerix_h[k] - zrerix_d[k])
          + abs(zrermi_h[k] - zrermi_d[k]) + abs(zrerrg_h[k] - zrerrg_d[k])
          + abs(zrerpk_h[k] - zrerpk_d[k])
          + abs(zimltb_h[k] - zimltb_d[k]) + abs(zimlix_h[k] - zimlix_d[k])
          + abs(zimlmi_h[k] - zimlmi_d[k]) + abs(zimlrg_h[k] - zimlrg_d[k])
          + abs(zimlpk_h[k] - zimlpk_d[k]);
   }
   if(verbose > 0)
   {
      cout << scientific << setprecision(16);
     cout << "the error : " << err << endl;
   }
   return err;
}

double test_dbl10_complex_random ( int deg, int verbose )
{
   double* xrertb = new double[deg+1];
   double* xrerix = new double[deg+1];
   double* xrermi = new double[deg+1];
   double* xrerrg = new double[deg+1];
   double* xrerpk = new double[deg+1];
   double* xreltb = new double[deg+1];
   double* xrelix = new double[deg+1];
   double* xrelmi = new double[deg+1];
   double* xrelrg = new double[deg+1];
   double* xrelpk = new double[deg+1];
   double* ximrtb = new double[deg+1];
   double* ximrix = new double[deg+1];
   double* ximrmi = new double[deg+1];
   double* ximrrg = new double[deg+1];
   double* ximrpk = new double[deg+1];
   double* ximltb = new double[deg+1];
   double* ximlix = new double[deg+1];
   double* ximlmi = new double[deg+1];
   double* ximlrg = new double[deg+1];
   double* ximlpk = new double[deg+1];
   double* yrertb = new double[deg+1];
   double* yrerix = new double[deg+1];
   double* yrermi = new double[deg+1];
   double* yrerrg = new double[deg+1];
   double* yrerpk = new double[deg+1];
   double* yreltb = new double[deg+1];
   double* yrelix = new double[deg+1];
   double* yrelmi = new double[deg+1];
   double* yrelrg = new double[deg+1];
   double* yrelpk = new double[deg+1];
   double* yimrtb = new double[deg+1];
   double* yimrix = new double[deg+1];
   double* yimrmi = new double[deg+1];
   double* yimrrg = new double[deg+1];
   double* yimrpk = new double[deg+1];
   double* yimltb = new double[deg+1];
   double* yimlix = new double[deg+1];
   double* yimlmi = new double[deg+1];
   double* yimlrg = new double[deg+1];
   double* yimlpk = new double[deg+1];
   double* zrertb_h = new double[deg+1];
   double* zrerix_h = new double[deg+1];
   double* zrermi_h = new double[deg+1];
   double* zrerrg_h = new double[deg+1];
   double* zrerpk_h = new double[deg+1];
   double* zreltb_h = new double[deg+1];
   double* zrelix_h = new double[deg+1];
   double* zrelmi_h = new double[deg+1];
   double* zrelrg_h = new double[deg+1];
   double* zrelpk_h = new double[deg+1];
   double* zimrtb_h = new double[deg+1];
   double* zimrix_h = new double[deg+1];
   double* zimrmi_h = new double[deg+1];
   double* zimrrg_h = new double[deg+1];
   double* zimrpk_h = new double[deg+1];
   double* zimltb_h = new double[deg+1];
   double* zimlix_h = new double[deg+1];
   double* zimlmi_h = new double[deg+1];
   double* zimlrg_h = new double[deg+1];
   double* zimlpk_h = new double[deg+1];
   double* zrertb_d = new double[deg+1];
   double* zrerix_d = new double[deg+1];
   double* zrermi_d = new double[deg+1];
   double* zrerrg_d = new double[deg+1];
   double* zrerpk_d = new double[deg+1];
   double* zreltb_d = new double[deg+1];
   double* zrelix_d = new double[deg+1];
   double* zrelmi_d = new double[deg+1];
   double* zrelrg_d = new double[deg+1];
   double* zrelpk_d = new double[deg+1];
   double* zimrtb_d = new double[deg+1];
   double* zimrix_d = new double[deg+1];
   double* zimrmi_d = new double[deg+1];
   double* zimrrg_d = new double[deg+1];
   double* zimrpk_d = new double[deg+1];
   double* zimltb_d = new double[deg+1];
   double* zimlix_d = new double[deg+1];
   double* zimlmi_d = new double[deg+1];
   double* zimlrg_d = new double[deg+1];
   double* zimlpk_d = new double[deg+1];

   for(int k=0; k<=deg; k++)
   {
      random_deca_double
         (&xrertb[k],&xrerix[k],&xrermi[k],&xrerrg[k],&xrerpk[k],
          &xreltb[k],&xrelix[k],&xrelmi[k],&xrelrg[k],&xrelpk[k]);
      random_deca_double
         (&ximrtb[k],&ximrix[k],&ximrmi[k],&ximrrg[k],&ximrpk[k],
          &ximltb[k],&ximlix[k],&ximlmi[k],&ximlrg[k],&ximlpk[k]);
      random_deca_double
         (&yrertb[k],&yrerix[k],&yrermi[k],&yrerrg[k],&yrerpk[k],
          &yreltb[k],&yrelix[k],&yrelmi[k],&yrelrg[k],&yrelpk[k]);
      random_deca_double
         (&yimrtb[k],&yimrix[k],&yimrmi[k],&yimrrg[k],&yimrpk[k],
          &yimltb[k],&yimlix[k],&yimlmi[k],&yimlrg[k],&yimlpk[k]);
   }
   CPU_cmplx10_product
      (deg,
       xrertb,xrerix,xrermi,xrerrg,xrerpk,xreltb,xrelix,xrelmi,xrelrg,xrelpk,
       ximrtb,ximrix,ximrmi,ximrrg,ximrpk,ximltb,ximlix,ximlmi,ximlrg,ximlpk,
       yrertb,yrerix,yrermi,yrerrg,yrerpk,yreltb,yrelix,yrelmi,yrelrg,yrelpk,
       yimrtb,yimrix,yimrmi,yimrrg,yimrpk,yimltb,yimlix,yimlmi,yimlrg,yimlpk,
       zrertb_h,zrerix_h,zrermi_h,zrerrg_h,zrerpk_h,
       zreltb_h,zrelix_h,zrelmi_h,zrelrg_h,zrelpk_h,
       zimrtb_h,zimrix_h,zimrmi_h,zimrrg_h,zimrpk_h,
       zimltb_h,zimlix_h,zimlmi_h,zimlrg_h,zimlpk_h);

   if(verbose > 0)
   {
      cout << "Product of two random complex series :" << endl;

      for(int k=0; k<=deg; k++)
      {
         cout << "zrertb[" << k << "] : " << zrertb_h[k];
         cout << "  zrerix[" << k << "] : " << zrerix_h[k] << endl;
         cout << "zrermi[" << k << "] : " << zrermi_h[k];
         cout << "  zrerrg[" << k << "] : " << zrerrg_h[k] << endl;
         cout << "zrerpk[" << k << "] : " << zrerpk_h[k] << endl;
         cout << "zreltb[" << k << "] : " << zreltb_h[k];
         cout << "  zrelix[" << k << "] : " << zrelix_h[k] << endl;
         cout << "zrelmi[" << k << "] : " << zrelmi_h[k];
         cout << "  zrelrg[" << k << "] : " << zrelrg_h[k] << endl;
         cout << "zrelpk[" << k << "] : " << zrelpk_h[k] << endl;
         cout << "zimrtb[" << k << "] : " << zimrtb_h[k];
         cout << "  zimrix[" << k << "] : " << zimrix_h[k] << endl;
         cout << "zimrmi[" << k << "] : " << zimrmi_h[k];
         cout << "  zimrrg[" << k << "] : " << zimrrg_h[k] << endl;
         cout << "zimrpk[" << k << "] : " << zimrpk_h[k] << endl;
         cout << "zimltb[" << k << "] : " << zimltb_h[k];
         cout << "  zimlix[" << k << "] : " << zimlix_h[k] << endl;
         cout << "zimlmi[" << k << "] : " << zimlmi_h[k];
         cout << "  zimlrg[" << k << "] : " << zimlrg_h[k] << endl;
         cout << "zimlpk[" << k << "] : " << zimlpk_h[k] << endl;
      }
   }
   GPU_cmplx10_product
      (xrertb,xrerix,xrermi,xrerrg,xrerpk,xreltb,xrelix,xrelmi,xrelrg,xrelpk,
       ximrtb,ximrix,ximrmi,ximrrg,ximrpk,ximltb,ximlix,ximlmi,ximlrg,ximlpk,
       yrertb,yrerix,yrermi,yrerrg,yrerpk,yreltb,yrelix,yrelmi,yrelrg,yrelpk,
       yimrtb,yimrix,yimrmi,yimrrg,yimrpk,yimltb,yimlix,yimlmi,yimlrg,yimlpk,
       zrertb_d,zrerix_d,zrermi_d,zrerrg_d,zrerpk_d,
       zreltb_d,zrelix_d,zrelmi_d,zrelrg_d,zrelpk_d,
       zimrtb_d,zimrix_d,zimrmi_d,zimrrg_d,zimrpk_d,
       zimltb_d,zimlix_d,zimlmi_d,zimlrg_d,zimlpk_d,deg,1,deg+1,2);

   if(verbose > 0) cout << "GPU computed product :" << endl;

   double err = 0.0;

   for(int k=0; k<=deg; k++)
   {
      if(verbose > 0)
      {
         cout << "zrertb[" << k << "] : " << zrertb_d[k];
         cout << "  zrerix[" << k << "] : " << zrerix_d[k] << endl;
         cout << "zrermi[" << k << "] : " << zrermi_d[k];
         cout << "  zrerrg[" << k << "] : " << zrerrg_d[k] << endl;
         cout << "zrerpk[" << k << "] : " << zrerpk_d[k] << endl;
         cout << "zreltb[" << k << "] : " << zreltb_d[k];
         cout << "  zrelix[" << k << "] : " << zrelix_d[k] << endl;
         cout << "zrelmi[" << k << "] : " << zrelmi_d[k];
         cout << "  zrelrg[" << k << "] : " << zrelrg_d[k] << endl;
         cout << "zrelpk[" << k << "] : " << zrelpk_d[k] << endl;
         cout << "zimrtb[" << k << "] : " << zimrtb_d[k];
         cout << "  zimrix[" << k << "] : " << zimrix_d[k] << endl;
         cout << "zimrmi[" << k << "] : " << zimrmi_d[k];
         cout << "  zimrrg[" << k << "] : " << zimrrg_d[k] << endl;
         cout << "zimrpk[" << k << "] : " << zimrpk_d[k] << endl;
         cout << "zimltb[" << k << "] : " << zimltb_d[k];
         cout << "  zimlix[" << k << "] : " << zimlix_d[k] << endl;
         cout << "zimlmi[" << k << "] : " << zimlmi_d[k];
         cout << "  zimlrg[" << k << "] : " << zimlrg_d[k] << endl;
         cout << "zimlpk[" << k << "] : " << zimlpk_d[k] << endl;
      }
      err = err
          + abs(zrertb_h[k] - zrertb_d[k]) + abs(zrerix_h[k] - zrerix_d[k])
          + abs(zrermi_h[k] - zrermi_d[k]) + abs(zrerrg_h[k] - zrerrg_d[k])
          + abs(zrerpk_h[k] - zrerpk_d[k])
          + abs(zimltb_h[k] - zimltb_d[k]) + abs(zimlix_h[k] - zimlix_d[k])
          + abs(zimlmi_h[k] - zimlmi_d[k]) + abs(zimlrg_h[k] - zimlrg_d[k])
          + abs(zimlpk_h[k] - zimlpk_d[k]);
   }
   if(verbose > 0)
   {
      cout << scientific << setprecision(16);
     cout << "the error : " << err << endl;
   }
   return err;
}

double test_dbl10_real_exponential ( int deg, int verbose )
{
   double *xrtb = new double[deg+1];
   double *xrix = new double[deg+1];
   double *xrmi = new double[deg+1];
   double *xrrg = new double[deg+1];
   double *xrpk = new double[deg+1];
   double *xltb = new double[deg+1];
   double *xlix = new double[deg+1];
   double *xlmi = new double[deg+1];
   double *xlrg = new double[deg+1];
   double *xlpk = new double[deg+1];
   double *yrtb = new double[deg+1];
   double *yrix = new double[deg+1];
   double *yrmi = new double[deg+1];
   double *yrrg = new double[deg+1];
   double *yrpk = new double[deg+1];
   double *yltb = new double[deg+1];
   double *ylix = new double[deg+1];
   double *ylmi = new double[deg+1];
   double *ylrg = new double[deg+1];
   double *ylpk = new double[deg+1];
   double *zrtb_h = new double[deg+1];
   double *zrix_h = new double[deg+1];
   double *zrmi_h = new double[deg+1];
   double *zrrg_h = new double[deg+1];
   double *zrpk_h = new double[deg+1];
   double *zltb_h = new double[deg+1];
   double *zlix_h = new double[deg+1];
   double *zlmi_h = new double[deg+1];
   double *zlrg_h = new double[deg+1];
   double *zlpk_h = new double[deg+1];
   double *zrtb_d = new double[deg+1];
   double *zrix_d = new double[deg+1];
   double *zrmi_d = new double[deg+1];
   double *zrrg_d = new double[deg+1];
   double *zrpk_d = new double[deg+1];
   double *zltb_d = new double[deg+1];
   double *zlix_d = new double[deg+1];
   double *zlmi_d = new double[deg+1];
   double *zlrg_d = new double[deg+1];
   double *zlpk_d = new double[deg+1];
   double rrtb,rrix,rrmi,rrrg,rrpk;
   double rltb,rlix,rlmi,rlrg,rlpk;
   double sumrtb,sumrix,sumrmi,sumrrg,sumrpk;
   double sumltb,sumlix,sumlmi,sumlrg,sumlpk;

   random_dbl10_exponentials
      (deg,&rrtb,&rrix,&rrmi,&rrrg,&rrpk,&rltb,&rlix,&rlmi,&rlrg,&rlpk,
           xrtb,xrix,xrmi,xrrg,xrpk,xltb,xlix,xlmi,xlrg,xlpk,
           yrtb,yrix,yrmi,yrrg,yrpk,yltb,ylix,ylmi,ylrg,ylpk);

   CPU_dbl10_product
      (deg,xrtb,xrix,xrmi,xrrg,xrpk,xltb,xlix,xlmi,xlrg,xlpk,
           yrtb,yrix,yrmi,yrrg,yrpk,yltb,ylix,ylmi,ylrg,ylpk,
           zrtb_h,zrix_h,zrmi_h,zrrg_h,zrpk_h,
           zltb_h,zlix_h,zlmi_h,zlrg_h,zlpk_h);

   if(verbose > 0)
   {
      cout << scientific << setprecision(16);

      cout << "Product of series of exp(x) with series of exp(-x)," << endl;
      cout << "  for xrtb = " << rrtb << endl;
      cout << "   xrix = " << rrix << endl;
      cout << "   xrmi = " << rrmi << endl;
      cout << "   xrrg = " << rrrg << endl;
      cout << "   xrpk = " << rrpk << endl;
      cout << "   xltb = " << rltb << endl;
      cout << "   xlix = " << rlix << endl;
      cout << "   xlmi = " << rlmi << endl;
      cout << "   xlrg = " << rlrg << endl;
      cout << "   xlpk = " << rlpk << endl;
   }
/*
   for(int k=0; k<=deg; k++)
   {
      cout << "zrtb[" << k << "] : " << zrtb_h[k];
      cout << "  zrix[" << k << "] : " << zrix_h[k];
      cout << "  zrmi[" << k << "] : " << zrmi_h[k];
      cout << "  zrrg[" << k << "] : " << zrrg_h[k];
      cout << "  zrpk[" << k << "] : " << zrpk_h[k] << endl;
      cout << "zltb[" << k << "] : " << zltb_h[k];
      cout << "  zlix[" << k << "] : " << zlix_h[k];
      cout << "  zlmi[" << k << "] : " << zlmi_h[k];
      cout << "  zlrg[" << k << "] : " << zlrg_h[k];
      cout << "  zlpk[" << k << "] : " << zlpk_h[k] << endl;
   }
 */
   if(verbose > 0)
   {
      sumrtb = 0.0; sumrix = 0.0; sumrmi = 0.0;
      sumrrg = 0.0; sumrpk = 0.0;
      sumltb = 0.0; sumlix = 0.0; sumlmi = 0.0;
      sumlrg = 0.0; sumlpk = 0.0;

      for(int k=0; k<=deg; k++)
         daf_inc(&sumrtb,&sumrix,&sumrmi,&sumrrg,&sumrpk,
                 &sumltb,&sumlix,&sumlmi,&sumlrg,&sumlpk,
                 zrtb_h[k],zrix_h[k],zrmi_h[k],zrrg_h[k],zrpk_h[k],
                 zltb_h[k],zlix_h[k],zlmi_h[k],zlrg_h[k],zlpk_h[k]);

      cout << "Summation of all coefficients in the product ..." << endl;
      cout << "       highest part of the sum : " << sumrtb << endl;
      cout << "second highest part of the sum : " << sumrix << endl;
      cout << " third highest part of the sum : " << sumrmi << endl;
      cout << "fourth highest part of the sum : " << sumrrg << endl;
      cout << " fifth highest part of the sum : " << sumrpk << endl;
      cout << "  fifth lowest part of the sum : " << sumltb << endl;
      cout << " fourth lowest part of the sum : " << sumlix << endl;
      cout << "  third lowest part of the sum : " << sumlmi << endl;
      cout << " second lowest part of the sum : " << sumlrg << endl;
      cout << "        lowest part of the sum : " << sumlpk << endl;
   }
   GPU_dbl10_product
      (xrtb,xrix,xrmi,xrrg,xrpk,xltb,xlix,xlmi,xlrg,xlpk,
       yrtb,yrix,yrmi,yrrg,yrpk,yltb,ylix,ylmi,ylrg,ylpk,
       zrtb_d,zrix_d,zrmi_d,zrrg_d,zrpk_d,
       zltb_d,zlix_d,zlmi_d,zlrg_d,zlpk_d,deg,1,deg+1,1);

   if(verbose > 0)
   {
      sumrtb = 0.0; sumrix = 0.0; sumrmi = 0.0; sumrrg = 0.0; sumrpk = 0.0;
      sumltb = 0.0; sumlix = 0.0; sumlmi = 0.0; sumlrg = 0.0; sumlpk = 0.0;
   }
   double err = 0.0;

   for(int k=0; k<=deg; k++)
   {
      if(verbose > 0)
      {
         daf_inc(&sumrtb,&sumrix,&sumrmi,&sumrrg,&sumrpk,
                 &sumltb,&sumlix,&sumlmi,&sumlrg,&sumlpk,
                 zrtb_d[k],zrix_d[k],zrmi_d[k],zrrg_d[k],zrpk_d[k],
                 zltb_d[k],zlix_d[k],zlmi_d[k],zlrg_d[k],zlpk_d[k]);
      }
      err = err
          + abs(zrtb_h[k] - zrtb_d[k]) + abs(zrix_h[k] - zrix_d[k])
          + abs(zrmi_h[k] - zrmi_d[k]) + abs(zrrg_h[k] - zrrg_d[k])
          + abs(zrpk_h[k] - zrpk_d[k])
          + abs(zltb_h[k] - zltb_d[k]) + abs(zlix_h[k] - zlix_d[k])
          + abs(zlmi_h[k] - zlmi_d[k]) + abs(zlrg_h[k] - zlrg_d[k])
          + abs(zlpk_h[k] - zlpk_d[k]);
   }
   if(verbose > 0)
   {
      cout << "Summation of all coefficients in the GPU computed product ..."
           << endl;
      cout << "       highest part of the sum : " << sumrtb << endl;
      cout << "second highest part of the sum : " << sumrix << endl;
      cout << " third highest part of the sum : " << sumrmi << endl;
      cout << "fourth highest part of the sum : " << sumrrg << endl;
      cout << " fifth highest part of the sum : " << sumrpk << endl;
      cout << "  fifth lowest part of the sum : " << sumltb << endl;
      cout << " fourth lowest part of the sum : " << sumlix << endl;
      cout << "  third lowest part of the sum : " << sumlmi << endl;
      cout << " second lowest part of the sum : " << sumlrg << endl;
      cout << "        lowest part of the sum : " << sumlpk << endl;

      cout << "the error : " << err << endl;
   }
   return err;
}

double test_dbl10_complex_exponential ( int deg, int verbose )
{
   double* xrertb = new double[deg+1];
   double* xrerix = new double[deg+1];
   double* xrermi = new double[deg+1];
   double* xrerrg = new double[deg+1];
   double* xrerpk = new double[deg+1];
   double* xreltb = new double[deg+1];
   double* xrelix = new double[deg+1];
   double* xrelmi = new double[deg+1];
   double* xrelrg = new double[deg+1];
   double* xrelpk = new double[deg+1];
   double* ximrtb = new double[deg+1];
   double* ximrix = new double[deg+1];
   double* ximrmi = new double[deg+1];
   double* ximrrg = new double[deg+1];
   double* ximrpk = new double[deg+1];
   double* ximltb = new double[deg+1];
   double* ximlix = new double[deg+1];
   double* ximlmi = new double[deg+1];
   double* ximlrg = new double[deg+1];
   double* ximlpk = new double[deg+1];
   double* yrertb = new double[deg+1];
   double* yrerix = new double[deg+1];
   double* yrermi = new double[deg+1];
   double* yrerrg = new double[deg+1];
   double* yrerpk = new double[deg+1];
   double* yreltb = new double[deg+1];
   double* yrelix = new double[deg+1];
   double* yrelmi = new double[deg+1];
   double* yrelrg = new double[deg+1];
   double* yrelpk = new double[deg+1];
   double* yimrtb = new double[deg+1];
   double* yimrix = new double[deg+1];
   double* yimrmi = new double[deg+1];
   double* yimrrg = new double[deg+1];
   double* yimrpk = new double[deg+1];
   double* yimltb = new double[deg+1];
   double* yimlix = new double[deg+1];
   double* yimlmi = new double[deg+1];
   double* yimlrg = new double[deg+1];
   double* yimlpk = new double[deg+1];
   double* zrertb_h = new double[deg+1];
   double* zrerix_h = new double[deg+1];
   double* zrermi_h = new double[deg+1];
   double* zrerrg_h = new double[deg+1];
   double* zrerpk_h = new double[deg+1];
   double* zreltb_h = new double[deg+1];
   double* zrelix_h = new double[deg+1];
   double* zrelmi_h = new double[deg+1];
   double* zrelrg_h = new double[deg+1];
   double* zrelpk_h = new double[deg+1];
   double* zimrtb_h = new double[deg+1];
   double* zimrix_h = new double[deg+1];
   double* zimrmi_h = new double[deg+1];
   double* zimrrg_h = new double[deg+1];
   double* zimrpk_h = new double[deg+1];
   double* zimltb_h = new double[deg+1];
   double* zimlix_h = new double[deg+1];
   double* zimlmi_h = new double[deg+1];
   double* zimlrg_h = new double[deg+1];
   double* zimlpk_h = new double[deg+1];
   double* zrertb_d = new double[deg+1];
   double* zrerix_d = new double[deg+1];
   double* zrermi_d = new double[deg+1];
   double* zrerrg_d = new double[deg+1];
   double* zrerpk_d = new double[deg+1];
   double* zreltb_d = new double[deg+1];
   double* zrelix_d = new double[deg+1];
   double* zrelmi_d = new double[deg+1];
   double* zrelrg_d = new double[deg+1];
   double* zrelpk_d = new double[deg+1];
   double* zimrtb_d = new double[deg+1];
   double* zimrix_d = new double[deg+1];
   double* zimrmi_d = new double[deg+1];
   double* zimrrg_d = new double[deg+1];
   double* zimrpk_d = new double[deg+1];
   double* zimltb_d = new double[deg+1];
   double* zimlix_d = new double[deg+1];
   double* zimlmi_d = new double[deg+1];
   double* zimlrg_d = new double[deg+1];
   double* zimlpk_d = new double[deg+1];
   double rndrertb,rndrerix,rndrermi,rndrerrg,rndrerpk;
   double rndreltb,rndrelix,rndrelmi,rndrelrg,rndrelpk;
   double rndimrtb,rndimrix,rndimrmi,rndimrrg,rndimrpk;
   double rndimltb,rndimlix,rndimlmi,rndimlrg,rndimlpk;
   double sumrertb,sumrerix,sumrermi,sumrerrg,sumrerpk;
   double sumreltb,sumrelix,sumrelmi,sumrelrg,sumrelpk;
   double sumimrtb,sumimrix,sumimrmi,sumimrrg,sumimrpk;
   double sumimltb,sumimlix,sumimlmi,sumimlrg,sumimlpk;

   random_cmplx10_exponentials
      (deg,&rndrertb,&rndrerix,&rndrermi,&rndrerrg,&rndrerpk,
           &rndreltb,&rndrelix,&rndrelmi,&rndrelrg,&rndrelpk,
           &rndimrtb,&rndimrix,&rndimrmi,&rndimrrg,&rndimrpk,
           &rndimltb,&rndimlix,&rndimlmi,&rndimlrg,&rndimlpk,
           xrertb,xrerix,xrermi,xrerrg,xrerpk,
           xreltb,xrelix,xrelmi,xrelrg,xrelpk,
           ximrtb,ximrix,ximrmi,ximrrg,ximrpk,
           ximltb,ximlix,ximlmi,ximlrg,ximlpk,
           yrertb,yrerix,yrermi,yrerrg,yrerpk,
           yreltb,yrelix,yrelmi,yrelrg,yrelpk,
           yimrtb,yimrix,yimrmi,yimrrg,yimrpk,
           yimltb,yimlix,yimlmi,yimlrg,yimlpk);

   CPU_cmplx10_product
      (deg,xrertb,xrerix,xrermi,xrerrg,xrerpk,
           xreltb,xrelix,xrelmi,xrelrg,xrelpk,
           ximrtb,ximrix,ximrmi,ximrrg,ximrpk,
           ximltb,ximlix,ximlmi,ximlrg,ximlpk,
           yrertb,yrerix,yrermi,yrerrg,yrerpk,
           yreltb,yrelix,yrelmi,yrelrg,yrelpk,
           yimrtb,yimrix,yimrmi,yimrrg,yimrpk,
           yimltb,yimlix,yimlmi,yimlrg,yimlpk,
           zrertb_h,zrerix_h,zrermi_h,zrerrg_h,zrerpk_h,
           zreltb_h,zrelix_h,zrelmi_h,zrelrg_h,zrelpk_h,
           zimrtb_h,zimrix_h,zimrmi_h,zimrrg_h,zimrpk_h,
           zimltb_h,zimlix_h,zimlmi_h,zimlrg_h,zimlpk_h);
  
   if(verbose > 0)
   {
      cout << scientific << setprecision(16);

      cout << "Product of series of exp(x) with series of exp(-x)," << endl;
      cout << "  for xrertb = " << rndrertb;
      cout << "    xrerix = " << rndrerix << endl;
      cout << "    xrermi = " << rndrermi;
      cout << "    xrerrg = " << rndrerrg << endl;
      cout << "    xrerpk = " << rndrerpk << endl;
      cout << "    xreltb = " << rndreltb;
      cout << "    xrelix = " << rndrelix << endl;
      cout << "    xrelmi = " << rndrelmi;
      cout << "    xrelrg = " << rndrelrg << endl;
      cout << "    xrelpk = " << rndrelpk << endl;
      cout << "  for ximrtb = " << rndimrtb;
      cout << "    ximrix = " << rndimrix << endl;
      cout << "    ximrmi = " << rndimrmi;
      cout << "    ximrrg = " << rndimrrg << endl;
      cout << "    ximrpk = " << rndimrpk << endl;
      cout << "    ximltb = " << rndimltb;
      cout << "    ximlix = " << rndimlix << endl;
      cout << "    ximlmi = " << rndimlmi;
      cout << "    ximlrg = " << rndimlrg << endl;
      cout << "    ximlpk = " << rndimlpk << endl;
   
      sumrertb = 0.0; sumrerix = 0.0; sumrermi = 0.0;
      sumrerrg = 0.0; sumrerpk = 0.0;
      sumreltb = 0.0; sumrelix = 0.0; sumrelmi = 0.0;
      sumrelrg = 0.0; sumrelpk = 0.0;
      sumimrtb = 0.0; sumimrix = 0.0; sumimrmi = 0.0;
      sumimrrg = 0.0; sumimrpk = 0.0;
      sumimltb = 0.0; sumimlix = 0.0; sumimlmi = 0.0;
      sumimlrg = 0.0; sumimlpk = 0.0;

      for(int k=0; k<=deg; k++) 
      {
         daf_inc(&sumrertb,&sumrerix,&sumrermi,&sumrerrg,&sumrerpk,
                 &sumreltb,&sumrelix,&sumrelmi,&sumrelrg,&sumrelpk,
                 zrertb_h[k],zrerix_h[k],zrermi_h[k],zrerrg_h[k],zrerpk_h[k],
                 zreltb_h[k],zrelix_h[k],zrelmi_h[k],zrelrg_h[k],zrelpk_h[k]);
         daf_inc(&sumimrtb,&sumimrix,&sumimrmi,&sumimrrg,&sumimrpk,
                 &sumimltb,&sumimlix,&sumimlmi,&sumimlrg,&sumimlpk,
                 zimrtb_h[k],zimrix_h[k],zimrmi_h[k],zimrrg_h[k],zimrpk_h[k],
                 zimltb_h[k],zimlix_h[k],zimlmi_h[k],zimlrg_h[k],zimlpk_h[k]);
      }
      cout << "Summation of all coefficients of the product ..." << endl;
      cout << "  sumrertb : " << sumrertb;
      cout << "  sumrerix : " << sumrerix << endl;
      cout << "  sumrermi : " << sumrermi;
      cout << "  sumrerrg : " << sumrerrg << endl;
      cout << "  sumrerpk : " << sumrerpk << endl;
      cout << "  sumreltb : " << sumreltb;
      cout << "  sumrelix : " << sumrelix << endl;
      cout << "  sumrelmi : " << sumrelmi;
      cout << "  sumrelrg : " << sumrelrg << endl;
      cout << "  sumrelpk : " << sumrelpk << endl;
      cout << "  sumimrtb : " << sumimrtb;
      cout << "  sumimrix : " << sumimrix << endl;
      cout << "  sumimrmi : " << sumimrmi;
      cout << "  sumimrrg : " << sumimrrg << endl;
      cout << "  sumimrpk : " << sumimrpk << endl;
      cout << "  sumimltb : " << sumimltb;
      cout << "  sumimlix : " << sumimlix << endl;
      cout << "  sumimlmi : " << sumimlmi;
      cout << "  sumimlrg : " << sumimlrg << endl;
      cout << "  sumimlpk : " << sumimlpk << endl;
   }
   GPU_cmplx10_product
      (xrertb,xrerix,xrermi,xrerrg,xrerpk,xreltb,xrelix,xrelmi,xrelrg,xrelpk,
       ximrtb,ximrix,ximrmi,ximrrg,ximrpk,ximltb,ximlix,ximlmi,ximlrg,ximlpk,
       yrertb,yrerix,yrermi,yrerrg,yrerpk,yreltb,yrelix,yrelmi,yrelrg,yrelpk,
       yimrtb,yimrix,yimrmi,yimrrg,yimrpk,yimltb,yimlix,yimlmi,yimlrg,yimlpk,
       zrertb_d,zrerix_d,zrermi_d,zrerrg_d,zrerpk_d,
       zreltb_d,zrelix_d,zrelmi_d,zrelrg_d,zrelpk_d,
       zimrtb_d,zimrix_d,zimrmi_d,zimrrg_d,zimrpk_d,
       zimltb_d,zimlix_d,zimlmi_d,zimlrg_d,zimlpk_d,deg,1,deg+1,2);

   if(verbose > 0)
   {
      sumrertb = 0.0; sumrerix = 0.0; sumrermi = 0.0;
      sumrerrg = 0.0; sumrerpk = 0.0;
      sumreltb = 0.0; sumrelix = 0.0; sumrelmi = 0.0;
      sumrelrg = 0.0; sumrelpk = 0.0;
      sumimrtb = 0.0; sumimrix = 0.0; sumimrmi = 0.0;
      sumimrrg = 0.0; sumimrpk = 0.0;
      sumimltb = 0.0; sumimlix = 0.0; sumimlmi = 0.0;
      sumimlrg = 0.0; sumimlpk = 0.0;
   }
   double err = 0.0;

   for(int k=0; k<=deg; k++) 
   {
      if(verbose > 0)
      {
         daf_inc(&sumrertb,&sumrerix,&sumrermi,&sumrerrg,&sumrerpk,
                 &sumreltb,&sumrelix,&sumrelmi,&sumrelrg,&sumrelpk,
                 zrertb_d[k],zrerix_d[k],zrermi_d[k],zrerrg_d[k],zrerpk_d[k],
                 zreltb_d[k],zrelix_d[k],zrelmi_d[k],zrelrg_d[k],zrelpk_d[k]);
         daf_inc(&sumimrtb,&sumimrix,&sumimrmi,&sumimrrg,&sumimrpk,
                 &sumimltb,&sumimlix,&sumimlmi,&sumimlrg,&sumimlpk,
                 zimrtb_d[k],zimrix_d[k],zimrmi_d[k],zimrrg_d[k],zimrpk_d[k],
                 zimltb_d[k],zimlix_d[k],zimlmi_d[k],zimlrg_d[k],zimlpk_d[k]);
      }
      err = err
          + abs(zrertb_h[k] - zrertb_d[k]) + abs(zrerix_h[k] - zrerix_d[k])
          + abs(zrermi_h[k] - zrermi_d[k]) + abs(zrerrg_h[k] - zrerrg_d[k])
          + abs(zrerpk_h[k] - zrerpk_d[k])
          + abs(zimltb_h[k] - zimltb_d[k]) + abs(zimlix_h[k] - zimlix_d[k])
          + abs(zimlmi_h[k] - zimlmi_d[k]) + abs(zimlrg_h[k] - zimlrg_d[k])
          + abs(zimlpk_h[k] - zimlpk_d[k]);
   }
   if(verbose > 0)
   {
      cout << "Summation of all coefficients of the GPU computed product ..."
           << endl;
      cout << "  sumrertb : " << sumrertb;
      cout << "  sumrerix : " << sumrerix << endl;
      cout << "  sumrermi : " << sumrermi;
      cout << "  sumrerrg : " << sumrerrg << endl;
      cout << "  sumrerpk : " << sumrerpk << endl;
      cout << "  sumreltb : " << sumreltb;
      cout << "  sumrelix : " << sumrelix << endl;
      cout << "  sumrelmi : " << sumrelmi;
      cout << "  sumrelrg : " << sumrelrg << endl;
      cout << "  sumrelpk : " << sumrelpk << endl;
      cout << "  sumimrtb : " << sumimrtb;
      cout << "  sumimrix : " << sumimrix << endl;
      cout << "  sumimrmi : " << sumimrmi;
      cout << "  sumimrrg : " << sumimrrg << endl;
      cout << "  sumimrpk : " << sumimrpk << endl;
      cout << "  sumimltb : " << sumimltb;
      cout << "  sumimlix : " << sumimlix << endl;
      cout << "  sumimlmi : " << sumimlmi;
      cout << "  sumimlrg : " << sumimlrg << endl;
      cout << "  sumimlpk : " << sumimlpk << endl;

      cout << scientific << setprecision(16);
      cout << "the error : " << err << endl;
   }
   return err;
}
