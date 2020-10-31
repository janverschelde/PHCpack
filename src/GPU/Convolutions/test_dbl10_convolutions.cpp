/* Tests the product of two series in deca double precision. */

#include <iostream>
#include <iomanip>
#include <vector_types.h>
#include "deca_double_functions.h"
#include "random10_vectors.h"
#include "dbl10_convolutions_host.h"
#include "dbl10_convolutions_kernels.h"

using namespace std;

void test_real ( int deg );
/*
 * DESCRIPTION :
 *   Multiplies the power series of 1/(1-x) with 1+x,
 *   truncated to degree deg, for real coefficients. */

void test_real_random ( int deg );
/*
 * DESCRIPTION :
 *   Multiplies two power series of degree deg
 *   with random real coefficients. */

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
   if(deg > 0) test_real_random(deg);
   if(deg > 0) test_complex(deg);
   if(deg > 0) test_real_exponential(deg);
   if(deg > 0) test_complex_exponential(deg);

   return 0;
}

void test_real ( int deg )
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

   GPU_dbl10_product
      (xrtb,xrix,xrmi,xrrg,xrpk,xltb,xlix,xlmi,xlrg,xlpk,
       yrtb,yrix,yrmi,yrrg,yrpk,yltb,ylix,ylmi,ylrg,ylpk,
       zrtb_d,zrix_d,zrmi_d,zrrg_d,zrpk_d,
       zltb_d,zlix_d,zlmi_d,zlrg_d,zlpk_d,deg,1,deg+1);

   cout << "GPU computed product :" << endl;

   for(int k=0; k<=deg; k++)
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
   cout << endl;
}

void test_real_random ( int deg )
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

   cout << "The product of two random series : " << endl;

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

   GPU_dbl10_product
      (xrtb,xrix,xrmi,xrrg,xrpk,xltb,xlix,xlmi,xlrg,xlpk,
       yrtb,yrix,yrmi,yrrg,yrpk,yltb,ylix,ylmi,ylrg,ylpk,
       zrtb_d,zrix_d,zrmi_d,zrrg_d,zrpk_d,
       zltb_d,zlix_d,zlmi_d,zlrg_d,zlpk_d,deg,1,deg+1);

   cout << "GPU computed product :" << endl;

   for(int k=0; k<=deg; k++)
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
   cout << endl;
}

void test_complex ( int deg )
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

   GPU_cmplx10_product
      (xrertb,xrerix,xrermi,xrerrg,xrerpk,xreltb,xrelix,xrelmi,xrelrg,xrelpk,
       ximrtb,ximrix,ximrmi,ximrrg,ximrpk,ximltb,ximlix,ximlmi,ximlrg,ximlpk,
       yrertb,yrerix,yrermi,yrerrg,yrerpk,yreltb,yrelix,yrelmi,yrelrg,yrelpk,
       yimrtb,yimrix,yimrmi,yimrrg,yimrpk,yimltb,yimlix,yimlmi,yimlrg,yimlpk,
       zrertb_d,zrerix_d,zrermi_d,zrerrg_d,zrerpk_d,
       zreltb_d,zrelix_d,zrelmi_d,zrelrg_d,zrelpk_d,
       zimrtb_d,zimrix_d,zimrmi_d,zimrrg_d,zimrpk_d,
       zimltb_d,zimlix_d,zimlmi_d,zimlrg_d,zimlpk_d,deg,1,deg+1);

   cout << "GPU computed product :" << endl;

   for(int k=0; k<=deg; k++)
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
}

void test_real_exponential ( int deg )
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
   double frtb,frix,frmi,frrg,frpk;
   double fltb,flix,flmi,flrg,flpk;

   random_deca_double
      (&rrtb,&rrix,&rrmi,&rrrg,&rrpk,&rltb,&rlix,&rlmi,&rlrg,&rlpk);

   xrtb[0] = 1.0; xrix[0] = 0.0; xrmi[0] = 0.0; xrrg[0] = 0.0; xrpk[0] = 0.0;
   xltb[0] = 0.0; xlix[0] = 0.0; xlmi[0] = 0.0; xlrg[0] = 0.0; xlpk[0] = 0.0;
   yrtb[0] = 1.0; yrix[0] = 0.0; yrmi[0] = 0.0; yrrg[0] = 0.0; yrpk[0] = 0.0;
   yltb[0] = 0.0; ylix[0] = 0.0; ylmi[0] = 0.0; ylrg[0] = 0.0; ylpk[0] = 0.0;
   xrtb[1] = rrtb; yrtb[1] = -rrtb;
   xrix[1] = rrix; yrix[1] = -rrix;
   xrmi[1] = rrmi; yrmi[1] = -rrmi;
   xrrg[1] = rrrg; yrrg[1] = -rrrg;
   xrpk[1] = rrpk; yrpk[1] = -rrpk;
   xltb[1] = rltb; yltb[1] = -rltb;
   xlix[1] = rlix; ylix[1] = -rlix;
   xlmi[1] = rlmi; ylmi[1] = -rlmi;
   xlrg[1] = rlrg; ylrg[1] = -rlrg;
   xlpk[1] = rlpk; ylpk[1] = -rlpk;

   for(int k=2; k<=deg; k++)
   {
      daf_mul(xrtb[k-1],xrix[k-1],xrmi[k-1],xrrg[k-1],xrpk[k-1],
              xltb[k-1],xlix[k-1],xlmi[k-1],xlrg[k-1],xlpk[k-1],
              rrtb,rrix,rrmi,rrrg,rrpk,rltb,rlix,rlmi,rlrg,rlpk,
              &xrtb[k],&xrix[k],&xrmi[k],&xrrg[k],&xrpk[k],
              &xltb[k],&xlix[k],&xlmi[k],&xlrg[k],&xlpk[k]);
      // x[k] = x[k-1]*r
      daf_mul(yrtb[k-1],yrix[k-1],yrmi[k-1],yrrg[k-1],yrpk[k-1],
              yltb[k-1],ylix[k-1],ylmi[k-1],ylrg[k-1],ylpk[k-1],
              -rrtb,-rrix,-rrmi,-rrrg,-rrpk,
              -rltb,-rlix,-rlmi,-rlrg,-rlpk,
              &yrtb[k],&yrix[k],&yrmi[k],&yrrg[k],&yrpk[k],
              &yltb[k],&ylix[k],&ylmi[k],&ylrg[k],&ylpk[k]); 
      // y[k] = y[k-1]*(-r);
      frtb = (double) k;
                  frix = 0.0; frmi = 0.0; frrg = 0.0; frpk = 0.0;
      fltb = 0.0; flix = 0.0; flmi = 0.0; flrg = 0.0; flpk = 0.0;
      daf_div(xrtb[k],xrix[k],xrmi[k],xrrg[k],xrpk[k],
              xltb[k],xlix[k],xlmi[k],xlrg[k],xlpk[k],
              frtb,frix,frmi,frrg,frpk,fltb,flix,flmi,flrg,flpk,
              &xrtb[k],&xrix[k],&xrmi[k],&xrrg[k],&xrpk[k],
              &xltb[k],&xlix[k],&xlmi[k],&xlrg[k],&xlpk[k]);
      daf_div(yrtb[k],yrix[k],yrmi[k],yrrg[k],yrpk[k],
              yltb[k],ylix[k],ylmi[k],ylrg[k],ylpk[k],
              frtb,frix,frmi,frrg,frpk,fltb,flix,flmi,flrg,flpk,
              &yrtb[k],&yrix[k],&yrmi[k],&yrrg[k],&yrpk[k],
              &yltb[k],&ylix[k],&ylmi[k],&ylrg[k],&ylpk[k]);
   }

   CPU_dbl10_product
      (deg,xrtb,xrix,xrmi,xrrg,xrpk,xltb,xlix,xlmi,xlrg,xlpk,
           yrtb,yrix,yrmi,yrrg,yrpk,yltb,ylix,ylmi,ylrg,ylpk,
           zrtb_h,zrix_h,zrmi_h,zrrg_h,zrpk_h,
           zltb_h,zlix_h,zlmi_h,zlrg_h,zlpk_h);

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
   double sumrtb = 0.0;
   double sumrix = 0.0;
   double sumrmi = 0.0;
   double sumrrg = 0.0;
   double sumrpk = 0.0;
   double sumltb = 0.0;
   double sumlix = 0.0;
   double sumlmi = 0.0;
   double sumlrg = 0.0;
   double sumlpk = 0.0;

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

   GPU_dbl10_product
      (xrtb,xrix,xrmi,xrrg,xrpk,xltb,xlix,xlmi,xlrg,xlpk,
       yrtb,yrix,yrmi,yrrg,yrpk,yltb,ylix,ylmi,ylrg,ylpk,
       zrtb_d,zrix_d,zrmi_d,zrrg_d,zrpk_d,
       zltb_d,zlix_d,zlmi_d,zlrg_d,zlpk_d,deg,1,deg+1);

   sumrtb = 0.0; sumrix = 0.0; sumrmi = 0.0; sumrrg = 0.0; sumrpk = 0.0;
   sumltb = 0.0; sumlix = 0.0; sumlmi = 0.0; sumlrg = 0.0; sumlpk = 0.0;

   for(int k=0; k<=deg; k++)
      daf_inc(&sumrtb,&sumrix,&sumrmi,&sumrrg,&sumrpk,
              &sumltb,&sumlix,&sumlmi,&sumlrg,&sumlpk,
              zrtb_d[k],zrix_d[k],zrmi_d[k],zrrg_d[k],zrpk_d[k],
              zltb_d[k],zlix_d[k],zlmi_d[k],zlrg_d[k],zlpk_d[k]);

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
}

void test_complex_exponential ( int deg )
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
   double tmprtb,tmprix,tmprmi,tmprrg,tmprpk;
   double tmpltb,tmplix,tmplmi,tmplrg,tmplpk;

   random_deca_double
      (&rndrertb,&rndrerix,&rndrermi,&rndrerrg,&rndrerpk,
       &rndreltb,&rndrelix,&rndrelmi,&rndrelrg,&rndrelpk);  // cos(a)

   daf_sqr(rndrertb,rndrerix,rndrermi,rndrerrg,rndrerpk,
           rndreltb,rndrelix,rndrelmi,rndrelrg,rndrelpk,
           &tmprtb,&tmprix,&tmprmi,&tmprrg,&tmprpk,
           &tmpltb,&tmplix,&tmplmi,&tmplrg,&tmplpk);        // cos^2(a)
   daf_minus(&tmprtb,&tmprix,&tmprmi,&tmprrg,&tmprpk,
             &tmpltb,&tmplix,&tmplmi,&tmplrg,&tmplpk);      // -cos^2(a)
   daf_inc_d(&tmprtb,&tmprix,&tmprmi,&tmprrg,&tmprpk,
             &tmpltb,&tmplix,&tmplmi,&tmplrg,&tmplpk,1.0);  // 1-cos^2(a)
   daf_sqrt(tmprtb,tmprix,tmprmi,tmprrg,tmprpk,
            tmpltb,tmplix,tmplmi,tmplrg,tmplpk,
            &rndimrtb,&rndimrix,&rndimrmi,&rndimrrg,&rndimrpk,
            &rndimltb,&rndimlix,&rndimlmi,&rndimlrg,&rndimlpk); // sin is sqrt

   xrertb[0] = 1.0; xrerix[0] = 0.0; xrermi[0] = 0.0;
   xrerrg[0] = 0.0; xrerpk[0] = 0.0;
   xreltb[0] = 0.0; xrelix[0] = 0.0; xrelmi[0] = 0.0;
   xrelrg[0] = 0.0; xrelpk[0] = 0.0;
   yrertb[0] = 1.0; yrerix[0] = 0.0; yrermi[0] = 0.0;
   yrerrg[0] = 0.0; yrerpk[0] = 0.0;
   yreltb[0] = 0.0; yrelix[0] = 0.0; yrelmi[0] = 0.0;
   yrelrg[0] = 0.0; yrelpk[0] = 0.0;
   ximrtb[0] = 0.0; ximrix[0] = 0.0; ximrmi[0] = 0.0;
   ximrrg[0] = 0.0; ximrpk[0] = 0.0;
   ximltb[0] = 0.0; ximlix[0] = 0.0; ximlmi[0] = 0.0;
   ximlrg[0] = 0.0; ximlpk[0] = 0.0;
   yimrtb[0] = 0.0; yimrix[0] = 0.0; yimrmi[0] = 0.0;
   ximrrg[0] = 0.0; yimrpk[0] = 0.0;
   yimltb[0] = 0.0; yimlix[0] = 0.0; yimlmi[0] = 0.0;
   ximlrg[0] = 0.0; yimlpk[0] = 0.0;
   xrertb[1] = rndrertb; xrerix[1] = rndrerix; xrermi[1] = rndrermi;
   xrerrg[1] = rndrerrg; xrerpk[1] = rndrerpk;
   xreltb[1] = rndreltb; xrelix[1] = rndrelix; xrelmi[1] = rndrelmi;
   xrelrg[1] = rndrelrg; xrelpk[1] = rndrelpk;
   ximrtb[1] = rndimrtb; ximrix[1] = rndimrix; ximrmi[1] = rndimrmi;
   ximrrg[1] = rndimrrg; ximrpk[1] = rndimrpk;
   ximltb[1] = rndimltb; ximlix[1] = rndimlix; ximlmi[1] = rndimlmi;
   ximlrg[1] = rndimlrg; ximlpk[1] = rndimlpk;
   yrertb[1] = -rndrertb; yrerix[1] = -rndrerix; yrermi[1] = -rndrermi;
   yrerrg[1] = -rndrerrg; yrerpk[1] = -rndrerpk;
   yreltb[1] = -rndreltb; yrelix[1] = -rndrelix; yrelmi[1] = -rndrelmi;
   yrelrg[1] = -rndrelrg; yrelpk[1] = -rndrelpk;
   yimrtb[1] = -rndimrtb; yimrix[1] = -rndimrix; yimrmi[1] = -rndimrmi;
   yimrrg[1] = -rndimrrg; yimrpk[1] = -rndimrpk;
   yimltb[1] = -rndimltb; yimlix[1] = -rndimlix; yimlmi[1] = -rndimlmi;
   yimlrg[1] = -rndimlrg; yimlpk[1] = -rndimlpk;

   for(int k=2; k<=deg; k++)
   {
      // xre[k] = (xre[k-1]*cr - xim[k-1]*sr)/k;
      daf_mul(xrertb[k-1],xrerix[k-1],xrermi[k-1],xrerrg[k-1],xrerpk[k-1],
              xreltb[k-1],xrelix[k-1],xrelmi[k-1],xrelrg[k-1],xrelpk[k-1],
              rndrertb,rndrerix,rndrermi,rndrerrg,rndrerpk,
              rndreltb,rndrelix,rndrelmi,rndrelrg,rndrelpk,
              &xrertb[k],&xrerix[k],&xrermi[k],&xrerrg[k],&xrerpk[k],
              &xreltb[k],&xrelix[k],&xrelmi[k],&xrelrg[k],&xrelpk[k]);
      daf_mul(ximrtb[k-1],ximrix[k-1],ximrmi[k-1],ximrrg[k-1],ximrpk[k-1],
              ximltb[k-1],ximlix[k-1],ximlmi[k-1],ximlrg[k-1],ximlpk[k-1],
              rndimrtb,rndimrix,rndimrmi,rndimrrg,rndimrpk,
              rndimltb,rndimlix,rndimlmi,rndimlrg,rndimlpk,
              &tmprtb,&tmprix,&tmprmi,&tmprrg,&tmprpk,
              &tmpltb,&tmplix,&tmplmi,&tmplrg,&tmplpk);
      daf_minus(&tmprtb,&tmprix,&tmprmi,&tmprrg,&tmprpk,
                &tmpltb,&tmplix,&tmplmi,&tmplrg,&tmplpk);
      daf_inc(&xrertb[k],&xrerix[k],&xrermi[k],&xrerrg[k],&xrerpk[k],
              &xreltb[k],&xrelix[k],&xrelmi[k],&xrelrg[k],&xrelpk[k],
              tmprtb,tmprix,tmprmi,tmprrg,tmprpk,
              tmpltb,tmplix,tmplmi,tmplrg,tmplpk);
      tmprtb = (double) k;
                    tmprix = 0.0; tmprmi = 0.0; tmprrg = 0.0; tmprpk = 0.0;
      tmpltb = 0.0; tmplix = 0.0; tmplmi = 0.0; tmplrg = 0.0; tmplpk = 0.0;
      daf_div(xrertb[k],xrerix[k],xrermi[k],xrerrg[k],xrerpk[k],
              xreltb[k],xrelix[k],xrelmi[k],xrelrg[k],xrelpk[k],
              tmprtb,tmprix,tmprmi,tmprrg,tmprpk,
              tmpltb,tmplix,tmplmi,tmplrg,tmplpk,
              &xrertb[k],&xrerix[k],&xrermi[k],&xrerrg[k],&xrerpk[k],
              &xreltb[k],&xrelix[k],&xrelmi[k],&xrelrg[k],&xrelpk[k]);
      // xim[k] = (xre[k-1]*sr + xim[k-1]*cr)/k;
      daf_mul(xrertb[k-1],xrerix[k-1],xrermi[k-1],xrerrg[k-1],xrerpk[k-1],
              xreltb[k-1],xrelix[k-1],xrelmi[k-1],xrelrg[k-1],xrelpk[k-1],
              rndimrtb,rndimrix,rndimrmi,rndimrrg,rndimrpk,
              rndimltb,rndimlix,rndimlmi,rndimlrg,rndimlpk,
              &ximrtb[k],&ximrix[k],&ximrmi[k],&ximrrg[k],&ximrpk[k],
              &ximltb[k],&ximlix[k],&ximlmi[k],&ximlrg[k],&ximlpk[k]);
      daf_mul(ximrtb[k-1],ximrix[k-1],ximrmi[k-1],ximrrg[k-1],ximrpk[k-1],
              ximltb[k-1],ximlix[k-1],ximlmi[k-1],ximlrg[k-1],ximlpk[k-1],
              rndrertb,rndrerix,rndrermi,rndrerrg,rndrerpk,
              rndreltb,rndrelix,rndrelmi,rndrelrg,rndrelpk,
              &tmprtb,&tmprix,&tmprmi,&tmprrg,&tmprpk,
              &tmpltb,&tmplix,&tmplmi,&tmplrg,&tmplpk);
      daf_inc(&ximrtb[k],&ximrix[k],&ximrmi[k],&ximrrg[k],&ximrpk[k],
              &ximltb[k],&ximlix[k],&ximlmi[k],&ximlrg[k],&ximlpk[k],
              tmprtb,tmprix,tmprmi,tmprrg,tmprpk,
              tmpltb,tmplix,tmplmi,tmplrg,tmplpk);
      tmprtb = (double) k;
                    tmprix = 0.0; tmprmi = 0.0; tmprrg = 0.0; tmprpk = 0.0;
      tmpltb = 0.0; tmplix = 0.0; tmplmi = 0.0; tmplrg = 0.0; tmplpk = 0.0;
      daf_div(ximrtb[k],ximrix[k],ximrmi[k],ximrrg[k],ximrpk[k],
              ximltb[k],ximlix[k],ximlmi[k],ximlrg[k],ximlpk[k],
              tmprtb,tmprix,tmprmi,tmprrg,tmprpk,
              tmpltb,tmplix,tmplmi,tmplrg,tmplpk,
              &ximrtb[k],&ximrix[k],&ximrmi[k],&ximrrg[k],&ximrpk[k],
              &ximltb[k],&ximlix[k],&ximlmi[k],&ximlrg[k],&ximlpk[k]);
      // yre[k] = (yre[k-1]*(-cr) - yim[k-1]*(-sr))/k;
      daf_mul(yrertb[k-1],yrerix[k-1],yrermi[k-1],yrerrg[k-1],yrerpk[k-1],
              yreltb[k-1],yrelix[k-1],yrelmi[k-1],yrelrg[k-1],yrelpk[k-1],
              -rndrertb,-rndrerix,-rndrermi,-rndrerrg,-rndrerpk,
              -rndreltb,-rndrelix,-rndrelmi,-rndrelrg,-rndrelpk,
              &yrertb[k],&yrerix[k],&yrermi[k],&yrerrg[k],&yrerpk[k],
              &yreltb[k],&yrelix[k],&yrelmi[k],&yrelrg[k],&yrelpk[k]);
      daf_mul(yimrtb[k-1],yimrix[k-1],yimrmi[k-1],yimrrg[k-1],yimrpk[k-1],
              yimltb[k-1],yimlix[k-1],yimlmi[k-1],yimlrg[k-1],yimlpk[k-1],
              -rndimrtb,-rndimrix,-rndimrmi,-rndimrrg,-rndimrpk,
              -rndimltb,-rndimlix,-rndimlmi,-rndimlrg,-rndimlpk,
              &tmprtb,&tmprix,&tmprmi,&tmprrg,&tmprpk,
              &tmpltb,&tmplix,&tmplmi,&tmplrg,&tmplpk);
      daf_minus(&tmprtb,&tmprix,&tmprmi,&tmprrg,&tmprpk,
                &tmpltb,&tmplix,&tmplmi,&tmplrg,&tmplpk);
      daf_inc(&yrertb[k],&yrerix[k],&yrermi[k],&yrerrg[k],&yrerpk[k],
              &yreltb[k],&yrelix[k],&yrelmi[k],&yrelrg[k],&yrelpk[k],
              tmprtb,tmprix,tmprmi,tmprrg,tmprpk,
              tmpltb,tmplix,tmplmi,tmplrg,tmplpk);
      tmprtb = (double) k;
                    tmprix = 0.0; tmprmi = 0.0; tmprrg = 0.0; tmprpk = 0.0;
      tmpltb = 0.0; tmplix = 0.0; tmplmi = 0.0; tmplrg = 0.0; tmplpk = 0.0;
      daf_div(yrertb[k],yrerix[k],yrermi[k],yrerrg[k],yrerpk[k],
              yreltb[k],yrelix[k],yrelmi[k],yrelrg[k],yrelpk[k],
              tmprtb,tmprix,tmprmi,tmprrg,tmprpk,
              tmpltb,tmplix,tmplmi,tmplrg,tmplpk,
              &yrertb[k],&yrerix[k],&yrermi[k],&yrerrg[k],&yrerpk[k],
              &yreltb[k],&yrelix[k],&yrelmi[k],&yrelrg[k],&yrelpk[k]);
      // yim[k] = (yre[k-1]*(-sr) + yim[k-1]*(-cr))/k;
      daf_mul(yrertb[k-1],yrerix[k-1],yrermi[k-1],yrerrg[k-1],yrerpk[k-1],
              yreltb[k-1],yrelix[k-1],yrelmi[k-1],yrelrg[k-1],yrelpk[k-1],
              -rndimrtb,-rndimrix,-rndimrmi,-rndimrrg,-rndimrpk,
              -rndimltb,-rndimlix,-rndimlmi,-rndimlrg,-rndimlpk,
              &yimrtb[k],&yimrix[k],&yimrmi[k],&yimrrg[k],&yimrpk[k],
              &yimltb[k],&yimlix[k],&yimlmi[k],&yimlrg[k],&yimlpk[k]);
      daf_mul(yimrtb[k-1],yimrix[k-1],yimrmi[k-1],yimrrg[k-1],yimrpk[k-1],
              yimltb[k-1],yimlix[k-1],yimlmi[k-1],yimlrg[k-1],yimlpk[k-1],
              -rndrertb,-rndrerix,-rndrermi,-rndrerrg,-rndrerpk,
              -rndreltb,-rndrelix,-rndrelmi,-rndrelrg,-rndrelpk,
              &tmprtb,&tmprix,&tmprmi,&tmprrg,&tmprpk,
              &tmpltb,&tmplix,&tmplmi,&tmplrg,&tmplpk);
      daf_inc(&yimrtb[k],&yimrix[k],&yimrmi[k],&yimrrg[k],&yimrpk[k],
              &yimltb[k],&yimlix[k],&yimlmi[k],&yimlrg[k],&yimlpk[k],
              tmprtb,tmprix,tmprmi,tmprrg,tmprpk,
              tmpltb,tmplix,tmplmi,tmplrg,tmplpk);
      tmprtb = (double) k;
                    tmprix = 0.0; tmprmi = 0.0; tmprrg = 0.0; tmprpk = 0.0;
      tmpltb = 0.0; tmplix = 0.0; tmplmi = 0.0; tmplrg = 0.0; tmplpk = 0.0;
      daf_div(yimrtb[k],yimrix[k],yimrmi[k],yimrrg[k],yimrpk[k],
              yimltb[k],yimlix[k],yimlmi[k],yimlrg[k],yimlpk[k],
              tmprtb,tmprix,tmprmi,tmprrg,tmprpk,
              tmpltb,tmplix,tmplmi,tmplrg,tmplpk,
              &yimrtb[k],&yimrix[k],&yimrmi[k],&yimrrg[k],&yimrpk[k],
              &yimltb[k],&yimlix[k],&yimlmi[k],&yimlrg[k],&yimlpk[k]);
   }

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

   double sumrertb = 0.0;
   double sumrerix = 0.0;
   double sumrermi = 0.0;
   double sumrerrg = 0.0;
   double sumrerpk = 0.0;
   double sumreltb = 0.0;
   double sumrelix = 0.0;
   double sumrelmi = 0.0;
   double sumrelrg = 0.0;
   double sumrelpk = 0.0;
   double sumimrtb = 0.0;
   double sumimrix = 0.0;
   double sumimrmi = 0.0;
   double sumimrrg = 0.0;
   double sumimrpk = 0.0;
   double sumimltb = 0.0;
   double sumimlix = 0.0;
   double sumimlmi = 0.0;
   double sumimlrg = 0.0;
   double sumimlpk = 0.0;

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

   GPU_cmplx10_product
      (xrertb,xrerix,xrermi,xrerrg,xrerpk,xreltb,xrelix,xrelmi,xrelrg,xrelpk,
       ximrtb,ximrix,ximrmi,ximrrg,ximrpk,ximltb,ximlix,ximlmi,ximlrg,ximlpk,
       yrertb,yrerix,yrermi,yrerrg,yrerpk,yreltb,yrelix,yrelmi,yrelrg,yrelpk,
       yimrtb,yimrix,yimrmi,yimrrg,yimrpk,yimltb,yimlix,yimlmi,yimlrg,yimlpk,
       zrertb_d,zrerix_d,zrermi_d,zrerrg_d,zrerpk_d,
       zreltb_d,zrelix_d,zrelmi_d,zrelrg_d,zrelpk_d,
       zimrtb_d,zimrix_d,zimrmi_d,zimrrg_d,zimrpk_d,
       zimltb_d,zimlix_d,zimlmi_d,zimlrg_d,zimlpk_d,deg,1,deg+1);

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
              zrertb_d[k],zrerix_d[k],zrermi_d[k],zrerrg_d[k],zrerpk_d[k],
              zreltb_d[k],zrelix_d[k],zrelmi_d[k],zrelrg_d[k],zrelpk_d[k]);
      daf_inc(&sumimrtb,&sumimrix,&sumimrmi,&sumimrrg,&sumimrpk,
              &sumimltb,&sumimlix,&sumimlmi,&sumimlrg,&sumimlpk,
              zimrtb_d[k],zimrix_d[k],zimrmi_d[k],zimrrg_d[k],zimrpk_d[k],
              zimltb_d[k],zimlix_d[k],zimlmi_d[k],zimlrg_d[k],zimlpk_d[k]);
   }
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
}

