/* Tests the product of two series in triple double precision. */

#include <iostream>
#include <iomanip>
#include <vector_types.h>
#include "triple_double_functions.h"
#include "random3_vectors.h"
#include "dbl3_convolutions_host.h"
#include "dbl3_convolutions_kernels.h"

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

   cout << "Series of 1/(1-x) multiplied with 1-x : " << endl;

   for(int k=0; k<=deg; k++)
   {
      cout << "zhi[" << k << "] : " << zhi_h[k];
      cout << "  zmi[" << k << "] : " << zmi_h[k];
      cout << "  zlo[" << k << "] : " << zlo_h[k] << endl;
   }

   GPU_dbl3_product(xhi,xmi,xlo,yhi,ymi,ylo,zhi_d,zmi_d,zlo_d,deg,1,deg+1);

   cout << "GPU computed product :" << endl;

   for(int k=0; k<=deg; k++)
   {
      cout << "zhi[" << k << "] : " << zhi_d[k];
      cout << "  zmi[" << k << "] : " << zmi_d[k];
      cout << "  zlo[" << k << "] : " << zlo_d[k] << endl;
   }
   cout << endl;
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
      cout << "zrehi[" << k << "] : " << zrehi[k];
      cout << "  zremi[" << k << "] : " << zremi[k];
      cout << "  zrelo[" << k << "] : " << zrelo[k] << endl;
      cout << "zimhi[" << k << "] : " << zimhi[k];
      cout << "  zimmi[" << k << "] : " << zimmi[k];
      cout << "  zimlo[" << k << "] : " << zimlo[k] << endl;
   }
}

void test_real_exponential ( int deg )
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
   double fhi,fmi,flo;

   random_triple_double(&rhi,&rmi,&rlo);

   xhi[0] = 1.0; xmi[0] = 0.0; xlo[0] = 0.0;
   yhi[0] = 1.0; ymi[0] = 0.0; ylo[0] = 0.0;
   xhi[1] = rhi; xmi[1] = rmi; xlo[1] = rlo;
   yhi[1] = -rhi; ymi[1] = -rmi; ylo[1] = -rlo;

   for(int k=2; k<=deg; k++)
   {
      tdf_mul(xhi[k-1],xmi[k-1],xlo[k-1],rhi,rmi,rlo,
              &xhi[k],&xmi[k],&xlo[k]);             // x[k] = x[k-1]*r;
      tdf_mul(yhi[k-1],ymi[k-1],ylo[k-1],-rhi,-rmi,-rlo,
              &yhi[k],&ymi[k],&ylo[k]); 
      // y[k] = y[k-1]*(-r);
      fhi = (double) k; fmi = 0.0; flo = 0.0;
      tdf_div(xhi[k],xmi[k],xlo[k],fhi,fmi,flo,&xhi[k],&xmi[k],&xlo[k]);
      tdf_div(yhi[k],ymi[k],ylo[k],fhi,fmi,flo,&yhi[k],&ymi[k],&ylo[k]);
   }

   CPU_dbl3_product(deg,xhi,xmi,xlo,yhi,ymi,ylo,zhi_h,zmi_h,zlo_h);

   cout << scientific << setprecision(16);

   cout << "Product of series of exp(x) with series of exp(-x)," << endl;
   cout << "  for xhi = " << rhi << endl;
   cout << "      xmi = " << rmi << endl;
   cout << "  and xlo = " << rlo << endl;

   double sumhi = 0.0;
   double summi = 0.0;
   double sumlo = 0.0;

   for(int k=0; k<=deg; k++)
      tdf_inc(&sumhi,&summi,&sumlo,zhi_h[k],zmi_h[k],zlo_h[k]);

   cout << "Summation of all coefficients in the product ..." << endl;
   cout << "  high part of sum : " << sumhi << endl;
   cout << "middle part of sum : " << summi << endl;
   cout << "   low part of sum : " << sumlo << endl;

   GPU_dbl3_product(xhi,xmi,xlo,yhi,ymi,ylo,zhi_d,zmi_d,zlo_d,deg,1,deg+1);

   sumhi = 0.0; sumlo = 0.0;

   for(int k=0; k<=deg; k++)
      tdf_inc(&sumhi,&summi,&sumlo,zhi_d[k],zmi_d[k],zlo_d[k]);

   cout << "Summation of all coefficients in the GPU computed product ..."
        << endl;
   cout << "  high part of sum : " << sumhi << endl;
   cout << "middle part of sum : " << summi << endl;
   cout << "   low part of sum : " << sumlo << endl;
}

void test_complex_exponential ( int deg )
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

   random_triple_double(&rndrehi,&rndremi,&rndrelo);       // cos(a)

   tdf_sqr(rndrehi,rndremi,rndrelo,&tmphi,&tmpmi,&tmplo);  // cos^2(a)
   tdf_minus(&tmphi,&tmpmi,&tmplo);                        // -cos^2(a)
   tdf_inc_d(&tmphi,&tmpmi,&tmplo,1.0);                    // 1-cos^2(a)
   tdf_sqrt(tmphi,tmpmi,tmplo,&rndimhi,&rndimmi,&rndimlo); // sin is sqrt

   xrehi[0] = 1.0; xremi[0] = 0.0; xrelo[0] = 0.0;
   yrehi[0] = 1.0; yremi[0] = 0.0; yrelo[0] = 0.0;
   ximhi[0] = 0.0; ximmi[0] = 0.0; ximlo[0] = 0.0;
   yimhi[0] = 0.0; yimmi[0] = 0.0; yimlo[0] = 0.0;
   xrehi[1] = rndrehi; xremi[1] = rndremi; xrelo[1] = rndrelo;
   ximhi[1] = rndimhi; ximmi[1] = rndimmi; ximlo[1] = rndimlo;
   yrehi[1] = -rndrehi; yremi[1] = -rndremi; yrelo[1] = -rndrelo;
   yimhi[1] = -rndimhi; yimmi[1] = -rndimmi; yimlo[1] = -rndimlo;

   for(int k=2; k<=deg; k++)
   {
      // xre[k] = (xre[k-1]*cr - xim[k-1]*sr)/k;
      tdf_mul(xrehi[k-1],xremi[k-1],xrelo[k-1],rndrehi,rndremi,rndrelo,
              &xrehi[k],&xremi[k],&xrelo[k]);
      tdf_mul(ximhi[k-1],ximmi[k-1],ximlo[k-1],rndimhi,rndimmi,rndimlo,
              &tmphi,&tmpmi,&tmplo);
      tdf_minus(&tmphi,&tmpmi,&tmplo);
      tdf_inc(&xrehi[k],&xremi[k],&xrelo[k],tmphi,tmpmi,tmplo);
      tmphi = (double) k; tmpmi = 0.0; tmplo = 0.0;
      tdf_div(xrehi[k],xremi[k],xrelo[k],tmphi,tmpmi,tmplo,
              &xrehi[k],&xremi[k],&xrelo[k]);
      // xim[k] = (xre[k-1]*sr + xim[k-1]*cr)/k;
      tdf_mul(xrehi[k-1],xremi[k-1],xrelo[k-1],rndimhi,rndimmi,rndimlo,
              &ximhi[k],&ximmi[k],&ximlo[k]);
      tdf_mul(ximhi[k-1],ximmi[k-1],ximlo[k-1],rndrehi,rndremi,rndrelo,
              &tmphi,&tmpmi,&tmplo);
      tdf_inc(&ximhi[k],&ximmi[k],&ximlo[k],tmphi,tmpmi,tmplo);
      tmphi = (double) k; tmpmi = 0.0; tmplo = 0.0;
      tdf_div(ximhi[k],ximmi[k],ximlo[k],tmphi,tmpmi,tmplo,
              &ximhi[k],&ximmi[k],&ximlo[k]);
      // yre[k] = (yre[k-1]*(-cr) - yim[k-1]*(-sr))/k;
      tdf_mul(yrehi[k-1],yremi[k-1],yrelo[k-1],-rndrehi,-rndremi,-rndrelo,
              &yrehi[k],&yremi[k],&yrelo[k]);
      tdf_mul(yimhi[k-1],yimmi[k-1],yimlo[k-1],-rndimhi,-rndimmi,-rndimlo,
              &tmphi,&tmpmi,&tmplo);
      tdf_minus(&tmphi,&tmpmi,&tmplo);
      tdf_inc(&yrehi[k],&yremi[k],&yrelo[k],tmphi,tmpmi,tmplo);
      tmphi = (double) k; tmpmi = 0.0; tmplo = 0.0;
      tdf_div(yrehi[k],yremi[k],yrelo[k],tmphi,tmpmi,tmplo,
              &yrehi[k],&yremi[k],&yrelo[k]);
      // yim[k] = (yre[k-1]*(-sr) + yim[k-1]*(-cr))/k;
      tdf_mul(yrehi[k-1],yremi[k-1],yrelo[k-1],-rndimhi,-rndimmi,-rndimlo,
              &yimhi[k],&yimmi[k],&yimlo[k]);
      tdf_mul(yimhi[k-1],yimmi[k-1],yimlo[k-1],-rndrehi,-rndremi,-rndrelo,
              &tmphi,&tmpmi,&tmplo);
      tdf_inc(&yimhi[k],&yimmi[k],&yimlo[k],tmphi,tmpmi,tmplo);
      tmphi = (double) k; tmpmi = 0.0; tmplo = 0.0;
      tdf_div(yimhi[k],yimmi[k],yimlo[k],tmphi,tmpmi,tmplo,
              &yimhi[k],&yimmi[k],&yimlo[k]);
   }

   CPU_cmplx3_product(deg,xrehi,xremi,xrelo,ximhi,ximmi,ximlo,
                          yrehi,yremi,yrelo,yimhi,yimmi,yimlo,
                          zrehi_h,zremi_h,zrelo_h,zimhi_h,zimmi_h,zimlo_h);

   cout << scientific << setprecision(16);

   cout << "Product of series of exp(x) with series of exp(-x)," << endl;
   cout << "  for xrehi = " << rndrehi;
   cout << "      xremi = " << rndremi << endl;
   cout << "  and xrelo = " << rndrelo << endl;
   cout << "  for ximhi = " << rndimhi;
   cout << "      ximmi = " << rndimmi << endl;
   cout << "  and ximlo = " << rndimlo << endl;

   double sumrehi = 0.0;
   double sumremi = 0.0;
   double sumrelo = 0.0;
   double sumimhi = 0.0;
   double sumimmi = 0.0;
   double sumimlo = 0.0;

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

   GPU_cmplx3_product
      (xrehi,xremi,xrelo,ximhi,ximmi,ximlo,
       yrehi,yremi,yrelo,yimhi,yimmi,yimlo,
       zrehi_d,zremi_d,zrelo_d,zimhi_d,zimmi_d,zimlo_d,deg,1,deg+1);

   sumrehi = 0.0; sumremi = 0.0; sumrelo = 0.0;
   sumimhi = 0.0; sumimmi = 0.0; sumimlo = 0.0;

   for(int k=0; k<=deg; k++) 
   {
      tdf_inc(&sumrehi,&sumremi,&sumrelo,zrehi_d[k],zremi_d[k],zrelo_d[k]);
      tdf_inc(&sumimhi,&sumimmi,&sumimlo,zimhi_d[k],zimmi_d[k],zimlo_d[k]);
   }
   cout << "Summation of all coefficients of the GPU computed product ..."
        << endl;
   cout << "  sumrehi : " << sumrehi;
   cout << "  sumremi : " << sumremi << endl;
   cout << "  sumrelo : " << sumrelo << endl;
   cout << "  sumimhi : " << sumimhi;
   cout << "  sumimmi : " << sumimmi << endl;
   cout << "  sumimlo : " << sumimlo << endl;
}
