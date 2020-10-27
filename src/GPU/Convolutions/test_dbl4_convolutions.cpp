/* Tests the product of two series in quad double precision. */

#include <iostream>
#include <iomanip>
#include "quad_double_functions.h"
#include "random4_vectors.h"
#include "dbl4_convolutions_host.h"

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
   double* xhihi = new double[deg+1];
   double* xlohi = new double[deg+1];
   double* xhilo = new double[deg+1];
   double* xlolo = new double[deg+1];
   double* yhihi = new double[deg+1];
   double* ylohi = new double[deg+1];
   double* yhilo = new double[deg+1];
   double* ylolo = new double[deg+1];
   double* zhihi = new double[deg+1];
   double* zlohi = new double[deg+1];
   double* zhilo = new double[deg+1];
   double* zlolo = new double[deg+1];

   for(int k=0; k<=deg; k++)
   {
      xhihi[k] = 1.0; xlohi[k] = 0.0;
      xhilo[k] = 0.0; xlolo[k] = 0.0;
      yhihi[k] = 0.0; ylohi[k] = 0.0;
      yhilo[k] = 0.0; ylolo[k] = 0.0;
   }
   yhihi[0] = 1.0; yhihi[1] = -1.0;

   CPU_dbl4_product(deg,xhihi,xlohi,xhilo,xlolo,yhihi,ylohi,yhilo,ylolo,
                        zhihi,zlohi,zhilo,zlolo);

   cout << "Series of 1/(1-x) multiplied with 1-x : " << endl;

   for(int k=0; k<=deg; k++)
   {
      cout << "zhihi[" << k << "] : " << zhihi[k];
      cout << "  zlohi[" << k << "] : " << zlohi[k];
      cout << "  zhilo[" << k << "] : " << zhilo[k];
      cout << "  zlolo[" << k << "] : " << zlolo[k] << endl;
   }
}

void test_complex ( int deg )
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
   double* zrehihi = new double[deg+1];
   double* zrelohi = new double[deg+1];
   double* zrehilo = new double[deg+1];
   double* zrelolo = new double[deg+1];
   double* zimhihi = new double[deg+1];
   double* zimlohi = new double[deg+1];
   double* zimhilo = new double[deg+1];
   double* zimlolo = new double[deg+1];

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
           zrehihi,zrelohi,zrehilo,zimlolo,zimhihi,zimlohi,zimhilo,zimlolo);

   cout << "Series of 1/(1-x) multiplied with 1-x : " << endl;

   for(int k=0; k<=deg; k++)
   {
      cout << "zrehihi[" << k << "] : " << zrehihi[k];
      cout << "  zrelohi[" << k << "] : " << zrelohi[k];
      cout << "  zrehilo[" << k << "] : " << zrehilo[k];
      cout << "  zrelolo[" << k << "] : " << zrelolo[k] << endl;
      cout << "zimhihi[" << k << "] : " << zimhihi[k];
      cout << "  zimlohi[" << k << "] : " << zimlohi[k];
      cout << "  zimhilo[" << k << "] : " << zimhilo[k];
      cout << "  zimlolo[" << k << "] : " << zimlolo[k] << endl;
   }
}

void test_real_exponential ( int deg )
{
   double *xhihi = new double[deg+1];
   double *xlohi = new double[deg+1];
   double *xhilo = new double[deg+1];
   double *xlolo = new double[deg+1];
   double *yhihi = new double[deg+1];
   double *ylohi = new double[deg+1];
   double *yhilo = new double[deg+1];
   double *ylolo = new double[deg+1];
   double *zhihi = new double[deg+1];
   double *zlohi = new double[deg+1];
   double *zhilo = new double[deg+1];
   double *zlolo = new double[deg+1];
   double rhihi,rlohi,rhilo,rlolo;
   double fhihi,flohi,fhilo,flolo;

   random_quad_double(&rhihi,&rlohi,&rhilo,&rlolo);

   xhihi[0] = 1.0; xlohi[0] = 0.0; xhilo[0] = 0.0; xlolo[0] = 0.0;
   yhihi[0] = 1.0; ylohi[0] = 0.0; yhilo[0] = 0.0; ylolo[0] = 0.0;
   xhihi[1] = rhihi; xlohi[1] = rlohi;
   xhilo[1] = rhilo; xlolo[1] = rlolo;
   yhihi[1] = -rhihi; ylohi[1] = -rlohi;
   yhilo[1] = -rhilo; ylolo[1] = -rlolo;

   for(int k=2; k<=deg; k++)
   {
      qdf_mul(xhihi[k-1],xlohi[k-1],xhilo[k-1],xlolo[k-1],
              rhihi,rlohi,rhilo,rlolo,
              &xhihi[k],&xlohi[k],&xhilo[k],&xlolo[k]);    // x[k] = x[k-1]*r;
      qdf_mul(yhihi[k-1],ylohi[k-1],yhilo[k-1],ylolo[k-1],
              -rhihi,-rlohi,-rhilo,-rlolo,
              &yhihi[k],&ylohi[k],&yhilo[k],&ylolo[k]); 
      // y[k] = y[k-1]*(-r);
      fhihi = (double) k; flohi = 0.0; fhilo = 0.0; flolo = 0.0;
      qdf_div(xhihi[k],xlohi[k],xhilo[k],xlolo[k],fhihi,flohi,fhilo,flolo,
              &xhihi[k],&xlohi[k],&xhilo[k],&xlolo[k]);
      qdf_div(yhihi[k],ylohi[k],yhilo[k],ylolo[k],fhihi,flohi,fhilo,flolo,
              &yhihi[k],&ylohi[k],&yhilo[k],&ylolo[k]);
   }

   CPU_dbl4_product(deg,xhihi,xlohi,xhilo,xlolo,yhihi,ylohi,yhilo,ylolo,
                        zhihi,zlohi,zhilo,zlolo);

   cout << scientific << setprecision(16);

   cout << "Product of series of exp(x) with series of exp(-x)," << endl;
   cout << "  for xhihi = " << rhihi;
   cout << "      xlohi = " << rlohi << endl;
   cout << "      xhilo = " << rhilo;
   cout << "  and xlolo = " << rlolo << endl;

   double sumhihi = 0.0;
   double sumlohi = 0.0;
   double sumhilo = 0.0;
   double sumlolo = 0.0;

   for(int k=0; k<=deg; k++)
      qdf_inc(&sumhihi,&sumlohi,&sumhilo,&sumlolo,
              zhihi[k],zlohi[k],zhilo[k],zlolo[k]);

   cout << "Summation of all coefficients in the product ..." << endl;
   cout << "    highest part of sum : " << sumhihi << endl;
   cout << "2nd highest part of sum : " << sumlohi << endl;
   cout << " 2nd lowest part of sum : " << sumhilo << endl;
   cout << "     lowest part of sum : " << sumlolo << endl;
}

void test_complex_exponential ( int deg )
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
   double* zrehihi = new double[deg+1];
   double* zrelohi = new double[deg+1];
   double* zrehilo = new double[deg+1];
   double* zrelolo = new double[deg+1];
   double* zimhihi = new double[deg+1];
   double* zimlohi = new double[deg+1];
   double* zimhilo = new double[deg+1];
   double* zimlolo = new double[deg+1];
   double rndrehihi,rndrelohi,rndrehilo,rndrelolo;
   double rndimhihi,rndimlohi,rndimhilo,rndimlolo;
   double tmphihi,tmplohi,tmphilo,tmplolo;

   random_quad_double(&rndrehihi,&rndrelohi,
                      &rndrehilo,&rndrelolo);              // cos(a)

   qdf_sqr(rndrehihi,rndrelohi,rndrehilo,rndrelolo,
           &tmphihi,&tmplohi,&tmphilo,&tmplolo);           // cos^2(a)
   qdf_minus(&tmphihi,&tmplohi,&tmphilo,&tmplolo);         // -cos^2(a)
   qdf_inc_d(&tmphihi,&tmplohi,&tmphilo,&tmplolo,1.0);     // 1-cos^2(a)
   qdf_sqrt(tmphihi,tmplohi,tmphilo,tmplolo,
            &rndimhihi,&rndimlohi,&rndimhilo,&rndimlolo);  // sin is sqrt

   xrehihi[0] = 1.0; xrelohi[0] = 0.0; xrehilo[0] = 0.0; xrelolo[0] = 0.0;
   yrehihi[0] = 1.0; yrelohi[0] = 0.0; yrehilo[0] = 0.0; yrelolo[0] = 0.0;
   ximhihi[0] = 0.0; ximlohi[0] = 0.0; ximhilo[0] = 0.0; ximlolo[0] = 0.0;
   yimhihi[0] = 0.0; yimlohi[0] = 0.0; yimhilo[0] = 0.0; yimlolo[0] = 0.0;
   xrehihi[1] = rndrehihi; xrelohi[1] = rndrelohi;
   xrehilo[1] = rndrehilo; xrelolo[1] = rndrelolo;
   ximhihi[1] = rndimhihi; ximlohi[1] = rndimlohi;
   ximhilo[1] = rndimhilo; ximlolo[1] = rndimlolo;
   yrehihi[1] = -rndrehihi; yrelohi[1] = -rndrelohi;
   yrehilo[1] = -rndrehilo; yrelolo[1] = -rndrelolo;
   yimhihi[1] = -rndimhihi; yimlohi[1] = -rndimlohi;
   yimhilo[1] = -rndimhilo; yimlolo[1] = -rndimlolo;

   for(int k=2; k<=deg; k++)
   {
      // xre[k] = (xre[k-1]*cr - xim[k-1]*sr)/k;
      qdf_mul(xrehihi[k-1],xrelohi[k-1],xrehilo[k-1],xrelolo[k-1],
              rndrehihi,rndrelohi,rndrehilo,rndrelolo,
              &xrehihi[k],&xrelohi[k],&xrehilo[k],&xrelolo[k]);
      qdf_mul(ximhihi[k-1],ximlohi[k-1],ximhilo[k-1],ximlolo[k-1],
              rndimhihi,rndimlohi,rndimhilo,rndimlolo,
              &tmphihi,&tmplohi,&tmphilo,&tmplolo);
      qdf_minus(&tmphihi,&tmplohi,&tmphilo,&tmplolo);
      qdf_inc(&xrehihi[k],&xrelohi[k],&xrehilo[k],&xrelolo[k],
              tmphihi,tmplohi,tmphilo,tmplolo);
      tmphihi = (double) k; tmplohi = 0.0; tmphilo = 0.0; tmplolo = 0.0;
      qdf_div(xrehihi[k],xrelohi[k],xrehilo[k],xrelolo[k],
              tmphihi,tmplohi,tmphilo,tmplolo,
              &xrehihi[k],&xrelohi[k],&xrehilo[k],&xrelolo[k]);
      // xim[k] = (xre[k-1]*sr + xim[k-1]*cr)/k;
      qdf_mul(xrehihi[k-1],xrelohi[k-1],xrehilo[k-1],xrelolo[k-1],
              rndimhihi,rndimlohi,rndimhilo,rndimlolo,
              &ximhihi[k],&ximlohi[k],&ximhilo[k],&ximlolo[k]);
      qdf_mul(ximhihi[k-1],ximlohi[k-1],ximhilo[k-1],ximlolo[k-1],
              rndrehihi,rndrelohi,rndrehilo,rndrelolo,
              &tmphihi,&tmplohi,&tmphilo,&tmplolo);
      qdf_inc(&ximhihi[k],&ximlohi[k],&ximhilo[k],&ximlolo[k],
              tmphihi,tmplohi,tmphilo,tmplolo);
      tmphihi = (double) k; tmplohi = 0.0; tmphilo = 0.0; tmplolo = 0.0;
      qdf_div(ximhihi[k],ximlohi[k],ximhilo[k],ximlolo[k],
              tmphihi,tmplohi,tmphilo,tmplolo,
              &ximhihi[k],&ximlohi[k],&ximhilo[k],&ximlolo[k]);
      // yre[k] = (yre[k-1]*(-cr) - yim[k-1]*(-sr))/k;
      qdf_mul(yrehihi[k-1],yrelohi[k-1],yrehilo[k-1],yrelolo[k-1],
              -rndrehihi,-rndrelohi,-rndrehilo,-rndrelolo,
              &yrehihi[k],&yrelohi[k],&yrehilo[k],&yrelolo[k]);
      qdf_mul(yimhihi[k-1],yimlohi[k-1],yimhilo[k-1],yimlolo[k-1],
              -rndimhihi,-rndimlohi,-rndimhilo,-rndimlolo,
              &tmphihi,&tmplohi,&tmphilo,&tmplolo);
      qdf_minus(&tmphihi,&tmplohi,&tmphilo,&tmplolo);
      qdf_inc(&yrehihi[k],&yrelohi[k],&yrehilo[k],&yrelolo[k],
              tmphihi,tmplohi,tmphilo,tmplolo);
      tmphihi = (double) k; tmplohi = 0.0; tmphilo = 0.0; tmplolo = 0.0;
      qdf_div(yrehihi[k],yrelohi[k],yrehilo[k],yrelolo[k],
              tmphihi,tmplohi,tmphilo,tmplolo,
              &yrehihi[k],&yrelohi[k],&yrehilo[k],&yrelolo[k]);
      // yim[k] = (yre[k-1]*(-sr) + yim[k-1]*(-cr))/k;
      qdf_mul(yrehihi[k-1],yrelohi[k-1],yrehilo[k-1],yrelolo[k-1],
              -rndimhihi,-rndimlohi,-rndimhilo,-rndimlolo,
              &yimhihi[k],&yimlohi[k],&yimhilo[k],&yimlolo[k]);
      qdf_mul(yimhihi[k-1],yimlohi[k-1],yimhilo[k-1],yimlolo[k-1],
              -rndrehihi,-rndrelohi,-rndrehilo,-rndrelolo,
              &tmphihi,&tmplohi,&tmphilo,&tmplolo);
      qdf_inc(&yimhihi[k],&yimlohi[k],&yimhilo[k],&yimlolo[k],
              tmphihi,tmplohi,tmphilo,tmplolo);
      tmphihi = (double) k; tmplohi = 0.0; tmphilo = 0.0; tmplolo = 0.0;
      qdf_div(yimhihi[k],yimlohi[k],yimhilo[k],yimlolo[k],
              tmphihi,tmplohi,tmphilo,tmplolo,
              &yimhihi[k],&yimlohi[k],&yimhilo[k],&yimlolo[k]);
   }

   CPU_cmplx4_product
      (deg,xrehihi,xrelohi,xrehilo,xrelolo,ximhihi,ximlohi,ximhilo,ximlolo,
           yrehihi,yrelohi,yrehilo,yrelolo,yimhihi,yimlohi,yimhilo,yimlolo,
           zrehihi,zrelohi,zrehilo,zrelolo,zimhihi,zimlohi,zimhilo,zimlolo);

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

   double sumrehihi = 0.0;
   double sumrelohi = 0.0;
   double sumrehilo = 0.0;
   double sumrelolo = 0.0;
   double sumimhihi = 0.0;
   double sumimlohi = 0.0;
   double sumimhilo = 0.0;
   double sumimlolo = 0.0;

   for(int k=0; k<=deg; k++) 
   {
      qdf_inc(&sumrehihi,&sumrelohi,&sumrehilo,&sumrelolo,
              zrehihi[k],zrelohi[k],zrehilo[k],zrelolo[k]);
      qdf_inc(&sumimhihi,&sumimlohi,&sumimhilo,&sumimlolo,
              zimhihi[k],zimlohi[k],zimhilo[k],zimlolo[k]);
   }
   cout << "Summation of all coefficients of the product ..." << endl;
   cout << "  sumrehihi : " << sumrehihi;
   cout << "  sumrelohi : " << sumrelohi << endl;
   cout << "  sumrehilo : " << sumrehilo;
   cout << "  sumrelolo : " << sumrelolo << endl;
   cout << "  sumimhihi : " << sumimhihi;
   cout << "  sumimloi : " << sumimlohi << endl;
   cout << "  sumimhilo : " << sumimhilo;
   cout << "  sumimlolo : " << sumimlolo << endl;
}
