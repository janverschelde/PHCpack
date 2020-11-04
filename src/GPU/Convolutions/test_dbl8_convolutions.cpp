/* Tests the product of two series in octo double precision. */

#include <iostream>
#include <iomanip>
#include <vector_types.h>
#include "octo_double_functions.h"
#include "random8_vectors.h"
#include "random8_series.h"
#include "dbl8_convolutions_host.h"
#include "dbl8_convolutions_kernels.h"

using namespace std;

void test_real ( int deg );
/*
 * DESCRIPTION :
 *   Multiplies the power series of 1/(1-x) with 1+x,
 *   truncated to degree deg, for real coefficients. */

void test_real_random ( int deg );
/*
 * DESCRIPTION :
 *   Multiplies two random power series. */

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
   double* xhihihi = new double[deg+1];
   double* xlohihi = new double[deg+1];
   double* xhilohi = new double[deg+1];
   double* xlolohi = new double[deg+1];
   double* xhihilo = new double[deg+1];
   double* xlohilo = new double[deg+1];
   double* xhilolo = new double[deg+1];
   double* xlololo = new double[deg+1];
   double* yhihihi = new double[deg+1];
   double* ylohihi = new double[deg+1];
   double* yhilohi = new double[deg+1];
   double* ylolohi = new double[deg+1];
   double* yhihilo = new double[deg+1];
   double* ylohilo = new double[deg+1];
   double* yhilolo = new double[deg+1];
   double* ylololo = new double[deg+1];
   double* zhihihi_h = new double[deg+1];
   double* zlohihi_h = new double[deg+1];
   double* zhilohi_h = new double[deg+1];
   double* zlolohi_h = new double[deg+1];
   double* zhihilo_h = new double[deg+1];
   double* zlohilo_h = new double[deg+1];
   double* zhilolo_h = new double[deg+1];
   double* zlololo_h = new double[deg+1];
   double* zhihihi_d = new double[deg+1];
   double* zlohihi_d = new double[deg+1];
   double* zhilohi_d = new double[deg+1];
   double* zlolohi_d = new double[deg+1];
   double* zhihilo_d = new double[deg+1];
   double* zlohilo_d = new double[deg+1];
   double* zhilolo_d = new double[deg+1];
   double* zlololo_d = new double[deg+1];

   for(int k=0; k<=deg; k++)
   {
      xhihihi[k] = 1.0; xlohihi[k] = 0.0;
      xhilohi[k] = 0.0; xlolohi[k] = 0.0;
      xhihilo[k] = 0.0; xlohilo[k] = 0.0;
      xhilolo[k] = 0.0; xlololo[k] = 0.0;
      yhihihi[k] = 0.0; ylohihi[k] = 0.0;
      yhilohi[k] = 0.0; ylolohi[k] = 0.0;
      yhihilo[k] = 0.0; ylohilo[k] = 0.0;
      yhilolo[k] = 0.0; ylololo[k] = 0.0;
   }
   yhihihi[0] = 1.0; yhihihi[1] = -1.0;

   CPU_dbl8_product
      (deg,xhihihi,xlohihi,xhilohi,xlolohi,xhihilo,xlohilo,xhilolo,xlololo,
           yhihihi,ylohihi,yhilohi,ylolohi,yhihilo,ylohilo,yhilolo,ylololo,
           zhihihi_h,zlohihi_h,zhilohi_h,zlolohi_h,
           zhihilo_h,zlohilo_h,zhilolo_h,zlololo_h);

   cout << "Series of 1/(1-x) multiplied with 1-x : " << endl;

   for(int k=0; k<=deg; k++)
   {
      cout << "zhihihi[" << k << "] : " << zhihihi_h[k];
      cout << "  zlohihi[" << k << "] : " << zlohihi_h[k];
      cout << "  zhilohi[" << k << "] : " << zhilohi_h[k];
      cout << "  zlolohi[" << k << "] : " << zlolohi_h[k] << endl;
      cout << "zhihilo[" << k << "] : " << zhihilo_h[k];
      cout << "  zlohilo[" << k << "] : " << zlohilo_h[k];
      cout << "  zhilolo[" << k << "] : " << zhilolo_h[k];
      cout << "  zlololo[" << k << "] : " << zlololo_h[k] << endl;
   }

   GPU_dbl8_product
      (xhihihi,xlohihi,xhilohi,xlolohi,xhihilo,xlohilo,xhilolo,xlololo,
       yhihihi,ylohihi,yhilohi,ylolohi,yhihilo,ylohilo,yhilolo,ylololo,
       zhihihi_d,zlohihi_d,zhilohi_d,zlolohi_d,
       zhihilo_d,zlohilo_d,zhilolo_d,zlololo_d,deg,1,deg+1);

   cout << "GPU computed product :" << endl;

   for(int k=0; k<=deg; k++)
   {
      cout << "zhihihi[" << k << "] : " << zhihihi_d[k];
      cout << "  zlohihi[" << k << "] : " << zlohihi_d[k];
      cout << "  zhilohi[" << k << "] : " << zhilohi_d[k];
      cout << "  zlolohi[" << k << "] : " << zlolohi_d[k] << endl;
      cout << "zhihilo[" << k << "] : " << zhihilo_d[k];
      cout << "  zlohilo[" << k << "] : " << zlohilo_d[k];
      cout << "  zhilolo[" << k << "] : " << zhilolo_d[k];
      cout << "  zlololo[" << k << "] : " << zlololo_d[k] << endl;
   }
}

void test_real_random ( int deg )
{
   double* xhihihi = new double[deg+1];
   double* xlohihi = new double[deg+1];
   double* xhilohi = new double[deg+1];
   double* xlolohi = new double[deg+1];
   double* xhihilo = new double[deg+1];
   double* xlohilo = new double[deg+1];
   double* xhilolo = new double[deg+1];
   double* xlololo = new double[deg+1];
   double* yhihihi = new double[deg+1];
   double* ylohihi = new double[deg+1];
   double* yhilohi = new double[deg+1];
   double* ylolohi = new double[deg+1];
   double* yhihilo = new double[deg+1];
   double* ylohilo = new double[deg+1];
   double* yhilolo = new double[deg+1];
   double* ylololo = new double[deg+1];
   double* zhihihi_h = new double[deg+1];
   double* zlohihi_h = new double[deg+1];
   double* zhilohi_h = new double[deg+1];
   double* zlolohi_h = new double[deg+1];
   double* zhihilo_h = new double[deg+1];
   double* zlohilo_h = new double[deg+1];
   double* zhilolo_h = new double[deg+1];
   double* zlololo_h = new double[deg+1];
   double* zhihihi_d = new double[deg+1];
   double* zlohihi_d = new double[deg+1];
   double* zhilohi_d = new double[deg+1];
   double* zlolohi_d = new double[deg+1];
   double* zhihilo_d = new double[deg+1];
   double* zlohilo_d = new double[deg+1];
   double* zhilolo_d = new double[deg+1];
   double* zlololo_d = new double[deg+1];

   for(int k=0; k<=deg; k++)
   {
      random_octo_double
         (&xhihihi[k],&xlohihi[k],&xhilohi[k],&xlolohi[k],
          &xhihilo[k],&xlohilo[k],&xhilolo[k],&xlololo[k]);
      random_octo_double
         (&yhihihi[k],&ylohihi[k],&yhilohi[k],&ylolohi[k],
          &yhihilo[k],&ylohilo[k],&yhilolo[k],&ylololo[k]);
   }

   CPU_dbl8_product
      (deg,xhihihi,xlohihi,xhilohi,xlolohi,xhihilo,xlohilo,xhilolo,xlololo,
           yhihihi,ylohihi,yhilohi,ylolohi,yhihilo,ylohilo,yhilolo,ylololo,
           zhihihi_h,zlohihi_h,zhilohi_h,zlolohi_h,
           zhihilo_h,zlohilo_h,zhilolo_h,zlololo_h);

   cout << "Series of 1/(1-x) multiplied with 1-x : " << endl;

   cout << scientific << setprecision(16);

   for(int k=0; k<=deg; k++)
   {
      cout << "zhihihi[" << k << "] : " << zhihihi_h[k];
      cout << "  zlohihi[" << k << "] : " << zlohihi_h[k];
      cout << "  zhilohi[" << k << "] : " << zhilohi_h[k];
      cout << "  zlolohi[" << k << "] : " << zlolohi_h[k] << endl;
      cout << "zhihilo[" << k << "] : " << zhihilo_h[k];
      cout << "  zlohilo[" << k << "] : " << zlohilo_h[k];
      cout << "  zhilolo[" << k << "] : " << zhilolo_h[k];
      cout << "  zlololo[" << k << "] : " << zlololo_h[k] << endl;
   }

   GPU_dbl8_product
      (xhihihi,xlohihi,xhilohi,xlolohi,xhihilo,xlohilo,xhilolo,xlololo,
       yhihihi,ylohihi,yhilohi,ylolohi,yhihilo,ylohilo,yhilolo,ylololo,
       zhihihi_d,zlohihi_d,zhilohi_d,zlolohi_d,
       zhihilo_d,zlohilo_d,zhilolo_d,zlololo_d,deg,1,deg+1);

   cout << "GPU computed product :" << endl;

   for(int k=0; k<=deg; k++)
   {
      cout << "zhihihi[" << k << "] : " << zhihihi_d[k];
      cout << "  zlohihi[" << k << "] : " << zlohihi_d[k];
      cout << "  zhilohi[" << k << "] : " << zhilohi_d[k];
      cout << "  zlolohi[" << k << "] : " << zlolohi_d[k] << endl;
      cout << "zhihilo[" << k << "] : " << zhihilo_d[k];
      cout << "  zlohilo[" << k << "] : " << zlohilo_d[k];
      cout << "  zhilolo[" << k << "] : " << zhilolo_d[k];
      cout << "  zlololo[" << k << "] : " << zlololo_d[k] << endl;
   }
}

void test_complex ( int deg )
{
   double* xrehihihi = new double[deg+1];
   double* xrelohihi = new double[deg+1];
   double* xrehilohi = new double[deg+1];
   double* xrelolohi = new double[deg+1];
   double* xrehihilo = new double[deg+1];
   double* xrelohilo = new double[deg+1];
   double* xrehilolo = new double[deg+1];
   double* xrelololo = new double[deg+1];
   double* ximhihihi = new double[deg+1];
   double* ximlohihi = new double[deg+1];
   double* ximhilohi = new double[deg+1];
   double* ximlolohi = new double[deg+1];
   double* ximhihilo = new double[deg+1];
   double* ximlohilo = new double[deg+1];
   double* ximhilolo = new double[deg+1];
   double* ximlololo = new double[deg+1];
   double* yrehihihi = new double[deg+1];
   double* yrelohihi = new double[deg+1];
   double* yrehilohi = new double[deg+1];
   double* yrelolohi = new double[deg+1];
   double* yrehihilo = new double[deg+1];
   double* yrelohilo = new double[deg+1];
   double* yrehilolo = new double[deg+1];
   double* yrelololo = new double[deg+1];
   double* yimhihihi = new double[deg+1];
   double* yimlohihi = new double[deg+1];
   double* yimhilohi = new double[deg+1];
   double* yimlolohi = new double[deg+1];
   double* yimhihilo = new double[deg+1];
   double* yimlohilo = new double[deg+1];
   double* yimhilolo = new double[deg+1];
   double* yimlololo = new double[deg+1];
   double* zrehihihi_h = new double[deg+1];
   double* zrelohihi_h = new double[deg+1];
   double* zrehilohi_h = new double[deg+1];
   double* zrelolohi_h = new double[deg+1];
   double* zrehihilo_h = new double[deg+1];
   double* zrelohilo_h = new double[deg+1];
   double* zrehilolo_h = new double[deg+1];
   double* zrelololo_h = new double[deg+1];
   double* zimhihihi_h = new double[deg+1];
   double* zimlohihi_h = new double[deg+1];
   double* zimhilohi_h = new double[deg+1];
   double* zimlolohi_h = new double[deg+1];
   double* zimhihilo_h = new double[deg+1];
   double* zimlohilo_h = new double[deg+1];
   double* zimhilolo_h = new double[deg+1];
   double* zimlololo_h = new double[deg+1];
   double* zrehihihi_d = new double[deg+1];
   double* zrelohihi_d = new double[deg+1];
   double* zrehilohi_d = new double[deg+1];
   double* zrelolohi_d = new double[deg+1];
   double* zrehihilo_d = new double[deg+1];
   double* zrelohilo_d = new double[deg+1];
   double* zrehilolo_d = new double[deg+1];
   double* zrelololo_d = new double[deg+1];
   double* zimhihihi_d = new double[deg+1];
   double* zimlohihi_d = new double[deg+1];
   double* zimhilohi_d = new double[deg+1];
   double* zimlolohi_d = new double[deg+1];
   double* zimhihilo_d = new double[deg+1];
   double* zimlohilo_d = new double[deg+1];
   double* zimhilolo_d = new double[deg+1];
   double* zimlololo_d = new double[deg+1];

   for(int k=0; k<=deg; k++)
   {
      xrehihihi[k] = 1.0; xrelohihi[k] = 0.0;
      xrehilohi[k] = 0.0; xrelolohi[k] = 0.0;
      xrehihilo[k] = 0.0; xrelohilo[k] = 0.0;
      xrehilolo[k] = 0.0; xrelololo[k] = 0.0;
      ximhihihi[k] = 0.0; ximlohihi[k] = 0.0;
      ximhilohi[k] = 0.0; ximlolohi[k] = 0.0;
      ximhihilo[k] = 0.0; ximlohilo[k] = 0.0;
      ximhilolo[k] = 0.0; ximlololo[k] = 0.0;
      yrehihihi[k] = 0.0; yrelohihi[k] = 0.0;
      yrehilohi[k] = 0.0; yrelolohi[k] = 0.0;
      yrehihilo[k] = 0.0; yrelohilo[k] = 0.0;
      yrehilolo[k] = 0.0; yrelololo[k] = 0.0;
      yimhihihi[k] = 0.0; yimlohihi[k] = 0.0;
      yimhilohi[k] = 0.0; yimlolohi[k] = 0.0;
      yimhihilo[k] = 0.0; yimlohilo[k] = 0.0;
      yimhilolo[k] = 0.0; yimlololo[k] = 0.0;
   }
   yrehihihi[0] = 1.0; yrehihihi[1] = -1.0;

   CPU_cmplx8_product
      (deg,xrehihihi,xrelohihi,xrehilohi,xrelolohi,
           xrehihilo,xrelohilo,xrehilolo,xrelololo,
           ximhihihi,ximlohihi,ximhilohi,ximlolohi,
           ximhihilo,ximlohilo,ximhilolo,ximlololo,
           yrehihihi,yrelohihi,yrehilohi,yimlolohi,
           yrehihilo,yrelohilo,yrehilolo,yimlololo,
           yimhihihi,yimlohihi,yimhilohi,yimlolohi,
           yimhihilo,yimlohilo,yimhilolo,yimlololo,
           zrehihihi_h,zrelohihi_h,zrehilohi_h,zrelolohi_h,
           zrehihilo_h,zrelohilo_h,zrehilolo_h,zrelololo_h,
           zimhihihi_h,zimlohihi_h,zimhilohi_h,zimlolohi_h,
           zimhihilo_h,zimlohilo_h,zimhilolo_h,zimlololo_h);

   cout << "Series of 1/(1-x) multiplied with 1-x : " << endl;

   for(int k=0; k<=deg; k++)
   {
      cout << "zrehihihi[" << k << "] : " << zrehihihi_h[k];
      cout << "  zrelohihi[" << k << "] : " << zrelohihi_h[k];
      cout << "  zrehilohi[" << k << "] : " << zrehilohi_h[k];
      cout << "  zrelolohi[" << k << "] : " << zrelolohi_h[k] << endl;
      cout << "zrehihilo[" << k << "] : " << zrehihilo_h[k];
      cout << "  zrelohilo[" << k << "] : " << zrelohilo_h[k];
      cout << "  zrehilolo[" << k << "] : " << zrehilolo_h[k];
      cout << "  zrelololo[" << k << "] : " << zrelololo_h[k] << endl;
      cout << "zimhihihi[" << k << "] : " << zimhihihi_h[k];
      cout << "  zimlohihi[" << k << "] : " << zimlohihi_h[k];
      cout << "  zimhilohi[" << k << "] : " << zimhilohi_h[k];
      cout << "  zimlolohi[" << k << "] : " << zimlolohi_h[k] << endl;
      cout << "zimhihilo[" << k << "] : " << zimhihilo_h[k];
      cout << "  zimlohilo[" << k << "] : " << zimlohilo_h[k];
      cout << "  zimhilolo[" << k << "] : " << zimhilolo_h[k];
      cout << "  zimlololo[" << k << "] : " << zimlololo_h[k] << endl;
   }

   GPU_cmplx8_product
      (xrehihihi,xrelohihi,xrehilohi,xrelolohi,
       xrehihilo,xrelohilo,xrehilolo,xrelololo,
       ximhihihi,ximlohihi,ximhilohi,ximlolohi,
       ximhihilo,ximlohilo,ximhilolo,ximlololo,
       yrehihihi,yrelohihi,yrehilohi,yrelolohi,
       yrehihilo,yrelohilo,yrehilolo,yrelololo,
       yimhihihi,yimlohihi,yimhilohi,yimlolohi,
       yimhihilo,yimlohilo,yimhilolo,yimlololo,
       zrehihihi_d,zrelohihi_d,zrehilohi_d,zrelolohi_d,
       zrehihilo_d,zrelohilo_d,zrehilolo_d,zrelololo_d,
       zimhihihi_d,zimlohihi_d,zimhilohi_d,zimlolohi_d,
       zimhihilo_d,zimlohilo_d,zimhilolo_d,zimlololo_d,deg,1,deg+1);

   cout << "GPU computed product :" << endl;

   for(int k=0; k<=deg; k++)
   {
      cout << "zrehihihi[" << k << "] : " << zrehihihi_d[k];
      cout << "  zrelohihi[" << k << "] : " << zrelohihi_d[k];
      cout << "  zrehilohi[" << k << "] : " << zrehilohi_d[k];
      cout << "  zrelolohi[" << k << "] : " << zrelolohi_d[k] << endl;
      cout << "zrehihilo[" << k << "] : " << zrehihilo_d[k];
      cout << "  zrelohilo[" << k << "] : " << zrelohilo_d[k];
      cout << "  zrehilolo[" << k << "] : " << zrehilolo_d[k];
      cout << "  zrelololo[" << k << "] : " << zrelololo_d[k] << endl;
      cout << "zimhihihi[" << k << "] : " << zimhihihi_d[k];
      cout << "  zimlohihi[" << k << "] : " << zimlohihi_d[k];
      cout << "  zimhilohi[" << k << "] : " << zimhilohi_d[k];
      cout << "  zimlolohi[" << k << "] : " << zimlolohi_d[k] << endl;
      cout << "zimhihilo[" << k << "] : " << zimhihilo_d[k];
      cout << "  zimlohilo[" << k << "] : " << zimlohilo_d[k];
      cout << "  zimhilolo[" << k << "] : " << zimhilolo_d[k];
      cout << "  zimlololo[" << k << "] : " << zimlololo_d[k] << endl;
   }
}

void test_real_exponential ( int deg )
{
   double* xhihihi = new double[deg+1];
   double* xlohihi = new double[deg+1];
   double* xhilohi = new double[deg+1];
   double* xlolohi = new double[deg+1];
   double* xhihilo = new double[deg+1];
   double* xlohilo = new double[deg+1];
   double* xhilolo = new double[deg+1];
   double* xlololo = new double[deg+1];
   double* yhihihi = new double[deg+1];
   double* ylohihi = new double[deg+1];
   double* yhilohi = new double[deg+1];
   double* ylolohi = new double[deg+1];
   double* yhihilo = new double[deg+1];
   double* ylohilo = new double[deg+1];
   double* yhilolo = new double[deg+1];
   double* ylololo = new double[deg+1];
   double* zhihihi_h = new double[deg+1];
   double* zlohihi_h = new double[deg+1];
   double* zhilohi_h = new double[deg+1];
   double* zlolohi_h = new double[deg+1];
   double* zhihilo_h = new double[deg+1];
   double* zlohilo_h = new double[deg+1];
   double* zhilolo_h = new double[deg+1];
   double* zlololo_h = new double[deg+1];
   double* zhihihi_d = new double[deg+1];
   double* zlohihi_d = new double[deg+1];
   double* zhilohi_d = new double[deg+1];
   double* zlolohi_d = new double[deg+1];
   double* zhihilo_d = new double[deg+1];
   double* zlohilo_d = new double[deg+1];
   double* zhilolo_d = new double[deg+1];
   double* zlololo_d = new double[deg+1];
   double rhihihi,rlohihi,rhilohi,rlolohi;
   double rhihilo,rlohilo,rhilolo,rlololo;

   random_dbl8_exponentials
      (deg,&rhihihi,&rlohihi,&rhilohi,&rlolohi,
           &rhihilo,&rlohilo,&rhilolo,&rlololo,
           xhihihi,xlohihi,xhilohi,xlolohi,xhihilo,xlohilo,xhilolo,xlololo,
           yhihihi,ylohihi,yhilohi,ylolohi,yhihilo,ylohilo,yhilolo,ylololo);

   CPU_dbl8_product
      (deg,xhihihi,xlohihi,xhilohi,xlolohi,xhihilo,xlohilo,xhilolo,xlololo,
           yhihihi,ylohihi,yhilohi,ylolohi,yhihilo,ylohilo,yhilolo,ylololo,
           zhihihi_h,zlohihi_h,zhilohi_h,zlolohi_h,
           zhihilo_h,zlohilo_h,zhilolo_h,zlololo_h);

   cout << scientific << setprecision(16);

   cout << "Product of series of exp(x) with series of exp(-x)," << endl;
   cout << "  for xhihihi = " << rhihihi;
   cout << "  xlohihi = " << rlohihi << endl;
   cout << "      xhilohi = " << rhilohi;
   cout << "  xlolohi = " << rlolohi << endl;
   cout << "      xhihilo = " << rhihilo;
   cout << "  xlohilo = " << rlohilo << endl;
   cout << "      xhilolo = " << rhilolo;
   cout << "  xlololo = " << rlololo << endl;

   double sumhihihi = 0.0;
   double sumlohihi = 0.0;
   double sumhilohi = 0.0;
   double sumlolohi = 0.0;
   double sumhihilo = 0.0;
   double sumlohilo = 0.0;
   double sumhilolo = 0.0;
   double sumlololo = 0.0;

   for(int k=0; k<=deg; k++)
      odf_inc(&sumhihihi,&sumlohihi,&sumhilohi,&sumlolohi,
              &sumhihilo,&sumlohilo,&sumhilolo,&sumlololo,
              zhihihi_h[k],zlohihi_h[k],zhilohi_h[k],zlolohi_h[k],
              zhihilo_h[k],zlohilo_h[k],zhilolo_h[k],zlololo_h[k]);

   cout << "Summation of all coefficients in the product ..." << endl;
   cout << "    highest part of sum : " << sumhihihi << endl;
   cout << "2nd highest part of sum : " << sumlohihi << endl;
   cout << "3rd highest part of sum : " << sumhilohi << endl;
   cout << "4th highest part of sum : " << sumlolohi << endl;
   cout << " 4th lowest part of sum : " << sumhihilo << endl;
   cout << " 3rd lowest part of sum : " << sumlohilo << endl;
   cout << " 2nd lowest part of sum : " << sumhilolo << endl;
   cout << "     lowest part of sum : " << sumlololo << endl;

   GPU_dbl8_product
      (xhihihi,xlohihi,xhilohi,xlolohi,xhihilo,xlohilo,xhilolo,xlololo,
       yhihihi,ylohihi,yhilohi,ylolohi,yhihilo,ylohilo,yhilolo,ylololo,
       zhihihi_d,zlohihi_d,zhilohi_d,zlolohi_d,
       zhihilo_d,zlohilo_d,zhilolo_d,zlololo_d,deg,1,deg+1);

   cout << "GPU computed product :" << endl;

   for(int k=0; k<=deg; k++)
   {
      cout << "zhihihi[" << k << "] : " << zhihihi_d[k];
      cout << "  zlohihi[" << k << "] : " << zlohihi_d[k];
      cout << "  zhilohi[" << k << "] : " << zhilohi_d[k];
      cout << "  zlolohi[" << k << "] : " << zlolohi_d[k] << endl;
      cout << "zhihilo[" << k << "] : " << zhihilo_d[k];
      cout << "  zlohilo[" << k << "] : " << zlohilo_d[k];
      cout << "  zhilolo[" << k << "] : " << zhilolo_d[k];
      cout << "  zlololo[" << k << "] : " << zlololo_d[k] << endl;
   }

   sumhihihi = 0.0; sumlohihi = 0.0; sumhilohi = 0.0; sumlolohi = 0.0;
   sumhihilo = 0.0; sumlohilo = 0.0; sumhilolo = 0.0; sumlololo = 0.0;

   for(int k=0; k<=deg; k++)
      odf_inc(&sumhihihi,&sumlohihi,&sumhilohi,&sumlolohi,
              &sumhihilo,&sumlohilo,&sumhilolo,&sumlololo,
              zhihihi_d[k],zlohihi_d[k],zhilohi_d[k],zlolohi_d[k],
              zhihilo_d[k],zlohilo_d[k],zhilolo_d[k],zlololo_d[k]);

   cout << "Summation of all coefficients in the GPU computed product ..."
        << endl;
   cout << "    highest part of sum : " << sumhihihi << endl;
   cout << "2nd highest part of sum : " << sumlohihi << endl;
   cout << "3rd highest part of sum : " << sumhilohi << endl;
   cout << "4th highest part of sum : " << sumlolohi << endl;
   cout << " 4th lowest part of sum : " << sumhihilo << endl;
   cout << " 3rd lowest part of sum : " << sumlohilo << endl;
   cout << " 2nd lowest part of sum : " << sumhilolo << endl;
   cout << "     lowest part of sum : " << sumlololo << endl;
}

void test_complex_exponential ( int deg )
{
   double* xrehihihi = new double[deg+1];
   double* xrelohihi = new double[deg+1];
   double* xrehilohi = new double[deg+1];
   double* xrelolohi = new double[deg+1];
   double* xrehihilo = new double[deg+1];
   double* xrelohilo = new double[deg+1];
   double* xrehilolo = new double[deg+1];
   double* xrelololo = new double[deg+1];
   double* ximhihihi = new double[deg+1];
   double* ximlohihi = new double[deg+1];
   double* ximhilohi = new double[deg+1];
   double* ximlolohi = new double[deg+1];
   double* ximhihilo = new double[deg+1];
   double* ximlohilo = new double[deg+1];
   double* ximhilolo = new double[deg+1];
   double* ximlololo = new double[deg+1];
   double* yrehihihi = new double[deg+1];
   double* yrelohihi = new double[deg+1];
   double* yrehilohi = new double[deg+1];
   double* yrelolohi = new double[deg+1];
   double* yrehihilo = new double[deg+1];
   double* yrelohilo = new double[deg+1];
   double* yrehilolo = new double[deg+1];
   double* yrelololo = new double[deg+1];
   double* yimhihihi = new double[deg+1];
   double* yimlohihi = new double[deg+1];
   double* yimhilohi = new double[deg+1];
   double* yimlolohi = new double[deg+1];
   double* yimhihilo = new double[deg+1];
   double* yimlohilo = new double[deg+1];
   double* yimhilolo = new double[deg+1];
   double* yimlololo = new double[deg+1];
   double* zrehihihi_h = new double[deg+1];
   double* zrelohihi_h = new double[deg+1];
   double* zrehilohi_h = new double[deg+1];
   double* zrelolohi_h = new double[deg+1];
   double* zrehihilo_h = new double[deg+1];
   double* zrelohilo_h = new double[deg+1];
   double* zrehilolo_h = new double[deg+1];
   double* zrelololo_h = new double[deg+1];
   double* zimhihihi_h = new double[deg+1];
   double* zimlohihi_h = new double[deg+1];
   double* zimhilohi_h = new double[deg+1];
   double* zimlolohi_h = new double[deg+1];
   double* zimhihilo_h = new double[deg+1];
   double* zimlohilo_h = new double[deg+1];
   double* zimhilolo_h = new double[deg+1];
   double* zimlololo_h = new double[deg+1];
   double* zrehihihi_d = new double[deg+1];
   double* zrelohihi_d = new double[deg+1];
   double* zrehilohi_d = new double[deg+1];
   double* zrelolohi_d = new double[deg+1];
   double* zrehihilo_d = new double[deg+1];
   double* zrelohilo_d = new double[deg+1];
   double* zrehilolo_d = new double[deg+1];
   double* zrelololo_d = new double[deg+1];
   double* zimhihihi_d = new double[deg+1];
   double* zimlohihi_d = new double[deg+1];
   double* zimhilohi_d = new double[deg+1];
   double* zimlolohi_d = new double[deg+1];
   double* zimhihilo_d = new double[deg+1];
   double* zimlohilo_d = new double[deg+1];
   double* zimhilolo_d = new double[deg+1];
   double* zimlololo_d = new double[deg+1];
   double rndrehihihi,rndrelohihi,rndrehilohi,rndrelolohi;
   double rndrehihilo,rndrelohilo,rndrehilolo,rndrelololo;
   double rndimhihihi,rndimlohihi,rndimhilohi,rndimlolohi;
   double rndimhihilo,rndimlohilo,rndimhilolo,rndimlololo;

   random_cmplx8_exponentials
      (deg,&rndrehihihi,&rndrelohihi,&rndrehilohi,&rndrelolohi,
           &rndrehihilo,&rndrelohilo,&rndrehilolo,&rndrelololo,
           &rndimhihihi,&rndimlohihi,&rndimhilohi,&rndimlolohi,
           &rndimhihilo,&rndimlohilo,&rndimhilolo,&rndimlololo,
           xrehihihi,xrelohihi,xrehilohi,xrelolohi,
           xrehihilo,xrelohilo,xrehilolo,xrelololo,
           ximhihihi,ximlohihi,ximhilohi,ximlolohi,
           ximhihilo,ximlohilo,ximhilolo,ximlololo,
           yrehihihi,yrelohihi,yrehilohi,yrelolohi,
           yrehihilo,yrelohilo,yrehilolo,yrelololo,
           yimhihihi,yimlohihi,yimhilohi,yimlolohi,
           yimhihilo,yimlohilo,yimhilolo,yimlololo);

   CPU_cmplx8_product
      (deg,xrehihihi,xrelohihi,xrehilohi,xrelolohi,
           xrehihilo,xrelohilo,xrehilolo,xrelololo,
           ximhihihi,ximlohihi,ximhilohi,ximlolohi,
           ximhihilo,ximlohilo,ximhilolo,ximlololo,
           yrehihihi,yrelohihi,yrehilohi,yrelolohi,
           yrehihilo,yrelohilo,yrehilolo,yrelololo,
           yimhihihi,yimlohihi,yimhilohi,yimlolohi,
           yimhihilo,yimlohilo,yimhilolo,yimlololo,
           zrehihihi_h,zrelohihi_h,zrehilohi_h,zrelolohi_h,
           zrehihilo_h,zrelohilo_h,zrehilolo_h,zrelololo_h,
           zimhihihi_h,zimlohihi_h,zimhilohi_h,zimlolohi_h,
           zimhihilo_h,zimlohilo_h,zimhilolo_h,zimlololo_h);

   cout << scientific << setprecision(16);

   cout << "Product of series of exp(x) with series of exp(-x)," << endl;
   cout << "  for xrehihihi = " << rndrehihihi;
   cout << "      xrelohihi = " << rndrelohihi << endl;
   cout << "      xrehilohi = " << rndrehilohi;
   cout << "      xrelolohi = " << rndrelolohi << endl;
   cout << "      xrehihilo = " << rndrehihilo;
   cout << "      xrelohilo = " << rndrelohilo << endl;
   cout << "      xrehilolo = " << rndrehilolo;
   cout << "  and xrelololo = " << rndrelololo << endl;
   cout << "  for ximhihihi = " << rndimhihihi;
   cout << "      ximlohihi = " << rndimlohihi << endl;
   cout << "      ximhilohi = " << rndimhilohi;
   cout << "      ximlolohi = " << rndimlolohi << endl;
   cout << "      ximhihilo = " << rndimhihilo;
   cout << "      ximlohilo = " << rndimlohilo << endl;
   cout << "      ximhilolo = " << rndimhilolo;
   cout << "  and ximlololo = " << rndimlololo << endl;

   double sumrehihihi = 0.0;
   double sumrelohihi = 0.0;
   double sumrehilohi = 0.0;
   double sumrelolohi = 0.0;
   double sumrehihilo = 0.0;
   double sumrelohilo = 0.0;
   double sumrehilolo = 0.0;
   double sumrelololo = 0.0;
   double sumimhihihi = 0.0;
   double sumimlohihi = 0.0;
   double sumimhilohi = 0.0;
   double sumimlolohi = 0.0;
   double sumimhihilo = 0.0;
   double sumimlohilo = 0.0;
   double sumimhilolo = 0.0;
   double sumimlololo = 0.0;

   for(int k=0; k<=deg; k++) 
   {
      odf_inc(&sumrehihihi,&sumrelohihi,&sumrehilohi,&sumrelolohi,
              &sumrehihilo,&sumrelohilo,&sumrehilolo,&sumrelololo,
              zrehihihi_h[k],zrelohihi_h[k],zrehilohi_h[k],zrelolohi_h[k],
              zrehihilo_h[k],zrelohilo_h[k],zrehilolo_h[k],zrelololo_h[k]);
      odf_inc(&sumimhihihi,&sumimlohihi,&sumimhilohi,&sumimlolohi,
              &sumimhihilo,&sumimlohilo,&sumimhilolo,&sumimlololo,
              zimhihihi_h[k],zimlohihi_h[k],zimhilohi_h[k],zimlolohi_h[k],
              zimhihilo_h[k],zimlohilo_h[k],zimhilolo_h[k],zimlololo_h[k]);
   }
   cout << "Summation of all coefficients of the product ..." << endl;
   cout << "  sumrehihihi : " << sumrehihihi;
   cout << "  sumrelohihi : " << sumrelohihi << endl;
   cout << "  sumrehilohi : " << sumrehilohi;
   cout << "  sumrelolohi : " << sumrelolohi << endl;
   cout << "  sumrehihilo : " << sumrehihilo;
   cout << "  sumrelohilo : " << sumrelohilo << endl;
   cout << "  sumrehilolo : " << sumrehilolo;
   cout << "  sumrelololo : " << sumrelololo << endl;
   cout << "  sumimhihihi : " << sumimhihihi;
   cout << "  sumimlohihi : " << sumimlohihi << endl;
   cout << "  sumimhilohi : " << sumimhilohi;
   cout << "  sumimlolohi : " << sumimlolohi << endl;
   cout << "  sumimhihilo : " << sumimhihilo;
   cout << "  sumimlohilo : " << sumimlohilo << endl;
   cout << "  sumimhilolo : " << sumimhilolo;
   cout << "  sumimlololo : " << sumimlololo << endl;

   GPU_cmplx8_product
      (xrehihihi,xrelohihi,xrehilohi,xrelolohi,
       xrehihilo,xrelohilo,xrehilolo,xrelololo,
       ximhihihi,ximlohihi,ximhilohi,ximlolohi,
       ximhihilo,ximlohilo,ximhilolo,ximlololo,
       yrehihihi,yrelohihi,yrehilohi,yrelolohi,
       yrehihilo,yrelohilo,yrehilolo,yrelololo,
       yimhihihi,yimlohihi,yimhilohi,yimlolohi,
       yimhihilo,yimlohilo,yimhilolo,yimlololo,
       zrehihihi_d,zrelohihi_d,zrehilohi_d,zrelolohi_d,
       zrehihilo_d,zrelohilo_d,zrehilolo_d,zrelololo_d,
       zimhihihi_d,zimlohihi_d,zimhilohi_d,zimlolohi_d,
       zimhihilo_d,zimlohilo_d,zimhilolo_d,zimlololo_d,deg,1,deg+1);

   sumrehihihi = 0.0; sumrelohihi = 0.0;
   sumrehilohi = 0.0; sumrelolohi = 0.0;
   sumrehihilo = 0.0; sumrelohilo = 0.0;
   sumrehilolo = 0.0; sumrelololo = 0.0;
   sumimhihihi = 0.0; sumimlohihi = 0.0;
   sumimhilohi = 0.0; sumimlolohi = 0.0;
   sumimhihilo = 0.0; sumimlohilo = 0.0;
   sumimhilolo = 0.0; sumimlololo = 0.0;

   for(int k=0; k<=deg; k++) 
   {
      odf_inc(&sumrehihihi,&sumrelohihi,&sumrehilohi,&sumrelolohi,
              &sumrehihilo,&sumrelohilo,&sumrehilolo,&sumrelololo,
              zrehihihi_d[k],zrelohihi_d[k],zrehilohi_d[k],zrelolohi_d[k],
              zrehihilo_d[k],zrelohilo_d[k],zrehilolo_d[k],zrelololo_d[k]);
      odf_inc(&sumimhihihi,&sumimlohihi,&sumimhilohi,&sumimlolohi,
              &sumimhihilo,&sumimlohilo,&sumimhilolo,&sumimlololo,
              zimhihihi_d[k],zimlohihi_d[k],zimhilohi_d[k],zimlolohi_d[k],
              zimhihilo_d[k],zimlohilo_d[k],zimhilolo_d[k],zimlololo_d[k]);
   }
   cout << "Summation of all coefficients of the GPU computed product ..."
        << endl;
   cout << "  sumrehihihi : " << sumrehihihi;
   cout << "  sumrelohihi : " << sumrelohihi << endl;
   cout << "  sumrehilohi : " << sumrehilohi;
   cout << "  sumrelolohi : " << sumrelolohi << endl;
   cout << "  sumrehihilo : " << sumrehihilo;
   cout << "  sumrelohilo : " << sumrelohilo << endl;
   cout << "  sumrehilolo : " << sumrehilolo;
   cout << "  sumrelololo : " << sumrelololo << endl;
   cout << "  sumimhihihi : " << sumimhihihi;
   cout << "  sumimlohihi : " << sumimlohihi << endl;
   cout << "  sumimhilohi : " << sumimhilohi;
   cout << "  sumimlolohi : " << sumimlolohi << endl;
   cout << "  sumimhihilo : " << sumimhihilo;
   cout << "  sumimlohilo : " << sumimlohilo << endl;
   cout << "  sumimhilolo : " << sumimhilolo;
   cout << "  sumimlololo : " << sumimlololo << endl;
}
