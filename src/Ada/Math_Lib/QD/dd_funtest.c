/* Test on the sin and cos functions in double double precision. */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "double_double.h"

int test_sincostaylor ( void );
/*
 * DESCRIPTION :
 *   Tests sin and cos on pi/32 using Taylor expansions. */

int test_sincos ( void );
/*
 * DESCRIPTION :
 *   Tests sin and cos on random doubles. */

int main ( void )
{
   int fail = test_sincostaylor();

   fail = test_sincos();

   return 0;
}

int test_sincostaylor ( void )
{
   double dd_pi16[2];
   double dd_pi32[2];
   double sinpi32a[2],sinpi32b[2];
   double cospi32a[2],cospi32b[2];
   double x[2];
   double y[2];
   const double pi16_hi = 1.963495408493620697e-01; // high part of pi/16
   const double pi16_lo = 7.654042494670957545e-18; // low part of pi/16

   dd_pi16[0] = pi16_hi;
   dd_pi16[1] = pi16_lo;

   dd_mul_pwr2(dd_pi16,0.5,dd_pi32);

   dd_sin_taylor(dd_pi32,sinpi32a);
   dd_cos_taylor(dd_pi32,cospi32a);
   dd_sincos_taylor(dd_pi32,sinpi32b,cospi32b);

   printf("test on sincostaylor ...\n");

   printf("sin(pi/32) : "); dd_write(sinpi32a,32); printf("\n");
   printf("sin(pi/32) : "); dd_write(sinpi32b,32); printf("\n");
   printf("cos(pi/32) : "); dd_write(cospi32a,32); printf("\n");
   printf("cos(pi/32) : "); dd_write(cospi32b,32); printf("\n");

   printf("cospi32b hi : %.16e\n", cospi32b[0]);
   printf("cospi32b lo : %.16e\n", cospi32b[1]);

   dd_sqr(sinpi32a,x);
   dd_sqr(cospi32a,y);
   dd_inc(x,y);
   printf("sin(pi/32)^2 + cos(pi/32)^2 : "); dd_write(x,32);
   printf("\n");

   return 0;
}

int test_sincos ( void )
{
   const double pi16_hi = 1.963495408493620697e-01; // high part of pi/16
   const double pi16_lo = 7.654042494670957545e-18; // low part of pi/16
   double rnd;

   srand(time(NULL));
   rnd = (double) rand();
   rnd = rnd/RAND_MAX;
   rnd = 6.0*rnd - 3.0;

   double ddrnd[2];
   ddrnd[0] = rnd; ddrnd[1] = 0.0;

   printf("test on sin and cos of a random number r ...\n");

   printf("r : "); dd_write(ddrnd,32); printf("\n");

   double sinrnd[2],cosrnd[2],x[2],y[2];

   printf("computing sin(r) ...\n");
   dd_sin(ddrnd,sinrnd);
   printf("sin(r) : "); dd_write(sinrnd,32); printf("\n");
   printf("computing cos(r) ...\n");
   dd_cos(ddrnd,cosrnd);
   printf("cos(r) : "); dd_write(cosrnd,32); printf("\n");

   dd_sqr(sinrnd,x);
   dd_sqr(cosrnd,y);
   dd_inc(x,y);
   printf("sin(r)^2 + cos(r)^2 : "); dd_write(x,32); printf("\n");

   ddrnd[0] = pi16_hi; ddrnd[1] = pi16_lo;

   printf("test on sin and cos of r = pi/32 ...\n");

   printf("r : "); dd_write(ddrnd,32); printf("\n");

   printf("computing sin(r) ...\n"); 
   dd_sin(ddrnd,sinrnd);
   printf("sin(r) : "); dd_write(sinrnd,32); printf("\n");
   printf("computing cos(r) ...\n");
   dd_cos(ddrnd,cosrnd);
   printf("cos(r) : "); dd_write(cosrnd,32); printf("\n");

   dd_sqr(sinrnd,x);
   dd_sqr(cosrnd,y);
   dd_inc(x,y);
   printf("sin(r)^2 + cos(r)^2 : "); dd_write(x,32); printf("\n");

   return 0;
}
