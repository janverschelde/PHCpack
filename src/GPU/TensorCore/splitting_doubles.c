/* Collection of functions to split doubles. */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "splitting_doubles.h"

void write_52bits ( int k, uint64 nbr )
{
   if(k > 0)
   {
      write_52bits(k-1,nbr/2);
      printf("%d", nbr%2);
      if(k % 4 == 0) printf(" ");
   }
}

void write_52double ( double nbr )
{
   int exponent;
   double fraction = frexp(nbr, &exponent );
   double shifted = ldexp(fraction, 52);
   uint64 int64fac = (uint64) shifted;

   write_52bits(52, int64fac);
   printf(" %d\n", exponent);
}

uint64 last_bits ( int k, uint64 nbr )
{
   uint64 mask = 1;
   uint64 pwr = 2;

   for(int i=1; i<k; i++)
   {
      mask = mask + pwr;
      pwr = 2*pwr;
   }
   return (nbr & mask); // bit masking
}

uint64 quarter_bits
 ( uint64 nbr, uint64 *b0, uint64 *b1, uint64 *b2, uint64 *b3, int vrblvl )
{
   uint64 mask = 1;
   uint64 pwr = 2;
   uint64 rest = nbr;

   for(int i=1; i<13; i++)
   {
      mask = mask + pwr;
      pwr = 2*pwr;
   }
   if(vrblvl > 0)
      printf(" m : "); write_52bits(52, mask); printf("\n");

   *b3 = (nbr & mask); // last 13 bits
   rest = rest - *b3;
   
   if(vrblvl > 0)
   {
      printf(" h : %llx\n", nbr);
      printf(" b : "); write_52bits(52, nbr); printf("\n");
      printf("b3 : "); write_52bits(52, *b3); printf("\n");
   }
   for(int i=0; i<13; i++)
   {
      mask = mask + pwr;
      pwr = 2*pwr;
   }
   if(vrblvl > 0)
      printf(" m : "); write_52bits(52, mask); printf("\n");

   *b2 = (rest & mask); // next 13 bits
   rest = rest - *b2;

   if(vrblvl > 0)
   {
      printf("b2 : "); write_52bits(52,*b2); printf("\n");
   }
   for(int i=0; i<13; i++)
   {
      mask = mask + pwr;
      pwr = 2*pwr;
   }
   if(vrblvl > 0)
      printf(" m : "); write_52bits(52, mask); printf("\n");

   *b1 = (rest & mask); // next 13 bits

   if(vrblvl > 0)
   {
      printf("b1 : "); write_52bits(52,*b1); printf("\n");
   }
   *b0 = nbr - *b1 - *b2 - *b3;

   if(vrblvl > 0)
   {
      printf("b0 : "); write_52bits(52, *b0); printf("\n");
      rest = nbr - *b0 - *b1 - *b2 - *b3;
      printf(" r : "); write_52bits(52, rest);
      if(rest == 0)
         printf(" OKAY\n");
      else
         printf(" BUG!\n");
   }
}

double first_half ( double x, int vrblvl )
{
   int exponent;
   double fraction = frexp(x, &exponent );
   double shifted = ldexp(fraction, 52);
   uint64 int64fac = (uint64) shifted;
   uint64 second_part = last_bits(26, int64fac);
   uint64 first_part = int64fac - second_part;
   double result = ldexp(first_part, exponent-52);

   if(vrblvl > 0)
   {
      printf("exponent : %d\n", exponent);
      printf("fraction : %lld\n", int64fac);
      printf(" h : %llx\n", int64fac);
      printf(" b : "); write_52bits(52, int64fac); printf("\n");
      printf("b0 : "); write_52bits(52, first_part); printf("\n");
      printf("b1 : "); write_52bits(52, second_part); printf("\n");
   }
   return result;
}

void half_split ( double x, double *x0, double *x1, int vrblvl )
{
   *x0 = first_half(x, vrblvl);
   *x1 = x - *x0;
}

int leading_zeros ( uint64 nbr, int idxpwr, int vrblvl )
{
   int result = 0;
   int idx = idxpwr;
   const uint64 two = 2;
   uint64 threshold = two << idx; // 2**idx

   if(vrblvl > 0)
       printf("threshold : %lld, nbr : %lld\n", threshold, nbr);

   while(nbr < threshold)
   {
      result = result + 1;
      threshold = threshold/2;
      if(vrblvl > 0)
         printf("threshold : %lld, cnt : %d\n", threshold, result);
   }
   
   return result;
}

void quarter_split
 ( double x, double *x0, double *x1, double *x2, double *x3, int vrblvl )
{
   int exponent;
   double fraction = frexp(x, &exponent);
   double shifted = ldexp(fraction, 52);
   uint64 int64fac = (uint64) shifted;
   uint64 f0,f1,f2,f3;
   double xf0,xf1,xf2,xf3;
   int cnt;

   if(vrblvl > 0) printf("exponent : %d\n", exponent);

   quarter_bits(int64fac,&f0,&f1,&f2,&f3,vrblvl);

   xf0 = (double) f0;
   *x0 = ldexp(xf0, exponent - 52);

   cnt = (f1 == 0) ? 0 : leading_zeros(f1, 37, vrblvl);
   if(vrblvl > 0) printf("#leading zeros in f1 : %d\n", cnt);
   xf1 = (double) f1;
   // *x1 = ldexp(xf1, exponent - 52 - cnt); // does not work ...
   if(cnt == 0)
      *x1 = ldexp(xf1, exponent - 52);
   else
   {
      for(int i=0; i<cnt; i++) f1 = 2*f1;
      xf1 = (double) f1;
      *x1 = ldexp(xf1, exponent - 52);
      for(int i=0; i<cnt; i++) *x1 = *x1/2.0;
   }
   cnt = (f2 == 0) ? 0 : leading_zeros(f2, 24, vrblvl);
   if(vrblvl > 0) printf("#leading zeros in f2 : %d\n", cnt);
   xf2 = (double) f2;
   if(cnt == 0)
      *x2 = ldexp(xf2, exponent - 52);
   else
   {
      for(int i=0; i<cnt; i++) f2 = 2*f2;
      xf2 = (double) f2;
      *x2 = ldexp(xf2, exponent - 52);
      for(int i=0; i<cnt; i++) *x2 = *x2/2.0;
   }
   if(vrblvl > 0) // just for testing purposes
   {
      cnt = (f3 == 0) ? 0 : leading_zeros(f3, 11, vrblvl);
      printf("#leading zeros in f3 : %d\n", cnt);
   }
   *x3 = x - *x0 - *x1 - *x2;
}

int test_half_split ( void )
{
   double x,x0,x1,s,e;

   srand(time(NULL));

   x = ((double) rand())/RAND_MAX;
   printf(" x : %.15e\n", x);

   half_split(x,&x0,&x1,1);

   printf("x0 : %.15e\n", x0);
   printf("x1 : %.15e\n", x1);
   s = x0 + x1;
   printf("x0 + x1 : %.15e\n", s);
   printf("      x : %.15e\n", x);
   e = abs(x - s);
   printf("  error : %.3e\n", e);

   return 0;
}

int test_quarter_split ( void )
{
   double x,x0,x1,x2,x3,s,e;

   srand(time(NULL));

   x = ((double) rand())/RAND_MAX;
   printf(" x : %.15e\n", x);

   quarter_split(x,&x0,&x1,&x2,&x3,1);

   printf("x0 : %.15e\n", x0);
   printf(" b : "); write_52double(x0);
   printf("x1 : %.15e\n", x1);
   printf(" b : "); write_52double(x1);
   printf("x2 : %.15e\n", x2);
   printf(" b : "); write_52double(x2);
   printf("x3 : %.15e\n", x3);
   printf(" b : "); write_52double(x3);

   printf("                x : %.15e\n", x);
   e = fabs(x - x0);
   printf("               x0 : %.15e, error : %.3e\n", x0, e);
   s = x0 + x1; e = fabs(x - s);
   printf("          x0 + x1 : %.15e, error : %.3e\n", s, e);
   s = s + x2; e = fabs(x - s);
   printf("     x0 + x1 + x2 : %.15e, error : %.3e\n", s, e);
   s = s + x3; e = fabs(x - s);
   printf("x0 + x1 + x2 + x3 : %.15e, error : %.3e\n", s, e);

   return 0;
}

int main ( void )
{
   // return test_half_split();
   return test_quarter_split();
}
