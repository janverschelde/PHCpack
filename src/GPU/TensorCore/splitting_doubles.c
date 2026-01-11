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

uint64 last_bits ( int k, uint64 nbr )
{
   uint64 mask = 1;
   uint64 pwr = 2;

   for(int i=1; i<k; i++)
   {
      mask = mask + pwr;
      pwr = 2*pwr;
   }
   return (nbr & mask);
}

double first_half ( double x, int vrblvl )
{
   int exponent;
   double fraction = frexp(x, &exponent );
   double shifted = ldexp(fraction, 52);
   uint64 int64fac = (uint64) shifted;
   uint64 second_part = last_bits(26,int64fac);
   uint64 first_part = int64fac - second_part;
   double result = ldexp(first_part,exponent-52);

   if(vrblvl > 0)
   {
      printf("exponent : %d\n", exponent);
      printf("fraction : %lld\n", int64fac);
      printf(" h : %llx\n", int64fac);
      printf(" b : "); write_52bits(52,int64fac); printf("\n");
      printf("b0 : "); write_52bits(52,first_part); printf("\n");
      printf("b1 : "); write_52bits(52,second_part); printf("\n");
   }

   return result;
}

void half_split ( double x, double *x0, double *x1, int vrblvl )
{
   *x0 = first_half(x,vrblvl);
   *x1 = x - *x0;
}

int test ( void )
{
   double x,x0,x1,s,e;

   srand(time(NULL));

   x = ((double) rand())/RAND_MAX;
   printf("x : %.15e\n", x);

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

int main ( void )
{
   return test();
}
