/* Tests the collection of functions to split doubles. */

#include <stdio.h>
#include <time.h>
#include <math.h>
#include "splitting_doubles.h"

int test_half_split ( void );
/*
 * Generates a random number, splits in two equal sized halves, and then
 * checks if adding the parts gives the original number. */

int test_quarter_split ( void );
/*
 * Generates a random number, splits in four equal sized halves, and then
 * checks if adding the parts gives the original number. */

int main ( void )
{
   int fail;

   fail = test_half_split();
   if(fail == 1)
      printf("\nTest on half split failed?!!!\n\n");
   else
      printf("\nTest on half split succeeded.\n\n");

   fail = test_quarter_split();
   if(fail == 1)
      printf("\nTest on quarter split failed?!!!\n\n");
   else
      printf("\nTest on quarter split succeeded.\n\n");

   return fail;
}

int test_half_split ( void )
{
   double x,x0,x1,s,e;

   srand(time(NULL));

   x = ((double) rand())/RAND_MAX;
   printf(" x : %.15e\n", x);

   half_split(x, &x0, &x1, 1);

   printf("x0 : %.15e\n", x0);
   printf("x1 : %.15e\n", x1);
   s = x0 + x1;
   printf("x0 + x1 : %.15e\n", s);
   printf("      x : %.15e\n", x);
   e = abs(x - s);
   printf("  error : %.3e\n", e);

   return not(e == 0);
}

int test_quarter_split ( void )
{
   double x,x0,x1,x2,x3,s,e;

   srand(time(NULL));

   x = ((double) rand())/RAND_MAX;
   printf(" x : %.15e\n", x);

   quarter_split(x, &x0, &x1, &x2, &x3, 1);

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

   return not(e == 0.0);
}

