/* Test on the operations declared in the triple_double.h file. */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include "triple_double.h"

void write ( const double* x );
/*
 * DESCRIPTION :
 *   Basic output of the three doubles in a given
 *   triple double number x. */

int basic_test ( void );
/* 
 * DESCRIPTION :
 *   A test on the basic functions of triple_double consists in
 *   generate two random doubles a, b and then calling the operations. */

int main ( void )
{
   int fail = basic_test();

   return 0;
}

void write ( const double* x )
{
   printf("  hi = %21.14e\n",x[0]);
   printf("  mi = %21.14e\n",x[1]);
   printf("  lo = %21.14e\n",x[2]);
}

int basic_test ( void )
{
   double a[3],b[3],c[3],d[3];

   printf("\ntesting some basic operations ...\n");
   srand(time(NULL));

   random_triple_double(a);
   printf("the first random number a :\n"); write(a);
   random_triple_double(b);
   printf("the second random number b :\n"); write(b);

   td_add(a,b,c); printf("a+b :\n"); write(c);
   td_sub(c,a,d); printf("a+b-a :\n"); write(d);

   td_mul(a,b,c); printf("a*b :\n"); write(c);
   td_div(c,a,d); printf("a*b/a :\n"); write(d);

   return 0;
}
