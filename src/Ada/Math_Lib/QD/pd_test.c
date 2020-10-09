/* Test on the operations declared in the penta_double.h file. */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include "penta_double.h"

void write ( const double* x );
/*
 * DESCRIPTION :
 *   Basic output of the three doubles in a penta double number x. */

int basic_test ( void );
/* 
 * DESCRIPTION :
 *   A test on the basic functions of penta_double consists in
 *   generate two random doubles a, b and then calling the operations. */

int main ( void )
{
   int fail = basic_test();

   return 0;
}

void write ( const double* x )
{
   printf("  thumb = %21.14e",x[0]);
   printf("  index = %21.14e\n",x[1]);
   printf("  middle = %21.14e\n",x[2]);
   printf("  ring = %21.14e",x[3]);
   printf("   pink = %21.14e\n",x[4]);
}

int basic_test ( void )
{
   double a[5],b[5],c[5],d[5];

   printf("Testing some basic operations ...\n");
   srand(time(NULL));

   random_penta_double(a);
   printf("the first random number a :\n"); write(a);
   random_penta_double(b);
   printf("the second random number b :\n"); write(b);

   pd_add(a,b,c); printf("a+b :\n"); write(c);
   pd_sub(c,a,d); printf("a+b-a :\n"); write(d);

   pd_mul(a,b,c); printf("a*b :\n"); write(c);
   pd_div(c,a,d); printf("a*b/a :\n"); write(d);

   return 0;
}
