/* Test on the operations declared in the hexa_double.h file. */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include "hexa_double.h"

int basic_test ( void );
/* 
 * DESCRIPTION :
 *   Tests basic arithmetic on randomly generated hexa doubles. */

int main ( void )
{
   int fail = basic_test();

   return 0;
}

int basic_test ( void )
{
   double a[16],b[16],c[16],d[16];

   printf("Testing some basic operations ...\n");
   srand(time(NULL));

   hd_random(a);
   printf("the first random number a :\n"); hd_write_doubles(a);
   hd_random(b);
   printf("the second random number b :\n"); hd_write_doubles(b);

   hd_add(a,b,c); printf("a+b :\n"); hd_write_doubles(c);
   hd_sub(c,a,d); printf("a+b-a :\n"); hd_write_doubles(d);

   hd_mul(a,b,c); printf("a*b :\n"); hd_write_doubles(c);
   hd_div(c,a,d); printf("a*b/a :\n"); hd_write_doubles(d);

   return 0;
}
