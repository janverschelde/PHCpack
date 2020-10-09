/* Test on the operations declared in the penta_double.h file. */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include "penta_double.h"

int basic_test ( void );
/* 
 * DESCRIPTION :
 *   Tests basic arithmetic on randomly generated penta doubles. */

int main ( void )
{
   int fail = basic_test();

   return 0;
}

int basic_test ( void )
{
   double a[5],b[5],c[5],d[5];

   printf("Testing some basic operations ...\n");
   srand(time(NULL));

   pd_random(a);
   printf("the first random number a :\n"); pd_write_doubles(a);
   pd_random(b);
   printf("the second random number b :\n"); pd_write_doubles(b);

   pd_add(a,b,c); printf("a+b :\n"); pd_write_doubles(c);
   pd_sub(c,a,d); printf("a+b-a :\n"); pd_write_doubles(d);

   pd_mul(a,b,c); printf("a*b :\n"); pd_write_doubles(c);
   pd_div(c,a,d); printf("a*b/a :\n"); pd_write_doubles(d);

   return 0;
}
