/* Test on the operations declared in the triple_double.h file. */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include "triple_double.h"

int basic_test ( void );
/* 
 * DESCRIPTION :
 *   Tests the basic arithmetic on randomly generated triple doubles. */

int main ( void )
{
   int fail = basic_test();

   return 0;
}

int basic_test ( void )
{
   double a[3],b[3],c[3],d[3];

   printf("Testing some basic operations ...\n");
   srand(time(NULL));

   td_random(a);
   printf("the first random number a :\n"); td_write_doubles(a);
   td_random(b);
   printf("the second random number b :\n"); td_write_doubles(b);

   td_add(a,b,c); printf("a+b :\n"); td_write_doubles(c);
   td_sub(c,a,d); printf("a+b-a :\n"); td_write_doubles(d);

   td_mul(a,b,c); printf("a*b :\n"); td_write_doubles(c);
   td_div(c,a,d); printf("a*b/a :\n"); td_write_doubles(d);

   return 0;
}
