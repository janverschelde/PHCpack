/* Test on the operations declared in the octo_double.h file. */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include "octo_double.h"

int basic_test ( void );
/* 
 * DESCRIPTION :
 *   Tests basic arithmetic on randomly generated octo doubles. */

int main ( void )
{
   int fail = basic_test();

   return 0;
}

int basic_test ( void )
{
   double a[8],b[8],c[8],d[8];

   printf("Testing some basic operations ...\n");
   srand(time(NULL));

   od_random(a);
   printf("the first random number a :\n"); od_write_doubles(a);
   od_random(b);
   printf("the second random number b :\n"); od_write_doubles(b);

   od_add(a,b,c); printf("a+b :\n"); od_write_doubles(c);
   od_sub(c,a,d); printf("a+b-a :\n"); od_write_doubles(d);

   od_mul(a,b,c); printf("a*b :\n"); od_write_doubles(c);
   od_div(c,a,d); printf("a*b/a :\n"); od_write_doubles(d);

   return 0;
}
