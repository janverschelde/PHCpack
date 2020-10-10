/* Test on the operations declared in the deca_double.h file. */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include "deca_double.h"

int basic_test ( void );
/* 
 * DESCRIPTION :
 *   Tests basic arithmetic on randomly generated deca doubles. */

int main ( void )
{
   int fail = basic_test();

   return 0;
}

int basic_test ( void )
{
   double a[10],b[10],c[10],d[10];

   printf("Testing some basic operations ...\n");
   srand(time(NULL));

   da_random(a);
   printf("the first random number a :\n"); da_write_doubles(a);
   da_random(b);
   printf("the second random number b :\n"); da_write_doubles(b);

   da_add(a,b,c); printf("a+b :\n"); da_write_doubles(c);
   da_sub(c,a,d); printf("a+b-a :\n"); da_write_doubles(d);

   da_mul(a,b,c); printf("a*b :\n"); da_write_doubles(c);
   da_div(c,a,d); printf("a*b/a :\n"); da_write_doubles(d);

   return 0;
}
