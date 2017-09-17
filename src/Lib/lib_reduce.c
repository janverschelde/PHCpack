/* Tests the reduction of systems through the C interface. */

#include <stdio.h>
#include <stdlib.h>
#include "syscon.h"
#include "phcpack.h"
#include "reducers.h"

int greetings( void );
/*
 * DESCRIPTION :
 *   Displays the version string of PHCpack. */

int standard_reducer ( void );
/*
 * DESCRIPTION :
 *   Prompts the user for a system and performs the scaling
 *   with standard double precision arithmetic. */

int main ( int argc, char *argv[] )
{
   int fail,choice;

   adainit();

   fail = greetings();

   fail = standard_reducer();

   adafinal();

   return fail;
}

int greetings( void )
{
   int fail,n;
   char s[40];

   fail = version_string(&n,s);

   if(fail == 0) printf("Testing the scaling in %s ...\n",s);

   return fail;
}

int standard_reducer ( void )
{
   int fail,dim,i;

   fail = syscon_read_standard_system();
   if(fail == 0)
   {
      printf("\nThe system in the container : \n");
      fail = syscon_write_standard_system();
      fail = standard_reduce_system(1);
      printf("\nThe system after reduction : \n");
      fail = syscon_write_standard_system();
   }
   return fail;
}
