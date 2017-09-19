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
 *   Prompts the user for a system and performs linear reduction
 *   with standard double precision arithmetic. */

int dobldobl_reducer ( void );
/*
 * DESCRIPTION :
 *   Prompts the user for a system and performs linear reduction
 *   with double double precision arithmetic. */

int quaddobl_reducer ( void );
/*
 * DESCRIPTION :
 *   Prompts the user for a system and performs linear reduction
 *   with quad double precision arithmetic. */

int main ( int argc, char *argv[] )
{
   int fail,choice;

   adainit();

   fail = greetings();

   printf("\nMENU for the precision of linear reduction :\n"); 
   printf("  0. use standard double precision arithmetic; or\n");
   printf("  1. use double double precision arithmetic; or\n");
   printf("  2. use quad double precision arithmetic.\n");
   printf("Type 0, 1, or 2 to make your choice : ");
   scanf("%d",&choice);

   if(choice == 0)
      fail = standard_reducer();
   else if(choice == 1)
      fail = dobldobl_reducer();
   else if(choice == 2)
      fail = quaddobl_reducer();

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
      fail = standard_row_reduce_system(1);
      printf("\nThe system after reduction : \n");
      fail = syscon_write_standard_system();
   }
   return fail;
}

int dobldobl_reducer ( void )
{
   int fail,dim,i;

   fail = syscon_read_dobldobl_system();
   if(fail == 0)
   {
      printf("\nThe system in the container : \n");
      fail = syscon_write_dobldobl_system();
      fail = dobldobl_row_reduce_system(1);
      printf("\nThe system after reduction : \n");
      fail = syscon_write_dobldobl_system();
   }
   return fail;
}

int quaddobl_reducer ( void )
{
   int fail,dim,i;

   fail = syscon_read_quaddobl_system();
   if(fail == 0)
   {
      printf("\nThe system in the container : \n");
      fail = syscon_write_quaddobl_system();
      fail = quaddobl_row_reduce_system(1);
      printf("\nThe system after reduction : \n");
      fail = syscon_write_quaddobl_system();
   }
   return fail;
}
