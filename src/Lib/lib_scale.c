/* Tests the scaling of systems and solutions through the C interface. */

#include <stdio.h>
#include <stdlib.h>
#include "syscon.h"
#include "solcon.h"
#include "phcpack.h"

int greetings( void );
/*
 * DESCRIPTION :
 *   Displays the version string of PHCpack. */

int standard_scaler ( void );
/*
 * DESCRIPTION :
 *   Prompts the user for a system and performs the scaling
 *   with standard double precision arithmetic. */

int dobldobl_scaler ( void );
/*
 * DESCRIPTION :
 *   Prompts the user for a system and performs the scaling
 *   with double double precision arithmetic. */

int quaddobl_scaler ( void );
/*
 * DESCRIPTION :
 *   Prompts the user for a system and performs the scaling
 *   with quad double precision arithmetic. */

int main ( int argc, char *argv[] )
{
   int fail,choice;

   adainit();

   fail = greetings();

   printf("\nMENU for the precision of the scalers :\n"); 
   printf("  0. use standard double precision arithmetic; or\n");
   printf("  1. use double double precision arithmetic; or\n");
   printf("  2. use quad double precision arithmetic.\n");
   printf("Type 1, 2, or 3 to make your choice : ");
   scanf("%d",&choice);

   if(choice == 0)
      fail = standard_scaler();
   else if(choice == 1)
      fail = dobldobl_scaler();
   else if(choice == 2)
      fail = quaddobl_scaler();

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

int standard_scaler ( void )
{
   int fail,dim;

   fail = syscon_read_system();
   if(fail == 0)
   {
      printf("\nThe system in the container : \n");
      fail = syscon_write_system();
      fail = syscon_number_of_polynomials(&dim);
      if(fail == 0)
      {
         double c[4*dim+2];
         fail = standard_scale_system(2,c);
         printf("The estimated inverse condition number : %.3e\n",c[4*dim]);
         printf("\nThe system in the container : \n");
         fail = syscon_write_system();

      }
   }
   return fail;
}

int dobldobl_scaler ( void )
{
   int fail,dim;

   fail = syscon_read_dobldobl_system();
   if(fail == 0)
   {
      printf("\nThe system in the container : \n");
      fail = syscon_write_dobldobl_system();
      fail = syscon_number_of_dobldobl_polynomials(&dim);
      if(fail == 0)
      {
         double c[8*dim+4];
         fail = dobldobl_scale_system(2,c);
         printf("The estimated inverse condition number : %.3e\n",c[8*dim]);
         printf("\nThe system in the container : \n");
         fail = syscon_write_dobldobl_system();

      }
   }
   return fail;
}

int quaddobl_scaler ( void )
{
   int fail,dim;

   fail = syscon_read_quaddobl_system();
   if(fail == 0)
   {
      printf("\nThe system in the container : \n");
      fail = syscon_write_quaddobl_system();
      fail = syscon_number_of_quaddobl_polynomials(&dim);
      if(fail == 0)
      {
         double c[16*dim+8];
         fail = quaddobl_scale_system(2,c);
         printf("The estimated inverse condition number : %.3e\n",c[16*dim]);
         printf("\nThe system in the container : \n");
         fail = syscon_write_quaddobl_system();

      }
   }
   return fail;
}
