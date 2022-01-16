/* Tests the scaling of systems and solutions through the C interface. */

#include <stdio.h>
#include <stdlib.h>
#include "syscon.h"
#include "solcon.h"
#include "phcpack.h"
#include "scalers.h"

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
   printf("Type 0, 1, or 2 to make your choice : ");
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
   int fail,dim,i;

   fail = syscon_read_standard_system();
   if(fail == 0)
   {
      printf("\nThe system in the container : \n");
      fail = syscon_write_standard_system();
      fail = syscon_number_of_standard_polynomials(&dim);
      if(fail == 0)
      {
         double c[4*dim+2];
         fail = standard_scale_system(2,c);
         printf("The estimated inverse condition number : %.3e\n",c[4*dim]);
         printf("The scaling coefficients : \n");
         for(i=0; i<4*dim; i=i+2)
            printf("%2d : %.15e  %.15e\n",i/2,c[i],c[i+1]);
         printf("\nThe system in the container : \n");
         fail = syscon_write_standard_system();

      }
   }
   return fail;
}

int dobldobl_scaler ( void )
{
   int fail,dim,i;

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
         printf("The scaling coefficients : \n");
         for(i=0; i<8*dim; i=i+4)
            printf("%2d : %.15e %.15e  %.15e %.15e\n",
                   i/4,c[i],c[i+1],c[i+2],c[i+3]);
         printf("\nThe system in the container : \n");
         fail = syscon_write_dobldobl_system();

      }
   }
   return fail;
}

int quaddobl_scaler ( void )
{
   int fail,dim,i;

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
         printf("The scaling coefficients : \n");
         for(i=0; i<16*dim; i=i+4)
            printf("%2d : %.15e %.15e %.15e %.15e  %.15e %.15e %.15e %.15e\n",
                   i/8,c[i],c[i+1],c[i+2],c[i+3],c[i+4],c[i+5],c[i+6],c[i+7]);
         printf("\nThe system in the container : \n");
         fail = syscon_write_quaddobl_system();

      }
   }
   return fail;
}
