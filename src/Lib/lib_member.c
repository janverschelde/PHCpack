/* Tests the membership test for a witness set and test point. */

#include <stdio.h>
#include "syscon.h"
#include "solcon.h"
#include "witset.h"

#define verbose 1 /* verbose flag */

int standard_membership_test ( void );
/*
 * DESCRIPTION :
 *   Prompts for a witness set and a test point.
 *   Runs the membership test in standard double precision. */

int dobldobl_membership_test ( void );
/*
 * DESCRIPTION :
 *   Prompts for a witness set and a test point.
 *   Runs the membership test in double double precision. */

int quaddobl_membership_test ( void );
/*
 * DESCRIPTION :
 *   Prompts for a witness set and a test point.
 *   Runs the membership test in quad double precision. */

int main ( int argc, char *argv[] )
{
   char ans;
   int precision,fail;

   adainit();

   printf("\nMENU for the working precision :\n");
   printf("  0. standard double precision;\n");
   printf("  1. double double precision;\n");
   printf("  2. quad double precision.\n");
   printf("Type 0, 1, or 2 to make a choice : ");
   scanf("%d",&precision);
   scanf("%c",&ans); /* skip end of line character */

   if(precision == 0)
      fail = standard_membership_test();
   else if(precision == 1)
      fail = dobldobl_membership_test();
   else if(precision == 2)
      fail = quaddobl_membership_test();
   else
      printf("Selected precision level is not supported.\n");

   adafinal();

   return 0;
}

int standard_membership_test ( void )
{
   int fail,n,dim,deg;
   char ans;

   printf("\nReading a witness set ...\n");
   fail = read_witness_set(&n,&dim,&deg);

   if(verbose>0)  /* only in verbose mode */
   {
      printf("\nThe ambient dimension : %d.\n",n);
      printf("The dimension of the solution set : %d.\n",dim);
      printf("The degree of the solution set : %d.\n",deg);
      printf("\nDo you want to see the embedded system ? (y/n) ");
      scanf("%c",&ans);
      if(ans == 'y')
      {
         printf("\nThe system read :\n");
         fail = syscon_write_standard_system();
      }
      scanf("%c",&ans); /* skip end of line character */
      printf("\nDo you want to see the solutions ? (y/n) ");
      scanf("%c",&ans);
      if(ans == 'y')
      {
         printf("\nThe solutions read :\n");
         fail = solcon_write_standard_solutions();
      }
      scanf("%c",&ans); /* skip end of line character */
   }

   return 0;
}

int dobldobl_membership_test ( void )
{
   int fail,n,dim,deg;
   char ans;

   printf("\nReading a witness set ...\n");
   fail = read_dobldobl_witness_set(&n,&dim,&deg);

   if(verbose>0)  /* only in verbose mode */
   {
      printf("\nThe ambient dimension : %d.\n",n);
      printf("The dimension of the solution set : %d.\n",dim);
      printf("The degree of the solution set : %d.\n",deg);
      printf("\nDo you want to see the embedded system ? (y/n) ");
      scanf("%c",&ans);
      if(ans == 'y')
      {
         printf("\nThe system read :\n");
         fail = syscon_write_dobldobl_system();
      }
      scanf("%c",&ans); /* skip end of line character */
      printf("\nDo you want to see the solutions ? (y/n) ");
      scanf("%c",&ans);
      if(ans == 'y')
      {
         printf("\nThe solutions read :\n");
         fail = solcon_write_dobldobl_solutions();
      }
      scanf("%c",&ans); /* skip end of line character */
   }

   return 0;
}

int quaddobl_membership_test ( void )
{
   int fail,n,dim,deg;
   char ans;

   printf("\nReading a witness set ...\n");
   fail = read_quaddobl_witness_set(&n,&dim,&deg);

   if(verbose>0)  /* only in verbose mode */
   {
      printf("\nThe ambient dimension : %d.\n",n);
      printf("The dimension of the solution set : %d.\n",dim);
      printf("The degree of the solution set : %d.\n",deg);
      printf("\nDo you want to see the embedded system ? (y/n) ");
      scanf("%c",&ans);
      if(ans == 'y')
      {
         printf("\nThe system read :\n");
         fail = syscon_write_quaddobl_system();
      }
      scanf("%c",&ans); /* skip end of line character */
      printf("\nDo you want to see the solutions ? (y/n) ");
      scanf("%c",&ans);
      if(ans == 'y')
      {
         printf("\nThe solutions read :\n");
         fail = solcon_write_quaddobl_solutions();
      }
      scanf("%c",&ans); /* skip end of line character */
   }

   return 0;
}
