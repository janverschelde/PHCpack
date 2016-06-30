/* Prompts the user for the level of precision, reads a system with solutions
   and then calls the functions to compute power series solutions. */ 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "phcpack.h"
#include "syscon.h"
#include "solcon.h"
#include "series.h"

void standard_test ( void ); /* test in standard double precision */

void dobldobl_test ( void ); /* test in double double precision */

void quaddobl_test ( void ); /* test in quad double precision */

void ask_options ( int *idx, int *nbr, int *verbose );
/*
 * DESCRIPTION :
 *   Prompts the user for the options of the power series method.
 *
 * ON RETURN :
 *   idx       index of the series parameter;
 *   nbr       number of steps with Newton's method;
 *   verbose   0 for false, 1 for true. */

int main ( int argc, char *argv[] )
{
   int choice;
   char ch;

   adainit();

   printf("\nMENU to select the working precision : \n");
   printf("  0. standard double precision;\n");
   printf("  1. double double precision;\n");
   printf("  2. quad double precision.\n");
   printf("Type 0, 1, or 2 to select the precision : ");
   scanf("%d",&choice);
   scanf("%c",&ch); /* skip newline symbol */

   if(choice == 0)
      standard_test();
   else if(choice == 1)
      dobldobl_test();
   else if(choice == 2)
      quaddobl_test();
   else
      printf("invalid selection, please try again\n");

   adafinal();

   return 0;
}

void standard_test ( void )
{
   int fail,idx,nbr,verbose;

   printf("\nNewton power series method in double precision ...\n");
   fail = read_standard_start_system();
   fail = copy_start_system_to_container();
   fail = copy_start_solutions_to_container();

   ask_options(&idx,&nbr,&verbose);
   fail = standard_Newton_series(idx,nbr,verbose);
}

void dobldobl_test ( void )
{
   int fail,idx,nbr,verbose;

   printf("\nNewton power series method in double double precision ...\n");
   fail = read_dobldobl_start_system();
   fail = copy_dobldobl_start_system_to_container();
   fail = copy_dobldobl_start_solutions_to_container();

   ask_options(&idx,&nbr,&verbose);
   fail = dobldobl_Newton_series(idx,nbr,verbose);
}

void quaddobl_test ( void )
{
   int fail,idx,nbr,verbose;

   printf("\nNewton power series method in quad double precision ...\n");
   fail = read_quaddobl_start_system();
   fail = copy_quaddobl_start_system_to_container();
   fail = copy_quaddobl_start_solutions_to_container();

   ask_options(&idx,&nbr,&verbose);
   fail = quaddobl_Newton_series(idx,nbr,verbose);
}

void ask_options ( int *idx, int *nbr, int *verbose )
{
   char ch;

   printf("\nReading the settings of the power series method ...\n");
   printf("  Give the index of the series parameter : "); scanf("%d",idx);
   printf("  Give the number of Newton steps : "); scanf("%d",nbr);
   printf("  Extra output during computations ? (y/n) ");
   scanf("%c",&ch); /* skip newline symbol */
   scanf("%c",&ch);
   *verbose = (ch == 'y' ? 1 : 0);
   scanf("%c",&ch); /* skip newline symbol */
}
