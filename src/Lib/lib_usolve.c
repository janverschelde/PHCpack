/* simple test on univariate root finding */

#include <stdio.h>
#include <stdlib.h>
#include "syscon.h"
#include "solcon.h"
#include "unisolvers.h"

void standard_test ( void );
/*
 * Prompts the user for a polynomial and then computes the roots
 * in standard double precision. */

void dobldobl_test ( void );
/*
 * Prompts the user for a polynomial and then computes the roots
 * in double double precision. */

void quaddobl_test ( void );
/*
 * Prompts the user for a polynomial and then computes the roots
 * in quad double precision. */

int main(void)
{
   int choice;

   printf("\nMENU for testing univariate root finding :\n");
   printf("  0. in standard double precision;\n");
   printf("  1. in double double precision;\n");
   printf("  2. in quad double precision.\n");
   printf("Type 0, 1, or 2 to select the precision : ");
   scanf("%d",&choice);

   adainit();

   if(choice == 0)
      standard_test();
   else if(choice == 1)
      dobldobl_test();
   else if(choice == 2)
      quaddobl_test();
   else
      printf("invalid choice, please try again...\n");

   adafinal();

   return 0;
}

void standard_test ( void )
{
   int fail,max,nit;
   const double eps = 1.0e-12;

   fail = syscon_read_system();
   fail = syscon_write_system();
   
   printf("\nGive the maximum number of iterations : ");
   scanf("%d",&max);
   
   fail = solve_with_standard_doubles(max,eps,&nit);

   printf("\nNumber of iterations : %d\n",nit);
   fail = solcon_write_solutions();
}

void dobldobl_test ( void )
{
   int fail,max,nit;
   const double eps = 1.0e-24;

   fail = syscon_read_dobldobl_system();
   fail = syscon_write_dobldobl_system();
   
   printf("\nGive the maximum number of iterations : ");
   scanf("%d",&max);
   
   fail = solve_with_double_doubles(max,eps,&nit);

   printf("\nNumber of iterations : %d\n",nit);
   fail = solcon_write_dobldobl_solutions();
}

void quaddobl_test ( void )
{
   int fail,max,nit;
   const double eps = 1.0e-48;

   fail = syscon_read_quaddobl_system();
   fail = syscon_write_quaddobl_system();
   
   printf("\nGive the maximum number of iterations : ");
   scanf("%d",&max);
   
   fail = solve_with_quad_doubles(max,eps,&nit);

   printf("\nNumber of iterations : %d\n",nit);
   fail = solcon_write_quaddobl_solutions();
}
