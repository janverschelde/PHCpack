/* Prompts the user for the level of precision, reads a system with solutions
   and then calls the functions to compute power series solutions. */ 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "phcpack.h"
#include "syscon.h"
#include "syspool.h"
#include "series.h"

void standard_query_pool ( void );
/*
 * Interactive function to query the pool for systems
 * in standard double precision. */

void dobldobl_query_pool ( void );
/*
 * Interactive function to query the pool for systems
 * in double double precision. */

void quaddobl_query_pool ( void );
/*
 * Interactive function to query the pool for systems
 * in quad double precision. */

void standard_series_test ( void );
/* 
 * Tests Newton's method for power series in standard double precision. */

void dobldobl_series_test ( void );
/*
 * Tests Newton's method for power series in double double precision. */

void quaddobl_series_test ( void );
/*
 * Tests Newton's method for power series in quad double precision. */

void standard_Pade_test ( void );
/* 
 * Tests the Pade approximant constructor in standard double precision. */

void dobldobl_Pade_test ( void );
/*
 * Tests the Pade approximant constructor in double double precision. */

void quaddobl_Pade_test ( void );
/*
 * Tests the Pade approximant constructor in quad double precision. */

void ask_options_for_series ( int *idx, int *maxdeg, int *nbr, int *verbose );
/*
 * DESCRIPTION :
 *   Prompts the user for the options of the power series method.
 *
 * ON RETURN :
 *   idx       index of the series parameter;
 *   nbr       number of steps with Newton's method;
 *   maxdeg    maximal degree of the series;
 *   verbose   0 for false, 1 for true. */

void ask_options_for_Pade
 ( int *idx, int *numdeg, int *dendeg, int *nbr, int *verbose );
/*
 * DESCRIPTION :
 *   Prompts the user for the options of the Pade approximant. 
 *
 * ON RETURN :
 *   idx       index of the series parameter;
 *   numdeg    degree of the numerator;
 *   dendeg    degree of the denominator;
 *   nbr       number of steps with Newton's method;
 *   verbose   0 for false, 1 for true. */

int main ( int argc, char *argv[] )
{
   int choice;
   char ch;

   adainit();

   printf("\nMENU to select the test and the working precision : \n");
   printf("  0. power series in standard double precision;\n");
   printf("  1. power series in double double precision;\n");
   printf("  2. power series in quad double precision.\n");
   printf("  3. Pade approximant in standard double precision;\n");
   printf("  4. Pade approximant in double double precision;\n");
   printf("  5. Pade approximant in quad double precision.\n");
   printf("Type 0, 1, 2, 3, 4, or 5 to select the precision : ");
   scanf("%d",&choice);
   scanf("%c",&ch); /* skip newline symbol */

   if(choice == 0)
      standard_series_test();
   else if(choice == 1)
      dobldobl_series_test();
   else if(choice == 2)
      quaddobl_series_test();
   else if(choice == 3)
      standard_Pade_test();
   else if(choice == 4)
      dobldobl_Pade_test();
   else if(choice == 5)
      quaddobl_Pade_test();
   else
      printf("invalid selection, please try again\n");

   adafinal();

   return 0;
}

void standard_query_pool ( void )
{
   int nbr,idx;

   syspool_standard_size(&nbr);
   printf("  Computed %d solution systems.\n",nbr);
   while(1)
   {
      printf("\nGive index for the system you want to see : ");
      scanf("%d",&idx);
      if(idx <= 0) break;
      printf("Copying system %d to container ...\n",idx);
      syspool_copy_to_standard_container(idx);
      printf("The system in the container :\n");
      syscon_write_standard_system();
   }
}

void dobldobl_query_pool ( void )
{
   int nbr,idx;

   syspool_dobldobl_size(&nbr);
   printf("  Computed %d solution systems.\n",nbr);
   while(1)
   {
      printf("\nGive index for the system you want to see (0 to exit) : ");
      scanf("%d",&idx);
      if(idx <= 0) break;
      printf("Copying system %d to container ...\n",idx);
      syspool_copy_to_dobldobl_container(idx);
      printf("The system in the container :\n");
      syscon_write_dobldobl_system();
   }
}

void quaddobl_query_pool ( void )
{
   int nbr,idx;

   syspool_quaddobl_size(&nbr);
   printf("  Computed %d solution systems.\n",nbr);
   while(1)
   {
      printf("\nGive index for the system you want to see : ");
      scanf("%d",&idx);
      if(idx <= 0) break;
      printf("Copying system %d to container ...\n",idx);
      syspool_copy_to_quaddobl_container(idx);
      printf("The system in the container :\n");
      syscon_write_quaddobl_system();
   }
}

void standard_series_test ( void )
{
   int fail,idx,maxdeg,nbr,verbose;

   printf("\nNewton power series method in double precision ...\n");
   fail = read_standard_start_system();
   fail = copy_start_system_to_container();
   fail = copy_start_solutions_to_container();

   ask_options_for_series(&idx,&maxdeg,&nbr,&verbose);
   fail = standard_Newton_series(idx,maxdeg,nbr,verbose);

   printf("\nDone with Newton's method.");
   standard_query_pool();
}

void dobldobl_series_test ( void )
{
   int fail,idx,maxdeg,nbr,verbose;

   printf("\nNewton power series method in double double precision ...\n");
   fail = read_dobldobl_start_system();
   fail = copy_dobldobl_start_system_to_container();
   fail = copy_dobldobl_start_solutions_to_container();

   ask_options_for_series(&idx,&maxdeg,&nbr,&verbose);
   fail = dobldobl_Newton_series(idx,maxdeg,nbr,verbose);

   printf("\nDone with Newton's method.");
   dobldobl_query_pool();
}

void quaddobl_series_test ( void )
{
   int fail,idx,maxdeg,nbr,verbose;

   printf("\nNewton power series method in quad double precision ...\n");
   fail = read_quaddobl_start_system();
   fail = copy_quaddobl_start_system_to_container();
   fail = copy_quaddobl_start_solutions_to_container();

   ask_options_for_series(&idx,&maxdeg,&nbr,&verbose);
   fail = quaddobl_Newton_series(idx,maxdeg,nbr,verbose);

   printf("\nDone with Newton's method.");
   quaddobl_query_pool();
}

void ask_options_for_series ( int *idx, int *maxdeg, int *nbr, int *verbose )
{
   char ch;

   printf("\nReading the settings of the power series method ...\n");
   printf("  Give the index of the series parameter : "); scanf("%d",idx);
   printf("  Give the maximal degree of the series : "); scanf("%d",maxdeg);
   printf("  Give the number of Newton steps : "); scanf("%d",nbr);
   printf("  Extra output during computations ? (y/n) ");
   scanf("%c",&ch); /* skip newline symbol */
   scanf("%c",&ch);
   *verbose = (ch == 'y' ? 1 : 0);
   scanf("%c",&ch); /* skip newline symbol */
}

void ask_options_for_Pade
 ( int *idx, int *numdeg, int *dendeg, int *nbr, int *verbose )
{
   char ch;

   printf("\nReading the settings of the Pade approximant constructor ...\n");
   printf("  Give the index of the series parameter : "); scanf("%d",idx);
   printf("  Give the degree of the numerator : "); scanf("%d",numdeg);
   printf("  Give the degree of the denominator : "); scanf("%d",dendeg);
   printf("  Give the number of Newton steps : "); scanf("%d",nbr);
   printf("  Extra output during computations ? (y/n) ");
   scanf("%c",&ch); /* skip newline symbol */
   scanf("%c",&ch);
   *verbose = (ch == 'y' ? 1 : 0);
   scanf("%c",&ch); /* skip newline symbol */
}

void standard_Pade_test ( void )
{
   int fail,idx,numdeg,dendeg,nbr,verbose;

   printf("\nPade approximant creator in double precision ...\n");
   fail = read_standard_start_system();
   fail = copy_start_system_to_container();
   fail = copy_start_solutions_to_container();

   ask_options_for_Pade(&idx,&numdeg,&dendeg,&nbr,&verbose);
   fail = standard_Pade_approximant(idx,numdeg,dendeg,nbr,verbose);

   printf("\nDone with the Pade approximant constructor.");
   standard_query_pool();
}

void dobldobl_Pade_test ( void )
{
   int fail,idx,numdeg,dendeg,nbr,verbose;

   printf("\nPade approximant creator in double double precision ...\n");
   fail = read_dobldobl_start_system();
   fail = copy_dobldobl_start_system_to_container();
   fail = copy_dobldobl_start_solutions_to_container();

   ask_options_for_Pade(&idx,&numdeg,&dendeg,&nbr,&verbose);
   fail = dobldobl_Pade_approximant(idx,numdeg,dendeg,nbr,verbose);

   printf("\nDone with the Pade approximant constructor.");
   dobldobl_query_pool();
}

void quaddobl_Pade_test ( void )
{
   int fail,idx,numdeg,dendeg,nbr,verbose;

   printf("\nPade approximant creator in quad double precision ...\n");
   fail = read_quaddobl_start_system();
   fail = copy_quaddobl_start_system_to_container();
   fail = copy_quaddobl_start_solutions_to_container();

   ask_options_for_Pade(&idx,&numdeg,&dendeg,&nbr,&verbose);
   fail = quaddobl_Pade_approximant(idx,numdeg,dendeg,nbr,verbose);

   printf("\nDone with the Pade approximant constructor.");
   quaddobl_query_pool();
}
