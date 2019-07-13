/* simple test program in C on the tableau format of a polynomial system */

#include <stdio.h>
#include <stdlib.h>
#include "syscon.h"
#include "tabform.h"

int make_standard_tableau_form ( int verbose );
/*
 * DESCRIPTION :
 *   Sets up the data structures for the tableau form of the system
 *   in the container with standard double precision coefficients.
 *   If verbose, then extra information is written to screen.
 *   Returns the failure code of the retrieval operations. */

int retrieve_standard_tableau_form ( int verbose );
/*
 * DESCRIPTION :
 *   Tests the retrieval of the system in the container in tableau format. */

int test_standard_container ( void );
/*
 * DESCRIPTION :
 *   Prompts the user for a polynomial system,
 *   computes and writes its tableau form. */

int main ( void )
{
   adainit();
   int fail = test_standard_container();
   adafinal();

   return 0;
}

int make_standard_tableau_form ( int verbose )
{
   int fail,neq,nvr,nbtsum;
   int *nbterms;
   int *exponents;
   double *coefficients;

   fail = syscon_number_of_standard_polynomials(&neq);
   if(fail != 0) return fail;
   fail = syscon_number_of_symbols(&nvr);
   if(fail != 0) return fail;

   printf("\nThe number of equations in the system : %d\n",neq);
   printf("The number of variables in the system : %d\n\n",nvr);

   nbterms = (int*)calloc(neq,sizeof(int));

   fail = number_of_standard_terms(neq,nbterms,&nbtsum,1);
   
   coefficients = (double*)calloc(2*nbtsum,sizeof(double));
   exponents = (int*)calloc(nvr*nbtsum,sizeof(int));

   fail = standard_tableau_form
             (neq,nvr,nbterms,coefficients,exponents,verbose);
   if(fail == 0)
   {
      printf("The tableau format :\n");
      write_standard_tableau_form(neq,nvr,nbterms,coefficients,exponents);

      fail = syscon_clear_standard_system();

      store_standard_tableau_form
         (neq,nvr,nbterms,coefficients,exponents,verbose);

      printf("After storing the system defined by the tableau format :\n");
      fail = syscon_write_standard_system();
   }
   return fail;
}

int retrieve_standard_tableau_form ( int verbose )
{
   int fail,neq,nvr,nbt;

   fail = load_standard_tableau_dimensions(&neq,&nvr,&nbt);

   if(verbose > 0)
   {
      printf("-> the number of equations : %d\n", neq);
      printf("-> the number of variables : %d\n", nvr);
      printf("-> total number of terms : %d\n", nbt);
   }
   {
      int nbterms[neq];
      int exponents[nvr*nbt];
      double coefficients[2*nbt];
   }
   return fail;
}

int test_standard_container ( void )
{
   int fail,verbose = 0;
   char ans;

   printf("Verbose ? (y/n) ");
   ans = getchar();
   if(ans == 'y') verbose = 1;

   fail = syscon_read_standard_system();
   if(fail == 0) fail = syscon_write_standard_system();

   if(fail == 0) fail = make_standard_tableau_form(verbose);

   if(fail == 0) fail = retrieve_standard_tableau_form(verbose);

   return fail;
}
