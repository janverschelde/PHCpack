/* simple test program in C on the tableau format of a polynomial system */

#include <stdio.h>
#include <stdlib.h>
#include "syscon.h"
#include "tabform.h"

int make_standard_tableau_form ( void );
/*
 * DESCRIPTION :
 *   Sets up the data structures for the tableau form of the system
 *   in the container with standard double precision coefficients.
 *   Returns the failure code of the retrieval operations. */

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

int make_standard_tableau_form ( void )
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

   fail = standard_tableau_form(neq,nvr,nbterms,coefficients,exponents,1);
   if(fail == 0)
   {
      printf("The tableau format :\n");
      write_standard_tableau_form(neq,nvr,nbterms,coefficients,exponents);
   }
   return fail;
}

int test_standard_container ( void )
{
   int n,fail,*d;
   double *c;

   fail = syscon_read_standard_system();
   if(fail == 0) fail = syscon_write_standard_system();

   if(fail == 0) fail = make_standard_tableau_form();

   return fail;
}
