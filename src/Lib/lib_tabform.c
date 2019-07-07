/* simple test program in C on the tableau format of a polynomial system */

#include <stdio.h>
#include <stdlib.h>
#include "syscon.h"

int total_number_of_standard_terms ( void );
/*
 * DESCRIPTION :
 *   Returns the total number of terms in the system stored in the
 *   container for systems with standard double precision coefficients. */

void standard_tableau_form
 ( int neq, int nvr, int *nbterms, double *coefficients, int *exponents,
   int verbose );
/*
 * DESCRIPTION :
 *   Given allocated data structures, defines the coefficients and exponents
 *   of the tableau form of the system stored in the container for systems
 *   with standard double precision coefficients.
 *
 * ON ENTRY :
 *   neq            number of equations of the system;
 *   nvr            number of variables of the system;
 *   nbterms        array of size neq, nbterms[k] contains the number
 *                  of terms in the (k+1)-th equation in the system;
 *   coefficients   allocated for 2 times the total number of terms;
 *   exponents      allocated for nvr times the total number of terms;
 *   verbose        if > 0, then the tableau form is written,
 *                  if = 0, then the function remains silent.
 *
 * ON RETURN :
 *   coefficients   coefficients of the system;
 *   exponents      exponents of the system. */

void write_standard_tableau_form
 ( int neq, int nvr, int *nbterms, double *coefficients, int *exponents );
/*
 * DESCRIPTION :
 *   Writes the tableau form of the system to screen.
 *
 * ON ENTRY :
 *   neq            number of equations of the system;
 *   nvr            number of variables of the system;
 *   nbterms        array of size neq, nbterms[k] contains the number
 *                  of terms in the (k+1)-th equation in the system;
 *   coefficients   coefficients of the system;
 *   exponents      exponents of the system. */

void make_standard_tableau_form ( void );
/*
 * DESCRIPTION :
 *   Sets up the data structures for the tableau form of the system
 *   in the container with standard double precision coefficients. */

void test_standard_container ( void );
/*
 * DESCRIPTION :
 *   Prompts the user for a polynomial system
 *   and writes its tableau form. */

int main ( void )
{
   adainit();
   test_standard_container();
   adafinal();

   return 0;
}

int total_number_of_standard_terms ( void )
{
   int result = 0;
   int fail,neq,idx,nbt;

   fail = syscon_number_of_standard_polynomials(&neq);

   for(idx=1; idx <= neq; idx++)
   {
      fail = syscon_number_of_standard_terms(idx,&nbt);
      result = result + nbt;
   }
   return result;
}

void standard_tableau_form
 ( int neq, int nvr, int *nbterms, double *coefficients, int *exponents,
   int verbose )
{
   int wrk[nvr];  // work space to hold the exponents of a term
   double cff[2]; // holds real and imaginary part of a coefficient
   int cffidx = 0;
   int expidx = 0;
   int fail,idx,idxtrm,k;

   if(verbose > 0) printf("%d\n", neq);
   for(idx=0; idx < neq; idx++)
   {
      if(verbose > 0) printf("%d\n", nbterms[idx]);
      for(idxtrm=1; idxtrm <= nbterms[idx]; idxtrm++)
      {
         fail = syscon_retrieve_standard_term(idx+1,idxtrm,neq,wrk,cff);
         coefficients[cffidx++] = cff[0];
         coefficients[cffidx++] = cff[1];
         if(verbose > 0)
         {
            printf("%22.15e  %22.15e", cff[0], cff[1]);
            for (k=0; k < nvr; k++) printf(" %d", wrk[k]);
            printf("\n");
         }
         for (k=0; k < nvr; k++) exponents[expidx++] = wrk[k];
      }
   }
}

void write_standard_tableau_form
 ( int neq, int nvr, int *nbterms, double *coefficients, int *exponents )
{
   int cffidx = 0;
   int expidx = 0;
   int idx,idxtrm,k;

   printf("%d\n", neq);
   for(idx=0; idx < neq; idx++)
   {
      printf("%d\n", nbterms[idx]);
      for(idxtrm=1; idxtrm <= nbterms[idx]; idxtrm++)
      {
         printf("%22.15e  %22.15e",
                coefficients[cffidx], coefficients[cffidx+1]);
         cffidx = cffidx + 2;
         for (k=0; k < nvr; k++) printf(" %d", exponents[expidx++]);
         printf("\n");
      }
   }
}

void make_standard_tableau_form ( void )
{
   int fail,neq,nvr,nbt,idx,nbtsum;
   int *nbterms;
   int *exponents;
   double *coefficients;

   fail = syscon_number_of_standard_polynomials(&neq);
   fail = syscon_number_of_symbols(&nvr);

   printf("\nThe number of equations in the system : %d\n",neq);
   printf("The number of variables in the system : %d\n\n",nvr);

   nbtsum = total_number_of_standard_terms();
   nbterms = (int*)calloc(neq,sizeof(int));
   
   for(idx=1; idx <= neq; idx++)
   {
      fail = syscon_number_of_standard_terms(idx,&nbt);
      printf("  #terms in polynomial %d : %d\n", idx, nbt);
      nbterms[idx-1] = nbt;
   }
   printf("Total number of terms : %d = %d", nbtsum, nbterms[0]);
   for(idx=1; idx<neq; idx++) printf(" + %d", nbterms[idx]);
   printf("\n");
   coefficients = (double*)calloc(2*nbtsum,sizeof(double));
   exponents = (int*)calloc(nvr*nbtsum,sizeof(int));

   standard_tableau_form(neq,nvr,nbterms,coefficients,exponents,1);

   printf("The tableau format :\n");
   write_standard_tableau_form(neq,nvr,nbterms,coefficients,exponents);
}

void test_standard_container ( void )
{
   int n,fail,*d;
   double *c;

   fail = syscon_read_standard_system();
   fail = syscon_write_standard_system();

   make_standard_tableau_form();
}
