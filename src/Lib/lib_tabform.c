/* simple test program in C on the tableau format of a polynomial system */

#include <stdio.h>
#include <stdlib.h>
#include "syscon.h"

void standard_tableau_form ( void );
/*
 * DESCRIPTION :
 *   Sets up the data structures for the tableau form of the system
 *   in the container with standard double precision coefficients. */

void write_standard_tableau_form ( int n );
/*
 * DESCRIPTION :
 *   Writes the tableau form of the system in the container
 *   with standard double precision coefficients.
 *   The dimension of the system is given in n. */

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

void write_standard_tableau_form ( int n )
{
   int i,j,k,nt,fail;
   double c[2];
   int d[n];

   for(i=1; i<=n; i++)
   {
      fail = syscon_number_of_standard_terms(i,&nt);
      printf("  #terms in polynomial %d : %d\n", i, nt);

      for(j=1; j<=nt; j++)
      {
         fail = syscon_retrieve_standard_term(i,j,n,d,c);
         printf("%22.15e  %22.15e", c[0], c[1]);
         for (k=0; k<n; k++) printf(" %d", d[k]);
         printf("\n");
      }
   }
}

void standard_tableau_form ( void )
{
   int fail,neq,nbt,idx,nbtsum,cffidx,expidx,idxtrm,k;
   double cff[2];
   double* coefficients;
   int* nbterms;
   int* exponents;

   fail = syscon_number_of_standard_polynomials(&neq);

   nbterms = (int*)calloc(neq,sizeof(int));
   
   nbtsum = 0;
   for(idx=1; idx <= neq; idx++)
   {
      fail = syscon_number_of_standard_terms(idx,&nbt);
      printf("  #terms in polynomial %d : %d\n", idx, nbt);
      nbterms[idx-1] = nbt;
      nbtsum = nbtsum + nbt;
   }
   printf("Total number of terms : %d = %d", nbtsum, nbterms[0]);
   for(idx=1; idx<neq; idx++) printf(" + %d", nbterms[idx]);
   printf("\n");
   coefficients = (double*)calloc(2*nbtsum,sizeof(double));
   exponents = (int*)calloc(neq*nbtsum,sizeof(int));
   {
      int wrk[neq];

      cffidx = 0; expidx = 0;
      printf("%d\n", neq);
      for(idx=0; idx < neq; idx++)
      {
         printf("%d\n", nbterms[idx]);
         for(idxtrm=1; idxtrm <= nbterms[idx]; idxtrm++)
         {
            fail = syscon_retrieve_standard_term(idx+1,idxtrm,neq,wrk,cff);
            coefficients[cffidx++] = cff[0];
            coefficients[cffidx++] = cff[1];
            printf("%22.15e  %22.15e", cff[0], cff[1]);
            for (k=0; k < neq; k++) printf(" %d", wrk[k]);
            for (k=0; k < neq; k++) exponents[expidx++] = wrk[k];
            printf("\n");
         }
      }
   }
   printf("The tableau format :\n");
   cffidx = 0; expidx = 0;
   printf("%d\n", neq);
   for(idx=0; idx < neq; idx++)
   {
      printf("%d\n", nbterms[idx]);
      for(idxtrm=1; idxtrm <= nbterms[idx]; idxtrm++)
      {
         printf("%22.15e  %22.15e",
                coefficients[cffidx], coefficients[cffidx+1]);
         cffidx = cffidx + 2;
         for (k=0; k < neq; k++) printf(" %d", exponents[expidx++]);
         printf("\n");
      }
   }
}

void test_standard_container ( void )
{
   int n,fail,*d;
   double *c;

   fail = syscon_read_standard_system();
   fail = syscon_write_standard_system();
   fail = syscon_number_of_standard_polynomials(&n);

   printf("\nThe dimension of the system : %d\n",n);

   write_standard_tableau_form(n);
   standard_tableau_form();
}
