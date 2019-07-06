/* simple test program in C on the tableau format of a polynomial system */

#include <stdio.h>
#include <stdlib.h>
#include "syscon.h"

void standard_tableau_form ( int n );
/*
 * DESCRIPTION :
 *   Writes the tableau form of the system in container
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

void standard_tableau_form ( int n )
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

void test_standard_container ( void )
{
   int n,fail,*d;
   double *c;

   fail = syscon_read_standard_system();
   fail = syscon_write_standard_system();
   fail = syscon_number_of_standard_polynomials(&n);

   printf("\nThe dimension of the system : %d\n",n);

   standard_tableau_form(n);
}
