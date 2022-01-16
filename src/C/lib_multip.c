/* tests the computation of the multiplicity structure */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "phcpack.h"
#include "syscon.h"
#include "solcon.h"
#include "multiplicity.h"

void test_standard_multiplicity_structure ( int order );
/* 
 * DESCRIPTION :
 *   Tests multiplicity structure computation with standard arithmetic,
 *   where order is the maximum differentation order. */

void test_dobldobl_multiplicity_structure ( int order );
/* 
 * DESCRIPTION :
 *   tests multiplicity structure computation with double double arithmetic,
 *   where order is the maximum differentiation order. */

void test_quaddobl_multiplicity_structure ( int order );
/*
 * DESCRIPTION :
 *   tests multiplicity structure computation with quad double arithmetic,
 *   where order is the maximum differentiation order. */

int main ( int argc, char *argv[] )
{
   int choice,order;
   char ch;

   adainit();

   printf("\nMENU to compute the multiplicity structure : \n");
   printf("  1. test with standard double arithmetic;\n");
   printf("  2. test with double double arithmetic;\n");
   printf("  3. test with quad double arithmetic.\n");
   printf("Type 0, 1, 2, or 3 to select : ");
   scanf("%d",&choice);
   scanf("%c",&ch); /* skip newline symbol */

   printf("Give the maximum differentiation order : ");
   scanf("%d",&order);

   if(choice == 1)
      test_standard_multiplicity_structure(order);
   else if(choice == 2)
      test_dobldobl_multiplicity_structure(order);
   else if(choice == 3)
      test_quaddobl_multiplicity_structure(order);
   else
      printf("invalid selection, please try again\n");

   adafinal();

   return 0;
}

void test_standard_multiplicity_structure ( int order )
{
   const double tol = 1.0e-8;
   int fail,dim,len,mult;
   int hilb[order+1];

   printf("\nTesting with standard arithmetic ...\n");
   fail = syscon_read_standard_system();
   fail = syscon_number_of_standard_polynomials(&dim);
   printf("The system container has %d polynomials.\n",dim);
   fail = solcon_read_standard_solutions();
   fail = solcon_number_of_standard_solutions(&len);
   printf("The solution container has size %d.\n",len);
   fail = solcon_dimension_of_standard_solutions(&dim);
   printf("The solutions in the container have dimension %d.\n",dim);
   fail = standard_multiplicity_structure(order,1,tol,&mult,hilb);
   printf("The multiplicity : %d.\n",mult);
   printf("Values of the Hilbert function :");
   int k;
   for(k=0; k<order+1; k++)
      printf(" %d",hilb[k]);
   printf("\n");
}

void test_dobldobl_multiplicity_structure ( int order )
{
   const double tol = 1.0e-8;
   int fail,dim,len,mult;
   int hilb[order+1];

   printf("\nTesting with dobldobl arithmetic ...\n");
   fail = syscon_read_dobldobl_system();
   fail = syscon_number_of_dobldobl_polynomials(&dim);
   printf("The system container has %d polynomials.\n",dim);
   fail = solcon_read_dobldobl_solutions();
   fail = solcon_number_of_dobldobl_solutions(&len);
   printf("The solution container has size %d.\n",len);
   fail = solcon_dimension_of_dobldobl_solutions(&dim);
   printf("The solutions in the container have dimension %d.\n",dim);
   fail = dobldobl_multiplicity_structure(order,1,tol,&mult,hilb);
   printf("The multiplicity : %d.\n",mult);
   printf("Values of the Hilbert function :");
   int k;
   for(k=0; k<order+1; k++)
      printf(" %d",hilb[k]);
   printf("\n");
}

void test_quaddobl_multiplicity_structure ( int order )
{
   const double tol = 1.0e-8;
   int fail,dim,len,mult;
   int hilb[order+1];

   printf("\nTesting with quaddobl arithmetic ...\n");
   fail = syscon_read_quaddobl_system();
   fail = syscon_number_of_quaddobl_polynomials(&dim);
   printf("The system container has %d polynomials.\n",dim);
   fail = solcon_read_quaddobl_solutions();
   fail = solcon_number_of_quaddobl_solutions(&len);
   printf("The solution container has size %d.\n",len);
   fail = solcon_dimension_of_quaddobl_solutions(&dim);
   printf("The solutions in the container have dimension %d.\n",dim);
   fail = quaddobl_multiplicity_structure(order,1,tol,&mult,hilb);
   printf("The multiplicity : %d.\n",mult);
   printf("Values of the Hilbert function :");
   int k;
   for(k=0; k<order+1; k++)
      printf(" %d",hilb[k]);
   printf("\n");
}
