/*  Get the inverse matrix of an complex matrix with LU-decomposition, test if A * A_inverse = I */
/*  where A is a square matrix and I is an identity matrix */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "dc_matrix.h"
#include "dc_inverse.h"

void Random_Test ( int n);
/* tests the computation of the inverse matrix on a random n-by-n matrix */

void Interactive_Test ( int n);
 /* tests the computation of the inverse matrix on a user-given n-by-n matrix */

int main(void)
{
   int ans,n;

   srand(time(NULL));

   printf("\nTesting inverse function for matrices of complex doubles.\n");

   for(;;)
   {
      printf("\nChoose one of the following options:\n");
      printf("  0. exit the loop; or\n");
      printf("  1. test on randomly generated matrices; or\n");
      printf("  2. test user-given linear system.\n");
      printf("Type 0, 1, or 2 to select : "); scanf("%d", &ans);
      if (ans == 0) break;
      if ((ans == 1) || (ans == 2))
      {
         printf("\nGive the dimension : ");
         scanf("%d", &n);
         if (ans == 1)
            Random_Test(n);
         else
            Interactive_Test(n);
      }
      else
         printf("Invalid choice.  Please try again...\n");
   }

   return 0;
}

void Interactive_Test ( int n )
{
   dcmplx a[n][n], a_inverse[n][n], I[n][n];

   printf("Reading a %d-by-%d matrix...\n", n, n);
   read_dcmatrix(n,n,a);
   printf("Writing a %d-by-%d matrix : \n", n, n);
   print_dcmatrix(n,n,a);

   copy_dcmatrix(n, n, a, a_inverse);
   dcinverse (n, a_inverse);
   multiply_dcmatrix(n, n, n, a, a_inverse, I);
   printf("The result  of a * a_inverse:\n");
   print_dcmatrix(n, n, I);
}

void Random_Test ( int n )
{
   dcmplx a[n][n], a_inverse[n][n], I[n][n];

   random_dcmatrix(n,n,a);
   printf("A random complex %d-by-%d matrix :\n", n, n);
   print_dcmatrix(n,n,a);

   copy_dcmatrix(n, n, a, a_inverse);
   dcinverse (n, a_inverse);
   multiply_dcmatrix(n, n, n, a, a_inverse, I);
   printf("The result  of a * a_inverse:\n");
   print_dcmatrix(n, n, I);
}


