#include <stdio.h>

int phc_sols ( int n, int r, int m, double s[m] )

/* This C function takes as input r solution vectors of length n,
   for a total of m double coefficients. */

{
   int i,j,k=0;

   printf("Dimension of solution vectors : %d.\n", n);
   printf("Number of solutions : %d.\n", r);

   for(i=0; i<r; i++)
   {
      printf("solution %d : \n", i+1);
      for(j=0; j<n; j++)
      {
         printf("  %.15lf", s[k++]);
         printf("  %.15lf\n", s[k++]);
      }
   }

   return 0;
}
