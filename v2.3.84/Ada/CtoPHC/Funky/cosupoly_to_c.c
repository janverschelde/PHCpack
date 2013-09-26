#include <stdio.h>

int cosupoly_to_c ( int n, int ns, int *s, int nc, double *c )
/* This C function takes as input the support and coefficients of
   a complex polynomial in several variables. */
{
   int i;

   printf("Entering the C function ...\n");

   printf("The number of variables : %d.\n", n);
   printf("Length of support : %d.\n", ns); 
   printf("The support : ");
   for(i=0; i<ns; i++) printf(" %d", s[i]);  
   printf(".\n");
   printf("Number of coefficients : %d.\n", nc);
   printf("The coefficients :\n");
   for(i=0; i<nc; i++) printf(" %.15lf\n", c[i]);

   printf("... leaving the C function.\n");

   return 0;
}
