#include <stdio.h>

int path_sols ( int n, int nbsols, int mul[nbsols],
                int nbcff, double cff[nbcff] )

/* This C function takes as input the solutions at the end of paths,
   obtained as output of the Ada path trackers.

   ON ENTRY :
     n          dimension of the solution vectors;
     nbsols     number of solution vectors;
     mul        multiplicities of each solution vector;
     nbcff      number of coefficients in the solution vectors;
     cff        all coefficients in the solution vectors.
*/

{
   int i;

   printf("Dimension of the solution vectors : %d.\n", n);
   printf("Number of solution maps : %d.\n", nbsols); 
   printf("The multiplicities : ");
   for(i=0; i<nbsols; i++)
      printf(" %d", mul[i]);  
   printf(".\n");
   printf("Number of coefficients : %d.\n", nbcff);
   printf("The coefficients :\n");
   for(i=0; i<nbcff; i++)
      printf(" %.15le\n", cff[i]);

   return 0;
}
