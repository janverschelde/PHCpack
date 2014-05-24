/* tests the C interface to the Littlewood-Richardson homotopies */

#include <stdio.h>
#include <stdlib.h>
#include "schubert.h"

extern void adainit();
extern void adafinal();

int main ( int argc, char *argv[] )
{
   int fail,n,k,c,i,j,r;
   int *brackets;
   int verbose=0;

   printf("\nresolving a general Schubert intersection condition ...\n");
   printf("  give the ambient dimension of the space : "); scanf("%d",&n);
   printf("  give the dimension of the solution planes : "); scanf("%d",&k);
   printf("  give the number of intersection conditions : "); scanf("%d",&c);
   printf("\nreading %d conditions on %d-planes in %d-space ...\n",c,k,n);

   brackets = (int*) calloc(c*k, sizeof(int));
   for(i=0; i<c; i++)
   {
      printf("  give %d integers for condition %d : ",k,i+1);
      for(j=0; j<k; j++) scanf("%d",&brackets[i*k+j]);
   }
   printf("\nthe intersection conditions :\n");
   for(i=0; i<c; i++)
   {
      printf("[%d", brackets[i*k]);
      for(j=1; j<k; j++) printf(" %d",brackets[i*k+j]);
      printf("]");
   }
   printf("\n");

   adainit();

   fail = resolve_Schubert_conditions(n,k,c,brackets,verbose,&r);
   printf("The formal root count : %d\n",r);

   adafinal();

   return 0;
}
