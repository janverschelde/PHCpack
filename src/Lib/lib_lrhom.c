/* tests the C interface to the Littlewood-Richardson homotopies */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "schubert.h"

extern void adainit();
extern void adafinal();

int main ( int argc, char *argv[] )
{
   int fail,n,k,c,i,j,r;
   int *brackets;
   int verbose=0;
   char ans;

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

   scanf("%c",&ans); /* read the new line character */
   printf("Run the Littlewood-Richardson homotopies ? (y/n) ");
   scanf("%c",&ans);

   adainit();

   if(ans != 'y')
      fail = resolve_Schubert_conditions(n,k,c,brackets,verbose,&r);
   else
   {
      double flags[2*(c-2)*n*n]; /* real + imaginary parts stored rowwise */
      char filename[80];
      int nbname;

      scanf("%c",&ans); /* skip newline symbol */
      printf("Give the name of the output file : ");
      scanf("%s",filename);
      scanf("%c",&ans); /* skip newline symbol */

      nbname = strlen(filename);
      printf("Number of characters in %s is %d\n",filename,nbname);

      fail = Littlewood_Richardson_homotopies
               (n,k,c,brackets,verbose,nbname,filename,&r,flags);
   }
   printf("The formal root count : %d\n",r);

   adafinal();

   return 0;
}
