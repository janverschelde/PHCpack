/* tests the C interface to the Littlewood-Richardson homotopies */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "schubert.h"

extern void adainit();
extern void adafinal();

int call_Littlewood_Richardson_homotopies
 ( int n, int k, int c, int *brackets, int verbose, int *r );
/*
 * DESCRIPTION :
 *   Tests the call to the Littlewood-Richardson homotopies.
 *
 * ON ENTRY :
 *   n         dimension of the ambient space;
 *   k         dimension of the solution planes;
 *   c         number of intersection conditions;
 *   brackets  contains c*k integer numbers with the intersection conditions;
 *   verbose   1 if intermediate output is needed, 0 to be silent.
 *
 * ON RETURN :
 *   r         the formal root count.  */

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
      fail = call_Littlewood_Richardson_homotopies(n,k,c,brackets,verbose,&r);

   printf("The formal root count : %d\n",r);

   adafinal();

   return 0;
}

int call_Littlewood_Richardson_homotopies
 ( int n, int k, int c, int *brackets, int verbose, int *r )
{
   const int size = 2*(c-2)*n*n;
   double flags[size]; /* real + imaginary parts stored rowwise */
   char ans,filename[80];
   int i,fail,nbname;

   scanf("%c",&ans); /* skip newline symbol */
   printf("Give the name of the output file : ");
   scanf("%s",filename);
   scanf("%c",&ans); /* skip newline symbol */

   nbname = strlen(filename);
   printf("Number of characters in %s is %d\n",filename,nbname);

   fail = Littlewood_Richardson_homotopies
            (n,k,c,brackets,verbose,nbname,filename,r,flags);

   printf("\nThe coefficients of the fixed flag :");
   for(i=0; i<size; i++)
   {
      if(i % (2*n) == 0) printf("\n");
      printf(" %+.1e", flags[i]);
   }
   printf("\n");
   return fail;
}
