/* tests the C interface to the Littlewood-Richardson homotopies */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "schubert.h"

extern void adainit();
extern void adafinal();

int call_Littlewood_Richardson_homotopies
 ( int n, int k, int c, int *brackets, int verbose, int verify, int size,
   int *r );
/*
 * DESCRIPTION :
 *   Tests the call to the Littlewood-Richardson homotopies.
 *
 * ON ENTRY :
 *   n         dimension of the ambient space;
 *   k         dimension of the solution planes;
 *   c         number of intersection conditions;
 *   brackets  contains c*k integer numbers with the intersection conditions;
 *   verify    1 if diagnostic verification is needed, 0 if no;
 *   verbose   1 if intermediate output is needed, 0 to be silent;
 *   size      number of doubles in the flags, depends on the precision.
 *
 * ON RETURN :
 *   r         the formal root count.  */

int main ( int argc, char *argv[] )
{
   int fail,n,k,c,i,j,r;
   int *brackets;
   int verbose=0;
   int verify=0;
   char ans;
   int precision,size;

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
      printf("\nMENU for the precision : \n");
      printf("  0. standard double precision; \n");
      printf("  1. double double precision; \n");
      printf("  2. quad double precision; \n");
      printf("Type 0, 1, or 2 to select the precision : ");
      scanf("%d",&precision);

      if(precision == 0) size = 2*(c-2)*n*n;
      if(precision == 1) size = 4*(c-2)*n*n;
      if(precision == 2) size = 8*(c-2)*n*n;

      fail = call_Littlewood_Richardson_homotopies
               (n,k,c,brackets,verbose,verify,size,&r);
   }

   printf("The formal root count : %d\n",r);

   adafinal();

   return 0;
}

int call_Littlewood_Richardson_homotopies
 ( int n, int k, int c, int *brackets, int verbose, int verify,
   int size, int *r )
{
   double flags[size]; /* real + imaginary parts stored rowwise */
   char ans,filename[80];
   int i,fail,nbname,prec;

   scanf("%c",&ans); /* skip newline symbol */
   printf("Give the name of the output file : ");
   scanf("%s",filename);
   /* scanf("%c",&ans);  skip newline symbol */

   nbname = strlen(filename);
   printf("Number of characters in %s is %d\n",filename,nbname);

   if(prec == 0)
      fail = standard_Littlewood_Richardson_homotopies
                (n,k,c,brackets,verbose,verify,nbname,filename,r,flags);
   if(prec == 1)
      fail = dobldobl_Littlewood_Richardson_homotopies
                (n,k,c,brackets,verbose,verify,nbname,filename,r,flags);
   if(prec == 2)
      fail = quaddobl_Littlewood_Richardson_homotopies
                (n,k,c,brackets,verbose,verify,nbname,filename,r,flags);

   printf("\nThe coefficients of the fixed flag :");
   for(i=0; i<size; i++)
   {
      if(i % (2*n) == 0) printf("\n");
      printf(" %+.1e", flags[i]);
   }
   printf("\n");
   return fail;
}
