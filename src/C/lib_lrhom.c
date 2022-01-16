/* tests the C interface to the Littlewood-Richardson homotopies */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "syscon.h"
#include "solcon.h"
#include "schubert.h"

extern void adainit();
extern void adafinal();

void prompt_for_options
 ( int *verify, int *minrep, int *tosquare );
/*
 * DESCRIPTION :
 *   Prompts the user for the options of the Littlewood-Richardson run.
 *
 * ON RETURN :
 *   verify    1 if diagnostic verification is needed, 0 if no;
 *   minrep    1 for a minimal problem formulation, 0 for all minors;
 *   tosquare  1 to square overdetermined systems, 0 for Gauss-Newton. */

int call_Littlewood_Richardson_homotopies
 ( int n, int k, int c, int *brackets, int precision,
   int verbose, int verify, int minrep, int tosquare,
   int size, int *r );
/*
 * DESCRIPTION :
 *   Tests the call to the Littlewood-Richardson homotopies.
 *
 * ON ENTRY :
 *   n         dimension of the ambient space;
 *   k         dimension of the solution planes;
 *   c         number of intersection conditions;
 *   brackets  contains c*k integer numbers with the intersection conditions;
 *   precision 0 for standard double, 1 for double double, 2 for quad double;
 *   verify    1 if diagnostic verification is needed, 0 if no;
 *   verbose   1 if intermediate output is needed, 0 to be silent;
 *   minrep    1 for a minimal problem formulation, 0 for all minors;
 *   tosquare  1 to square overdetermined systems, 0 for Gauss-Newton;
 *   size      number of doubles in the flags, depends on the precision.
 *
 * ON RETURN :
 *   r         the formal root count.  */

void write_results ( int prec );
/*
 * DESCRIPTION :
 *   Writes the contents of the systems and solutions containers,
 *   for precision level 0, 1, or 2. */

int main ( int argc, char *argv[] )
{
   int fail,n,k,c,i,j,r;
   int *brackets;
   int verbose=0;
   int verify=0;
   int minrep=1;
   int tosquare=0;
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

      prompt_for_options(&verify,&minrep,&tosquare);

      fail = call_Littlewood_Richardson_homotopies
               (n,k,c,brackets,precision,verbose,verify,minrep,tosquare,
                size,&r);
   }

   printf("The formal root count : %d\n",r);

   adafinal();

   return 0;
}

void prompt_for_options
 ( int *verify, int *minrep, int *tosquare )
{
   char ans;

   printf("\nVerify the intersection conditions ? (y/n) ");
   scanf("%c",&ans); /* skip new line character */
   scanf("%c",&ans);
   if(ans == 'y')
      *verify = 1;
   else
      *verify = 0;

   printf("Use a minimal problem formulation ? (y/n) ");
   scanf("%c",&ans); /* skip new line character */
   scanf("%c",&ans);
   if(ans == 'y')
      *minrep = 1;
   else
      *minrep = 0;

   printf("Square the overdetermined systems ? (y/n) ");
   scanf("%c",&ans); /* skip new line character */
   scanf("%c",&ans);
   if(ans == 'y')
      *tosquare = 1;
   else
      *tosquare = 0;
}

int call_Littlewood_Richardson_homotopies
 ( int n, int k, int c, int *brackets, int precision,
   int verbose, int verify, int minrep, int tosquare,
   int size, int *r )
{
   double flags[size]; /* real + imaginary parts stored rowwise */
   char ans,filename[80];
   int i,fail,nbname;

   scanf("%c",&ans); /* skip newline symbol */

   if(verify == 1)
   {
      ans = 'y';
      printf("\n");
   }
   else
   {
      printf("\nOutput to separate file ? (y/n) ? ");
      scanf("%c",&ans);
      scanf("%c",&ans); /* skip newline symbol */
   }
   if(ans != 'y')
   {
      nbname = 0;
      filename[0] = '\0';  
   }
   else
   {
      printf("Give the name of the output file : ");
      scanf("%s",filename);
      /* scanf("%c",&ans);  skip newline symbol */
      nbname = strlen(filename);
   }
   printf("Number of characters in %s is %d\n",filename,nbname);

   if(precision == 0)
      fail = standard_Littlewood_Richardson_homotopies
                (n,k,c,brackets,verbose,verify,minrep,tosquare,
                 nbname,filename,r,flags);
   if(precision == 1)
      fail = dobldobl_Littlewood_Richardson_homotopies
                (n,k,c,brackets,verbose,verify,minrep,tosquare,
                 nbname,filename,r,flags);
   if(precision == 2)
      fail = quaddobl_Littlewood_Richardson_homotopies
                (n,k,c,brackets,verbose,verify,minrep,tosquare,
                 nbname,filename,r,flags);

   printf("\nThe coefficients of the fixed flag :");
   for(i=0; i<size; i++)
   {
      if(i % (2*n) == 0) printf("\n");
      printf(" %+.1e", flags[i]);
   }
   printf("\n");

   if(nbname > 0)
      printf("See the file %s for the results.\n",filename);
   else
      write_results(precision);

   return fail;
}

void write_results ( int prec ) 
{
   if(prec==0)
   {
      printf("THE SYSTEM SOLVED :\n");
      syscon_write_standard_system();
      printf("THE SOLUTIONS : \n");
      solcon_write_standard_solutions();
   }
   if(prec==1)
   {
      printf("THE SYSTEM SOLVED :\n");
      syscon_write_dobldobl_system();
      printf("THE SOLUTIONS : \n");
      solcon_write_dobldobl_solutions();
   }
   if(prec==2)
   {
      printf("THE SYSTEM SOLVED :\n");
      syscon_write_quaddobl_system();
      printf("THE SOLUTIONS : \n");
      solcon_write_quaddobl_solutions();
   }
}
