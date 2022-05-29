/* test on some of the operations in witset */

#include <stdio.h>
#include "syscon.h"
#include "solcon.h"
#include "witset.h"

int main ( int argc, char *argv[] )
{
   int i,n,fail;
   const int max = 80;
   char p[max];
   char nl;

   adainit();

   printf("\nMaking a witness set for a polynomial ...\n");
   printf("  give the number of variables : ");
   scanf("%d",&n); scanf("%c",&nl); /* skip new line */
   printf("  give a polynomial (end with \';\'): ");
   for(i=0; i<max; i++)
   {
      scanf("%c",&p[i]);
      if(p[i] == ';')
      {
         p[++i] = '\0';
         break;
      }
   }
   printf("Your polynomial in %d variables is %s\n",n,p);

   fail = standard_witset_of_hypersurface(n,i+1,p);

   printf("\nThe embedded system :\n");
   syscon_write_standard_system();
   printf("\nThe witness points :\n");
   solcon_write_standard_solutions();

   adafinal();

   return 0;
}
