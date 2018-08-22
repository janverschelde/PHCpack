/* reads a system and then computes a numerical irreducible decomposition */

#include <stdio.h>
#include "witsols.h"

void read_solver_options
 ( int *nbtasks, int *topdim, int *filter, int *factor, int *verbose );
/*
 * DESCRIPTION :
 *   Prompts the user for the number of tasks, top dimension,
 *   asks if the witness supersets need filtering and factoring,
 *   and sets the verbose flag.
 *
 * ON RETURN :
 *   nbtasks   equals the number of tasks for multitasking,
 *   topdim    the top dimension to start the homotopy cascades,
 *   filter    0 or 1 flag to filter the witness supersets, 
 *   factor    0 or 1 flag to factor the witness sets,
 *   verbose   for intermediate output. */

int main ( int argc, char *argv[] )
{
   int nbtasks,topdim,filter,factor,verbose,fail;

   adainit();

   read_solver_options(&nbtasks,&topdim,&filter,&factor,&verbose);

   fail = standard_polysys_solve(nbtasks,topdim,filter,factor,verbose);

   adafinal();

   return 0;
}

void read_solver_options
 ( int *nbtasks, int *topdim, int *filter, int *factor, int *verbose )
{
   printf("Give the number of tasks : ");
   scanf("%d", nbtasks);

   printf("Give the top dimension : ");
   scanf("%d", topdim);

   printf("Filter witness supersets ? (1 = yes, 0 = no) ");
   scanf("%d", filter);
   if(*filter == 1)
   {
      printf("Factor witness sets ? (1 = yes, 0 = no) ");
      scanf("%d", factor);
   }
   else
      *factor = 0;

   printf("Verbose mode ? (1 = yes, 0 = no) ");
   scanf("%d", verbose);
}
