/* Tests the funcions defined in sweep.h */

#include <stdio.h>
#include <stdlib.h>
#include "sweep.h"

int main ( int argc, char *argv[] )
{
   int fail,nbequ,nbvar,nbpar,nbr;

   adainit();

   printf("Give the number of equations : "); scanf("%d",&nbequ);
   printf("Give the number of variables : "); scanf("%d",&nbvar);
   printf("Give the number of parameters : "); scanf("%d",&nbpar);
   {
      int pars[nbpar],k;
      for(k=0; k<nbpar; k++)
      {
         printf("Give the index of parameter %d : ",k);
         scanf("%d",&pars[k]);
      }
      sweep_define_parameters_numerically(nbequ,nbvar,nbpar,pars);
   }
   fail = sweep_get_number_of_equations(&nbr);
   printf("The number of equations : %d.\n",nbr);
   fail = sweep_get_number_of_variables(&nbr);
   printf("The number of variables : %d.\n",nbr);
   fail = sweep_get_number_of_parameters(&nbr);
   printf("The number of parameters : %d.\n",nbr);
   {
      int idxpars[nbpar],k;
      fail = sweep_get_indices_numerically(idxpars);
      printf("THe indices of the parameters :");
      for(k=0; k<nbpar; k++) printf(" %d",idxpars[k]);
      printf("\n");
   }

   adafinal();

   return fail;
}
