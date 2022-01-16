/* reads a system and then computes a numerical irreducible decomposition */

#include <stdio.h>
#include "syscon.h"
#include "solcon.h"
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

int standard_polysys_solver ( void );
/*
 * DESCRIPTION :
 *   Calls the solver in standard double precision on a polynomial system. */

int standard_laursys_solver ( void );
/*
 * DESCRIPTION :
 *   Calls the solver in standard double precision on a Laurent system. */

int dobldobl_polysys_solver ( void );
/*
 * DESCRIPTION :
 *   Calls the solver in double double precision on a polynomial system. */

int dobldobl_laursys_solver ( void );
/*
 * DESCRIPTION :
 *   Calls the solver in double double precision on a Laurent system. */

int quaddobl_polysys_solver ( void );
/*
 * DESCRIPTION :
 *   Calls the solver in quad double precision on a polynomial system. */

int quaddobl_laursys_solver ( void );
/*
 * DESCRIPTION :
 *   Calls the solver in quad double precision on a Laurent system. */

int standard_polysys_write ( int topdim );
/*
 * DESCRIPTION :
 *   Writes the number of points for dimensions 0 to topdim,
 *   for the polynomial system solved in standard double precision. */

int standard_laursys_write ( int topdim );
/*
 * DESCRIPTION :
 *   Writes the number of points for dimensions 0 to topdim,
 *   for the Laurent polynomial system solved in standard double precision. */

int dobldobl_polysys_write ( int topdim );
/*
 * DESCRIPTION :
 *   Writes the number of points for dimensions 0 to topdim,
 *   for the polynomial system solved in double double precision. */

int dobldobl_laursys_write ( int topdim );
/*
 * DESCRIPTION :
 *   Writes the number of points for dimensions 0 to topdim,
 *   for the Laurent polynomial system solved in double double precision. */

int quaddobl_polysys_write ( int topdim );
/*
 * DESCRIPTION :
 *   Writes the number of points for dimensions 0 to topdim,
 *   for the polynomial system solved in quad double precision. */

int quaddobl_laursys_write ( int topdim );
/*
 * DESCRIPTION :
 *   Writes the number of points for dimensions 0 to topdim,
 *   for the Laurent polynomial system solved in quad double precision. */

int main ( int argc, char *argv[] )
{
   int fail,choice;

   adainit();

   printf("\nMENU for the precision and type of system :\n");
   printf("  0. standard double precision on a polynomial system\n");
   printf("  1. standard double precision on a Laurent polynomial system\n");
   printf("  2. double double precision on a polynomial system\n");
   printf("  3. double double precision on a Laurent polynomial system\n");
   printf("  4. quad double precision on a polynomial system\n");
   printf("  5. quad double precision on a Laurent polynomial system\n");
   printf("Type 0, 1, 2, 3, 4, or 5 : "); scanf("%d", &choice);

   if(choice == 0)
      fail = standard_polysys_solver();
   else if(choice == 1)
      fail = standard_laursys_solver();
   else if(choice == 2)
      fail = dobldobl_polysys_solver();
   else if(choice == 3)
      fail = dobldobl_laursys_solver();
   else if(choice == 4)
      fail = quaddobl_polysys_solver();
   else if(choice == 5)
      fail = quaddobl_laursys_solver();
   else
      printf("invalid choice\n");

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

int standard_polysys_solver ( void )
{
   int fail,nbq,nbtasks,topdim,filter,factor,verbose;

   fail = syscon_read_standard_system();
   fail = syscon_number_of_standard_polynomials(&nbq);
   printf("-> read %d polynomials\n", nbq);
   printf("\n");
   read_solver_options(&nbtasks,&topdim,&filter,&factor,&verbose);

   printf("\nCalling the solver ...\n\n");
   fail = standard_polysys_solve(nbtasks,topdim,filter,factor,verbose);
   fail = standard_polysys_write(topdim);

   return fail;
}

int standard_laursys_solver ( void )
{
   int fail,nbq,nbtasks,topdim,filter,factor,verbose;

   fail = syscon_read_standard_Laurent_system();
   fail = syscon_number_of_standard_Laurentials(&nbq);
   printf("-> read %d polynomials\n", nbq);
   printf("\n");
   read_solver_options(&nbtasks,&topdim,&filter,&factor,&verbose);

   printf("\nCalling the solver ...\n\n");
   fail = standard_laursys_solve(nbtasks,topdim,filter,factor,verbose);
   fail = standard_laursys_write(topdim);

   return fail;
}

int dobldobl_polysys_solver ( void )
{
   int fail,nbq,nbtasks,topdim,filter,factor,verbose;

   fail = syscon_read_dobldobl_system();
   fail = syscon_number_of_dobldobl_polynomials(&nbq);
   printf("-> read %d polynomials\n", nbq);
   printf("\n");
   read_solver_options(&nbtasks,&topdim,&filter,&factor,&verbose);

   printf("\nCalling the solver ...\n\n");
   fail = dobldobl_polysys_solve(nbtasks,topdim,filter,factor,verbose);
   fail = dobldobl_polysys_write(topdim);

   return fail;
}

int dobldobl_laursys_solver ( void )
{
   int fail,nbq,nbtasks,topdim,filter,factor,verbose;

   fail = syscon_read_dobldobl_Laurent_system();
   fail = syscon_number_of_dobldobl_Laurentials(&nbq);
   printf("-> read %d polynomials\n", nbq);
   printf("\n");
   read_solver_options(&nbtasks,&topdim,&filter,&factor,&verbose);

   printf("\nCalling the solver ...\n\n");
   fail = dobldobl_laursys_solve(nbtasks,topdim,filter,factor,verbose);
   fail = dobldobl_laursys_write(topdim);

   return fail;
}

int quaddobl_polysys_solver ( void )
{
   int fail,nbq,nbtasks,topdim,filter,factor,verbose;

   fail = syscon_read_quaddobl_system();
   fail = syscon_number_of_quaddobl_polynomials(&nbq);
   printf("-> read %d polynomials\n", nbq);
   printf("\n");
   read_solver_options(&nbtasks,&topdim,&filter,&factor,&verbose);

   printf("\nCalling the solver ...\n\n");
   fail = quaddobl_polysys_solve(nbtasks,topdim,filter,factor,verbose);
   fail = quaddobl_polysys_write(topdim);

   return fail;
}

int quaddobl_laursys_solver ( void )
{
   int fail,nbq,nbtasks,topdim,filter,factor,verbose;

   fail = syscon_read_quaddobl_Laurent_system();
   fail = syscon_number_of_quaddobl_Laurentials(&nbq);
   printf("-> read %d polynomials\n", nbq);
   printf("\n");
   read_solver_options(&nbtasks,&topdim,&filter,&factor,&verbose);

   printf("\nCalling the solver ...\n\n");
   fail = quaddobl_laursys_solve(nbtasks,topdim,filter,factor,verbose);
   fail = quaddobl_laursys_write(topdim);

   return fail;
}

int standard_polysys_write ( int topdim )
{
   int fail,dim,len;

   for(dim=0; dim<=topdim; dim++)
   {
      fail = copy_standard_polysys_witset(dim);
      fail = solcon_number_of_standard_solutions(&len);
      printf("-> number of points at dimension %d : %d\n",dim,len);
   }

   return fail;
}

int standard_laursys_write ( int topdim )
{
   int fail,dim,len;

   for(dim=0; dim<=topdim; dim++)
   {
      fail = copy_standard_laursys_witset(dim);
      fail = solcon_number_of_standard_solutions(&len);
      printf("-> number of points at dimension %d : %d\n",dim,len);
   }

   return fail;
}

int dobldobl_polysys_write ( int topdim )
{
   int fail,dim,len;

   for(dim=0; dim<=topdim; dim++)
   {
      fail = copy_dobldobl_polysys_witset(dim);
      fail = solcon_number_of_dobldobl_solutions(&len);
      printf("-> number of points at dimension %d : %d\n",dim,len);
   }

   return fail;
}

int dobldobl_laursys_write ( int topdim )
{
   int fail,dim,len;

   for(dim=0; dim<=topdim; dim++)
   {
      fail = copy_dobldobl_laursys_witset(dim);
      fail = solcon_number_of_dobldobl_solutions(&len);
      printf("-> number of points at dimension %d : %d\n",dim,len);
   }

   return fail;
}

int quaddobl_polysys_write ( int topdim )
{
   int fail,dim,len;

   for(dim=0; dim<=topdim; dim++)
   {
      fail = copy_quaddobl_polysys_witset(dim);
      fail = solcon_number_of_quaddobl_solutions(&len);
      printf("-> number of points at dimension %d : %d\n",dim,len);
   }

   return fail;
}

int quaddobl_laursys_write ( int topdim )
{
   int fail,dim,len;

   for(dim=0; dim<=topdim; dim++)
   {
      fail = copy_quaddobl_laursys_witset(dim);
      fail = solcon_number_of_quaddobl_solutions(&len);
      printf("-> number of points at dimension %d : %d\n",dim,len);
   }

   return fail;
}
