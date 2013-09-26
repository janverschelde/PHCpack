/* testing linear-product root counts and random linear-product systems */

#include <stdio.h>
#include <stdlib.h>
#include "phcpack.h"
#include "product.h"

int compute_root_count ( int *r );
/*
 * DESCRIPTION :
 *   prompts the user for a polynomial system, constructs a supporting
 *   set structure, and returns the corresponding root count in r. */ 

int construct_start_system ( int r );
/*
 * DESCRPTION :
 *   constructs a random linear-product system based on the supporting
 *   set structure and with root count in r.
 *
 * REQUIRED : compute_root_count() was executed. */

int main ( void )
{
   int r;

   printf("\nTesting linear-product root counts and systems...\n");

   adainit();

   compute_root_count(&r);
   construct_start_system(r);
   clear_set_structure();

   adafinal();

   return 0;
}

int compute_root_count ( int *r )
{
   int fail;

   fail = syscon_read_system();
   printf("\nThe system in the container : \n");
   fail = syscon_write_system();
   fail = supporting_set_structure();
   printf("\nA supporting set structure : \n");
   fail = write_set_structure();
   fail = linear_product_root_count(r);
   printf("\nThe linear-product root count : %d\n",*r);

   return fail;
}

int construct_start_system ( int r )
{
   int fail;
   int nbsols;

   fail = random_linear_product_system();
   printf("\nA random linear-product system :\n");
   fail = syscon_write_system();
   fail = solve_linear_product_system();
   printf("\nThe solutions : \n");
   fail = solcon_write_solutions();
   fail = solcon_number_of_solutions(&nbsols);
   if(r == nbsols)
      printf("\nComputed %d solutions, as many as the root count.\n",nbsols);
   else
      printf("\nNumber of solutions computed %d /= %d, the root count ?!!\n",
             nbsols,r);

   return fail;
}
