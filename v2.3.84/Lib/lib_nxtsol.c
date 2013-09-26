/* reads a target and start system and then solves the target system,
   using the start system in an artificial-parameter homotopy */

#include <stdio.h>
#include "solcon.h"
#include "phcpack.h"
#include "jump_track.h"
#include "next_track.h"

int prompt_for_precision ( void );
/*
 * DESCRIPTION :
 *   Prompts the user for the level of precision and returns
 *   0 for standard double precision, 1 for double double precision,
 *   and 2 for quad double precision. */

int call_initialize_standard_homotopy ( int *index );
/*
 * DESCRIPTION :
 *   Prepares the containers to initialize the homotopy to track a path
 *   in standard double precision.  Returns in index the number of the
 *   solution in the container selected as start solution. */

int call_initialize_dobldobl_homotopy ( int *index );
/*
 * DESCRIPTION :
 *   Prepares the containers to initialize the homotopy to track a path
 *   in double double precision.  Returns in index the number of the
 *   solution in the container selected as start solution. */

int call_initialize_quaddobl_homotopy ( int *index );
/*
 * DESCRIPTION :
 *   Prepares the containers to initialize the homotopy to track a path
 *   in quad double precision.  Returns in index the number of the
 *   solution in the container selected as start solution. */

int write_standard_solution ( int index );
/*
 * DESCRIPTION :
 *   Writes the solution in the standard solutions container
 *   at position equal to the value of the given index to screen. */

int write_dobldobl_solution ( int index );
/*
 * DESCRIPTION :
 *   Writes the solution in the double double solutions container
 *   at position equal to the value of the given index to screen. */

int write_quaddobl_solution ( int index );
/*
 * DESCRIPTION :
 *   Writes the solution in the quad double solutions container
 *   at position equal to the value of the given index to screen. */

int call_standard_path_tracker ( int index );
/*
 * DESCRIPTION :
 *   Calls the path tracker in standard double precision,
 *   starting a solution with the numbers in the index. */

int call_dobldobl_path_tracker ( int index );
/*
 * DESCRIPTION :
 *   Calls the path tracker in double double precision,
 *   starting a solution with the numbers in the index. */

int call_quaddobl_path_tracker ( int index );
/*
 * DESCRIPTION :
 *   Calls the path tracker in quad double precision,
 *   starting a solution with the numbers in the index. */

int main ( int argc, char *argv[] )
{
   int fail,nbsol;

   adainit();

   int level = prompt_for_precision();
   if(level == 0)
   {
      fail = call_initialize_standard_homotopy(&nbsol);
      fail = call_standard_path_tracker(nbsol);
      fail = clear_standard_tracker();
   }
   else if(level == 1)
   {
      fail = call_initialize_dobldobl_homotopy(&nbsol);
      fail = call_dobldobl_path_tracker(nbsol);
      fail = clear_dobldobl_tracker();
   }
   else
   {
      fail = call_initialize_quaddobl_homotopy(&nbsol);
      fail = call_quaddobl_path_tracker(nbsol);
      fail = clear_quaddobl_tracker();
   }

   adafinal();

   return 0;
}

int prompt_for_precision ( void )
{
   int answer = 0;
   char nlc;

   printf("\n");
   printf("Welcome to the path tracking with generators ...\n");
   printf("  0. run in standard double precision arithmetic;\n");
   printf("  1. run in double double precision arithmetic;\n");
   printf("  2. run in quad double precision arithmetic.\n");
   printf("Type 0, 1, or 2 to select precision : ");
   scanf("%d",&answer);
   scanf("%c",&nlc);     /* skip new line symbol */

   return answer;
}

int call_initialize_standard_homotopy ( int *index )
{
   int fail,len;

   fail = read_target_system_without_solutions();
   fail = read_start_system();
   fail = copy_start_solutions_to_container();
   fail = solcon_number_of_solutions(&len);
   printf("number of start solutions : %d\n",len);
   printf("-> give index of solution : "); scanf("%d",index);
   fail = initialize_standard_homotopy();
   fail = initialize_standard_solution(*index);

   return fail;
}

int call_initialize_dobldobl_homotopy ( int *index )
{
   int fail,len;

   fail = read_dobldobl_target_system(); /* no _without_solutions ! */
   fail = read_dobldobl_start_system();
   fail = copy_dobldobl_start_solutions_to_container();
   fail = solcon_number_of_dobldobl_solutions(&len);
   printf("number of start solutions : %d\n",len);
   printf("-> give index of solution : "); scanf("%d",index);
   fail = initialize_dobldobl_homotopy();
   fail = initialize_dobldobl_solution(*index);

   return fail;
}

int call_initialize_quaddobl_homotopy ( int *index )
{
   int fail,len;

   fail = read_quaddobl_target_system(); /* no _without_solutions ! */
   fail = read_quaddobl_start_system();
   fail = copy_quaddobl_start_solutions_to_container();
   fail = solcon_number_of_quaddobl_solutions(&len);
   printf("number of start solutions : %d\n",len);
   printf("-> give index of solution : "); scanf("%d",index);
   fail = initialize_quaddobl_homotopy();
   fail = initialize_quaddobl_solution(*index);

   return fail;
}

int write_standard_solution ( int index )
{
   int fail,nb;

   fail = solcon_length_solution_string(index,&nb);
   {
      char solution[nb];

      fail = solcon_write_solution_string(index,nb,solution);
      printf("\nsolution %d :\n%s\n",index,solution);
   }

   return fail;
}

int write_dobldobl_solution ( int index )
{
   int fail,nb;

   fail = solcon_length_dobldobl_solution_string(index,&nb);
   {
      char solution[nb];

      fail = solcon_write_dobldobl_solution_string(index,nb,solution);
      printf("\nsolution %d :\n%s\n",index,solution);
   }

   return fail;
}

int write_quaddobl_solution ( int index )
{
   int fail,nb;

   fail = solcon_length_quaddobl_solution_string(index,&nb);
   {
      char solution[nb];

      fail = solcon_write_quaddobl_solution_string(index,nb,solution);
      printf("\nsolution %d :\n%s\n",index,solution);
   }

   return fail;
}

int call_standard_path_tracker ( int index )
{
   int fail;
   char answer;

   fail = write_standard_solution(index);
   do
   {
      fail = next_standard_solution(index);
      fail = write_standard_solution(index);
      printf("Continue to next step ? (y/n) ");
      scanf("%c",&answer); /* get trailing new line...*/
      scanf("%c",&answer);
   }
   while(answer == 'y');

   return fail;
}

int call_dobldobl_path_tracker ( int index )
{
   int fail;
   char answer;

   fail = write_dobldobl_solution(index);
   do
   {
      fail = next_dobldobl_solution(index);
      fail = write_dobldobl_solution(index);
      printf("Continue to next step ? (y/n) ");
      scanf("%c",&answer); /* get trailing new line...*/
      scanf("%c",&answer);
   }
   while(answer == 'y');

   return fail;
}

int call_quaddobl_path_tracker ( int index )
{
   int fail;
   char answer;

   fail = write_quaddobl_solution(index);
   do
   {
      fail = next_quaddobl_solution(index);
      fail = write_quaddobl_solution(index);
      printf("Continue to next step ? (y/n) ");
      scanf("%c",&answer); /* get trailing new line...*/
      scanf("%c",&answer);
   }
   while(answer == 'y');

   return fail;
}
