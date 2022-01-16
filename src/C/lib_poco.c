/* reads a target and start system and then solves the target system,
   using the start system in an artificial-parameter homotopy */

#include <stdio.h>
#include <stdlib.h>
#include "phcpack.h"

int run_standard_continuation ( void );
/* runs in standard floating-point arithmetic */

int run_dobldobl_continuation ( void );
/* runs in double double arithmetic */

int run_quaddobl_continuation ( void );
/* runs in quad double arithmetic */

int run_multprec_continuation ( void );
/* runs in multiprecision arithmetic */

int run_standard_Laurent_continuation ( void );
/* solves Laurent system in standard floating-point arithmetic */

int run_dobldobl_Laurent_continuation ( void );
/* solves Laurent system in double double arithmetic */

int run_quaddobl_Laurent_continuation ( void );
/* solves Laurent system in quad double arithmetic */

int main ( int argc, char *argv[] )
{
   int fail,choice;
   char ch;

   adainit();

   printf("\nMENU for running increment-and-fix continuation :\n");
   printf("  1. run in standard double arithmetic;\n");
   printf("  2. run in double double arithmetic;\n");
   printf("  3. run in quad double arithmetic.\n");
   printf("  4. run in multiprecision arithmetic;\n");
   printf("  5. solve a Laurent system in standard double precision;\n");
   printf("  6. solve a Laurent system in double double precision;\n");
   printf("  7. solve a Laurent system in quad double precision.\n");
   printf("Type 1, 2, 3, 4, 5, 6, or 7 to select a run : ");
   scanf("%d",&choice);
   scanf("%c",&ch);     /* skip new line */

   if(choice == 1)
      fail = run_standard_continuation();
   else if(choice == 2)
      fail = run_dobldobl_continuation();
   else if(choice == 3)
      fail = run_quaddobl_continuation();
   else if(choice == 4)
      fail = run_multprec_continuation();
   else if(choice == 5)
      fail = run_standard_Laurent_continuation();
   else if(choice == 6)
      fail = run_dobldobl_Laurent_continuation();
   else if(choice == 7)
      fail = run_quaddobl_Laurent_continuation();
   else
      printf("Invalid choice, please try again.\n");

   adafinal();

   return 0;
}

int run_standard_continuation ( void )
{
   int fail, nbtasks;

   printf("\nCalling the standard path trackers in PHCpack...\n");
   fail = read_standard_target_system();
   fail = read_standard_start_system();
   fail = define_output_file();
   fail = write_standard_target_system();
   fail = write_standard_start_system();
   fail = write_start_solutions();
   printf("\n");
   fail = tune_continuation_parameters();
   printf("\n");
   fail = determine_output_during_continuation();
   printf("\nGive the number of tasks (0 for no multitasking ) : ");
   scanf("%d", &nbtasks);
   printf("\nSee the output file for results ...\n");
   fail = solve_by_standard_homotopy_continuation(nbtasks);
   // fail = write_target_solutions(); // root refiner

   return fail;
}

int run_dobldobl_continuation ( void )
{
   int fail, nbtasks;

   printf("\nCalling the double double path trackers in PHCpack...\n");
   fail = read_dobldobl_target_system();
   fail = read_dobldobl_start_system();
   fail = define_output_file();
   fail = write_dobldobl_target_system();
   fail = write_dobldobl_start_system();
   fail = write_dobldobl_start_solutions();
   printf("\n");
   fail = autotune_continuation_parameters(0,14); // (0, 24) is too severe
   fail = tune_continuation_parameters();
   printf("\n");
   fail = determine_output_during_continuation();
   printf("\nGive the number of tasks (0 for no multitasking ) : ");
   scanf("%d", &nbtasks);
   printf("\nSee the output file for results ...\n");
   fail = solve_by_dobldobl_homotopy_continuation(nbtasks);
   fail = write_dobldobl_target_solutions();

   return fail;
}

int run_quaddobl_continuation ( void )
{
   int fail, nbtasks;

   printf("\nCalling the quad double path trackers in PHCpack...\n");
   fail = read_quaddobl_target_system();
   fail = read_quaddobl_start_system();
   fail = define_output_file();
   fail = write_quaddobl_target_system();
   fail = write_quaddobl_start_system();
   fail = write_quaddobl_start_solutions();
   printf("\n");
   fail = autotune_continuation_parameters(0,14); // (0, 48) is too severe
   fail = tune_continuation_parameters();
   printf("\n");
   fail = determine_output_during_continuation();
   printf("\nGive the number of tasks (0 for no multitasking ) : ");
   scanf("%d", &nbtasks);
   printf("\nSee the output file for results ...\n");
   fail = solve_by_quaddobl_homotopy_continuation(nbtasks);
   fail = write_quaddobl_target_solutions();

   return fail;
}

int run_multprec_continuation ( void )
{
   int fail,deci;
   char ch;

   printf("\nCalling the multiprecision path trackers in PHCpack...\n");

   printf("-> give the number of decimal places : ");
   scanf("%d", &deci);
   scanf("%c",&ch);     /* skip new line */

   fail = read_multprec_target_system(deci);
   fail = read_multprec_start_system(deci);
   fail = define_output_file();
   fail = write_multprec_target_system();
   fail = write_multprec_start_system();
   fail = write_multprec_start_solutions();
   printf("\n");
   fail = autotune_continuation_parameters(0,deci);
   fail = tune_continuation_parameters();
   printf("\n");
   fail = determine_output_during_continuation();
   fail = solve_by_multprec_homotopy_continuation(deci);

   return fail;
}

int run_standard_Laurent_continuation ( void )
{
   int fail, nbtasks;

   printf("\nCalling the standard path trackers in PHCpack...\n");
   fail = read_standard_target_Laurent_system();
   fail = read_standard_start_Laurent_system();
   fail = define_output_file();
   fail = write_standard_target_Laurent_system();
   fail = write_standard_start_Laurent_system();
   fail = write_start_solutions();
   printf("\n");
   fail = tune_continuation_parameters();
   printf("\n");
   fail = determine_output_during_continuation();
   printf("\nGive the number of tasks (0 for no multitasking ) : ");
   scanf("%d", &nbtasks);
   printf("\nSee the output file for results ...\n");
   fail = solve_by_standard_Laurent_homotopy_continuation(nbtasks);
   fail = write_target_solutions();

   return fail;
}

int run_dobldobl_Laurent_continuation ( void )
{
   int fail, nbtasks;

   printf("\nCalling the double double path trackers in PHCpack...\n");
   fail = read_dobldobl_target_Laurent_system();
   fail = read_dobldobl_start_Laurent_system();
   fail = define_output_file();
   fail = write_dobldobl_target_Laurent_system();
   fail = write_dobldobl_start_Laurent_system();
   fail = write_dobldobl_start_solutions();
   printf("\n");
   fail = autotune_continuation_parameters(0,14); // (0, 24) is too severe
   fail = tune_continuation_parameters();
   printf("\n");
   fail = determine_output_during_continuation();
   printf("\nGive the number of tasks (0 for no multitasking ) : ");
   scanf("%d", &nbtasks);
   printf("\nSee the output file for results ...\n");
   fail = solve_by_dobldobl_Laurent_homotopy_continuation(nbtasks);
   fail = write_dobldobl_target_solutions();

   return fail;
}

int run_quaddobl_Laurent_continuation ( void )
{
   int fail, nbtasks;

   printf("\nCalling the quad double path trackers in PHCpack...\n");
   fail = read_quaddobl_target_Laurent_system();
   fail = read_quaddobl_start_Laurent_system();
   fail = define_output_file();
   fail = write_quaddobl_target_Laurent_system();
   fail = write_quaddobl_start_Laurent_system();
   fail = write_quaddobl_start_solutions();
   printf("\n");
   fail = autotune_continuation_parameters(0,14); // (0, 48) is too severe
   fail = tune_continuation_parameters();
   printf("\n");
   fail = determine_output_during_continuation();
   printf("\nGive the number of tasks (0 for no multitasking ) : ");
   scanf("%d", &nbtasks);
   printf("\nSee the output file for results ...\n");
   fail = solve_by_quaddobl_Laurent_homotopy_continuation(nbtasks);
   fail = write_quaddobl_target_solutions();

   return fail;
}
