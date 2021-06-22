/* test on the functions in gpupath_qd */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "phcpack.h"
#include "syscon.h"
#include "solcon.h"
#include "gpupath_qd.h"

void test_newton ( void );
/* 
 * DESCRIPTION :
 *   The user is prompted for a start system
 *   and Newton's method is called. */

void test_onepath ( void );
/* 
 * DESCRIPTION :
 *   The user is prompted for a target and start system
 *   and one solution path is tracked. */

void test_manypaths ( void );
/* 
 * DESCRIPTION :
 *   The user is prompted for a target and start system
 *   and many solution paths are tracked. */

int main ( int argc, char *argv[] )
{
   int choice;
   char ch;

   adainit();

   printf("\nMENU to apply algorithmic differentiation on the GPU: \n");
   printf("  0. in Newton's method to correct one solution;\n");
   printf("  1. to track one solution path;\n");
   printf("  2. to track many solution paths.\n");
   printf("Type 0, 1, or 2 to select : ");
   scanf("%d",&choice);
   scanf("%c",&ch); /* skip newline symbol */

   if(choice == 0)
      test_newton();
   else if(choice == 1)
      test_onepath();
   else if(choice == 2)
      test_manypaths();
   else
      printf("invalid selection, please try again\n");

   adafinal();

   return 0;
}

void test_newton ( void )
{
   int fail,dim,len;

   printf("\nRunning Newton's method ...\n");
   fail = read_quaddobl_start_system();
   fail = copy_quaddobl_start_system_to_container();
   fail = copy_quaddobl_start_solutions_to_container();
   fail = syscon_number_of_quaddobl_polynomials(&dim);
   printf("The system container has %d polynomials.\n",dim);
   fail = solcon_number_of_quaddobl_solutions(&len);
   printf("The solution container has size %d.\n",len);
   fail = gpu_newton_qd(2,1);
   printf("The solutions after Newton's method :\n");
   fail = solcon_write_quaddobl_solutions();
}

void test_onepath ( void )
{
   int fail,len;

   printf("\nTracking one solution path ...\n");
   fail = read_quaddobl_target_system();
   fail = read_quaddobl_start_system();
   fail = copy_quaddobl_start_solutions_to_container();
   fail = solcon_number_of_quaddobl_solutions(&len);
   printf("The solution container has size %d.\n",len);
   fail = gpu_onepath_qd(2,1,1.0,0.0);
   printf("The solutions after Newton's method :\n");
   fail = solcon_write_quaddobl_solutions();
}

void test_manypaths ( void )
{
   int fail,len;

   printf("\nTracking many solution paths ...\n");
   fail = read_quaddobl_target_system();
   fail = read_quaddobl_start_system();
   fail = copy_quaddobl_start_solutions_to_container();
   fail = solcon_number_of_quaddobl_solutions(&len);
   printf("The solution container has size %d.\n",len);
   fail = gpu_manypaths_qd(2,1,1.0,0.0);
   printf("The solutions after Newton's method :\n");
   fail = solcon_write_quaddobl_solutions();
}
