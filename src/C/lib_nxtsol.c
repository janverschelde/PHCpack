/* reads a target and start system and then solves the target system,
   using the start system in an artificial-parameter homotopy */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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

int prompt_for_gamma ( double *regamma, double *imgamma );
/*
 * DESCRIPTION :
 *   Prompts the user for the real and imaginary part of a complex number,
 *   returned in the parameters regamma and imgamma. */

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

int call_initialize_multprec_homotopy ( int *index );
/*
 * DESCRIPTION :
 *   Prepares the containers to initialize the homotopy to track a path
 *   in multiprecision.  Returns in index the number of the
 *   solution in the container selected as start solution. */

int call_initialize_varbprec_solution ( char *name, int *index );
/*
 * DESCRIPTION :
 *   Given the name of the file with start system and start solutions,
 *   the user is prompted for an index and then the solution is read
 *   from file and used to initialize the variable precision homotopy. */

int call_initialize_varbprec_homotopy ( int *index );
/*
 * DESCRIPTION :
 *   Prompts the user for data to initialize the variable precision
 *   path tracker with a homotopy and an initial solution.
 *   On return is the index of the start solution given on file
 *   with the start system that was selected as start solution. */

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

int write_multprec_solution ( int index );
/*
 * DESCRIPTION :
 *   Writes the solution in the multiprecision solutions container
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

int call_multprec_path_tracker ( int index );
/*
 * DESCRIPTION :
 *   Calls the path tracker in multiprecision,
 *   starting a solution with the numbers in the index. */

int call_varbprec_path_tracker ( void );
/*
 * DESCRIPTION :
 *   Calls the path tracker in variable precision,
 *   with the solution the path tracker was initialized with. */ 

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
   else if(level == 2)
   {
      fail = call_initialize_quaddobl_homotopy(&nbsol);
      fail = call_quaddobl_path_tracker(nbsol);
      fail = clear_quaddobl_tracker();
   }
   else if(level == 3)
   {
      fail = call_initialize_multprec_homotopy(&nbsol);
      fail = call_multprec_path_tracker(nbsol);
      fail = clear_multprec_tracker();
   }
   else
   {
      fail = call_initialize_varbprec_homotopy(&nbsol);
      fail = call_varbprec_path_tracker();
      fail = clear_varbprec_tracker();
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
   printf("  2. run in quad double precision arithmetic;\n");
   printf("  3. run in multiprecision arithmetic;\n");
   printf("  4. run in variable precision arithmetic.\n");
   printf("Type 0, 1, 2, 3, or 4 to select precision : ");
   scanf("%d",&answer);
   scanf("%c",&nlc);     /* skip new line symbol */

   return answer;
}

int prompt_for_gamma ( double *regamma, double *imgamma )
{
   printf("Give the real part : "); scanf("%lf", regamma);
   printf("Give the imaginary part : "); scanf("%lf", imgamma);

   return 0;
}

int call_initialize_standard_homotopy ( int *index )
{
   int fail,len,fixed;
   double regamma = 0.0;
   double imgamma = 0.0;

   fail = read_target_system_without_solutions();
   fail = read_standard_start_system();
   fail = copy_start_solutions_to_container();
   fail = solcon_number_of_standard_solutions(&len);
   printf("Number of start solutions : %d\n",len);
   printf("-> give index of solution : "); scanf("%d",index);
   printf("Fixed gamma constant ? (1 = yes/0 = no) "); scanf("%d",&fixed);
   if(fixed != 1) fail = prompt_for_gamma(&regamma,&imgamma);
   fail = initialize_standard_homotopy(fixed,regamma,imgamma);
   fail = initialize_standard_solution(*index);

   return fail;
}

int call_initialize_dobldobl_homotopy ( int *index )
{
   int fail,len,fixed;
   double regamma = 0.0;
   double imgamma = 0.0;

   fail = read_dobldobl_target_system(); /* no _without_solutions ! */
   fail = read_dobldobl_start_system();
   fail = copy_dobldobl_start_solutions_to_container();
   fail = solcon_number_of_dobldobl_solutions(&len);
   printf("Number of start solutions : %d\n",len);
   printf("-> give index of solution : "); scanf("%d",index);
   printf("Fixed gamma constant ? (1 = yes/0 = no) "); scanf("%d",&fixed);
   if(fixed != 1) fail = prompt_for_gamma(&regamma,&imgamma);
   fail = initialize_dobldobl_homotopy(fixed,regamma,imgamma);
   fail = initialize_dobldobl_solution(*index);

   return fail;
}

int call_initialize_quaddobl_homotopy ( int *index )
{
   int fail,len,fixed;
   double regamma = 0.0;
   double imgamma = 0.0;

   fail = read_quaddobl_target_system(); /* no _without_solutions ! */
   fail = read_quaddobl_start_system();
   fail = copy_quaddobl_start_solutions_to_container();
   fail = solcon_number_of_quaddobl_solutions(&len);
   printf("Number of start solutions : %d\n",len);
   printf("-> give index of solution : "); scanf("%d",index);
   printf("Fixed gamma constant ? (1 = yes/0 = no) "); scanf("%d",&fixed);
   if(fixed != 1) fail = prompt_for_gamma(&regamma,&imgamma);
   fail = initialize_quaddobl_homotopy(fixed,regamma,imgamma);
   fail = initialize_quaddobl_solution(*index);

   return fail;
}

int call_initialize_multprec_homotopy ( int *index )
{
   int fail,len,deci,fixed;
   char nlc;

   printf("\ngive the number of decimal places in the working precision : ");
   scanf("%d",&deci);
   scanf("%c",&nlc);     /* skip new line symbol */

   fail = read_multprec_target_system(deci); /* no _without_solutions ! */
   fail = read_multprec_start_system(deci);
   fail = copy_multprec_start_solutions_to_container();
   fail = solcon_number_of_multprec_solutions(&len);
   printf("Number of start solutions : %d\n",len);
   printf("-> give index of solution : "); scanf("%d",index);
   printf("Fixed gamma constant ? (1 = yes/0 = no) "); scanf("%d",&fixed);
   fail = initialize_multprec_homotopy(fixed,deci);
   fail = initialize_multprec_solution(*index);

   return fail;
}

int call_initialize_varbprec_solution ( char *name, int *index )
{
   FILE *fp;
   char *sol;
   int fail,len,nc,nv;

   fp = fopen(name,"r");
   if(fp == NULL)
   {
      printf("File with name %s could not be opened for reading!\n",name);
      *index = -1;
      fail = -1;
   }
   else
   {
      printf("Give the index of the start solution : ");
      scanf("%d",index);
      sol = read_solution_banner_and_string(fp,*index,&len,&nv);
      printf("Solution %d :\n%s\n",*index,sol);
      nc = strlen(sol);
      fail = initialize_varbprec_solution(nv,nc,sol);
   }
   fclose(fp);
   return fail;
}

int call_initialize_varbprec_homotopy ( int *index )
{
   int fail,nc,lentar,lensta,nq,nv,fix,i,idx;
   char name[80];
   char *target,*start;

   printf("\nGive the name of a file to read the target system : ");
   scanf("%s",name);
   nc = strlen(name);
   
   target = read_polynomials_from_file(nc,name,&lentar,&nq,&nv,&fail);
   if(fail != 0)
      printf("Some failure occurred.\n");
   else
   {
      printf("Read %d polynomials in %d variables : \n%s\n",nq,nv,target);

      printf("\nGive the name of a file to read the start system : ");
      scanf("%s",name);
      nc = strlen(name);
   
      start = read_polynomials_from_file(nc,name,&lensta,&nq,&nv,&fail);
      if(fail != 0)
         printf("Some failure occurred.\n");
      else
      {
         printf("Read %d polynomials in %d variables : \n%s\n",nq,nv,start);
         printf("Fixed gamma constant ? (1 = yes/0 = no) ");
         scanf("%d",&fix);
         fail = initialize_varbprec_homotopy(fix,lentar,target,lensta,start);
         if(fail == 0) call_initialize_varbprec_solution(name,&idx);
      }
   }

   return fail;
}

int write_standard_solution ( int index )
{
   int fail,nb;

   fail = solcon_length_standard_solution_string(index,&nb);
   {
      char solution[nb];

      fail = solcon_write_standard_solution_string(index,nb,solution);
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

int write_multprec_solution ( int index )
{
   int fail,nb;

   fail = solcon_length_multprec_solution_string(index,&nb);
   {
      char solution[nb];

      fail = solcon_write_multprec_solution_string(index,nb,solution);
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

int call_multprec_path_tracker ( int index )
{
   int fail;
   char answer;

   fail = write_multprec_solution(index);
   do
   {
      fail = next_multprec_solution(index);
      fail = write_multprec_solution(index);
      printf("Continue to next step ? (y/n) ");
      scanf("%c",&answer); /* get trailing new line...*/
      scanf("%c",&answer);
   }
   while(answer == 'y');

   return fail;
}

int call_varbprec_path_tracker ( void )
{
   int fail,len;
   int want = 8;
   int maxprc = 256;
   int maxitr = 3;
   int vrb = 1;
   char *sol,answer;

   do
   {
      sol = next_varbprec_solution(want,maxprc,maxitr,vrb,&len,&fail);
      if(fail == 0)
      {
         printf("The next solution :\n%s\n",sol);
         free(sol); /* only free when no failure! */
      }
      printf("Continue to next step ? (y/n) ");
      scanf("%c",&answer); /* get trailing new line...*/
      scanf("%c",&answer);
   }
   while(answer == 'y');

   return fail;
}
