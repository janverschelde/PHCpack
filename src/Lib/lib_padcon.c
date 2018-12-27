/* simple test on running a Pade continuation */

#include <stdio.h>
#include <string.h>
#include "solcon.h"
#include "phcpack.h"
#include "padcon.h"

void write_homotopy_continuation_parameters ( void );
/*
 * DESCRIPTION :
 *   Writes the current values of the homotopy continuation parameters. */

void tune_homotopy_continuation_parameters ( void );
/*
 * DESCRIPTION :
 *   Interactive loop to tune the homotopy continuation parameters. */

void prompt_for_output_file ( int* nbc, char* name, int *verbose );
/*
 * DESCRIPTION :
 *   Prompts the user for the verbosity, returned in verbose.
 *   Prompts the user for an output file, returned in name,
 *   with in nbc the number of characters in name. */

void standard_track ( int nbc, char* name, int verbose );
/*
 * DESCRIPTION :
 *   Tracks in standard double precision.  On input is the name of the
 *   output file, with its number of characters in nbc.
 *   If nbc = 0, then no output is written to file,
 *   otherwise, the output is written to the name with file name.
 *   If verbose > 0, then more output is written. */

void dobldobl_track ( int nbc, char* name, int verbose );
/*
 * DESCRIPTION :
 *   Tracks in double double precision.  On input is the name of the
 *   output file, with its number of characters in nbc.
 *   If nbc = 0, then no output is written to file,
 *   otherwise, the output is written to the name with file name.
 *   If verbose > 0, then more output is written. */

void quaddobl_track ( int nbc, char* name, int verbose );
/*
 * DESCRIPTION :
 *   Tracks in quad double precision.  On input is the name of the
 *   output file, with its number of characters in nbc.
 *   If nbc = 0, then no output is written to file,
 *   otherwise, the output is written to the name with file name.
 *   If verbose > 0, then more output is written. */

void standard_next_step ( void );
/*
 *  DESCRIPTION :
 *    Runs the series-Pade tracker step by step in double precision. */

void dobldobl_next_step ( void );
/*
 *  DESCRIPTION :
 *    Runs the series-Pade tracker step by step in double double precision. */

void quaddobl_next_step ( void );
/*
 *  DESCRIPTION :
 *    Runs the series-Pade tracker step by step in quad double precision. */

int main ( void )
{
   adainit();

   char nlsb,ans;
   int fail,precision,nbchar,verbose;
   char filename[80];

   printf("\nMENU for the precision :\n");
   printf("  0. double precision\n");
   printf("  1. double double precision\n");
   printf("  2. quad double precision\n");
   printf("Type 0, 1, or 2 to select the precision : ");
   scanf("%d", &precision);

   printf("\nTuning the homotopy continuation parameters ...\n");

   fail = padcon_set_default_parameters();
   tune_homotopy_continuation_parameters();

   scanf("%c", &nlsb); // swallow new line symbol

   printf("\nStep-by-step run ? (y/n) ");
   ans = getchar();
   scanf("%c", &nlsb); // swallow new line symbol
   if(ans == 'y')
   {
      if(precision == 0) standard_next_step();
      if(precision == 1) dobldobl_next_step();
      if(precision == 2) quaddobl_next_step();
   }
   else
   {
      prompt_for_output_file(&nbchar,filename,&verbose);

      if(nbchar > 0)
         printf("\nThe name of the output file is %s.\n", filename);

      if(precision == 0) standard_track(nbchar,filename,verbose);
      if(precision == 1) dobldobl_track(nbchar,filename,verbose);
      if(precision == 2) quaddobl_track(nbchar,filename,verbose);
   }
   adafinal();

   return 0;
}

void write_homotopy_continuation_parameters ( void )
{
   double fltval[2];
   int fail,intval;

   printf("Values of the HOMOTOPY CONTINUATION PARAMETERS :\n");
   fail = padcon_get_homotopy_continuation_parameter(1,fltval);
   printf(" 1. gamma : ");
   printf("%+.14E", fltval[0]); printf("  %+.14E\n", fltval[1]);
   printf(" 2. degree of numerator of Pade approximant    : ");
   fail = padcon_get_homotopy_continuation_parameter(2,&fltval[0]);
   printf("%d\n", (int) fltval[0]);
   printf(" 3. degree of denominator of Pade approximant  : ");
   fail = padcon_get_homotopy_continuation_parameter(3,&fltval[0]);
   printf("%d\n", (int) fltval[0]);
   printf(" 4. maximum step size                          : ");
   fail = padcon_get_homotopy_continuation_parameter(4,&fltval[0]);
   printf("%.2E\n", fltval[0]);
   printf(" 5. minimum step size                          : ");
   fail = padcon_get_homotopy_continuation_parameter(5,&fltval[0]);
   printf("%.2E\n", fltval[0]);
   printf(" 6. multiplication factor of the series step   : ");
   fail = padcon_get_homotopy_continuation_parameter(6,&fltval[0]);
   printf("%.2E\n", fltval[0]);
   printf(" 7. multiplication factor of the pole radius   : ");
   fail = padcon_get_homotopy_continuation_parameter(7,&fltval[0]);
   printf("%.2E\n", fltval[0]);
   printf(" 8. tolerance on the residual of the predictor : ");
   fail = padcon_get_homotopy_continuation_parameter(8,&fltval[0]);
   printf("%.2E\n", fltval[0]);
   printf(" 9. tolerance on the residual of the corrector : ");
   fail = padcon_get_homotopy_continuation_parameter(9,&fltval[0]);
   printf("%.2E\n", fltval[0]);
   printf("10. tolerance on zero series coefficients      : ");
   fail = padcon_get_homotopy_continuation_parameter(10,&fltval[0]);
   printf("%.2E\n", fltval[0]);
   printf("11. maximum number of corrector steps          : ");
   fail = padcon_get_homotopy_continuation_parameter(11,&fltval[0]);
   printf("%d\n", (int) fltval[0]);
   printf("12. maximum steps on a path                    : ");
   fail = padcon_get_homotopy_continuation_parameter(12,&fltval[0]);
   printf("%d\n", (int) fltval[0]);
}

void tune_homotopy_continuation_parameters ( void )
{
   int choice,intpar,fail;
   double fltpar;
   double gamma[2];

   do
   {
      write_homotopy_continuation_parameters();
      printf("Type a number to change a value, or 0 to exit : ");
      scanf("%d", &choice);
      if(choice == 1)
      {
         printf("-> give the real part of the new gamma : ");
         scanf("%lf", &gamma[0]);
         printf("-> give the imaginary part of the new gamma : ");
         scanf("%lf", &gamma[1]);
         fail = padcon_set_homotopy_continuation_parameter(1,gamma);
      }
      else if(choice == 2)
      {
         printf("-> give a new numerator degree for the Pade approximant : ");
         scanf("%d", &intpar); fltpar = (double) intpar;
      }
      else if(choice == 3)
      {
         printf("-> give a new denominator degree for the Pade approximant : ");
         scanf("%d", &intpar); fltpar = (double) intpar;
      }
      else if(choice == 4)
      {
         printf("-> give a new value for the maximum step size : ");
         scanf("%lf", &fltpar);
      }
      else if(choice == 5)
      {
         printf("-> give a new value for the minimum step size  : ");
         scanf("%lf", &fltpar);
      }
      else if(choice == 6)
      {
         printf("-> give a new multiplication factor for the series step : ");
         scanf("%lf", &fltpar);
      }
      else if(choice == 7)
      {
         printf("-> give a new multiplication factor for the pole radius : ");
         scanf("%lf", &fltpar);
      }
      else if(choice == 8)
      {
         printf("-> give a new tolerance on the predictor residual : ");
         scanf("%lf", &fltpar);
      }
      else if(choice == 9)
      {
         printf("-> give a new tolerance on the corrector residual : ");
         scanf("%lf", &fltpar);
      }
      else if(choice == 10)
      {
         printf("-> give a new tolerance on a zero series coefficient : ");
         scanf("%lf", &fltpar);
      }
      else if(choice == 11)
      {
         printf("-> give a new maximum number of corrector steps : ");
         scanf("%d", &intpar); fltpar = (double) intpar;
      }
      else if(choice == 12)
      {
         printf("-> give a new maximum number of steps on a path : ");
         scanf("%d", &intpar); fltpar = (double) intpar;
      }
      if((choice > 1) && (choice < 13))
      {
         // printf("setting parameter %d to %.3e ...\n", choice, fltpar);
         fail = padcon_set_homotopy_continuation_parameter(choice,&fltpar);
      }
   }
   while (choice != 0);
}

void prompt_for_output_file ( int* nbc, char* name, int *verbose )
{
   char ans,nlsb;

   printf("\nVerbose?  Extra output desired ? (y/n) ? ");
   scanf("%c",&ans);
   scanf("%c",&nlsb); /* skip newline symbol */
   *verbose = (ans == 'y');

   printf("\nOutput to separate file ? (y/n) ? ");
   scanf("%c",&ans);
   scanf("%c",&nlsb); /* skip newline symbol */

   if(ans != 'y')
      *nbc = 0;
   else
   {
      printf("-> give the name of the output file : ");
      scanf("%s",name);
      scanf("%c",&nlsb); /* skip newline symbol */
      *nbc = strlen(name);
   }
}

void standard_track ( int nbc, char* name, int verbose )
{
   int fail,length;

   fail = read_standard_target_system();
   fail = read_standard_start_system();
   fail = copy_start_solutions_to_container();
   fail = solcon_number_of_standard_solutions(&length);
   printf("Read %d start solutions.\n", length);

   if(nbc > 0) printf("\nSee the output file %s ...\n", name);

   fail = padcon_standard_track(nbc,name,verbose);

   if(nbc == 0) fail = solcon_write_standard_solutions();
}

void dobldobl_track ( int nbc, char* name, int verbose )
{
   int fail,length;

   fail = read_dobldobl_target_system();
   fail = read_dobldobl_start_system();
   fail = copy_dobldobl_start_solutions_to_container();
   fail = solcon_number_of_dobldobl_solutions(&length);
   printf("Read %d start solutions.\n", length);

   if(nbc > 0) printf("\nSee the output file %s ...\n", name);

   fail = padcon_dobldobl_track(nbc,name,verbose);

   if(nbc == 0) fail = solcon_write_dobldobl_solutions();
}

void quaddobl_track ( int nbc, char* name, int verbose )
{
   int fail,length;

   fail = read_quaddobl_target_system();
   fail = read_quaddobl_start_system();
   fail = copy_quaddobl_start_solutions_to_container();
   fail = solcon_number_of_quaddobl_solutions(&length);
   printf("Read %d start solutions.\n", length);

   if(nbc > 0) printf("\nSee the output file %s ...\n", name);

   fail = padcon_quaddobl_track(nbc,name,verbose);

   if(nbc == 0) fail = solcon_write_quaddobl_solutions();
}

void standard_next_step ( void )
{
   int fail,length,index,failed,strlensol;
   char contstep,nlsb;
   double frp,re_pole,im_pole;

   fail = read_standard_target_system();
   fail = read_standard_start_system();
   fail = copy_start_solutions_to_container();
   fail = solcon_number_of_standard_solutions(&length);
   printf("Read %d start solutions.\n", length);

   fail = padcon_standard_initialize_homotopy(1);
   for(index=1; index <= length; index++)
   {
      printf("\nTracking path %d ...\n", index);
      fail = padcon_initialize_standard_solution(index,1);
      do
      {
         fail = padcon_standard_predict_correct(&failed,1);
         if(failed != 0)
            contstep = 'n';
         else
         {
            fail = padcon_get_standard_forward_pole_radius(&frp);
            fail = padcon_get_standard_closest_pole(&re_pole,&im_pole);
            printf("Smallest forward pole radius : %.3e\n", frp);
            if(re_pole >= 0.0)
               printf("Closest pole : %.14e  %.14e\n", re_pole, im_pole);
            printf("continue ? (y/n) ");
            contstep = getchar();
            scanf("%c",&nlsb);
         }
      }
      while(contstep == 'y');

      fail = padcon_get_standard_solution(index,1);
      fail = solcon_length_standard_solution_string(index,&strlensol);
      {
         char idxsol[strlensol+1];

         fail = solcon_write_standard_solution_string
                  (index,strlensol,idxsol);
         printf("Solution %d :\n%s\n", index, idxsol);
      }
      if(index == length) break;
      printf("Continue with the next solution ? (y/n) ");
      contstep = getchar();
      scanf("%c",&nlsb);
      if(contstep != 'y') break;
   }
}

void dobldobl_next_step ( void )
{
   int fail,length,index,failed,strlensol;
   char contstep,nlsb;
   double frp,re_pole,im_pole;

   fail = read_dobldobl_target_system();
   fail = read_dobldobl_start_system();
   fail = copy_dobldobl_start_solutions_to_container();
   fail = solcon_number_of_dobldobl_solutions(&length);
   printf("Read %d start solutions.\n", length);

   fail = padcon_dobldobl_initialize_homotopy(1);
   for(index = 1; index <= length; index++)
   {
      printf("\nTracking path %d ...\n", index);
      fail = padcon_initialize_dobldobl_solution(index,1);
      do
      {
         fail = padcon_dobldobl_predict_correct(&failed,1);
         if(failed != 0)
            contstep = 'n';
         else
         {
            fail = padcon_get_dobldobl_forward_pole_radius(&frp);
            fail = padcon_get_dobldobl_closest_pole(&re_pole,&im_pole);
            printf("Smallest forward pole radius : %.3e\n", frp);
            if(re_pole >= 0.0)
               printf("Closest pole : %.14e  %.14e\n", re_pole, im_pole);
            printf("continue ? (y/n) ");
            contstep = getchar();
            scanf("%c",&nlsb);
         }
      }
      while(contstep == 'y');

      fail = padcon_get_dobldobl_solution(index,1);
      fail = solcon_length_dobldobl_solution_string(index,&strlensol);
      {
         char idxsol[strlensol+1];

         fail = solcon_write_dobldobl_solution_string
                  (index,strlensol,idxsol);
         printf("Solution %d :\n%s\n", index, idxsol);
      }
      if(index == length) break;
      printf("Continue with the next solution ? (y/n) ");
      contstep = getchar();
      scanf("%c",&nlsb);
      if(contstep != 'y') break;
   }
}

void quaddobl_next_step ( void )
{
   int fail,length,index,failed,strlensol;
   char contstep,nlsb;
   double frp,re_pole,im_pole;

   fail = read_quaddobl_target_system();
   fail = read_quaddobl_start_system();
   fail = copy_quaddobl_start_solutions_to_container();
   fail = solcon_number_of_quaddobl_solutions(&length);
   printf("Read %d start solutions.\n", length);

   fail = padcon_quaddobl_initialize_homotopy(1);
   for(index = 1; index <= length; index++)
   {
      printf("\nTracking path %d ...\n", index);
      fail = padcon_initialize_quaddobl_solution(index,1);
      do
      {
         fail = padcon_quaddobl_predict_correct(&failed,1);
         if(failed != 0)
            contstep = 'n';
         else
         {
            fail = padcon_get_quaddobl_forward_pole_radius(&frp);
            fail = padcon_get_quaddobl_closest_pole(&re_pole,&im_pole);
            printf("Smallest forward pole radius : %.3e\n", frp);
            if(re_pole >= 0.0)
               printf("Closest pole : %.14e  %.14e\n", re_pole, im_pole);
            printf("continue ? (y/n) ");
            contstep = getchar();
            scanf("%c",&nlsb);
         }
      }
      while(contstep == 'y');

      fail = padcon_get_quaddobl_solution(index,1);
      fail = solcon_length_quaddobl_solution_string(index,&strlensol);
      {
         char idxsol[strlensol+1];

         fail = solcon_write_quaddobl_solution_string
                  (index,strlensol,idxsol);
         printf("Solution %d :\n%s\n", index, idxsol);
      }
      if(index == length) break;
      printf("Continue with the next solution ? (y/n) ");
      contstep = getchar();
      scanf("%c",&nlsb);
      if(contstep != 'y') break;
   }
}
