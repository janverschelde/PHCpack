/* simple test on running a Pade continuation */

#include <stdio.h>
#include <string.h>
#include "syscon.h"
#include "solcon.h"
#include "phcpack.h"
#include "padcon.h"

int prompt_for_artificial ( void );
/*
 * DESCRIPTION :
 *   Ask the user if the homotopy is artificial or not.
 *   Returns 1 for an artificial parameter homotopy, 0 otherwise. */

int prompt_for_parameter ( void );
/*
 * DESCRIPTION :
 *   Displays the variable names and prompts the user for the
 *   index of the parameter, which is returned. */

void prompt_for_output_file ( int* nbc, char* name, int *verbose );
/*
 * DESCRIPTION :
 *   Prompts the user for the verbosity, returned in verbose.
 *   Prompts the user for an output file, returned in name,
 *   with in nbc the number of characters in name. */

void standard_track ( int nbc, char *name, int verbose );
/*
 * DESCRIPTION :
 *   Tracks in standard double precision.  On input is the name of the
 *   output file, with its number of characters in nbc.
 *   If nbc = 0, then no output is written to file,
 *   otherwise, the output is written to the name with file name.
 *   If verbose > 0, then more output is written. */

void dobldobl_track ( int nbc, char *name, int verbose );
/*
 * DESCRIPTION :
 *   Tracks in double double precision.  On input is the name of the
 *   output file, with its number of characters in nbc.
 *   If nbc = 0, then no output is written to file,
 *   otherwise, the output is written to the name with file name.
 *   If verbose > 0, then more output is written. */

void quaddobl_track ( int nbc, char *name, int verbose );
/*
 * DESCRIPTION :
 *   Tracks in quad double precision.  On input is the name of the
 *   output file, with its number of characters in nbc.
 *   If nbc = 0, then no output is written to file,
 *   otherwise, the output is written to the name with file name.
 *   If verbose > 0, then more output is written. */

int maximal_series_degree ( void );
/*
 * DESCRIPTION :
 *   After the homotopy continuation parameters are set,
 *   the maximal degree of the series is determined as well.
 *   This function returns the maximal series degree. */

void show_standard_series_coefficients ( int dim );
/*
 * DESCRIPTION :
 *   Shows all coefficients of the power series computed by
 *   the predictor in standard double precision for a system with
 *   as many variables as the value of dim. */

void show_dobldobl_series_coefficients ( int dim );
/*
 * DESCRIPTION :
 *   Shows all coefficients of the power series computed by
 *   the predictor in double double precision for a system with
 *   as many variables as the value of dim. */

void show_quaddobl_series_coefficients ( int dim );
/*
 * DESCRIPTION :
 *   Shows all coefficients of the power series computed by
 *   the predictor in quad double precision for a system with
 *   as many variables as the value of dim. */

void show_standard_pade_coefficients ( int dim );
/*
 * DESCRIPTION :
 *   Shows all coefficients of the Pade approximants computed by
 *   the predictor in standard double precision for a system with
 *   as many variables as the value of dim. */

void show_dobldobl_pade_coefficients ( int dim );
/*
 * DESCRIPTION :
 *   Shows all coefficients of the Pade approximants computed by
 *   the predictor in double double precision for a system with
 *   as many variables as the value of dim. */

void show_quaddobl_pade_coefficients ( int dim );
/*
 * DESCRIPTION :
 *   Shows all coefficients of the Pade approximants computed by
 *   the predictor in quad double precision for a system with
 *   as many variables as the value of dim. */

void show_standard_poles ( int dim );
/*
 * DESCRIPTION :
 *   Shows all poles computed by the predictor in double precision.
 *   The value of dim equals the number of variables. */

void show_dobldobl_poles ( int dim );
/*
 * DESCRIPTION :
 *   Shows all poles computed by the predictor in double double precision.
 *   The value of dim equals the number of variables. */

void show_quaddobl_poles ( int dim );
/*
 * DESCRIPTION :
 *   Shows all poles computed by the predictor in quad double precision.
 *   The value of dim equals the number of variables. */

void show_standard_solution ( int index );
/*
 * DESCRIPTION :
 *   Shows the current solution in standard double precision,
 *   at position index in the solutions container. */

void show_dobldobl_solution ( int index );
/*
 * DESCRIPTION :
 *   Shows the current solution in double double precision,
 *   at position index in the solutions container. */

void show_quaddobl_solution ( int index );
/*
 * DESCRIPTION :
 *   Shows the current solution in quad double precision,
 *   at position index in the solutions container. */

void standard_next_loop ( int index );
/*
 *  DESCRIPTION :
 *    Runs the series-Pade tracker step by step in double precision
 *    on the solution with the given index in the double precision
 *    solutions container. */

void dobldobl_next_loop ( int index );
/*
 *  DESCRIPTION :
 *    Runs the series-Pade tracker step by step in double double precision
 *    on the solution with the given index in the double double precision
 *    solutions container. */

void quaddobl_next_loop ( int index );
/*
 *  DESCRIPTION :
 *    Runs the series-Pade tracker step by step in quad double precision
 *    on the solution with the given index in the quad double precision
 *    solutions container. */

void standard_next_step ( int homo );
/*
 *  DESCRIPTION :
 *    Runs the series-Pade tracker step by step in double precision.
 *    If homo = 0, then tracking happens in affine coordinates,
 *    otherwise projective coordinate transformations are applied. */

void dobldobl_next_step ( int homo );
/*
 *  DESCRIPTION :
 *    Runs the series-Pade tracker step by step in double double precision.
 *    If homo = 0, then tracking happens in affine coordinates,
 *    otherwise projective coordinate transformations are applied. */

void quaddobl_next_step ( int homo );
/*
 *  DESCRIPTION :
 *    Runs the series-Pade tracker step by step in quad double precision.
 *    If homo = 0, then tracking happens in affine coordinates,
 *    otherwise projective coordinate transformations are applied. */

int main ( void )
{
   adainit();

   char nlsb,ans;
   int fail,precision,nbchar,homo,verbose;
   char filename[80];

   printf("\nMENU for the precision :\n");
   printf("  0. double precision\n");
   printf("  1. double double precision\n");
   printf("  2. quad double precision\n");
   printf("Type 0, 1, or 2 to select the precision : ");
   scanf("%d", &precision);

   printf("\nTuning the homotopy continuation parameters ...\n");

   fail = padcon_set_default_parameters();
   padcon_tune_homotopy_continuation_parameters();

   scanf("%c", &nlsb); // swallow new line symbol

   printf("\nStep-by-step run ? (y/n) ");
   ans = getchar();
   scanf("%c", &nlsb); // swallow new line symbol

   if(ans == 'y')
   {
      homo = padcon_prompt_for_homogenization();

      if(precision == 0) standard_next_step(homo);
      if(precision == 1) dobldobl_next_step(homo);
      if(precision == 2) quaddobl_next_step(homo);
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

int prompt_for_artificial ( void )
{
   char answer,nlsb;

   printf("\nEither a homotopy has a parameter among its variables,\n");
   printf("or the parameter is artificial and the homotopy connects\n");
   printf("a target system to a start system with known solutions.\n");
   printf("Is the homotopy an artificial parameter homotopy ? (y/n) ");

   answer = getchar();
   scanf("%c", &nlsb); // swallow new line symbol

   if(answer == 'y')
      return 1;
   else
      return 0;
}

int prompt_for_parameter ( void )
{
   int fail,nbr,idx;
   char nlsb;

   fail = syscon_number_of_symbols(&nbr);

   printf("\nThe names of the variables :");
   fail = syscon_write_symbols();

   printf("\nGive the index of the continuation parameter : ");
   scanf("%d",&idx);
   scanf("%c",&nlsb); // swallow the newline symbol

   if(idx < 1) printf("Index is too small!\n");
   if(idx > nbr) printf("Index is too large!\n");

   return idx;
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

void standard_track ( int nbc, char *name, int verbose )
{
   int fail,length,mhom,dim;

   fail = read_standard_target_system();
   fail = read_standard_start_system();
   fail = copy_start_solutions_to_container();
   fail = solcon_number_of_standard_solutions(&length);
   printf("Read %d start solutions.\n", length);

   fail = syscon_number_of_symbols(&dim);
   mhom = padcon_prompt_for_multi_homogenization(dim);

   if(mhom > 1) 
   {
      int idz[dim];
      padcon_define_partition(mhom,dim,idz);
      padcon_standard_multi_projective_transformation(dim,mhom,idz);
      if(nbc > 0) printf("\nSee the output file %s ...\n", name);
      fail = padcon_standard_track(nbc,name,0,verbose,mhom,dim,idz);
      fail = solcon_standard_multi_affinization(dim,mhom,idz);
      fail = copy_target_system_to_container();
      fail = syscon_standard_multi_affinization(mhom);
   }
   else
   {
      if(mhom == 1) padcon_standard_projective_transformation();
      if(nbc > 0) printf("\nSee the output file %s ...\n", name);
      fail = padcon_standard_track(nbc,name,0,verbose,mhom,0,&fail);
      if(mhom == 1)
      {
         fail = solcon_standard_one_affinization();
         fail = copy_target_system_to_container();
         fail = syscon_standard_one_affinization();
      }
   }
   fail = syscon_write_standard_system();
   if(nbc != 0) fail = solcon_write_solution_banner_to_defined_output_file();
   fail = solcon_write_standard_solutions();
}

void dobldobl_track ( int nbc, char *name, int verbose )
{
   int fail,length,mhom,dim;

   fail = read_dobldobl_target_system();
   fail = read_dobldobl_start_system();
   fail = copy_dobldobl_start_solutions_to_container();
   fail = solcon_number_of_dobldobl_solutions(&length);
   printf("Read %d start solutions.\n", length);

   fail = syscon_number_of_symbols(&dim);
   mhom = padcon_prompt_for_multi_homogenization(dim);

   if(mhom > 1) 
   {
      int idz[dim];
      padcon_define_partition(mhom,dim,idz);
      padcon_dobldobl_multi_projective_transformation(dim,mhom,idz);
      if(nbc > 0) printf("\nSee the output file %s ...\n", name);
      fail = padcon_dobldobl_track(nbc,name,0,verbose,mhom,dim,idz);
      fail = solcon_dobldobl_multi_affinization(dim,mhom,idz);
      fail = copy_dobldobl_target_system_to_container();
      fail = syscon_dobldobl_multi_affinization(mhom);
   }
   else
   {
      if(mhom == 1) padcon_dobldobl_projective_transformation();
      if(nbc > 0) printf("\nSee the output file %s ...\n", name);
      fail = padcon_dobldobl_track(nbc,name,0,verbose,mhom,0,&fail);
      if(mhom == 1)
      {
         fail = solcon_dobldobl_one_affinization();
         fail = copy_dobldobl_target_system_to_container();
         fail = syscon_dobldobl_one_affinization();
      }
   }
   fail = syscon_write_dobldobl_system();
   if(nbc != 0) fail = solcon_write_solution_banner_to_defined_output_file();
   fail = solcon_write_dobldobl_solutions();
}

void quaddobl_track ( int nbc, char *name, int verbose )
{
   int fail,length,mhom,dim;

   fail = read_quaddobl_target_system();
   fail = read_quaddobl_start_system();
   fail = copy_quaddobl_start_solutions_to_container();
   fail = solcon_number_of_quaddobl_solutions(&length);
   printf("Read %d start solutions.\n", length);

   fail = syscon_number_of_symbols(&dim);
   mhom = padcon_prompt_for_multi_homogenization(dim);

   if(mhom > 1) 
   {
      int idz[dim];
      padcon_define_partition(mhom,dim,idz);
      padcon_quaddobl_multi_projective_transformation(dim,mhom,idz);
      if(nbc > 0) printf("\nSee the output file %s ...\n", name);
      fail = padcon_quaddobl_track(nbc,name,0,verbose,mhom,dim,idz);
      fail = solcon_quaddobl_multi_affinization(dim,mhom,idz);
      fail = copy_quaddobl_target_system_to_container();
      fail = syscon_quaddobl_multi_affinization(mhom);
   }
   else
   {
      if(mhom == 1) padcon_quaddobl_projective_transformation();
      if(nbc > 0) printf("\nSee the output file %s ...\n", name);
      fail = padcon_quaddobl_track(nbc,name,0,verbose,mhom,0,&fail);
      if(mhom == 1) 
      {
         fail = solcon_quaddobl_one_affinization();
         fail = copy_quaddobl_target_system_to_container();
         fail = syscon_quaddobl_one_affinization();
      }
   }
   fail = syscon_write_quaddobl_system();
   if(nbc != 0) fail = solcon_write_solution_banner_to_defined_output_file();
   fail = solcon_write_quaddobl_solutions();
}

int maximal_series_degree ( void )
{
   int fail,numdeg,dendeg;
   double val;

   fail = padcon_get_homotopy_continuation_parameter(2,&val);
   numdeg = (int) val;
   fail = padcon_get_homotopy_continuation_parameter(3,&val);
   dendeg = (int) val;

   return numdeg + dendeg;
}

void show_standard_series_coefficients ( int dim )
{
   int fail,leadidx,idx;
   const int maxdeg = maximal_series_degree();
   double re,im;
  
   for(leadidx=1; leadidx<=dim; leadidx++)
   {
      printf("Coefficient of series at component %d :\n", leadidx);
      for(idx=0; idx<=maxdeg; idx++)
      {
         fail = padcon_get_standard_series_coefficient(leadidx,idx,0,&re,&im);
         printf(" %d : %.14e  %.14e\n", idx, re, im);
      }
   }
}

void show_dobldobl_series_coefficients ( int dim )
{
   int fail,leadidx,idx;
   const int maxdeg = maximal_series_degree();
   double re,im;
  
   for(leadidx=1; leadidx<=dim; leadidx++)
   {
      printf("Coefficient of series at component %d :\n", leadidx);
      for(idx=0; idx<=maxdeg; idx++)
      {
         fail = padcon_get_dobldobl_series_coefficient(leadidx,idx,0,&re,&im);
         printf(" %d : %.14e  %.14e\n", idx, re, im);
      }
   }
}

void show_quaddobl_series_coefficients ( int dim )
{
   int fail,leadidx,idx;
   const int maxdeg = maximal_series_degree();
   double re,im;
  
   for(leadidx=1; leadidx<=dim; leadidx++)
   {
      printf("Coefficient of series at component %d :\n", leadidx);
      for(idx=0; idx<=maxdeg; idx++)
      {
         fail = padcon_get_quaddobl_series_coefficient(leadidx,idx,0,&re,&im);
         printf(" %d : %.14e  %.14e\n", idx, re, im);
      }
   }
}

void show_standard_pade_coefficients ( int dim )
{
   int fail,lead,idx,numdeg,dendeg;
   double val,re,im;
  
   fail = padcon_get_homotopy_continuation_parameter(2,&val);
   numdeg = (int) val;
   fail = padcon_get_homotopy_continuation_parameter(3,&val);
   dendeg = (int) val;

   for(lead=1; lead<=dim; lead++)
   {
      printf("Numerator coefficient of Pade approximant %d :\n", lead);
      for(idx=0; idx<=numdeg; idx++)
      {
         fail = padcon_get_standard_numerator_coefficient(lead,idx,0,&re,&im);
         printf(" %d : %.14e  %.14e\n", idx, re, im);
      }
   }
   for(lead=1; lead<=dim; lead++)
   {
      printf("Denominator coefficient of Pade approximant %d :\n", lead);
      for(idx=0; idx<=dendeg; idx++)
      {
         fail = padcon_get_standard_denominator_coefficient(lead,idx,0,&re,&im);
         printf(" %d : %.14e  %.14e\n", idx, re, im);
      }
   }
}

void show_dobldobl_pade_coefficients ( int dim )
{
   int fail,lead,idx,numdeg,dendeg;
   double val,re,im;
  
   fail = padcon_get_homotopy_continuation_parameter(2,&val);
   numdeg = (int) val;
   fail = padcon_get_homotopy_continuation_parameter(3,&val);
   dendeg = (int) val;

   for(lead=1; lead<=dim; lead++)
   {
      printf("Numerator coefficient of Pade approximant %d :\n", lead);
      for(idx=0; idx<=numdeg; idx++)
      {
         fail = padcon_get_dobldobl_numerator_coefficient(lead,idx,0,&re,&im);
         printf(" %d : %.14e  %.14e\n", idx, re, im);
      }
   }
   for(lead=1; lead<=dim; lead++)
   {
      printf("Denominator coefficient of Pade approximant %d :\n", lead);
      for(idx=0; idx<=dendeg; idx++)
      {
         fail = padcon_get_dobldobl_denominator_coefficient(lead,idx,0,&re,&im);
         printf(" %d : %.14e  %.14e\n", idx, re, im);
      }
   }
}

void show_quaddobl_pade_coefficients ( int dim )
{
   int fail,lead,idx,numdeg,dendeg;
   double val,re,im;
  
   fail = padcon_get_homotopy_continuation_parameter(2,&val);
   numdeg = (int) val;
   fail = padcon_get_homotopy_continuation_parameter(3,&val);
   dendeg = (int) val;

   for(lead=1; lead<=dim; lead++)
   {
      printf("Numerator coefficient of Pade approximant %d :\n", lead);
      for(idx=0; idx<=numdeg; idx++)
      {
         fail = padcon_get_quaddobl_numerator_coefficient(lead,idx,0,&re,&im);
         printf(" %d : %.14e  %.14e\n", idx, re, im);
      }
   }
   for(lead=1; lead<=dim; lead++)
   {
      printf("Denominator coefficient of Pade approximant %d :\n", lead);
      for(idx=0; idx<=dendeg; idx++)
      {
         fail = padcon_get_quaddobl_denominator_coefficient(lead,idx,0,&re,&im);
         printf(" %d : %.14e  %.14e\n", idx, re, im);
      }
   }
}

void show_standard_poles ( int dim )
{
   int fail,dendeg,leadidx,poleidx;
   double val,re,im;

   fail = padcon_get_homotopy_continuation_parameter(3,&val);
   dendeg = (int) val;

   for(leadidx=1; leadidx<=dim; leadidx++)
   {
      printf("poles for Pade approximant %d :\n", leadidx);
      for(poleidx=1; poleidx<=dendeg; poleidx++)
      {
         fail = padcon_get_standard_pole(leadidx,poleidx,0,&re,&im);
         printf(" %d : %.14e  %.14e\n", poleidx, re, im);
      }
   }
}

void show_dobldobl_poles ( int dim )
{
   int fail,dendeg,leadidx,poleidx;
   double val,re,im;

   fail = padcon_get_homotopy_continuation_parameter(3,&val);
   dendeg = (int) val;

   for(leadidx=1; leadidx<=dim; leadidx++)
   {
      printf("poles for Pade approximant %d :\n", leadidx);
      for(poleidx=1; poleidx<=dendeg; poleidx++)
      {
         fail = padcon_get_dobldobl_pole(leadidx,poleidx,0,&re,&im);
         printf(" %d : %.14e  %.14e\n", poleidx, re, im);
      }
   }
}

void show_quaddobl_poles ( int dim )
{
   int fail,dendeg,leadidx,poleidx;
   double val,re,im;

   fail = padcon_get_homotopy_continuation_parameter(3,&val);
   dendeg = (int) val;

   for(leadidx=1; leadidx<=dim; leadidx++)
   {
      printf("poles for Pade approximant %d :\n", leadidx);
      for(poleidx=1; poleidx<=dendeg; poleidx++)
      {
         fail = padcon_get_quaddobl_pole(leadidx,poleidx,0,&re,&im);
         printf(" %d : %.14e  %.14e\n", poleidx, re, im);
      }
   }
}

void show_standard_solution ( int index )
{
   int fail,strlensol;

   fail = padcon_get_standard_solution(index,0); // 1);
   fail = solcon_length_standard_solution_string(index,&strlensol);
   {
      char idxsol[strlensol+1];

      fail = solcon_write_standard_solution_string(index,strlensol,idxsol);
      printf("Solution %d :\n%s\n", index, idxsol);
   }
}

void show_dobldobl_solution ( int index )
{
   int fail,strlensol;

   fail = padcon_get_dobldobl_solution(index,0); // 1);
   fail = solcon_length_dobldobl_solution_string(index,&strlensol);
   {
      char idxsol[strlensol+1];

      fail = solcon_write_dobldobl_solution_string(index,strlensol,idxsol);
      printf("Solution %d :\n%s\n", index, idxsol);
   }
}

void show_quaddobl_solution ( int index )
{
   int fail,strlensol;

   fail = padcon_get_quaddobl_solution(index,0); // 1);
   fail = solcon_length_quaddobl_solution_string(index,&strlensol);
   {
      char idxsol[strlensol+1];

      fail = solcon_write_quaddobl_solution_string(index,strlensol,idxsol);
      printf("Solution %d :\n%s\n", index, idxsol);
   }
}

void standard_next_loop ( int index )
{
   int fail,length,failed,dim;
   char contstep,nlsb;
   double frp,re_pole,im_pole,tval,step,pstep,dstep,eta;

   fail = padcon_initialize_standard_solution(index,1);
   fail = solcon_dimension_of_standard_solutions(&dim);
   do
   {
      fail = padcon_standard_predict_correct(&failed,1);
      show_standard_series_coefficients(dim);
      show_standard_pade_coefficients(dim);
      show_standard_poles(dim);
      if(failed != 0)
         contstep = 'n';
      else
      {
         fail = padcon_get_standard_pole_radius(&frp);
         fail = padcon_get_standard_closest_pole(&re_pole,&im_pole);
         fail = padcon_get_standard_t_value(&tval);
         fail = padcon_get_standard_step_size(&step);
         fail = padcon_get_standard_pole_step(&pstep);
         fail = padcon_get_standard_hessian_step(&dstep);
         fail = padcon_get_standard_estimated_distance(&eta);
         printf("pole step : %.3e  Hessian step : %.3e\n", pstep, dstep);
         printf("t : %.3e  step : %.3e\n", tval, step);
         printf("Smallest pole radius : %.3e  eta : %.3e\n", frp, eta);
         printf("Closest pole : %.14e  %.14e\n", re_pole, im_pole);
         show_standard_solution(index);
         printf("continue ? (y/n) ");
         contstep = getchar();
         scanf("%c",&nlsb);
      }
   }
   while(contstep == 'y');
}

void dobldobl_next_loop ( int index )
{
   int fail,length,failed,dim;
   char contstep,nlsb;
   double frp,re_pole,im_pole,tval,step,pstep,dstep,eta;

   fail = padcon_initialize_dobldobl_solution(index,1);
   fail = solcon_dimension_of_dobldobl_solutions(&dim);
   do
   {
      fail = padcon_dobldobl_predict_correct(&failed,1);
      show_dobldobl_series_coefficients(dim);
      show_dobldobl_pade_coefficients(dim);
      show_dobldobl_poles(dim);
      if(failed != 0)
         contstep = 'n';
      else
      {
         fail = padcon_get_dobldobl_pole_radius(&frp);
         fail = padcon_get_dobldobl_closest_pole(&re_pole,&im_pole);
         fail = padcon_get_dobldobl_t_value(&tval);
         fail = padcon_get_dobldobl_step_size(&step);
         fail = padcon_get_dobldobl_pole_step(&pstep);
         fail = padcon_get_dobldobl_hessian_step(&dstep);
         fail = padcon_get_dobldobl_estimated_distance(&eta);
         printf("pole step : %.3e  Hessian step : %.3e\n", pstep, dstep);
         printf("t : %.3e  step : %.3e\n", tval, step);
         printf("Smallest pole radius : %.3e  eta : %.3e\n", frp, eta);
         printf("Closest pole : %.14e  %.14e\n", re_pole, im_pole);
         show_dobldobl_solution(index);
         printf("continue ? (y/n) ");
         contstep = getchar();
         scanf("%c",&nlsb);
      }
   }
   while(contstep == 'y');
}

void quaddobl_next_loop ( int index )
{
   int fail,length,failed,dim;
   char contstep,nlsb;
   double frp,re_pole,im_pole,tval,step,pstep,dstep,eta;

   fail = padcon_initialize_quaddobl_solution(index,1);
   fail = solcon_dimension_of_quaddobl_solutions(&dim);
   do
   {
      fail = padcon_quaddobl_predict_correct(&failed,1);
      show_quaddobl_series_coefficients(dim);
      show_quaddobl_pade_coefficients(dim);
      show_quaddobl_poles(dim);
      if(failed != 0)
         contstep = 'n';
      else
      {
         fail = padcon_get_quaddobl_pole_radius(&frp);
         fail = padcon_get_quaddobl_closest_pole(&re_pole,&im_pole);
         fail = padcon_get_quaddobl_t_value(&tval);
         fail = padcon_get_quaddobl_step_size(&step);
         fail = padcon_get_quaddobl_pole_step(&pstep);
         fail = padcon_get_quaddobl_hessian_step(&dstep);
         fail = padcon_get_quaddobl_estimated_distance(&eta);
         printf("pole step : %.3e  Hessian step : %.3e\n", pstep, dstep);
         printf("t : %.3e  step : %.3e\n", tval, step);
         printf("Smallest pole radius : %.3e  eta : %.3e\n", frp, eta);
         printf("Closest pole : %.14e  %.14e\n", re_pole, im_pole);
         show_quaddobl_solution(index);
         printf("continue ? (y/n) ");
         contstep = getchar();
         scanf("%c",&nlsb);
      }
   }
   while(contstep == 'y');
}

void standard_next_step ( int homo )
{
   int fail,length,index,strlensol,dim,arthom;
   int idx=0;
   char contstep,nlsb;

   arthom = prompt_for_artificial();

   if(arthom == 1)
   {
      fail = read_standard_target_system();
      fail = read_standard_start_system();
      fail = copy_start_solutions_to_container();
   }
   else
   {
      fail = read_standard_target_system();
      fail = copy_target_solutions_to_container();

      idx = prompt_for_parameter();
      fail = solcon_standard_drop_coordinate_by_index(idx);
      fail = solcon_standard_set_continuation_parameter();
   }
   fail = solcon_number_of_standard_solutions(&length);
   fail = solcon_dimension_of_standard_solutions(&dim);
   printf("Read %d start solutions in dimension %d.\n", length, dim);

   if(idx == 0)
      fail = padcon_standard_initialize_homotopy(1,homo);
   else
      fail = padcon_standard_initialize_parameter_homotopy(idx,1);

   for(index=1; index <= length; index++)
   {
      printf("\nTracking path %d ...\n", index);
      standard_next_loop(index);
      if(index == length) break;
      printf("Continue with the next solution ? (y/n) ");
      contstep = getchar();
      scanf("%c",&nlsb);
      if(contstep != 'y') break;
   }
}

void dobldobl_next_step ( int homo )
{
   int fail,length,index,strlensol,dim,arthom;
   int idx=0;
   char contstep,nlsb;

   arthom = prompt_for_artificial();

   if(arthom == 1)
   {
      fail = read_dobldobl_target_system();
      fail = read_dobldobl_start_system();
      fail = copy_dobldobl_start_solutions_to_container();
   }
   else
   {
      fail = read_dobldobl_target_system();
      fail = copy_dobldobl_target_solutions_to_container();

      idx = prompt_for_parameter();
      fail = solcon_dobldobl_drop_coordinate_by_index(idx);
      fail = solcon_dobldobl_set_continuation_parameter();
   }
   fail = solcon_number_of_dobldobl_solutions(&length);
   fail = solcon_dimension_of_dobldobl_solutions(&dim);
   printf("Read %d start solutions of dimension %d.\n", length, dim);

   if(idx == 0)
      fail = padcon_dobldobl_initialize_homotopy(1,homo);
   else
      fail = padcon_dobldobl_initialize_parameter_homotopy(idx,1);

   for(index = 1; index <= length; index++)
   {
      printf("\nTracking path %d ...\n", index);
      dobldobl_next_loop(index);
      if(index == length) break;
      printf("Continue with the next solution ? (y/n) ");
      contstep = getchar();
      scanf("%c",&nlsb);
      if(contstep != 'y') break;
   }
}

void quaddobl_next_step ( int homo )
{
   int fail,length,index,failed,strlensol,dim,arthom;
   int idx=0;
   char contstep,nlsb;
   double frp,re_pole,im_pole,tval,step;

   arthom = prompt_for_artificial();

   if(arthom == 1)
   {
      fail = read_quaddobl_target_system();
      fail = read_quaddobl_start_system();
      fail = copy_quaddobl_start_solutions_to_container();
   }
   else
   {
      fail = read_quaddobl_target_system();
      fail = copy_quaddobl_target_solutions_to_container();

      idx = prompt_for_parameter();
      fail = solcon_quaddobl_drop_coordinate_by_index(idx);
      fail = solcon_quaddobl_set_continuation_parameter();
   }
   fail = solcon_number_of_quaddobl_solutions(&length);
   fail = solcon_dimension_of_quaddobl_solutions(&dim);
   printf("Read %d start solutions of dimension %d.\n", length, dim);

   if(idx == 0)
      fail = padcon_quaddobl_initialize_homotopy(1,homo);
   else
      fail = padcon_quaddobl_initialize_parameter_homotopy(idx,1);

   for(index = 1; index <= length; index++)
   {
      printf("\nTracking path %d ...\n", index);
      quaddobl_next_loop(index);
      if(index == length) break;
      printf("Continue with the next solution ? (y/n) ");
      contstep = getchar();
      scanf("%c",&nlsb);
      if(contstep != 'y') break;
   }
}
