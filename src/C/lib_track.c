/* Reads a target and start system and then solves the target system,
   using the start system in an artificial-parameter homotopy. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "syscon.h"
#include "solcon.h"
#include "phcpack.h"
#include "jump_track.h"
#include "witset.h"

void ask_menu_selection ( int *kind, int *monitor );
/*
 * DESCRIPTION :
 *   Presents the user with a menu asking about the type of start system
 *   and asks whether the user wants to monitor the progress on screen.
 *
 * ON RETURN :
 *   kind      represents the type of start system:
 *             = 1 for total degree start system,
 *             = 2 for linear-product start system,
 *             = 3 for cheater's homotopy start system,
 *             = 4 for one-level-down cascade homotopy,
 *             = 5 for starting a diagonal homotopy,
 *             = 6 for computing a witness set for a hypersurface,
 *             = 7 for collapsing the extrinsic diagonal;
 *             = 8 for removing the last slack variable;
 *   monitor   indicates whether the user wants output to screen:
 *             = 1 if the user want to monitor progress on screen,
 *             = 0 if not output to screen is wanted. */

int read_file_names ( int kind );
/*
 * DESCRIPTION :
 *   The user is asked about the location of the input files
 *   and requested to make choices about the path tracking process.
 *
 * ON ENTRY :
 *   kind      a number between 1 and 4, indicating the type of homotopy.
 *
 * ON RETURN :
 *   kind and monitor are returned by ask_menu_selection above; 
 *   the function returns the value of the fail parameter. */

int create_the_homotopy ( int kind );
/*
 * DESCRIPTION :
 *   Depending on the kind of homotopy, the start and target
 *   system are defined and written to the output file. 
 *
 * ON ENTRY :
 *   kind      must be a natural number in the range 1..4;
 *
 * REQUIRED :
 *   The user has given the appropriate data corresponding to
 *   the kind of homotopy and the output file has been defined. */

int create_the_diagonal_homotopy ( int a, int b );
/*
 * DESCRIPTION :
 *   Creates the diagonal homotopy to start the cascade to intersect
 *   two witness sets of dimensions a and b respectively.
 *   The target and start system are written to file. */

int tune_parameters_and_output ( void );
/*
 * DESCRIPTION :
 *   Allows the user to tune the continuation parameters
 *   and to determine the level of output during path tracking. */

int run_path_trackers ( int kind, int monitor );
/*
 * DESCRIPTION :
 *   Creates a homotopy based on the available start and target system,
 *   scans the dimensions from the solution file and tracks paths.
 *
 * ON ENTRY :
 *   kind and monitor are obtained via ask_menu_selection above.
 *
 * ON RETURN :
 *   Returns the value of the returning fail parameter. */

int prepare_for_the_next_level_down ( int dimension );
/*
 * DESCRIPTION :
 *   Interactive routine that prompts the user for the system
 *   that was just solved, removes then one embedding variable,
 *   and sets up the cascade homotopy for the next level down.
 *   The dimension is the current value for the ambient dimension. */

int run_cascade ( int monitor );
/*
 * DESCRIPTION :
 *   Runs the cascade one level down.
 *
 * ON ENTRY :
 *   monitor is flag to indicate if user wants to monitor
 *           the progress of the computations on screen.
 *
 * ON RETURN :
 *   the value of the fail parameter. */

int read_two_witness_sets
      ( int *n1, int *n2, int *dim1, int *dim2, int *deg1, int *deg2,
        int *cd );
/*
 * DESCRIPTION :
 *   The user is prompted for two witness sets.
 *
 * ON RETURN :
 *   n1        ambient dimension for the first witness set;
 *   n2        ambient dimension for the second witness set;
 *   dim1      dimension of the solutions represented by the 1st witness set;
 *   dim2      dimension of the solutions represented by the 2nd witness set;
 *   deg1      degree of the 1st witness set, i.e.: #solutions in 1st set;
 *   deg2      degree of the 2nd witness set, i.e.: #solutions in 2nd set;
 *   cd        cascade dimension: #vars and #eqs in the homotopy. */

int intersect_two_witness_sets
      ( int n1, int n2, int dim1, int dim2, int deg1, int deg2,
        int cd, int monitor );
/*
 * DESCRIPTION :
 *   Runs the path trackers to start the extrinsic cascade to
 *   intersect two witness sets.
 *
 * ON ENTRY :
 *   n1       ambient dimension of the first witness set;
 *   n2       ambient dimension of the second witness set;
 *   dim1     dimension of the first solution component;
 *   dim2     dimension of the second solution component;
 *   deg1     degree of the first solution set;
 *   deg2     degree of the second solution set;
 *   cd       cascade dimension, size of the homotopy;
 *   monitor  if 0, then no extra output will be written to screen,
 *            otherwise, the user can monitor the progress. */

int witness_set_for_hypersurface ( void );
/*
 * DESCRIPTION :
 *   Computes the witness set for one hypersurface. */

int collapse_extrinsic_diagonal ( void );
/*
 * DESCRIPTION :
 *   Collapses the diagonal in a witness set obtained
 *   via an extrinsic diagonal homotopy. */

int remove_last_slack_variable ( void );
/*
 * DESCRIPTION :
 *   Removes the last slack variable of a witness set.
 *   This operation is typically done before moving to
 *   the next level down in a cascade. */

int test_standard_crude_tracker ( int verbose );
/*
 * DESCRIPTION :
 *   Tests the crude path tracker in standard double precision. */

int test_dobldobl_crude_tracker ( int verbose );
/*
 * DESCRIPTION :
 *   Tests the crude path tracker in double double precision. */

int test_quaddobl_crude_tracker ( int verbose );
/*
 * DESCRIPTION :
 *   Tests the crude path tracker in quad double precision. */

int test_crude_trackers ( int verbose );
/*
 * DESCRIPTION :
 *   Prompts the precision and then calls the crude path trackers
 *   in double, double double, or quad double precision. */

int main ( int argc, char *argv[] )
{
   int fail,kind,monitor,n1,n2,a,b,d1,d2,cd;

   adainit();

   printf("\nCalling the path trackers in PHCpack...\n");

   ask_menu_selection(&kind,&monitor);

   if(kind == 9)
   {
      fail = test_crude_trackers(monitor);
      adafinal(); return fail;
   }
   if(kind == 6)
   {
      fail = witness_set_for_hypersurface();
      adafinal(); return fail;
   }
   if(kind == 7)
   {
      fail = collapse_extrinsic_diagonal();
      adafinal(); return fail;
   }
   if(kind == 8)
   {
      fail = remove_last_slack_variable();
      adafinal(); return fail;
   }

   if(kind < 5)
   {
      fail = read_file_names(kind);
      fail = create_the_homotopy(kind);
   }
   else
   {
      fail = read_two_witness_sets(&n1,&n2,&a,&b,&d1,&d2,&cd);
      fail = create_the_diagonal_homotopy(a,b); 
   }

   if(fail > 0)
      printf("creation of homotopy fail = %d, no path tracking ...\n",fail);
   else
   {
      fail = tune_parameters_and_output();
      if(kind == 4)
         fail = run_cascade(monitor);
      else if(kind < 5)
         fail = run_path_trackers(kind,monitor);
      else
         fail = intersect_two_witness_sets(n1,n2,a,b,d1,d2,cd,monitor);
   }

   adafinal();

   return 0;
}

void ask_menu_selection ( int *kind, int *monitor )
{
   char ans;

   printf("\nMENU for type of start system and homotopy :\n");
   printf("  1. start system is based on total degree;\n");
   printf("  2. a linear-product start system will be given;\n");
   printf("  3. start system and start solutions are provided;\n");
   printf("  4. the homotopy is a cascade to go one level down;\n");
   printf("  5. diagonal homotopy to start a cascade;\n");
   printf("  6. compute witness set for a hypersurface;\n");
   printf("  7. collapse extrinsic diagonal of a witness set;\n");
   printf("  8. remove the last slack variable of a witness set;\n");
   printf("  9. test the crude path trackers.\n");
   printf("Type 1, 2, 3, 4, 5, 6, 7, 8, or 9 to select : ");
   scanf("%d",kind);
   if(*kind < 6)
   {
      printf("\nDo you want to monitor progress of solver on screen ? (y/n) ");
      scanf("%c",&ans);  /* skip previous new line symbol */
      scanf("%c",&ans);  /* get answer of the user */
      if(ans = 'y')
         *monitor = 1;
      else
         *monitor = 0;
      scanf("%c",&ans);  /* skip current new line symbol */
   }
   if(*kind == 9) *monitor = 1; /* verbose mode default when testing */
}

int read_file_names ( int kind )
{
   int fail,n,dim;
   char start[80],ch;

   if(kind < 4) fail = read_target_system_without_solutions();

   printf("\nReading the name of the file for the start system...");
   printf("\nGive a string of characters : "); scanf("%s",start);

   scanf("%c",&ch); /* skip previous newline symbol */

   fail = define_output_file();

   n = (int) strlen(start);
   if(kind == 2)
      fail = read_named_linear_product_start_system(n,start);
   else 
      fail = read_named_start_without_solutions(n,start);

   if(kind >= 3) fail = solcon_scan_solution_banner();

   if(kind == 4)
   {
      fail = copy_start_system_to_container();
      fail = syscon_sort_embed_symbols(&dim);
      printf("the top dimension is %d\n",dim);
      fail = copy_container_to_start_system();
   }

   return fail;
}

int create_the_homotopy ( int kind )
{
   int fail;

   if(kind < 4)
      fail = create_homotopy();
   else if(kind == 4)
      fail = create_cascade_homotopy();
   else
      return 1;

   if(fail == 0)
   {
      fail = write_standard_target_system();
      fail = write_standard_start_system();
   }

   return fail;
}

int create_the_diagonal_homotopy ( int a, int b )
{
   int fail;

   fail = standard_diagonal_homotopy(a,b);
   if(fail == 0)
   {
      fail = write_standard_target_system();
      fail = write_standard_start_system();
   }
   return fail;
}

int tune_parameters_and_output ( void )
{
   int fail;

   printf("\n");
   fail = tune_continuation_parameters();
   printf("\n");
   fail = determine_output_during_continuation();
   printf("\nSee the output file for results ...\n\n");

   return fail;
}

int run_path_trackers ( int kind, int monitor )
{
   int fail,len,dim,m,i,nbstep,nbfail,nbiter,nbsyst,cnt,ind;
   double *sol;

   if(kind >= 3)
      fail = solcon_read_solution_dimensions(&len,&dim);
   else
   {
      fail = copy_target_system_to_container();
      fail = syscon_number_of_standard_polynomials(&dim);
      fail = syscon_total_degree(&len);
   }

   fail = solcon_write_solution_banner_to_defined_output_file();
   fail = solcon_write_solution_dimensions_to_defined_output_file(len,dim);

   sol = (double*)calloc(2*dim+5,sizeof(double));
   cnt = 0;
   for(i=1; i<=len; i++)
   {
      if(kind == 1)
         fail = solcon_compute_total_degree_solution(dim,i,&m,sol);
      else if(kind == 2)
      {
         i = i-1;
         fail = solcon_next_linear_product_solution(dim,&i,&m,sol);
      }
      else if(kind >= 3)
         fail = solcon_read_next_solution(dim,&m,sol);
      else
         fail = 1;

      if(fail == 0)
      {
         fail = reporting_path_tracker
                  (dim,&m,sol,&nbstep,&nbfail,&nbiter,&nbsyst);
         if(monitor == 1)
            printf("%d : #step : %3d #fail : %3d #iter : %3d #syst : %3d\n",
                   i,nbstep,nbfail,nbiter,nbsyst);
         ind = i;
         fail = write_next_solution_with_diagnostics
                  (&ind,dim,m,sol,nbstep,nbfail,nbiter,nbsyst);
         cnt++;
      }
      else
      {
         if(monitor == 1) printf("%d : no path\n",i); break;
      }
   }
   printf("Wrote %d solutions to file.\n",cnt);
   return fail;
}

int prepare_for_next_level_down ( int dimension )
{
   char start[80],ch;
   int n,fail;

   printf("\nremoving last slack variable first ...\n");
   remove_last_slack_variable();
   fail = syscon_remove_symbol_from_table(dimension);
   printf("\nnow ready for the next level down ...\n");
   printf("\nReading the name of the file for the start system...");
   printf("\nGive a string of characters : "); scanf("%s",start);
   scanf("%c",&ch); /* skip previous newline symbol */
   fail = define_output_file();
   n = (int) strlen(start);
 /*  fail = syscon_clear_symbol_table(); */
   fail = read_named_start_without_solutions(n,start);
   fail = solcon_scan_solution_banner();
   printf("\ncreating the cascade homotopy ...\n\n");
   fail = create_cascade_homotopy();
   fail = write_standard_target_system();
   fail = write_standard_start_system();
}

int run_cascade ( int monitor )
{
   int fail,len,dim,m,i,nbstep,nbfail,nbiter,nbsyst,cnt,ind;
   char ans = 'y';
   double *sol;
   while(ans == 'y')
   {
      fail = solcon_read_solution_dimensions(&len,&dim);
      fail = solcon_write_solution_banner_to_defined_output_file();
      fail = solcon_write_solution_dimensions_to_defined_output_file(len,dim);
      sol = (double*)calloc(2*dim+5,sizeof(double));
      cnt = 0;
      for(i=1; i<=len; i++)
      {
         fail = solcon_read_next_solution(dim,&m,sol);
         if(fail == 0)
         {
            fail = reporting_path_tracker
                     (dim,&m,sol,&nbstep,&nbfail,&nbiter,&nbsyst);
            if(monitor == 1)
               printf("%d : #step : %3d #fail : %3d #iter : %3d #syst : %3d\n",
                      i,nbstep,nbfail,nbiter,nbsyst);
            ind = i;
            fail = write_next_solution_with_diagnostics
                     (&ind,dim,m,sol,nbstep,nbfail,nbiter,nbsyst);
            cnt++;
         }
         else
         {
            if(monitor == 1) printf("%d : no path\n",i); break;
         }
      }
      printf("Wrote %d solutions to file.\n",cnt);
      close_output_file();
      solcon_close_solution_input_file(0);
      printf("\ncontinue ? (y/n) "); ans = getchar();
      if(ans == 'y') prepare_for_next_level_down(dim);
   }
   return fail;
}

void check_on_two_witness_sets
      ( int n1, int n2, int dim1, int dim2, int deg1, int deg2, int cd )
/* this is only for testing purposes: a sanity check on the inputs ... */
{
   double sol1[2*n1+5];
   double sol2[2*n2+5];
   double ps[2*cd+5];
   int fail,m,deg,dim,i,j;

   printf("\nDoing a sanity check on the solutions on file...\n");

   fail = solcon_read_next_witness_point(1,n1,&m,sol1);
   fail = solcon_read_next_witness_point(2,n2,&m,sol2);

   fail = solcon_extrinsic_product(dim1,dim2,n1,n2,sol1,sol2,cd,ps);
   printf("\nThe product of the first two solutions :");
   for(i=0; i<2*cd+2; i++)
   {
      if(i%2 == 0) printf("\n");
      printf("  %.14lf",ps[i]);
   }
   printf("\n");
   fail = solcon_reset_input_file(1,&deg,&dim);
   printf("after resetting 1 : n = %d  degree = %d\n",dim,deg);
   fail = solcon_reset_input_file(2,&deg,&dim);
   printf("after resetting 2 : n = %d  degree = %d\n",dim,deg);

   for(i=0; i<deg1; i++)
   {
      fail = solcon_read_next_witness_point(1,n1,&m,sol1);
      for(j=0; j<deg2; j++)
      {
         fail = solcon_read_next_witness_point(2,n2,&m,sol2);
         printf("%d x %d starts with (%.3lf , %.3lf) \n",
                i,j,sol1[2],sol2[2]);
      }
      if(i<deg1-1) 
      {
         fail = solcon_reset_input_file(2,&deg,&dim);
         printf("after resetting 2 : n = %d  degree = %d\n",dim,deg);
      }
   }
   printf("\n");
   fail = solcon_reset_input_file(1,&deg,&dim);
   printf("after resetting 1 : n = %d  degree = %d\n",dim,deg);
   fail = solcon_reset_input_file(2,&deg,&dim);
   printf("after resetting 2 : n = %d  degree = %d\n",dim,deg);
}

int read_two_witness_sets 
      ( int *n1, int *n2, int *dim1, int *dim2, int *deg1, int *deg2,
        int *cd )
{
   int fail;

   printf("\n");
   fail = read_a_witness_set(1,n1,dim1,deg1);
   printf("  n = %d  dimension = %d  degree = %d\n",*n1,*dim1,*deg1);
   printf("\n");
   fail = read_a_witness_set(2,n2,dim2,deg2);
   printf("  n = %d  dimension = %d  degree = %d\n",*n2,*dim2,*deg2);

   fail = define_output_file();

   fail = extrinsic_top_diagonal_dimension(*n1,*n2,*dim1,*dim2,cd);
   printf("The top dimension of the extrinsic diagonal cascade : %d.\n",*cd);

   /* check_on_two_witness_sets(*n1,*n2,*dim1,*dim2,*deg1,*deg2,*cd); */

   return fail;
}

int intersect_two_witness_sets 
      ( int n1, int n2, int dim1, int dim2, int deg1, int deg2,
	int cd, int monitor )
{
   int fail,m,i,j,k,len,nbstep,nbfail,nbiter,nbsyst,ind;
   double sol1[2*n1+5];
   double sol2[2*n2+5];
   double ps[2*cd+5];

   fail = solcon_write_solution_banner_to_defined_output_file();
   len = deg1*deg2;
   fail = solcon_write_solution_dimensions_to_defined_output_file(len,cd);

   for(k=0,i=0,j=0; k<len; k++)
   {
      if(monitor == 1) printf("%d x %d = %d ... ",i,j,k);
     
      fail = get_next_start_product
               (&i,&j,monitor,n1,n2,dim1,dim2,deg1,deg2,cd,sol1,sol2,ps);
      m = 1;
      fail = reporting_path_tracker(cd,&m,ps,&nbstep,&nbfail,&nbiter,&nbsyst);

      if(monitor == 1)
         printf(" #step : %3d #fail : %3d #iter : %3d #syst : %3d\n",
                nbstep,nbfail,nbiter,nbsyst);

      ind = k;
      fail = write_next_solution_with_diagnostics
               (&ind,cd,m,ps,nbstep,nbfail,nbiter,nbsyst);
   }
   return fail;
}

int witness_set_for_hypersurface ( void )
{
   int fail,k,n;
   char outfile[80];

   fail = syscon_read_standard_system();

   printf("\nReading the name of the output file...");
   printf("\nGive a string of characters : "); scanf("%s",outfile);

   printf("\nThe system in the container :\n");
   fail = syscon_write_standard_system();
   fail = syscon_number_of_standard_polynomials(&n);
   printf("Number of polynomials in container : %d\n",n);
   printf("Give the number of the equation : ");
   scanf("%d",&k);

   n = (int) strlen(outfile);
   fail = hypersurface_witness_set(k,n,outfile);

   return fail;
}

int collapse_extrinsic_diagonal ( void )
{
   int m1,m2,n,dim,deg,add,fail;
   char ch,infile[80],outfile[80];

   scanf("%c",&ch); /* skip newline symbol from input */

   printf("\nReading the name of the file for a witness set...");
   printf("\nGive a string of characters : "); scanf("%s",infile);
   m1 = (int) strlen(infile);
   printf("\nReading the name of the output file...");
   printf("\nGive a string of characters : "); scanf("%s",outfile);
   m2 = (int) strlen(outfile);

   fail = read_witness_set_from_file(m1,infile,&n,&dim,&deg);
   printf("\nThe current number of slack variables : %d\n",dim);
   printf("Give number of slack variables to add : ");
   scanf("%d",&add);
   fail = standard_collapse_diagonal(dim,add);
   fail = write_witness_set_to_file(m2,outfile);

   return fail;
}

int remove_last_slack_variable ( void )
{
   int m1,m2,n,dim,deg,add,fail;
   char ch,infile[80],outfile[80];

   scanf("%c",&ch); /* skip newline symbol from input */

   printf("\nReading the name of the file for a witness set...");
   printf("\nGive a string of characters : "); scanf("%s",infile);
   m1 = (int) strlen(infile);
   printf("\nReading the name of the output file...");
   printf("\nGive a string of characters : "); scanf("%s",outfile);
   m2 = (int) strlen(outfile);

   fail = read_witness_set_from_file(m1,infile,&n,&dim,&deg);
   printf("\nThe current number of slack variables : %d\n",dim);
   fail = remove_last_slack(dim);
   fail = write_witness_set_to_file(m2,outfile);
}

int test_standard_crude_tracker ( int verbose )
{
   int fail;

   printf("\ntesting the crude trackers in double precision ...\n");

   fail = read_standard_target_system();
   fail = copy_container_to_target_system();
   fail = read_standard_start_system();
   fail = copy_container_to_start_system();
   fail = copy_container_to_start_solutions();
   fail = create_homotopy();
   fail = standard_crude_tracker(verbose);

   return 0;
}

int test_dobldobl_crude_tracker ( int verbose )
{
   int fail;

   printf("\ntesting the crude trackers in double double precision ...\n");

   fail = read_dobldobl_target_system();
   fail = copy_dobldobl_container_to_target_system();
   fail = read_dobldobl_start_system();
   fail = copy_dobldobl_container_to_start_system();
   fail = copy_dobldobl_container_to_start_solutions();
   fail = create_dobldobl_homotopy();
   fail = dobldobl_crude_tracker(verbose);

   return 0;
}

int test_quaddobl_crude_tracker ( int verbose )
{
   int fail;

   printf("\ntesting the crude trackers in quad double precision ...\n");

   fail = read_quaddobl_target_system();
   fail = copy_quaddobl_container_to_target_system();
   fail = read_quaddobl_start_system();
   fail = copy_quaddobl_container_to_start_system();
   fail = copy_quaddobl_container_to_start_solutions();
   fail = create_quaddobl_homotopy();
   fail = quaddobl_crude_tracker(verbose);

   return 0;
}

int test_crude_trackers ( int verbose )
{
   int precision = 0;
   char nlch;

   printf("\nMENU for the working precision :\n");
   printf("  0. double precision,\n");
   printf("  1. double double precision, or\n");
   printf("  2. quad double precision.\n");
   printf("Type 0, 1, or 2 to select the precision : ");
   scanf("%d",&precision);
   scanf("%c",&nlch); /* swallow new line character */

   if(precision == 0)
      return test_standard_crude_tracker(verbose);
   else if(precision == 1)
      return test_dobldobl_crude_tracker(verbose);
   else if(precision == 2)
      return test_quaddobl_crude_tracker(verbose);
   else
      printf("Invalid value for the precision.");

   return 0;
}
