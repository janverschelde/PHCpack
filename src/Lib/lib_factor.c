/* The main program prompts the user for a witness set and then runs the
   basic version of the monodromy loop algorithm to decompose the witness
   set into irreducible factors.  Double double and quad double versions
   of the monodromy breakup method are provided as well. */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "phcpack.h"
#include "syscon.h"
#include "solcon.h"
#include "witset.h"

#define verbose 1 /* verbose flag */

int set_standard_trace_slice ( int first );
/*
 * DESCRIPTION :
 *   Sets the constant coefficient of the slice used in the linear trace,
 *   in standard double precision,
 *   when called for the first time, first must have the value 1. */

int set_dobldobl_trace_slice ( int first );
/*
 * DESCRIPTION :
 *   Sets the constant coefficient of the slice used in the linear trace,
 *   in standard double precision,
 *   when called for the first time, first must have the value 1. */

int set_quaddobl_trace_slice ( int first );
/*
 * DESCRIPTION :
 *   Sets the constant coefficient of the slice used in the linear trace,
 *   in quad double precision,
 *   when called for the first time, first must have the value 1. */

int store_standard_gammas ( int n );
/*
 * DESCRIPTION :
 *   Generates n random complex coefficients in standard double precision
 *   and stores these coefficients into the monodromy breakup machine. */

int standard_monodromy_breakup
 ( int islaurent, int nbloops, int n, int k, int d );
/*
 * DESCRIPTION :
 *   Uses at most nbloops to break up a k-dimensional set of degree d
 *   in n-space, in standard double precision.
 *   If the witness set is defined by a Laurent polynomial system,
 *   then islaurent must be equal to one. */

int dobldobl_monodromy_breakup
 ( int islaurent, int nbloops, int n, int k, int d );
/*
 * DESCRIPTION :
 *   Uses at most nbloops to break up a k-dimensional set of degree d 
 *   in n-space, in double double precision.
 *   If the witness set is defined by a Laurent polynomial system,
 *   then islaurent must be equal to one. */ 

int quaddobl_monodromy_breakup
 ( int islaurent, int nbloops, int n, int k, int d );
/*
 * DESCRIPTION :
 *   Uses at most nbloops to break up a k-dimensional set of degree d
 *   in n-space, in quad double precision.
 *   If the witness set is defined by a Laurent polynomial system,
 *   then islaurent must be equal to one. */ 

int standard_test ( int islaurent );
/*
 * DESCRIPTION :
 *   Performs the monodromy breakup in standard double precision.
 *   If islaurent equals one , then the witness set is defined by
 *   a Laurent polynomial system. */

int dobldobl_test ( int islaurent );
/*
 * DESCRIPTION :
 *   Performs the monodromy breakup in double double precision.
 *   If islaurent equals one, then the witness set is defined by
 *   a Laurent polynomial system. */

int quaddobl_test ( int isluarent );
/*
 * DESCRIPTION :
 *   Performs the monodromy breakup in quad double precision.
 *   If islaurent equals one, then the witness set is defined by
 *   a Laurent polynomial system. */

int main ( int argc, char *argv[] )
{
   char ans;
   int precision,fail,laurent;
   int seed = time(NULL);

   adainit();
   srand(seed); /* monodromy loops are probabilistic ! */

   printf("Seed for the random number generators : %d.\n",seed);

   printf("\nMENU for the working precision :\n");
   printf("  0. standard double precision;\n");
   printf("  1. double double precision;\n");
   printf("  2. quad double precision.\n");
   printf("Type 0, 1, or 2 to make a choice : ");
   scanf("%d",&precision);

   scanf("%c",&ans); /* skip end of line character */

   printf("\nLaurent polynomial system ? (y/n) ");
   scanf("%c",&ans);
   laurent = (ans == 'y');
   scanf("%c",&ans); /* skip end of line character */

   if(precision == 0)
      fail = standard_test(laurent);
   else if(precision == 1)
      fail = dobldobl_test(laurent);
   else if(precision == 2)
      fail = quaddobl_test(laurent);
   else
      printf("Selected precision level is not supported.\n");

   adafinal();

   return 0;
}

int set_standard_trace_slice ( int first )
{
   int fail;
   double r[2];

   r[1] = 0.0;
   if(first == 1)                  /* determine constant coefficient */
      r[0] = -1.0;
   else
      r[0] = +1.0;

   fail = assign_standard_coefficient_of_slice(0,0,r);

   return fail;
}

int set_dobldobl_trace_slice ( int first )
{
   int fail;
   double r[4];

   r[1] = 0.0; r[2] = 0.0; r[3] = 0.0;
   if(first == 1)                  /* determine constant coefficient */
      r[0] = -1.0;
   else
      r[0] = +1.0;

   fail = assign_dobldobl_coefficient_of_slice(0,0,r);

   return fail;
}

int set_quaddobl_trace_slice ( int first )
{
   int fail;
   double r[8];

   r[1] = 0.0; r[2] = 0.0; r[3] = 0.0; r[4] = 0.0;
   r[5] = 0.0; r[6] = 0.0; r[7] = 0.0;
   if(first == 1)                  /* determine constant coefficient */
      r[0] = -1.0;
   else
      r[0] = +1.0;

   fail = assign_quaddobl_coefficient_of_slice(0,0,r);

   return fail;
}

int store_standard_gammas ( int n )
{
   double re_gamma[n];
   double im_gamma[n];
   int i;
    
   for(i=0; i<n; i++)
      random_complex(&re_gamma[i],&im_gamma[i]);
   store_standard_gamma(n,re_gamma,im_gamma);
}

int store_dobldobl_gammas ( int n )
{
   double re_gamma[2*n];
   double im_gamma[2*n];
   int i;
    
   for(i=0; i<n; i++)
      random_dobldobl_complex(&re_gamma[2*i],&im_gamma[2*i]);
   store_dobldobl_gamma(n,re_gamma,im_gamma);
}

int store_quaddobl_gammas ( int n )
{
   double re_gamma[4*n];
   double im_gamma[4*n];
   int i;
    
   for(i=0; i<n; i++)
      random_quaddobl_complex(&re_gamma[4*i],&im_gamma[4*i]);
   store_quaddobl_gamma(n,re_gamma,im_gamma);
}

int standard_monodromy_breakup 
 ( int islaurent, int nbloops, int n, int k, int d )
{
   int fail,i,j,done;
   double err,dis;

   if(islaurent == 1)
      fail = initialize_standard_Laurent_sampler(k);
   else
      fail = initialize_standard_sampler(k);
   fail = initialize_standard_monodromy(nbloops,d,k);
   if(verbose>0)
      printf("\nInitialized sampler and monodromy permutations ...\n");

   fail = store_standard_solutions();  /* store in monodromy permutations */

   printf("... initializing the grid in standard double precision ...\n");
   for(i=1; i<=2; i++)        /* initialize grid for trace validation */
   {
      fail = set_standard_trace_slice(i);   /* fix constant of slice */
      fail = store_standard_gammas(n);  /* generate random gamma constants */
      fail = standard_track_paths(islaurent);
      fail = store_standard_solutions();   /* store solutions in the grid */
      fail = restore_standard_solutions(); /* use original solutions at start */
      fail = swap_standard_slices();       /* go back to original slices */
   }
   fail = standard_trace_grid_diagnostics(&err,&dis);
   printf("Trace grid diagnostics : \n");
   printf("  largest error of the samples : %.3e\n", err);
   printf("  smallest distance between samples : %.3e\n", dis);
   done = 0;
   for(i=1; (i<=nbloops) && (done==0); i++)  /* perform monodromy loops */
   {
      printf("... starting loop #\%d in standard double precision ...\n",i);
      fail = new_standard_slices(k,n);
      fail = store_standard_gammas(n);
      fail = standard_track_paths(islaurent);  /* swapping slices happens */
      fail = solcon_clear_standard_solutions();
      fail = store_standard_gammas(n);
      fail = standard_track_paths(islaurent);
      fail = store_standard_solutions();
      {
         int permutation[d];
         fail = permutation_after_standard_loop(d,permutation);
         printf("the permutation at step %d:",i);
         for(j=0; j<d; j++) printf(" %d",permutation[j]);
         printf("\n");
      }
      fail = standard_monodromy_permutation(d,&done);
      {
         int nf,deg,jp,w[d];
         double trace_difference;
         fail = number_of_standard_factors(&nf);
         printf("number of irreducible factors : %d\n",nf);
         for(j=1; j<=nf; j++)
         {
            fail = witness_points_of_standard_factor(j,&deg,w);
            printf(" %d :",j);
            for(jp=0; jp<deg; jp++) printf(" %d",w[jp]);
            fail = standard_trace_sum_difference(deg,w,&trace_difference);
            printf(" : %.3e\n",trace_difference);
         }
      }
      fail = restore_standard_solutions();
   }
   if(done==1)
      printf("Found factorization using %d loops.\n", i-1);
   else
      printf("Failed to factor using %d loops.\n", i-1);

   return fail;
}

int dobldobl_monodromy_breakup
 ( int islaurent, int nbloops, int n, int k, int d )
{
   int fail,i,j,done;
   double err,dis;

   if(islaurent == 1)
      fail = initialize_dobldobl_Laurent_sampler(k);
   else
      fail = initialize_dobldobl_sampler(k);
   fail = initialize_dobldobl_monodromy(nbloops,d,k);
   if(verbose>0)
      printf("\nInitialized sampler and monodromy permutations ...\n");

   fail = store_dobldobl_solutions();  /* store in monodromy permutations */

   printf("... initializing the grid in double double precision...\n");
   for(i=1; i<=2; i++)        /* initialize grid for trace validation */
   {
      fail = set_dobldobl_trace_slice(i);   /* fix constant of slice */
      fail = store_dobldobl_gammas(n);  /* generate random gamma constants */
      fail = dobldobl_track_paths(islaurent);
      fail = store_dobldobl_solutions();   /* store solutions in the grid */
      fail = restore_dobldobl_solutions(); /* use original sols at start */
      fail = swap_dobldobl_slices();       /* go back to original slices */
   }
   fail = dobldobl_trace_grid_diagnostics(&err,&dis);
   printf("Trace grid diagnostics : \n");
   printf("  largest error of the samples : %.3e\n", err);
   printf("  smallest distance between samples : %.3e\n", dis);
   done = 0;
   for(i=1; (i<=nbloops) && (done==0); i++)  /* perform monodromy loops */
   {
      printf("... starting loop #\%d in double double precision ...\n",i);
      fail = new_dobldobl_slices(k,n);
      fail = store_dobldobl_gammas(n);
      fail = dobldobl_track_paths(islaurent); /* swapping slices happens */
      fail = solcon_clear_dobldobl_solutions();
      fail = store_dobldobl_gammas(n);
      fail = dobldobl_track_paths(islaurent);
      fail = store_dobldobl_solutions();
      {
         int permutation[d];

         fail = permutation_after_dobldobl_loop(d,permutation);
         printf("the permutation at step %d:",i);
         for(j=0; j<d; j++) printf(" %d",permutation[j]);
         printf("\n");
      }
      fail = dobldobl_monodromy_permutation(d,&done);
      {
         int nf,deg,jp,w[d];
         double trace_difference;

         fail = number_of_dobldobl_factors(&nf);
         printf("number of irreducible factors : %d\n",nf);
         for(j=1; j<=nf; j++)
         {
            fail = witness_points_of_dobldobl_factor(j,&deg,w);
            printf(" %d :",j);
            for(jp=0; jp<deg; jp++) printf(" %d",w[jp]);
            fail = dobldobl_trace_sum_difference(deg,w,&trace_difference);
            printf(" : %.3e\n",trace_difference);
         }
      }
      fail = restore_dobldobl_solutions();
   }
   if(done==1)
      printf("Found factorization using %d loops.\n", i-1);
   else
      printf("Failed to factor using %d loops.\n", i-1);

   return fail;
}

int quaddobl_monodromy_breakup
 ( int islaurent, int nbloops, int n, int k, int d )
{
   int fail,i,j,done;
   double err,dis;

   if(islaurent == 1)
      fail = initialize_quaddobl_Laurent_sampler(k);
   else
      fail = initialize_quaddobl_sampler(k);
   fail = initialize_quaddobl_monodromy(nbloops,d,k);
   if(verbose>0)
      printf("\nInitialized sampler and monodromy permutations ...\n");

   fail = store_quaddobl_solutions();  /* store in monodromy permutations */

   printf("... initializing the grid in quad double precision ...\n");
   for(i=1; i<=2; i++)        /* initialize grid for trace validation */
   {
      fail = set_quaddobl_trace_slice(i);   /* fix constant of slice */
      fail = store_quaddobl_gammas(n);  /* generate random gamma constants */
      fail = quaddobl_track_paths(islaurent);
      fail = store_quaddobl_solutions();   /* store solutions in the grid */
      fail = restore_quaddobl_solutions(); /* use original sols at start */
      fail = swap_quaddobl_slices();       /* go back to original slices */
   }
   fail = quaddobl_trace_grid_diagnostics(&err,&dis);
   printf("Trace grid diagnostics : \n");
   printf("  largest error of the samples : %.3e\n", err);
   printf("  smallest distance between samples : %.3e\n", dis);
   done = 0;
   for(i=1; (i<=nbloops) && (done==0); i++)  /* perform monodromy loops */
   {
      printf("... starting loop #\%d in quad double precision ...\n",i);
      fail = new_quaddobl_slices(k,n);
      fail = store_quaddobl_gammas(n);
      fail = quaddobl_track_paths(islaurent); /* swapping slices happens */
      fail = solcon_clear_quaddobl_solutions();
      fail = store_quaddobl_gammas(n);
      fail = quaddobl_track_paths(islaurent);
      fail = store_quaddobl_solutions();
      {
         int permutation[d];

         fail = permutation_after_quaddobl_loop(d,permutation);
         printf("the permutation at step %d:",i);
         for(j=0; j<d; j++) printf(" %d",permutation[j]);
         printf("\n");
      }
      fail = quaddobl_monodromy_permutation(d,&done);
      {
         int nf,deg,jp,w[d];
         double trace_difference;

         fail = number_of_quaddobl_factors(&nf);
         printf("number of irreducible factors : %d\n",nf);
         for(j=1; j<=nf; j++)
         {
            fail = witness_points_of_quaddobl_factor(j,&deg,w);
            printf(" %d :",j);
            for(jp=0; jp<deg; jp++) printf(" %d",w[jp]);
            fail = quaddobl_trace_sum_difference(deg,w,&trace_difference);
            printf(" : %.3e\n",trace_difference);
         }
      }
      fail = restore_quaddobl_solutions();
   }
   if(done==1)
      printf("Found factorization using %d loops.\n", i-1);
   else
      printf("Failed to factor using %d loops.\n", i-1);

   return fail;
}

int standard_test ( int islaurent )
{
   int fail,n,dim,deg,nbloops,kind;
   char ans;

   printf("\nReading a witness set ...\n");
   if(islaurent == 1)
      fail = read_standard_Laurent_witness_set(&n,&dim,&deg);
   else
      fail = read_witness_set(&n,&dim,&deg);

   printf("\nMENU for the kind of output :\n");
   printf("  0. remain silent with no intermediate output;\n");
   printf("  1. all intermediate output goes to screen;\n");
   printf("  2. give a file name for all intermediate output.\n");
   printf("Type 0, 1, or 2 to make a choice : "); scanf("%d",&kind);

   scanf("%c",&ans); /* skip end of line character */

   if(kind == 0) fail = set_standard_state_to_silent();
   if(kind == 2) fail = define_output_file();

   if(verbose>0)  /* only in verbose mode */
   {
      printf("\nThe ambient dimension : %d.\n",n);
      printf("The dimension of the solution set : %d.\n",dim);
      printf("The degree of the solution set : %d.\n",deg);
      printf("\nDo you wish to see the embedded system ? (y/n) ");
      scanf("%c",&ans);
      if(ans == 'y')
      {
         printf("\nThe system read :\n");
         if(islaurent == 1)
            fail = syscon_write_standard_Laurent_system();
         else
            fail = syscon_write_standard_system();
      }
   }

   fail = assign_labels(n,deg,0);

   printf("\nGive the maximum number of loops : ");
   scanf("%d",&nbloops);

   fail = standard_monodromy_breakup(islaurent,nbloops,n,dim,deg);

   if(kind == 2) printf("See the output file for results.\n");

   return 0;
}

int dobldobl_test ( int islaurent )
{
   int fail,n,dim,deg,nbloops,kind;
   char ans;

   printf("\nReading a witness set ...\n");
   if(islaurent == 1)
      fail = read_dobldobl_Laurent_witness_set(&n,&dim,&deg);
   else
      fail = read_dobldobl_witness_set(&n,&dim,&deg);

   printf("\nMENU for the kind of output :\n");
   printf("  0. remain silent with no intermediate output;\n");
   printf("  1. all intermediate output goes to screen;\n");
   printf("  2. give a file name for all intermediate output.\n");
   printf("Type 0, 1, or 2 to make a choice : "); scanf("%d",&kind);

   scanf("%c",&ans); /* skip end of line character */

   if(kind == 0) fail = set_dobldobl_state_to_silent();
   if(kind == 2) fail = define_output_file();

   if(verbose>0)  /* only in verbose mode */
   {
      printf("\nThe ambient dimension : %d.\n",n);
      printf("The dimension of the solution set : %d.\n",dim);
      printf("The degree of the solution set : %d.\n",deg);
      printf("\nDo you wish to see the embedded system ? (y/n) ");
      scanf("%c",&ans);
      if(ans == 'y')
      {
         printf("\nThe system read :\n");
         if(islaurent == 1)
            fail = syscon_write_dobldobl_Laurent_system();
         else
            fail = syscon_write_dobldobl_system();
      }
   }
   fail = dobldobl_assign_labels(n,deg);

   printf("\nGive the maximum number of loops : ");
   scanf("%d",&nbloops);

   fail = dobldobl_monodromy_breakup(islaurent,nbloops,n,dim,deg);

   if(kind == 2) printf("See the output file for results.\n");

   return 0;
}

int quaddobl_test ( int islaurent )
{
   int fail,n,dim,deg,nbloops,kind;
   char ans;

   printf("\nReading a witness set ...\n");
   if(islaurent == 1)
      fail = read_quaddobl_Laurent_witness_set(&n,&dim,&deg);
   else
      fail = read_quaddobl_witness_set(&n,&dim,&deg);

   printf("\nMENU for the kind of output :\n");
   printf("  0. remain silent with no intermediate output;\n");
   printf("  1. all intermediate output goes to screen;\n");
   printf("  2. give a file name for all intermediate output.\n");
   printf("Type 0, 1, or 2 to make a choice : "); scanf("%d",&kind);

   scanf("%c",&ans); /* skip end of line character */

   if(kind == 0) fail = set_quaddobl_state_to_silent();
   if(kind == 2) fail = define_output_file();

   if(verbose>0)  /* only in verbose mode */
   {
      printf("\nThe ambient dimension : %d.\n",n);
      printf("The dimension of the solution set : %d.\n",dim);
      printf("The degree of the solution set : %d.\n",deg);
      printf("\nDo you wish to see the embedded system ? (y/n) ");
      scanf("%c",&ans);
      if(ans == 'y')
      {
         printf("\nThe system read :\n");
         if(islaurent == 1)
            fail = syscon_write_quaddobl_Laurent_system();
         else
            fail = syscon_write_quaddobl_system();
      }
   }
   fail = quaddobl_assign_labels(n,deg);

   printf("\nGive the maximum number of loops : ");
   scanf("%d",&nbloops);

   fail = quaddobl_monodromy_breakup(islaurent,nbloops,n,dim,deg);

   if(kind == 2) printf("See the output file for results.\n");

   return 0;
}
