/* The main program prompts the user for a witness set and then runs the
   basic version of the monodromy loop algorithm to decompose the witness
   set into irreducible factors. */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "phcpack.h"
#include "solcon.h"
#include "witset.h"

#define v 0 /* verbose flag */

int assign_labels ( int n, int nbsols );
/* 
 * DESCRIPTION :
 *   Assigns a unique label between 1 and nbsols for each solution in the
 *   solutions container, using the multiplicity field of the solution. */

int set_trace_slice ( int first );
/*
 * DESCRIPTION :
 *   Sets the coefficient of the slice used in the linear trace,
 *   when called for the first time, first must have the value 1. */

int store_gammas ( int n );
/*
 * DESCRIPTION :
 *   Generates n random complex coefficients
 *   and stores these coefficients into the monodromy breakup machine. */

int monodromy_breakup ( int nbloops, int n, int k, int d );
/*
 * DESCRIPTION :
 *   Uses nbloops to break up a k-dimensional set of degree d in n-space. */ 

int main ( int argc, char *argv[] )
{
   int fail,n,dim,deg,nbloops,kind;
   char ans;

   adainit();
   srand(time(NULL));

   printf("\nReading a witness set ...\n");
   fail = read_witness_set(&n,&dim,&deg);
   printf("\nMENU for the kind of output :\n");
   printf("  0. remain silent with no intermediate output;\n");
   printf("  1. all intermediate output goes to screen;\n");
   printf("  2. give a file name for all intermediate output.\n");
   printf("Type 0, 1, or 2 to make a choice : "); scanf("%d",&kind);

   scanf("%c",&ans); /* skip end of line character */

   if(kind == 0) fail = set_state_to_silent();
   if(kind == 2) fail = define_output_file();

   if(v>0)  /* only in verbose mode */
   {
      printf("\nThe ambient dimension : %d.\n",n);
      printf("The dimension of the solution set : %d.\n",dim);
      printf("The degree of the solution set : %d.\n",deg);
      printf("\nDo you wish to see the embedded system ? (y/n) ");
      scanf("%c",&ans);
      if(ans == 'y')
      {
         printf("\nThe system read :\n");
         fail = print_system();
      }
   }

   fail = assign_labels(n,deg);

   printf("\nGive the number of loops : ");
   scanf("%d",&nbloops);

   fail = monodromy_breakup(nbloops,n,dim,deg);

   if(kind == 2) printf("See the output file for results.\n");

   adafinal();

   return 0;
}

int assign_labels ( int n, int nbsols )
{
   int i,j,m,fail;
   double x[2*n+5];

   for(i=1; i<=nbsols; i++)
   {
      fail = solcon_retrieve_solution(n,i,&m,x);
      m = i;
      fail = solcon_replace_solution(n,i,m,x);
   }

   return fail;
}

int set_trace_slice ( int first )
{
   int fail;
   double r[2];

   if(first == 1)                  /* determine constant coefficient */
   { 
      r[0] = -1.0; r[1] = 0.0;
   }
   else
   {
      r[0] = +1.0; r[1] = 0.0;
   }
   fail = assign_coefficient_of_slice(0,0,r);

   return fail;
}

int store_gammas ( int n )
{
   double re_gamma[n];
   double im_gamma[n];
   int i;
    
   for(i=0; i<n; i++)
      random_complex(&re_gamma[i],&im_gamma[i]);
   store_gamma(n,re_gamma,im_gamma);
}

int monodromy_breakup ( int nbloops, int n, int k, int d )
{
   int fail,i,j,done;

   fail = initialize_sampler(k);
   fail = initialize_monodromy(nbloops,d,k);
   if(v>0)
      printf("\nInitialized sampler and monodromy permutations ...\n");

   fail = store_solutions();  /* store in monodromy permutations */

   printf("... initializing the grid ...\n");
   for(i=1; i<=2; i++)        /* initialize grid for trace validation */
   {
      fail = set_trace_slice(i);  /* fix constant coefficient of slice */
      fail = store_gammas(n);     /* generate random gamma constants */
      fail = track_paths();
      fail = store_solutions();   /* store solutions in the grid */
      fail = restore_solutions(); /* use original solutions at start */
      fail = swap_slices();       /* go back to original slices */
   }

   done = 0;
   for(i=1;(i<=nbloops) && (done==0); i++)  /* perform monodromy loops */
   {
      printf("... starting loop #\%d ...\n",i);
      fail = new_slices(k,n);
      fail = store_gammas(n);
      fail = track_paths();     /* swapping slices happens here */
      fail = solcon_clear_solutions();
      fail = store_gammas(n);
      fail = track_paths();
      fail = store_solutions();
      {
         int permutation[d];
         fail = permutation_after_loop(d,permutation);
         printf("the permutation at step %d:",i);
         for(j=0; j<d; j++) printf(" %d",permutation[j]);
         printf("\n");
      }
      fail = monodromy_permutation(d,&done);
      {
         int nf,deg,jp,w[d];
         double trace_difference;
         fail = number_of_irreducible_factors(&nf);
         printf("number of irreducible factors : %d\n",nf);
         for(j=1; j<=nf; j++)
         {
            fail = witness_points_of_irreducible_factor(j,&deg,w);
            printf(" %d :",j);
            for(jp=0; jp<deg; jp++) printf(" %d",w[jp]);
            fail = trace_sum_difference(deg,w,&trace_difference);
            printf(" : %.3e\n",trace_difference);
         }
      }
      fail = restore_solutions();
   }
   if(done==1)
     printf("Found factorization using %d loops.\n", i-1);
   else
     printf("Failed to factor using %d loops.\n", i-1);

   return fail;
}
