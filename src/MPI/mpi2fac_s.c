/* The program "mpi2fac_s" is a parallel implementation of the factorization
 * of a pure positive dimensional solution set into irreducible components,
 * using static load distribution in master/slave programming model.
 * The program prompts the user for an embedded system (system, slices, 
 * and witness set) and for the maximal number of monodromy loops.  */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "mpi.h"
#include "../Lib/phcpack.h"
#include "../Lib/syscon.h"
#include "../Lib/solcon.h"
#include "parallel_phcpack.h"
#include "parallel_monodromy.h"

#define v 0  /* verbose flag:
                  0 no output during the computations,
                  1 only one line messages on progress of computations
                  2 display of data, like systems and solutions */

int monodromy_breakup ( int myid, int n, int dim, int deg, int nbloops,
                        int numprocs, double *trackwtime );
/* 
 * DESCRIPTION :
 *   Does at most as many monodromy loops as the value in nbloops,
 *   to factor a set of dimension dim and degree deg in n-space.
 *
 * ON ENTRY :
 *   myid          identification of the node;
 *   n             ambient dimension;
 *   dim           dimension of the solution set;
 *   deg           degree of the solution set;
 *   nbloops       maximal number of allowed loops;
 *   numprocs      number of processors.
 *
 * ON RETURN :
 *   trackwtime    time node myid has been tracking paths;
 *   the return value of the function is either 0 or 1:
 *     0           if nbloops sufficed to factor the solution set,
 *     1           otherwise: nbloops was set too small. */

int main ( int argc, char *argv[] )
{
   int myid,numprocs,n,dim,deg,nbloops,fail;
   double startwtime,trackwtime,wtime,*mytime;
   MPI_Status status;

   adainit();
   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
   MPI_Comm_rank(MPI_COMM_WORLD,&myid);
   srand(time(NULL));   /* for different random gamma constant */
   if(myid == 0)
   {
      mytime = (double*) calloc(numprocs, sizeof(double));
      startwtime = MPI_Wtime();
      fail = read_witness_set(&n,&dim,&deg);
      fail = define_output_file();
      printf("Give the number of loops : "); scanf("%d",&nbloops);
   }
   else
      trackwtime = MPI_Wtime();

   MPI_Barrier(MPI_COMM_WORLD);  /* wait for node 0 */
   dimension_broadcast(myid,&n);
   MPI_Bcast(&nbloops,1,MPI_INT,0,MPI_COMM_WORLD);

   fail = monodromy_breakup(myid,n,dim,deg,nbloops,numprocs,&trackwtime);

   MPI_Barrier(MPI_COMM_WORLD);
   if(myid == 0)
     wtime = MPI_Wtime() - startwtime;
   else
     wtime = trackwtime;
   MPI_Gather(&wtime,1,MPI_DOUBLE,mytime,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   if(myid == 0) print_time(mytime,numprocs);
   MPI_Finalize();
   adafinal();

   return 0;
}

int monodromy_breakup ( int myid, int n, int dim, int deg, int nbloops,
                        int numprocs, double *trackwtime )
{
   int mysolnum,i,done,fail;

   if(v>0) printf("Node %d knows n = %d, #loops = %d.\n",myid,n,nbloops);

   monomials_broadcast(myid,n);
   MPI_Bcast(&deg,1,MPI_INT,0,MPI_COMM_WORLD);
   if(v>0) printf("Node %d went through bcast of degree %d.\n",myid,deg);
   solutions_distribute(myid,deg,n,numprocs,&mysolnum);
   
   MPI_Bcast(&dim,1,MPI_INT,0,MPI_COMM_WORLD);
   fail = initialize_standard_sampler(dim);

   if(myid == 0)
   {
      fail = validate_solutions();                   /* sanity check */
      fail = initialize_standard_monodromy(nbloops,deg,dim);
      /* initialize traces */
      if(v>1) print_solutions(myid);
      fail = store_standard_solutions();
   }

   for(i=1; i<=2; i++)
   {
     if((v>0)&&(myid == 0)) printf("slice %d for linear trace...\n",i);
     trace_slice_broadcast(myid,i);
     gamma_broadcast(myid,n);
     f_track_paths(myid,deg,n,numprocs,mysolnum,trackwtime);
     if(myid == 0)
     {
        fail = store_standard_solutions();
        fail = restore_standard_solutions();
     }
     solutions_distribute(myid,deg,n,numprocs,&mysolnum);
     if(myid == 0) f_swap_slices(myid);
   }

   slices_broadcast(myid,dim,n);

   for(i=0; i<nbloops; i++)
   {
      if((v>0)&&(myid == 0)) printf("Starting loop #\%d...\n",i);
      gamma_broadcast(myid,n);
      f_track_paths(myid,deg,n,numprocs,mysolnum,trackwtime);
      if(myid != 0) f_swap_slices(myid);
      if(myid == 0)
      {
         if(v>1) print_solutions(myid);
         solcon_clear_standard_solutions();
      }
      gamma_broadcast(myid,n);
      f_track_paths(myid,deg,n,numprocs,mysolnum,trackwtime);
      if(myid == 0)
      {
         if(v>1) print_solutions(myid);
         fail = store_standard_solutions();
         fail = standard_monodromy_permutation(deg,&done);
         solcon_clear_standard_solutions();
      }

      MPI_Bcast(&done,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);   /* every node must know if done */

      if(done == 1)
      {
         if (v>0) printf("Node %d leaves the loop at step %d.\n", myid,i+1);
         if(myid == 0) printf("Executed %d monodromy loops.\n", i+1);
         return 0;
      }
      else
      {
         if(myid != 0) f_swap_slices(myid);
         slices_broadcast(myid,dim,n);               /* new slices */
         if(myid == 0) fail = restore_standard_solutions();
         MPI_Bcast(&deg,1,MPI_INT,0,MPI_COMM_WORLD);  
         solutions_distribute(myid,deg,n,numprocs,&mysolnum);
      }
   }
   if(myid == 0) printf("Executed %d monodromy loops.\n", i+1);
   return 1;
}
