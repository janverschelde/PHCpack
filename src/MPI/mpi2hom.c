/* The following is a simple interface of MPI with PHCpack:
 * after prompting the user for a target and start system,
 * the homotopy is then broadcasted to all the nodes. */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "mpi.h"
#include "../Lib/syscon.h"
#include "../Lib/solcon.h"
#include "../Lib/phcpack.h"
#include "../Lib/jump_track.h"
#include "parallel_phcpack.h"

#define SEND_SOL 100 /* message tag for sending start solution */
#define SEND_MUL 101 /* message tag for sending multiplicity */

extern void adainit( void );
extern int _ada_use_c2phc4c ( int task, int *a, int *b, double *c );
extern void adafinal( void );

void test_broadcast ( int myid, int start );
/*
 * Writes the target system if start == 0,
 * otherwise writes the start system. */

int main ( int argc, char *argv[] )
{
   int myid,numprocs,n,nbsols,mysolnum,fail,len;
   MPI_Status status; 

   adainit();
   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
   MPI_Comm_rank(MPI_COMM_WORLD,&myid);

   if(myid == 0) 
   {
      fail = read_target_system_without_solutions();
      fail = copy_target_system_to_container();
      fail = syscon_number_of_standard_polynomials(&n);
      printf("The dimension is %d.\n",n);
   }

   dimension_broadcast(myid,&n);
   monomials_broadcast(myid,n);
   test_broadcast(myid,0); 

   MPI_Barrier(MPI_COMM_WORLD);

   fail = start_system_broadcast(myid,n,&nbsols);

   test_broadcast(myid,1);

   MPI_Barrier(MPI_COMM_WORLD);

   MPI_Bcast(&nbsols,1,MPI_INT,0,MPI_COMM_WORLD);
   solutions_distribute(myid,nbsols,n,numprocs,&mysolnum);

   fail = solcon_number_of_standard_solutions(&len);
   printf("Node %d has %d solutions.\n",myid,len);

   MPI_Finalize();
   adafinal();

   return 0;
}

void test_broadcast ( int myid, int start )
{
   int fail,dim;

   fail = syscon_number_of_standard_polynomials(&dim);

   printf("Dimension %d at node %d, the system:\n",dim,myid);
   if(start == 0)
   {
      fail = copy_container_to_target_system();
      fail = write_standard_target_system();
   }
   if(start == 1) fail = write_standard_start_system();
}
