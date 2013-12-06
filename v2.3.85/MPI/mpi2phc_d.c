/* The following is a simple interface of MPI with PHCpack:
 * after prompting the user for a target and start system,
 * the homotopy is then broadcasted to all the nodes and
 * the paths are distributed dynamically over the nodes. */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "mpi.h"
#include "../Lib/phcpack.h"
#include "../Lib/syscon.h"
#include "../Lib/solcon.h"
#include "parallel_phcpack.h"

#define SEND_SSOL 100 /* tag for sending start solution */
#define SEND_SMUL 101 /* tag for sending multiplicity of start solution */
#define SEND_TSOL 102 /* tag for sending target solution */
#define SEND_TMUL 103 /* tag for sending multiplicity of target solution */

#define v 0  /* verbose flag */

extern void adainit( void );
extern void adafinal( void );

int track_one_path ( int n, int *m, double *sol );
/*
 * DESCRIPTION :
 *   Tracks one path starting at the given start solution.
 *
 * ON ENTRY :
 *   n        dimension of the solution;
 *   m        label of the solution, multiplicity flag;
 *   sol      coordinates of the solution vector. */

void dynamic_load ( int myid, int n, int numprocs, int nbsols, int *mysolnb );
/*
 * DESCRIPTION :
 *   The manager distributes start solutions to the node for tracking.
 *   After path tracking, the nodes return their solutions to the manager.
 *   The manager writes the solution to file and distributes a new jobs. */ 

int copy_broadcast ( int myid );
/*
 * DESCRIPTION :
 *   All worker nodes copy the result of the broadcast from the
 *   system container as the target system in PHCpack.
 *   This routine offers a first test on the communication network. */

void print_monomials ( void );
/*
 * DESCRIPTION :
 *   Writes the monomials in the container to screen.
 *   Only used for testing purposes. */

int print_homotopy ( void );
/* 
 * DESCRIPTION :
 *   Writes the target and start system with all start solutions
 *   to the defined output file. */

void print_time_summary ( double *time, int numprocs, int *mysolnb );
/*
 * DESCRIPTION :
 *   Prints the total wall time spent and the number of paths tracked
 *   on every processor. */

int main ( int argc, char *argv[] )
{
   int myid,numprocs,n,nbsols,mysolnb,*mysol,fail;
   double startwtime,endwtime,wtime,*time,vcp[34];
   MPI_Status status;
 
   adainit();
   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
   MPI_Comm_rank(MPI_COMM_WORLD,&myid);

   if(myid == 0)
   {
      time = (double*)calloc(numprocs,sizeof(double));
      mysol = (int*)calloc(numprocs,sizeof(int));
      startwtime = MPI_Wtime();
      fail = read_target_system();
      fail = copy_target_system_to_container();
      fail = syscon_number_of_polynomials(&n);
   }

   dimension_broadcast(myid,&n);
   monomials_broadcast(myid,n);
   fail = copy_broadcast(myid); 
   fail = start_system_broadcast(myid,n,&nbsols);

   if(myid == 0)
   {
      fail = define_output_file(); printf("\n");
      fail = print_homotopy();
      fail = tune_continuation_parameters(); printf("\n");
     /* fail = determine_output_during_continuation(); */
      fail = retrieve_continuation_parameters(vcp);
      printf("\nSee the output file for results...\n\n");
   }

   MPI_Bcast(&nbsols,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(vcp,34,MPI_DOUBLE,0,MPI_COMM_WORLD);
   
   if(myid != 0) 
   {
      fail = set_continuation_parameters(vcp);
      startwtime = MPI_Wtime();
   }
   dynamic_load(myid,n,numprocs,nbsols,&mysolnb);
   if(myid != 0) endwtime = MPI_Wtime();

   if(myid == 0)
   {
      fail = solcon_write_solutions();
      fail = solcon_clear_solutions();
   }

   MPI_Barrier(MPI_COMM_WORLD);
   if(myid == 0) endwtime = MPI_Wtime();
   wtime = endwtime-startwtime;
   MPI_Gather(&wtime,1,MPI_DOUBLE,time,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Gather(&mysolnb,1,MPI_INT,mysol,1,MPI_INT,0,MPI_COMM_WORLD);
   if(myid == 0) print_time_summary(time,numprocs,mysol);
   MPI_Finalize();
   adafinal();

   return 0;
}

void dynamic_load ( int myid, int n, int numprocs, int nbsols, int *mysolnb )
{
   int dest, m[2], send[2], i, k, fail;
   double sol[2*n+5];
   MPI_Status status;

   m[0] = n;
   *mysolnb = 0;
   if(myid == 0)
   {
      /* distributes one start root for each client at the beginning */
      for(dest=1; dest<=numprocs-1; dest++)
      {
         fail = solcon_retrieve_solution(n,dest,&m[1],sol);
         send[0] = dest; send[1] = m[1];
         MPI_Send(send,2,MPI_INT,dest,SEND_SMUL,MPI_COMM_WORLD);
         MPI_Send(sol,2*n+5,MPI_DOUBLE,dest,SEND_SSOL,MPI_COMM_WORLD);  
         (*mysolnb)++;
      }
      /* collects target roots and distributes remaining start roots */ 
      for(k=numprocs; k<=nbsols+numprocs-1; k++)
      {
         MPI_Recv(&m[1],1,MPI_INT,MPI_ANY_SOURCE,SEND_TMUL,
                  MPI_COMM_WORLD,&status);
         MPI_Recv(sol,2*n+5,MPI_DOUBLE,MPI_ANY_SOURCE,SEND_TSOL,
                  MPI_COMM_WORLD,&status);
         fail = solcon_replace_solution(n,k-numprocs+1,m[1],sol);
         fail = solcon_retrieve_solution(n,k,&m[1],sol);
         send[0] = k; send[1] = m[1];
         dest = status.MPI_SOURCE;
         MPI_Send(send,2,MPI_INT,dest,SEND_SMUL,MPI_COMM_WORLD);
         if(k<=nbsols)
         {
            MPI_Send(sol,2*n+5,MPI_DOUBLE,dest,SEND_SSOL,MPI_COMM_WORLD);
            (*mysolnb)++;
         }
      }
   }
   else
   {
      while(1)
      {
         MPI_Recv(send,2,MPI_INT,0,SEND_SMUL,MPI_COMM_WORLD,&status);
         if (send[0]>nbsols) break;
         m[1] = send[1];
         MPI_Recv(sol,2*n+5,MPI_DOUBLE,0,SEND_SSOL,MPI_COMM_WORLD,&status); 
         fail = track_one_path(m[0],&m[1],sol);
         MPI_Send(&m[1],1,MPI_INT,0,SEND_TMUL,MPI_COMM_WORLD);
         MPI_Send(sol,2*n+5,MPI_DOUBLE,0,SEND_TSOL,MPI_COMM_WORLD);
         (*mysolnb)++;
      }
      if(v>0) printf("Node %d done.\n", myid);
   }
}

int track_one_path ( int n, int *m, double *sol )
{
   int fail;

   fail = solcon_append_solution(n,*m,sol);
   fail = copy_container_to_start_solutions();
   fail = solcon_clear_solutions();
   fail = solve_by_homotopy_continuation();
   fail = copy_target_solutions_to_container();
   fail = solcon_retrieve_solution(n,1,m,sol);
   fail = solcon_clear_solutions();

   return fail;
}

int copy_broadcast ( int myid )
{
   int fail = 0;

/* the code below is a sanity check on the communication network :
 
   if(myid == 1) 
   {
      printf("The monomials at node %d:\n", myid);
      print_monomials();
      printf("The polynomial system at node %d:\n", myid);
      fail = print_system();
   }

   end of sanity check on the communication network. */

   if(myid != 0)
   {
      fail = copy_container_to_target_system();
      fail = syscon_clear_system();
   }

   return fail;
}

void print_monomials ( void ) 
{
   int *d,i,j,k,n,mm,fail;
   double c[2];

   fail = syscon_number_of_polynomials(&n);

   d = (int*)calloc(n,sizeof(int));

   for(i=1; i<=n; i++)
   {
      fail = syscon_number_of_terms(i,&mm);
      printf("Polynomial %d has %d monomials :\n",i,mm);
      for(j=1; j<=mm; j++)
      {
         fail = syscon_retrieve_term(i,j,n,d,c);
         printf(" %.15e  %.15e",c[0],c[1]);
         for(k=0; k<n; k++) printf(" %d",d[k]);
         printf("\n");
      }
   }
}

int print_homotopy ( void )
{
   int fail;

   fail = write_target_system();
   fail = write_start_system();
   fail = write_start_solutions();

   return fail;
}

void print_time_summary ( double *time, int numprocs, int *mysolnb )
{
   int i;
   printf("\nTotal wall time = %lf seconds on %d processors\n",
          time[0],numprocs);
   printf("Total number of paths tracked: %d.\n", mysolnb[0]);

   for(i=1; i<numprocs; i++)
   {
      printf("The wall time spent on NO. %d processor = %lf seconds.\n",
             i,time[i]);
      printf("The number of paths tracked on NO. %d processor: %d.\n",
             i,mysolnb[i]);
   }
   free(time);
   free(mysolnb);
}
