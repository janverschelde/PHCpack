/* The following is a simple interface of MPI with PHCpack:
 * after prompting the user for a target and start system,
 * the homotopy is then broadcasted to all the nodes and
 * the paths are distributed dynamically over the nodes. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "mpi.h"
#include "../Lib/phcpack.h"
#include "../Lib/syscon.h"
#include "../Lib/solcon.h"
#include "../Lib/jump_track.h"
#include "../Lib/witset.h"
#include "parallel_phcpack.h"

#define SEND_SSOL 100 /* tag for sending start solution */
#define SEND_SMUL 101 /* tag for sending multiplicity of start solution */
#define SEND_TSOL 102 /* tag for sending target solution */
#define SEND_TMUL 103 /* tag for sending multiplicity of target solution */

#define v 0  /* verbose flag */

extern void adainit ( void );
extern void adafinal ( void );

void ask_menu_selection ( int *kind );
/*
 * DESCRIPTION :
 *   Presents the user with a menu asking about the type of start system.
 *
 * ON RETURN :
 *   kind     represents the type of start system:
 *             = 1 for total degree start system,
 *             = 2 for linear-product start system,
 *             = 3 for cheater's homotopy start system,
 *             = 4 for cascade homotopy to go one level down,
 *             = 5 for starting an extrinsic diagonal homotopy. */

void read_dimension_of_system ( int nc, char name[nc], int *n );
/*
 * DESCRIPTION :
 *   Opens the file with the given name and reads the first number
 *   which is returned in n.
 *
 * ON ENTRY :
 *   nc       number of characters in the name;
 *   name     name of an input file which contains a system.
 *
 * ON RETURN :
 *   n        the dimension of the system on the file. */

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

int track_one_path ( int n, int i, int *m, double *sol );
/*
 * DESCRIPTION :
 *   Tracks one path starting at the given start solution.
 *
 * ON ENTRY :
 *   n        dimension of the solution;
 *   i        number of the solution;
 *   m        label of the solution, multiplicity flag;
 *   sol      coordinates of the solution vector. */

int get_start ( int n, int kind, int k, int *m, double sol[2*n+6] );
/*
 * DESCRIPTION :
 *   Returns in sol the k-th start solution, depending on kind. */

void track_all_paths ( int myid, int n, int p, int kind, int nbp, int *mynbp );
/*
 * DESCRIPTION :
 *   Distribution of path tracking jobs using dynamic load balancing.
 *   The manager distributes start solutions to the workers for tracking.
 *   After tracking a path, the worker sends the solution to the manager.
 *   The manager writes the solution to file and distributes a new job.
 *
 * REQUIRED :
 *   All processors have received the homotopy.
 *
 * ON ENTRY :
 *   myid     id of the processor executing this routine;
 *   n        is the ambient dimension;
 *   p        the number of processors;
 *   kind     type of start system used:
 *            = 1 for total degree start system,
 *            = 2 for linear-product start system,
 *            = 3 for cheater's homotopy start system,
 *            = 4 for cascade one-level-down homotopy;
 *   nbp      total number of solution paths to be tracked.
 *
 * ON RETURN :
 *   mynbp    is number of paths tracked by executing processor. */ 

void distribute_paths ( int n, int p, int kind, int nbp, int *mynbp );
/*
 * DESCRIPTION :
 *   Code executed by the manager, distributing the paths,
 *   called by track_all_paths, using the same parameters. */

void intersect_two_witness_sets
        ( int myid, int n, int p, int nbp, int *mynbp,
          int n1, int n2, int dim1, int dim2, int deg1, int deg2 );
/*
 * DESCRIPTION :
 *   Distributed path tracking to intersect two witness sets.
 *
 * ON ENTRY :
 *   myid     id of the processor executing this routine;
 *   n        is the ambient dimension of the diagonal homotopy;
 *   p        the number of processors;
 *   nbp      total number of solution paths to be tracked;
 *   n1        ambient dimension for the first witness set;
 *   n2        ambient dimension for the second witness set;
 *   dim1      dimension of the solutions represented by the 1st witness set;
 *   dim2      dimension of the solutions represented by the 2nd witness set;
 *   deg1      degree of the 1st witness set, i.e.: #solutions in 1st set;
 *   deg2      degree of the 2nd witness set, i.e.: #solutions in 2nd set.
 *
 * ON RETURN :
 *   mynbp    is number of paths tracked by executing processor. */

void distribute_to_intersect
        ( int n, int p, int nbp, int *mynbp,
          int n1, int n2, int dim1, int dim2, int deg1, int deg2 );
/*
 * DESCRIPTION :
 *   Code executed by the manager for distributing deg1*deg2 paths,
 *   with the same parameters as in intersect_two_witness_sets. */

int copy_broadcast ( int myid );
/*
 * DESCRIPTION :
 *   All worker nodes copy the result of the broadcast from the
 *   system container as the target system in PHCpack.
 *   This routine offers a first test on the communication network. */

int print_homotopy ( int kind );
/* 
 * DESCRIPTION :
 *   Writes the target and start system with all start solutions
 *   to the defined output file. */

int main ( int argc, char *argv[] )
{
   int myid,numprocs,n,nbsols,mysolnb,*mysol,fail,kind,nc;
   int n1,n2,dim1,dim2,deg1,deg2,cd;
   double startwtime,endwtime,wtime,*time,vcp[34];
   char *name,nl;
   MPI_Status status;
 
   adainit();
   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
   MPI_Comm_rank(MPI_COMM_WORLD,&myid);

   if(myid == 0)
   {
      time = (double*)calloc(numprocs,sizeof(double));
      mysol = (int*)calloc(numprocs,sizeof(int));
      name = (char*)calloc(80,sizeof(char));
      startwtime = MPI_Wtime();
      ask_menu_selection(&kind);
      if(kind < 4)
      {
         fail = read_target_system_without_solutions();
         fail = copy_target_system_to_container();
         fail = syscon_number_of_standard_polynomials(&n);
      }
      if(kind < 5)
      {
         printf("\nReading the name of the file for the start system...");
         printf("\nGive a string of characters : "); scanf("%s",name);
         scanf("%c",&nl); /* skip newline symbol for next reading ...*/
         nc = (int) strlen(name);
         if(kind == 4) read_dimension_of_system(nc,name,&n);
      }
      if(kind == 5)
      {
         fail = read_two_witness_sets(&n1,&n2,&dim1,&dim2,&deg1,&deg2,&cd);
         fail = standard_diagonal_homotopy(dim1,dim2);
         n = cd;
         nbsols = deg1*deg2;
      }
   }
   MPI_Bcast(&kind,1,MPI_INT,0,MPI_COMM_WORLD);
   dimension_broadcast(myid,&n);
   if(kind < 4)
   {
      monomials_broadcast(myid,n);
      fail = copy_broadcast(myid); 
   }
   if(kind < 5)
      fail = named_start_system_broadcast(myid,n,&nbsols,kind,nc,name);

   if(myid == 0)
   {
      if(v>0) printf("# paths to track : %d\n",nbsols);
      fail = define_output_file(); printf("\n");
      fail = print_homotopy(kind);
      fail = tune_continuation_parameters(); printf("\n");
     /* fail = determine_output_during_continuation(); */
      fail = retrieve_continuation_parameters(vcp);
      write_solution_banner_to_defined_output_file(nbsols,n);
      printf("\nSee the output file for results...\n\n");
   }

   MPI_Bcast(&nbsols,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(vcp,34,MPI_DOUBLE,0,MPI_COMM_WORLD);
   
   if(myid != 0) 
   {
      fail = set_continuation_parameters(vcp);
      startwtime = MPI_Wtime();
   }
   if(kind == 4)
      fail = create_cascade_homotopy();
   else
   {
      if(kind == 5) fail = homotopy_broadcast(myid,n);
      fail = create_homotopy_with_given_gamma
                (-0.824263733224601,0.566206056211557);
   }
   if(kind < 5) 
      track_all_paths(myid,n,numprocs,kind,nbsols,&mysolnb);
   else
      intersect_two_witness_sets
        (myid,n,numprocs,nbsols,&mysolnb,n1,n2,dim1,dim2,deg1,deg2);

   if(myid != 0) endwtime = MPI_Wtime();

   MPI_Barrier(MPI_COMM_WORLD);
   if(myid == 0) endwtime = MPI_Wtime();
   wtime = endwtime-startwtime;
   MPI_Gather(&wtime,1,MPI_DOUBLE,time,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Gather(&mysolnb,1,MPI_INT,mysol,1,MPI_INT,0,MPI_COMM_WORLD);
   if(myid == 0)
   {
      write_time_and_paths_to_defined_output_file(numprocs,time,mysol);
      free(time); free(mysol);
   }
   MPI_Finalize();
   adafinal();

   return 0;
}

void ask_menu_selection ( int *kind )
{
   char nl;

   printf("\nMENU for type of start system or homotopy :\n");
   printf("  1. start system is based on total degree;\n");
   printf("  2. a linear-product start system will be given;\n");
   printf("  3. start system and start solutions are provided;\n");
   printf("  4. the homotopy is a cascade to go one level down;\n");
   printf("  5. start extrinsic diagonal homotopy to intersect 2 sets.\n");
   printf("Type 1, 2, 3, 4, or 5 to select type of start system : ");
   scanf("%d",kind);

   scanf("%c",&nl); /* skip newline symbol for next reading ...*/
}

void read_dimension_of_system ( int nc, char name[nc], int *n )
{
   FILE *fp;

   *n = 0;
   
   if(v>0) printf("opening file %s\nto read the dimension ...\n",name);

   fp = fopen(name,"r");
   fscanf(fp,"%d",n);

   if(v>0) printf("found the dimension n = %d\n",*n);

   fclose(fp);
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

   fail = extrinsic_top_diagonal_dimension(*n1,*n2,*dim1,*dim2,cd);
   printf("The top dimension of the extrinsic diagonal cascade : %d.\n",*cd);

   return fail;
}

int get_start ( int n, int kind, int k, int *m, double sol[2*n+6] )
{
   int fail;

   if(kind == 1)
      fail = solcon_compute_total_degree_solution(n,k,m,sol);
   else if(kind == 2)
   {
      int kk = k-1;
      fail = solcon_next_linear_product_solution(n,&kk,m,sol);
   }
   else
      fail = solcon_read_next_solution(n,m,sol);

   return fail;
}

void distribute_paths ( int n, int p, int kind, int nbp, int *mynbp )
{
   int dest,m[2],send[2],fail,label_sol;
   int cnt_sol = 0;     /* counts #solutions written to file */
   int cnt_nbp = 0;     /* counts #solutions to start a path */
   int finished = 0;    /* counts #workers who are finished */
   int done = 0;        /* flag to indicate there are no more paths */
   double sol[2*n+7];
   MPI_Status status;

   m[0] = n;
   *mynbp = 0;

   for(dest=1; dest <= p-1; dest++)  /* one path to each worker at start */
   {
      fail = 1;
      while((fail > 0) && (cnt_nbp < nbp))
      {
         cnt_nbp++;
         fail = get_start(n,kind,cnt_nbp,&m[1],sol);
      }

      if(fail > 0)
         break;        /* cnt_nbp == nbp and no paths left to distribute */
      else
      {
         send[0] = dest; send[1] = m[1];
         if(v>0) printf("sending path %d to node %d.\n",cnt_nbp,dest);
         MPI_Send(send,2,MPI_INT,dest,SEND_SMUL,MPI_COMM_WORLD);
         MPI_Send(sol,2*n+5,MPI_DOUBLE,dest,SEND_SSOL,MPI_COMM_WORLD);  
         (*mynbp)++;
      }
   }                             /* collect results and distribute paths */
   while((cnt_sol < *mynbp) && (finished < p-1)) 
   {
      if(v>0)
      {
         printf("cnt_nbp = %d  cnt_sol = %d  mynbp = %d\n",
                 cnt_nbp,cnt_sol,*mynbp);
         printf("receiving path %d...\n",cnt_sol+1);
      }
      MPI_Recv(sol,2*n+7,MPI_DOUBLE,MPI_ANY_SOURCE,SEND_TSOL,
               MPI_COMM_WORLD,&status);
      m[1] = (int) sol[2*n+5];
      label_sol = (int) sol[2*n+6] - 1;
      fail = solcon_write_next_solution_to_defined_output_file
               (&label_sol,n,m[1],sol);
      cnt_sol++;

      dest = status.MPI_SOURCE;
      if(v>0) printf("looking to give to node %d a job...\n",dest);

      fail = 1;
      if(done == 0)
      {
         if(cnt_nbp < nbp)
         {
            cnt_nbp++;
            fail = get_start(n,kind,cnt_nbp,&m[1],sol);
         }
         if(fail > 0) done = 1;
         if((v>0) && (done == 1)) printf("no more paths to distribute\n");
      }

      send[1] = m[1];
      if(fail > 0)  /* there is no path left to track */
      {
         send[0] = nbp+1; send[1] = m[1];
         if(v>0) printf("sending termination %d to node %d\n",send[0],dest);
         MPI_Send(send,2,MPI_INT,dest,SEND_SMUL,MPI_COMM_WORLD);
         finished++;
         if(v>0) printf("number of finished workers : %d\n",finished);
      }
      else
      {
         send[0] = cnt_nbp;
         if(v>0) printf("sending path %d to node %d\n",cnt_nbp,dest);
         MPI_Send(send,2,MPI_INT,dest,SEND_SMUL,MPI_COMM_WORLD);
         MPI_Send(sol,2*n+5,MPI_DOUBLE,dest,SEND_SSOL,MPI_COMM_WORLD);
         (*mynbp)++;
      }
   }
}

void track_all_paths ( int myid, int n, int p, int kind, int nbp, int *mynbp )
{
   int m[2],send[2],i,k,fail;
   int cnt = 0;
   double sol[2*n+7];
   MPI_Status status;

   m[0] = n;
   *mynbp = 0;

   if(myid == 0)
      distribute_paths(n,p,kind,nbp,mynbp);
   else
   {
      while(1)
      {
         MPI_Recv(send,2,MPI_INT,0,SEND_SMUL,MPI_COMM_WORLD,&status);
         if(v>0) 
	 {  
            if(send[0] <= nbp)
               printf("Node %d receives path %d\n",myid,send[0]);
            else
               printf("Node %d receives %d and breaks.\n",myid,send[0]);
         }
         if(send[0] > nbp) break;
         m[1] = send[1];
         MPI_Recv(sol,2*n+5,MPI_DOUBLE,0,SEND_SSOL,MPI_COMM_WORLD,&status); 
         fail = track_one_path(m[0],send[0],&m[1],sol);
         sol[2*n+5] = (double) m[1];
         sol[2*n+6] = (double) send[0];
         MPI_Send(sol,2*n+7,MPI_DOUBLE,0,SEND_TSOL,MPI_COMM_WORLD);
         (*mynbp)++;
      }
      if(v>0) printf("Node %d done.\n", myid);
   }
}

void distribute_to_intersect
        ( int n, int p, int nbp, int *mynbp,
          int n1, int n2, int dim1, int dim2, int deg1, int deg2 )
{
   const int cd = n;
   double sol1[2*n1+5];
   double sol2[2*n2+5];
   double ps[2*cd+5];
   int dest,m[2],send[2],fail,label_sol,i,j;
   int cnt_sol = 0;     /* counts #solutions written to file */
   int cnt_nbp = 0;     /* counts #solutions to start a path */
   int finished = 0;    /* counts #workers who are finished */
   int done = 0;        /* flag to indicate there are no more paths */
   MPI_Status status;

   m[0] = n;
   *mynbp = 0;
   i=0; j=0;

   for(dest=1; dest <= p-1; dest++)  /* one path to each worker at start */
   {
      fail = 1;
      while((fail > 0) && (cnt_nbp < nbp))
      {
         cnt_nbp++;
         if(v > 0)
            fail = get_next_start_product
                     (&i,&j,1,n1,n2,dim1,dim2,deg1,deg2,cd,sol1,sol2,ps);
         else
            fail = get_next_start_product
                     (&i,&j,0,n1,n2,dim1,dim2,deg1,deg2,cd,sol1,sol2,ps);
         m[1] = 1;
      }

      if(fail > 0)
         break;        /* cnt_nbp == nbp and no paths left to distribute */
      else
      {
         send[0] = dest; send[1] = m[1];
         if(v>0) printf("sending path %d to node %d.\n",cnt_nbp,dest);
         MPI_Send(send,2,MPI_INT,dest,SEND_SMUL,MPI_COMM_WORLD);
         MPI_Send(ps,2*n+5,MPI_DOUBLE,dest,SEND_SSOL,MPI_COMM_WORLD);  
         (*mynbp)++;
      }
   }                             /* collect results and distribute paths */
   while((cnt_sol < *mynbp) && (finished < p-1)) 
   {
      if(v>0)
      {
         printf("cnt_nbp = %d  cnt_sol = %d  mynbp = %d\n",
                 cnt_nbp,cnt_sol,*mynbp);
         printf("receiving path %d...\n",cnt_sol+1);
      }
      MPI_Recv(ps,2*n+7,MPI_DOUBLE,MPI_ANY_SOURCE,SEND_TSOL,
               MPI_COMM_WORLD,&status);
      m[1] = (int) ps[2*n+5];
      label_sol = (int) ps[2*n+6] - 1;
      fail = solcon_write_next_solution_to_defined_output_file
               (&label_sol,n,m[1],ps);
      cnt_sol++;

      dest = status.MPI_SOURCE;
      if(v>0) printf("looking to give to node %d a job...\n",dest);

      fail = 1;
      if(done == 0)
      {
         if(cnt_nbp < nbp)
         {
            cnt_nbp++;
            if(v > 0)
               fail = get_next_start_product
                        (&i,&j,1,n1,n2,dim1,dim2,deg1,deg2,cd,sol1,sol2,ps);
            else
               fail = get_next_start_product
                        (&i,&j,0,n1,n2,dim1,dim2,deg1,deg2,cd,sol1,sol2,ps);
            m[1] = 1;
         }
         if(fail > 0) done = 1;
         if((v>0) && (done == 1)) printf("no more paths to distribute\n");
      }

      send[1] = m[1];
      if(fail > 0)  /* there is no path left to track */
      {
         send[0] = nbp+1; send[1] = m[1];
         if(v>0) printf("sending termination %d to node %d\n",send[0],dest);
         MPI_Send(send,2,MPI_INT,dest,SEND_SMUL,MPI_COMM_WORLD);
         finished++;
         if(v>0) printf("number of finished workers : %d\n",finished);
      }
      else
      {
         send[0] = cnt_nbp;
         if(v>0) printf("sending path %d to node %d\n",cnt_nbp,dest);
         MPI_Send(send,2,MPI_INT,dest,SEND_SMUL,MPI_COMM_WORLD);
         MPI_Send(ps,2*n+5,MPI_DOUBLE,dest,SEND_SSOL,MPI_COMM_WORLD);
         (*mynbp)++;
      }
   }
}

void intersect_two_witness_sets
        ( int myid, int n, int p, int nbp, int *mynbp,
          int n1, int n2, int dim1, int dim2, int deg1, int deg2 )
{
   int m[2],send[2],i,k,fail;
   int cnt = 0;
   double sol[2*n+7];
   MPI_Status status;

   m[0] = n;
   *mynbp = 0;

   if(myid == 0)
      distribute_to_intersect(n,p,nbp,mynbp,n1,n2,dim1,dim2,deg1,deg2);
   else
   {
      while(1)
      {
         MPI_Recv(send,2,MPI_INT,0,SEND_SMUL,MPI_COMM_WORLD,&status);
         if(v>0) 
	 {  
            if(send[0] <= nbp)
               printf("Node %d receives path %d\n",myid,send[0]);
            else
               printf("Node %d receives %d and breaks.\n",myid,send[0]);
         }
         if(send[0] > nbp) break;
         m[1] = send[1];
         MPI_Recv(sol,2*n+5,MPI_DOUBLE,0,SEND_SSOL,MPI_COMM_WORLD,&status); 
         fail = track_one_path(m[0],send[0],&m[1],sol);
         sol[2*n+5] = (double) m[1];
         sol[2*n+6] = (double) send[0];
         MPI_Send(sol,2*n+7,MPI_DOUBLE,0,SEND_TSOL,MPI_COMM_WORLD);
         (*mynbp)++;
      }
      if(v>0) printf("Node %d done.\n", myid);
   }
}

int track_one_path ( int n, int i, int *m, double *sol )
{
   int fail,nbstep,nbfail,nbiter,nbsyst;

   fail = silent_path_tracker
            (n,m,sol,&nbstep,&nbfail,&nbiter,&nbsyst);
   if(v>0) printf("%d : #step : %3d #fail : %3d #iter : %3d #syst : %3d\n",
                  i,nbstep,nbfail,nbiter,nbsyst);

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
      fail = syscon_clear_standard_system();
   }

   return fail;
}

int print_homotopy ( int kind )
{
   int fail;

   if(kind != 4) fail = write_standard_target_system();
   fail = write_standard_start_system();
   fail = write_start_solutions();

   return fail;
}
