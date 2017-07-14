/* The program "mpi2fac_d" is a parallel implementation of the factorization
 * of a pure positive dimensional solution set into irreducible components,
 * using a dynamic load distribution.
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
#include "manage_components.h"

#define v 1  /* verbose flag:
                  0 no output during the computations,
                  1 only one line messages on progress of computations
                  2 display of data, like systems and solutions */

#define tf 0  /* testing flag:
                  0 no interactive test is necessary;
                  1 some interactive tests will be run */

int initialize_slices ( int myid, int n, int dim, int deg, int nb_slices,
                        int numprocs, double *trackwtime );
/* 
 * DESCRIPTION :
 *   Computes as many witness sets as the value of nb_slices + 2,
 *   to decompose a solution set of degree deg and dimension dim 
 *   in n-space  into irreducible factors.
 *
 * ON ENTRY :
 *   myid          identification of the node;
 *   n             ambient dimension;
 *   dim           dimension of the solution set;
 *   deg           degree of the solution set;
 *   nb_slices     maximal number of slices to be used;
 *   numprocs      number of processors.
 *
 * ON RETURN :
 *   trackwtime    time node myid has been tracking paths;
 *   the return value of the function is either 0 or 1:
 *     0           if the trace grid and all slices are well,
 *     1           otherwise: if the trace grid is not good. */

void start_walk(int sliceA, int start_label, int sliceB, SLAVE_STATE* state);
/*
 * DESCRIPTION :
 *   This is the main job scheduling function. */

void run_master(int n, int s, int max_s, int np, double* breakup_time);
/*
 * DESCRIPTION :
 *   This function is executed by the root node 0. */

void run_slave(int n, double* trackwtime);
/*
 * DESCRIPTION :
 *   This function is executed by the computational nodes. */

void pick_slices_random(int s, int* i, int* j);
/*
 * DESCRIPTION :
 *   Returns two random numbers "i","j" from 0 to "s"-1. */

void pick_slices_by_rejection_rate(int s, int max_s, ISLICE_STATS* stats, int* i, int* j);
/*
 * DESCRIPTION :
 *   Picks a pair of slices with the total lowest rejection rate. */

void pick_slices_by_truncated_rejection_history(int s, int max_s, ISLICE_STATS* stats, 
						int* i, int* j); 
/*
 * DESCRIPTION :
 *   Picks a pair of slices with the lowest rejection rate in the short past. */

int main ( int argc, char *argv[] )
{
   int myid,numprocs,n,dim,deg,nb_slices,fail;
   double startwtime,trackwtime,wtime,*mytime,init_time,breakup_time;
   MPI_Status status;

   adainit();
   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
   MPI_Comm_rank(MPI_COMM_WORLD,&myid);
   srand(time(NULL));   /* for different random gamma constant */
   if(myid == 0)
   {
      mytime = (double*) calloc(numprocs, sizeof(double));
      fail = read_witness_set(&n,&dim,&deg);
      fail = define_output_file();
      printf("\nGive the number of slices : "); scanf("%d",&nb_slices);
      startwtime = MPI_Wtime();
   }
   else
      trackwtime = MPI_Wtime();

   MPI_Barrier(MPI_COMM_WORLD);  /* everybody must wait for node 0 */
   dimension_broadcast(myid,&n);
   MPI_Bcast(&nb_slices,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&dim,1,MPI_INT,0,MPI_COMM_WORLD); 
   MPI_Bcast(&deg,1,MPI_INT,0,MPI_COMM_WORLD);
   
   /* initialize slices */
   init_time = MPI_Wtime();
   fail = initialize_slices(myid,n,dim,deg,nb_slices,numprocs,&trackwtime);
   init_time = MPI_Wtime()-init_time;
   
   MPI_Barrier(MPI_COMM_WORLD);
   if(myid == 0)
      run_master(deg,nb_slices+1,nb_slices+1,numprocs, &breakup_time);
   else {
     trackwtime = 0; //counting slaves time only after this point 
     run_slave(deg,&trackwtime);
   }
   MPI_Barrier(MPI_COMM_WORLD);
   if(myid == 0)
     wtime = MPI_Wtime() - startwtime;
   else
     wtime = trackwtime;
   MPI_Gather(&wtime,1,MPI_DOUBLE,mytime,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   if(myid == 0) {
     printf("Initialization time = %f seconds.\n", init_time);
     printf("Breakup time = %f seconds.\n", breakup_time);     
     print_time(mytime,numprocs);
   }
   MPI_Finalize();
   adafinal();

   return 0;
}

int initialize_slices ( int myid, int n, int dim, int deg, int nb_slices,
                        int numprocs, double *trackwtime )
{
   int fail;

   fail = initialize_standard_monodromy(nb_slices,deg,dim);
   fail = witness_set_distribute(myid,n,dim,deg,numprocs);
   MPI_Barrier(MPI_COMM_WORLD);  /* initialization of sampler done */
   all_slices_broadcast(myid,nb_slices,dim,n);

   fail = build_trace_grid(myid,n,dim,deg,numprocs,trackwtime);

   MPI_Barrier(MPI_COMM_WORLD);  /* wait for finish build_trace_grid */

   trace_grid_broadcast(myid,deg,n);
   MPI_Barrier(MPI_COMM_WORLD);  /* everybody waits for trace grid */
   if(tf > 0)
   {
      if(myid == 0) printf("\nTesting linear trace operations ...\n");
      interactive_trace_test(myid);
   }

   fail = track_paths_to_new_slices
            (myid,nb_slices,deg,dim,n,numprocs,trackwtime);

   MPI_Barrier(MPI_COMM_WORLD);
   if(tf > 0)
   {
      if(myid == 0) printf("\nRunning the in_slice_test ...\n");
      interactive_in_slice_test(myid);
   }

   MPI_Barrier(MPI_COMM_WORLD);

   if(tf > 0)
   {
      if(myid == 0) printf("\nTesting sample_loop interactively ...\n");
      interactive_track_to_make_loops(myid);
   }

   return fail;
}

void start_walk (int sliceA, int start_label, int sliceB, SLAVE_STATE* state) 
{ 
   /* change state */
   state->status = TRACKING;
   //(is_in_slice(start_label, sliceB)) ? TRACKING : NEW_POINT;
   state->sliceA = sliceA;
   state->start_label = start_label;
   state->sliceB = sliceB;        

   /* send the job to slave */
    if (v>=3) printf("M: sending a job to slave %d\n", state->num);
   MPI_Send(state, sizeof(SLAVE_STATE), MPI_CHAR, 
            state->num, /* the number of the slave */ 
 	    SEND_TRACK_JOB, MPI_COMM_WORLD);
}

void run_master(int n, int s, int max_s, int np, double* breakup_time)
{
  int i,j,k,sl,t,c,c_nodes,pair_index;
  double endwtime, startwtime;

  POINT pt;
  
  int n_paths = 0;
  int n_rejected = 0;

  ISLICE_STATS* stats = (ISLICE_STATS*)calloc(max_s*max_s,sizeof(ISLICE_STATS));
  SLAVE_STATE* state;
  int* ind; 

  /* ind[i] = the number of the component that i-th point belongs to */
  int n_chunks;
  COMPONENT* part;

  MPI_Status status;
  int flag, ept, spt; 

  startwtime = MPI_Wtime();
  if (v>1) printf("runmaster: n=%d, s=%d, max_s=%d, np=%d\n", n,s,max_s,np); 

  for (i=0; i<max_s; i++) for (j=i+1; j<max_s; j++)
    stats[i*max_s+j].n_rejected = stats[i*max_s+j].n_paths = 0;

  /* records what slave is doing */
  state = (SLAVE_STATE*) calloc(np, sizeof(SLAVE_STATE)); 
  for (sl=1; sl<np; sl++) {
    state[sl].num = sl;
    state[sl].status = FREE; 
  }

  ind = (int*)calloc(n+1,sizeof(int)); 
  
  n_chunks = n;
  part = (COMPONENT*)calloc(n,sizeof(COMPONENT)); 
  /* partition: list of lists */
  
  for(i=0; i<n; i++){
    ind[i+1] = i; /* each point is a singleton in the beginning, 
		   point labels start from 1 */
    makeCOMPONENT(&part[i], i, 1); 
    part[i].pts[0] = i+1;
    if (is_irreducible(&part[i])) {
      part[i].is_irred = 1;
      n_chunks--;
    }
    part[i].n_nodes_tracking = 0;
  } 

  /* MAIN LOOP */
#define MAX_PATHS 10000 /*!!!*/
  while (n_chunks>0 
	 && n_paths<MAX_PATHS) {    
    /* make busy the slaves */ 
    while(free_slave(state, np)>0) {
      /* pick a pair of slices */ 
      //pick_slices_random(s,&i,&j);
      //pick_slices_by_rejection_rate(s, max_s, stats,&i,&j);
      pick_slices_by_truncated_rejection_history(s, max_s, stats,&i,&j);
	      
      /* find the smallest reducible component 
	 with the smallest number of nodes assigned to it
	 (!!! we should implement sorting instead!!!)*/
      c = -1; c_nodes = 0;
      for (t=1; t<n; t++){
	if (part[t].size!=0 
	    && !part[t].is_irred 
	    && (c==-1 
		|| part[t].size<part[c].size 
		|| (part[t].size==part[c].size 
		    && part[c].n_nodes_tracking < c_nodes)
		)) 
	  {
	    c_nodes = part[c].n_nodes_tracking;
	    c = t;
	  }
      }
      start_walk(i, part[c].pts[0], j, &state[free_slave(state, np)]);
      part[c].n_nodes_tracking++; //one more slave is assigned to this component
      n_paths++; 
    }
    
    /* collect results */
    MPI_Iprobe(MPI_ANY_SOURCE, SEND_END_POINT, MPI_COMM_WORLD, &flag, &status); 
    while (flag) {
      sl = status.MPI_SOURCE; /* the number of slave who sent this */
      if (v>=3) printf("M: received a message from slave %d\n", sl);
      MPI_Recv(&ept, /* end point */ 
	       1, MPI_INT, sl,
	       SEND_END_POINT, MPI_COMM_WORLD, &status);
      spt = state[sl].start_label; /* start point */
      part[ind[spt]].n_nodes_tracking--; // one fewer slave is assigned to this component
      pair_index = state[sl].sliceA*max_s+state[sl].sliceB; /* index for "stats"[][] */
      
      if (v>=2) printf("M: there is a path from %d to %d\n", spt, ept);
      if (ept>0 && ept<=n 
	  && ind[ept] != ind[spt] /* if not in the same component... */ 
	  && !part[ind[ept]].is_irred) 
	{ 
	  record_path_outcome(&stats[pair_index], PATH_LED_TO_MERGER); 
	  mergeCOMPONENT(&part[ind[ept]],&part[ind[spt]],ind);
	  n_chunks--;
	  if (is_irreducible(&part[ind[ept]])) {
	    part[ind[ept]].is_irred = 1;
	    n_chunks--;
	  }
	  if (v>=3) {
	    printf("#chunks = %d, index: ", n_chunks);
	    for (i=1; i<=n; i++)
	      printf("[%d]", ind[i]);
	    printf("\n");
	  }
	} else {
	  n_rejected++;
	  if (ept<=0 || ept>n) t=PATH_DIVERGED; 
	  else if (ind[ept] != ind[spt] && part[ind[ept]].is_irred) t=PATH_CROSSED; 
	  else t=PATH_STAYED_WITHIN_COMPONENT;
	  record_path_outcome(&stats[pair_index],t); 
	}
      if (v>=2) print_partition(n, part); /* print the updated partition */
      state[sl].status = FREE;
    
    MPI_Iprobe(MPI_ANY_SOURCE, SEND_END_POINT, MPI_COMM_WORLD, &flag, &status); 
  }// end "collect results"
}

  /* TERMINATE SLAVES */
  for (sl=1; sl<np; sl++)
    MPI_Send(&sl, 1, MPI_INT, 
	     sl, /* the number of the slave */ 
	     SEND_QUIT, MPI_COMM_WORLD);  

  /* PRINT STATISTICS */
  printf(
"*** Statistics ***\n\
Summary: \n\
paths tracked\n\
  %d -- for trace grid,\n\
  %d -- for slices initialization;\n\
after that \n\
  %d more paths tracked,\n\
  %d rejected.\n", 2*n, n*(s-1), n_paths, n_rejected);

  printf("For pairs of slices:\n");
  for (i=0; i<s; i++) for (j=i+1; j<s; j++){
    int total;
    printf("  (%d)<->(%d): \n\
    history [", i,j);
    total = MIN(stats[i*max_s+j].n_paths, N_PATHS_TO_REMEMBER);
    for(k=0; k<total; k++){
      switch(stats[i*max_s+j].path_outcome[k]) {
      case PATH_LED_TO_MERGER: printf("M"); break;
      case PATH_DIVERGED: printf("d"); break;
      case PATH_STAYED_WITHIN_COMPONENT: printf("c"); break;
      case PATH_CROSSED: printf("x"); break;
      };
    }
    printf("]\n\
    %d paths tracked,\n\
    %d rejected.\n", 
	   stats[i*max_s+j].n_paths, stats[i*max_s+j].n_rejected);
  }  
  print_partition(n,part);

  /* CLEAN-UP */
  free(ind);
  for(i=0; i<n; i++){
    killCOMPONENT(&part[i]); 
  } 
  free(part);
  free(stats);

  endwtime = MPI_Wtime();
  *breakup_time += (endwtime - startwtime);
}//end{run_master}

void run_slave(int n, double* trackwtime) 
{
  MPI_Status status;
  SLAVE_STATE st;
  int target_label,fail,flag = 0;
  int sl; 
  double endwtime, startwtime;
  while(!flag)   /* while there is no "quit" message */
  {
    MPI_Iprobe(MPI_ANY_SOURCE, SEND_TRACK_JOB, MPI_COMM_WORLD, &flag, &status);
    if (flag) {
      MPI_Recv(&st, sizeof(SLAVE_STATE), MPI_CHAR, status.MPI_SOURCE,
	       SEND_TRACK_JOB, MPI_COMM_WORLD, &status);
      if (v>=3) printf("S_%d: tracking label %d from slice %d to slice %d\n", 
	     st.num, st.start_label, st.sliceA, st.sliceB);

      startwtime = MPI_Wtime();
      fail = standard_sample_loop
               (st.sliceA,st.sliceB,st.start_label,&target_label);
      endwtime = MPI_Wtime();
      *trackwtime += (endwtime - startwtime);
      
      // end_point = (rand() % n) + 1; // simulate the discovery of the endpoint

      if (v>=3) printf("S_%d: done tracking, returning label %d\n", st.num, target_label);
      MPI_Send(&target_label, 1, MPI_INT, 
	   0, /*status.MPI_SOURCE, /* send to where the job was received from */ 
	   SEND_END_POINT, MPI_COMM_WORLD);
    }
    MPI_Iprobe(MPI_ANY_SOURCE, SEND_QUIT, MPI_COMM_WORLD, &flag, &status);
  }
  
  /* slave receives its number and quits */ 
  MPI_Recv(&sl, 1, MPI_INT, status.MPI_SOURCE,
	   SEND_QUIT, MPI_COMM_WORLD, &status);
  if (v>=3) printf("S_%d: finita la comedia\n", sl);
}

void pick_slices_random(int s, int* i, int* j)
{
  int a,b;
  a = rand()%s;
  b = (a + rand()%(s-1)+1)%s;
  if (a<b) {*i=a; *j=b;} else {*i=b; *j=a;};  
}

void pick_slices_by_rejection_rate(int s, int max_s, ISLICE_STATS* stats, int* i, int* j) 
{
  int ii,jj; double rr = 2, tr;
  while (rr>1) {
    for (ii=0; ii<s; ii++)   
      for (jj=ii+1; jj<s; jj++)  
	if (isGood(&stats[ii*max_s+jj])
	    // && ii>0 // avoid the base slice for now !!!
	    ) {
	  tr = (stats[ii*max_s+jj].n_paths==0) ? 0: 
	    (double) stats[ii*max_s+jj].n_rejected / stats[ii*max_s+jj].n_paths;
	  if (tr<rr) 
	    { rr = tr; *i = ii; *j = jj; }
	}
  };
}

void pick_slices_by_truncated_rejection_history(int s, int max_s, ISLICE_STATS* stats, 
						int* i, int* j) 
{
  int ii,jj,t,total,c;
  double rate,top = -1;
  for (ii=0; ii<s; ii++)   
    for (jj=ii+1; jj<s; jj++)  
      {	  
	c = 0;
	total = MIN(stats[ii*max_s+jj].n_paths, N_PATHS_TO_REMEMBER);
	for (t=0; t<total; t++) 
	  if (stats[ii*max_s+jj].path_outcome[t] == PATH_LED_TO_MERGER)
	    c++; // count the number of successful paths

	rate = (2*total<N_PATHS_TO_REMEMBER) ? 1 : (double)c/total;
	/* there should be at least a half of history record, 
	   before a heuristic choice is made */ 

	if (rate>top) 
	  { top = rate; *i = ii; *j = jj; }
      }
}
