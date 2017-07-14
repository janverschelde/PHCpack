/* The program "monodromyMPI" is a parallel implementation of the
 * factorization of a pure positive dimensional solution set
 * into irreducible components.  The program prompts the user
 * for an embedded system (system, slices, and witness set)
 * and the dimension of the solution set.  */

#define v 2  /* verbose flag:
                  0 no output during the computations,
                  1 only one line messages on progress of computations
                  2 display of data, like systems and solutions 
		  3 low-level debugging messages 
	     */


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

#define SEND_SSOL 100 /* tag for sending start solution */
#define SEND_SMUL 101 /* tag for sending multiplicity of start solution */
#define SEND_TSOL 102 /* tag for sending target solution */
#define SEND_TMUL 103 /* tag for sending multiplicity of target solution */

void start_walk(int sliceA, int start_label, int sliceB, SLAVE_STATE* state);
/*
 * DESCRIPTION :
 *   This is the main job scheduling function. */

void run_master(int n, int s, int max_s, int np);
/*
 * DESCRIPTION :
 *   This function is executed by the root node 0. */

void run_slave(int n);
/*
 * DESCRIPTION :
 *   This function is executed by the computational nodes. */

int main( int argc, char *argv[] )
{
  int fail,myid,numprocs,
    n/*ambient dim*/, 
    dim/*dim(solution set)*/, 
    nb_slices,
    max_nb_slices, 
    deg/*deg(witness set)*/;
  double startwtime,trackwtime,wtime,*mytime;
  MPI_Status status;

  adainit();
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  srand(time(NULL));   /* for different random gamma constant */
  if(myid == 0) {
    mytime = (double*) calloc(numprocs, sizeof(double));
    startwtime = MPI_Wtime(); 
    
    fail = read_witness_set(&n,&dim,&deg);
    fail = define_output_file();

    printf("Give the number of slices: "); scanf("%d",&nb_slices);
    printf("Give the maximal number of slices: "); scanf("%d",&max_nb_slices);
  } else {
    trackwtime = MPI_Wtime(); 
    if (v>=3) printf("\nLaunching slave %d\n", myid);
  }
  
  MPI_Barrier(MPI_COMM_WORLD);  /* everybody must wait for node 0 */
  dimension_broadcast(myid,&n);
  MPI_Bcast(&nb_slices,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&max_nb_slices,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&dim,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&deg,1,MPI_INT,0,MPI_COMM_WORLD);
  if (v>=1) printf("Node %d: went through bcast of degree,dim,#slices.\n",myid);

  fail = initialize_standard_monodromy(max_nb_slices,deg,dim);
  fail = build_trace_grid(myid,n,dim,deg,numprocs,&trackwtime);

  MPI_Barrier(MPI_COMM_WORLD);  /* wait for finish build_trace_grid */
  trace_grid_broadcast(myid,deg,n);
  all_slices_broadcast(myid,max_nb_slices,dim,n);
  if (v>=1) printf("Node %d: construction and broadcast of trace grid and slices... done.\n",myid);
  MPI_Barrier(MPI_COMM_WORLD);  /* everybody waits for trace grid */
    
  fail = track_paths_to_new_slices(myid,nb_slices,deg,dim,n,numprocs,&trackwtime);

  MPI_Barrier(MPI_COMM_WORLD);


  if (myid ==0) run_master(deg, nb_slices+1, max_nb_slices+1, numprocs);
  else run_slave(deg);
   
  MPI_Barrier(MPI_COMM_WORLD);
  
  if(myid == 0)
    wtime = MPI_Wtime() - startwtime;
  else
    wtime = trackwtime;
  MPI_Gather(&wtime,1,MPI_DOUBLE,mytime,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  if(myid == 0) {
    print_time(mytime,numprocs);
  }
   
  MPI_Finalize();  
  adafinal();
  return 0;
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

void run_master(int n, int s, int max_s, int np)
{
  int i,j,sl,t,c;
  POINT pt;
  
  int n_paths = n*(s-1);
  int n_rejected = 0;
  ISLICE_STATS* stats = (ISLICE_STATS*) calloc(max_s*max_s,sizeof(ISLICE_STATS));
  SLAVE_STATE* state;
  int* ind;  /* ind[i] = the number of the component that i-th point belongs to */
  int n_chunks;
  COMPONENT* part;

  MPI_Status status;
  int flag, ept, spt; 

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
  } 

  /* MAIN LOOP */
#define MAX_PATHS 1000 /*!!!*/
  while (n_chunks>0 
	 && n_paths<MAX_PATHS) {    
    /* make busy the slaves */ 
    while(free_slave(state, np)>0) {
      /* pick a pair of slices such that the rejection rate is minimal */ 
      double tr, rr = 2; int ii,jj;
      
      while (rr>1) {
	for (ii=0; ii<s; ii++)   
	  for (jj=ii+1; jj<s; jj++)  
	    if (isGood(&stats[ii*max_s+jj])
		&& ii>0 // avoid the base slice for now !!!
		) {
	      tr = (stats[ii*max_s+jj].n_paths==0) ? 0: 
		(double) stats[ii*max_s+jj].n_rejected / stats[ii*max_s+jj].n_paths;
	      if (tr<rr) 
		{ rr = tr; i = ii; j = jj; }
	    }
	if (rr>1) {
	  s = s+1;
	  if (s>max_s) {
	    perror("maximum number of slices exceeded");
	    exit(max_s);
	  }
	}
      };
	      
      /* find the smallest reducible component 
	 (!!! we should implement sorting instead!!!)*/
      c = -1;
      for (t=1; t<n; t++){
	if (part[t].size!=0 
	    && !part[t].is_irred 
	    && (c==-1 || part[t].size<part[c].size))
	  c = t;
      }
      start_walk(i, part[c].pts[0], j, &state[free_slave(state, np)]);
      n_paths++; stats[i*max_s+j].n_paths++;
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

      /* if the endpoint was not labelled before...
      if (state[sl].status == NEW_POINT) {
	if (ept != spt) { // end point label should equal the start label
	  perror("new label got matched with an old point");
	  exit(NEW_POINT);
	}
	} else {*/
      if (v>=2) printf("M: there is a path from %d to %d\n", spt, ept);
      if (ept>0 && ept<=n && ind[ept] != ind[spt]) { /* if not in the same component... */
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
	stats[state[sl].sliceA*max_s+state[sl].sliceB].n_rejected++;
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
Summary: the total of %d paths was tracked, %d rejected.\n ", 
	 n_paths, n_rejected);
  for (i=0; i<s; i++) for (j=i+1; j<s; j++)
    printf("(%d)<->(%d): the total of %d paths was tracked, %d rejected.\n ", 
	   i,j,stats[i*max_s+j].n_paths, stats[i*max_s+j].n_rejected);  

  /* CLEAN-UP */
  free(ind);
  for(i=0; i<n; i++){
    killCOMPONENT(&part[i]); 
  } 
  free(part);
  free(stats);
}

void run_slave(int n) 
{
  MPI_Status status;
  SLAVE_STATE st;
  int target_label,fail,flag = 0;
  int sl; 

  while(!flag)   /* while there is no "quit" message */
  {
    MPI_Iprobe(MPI_ANY_SOURCE, SEND_TRACK_JOB, MPI_COMM_WORLD, &flag, &status);
    if (flag) {
      MPI_Recv(&st, sizeof(SLAVE_STATE), MPI_CHAR, status.MPI_SOURCE,
	       SEND_TRACK_JOB, MPI_COMM_WORLD, &status);
      if (v>=3) printf("S_%d: tracking label %d from slice %d to slice %d\n", 
	     st.num, st.start_label, st.sliceA, st.sliceB);
      
      fail = standard_sample_loop
               (st.sliceA,st.sliceB,st.start_label,&target_label);

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
