/* The program "mpi2fac" is a parallel implementation of the
 * factorization of a pure positive dimensional solution set
 * into irreducible components.  The program prompts the user
 * for an embedded system (system, slices, and witness set)
 * and the dimension of the solution set.  */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "mpi.h"

#define SEND_SSOL 100 /* tag for sending start solution */
#define SEND_SMUL 101 /* tag for sending multiplicity of start solution */
#define SEND_TSOL 102 /* tag for sending target solution */
#define SEND_TMUL 103 /* tag for sending multiplicity of target solution */

#define v 0  /* verbose flag:
                  0 no output during the computations,
                  1 only one line messages on progress of computations
                  2 display of data, like systems and solutions */

extern void adainit( void );
extern int _ada_use_c2phc ( int task, int *a, int *b, double *c );
extern void adafinal( void );

void read_embedding ( int myid, int *n, int *dim, int *deg );
/*
 * DESCRIPTION :
 *   The root node read the embedding, after broadcasting the ambient
 *   dimension n, every node intializes its system container with n.
 *
 * ON ENTRY :
 *   myid     number of the node.
 *
 * ON RETURN :
 *   n        ambient dimension, i.e.: total number of variables;
 *   dim      dimension of the solution set to factor;
 *   deg      degree of the solution set.  */

void monomials_broadcast ( int myid, int n );
/*
 * DESCRIPTION :
 *   The system container at the root node is broadcasted to all nodes.
 *
 * ON ENTRY :
 *   myid     number of the node.
 *   n        ambient dimension.  */

void solutions_distribute ( int myid, int nbsols, int n, int nprocs,
                            int *solnum );
/*
 * DESCRIPTION :
 *   Distributes the solutions in the container at the root
 *   to the solutions container at each node.
 *
 * ON ENTRY : 
 *   myid     number of the node;
 *   nbsols   total number of solutions to distribute;
 *   n        length of the solution vectors;
 *   nprocs   the number of nodes;
 *
 * ON RETURN : 
 *   solnum   number of solutions assigned to the node.  */

void initialize_sampler ( int myid, int dim );
/*
 * DESCRIPTION :
 *   Every node (except for the root) initializes the sampling machine
 *   with the system and solutions in their container.
 *
 * ON ENTRY :
 *   myid     number of the node;
 *   dim      dimension of the witness set.  */

void random_complex ( double *re, double *im );
/*
 * DESCRIPTION :
 *   Returns a random complex number on the complex unit circle,
 *   re is the real and im the imaginary part of the random number.  */

void slices_broadcast ( int myid, int k, int n );
/*
 * DESCRIPTION :
 *   The root node generates k random slides in n-space,
 *   which are then broadcasted to all nodes.  */

void trace_slice_broadcast ( int myid, int first );
/* 
 * DESCRIPTION :
 *   Sets the coefficient of the slice used in the linear trace,
 *   for the first time, first must have the value 1. */

void gamma_broadcast ( int myid, int n );
/* 
 * DESCRIPTION :
 *   The root node generates n random complex numbers,
 *   which are broadcasted to all nodes.  */

void initialize_monodromy ( int n, int d, int k );
/* 
 * DESCRIPTION :
 *   Initialize the package Monodromy_Permutations for n loops,
 *   to factor a k-dimensional solution component of degree d.  */

void swap_slices ( int myid );
/* 
 * DESCRIPTION :
 *   Swaps the current slices with new slices and takes new solutions
 *   as start to turn back. */

void solutions_collect ( int myid, int nbsols, int n, 
                         int numprocs, int mysolnum );
/* collects all target solutions in the container at the root
 * from the solutions container at each node */

void broadcast_gamma ( int myid );
/* root generates a random constant which then gets broadcasted
 * and stored in the PHCpack_Operations module of all nodes */

void track_paths ( int myid, int deg, int n, int numprocs, int mysolnum,
                   double *trackwtime );
/* nodes track the paths defined by the homotopy and the start solutions,
 * the root node collects the target solutions */

void monodromy_breakup ( int myid, int n, int dim, int deg, int nbloops,
                         int numprocs, double *trackwtime );
/* does as many monodromy loops as the value in nbloops, to factor a set
 * of dimension dim and degree deg in n-space */

void monodromy_permutation ( int d, int *done );
/* computes the permutation by last stored solution,
 * d is the number of solutions or the degree of the set,
 * if *done == 1 on return, then the linear trace test has certified
 * the current monodromy breakup, otherwise we are not done yet */

void print_monomials ( void );
/* prints the monomials in the container */

void print_system ( void );
/* prints the polynomial system in the container */

void print_solutions ( int myid );
/* prints the solutions in the container */

void refine_solutions ( void );
/* applies a root refiner to the solutions in the container */

void store_solutions ( void );
/* store the solutions in the container to the Monodromy_Permutations */

void restore_solutions ( void );
/* restores the original first list of solutions to the container */

void clear_solutions ( void );
/* clears the solutions in the container */

void print_time ( double *time, int numprocs );
/* prints the time spent on each of the processor */

int main ( int argc, char *argv[] )
{
   int myid,numprocs,n,dim,deg,nbloops;
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
   }
   else
      trackwtime = MPI_Wtime();

   read_embedding(myid,&n,&dim,&deg);
   if(myid == 0) printf("Give the number of loops : "); scanf("%d",&nbloops);
   MPI_Bcast(&nbloops,1,MPI_INT,0,MPI_COMM_WORLD);

   monodromy_breakup(myid,n,dim,deg,nbloops,numprocs,&trackwtime);

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

void read_embedding ( int myid, int *n, int *dim, int *deg )
{
   int fail,m;
   double *c;   
   int d[2];

/* if(v>0) printf("Node %d has entered read_embedding.\n", myid); */

   if(myid == 0)
   {
      fail = _ada_use_c2phc(41,n,d,c);    /* read the witness set */
      *dim = d[0];
      *deg = d[1];
      if(v>0) printf("The ambient dimension is %d.\n", *n); 
      if(v>0) printf("The dimension of the solution set is %d.\n", d[0]); 
      if(v>0) printf("The degree of the solution set is %d.\n", d[1]); 
      if(v>1)
      {
         printf("The embedded system at the root :\n");
         fail = _ada_use_c2phc(21,n,d,c); /* write system in container */
      }
      if(v>1)
      {
         printf("Root node has the following solutions in container :\n");
         fail = _ada_use_c2phc(31,n,d,c); /* write solutions in container */
      }
   }
 
   MPI_Bcast(n,1,MPI_INT,0,MPI_COMM_WORLD);  

   if(myid != 0)
   {
      if(v>0) printf("Node %d knowns that the dimension is %d\n", myid, *n);
      fail = _ada_use_c2phc(23,n,d,c);    /* initialize container */
      if(v>0) fail = _ada_use_c2phc(22,&m,d,c);  /* get dimension as test */
      if(v>0) printf("  and initialized container with dimension %d.\n", m);
   }

   if(v>0) printf("Node %d is about to leave read_embedding.\n", myid);
}

void monomials_broadcast ( int myid, int n )
{
   double cff[2];
   int i,j,exp[n],m[3],mm,fail;

   if(v>0) printf("Node %d has entered monomials_broadcast.\n", myid);

   for(i=1; i<=n; i++)
   {
      m[1] = i;
      if(myid == 0)
      {
         fail = _ada_use_c2phc(24,m,exp,cff);    /* get #monomials */
         mm = m[0];
         if(v>1) printf("Polynomial %d has %d monomials.\n",i,mm);
      }

      MPI_Bcast(&mm,1,MPI_INT,0,MPI_COMM_WORLD);

      m[0] = n;
      for(j=1; j<=mm; j++)    /* broadcast j-th term of i-th polynomial */
      {
         m[2] = j;
         if(myid == 0) fail = _ada_use_c2phc(25,m,exp,cff); 
         MPI_Bcast(cff,2,MPI_DOUBLE,0,MPI_COMM_WORLD); 
         MPI_Bcast(exp,n,MPI_INT,0,MPI_COMM_WORLD); 
         if(myid != 0) fail = _ada_use_c2phc(26,m,exp,cff); 
      }
   }

   if(v>0) printf("Node %d is about to leave monomials_broadcast.\n", myid);
}

void swap_slices ( int myid )
{
   double *c;
   int fail,*a,*b;

   if(v>0) printf("Node %d is swapping slices ...\n",myid);
   fail = _ada_use_c2phc(46,a,b,c); /* swap start with new slices */
}

void solutions_distribute ( int myid, int nbsols, int n, int nprocs,
                            int *solnum )
{
   int i,k,m[2],fail,dest;
   double sol[2*n+5];
   MPI_Status status;

   m[0] = n;
   if(myid == 0)
   {
      for(k=1; k<=nbsols; k++)
      {
         dest = k%(nprocs-1);
         if(dest == 0) dest=nprocs-1;    
         fail = _ada_use_c2phc(34,&k,&m[1],sol);
         m[1] = k;  /* multiplicity field is solution ID label */
         MPI_Send(&m[1],1,MPI_INT,dest,SEND_SMUL,MPI_COMM_WORLD);
         MPI_Send(sol,2*n+5,MPI_DOUBLE,dest,SEND_SSOL,MPI_COMM_WORLD);
         if(v>0) printf("Root sends solution %d to node %d.\n",k,dest);
      }
   }
   else
   {
      *solnum = nbsols/(nprocs-1);
      if(myid<=(nbsols%(nprocs-1))) *solnum = *solnum+1;
      for(k=1; k<=*solnum; k++)
      {
         MPI_Recv(&m[1],1,MPI_INT,0,SEND_SMUL,MPI_COMM_WORLD,&status);
         MPI_Recv(sol,2*n+5,MPI_DOUBLE,0,SEND_SSOL,MPI_COMM_WORLD,&status); 
         fail = _ada_use_c2phc(36,&k,m,sol);
         if(v>0) printf("Node %d receives solution %d.\n",myid,m[1]); 
      }
   }
}

void solutions_collect ( int myid, int nbsols, int n,
                         int numprocs, int mysolnum )
{
   int k,m[2],fail,dest, *a, *b;
  /* double sol[2*n+5], *c; */
   double sol[2*n+6], *c;
   MPI_Status status;

   if(myid == 0) fail = _ada_use_c2phc(37,a,b,c); /* clear solution container */

   m[0] = n;
   if(myid != 0)
   {
     for(k=1; k<=mysolnum; k++)
     {
        fail = _ada_use_c2phc(34,&k,&m[1],sol);
      /*  MPI_Send(&m[1],1,MPI_INT,0,SEND_TMUL,MPI_COMM_WORLD);
          MPI_Send(sol,2*n+5,MPI_DOUBLE,0,SEND_TSOL,MPI_COMM_WORLD); */
        sol[2*n+5] = (double) m[1];
        MPI_Send(sol,2*n+6,MPI_DOUBLE,0,SEND_TSOL,MPI_COMM_WORLD);
     }
     /* printf("Number of iterations done by node %d: %d\n", myid, k-1); */
   }
   else
   {
     for(k=1; k<=nbsols; k++)
     {
     /*  MPI_Recv
        (&m[1],1,MPI_INT,MPI_ANY_SOURCE,SEND_TMUL,MPI_COMM_WORLD,&status);
       MPI_Recv
        (sol,2*n+5,MPI_DOUBLE,MPI_ANY_SOURCE,SEND_TSOL,MPI_COMM_WORLD,&status);
     */
       MPI_Recv
        (sol,2*n+6,MPI_DOUBLE,MPI_ANY_SOURCE,SEND_TSOL,MPI_COMM_WORLD,&status);
       m[1] = (int) sol[2*n+5];
       fail = _ada_use_c2phc(36,&k,m,sol);
     }
     /* printf("Number of iterations done by root : %d\n", k-1); */
   }
  
   if(v>1)
   {
      printf("After collect, node %d has in its container:\n",myid);
      fail = _ada_use_c2phc(31,a,b,c);        /* write the target solutions */
   }
   if(myid != 0)
      fail = _ada_use_c2phc(37,a,b,c);        /* clear solution container */

   /* printf("Mysolnum is %d, Tolsolnum is %d\n", mysolnum, nbsols); */
}

void initialize_sampler ( int myid, int dim )
{
   int *b,fail;
   double *c;

   if(v>0) printf("Node %d entering initialization of sampler.\n",myid);
   fail = _ada_use_c2phc(42,&dim,b,c);        /* initialize sampler */
}

void random_complex ( double *re, double *im )
{
   double angle = ((double) rand())/RAND_MAX;
   *re = cos(angle);
   *im = sin(angle);
}

void slices_broadcast ( int myid, int k, int n )
{
   int i,j,fail;
   double r[2];

   for(i=0; i<k; i++)
      for(j=0; j<n+1; j++)
      {
         if(myid == 0)
         {
            random_complex(&r[0],&r[1]);
            if(v>0) printf("Generated r = %.15lf + %.15lf*I\n",r[0],r[1]);
         }
         MPI_Bcast(r,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
         /*if(myid != 0)*/  /* also root must know the slices! */
         {
            if(v>0) printf("Node %d received %.15lf + %.15lf*I\n",
                           myid,r[0],r[1]);
            fail = _ada_use_c2phc(43,&i,&j,r);  /* assign coefficient */
         }
      }
}

void trace_slice_broadcast ( int myid, int first )
{
   int i,j,fail;
   double r[2];

   i = 0; j = 0;

   if(myid == 0)
   {
      if(first == 1)
      {   
         r[0] = -1.0; r[1] = 0.0;
      }
      else
      {
         r[0] = +1.0; r[1] = 0.0;
      }
   }
   MPI_Bcast(r,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
   fail = _ada_use_c2phc(43,&i,&j,r);
}

void gamma_broadcast ( int myid, int n )
{
   int i,*b,fail;
   double r[2];

   for(i=0; i<n; i++)
   {
      if(myid == 0)
      {
         random_complex(&r[0],&r[1]);
         if(v>2) printf("Generated r = %.15lf + %.15lf*I\n",r[0],r[1]);
      }
      MPI_Bcast(r,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
      if(myid != 0)
      {
         if(v>2) printf("Node %d received %.15lf + %.15lf*I\n",
                        myid,r[0],r[1]);
         fail = _ada_use_c2phc(44,&i,b,r);  /* store gamma */
      }
   }
}

void broadcast_gamma ( int myid )
{
   int *a,*b,fail;
   double gamma[2];

   if(myid == 0)
   {
      double angle = ((double) rand())/RAND_MAX;
      gamma[0] = cos(angle);
      gamma[1] = sin(angle);
      if(v>0)
         printf("generated gamma = %.15lf + %.15lf*I\n",gamma[0],gamma[1]);
   }
   MPI_Bcast(gamma,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
   if(v>0)
      printf("gamma at node %d = %.15lf + %.15lf*I\n",myid,gamma[0],gamma[1]);
   fail = _ada_use_c2phc(44,a,b,gamma);       /* store gamma in phc */
}

void track_paths ( int myid, int deg, int n, int numprocs, int mysolnum,
                   double *trackwtime )
{
   int *a,*b,fail;
   double *c,startwtime,endwtime;

   if(myid != 0) 
   {
      startwtime = MPI_Wtime();
      fail = _ada_use_c2phc(45,a,b,c);          /* do path tracking */
      if(v>0) printf("Node %d done tracking %d paths.\n",myid,mysolnum);
      if(v>1) printf("Solutions computed by node %d : \n",myid);
      if(v>1) fail = _ada_use_c2phc(31,a,b,c);  /* write solutions container */
      endwtime = MPI_Wtime();
      *trackwtime += (endwtime - startwtime);
   }
   solutions_collect(myid,deg,n,numprocs,mysolnum);
   if(myid == 0)
   {
      if(v>1)
      {
         printf("Root node has solutions in its container :\n");
         fail = _ada_use_c2phc(31,a,b,c);    /* write the target solutions */
      }
      fail = _ada_use_c2phc(46,a,b,c);       /* swap start with new slices */
      fail = _ada_use_c2phc(47,a,b,c);       /* copy target system */
      fail = _ada_use_c2phc(9,a,b,c);        /* validate the solutions */
   }
}

void monodromy_breakup( int myid, int n, int dim, int deg, int nbloops,
                        int numprocs, double *trackwtime )
{
   int mysolnum,i,done;

   if(v>0) printf("Node %d knows n = %d, #loops = %d.\n",myid,n,nbloops);

   monomials_broadcast(myid,n);
   MPI_Bcast(&deg,1,MPI_INT,0,MPI_COMM_WORLD);
   if(v>0) printf("Node %d went through bcast of degree %d.\n",myid,deg);
   solutions_distribute(myid,deg,n,numprocs,&mysolnum);
   
   MPI_Bcast(&dim,1,MPI_INT,0,MPI_COMM_WORLD);
   initialize_sampler(myid,dim);

   if(myid == 0)
   {
      printf("doing sanity check...\n");
      refine_solutions();                     /* sanity check */
      initialize_monodromy(nbloops,deg,dim);  /* initialize traces */
      if(v>0) print_solutions(myid);
      store_solutions();
   }

   for(i=1; i<=2; i++)
   {
     if((v>0)&&(myid == 0)) printf("slice %d for linear trace...\n",i);
     trace_slice_broadcast(myid,i);
     gamma_broadcast(myid,n);
     track_paths(myid,deg,n,numprocs,mysolnum,trackwtime);
     if(myid == 0) { store_solutions(); restore_solutions(); }
     solutions_distribute(myid,deg,n,numprocs,&mysolnum);
     if(myid == 0) swap_slices(myid);
   }

   slices_broadcast(myid,dim,n);

   for(i=0; i<nbloops; i++)
   {
      if((v>0)&&(myid == 0)) printf("Starting loop #\%d...\n",i);
      gamma_broadcast(myid,n);
      track_paths(myid,deg,n,numprocs,mysolnum,trackwtime);
      if(myid != 0) swap_slices(myid);
      if(myid == 0)
      {
         if(v>0) print_solutions(myid);
         clear_solutions();
      }
      gamma_broadcast(myid,n);
      track_paths(myid,deg,n,numprocs,mysolnum,trackwtime);
      if(myid == 0)
      {
         if(v>0) print_solutions(myid);
         store_solutions();
         monodromy_permutation(deg,&done);
         clear_solutions();
      }

      MPI_Bcast(&done,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);   /* every node must know if done */

      if(done == 1)
      {
         if (v>0) printf("Node %d leaves the loop at step %d.\n", myid,i+1);
         if(myid == 0) printf("Executed %d monodromy loops.\n", i+1);
         return;
      }
      else
      {
         if(myid != 0) swap_slices(myid);
         slices_broadcast(myid,dim,n);               /* new slices */
         if(myid == 0) restore_solutions();
         MPI_Bcast(&deg,1,MPI_INT,0,MPI_COMM_WORLD);  
         solutions_distribute(myid,deg,n,numprocs,&mysolnum);
      }
   }
   if(myid == 0) printf("Executed %d monodromy loops.\n", i+1);
}

void initialize_monodromy ( int n, int d, int k )
{
   int fail;
   double *c;
   int dk[2];

   dk[0] = d;
   dk[1] = k;
   fail = _ada_use_c2phc(50,&n,dk,c); /* initialize Monodromy_Permutations */
}

void monodromy_permutation ( int d, int *done )
{
   int *a,*b,fail,i;
   int permutation[d];
   int nf[2];
   double *c;

   fail = _ada_use_c2phc(52,a,permutation,c);   /* compute permutation */
   printf("the permutation :");
   for (i=0; i<d; i++)
     printf(" %d",permutation[i]);
   nf[0] = d;
   nf[1] = 0;
   fail = _ada_use_c2phc(53,nf,permutation,c);  /* update decomposition */
   printf(" : %d -> %d\n",nf[0],nf[1]);
   fail = _ada_use_c2phc(55,done,b,c);          /* apply linear trace */
   fail = _ada_use_c2phc(54,a,b,c);             /* write decomposition */
}

void print_monomials ( void ) 
{
   int *d,i,j,k,n,m[3],mm,fail;
   double c[2];

   fail = _ada_use_c2phc(22,&n,d,c);

   d = (int*)calloc(n,sizeof(int));

   for(i=1; i<=n; i++)
   {
      m[1] = i;
      fail = _ada_use_c2phc(24,m,d,c);
      mm = m[0];
      m[0] = n;
      printf("Polynomial %d has %d monomials :\n",i,mm);
      for(j=1; j<=mm; j++)
      {
         m[2] = j;
         fail = _ada_use_c2phc(25,m,d,c);
         printf(" %.15e  %.15e",c[0],c[1]);
         for(k=0; k<n; k++) printf(" %d",d[k]);
         printf("\n");
      }
   }
}

void print_system ( void ) 
{
   int *a,*b,fail;
   double *c;

   fail = _ada_use_c2phc(21,a,b,c);
}

void print_solutions ( int myid )
{
   int *a,*b,fail;
   double *c;

   printf("Node %d has in its container the solutions :\n",myid);
   fail = _ada_use_c2phc(31,a,b,c); /* write solutions in container */
}

void refine_solutions ( void )
{
   int *a,*b,fail;
   double *c;

   fail = _ada_use_c2phc(9,a,b,c);       /* validate the solutions */
}

void store_solutions ( void )
{
   int *a,*b,fail;
   double *c;

   fail = _ada_use_c2phc(51,a,b,c);   /* move solutions to permutations */
}

void restore_solutions ( void )
{
   int *a,*b,fail;
   double *c;

   fail = _ada_use_c2phc(48,a,b,c);   /* first solutions in container */
}

void clear_solutions ( void )
{
   int *a,*b,fail;
   double *c;

   fail = _ada_use_c2phc(37,a,b,c); /* clear solutions in container */
}

void print_time ( double *time, int numprocs )
{
   int i;
   printf("\nTotal wall time = %lf seconds on %d processors\n",
          time[0],numprocs);
   for(i=1; i<numprocs; i++)
      printf("The wall time spent on processor %d = %lf seconds.\n",
             i,time[i]);
   free(time);
}
