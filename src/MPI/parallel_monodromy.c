/* The file "parallel_monodromy.c" contains the definitions of the prototypes
 * in "parallel_monodromy.h". */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "mpi.h"
#include "../Lib/phcpack.h"
#include "../Lib/syscon.h"
#include "../Lib/solcon.h"
#include "../Lib/witset.h"

#define SEND_SSOL 100 /* tag for sending start solution */
#define SEND_SMUL 101 /* tag for sending multiplicity of start solution */
#define SEND_TSOL 102 /* tag for sending target solution */
#define SEND_TMUL 103 /* tag for sending multiplicity of target solution */

#define v 0  /* verbose flag:
                  0 no output during the computations,
                  1 only one line messages on progress of computations
                  2 display of data, like systems and solutions */

void f_swap_slices ( int myid )
{
   int fail;

   if(v>0) printf("Node %d is swapping slices ...\n",myid);
   fail = swap_standard_slices();
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
         if(v>0) printf("Node %d received %.15lf + %.15lf*I\n",myid,r[0],r[1]);
         fail = assign_standard_coefficient_of_slice(i,j,r);
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
   fail = assign_standard_coefficient_of_slice(i,j,r);
}

void gamma_broadcast ( int myid, int n )
{
   int i,*b,fail;
   double re_gamma[n],im_gamma[n];

   if(myid == 0)
      for(i=0; i<n; i++)
      {
         random_complex(&re_gamma[i],&im_gamma[i]);
         if(v>2) printf("Generated r = %.15lf + %.15lf*I\n",
                        re_gamma[i],im_gamma[i]);
      }
   MPI_Bcast(re_gamma,n,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(im_gamma,n,MPI_DOUBLE,0,MPI_COMM_WORLD);
   if(myid != 0)
   {
      for(i=0; i<n; i++)
         if(v>2) printf("Node %d received %.15lf + %.15lf*I\n",
                        myid,re_gamma[i],im_gamma[i]);
      fail = store_standard_gamma(n,re_gamma,im_gamma);
   }
}

int witness_set_distribute ( int myid, int n, int dim, int deg, int np )
{
   int fail,mysolnum;

   if(v>0) printf("Node %d is in witness_set_distribute.\n",myid);

   monomials_broadcast(myid,n);
   if(v>1) if(myid == 0) print_solutions(myid);
   solutions_distribute(myid,deg,n,np,&mysolnum);
   fail = initialize_standard_sampler(dim);

   if(v>0) printf("Node %d is leaving witness_set_distribute.\n",myid);

   return fail;
}

void f_track_paths ( int myid, int deg, int n, int numprocs, int mysolnum,
                     double *trackwtime )
{
   int fail;
   double startwtime,endwtime;

   if(myid != 0) 
   {
      startwtime = MPI_Wtime();
      fail = standard_sample_to_new_slices();      /* do path tracking */
      if(v>0) printf("Node %d done tracking %d paths.\n",myid,mysolnum);
      if(v>1) printf("Solutions computed by node %d : \n",myid);
      if(v>1) fail = solcon_write_solutions();
      endwtime = MPI_Wtime();
      *trackwtime += (endwtime - startwtime);
   }
   solutions_collect(myid,deg,n,numprocs,mysolnum);
   if(myid == 0)
   {
      if(v>1)
      {
         printf("Root node has solutions in its container :\n");
         fail = solcon_write_solutions();
      }
      fail = swap_standard_slices();
      fail = standard_witness_set_to_system_container();
      fail = validate_solutions();
   }
}

int build_trace_grid ( int myid, int n, int dim, int deg,
                       int numprocs, double *trackwtime )
{
   int mysolnum,i,done,fail;
   double err,dis;
   double tol = 1.0e-8;

   if(v>0) printf("Node %d knows n = %d.\n",myid,n);

   if(myid == 0)
   {
      fail = validate_solutions();               /* sanity check */
      if(v>1) print_solutions(myid);
      fail = store_standard_solutions();
   }

   solutions_distribute(myid,deg,n,numprocs,&mysolnum);
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

   if(myid == 0)
   {
      fail = standard_trace_grid_diagnostics(&err,&dis);
      printf("maximal error on samples in trace grid : %.2e\n",err);
      printf("minimal distance between samples in grid : %.2e\n",dis);
      if((err > tol)||(dis < tol)) return 1;
   }

   return 0;
}

void trace_grid_broadcast ( int myid, int d, int n )
{
   int i,fail;

   if(v > 0) printf("Node %d is entering trace_grid_broadcast.\n",myid);

   for(i=0; i<3; i++)
   {
      if(myid == 0) fail = retrieve_solutions_on_grid(i);
      MPI_Barrier(MPI_COMM_WORLD);  /* everybody must wait for node 0 */
      if(v>0) printf("Node %d before broadcasting solutions %d.\n",myid,i);
      solutions_broadcast(myid,d,n);
      if(v>0) printf("Node %d after broadcasting solutions %d.\n",myid,i);

      if(i == 0)
         f_swap_slices(myid);
      else
         trace_slice_broadcast(myid,i);

      if(myid != 0) 
      {
         f_swap_slices(myid);
	 fail = store_standard_solutions();
      }
   }
}

void all_slices_broadcast ( int myid, int nb_slices, int dim, int n )
{
   int i,j,k,fail;
   int nb_cff = (n+1)*dim; /* number of complex coefficients in one slice */
   int nb_cff2 = 2*nb_cff; /* a complex number is stored as two doubles */
   double cff[nb_cff2];    /* so a slice has twice as many doubles */
   double re,im;           /* real and imaginary part of a complex */

   fail = initialize_hyperplane_sections(nb_slices);

   for(i=0; i<nb_slices; i++)
   {
      if(myid == 0)
      {
         for(j=0,k=0; j<nb_cff; j++,k+=2)
         {
            random_complex(&re,&im);
            cff[k] = re;
            cff[k+1] = im;
         }
       /* printf("Node 0 has stored slice %d with cff1 = %.2le\n",i,cff[0]); */
      }

      MPI_Bcast(cff,nb_cff2,MPI_DOUBLE,0,MPI_COMM_WORLD);

      fail = store_new_hyperplane_sections(nb_cff2,dim,n,cff);

   /* if(myid == 1) */ /* just for testing purposes */
    /*     printf("Node 1 has stored slice %d with cff1 = %.2le\n",i,cff[0]); */
   }
}

int track_paths_to_new_slices ( int myid, int nb_slices, int d, int k, int n,
                                int numprocs, double *trackwtime )
{
   int i,fail,mysolnum;

   for(i=1; i<=nb_slices; i++)
   {
/* if(v>0) */ if(myid == 0) printf("Tracking to slice #\%d...\n",i);

      if(myid == 0) fail = restore_standard_solutions();
      /* gamma_broadcast(myid,n); */
      fail = set_target_hyperplane_sections(i);
      solutions_distribute(myid,d,n,numprocs,&mysolnum);
      f_track_paths(myid,d,n,numprocs,mysolnum,trackwtime);
      MPI_Barrier(MPI_COMM_WORLD);  /* every body waits for broadcast */

/* if(v>0) */ if(myid == 0) printf("Broadcasting the solutions ...\n");

      if(myid == 0) fail = store_standard_solutions();
      solutions_broadcast(myid,d,n);
      if(myid != 0) fail = store_standard_solutions();
   }
}

void interactive_trace_test ( int myid )
{
   int k,*f,i,fail,dest;
   double d,d1;
   char ans = 'y';
   MPI_Status status;

   if(v>0) printf("Node %d has entered interactive trace test...\n",myid);
   if(myid == 0)
   {
      printf("Give destination node for test : ");
      scanf("%d",&dest);
   }
   MPI_Bcast(&dest,1,MPI_INT,0,MPI_COMM_WORLD);

   while(ans == 'y')
   {
      if(myid == 0)
      {
         printf("Give the degree of the test factor : ");
         scanf("%d",&k);
         f = (int*)calloc(k,sizeof(int));
         for(i=0; i<k; i++)
         {
            printf("Give label %d of the factor : ",i+1);
            scanf("%d",&f[i]);
         }
         printf("The labels");
         for(i=0; i<k; i++) printf(" %d",f[i]);

         fail = standard_trace_sum_difference(k,f,&d);
         printf(" have trace sum difference %.2e\n",d);

         if(dest != 0)
         {
            MPI_Send(&k,1,MPI_INT,dest,SEND_SMUL,MPI_COMM_WORLD);
            MPI_Send(f,k,MPI_INT,dest,SEND_SMUL,MPI_COMM_WORLD);
            MPI_Recv(&d1,1,MPI_DOUBLE,dest,0,MPI_COMM_WORLD,&status);
            printf("Node 0 receives %.2e from Node %d.\n",d1,dest);
         }
      }
      if(dest != 0) if(myid == dest)
      {
         MPI_Recv(&k,1,MPI_INT,0,SEND_SMUL,MPI_COMM_WORLD,&status);
         f = (int*)calloc(k,sizeof(int));
         MPI_Recv(f,k,MPI_INT,0,SEND_SMUL,MPI_COMM_WORLD,&status);
         fail = standard_trace_sum_difference(k,f,&d);
	 printf("Node %d computes trace sum difference %.2e\n",myid,d);
         MPI_Send(&d,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
      }
      if(myid == 0)
      {
         printf("Do you want more tests ? (y/n) ");
         scanf("%c",&ans);  /* skip new line symbol */
         scanf("%c",&ans);
      }
      MPI_Bcast(&ans,1,MPI_CHAR,0,MPI_COMM_WORLD);
   }
}

void interactive_in_slice_test ( int myid )
{
   int label,slice,p,p1,dest,fail;
   char ans = 'y';
   MPI_Status status;

   if(myid == 0)
   {
      printf("Give destination node for test : ");
      scanf("%d",&dest);
   }
   MPI_Bcast(&dest,1,MPI_INT,0,MPI_COMM_WORLD);

   while(ans == 'y')
   {
      if(myid == 0)
      {
         printf("Give a label to a solution : "); scanf("%d",&label);
         printf("Give a slice number : "); scanf("%d",&slice);
         fail = in_slice(label,slice,&p);
         printf("The position of %d in slice %d is %d.\n",label,slice,p);

         MPI_Send(&label,1,MPI_INT,dest,SEND_SMUL,MPI_COMM_WORLD);
         MPI_Send(&slice,1,MPI_INT,dest,SEND_SMUL,MPI_COMM_WORLD);

         MPI_Recv(&p1,1,MPI_INT,dest,SEND_SMUL,MPI_COMM_WORLD,&status);
         printf("Node 0 receives position %d from Node %d.\n",p1,dest);
      }
      if(myid == dest)
      {
         MPI_Recv(&label,1,MPI_INT,0,SEND_SMUL,MPI_COMM_WORLD,&status);
         MPI_Recv(&slice,1,MPI_INT,0,SEND_SMUL,MPI_COMM_WORLD,&status);
         fail = in_slice(label,slice,&p);
         printf("Node %d finds position of %d in slice %d to be %d.\n",
                dest,label,slice,p);
         MPI_Send(&p,1,MPI_INT,0,SEND_SMUL,MPI_COMM_WORLD);
      }
      if(myid == 0)
      {
         printf("Do you want more tests ? (y/n) ");
         scanf("%c",&ans); /* skip new line symbol */
         scanf("%c",&ans);
      }
      MPI_Bcast(&ans,1,MPI_CHAR,0,MPI_COMM_WORLD);
   }
}

void interactive_track_to_make_loops ( int myid )
{
   int start_slice,target_slice,start_label,target_label,dest,fail;
   char ans = 'y';
   MPI_Status status;

   if(myid == 0)
   {
      printf("Give destination node for test : ");
      scanf("%d",&dest);
   }
   MPI_Bcast(&dest,1,MPI_INT,0,MPI_COMM_WORLD);

   while(ans == 'y')
   {
      if(myid == 0)
      {
         printf("Give the index of the start slice : ");
         scanf("%d",&start_slice);
         printf("Give the index of the target slice : ");
         scanf("%d",&target_slice);
         printf("Give the label of the start solution : ");
         scanf("%d",&start_label);
         printf("Node 0 is sending (%d,%d,%d) to Node %d.\n",
                start_slice,target_slice,start_label,dest);

         MPI_Send(&start_slice,1,MPI_INT,dest,SEND_SMUL,MPI_COMM_WORLD);
         MPI_Send(&target_slice,1,MPI_INT,dest,SEND_SMUL,MPI_COMM_WORLD);
         MPI_Send(&start_label,1,MPI_INT,dest,SEND_SMUL,MPI_COMM_WORLD);

         printf("Node 0 has sent (%d,%d,%d) to Node %d.\n",
                start_slice,target_slice,start_label,dest);

         MPI_Recv(&target_label,1,MPI_INT,dest,SEND_SMUL,
                  MPI_COMM_WORLD,&status);

         printf("Node 0 receives label of target %d from Node %d.\n",
                target_label,dest);
      }
      if(myid == dest)
      {
         MPI_Recv(&start_slice,1,MPI_INT,0,SEND_SMUL,MPI_COMM_WORLD,&status);
         MPI_Recv(&target_slice,1,MPI_INT,0,SEND_SMUL,MPI_COMM_WORLD,&status);
         MPI_Recv(&start_label,1,MPI_INT,0,SEND_SMUL,MPI_COMM_WORLD,&status);

         printf("Node %d received (%d,%d,%d) from Node 0.\n",
                dest,start_slice,target_slice,start_label);

         fail = standard_sample_loop
                  (start_slice,target_slice,start_label,&target_label);
         printf("Node %d computes label of the target solution : %d.\n",
                dest,target_label);

         MPI_Send(&target_label,1,MPI_INT,0,SEND_SMUL,MPI_COMM_WORLD);
      }
      if(myid == 0)
      {
         printf("Do you want more tests ? (y/n) ");
         scanf("%c",&ans); /* skip new line symbol */
         scanf("%c",&ans);
      }
      MPI_Bcast(&ans,1,MPI_CHAR,0,MPI_COMM_WORLD);
   }
}
