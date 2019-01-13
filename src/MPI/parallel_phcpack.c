/* The file "parallel_phcpack.c" collects the definition of the functions
 * with prototypes in "parallel_phcpack.h". */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "mpi.h"
#include "../Lib/phcpack.h"
#include "../Lib/syscon.h"
#include "../Lib/solcon.h"
#include "../Lib/jump_track.h"
#include "parallel_phcpack.h"

#define v 0  /* verbose flag:
                  0 no output during the computations,
                  1 only one line messages on progress of computations
                  2 display of data, like systems and solutions */

void dimension_broadcast ( int myid, int *n )
{
   int fail,m;

   if(v>0) printf("Node %d has entered broadcast_dimension.\n", myid);

   MPI_Bcast(n,1,MPI_INT,0,MPI_COMM_WORLD);  

   if(myid != 0)
   {
      if(v>0) printf("Node %d knowns that the dimension is %d\n", myid, *n);
      /* initialize container */
      fail = syscon_initialize_number_of_standard_polynomials(*n);
      if(v>0)
      {  /* get dimension as test */
         fail = syscon_number_of_standard_polynomials(&m);
         printf("  and initialized container with dimension %d.\n", m);
      }
   }

   if(v>0) printf("Node %d is leaving broadcast_dimension.\n", myid);
}

void standard_dimensions_broadcast ( int myid, int *nbequ, int *nbvar )
{
   int fail,m;

   if(v>0) printf("Node %d has entered dimensions_broadcast.\n", myid);

   MPI_Bcast(nbequ,1,MPI_INT,0,MPI_COMM_WORLD);  
   MPI_Bcast(nbvar,1,MPI_INT,0,MPI_COMM_WORLD);  

   if(myid != 0)
   {
      if(v>0) 
      {
          printf("Node %d has the number of equations as %d.\n", myid, *nbequ);
          printf("Node %d has the number of variables as %d.\n", myid, *nbvar);
      }
      /* initialize container */
      fail = syscon_initialize_number_of_standard_polynomials(*nbequ);
      if(v>0)
      {  /* get dimension as test */
         fail = syscon_number_of_standard_polynomials(&m);
         printf("  and initialized container with dimension %d.\n", m);
      }
   }

   if(v>0) printf("Node %d is leaving dimensions_broadcast.\n", myid);
}

void dobldobl_dimensions_broadcast ( int myid, int *nbequ, int *nbvar )
{
   int fail,m;

   if(v>0) printf("Node %d has entered dimensions_broadcast.\n", myid);

   MPI_Bcast(nbequ,1,MPI_INT,0,MPI_COMM_WORLD);  
   MPI_Bcast(nbvar,1,MPI_INT,0,MPI_COMM_WORLD);  

   if(myid != 0)
   {
      if(v>0) 
      {
          printf("Node %d has the number of equations as %d.\n", myid, *nbequ);
          printf("Node %d has the number of variables as %d.\n", myid, *nbvar);
      }
      /* initialize container */
      fail = syscon_initialize_number_of_dobldobl_polynomials(*nbequ);
      if(v>0)
      {  /* get dimension as test */
         fail = syscon_number_of_dobldobl_polynomials(&m);
         printf("  and initialized container with dimension %d.\n", m);
      }
   }

   if(v>0) printf("Node %d is leaving dimensions_broadcast.\n", myid);
}

void quaddobl_dimensions_broadcast ( int myid, int *nbequ, int *nbvar )
{
   int fail,m;

   if(v>0) printf("Node %d has entered dimensions_broadcast.\n", myid);

   MPI_Bcast(nbequ,1,MPI_INT,0,MPI_COMM_WORLD);  
   MPI_Bcast(nbvar,1,MPI_INT,0,MPI_COMM_WORLD);  

   if(myid != 0)
   {
      if(v>0) 
      {
          printf("Node %d has the number of equations as %d.\n", myid, *nbequ);
          printf("Node %d has the number of variables as %d.\n", myid, *nbvar);
      }
      /* initialize container */
      fail = syscon_initialize_number_of_quaddobl_polynomials(*nbequ);
      if(v>0)
      {  /* get dimension as test */
         fail = syscon_number_of_quaddobl_polynomials(&m);
         printf("  and initialized container with dimension %d.\n", m);
      }
   }

   if(v>0) printf("Node %d is leaving dimensions_broadcast.\n", myid);
}

void monomials_broadcast ( int myid, int n )
{
   double cff[2];
   int i,j,exp[n],mm,fail;

   if(v>0) printf("Node %d has entered monomials_broadcast.\n", myid);

   for(i=1; i<=n; i++)
   {
      if(myid == 0)
      {
         fail = syscon_number_of_standard_terms(i,&mm); /* get #monomials */
         if(v>1) printf("Polynomial %d has %d monomials.\n",i,mm);
      }

      MPI_Bcast(&mm,1,MPI_INT,0,MPI_COMM_WORLD);

      for(j=1; j<=mm; j++)    /* broadcast j-th term of i-th polynomial */
      {
         if(myid == 0) fail = syscon_retrieve_standard_term(i,j,n,exp,cff); 
         MPI_Bcast(cff,2,MPI_DOUBLE,0,MPI_COMM_WORLD); 
         MPI_Bcast(exp,n,MPI_INT,0,MPI_COMM_WORLD); 
         if(myid != 0) fail = syscon_add_standard_term(i,n,exp,cff);
      }
   }

   if(v>0) printf("Node %d is about to leave monomials_broadcast.\n", myid);
}

void dobldobl_monomials_broadcast ( int myid, int n )
{
   double cff[4];
   int i,j,exp[n],mm,fail;

   if(v>0) printf("Node %d has entered monomials_broadcast.\n", myid);

   for(i=1; i<=n; i++)
   {
      if(myid == 0)
      {
         fail = syscon_number_of_dobldobl_terms(i,&mm); /* get #monomials */
         if(v>1) printf("Polynomial %d has %d monomials.\n",i,mm);
      }

      MPI_Bcast(&mm,1,MPI_INT,0,MPI_COMM_WORLD);

      for(j=1; j<=mm; j++)    /* broadcast j-th term of i-th polynomial */
      {
         if(myid == 0) fail = syscon_retrieve_dobldobl_term(i,j,n,exp,cff); 
         MPI_Bcast(cff,4,MPI_DOUBLE,0,MPI_COMM_WORLD); 
         MPI_Bcast(exp,n,MPI_INT,0,MPI_COMM_WORLD); 
         if(myid != 0) fail = syscon_add_dobldobl_term(i,n,exp,cff);
      }
   }

   if(v>0) printf("Node %d is about to leave monomials_broadcast.\n", myid);
}

void quaddobl_monomials_broadcast ( int myid, int n )
{
   double cff[8];
   int i,j,exp[n],mm,fail;

   if(v>0) printf("Node %d has entered monomials_broadcast.\n", myid);

   for(i=1; i<=n; i++)
   {
      if(myid == 0)
      {
         fail = syscon_number_of_quaddobl_terms(i,&mm); /* get #monomials */
         if(v>1) printf("Polynomial %d has %d monomials.\n",i,mm);
      }

      MPI_Bcast(&mm,1,MPI_INT,0,MPI_COMM_WORLD);

      for(j=1; j<=mm; j++)    /* broadcast j-th term of i-th polynomial */
      {
         if(myid == 0) fail = syscon_retrieve_quaddobl_term(i,j,n,exp,cff); 
         MPI_Bcast(cff,8,MPI_DOUBLE,0,MPI_COMM_WORLD); 
         MPI_Bcast(exp,n,MPI_INT,0,MPI_COMM_WORLD); 
         if(myid != 0) fail = syscon_add_quaddobl_term(i,n,exp,cff);
      }
   }

   if(v>0) printf("Node %d is about to leave monomials_broadcast.\n", myid);
}

int start_system_broadcast ( int myid, int n, int *nbsols )
{
   int fail;

   if(myid == 0)                           /* manager reads start system */
   {
      fail = read_standard_start_system();
      fail = copy_start_system_to_container();
      fail = copy_start_solutions_to_container();
   }
   else
      /* initialize system container */
      fail = syscon_initialize_number_of_standard_polynomials(n);

   monomials_broadcast(myid,n);            /* broadcast container */

   if(myid != 0)                           /* copy result of broadcast */
   {
      fail = copy_container_to_start_system();
      fail = syscon_clear_standard_system(); /* clear system container */
   }
   else
      fail = solcon_number_of_standard_solutions(nbsols);

   return fail;
}

int start_system_broadcast_without_solutions ( int myid, int n, int *nbsols )
{
   int fail,dim;

   if(myid == 0)                           /* manager reads start system */
   {
      fail = read_start_system_without_solutions();
      fail = copy_start_system_to_container();
     /* fail = solcon_open_solution_input_file(); */
      fail = solcon_scan_solution_banner();
      fail = solcon_read_solution_dimensions(nbsols,&dim);
      printf("read the solution dimensions : %d\n",*nbsols);
   }
   else
      /* initialize system container */
      fail = syscon_initialize_number_of_standard_polynomials(n);

   monomials_broadcast(myid,n);            /* broadcast container */

   if(myid != 0)                           /* copy result of broadcast */
   {
      fail = copy_container_to_start_system();
      fail = syscon_clear_standard_system();  /* clear system container */
   }

   return fail;
}

int start_in_container_broadcast ( int myid, int n )
{
   int fail;

   if(myid == 0)
   {
      if(v>0) printf("Manager copies target system to container.\n");
      fail = copy_start_system_to_container();
   }
   else
   {
      if(v>0) printf("Node %d initializes container with n = %d.\n",myid,n);
      fail = syscon_initialize_number_of_standard_polynomials(n);
   }

   monomials_broadcast(myid,n);

   if(myid != 0)
   {
      if(v>0) printf("Node %d copies container to target.\n",myid);
      fail = copy_container_to_start_system();
      fail = syscon_clear_standard_system();
   }

   return fail;
}

int target_in_container_broadcast ( int myid, int n )
{
   int fail;

   if(myid == 0)
   {
      if(v>0) printf("Manager copies target system to container.\n");
      fail = copy_target_system_to_container();
   }
   else
   {
      if(v>0) printf("Node %d initializes container with n = %d.\n",myid,n);
      fail = syscon_initialize_number_of_standard_polynomials(n);
   }

   monomials_broadcast(myid,n);

   if(myid != 0)
   {
      if(v>0) printf("Node %d copies container to target.\n",myid);
      fail = copy_container_to_target_system();
      fail = syscon_clear_standard_system();
   }

   return fail;
}

int homotopy_broadcast ( int myid, int n )
{
   int fail;

   fail = start_in_container_broadcast(myid,n);
   fail = target_in_container_broadcast(myid,n);

   return fail;
}

int named_start_system_broadcast
    ( int myid, int n, int *nbsols, int kind, int nc, char name[nc] )
{
   int fail,dim;

   if(myid == 0)                           /* manager reads start system */
   {
      if(kind == 2)
         fail = read_named_linear_product_start_system(nc,name);
      else
         fail = read_named_start_without_solutions(nc,name);

      fail = copy_start_system_to_container();

      if(kind >= 3)                        /* cheater's homotopy */
      {
         fail = solcon_scan_solution_banner();
         fail = solcon_read_solution_dimensions(nbsols,&dim);
      }
      else
         fail = syscon_total_degree(nbsols);

      if(v>0) printf("the number of solution paths : %d\n",*nbsols);
   }
   else
      /* initialize system container */
      fail = syscon_initialize_number_of_standard_polynomials(n);

   monomials_broadcast(myid,n);            /* broadcast container */

   if(myid != 0)                           /* copy result of broadcast */
   {
      fail = copy_container_to_start_system();
      fail = syscon_clear_standard_system(); /* clear system container */
   }

   return fail;
}

void solutions_broadcast ( int myid, int nbsols, int n )
{
   int i,m,fail;
   int s = 2*n+5;
   double sol[s];

   if(myid != 0) fail = solcon_clear_standard_solutions();

   if(myid == 0) printf("Broadcasting");
   for(i=1; i<= nbsols; i++)
   {
      if(myid == 0) fail = solcon_retrieve_standard_solution(n,i,&m,sol);
      MPI_Bcast(&m,1,MPI_INT,0,MPI_COMM_WORLD);
      if(myid == 0) printf(" %d",m);
      MPI_Bcast(sol,s,MPI_DOUBLE,0,MPI_COMM_WORLD);
      if(myid != 0) fail = solcon_append_standard_solution(n,m,sol);
   }
   if(myid == 0) printf("\n");
}

void solutions_distribute
 ( int myid, int nbsols, int n, int nprocs, int *solnum )
{
   const int lensol = 2*n+5;
   int i,k,m,fail,dest,idx;
   double sol[lensol];
   MPI_Status status;

   if(myid == 0)
   {
      fail = solcon_retrieve_next_standard_initialize();
      for(k=1; k<=nbsols; k++)
      {
         dest = k%(nprocs-1);
         if(dest == 0) dest=nprocs-1;    
         // fail = solcon_retrieve_standard_solution(n,k,&m,sol);
         fail = solcon_retrieve_next_standard_solution(n,&idx,&m,sol);
         m = k;       /* multiplicity field labels solution ID */
         MPI_Send(&m,1,MPI_INT,dest,SEND_SMUL,MPI_COMM_WORLD);
         MPI_Send(sol,lensol,MPI_DOUBLE,dest,SEND_SSOL,MPI_COMM_WORLD);
         if(v>0) printf("Root sends solution %d to node %d.\n",k,dest);
      }
   }
   else
   {
      *solnum = nbsols/(nprocs-1);
      if(myid<=(nbsols%(nprocs-1))) *solnum = *solnum+1;
      for(k=1; k<=*solnum; k++)
      {
         MPI_Recv(&m,1,MPI_INT,0,SEND_SMUL,MPI_COMM_WORLD,&status);
         MPI_Recv(sol,lensol,MPI_DOUBLE,0,SEND_SSOL,MPI_COMM_WORLD,&status); 
         fail = solcon_append_standard_solution(n,m,sol);
         if(v>0) printf("Node %d receives solution %d.\n",myid,m); 
      }
   }
}

void dobldobl_solutions_distribute
 ( int myid, int nbsols, int n, int nprocs, int *solnum, int verbose )
{
   const int lensol = 4*n+10;
   int i,k,m,fail,dest,idx;
   double sol[lensol];
   MPI_Status status;

   if(myid == 0)
   {
      fail = solcon_retrieve_next_dobldobl_initialize();
      for(k=1; k<=nbsols; k++)
      {
         dest = k%(nprocs-1);
         if(dest == 0) dest=nprocs-1;    
         // fail = solcon_retrieve_dobldobl_solution(n,k,&m,sol);
         fail = solcon_retrieve_next_dobldobl_solution(n,&idx,&m,sol);
         m = k;       /* multiplicity field labels solution ID */
         MPI_Send(&m,1,MPI_INT,dest,SEND_SMUL,MPI_COMM_WORLD);
         MPI_Send(sol,lensol,MPI_DOUBLE,dest,SEND_SSOL,MPI_COMM_WORLD);
         if(verbose>0) printf("Root sends solution %d to node %d.\n",k,dest);
      }
   }
   else
   {
      *solnum = nbsols/(nprocs-1);
      if(myid<=(nbsols%(nprocs-1))) *solnum = *solnum+1;
      for(k=1; k<=*solnum; k++)
      {
         MPI_Recv(&m,1,MPI_INT,0,SEND_SMUL,MPI_COMM_WORLD,&status);
         MPI_Recv(sol,lensol,MPI_DOUBLE,0,SEND_SSOL,MPI_COMM_WORLD,&status); 
         fail = solcon_append_dobldobl_solution(n,m,sol);
         if(verbose>0) printf("Node %d receives solution %d.\n",myid,m); 
      }
   }
}

void quaddobl_solutions_distribute
 ( int myid, int nbsols, int n, int nprocs, int *solnum, int verbose )
{
   const int lensol = 8*n+20;
   int i,k,m,fail,dest,idx;
   double sol[lensol];
   MPI_Status status;

   if(myid == 0)
   {
      fail = solcon_retrieve_next_quaddobl_initialize();
      for(k=1; k<=nbsols; k++)
      {
         dest = k%(nprocs-1);
         if(dest == 0) dest=nprocs-1;    
         // fail = solcon_retrieve_quaddobl_solution(n,k,&m,sol);
         fail = solcon_retrieve_next_quaddobl_solution(n,&idx,&m,sol);
         m = k;       /* multiplicity field labels solution ID */
         MPI_Send(&m,1,MPI_INT,dest,SEND_SMUL,MPI_COMM_WORLD);
         MPI_Send(sol,lensol,MPI_DOUBLE,dest,SEND_SSOL,MPI_COMM_WORLD);
         if(verbose>0) printf("Root sends solution %d to node %d.\n",k,dest);
      }
   }
   else
   {
      *solnum = nbsols/(nprocs-1);
      if(myid<=(nbsols%(nprocs-1))) *solnum = *solnum+1;
      for(k=1; k<=*solnum; k++)
      {
         MPI_Recv(&m,1,MPI_INT,0,SEND_SMUL,MPI_COMM_WORLD,&status);
         MPI_Recv(sol,lensol,MPI_DOUBLE,0,SEND_SSOL,MPI_COMM_WORLD,&status); 
         fail = solcon_append_quaddobl_solution(n,m,sol);
         if(verbose>0) printf("Node %d receives solution %d.\n",myid,m); 
      }
   }
}

void solutions_collect
 ( int myid, int nbsols, int n, int numprocs, int mysolnum )
{
   const int lensol = 2*n+6;
   int k,m,fail,dest;
   double sol[lensol];  /* very important to send "m" along with "sol" */
   MPI_Status status;

   if(myid == 0) fail = solcon_clear_standard_solutions(); 

   if(myid != 0)
   {
      for(k=1; k<=mysolnum; k++)
      {
         fail = solcon_retrieve_standard_solution(n,k,&m,sol);
         sol[lensol-1] = (double) m;  /* m is last entry of array sol */
         MPI_Send(sol,lensol,MPI_DOUBLE,0,SEND_TSOL,MPI_COMM_WORLD);
      }
   }
   else
   {
      for(k=1; k<=nbsols; k++)
      {
         MPI_Recv(sol,lensol,MPI_DOUBLE,MPI_ANY_SOURCE,SEND_TSOL,
                  MPI_COMM_WORLD,&status);
         m = (int) sol[lensol-1];
         fail = solcon_append_standard_solution(n,m,sol);
      }
   }
  
   if(v>1)
   {
      printf("After collect, node %d has in its container:\n",myid);
      fail = solcon_write_standard_solutions();
   }

   if(myid != 0) fail = solcon_clear_standard_solutions();
}

void dobldobl_solutions_collect
 ( int myid, int nbsols, int n, int numprocs, int mysolnum )
{
   const int lensol = 4*n+11;
   int k,m,fail,dest;
   double sol[lensol];  /* very important to send "m" along with "sol" */
   MPI_Status status;

   if(myid == 0) fail = solcon_clear_dobldobl_solutions(); 

   if(myid != 0)
   {
      for(k=1; k<=mysolnum; k++)
      {
         fail = solcon_retrieve_dobldobl_solution(n,k,&m,sol);
         sol[lensol-1] = (double) m;  /* m is last entry of array sol */
         MPI_Send(sol,lensol,MPI_DOUBLE,0,SEND_TSOL,MPI_COMM_WORLD);
      }
   }
   else
   {
      for(k=1; k<=nbsols; k++)
      {
         MPI_Recv(sol,lensol,MPI_DOUBLE,MPI_ANY_SOURCE,SEND_TSOL,
                  MPI_COMM_WORLD,&status);
         m = (int) sol[lensol-1];
         fail = solcon_append_dobldobl_solution(n,m,sol);
      }
   }
  
   if(v>1)
   {
      printf("After collect, node %d has in its container:\n",myid);
      fail = solcon_write_dobldobl_solutions();
   }

   if(myid != 0) fail = solcon_clear_dobldobl_solutions();
}

void quaddobl_solutions_collect
 ( int myid, int nbsols, int n, int numprocs, int mysolnum )
{
   const int lensol = 8*n+21;
   int k,m,fail,dest;
   double sol[lensol];  /* very important to send "m" along with "sol" */
   MPI_Status status;

   if(myid == 0) fail = solcon_clear_quaddobl_solutions(); 

   if(myid != 0)
   {
      for(k=1; k<=mysolnum; k++)
      {
         fail = solcon_retrieve_quaddobl_solution(n,k,&m,sol);
         sol[lensol-1] = (double) m;  /* m is last entry of array sol */
         MPI_Send(sol,lensol,MPI_DOUBLE,0,SEND_TSOL,MPI_COMM_WORLD);
      }
   }
   else
   {
      for(k=1; k<=nbsols; k++)
      {
         MPI_Recv(sol,lensol,MPI_DOUBLE,MPI_ANY_SOURCE,SEND_TSOL,
                  MPI_COMM_WORLD,&status);
         m = (int) sol[lensol-1];
         fail = solcon_append_quaddobl_solution(n,m,sol);
      }
   }
  
   if(v>1)
   {
      printf("After collect, node %d has in its container:\n",myid);
      fail = solcon_write_quaddobl_solutions();
   }

   if(myid != 0) fail = solcon_clear_quaddobl_solutions();
}

void print_monomials ( void ) 
{
   int *d,i,j,k,n,mm,fail;
   double c[2];

   fail = syscon_number_of_standard_polynomials(&n);

   d = (int*)calloc(n,sizeof(int));

   for(i=1; i<=n; i++)
   {
      fail = syscon_number_of_standard_terms(i,&mm);
      printf("Polynomial %d has %d monomials :\n",i,mm);
      for(j=1; j<=mm; j++)
      {
         fail = syscon_retrieve_standard_term(i,j,n,d,c);
         printf(" %.15e  %.15e",c[0],c[1]);
         for(k=0; k<n; k++) printf(" %d",d[k]);
         printf("\n");
      }
   }
}

void write_solution_banner_to_defined_output_file ( int nbsols, int n )
{
   int fail;

   fail = write_string_to_defined_output_file(17,"\nTHE SOLUTIONS :\n");
   fail = write_integers_to_defined_output_file(1,&nbsols);
   fail = write_string_to_defined_output_file(1," ");
   fail = write_integers_to_defined_output_file(1,&n);
   fail = write_string_to_defined_output_file
     (61,"\n===========================================================\n");
}

void print_solutions ( int myid )
{
   int *a,*b,fail;
   double *c;

   printf("Node %d has in its container the solutions :\n",myid);
   fail = solcon_write_standard_solutions();
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

void write_time_and_paths_to_defined_output_file
       ( int p, double time[p], int paths[p] )
{
   int i,fail;

   fail = write_string_to_defined_output_file(19,"\nTotal wall time = ");
   fail = write_doubles_to_defined_output_file(1,&time[0]);
   fail = write_string_to_defined_output_file(12," seconds on ");
   fail = write_integers_to_defined_output_file(1,&p);
   fail = write_string_to_defined_output_file(13," processors.\n"),

   fail = write_string_to_defined_output_file
            (32,"Total number of paths tracked : ");
   fail = write_integers_to_defined_output_file(1,&paths[0]);
   fail = write_string_to_defined_output_file(2,".\n");

   for(i=1; i<p; i++)
   {
      fail = write_string_to_defined_output_file
               (33,"The wall time spent on processor ");
      fail = write_integers_to_defined_output_file(1,&i);
      fail = write_string_to_defined_output_file(3," = ");
      fail = write_doubles_to_defined_output_file(1,&time[i]);
      fail = write_string_to_defined_output_file(10," seconds,\n");

      fail = write_string_to_defined_output_file
               (41,"    number of paths tracked by processor ");
      fail = write_integers_to_defined_output_file(1,&i);
      fail = write_string_to_defined_output_file(3," : ");
      fail = write_integers_to_defined_output_file(1,&paths[i]);
      fail = write_string_to_defined_output_file(2,".\n");
   }
}
