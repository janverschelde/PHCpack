/* The following is a simple interface of MPI with PHCpack:
 * after prompting the user for a list of solutions,
 * the solution list is broadcasted to all the nodes. */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "mpi.h"

extern void adainit( void );
extern int _ada_use_solcon ( int task, int *a, int *b, double *c );
extern void adafinal( void );

void standard_dimension_broadcast ( int myid, int *len, int *n );
/* 
 * DESCRIPTION :
 *   gets the length of the solution list and the dimension n
 *   from the root node, and broadcasts then to all the nodes,
 *   for the standard double precision. */

void dobldobl_dimension_broadcast ( int myid, int *len, int *n );
/* 
 * DESCRIPTION :
 *   gets the length of the solution list and the dimension n
 *   from the root node, and broadcasts then to all the nodes,
 *   for the double double precision. */

void quaddobl_dimension_broadcast ( int myid, int *len, int *n );
/* 
 * DESCRIPTION :
 *   gets the length of the solution list and the dimension n
 *   from the root node, and broadcasts then to all the nodes,
 *   for the quad double precision. */

void standard_solutions_broadcast ( int myid, int len, int n );
/* 
 * DESCRIPTION :
 *   broadcasts all solutions in the container at the root
 *   to the solutions container at each node,
 *   for the standard double precision. */

void dobldobl_solutions_broadcast ( int myid, int len, int n );
/* 
 * DESCRIPTION :
 *   broadcasts all solutions in the container at the root
 *   to the solutions container at each node,
 *   for the double double precision. */

void quaddobl_solutions_broadcast ( int myid, int len, int n );
/* 
 * DESCRIPTION :
 *   broadcasts all solutions in the container at the root
 *   to the solutions container at each node,
 *   for the quad double precision. */

void standard_print_broadcast ( int myid );
/* 
 * DESCRIPTION :
 *   prints the result of the broadcast at the node 1,
 *   for the standard double precision. */

void dobldobl_print_broadcast ( int myid );
/* 
 * DESCRIPTION :
 *   prints the result of the broadcast at the node 1,
 *   for the double double precision. */

void quaddobl_print_broadcast ( int myid );
/* 
 * DESCRIPTION :
 *   prints the result of the broadcast at the node 1,
 *   for the double double precision. */

void standard_print_solutions ( void );
/* 
 * DESCRIPTION :
 *   prints the solutions in the container,
 *   for the standard double precision. */

void dobldobl_print_solutions ( void );
/* 
 * DESCRIPTION :
 *   prints the solutions in the container,
 *   for the double double precision. */

void quaddobl_print_solutions ( void );
/* 
 * DESCRIPTION :
 *   prints the solutions in the container,
 *   for the quad double precision. */

int prompt_for_precision ( void );
/*
 * DESCRIPTION :
 *   Shows the user the available precision level
 *   and prompts the user to make a choice.
 *   The number on return is either 
 *     0 for standard double precision,
 *     1 for double double precision, or
 *     2 for quad double precision. */

int main ( int argc, char *argv[] )
{
   int myid,numprocs,len,n,choice;
   MPI_Status status; 

   adainit();
   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
   MPI_Comm_rank(MPI_COMM_WORLD,&myid);

   if(myid == 0) choice = prompt_for_precision();
   /* every node must know the precision level */
   MPI_Bcast(&choice,1,MPI_INT,0,MPI_COMM_WORLD);  

   if(choice == 0)
   {
      standard_dimension_broadcast(myid,&len,&n);
      standard_solutions_broadcast(myid,len,n);
      standard_print_broadcast(myid);
   }
   else if(choice == 1)
   {
      dobldobl_dimension_broadcast(myid,&len,&n);
      dobldobl_solutions_broadcast(myid,len,n);
      dobldobl_print_broadcast(myid);
   }
   else if(choice == 2)
   {
      quaddobl_dimension_broadcast(myid,&len,&n);
      quaddobl_solutions_broadcast(myid,len,n);
      quaddobl_print_broadcast(myid);
   }
   else
      if(myid == 0) printf("Invalid selection.\n");

   MPI_Barrier(MPI_COMM_WORLD);

   MPI_Finalize();
   adafinal();

   return 0;
}

int prompt_for_precision ( void )
{
   char newline;
   int choice;

   printf("\nMENU for selecting the precision level :\n");
   printf("  0. broadcast system in standard double precision;\n");
   printf("  1. broadcast system in double double precision;\n");
   printf("  2. broadcast system in quad double precision.\n");
   printf("Type 0, 1, or 2 to select precision level : ");
   scanf("%d",&choice);
   scanf("%c",&newline);

   return choice;
}

void standard_dimension_broadcast ( int myid, int *len, int *n )
{
   int k,fail;
   double *c;

   if(myid == 0)
   {
      printf("\n");
      fail = _ada_use_solcon(0,&k,&k,c);   /* prompts user for solutions */
      /* printf("The solution list in the container :\n"); */
      /* prints the solution list on screen */
      /* fail = _ada_use_solcon(1,&k,&k,c); */
      fail = _ada_use_solcon(2,&k,len,c);  /* get #solutions */
      fail = _ada_use_solcon(3,&k,n,c);    /* get dimension */
      printf("Root read %d solutions of dimension %d.\n",*len,*n); 
   }
   MPI_Bcast(n,1,MPI_INT,0,MPI_COMM_WORLD);  
   MPI_Bcast(len,1,MPI_INT,0,MPI_COMM_WORLD);  

   if(myid != 0)
      printf("Node %d expects %d solutions of dimension %d\n",myid,*len,*n);
}

void dobldobl_dimension_broadcast ( int myid, int *len, int *n )
{
   int k,fail;
   double *c;

   if(myid == 0)
   {
      printf("\n");
      fail = _ada_use_solcon(40,&k,&k,c);   /* prompts user for solutions */
      /* printf("The solution list in the container :\n"); */
      /* prints the solution list on screen */
      /* fail = _ada_use_solcon(41,&k,&k,c); */
      fail = _ada_use_solcon(42,&k,len,c);  /* get #solutions */
      fail = _ada_use_solcon(43,&k,n,c);    /* get dimension */
      printf("Root read %d solutions of dimension %d.\n",*len,*n); 
   }
   MPI_Bcast(n,1,MPI_INT,0,MPI_COMM_WORLD);  
   MPI_Bcast(len,1,MPI_INT,0,MPI_COMM_WORLD);  

   if(myid != 0)
      printf("Node %d expects %d solutions of dimension %d\n",myid,*len,*n);
}

void quaddobl_dimension_broadcast ( int myid, int *len, int *n )
{
   int k,fail;
   double *c;

   if(myid == 0)
   {
      printf("\n");
      fail = _ada_use_solcon(80,&k,&k,c);   /* prompts user for solutions */
      /* printf("The solution list in the container :\n"); */
      /* prints the solution list on screen */
      /* fail = _ada_use_solcon(81,&k,&k,c); */
      fail = _ada_use_solcon(82,&k,len,c);  /* get #solutions */
      fail = _ada_use_solcon(83,&k,n,c);    /* get dimension */
      printf("Root read %d solutions of dimension %d.\n",*len,*n); 
   }
   MPI_Bcast(n,1,MPI_INT,0,MPI_COMM_WORLD);  
   MPI_Bcast(len,1,MPI_INT,0,MPI_COMM_WORLD);  

   if(myid != 0)
      printf("Node %d expects %d solutions of dimension %d\n",myid,*len,*n);
}

void standard_solutions_broadcast ( int myid, int len, int n )
{
   int i,k,m[2],fail;
   double sol[2*n+5];

   m[0] = n;
   for(k=1; k<=len; k++)
   {
      if(myid == 0) fail = _ada_use_solcon(4,&k,&m[1],sol);
      MPI_Bcast(&m[1],1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(sol,2*n+5,MPI_DOUBLE,0,MPI_COMM_WORLD);
      if(myid != 0) fail = _ada_use_solcon(6,&k,m,sol);
   }
}

void dobldobl_solutions_broadcast ( int myid, int len, int n )
{
   int i,k,m[2],fail;
   double sol[4*n+10];

   m[0] = n;
   for(k=1; k<=len; k++)
   {
      if(myid == 0) fail = _ada_use_solcon(44,&k,&m[1],sol);
      MPI_Bcast(&m[1],1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(sol,4*n+10,MPI_DOUBLE,0,MPI_COMM_WORLD);
      if(myid != 0) fail = _ada_use_solcon(46,&k,m,sol);
   }
}

void quaddobl_solutions_broadcast ( int myid, int len, int n )
{
   int i,k,m[2],fail;
   double sol[8*n+20];

   m[0] = n;
   for(k=1; k<=len; k++)
   {
      if(myid == 0) fail = _ada_use_solcon(84,&k,&m[1],sol);
      MPI_Bcast(&m[1],1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(sol,8*n+20,MPI_DOUBLE,0,MPI_COMM_WORLD);
      if(myid != 0) fail = _ada_use_solcon(86,&k,m,sol);
   }
}

void standard_print_broadcast ( int myid )
{
   if(myid == 1)
   {
      printf("The solutions at node %d:\n",myid);
      standard_print_solutions();
   }
}

void dobldobl_print_broadcast ( int myid )
{
   if(myid == 1)
   {
      printf("The solutions at node %d:\n",myid);
      dobldobl_print_solutions();
   }
}

void quaddobl_print_broadcast ( int myid )
{
   if(myid == 1)
   {
      printf("The solutions at node %d:\n",myid);
      quaddobl_print_solutions();
   }
}

void standard_print_solution
 ( int n, double rt, double it, int m, double *v,
   double err, double rco, double res )
{
   int i;

   printf("t : %.15e  %.15e\n",rt,it);
   printf("m : %d\n", m);
   printf("The solution for t :\n");
   for(i=0; i<n; i++)
      printf(" %.15e  %.15e\n", v[2*i], v[2*i+1]);
   printf("== err : %.3e = rco : %.3e = res : %.3e ==\n",err,rco,res);
}

void dobldobl_print_solution
 ( int n, double *rt, double *it, int m, double *v,
   double err, double rco, double res )
{
   int i;

   printf("t : %.15e  %.15e  %.15e  %.15e\n",rt[0],rt[1],it[0],it[1]);
   printf("m : %d\n", m);
   printf("The solution for t :\n");
   for(i=0; i<n; i++)
      printf(" %.15e  %.15e  %.15e  %.15e\n",
             v[4*i],v[4*i+1],v[4*i+2],v[4*i+3]);
   printf("== err : %.3e = rco : %.3e = res : %.3e ==\n",err,rco,res);
}

void quaddobl_print_solution
 ( int n, double *rt, double *it, int m, double *v,
   double err, double rco, double res )
{
   int i;

   printf("t : %.15e  %.15e  %.15e  %.15e  %.15e  %.15e  %.15e  %.15e\n",
          rt[0],rt[1],rt[2],rt[3],it[0],it[1],it[2],it[3]);
   printf("m : %d\n", m);
   printf("The solution for t :\n");
   for(i=0; i<n; i++)
      printf(" %.15e  %.15e  %.15e  %.15e  %.15e  %.15e  %.15e  %.15e\n",
             v[4*i],v[4*i+1],v[4*i+2],v[4*i+3],
             v[4*i+4],v[4*i+5],v[4*i+6],v[4*i+7]);
   printf("== err : %.3e = rco : %.3e = res : %.3e ==\n",err,rco,res);
}

void standard_print_solutions ( void ) 
{
   int k,len,n,m,fail;
   double *c;

   fail = _ada_use_solcon(2,&k,&len,c);
   fail = _ada_use_solcon(3,&k,&n,c);

   c = (double*)calloc(2*n+5,sizeof(double));

   for(k=1; k<=len; k++)
   {
      fail = _ada_use_solcon(4,&k,&m,c);
      printf("Solution %d : \n", k);
      standard_print_solution(n,c[0],c[1],m,c+2,c[2*n+2],c[2*n+3],c[2*n+4]);
   }
}

void dobldobl_print_solutions ( void ) 
{
   int k,len,n,m,fail;
   double *c;

   fail = _ada_use_solcon(42,&k,&len,c);
   fail = _ada_use_solcon(43,&k,&n,c);

   c = (double*)calloc(4*n+10,sizeof(double));

   for(k=1; k<=len; k++)
   {
      fail = _ada_use_solcon(44,&k,&m,c);
      printf("Solution %d : \n", k);
      dobldobl_print_solution(n,c,c+2,m,c+4,c[4*n+4],c[4*n+6],c[4*n+8]);
   }
}

void quaddobl_print_solutions ( void ) 
{
   int k,len,n,m,fail;
   double *c;

   fail = _ada_use_solcon(82,&k,&len,c);
   fail = _ada_use_solcon(83,&k,&n,c);

   c = (double*)calloc(8*n+20,sizeof(double));

   for(k=1; k<=len; k++)
   {
      fail = _ada_use_solcon(84,&k,&m,c);
      printf("Solution %d : \n", k);
      quaddobl_print_solution(n,c,c+4,m,c+8,c[8*n+8],c[4*n+12],c[8*n+16]);
   }
}
