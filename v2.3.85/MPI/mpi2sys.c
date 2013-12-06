/* The following is a simple interface of MPI with PHCpack:
 * after prompting the user for a polynomial system,
 * the system is broadcasted to all the nodes. */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "mpi.h"

extern void adainit( void );
extern int _ada_use_syscon ( int task, int *a, int *b, double *c );
extern void adafinal( void );

void standard_dimension_broadcast ( int myid, int *n );
/*
 * DESCRIPTION :
 *   Gets the dimension n from the root node,
 *   broadcasts the dimension then to all the nodes, and
 *   initializes the systems container at each node,
 *   for the standard double precision. */

void dobldobl_dimension_broadcast ( int myid, int *n );
/*
 * DESCRIPTION :
 *   Gets the dimension n from the root node,
 *   broadcasts the dimension then to all the nodes, and
 *   initializes the systems container at each node,
 *   for the double double precision. */

void quaddobl_dimension_broadcast ( int myid, int *n );
/*
 * DESCRIPTION :
 *   Gets the dimension n from the root node,
 *   broadcasts the dimension then to all the nodes, and
 *   initializes the systems container at each node,
 *   for the quad double precision. */

void standard_monomials_broadcast ( int myid, int n );
/* 
 * DESCRIPTION :
 *   Broadcasts all monomials in the system to the nodes,
 *   adding the monomials to the containers at each node,
 *   for the standard double precision. */

void dobldobl_monomials_broadcast ( int myid, int n );
/* 
 * DESCRIPTION :
 *   Broadcasts all monomials in the system to the nodes,
 *   adding the monomials to the containers at each node,
 *   for the double double precision. */

void quaddobl_monomials_broadcast ( int myid, int n );
/* 
 * DESCRIPTION :
 *   Broadcasts all monomials in the system to the nodes,
 *   adding the monomials to the containers at each node,
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
 *   for the quad double precision. */

void standard_print_monomials ( void );
/* 
 * DESCRIPTION :
 *   prints the monomials in the container,
 *   for the standard double precision. */

void dobldobl_print_monomials ( void );
/* 
 * DESCRIPTION :
 *   prints the monomials in the container,
 *   for the double double precision. */

void quaddobl_print_monomials ( void );
/* 
 * DESCRIPTION :
 *   prints the monomials in the container,
 *   for the quad double precision. */

void standard_print_system ( void );
/*
 * DESCRIPTION :
 *   prints the polynomial system in the container
 *   for the standard double precision. */

void dobldobl_print_system ( void );
/*
 * DESCRIPTION :
 *   prints the polynomial system in the container
 *   for the double double precision. */

void quaddobl_print_system ( void );
/*
 * DESCRIPTION :
 *   prints the polynomial system in the container
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
   int myid,numprocs,n,choice;
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
      standard_dimension_broadcast(myid,&n);
      standard_monomials_broadcast(myid,n);
      standard_print_broadcast(myid); 
   }
   else if(choice == 1)
   {
      dobldobl_dimension_broadcast(myid,&n);
      dobldobl_monomials_broadcast(myid,n);
      dobldobl_print_broadcast(myid); 
   }
   else if(choice == 2)
   {
      quaddobl_dimension_broadcast(myid,&n);
      quaddobl_monomials_broadcast(myid,n);
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
   /* scanf("%c",&newline); */

   return choice;
}

void standard_dimension_broadcast ( int myid, int *n )
{
   int *d,fail,m;
   double *c;

   if(myid == 0)
   {
      fail = _ada_use_syscon(0,n,d,c);    /* read system */
      fail = _ada_use_syscon(1,n,d,c);    /* write system */
      fail = _ada_use_syscon(2,n,d,c);    /* get dimension */
      printf("The dimension is %d.\n", *n); 
   }
 
   MPI_Bcast(n,1,MPI_INT,0,MPI_COMM_WORLD);  

   if(myid != 0)
   {
      printf("Node %d knowns that the dimension is %d\n", myid, *n);
      fail = _ada_use_syscon(3,n,d,c);    /* initialize container */
      fail = _ada_use_syscon(2,&m,d,c);   /* get dimension as test */
      printf("  and initialized container with dimension %d.\n", m);
   }
}

void dobldobl_dimension_broadcast ( int myid, int *n )
{
   int *d,fail,m;
   double *c;

   if(myid == 0)
   {
      fail = _ada_use_syscon(200,n,d,c);    /* read system */
      fail = _ada_use_syscon(201,n,d,c);    /* write system */
      fail = _ada_use_syscon(202,n,d,c);    /* get dimension */
      printf("The dimension is %d.\n", *n); 
   }
 
   MPI_Bcast(n,1,MPI_INT,0,MPI_COMM_WORLD);  

   if(myid != 0)
   {
      printf("Node %d knowns that the dimension is %d\n", myid, *n);
      fail = _ada_use_syscon(203,n,d,c);    /* initialize container */
      fail = _ada_use_syscon(202,&m,d,c);   /* get dimension as test */
      printf("  and initialized container with dimension %d.\n", m);
   }
}

void quaddobl_dimension_broadcast ( int myid, int *n )
{
   int *d,fail,m;
   double *c;

   if(myid == 0)
   {
      fail = _ada_use_syscon(210,n,d,c);    /* read system */
      fail = _ada_use_syscon(211,n,d,c);    /* write system */
      fail = _ada_use_syscon(212,n,d,c);    /* get dimension */
      printf("The dimension is %d.\n", *n); 
   }
 
   MPI_Bcast(n,1,MPI_INT,0,MPI_COMM_WORLD);  

   if(myid != 0)
   {
      printf("Node %d knowns that the dimension is %d\n", myid, *n);
      fail = _ada_use_syscon(213,n,d,c);    /* initialize container */
      fail = _ada_use_syscon(212,&m,d,c);   /* get dimension as test */
      printf("  and initialized container with dimension %d.\n", m);
   }
}

void standard_monomials_broadcast ( int myid, int n )
{
   double cff[2];
   int i,j,exp[n],m[3],mm,fail;

   for(i=1; i<=n; i++)
   {
      m[1] = i;
      if(myid == 0)
      {
         fail = _ada_use_syscon(4,m,exp,cff);    /* get #monomials */
         mm = m[0];
         printf("Polynomial %d has %d monomials.\n",i,mm);
      }
      MPI_Bcast(&mm,1,MPI_INT,0,MPI_COMM_WORLD);

      m[0] = n;
      for(j=1; j<=mm; j++)    /* broadcast j-th term of i-th polynomial */
      {
         m[2] = j;
         if(myid == 0) fail = _ada_use_syscon(5,m,exp,cff); 
         MPI_Bcast(cff,2,MPI_DOUBLE,0,MPI_COMM_WORLD); 
         MPI_Bcast(exp,n,MPI_INT,0,MPI_COMM_WORLD); 
         if(myid != 0) fail = _ada_use_syscon(6,m,exp,cff); 
      }
   }
}

void dobldobl_monomials_broadcast ( int myid, int n )
{
   double cff[4];
   int i,j,exp[n],m[3],mm,fail;

   for(i=1; i<=n; i++)
   {
      m[1] = i;
      if(myid == 0)
      {
         fail = _ada_use_syscon(204,m,exp,cff);    /* get #monomials */
         mm = m[0];
         printf("Polynomial %d has %d monomials.\n",i,mm);
      }
      MPI_Bcast(&mm,1,MPI_INT,0,MPI_COMM_WORLD);

      m[0] = n;
      for(j=1; j<=mm; j++)    /* broadcast j-th term of i-th polynomial */
      {
         m[2] = j;
         if(myid == 0) fail = _ada_use_syscon(205,m,exp,cff); 
         MPI_Bcast(cff,4,MPI_DOUBLE,0,MPI_COMM_WORLD); 
         MPI_Bcast(exp,n,MPI_INT,0,MPI_COMM_WORLD); 
         if(myid != 0) fail = _ada_use_syscon(206,m,exp,cff); 
      }
   }
}

void quaddobl_monomials_broadcast ( int myid, int n )
{
   double cff[8];
   int i,j,exp[n],m[3],mm,fail;

   for(i=1; i<=n; i++)
   {
      m[1] = i;
      if(myid == 0)
      {
         fail = _ada_use_syscon(214,m,exp,cff);    /* get #monomials */
         mm = m[0];
         printf("Polynomial %d has %d monomials.\n",i,mm);
      }
      MPI_Bcast(&mm,1,MPI_INT,0,MPI_COMM_WORLD);

      m[0] = n;
      for(j=1; j<=mm; j++)    /* broadcast j-th term of i-th polynomial */
      {
         m[2] = j;
         if(myid == 0) fail = _ada_use_syscon(215,m,exp,cff); 
         MPI_Bcast(cff,8,MPI_DOUBLE,0,MPI_COMM_WORLD); 
         MPI_Bcast(exp,n,MPI_INT,0,MPI_COMM_WORLD); 
         if(myid != 0) fail = _ada_use_syscon(216,m,exp,cff); 
      }
   }
}

void standard_print_broadcast ( int myid )
{
   if(myid == 1)
   {
      printf("The monomials at node %d:\n", myid);
      standard_print_monomials();
      printf("The polynomial system at node %d:\n", myid);
      standard_print_system();
   }
}

void dobldobl_print_broadcast ( int myid )
{
   if(myid == 1)
   {
      printf("The monomials at node %d:\n", myid);
      dobldobl_print_monomials();
      printf("The polynomial system at node %d:\n", myid);
      dobldobl_print_system();
   }
}

void quaddobl_print_broadcast ( int myid )
{
   if(myid == 1)
   {
      printf("The monomials at node %d:\n", myid);
      quaddobl_print_monomials();
      printf("The polynomial system at node %d:\n", myid);
      quaddobl_print_system();
   }
}

void standard_print_monomials ( void ) 
{
   int *d,i,j,k,n,m[3],mm,fail;
   double c[2];

   fail = _ada_use_syscon(2,&n,d,c);

   d = (int*)calloc(n,sizeof(int));

   for(i=1; i<=n; i++)
   {
      m[1] = i;
      fail = _ada_use_syscon(4,m,d,c);
      mm = m[0];
      m[0] = n;
      printf("Polynomial %d has %d monomials :\n",i,mm);
      for(j=1; j<=mm; j++)
      {
         m[2] = j;
         fail = _ada_use_syscon(5,m,d,c);
         printf(" %.15e  %.15e",c[0],c[1]);
         for(k=0; k<n; k++) printf(" %d",d[k]);
         printf("\n");
      }
   }
}

void dobldobl_print_monomials ( void ) 
{
   int *d,i,j,k,n,m[3],mm,fail;
   double c[4];

   fail = _ada_use_syscon(202,&n,d,c);

   d = (int*)calloc(n,sizeof(int));

   for(i=1; i<=n; i++)
   {
      m[1] = i;
      fail = _ada_use_syscon(204,m,d,c);
      mm = m[0];
      m[0] = n;
      printf("Polynomial %d has %d monomials :\n",i,mm);
      for(j=1; j<=mm; j++)
      {
         m[2] = j;
         fail = _ada_use_syscon(205,m,d,c);
         printf(" %.15e  %.15e  %.15e  %.15e",c[0],c[1],c[2],c[3]);
         for(k=0; k<n; k++) printf(" %d",d[k]);
         printf("\n");
      }
   }
}

void quaddobl_print_monomials ( void ) 
{
   int *d,i,j,k,n,m[3],mm,fail;
   double c[8];

   fail = _ada_use_syscon(212,&n,d,c);

   d = (int*)calloc(n,sizeof(int));

   for(i=1; i<=n; i++)
   {
      m[1] = i;
      fail = _ada_use_syscon(214,m,d,c);
      mm = m[0];
      m[0] = n;
      printf("Polynomial %d has %d monomials :\n",i,mm);
      for(j=1; j<=mm; j++)
      {
         m[2] = j;
         fail = _ada_use_syscon(215,m,d,c);
         printf(" %.15e  %.15e  %.15e  %.15e  %.15e  %.15e  %.15e  %.15e",
                c[0],c[1],c[2],c[3],c[4],c[5],c[6],c[7]);
         for(k=0; k<n; k++) printf(" %d",d[k]);
         printf("\n");
      }
   }
}

void standard_print_system ( void ) 
{
   int *d,*n,fail;
   double *c;

   fail = _ada_use_syscon(1,n,d,c);
}

void dobldobl_print_system ( void ) 
{
   int *d,*n,fail;
   double *c;

   fail = _ada_use_syscon(201,n,d,c);
}

void quaddobl_print_system ( void ) 
{
   int *d,*n,fail;
   double *c;

   fail = _ada_use_syscon(211,n,d,c);
}
