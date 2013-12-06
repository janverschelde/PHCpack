/* The following is a simple interface of MPI with PHCpack:
 * after prompting the user for a target and start system,
 * the homotopy is then broadcasted to all the nodes. */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "mpi.h"

#define SEND_SOL 100 /* message tag for sending start solution */
#define SEND_MUL 101 /* message tag for sending multiplicity */

extern void adainit( void );
extern int _ada_use_c2phc ( int task, int *a, int *b, double *c );
extern void adafinal( void );

void dimension_broadcast ( int myid, int *n );
/* gets the dimension n from the root node,
 * broadcasts the dimension then to all the nodes, and
 * initializes the systems container at each node */

void monomials_broadcast ( int myid, int n );
/* broadcasts all monomials in the system to the nodes,
 * adding the monomials to the containers at each node */

void solutions_distribute ( int myid, int nbsols, int n, int numprocs, int *mysolnum );
/* distributes all solutions in the container at the root
 * to the solutions container at each node */

void copy_broadcast ( int myid );
/* copies the result of the broadcast to PHCpack data and
 * prints the result of the broadcast at the node */

void start_system_broadcast ( int myid, int n, int *nbsols );
/* broadcast of start system and start solutions to nodes,
 * n is the dimension and nbsols the number of solutions */

void print_monomials ( void );
/* prints the monomials in the container */

void print_system ( void );
/* prints the polynomial system in the container */

void print_homotopy ( void );
/* prints the target and start system in the homotopy */

int main ( int argc, char *argv[] )
{
   int myid,numprocs,n,nbsols,mysolnum;
   MPI_Status status; 

   adainit();
   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
   MPI_Comm_rank(MPI_COMM_WORLD,&myid);

   dimension_broadcast(myid,&n);
   monomials_broadcast(myid,n);
   copy_broadcast(myid); 
   start_system_broadcast(myid,n,&nbsols);

   MPI_Bcast(&nbsols,1,MPI_INT,0,MPI_COMM_WORLD);
   solutions_distribute(myid,nbsols,n,numprocs,&mysolnum);

  // if(myid == 1) print_homotopy();  

   MPI_Finalize();
   adafinal();

   return 0;
}

void dimension_broadcast ( int myid, int *n )
{
   int *d,fail,m;
   double *c;

   if(myid == 0)
   {
      fail = _ada_use_c2phc(11,n,d,c);    /* read target system */
     /* fail = _ada_use_c2phc(12,n,d,c); */   /* write target system */
      fail = _ada_use_c2phc(1,n,d,c);     /* copy target to container */
     /* fail = _ada_use_c2phc(21,n,d,c); */   /* write system in container */
      fail = _ada_use_c2phc(22,n,d,c);    /* get dimension */
     /* printf("The dimension is %d.\n", *n); */
   }
 
   MPI_Bcast(n, 1, MPI_INT, 0, MPI_COMM_WORLD);  

   if(myid != 0)
   {
     /* printf("Node %d knowns that the dimension is %d\n", myid, *n); */
      fail = _ada_use_c2phc(23,n,d,c);    /* initialize container */
     /* fail = _ada_use_c2phc(22,&m,d,c); */  /* get dimension as test */
     /* printf("  and initialized container with dimension %d.\n", m); */
   }
}

void monomials_broadcast ( int myid, int n )
{
   double cff[2];
   int i,j,exp[n],m[3],mm,fail;

   for(i=1; i<=n; i++)
   {
      m[1] = i;
      if(myid == 0)
      {
         fail = _ada_use_c2phc(24,m,exp,cff);    /* get #monomials */
         mm = m[0];
        /* printf("Polynomial %d has %d monomials.\n",i,mm); */
      }

      MPI_Bcast(&mm, 1, MPI_INT, 0, MPI_COMM_WORLD);

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
}

void solutions_distribute ( int myid, int nbsols, int n, int numprocs, int *mysolnum )
{
   int i,k,m[2],fail,dest, *a, *b;
   double sol[2*n+5], *c;
   MPI_Status status;

   m[0] = n;
   if(myid==0)
   {
     for(k=1; k<=nbsols; k++)
     {
       dest = k%(numprocs-1);
       if(dest==0)  dest=numprocs-1;    
       fail = _ada_use_c2phc(34,&k,&m[1],sol);
       MPI_Send(&m[1],1,MPI_INT,dest,SEND_MUL,MPI_COMM_WORLD);
       MPI_Send(sol,2*n+5,MPI_DOUBLE,dest,SEND_SOL,MPI_COMM_WORLD);
printf("%dth root is sending to %d\n", k, dest);
     }
   }
   else
   {
     *mysolnum=nbsols/(numprocs-1);
     if(myid<=(nbsols%(numprocs-1)))
       *mysolnum++;
     for(k=1; k<=*mysolnum; k++)
     {
       MPI_Recv(&m[1],1,MPI_INT,0,SEND_MUL,MPI_COMM_WORLD,&status);
       MPI_Recv(sol,2*n+5,MPI_DOUBLE,0,SEND_SOL,MPI_COMM_WORLD,&status); 
       fail = _ada_use_c2phc(36,&k,m,sol);
printf("%d is receiving %dth root.\n", myid,k);
     }
   }

   if(myid !=0 )
   {
      fail = _ada_use_c2phc(8,a,b,c);     /* copy container to start sols */
   }
   fail = _ada_use_c2phc(37,a,b,c);       /* clear solution container */
}

void copy_broadcast ( int myid )
{
   int n,*a,*b,fail;
   double *c;

/*   if(myid == 1)
   {
      printf("The monomials at node %d:\n", myid);
      print_monomials();
      printf("The polynomial system at node %d:\n", myid);
      print_system();
   }
*/
   if(myid != 0)
   {
      fail = _ada_use_c2phc(2,a,b,c);     /* copy container to target */
      fail = _ada_use_c2phc(27,a,b,c);    /* clear systems container */
     /* fail = _ada_use_c2phc(12,a,b,c); */  /* write target system */
   }
}

void start_system_broadcast ( int myid, int n, int *nbsols )
{
   int *a,*b,fail;
   double *c;

   if(myid == 0)
   {
      fail = _ada_use_c2phc(13,a,b,c);    /* read start system */
     /* fail = _ada_use_c2phc(14,a,b,c); */  /* write start system */
      fail = _ada_use_c2phc(3,a,b,c);     /* copy start to container */
     /* fail = _ada_use_c2phc(21,a,b,c); */  /* write system in container */
      fail = _ada_use_c2phc(7,a,b,c);     /* copy start sols to container */
     /* fail = _ada_use_c2phc(15,a,b,c); */  /* write start solutions */
   }
   if(myid != 0)
      fail = _ada_use_c2phc(23,&n,b,c);   /* initialize container */

   monomials_broadcast(myid,n);           /* broadcast container */

   if(myid != 0)                          /* copy result of broadcast */
   {
      fail = _ada_use_c2phc(4,a,b,c);     /* copy container to target */
      fail = _ada_use_c2phc(27,a,b,c);    /* clear systems container */
     /* fail = _ada_use_c2phc(14,a,b,c); */  /* write start system */
   }

   if(myid == 0)
   {
      fail = _ada_use_c2phc(32,a,nbsols,c);  /* get #solutions */
     /* printf("There are %d start solutions.\n", *nbsols); */
   }
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

void print_homotopy ( void )
{
   int *a,*b,fail;
   double *c;

   printf("\n");
   fail = _ada_use_c2phc(12,a,b,c);
   fail = _ada_use_c2phc(14,a,b,c);
   printf("The start solutions :\n");
   fail = _ada_use_c2phc(15,a,b,c);
}
