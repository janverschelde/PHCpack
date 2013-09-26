/* The following program tries to send an integer vector around... */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "mpi.h"
#include "call_phc_rw.h"

extern void adainit();
extern void adafinal();

#define SEND 100

void print_vector ( int n, double a[n] );
/* writes the vector to screen */

int main ( int argc, char *argv[] )
{
    int myid,n,i,numprocs;
    double *x;

    MPI_Status status; 

    adainit();

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);

    if(myid==0)
    {
        srand( (unsigned)time( NULL ) );
        printf("Please give the number of doubles : ");
        scanf("%d", &n);
    }

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);    
    /* printf("\nallocating memory for %d doubles...\n", n); */
    x = (double*) calloc(n, sizeof(double)); 

    if(myid==0)
    {
        call_phc_rw(n,0,x);
    }

    MPI_Bcast(x, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if(myid==0)
    {
        for(i=1; i<numprocs; i++)
        {
           MPI_Recv(x, n, MPI_DOUBLE, i, SEND, MPI_COMM_WORLD, &status);
           printf("received from node %d :\n", i);
           /* print_vector(n,x); */
           call_phc_rw(n,1,x);
        }
    }
    else
    {
        printf("Node %d received vector\n", myid);
        call_phc_rw(n,1,x);
        printf("Sending vector... ");
        /* print_vector(n,x); */
        MPI_Send(x, n, MPI_DOUBLE, 0, SEND, MPI_COMM_WORLD);
    }

    free(x);

    MPI_Finalize();

    adafinal();

    return 0;
}

void print_vector ( int n, double a[n] )
{
   int i;

   for(i=0; i<n; i++)
      printf(" %2.1f", a[i]);
   printf("\n");
}
