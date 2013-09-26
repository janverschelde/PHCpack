/* The following is a simple interface of MPI with PHCpack:
 * after prompting the user for a target and start system,
 * the start roots are distributed to all the nodes
 * dynamically, overlapping computation and communication
 * to improve the efficiency   */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "mpi.h"

#define SEND_SSOL 100    /* tag for sending start solution */
#define SEND_SMUL 101    /* tag for sending multiplicity of start solution */
#define SEND_TSOL 102    /* tag for sending target solution */
#define SEND_TMUL 103    /* tag for sending multiplicity of target solution */
#define SEND_FIN_M 104   /* tag for sending multiplicity of last solution */
#define SEND_FIN_SOL 105 /* tag for sending the final solution of each node */

#define v 0  /* verbose flag: 0 is silent, 1 is verbose */

extern void adainit( void );
extern int _ada_use_c2phc ( int task, int *a, int *b, double *c );
extern void adafinal( void );

void dimension_broadcast ( int myid, int *n );
/*
 * DESCRIPTION :
 *   Every node receives the dimension n from the root node,
 *   and initializes the systems container at each node. */

void monomials_broadcast ( int myid, int n );
/*
 * DESCRIPTION :
 *   Broadcasts all monomials in the system to the nodes,
 *   adding the monomials to the containers at each node. */

void track_path ( int m[2], double *sol );
/* 
 * DESCRIPTION :
 *   Tracks a path starting at the given solution. */

void dynamic_load
 ( int myid, FILE *fp, int szsymb, char *symbols,
   int n, int numprocs, int nbsols, int *mysolnb );
/* 
 * DESCRIPTION :
 *   The root node distributes one start solution to other nodes each time.
 *   After finishing path tracking, the nodes return their target
 *   solutions to the root node.  The root node prints out the solution 
 *   to the file with point fp and distributes a new job to the node 
 *   which finishes its job.
 *
 * ON ENTRY :
 *   myid       rank of the MPI process;
 *   fp         pointer to a file which must be opened for output;
 *   szsymb     number of characters in the string symbols;
 *   symbols    string of length nc which contains the symbols of 
 *              the variables, each variable is separated by a space;
 *   n          dimension of the system;
 *   numprocs   number of processes;
 *   nbsols     number of solution paths to be tracked;
 *   mysolnb    array of size numprocs.
 *
 * ON RETURN :
 *   mysolnb    entry k contains the number of paths tracked
 *              by process with rank k. */

void print_solution
 ( FILE *fp, int szsymb, char *symbols, int n, int m, double *sol );
/*
 * DESCRIPTION :
 *   The root node prints out the target solutions to file.
 *   The information of the symbol table is in the string symbols
 *   of size szsymb. */

void copy_broadcast ( int myid );
/*
 * DESCRIPTION :
 *   Copies the result of the broadcast to PHCpack data and
 *   prints the result of the broadcast at the node. */

void start_system_broadcast ( int myid, int n, int *nbsols );
/*
 * DESCRIPTION :
 *   broadcast of start system and start solutions to nodes,
 *   n is the dimension and nbsols the number of solutions. */

void print_monomials ( void );
/*
 * DESCRIPTION :
 *   Prints the monomials in the container. */

/* void print_system ( void ); */
/* prints the polynomial system in the container */

void print_homotopy ( void );
/*
 * DESCRIPTION :
 *   Prints the target and start system in the homotopy. */

void print_time ( FILE *fp, double *time, int numprocs, int *mysolnb );
/*
 * DESCRIPTION :
 *   Prints to file with point fp the time spent and the number of paths
 *   tracked on each of the processors. */

int main ( int argc, char *argv[] )
{
   int myid,numprocs,n,nbsols,mysolnb,*mysol;
   double startwtime,endwtime,wtime,*time;
   char name[80];
   FILE *fp;
   int size_of_symbols = 1024;     /* total size of all symbols */
   int isymbols[size_of_symbols];
   char symbols[size_of_symbols];  /* symbols for the variables */
   MPI_Status status;
 
   adainit();
   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
   MPI_Comm_rank(MPI_COMM_WORLD,&myid);

   if(myid == 0)
   {
      char nl;

      time = (double*) calloc(numprocs, sizeof(double));
      mysol = (int*) calloc(numprocs, sizeof(int));
      startwtime = MPI_Wtime();

      printf("\nGive the name of the output file : ");
      scanf("%s",name); scanf("%c",&nl); /* skip newline symbol */
      if(v) printf("opening %s for writing ...\n",name);
      fp = fopen(name,"w");
   }
   dimension_broadcast(myid,&n);
   if(v) if(myid == 0) printf("\nbroadcasting the monomials ...\n");
   monomials_broadcast(myid,n);
   copy_broadcast(myid); 
   if(v) if(myid == 0) printf("\nbroadcasting the start system ...\n");
   start_system_broadcast(myid,n,&nbsols);

   if(myid == 0)  /* reading the symbol table into a string */
   {
      int fail,k;
      double *c;

      fail = _ada_use_c2phc(295,&size_of_symbols,isymbols,c);
      for(k=0; k<size_of_symbols; k++) symbols[k] = (char) isymbols[k];
      symbols[size_of_symbols] = '\0';
      if(v) printf("the symbols : %s\n",symbols);
   }
   MPI_Bcast(&nbsols,1,MPI_INT,0,MPI_COMM_WORLD);
   
   if(myid != 0) startwtime = MPI_Wtime();
   if(myid == 0) if(v) printf("\nLaunching path tracking jobs ...\n\n");
   dynamic_load(myid,fp,size_of_symbols,symbols,n,numprocs,nbsols,&mysolnb);
   if(myid != 0) endwtime = MPI_Wtime();

   MPI_Barrier(MPI_COMM_WORLD);
   if(myid == 0) endwtime = MPI_Wtime();
   wtime = endwtime-startwtime;
   MPI_Gather(&wtime,1,MPI_DOUBLE,time,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Gather(&mysolnb,1,MPI_INT,mysol,1,MPI_INT,0,MPI_COMM_WORLD);
   if(myid == 0) print_time(fp,time,numprocs,mysol);
   if(myid == 0) fclose(fp);
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
   MPI_Bcast(n,1,MPI_INT,0,MPI_COMM_WORLD);  

   if(myid != 0)
   {
      if(v) printf("Node %d knows that the dimension is %d\n", myid, *n);
      fail = _ada_use_c2phc(23,n,d,c);    /* initialize container */
      if(v) fail = _ada_use_c2phc(22,&m,d,c);  /* get dimension as test */
      if(v) printf("  and initialized container with dimension %d.\n", m);
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
         mm = m[0]; /* printf("Polynomial %d has %d monomials.\n",i,mm); */
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

void dynamic_load
 ( int myid, FILE *fp, int szsymb, char *symbols,
   int n, int numprocs, int nbsols, int *mysolnb )
{
   int dest, m[2], send[2], i, k, fail, len=2*n+5, send_m;
   double sol[len],t_sol[len], send_sol[len];
   MPI_Status status, stat[4];
   MPI_Request request[4];

   m[0] = n;
   *mysolnb = 0;
   if(myid == 0)
   {
      /*distributes 2 start roots for each client at the beginning */
      for(k=1; k<=2*(numprocs-1); k++)
      {
         fail = _ada_use_c2phc(34,&k,&m[1],sol); /* retrieve solution */
         send[0] = k; send[1] = m[1]; 
         if(k<=numprocs-1)
            dest = k;
         else
            dest=k-numprocs+1;

         MPI_Isend(send,2,MPI_INT,dest,SEND_SMUL,MPI_COMM_WORLD,&(request[0]));
         MPI_Isend(sol,2*n+5,MPI_DOUBLE,dest,SEND_SSOL,
                   MPI_COMM_WORLD,&(request[1]));  
      }
      /* fprintf(fp,"THE SOLUTIONS :\n\n"); */
      fprintf(fp,"%d %d\n", nbsols, n);
      fprintf
         (fp,"===========================================================\n");

      /*collects the target roots and distributes the remaining start roots */ 
      for(k=2*numprocs-1; k<=nbsols+numprocs-1; k++)
      {
         MPI_Recv(&m[1],1,MPI_INT,MPI_ANY_SOURCE,SEND_TMUL,
                  MPI_COMM_WORLD,&status);
         dest = status.MPI_SOURCE;
         MPI_Recv(sol,2*n+5,MPI_DOUBLE,dest,SEND_TSOL,MPI_COMM_WORLD,&status);

         fprintf(fp,"solution %d :\n", k-2*numprocs+2);
         print_solution(fp,szsymb,symbols,n,m[1],sol); 

         fail = _ada_use_c2phc(34,&k,&m[1],sol); /* retrieve solution */
         send[0] = k; send[1] = m[1];
         MPI_Send(send,2,MPI_INT,dest,SEND_SMUL,MPI_COMM_WORLD);
         MPI_Send(sol,2*n+5,MPI_DOUBLE,dest,SEND_SSOL,MPI_COMM_WORLD);
         (*mysolnb)++;
      }
      for(k=nbsols-numprocs+2; k<=nbsols; k++)
      /* receives last root from each node */
      {
         MPI_Recv(&m[1],1,MPI_INT,MPI_ANY_SOURCE,SEND_FIN_M,
                  MPI_COMM_WORLD,&status);
         MPI_Recv(sol,2*n+5,MPI_DOUBLE,MPI_ANY_SOURCE,SEND_FIN_SOL,
                  MPI_COMM_WORLD,&status);
         fprintf(fp,"solution %d :\n", k);
         print_solution(fp,szsymb,symbols,n,m[1],sol);
         (*mysolnb)++;
      }
   }
   else
   {
      /*receives 2 jobs at the beginning, 1 will be sent back after finished,
        the other will be computed while communicating */
      MPI_Recv(send,2,MPI_INT,0,SEND_SMUL,MPI_COMM_WORLD,&status);
      m[1] = send[1];
      MPI_Recv(sol,2*n+5,MPI_DOUBLE,0,SEND_SSOL,MPI_COMM_WORLD,&status);
      track_path(m,sol);
      send_m = m[1];
      for(i=0; i<len; i++) send_sol[i] = sol[i];
      (*mysolnb)++;

      MPI_Recv(send,2,MPI_INT,0,SEND_SMUL,MPI_COMM_WORLD,&status);
      m[1] = send[1];
      MPI_Recv(sol,2*n+5,MPI_DOUBLE,0,SEND_SSOL,MPI_COMM_WORLD,&status);

      while(1)
      {
         (*mysolnb)++;
         MPI_Isend(&send_m,1,MPI_INT,0,SEND_TMUL,MPI_COMM_WORLD,&(request[0]));
         MPI_Isend(send_sol,2*n+5,MPI_DOUBLE,0,SEND_TSOL,
                   MPI_COMM_WORLD,&(request[1]));

         MPI_Irecv(send,2,MPI_INT,0,SEND_SMUL,MPI_COMM_WORLD,&(request[2]));
         MPI_Irecv(t_sol,2*n+5,MPI_DOUBLE,0,SEND_SSOL,
                   MPI_COMM_WORLD,&(request[3])); 
         track_path(m,sol);

         MPI_Waitall(4,request,stat); /*waits for a new job*/

         send_m = m[1];
         for(i=0; i<len; i++) send_sol[i] = sol[i];

         if(send[0]>nbsols) break;
         for(i=0; i<len; i++) sol[i] = t_sol[i];
         m[1] = send[1];
      }
      MPI_Send(&send_m,1,MPI_INT,0,SEND_FIN_M,MPI_COMM_WORLD);
      MPI_Send(send_sol,2*n+5,MPI_DOUBLE,0,SEND_FIN_SOL,MPI_COMM_WORLD);

      if(v) printf("Node %d done.\n", myid);
   }
}

void print_solution
 ( FILE *fp, int szsymb, char *symbols, int n, int m, double *sol )
{
   int i,k; 

   fprintf(fp,"t : %.15E  %.15E\n",sol[0],sol[1]);
   fprintf(fp,"m : %d\n",m);
   fprintf(fp,"The solution for t :\n");
   k = 0;
   while((k < szsymb) && (symbols[k] == ' ')) k++;
   for(i=0; i<n; i++)
   {
      fprintf(fp," ");
      while((k < szsymb) && (symbols[k] != ' ')) /* write symbol */
         fprintf(fp,"%c",symbols[k++]);
        
      fprintf(fp," : %.15E  %.15E\n",sol[2*i+2],sol[2*i+3]);
      if(k < szsymb) k++;  /* skip space between symbols */
   }
   fprintf(fp,"== err : %.3E = rco : %.3E = res : %.3E ==\n",
           sol[2*n+2],sol[2*n+3],sol[2*n+4]);
}

void track_path ( int m[2], double *sol )
{
   int *a,*b,fail,k=1;
   double *c;
 
   fail = _ada_use_c2phc(36,&k,m,sol); /* append a solution to the container */ 
   fail = _ada_use_c2phc(8,a,b,c);     /* copy container to start sols */
   fail = _ada_use_c2phc(37,a,b,c);    /* clear solution container */
   fail = _ada_use_c2phc(16,a,b,c);    /* do path tracking */
   fail = _ada_use_c2phc(5,a,b,c);     /* copy target sols to container */
   fail = _ada_use_c2phc(34,&k,&m[1],sol); 
   /* return a solution from the container */
   fail = _ada_use_c2phc(37,a,b,c);    /* clear solution container */
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
      if(v) printf("There are %d start solutions.\n",*nbsols);
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

/*
void print_system ( void ) 
{
   int *a,*b,fail;
   double *c;

   fail = _ada_use_c2phc(21,a,b,c);
}
*/

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

void print_time ( FILE *fp, double *time, int numprocs, int *mysolnb )
{
   int i;
   fprintf(fp,"\nTotal wall time = %lf seconds on %d processors\n",
           time[0],numprocs);
   fprintf(fp,"Total number of paths tracked: %d.\n", mysolnb[0]);

   for(i=1; i<numprocs; i++)
   {
      fprintf(fp,"The wall time spent on NO. %d processor = %lf seconds.\n",
              i,time[i]);
      fprintf(fp,"The number of paths tracked on NO. %d processor: %d.\n",
              i,mysolnb[i]);
   }
   free(time);
   free(mysolnb);
}
