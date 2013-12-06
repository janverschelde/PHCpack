#include<stdio.h>
#include<stdlib.h>
#include "mpi.h"
#include "parallel_tree.h"
#include "queue.h"

#define SEND_PIVOT 100   /* tag for sending pivots */
#define SEND_SSOL  101   /* tag for sending start solutions */
#define SEND_PIVOT_B 102 /* tag for sending pivot back to the root node */
#define SEND_TSOL  103   /* tag for sending target solutions */
#define SEND_RES 104     /* tag for sending residuals */

#define v 1 /* verbose flag */

extern void adainit();
extern int _ada_use_c2pieri ( int job, int *a, int *b, double *c );
extern void adafinal();

void dimension_broadcast ( int myid, int *m, int *p, int *q, int *nd_tp );
/* Root node reads the dimension information and the pivot method to use, then
   broadcasts to all the nodes */

void input_planes_broadcast( int myid, int m, int p, int q );
/* Root node initializes the machine with m*p + q*(m+p) random input m-planes,
 * then broadcasts to all the nodes */

void interpolation_points_broadcast(int myid, int m, int p, int q);
/* Root node initializes the machine with m*p + q*(m+p) random interpolation 
   points, then broadcasts to all the nodes  */

int sever_distribute ( int m, int p, int q, int numprocs, int nd_tp );
/* Server generates jobs with a tree structure, then distributes jobs to other
   nodes with a queue. Returns the total number of solutions */ 

void client_compute ( int m, int p, int q );
/* Clients do the computation jobs received from the server with Ada program */


int main(int argc, char* argv[])
{
  int numprocs,myid,fail=-1;
  int m,p,q,sum,nd_tp;
  double startwtime, endwtime;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

  srand(time(NULL));
  adainit(); 

  dimension_broadcast(myid,&m,&p,&q,&nd_tp); 
  input_planes_broadcast(myid,m,p,q);
  interpolation_points_broadcast(myid,m,p,q);

  if(myid==0)
  {
    startwtime = MPI_Wtime();
    sum = server_distribute(m,p,q,numprocs,nd_tp);
    endwtime = MPI_Wtime();
    printf("\nTotal wall time = %lf seconds on %d processors\n",
           endwtime-startwtime, numprocs);
  }
  else 
    client_compute(m,p,q);  

  adafinal();      
  MPI_Finalize();

  return 0;
}

void dimension_broadcast ( int myid, int *m, int *p, int *q, int *nd_tp )
{
  int fail=-1,mpq[3];
  int *b;
  double *c;

  if(myid==0)
  {
    printf("Give p, dimension of the solution planes :");
    scanf("%d", p);
    printf("Give m, the co-dimension so that n = m+p :");
    scanf("%d", m);
    printf("Give q, the degree of the maps :");
    scanf("%d", q);
    while(1)
    {
      printf("Menu of pivot methods:
      1. top pivots.
      2. bottom pivots.
      3. mixed pivots.\n");
      printf("Type 1, 2, 3 to choose: ");
      scanf("%d", nd_tp);
      if(*nd_tp<1 || *nd_tp>3)
        printf("Please choose 1~3.\n");
      else
        break;
    }
  }
  MPI_Bcast(p, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(m, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(q, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(nd_tp, 1, MPI_INT, 0, MPI_COMM_WORLD);


  mpq[0] = *m;
  mpq[1] = *p;
  mpq[2] = *q;

  fail = _ada_use_c2pieri(1,mpq,b,c);

  if(fail!=0)
  {
    if(myid==0)
      printf("Invalid choice. Please quit the program and try again...\n");
  }
  if(v==1)
    printf("No. %d processor:p=%d,m=%d,q=%d\n",myid,*p,*m,*q);   

} 

void input_planes_broadcast( int myid, int m, int p, int q )
{
   int n = m*p + q*(m+p);
   int i,nbcff = 2*(m+p)*m*n;
   int fail,mpq[3];
   double c[nbcff];

   mpq[0] = m; mpq[1] = p; mpq[2] = q;     

   if(myid==0)
   {
     for (i=0; i<nbcff;i++) c[i] = ((double) rand())/RAND_MAX;
   }
   MPI_Bcast(c, nbcff, MPI_DOUBLE, 0, MPI_COMM_WORLD);

   fail = _ada_use_c2pieri(2,mpq,&n,c);

}

void interpolation_points_broadcast(myid,m,p,q)
{
   int n = m*p + q*(m+p);
   int i,nbcff = 2*(m+p)*m*n;
   int fail,mpq[3];
   double c[nbcff];

   mpq[0] = m; mpq[1] = p; mpq[2] = q;
   nbcff = 2*n;

   if(myid==0)
     for (i=0; i<nbcff;i++) c[i] = ((double) rand())/RAND_MAX;

   MPI_Bcast(c, nbcff, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   fail = _ada_use_c2pieri(3,mpq,&n,c);
}

void package(int p,int *start,int *target,int length,int address,int pack[4*p+2])
/* Saves top and bottom pivots of the start solution in the first 2*p elements,
   appends top and bottom pivots of the target solution in the second 2*p 
   elements,length of solution in the (4p+1)th position,address of the node in
   the (4p+2)th position      */
{
  int i;
  
  for(i=0; i<2*p; i++)
    pack[i] = start[i];

  for(i=0; i<2*p; i++)
    pack[i+2*p] = target[i];

  pack[4*p] = length;
  pack[4*p+1] = address;
}

void print_solution(int length, double *sol)
/* Prints out the solutions */
{
   int i;

   for(i=0; i<length; i=i+2)
   {
     if(sol[i]>0) printf(" ");
     printf("%.15e\t",sol[i]);
     if(sol[i+1]>0) printf(" ");
     printf("%.15e\n",sol[i+1]);
   }
}

int server_distribute ( int m, int p, int q, int numprocs, int nd_tp )
{
  int length,buffer[4*p+2],counter,sol_num=0,/* the total number of solutions*/
      dest=1,source,client_number = numprocs-1,finish=-1;
  int i,plength=4*p+2,pack[plength],dim = m*p+q*(m+p),wait_cpu[numprocs-1]; 

  Node *leave,*nd; 
  double *sol=NULL,res;
  queue *front=NULL, *end=NULL; 
  /* q points the front of the queue, end points the end of the queue */
  job *j;  
  MPI_Status status;

  if(nd_tp==1)
  {
    nd=Top_Root(m,p,q);
    leave=Top_Leave(m,p,q);
  }
  else if(nd_tp==2)
  {
    nd=Bottom_Root(p);
    leave=Bottom_Leave(m,p,q);
  }
  else
  {
    leave=Mixed_Leave(m,p,q);
    Mixed_Root(m,p,q,leave,&front,&end);     
  }

  for(i=0; i<numprocs-1; i++)
    wait_cpu[i] = i+1;
  /* At the beginning all the client processors are waiting for jobs */

  //////////// Check nd for mixed!!!!!!!!!!!!!!!!!

  while(1)
  { 
      if(nd->level==dim)
      /* if nd is a leave,prints out the solution and put the 
         client ID in an array */
      {
        sol_num++;
        printf("Soltion No. %d\n",sol_num);
        print_solution(dim*2,sol);
        MPI_Recv(&res,1,MPI_DOUBLE,MPI_ANY_SOURCE,SEND_RES,
                 MPI_COMM_WORLD,&status);
        printf("The residual is : %.15e\n\n",res);
        Deallocate(nd,p);
	/*deallocate memory for the leave*/

        if(v==1)
          printf("client_number=%d\n",client_number);
        if(client_number==(numprocs-1)) /* All of the clients returned leaves */
          break;
      }
      else
      { 
        if(nd_tp==1)
        {
          if(v==1)
          {
            printf("Create top parent for node ");
            Print_Localization(p,nd->top);
          }
          Q_Create_Top_Parent(nd,m,p,q,leave,sol,&front,&end);
        }
        else if(nd_tp==2)
        {
          if(v==1)
          {
            printf("Create bottom parent for node ");
            Print_Localization(p,nd->bottom);
          }
          Q_Create_Bottom_Parent(nd,m,p,q,leave,sol,&front,&end);
        }
        else
          printf("haha\n");         
      }
      if(sol!=NULL)
        free(sol);
         
      while((front!=NULL)&&(client_number>0))
    /* Distributes one job to each of the client processors at the beginning */
      {       
         front = pop(front, &j);
         package(p,j->start_pivot,j->target_pivot,j->length,j->address,pack);
         MPI_Send(pack,plength,MPI_INT,wait_cpu[client_number-1],
                  SEND_PIVOT,MPI_COMM_WORLD);
         MPI_Send(j->start_sol,j->length,MPI_DOUBLE,wait_cpu[client_number-1],
                  SEND_SSOL,MPI_COMM_WORLD);
	 if(v==1) 
         {
           printf("Sending a job to No. %d processor: ",wait_cpu[client_number-1]);
           Print_Localization(plength, pack);
         }

         Deallocate_job(j);
         client_number--; 
      }

      MPI_Recv(buffer,plength,MPI_INT,MPI_ANY_SOURCE,SEND_PIVOT_B,
               MPI_COMM_WORLD,&status);
      length = buffer[4*p];
      sol = (double*)calloc(length, sizeof(double));
      source = status.MPI_SOURCE;

      MPI_Recv(sol,length,MPI_DOUBLE,source,SEND_TSOL,
               MPI_COMM_WORLD,&status);
      wait_cpu[client_number++] = source;

      if(v==1)
      { 
        printf("Received pivot array from No. %d processor: ", source);
        Print_Localization(4*p+2, buffer); 
        printf("Received solution from No. %d processor:\n", source);
        for(i=0; i<length; i=i+2)
        {
          if(sol[i]>0) printf(" ");
          printf("%.14e\t",sol[i]);
          if(sol[i+1]>0) printf(" ");
          printf("%.14e\n",sol[i+1]);
        }
      }
      nd = (Node*) buffer[4*p+1];  /* Assigns the node address to nd */
    }  

  /* Tell all the clients no more jobs to be done */
  for(i=1; i<numprocs; i++)
     MPI_Send(&finish,1,MPI_INT,i,SEND_PIVOT,MPI_COMM_WORLD);
  
  Deallocate(leave,p);
  printf("The number of roots : %d\n", sol_num);
  return sol_num;
}

void client_compute ( int m, int p, int q )
{
  int length,plength=4*p+2,buffer[plength],
      i,dim=m*p+q*(m+p),fail,*a,*b,deg_freedom;
  double ssol[dim*2],tsol[dim*2],*c,res; 
  MPI_Status status;

  while(1)
  { 
    MPI_Recv(buffer,plength,MPI_INT,0,SEND_PIVOT,MPI_COMM_WORLD,&status);
    if(buffer[0]==-1) break;
    length = buffer[4*p];
    MPI_Recv(ssol,length,MPI_DOUBLE,0,SEND_SSOL,MPI_COMM_WORLD,&status);

    fail = _ada_use_c2pieri(4,&p,buffer,c); 
    /* store pivots of pattern at start solution curve */

    fail = _ada_use_c2pieri(5,&p,&buffer[2*p],c);
    /* store pivots of pattern at target solution curve */

    deg_freedom = length/2; 
    fail = _ada_use_c2pieri(6,&deg_freedom,b,ssol);
    
    /* store coefficients of start solution curve */
    if(v==1)
    {
      printf("The start curve is:\n");
      for(i=0;i<2*deg_freedom;i=i+2)
      {
        if(ssol[i]>0) printf(" ");
        printf("%.14e\t",ssol[i]);
        if(ssol[i+1]>0) printf(" ");
        printf("%.14e\n",ssol[i+1]);
      }
    }
    fail = _ada_use_c2pieri(9,a,b,tsol);
    /* track solution path {8: without or 9: with} intermediate output */

    deg_freedom++;
    fail = _ada_use_c2pieri(7,&deg_freedom,b,tsol);
    /* retrieve coefficients of target solution curve */
    
    buffer[4*p] = buffer[4*p]+2;
    MPI_Send(buffer,plength,MPI_INT,0,SEND_PIVOT_B,MPI_COMM_WORLD);
    MPI_Send(tsol,length+2,MPI_DOUBLE,0,SEND_TSOL,MPI_COMM_WORLD);

    fail = _ada_use_c2pieri(11,a,b,&res);  
    /* verify intersection conditions {10: without or 11: with} output */      
    if(deg_freedom==dim)
    {
      MPI_Send(&res,1,MPI_DOUBLE,0,SEND_RES,MPI_COMM_WORLD);
    }
  }
}
