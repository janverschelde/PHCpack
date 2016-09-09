#include <stdio.h> 
#include <stdlib.h> 
#include <time.h>
#include <math.h>
#include "mpi.h"
#include "parallel_phcpack.h"
#include "parallel_cells.h"

#include "../Lib/celcon.h"
#include "../Lib/solcon.h"

#define SEND_CELL 100 /* tag for sending cell number and #solutions */
#define SEND_SUPP 101 /* tag for sending  cell support */
#define SEND_NORMAL 102 /* tag for sending cell normal */
#define SEND_MUL 103 /* tag for sending mulplicity of the solution */
#define SEND_SOL 104 /* tag for sending solution */

#define v 2  /* verbose flag:
                  0 no output during the computations,
                  1 only one line messages on progress of computations
                  2 display of data, like systems and solutions */

/* extern void adainit( void );
   extern int _ada_use_c2phc ( int task, int *a, int *b, double *c );
   extern void adafinal( void ); */

typedef struct lis lis;
typedef struct lisStack  lisStack;

struct lis
{
   int* data;
   lis* next;
};

struct lisStack
{ 
   int count;
   lis* top;
   lis* cur;
};

void lis_init ( lis *c, int *m, lis *ptr );

void ls_init ( lisStack *cs );

void ls_push ( lisStack *cs, int *m );

void ls_pop ( lisStack *cs );

int ls_next ( lisStack *cs );

int* ls_cur ( lisStack *cs );

int ls_isempty ( lisStack *cs );

void Print_Integer_Array ( int n, int *a );
/*
 * DESCRIPTION :
 *   Writes n integer numbers in a to screen. */

void Print_Double_Array ( int n, double *a );
/*
 * DESCRIPTION :
 *   Writes n doubles stored in a to screen. */

void print_cell ( int n, int *cell );

void distribute_cells ( int myid, int np, int nspt, int dim, int *nbpaths );
/*
 * DESCRIPTION :
 *   Static distribution of the cells among the computing nodes. 
 *
 * ON ENTRY :
 *   myid     identification number of the node;
 *   np       total number of processors - 1;
 *   nspt     number of distinct supports;
 *   dim      ambient dimension.
 *
 * ON RETURN :
 *   nbpaths  the number of paths tracked by the computing node. */

int main ( int argc, char *argv[] )
{
   int np,myid,dim,nspt,nbpaths,*mysol;
   double startwtime,endwtime,wtime,*time;
   MPI_Status status;

   adainit();
   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD,&np);
   MPI_Comm_rank(MPI_COMM_WORLD,&myid);
   
   if(myid == 0)
   {
      time = (double*)calloc(np,sizeof(double));
      mysol = (int*)calloc(np,sizeof(int));
      startwtime = MPI_Wtime();
   }
   else
      startwtime = MPI_Wtime();
 
   retrieve_dimensions(myid,&nspt,&dim);
   supports_broadcast(myid,nspt,dim);
   system_broadcast(myid,dim-1);
   distribute_cells(myid,np,nspt,dim,&nbpaths);
   
   endwtime = MPI_Wtime();
   wtime = endwtime-startwtime;
   MPI_Gather(&wtime,1,MPI_DOUBLE,time,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Gather(&nbpaths,1,MPI_INT,mysol,1,MPI_INT,0,MPI_COMM_WORLD);
   if(myid == 0)
   {
      printf("\nTotal wall time = %lf seconds on %d processors\n",
             time[0],np);
      write_time_and_paths_to_defined_output_file(np,time,mysol);
      free(time); free(mysol);
   }
    
   MPI_Finalize();
   adafinal();
   return 0;
}

void distribute_cells ( int myid, int np, int nspt, int dim, int *nbpaths )
{
   int nbcell,*cell,nbsols,i,j,k,left[np],L,R,*m,sn,mysolnum;
   int count,labSize,cnt,dest;
   MPI_Status status;
   MPI_Comm com = MPI_COMM_WORLD;
   int length[nspt],*labels ;
   int fail,A[2],n;
   double normal[dim],sol[2*(dim-1)+5],*c;
   lisStack *s;
   int *a,*b;
   
   n=dim-1;A[0]=n;
   if(myid == 0)
   { 
      fail=celcon_number_of_cells(&nbcell);  /* get the number of cells */
      cell=(int*)calloc(nbcell,sizeof(int));
      for(i=0; i<nbcell; i++)
         fail = celcon_mixed_volume(i+1,&cell[i]);

      nbsols=0;
      for(i=0; i<nbcell; i++) nbsols=nbsols+cell[i];   
      if(v>0) printf("The number cells are %d\n",nbcell);
      if(v>0) print_cell(nbcell,cell);
      if(v>0) printf("The total solutions are %d\n",nbsols);
   }
   MPI_Bcast(&nbsols,1,MPI_INT,0,com);
   
   if(myid == 0)
   {
      if(v>0) printf("\n");
      left[0] = 0;
      for(i=1; i<np; i++)
         if(i <= (nbsols%(np-1))) 
            left[i] = nbsols/(np-1)+1; 
         else
            left[i] = nbsols/(np-1);
      if(v>0) printf("left:");
      if(v>0) Print_Integer_Array(np,left);
   }

   MPI_Scatter(left,1,MPI_INT,&mysolnum,1,MPI_INT,0,com);

   if(myid == 0)
   {
      fail = celcon_number_of_points_in_cell(1,nspt,length);
      labSize=0;
      for(i=0; i<nspt; i++) labSize = labSize+length[i];
      labSize = 1+nspt+labSize;
   }
   MPI_Bcast(&labSize,1,MPI_INT,0,MPI_COMM_WORLD);
   labels = (int*)calloc(labSize,sizeof(int));
   m=(int*)calloc(3,sizeof(int));

   if(myid==0)
   {
      L=1; R=np-1;
      for(i=0; i<nbcell; i++)
      { 
         m[0] = i+1;
         m[2] = labSize;
         celcon_retrieve_mixed_cell(dim,nspt,i+1,labels,normal);
         
         sn=1;
         while(cell[i]!=0)
         { 
            if(cell[i]>=left[L])
            {   
               m[1] = left[L]; m[2] = sn;      
	       MPI_Send(&m[1],2,MPI_INT,L,SEND_CELL,com);
               /*
               if(v>0)
                  printf("%2d paths from cell %d is sending to node %d\n",
                         m[1],m[0],L);
                */
               MPI_Send(labels,labSize,MPI_INT,L,SEND_SUPP,com);
               MPI_Send(normal,dim,MPI_DOUBLE,L,SEND_NORMAL,com) ;
               cell[i] = cell[i]-left[L]; sn = sn+left[L];
               L++;  
            }
            else if(cell[i]>=left[R])
            { 
               m[1] = left[R]; m[2] = sn;      
               MPI_Send(&m[1],2,MPI_INT,R,SEND_CELL,com);
               /*
               if(v>0)
                  printf("%2d paths from cell %d is sending to node %d\n",
                         m[1],m[0],R);
               */
               MPI_Send(labels,labSize,MPI_INT,R,SEND_SUPP,com);
               MPI_Send(normal,dim,MPI_DOUBLE,R,SEND_NORMAL,com) ;
               cell[i] = cell[i]-left[R]; sn = sn+left[R];
               R--;   
            }
            else
            {              
               m[1] = cell[i]; m[2] = sn;     
               MPI_Send(&m[1],2,MPI_INT,R,SEND_CELL,com);
               /*
               if(v>0)
                  printf("%2d paths from cell %d is sending to node %d\n",
                         m[1],m[0],R);
               */
               MPI_Send(labels,labSize,MPI_INT,R,SEND_SUPP,com);
               MPI_Send(normal,dim,MPI_DOUBLE,R,SEND_NORMAL,com);
               left[R]=left[R]-cell[i]; sn = sn+cell[i];
               cell[i]=0;  
            }
         }
      }
      if(v>0) printf("****************************************************\n");
      printf("writing random coefficient system and its solutions to file\n");
      fail = celcon_write_standard_random_coefficient_system();
      fail = solcon_write_solution_dimensions_to_defined_output_file(nbsols,n);
      cnt = 0;
      for(k=1; k<=nbsols; k++)
      {
         MPI_Recv(&A[1],1,MPI_INT,MPI_ANY_SOURCE,SEND_MUL,
                  MPI_COMM_WORLD,&status);
         MPI_Recv(sol,2*n+5,MPI_DOUBLE,MPI_ANY_SOURCE,SEND_SOL,
                  MPI_COMM_WORLD,&status);
         fail = solcon_write_next_solution_to_defined_output_file
                  (&cnt,n,A[1],sol);
      }
      *nbpaths = nbsols;

   } /* myid=0 finish */
   else
   {  
      *nbpaths = mysolnum;
      if(v>0) printf("Node %d has %d paths\n",myid,mysolnum);
      s=(lisStack*)calloc(1,sizeof(lisStack));
      ls_init(s);
      count = 0; sn = 0;
      while(count < mysolnum)
      {  
         MPI_Recv(&m[1],2,MPI_INT,0,SEND_CELL,com,&status);
         sn++;
         m[0] = sn;
         ls_push(s,m);
         count = count+m[1];
         /*
     	 if(v>0) printf("Node %d is receving %2d paths from cell %d\n",
		        myid,m[1],m[0]);
         */
	 MPI_Recv(labels,labSize,MPI_INT,0,SEND_SUPP,com,&status) ;
	 MPI_Recv(normal,dim,MPI_DOUBLE,0,SEND_NORMAL,com,&status) ;
      	 fail = celcon_append_mixed_cell(dim,nspt,labSize,labels,normal);
      }
      fail = celcon_standard_polyhedral_homotopy();
      
      for(i=1; i<=sn; i++)
      {
         m = ls_cur(s);
         fail = celcon_solve_standard_start_system(m[0],&R);

         if(fail == 1)
         {
            printf("Solving start system failed.\n");
            printf("Node %d skips cell %d with volume %d...\n",myid,m[0],R);
         }
         else
         {
    /* printf("found %d start solutions from cell %d\n",R,m[0]); */

            fail=celcon_mixed_volume(m[0],&L);

    /* if(R==L)  printf("#start solutions equals mixed volume %d, OK\n",L);
       else  printf("#start solutions not equals mixed volume %d!!!, \n",L);
    */
            for(j=m[2]; j<m[1]+m[2]; j++)
            {          
               fail = celcon_track_standard_solution_path(m[0],j,0);
               fail = solcon_clear_standard_solutions();
               fail = celcon_copy_target_standard_solution_to_container
                        (m[0],j-m[2]+1);
               fail = solcon_retrieve_standard_solution(n,1,&A[1],sol);
               MPI_Send(&A[1],1,MPI_INT,0,SEND_MUL,MPI_COMM_WORLD);
               MPI_Send(sol,2*n+5,MPI_DOUBLE,0,SEND_SOL,MPI_COMM_WORLD);
            }
         }
         ls_pop(s);
      }
   } /* end else */

 /*  MPI_Barrier(com);  */  
   if(myid == 0)
      free(cell);
   else
      free(s);
   free(labels);
   free(m);
}

void Print_Integer_Array ( int n, int* a )
{
   int i;
   for(i=0;i<n;i++) printf("  %d",a[i]);
   printf("\n");
}

void Print_Double_Array ( int n, double *a )
{
   int i;
   for(i=0;i<n;i++) printf("  %f",a[i]);
   printf("\n");
}

void print_cell ( int n, int *cell )
{
  int i;
  printf("The start solution in each cell is as follows:\n");
  for(i=0;i<n;i++)
     printf("Cell %d:   %d\n",i+1,cell[i]);
}

void lis_init ( lis *c, int *m, lis *ptr )
{
   int i;
   c->data=(int*)calloc(3,sizeof(int));
   for(i=0;i<3;i++) c->data[i]=m[i];
   c->next=ptr;
}

void ls_init ( lisStack* cs )
{
   cs->count=0;
   cs->cur=cs->top=0;
}

void ls_push ( lisStack *cs, int *m )
{
   lis *c;
   c=(lis*)calloc(1,sizeof(lis));
   lis_init(c,m,cs->top);
   cs->cur=cs->top=c;
   ++cs->count;
}

void ls_pop ( lisStack *cs )
{
   lis *ptr=cs->top;
   cs->cur=cs->top=(cs->top)->next;
   --cs->count;
}

int ls_next ( lisStack *cs )
{
   if(cs->cur->next)
   {
      cs->cur=cs->cur->next;
      return 1;
   }
   else
      return 0;
}

int *ls_cur ( lisStack *cs )
{ 
   return cs->cur->data;
}

int ls_isempty ( lisStack *cs )
{
   if(cs->top==0)
      return 1;
   else
      return 0;
}
