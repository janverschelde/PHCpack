#include <stdio.h> 
#include <stdlib.h> 
#include <time.h>
#include <math.h>
#include "mpi.h"
#include "parallel_phcpack.h"
#include "parallel_cells.h"

#include "../Lib/celcon.h"
#include "../Lib/solcon.h"

#define SEND_CELL1 105
#define SEND_CELL 100 /* tag for sending cell number and #solutions */
#define SEND_SUPP 101 /* tag for sending  cell support */
#define SEND_NORMAL 102 /* tag for sending cell normal */
#define SEND_MUL 103 /* tag for sending mulplicity of the solution */
#define SEND_SOL 104 /* tag for sending solution */

#define v 0  /* verbose flag:
                  0 no output during the computations,
                  1 only one line messages on progress of computations
                  2 display of data, like systems and solutions */

/* extern void adainit( void );
   extern int _ada_use_c2phc ( int task, int *a, int *b, double *c );
   extern void adafinal( void ); */

void Print_Integer_Array ( int n, int *a );
/*
 * DESCRIPTION :
 *   Writes n integer numbers in a to screen. */

void Print_Double_Array ( int n, double *a );
/*
 * DESCRIPTION :
 *   Writes n doubles stored in a to screen. */

void print_cell ( int n, int *cell );

void distribute_cells ( int myid, int np, int nspt, int dim, int *mysolnb ); 
/*
 * DESCRIPTION :
 *   Dynamic distribution of the cells among the computing nodes. */

int main ( int argc, char *argv[] )
{
   int np,myid,dim,nspt,mysolnb,*mysol;
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
   distribute_cells(myid,np,nspt,dim,&mysolnb);
   
   endwtime = MPI_Wtime();
   wtime = endwtime-startwtime;
   MPI_Gather(&wtime,1,MPI_DOUBLE,time,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Gather(&mysolnb,1,MPI_INT,mysol,1,MPI_INT,0,MPI_COMM_WORLD);
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

void distribute_cells ( int myid, int np, int nspt, int dim, int *mysolnb )
{
   int nbcell,*cell,nbsols,i,j,R,m[2],sn,dest;
   int count,labSize,cnt;
   MPI_Status status;
   MPI_Request request;
   MPI_Comm com = MPI_COMM_WORLD;
   int length[nspt],*labels ;
   int fail,k,mult,n,mold,st,lastcellpath;
   double normal[dim],sol[2*(dim-1)+5];
   
   n=dim-1;
   if(myid == 0)
   { 
      fail=celcon_number_of_cells(&nbcell);  /* get the number of cells */
      cell=(int*)calloc(nbcell,sizeof(int));
      for(i=0;i<nbcell;i++)
         fail = celcon_mixed_volume(i+1,&cell[i]);

      nbsols = 0;
      for(i=0;i<nbcell;i++) nbsols=nbsols+cell[i];   
      if(v>0) printf("The total number of cells is %d.\n",nbcell);
      if(v>0) print_cell(nbcell,cell);
      if(v>0) printf("The total #solutions is %d.\n",nbsols);
   }
   MPI_Bcast(&nbcell,1,MPI_INT,0,com);
   *mysolnb=0;
   
   if(myid == 0)
   {
      fail=celcon_number_of_points_in_cell(1,nspt,length);
      labSize=0;
      for(i=0;i<nspt;i++) labSize=labSize+length[i];
      labSize=1+nspt+labSize;
   }
   MPI_Bcast(&labSize,1,MPI_INT,0,MPI_COMM_WORLD);
   labels = (int*)calloc(labSize,sizeof(int));
   if(myid == 0)
   {  
      printf("writing random coefficient system and its solutions to file\n");
      fail=celcon_write_standard_random_coefficient_system();
      fail=solcon_write_solution_dimensions_to_defined_output_file(nbsols,n);
      cnt=0;
      if(nbcell<np)
      {
         k=1;
         for(i=1;i<=nbcell;i++)
         {
            celcon_retrieve_mixed_cell(dim,nspt,i,labels,normal);
            *mysolnb=*mysolnb+cell[i-1];
            for(j=1;j<=cell[i-1];j++)
            {
               m[0]=i; m[1]=j;
               if(k<np)
               {
                  MPI_Send(m,2,MPI_INT,k,SEND_CELL,com);
         	  MPI_Send(labels,labSize,MPI_INT,k,SEND_SUPP,com); 
                  MPI_Send(normal,dim,MPI_DOUBLE,k,SEND_NORMAL,com);
         	  k++;
	       }
               else
               {
                  MPI_Recv(&R,1,MPI_INT,MPI_ANY_SOURCE,SEND_CELL1,com,&status);
                  dest=status.MPI_SOURCE;
                  MPI_Send(m,2,MPI_INT,dest,SEND_CELL,com);
                  if(R!=i)
                  {
                     MPI_Send(labels,labSize,MPI_INT,dest,SEND_SUPP,com);
                     MPI_Send(normal,dim,MPI_DOUBLE,dest,SEND_NORMAL,com); 
                  }
                  MPI_Recv(&mult,1,MPI_INT,MPI_ANY_SOURCE,SEND_MUL,
                           com,&status);
                  MPI_Recv(sol,2*n+5,MPI_DOUBLE,MPI_ANY_SOURCE,SEND_SOL,
                           com,&status);
                  fail = solcon_write_next_solution_to_defined_output_file
                           (&cnt,n,mult,sol);  
               }   
            } /* end for j */
         } /* end for i */
         m[0] = 0;
         for(i=1;i<np;i++)
	 {
            MPI_Send(m,2,MPI_INT,i,SEND_CELL,com);
            if(i<k)
	    {
               MPI_Recv(&mult,1,MPI_INT,MPI_ANY_SOURCE,SEND_MUL,com,&status);
               MPI_Recv(sol,2*n+5,MPI_DOUBLE,MPI_ANY_SOURCE,SEND_SOL,
                        com,&status);
               fail = solcon_write_next_solution_to_defined_output_file
                        (&cnt,n,mult,sol); 
            }
         }
      } /* end big if */
      else
      {
         for(i=1;i<=np-1;i++)
         {
            m[0]=i;
            m[1]=cell[i-1];
            celcon_retrieve_mixed_cell(dim,nspt,i,labels,normal);
            MPI_Send(m,2,MPI_INT,i,SEND_CELL,com);
            MPI_Send(labels,labSize,MPI_INT,i,SEND_SUPP,com);
            MPI_Send(normal,dim,MPI_DOUBLE,i,SEND_NORMAL,com);
            *mysolnb=*mysolnb+m[1];
         }   
         lastcellpath=cell[nbcell-1];
      
         for(i=np;i<nbcell+np-1+lastcellpath;i++)
         {
            MPI_Recv(&R,1,MPI_INT,MPI_ANY_SOURCE,SEND_CELL1,com,&status);
            dest=status.MPI_SOURCE;
            if(i<nbcell)
            {
               m[0]=i; m[1]=cell[i-1];
            }
            else  if(i<nbcell+lastcellpath)
                    {  m[0]=nbcell;  m[1]=i-nbcell+1; }
              else {  m[0]=0; m[1]=0; }
              MPI_Send(m,2,MPI_INT,dest,SEND_CELL,com);
	      if(v>0) printf("root sends %d paths of cell %d to worker %d\n",
                              m[1],m[0],dest);
              if(i<nbcell+lastcellpath && R!=0)
              {
                 celcon_retrieve_mixed_cell(dim,nspt,i,labels,normal);
                 MPI_Send(labels,labSize,MPI_INT,dest,SEND_SUPP,com);
                 MPI_Send(normal,dim,MPI_DOUBLE,dest,SEND_NORMAL,com);
		 if(v>0) printf("root sends label to worker %d\n",dest);
              }
              if(R==0) R=1;
              for(j=0; j<R; j++)
              {
                 MPI_Recv(&mult,1,MPI_INT,MPI_ANY_SOURCE,SEND_MUL,
                          MPI_COMM_WORLD,&status);
                 MPI_Recv(sol,2*n+5,MPI_DOUBLE,MPI_ANY_SOURCE,SEND_SOL,
                          MPI_COMM_WORLD,&status);
                 fail = solcon_write_next_solution_to_defined_output_file
                          (&cnt,n,mult,sol);
                 *mysolnb=*mysolnb+R;
              }
          
           } 
      } /* end big else */
   } /* myid=0 finish */
   else
   { 
     if(nbcell<np)
     {
        sn=0; mold=0;
        while(1)
	{
           MPI_Recv(m,2,MPI_INT,0,SEND_CELL,com,&status);
           if(m[0] == 0) break;
            
           if(m[0]!=mold)
	   {
              sn++; st = 1;
              mold = m[0];
              MPI_Recv(labels,labSize,MPI_INT,0,SEND_SUPP,com,&status) ;
              MPI_Recv(normal,dim,MPI_DOUBLE,0,SEND_NORMAL,com,&status) ;
              fail = celcon_append_mixed_cell(dim,nspt,labSize,labels,normal);
              fail = celcon_standard_polyhedral_homotopy();
              fail = celcon_solve_standard_start_system(sn,&R);
           }
           else
              st++;
            
           fail = celcon_track_standard_solution_path(sn,m[1],0);
           fail = solcon_clear_standard_solutions();
           fail = celcon_copy_target_standard_solution_to_container(sn,st);
           fail = solcon_retrieve_standard_solution(n,1,&mult,sol);
           MPI_Send(&mult,1,MPI_INT,0,SEND_MUL,com);
           MPI_Send(sol,2*n+5,MPI_DOUBLE,0,SEND_SOL,com);
            
           *mysolnb=*mysolnb+1;
           MPI_Send(&m[0],1,MPI_INT,0,SEND_CELL1,com);                    
        } /* end while */
     } /* end if */
     else
     {
        sn=0; mold=0;st=0;
        while(1) 
        {
           MPI_Recv(m,2,MPI_INT,0,SEND_CELL,com,&status);
           if(v>0) printf("worker %d receives cell %d\n",myid,m[0]);
           if(m[0]==0) break;
           if(m[0]!=mold)
           {
              sn++; mold=m[0];
              MPI_Recv(labels,labSize,MPI_INT,0,SEND_SUPP,com,&status) ;
              MPI_Recv(normal,dim,MPI_DOUBLE,0,SEND_NORMAL,com,&status) ;
              if(v>0) printf("worker %d receives label\n",myid);
              fail = celcon_append_mixed_cell(dim,nspt,labSize,labels,normal);
              fail = celcon_standard_polyhedral_homotopy();
              fail = celcon_solve_standard_start_system(sn,&R);
           }
           if(m[0] != nbcell)
           {
              R=m[1];
              for(j=1;j<=m[1];j++)
              {
                 fail = celcon_track_standard_solution_path(sn,j,0);
                 fail = solcon_clear_standard_solutions();  
                 fail = celcon_copy_target_standard_solution_to_container(sn,j);
                 fail = solcon_retrieve_standard_solution(n,1,&mult,sol);
                 MPI_Send(&mult,1,MPI_INT,0,SEND_MUL,MPI_COMM_WORLD);
                 MPI_Send(sol,2*n+5,MPI_DOUBLE,0,SEND_SOL,MPI_COMM_WORLD);
              } 
           }
           else
	   {
              R=0; st++;
              fail = celcon_track_standard_solution_path(sn,m[1],0);
              fail = solcon_clear_standard_solutions();
              fail = celcon_copy_target_standard_solution_to_container(sn,st);
              fail = solcon_retrieve_standard_solution(n,1,&mult,sol);
              MPI_Send(&mult,1,MPI_INT,0,SEND_MUL,MPI_COMM_WORLD);
              MPI_Send(sol,2*n+5,MPI_DOUBLE,0,SEND_SOL,MPI_COMM_WORLD);
           }
           *mysolnb=*mysolnb+R+1;
           MPI_Send(&R,1,MPI_INT,0,SEND_CELL1,com);
         } /* end while */
      } /* end big else */
   } /* end worknode */

 /* MPI_Barrier(com); */
   if(myid == 0)
   {
      free(cell);
      *mysolnb = nbsols;  /* a patch ... */
   }
   free(labels);
}

void Print_Integer_Array ( int n, int* a )
{
   int i;
   for(i=0;i<n;i++) printf("  %d",a[i]);
   printf("\n");
}

void Print_Double_Array ( int n, double* a )
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
