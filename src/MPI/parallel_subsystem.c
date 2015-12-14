/* parallel subsystem-by-subsystem solver                */
/* first working version for katsura8 on Nov. 29 2007    */
/* on Dec. 1, only one job queue lp is used for all      */
/* TaskCount works fine, also clean up part              */
/* on Dec. 3, taking care of remainders case katsura10,  */
/*            but has problem with katsura7              */
/* Dec. 6 fixed file name(name_1, name_2) problems       */
/*        NUM_GROUP and NUM_VARS are imported from input */
/* Dec. 14 added function compute_witness_sets           */
/*         works for katsura6, 7, 8, 9, 10, 12           */
/* Dec. 17 fixed bug, works for katsura5, 11             */
/* Dec. 29 works for adjacent minors(24,28,22)           */
/*                   redeco 5,7,8                        */
/*         no first_cascade flag                         */
/* March 12 adding timings for each stage                */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include "mpi.h"
#include "../Lib/syscon.h"
#include "../Lib/phcpack.h"
#include "../Lib/solcon.h"
#include "../Lib/witset.h"
#include "parallel_phcpack.h"
#include "from_mpi2track.h"
#include "job_queue.h"
#include "idle_queue.h"

#define JOB_TAG 90
#define R_TAG 91
#define COUNT_TAG 92
#define EQN_TAG 93
#define DIM_TAG 94

#define v 1 /* verbose flag:
             0 only major output 
             1 display job queue info. 
             2 details after each step 
             3 detailed info.          */

extern void adainit();
extern void adafinal();

void combine_points(int *buffer, int *newbuffer);
/* Combine two witness points contained in buffer, then put in newbuffer */

void pack(int length, JOB *j, int package[length]);
/* Pack JOB */

LISTELEMENT* add_jobs (WSET *ws1, WSET *ws2, LISTELEMENT *listpointer);
/* Manager create & add jobs to the queue */

void manager_initialize(WSET *ws, int NUM_GROUP, int stages, 
                        int numbprocs, char filename[50]);
/* 1st part: Initialization phase. 
   The manager send eqn. index, then do synchronization */
/* 2nd part: mainloop initialization. 
   The manager read witness sets info. into the state table, */

void mainloop_initialize(WSET *ws, int NUM_GROUP, int stages, 
                         int numbprocs, char filename[50]);               
/* Called by func. manager_initialize. */ 
/* The manager read witness sets info. into the state table, */

void slave_initialize( int id, int NUM_GROUP, char filename[50]);
/* Corresponds to the 1st part of manager_initialize. */
/* Solving eqns. indicated by the indices */

void send_collect(WSET *ws, int num_group, int stages, int cnt_stage, 
                  int numbprocs, int cnt_step, IDLE_ELEMENT **ie, 
                  LISTELEMENT **listpointer, int *TaskCount, int NUM_GROUP, 
                  int *cnt_index, int update_flag, int n1, int n2, int cd);
/* Executed by the manager. Send jobs and collect results */
 
void run_slave(int id);
/* Code executed by the work with id */

void obtain_namebase(char filename[50], JOB *temp_j, 
                     char *namebase, char *name_1, char *name_2);
/* Obtain the filename for later input & output 
   based on info. contained in filename and temp_j. 
   namebase, name_1, name_2 are returned.
   name_1 and name_2 indicate the file contain witness set 
   for starting a diagonal homotopy. */

void source_to_string(int* source, int num, char* str);
/* Generate a string from source */

void obtain_string_from_source(JOB *temp_j, char filename[50], char *s_1, 
                               char *s_2, char *name_1, char *name_2);
/* Obtain two strings from source_1 and source_2 contained in JOB temp_j */

int read_two_witness_sets_from_file(int *n1, int *n2, int *dim1, 
                                    int *dim2, int *deg1, int *deg2,
                                    char *name_1, char *name_2, int *cd );
/* Read two witness sets from file indicated by name_1 and name_2 */

void compute_witness_sets(int id, int stages, int cnt_stage, int numbprocs,
                         int NUM_GROUP, int NUM_VARS, int num_group, WSET *ws,
                         int *TaskCount, int *cnt_index, char filename[50], 
                         int *remainder_pos, int remainder[2][2]); 
/* Compute witness sets for all eqn. groups at one stage 
   If there are two remainders, taking care of remainders 
   at the end of each stage.                             */                      
                       
void start_a_diagonal_homotopy(int myid, int numbprocs, char *name_1, 
                               char *name_2, char *outfile, WSET *ws,
                               int num_group, int stages, int cnt_stage, 
                               int cnt_step,
                               IDLE_ELEMENT **ie, LISTELEMENT **listpointer,
                               int *TaskCount, int NUM_GROUP, int NUM_VARS, 
                               int *cnt_index, int *expected,
                               int *dim, int update_flag);     
                                                     
/* First step: diagonal homotopy to start a cascade 
   returned is expected dimension & dimension of the witness set. */

void cascade_one_level_down(int myid, int numbprocs, char *infile, 
                            char *outfile, WSET *ws, int num_group, 
                            int stages, int cnt_stage, int cnt_step, IDLE_ELEMENT **ie, 
                            LISTELEMENT **listpointer, int *TaskCount,
                            int NUM_GROUP, int *dimension, int *cnt_index, 
                            int update_flag);                           
/* Second step: cascade to go one level down */
 
int remove_last_slack_variable(char *infile, char *outfile);
/* Removes the last slack variable of a witness set. */
                                                                 
void collapse(int cnt_stage, int stages, 
              char *infile, char *outfile, 
              WSET *ws, int *cnt_index,
              int *expected);                 
/* intersection step 3. executed by the manager */

int main ( int argc, char *argv[] )
{
   int myid, stages, numbprocs, *TaskCount, total=0,
       remainder[2][2], remainder_pos=-1, 
       remainderjob=0, cnt_index=0,
       NUM_GROUP, NUM_VARS, num_group=0, 
       cnt_stage=1, i, slv, count, r, c;   
   double startwtime,slavewtime,wtime,*mytime,
          stagewtime;
   char filename[50], ans;
   WSET *ws; 
   
   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD,&numbprocs);
   MPI_Comm_rank(MPI_COMM_WORLD,&myid);
   
   adainit(); 
   
   if(myid == 0)    
   {
     printf("Give a file name (less than 50 characters) : "); 
     scanf("%s",filename); 
     scanf("%c",&ans);  /* skip current new line symbol */ 
     
     printf("\n# of eqn. group :");
     scanf("%d",&NUM_GROUP);
     scanf("%c",&ans);
     
     printf("\n# of vars :");
     scanf("%d",&NUM_VARS);
     scanf("%c",&ans);
     
   }
   /* broadcast data */
   MPI_Barrier(MPI_COMM_WORLD);
   if(myid == 0) printf("\nManager broadcast data. \n");
   MPI_Bcast(filename,50*sizeof(char),MPI_CHAR,0,MPI_COMM_WORLD); 
   MPI_Bcast(&NUM_GROUP,sizeof(int),MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&NUM_VARS,sizeof(int),MPI_INT,0,MPI_COMM_WORLD); 
   MPI_Barrier(MPI_COMM_WORLD); 
   
   num_group = NUM_GROUP;
   stages = (int)ceil(log(NUM_GROUP)/log(2));
   if(myid == 0)
   {  
     /* TaskCount = (int*) calloc(numbprocs-1, sizeof(int));*/
     TaskCount = (int*) malloc(sizeof(int)*(numbprocs-1));   
     ws = (WSET*) calloc(NUM_GROUP*(stages+1), sizeof(WSET)); 
    /*  mytime = (double*) calloc(numbprocs, sizeof(double)); */   
     mytime = (double*) malloc(sizeof(double)*numbprocs);

     /* initialize TaskCount & remainder array */
     cnt_index=0;
     for(slv=1;slv<=numbprocs-1;slv++)
     {
     *(TaskCount+slv-1)=0;         
     }
     remainder_pos = -1;
     for(i=1;i<=2;i++)
     {
     remainder[i-1][0] = 0;
     remainder[i-1][1] = 0;
     }
   }  
   /* initialization */       
   if(myid == 0)
   {  
     startwtime = MPI_Wtime();
     manager_initialize(ws, NUM_GROUP, stages, numbprocs, filename);
   }
   else 
   {
     slavewtime = MPI_Wtime();
     slave_initialize(myid, NUM_GROUP, filename);
   }
   /* Initialization phase and mainloop initialization are DONE */
   if(myid == 0)
   {
     if((int)fmod(num_group,2)==1) /* there is a remainder */
     { 
       remainder_pos++; 
       remainder[remainder_pos][0] = num_group-1;
       remainder[remainder_pos][1] = 0; 
       printf("record remainder. ws[%d][%d]\n", 
               remainder[remainder_pos][0], remainder[remainder_pos][1]); 
     }
   }
   num_group = (int)floor(num_group/2);
   /* stages */ 
   for(cnt_stage=1;cnt_stage<=stages;cnt_stage++)
   {
     if(myid == 0)
      {
       printf("*********** stage %d **************** \n", cnt_stage); 
       stagewtime = MPI_Wtime();
      } 
     compute_witness_sets(myid,stages,cnt_stage,numbprocs,NUM_GROUP,NUM_VARS,
                          num_group,ws,TaskCount,&cnt_index,filename,
                          &remainder_pos,remainder); 
     MPI_Barrier(MPI_COMM_WORLD);                     
                 
     /* once all groups are finished, clear cnt_index */
     if(myid == 0)
     {
      printf("outside of compute_witness_sets\n");
      printf("remainder_pos=%d\n", remainder_pos);     
      cnt_index = 0; 
      if((int)fmod(num_group,2)==1 && cnt_stage!=stages) 
      /* there is a remainder */
      {  
       remainder_pos++; 
       remainder[remainder_pos][0] = num_group-1;
       remainder[remainder_pos][1] = cnt_stage; 
       printf("record remainder: ws[%d][%d]\n", 
               remainder[remainder_pos][0], 
               remainder[remainder_pos][1]);
       printf("remainder_pos=%d\n", remainder_pos); 
      }
     }  
     num_group = (int)floor(num_group/2); 
     if(num_group==0 && cnt_stage!=stages) num_group=1; 
     /* print out remainders info. */
     if(v>0)
     {
      if(myid == 0)
      {
      if(remainder_pos==1) /* two remainders exist */
       {printf("remainders:\n");
        printf("ws[%d][%d] and ws[%d][%d]\n", 
               remainder[0][0],remainder[0][1],
               remainder[1][0],remainder[1][1]);
       }
      if(remainder_pos==0 && cnt_stage!=stages)
       {
        printf("remainder:\n");
        printf("ws[%d][%d] \n", 
               remainder[0][0],remainder[0][1]);     
       }
      }
     }
     if(myid == 0)
     {  printf("stage %d: %lf seconds\n", cnt_stage, MPI_Wtime()-stagewtime);
        
     }
   }/* cnt_stage, stages */
      
   if(myid == 0)
     wtime = MPI_Wtime()-startwtime;
   else 
     wtime = MPI_Wtime()-slavewtime;
   MPI_Gather(&wtime,1,MPI_DOUBLE,mytime,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      
   if(myid == 0)
    {
     
     print_time(mytime,numbprocs);
     printf("\nTask tallies:\n");
     for(slv=1;slv<=numbprocs-1;slv++)
     {
     total = total + *(TaskCount+slv-1);
     printf("Node %d ----- %d jobs \n",slv, *(TaskCount+slv-1));
     }
     printf("----------------------\n");
     printf("             %d jobs! \n", total); 
   
     /* clean up */
     free(TaskCount);
     free(mytime);
     for(r=0;r<NUM_GROUP;r++)
       for(c=0;c<=stages;c++)
         kill_ws(ws+r*(stages+1)+c);    
    }
     
   adafinal();      
   MPI_Finalize();
   return 0;
}

void combine_points(int *buffer, int *newbuffer)
{
   int eqn_1, dim_1, num_vars_1, solution_1, eqn_2, i;
   
   eqn_1 = buffer[0];
   dim_1 = buffer[buffer[0]+1];
   num_vars_1 = buffer[buffer[0]+2];
   solution_1 = buffer[buffer[0]+3];
   eqn_2 = buffer[buffer[0]+4];
 
   newbuffer[0] = eqn_1 + eqn_2;
   for(i=0;i<=eqn_1-1;i++)
     newbuffer[i+1] = buffer[i+1]; /* wp1 source */
   for(i=0;i<=eqn_2-1;i++)
     newbuffer[eqn_1+i+1] = buffer[eqn_1+5+i]; /* wp2 source */
   newbuffer[eqn_1+eqn_2+1] = dim_1;
   newbuffer[eqn_1+eqn_2+2] = num_vars_1;
   /*newbuffer[eqn_1+eqn_2+3] = solution_1; */
   newbuffer[eqn_1+eqn_2+3] = 0;
}

void pack(int length, JOB *j, int package[length])
{
   int i,num,*newpackage;
   
   newpackage = (int*) malloc((length-3)*sizeof(int));
   
   /* pack witness point 1 */
   package[0] = j->num_eqns_1;
   for(i=1;i<=j->num_eqns_1;i++)
     package[i] = j->source_1[i-1];
   package[j->num_eqns_1+1] = j->dim_1;
   package[j->num_eqns_1+2] = j->num_vars_1;
   /*package[j->num_eqns_1+3] = j->solution_1;*/
   package[j->num_eqns_1+3] = 0;
  
   /* pack witness point 2 */
   if(j->num_eqns_2!=0)
   { 
     num = j->num_eqns_1+4;
     package[num] = j->num_eqns_2;
     for(i=1;i<=j->num_eqns_2;i++)
       package[num+i] = j->source_2[i-1];
     package[num+j->num_eqns_2+1] = j->dim_2;
     package[num+j->num_eqns_2+2] = j->num_vars_2;
     /*package[num+j->num_eqns_2+3] = j->solution_2;*/
     package[num+j->num_eqns_2+3] = 0;
   }
   combine_points(package, newpackage);
   for(i=1;i<=length-3;i++)
     package[i-1]=newpackage[i-1];
   package[length-3] = 0;
   package[length-2] = 0;
   package[length-1] = 0; 
   free(newpackage);
}

LISTELEMENT*  add_jobs (WSET *ws1, WSET *ws2, LISTELEMENT *listpointer)
{
   int i,k,kk;
   JOB *j;
   for(k=0;k<=ws1->deg-1;k++)
     for(kk=0;kk<=ws2->deg-1;kk++)
     {
     j = (JOB *) calloc (1, sizeof(JOB));
     make_job(j,ws1->num_eqns,ws2->num_eqns); 
     /* witness point 1 */
     for(i=0;i<=ws1->num_eqns-1;i++)
       j->source_1[i]=ws1->source[i];
     j->num_eqns_1 = ws1->num_eqns;
     j->num_vars_1 = ws1->num_vars;
     j->dim_1 = ws1->dim;
     j->deg_1 = ws1->deg;
     /*j->solution_1 = ws1->sols[k]; */
     /* witness point 2 */
     for(i=0;i<=ws2->num_eqns-1;i++)
       j->source_2[i]=ws2->source[i];
     j->num_eqns_2 = ws2->num_eqns;
     j->num_vars_2 = ws2->num_vars;
     j->dim_2 = ws2->dim;
     j->deg_2 = ws2->deg;
     /*j->solution_2 = ws2->sols[kk]; */
 
     listpointer = additem(listpointer, j);
     }
   if(v>3) printf("adding jobs is DONE\n");
   return listpointer;
}

void manager_initialize(WSET *ws, int NUM_GROUP, int stages, 
                        int numbprocs, char filename[50])                        
{
   int i,slv,eqnid=0,eqnpack[NUM_GROUP],new_eqnpack[NUM_GROUP],
       quotient,remain,eqnsize,flag,quit_collect=0,cnt_recv;
   MPI_Status status;
       
   /* Initialization phase: sending eqn. index to workers */
   if(NUM_GROUP<=numbprocs-1) /* worker with id<=NUM_GROUP will get one eqn. */
   {
   if(v>3) printf("NUM_GROUP<=numbprocs-1\n");
   for(slv=1;slv<=NUM_GROUP;slv++)
    { 
     eqnid = slv;
     MPI_Send(&eqnid,1,MPI_INT,slv,EQN_TAG,MPI_COMM_WORLD);
     if(v>3 ) printf("manager sent eqn. index %d to %d\n", eqnid, slv);
    }   
   }
   else /* NUM_GROUP>numbprocs-1 */
   {
   if(v>3) printf("NUM_GROUP>numbprocs-1\n");
   quotient = NUM_GROUP/(numbprocs-1);
   remain = (int)fmod(NUM_GROUP,numbprocs-1);   
   if(v>3) printf("quotient=%d,remain=%d\n",quotient,remain);
   eqnid = 1;
   for(slv=1;slv<=remain;slv++)
    {
     for(i=1;i<=quotient+1;i++)
     {
       MPI_Send(&eqnid,1,MPI_INT,slv,EQN_TAG,MPI_COMM_WORLD);
       eqnid++;
     }
    }
   for(slv=remain+1;slv<=numbprocs-1;slv++)
    {
     for(i=1;i<=quotient;i++)
     {
      MPI_Send(&eqnid,1,MPI_INT,slv,EQN_TAG,MPI_COMM_WORLD);
      eqnid++;
     }   
    }
    if(v>3) printf("the last eqnid=%d\n", eqnid-1);
    assert((eqnid-1)<=NUM_GROUP);
   }
   /* send 0 to everyone */
   if(v>3) printf("manager sent 0 to everyone.\n");
   eqnid = 0;
   for(slv=1;slv<=numbprocs-1;slv++)
     MPI_Send(&eqnid,1,MPI_INT,slv,EQN_TAG,MPI_COMM_WORLD);  
 
   /* initial collect */
   cnt_recv = 0;
   while(!quit_collect)
   {
     MPI_Iprobe(MPI_ANY_SOURCE,EQN_TAG,MPI_COMM_WORLD,&flag,&status); 
     while(flag) /* while flag -- pending recv */
     {  
     slv = status.MPI_SOURCE;    /* which slave sent this JOB back */
     MPI_Get_count(&status, MPI_INT, &eqnsize);
     MPI_Recv(new_eqnpack,eqnsize,MPI_INT,slv,EQN_TAG,MPI_COMM_WORLD,&status);
     if(v>3)
     {
      printf("manager recv returned eqn index from node %d with eqnsize=%d\n",
              slv, eqnsize); fflush;
      for(i=0;i<eqnsize;i++)
        printf("%d\n", new_eqnpack[i]);
      fflush;
     }
     cnt_recv++;
     if(cnt_recv==numbprocs-1)  /* all workers have sent back */
       { quit_collect=1; break;}
     MPI_Iprobe(MPI_ANY_SOURCE,EQN_TAG,MPI_COMM_WORLD,&flag,&status);   
     }
   }
   printf("Initialization phase is DONE!\n"); fflush;
   /* End of Initialization phase */
   
   /* mainloop initialization */
   mainloop_initialize(ws, NUM_GROUP, stages, numbprocs, filename);
}

void mainloop_initialize(WSET* ws, int NUM_GROUP, int stages, int numbprocs, char filename[50])
{
   int r,c,i,k,slv,fail,n,dim,deg;
   char *name;
   WSET *ws_p;
   
   /* obtain the witness sets info. from files */
   name = (char *) calloc (50, sizeof(char));
   for(r=0;r<NUM_GROUP;r++)
     for(c=0;c<=stages;c++)
     {   
     ws_p = ws+r*(stages+1)+c;
     if(c==0)
     {
        sprintf(name, "%s_%d", filename, r+1);
        fail = read_witness_set_from_file((int)strlen(name),name,&n,&dim,&deg);
        if(fail>0) printf("failed to read a witness set.\n");
        printf("  n = %d  dimension = %d  degree = %d\n",n,dim,deg);
        fail = solcon_clear_standard_solutions();
        fail = syscon_clear_standard_system( ); 
        make_ws(ws_p,1,n,deg,dim);
        ws_p->source[0] = r+1;
        /*for(k=0;k<ws_p->deg;k++)
        (ws+r*(stages+1)+c)->sols[k] = k+1;*/
        if(v>2) 
        {printf("ws[%d][%d]\n", r,c);}
        }
       else make_ws(ws_p,0,0,0,0);
     }
   free(name);
   printf("Mainloop initialization is DONE!\n"); fflush;
}

void slave_initialize(int id, int NUM_GROUP, char filename[50])
{
   int flag,*buffer,*newbuffer,count,i,eqn_1,dim_1,solution_1,eqn_2,
       cnt_eqns,eqnid,eqn_pack[NUM_GROUP],quit_flag=0,
       fail,k,n;
   char outfile[80];
   MPI_Status status;
   
   if(v>2) printf("Node %d knows the filename: %s\n",id, filename);fflush;
   fail = read_named_target_without_solutions ((int) strlen(filename),filename);
   fail = copy_target_system_to_container();
   
   cnt_eqns = 0;
   while(!quit_flag)
   {
     MPI_Iprobe(0,MPI_ANY_TAG,MPI_COMM_WORLD,&flag,&status); 
     while(flag) /* while flag -- pending recv */
     {  
     if(status.MPI_TAG==EQN_TAG)
     {
     MPI_Recv(&eqnid,1,MPI_INT,0,EQN_TAG,MPI_COMM_WORLD,&status);
     if(eqnid!=0) 
     {
      eqn_pack[cnt_eqns] = eqnid; 
      cnt_eqns++;
      if(v>2) printf("node %d solving eqn. %d; cnt_eqns=%d \n",
                      id, eqnid,  cnt_eqns); fflush;

     /* solving eqn. */
     sprintf(outfile, "%s_%d", filename, eqnid );
     n = (int) strlen(outfile);
     fail = hypersurface_witness_set(eqnid,n,outfile);
     /* end of solving eqn. */
     } 
     else 
     {
     if(v>2) printf("node %d recv end signal and send back: \n", id);
     /* send the index of the eqns for which the worker has computed back */
     eqn_pack[cnt_eqns] = 0; /* the last number is 0 */
     cnt_eqns++;
     if(v>2) 
     {
      printf("cnt_eqns=%d\n", cnt_eqns);
      for(i=0;i<cnt_eqns;i++)
        printf("%d\n", eqn_pack[i]);
     }
     MPI_Send(eqn_pack,cnt_eqns,MPI_INT,0,EQN_TAG,MPI_COMM_WORLD);
     /* clean the system and solution container */
     fail = solcon_clear_standard_solutions( );  
     if(fail>0) printf("fail to clear solution container.\n");
     fail = syscon_clear_standard_system( );
     if(fail>0) printf("fail to clear system container.\n");
     /* end of cleaning the system and solution container */
     quit_flag =1; 
     }
     fflush; 
     } /* status.MPI_TAG */
     MPI_Iprobe(0,MPI_ANY_TAG,MPI_COMM_WORLD,&flag,&status);       
     } /* flag */
   } /* quit_flag */
}

void send_collect(WSET *ws, int num_group, int stages, int cnt_stage, 
                  int numbprocs, int cnt_step, IDLE_ELEMENT **ie, 
                  LISTELEMENT **listpointer, int *TaskCount, int NUM_GROUP, 
                  int *cnt_index, int update_flag, int n1, int n2, int cd)
{
   int slv,flght,flag,*buffer,*package,count,cnt_count,
       idle_procs[numbprocs],idle_pos,cnt_idle,
       r,c,s_max,s_min,
       deg,n=cd,m[2],send[3],label_sol,fail,i=0,j=0;
   double sol1[2*n1+5], sol2[2*n2+5], ps[2*cd+5];
   WSET *ws_p;
   JOB *temp_j;
   MPI_Status status;
   
   if(v>3) printf("inside of send_collect\n"); fflush;
   temp_j = (JOB *)((*listpointer)->dataitem); 
   count = temp_j->num_eqns_1 + temp_j->num_eqns_2 + 8; 
   if(cnt_step==2 && update_flag) 
   {   
     package = (int*) malloc(count*sizeof(int));
     pack(count, temp_j, package);    
   } 
   
   cnt_idle = num_idle(*ie);
   /* Initial JOB distribution  */
   while(*ie != NULL)                 /* there are idle processors */
   {                                     
   if(*listpointer!=NULL)            /* JOB queue is not empty */
     {
     /* pop one JOB from the queue */
     *listpointer = removeitem (*listpointer);
     slv = (*ie)->data;
     *ie = removeslv(*ie);
     /* obtain start solution & send */
     if(v>3) printf("sent a job & path to node %d.\n",slv);
     if(cnt_step==1)
     {   
     if(v>3)
         fail = get_next_start_product
                     (&i,&j,1,temp_j->num_vars_1,temp_j->num_vars_2,
                      temp_j->dim_1,temp_j->dim_2,
                      temp_j->deg_1,temp_j->deg_2,cd,sol1,sol2,ps);
     else
         fail = get_next_start_product
                     (&i,&j,0,temp_j->num_vars_1,temp_j->num_vars_2,
                      temp_j->dim_1,temp_j->dim_2,
                      temp_j->deg_1,temp_j->deg_2,cd,sol1,sol2,ps);
     m[0] = n;  m[1] = 1;
     send[0] = slv; send[1] = m[1]; send[2] = m[0];
     MPI_Send(send,3,MPI_INT,slv,SEND_SMUL,MPI_COMM_WORLD);
     MPI_Send(ps,2*n+5,MPI_DOUBLE,slv,SEND_SSOL,MPI_COMM_WORLD);
     }
     if(cnt_step==2)
     {
      fail = solcon_read_next_solution(n,&m[1],ps);
      if(fail>0) printf("solcon_read_next_solution fail!\n");
      m[0] = n; 
      send[0] = slv; send[1] = m[1]; send[2] = m[0];
      MPI_Send(send,3,MPI_INT,slv,SEND_SMUL,MPI_COMM_WORLD);
      MPI_Send(ps,2*n+5,MPI_DOUBLE,slv,SEND_SSOL,MPI_COMM_WORLD); 
     }
     /* end of obtaining start solution & sending */
     *(TaskCount+slv-1)=*(TaskCount+slv-1)+1;
     } 
   else
     break;
   }
   flght = cnt_idle - num_idle(*ie);
   
   flag = 0;
   while(flght>0 || *listpointer!=NULL)    /* JOB queue loop */
   {    
    if(flght>0)
    {
     MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&flag,&status); 
     while(flag) /* while flag -- pending recv */
       {  
        if(v>3) printf("manager starting recv... \n");   
        slv = status.MPI_SOURCE;    /* which slave sent this JOB back */
        /* recv end solution */   
        MPI_Recv(ps,2*n+7,MPI_DOUBLE,MPI_ANY_SOURCE,SEND_TSOL,
                 MPI_COMM_WORLD,&status);
        m[1] = (int) ps[2*n+5];
        label_sol = (int) ps[2*n+6] - 1;
        fail = solcon_write_next_solution_to_defined_output_file
               (&label_sol,n,m[1],ps);
        /* end of recv. end solution */
        
        /* update idle processor list */
        *ie = addslv(*ie, slv);
        
        /* update corresponding cell when cnt_step==2 && update_flag is TRUE */
        if(cnt_step==2 && update_flag) 
        {
         c = cnt_stage;
         r = *cnt_index;
                  
         ws_p = ws+r*(stages+1)+c; 
        
         if(ws_p->count==0)    /* first returned JOB for current witness set */
         { 
          printf("update num_eqns, source, dim, deg\n");   
          ws_p->num_eqns = package[0];
          ws_p->source = (int*) malloc(package[0]*sizeof(int));
          deg = 1;
          for(cnt_count=0;cnt_count<package[0];cnt_count++)
          {
           ws_p->source[cnt_count] = package[cnt_count+1];
           printf("ws[%d][%d].source[%d]=%d\n",r,c,cnt_count,ws_p->source[cnt_count]);
           deg = (ws+(package[cnt_count+1]-1)*(stages+1))->deg*deg;
          }
          ws_p->dim = package[package[0]+1];
          /*ws_p->sols = (int*) malloc(deg*sizeof(int));*/
          ws_p->deg = deg;
         }  
         /*ws_p->sols[ws_p->count]=ws_p->count+1; */
         ws_p->count++;
        
         if(ws_p->count==ws_p->deg)     /* this witness set is complete */
         {          
          if(ws_p->num_eqns==NUM_GROUP) 
           {
            printf("\nrecord [%d][%d]\n", r,c);
            printf("ALL DONE! aha........\n");
            print_ws(ws_p); 
           }  
         }      /* if: count==deg */
        }       /* cnt_step == 2 && update_flag */               
                  
        if(*listpointer!=NULL)              /* JOB queue is not empty */
        {
          /* pop one JOB from the queue */
          *listpointer = removeitem (*listpointer);
          slv = (*ie)->data;
          *ie = removeslv(*ie); 
          if(v>3) printf("sending a job & path to node %d.\n",slv);         
          /* obtain start solution & send */
          if(cnt_step==1)
          {
           if(v>3)
             fail = get_next_start_product
                     (&i,&j,1,temp_j->num_vars_1,temp_j->num_vars_2,
                      temp_j->dim_1,temp_j->dim_2,
                      temp_j->deg_1,temp_j->deg_2,cd,sol1,sol2,ps);
           else
             fail = get_next_start_product
                     (&i,&j,0,temp_j->num_vars_1,temp_j->num_vars_2,
                      temp_j->dim_1,temp_j->dim_2,
                      temp_j->deg_1,temp_j->deg_2,cd,sol1,sol2,ps);
             m[0] = n;  m[1] = 1;
             send[0] = slv; send[1] = m[1]; send[2] = m[0];
             MPI_Send(send,3,MPI_INT,slv,SEND_SMUL,MPI_COMM_WORLD);
             MPI_Send(ps,2*n+5,MPI_DOUBLE,slv,SEND_SSOL,MPI_COMM_WORLD);
          }
          if(cnt_step==2)
          {
            fail = solcon_read_next_solution(n,&m[1],ps);
            m[0] = n; 
            send[0] = slv; send[1] = m[1]; send[2] = m[0];
            MPI_Send(send,3,MPI_INT,slv,SEND_SMUL,MPI_COMM_WORLD);
            MPI_Send(ps,2*n+5,MPI_DOUBLE,slv,SEND_SSOL,MPI_COMM_WORLD); 
          }   
          /* end of obtaining start solution & sending */
          *(TaskCount+slv-1)=*(TaskCount+slv-1)+1; 
        }
        else                       /* JOB queue is empty */
        {
          flght=flght-1;           /* one in-flight task less */         
          if(v>3)
          {printf ("Job queue empty!\n");        
           printf("flght=%d\n",flght);
          }  
        }
        MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&flag,&status);
       } /* while flag */
    } /* if flght */  
   }  /* while flght listpointer */
   /* send termination to all workers */
   if(v>3) printf("\nmanager sent termination\n");
   for(slv=1;slv<=numbprocs-1;slv++)
   {
    count = -1; 
    MPI_Send(&count,1,MPI_INT,slv,COUNT_TAG,MPI_COMM_WORLD);
   } 
}

void run_slave(int id)
{
   int flag=0,count,i,
       eqn_1,dim_1,solution_1,eqn_2,quit_flag=0,
       m[2],send[2],fail,n;
   double *sol;
   MPI_Status status;
   
   if(v>3) printf("node %d is in run_slave ... \n", id); fflush;
   while(!quit_flag)
   {
    MPI_Iprobe(0,MPI_ANY_TAG,MPI_COMM_WORLD,&flag,&status);
    while(flag) /* while flag -- pending recv */
    {  
     if(status.MPI_TAG==COUNT_TAG)
     {
     MPI_Recv(&count,1,MPI_INT,0,COUNT_TAG,MPI_COMM_WORLD,&status);
     if(v>3) {printf("node %d received count=%d\n", id, count); fflush;} 
     if(count==-1) 
      {
       quit_flag=1; 
       fail = solcon_clear_standard_solutions( );  
       if(fail>0) printf("fail to clear solution container.\n");
       fail = syscon_clear_standard_system( );
       if(fail>0) printf("fail to clear system container.\n");
       if(v>3) printf("node %d terminate!! \n", id); 
       fflush; 
       break;
      }     
     }
     else    /* start receiving jobs */
     {
     if(v>3) printf("node %d starting recv... \n", id); fflush;       
     /* start recv start solution & do path tracking */
     MPI_Recv(send,3,MPI_INT,0,SEND_SMUL,MPI_COMM_WORLD,&status);
     m[1] = send[1]; m[0] = send[2], n = send[2];    
     sol = (double*)calloc(2*n+7,sizeof(double));
     MPI_Recv(sol,2*n+5,MPI_DOUBLE,0,SEND_SSOL,MPI_COMM_WORLD,&status);
     fail = track_one_path(m[0],send[0],&m[1],sol);
     *(sol+2*n+5) = (double) m[1];
     *(sol+2*n+6) = (double) send[0]; 
     /* send result back */
     MPI_Send(sol,2*n+7,MPI_DOUBLE,0,SEND_TSOL,MPI_COMM_WORLD);
     free(sol);
     /* end of recv. start solution & path tracking */
     } /* status.MPI_TAG */
     MPI_Iprobe(0,MPI_ANY_TAG,MPI_COMM_WORLD,&flag,&status);       
    } /* flag */
   } /* while quit_flag */
}

void obtain_namebase(char filename[50], JOB *temp_j, char *namebase, 
                     char *name_1, char*name_2)
{
   char *s_1, *s_2; /* the converted string of numbers */
    
   s_1 = (char *) calloc (3*temp_j->num_eqns_1, sizeof(char));
   s_2 = (char *) calloc (3*temp_j->num_eqns_2, sizeof(char));
   /*printf("in obtain_namebase\n");      */
   obtain_string_from_source(temp_j, filename, s_1, s_2, name_1, name_2);
   sprintf(namebase, "%s_%s%s", filename, s_1, s_2);
   /* clean up */
   free(s_1); free(s_2);   
}

void source_to_string(int* source, int num, char* str) 
{
   int i;
   char temp[3]; 
  
   /*printf("in source_to_string\n");*/
   for(i=0;i<num;i++)
   {
    sprintf(temp, "%d", source[i]);
    str = strcat(str, temp);
   }
   /*printf("str=%s\n", str);*/
}

void obtain_string_from_source(JOB *temp_j, char filename[50], 
                               char *s_1, char *s_2, 
                               char *name_1, char *name_2)
{   
   /*printf("in obtain_string_from_source\n");          */
   source_to_string(temp_j->source_1, temp_j->num_eqns_1, s_1); 
   sprintf(name_1, "%s_%s", filename, s_1);
   /*printf("name_1=%s\n", name_1);*/
   source_to_string(temp_j->source_2, temp_j->num_eqns_2, s_2); 
   sprintf(name_2, "%s_%s", filename, s_2);
   /*printf("name_2=%s\n", name_2);*/
}

/* in fact, user will be prompted for the filename. same as read_two_witness_sets in mpi2track.c */
int read_two_witness_sets_from_file(int *n1, int *n2, int *dim1, int *dim2, 
                                    int *deg1, int *deg2, char *name_1, 
                                    char *name_2, int *cd ) 
{     
   int fail;
         
   /*fail=read_witness_set_from_file((int)strlen(name_1),name_1,n1,dim1,deg1);*/
   fail=read_a_witness_set(1,n1,dim1,deg1);
   printf("  n = %d  dimension = %d  degree = %d\n",*n1,*dim1,*deg1);
         
   /*fail=read_witness_set_from_file((int)strlen(name_2),name_2,n2,dim2,deg2);*/
   fail=read_a_witness_set(2,n2,dim2,deg2);
   printf("  n = %d  dimension = %d  degree = %d\n",*n2,*dim2,*deg2);
   
   fail=extrinsic_top_diagonal_dimension(*n1,*n2,*dim1,*dim2,cd);
   printf("The top dimension of the extrinsic diagonal cascade : %d.\n",*cd);
             
   return fail; 
}

void compute_witness_sets(int myid, int stages, int cnt_stage, int numbprocs,
                         int NUM_GROUP, int NUM_VARS, int num_group, WSET *ws,
                         int *TaskCount, int *cnt_index, char filename[50], 
                         int *remainder_pos, int remainder[2][2])
{
   int cnt_step, expected=0, dim=0, cnt_dim=0, update_flag=0, fail,
       cnt_numgroup, slv, i, remainderjob=0, quit=1, record_flag=0,
       dimension;
   WSET *wspointer_1, *wspointer_2; 
   LISTELEMENT *lp=NULL;
   IDLE_ELEMENT *ie=NULL;
   char *name_1, *name_2, *namebase, *infile, *outfile;
   JOB *temp_j;
   MPI_Status status;
   double groupwtime;
   
   if(myid == 0)
   {
     for(slv=1;slv<=numbprocs-1;slv++)     
        ie = addslv(ie, slv);
     if((int)fmod(num_group,2)==0 && *remainder_pos==1) record_flag=1;  
   }
   while(quit)
   {
   for(cnt_numgroup=1;cnt_numgroup<=num_group;cnt_numgroup++)
   { 
   /* 3 steps: starting of a diagonal homotopy; cascading; collapse; */
   for(cnt_step=1;cnt_step<=3;cnt_step++)
   {
     if(myid == 0) {groupwtime=MPI_Wtime();
                    printf("\n -------- cnt_step=%d --------------\n", cnt_step); }
     else if(v>3) printf("node %d is ready! \n", myid);
     /* for starting a diagonal homotopy */
     
     if(cnt_step==1)
     { if(myid==0) 
       {  
         if(remainderjob!=1)
         {
          printf("start the JOBs for GROUP %d\n", cnt_numgroup); 
          printf("ws[%d][%d]:\n", cnt_numgroup*2-2, cnt_stage-1);
          print_ws(ws+(cnt_numgroup*2-2)*(stages+1)+(cnt_stage-1));
          printf("ws[%d][%d]:\n", cnt_numgroup*2-1, cnt_stage-1);
          print_ws(ws+(cnt_numgroup*2-1)*(stages+1)+(cnt_stage-1));
          wspointer_1 = ws+(cnt_numgroup*2-2)*(stages+1)+(cnt_stage-1);
          wspointer_2 = ws+(cnt_numgroup*2-1)*(stages+1)+(cnt_stage-1);
         }
         if(remainderjob==1 || (cnt_stage==stages && *remainder_pos==1))
         {
          printf("start the JOBs for remainders:\n");
          printf("ws[%d][%d]:\n", remainder[0][0], remainder[0][1]);
          print_ws(ws+remainder[0][0]*(stages+1)+remainder[0][1]);
          printf("ws[%d][%d]:\n", remainder[1][0], remainder[1][1]);
          print_ws(ws+remainder[1][0]*(stages+1)+remainder[1][1]);
          wspointer_1 = ws+remainder[0][0]*(stages+1)+remainder[0][1];
          wspointer_2 = ws+remainder[1][0]*(stages+1)+remainder[1][1]; 
          *remainder_pos = -1;
          remainderjob = 1;
         }
         lp = add_jobs(wspointer_1, wspointer_2, lp); 
         temp_j = (JOB *) (lp->dataitem);  
         namebase = (char *) calloc(100+temp_j->num_eqns_1+temp_j->num_eqns_2,sizeof(char));
         name_1 = (char *) calloc(100+temp_j->num_eqns_1+temp_j->num_eqns_2, sizeof(char));
         name_2 = (char *) calloc(100+temp_j->num_eqns_1+temp_j->num_eqns_2, sizeof(char));                           
         infile = (char *) calloc(100+temp_j->num_eqns_1+temp_j->num_eqns_2, sizeof(char));
         outfile = (char *) calloc(100+temp_j->num_eqns_1+temp_j->num_eqns_2, sizeof(char));
         obtain_namebase(filename, temp_j, namebase, name_1, name_2);
            
         sprintf(outfile, "%s_%d", namebase, cnt_step);
         if(cnt_stage!=1 && remainderjob!=1)
         { name_1 = strcat(name_1, "_3"); name_2 = strcat(name_2, "_3");}
         else
         {
         if(remainderjob==1)
         {
           if(remainder[0][1]!=0)  /* 1st remainder is not from ws[x][0] */
             name_1 = strcat(name_1, "_3");
           /* 2nd remainder is not from ws[xx][0] */
           name_2 = strcat(name_2, "_3"); 
           printf("change remainderjob=0\n");
           remainderjob=0;
         }
         }
         printf("name_1=%s\n", name_1);
         printf("name_2=%s\n", name_2);
         printf("outfile=%s\n", outfile);    
       }
       
       MPI_Barrier(MPI_COMM_WORLD);               
       start_a_diagonal_homotopy(myid, numbprocs, name_1, name_2, outfile, 
                                 ws, num_group, stages, cnt_stage, cnt_step,
                                 &ie, &lp, TaskCount, NUM_GROUP, NUM_VARS, 
                                 cnt_index, &expected, &dim, update_flag);                                
       MPI_Barrier(MPI_COMM_WORLD);       
       /*stack overflow, erroneous access error! WHY????? 
       /*MPI_Bcast(&dim,sizeof(int),MPI_INT,0,MPI_COMM_WORLD);   */
       
       if(myid==0) 
       {
       for(slv=1;slv<=numbprocs-1;slv++)
         MPI_Send(&dim,sizeof(int),MPI_INT,slv,DIM_TAG,MPI_COMM_WORLD); 
       }
       else 
         MPI_Recv(&dim,sizeof(int),MPI_INT,0,DIM_TAG,MPI_COMM_WORLD,&status);  
       MPI_Barrier(MPI_COMM_WORLD);  
       if(v>3) printf("node %d knows dim=%d\n", myid, dim);  
       
       if(myid==0) 
       { printf("dim=%d\n", dim);
         if(v>1)
         {            
         printf("outside of start_a_diagonal_homotopy\n");
         printf("expected dimension=%d\n", expected);
         printf("dim=%d\n", dim);
         printf("%d JOBs in lp\n", length(lp));
         printf("%d idles\n", num_idle(ie));
         printf("indexing........ with cnt_index=%d\n", *cnt_index);
         }
       }
       else if(v>1) {printf("node %d finishes step 1 \n", myid); fflush;}
     } /* cnt_step==1 */
          
     MPI_Barrier(MPI_COMM_WORLD);
     /* for one-level-down cascade homotopy */
     if(cnt_step==2)
     {
       if(myid==0) 
       {
        /* the state table is updated only at cnt_step==2 */
        if(dim==1) update_flag=1;
        lp = add_jobs(wspointer_1, wspointer_2, lp);   
        printf("namebase=%s\n", namebase);                         
        sprintf(infile, "%s_%d", namebase, cnt_step-1);
        sprintf(outfile, "%s_%d", namebase, cnt_step);
        printf("infile=%s\n", infile);
        printf("outfile=%s\n", outfile);
       }    
       MPI_Barrier(MPI_COMM_WORLD);
       
       cascade_one_level_down(myid, numbprocs, infile, outfile, 
                                ws, num_group, stages, cnt_stage, cnt_step,
                                &ie, &lp, TaskCount, NUM_GROUP, &dimension,
                                cnt_index, update_flag);  
                                                              
       
       /* dim => (dim-1) times of {remove last slack variable and cascade} */
       if(dim>1)
       {
        for(cnt_dim=1;cnt_dim<=dim-1;cnt_dim++)
        {
         if(myid==0)
         {
          /*the state table need to be updated when the last cascade is done. */
          if(cnt_dim==dim-1) update_flag=1; 
          printf("removing last slack variable \n");
          sprintf(infile, "%s", outfile);
          sprintf(outfile, "%s%s", outfile, "r"); 
          printf("infile=%s\n", infile);
          printf("outfile=%s\n", outfile);
          
          remove_last_slack_variable(infile, outfile); 
          printf("dimension=%d\n", dimension);
          fail = syscon_remove_symbol_from_table(dimension);
          
          sprintf(infile, "%s", outfile);
          sprintf(outfile, "%s%d", outfile, cnt_step);
          printf("%d times cascade_one_level_down\n", cnt_dim+1);
          printf("infile=%s\n", infile);
          printf("outfile=%s\n", outfile);  
           
          lp = add_jobs(wspointer_1, wspointer_2, lp);                             
         }
         MPI_Barrier(MPI_COMM_WORLD);
         cascade_one_level_down(myid, numbprocs, infile, outfile, 
                                ws, num_group, stages, cnt_stage, cnt_step,
                                &ie, &lp, TaskCount, NUM_GROUP, &dimension,
                                cnt_index, update_flag);                                                            
        }
       } /* dim>1 */
     } /* cnt_step==2 */
     MPI_Barrier(MPI_COMM_WORLD);
     
     /* for collapsing the extrinsic diagonal */
     if(cnt_step==3)
     {  /* only manager does collapsing */
      if(myid==0) 
        {           
         sprintf(infile, "%s", outfile);
         sprintf(outfile, "%s_%d", namebase, cnt_step);
         printf("infile=%s\n", infile);
         printf("outfile=%s\n", outfile);
         
         collapse(cnt_stage, stages, 
                  infile, outfile, 
                  ws, cnt_index, 
                  &expected);
          
         if(v>1)
         {            
          printf("outside of collapse\n");
          printf("%d idles\n", num_idle(ie));
          printf("indexing........ with cnt_index=%d\n", *cnt_index);
         }
        } /* myid==0 */
     } /* cnt_step==3 */ 
     
     MPI_Barrier(MPI_COMM_WORLD);
     if(myid==0) update_flag=0;
                                
   } /* for cnt_step */
   if(myid == 0) 
    {
     free(namebase); free(name_1); free(name_2); free(infile); free(outfile);
     printf("group %d: %lf seconds ", cnt_numgroup, MPI_Wtime()-groupwtime); 
    }
   }/* for cnt_numgroup */ 
   if(myid == 0)
   {
   if(*remainder_pos==1)
    {
     num_group=1;
     remainderjob=1;
    }
   else 
    {
     quit=0; 
     if(record_flag)
     {
       printf("in compute_witness_sets\n");
       printf("*remainder_pos=%d\n", *remainder_pos);
       (*remainder_pos)++;                          
       remainder[*remainder_pos][0] = *cnt_index-1;
       remainder[*remainder_pos][1] = cnt_stage;  
       printf("record remainder: ws[%d][%d]\n",  
               remainder[*remainder_pos][0],      
               remainder[*remainder_pos][1]);  
       printf("*remainder_pos=%d\n", *remainder_pos);
     }
    }   
   }
   MPI_Bcast(&num_group,sizeof(int),MPI_INT,0,MPI_COMM_WORLD); 
   MPI_Bcast(&quit,sizeof(int),MPI_INT,0,MPI_COMM_WORLD); 
   MPI_Barrier(MPI_COMM_WORLD);
   }
}

void start_a_diagonal_homotopy(int myid, int numbprocs, char *name_1, 
                               char *name_2, char *outfile, WSET *ws,
                               int num_group, int stages, int cnt_stage, 
                               int cnt_step,
                               IDLE_ELEMENT **ie, LISTELEMENT **listpointer,
                               int *TaskCount, int NUM_GROUP, int NUM_VARS, 
                               int *cnt_index, int *expected,
                               int *dim, int update_flag)                           
{
   int n, deg, n1, n2, dim1, dim2, deg1, deg2, cd, nbsols,
       fail, i, slv, count;
   double vcp[34];
  
   if(myid==0)
   {
     printf("manager is in inside of start_a_diagonal_homotopy\n"); 
                           
     fail = read_two_witness_sets_from_file(&n1,&n2,&dim1,&dim2,
                                            &deg1,&deg2,name_1,name_2,&cd);
     n = cd; /* dimension */
     nbsols = deg1*deg2;
     printf("#paths to track: %d\n", nbsols);
     
     fail = define_output_file_with_string ((int)strlen(outfile), outfile);
     
     fail = standard_diagonal_homotopy(dim1,dim2);
     if(fail == 0)
     {
      fail = write_standard_target_system();
      fail = write_standard_start_system();
     }
     /*fail = tune_continuation_parameters(); printf("\n");*/
     fail = retrieve_continuation_parameters(vcp);
     write_solution_banner_to_defined_output_file(nbsols,n);
     printf("\nSee the output file %s for results...\n\n", outfile);  
     *expected = dim1+dim2-NUM_VARS;  
     printf("expected dimension=%d\n", *expected);
   }       
   
   dimension_broadcast(myid,&n);
  
   MPI_Bcast(&nbsols,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(vcp,34,MPI_DOUBLE,0,MPI_COMM_WORLD);
   
   if(myid != 0) 
   {
      fail = set_continuation_parameters(vcp);
   }
   
   fail = homotopy_broadcast(myid,n);
   fail = create_homotopy_with_given_gamma
                (-0.824263733224601,0.566206056211557);
   
   if(myid == 0)
   {
     send_collect(ws, num_group, stages, cnt_stage, numbprocs, cnt_step, 
                  ie, listpointer, TaskCount, NUM_GROUP,
                  cnt_index, update_flag, 
                  n1, n2, cd);
   
     if(v>3)
     {   
     printf("after send_collect \n");
     printf("%d JOBs in listpointer. \n", length(*listpointer));
     printf("%d idles\n", num_idle(*ie));
     printf("indexing........ with cnt_index=%d\n", *cnt_index);    
     }                 
   }
   else
   {   
     if(v>3) printf("node %d will run run_slave \n", myid); fflush;
     run_slave(myid); 
   }
   
   MPI_Barrier(MPI_COMM_WORLD);
                            
   if(myid == 0)
   { 
     if(v>3) printf("manager clear data \n"); 
     fail = clear_data(); 
     if(fail>0) printf("manager fail to clear data.\n"); 
     fail = clear_homotopy(); 
     if(fail>0) printf("manager fail to clear homotopy.\n"); 
     fail = solcon_close_solution_input_file(1);
     if(fail>0) printf("fail to close witness set 1.\n");
     fail = solcon_close_solution_input_file(2); 
     if(fail>0) printf("fail to close witness set 2.\n");
     fail = solcon_clear_standard_solutions( );  
     if(fail>0) printf("fail to clear solution container.\n");
     fail = close_output_file();
     if(fail>0) printf("fail to close output file. \n"); 
     fail = read_witness_set_from_file((int)strlen(outfile),
                                       outfile,&n,dim,&deg); 
     fail = solcon_clear_standard_solutions();
     fail = syscon_clear_standard_system( );
     printf("end of start_a_diagonal_homotopy\n");
     fflush;
   }
   else 
   {
     if(v>3) printf("node %d clear data.\n", myid); 
     fail = clear_data();
     if(fail>0) printf("node %d fail to clear data.\n", myid);
     fail = clear_homotopy();
     if(fail>0) printf("node %d fail to clear homotopy.\n", myid);
   }
}

void cascade_one_level_down(int myid, int numbprocs, char *infile, 
                            char *outfile, WSET *ws, int num_group, 
                            int stages, int cnt_stage, int cnt_step, IDLE_ELEMENT **ie, 
                            LISTELEMENT **listpointer, int *TaskCount,
                            int NUM_GROUP, int *dimension, int *cnt_index, 
                            int update_flag)
{
   int n,/* dimension */
       nbsols, dim, deg, fail, i,
       n1, n2, cd;   
   double vcp[34];
  
   if(myid == 0)
   {
     printf("manager is in inside of cascade_one_level_down\n");        
     read_dimension_of_system((int)strlen(infile),infile,&n);
   }
   
   dimension_broadcast(myid,&n);
   
   if(myid == 0)
   {
     fail = syscon_clear_symbol_table();
     fail = read_named_start_without_solutions((int)strlen(infile),infile);  
     /*fail = copy_start_system_to_container();*/
     fail = copy_start_system_to_container();
     fail = syscon_sort_embed_symbols(&dim);
     printf("the top dimension is %d\n",dim);
     fail = copy_container_to_start_system();
      
     fail = solcon_scan_solution_banner();
     fail = solcon_read_solution_dimensions(&nbsols,dimension); 
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
   
   fail = create_cascade_homotopy(); 
   
   MPI_Barrier(MPI_COMM_WORLD);
   if(myid == 0)
   {
      if(v>3) printf("# paths to track : %d\n",nbsols);
      fail = define_output_file_with_string ((int)strlen(outfile), outfile);
      
      fail = write_standard_target_system();
      fail = write_standard_start_system(); 
      /*fail = tune_continuation_parameters(); printf("\n");*/
      fail = retrieve_continuation_parameters(vcp);     
      write_solution_banner_to_defined_output_file(nbsols,n);
      printf("\nSee the output file for results...\n\n");
   }
   
   MPI_Bcast(&nbsols,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(vcp,34,MPI_DOUBLE,0,MPI_COMM_WORLD);
    
   if(myid != 0) 
   {
      fail = set_continuation_parameters(vcp); 
   } 
   if(myid == 0)
   {
      n1=0; n2=0; cd=n;
      send_collect(ws, num_group, stages, cnt_stage, numbprocs, cnt_step, 
                   ie, listpointer, TaskCount, NUM_GROUP,
                   cnt_index, update_flag, 
                   n1, n2, cd);  
      if(update_flag) *cnt_index=*cnt_index+1;
      
      if(v>3)
      {                         
       printf("after send_collect \n");
       printf("%d JOBs in listpointer. \n", length(*listpointer));
       printf("%d idles\n", num_idle(*ie));
       printf("indexing........ with cnt_index=%d\n", *cnt_index);  
      }                             
   }
   else
   {  
      if(v>3) printf("node %d will run run_slave \n", myid); fflush;
      run_slave(myid); 
   }
      
   MPI_Barrier(MPI_COMM_WORLD);
                          
   if(myid == 0)
   { 
      if(v>3) printf("manager clear data \n"); 
      fail = clear_data(); 
      if(fail>0) printf("manager fail to clear data.\n"); fflush;
      fail = solcon_close_solution_input_file(0);
      if(fail>0) printf("fail to close solution input file.\n");
      fail = solcon_clear_standard_solutions( );  
      if(fail>0) printf("fail to clear solution container.\n");
      fail = close_output_file();
      if(fail>0) printf("fail to close output file. \n");
      fflush;	 
   } 
   else 
   {
      if(v>3) printf("node %d clear data \n", myid); 
      fail = clear_data(); 
      if(fail>0) printf("node %d fail to clear data.\n", myid); 
      fflush; 
   }
}

int remove_last_slack_variable(char *infile, char *outfile)
{
   int n,dim,deg,fail;
   
   fail = read_witness_set_from_file((int)strlen(infile),infile,&n,&dim,&deg);
   printf("\nThe current number of slack variables : %d\n",dim);
   fail = remove_last_slack(dim);
   fail = write_witness_set_to_file((int)strlen(outfile),outfile);
   fail = solcon_clear_standard_solutions();
   fail = syscon_clear_standard_system( ); 
   return fail;
}

void collapse(int cnt_stage, int stages, 
              char *infile, char *outfile,
              WSET *ws, int *cnt_index, 
              int *expected)
{
   int n, dim, deg, add, fail, i, cnt_dim;
     
   /* check which input file to use */
   fail = read_witness_set_from_file((int)strlen(infile),infile,&n,&dim,&deg);
   printf("n=%d, dim=%d, deg=%d\n", n, dim, deg);
         
   add = *expected-dim; 
   printf("add=%d\n", add);
   fail = standard_collapse_diagonal(dim,add);
   fail = write_witness_set_to_file((int)strlen(outfile),outfile);
     
   /* update the value of num_vars in the state table */
   fail = read_witness_set_from_file((int)strlen(outfile),outfile,&n,&dim,&deg);
   printf("update ws[%d][%d]\n", *cnt_index-1, cnt_stage);
   (ws+(*cnt_index-1)*(stages+1)+(cnt_stage))->num_vars = n;
   (ws+(*cnt_index-1)*(stages+1)+(cnt_stage))->dim = dim; 
      
   fail = solcon_close_solution_input_file(0);
}
