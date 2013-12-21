/* Implementation of job_queue.h  */

#include "job_queue.h"
#include <stdlib.h>
#include <stdio.h>

void make_ws(WSET *ws, int num_eqns, int num_vars, int deg, int dim) 
{
    if(num_eqns != 0)
      ws->source = (int*) calloc(num_eqns,sizeof(int));
    /*if(deg != 0)
       ws->sols = (int*) calloc(deg,sizeof(int));*/
    ws->num_eqns = num_eqns; 
    ws->num_vars = num_vars;
    ws->deg = deg;
    ws->dim = dim; 
    ws->count = 0;
}

void print_ws(WSET *ws)
/* prints out the witness sets */
{
   int i;

   printf("Source:");
   for(i=0;i<=ws->num_eqns-1;i++)
      printf("%d ",ws->source[i]);
   printf("\ndim=%d",ws->dim);
   printf("\nnum_vars=%d",ws->num_vars);
   /*printf("\nSolutions:");
   for(i=0;i<=ws->deg-1;i++)
      printf("%d ",ws->sols[i]); */
   printf("\ncount=%d\n",ws->count);
}

void kill_ws(WSET *ws) 
{
   if(ws->num_eqns != 0) 
   {
     free(ws->source);
     ws->source = NULL;    
     ws->num_eqns = 0;
   }
   if(ws->deg != 0)
   {
   /*free(ws->sols);
   ws->sols = NULL;  */
   ws->deg = 0; 
   }
   ws->dim = 0;
   ws->num_vars = 0;
}

void make_job(JOB *j, int num_eqns_1, int num_eqns_2) 
{
   if(num_eqns_1!=0) 
   j->source_1 = (int*) calloc(num_eqns_1,sizeof(int));
   j->num_eqns_1 = num_eqns_1;
   j->num_vars_1 = 0;
   j->dim_1 = 0;
   j->deg_1 = 0;
   /*j->solution_1 = 0; */
   if(num_eqns_2!=0)
     j->source_2 = (int*) calloc(num_eqns_2,sizeof(int));
   j->num_eqns_2 = num_eqns_2;
   j->num_vars_2 = 0;
   j->dim_2 = 0;
   j->deg_2 = 0;
   /*j->solution_2 = 0;*/
}

void kill_job(JOB *j) 
{
   free(j->source_1);
   j->source_1 = NULL;
   free(j->source_2);
   j->source_2 = NULL;
   j->num_eqns_1 = 0;
   j->num_eqns_2 = 0;
   j->num_vars_1 = 0;
   j->num_vars_2 = 0;	
}

LISTELEMENT* additem (LISTELEMENT *listpointer, JOB *data) 
{
   LISTELEMENT *lp = listpointer;

   if (listpointer != NULL) 
   {
     while (listpointer -> link != NULL)
       listpointer = (LISTELEMENT *)listpointer -> link;
     listpointer -> link = (LISTELEMENT  *) malloc (sizeof (LISTELEMENT));
     listpointer = (LISTELEMENT *)listpointer -> link;
     listpointer -> link = NULL;
     listpointer -> dataitem = (JOB *)data;
     return lp;
   }
   else 
   {
     listpointer = (LISTELEMENT  *) malloc (sizeof(LISTELEMENT));
     listpointer -> link = NULL;
     listpointer -> dataitem = (JOB *)data;
     return listpointer;
   }
}

LISTELEMENT* removeitem (LISTELEMENT *listpointer) 
{
   LISTELEMENT * tempp;
   /*JOB *temp_j = (JOB *)(listpointer->dataitem);
   printf ("Element removed is solution_1=%d, solution_2=%d\n", temp_j->solution_1, temp_j->solution_2); */
   tempp = (LISTELEMENT *)listpointer -> link;
   free (listpointer);
   return tempp;
}

int length (LISTELEMENT *listpointer)
{
   if(listpointer != NULL)
     return 1+length((LISTELEMENT *)listpointer->link);
   else
     return 0;
}

void printqueue (LISTELEMENT *listpointer) {
   int i,j;
   JOB *temp_j;
   if (listpointer == NULL)
     printf ("queue is empty!\n");
   else
   {
     j=1;
     while (listpointer != NULL) 
     {
       temp_j = (JOB *)(listpointer->dataitem);
     printf("\nNo. %d JOB\n", j);
     printf("witness point 1:\n");
     for(i=0;i<=temp_j->num_eqns_1-1;i++)
     {
       printf("%d ",temp_j->source_1[i]);
     }
     printf ("\n%d", temp_j->dim_1);
     printf ("\n%d", temp_j->num_vars_1);
     /*printf ("\n%d", temp_j->solution_1);*/
     if(temp_j->num_eqns_2 != 0)
     {
       printf("\nwitness point 2:\n");
       for(i=0;i<=temp_j->num_eqns_2-1;i++)
       {
        printf("%d ",temp_j->source_2[i]);
       }
       printf ("\n%d", temp_j->num_vars_2);
       printf ("\n%d", temp_j->dim_2);
       /*printf ("\n%d", temp_j->solution_2);*/
     }
        
     listpointer = (LISTELEMENT *)listpointer -> link;
     j++;
     }
   }
   printf ("\n");
}

void clearqueue (LISTELEMENT **listpointer) 
{
   while (*listpointer != NULL) 
   {
     *listpointer = removeitem (*listpointer);
   }
}
