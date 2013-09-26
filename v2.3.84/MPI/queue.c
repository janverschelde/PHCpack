/* Implementation of a queue as a list */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "queue.h"
#include "parallel_tree.h"

void Allocate_job ( job *j, int p, int length )
{

  j->start_pivot = (int*) calloc(2*p, sizeof(int));
  j->target_pivot = (int*) calloc(2*p,sizeof(int));
  j->start_sol = (double*) calloc(length, sizeof(double));
}

void Deallocate_job ( job *j )
{
  free(j->start_pivot);
  free(j->target_pivot);
  free(j->start_sol);
  free(j);
}

int length ( queue *q )
{
   if (q==NULL)
      return 0;
   else 
      return 1+length(q->next);
} 

void append( queue **q, queue **end, job *j )
{
   queue* nq = (queue*) calloc(1, sizeof(queue));
 
   nq->item = j;

   if(*q==NULL)
   {
      nq->next = NULL;
      *end = nq;
      *q = nq;
   }
   else
   {
     (*end)->next = nq;  
     nq->next = NULL;
     *end = nq;    
   }
}

queue* pop( queue*q, job **job_ptr )
{
   queue *nq;

   assert(q!=NULL);
   *job_ptr = q->item; 
   nq = q->next;
   free(q);  

   return nq;
}

