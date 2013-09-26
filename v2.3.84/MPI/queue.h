
/* Declaration of a queue as a list */

#ifndef QUEUE_H
#define QUEUE_H

typedef struct 
{
  int *start_pivot;
  int *target_pivot;
  int length;
  int address;
  double *start_sol;
}job;

typedef struct Queue queue;

struct Queue 
{
   job   *item;  /* all of the roots' information  */
   queue *next;  /* points to next node */
}; 

void Allocate_job ( job *j, int p, int length );
/* Allocates memory for a job in the queue */

void Deallocate_job ( job *j );
/* Deallocates memory for a job in the queue */

int length ( queue *q ); 
/* Returns the number of elements in the queue */

void append( queue **front, queue **end, job* item );
/* Appends a job to the end of the queue */

queue* pop( queue *q, job **job_ptr );
/* Returns a job from the front of the queue*/

#endif
