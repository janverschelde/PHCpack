#ifndef JOB_QUEUE_H
#define JOB_QUEUE_H

typedef struct w_set WSET;
typedef struct job JOB;
typedef struct listelement LISTELEMENT;

struct w_set {
   int *source;
   int num_eqns;
   int num_vars;
   int dim;
   /*int *sols;*/
   int deg;
   int count;     /* number of sols has been finished */
};

struct job {
   int *source_1;
   int num_vars_1;
   int num_eqns_1;
   int dim_1;
   int deg_1;
   /*int solution_1;*/
   int *source_2;
   int num_vars_2;
   int num_eqns_2;
   int dim_2;
   int deg_2;
   /*int solution_2;*/
};

struct listelement {
   JOB *dataitem;
   LISTELEMENT *link;
};

void make_ws(WSET *ws, int num_eqns, int num_vars, int deg, int dim);
/* Create a witness set with the content of the parameters. */

void print_ws(WSET *ws);
/* prints out the witness sets */

void kill_ws(WSET *ws);
/* Destroy the witness set ws. */

void make_job(JOB *j, int num_eqns_1, int num_eqns_2);
/* Create a JOB with the content of the parameters. */

void kill_job(JOB *j);
/* Destroy the JOB j. */

LISTELEMENT* additem (LISTELEMENT *listpointer, JOB *data);
/* Adds a JOB to the queue. */

LISTELEMENT* removeitem (LISTELEMENT *listpointer);
/* Returns the pointer to the next available JOB in the queue. */

int length (LISTELEMENT *listpointer);
/* Returns the number of JOBs in the queue. */

void printqueue (LISTELEMENT *listpointer);
/* Print out all JOBs in the queue. */

void clearqueue (LISTELEMENT **listpointer);
/* Clear the queue. */

#endif


