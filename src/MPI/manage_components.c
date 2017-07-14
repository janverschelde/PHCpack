/* The file "manage_components.c" contains the definitions of the
 * protypes in manage_components.h". */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "manage_components.h"

#define v 1  /* verbose flag:
                  0 no output during the computations,
                  1 only one line messages on progress of computations
                  2 display of data, like systems and solutions 
		  3 low-level debugging messages 
	     */

int SWAP(int* X, int* Y)
{ 
   int T = *X;
   *X = *Y;
   *Y = T;
}  

void makeCOMPONENT(COMPONENT* c, int n, int s)
{
   c->num = n;
   c->size = s;
   c->pts = (int*)calloc(s,sizeof(int));
   c->is_irred = 0;
}

void killCOMPONENT(COMPONENT* c)
{
   free(c->pts);
   c->pts = NULL;
   c->size = 0;
}

bool is_in_COMPONENT(int a, COMPONENT* c)
{
   int i;
   for (i=0; i<c->size; i++)
   {
      if (c->pts[i]==a) return 1;
   }
   return 0;
}

void mergeCOMPONENT(COMPONENT* a, COMPONENT* b, int* ind)
{
   COMPONENT t = *a;
   int i;
   if (v==2) printf("Merging components %d and %d\n", a->num, b->num);
  
  /* the elements of b now belong to a */
   for (i=0; i<b->size; i++)
   {
      ind[b->pts[i]] = a->num;
   }

   a->size = t.size+b->size;
   a->n_nodes_tracking = t.n_nodes_tracking + b->n_nodes_tracking;
   a->pts = (int*)calloc(a->size, sizeof(int));
   memcpy(a->pts, t.pts, t.size*sizeof(int));
   memcpy(a->pts+t.size, b->pts, b->size*sizeof(int));
   killCOMPONENT(&t);
   killCOMPONENT(b);
}

int free_slave(SLAVE_STATE* state, int np)
{
   int sl;
   for (sl=1; sl<np; sl++)
   {
      if (state[sl].status == FREE) return sl;
   }
   return 0;
}

bool isGood(ISLICE_STATS* s) 
{
   return s->n_paths < MIN_ATTEMPTS || 
            100*s->n_rejected < MAX_PERCENT_REJECTED*s->n_paths;
}   

void print_component(COMPONENT* c)
{
  int j;
  if(c->is_irred) printf("["); else printf("(");
  for (j=0; j<c->size; j++)
    {
      printf("%d", c->pts[j]);
      if(j!=c->size-1) printf(" ");
    }
  if(c->is_irred) printf("]"); else printf(")");
}
 
void print_partition(int n, COMPONENT* part)
{
   int i,c = 0;
   for(i=0; i<n; i++)
   {
      if(part[i].size != 0)
      {
         c++;
	 print_component(&part[i]);
      }
   }
   printf(", #comp = %d\n", c);
}

bool is_in_slice(int pt, int slice)
{ 
   return 1;
}

bool is_irreducible(COMPONENT* c) 
{ 
   int fail;
   double d;
   fail = standard_trace_sum_difference(c->size, c->pts, &d);
   if (v>=3) printf("Comp #%d: has trace sum difference %.2e\n",c->num,d);
   return d<TRACE_TOLERANCE;
}

void record_path_outcome(ISLICE_STATS* s, int o)
{
  s->path_outcome[(s->n_paths++)%N_PATHS_TO_REMEMBER] = o;
  if (o == PATH_DIVERGED || o == PATH_STAYED_WITHIN_COMPONENT)
    s->n_rejected++;
}

int MIN(int X, int Y) {return (X>Y) ? Y : X;}
