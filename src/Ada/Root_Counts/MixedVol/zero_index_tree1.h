/* The file "zero_index_tree1.h" defines the data structures for the specific
 * case of level 0 in the index tree of LP problems. */

#ifndef _L0_IML1_H_
#define _L0_IML1_H_

#include <memory.h>
#include "index_tree_lp1.h"

typedef struct L0IdxNode L0IdxNode;
typedef struct L0_IML L0_IML;

struct L0IdxNode      /* structure to define linked list */
{
   int idx;           /* index of the node in the list is just a number */
   IndexNode *R;      /* linked list of LPdata */
   L0IdxNode *D;      /* points to the next item in the list */
};

struct L0_IML         /* index tree of LP problems at level zero */
{
   L0IdxNode *L0head; 
   L0IdxNode *L0curr; /* Current IT node ptr, for Add(-). */
   L0IdxNode *L0prev; /* parent of curr, for Add(-). */
   IndexNode *curr;   /* Current IT node ptr, for Add(-). */
   IndexNode *prev;   /* parent of curr, for Add(-). */
   LPPL *LP1;         /* dummy head of the LP address link (used in level 1)*/
};

void L0IdxNode_Init ( L0IdxNode *p, int i );
/*
 * DESCRIPTION :
 *   Initializes the linked list of nodes with LPdata.
 *
 * ON ENTRY :
 *   p         memory allocated for one structure of type L0IdxNode;
 *   i         value for p->idx.
 *
 * ON RETURN :
 *   p         p->idx == i and the pointers p->R and p->D are 0. */

void L0_IML_Init ( L0_IML *p );
/*
 * DESCRIPTION :
 *   Allocates memory for p->L0head and p->LP1, initializing both.
 *   The pointers p0->L0curr and p0->L0prev are set to p->L0head. */

int L0_Migrate ( L0_IML *p, IndexNode *inp );
/*
 * DESCRIPTION :
 *   Migrates the data of the index tree at level zero with the node inp.
 *
 * ON RETURN :
 *   p          p->L0head->D is updated and p->L0prev is cleared;
 *   inp        inp->S now has the value of p->L0prev->R;
 *   0 if p->L0head->D is empty, 1 otherwise.
 *   the list in p->L0prev->R is assigned to inp->S. */

int L0_FindInR ( L0_IML *p, int IDX );
/*
 * DESCRIPTION :
 *   Does horizontal search in p for a node with label idx equal to IDX.
 *
 * ON ENTRY :
 *   p         index tree with LP problems at level zero;
 *   IDX       label to an index node.
 *
 * ON RETURN :
 *   0 if there is no index node in the list starting with p->prev->S
 *     with the given label IDX, otherwise
 *   1 is returned if there is an index node with this index
 *     and then p->curr->idx == IDX. */

int L0_FindInD ( L0_IML *p, int IDX );
/*
 * DESCRIPTION :
 *   Does vertical search in p for a node with label idx equal to IDX.
 *
 * ON ENTRY :
 *   p         index tree with LP problems at level zero;
 *   IDX       label to an index node.
 *
 * ON RETURN :
 *   0 if there is no index node in the list starting with p->L0prev->D
 *     with the given label IDX, otherwise
 *   1 is returned if there is an index node with this index
 *     and then p->L0curr->idx == IDX. */

void L0_Add1
 ( L0_IML *p, int n, int *J, int N, int *I, double *X, double **A );
/*
 * DESCRIPTION :
 *   Updates the index tree at level zero with data for a LP problem.
 *
 * REQUIRED :
 *   The n entries in J are sorted in ascending order.
 *
 * ON ENTRY :
 *   p         current index tree at level zero of LP problems.
 *   n         dimension of the vector J;
 *   J         index vector of n entries, sorted in ascending order;
 *   N         dimension of the LP problem;
 *   I         index vector of dimension N with constraints involved;
 *   X         solution vector of dimension N;
 *   A         basis inverse, matrix of dimension N.
 *
 * ON RETURN :
 *   p         updated index tree at level zero of LP problems. */

void L0_Add2 
 ( L0_IML *p, int *J, int N, int *I, double *X, double **A );
/*
 * DESCRIPTION :
 *   Updates the index tree at level zero with data for a LP problem,
 *   for only two points.
 *
 * REQUIRED :
 *   The two entries in J are sorted in ascending order.
 *
 * ON ENTRY :
 *   p         current index tree at level zero of LP problems;
 *   J         two number sorted in ascending order;
 *   N         dimension of the LP problem;
 *   I         index vector of dimension N with constraints involved;
 *   X         solution vector of dimension N;
 *   A         basis inverse, matrix of dimension N.
 *
 * ON RETURN :
 *   p         updated index tree at level zero of LP problems. */

void L0_IML_Del ( L0_IML *p );
/*
 * DESCRIPTION :
 *   Deletes the index tree at level zero,
 *   calls the routine L0_free below. */

void L0_Free ( L0_IML *li );
/*
 * DESCRIPTION :
 *   Deallocation of all the LPdata in the index tree,
 *   called by the other routine L0_IML_Del. */

#endif
