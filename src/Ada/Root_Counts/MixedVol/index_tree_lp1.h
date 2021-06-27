/* The file "index_tree_lp1.h" defines the general data structures for an
 * index tree of LP problems, with unused level 0, treated separtedly by
 * the data structures and prototypes in zero_index_tree1.h. */

#ifndef _IT_LP1_H_
#define _IT_LP1_H_

#include <memory.h>

typedef struct LPdata LPdata;
typedef struct LPPL LPPL;
typedef struct IndexNode IndexNode;
typedef struct IT_LP IT_LP;

struct LPdata        /* LP data = indices to active constraints, */
{                    /*           feasible solution, and basis inverse */
   int dim;          /* size of the vectors xxx, JJJ, and matrix INV */
   double *xxx;      /* vector of range 0..dim-1, solution vector */
   double **INV;     /* matrix of dimension dim, basis inverse */
   int *JJJ;         /* index vector for the solution to the LP problem */
};   
  
struct LPPL          /* LP pointer linked list */
{
   LPdata *addr;     /* information of one node in the list */
   LPPL *next;       /* pointer to the next item in the list */
};

struct IndexNode     /* node in linked list of LPdata */
{
   int idx;          /* index of the node is just a number */
   LPdata *info;     /* information for an LP problem */
   IndexNode *S;     /* points to next element in linked list */
};

struct IT_LP         /* index tree for LP problems */
{
   int MaxLevels;    /* maximal #levels for IT[] */
   int CurLevel;     /* index to the current level */
                     /* the next 5 items are arrays of dimension MaxLevels */
   int *DIM;         /* dimension of LP in each level */
   int *NP;          /* NP[L] = #nodes in IT[L], including fixed node. */
   int *cell;        /* indices in horizontal branch of tree when full */
   int *InSpt;       /* indices to the supports */
   int *minNP;       /* minimal #nodes in IT[L], excluding fixed node. */

   LPPL **LP;        /* array of dimension MaxLevels of LP problems */
   LPPL *LPlast;     /* points to the last valid node in LP[CurLevel] */

   IndexNode **IT;   /* array of size MaxLevels, IT[0] is not used */
   IndexNode **last; /* last[] points to the last valid node in IT[] */
   IndexNode *curr;  /* pointer to the current IT node */
   IndexNode *prev;  /* pointer to the previous IT node */
};

void LPdata_Init ( LPdata *p, int n, int *J, double *x, double **A );
/*
 * DESCRIPTION :
 *   Initializes the data in p with the given parameters.
 *
 * ON ENTRY :
 *   p         memory allocated to hold one structure of type LPdata;
 *   n         dimension of the LP problem;
 *   J         index vector for the solution of the LP problem;
 *   x         vector of range 0..n-1, a solution vector;
 *   A         matrix of dimension n, basis inverse used in LP.
 *
 * ON RETURN :
 *   p         copies of J, x, and A are stored in p. */

void IndexNode_Init ( IndexNode *p, int i, LPdata *ptr );
/*
 * DESCRIPTION :
 *   Initializes the index node p with i and the LP data.
 *
 * ON ENTRY :
 *   p         memory allocated to hold one IndexNode structure;
 *   i         value for p->idx;
 *   ptr       data of the LP problem.
 *
 * ON RETURN :
 *   p         the pointer fields in p are assigned to i and ptr,
 *             no deep copies of the fields in ptr are made. */

void LPPL_Init ( LPPL *p, LPdata *A, LPPL *N );
/*
 * DESCRIPTION :
 *   Initializes the linked list of LP problems with the LP data in A
 *   and pointer to the next item in the list N.
 *
 * ON ENTRY :
 *   p         memory allocated for one structure;
 *   A         data for one LP problem;
 *   N         pointer information to the next item in the list.
 *
 * ON RETURN :
 *   p         structure with its two fields assigned to using A and N,
 *             no deep copies of the fields in p are made. */

void IT_LP_Init ( IT_LP *p, int nSpt, int *type );
/*
 * DESCRIPTION :
 *   Initialization of the index tree of linear programs p.
 *
 * ON ENTRY :
 *   p         memory allocated for one structure;
 *   nSpt      number of different supports;
 *   type      type of mixture of the supports,
 *             given by an array of range 0..nSpt-1.
 *
 * ON RETURN :
 *   p         allocated memory for all its members. */

int IT_IsEmpty ( IT_LP *p );
/*
 * DESCRIPTION :
 *   Returns the outcome of the test p->CurLevel < 1. */

int IT_IsFull ( IT_LP *p );  
/*
 * DESCRIPTION :
 *   Returns the outcome of the test p->CurLevel+1 >= p->MaxLevels. */

int IT_Level ( IT_LP *p );
/*
 * DESCRIPTION :
 *   Returns the value of p->CurLevel,
 *   i.e.: the value of the current level in the index tree. */

int IT_CurLPdim ( IT_LP *p );
/*
 * DESCRIPTION :
 *   Returns the value of p->DIM[p->CurLevel],
 *   i.e.: the dimension of the LP problems at the current level. */

int *IT_Cell ( IT_LP *p );
/*
 * DESCRIPTION :
 *   Returns the address of p->cell + 1,
 *   i.e.: the pointer containing the indices of a horizontal branch. */

int IT_CurSptIdx ( IT_LP *p );
/*
 * DESCRIPTION :
 *   Returns the value of p->InSpt[p->CurLevel],
 *   i.e.: the index of the support concerning the current node. */

int IT_MinNumPt ( IT_LP *p );
/*
 * DESCRIPTION :
 *   Returns the value of p->minNP[p->CurLevel],
 *   i.e.: the minimal #points required for moving to the next level. */

int IT_NumPt ( IT_LP *p );
/*
 * DESCRIPTION :
 *   Returns the value of p->NP[p->CurLevel],
 *   i.e.: the #points in the current vertical branch of the index tree. */

IndexNode *IT_FixedIdxNdPtr ( IT_LP *p );
/*
 * DESCRIPTION :
 *   For p nonempty, returns p->IT[p->CurLevel],
 *   i.e.: the current head of the vertical branch in the tree. */

void IT_StepBack ( IT_LP *p );
/*
 * DESCRIPTION :
 *   Sets the number of points at the current level to zero
 *   and decreases the value of p->CurLevel by one. */

int IT_Find ( IT_LP *p, int IDX );
/*
 * DESCRIPTION :
 *   Returns 0 if there is no index node in p with index equal to IDX,
 *   otherwise, 1 is returned and the p->curr->idx == IDX.
 *
 * ON ENTRY :
 *   p         index tree with LP problems;
 *   IDX       label to an index node.
 *
 * ON RETURN :
 *   0 if there is no index node with index equal to IDX,
 *   1 otherwise, in which case: p->curr->idx == IDX. */

IndexNode *IT_ResetCurLevelTo1 ( IT_LP *p );
/*
 * DESCRIPTION :
 *   Resets the value of p->CurLevel to 1, sets p->NP[1] to 1,
 *   sets p->prev to p->IT[1] and deallocates other entries;
 *   returns p->IT[1]. */

void IT_RenewNP1 ( IT_LP *p );
/*
 * DESCRIPTION :
 *   Increases the count of p->NP[1] by the number of elements in p->IT[1].
 *   The values of p->last[1] and p->cell[1] are changed as well. */

int IT_NextLevel ( IT_LP *p );
/*
 * DESCRIPTION :
 *   Returns 0 if p->CurLevel is equal to the p->MaxLevels, or
 *             if the number of nodes at the current level is fewer
 *                than the minimal number of nodes at that level,
 *   otherwise 1 is returned and pointers are swapped to move through
 *   to the next level in the index tree of LP problems. */

void IT_Add1
 ( IT_LP *p, int n, int *J, int nn, int *JJ, double *X, double **A );
/*
 * DESCRIPTION :
 *   Updates the index tree of LP problems with data for one LP problem,
 *   for an array of indices to points.
 *
 * ON ENTRY :
 *   p         index tree of LP problems;
 *   n         dimension of the array J;
 *   J         labels of points;
 *   nn        dimension of the LP problem;
 *   JJ        array of dimension nn with indices of the constraints involved;
 *   X         solution vector of dimension nn;
 *   A         basis inverse, matrix of dimension nn.
 *
 * ON RETURN :
 *   p         updated index tree of LP problems. */

void IT_Add2
 ( IT_LP *p, int oneidx, int nn, int *JJ, double *X, double **A );
/*
 * DESCRIPTION :
 *   Updates the index tree of LP problems with data for one LP problem,
 *   for just one point.
 *
 * ON ENTRY :
 *   p         index tree of LP problems;
 *   oneidx    index to the point in the support;
 *   nn        dimension of the LP problem;
 *   JJ        array of dimension nn with indices of the constraints involved;
 *   X         solution vector of dimension nn;
 *   A         basis inverse, matrix of dimension nn.
 *
 * ON RETURN :
 *   p         updated index tree of LP problems. */

void IT_LP_DEL ( IT_LP *p );
/*
 * DESCRIPTION :
 *   This is the main function to deallocate all memory occupied by p,
 *   it calls the other two IT_Free routines below. */

void IT_FreeIT ( IT_LP *p );
/*
 * DESCRIPTION :
 *   Deallocation of the index tree structures in p->IT. */

void IT_FreeLP ( IT_LP *p );
/*
 * DESCRIPTION :
 *   Deallocation of all the LP data in p->LP. */

#endif
