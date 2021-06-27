/* The file "cell_stack1.h" collects the definition of a stack of mixed cells
 * and the prototypes of the operations on the stack. */

#ifndef _CELL_STACK1_H_
#define _CELL_STACK1_H_

typedef struct cell cell;
typedef struct CellStack CellStack;

struct cell
{
   int *idx;     /* labels of the points which span the mixed cell */
   cell *next;   /* pointer to the next cell */
};

struct CellStack
{
   int size;     /* total #points which span a mixed cell */
   int count;    /* number of elements in the stack */
   cell *top;    /* top of the stack */
   cell *cur;    /* pointer to the current cell in the stack */
};   

void Cell_Init1 ( cell *c );
/*
 * DESCRIPTION :
 *   Sets c->next to zero. */
         
void Cell_Init2 ( cell *c, int n, int *J, cell *ptr );
/*
 * DESCRIPTION :
 *   Defines the content of the cell with the indices in J.
 *
 * ON ENTRY :
 *   c        memory allocated for one cell;
 *   n        number of indices of J;
 *   J        labels to the points spanning the mixed cell;
 *   ptr      pointer for the the next cell.
 *
 * ON RETURN :
 *   c        c->idx contains the values of J and
 *            c->next has the same value as ptr. */

void Cs_Init ( CellStack *cs, int n );
/*
 * DESCRIPTION :
 *   Initializes the stack cs to contain mixed cells spanned by n points. */
     
void Cs_Push ( CellStack *cs, int *J );
/*
 * DESCRIPTION :
 *   Pushes a mixed cell defined by the indices in J to the stack. */

void Cs_Pop ( CellStack *cs );
/*
 * DESCRIPTION :
 *   Pops the top element from the cell stack. */
      
int Cs_Next ( CellStack *cs );
/*
 * DESCRIPTION :
 *   Returns 0 if the next cell to the current cell is empty,
 *   otherwise returns 1 and puts the pointer to the current cell
 *   to the next cell. */

int *Cs_Cur ( CellStack *cs );
/*
 * DESCRIPTION :
 *   Returns the content of the current cell. */

void Cs_Top ( CellStack *cs );
/*
 * DESCRIPTION :
 *   Assigns the current cell pointer to the top of the stack. */
      
int Cs_IsEmpty ( CellStack *cs );
/*
 * DESCRIPTION :
 *   Returns 1 if the cell stack is empty, 0 otherwise. */

int Cs_Count ( CellStack *cs );
/*
 * DESCRIPTION :
 *   Returns the value of cs->count, the #cells in the stack. */
      
int *Csi ( CellStack *cs, int i );
/*
 * DESCRIPTION :
 *   Returns the content of the i-th cell in the stack. */

void Cs_Del ( CellStack *cs );
/*
 * DESCRIPTION :
 *   Pops all the cells of the stack, deallocating all memory. */

#endif
