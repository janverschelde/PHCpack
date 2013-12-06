#ifndef PIERI_TREE_H
#define PIERI_TREE_H

typedef struct Node Node;
typedef enum {top,bottom,mixed} Node_Type; 
 
struct Node
{
  Node_Type tp;           /* type of node */
  int level;
/*The level of a node is a natural number between 0 and m*p+q*(m+p),*/  
  int *top, *bottom;
  Node ***parents;        /* a matrix of pointers to node */ 
  
};

int Build_Bottom_Tree(int m, int p);
/* Creates a tree by incrementing the bottom pivots recursively and returns */
/* the total number of roots for q=0                                        */

int Q_Build_Bottom_Tree(int m, int p, int q);
/* Creates a tree by incrementing the bottom pivots recursively and returns */
/* the total number of roots for q>0                                        */

#endif
