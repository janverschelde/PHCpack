
#ifndef PARALLEL_TREE_H
#define PARALLEL_TREE_H

#include "queue.h"

typedef struct Node Node;
typedef enum {top,bottom,mixed} Node_Type; 
 
struct Node
{
  Node_Type tp;           /* type of node */
  int level;
/*The level of a node is a natural number between 0 and m*p+q*(m+p),*/  
  int *top, *bottom;
  Node ***parents;        /* a matrix of pointers to node */
//Node *next;  /* maintains a list for the roots of the mixed pivot method */  
};

typedef struct Queue_Node queue_node;

struct Queue_Node
{
   Node   *item;      /* a pointer to the node */
   queue_node *next;  /* points to next node   */
};


Node* Top_Root ( int m, int p, int q );
/* Returns a pointer for the root of the localization tree with top pivots, 
   with top and bottom pivots for the trivial localization pattern for any 
   p-plane in (m+p)-space.                                             */

Node* Bottom_Root ( int p );
/* Returns a pointer for the root of the localization tree with bottom pivots, 
   with top and bottom pivots for the trivial localization pattern for any 
   p-plane in (m+p)-space.                                             */

void Mixed_Root ( int m, int p, int q, Node *leave, queue_node **front, 
                  queue_node **end );
/* Finds all of the roots for the mixed pivots method and put them in the queue */

Node* Top_Leave ( int m, int p, int q );
/* Returns a pointer for the leaf of the localization tree with top pivots */        

Node* Bottom_Leave ( int m, int p, int q );
/* Returns a pointer for the leaf of the localization tree with bottom pivots */

Node* Mixed_Leave ( int m, int p, int q );
/* Returns a pointer for the leaf of the localization tree with mixed pivots */

void Create_Top_Parent ( Node *nd, int m, int p, double *sol, 
			    queue **q, queue **end );
/* Creates parents of the node by decrementing the top pivot            */


void Create_Bottom_Parent ( Node *nd, int m, int p, double *sol, 
			    queue **q, queue **end );
/* Creates parents of the node by incrementing the bottom pivot            */

void Q_Create_Top_Parent ( Node *nd, int m, int p, int q, Node *leave,
			    double *sol, queue **front, queue **end );
/* Creates parents of the node by decrementing the top pivot            */

void Q_Create_Bottom_Parent ( Node *nd, int m, int p, int q, Node *leave,
			    double *sol, queue **front, queue **end );
/* Creates parents of the node by incrementing the bottom pivot            */

void Q_Create_Mixed_Parent ( Node *nd, int m, int p, int q, Node *top_leave,
       Node *bottom_leave, int diff, double *sol, queue **front, queue **end);
/* Creates a parent of the node by decrementing the top pivot and incrementing
   the bottom pivot recursively */

void append_node( queue_node **q, queue_node **end, Node *nd );
/* Appends a node to the end of the queue for node */

queue_node* pop_node( queue_node*q, Node **nd_ptr );
/* Returns a job from the front of the queue for node */

#endif
