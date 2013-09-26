#include <stdio.h>
#include <stdlib.h>
#include "pieri_tree.h"

#define v_f 0  /* verbose flag */

void Allocate ( Node *nd, int p )
/* Allocate memory for a new node */
{
  int i;

  nd->top = (int*) calloc(p, sizeof(int));
  nd->bottom = (int*) calloc(p, sizeof(int));
  nd->parents = (Node***) calloc(p, sizeof(Node**));
  for(i=0; i<p; i++)
    nd->parents[i] = (Node**) calloc(p, sizeof(Node*));
}

void Deallocate ( Node *nd, int p )
/* Deallocate memory for a node */
{
  int i;
  
  free(nd->top);
  free(nd->bottom);
  for(i=0; i<p; i++)
    free(nd->parents[i]);
  free(nd->parents);
  free(nd);    
}

void Print_Localization(int p, int *a)
{
  int i;

  printf("[");
  for(i=0; i<p; i++)
  {
    if(i!=p-1)
      printf("%d ", a[i]);
    else
      printf("%d]\n", a[i]);
  }  
}

Node* Root ( int p )
/* Returns a pointer for the root of the localization tree, with top and     */
/* bottom pivots for the trivial localization pattern for any p-plane in     */
/* (m+p)-space.                                                              */
{
  Node *nd = (Node*) malloc(sizeof(Node));
  int i;  

  Allocate(nd,p);
  nd->tp = bottom;
  nd->level = 0;

  for(i=0; i<p; i++)
  {
    nd->top[i] = i+1;
    nd->bottom[i] = i+1;
  }

  if(v_f==1)
    printf("Build the tree:\n");
  return nd;
}

void Copy_Pivot ( int *a, int *b, int p )
/* Copy the values in array a to array b */
{
  int i;

  for(i=0; i<p; i++)
    b[i] = a[i];
}

int Bottom_Creatable ( Node *nd, int i, int m, int p )
/* Return true if the i-th bottom pivot can be incremented */
{
  if(nd->bottom[i]>(m+i+1))
    return 0;
  else if(i==(p-1))
    return (nd->bottom[i]<=(m+p-1));
  else
    return (nd->bottom[i]+1<nd->bottom[i+1]);
}					

Node* Leave ( int m, int p, int q )
/* Find the bottom pivot of the leaf */
{
  Node *nd = (Node*) malloc(sizeof(Node));
  int i, r, d;

  d = q/p;
  r = q-d*p;

  Allocate(nd,p);
  nd->tp = bottom;
  nd->level = m*p+q*(m+p);

  for(i=0; i<p; i++)
  {
    nd->top[i] = i+1;
    if(i<p-r) 
      nd->bottom[i] = d*(m+p)+m+r+i+1;
    else
      nd->bottom[i] = (d+1)*(m+p)+m+(i-p+r+1);
  }

  if(v_f==1)
  {
    printf("The leaf of the tree is:");
    Print_Localization(p, nd->bottom);
    printf("The level of the leaf: %d\n", nd->level);
  }
  return nd;
}

int Q_Bottom_Creatable ( Node *nd, int i, int m, int p, int q, Node *leave)
/* Return true if the i-th bottom pivot can be incremented */
{
    if(nd->bottom[i]>=leave->bottom[i])
      return 0;
    else if(i==(p-1) && nd->bottom[p-1]+1-nd->bottom[0]>=m+p)
      return 0;
    else if(i<p-1) 
      return(nd->bottom[i]+1<nd->bottom[i+1]);
}

Node* Create_Parent ( Node *nd, int i, int p )
/* The link on return is a newly created link to node with contents parent. */
{
  Node *parent = (Node*) malloc(sizeof(Node));

  Allocate(parent, p);
  parent->tp = bottom; 
  parent->level = nd->level+1;
  Copy_Pivot(nd->bottom,parent->bottom,p);
  parent->bottom[i] = nd->bottom[i]+1;
  Copy_Pivot(nd->top,parent->top,p);

  return parent;
}
			
int Create_Bottom_Parent ( Node *nd, int m, int p )
/* Creates a parent of the node by incrementing the bottom pivot recursively */
{
  static sum = 0;
  int i;
  Node *tmp;

  for(i=0; i<p; i++)
  {
    if(Bottom_Creatable(nd,i, m, p))
    {
      nd->parents[0][i] = Create_Parent(nd,i,p);
      if(v_f==1)
        Print_Localization(p,nd->parents[0][i]->bottom);
      tmp = nd->parents[0][i];
      if(i==p-1)       /*  All parents for node nd have been found,delete nd */
      {
        Deallocate(nd,p);
      }
      Create_Bottom_Parent(tmp,m,p);
    }
    else
    {
      if(nd->level==m*p)
      {
        sum++;
        if(v_f==1) printf("\n");
      }

      if(i==p-1)       /*  All parents for node nd have been found,delete nd */
      {
        Deallocate(nd,p);
      }
    }
  }
  return sum/p;
}

int Build_Bottom_Tree(int m, int p)
{
  int sum;
  Node* root;

  root = Root(p);
  sum = Create_Bottom_Parent(root, m, p);

  return sum;
}


int Q_Create_Bottom_Parent ( Node *nd, int m, int p, int q, Node *leave )
/* Creates a parent of the node by incrementing the bottom pivot recursively */
{
  static sum = 0;
  int i;
  Node *tmp;

  for(i=0; i<p; i++)
  {
    if(Q_Bottom_Creatable(nd,i,m,p,q,leave))
    { 
      nd->parents[0][i] = Create_Parent(nd,i,p);
      if(v_f==1)
        Print_Localization(p,nd->parents[0][i]->bottom);
      tmp = nd->parents[0][i];
      if(i==p-1)       /*  All parents for node nd have been found,delete nd */
      {
        Deallocate(nd,p);
      }
      Q_Create_Bottom_Parent(tmp,m,p,q,leave);
    }
    else
    {
      if(nd->level==leave->level)
      { 
        sum++;
        if(v_f==1) printf("\n");
      }

      if(i==p-1)       /*  All parents for node nd have been found,delete nd */
      {
        Deallocate(nd,p);
      }
    }
  }
  return sum/p;
}

int Q_Build_Bottom_Tree(int m, int p, int q)
{
  int sum;
  Node *root, *leave;

  if(q==0)
    return Build_Bottom_Tree(m,p); 
  else
  {
    leave = Leave(m,p,q);
    root = Root(p);
    sum = Q_Create_Bottom_Parent(root,m,p,q,leave);
    Deallocate(leave,p);
    return sum;
  }
}
