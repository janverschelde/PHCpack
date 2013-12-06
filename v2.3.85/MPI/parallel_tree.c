#include <stdio.h>
#include <stdlib.h>
#include "parallel_tree.h"
#include <assert.h>

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

Node* Top_Root ( int m, int p, int q )
{
  Node *nd = (Node*) malloc(sizeof(Node));
  int i, r, d;

  d = q/p;
  r = q-d*p;

  Allocate(nd,p);
  nd->tp = top;
  nd->level = 0;

  for(i=0; i<p; i++)
  {
    if(i<p-r)
    {
      nd->top[i] = d*(m+p)+m+r+i+1;
      nd->bottom[i] = d*(m+p)+m+r+i+1;
    } 
    else
    { 
      nd->top[i] = (d+1)*(m+p)+m+(i-p+r+1); 
      nd->bottom[i] = (d+1)*(m+p)+m+(i-p+r+1);
    } 
  }

  if(v_f==1)
  {
    printf("The root of the tree is:");
    Print_Localization(p, nd->top);
    printf("The level of the root: %d\n", nd->level);
    printf("Build the tree:\n");
  }
  return nd;
}

Node* Bottom_Root ( int p )
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

int Q_Roots_Valid ( int *pivot, int m, int p, int sum)
/* Return true if it is a valid root for the forest of the mixed pivots*/
{
  int i,pivot_sum=0;

  if((pivot[p-1]-pivot[0])>=(m+p))
    return 0;

  for(i=0;i<p-1;i++)
    if(pivot[i]>=pivot[i+1])
      return 0;

  for(i=0;i<p;i++)
    pivot_sum+=pivot[i];

  if(pivot_sum==sum)
    return 1;
  else
    return 0;
}

void Mixed_Root ( int m, int p, int q, Node *mixed_leave, 
                  queue_node **nq_front, queue_node **nq_end )
{
  int i,j,k,pivot[2*p],sum=0,
      *low_bound=mixed_leave->top,
      *up_bound=mixed_leave->bottom; 
  int counter=0,
      p_first=mixed_leave->bottom[0],
      p_last=mixed_leave->bottom[p-1];
  Node *nd;

  /* find the sum of the pivots of the roots */
  for(i=0; i<p; i++)
    sum=sum+up_bound[i]-low_bound[i];
  sum=sum/2+(1+p)*p/2;
  if(v_f==1)
     printf("The sum of the pivots of the roots: %d\n", sum);

  /* find the up bound for the first pivot */
  for(i=1; i<p_first; i++)
    if((i+p+i-1)*p/2>sum)
      break;
  up_bound[0]=i-1;
  if(v_f==1)
    printf("The up bound for the first pivot : %d\n", up_bound[0]);

  /* fine the low bound for the last pivot */
  for(i=p_last; i>p; i--)
    if((i+i-p+1)*p/2<sum)
      break; 
  low_bound[p-1]=i+1;
  if(v_f==1)
     printf("The low bound for the last pivot : %d\n", low_bound[p-1]);

  for(i=0; i<p; i++)
    pivot[i]=mixed_leave->top[i];
  
  if(v_f==1)
    printf("The roots for the mixed pivot methods : \n");
  i = p-2;
  while(i>=0)
  {
    while(pivot[p-1]<=p_last)
    {
      if(Q_Roots_Valid(pivot,m,p,sum))
      {
	 if(v_f==1)
           Print_Localization(p,pivot);
         counter++;
   
         nd = (Node*) malloc(sizeof(Node));
         Allocate(nd,p);
         nd->tp=mixed;
         nd->level=0;
         for(j=0; j<p; j++)
         {
           nd->top[j] = pivot[j];
           nd->bottom[j] = pivot[j];
         }
        append_node(nq_front,nq_end,nd);       
      }
      pivot[p-1]++;
    }

    for(k=p-2; k>i; k--)
    {
      if(pivot[k]<up_bound[k])
        break;
    }
       
    if(k>=i && pivot[k]<up_bound[k])    
    {
      pivot[k]++;
      for(j=p-1; j>k; j--)
      {
        pivot[j]=low_bound[j];
      }  
      
    }
    else
    {
      pivot[--i]++;
      for(j=p-1; j>i; j--)
      {
        pivot[j]=low_bound[j];
      }  
    }
  }
  if(v_f==1)
    printf("Total combination number: %d\n", counter);
}

void Copy_Pivot ( int *a, int *b, int *c, int p )
/* Appends the values in array b to array a and saves them in array c */
{
  int i;

  /* Copies the top pivot to array c */
  for(i=0; i<p; i++)
    c[i] = a[i];

  /* Copies the bottom pivot to array c */
  for(i=p; i<2*p; i++)
    c[i] = b[i-p];
}

void Copy_Pivot1 ( int *a, int *b, int p )
/* Copies the values in array a to array b */
{
  int i;

  for(i=0; i<p; i++)
    b[i] = a[i];
}

int Top_Creatable ( Node *nd, int i, int m, int p )
/* Return true if the i-th top pivot can be decremented */
{
  if(nd->top[i]<(i+1))
    return 0;
  else if(i==0)
    return (nd->top[i]>=2);
  else
    return (nd->top[i]-1>nd->top[i-1]);
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

Node* Top_Leave ( int m, int p, int q )
{
  Node *nd = (Node*) malloc(sizeof(Node));
  int i, r, d;

  d = q/p;
  r = q-d*p;

  Allocate(nd,p);
  nd->tp = top;
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
    Print_Localization(p, nd->top);
    printf("The level of the leaf: %d\n", nd->level);
  }
  return nd;
}

Node* Bottom_Leave ( int m, int p, int q )
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

Node* Mixed_Leave ( int m, int p, int q )
{
  Node *nd = (Node*) malloc(sizeof(Node));
  int i, r, d;

  d = q/p;
  r = q-d*p;

  Allocate(nd,p);
  nd->tp = mixed;
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
    printf("The leaf of the tree is:\n");
    Print_Localization(p, nd->top);
    Print_Localization(p, nd->bottom);
    printf("The level of the leaf: %d\n", nd->level);
  }
  return nd;
}

int Q_Top_Creatable ( Node *nd, int i, int m, int p, int q, Node *leave)
/* Return true if the i-th top pivot can be decremented */
{
    if(nd->top[i]<=leave->top[i])
      return 0;
    else if(i==0 && nd->top[p-1]+1-nd->top[0]>=m+p)
      return 0;
    else if(i>0)
      return(nd->top[i]-1>nd->top[i-1]);
    else
      return 1;
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
    else 
      return 1;
}

int Q_Mixed_Creatable ( Node *nd, int j, int k, int m, int p, int q,
                        Node *top_leave, Node *bottom_leave )
/* Return true if j-th top pivot can be decremented and k-th bottom
   pivot can be incremented */
{
  if(Q_Top_Creatable(nd,j,m,p,q,top_leave) && 
     Q_Bottom_Creatable(nd,k,m,p,q,bottom_leave))
    return 1;
  else
    return 0;  
}


Node* Create_Parent ( Node *nd, int i, int p )
/* The link on return is a newly created link to node with contents parent. */
{
  Node *parent = (Node*) malloc(sizeof(Node));

  Allocate(parent, p); 
  if(nd->tp==top)
  {
    parent->tp = top; 
    parent->level = nd->level+1;
    Copy_Pivot1(nd->top,parent->top,p);
    parent->top[i] = nd->top[i]-1;
    Copy_Pivot1(nd->bottom,parent->bottom,p);
    /* keep the bottom pivots unchanged */    
  }
  else if(nd->tp==bottom)
  {
    parent->tp = bottom; 
    parent->level = nd->level+1;
    Copy_Pivot1(nd->bottom,parent->bottom,p);
    parent->bottom[i] = nd->bottom[i]+1;
    Copy_Pivot1(nd->top,parent->top,p);
    /* keep the top pivots unchanged */
  }
  return parent;
}

Node* Create_Parent_Mixed ( Node *nd, int i, int j, int p, int diff )
/* The link on return is a newly created link to node with contents parent
   for mixed pivot method */
{
  Node *parent = (Node*) malloc(sizeof(Node));
 
  Allocate(parent, p);
  parent->tp = mixed;
  parent->level = nd->level+diff;
  Copy_Pivot1(nd->top,parent->top,p);
  Copy_Pivot1(nd->bottom,parent->bottom,p);
  if(diff!=1)   
    parent->top[i] = nd->top[i]-1;
  /* else the difference of levels between parent and child is one, leave 
     the top pivot unchanged */

  parent->bottom[j] = nd->bottom[j]+1;
  return parent;
}

void Copy_Sol ( double *a, double *b, int len )
/* Copy the values in array a to array b */
{
  int i;

  for(i=0; i<len; i++)
    b[i] = a[i];
}

void Create_Top_Parent ( Node *nd, int m, int p, double *sol, 
                           queue **q, queue **end )
{
  int i;
  Node *tmp;
  job *j;

  for(i=0; i<p; i++)
  {
    if(Top_Creatable(nd,i, m, p))
    {
      nd->parents[0][i] = Create_Parent(nd,i,p);
      if(v_f==1)
        Print_Localization(p,nd->parents[0][i]->top);

      tmp = nd->parents[0][i];
      j = (job*) calloc(1, sizeof(job));
      Allocate_job(j,p,nd->level*2);

      Copy_Pivot(nd->top,nd->bottom,j->start_pivot,p);
      Copy_Pivot(tmp->top,tmp->bottom,j->target_pivot,p);
      j->length = nd->level*2; 
      j->address = (int) tmp; /* Passes the address of the node */  
      Copy_Sol(sol,j->start_sol,j->length);
  
      append(q,end,j);
    }
    if(i==p-1)       /*  All parents for node nd have been found,delete nd */ 
      Deallocate(nd,p);     
  }
}
			
void Create_Bottom_Parent ( Node *nd, int m, int p, double *sol, 
                           queue **q, queue **end)
{
  int i;
  Node *tmp;
  job *j;

  for(i=0; i<p; i++)
  {
    if(Bottom_Creatable(nd,i, m, p))
    {
      nd->parents[0][i] = Create_Parent(nd,i,p);
      if(v_f==1)
        Print_Localization(p,nd->parents[0][i]->bottom);

      tmp = nd->parents[0][i];
      j = (job*) calloc(1, sizeof(job));
      Allocate_job(j,p,nd->level*2);

      Copy_Pivot(nd->top,nd->bottom,j->start_pivot,p);
      Copy_Pivot(tmp->top,tmp->bottom,j->target_pivot,p);
      j->length = nd->level*2; 
      j->address = (int) tmp; /* Passes the address of the node */  
      Copy_Sol(sol,j->start_sol,j->length);
  
      append(q,end,j);
    }
    if(i==p-1)     /*  All parents for node nd have been found,delete nd */ 
      Deallocate(nd,p);    
  }
}

void Q_Create_Top_Parent ( Node *nd, int m, int p, int q, Node *leave,
                           double *sol, queue **front, queue **end )
/* Creates a parent of the node by decrementing the top pivot recursively */
{
  int i;
  Node *tmp;
  job *j;

  if(q==0)
  {
    Create_Top_Parent(nd,m,p,sol,front,end); 
  }
  else
  {
    for(i=0; i<p; i++)
    {
      if(Q_Top_Creatable(nd,i,m,p,q,leave))
      {
        nd->parents[0][i] = Create_Parent(nd,i,p);
        if(v_f==1)
          Print_Localization(p,nd->parents[0][i]->top);
        tmp = nd->parents[0][i];
        j = (job*) calloc(1, sizeof(job));
        Allocate_job(j,p,nd->level*2);

        Copy_Pivot(nd->top,nd->bottom,j->start_pivot,p);
        Copy_Pivot(tmp->top,tmp->bottom,j->target_pivot,p);
        j->length = nd->level*2;
        j->address = (int) tmp; /* Passes the address of the node */
        Copy_Sol(sol,j->start_sol,j->length);

        append(front,end,j);
      }
      if(i==p-1)       /*  All parents for node nd have been found,delete nd */
        Deallocate(nd,p);
    }
  }
}

void Q_Create_Bottom_Parent ( Node *nd, int m, int p, int q, Node *leave,
                            double *sol, queue **front, queue **end)
/* Creates a parent of the node by incrementing the bottom pivot recursively */
{
  int i;
  Node *tmp;
  job *j;

  if(q==0)
  {
    Create_Bottom_Parent(nd,m,p,sol,front,end); 
  }
  else
  {
    for(i=0; i<p; i++)
    {
      if(Q_Bottom_Creatable(nd,i,m,p,q,leave))
      { 
        nd->parents[0][i] = Create_Parent(nd,i,p);
        if(v_f==1)
          Print_Localization(p,nd->parents[0][i]->bottom);
        tmp = nd->parents[0][i];
        j = (job*) calloc(1, sizeof(job));
        Allocate_job(j,p,nd->level*2);

        Copy_Pivot(nd->top,nd->bottom,j->start_pivot,p);
        Copy_Pivot(tmp->top,tmp->bottom,j->target_pivot,p);
        j->length = nd->level*2; 
        j->address = (int) tmp; /* Passes the address of the node */
        Copy_Sol(sol,j->start_sol,j->length);
    
        append(front,end,j);
      }
    }
    if(i==p-1)       /*  All parents for node nd have been found,delete nd */
      Deallocate(nd,p);
  }
}


void Q_Create_Mixed_Parent ( Node *nd, int m, int p, int q, Node *top_leave,
        Node *bottom_leave, int diff, double *sol, queue **front, queue **end)
/* Creates a parent of the node by decrementing the top pivot and incrementing
   the bottom pivot recursively */
{
  int j,k;
  Node *parent;
  job *J;
  
  if(diff!=1)
  {
    for(j=0; j<p; j++)
      for(k=0; k<p; k++)
      {              
        if(Q_Mixed_Creatable(nd,j,k,m,p,q,top_leave,bottom_leave))
        {    
          nd->parents[j][k] = Create_Parent_Mixed(nd,j,k,p,diff);
          parent = nd->parents[j][k];
          if(v_f==1)
          {
            printf("parent[%d][%d]:\n",j,k);
            Print_Localization(p,parent->top);
            Print_Localization(p,parent->bottom);
          }
                
          J = (job*) calloc(1, sizeof(job));
          Allocate_job(J,p,nd->level*2);
          Copy_Pivot(nd->top,nd->bottom,J->start_pivot,p);
          Copy_Pivot(parent->top,parent->bottom,J->target_pivot,p);
          J->length = nd->level*2;
          J->address = (int) parent;
          Copy_Sol(sol,J->start_sol,J->length);
          append(front,end,J);
        }
        if((j==p-1)&&(k==p-1)) /*All parents for node nd have been found*/
          Deallocate(nd,p);
      }     
  }
  else /* updates bottom pivot only for diff=1 */
  {
    j=0;
    for(k=0; k<p; k++)
    {              
      if(Q_Bottom_Creatable(nd,k,m,p,q,bottom_leave))
      {    
        nd->parents[j][k] = Create_Parent_Mixed(nd,j,k,p,diff);
        parent = nd->parents[j][k];
        if(v_f==1)
        {
          printf("parent[%d][%d]:\n",j,k);
          Print_Localization(p,parent->top);
          Print_Localization(p,parent->bottom);
        }
                
        J = (job*) calloc(1, sizeof(job));
        Allocate_job(J,p,nd->level*2);
        Copy_Pivot(nd->top,nd->bottom,J->start_pivot,p);
        Copy_Pivot(parent->top,parent->bottom,J->target_pivot,p);
        J->length = nd->level*2;
        J->address = (int) parent;
        Copy_Sol(sol,J->start_sol,J->length);
        append(front,end,J);
      }
      if((j==p-1)&&(k==p-1)) /*All parents for node nd have been found*/
        Deallocate(nd,p);
    }     
  }
}

void append_node( queue_node **front, queue_node **end, Node *nd )
{
   queue_node* nq = (queue_node*) calloc(1, sizeof(queue_node));

   nq->item = nd;

   if(*front==NULL)
   {
      nq->next = NULL;
      *end = nq;
      *front = nq;
   }
   else
   {
     (*end)->next = nq;
     nq->next = NULL;
     *end = nq;
   }
}

queue_node* pop_node( queue_node*q, Node **nd_ptr )
{
   queue_node *nq;

   assert(q!=NULL);
   *nd_ptr = q->item;
   nq = q->next;
   free(q);

   return nq;
}

