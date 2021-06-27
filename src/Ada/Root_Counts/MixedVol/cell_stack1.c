/* The file "cell_stack1.c" contains the definitions of the prototypes
 * in the file "cell_stack1.h". */

#include <stdio.h>
#include <stdlib.h>
#include "cell_stack1.h"

#define verbose_stack 0

void Cell_Init1 ( cell *c )
{
   c->next=0;
}
         
void Cell_Init2 ( cell *c, int n, int *J, cell *ptr )
{
   int i;

   c->idx = (int*)calloc(n,sizeof(int));
   for(i=0; i< n; i++) c->idx[i] = J[i];
   c->next=ptr; 
}
   
void Cs_Init ( CellStack *cs, int n )
{  
   cs->size=n;
   cs->count=0;
   cs->cur=cs->top=0; 
}
     
void Cs_Push ( CellStack *cs, int *J )
{
   cell *c;

   if(verbose_stack > 0)
   {
      int i;

      printf("Pushing indices of size %d on the stack :",cs->size);
      for(i=0; i< cs->size; i++) printf(" %d",J[i]); printf("\n");
   }

   c = (cell*)calloc(1,sizeof(cell));
         
   Cell_Init2(c,cs->size,J,cs->top);

   cs->cur = cs->top=c;

   ++(cs->count);
}

void Cs_Pop ( CellStack *cs )
{
   cell *ptr = cs->top;

   cs->cur = cs->top = (cs->top)->next;
   free(ptr->idx);
   free(ptr);
   --(cs->count);
}
      
int Cs_Next ( CellStack *cs )
{
   if(cs->cur->next)
   {
      cs->cur = cs->cur->next;
      return 1;
   }
   else
      return 0;
}

int *Cs_Cur ( CellStack *cs )
{
   return cs->cur->idx;
}

void Cs_Top ( CellStack *cs )
{ 
   cs->cur = cs->top;
}
      
int Cs_IsEmpty ( CellStack *cs )
{ 
   if(cs->top==0)
      return 1;
   else
      return 0; 
}

int Cs_Count ( CellStack *cs )
{
   return cs->count;
}
      
int *Csi ( CellStack *cs, int i )
{
   if(-1< i && i< cs->count)
   {
      int c=0;

      cs->cur = cs->top;
      while(c<i)
      {
         cs->cur=cs->cur->next;
         c++;
      }
      return cs->cur->idx;
   }
   else
      return 0;
}

void Cs_Del ( CellStack *cs )
{
   while(!Cs_IsEmpty(cs)) Cs_Pop(cs);
   cs->count = 0;
}
