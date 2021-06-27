/* This file "zero_index_tree1.c" defines the functions whose prototypes
 * are documented in the file "zero_index_tree1.h". */

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include "zero_index_tree1.h"

void L0IdxNode_Init ( L0IdxNode *p, int i )
{
   p->idx = i;
   p->R = 0;
   p->D = 0; 
}

void L0_IML_Init ( L0_IML *p ) 
{
   p->L0head = (L0IdxNode*)calloc(1,sizeof(L0IdxNode));
   L0IdxNode_Init(p->L0head,-1);
   p->L0prev = p->L0curr = p->L0head;
   p->LP1 = (LPPL*)calloc(1,sizeof(LPPL));
   LPPL_Init(p->LP1,0,0);
}

int L0_Migrate ( L0_IML *p, IndexNode *inp )
{
   if(!inp)
   {
      printf("\aMigrate: inp = 0!\n");
      exit(1);
   }
   if(p->L0head->D)
   {
      p->L0prev = p->L0head->D;
      inp->idx  = p->L0prev->idx;
      inp->S    = p->L0prev->R;

      p->L0head->D = p->L0prev->D;
      free(p->L0prev);
      return 1;
   }
   else
   {
      free(p->L0head);
      return 0;
   }
}

int L0_FindInR ( L0_IML *p, int IDX )
{
   for(p->curr = p->prev->S; p->curr; p->curr = p->curr->S)
   {
      if(IDX <= p->curr->idx)
      {
         if(IDX == p->curr->idx)
            return 1;              /* curr->idx = IDX */
         else
            return 0;              /* prev->idx < IDX < curr->idx */
      }
      p->prev = p->prev->S;
   }
   return 0; /* => all indices are smaller than IDX, */
             /* i.e. *** < prev->idx < IDX ( prev->S = curr = 0 ) */
} 

int L0_FindInD ( L0_IML *p, int IDX )
{
   for(p->L0curr = p->L0prev->D; p->L0curr; p->L0curr = p->L0curr->D)
   {
      if(IDX <= p->L0curr->idx)
      {
         if(IDX == p->L0curr->idx)
            return 1;              /* L0curr->idx = IDX */
         else
            return 0;              /* L0prev->idx < IDX < L0curr->idx */
      }
      p->L0prev = p->L0prev->D;
   }
   return 0; /* => all indices are smaller than IDX, */
             /* i.e. *** < L0prev->idx < IDX ( L0prev->D = L0curr = 0 ) */
} 

void L0_Add1 
 ( L0_IML *p, int n, int *J, int N, int *I, double *X, double **A )
{
   int LPnotused = 1;
   LPdata *lpp;
   int i,j;

   lpp = (LPdata*)calloc(1,sizeof(LPdata));
   LPdata_Init(lpp,N,I,X,A);
   p->L0prev = p->L0head;
   
   for(i=0; i<n; i++)                      /* loop thru all pts in J*/
   {
      j = i+2;
      if(L0_FindInD(p,J[i]))                    /* then L0curr != 0 */
      {
         if(i+1<n)
         {
            if(p->L0curr->R)
            {
               p->L0prev = p->L0curr;
                 p->curr = p->L0curr->R;
               if(p->curr->idx > J[i+1])
               {  
                  p->curr = (IndexNode*)calloc(1,sizeof(IndexNode));
                  IndexNode_Init(p->curr,J[i+1],lpp);
                  p->curr->S = p->L0prev->R;
                  p->L0prev->R = p->curr;
                  LPnotused = 0;
               }
               else if(p->curr->idx < J[i+1]) j = i+1;
               for( ; j<n; j++)
               {
                  p->prev = p->curr;
                  if(!L0_FindInR(p,J[j]))
                  {
	             p->curr = (IndexNode*)calloc(1,sizeof(IndexNode));
                     IndexNode_Init(p->curr,J[j],lpp);
                     p->curr->S = p->prev->S;
                     p->prev->S = p->curr;
                     LPnotused = 0;
                  }
               }
            }
            else
            {  
               p->curr = (IndexNode*)calloc(1,sizeof(IndexNode));
               IndexNode_Init(p->curr,J[i+1],lpp); 
               p->L0curr->R = p->curr;
               LPnotused = 0;
               for( ; j<n; j++) 
               {
                  p->prev = p->curr;               /* Add2S(J[j],lpp); */
                  p->curr = (IndexNode*)calloc(1,sizeof(IndexNode));
                  IndexNode_Init(p->curr,J[j],lpp);
                  p->curr->S = p->prev->S;
                  p->prev->S = p->curr;
               }
            }
         }
      }
      else                                    /* add J[i] after L0prev */
      {
         p->L0curr = (L0IdxNode*)calloc(1,sizeof(L0IdxNode));
         L0IdxNode_Init(p->L0curr,J[i]);
         LPnotused = 0;
         p->L0curr->D = p->L0prev->D;
         p->L0prev->D = p->L0curr; 
         if(i+1<n)
         {
            p->curr = (IndexNode*)calloc(1,sizeof(IndexNode));
            IndexNode_Init(p->curr,J[i+1],lpp);
            p->L0curr->R = p->curr;
            for( ; j<n; j++)
            {
               p->prev = p->curr;                 /* Add2S(J[j],lpp); */
               p->curr = (IndexNode*)calloc(1,sizeof(IndexNode));
               IndexNode_Init(p->curr,J[j],lpp);
               p->curr->S = p->prev->S;
               p->prev->S = p->curr;
            }
         }
      }
      p->L0prev = p->L0curr;
   }
   if(LPnotused)
   {
      free(lpp->JJJ);
      free(lpp->xxx);
      for(i=0; i<lpp->dim; i++) free((lpp->INV)[i]);
      free(lpp->INV);
      free(lpp);
   }
   else
   {
      p->LP1->next = (LPPL*)calloc(1,sizeof(LPPL));
      LPPL_Init(p->LP1->next,lpp,p->LP1->next);
   }
}

void L0_Add2 ( L0_IML *p, int *J, int N, int *I, double *X, double **A )
{
   LPdata *lpp;
   int i;

   p->L0prev = p->L0head;
   for(i=0; i<2; i++)                         /* loop thru all pts in J */
   {
      if(L0_FindInD(p,J[i]))                        /* then L0curr != 0 */
      {
         if(i==0)
         {
            if(p->L0curr->R)
            {
               p->L0prev = p->L0curr;
                 p->curr = p->L0curr->R;
               if(p->curr->idx > J[1])
               {
                  lpp = (LPdata*)calloc(1,sizeof(LPdata));
                  LPdata_Init(lpp,N,I,X,A);
                  p->LP1->next = (LPPL*)calloc(1,sizeof(LPPL));
                  LPPL_Init(p->LP1->next,lpp,p->LP1->next);
                  p->curr = (IndexNode*)calloc(1,sizeof(IndexNode));
                  IndexNode_Init(p->curr,J[1],lpp);
                  p->curr->S = p->L0prev->R;
                  p->L0prev->R = p->curr;
               }
               else if(p->curr->idx < J[1])
               {
                  p->prev = p->curr;
                  if(!L0_FindInR(p,J[1]))
                  {                                   /* Add2S(J[1],lpp); */
                     lpp = (LPdata*)calloc(1,sizeof(LPdata));
                     LPdata_Init(lpp,N,I,X,A);
                     p->LP1->next = (LPPL*)calloc(1,sizeof(LPPL));
                     LPPL_Init(p->LP1->next,lpp,p->LP1->next);
                     p->curr = (IndexNode*)calloc(1,sizeof(IndexNode));
                     IndexNode_Init(p->curr,J[1],lpp);
                     p->curr->S = p->prev->S;
                     p->prev->S = p->curr;
                  }
                  else return;     /* J[1] is in the right branch of J[0] */
               }
            }
            else
            {
               lpp = (LPdata*)calloc(1,sizeof(LPdata));
               LPdata_Init(lpp,N,I,X,A);
               p->LP1->next = (LPPL*)calloc(1,sizeof(LPPL));                  
               LPPL_Init(p->LP1->next,lpp,p->LP1->next);
               p->L0curr->R = (IndexNode*)calloc(1,sizeof(IndexNode));
               IndexNode_Init(p->L0curr->R,J[1],lpp); /* Add2E(J[1],lpp); */
            }
         }
      }
      else                                       /* add J[i] after L0prev */
      {
         if(i==0)
         {
            lpp = (LPdata*)calloc(1,sizeof(LPdata));
            LPdata_Init(lpp,N,I,X,A);
            p->LP1->next = (LPPL*)calloc(1,sizeof(LPPL));                  
            LPPL_Init(p->LP1->next,lpp,p->LP1->next);
            p->L0curr = (L0IdxNode*)calloc(1,sizeof(L0IdxNode));
            L0IdxNode_Init(p->L0curr,J[i]);
            p->L0curr->D = p->L0prev->D;
            p->L0prev->D = p->L0curr; 
            p->L0curr->R = (IndexNode*)calloc(1,sizeof(IndexNode));
            IndexNode_Init(p->L0curr->R,J[1],lpp);    /* Add2E(J[1],lpp); */
         }
         else
         {
            p->L0curr = (L0IdxNode*)calloc(1,sizeof(L0IdxNode));
            L0IdxNode_Init(p->L0curr,J[i]);
            p->L0curr->D = p->L0prev->D;
            p->L0prev->D = p->L0curr; 
         }
      }
      p->L0prev = p->L0curr;
   }
}

void L0_IML_Del ( L0_IML *p )
{
   L0_Free(p);
   free(p->LP1);             /* L0head is removed by Migrate when empty. */
}

void L0_Free ( L0_IML *li )
{
   LPPL *P = li->LP1->next;                            /* Clean LP1 link */
   int i;

   while(P)
   {
      li->LP1->next = P->next;
      free(P->addr->JJJ);
      free(P->addr->xxx);
      for(i=0; i<P->addr->dim; i++) free((P->addr->INV)[i]);
      free(P->addr->INV);                              /* delete P->addr */
      free(P);
      P = li->LP1->next;
   }    /* Multi index links and L0head would be removed by Migrate(-) ! */
}
