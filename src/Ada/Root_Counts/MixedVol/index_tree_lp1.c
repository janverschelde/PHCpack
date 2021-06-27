/* The file "index_tree_lp1.c" contains the definitions of the functions
 * whose prototypes are documented in the file "index_tree_lp1.h". */

#include <stdio.h>
#include <stdlib.h>
#include "index_tree_lp1.h"

void LPdata_Init ( LPdata *p, int n, int *J, double *x, double **A )
{
   int i,j;

   p->dim = n;
   p->JJJ = (int*)calloc(n,sizeof(int));
   p->xxx = (double*)calloc(n,sizeof(double));
   p->INV = (double**)calloc(n,sizeof(double*));

   for(i=0; i<n; i++)
   {
      p->JJJ[i] = J[i];
      p->xxx[i] = x[i];
      p->INV[i] = (double*)calloc(n,sizeof(double));
      for(j=0; j<n; j++) p->INV[i][j] = A[i][j];
   }
}

void IndexNode_Init ( IndexNode *p, int i, LPdata *ptr )
{
   p->idx = i;
   p->info = ptr;
   p->S = 0; 
}

void LPPL_Init ( LPPL *p, LPdata *A, LPPL *N )
{
   p->addr = A;
   p->next = N;
}

void IT_LP_Init ( IT_LP *p, int nSpt, int *type )
{
   int i,j,sum,itmp=0;

   for(i=0; i<nSpt; i++) itmp += type[i];
   p->MaxLevels = itmp + nSpt + 1;        /* "+1" is for unused level 0 */
   p->CurLevel = 1;

   p->NP    = (int*)calloc(p->MaxLevels,sizeof(int));
   p->DIM   = (int*)calloc(p->MaxLevels,sizeof(int));
   p->cell  = (int*)calloc(p->MaxLevels,sizeof(int));
   p->InSpt = (int*)calloc(p->MaxLevels,sizeof(int));
   p->minNP = (int*)calloc(p->MaxLevels,sizeof(int));

   p->LP = (LPPL**)calloc(p->MaxLevels,sizeof(LPPL*));
   for(i=0; i<p->MaxLevels; i++)
   {
      p->LP[i]=(LPPL*)calloc(1,sizeof(LPPL));
      LPPL_Init(p->LP[i],0,0);  
   }

   p->IT   = (IndexNode**)calloc(p->MaxLevels,sizeof(IndexNode*));
   p->last = (IndexNode**)calloc(p->MaxLevels,sizeof(IndexNode*));
   
   memset(p->IT,    0, p->MaxLevels * sizeof(IndexNode*));
   memset(p->last,  0, p->MaxLevels * sizeof(IndexNode*));
   memset(p->NP,    0, p->MaxLevels * sizeof(int));
   memset(p->DIM,   0, p->MaxLevels * sizeof(int)); 
   memset(p->cell,  0, p->MaxLevels * sizeof(int)); 
   memset(p->InSpt, 0, p->MaxLevels * sizeof(int)); 
   memset(p->minNP, 0, p->MaxLevels * sizeof(int));

   sum = 0;
   p->DIM[sum] = itmp++;
   for(i=0; i<nSpt; i++)
   {
      p->minNP[sum] = type[i]+1;
      p->InSpt[sum] = i;
      for(j=1; j<=type[i]; j++)
      {
         p->DIM[sum+j] = --itmp;
         p->minNP[sum+j] = type[i]+1-j;
      }
      sum += type[i]+1;
      if(sum < p->MaxLevels) p->DIM[sum] = itmp;
   }
   
   p->IT[1]=(IndexNode*)calloc(1,sizeof(IndexNode));
   IndexNode_Init(p->IT[1],-1,0);
   p->last[1] = p->prev = p->curr = p->IT[1]; /* points to dummy node */
   p->NP[1] = 1;
}

int IT_IsEmpty ( IT_LP *p ) 
{ 
   return (p->CurLevel<1);
}

int IT_IsFull ( IT_LP *p )  
{
   return (p->CurLevel+1 >= p->MaxLevels);
}
     
int IT_Level ( IT_LP *p )
{
   return p->CurLevel;
}

int IT_CurLPdim ( IT_LP *p )
{
   return p->DIM[p->CurLevel];
} 

int *IT_Cell ( IT_LP *p )
{
   return p->cell + 1;
}

int IT_CurSptIdx ( IT_LP *p )
{
   return p->InSpt[p->CurLevel];
}

int IT_MinNumPt ( IT_LP *p )
{
   return p->minNP[p->CurLevel];
}

int IT_NumPt ( IT_LP *p )
{
   return p->NP[p->CurLevel];
}

IndexNode *IT_FixedIdxNdPtr ( IT_LP *p )
{
   return (p->CurLevel>=0 ? p->IT[p->CurLevel] : 0);
}

void IT_StepBack ( IT_LP *p )
{
   (p->NP)[p->CurLevel--] = 0;
}

int IT_Find ( IT_LP *p, int IDX )
{
   for(p->curr=p->prev->S;
       p->curr != p->last[p->CurLevel]->S; p->curr=p->curr->S)
   {
      if(IDX <= p->curr->idx)
      {
         if(IDX == p->curr->idx)
            return 1; 
         else
            return 0; 
      }
      p->prev = p->prev->S;
   }
   return 0; 
} 

IndexNode *IT_ResetCurLevelTo1 ( IT_LP *p )
{
   p->prev = p->IT[1];
   p->curr = p->prev->S;

   while(p->curr)
   {
      p->prev->S = p->curr->S;
      free(p->curr);
      p->curr = p->prev->S;
   }
   p->NP[1]=1;
   p->CurLevel=1;

   return p->IT[1];
}

void IT_RenewNP1 ( IT_LP *p ) 
{
   for(p->prev=p->IT[1]; p->prev->S; p->prev=p->prev->S)
      ++(p->NP[1]);

   p->last[1] = p->prev;
   p->cell[1] = p->IT[1]->idx;
}

int IT_NextLevel ( IT_LP *p ) 
{
   IndexNode *tmp;

   if(p->CurLevel+1 >= p->MaxLevels) return 0;                  /* case 0 */
   else
   {
      if(p->NP[p->CurLevel] <= p->minNP[p->CurLevel]) return 0; /* case 1 */
                           /* now IT[CurLevel] has enough points to go on */
      if(p->IT[p->CurLevel+1])     /* backtracking (next level non-empty) */
      {
         tmp = p->IT[p->CurLevel+1];
               p->IT[p->CurLevel+1] = tmp->S;
         tmp->S = (p->last[p->CurLevel])->S;
                  (p->last[p->CurLevel])->S = tmp;
         tmp = p->IT[p->CurLevel]->S;
               p->IT[p->CurLevel]->S = tmp->S;
         tmp->S = p->IT[p->CurLevel+1];
                  p->IT[p->CurLevel+1] = tmp;
      }
      else 
      {
         tmp = p->IT[p->CurLevel]->S;
               p->IT[p->CurLevel+1] = tmp;
         p->IT[p->CurLevel]->S = tmp->S;
                                 tmp->S = 0;
      }
      if(p->NP[p->CurLevel] == 2) p->last[p->CurLevel] = p->IT[p->CurLevel];
      --p->NP[p->CurLevel++];
      ++p->NP[p->CurLevel];
      p->last[p->CurLevel] = p->IT[p->CurLevel];
      p->curr = p->IT[p->CurLevel];
      p->cell[p->CurLevel] = p->curr->idx;
      p->LPlast = p->LP[p->CurLevel];
      return 1;                                                 /* case 2 */
   }
}

void IT_Add1
 ( IT_LP *p, int n, int *J, int nn, int *JJ, double *X, double **A )
{
   int i,j,k;
   int AddIdx = 0;
   LPdata *lpp;
   IndexNode *ptr;

   p->prev = p->IT[p->CurLevel];                /* moved out of for loop */
                                              /* Create or Overwrite LP: */
   for(i=0; i<n; i++)      /* search for the 1st idx of J not in level 0 */
   {
      if(!IT_Find(p,J[i]))                        /* add J[i] after prev */
      {
         if(p->LPlast->next)
         {
            lpp = p->LPlast->next->addr;
            for(k=0; k<nn; k++)
            {
               lpp->JJJ[k] = JJ[k];
               lpp->xxx[k] = X[k];
               for(j=0; j<nn; j++) lpp->INV[k][j] = A[k][j];
            }
         }
         else
         {
            lpp = (LPdata*)calloc(1,sizeof(LPdata));
            LPdata_Init(lpp,nn,JJ,X,A);
            p->LPlast->next = (LPPL*)calloc(1,sizeof(LPPL));
            LPPL_Init(p->LPlast->next,lpp,0);
         }
         p->LPlast = p->LPlast->next;
         AddIdx = 1;
         break;
      }
   }
   while(AddIdx)
   {
      if(!p->last[p->CurLevel]->S)              /* IT[CurLevel] is full */
      {
         p->curr = (IndexNode*)calloc(1,sizeof(IndexNode));
         IndexNode_Init(p->curr,J[i],lpp);
         p->curr->S = p->prev->S;
         p->prev = p->prev->S = p->curr;
         if( (p->last[p->CurLevel])==(p->IT[p->CurLevel])
             || p->last[p->CurLevel]->idx < J[i] )
            p->last[p->CurLevel] = p->curr;
      }                                    /*  spare slots for LP data  */
      else if(p->last[p->CurLevel] == p->prev) /* after *last[CurLevel] */
      {
         ptr = p->prev->S;
         ptr->idx = J[i];
         ptr->info = lpp;
         p->last[p->CurLevel] = ptr;
         p->prev = ptr;
      }
      else                                     /* intermediate position */
      {
         ptr = p->last[p->CurLevel]->S;
         ptr->idx = J[i];
         ptr->info = lpp;
         p->prev->S = ptr;
         p->last[p->CurLevel]->S = ptr->S;
         ptr->S = p->curr;
         p->prev = ptr;
      }
      ++p->NP[p->CurLevel];
      AddIdx = 0;
      for(++i; i<n; i++) 
      {
         if( !IT_Find(p,J[i]) ) 
         { 
            AddIdx = 1;
            break;
         }
      }
   }
}

void IT_Add2 
 ( IT_LP *p, int oneidx, int nn, int *JJ, double *X, double **A ) 
{
   int i,j;
   LPdata *lpp;
   IndexNode *ptr ;

   p->prev = p->IT[p->CurLevel]; 
   if(!IT_Find(p,oneidx))                       /* add oneidx after prev */ 
   {
      if(p->LPlast->next)
      {
         lpp = p->LPlast->next->addr;
         for(i=0; i<nn; i++)
         {
            lpp->JJJ[i] = JJ[i];
            lpp->xxx[i] = X[i];
            for(j=0; j<nn; j++) lpp->INV[i][j] = A[i][j];
         }
      }
      else
      {
         lpp = (LPdata*)calloc(1,sizeof(LPdata));
         LPdata_Init(lpp,nn,JJ,X,A);
         p->LPlast->next = (LPPL*)calloc(1,sizeof(LPPL)); 
         LPPL_Init(p->LPlast->next,lpp,0);
      }
      p->LPlast = p->LPlast->next;

      if(!p->last[p->CurLevel]->S)               /* IT[CurLevel] is full */
      {
         /* LPdata* LP = new LPdata(nn,JJ,X,A);*/
         p->curr=(IndexNode*)calloc(1,sizeof(IndexNode));
         IndexNode_Init(p->curr,oneidx,lpp);
         p->curr->S = p->prev->S;
         p->prev = p->prev->S = p->curr;
         if( p->last[p->CurLevel]->idx < oneidx ) 
            p->last[p->CurLevel] = p->curr;
      }                               /* have spare slots for LP data .... */
      else if(p->last[p->CurLevel] == p->prev) /* after *last[CurLevel] */
      {
         ptr = p->prev->S;
         ptr->idx = oneidx;
         ptr->info = lpp;
         
         p->last[p->CurLevel] = ptr;
         p->prev = ptr;
      }
      else       /* intermediate position before last, though last->S !=0 */
      {
         ptr = p->last[p->CurLevel]->S;
         ptr->idx = oneidx;
         ptr->info = lpp;
         p->prev->S = ptr;
         p->last[p->CurLevel]->S = ptr->S;
         ptr->S = p->curr;
         p->prev = ptr;
      }
      ++p->NP[p->CurLevel];
   }
}

void IT_LP_DEL ( IT_LP *p )
{
   IT_FreeIT(p);
   IT_FreeLP(p); 
   free(p->IT);
   free(p->last);
   free(p->LP);
   free(p->NP);
   free(p->DIM);
   free(p->cell);
   free(p->InSpt);
   free(p->minNP);
}

void IT_FreeIT ( IT_LP *p ) 
{
   p->CurLevel = p->MaxLevels-1;
   while( p->CurLevel > 1 )
   {
      p->prev = p->IT[p->CurLevel];
      p->curr = p->prev->S;
      while(p->curr)
      {
         p->prev->S = p->curr->S;
         free(p->curr);
         p->curr = p->prev->S;
      }
      --p->CurLevel;
   }
   for(p->CurLevel=0; p->CurLevel<p->MaxLevels; p->CurLevel++)
      free(p->IT[p->CurLevel]);
}

void IT_FreeLP ( IT_LP *p ) 
{
   int i;
   LPdata *lpp;

   p->CurLevel = p->MaxLevels-1;
   while(p->CurLevel > 1)
   {
      p->LPlast = (p->LP[p->CurLevel])->next;
      while(p->LPlast)
      {
         p->LP[p->CurLevel]->next = p->LPlast->next;
         lpp = p->LPlast->addr;
         if(lpp)
         {
            free(lpp->JJJ);
            free(lpp->xxx);
            for(i=0; i<lpp->dim; i++) free(lpp->INV[i]);
            free(lpp->INV);
         }
         free(p->LPlast);
         p->LPlast = p->LP[p->CurLevel]->next;
      }
      --p->CurLevel;
   }
   for(p->CurLevel=0; p->CurLevel<p->MaxLevels; p->CurLevel++)
      free(p->LP[p->CurLevel]);
}
