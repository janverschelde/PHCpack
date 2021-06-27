/* The file "form_lp1.c" contains the definitions of the two functions
 * whose prototypes are documented in the file "form_lp1.h". */

#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "form_lp1.h"

void form_LP
 ( int nVar, int nSpt, int *SptType, int *SptIdx, int **RelTab, 
   double **ElmMtrx, int *NonZero, int *Lvl2CoDim, double ***A,
   int Lvl, LPdata *ptr, int *Cell, int *Cell_Orig, int *Lvl2LPdim, 
   int *Lvl2Spt, int *FixFrstPt, int *FixLstPt, int *MinNumPt,
   int **PtIn, int *Strt1Pt, int *End1Pt, int **ToOrig,
   int **Pre2Cur, int *LPdim, int *FrstPt, int *LstPt, 
   int *Cur2Pre, double *x, double **Binv, int *Bidx, int *info )
{
   int ibrk,ell,i,j,jj,i1,j1,k,kout,EndOldPart,FixPt_Orig,FixPt,nPtIn;
   double smallnum=1.0e-12,dtmp;
   double Elm1Var[nVar+1];
   const int Lvl1 = Lvl-1;
   const int ynFixLstPt = FixLstPt[Lvl];
   const int ynFixFrstPt = FixFrstPt[Lvl];
   const int CurSpt = Lvl2Spt[Lvl];         /* index of current support */
   const int CoDim = Lvl2CoDim[Lvl];
   const int PreLPdim = Lvl2LPdim[Lvl-1];

   *LPdim = Lvl2LPdim[Lvl];       /* Step 1 : search for LP constraints */
   FixPt = Cell[Lvl];                             /* last point fixed */
   if(Lvl == 2)
   {
      FixPt_Orig = FixPt;
      FixPt = Pre2Cur[Lvl1][FixPt];
   }
   else
      FixPt_Orig = ToOrig[Lvl1][FixPt];

   nPtIn = -1;                /* #points involved in the constraints - 1 */
   for(i1=0; i1<FrstPt[Lvl1]; i1++)
   {
      i = ToOrig[Lvl1][i1];
      /* if ( RelTab[FixPt_Orig][i] && PtIn[Lvl1][i1] ) */ 
      if(RelTab[FixPt_Orig][i])  /* always in */
      {
         ToOrig[Lvl][++nPtIn] = i;
         Cur2Pre[nPtIn] = i1;
         Pre2Cur[Lvl][i1] = nPtIn;
      }
      else
         Pre2Cur[Lvl][i1] = -1;   
   }
   FrstPt[Lvl] = nPtIn + 1;
   for(i1=FrstPt[Lvl1]; i1<FixPt; i1++)
   {
      i = ToOrig[Lvl1][i1];
      if(RelTab[FixPt_Orig][i] && PtIn[Lvl1][i1])
      {
         ToOrig[Lvl][++nPtIn] = i;
         Cur2Pre[nPtIn] = i1;
         Pre2Cur[Lvl][i1] = nPtIn;
      }
      else
         Pre2Cur[Lvl][i1] = -1;
   }

   *Strt1Pt = nPtIn + 1;

   for(i1=FixPt+1; i1<=LstPt[Lvl1]; i1++)
   {
      i = ToOrig[Lvl1][i1];
      if(RelTab[FixPt_Orig][i] && PtIn[Lvl1][i1])
      {
         ToOrig[Lvl][++nPtIn] = i;
         Cur2Pre[nPtIn] = i1;
         Pre2Cur[Lvl][i1] = nPtIn;
      }
      else
         Pre2Cur[Lvl][i1] = -1;
   }

   if(!ynFixLstPt && nPtIn-*Strt1Pt < MinNumPt[Lvl]-1)
   {
      *info = 1;
      return;    /* not enough pts to extend; */
   }

   *End1Pt = EndOldPart = LstPt[Lvl] = nPtIn;

   if(Lvl == 2)
      Cell_Orig[2] = Cell[2];
   else
      Cell_Orig[Lvl] = ToOrig[Lvl-1][Cell[Lvl]];

   if(ynFixLstPt)       /* The last pt fixed is last pt needed; */
   {                    /* Extend to next spt.*/

      *Strt1Pt = FrstPt[Lvl] = nPtIn + 1;
      for(i=SptIdx[CurSpt+1]; i<SptIdx[CurSpt+2]; i++)
      {
         ibrk = 0;
         for(j=1; j<=Lvl; j++)
         {
            if(!RelTab[i][Cell_Orig[j]])
            {
               ibrk = 1;
               break;
            }
         }
         if(!ibrk)
         {
            ToOrig[Lvl][++nPtIn] = i;
            Cur2Pre[nPtIn] = i;
            Pre2Cur[Lvl][i] = nPtIn;
         }
         else
            Pre2Cur[Lvl][i] = -1;
      }
      *End1Pt = LstPt[Lvl] = nPtIn;

      if(*End1Pt-*Strt1Pt < SptType[CurSpt+1])   
      { 
         *info = 1;
         return;    /* not enough pts; need irep(j)+1 pts */
      }
   }
   /* Step 2 : eliminate a variable by the eq. of the last fixed point */
   /*          form the equation used for elimination */

   for(i=0; i<PreLPdim; i++) Elm1Var[i] = A[Lvl1][FixPt][i];
   Elm1Var[nVar] = A[Lvl1][FixPt][nVar];

   /* eliminate the kout-th variable and shift columns of a to left */

   if(!ynFixFrstPt)    /* still in same support */
   {
      ibrk = 0;
      for(kout=PreLPdim-1; kout>=0; kout--)
         if(fabs(Elm1Var[kout]) > smallnum)
         {
            ibrk=1;
            break;
         }
      if(!ibrk) 
      {
         printf("form_LP : cannot find a variable to eliminate.\a\n");
         abort();
      }
      for(i=0; i<kout; i++) Elm1Var[i] /= Elm1Var[kout];
      Elm1Var[nVar] /= Elm1Var[kout];
      NonZero[CoDim] = kout;
      for(i=0; i<kout; i++) ElmMtrx[CoDim][i] = Elm1Var[i];
      ElmMtrx[CoDim][nVar] = Elm1Var[nVar];

      for(i=0; i<=EndOldPart; i++)
      {
         j = Cur2Pre[i];
         if(fabs(A[Lvl1][j][kout]) > smallnum)  /* to eliminate */
         {
            for(j1=0; j1<kout; j1++)
               A[Lvl][i][j1] = A[Lvl1][j][j1] 
                                  - A[Lvl1][j][kout]*Elm1Var[j1];
            A[Lvl][i][nVar] = A[Lvl1][j][nVar] 
                                  - A[Lvl1][j][kout]*Elm1Var[nVar];
         } 
         else
         {
            for(j1=0; j1<kout; j1++) A[Lvl][i][j1] = A[Lvl1][j][j1];
            A[Lvl][i][nVar] = A[Lvl1][j][nVar];
         }
 
         for(j1=kout; j1<PreLPdim-1; j1++)  /* to shift; */
            A[Lvl][i][j1] = A[Lvl1][j][j1+1];
      }
      if(ynFixLstPt)   /* set the values for the variable alpha0 */
      {
         for(i=0; i<=EndOldPart; i++)   /* constraints without alpha0 */
            A[Lvl][i][PreLPdim-1] = 0.0;
         for(i=FrstPt[Lvl]; i<=nPtIn; i++)
	 {                              /* constraints with alpha0 */
            j = ToOrig[Lvl][i];
            for(k=0; k<nVar; k++) A[Lvl][i][k] = A[0][j][k];
            A[Lvl][i][nVar] = A[0][j][nVar];
         }
         for(i1=1; i1<=CoDim; i1++)
         {
            kout = NonZero[i1];
            for(jj=0; jj<nVar+1; jj++) Elm1Var[jj] = ElmMtrx[i1][jj];
            for(i=FrstPt[Lvl]; i<=nPtIn; i++)
            {
               if(fabs(A[Lvl][i][kout]) > smallnum)  /* to eliminate */
               {
                  for(j1=0; j1<kout; j1++)
                     A[Lvl][i][j1] -= A[Lvl][i][kout]*Elm1Var[j1];
                  A[Lvl][i][nVar] -= A[Lvl][i][kout]*Elm1Var[nVar];
               }
               for(j1=kout; j1<nVar-i1; j1++)  /* to shift; */
                  A[Lvl][i][j1] = A[Lvl][i][j1+1];
            }
         }  /* kout = NonZero[CoDim] will be used below! */
         for(i=FrstPt[Lvl]; i<=*End1Pt; i++)
                                              /* constraints with alpha0 */
            A[Lvl][i][PreLPdim-1] = 1.0;   /* For added variable alpha0 */
      }
   }
   else  /* fixing the first point from a support, eliminate alpha0 */
   {
      kout = PreLPdim - 1;
      if(fabs(Elm1Var[kout]) > smallnum)  /* eliminate alpha0 */
      {
         for(i=0; i<kout; i++) Elm1Var[i] /= Elm1Var[kout];
         Elm1Var[nVar] /= Elm1Var[kout];
         for(i=0; i<FrstPt[Lvl]; i++)
	 {                               /* copy the previous a, b saved */
            j = Cur2Pre[i];
            for(j1=0; j1<PreLPdim-1; j1++)
               A[Lvl][i][j1] = A[Lvl1][j][j1];
            A[Lvl][i][nVar] = A[Lvl1][j][nVar];
         }
         for(i=FrstPt[Lvl]; i<=*End1Pt; i++)
         {                                       /* eliminate alpha0 */ 
            j = Cur2Pre[i];
            for(j1=0; j1<PreLPdim-1; j1++)
               A[Lvl][i][j1] = A[Lvl1][j][j1] - Elm1Var[j1];
            A[Lvl][i][nVar] = A[Lvl1][j][nVar] - Elm1Var[nVar];
         }
      }
   }
   ibrk = 0;
   if(Lvl == 2)
      for(i=0; i<PreLPdim; i++)
      {
         Bidx[i] = ptr->JJJ[i];
         if(Bidx[i] == FixPt_Orig)
         {
            k = i;
            ibrk = 1;
            break;
         }
      }
   else
      for(i=0; i<PreLPdim; i++)
      {
         Bidx[i] = ptr->JJJ[i];
         if(Bidx[i] == FixPt)
         {
            k = i;
            ibrk = 1;
            break;
         }
      }
   if(!ibrk)
   {
      printf("form_LP : no index match for reused info\a\n");
      abort();
   }
   for(i=k+1; i<PreLPdim; i++) Bidx[i-1] = ptr->JJJ[i];

   if(Lvl > 2)
   {
      for(i=0; i<PreLPdim-1; i++)
         if(Bidx[i]>-1) 
            Bidx[i] = Pre2Cur[Lvl][Bidx[i]];  /* may = -1 if eliminated */
   }
   else 
   {
      for(i=0; i<PreLPdim-1; i++)
         if(Bidx[i]>-1) 
         {
            Bidx[i] = Pre2Cur[1][Bidx[i]]; /* from level 0 to level 1 */
            if(Bidx[i]>-1) 
               Bidx[i] = Pre2Cur[Lvl][Bidx[i]]; /* level 1 to level 2 */
         }
   }
   for(i=0; i<kout; i++) x[i] = ptr->xxx[i];
   for(i=kout+1; i<PreLPdim; i++) x[i-1] = ptr->xxx[i];

   for(j=0; j<k; j++)
   {
      for(i=0; i<kout; i++) Binv[j][i] = ptr->INV[j][i];
      for(i=kout+1; i<PreLPdim; i++) Binv[j][i-1] = ptr->INV[j][i];
   }
   for(j=k+1; j<PreLPdim; j++)
   {
      for(i=0; i<kout; i++) Binv[j-1][i] = ptr->INV[j][i]; 
      for(i=kout+1; i<PreLPdim; i++) Binv[j-1][i-1] = ptr->INV[j][i];
   }
   if(ynFixLstPt)   /* the 1st round 1-pt test for next supp */
   {
      x[*LPdim-1] = DBL_MAX;  /* The variable alpha0 */
      ell = -1;
      for(i=*Strt1Pt; i<=*End1Pt; i++)
      {
         dtmp = A[Lvl][i][nVar];
         for(i1=0; i1<*LPdim-1; i1++) dtmp -= A[Lvl][i][i1]*x[i1];

         if(x[*LPdim-1] > dtmp) 
         {
            ell = i;
            x[*LPdim-1] = dtmp;
         }
      }
      if(ell < 0)
      {
         printf("FormLP : no ell\a\n");
         abort();
      }
      Bidx[*LPdim-1] = ell; /* constraint becomes last element of base */
      for(i=0; i<*LPdim; i++)            /*  update the inverse matrix */
         Binv[*LPdim-1][i] = 0.0;

      for(j=0; j<*LPdim-1; j++)   
      {
         dtmp = A[Lvl][ell][0]*Binv[j][0];
         for(i=1; i<*LPdim-1; i++) dtmp += A[Lvl][ell][i]*Binv[j][i];

         Binv[j][*LPdim-1] = -dtmp;
      }
      Binv[*LPdim-1][*LPdim-1] = 1.0;
   }
   for(i=FrstPt[Lvl]; i<*Strt1Pt; i++) PtIn[Lvl][i] = 1;
   for(i=*Strt1Pt; i<=*End1Pt; i++) PtIn[Lvl][i] = 0;
   *info = 0;
}

void form_LP1
 ( int nVar, int nSpt, int *SptType, int *SptIdx, int **RelTab,
   double ***A, int Lvl, int *Cell, int *Lvl2LPdim, int *FixLstPt,
   int *MinNumPt, int **PtIn, int **ToOrig, int **Pre2Cur,
   int *FrstPt, int *LstPt, int *info )
{
   int i,j,j1,FixPt,nPtIn,Strt1Pt;
   const int Lvl1 = Lvl - 1;
   const int ynFixLstPt = FixLstPt[Lvl];
   const int PreLPdim = Lvl2LPdim[Lvl-1];

   FixPt = Cell[Lvl];                              /* last point fixed */
   nPtIn = -1;         /* #points involved to form the constraints - 1 */
   FrstPt[Lvl] = 0;
            
   for(i=SptIdx[0]; i<FixPt; i++)
      if(RelTab[FixPt][i])
      {
         ToOrig[Lvl][++nPtIn] = i;
         Pre2Cur[Lvl][i] = nPtIn;  /* need for saved Bidx for next lvl */
      }
      else
         Pre2Cur[Lvl][i] = -1;
   Pre2Cur[Lvl][FixPt] = -1;
   Strt1Pt = nPtIn + 1;
   for(i=FixPt+1; i<SptIdx[1]; i++)
      if(RelTab[FixPt][i])
      {
         ToOrig[Lvl][++nPtIn] = i;
         Pre2Cur[Lvl][i] = nPtIn;   /* need for saved Bidx for next lvl */
      }
      else
         Pre2Cur[Lvl][i] = -1;
     
   if(!ynFixLstPt && nPtIn-Strt1Pt < MinNumPt[Lvl]-1)
   {
      *info = 1;
      return;   /* not enough pts to extend; */
   }
   LstPt[Lvl] = nPtIn;   
                                                 /* To form the matrix */
   for(i=0; i<=LstPt[Lvl]; i++)         /* copy the previous a, b saved */
   {
      j = ToOrig[Lvl][i];
      for(j1=0; j1<PreLPdim; j1++)
         A[Lvl][i][j1] = A[Lvl1][j][j1] - A[Lvl1][FixPt][j1];
      A[Lvl][i][nVar] = A[Lvl1][j][nVar] - A[Lvl1][FixPt][nVar];
   }
   for(i=0; i<=LstPt[Lvl]; i++) PtIn[Lvl][i] = 1;
   *info = 0;
}
