/* The file "mixed_volume1.c" contains the definitions of the prototypes
 * in the file "mixed_volume1.h". */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>  
#include "mixed_volume1.h"
#include "index_tree_lp1.h"
#include "form_lp1.h"
#include "relation_table1.h"
#include "one_level_lp1.h"

void write_relation_table ( int n, int **RelTab )
{
   int i,j;

   for(i=0; i<n; i++)
   {
      for(j=0; j<n; j++)
         printf(" %d",RelTab[i][j]);
      printf("\n");
   }
}

void MixedVol
 ( int nVar, int nSpt, int CellSize, int *SptType, int *SptIdx, int **Spt,
   double *lft, int *nbCells, CellStack *MCells, int *MVol )
{  
   int Vol,tolVol,i,j,k,m,info,Strt1Pt,End1Pt,nMCells,LPdim,Lvl;

   int *labels = (int*)calloc(CellSize,sizeof(int));
   int *NonZero = (int*)calloc(CellSize,sizeof(int));
   int *Lvl2CoDim = (int*)calloc(CellSize,sizeof(int));
   int *Cell = (int*)calloc(CellSize,sizeof(int));
   int *Cell_Orig = (int*)calloc(CellSize,sizeof(int));
   int *Lvl2LPdim = (int*)calloc(CellSize,sizeof(int));
   int *Lvl2Spt = (int*)calloc(CellSize,sizeof(int));
   int *MinNumPt = (int*)calloc(CellSize,sizeof(int));
   int *FixFrstPt = (int*)calloc(CellSize,sizeof(int));
   int *FixLstPt = (int*)calloc(CellSize,sizeof(int));
   int *FrstPt = (int*)calloc(CellSize,sizeof(int));
   int *LstPt = (int*)calloc(CellSize,sizeof(int));
   int *Cur2Pre = (int*)calloc(SptIdx[nSpt],sizeof(int));
   int *Bidx = (int*)calloc(nVar,sizeof(int));
   int **PtIn = (int**)calloc(CellSize,sizeof(int*));
   int **ToOrig = (int**)calloc(CellSize,sizeof(int*));
   int **Pre2Cur = (int**)calloc(CellSize,sizeof(int*));
   int **RelTab = (int**)calloc(SptIdx[nSpt],sizeof(int*));
   double *x = (double*)calloc(nVar,sizeof(double));
   double **ElmMtrx = (double**)calloc(CellSize,sizeof(double*));
   double **Binv = (double**)calloc(nVar,sizeof(double*));
   double ***A = (double***)calloc(CellSize,sizeof(double**));

   IT_LP *ItLp = (IT_LP*)calloc(1,sizeof(IT_LP));
   L0_IML *L0 = (L0_IML*)calloc(1,sizeof(L0_IML));
   LPdata *ptr;

   *nbCells = 0; *MVol = 0; Vol = 0; tolVol = 0;

   IT_LP_Init(ItLp,nSpt,SptType);
   L0_IML_Init(L0);   
         
   for(i=0; i<CellSize; i++)
   {
      PtIn[i] = (int*)calloc(SptIdx[nSpt],sizeof(int));
      ToOrig[i] = (int*)calloc(SptIdx[nSpt],sizeof(int));
      Pre2Cur[i] = (int*)calloc(SptIdx[nSpt],sizeof(int));
      ElmMtrx[i] = (double*)calloc(nVar+1,sizeof(double));
   }
   for(i=0; i<nVar; i++) Binv[i] = (double*)calloc(nVar,sizeof(double));

   for(i=0; i<SptIdx[nSpt]; i++)
      RelTab[i] = (int*)calloc(SptIdx[nSpt],sizeof(int));
   RelTable(nVar,nSpt,Spt,SptIdx,lft,RelTab,L0);
  /* write_relation_table(SptIdx[nSpt],RelTab); */
   
   info = SptIdx[nSpt];                                  /* #rows in A */
   for(i=0; i<nSpt-1; i++)
   {
      for(j=0; j<SptType[i]; j++) info += SptIdx[i+1];
      info += SptIdx[i+2];
   }
   for(j=0; j<SptType[nSpt-1]; j++) info += SptIdx[nSpt];

   for(i=0; i<CellSize; i++)
   {
      A[i] = (double**)calloc(info,sizeof(double*));
      for(j=0; j<info; j++)
         A[i][j] = (double*)calloc((nVar+1),sizeof(double));
   }
   info = 0;

   /* Note that the level 0 is the artificial head */

   Lvl2LPdim[0] = nVar + 1;  /* for conveniece, should = nVar, see below */
   Lvl2CoDim[0] = 0;         /* #variables eliminated. [0] not used */
   k = 0;
   for(i=0; i<nSpt; i++)
   {
      ++k;
      Lvl2LPdim[k] = Lvl2LPdim[k-1] - 1;
      Lvl2CoDim[k] = Lvl2CoDim[k-1];
      Lvl2Spt[k] = i;
      MinNumPt[k] = SptType[i];
      FixFrstPt[k] = 1;
      FixLstPt[k] = 0;
      for(j=1; j<=SptType[i]; j++)
      {
         ++k;
	 if(k == CellSize) break;
         Lvl2LPdim[k] = Lvl2LPdim[k-1] - 1;
         Lvl2CoDim[k] = Lvl2CoDim[k-1] + 1;
         Lvl2Spt[k] = i;
         MinNumPt[k] = MinNumPt[k-1] - 1;
         FixFrstPt[k] = 0;
         FixLstPt[k] = 0;
      }
      if(k == CellSize) break;
      Lvl2LPdim[k] += 1;                     /* add alpha0 for next spt */
      if(i < nSpt-1) MinNumPt[k] = SptType[i+1] + 1;
      FixLstPt[k] = 1;
   }      
   Lvl2LPdim[0] = nVar; /* =nVar+1-1, in RelTab, add alpha0, fix 1 point */

   for(i=SptIdx[0]; i<SptIdx[nSpt]; i++)               /* define A[0] */
   {
      for(j=0; j<nVar; j++) A[0][i][j] = -Spt[i][j];
      A[0][i][nVar] = lft[i];                          /* bIdx = nVar */
   }
   for(i=SptIdx[0]; i<SptIdx[1]; i++) PtIn[0][i] = 1;  /* define PtIn[0] */

   while(L0_Migrate(L0,IT_ResetCurLevelTo1(ItLp)))
   {
      IT_RenewNP1(ItLp);
      while(!IT_IsEmpty(ItLp))
      {
         if(!IT_IsFull(ItLp))
         {
            Lvl = IT_Level(ItLp);  /* level-0 is the artificial head */
            Cell[Lvl] = IT_FixedIdxNdPtr(ItLp)->idx; 
               /* The index of the point to be fixed */
            ptr = IT_FixedIdxNdPtr(ItLp)->info;
               /* The ptr to saved x, Binv and jj */

            if(Lvl <= 1)   /* eliminate and form new system Ax <= b */
            {
               Cell_Orig[Lvl] = Cell[Lvl];
               form_LP1(nVar,nSpt,SptType,SptIdx,RelTab,A,Lvl,Cell,
                        Lvl2LPdim,FixLstPt,MinNumPt,PtIn,
                        ToOrig,Pre2Cur,FrstPt,LstPt,&info);
            }
            else  /* eliminate, form new Ax <= b and solve the LP */
            {
               form_LP(nVar,nSpt,SptType,SptIdx,RelTab,ElmMtrx,NonZero,
                       Lvl2CoDim,A,Lvl,ptr,Cell,Cell_Orig,Lvl2LPdim,Lvl2Spt,
                       FixFrstPt,FixLstPt,MinNumPt,PtIn,
		       &Strt1Pt,&End1Pt,ToOrig,Pre2Cur,&LPdim,
                       FrstPt,LstPt,Cur2Pre,x,Binv,Bidx,&info);

               if(info == 0)
                  one_level_LP(Strt1Pt,End1Pt,PtIn[Lvl],
                               LPdim,A[Lvl],nVar,x,Binv,Bidx,ItLp);
            }
         } /* This level finished */
         while(!IT_NextLevel(ItLp) && !IT_IsEmpty(ItLp))
         {
            if(IT_IsFull(ItLp))
            {
               labels[0] = IT_Cell(ItLp)[0];
               labels[1] = IT_Cell(ItLp)[1];
               for(i=2; i<CellSize; i++)
                  labels[i] = ToOrig[i][IT_Cell(ItLp)[i]];
               Cs_Push(MCells,labels);   /* store the cell */
               CellVol(nVar,nSpt,Spt,SptType,labels,MVol);
               ++(*nbCells); 
               Vol = *MVol-tolVol;
               tolVol = *MVol;
            }
            IT_StepBack(ItLp);
         }
      } /* End of while ( !IT_IsEmpty(ItLp) ) */

   } /* End of while ( L0.Migrate(inp) ) */

   for(i=0; i<nVar; i++) free(Binv[i]); free(Binv);
   for(i=0; i<SptIdx[nSpt]; i++) free(RelTab[i]); free(RelTab);

   for(i=0; i<CellSize; i++)
   {
       for(j=0; j<info; j++) free(A[i][j]);
       free(A[i]);
   }
   free(A);

   for(i=0; i<CellSize; i++)
   {
      free(PtIn[i]); free(ElmMtrx[i]);
      free(ToOrig[i]); free(Pre2Cur[i]);
   }
   free(PtIn); free(ElmMtrx); free(ToOrig); free(Pre2Cur);
   free(Bidx); free(x); free(labels);
   free(Cell); free(Cell_Orig); free(Lvl2LPdim); free(Lvl2Spt);
   free(MinNumPt); free(FixFrstPt); free(FixLstPt); free(NonZero);
   free(Lvl2CoDim); free(FrstPt); free(LstPt); free(Cur2Pre);
}

int gcd ( int r1, int r2, int *k, int *l )
{
   int tempxx,tempyy,xx,yy,q,r3;

   *k = 0; *l = 1;
   xx = 1; yy = 0; q = r1/r2; r3 = r1%r2;

   while(r3)
   {
      tempxx = xx - q*(*k);
      tempyy = yy - q*(*l);
      xx =*k;
      yy = *l;
      *k = tempxx;
      *l = tempyy;
      r1 = r2;
      r2 = r3;
      q = r1/r2;
      r3 = r1%r2;
   }
   return ( r2 < 0 )? (*k = -(*k), *l = -(*l), r2 = -r2) : r2;
}

int cell_size ( int nSpt, int *SptType )
{
   int CellSize = nSpt;
   int i;

   for(i=0; i<nSpt; i++) CellSize += SptType[i];

   return CellSize;
}

void CellVol
 ( int nVar, int nSpt, int **Spt, int *SptType, int *Cell, int *Vol )
{
   int i,i1,j,j1,k,l,m,n,d,w,tmp;
   int iwork[nVar][nVar];

   d = -1;
   n = -1;
   for(i=0; i<nSpt; i++)
   {
      l = Cell[++d];
      for(j=0; j<SptType[i]; j++)  
      {
         m = Cell[++d];
         ++n;
         for(k=0; k<nVar; k++) iwork[n][k] = Spt[m][k] - Spt[l][k];
      }
   }
   for(i=0; i<nVar-1; i++)
   {
      j = i;
      while(j<nVar && iwork[j][i] == 0) ++j;

      if(j == nVar)
         printf( "Error code V001\n"); /* volume of a cell = 0 */
      else
      {
         if(j != i)
            for(w=0; w<nVar; w++)
            {
               tmp = iwork[i][w];
               iwork[i][w] = iwork[j][w];
               iwork[j][w] = tmp;
            }
   
         for(j=i+1; j<nVar; ++j)
            if(iwork[j][i] != 0)
            {
               d = gcd(iwork[i][i],iwork[j][i],&k,&l);
               m = -iwork[j][i]/d;
               n = iwork[i][i]/d;
               for(j1=i; j1<nVar; j1++)
               {
                  i1 = iwork[i][j1];  
                  iwork[i][j1] = (k)*i1 + l*iwork[j][j1];
                  iwork[j][j1] = m*i1 + n*iwork[j][j1];
               }
            }
      }
   }
   d = iwork[nVar-1][nVar-1];
   for(i=0; i<nVar-1; i++) d *= iwork[i][i];
      
   if (d < 0)
      *Vol -= d;
   else
      *Vol += d;  
}
         
int solve_linear_system ( int n, double A[n][n], double b[n] )
{
   int i,j,pivot,k;
   double temp,scale,*p,*q,maxvalue,tmp;

   for(i=0; i < n-1; ++i)
   {
      pivot = i;
      maxvalue = fabs(A[i][i]);
 
      for(j = i+1; j < n; ++j)   /* find pivot element */
      {
         temp = fabs(A[j][i]);
         if(temp > maxvalue)
         {
            maxvalue = temp;
            pivot = j;
         }
      }
      if(maxvalue == 0.0) return 0;
  
      if(pivot != i) /* swap two rows */
      {
         for(k=0; k<n; k++)
         { 
            tmp = A[i][k];
	    A[i][k] = A[pivot][k];
	    A[pivot][k] = tmp;
         }
         tmp = b[i];
	 b[i] = b[pivot];
	 b[pivot] = tmp;
      }
      for(j = i+1; j < n; ++j)   /* do elementary transformation */
      {
          scale = -A[j][i]/A[i][i];
          for(p = A[j]+i+1, q=A[i]+i+1; p < A[j] + n; )
            (*p++) += scale*(*q++);
  
          b[j] += scale*b[i];
      }
   }
   if(fabs(A[n-1][n-1]) == 0.0) return 0;
         
   b[n-1] /= A[n-1][n-1];      /* solve upper triangular linear system */
   for(i = n-2; i>=0; --i )
   {
      for(p = b+i+1, q=A[i]+i+1; p < b + n; )
         b[i] -= (*p++)*(*q++);
  
      b[i] /= A[i][i];
   }
   return 1;
}

void write_cells
 ( int ncfn, char output_file[ncfn],
   int nVar, int nSpt, int *SptType, int **Spt, double *lft,
   int CellSize, int nMCells, CellStack *MCells )
{
   FILE *fp;
   int i,k,m,kk,d,dd;
   double A[nVar][nVar];
   double rhs[nVar];
   int *labels = (int*)calloc(CellSize,sizeof(int));

   if((fp=fopen(output_file,"w")) == NULL)
   {
      printf("Opening the file %s for writing failed.\n",output_file);
      exit(1);
   }
   
   fprintf(fp,"%d\n",nVar);   
   fprintf(fp,"%d\n",nSpt);
   
   for(i=0; i<nSpt; i++) fprintf(fp,"%d  ",SptType[i]);
   fprintf(fp,"\n");

   fprintf(fp,"%d\n",nMCells);

   for(kk=0; kk<nMCells; kk++)
   {
      labels = Cs_Cur(MCells);
          
      d = 0; dd = -1;
      for(i=0; i<nSpt; i++)
      {
         for(m=0; m<SptType[i]; m++ )
         {
            dd = dd+1;
            for(k=0; k<nVar; k++)
               A[dd][k] = Spt[labels[d+m+1]][k] - Spt[labels[d+m]][k];

            rhs[dd] = lft[labels[d+m]] - lft[labels[d+m+1]];
         }
         d = d + SptType[i] + 1;
      }
      solve_linear_system(nVar,A,rhs);
      for(k=0; k<nVar; k++) fprintf(fp,"%2.14E\n",rhs[k]);
      fprintf(fp,"%2.14E\n",1.000);
      
      d = -1;
      for(i=0; i<nSpt; i++)
      {
         fprintf(fp,"%d\n",SptType[i]+1);
         for(m=0; m<SptType[i]+1; m++)
         {
            d++;
            for(k=0; k<nVar; k++) fprintf(fp,"%d  ",Spt[labels[d]][k]);
            fprintf(fp,"  %2.14E\n",lft[labels[d]]);
         }
      }
      fprintf(fp,"%d\n",0);
                 
      if(kk != nMCells-1) Cs_Pop(MCells); 
   }     
   fclose(fp);
   free(labels);
}
