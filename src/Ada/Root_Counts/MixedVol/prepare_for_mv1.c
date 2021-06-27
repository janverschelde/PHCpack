/* The file "prepare_for_mv1.c" contains the definitions of the prototypes
 * in the file "prepare_for_mv1.h". */

#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <stdio.h>
#include <math.h>
#include "prepare_for_mv1.h"

#define verbose_prepare 0
 /* 0 : no output; 1 : one line messages; 2 : extra data */

void write_coordinates ( int n, int p[n] )
/*
 * DESCRIPTION :
 *   Writes the coordinates of the point to screen.
 *   Auxilary routine only useful in verbose mode. */
{
   int i;

   for(i=0; i<n; i++) printf(" %d",p[i]); printf("\n");
}

void write_supports ( int nVar, int nSpt, int **Spt, int *SptIdx )
/*
 * DESCRIPTION :
 *   Auxiliary routine to write the supports when in verbose mode.
 *   The parameters are the same as in Pre4MV. */
{
   int i,j;

   for(i=0; i<nSpt; i++)
   {
      printf("Support %d :\n",i);
      for(j=SptIdx[i]; j<SptIdx[i+1]; j++)
         write_coordinates(nVar,Spt[j]);
   }
}

void Pre4MV ( int nVar, int nSpt, int *nS, int *SptType,
              int **Spt, int *SptIdx, int **Vtx, int *VtxIdx, int *OldIdx )
{
   int i,j,k,l,icnt,itmp;
   int ynReOrder,*OrderSpt,*iTmp,*iTmp1,**SptCopy,*NuIdx2OldIdxCopy ;
   int ibreak;   

   /* To sort the points of each support in lexicographicaly decending order */

   if(verbose_prepare>0)
   {
      printf("entering Pre4MV with nVar = %d and nSpt = %d\n",nVar,nSpt);
      printf("indices to supports : ");
      for(i=0; i<=nVar; i++) printf(" %d", SptIdx[i]); printf("\n");
      if(verbose_prepare>1) write_supports(nVar,nSpt,Spt,SptIdx);
   }

   SortSpt(nVar,nSpt,Spt,SptIdx);

   for(i=0; i<nVar; i++)
   {
      for(j=SptIdx[i]; j<SptIdx[i+1]-1; j++)
      {
         itmp = 0;
         for(k=0; k<nVar; k++) itmp += abs(Spt[j][k]-Spt[j+1][k]);
	 if(itmp == 0)
         {
	    printf("%d-th and %d-th points in %d-th support are same!\n",
                   j+1,j+2,i+1);
            exit(0);
	 }
      }
   }

   /* To locate all non-vertex points in each support */
   NonVertex(nVar,nSpt,SptIdx,Spt,OldIdx);

   icnt = -1;
   for(i=0; i<nVar; i++)
   {
      for(j=SptIdx[i]; j<SptIdx[i+1]; j++)
      {
         if(OldIdx[j] == 1)
         {
            icnt = icnt + 1;
            for (k=0; k<nVar; k++) Vtx[icnt][k] = Spt[j][k];
            OldIdx[icnt] = j;
        /*  For the vertex system, OldIdx[i] = j means that */ 
        /*  the current i-th point is the j-th point in the original system */
	 }
      }
      VtxIdx[i+1] = icnt+1;
   }
   VtxIdx[0] = 0;

   /* To reorder the supports if necessary */ 

   ynReOrder = 0;
   OrderSpt = (int*)calloc(nVar,sizeof(int));
   iTmp = (int*)calloc(nVar+1,sizeof(int));
   iTmp1 = (int*)calloc(nVar+1,sizeof(int));
   SptCopy = (int**)calloc(SptIdx[nVar],sizeof(int*));
   for(i=0; i<SptIdx[nVar]; i++)
      SptCopy[i] = (int*)calloc(nVar,sizeof(int));
   NuIdx2OldIdxCopy = (int*)calloc(SptIdx[nVar],sizeof(int));

   for(i=0; i<nVar; i++) iTmp[i] = 0;
   for(i=0; i<nVar; i++)
   {
      itmp = 1000000;
      for(j=0; j<nVar; j++)
      {
         if(iTmp[j] == 0 && VtxIdx[j+1]-VtxIdx[j] < itmp)
	 {
            itmp = VtxIdx[j+1]-VtxIdx[j];
            OrderSpt[i] = j;
	 }
      }
      if(OrderSpt[i] != i) ynReOrder = 1;
      iTmp[OrderSpt[i]] = 1;
   }

   /* The order of the supports : (OrderSpt[i],i=0,nVar-1) */
 
   if(ynReOrder == 1)
   {
      /* To change Spt, Vtx, SptIdx, VtxIdx  */
      iTmp[0] = 0;
      iTmp1[0] = 0;
      icnt = -1;
      for (i=0; i<nVar; i++)
      {
         for(j=SptIdx[OrderSpt[i]]; j<SptIdx[OrderSpt[i]+1]; j++)
	 {
	    icnt = icnt + 1;
            for(k=0; k<nVar; k++) SptCopy[icnt][k] = Spt[j][k];
	 }
         iTmp[i+1] = icnt+1;
      }
      for(j=0; j<=icnt; j++)
         for(i=0; i<nVar; i++) Spt[j][i] = SptCopy[j][i];

      icnt = -1;
      for(i=0; i<nVar; i++)
      {
         for(j=VtxIdx[OrderSpt[i]]; j<VtxIdx[OrderSpt[i]+1]; j++)
	 {
            icnt = icnt + 1;
            for(k=0; k<nVar; k++) SptCopy[icnt][k] = Vtx[j][k];
            NuIdx2OldIdxCopy[icnt]
               = OldIdx[j] - (SptIdx[OrderSpt[i]]-iTmp[i]);
	 }
         iTmp1[i+1] = icnt+1;
      }
      for(j=0; j<=icnt; j++)
      {
         for(i=0; i<nVar; i++) Vtx[j][i] = SptCopy[j][i];
	 OldIdx[j] = NuIdx2OldIdxCopy[j];
      }

      for(i=0; i<=nVar; i++) SptIdx[i] = iTmp[i];
      for(i=0; i<=nVar; i++) VtxIdx[i] = iTmp1[i];
   }
   
   /* To find whether there are some supports same in the vertex system */

   for(i=0; i<nVar; i++) SptType[i] = -1;

   nSpt = 0;
   for(i=0; i<nVar-1; i++)
   {
      if(SptType[i] < 0)
      {
         nSpt = nSpt + 1;
         for(j=i+1; j<nVar; j++)
	 {
            if(SptType[j] < 0)
	    { 
               if(VtxIdx[j+1]-VtxIdx[j] != VtxIdx[i+1]-VtxIdx[i]) continue;
	       ibreak = 0;
               for(k=VtxIdx[i]; k<VtxIdx[i+1]; k++)
	       {
                  for(l=0; l<nVar; l++)
		  {
                     if(Vtx[k][l] != Vtx[k+VtxIdx[j]-VtxIdx[i]][l])
		     {
			ibreak = 1;
		        break;
		     }
		  }
	          if(ibreak) break;
	       }
	       if(ibreak) continue;
               SptType[j] = i;
	    }
         }
      }
   }

   if(SptType[nVar-1] < 0) nSpt = nSpt + 1;

   /*  StpType[i] = -1 : the i-th support is the base support */
   /*          = j  : the i-th support is the same as the j-th base support */

   if(nSpt<nVar)
   {
      /* For the unmixed case, rearrange the supports according to the type */
      ReArrangeSpt(nVar,Spt,SptIdx,Vtx,VtxIdx,SptType,OldIdx);
   }
   else
   {
      nSpt = nVar;
      for(i=0; i<nVar; i++) SptType[i] = 1;
   }

   *nS=nSpt;

   if(verbose_prepare > 0)
   {
      printf("leaving Pre4MV with nS = %d\n",*nS);
      printf("Type of mixture : ");
      for(i=0; i<*nS; i++) printf(" %d",SptType[i]); printf("\n");
      printf("indices to supports : ");
      for(i=0; i<=*nS; i++) printf(" %d", VtxIdx[i]); printf("\n");
      write_supports(nVar,*nS,Vtx,VtxIdx);
      printf("The old indices : ");
      for(i=0; i<SptIdx[nVar]; i++) printf(" %d",OldIdx[i]); printf("\n");
   }

   free(OrderSpt); free(iTmp); free(iTmp1);
   for(i=0; i<SptIdx[nVar]; i++) free(SptCopy[i]);
   free(SptCopy);
   free(NuIdx2OldIdxCopy);
}

void SortSpt ( int nVar, int nSpt, int **Spt, int *SptIdx )
{   
   int iSpt,j,k,idx,i;
   
   for(iSpt=0; iSpt<nSpt; iSpt++)
   {
      for(i=SptIdx[iSpt]; i<SptIdx[iSpt+1]-1; i++)
      {
         idx = i;
         for(j=i+1; j<SptIdx[iSpt+1]; j++)
         {
            for(k=0; k<nVar; k++)
            {
               if(Spt[j][k] < Spt[idx][k])
               {
                  break;    /* LexiGt = 0 */
               }
               else
               {   
                  if(Spt[j][k] > Spt[idx][k])
                  {
                     idx = j;
                     break; /* LexiGt = 1; */
                  }
               }   
            }      
         }         
         if(idx>i)
         {
            for(j=0; j<nVar; j++)
            {
               k = Spt[i][j];
               Spt[i][j] = Spt[idx][j];
               Spt[idx][j] = k;
            }
         }
      }
   }
}

void LowerTriangular
 ( double **A, int jStrt, int jEnd, int na, int *rnk, int* ib )
{   
   int i, j, k, itmp;
   double dtmp;
   
   ib[0] = jStrt;                 /* the jStrt-th row of A is (1,0,...,0) */
   for(i=1; i<na; i++) ib[i] = -1;
   *rnk = 1;
   
   k = jStrt+1;
   while(*rnk < na && k <= jEnd)
   {
      dtmp = 1.0e-13;                      /* search for largest component */
      itmp = -1;                          /* from A[k][*rnk] to A[k][na-1] */
      for(i = *rnk; i<na; i++)
         if(fabs(A[k][i]) > dtmp)
         {
            dtmp = fabs(A[k][i]);
            itmp = i;
         }
      if(itmp >= 0)                                         /* nonzero row */
      {                                   /* To change aa(k,*) into e_itmp */
         for(i=0; i<itmp; i++) A[k][i] /= A[k][itmp];
         for(i=itmp+1; i<na; i++) A[k][i] /= A[k][itmp];
         for(j=k+1; j<=jEnd; j++)
         {
            for(i=0; i<itmp; i++) A[j][i] -= A[k][i]*A[j][itmp];
            for(i=itmp+1; i<na; i++) A[j][i] -= A[k][i]*A[j][itmp];
            A[j][itmp] /= A[k][itmp];
         }
         if(itmp != *rnk)         /* interchange columns rnk and itmp of A */
         {
            for(j=jStrt; j<=jEnd; j++)
            {  
               dtmp = A[j][*rnk];
               A[j][*rnk] = A[j][itmp];
               A[j][itmp] = dtmp;
            }
         }
         
         for(i=0; i<*rnk; i++) A[k][i] = 0.0e0;
         A[k][*rnk] = 1.0e0;
         for(i=*rnk+1; i<na; i++) A[k][i] = 0.0e0;
            
         ++(*rnk);
         ib[*rnk-1]=k;
      }   
      ++k;  
   }         
}              

void RSimplex
 ( double **A, double *b, int jStrt, int jEnd, int m,
   double *xb, double *sb, int *ib, double **binv, int *info )
{
   int i,j,k,ell;
   int ibrk;
   double sigb, sum, vkmin;
   double eps = 1.0e-6;
   double tol = 1.0e-10;
                
   *info = -1;  
   for(i=0; i<m; i++)
   {
      if(ib[i] == jStrt)
      {
         *info = i;
         break;
      }
   }
   while(1)
   {                                              /* compute u and find k */
      if(*info == -1)                    /* $e_{n+1}$ is not in the basis */
      {
         k = jStrt;
         vkmin = -1.0e0;
      }
      else                                   /* $e_{n+1}$ is in the basis */
      {                   /* find vkmin=v_k=min{c_i+a_i^Tu | i\notin I_b} */
         k = -1;
         vkmin = DBL_MAX;
         for(i=jStrt; i<=jEnd; i++)
         {                                         /* check if i is in ib */
            ibrk = 0;
            for(j=0; j<m; j++)
            {
               if(i == ib[j])
               {
                  ibrk = 1;
                  break;
               }
            }
            if(ibrk) continue;
                                                    /* if i is not in ib */
            sum = 0.0e0;
            for(j=0; j<m; j++) sum += A[i][j]*binv[*info][j];
            if(sum < vkmin)
            {
               vkmin = sum;
               k = i;
            }
         }
         if(vkmin > -eps) return;          /* found the optimal solution  */
      }              /* k-th column will get into the basis -- to form sb */
      for(i=0; i<m; i++)
      {
         sum = 0.0e0;
         for(j=0; j<m; j++) sum += binv[i][j]*A[k][j];
         sb[i] = sum;
      }
      ell = 0;
      sigb = DBL_MAX;
      for(i=0; i<m; i++)
      {
         if(sb[i] > eps && (xb[i]/sb[i]) < sigb)
         {
            sigb = (xb[i]/sb[i]);
            ell = i;
         }
      } /* ell-th row gets out of the basis and k-th row gets into the basis.
                                                    to find the new B^{-1} */
      sum = 0.0e0;
      for(i=0; i<m; i++) sum += binv[ell][i]*A[k][i];
            
      if(fabs(sum) <= eps) abort();
          
      sum = 1.0e0/sum;
      for(i=0; i<m; i++) binv[ell][i] *= sum;
         
      for(j=0; j<m; j++)
      {
         if(j != ell)   
         {
            sum = 0.0e0;
            for(i=0; i<m; i++) sum += binv[j][i]*A[k][i];
            for(i=0; i<m; i++) binv[j][i] -= sum*binv[ell][i];
         }
      }
      for(i=0; i<m; i++)           /* form the new basic feasible solution */
      {
         sum = 0.0e0;
         for(j=0; j<m; j++) sum += binv[i][j]*b[j];
         xb[i] = sum;
      }   
      if(ib[ell] == jStrt) *info = 0;      /* column jStrt is out of basis */
      ib[ell] = k;
      if(k == jStrt)
      {     
         *info = ell;
         if(xb[*info] > tol) return;        /* epsilon is already positive */
      }
   }   
}         

void ReArrangeSpt
 ( int nVar, int **Spt, int *SptIdx, int **Vtx, int *VtxIdx,
   int *SptType, int *OldIdx )
{ 
   int i,i1,icnt,iSpt,ith,j;
   int *iTmp,**SptCopy,*SptIdxCopy,*Nu2OldCopy;
   
   iTmp = (int*)calloc(nVar,sizeof(int));
   
   i = 0;
   iSpt = -1;
   while(i<nVar)
   {
      if(SptType[i] < 0)
      {
         ++iSpt;
         iTmp[iSpt] = i;
      }
      ++i;
   }
   i = 0;
   ith = -1;
   while(i<nVar)
   {
      if(SptType[i] < 0)
      {
         ++ith;
         icnt = 1;
         for(j=i+1; j<nVar; j++)
         {
            if(SptType[j] == i)
            {
               ++iSpt;
               ++icnt;
               iTmp[iSpt] = j;
            }
         }
         SptType[ith] = icnt;
      }
      i = i + 1;
   }
   /* SptType[i] = k : #supports which are the same as the i-th base support.
      To rearrange Spt, SptIdx */
    
   SptCopy = (int**)calloc(SptIdx[nVar],sizeof(int*));
   for(i=0; i<SptIdx[nVar]; i++)
      SptCopy[i] = (int*)calloc(nVar,sizeof(int));
   SptIdxCopy = (int*)calloc(nVar+1,sizeof(int));
       
   for(j=0; j<SptIdx[nVar]; j++)
      for(i=0; i<nVar; i++) SptCopy[j][i] = Spt[j][i];
   
   icnt = -1;
   for(i=0; i<nVar; i++)
   {   
      for(j=SptIdx[iTmp[i]]; j<SptIdx[iTmp[i]+1]; j++)
      {
         ++icnt;
         for (i1=0; i1<nVar; i1++) Spt[icnt][i1] = SptCopy[j][i1];
      }
      SptIdxCopy[i+1] = icnt+1;
   }
   SptIdxCopy[0] = 0;
               
   /* To rearrange OldIdx */
          
   Nu2OldCopy = (int*)calloc(SptIdx[nVar],sizeof(int));
       
   for(j=0; j<VtxIdx[nVar]; j++) Nu2OldCopy[j] = OldIdx[j];
    
   icnt = -1;
   for(i=0; i<nVar; i++)
      for(j=VtxIdx[iTmp[i]]; j<VtxIdx[iTmp[i]+1]; j++)
      {
         ++icnt;
         OldIdx[icnt]
           = Nu2OldCopy[j] - (SptIdx[iTmp[i]] - SptIdxCopy[i]);
      }
       
   for(i=0; i<=nVar; i++) SptIdx[i] = SptIdxCopy[i];
    
   /* To rearrange SptVtx, VtxIdx */
    
   for(j=0; j<VtxIdx[nVar]; j++)
      for(i=0; i<nVar; i++) SptCopy[j][i] = Vtx[j][i];
   
   icnt = -1;
   for(i=0; i<nVar; i++)
   {   
      for(j=VtxIdx[iTmp[i]]; j<VtxIdx[iTmp[i]+1]; j++)
      {
         ++icnt;
         for(i1=0; i1<nVar; i1++) Vtx[icnt][i1] = SptCopy[j][i1];
      }
      SptIdxCopy[i+1] = icnt+1;
   }           
   VtxIdx[0] = 0;
   for(i=1; i<=nVar; i++) VtxIdx[i] = SptIdxCopy[i];
   
   free(iTmp); free(Nu2OldCopy); free(SptIdxCopy);
   for(i=0; i<SptIdx[nVar]; i++) free(SptCopy[i]);
   free(SptCopy);
}

void NonVertex
 ( int nVar, int nSpt, int *SptIdx, int **Spt, int *ynVtx )
{
   int i,i1,info,itmp,j,jStrt,jEnd,k,PtIn,rnk;
   double tmp,**A,*b,*xb,*sb,**binv;
   int *ib,*ib0;
   
   sb = (double*)calloc(nVar+2,sizeof(double));
   A = (double**)calloc(SptIdx[nSpt]+1,sizeof(double*));
   for(i=0; i<=SptIdx[nSpt]; i++)
      A[i] = (double*)calloc(nVar+2,sizeof(double));
   b = (double*)calloc(SptIdx[nSpt]+1,sizeof(double));
   ib = (int*)calloc(nVar+2,sizeof(int));
   ib0 = (int*)calloc(nVar+2,sizeof(int));
   xb = (double*)calloc(nVar+2,sizeof(double));
   binv = (double**)calloc(nVar+2,sizeof(double*));
   for(i=0; i<nVar+2; i++)
      binv[i] = (double*)calloc(nVar+2,sizeof(double));
      
   for(i=0; i<SptIdx[nSpt]; i++) ynVtx[i] = 1;
   /* = 0: nonvertex point; to form the matrix for simplex method
           to choose a random point */
    
   srand((unsigned)(time(0)));
   for(i=0; i<nVar; i++) b[i] = (double)(rand())/RAND_MAX;
   
   for(j=0; j<SptIdx[nSpt]; j++)  /* lift the points */
   {
      tmp = 0.0e0;
      for(i1=0; i1<nVar; i1++) tmp=tmp+(Spt[j][i1]-b[i1])*(Spt[j][i1]-b[i1]);
      A[j+1][0] = -sqrt(tmp);  /* j+1: shift a slot for the added 1st column */
      for(i1=0; i1<nVar; i1++) A[j+1][i1+1] = Spt[j][i1];
      A[j+1][nVar+1] = 1.0e0;
   }   
   for(i=0; i<nSpt; i++)                                /* for each support */
   {   
      jStrt = SptIdx[i];          /* The first row is the lifting direction */
      jEnd = SptIdx[i+1];      /* +1: shift a slot for the added 1st column */
         
      A[jStrt][0] = 1.0e0;
      for(i1=1; i1<nVar+2; i1++) A[jStrt][i1] = 0.0e0;
             
      rnk = 0;                          /* To lower-trianglize the matrix A */

      LowerTriangular(A,jStrt,jEnd,nVar+2,&rnk,ib0);
             
      for(j=0; j<rnk; j++)
      {
         for (i1=0; i1<rnk; i1++) binv[j][i1] = 0.0e0;
         binv[j][j] = 1.0e0;
      }
            
      for(j=jStrt+1; j<=jEnd; j++) /* for each point */
      {
         PtIn = j;
         for(i1=0; i1<rnk; i1++) ib[i1] = -1;
         for(i1=0; i1<rnk; i1++)
            if(ib0[i1] == j)
            {
               ib[i1] = PtIn;
               PtIn = -PtIn;
               break;
            }   
         if(PtIn >= 0)
         {
            /* The j-th row will get into the basis. */
            /* To decide which column will get out of the basis */
            for(k=0; k<rnk; k++)
            {
               xb[k] = 0.0e0;
               for(i1=0; i1<rnk; i1++) xb[k] += binv[k][i1]*A[j][i1];
            }
            
            /* To find the biggest component from xb[0] to xb[*rnk-1] */
            tmp = 0.0e0;
            itmp = -1;
            for(k=1; k<rnk; k++)
               if(fabs(xb[k]) > tmp && ib[k] == -1)
               {
                  tmp = fabs(xb[k]);
                  itmp = k;
               }
            if(itmp == -1)
            {
               info = -1;
               abort();  /* It is impossible in our problem */
            }
            
            /* itmp-th row will get out of the basis */
            /* To find the new B^{-1} after itmp out and j in */
               
            for(k=0; k<itmp; k++) xb[k] /= xb[itmp];
            for(k=itmp+1; k<rnk; k++) xb[k] /= xb[itmp];
            for(i1=0; i1<rnk; i1++)
            {
               tmp = binv[itmp][i1];
               for (k=0; k<itmp; k++) binv[k][i1] -= xb[k]*tmp;   
               for (k=itmp+1; k<rnk; k++) binv[k][i1] -= xb[k]*tmp;
               binv[itmp][i1] = tmp/xb[itmp];
            }
            ib[itmp] = j;
         }   
            
         for(i1=0; i1<rnk; i1++)
         {
            if(ib[i1] != -1)
            {
               xb[i1] = 1.0e0;
            }
            else
            {
               xb[i1] = 0.0e0;
               ib[i1] = ib0[i1];
            }
         }
         for(i1=0; i1<rnk; i1++) b[i1] = A[j][i1];
               
         RSimplex(A,b,jStrt,jEnd,rnk,xb,sb,ib,binv,&info);
         if(info > -1)
         {     
            if(fabs(xb[info]) > 1.0e-10)
            {  
                                /* It is a non-vertex */
               ynVtx[j-1] = 0;  /* -1: since the added 1st column */
            }
         }   
         for(i1=0; i1<rnk; i1++) ib0[i1] = ib[i1];
      }
   }
               
   for(i=0; i<=SptIdx[nSpt]; i++) free(A[i]);
   free(A);
   free(b); free(ib); free(ib0); free(xb); free(sb);
   for(i=0; i<nVar+2; i++) free(binv[i]);
   free(binv); 
}
