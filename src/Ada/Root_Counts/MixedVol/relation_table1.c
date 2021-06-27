#include <string.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "relation_table1.h"

#define verbose_reltab 0

void RelTable ( int nVar, int nSpt, int **Spt, int *SptIdx,
                double *lft, int **RelTab, L0_IML *L0 )
{
   const int nv1 = nVar+1;
   int i,j,j1,info,ell,ynNuPt;
   int CurSpt,NxSpt,FixPt,TstPt,Strt1PtTst,nPtIn;
   int *Bidx_s,*LPidx,*Bidx;
   double dtmp;
   double *x_s,**a,*x,**Binv;

   double Binv_s[nVar][nVar];
   double *c = (double*)calloc(nVar+2,sizeof(double));
   int *J = (int*)calloc(2,sizeof(int));

   if(verbose_reltab > 0) printf("entering RelTable...\n");

   x_s = (double*)calloc(nVar,sizeof(double));
   Bidx_s = (int*)calloc(nVar,sizeof(int));

   a = (double**)calloc(SptIdx[nSpt]+1,sizeof(double*));
   for(i=0; i<=SptIdx[nSpt]; i++)
      a[i] = (double*)calloc(nVar+2,sizeof(double));

   LPidx = (int*)calloc(SptIdx[nSpt]+1,sizeof(int));

   x = (double*)calloc(nVar+2,sizeof(double));
   Bidx = (int*)calloc(nVar+2,sizeof(int));
   Binv = (double**)calloc(nVar+2,sizeof(double*));
   for(i=0; i<nVar+2; i++)
      Binv[i] = (double*)calloc(nVar+2,sizeof(double));

   for(i=0; i<SptIdx[nSpt]; i++)
      for(j=0; j<SptIdx[nSpt]; j++) RelTab[i][j] = 0;
 
   for(j=SptIdx[1]; j<SptIdx[nSpt]; j++)           /* [0] changing */
   {
      for(i=0; i<nVar; i++) a[j][i] = -Spt[j][i];
      a[j][nVar] = 1.0;                            /* WRT alpha0 */
      a[j][nVar+1] = lft[j];
   }
   
   for(CurSpt=0; CurSpt<nSpt; CurSpt++)
   {
      for(FixPt=SptIdx[CurSpt]; FixPt<SptIdx[CurSpt+1]; FixPt++)
      {
         nPtIn = -1;                /* ! # of constraints involved - 1 */
         for(j=SptIdx[CurSpt]; j<FixPt; j++)
            if(RelTab[FixPt][j])
            {
               LPidx[++nPtIn] = j;
               for(i=0; i<nVar; i++) a[j][i] = Spt[FixPt][i] - Spt[j][i];
               a[j][nVar] = 1.0;   /* for alpha0 */
               a[j][nv1] = lft[j] - lft[FixPt];  /* the constant term */
            }
   
         Strt1PtTst = nPtIn + 1;    /* ! starting index of pts for 1-pt test */
         for(j=FixPt+1; j<SptIdx[CurSpt+1]; j++)
         {
            LPidx[++nPtIn] = j;
            for(i=0; i<nVar; i++) a[j][i] = Spt[FixPt][i] - Spt[j][i];
            a[j][nVar] = 1.0;
            a[j][nv1] = lft[j] - lft[FixPt];
         }
         /* To extend the point by using 1-point test
            To find a nondegenerate solution for 1-point tests of the points */
    
         RlTbLP2_e(nPtIn+1,nVar,SptIdx[nSpt],a,nv1,LPidx,Bidx,x,Binv,&info);
      
         if(FixPt < SptIdx[CurSpt+1]-1)
         {
            /* To record any possible points in Bidx passed 1-point test */
            if(CurSpt == 0)
            {
               ynNuPt = 0;
               for(i=0; i<nVar; i++)
                  if(Bidx[i] > FixPt && !RelTab[FixPt][Bidx[i]])  
                  {
                     ynNuPt = 1;
                     RelTab[FixPt][Bidx[i]] = RelTab[Bidx[i]][FixPt] = 1;
                  }
               if(ynNuPt)
               {
                  J[0] = FixPt;
                  for(i=0; i<nVar; i++)
                     if(Bidx[i] > FixPt)
                     {
                        J[1] = Bidx[i];
                        L0_Add2(L0,J,nVar,Bidx,x,Binv);  /* add 2 points */
                     }
               }
            }
            else
               for(i=0; i<nVar; i++)
                  if(Bidx[i] > FixPt)
                  {
                     RelTab[Bidx[i]][FixPt] = RelTab[FixPt][Bidx[i]] = 1;
                     for(j=i+1; j<nVar; j++)
                        if(Bidx[j] > FixPt)
                           RelTab[Bidx[i]][Bidx[j]] =
                              RelTab[Bidx[j]][Bidx[i]] = 1;
                  }
            /* x, Bidx, Binv are available for 1-point test on other points */

            if(info >= 0)
               for(j=Strt1PtTst; j<=nPtIn; j++)
               {
                  TstPt = LPidx[j];
                  if(!RelTab[FixPt][TstPt])
                  {
                     for(i=0; i<nVar; i++) c[i] = -a[TstPt][i];
                     if(CurSpt == 0)
                        dlp1_1pt_s(nPtIn+1,nVar,a,nv1,c,LPidx,FixPt,TstPt,
                                   Bidx,x,Binv,RelTab,L0);
                     else
                        dlp1_1pt_i(nPtIn+1,nVar,a,nv1,c,LPidx,FixPt,TstPt,
                                   Bidx,x,Binv,RelTab);
                  }   
               }
            else
               for(j=Strt1PtTst; j<=nPtIn; j++)
               {
                  TstPt = LPidx[j];
                  if(!RelTab[FixPt][TstPt])
                  {
                     for(i=0; i<nVar; i++) c[i] = -a[TstPt][i];
                     if(CurSpt == 0)
                        dlp2_1pt_s(nPtIn+1,nVar,a,nv1,c,LPidx,FixPt,TstPt, 
                                   Bidx,x,Binv,RelTab,L0);
                     else
                        dlp2_1pt_i(nPtIn+1,nVar,a,nv1,c,LPidx,FixPt,TstPt,
                                   Bidx,x,Binv,RelTab);  
                  }
               }
         }
         for(i=0; i<nVar; i++)   /* save the starting point for next LP */
         {
            for(j=0; j<nVar; j++) Binv_s[i][j] = Binv[i][j];
            Bidx_s[i] = Bidx[i];
            x_s[i] = x[i];
         }
                   
         /* To check point FixPt in S_CurSpt and each point in S_i (i>CurSpt)*/
                     
         nPtIn = -1;  /* To get all constraints of pts related to FixPt */ 
         for(i=SptIdx[CurSpt]; i<SptIdx[CurSpt+1]; i++)  
            if(RelTab[FixPt][i])
            {
               LPidx[++nPtIn] = i;
               a[i][nVar] = 0.0;    /* no alpha0 */
            }   
         Strt1PtTst = nPtIn + 1;    /* starting index of pts for 1-pt test */
                
         for(NxSpt=CurSpt+1; NxSpt<nSpt; NxSpt++)
         {
            nPtIn = Strt1PtTst - 1;
            for(i=SptIdx[NxSpt]; i<SptIdx[NxSpt+1]; i++) LPidx[++nPtIn] = i;
            /*
              The part of Ax<=B for this support was formed at the beginning
              to do the 1-point test for the points in S_(NxSpt)
              to form the starting point from saved starting Bidx, x, Binv
            */
            for(i=0; i<nVar; i++)
            {
               Bidx[i] = Bidx_s[i];
               x[i] = x_s[i];
            }
            ell = -1;
            x[nVar] = DBL_MAX;  
            for(j1=Strt1PtTst; j1<=nPtIn; j1++)
            {
               j = LPidx[j1];
               dtmp = a[j][nv1];
               for(i=0; i<nVar; i++) dtmp -= a[j][i]*x_s[i];
               if(dtmp < x[nVar])
               {
                  x[nVar] = dtmp;  
                  ell = j;
               }
            }
            Bidx[nVar] = ell;
            for(i=0; i<nVar; i++)
            {
               for(j=0; j<nVar; j++) Binv[i][j] = Binv_s[i][j];
               Binv[i][nVar] = 0.0;
            }
            for(i=0; i<nVar; i++) Binv[nVar][i] = 0.0;
            for(j=0; j<nVar; j++)
            {
               Binv[j][nVar] = 0.0;
               for(i=0; i<nVar; i++) Binv[j][nVar] -= Binv[j][i]*a[ell][i];
            }        
            Binv[nVar][nVar] = 1.0;
            
            /* find a nondegenerate solution for 1-point tests of the points */
            
            TstPt = LPidx[Strt1PtTst];
            for(i=0; i<=nVar; i++) c[i] = -a[TstPt][i];
               
            RlTbLP2_a(nPtIn+1,nVar+1,a,nv1,c,LPidx,Bidx,x,Binv,&info);
            
            /* record any possible points in Bidx passed 1-point test */
            
            for(i=0; i<=nVar; i++)
            {
               if(Bidx[i] >= SptIdx[NxSpt])
               {
                  if(Bidx[i] >= TstPt)
                     RelTab[Bidx[i]][FixPt] = RelTab[FixPt][Bidx[i]] = 1;
                  for(j=i+1; j<=nVar; j++)
                     if(Bidx[j] >= SptIdx[NxSpt])
                        RelTab[Bidx[i]][Bidx[j]] = RelTab[Bidx[j]][Bidx[i]] = 1;
               }
            }
            /* do the 1-point test for other points */
               
            if(info >= 0)
               for(j=Strt1PtTst+1; j<=nPtIn; j++)
               {
                  TstPt = LPidx[j];
                  if(!RelTab[FixPt][TstPt])
                  {
                     for(i=0; i<=nVar; i++) c[i] = -a[TstPt][i];
                     dlp1_1pt_i(nPtIn+1,nVar+1,a,nv1,c,LPidx,FixPt,TstPt,
                                Bidx,x,Binv,RelTab);
                  }
               }
            else
               for(j=Strt1PtTst+1; j<=nPtIn; j++)
               {
                  TstPt = LPidx[j];
                  if(!RelTab[FixPt][TstPt])
                  {
                     for(i=0; i<=nVar; i++) c[i] = -a[TstPt][i];
                     dlp2_1pt_i(nPtIn+1,nVar+1,a,nv1,c,LPidx,FixPt,TstPt,
                                Bidx,x,Binv,RelTab);
                  }
               }
         }
      }
   }
   for(i=0; i<=SptIdx[nSpt]; i++) free(a[i]);
   free(a);
   free(J); free(x); free(Bidx); free(x_s); free(Bidx_s);
   for(i=0; i<nVar+2; i++) free(Binv[i]);
   free(Binv);
   free(LPidx); free(c);
                     
   if(verbose_reltab > 0)
   {
      printf("leaving RelTable, the relation table :\n");
      for(i=0; i<SptIdx[nSpt]; i++)
      {
         for(j=0; j<SptIdx[nSpt]; j++) printf(" %d",RelTab[i][j]);
         printf("\n");
      }
   }
}

void RlTbLP2_a 
 ( int ma, int na, double **a, int nv1, double *c,
   int *LPidx, int *Bidx, double *x, double **Binv, int *info )
{
   int ibrk,ell,k,k1,i,i1,j,flag;
   double eps=1.0e-6,sigj,vmax,bot,top,dtmp;
   double v[na];

   if(verbose_reltab > 0)
      printf("entered RlTbLP2_a with ma = %d and na = %d\n",ma,na);

   while(1)
   { 
      ibrk = 0;             /* step 1 : search for outgoing constraint */
      for(i=0; i<na; i++)
         if(Bidx[i] == -1)
         {
            ibrk = 1;
            break;
         }
      if(ibrk)
      {
         vmax = -1.0;
         for(i=0; i<na; i++)
            if(Bidx[i] == -1)
            {
               dtmp = c[0] * Binv[i][0];
               for(j=1; j<na; j++) dtmp += c[j] * Binv[i][j];
               if(fabs(dtmp) > vmax)
	       {
                  bot = dtmp;
                  vmax = fabs(bot);
                  k = i;
               }
            }
         if(vmax > eps)
            if(bot >= 0.0)
               for(i=0; i<na; i++) v[i] = Binv[k][i];
            else
               for(i=0; i<na; i++) v[i] = -Binv[k][i];
         else
            ibrk = 0;
      }
      if(!ibrk)
      {
         vmax = -DBL_MAX;
         for(i=0; i<na; i++)
            if(Bidx[i] != -1)
            {
               dtmp = c[0] * Binv[i][0];
               for(j=1; j<na; j++) dtmp += c[j] * Binv[i][j];
               if(dtmp > vmax)
               {
                  vmax = dtmp;
                  k = i;
               }
            }
         if(vmax < eps) break;  /* found optimal solution, solve smaller LP */

         for(i=0; i<na; i++) v[i] = Binv[k][i];
      }
      sigj = DBL_MAX;           /* step 2 : search for incoming contraint */
      ell = -1;
      for(i1=0; i1<ma; i1++)
      {
         i = LPidx[i1];
         ibrk = 0;
         for(j=0; j<na; j++)
            if(Bidx[j] == i)
            {
               ibrk = 1;
               break;
            }
         if(ibrk) continue;

         bot = a[i][0] * v[0];
         for(j=1; j<na; j++) bot += a[i][j] * v[j];
         if(bot >= -eps) continue;

         top = -a[i][nv1];
         for(j=0; j<na; j++) top += a[i][j] * x[j];
         top /= bot;
         if(top >= sigj) continue;

         sigj = top;
         ell = i;
      }
      if(ell < 0)
      {
         printf("RlTbLP2_a : LP unbounded\n\a");
         abort();
      }                                  /* step 3 : update x, Bidx, Binv */
      for(i=0; i<na; i++) x[i] -= sigj * v[i];
      top = a[ell][0] * Binv[k][0];
      for(i=1; i<na; i++) top += a[ell][i] * Binv[k][i];
      if(fabs(top) < eps)
      {
         printf("RlTbLP2_a : Base matrix is singular\a\n");
         abort();
      }
      top = 1.0/top;
      for(i=0; i<na; i++) Binv[k][i] *= top;
      for(j=0; j<na; j++)
         if(j != k)
         {
            top = a[ell][0] * Binv[j][0];
            for(i=1; i<na; i++) top += a[ell][i] * Binv[j][i];
            for(i=0; i<na; i++) Binv[j][i] -= top * Binv[k][i];
         }

      Bidx[k] = ell;
   }
   k1 = -1;
   *info = 0;
   for(i=0; i<na; i++)
      if(Bidx[i] > -1) ++(*info);

   if(*info == na) return;   /* case of full rank */

   k = 0;
   while(k < na)
   {
      if(Bidx[k] > -1 || k == k1 )
      {
         ++k;
         continue;
      }

      for(i=0; i<na; i++) v[i] = Binv[k][i];

      flag = 1;
      ell = -1;               /* step 2 : search incoming constraint */
      while(ell == -1)
      {
         sigj = DBL_MAX;
         for(i1=0; i1<ma; i1++)
         {
            i = LPidx[i1];
            ibrk = 0;
            for(j=0; j<na; j++)
               if(Bidx[j] == i)
               {
                  ibrk = 1;
                  break;
               }
            if(ibrk) continue;

            bot = a[i][0] * v[0];
            for(j=1; j<na; j++) bot += a[i][j] * v[j];
            if(bot >= -eps) continue;

            top = -a[i][nv1];
            for(j=0; j<na; j++) top += a[i][j] * x[j];
            top /= bot;
            if(top >= sigj) continue;

            sigj = top;
            ell = i;
         }

         ibrk = 0;
         if(ell < 0)
            if(flag == 1)
            {
               for(i=0; i<na; i++) v[i] = -Binv[k][i];
               flag = 0;
            }
            else
            {
               k = k + 1;
               ibrk = 1;
               break;
            }
      }
 
      if(ibrk) continue;
                                        /* step 3 : update x, Bidx, Binv */
      for(i=0; i<na; i++) x[i] -= sigj * v[i];

      top = a[ell][0] * Binv[k][0];
      for(i=1; i<na; i++) top += a[ell][i] * Binv[k][i];
      if(fabs(top) < eps)
      {
         printf("RlTbLP2_a : Base matrix is singular\a\n");
         abort();
      }
      top = 1.0/top;
      for(i=0; i<na; i++) Binv[k][i] *= top;
      for(j=0; j<na; j++)
         if(j != k)
         {
            top = a[ell][0] * Binv[j][0];
            for(i=1; i<na; i++) top += a[ell][i] * Binv[j][i];
            for(i=0; i<na; i++) Binv[j][i] -= top * Binv[k][i];
         }
      Bidx[k] = ell;
      *info = *info + 1;

      if(*info == na) 
         return;
      else
      { 
         k1 = k;
         k = 0;
      }   
   }
   *info = - *info;   /* the case of not full rank */
}

void RlTbLP2_e
 ( int ma, int na, int NumCol, double **a, int nv1, int *LPidx,
   int *Bidx, double *x, double **Binv, int *info )
{
   int ibrk,ell,k,k1,i,i1,j,flag;
   double eps=1.0e-6,sigj,vmax,bot,top,dtmp;
   const int nap1 = na + 1;
   const int map1 = ma + 1;
   double v[nap1];

   if(verbose_reltab > 0)
      printf("entered RlTbLP2_e with ma = %d and na = %d\n",ma,na);

/* Expand the system (add variable epsilon)
   Assumed a, b, Binv, x have extra 1 dimension to hold epsilon */

   for(i=0; i<na; i++) a[NumCol][i] = 0.0;  /* -epsilon<=0 */
   a[NumCol][na] = -1.0;
   a[NumCol][nv1] = 0.0;
   LPidx[ma] = NumCol;
   dtmp = 0.0;
   ell = NumCol;             /* the index of the equality equation */
   for(i=0; i<ma; i++)
   {
      j = LPidx[i];
      if(a[j][nv1] < -eps)
      {
         a[j][nap1-1] = -1;
         if(-a[j][nv1] > dtmp)
         {
            dtmp = -a[j][nv1];
            ell = j;
         }
      }
      else
         a[j][nap1-1] = 0.0;
   } 
   for(i=0; i<na; i++) Bidx[i] = -1;  /* form the first x, Bidx, Binv */
   Bidx[nap1-1] = ell;
   for(i=0; i<na; i++) x[i] = 0.0;
   x[na] = dtmp;
   for(i=0; i<nap1; i++)
   {
      for(j=0; j<nap1; j++) Binv[i][j] = 0.0;
      Binv[i][i] = 1.0;
   }
   for(j=0; j<na; j++) Binv[j][nap1-1] = a[ell][j];
   Binv[nap1-1][nap1-1] = -1.0;

    /* apply the LP algorithm to the larger LP
       step 1: To find which constraint to kick out */

   while(1)
   {
      ibrk = 0;
      for(i=0; i<nap1; i++)
         if(Bidx[i] == -1)
	 {
            ibrk = 1;
            break;
         }
      if(ibrk)
      {
         vmax = -1.0;
         for(i=0; i<nap1; i++)
            if(Bidx[i] == -1) 
            {
               dtmp = Binv[i][nap1-1];
               if(fabs(dtmp) > vmax)
               {
                  bot = dtmp;
                  vmax = fabs(bot);
                  k = i;
               }
            }
         if(vmax > eps)
            if(bot >= 0.0)
               for(i=0; i<nap1; i++) v[i] = Binv[k][i];
            else
               for(i=0; i<nap1; i++) v[i] = -Binv[k][i];
         else
            ibrk = 0;
      }
      if(!ibrk)
      {
         vmax = -DBL_MAX;
         for(i=0; i<nap1; i++)
            if(Bidx[i] != -1) 
            {
               dtmp = Binv[i][nap1-1];
               if(dtmp > vmax)
               {
                  vmax = dtmp;
                  k = i;
               }
            }
         if(vmax < eps) break;  /* found optimal solution, solve smaller LP */

         for(i=0; i<nap1; i++) v[i] = Binv[k][i];
      }
      sigj = DBL_MAX;           /* search for incoming constraint */
      ell = -1;
      for(i1=0; i1<map1; i1++)
      {
         i = LPidx[i1];
         ibrk = 0;
         for(j=0; j<nap1; j++)
            if(Bidx[j] == i)
            {
               ibrk = 1;
               break;
            }
         if(ibrk) continue;

         bot = a[i][0] * v[0];
         for(j=1; j<nap1; j++) bot += a[i][j] * v[j];
         if(bot >= -eps)  continue;

         top = -a[i][nv1];
         for(j=0; j<nap1; j++) top += a[i][j] * x[j];
         top /= bot;
         if(top >= sigj) continue;

         sigj = top;
         ell = i;
      }
      if(ell < 0)
      {
         printf("RlTbLP2_e : LP unbounded??\n\a");
         abort();
      }                                      /* step 3 : update x, Bidx, Binv */
      for(i=0; i<nap1; i++) x[i] -= sigj * v[i];
      top = a[ell][0] * Binv[k][0];
      for(i=1; i<nap1; i++) top += a[ell][i] * Binv[k][i];
      if(fabs(top) < eps)
      {
         printf("RlTbLP2_e : Base matrix is singular\a\n");
         abort();
      }
      top = 1.0/top;
      for(i=0; i<nap1; i++) Binv[k][i] *= top;
      for(j=0; j<nap1; j++)
         if(j != k)
         {
            top = a[ell][0] * Binv[j][0];
            for(i=1; i<nap1; i++) top += a[ell][i] * Binv[j][i];
            for(i=0; i<nap1; i++) Binv[j][i] -= top * Binv[k][i];
         }
      Bidx[k] = ell;
   }
   if(Bidx[nap1-1] != NumCol)  /* change x, Bidx, Binv into the ones */
   {                           /* without the x(na) = epsilon variable */
      k = -1;
      j = -1;
      while(j < nap1-1)
      {
         j = j + 1;
         if(Bidx[j] == NumCol)
	 { 
            k = j;
            j = nap1-1;
            break;
         }
      }
      if(k == -1)   /* for testing */
      {
         printf("RlTbLP2_e : no index ma in Bidx \a\n");
         abort();       /* not happen in our case */
      }
      Bidx[k] = Bidx[nap1-1];
      for(i=0; i<na; i++) Binv[k][i] = Binv[nap1-1][i];
   }
   for(i=0; i<ma; i++) a[LPidx[i]][nap1-1] = 0.0; /* reset values WRT epsilon */

   *info = 0;
   for(i=0; i<na; i++)
      if(Bidx[i] > -1) ++(*info);

   if(*info == na) return;   /* case of full rank */

   k1 = -1;
   k = 0;
   while(k < na)
   { 
      if(Bidx[k] > -1 || k == k1)
      {
         ++k;
         continue;
      }

      for(i=0; i<na; i++) v[i] = Binv[k][i];

      flag = 1;
      ell = -1;                   /* step 2 : search incoming constraint */
      while(ell == -1)
      {
         sigj = DBL_MAX;
         for(i1=0; i1<ma; i1++)
         {
            i = LPidx[i1];
            ibrk = 0;
            for(j=0; j<na; j++)
               if(Bidx[j] == i)
               {
                  ibrk = 1;
                  break;
               }
            if(ibrk) continue;

            bot = a[i][0] * v[0];
            for(j=1; j<na; j++) bot += a[i][j] * v[j];
            if(bot >= -eps) continue;

            top = -a[i][nv1];
            for(j=0; j<na; j++) top += a[i][j] * x[j];
            top /= bot;
            if(top >= sigj) continue;

            sigj = top;
            ell = i;
         }
         ibrk = 0;
         if(ell < 0) 
            if(flag == 1)
            {
               for(i=0; i<na; i++) v[i] = -Binv[k][i];
               flag = 0;
            }
            else
            {
               k = k + 1;
               ibrk = 1;
               break;
            }
      }

      if(ibrk) continue;
                                           /* step 3 : update x, Bidx, Binv */
      for(i=0; i<na; i++) x[i] -= sigj * v[i]; 

      top = a[ell][0] * Binv[k][0];
      for(i=1; i<na; i++) top += a[ell][i] * Binv[k][i];
      if(fabs(top) < eps)
      {
         printf("RlTbLP2_e : Base matrix is singular\a\n");
         abort();
      } 
      top = 1.0 / top;
      for(i=0; i<na; i++) Binv[k][i] *= top;
      for(j=0; j<na; j++)
         if(j != k)
         {
            top = a[ell][0] * Binv[j][0];
            for(i=1; i <na; i++) top += a[ell][i] * Binv[j][i];
            for(i=0; i <na; i++) Binv[j][i] -= top * Binv[k][i];
         }
      Bidx[k] = ell;
      ++(*info);

      if(*info == na)  
         return;
      else 
      {
         k1 = k;
         k = 0;
      }   
   }
   *info = - *info;   /* the case of not full rank */
}

void dlp2_1pt_i
 ( int ma, int na, double **a, int nv1, double *c, int *LPidx, int FixPt,
   int TstPt, int *Bidx, double *x, double **Binv, int **RelTab )
{
   int ibrk;
   int i,i1,j,k,ell;
   double eps=1.0e-6,vmax,top,bot,sigj,dtmp;

   if(verbose_reltab > 0)
      printf("entered dlp2_1pt_i with ma = %d and na = %d\n",ma,na);
 
   while(1)
   {
      vmax = -DBL_MAX;           /* step 1 : search for outgoing constraint */
      for(i=0; i<na; i++)
         if(Bidx[i] != -1)
         {
            dtmp = c[0] * Binv[i][0];
            for(j=1; j<na; j++) dtmp += c[j] * Binv[i][j];
            if(dtmp > vmax)
            {
               vmax = dtmp;
               k = i;
            }
         }
      if(vmax < eps) return;     /* leave with optimal solution */

      sigj = DBL_MAX;            /* step 2 : search for incoming constraint */
      ell = -1;
      for(i1=0; i1<ma; i1++)
      {
         i = LPidx[i1];
         ibrk = 0;
         for(j=0; j<na; j++)
            if(Bidx[j] == i)
            {
               ibrk = 1;
               break;
            }
         if(ibrk) continue;

         bot = a[i][0] * Binv[k][0];
         for(j=1; j<na; j++) bot += a[i][j] * Binv[k][j];
         if(bot >= -eps) continue;

         top = -a[i][nv1];
         for(j=0; j<na; j++) top += a[i][j] * x[j];
         top /= bot;
         if(top >= sigj) continue;

         sigj = top;
         ell = i;
      }

      if(ell < 0)
      {
         printf("dlp2_1pt_i : LP unbounded\n\a");
         abort();       
      }                                    /* step 3 : update x, Bidx, Binv */

      for(i=0; i<na; i++) x[i] -= sigj * Binv[k][i];

      top = a[ell][0] * Binv[k][0];
      for(i=1; i<na; i++) top += a[ell][i] * Binv[k][i];
      if(fabs(top) < eps)
      {
         printf("dlp2_1pt_i : Base matrix is singular\a\n");
         abort();
      }
      top = 1.0/top;
      for(i=0; i<na; i++) Binv[k][i] *= top;
      for(j=0; j<na; j++)
         if(j != k)
         {
            top = a[ell][0] * Binv[j][0];
            for(i=1; i<na; i++) top += a[ell][i] * Binv[j][i];
            for(i=0; i<na; i++) Binv[j][i] -= top * Binv[k][i];
         }

      Bidx[k] = ell;

      if(ell >= TstPt) 
      {
         for(i=0; i< na; i++)
            if(i != k && Bidx[i] > -1 && !RelTab[ell][Bidx[i]])
               RelTab[ell][Bidx[i]] = RelTab[Bidx[i]][ell] = 1;

         if(!RelTab[ell][FixPt]) 
            RelTab[ell][FixPt] = RelTab[FixPt][ell] = 1;
      }
   }
}

void dlp2_1pt_s
 ( int ma, int na, double **a, int nv1, double *c, int *LPidx, int FixPt,
   int TstPt, int *Bidx, double *x, double **Binv, int **RelTab, L0_IML *L0 )
{
   int ibrk;
   int i,i1,j,k,ell;
   double eps=1.0e-6,vmax,top,bot,sigj,dtmp;
   int *J = (int*) calloc(2,sizeof(int));

   if(verbose_reltab > 0)
      printf("entered dlp2_1pt_s with ma = %d and na = %d\n",ma,na);
 
   while(1)
   {
      vmax = -DBL_MAX;            /* step 1 : search for outgoing constraint */
      for(i=0; i<na; i++) 
         if(Bidx[i] != -1)
         {
            dtmp = c[0] * Binv[i][0];
            for(j=1; j<na; j++) dtmp += c[j] * Binv[i][j];
            if(dtmp > vmax)
            {
               vmax = dtmp;
               k = i;
            }
         }
 
      if(vmax < eps) return;      /* leave with optimal solution */

      sigj = DBL_MAX;             /* step 2 : search for incoming constraint */
      ell = -1;
      for(i1=0; i1<ma; i1++)
      {
         i = LPidx[i1];
         ibrk = 0;
         for(j=0; j<na; j++)  
            if(Bidx[j] == i)
            {
               ibrk = 1;
               break;
            }
         if(ibrk) continue;

         bot = a[i][0] * Binv[k][0];
         for(j=1; j<na; j++) bot += a[i][j] * Binv[k][j];
         if(bot >= -eps) continue;

         top = -a[i][nv1];
         for(j=0; j<na; j++) top += a[i][j] * x[j];
         top /= bot;
         if(top >= sigj) continue;

         sigj = top;
         ell = i;
      }

      if(ell < 0)
      {
         printf("dlp2_1pt_s : LP unbounded\n\a");
         abort();       
      }                               /* step 3 : update x, Bidx, Binv */

      for(i=0; i<na; i++) x[i] -= sigj * Binv[k][i];

      top = a[ell][0] * Binv[k][0];
      for(i=1; i<na; i++) top += a[ell][i] * Binv[k][i];
      if(fabs(top) < eps)
      {
         printf("dlp2_1pt_s : Base matrix is singular\a\n");
         abort();
      }
      top = 1.0/top;
      for(i=0; i<na; i++) Binv[k][i] *= top;
      for(j=0; j<na; j++)
         if(j != k)
         {
            top = a[ell][0] * Binv[j][0];
            for(i=1; i<na; i++) top += a[ell][i] * Binv[j][i];
            for(i=0; i<na; i++) Binv[j][i] -= top * Binv[k][i];
         }

      Bidx[k] = ell;
      if(ell >= TstPt && !RelTab[ell][FixPt])   /* save Bidx, x, Binv */
      {
         RelTab[ell][FixPt] = RelTab[FixPt][ell] = 1;
         J[0] = FixPt;
         J[1] = ell;                            /* J[0] < J[1] is assumed. */
         L0_Add2(L0,J,na,Bidx,x,Binv);          /* add 2 points */
      }
   }
   free(J);
}

void dlp1_1pt_i
 ( int ma, int na, double **a, int nv1, double *c, int *LPidx, int FixPt, 
   int TstPt, int *Bidx, double *x, double **Binv, int **RelTab )
{
   int ibrk,i,i1,j,k,ell;
   double eps=1.0e-6,vmax,top,bot,sigj,dtmp;

   if(verbose_reltab > 0)
      printf("entered dlp1_1pt_i with ma = %d and na = %d\n",ma,na);
 
   while(1)
   {
      vmax = -DBL_MAX;           /* step 1 : search for outgoing constraint */
      for(i=0; i<na; i++) 
      {
         dtmp = c[0] * Binv[i][0];
         for(j=1; j<na; j++) dtmp += c[j] * Binv[i][j];
         if(dtmp > vmax)
	 {
            vmax = dtmp;
            k = i;
         }
      }
      if(vmax < eps) return;     /* out with optimal solution */

      sigj = DBL_MAX;            /* step 2 : find incoming constraint */
      ell = -1;
      for(i1=0; i1<ma; i1++)
      {
         i = LPidx[i1];
         ibrk = 0;
         for(j=0; j<na; j++)  
            if(Bidx[j] == i)
            {
               ibrk = 1;
               break;
            }
         if(ibrk) continue;

         bot = a[i][0] * Binv[k][0];
         for(j=1; j<na; j++) bot += a[i][j] * Binv[k][j];
         if(bot >= -eps) continue;

         top = -a[i][nv1];
         for(j=0; j<na; j++) top += a[i][j] * x[j];
         top /= bot;
         if(top >= sigj) continue;

         sigj = top;
         ell = i;
      }
      if(ell < 0)
      {
         printf("dlp1_1pt_i : LP unbounded\n\a");
         abort();       
      }                                     /* step 3 : update x, Bidx, Binv */

      for(i=0; i<na; i++) x[i] -= sigj * Binv[k][i];

      top = a[ell][0] * Binv[k][0];
      for(i=1; i<na; i++)  top += a[ell][i] * Binv[k][i];
      if(fabs(top) < eps)
      {
         printf("dlp1_1pt_i : Base matrix is singular\a\n");
         abort();
      }
      top = 1.0/top;
      for(i=0; i<na; i++) Binv[k][i] *= top;
      for(j=0; j<na; j++)
         if(j != k)
         {
            top = a[ell][0]* Binv[j][0];
            for(i=1; i<na; i++) top += a[ell][i] * Binv[j][i];
            for(i=0; i<na; i++) Binv[j][i] -= top * Binv[k][i];
         }
      Bidx[k] = ell;

      if(ell >= TstPt) 
      {
         for(i=0; i<na; i++)
            if(i != k && !RelTab[ell][Bidx[i]])
               RelTab[ell][Bidx[i]] = RelTab[Bidx[i]][ell] = 1;

         if(!RelTab[ell][FixPt]) 
            RelTab[ell][FixPt] = RelTab[FixPt][ell] = 1;
      }
   }
}

void dlp1_1pt_s
 ( int ma, int na, double **a, int nv1, double *c, int *LPidx, int FixPt,
   int TstPt, int *Bidx, double *x, double **Binv, int **RelTab, L0_IML *L0 )
{
   int ibrk,i,i1,j,k,ell;
   double eps=1.0e-6,vmax,top,bot,sigj,dtmp;
   int *J = (int*) calloc(2,sizeof(int));

   if(verbose_reltab > 0)
      printf("entered dlp1_1pt_s with ma = %d and na = %d\n",ma,na);
 
   while(1)
   {
      vmax = -DBL_MAX;           /* step 1 : search outgoing constraint */
      for(i=0; i<na; i++) 
      {
         dtmp = c[0] * Binv[i][0];
         for(j=1; j<na; j++) dtmp += c[j] * Binv[i][j];
         if(dtmp > vmax)
         {
            vmax = dtmp;
            k = i;
         }
      }
 
      if(vmax < eps) return;     /* out with optimal solution */

      sigj = DBL_MAX;            /* step 2 : search incoming constraint */
      ell = -1;
      for(i1=0; i1<ma; i1++)
      {
         i = LPidx[i1];
         ibrk = 0;
         for(j=0; j<na; j++)
            if(Bidx[j] == i)
            {
               ibrk = 1;
               break;
            }
         if(ibrk) continue;

         bot = a[i][0] * Binv[k][0];
         for(j=1; j<na; j++) bot += a[i][j] * Binv[k][j];
         if(bot >= -eps) continue;

         top = -a[i][nv1];
         for(j=0; j<na; j++) top += a[i][j] * x[j];
         top /= bot;
         if(top >= sigj) continue;

         sigj = top;
         ell = i;
      }

      if(ell < 0)
      {
         printf("dlp1_1pt_s : LP unbounded\n\a");
         abort();       
      }                          /* step 3 : update x, Bidx, Binv */

      for(i=0; i<na; i++) x[i] -= sigj * Binv[k][i];

      top = a[ell][0] * Binv[k][0];
      for(i=1; i<na; i++) top += a[ell][i] * Binv[k][i];
      if(fabs(top) < eps)
      {
         printf("dlp1_1pt_s : Base matrix is singular\a\n");
         abort();
      }
      top = 1.0/top;
      for(i=0; i<na; i++) Binv[k][i] *= top;
      for(j=0; j<na; j++)
         if(j != k)
         {
            top = a[ell][0] * Binv[j][0];
            for(i=1; i<na; i++) top += a[ell][i] * Binv[j][i];
            for(i=0; i<na; i++) Binv[j][i] -= top * Binv[k][i];
         }

      Bidx[k] = ell;

      if(ell >= TstPt && !RelTab[ell][FixPt])   /* save Bidx, x, Binv */
      {
         RelTab[ell][FixPt] = RelTab[FixPt][ell] = 1;
         J[0] = FixPt;
         J[1] = ell;                            /* J[0] < J[1] is assumed .*/
         L0_Add2(L0,J,na,Bidx,x,Binv);          /* add 2 points */
      }
   }
   free(J);
}
