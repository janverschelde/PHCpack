/* The file "one_level_lp1.c" collects the definitions of the functions 
 * with prototypes in "one_level_lp1.h". */

#include <string.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "one_level_lp1.h"

#define verbose_sort 0
#define verbose_one_level_lp 0
#define verbose_dnulp2_a 0

void one_level_LP
 ( int Strt1Pt, int End1Pt, int *PtIn, int LPdim, double **A,
   int nVar, double *x, double **Binv, int *Bidx, IT_LP *ItLp )
{
   int i,j,info,TstPt;
   int *labels = (int*)calloc(nVar,sizeof(int));
   double *c = (double*)calloc(nVar+2,sizeof(double));

   if(verbose_one_level_lp > 0)
   {
      printf("entered one_level_LP ");
      printf("with Strt1Pt = %d and End1Pt = %d\n",Strt1Pt,End1Pt); 
   }
                           /* extend the cells by using 1-point test */
   TstPt = Strt1Pt;
   for(i=0; i<LPdim; i++) c[i] = -A[TstPt][i];
   dnulp2_a(End1Pt+1,LPdim,A,nVar,c,Bidx,x,Binv,&info);
   j = -1;     /* record any possible points in Bidx passed 1-point test */
   for(i=0; i<LPdim; i++)
      if(Bidx[i] >= TstPt)
      {
         PtIn[Bidx[i]] = 1;
         labels[++j] = Bidx[i];
      }
   if(++j > 0)
   {
      Sort(j,labels);
      IT_Add1(ItLp,j,labels,LPdim,Bidx,x,Binv);
   }
   if(info >= 0)             /* perform 1-point test for other points */
   {
      for(TstPt=Strt1Pt+1; TstPt<=End1Pt; TstPt++)
      {  
         if(!PtIn[TstPt])
         {
            for(i=0; i<LPdim; i++) c[i] = -A[TstPt][i];
            dlp1_1pts(End1Pt+1,LPdim,A,nVar,c,
                      TstPt,Bidx,x,Binv,PtIn,ItLp);
         }
      }
   }
   else
   {
      for(TstPt=Strt1Pt+1; TstPt<=End1Pt; TstPt++)
      {   
         if(!PtIn[TstPt])   
         {
            for(i=0; i<LPdim; i++) c[i] = -A[TstPt][i];
            dlp2_1pts(End1Pt+1,LPdim,A,nVar,c,
                      TstPt,Bidx,x,Binv,PtIn,ItLp);
         }
      }
   }
   free(labels); free(c);
}  

void dnulp2_a 
 ( int ma, int na, double **a, int nVar, double *c, int *Bidx,
   double *x, double **Binv, int *info )
{
   int ibrk,ell,k,k1,i,j,flag;
   double eps=1.0e-6,sigj,vmax,bot,top,dtmp;
   double v[na];

   if(verbose_dnulp2_a > 0)
   {
      printf("entered dnulp2_a with ma = %d and na = %d\n",ma,na);
   }

   while(1)
   { 
      ibrk = 0;                 /* step 1 : search outgoing constraint */
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
               dtmp = c[0]*Binv[i][0];
               for(j=1; j<na; j++) dtmp += c[j]*Binv[i][j];
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
               dtmp = c[0]*Binv[i][0];
               for(j=1; j<na; j++) dtmp += c[j]*Binv[i][j];
               if(dtmp > vmax)
               {
                  vmax = dtmp;
                  k = i;
               }
            }
         if(vmax < eps) break; /* found optimal solution, solve smaller LP */
         for(i=0; i<na; i++) v[i] = Binv[k][i];
      }
      sigj = DBL_MAX;               /* step 2 : search incoming constraint */
      ell = -1;
      for(i=0; i<ma; i++)
      {
         ibrk = 0;
         for(j=0; j<na; j++)
            if(Bidx[j] == i)
            {
               ibrk = 1;
               break;
            }
         if(ibrk) continue;
         bot = a[i][0]*v[0];
         for(j=1; j<na; j++) bot += a[i][j]*v[j];
         if(bot >= -eps) continue;
         top = -a[i][nVar];
         for(j=0; j<na; j++) top += a[i][j]*x[j];
         top /= bot;
         if(top >= sigj) continue;
         sigj = top;
         ell = i;
      }
      if(ell < 0)
      {
         printf("dnulp2_a : LP unbounded\n\a");
         abort();
      }                                 /* step 3: update x, Bidx, Binv */
      for(i=0; i<na; i++) x[i] -= sigj * v[i];
      top = a[ell][0] * Binv[k][0];
      for(i=1; i<na; i++) top += a[ell][i] * Binv[k][i];
      if(fabs(top) < eps)
      {
         printf("dnulp2_a : Base matrix is singular\a\n");
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
      if(Bidx[i] > -1) ++*info;

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
         for(i=0; i<ma; i++)
         {
            ibrk = 0;
            for(j=0; j<na; j++)
               if(Bidx[j] == i)
               {
                  ibrk = 1;
                  break;
               }
            if(ibrk) continue;
            bot = a[i][0]*v[0];
            for(j=1; j<na; j++) bot += a[i][j]*v[j];
            if(bot >= -eps) continue;

            top = -a[i][nVar];
            for(j=0; j<na; j++) top += a[i][j]*x[j];
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
                                       /* step 3: update x, Bidx, Binv */
      for(i=0; i<na; i++) x[i] -= sigj * v[i];
      top = a[ell][0] * Binv[k][0];
      for(i=1; i<na; i++) top += a[ell][i] * Binv[k][i];
      if(fabs(top) < eps)
      {
         printf("dnulp2_a : Base matrix is singular\a\n");
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

void dlp2_1pts
 ( int ma, int na, double **a, int nVar, double *c, int TstPt,
   int *Bidx, double *x, double **Binv, int *PtIn, IT_LP *ItLp )
{
   int  ibrk,ell,k,i,j;
   double eps=1.0e-6,sigj,vmax,bot,top,dtmp;

   while(1)
   {
      vmax = -DBL_MAX;              /* step 1 : find outgoing constraint */
      for(i=0; i<na; i++)
         if(Bidx[i] != -1)
         {
            dtmp = c[0]*Binv[i][0];
            for(j=1; j<na; j++) dtmp += c[j]*Binv[i][j];
            if(dtmp > vmax)
            {
               vmax = dtmp;
               k = i;
            }
         }
      if(vmax < eps) return;              /* leave with optimal solution */
      sigj = DBL_MAX;             /* step 2 : search incoming constraint */
      ell = -1;
      for(i=0; i<ma; i++)
      {
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
         top = -a[i][nVar];
         for(j=0; j<na; j++) top += a[i][j] * x[j];
         top /= bot;
         if(top >= sigj) continue;
         sigj = top;
         ell = i;
      }
      if(ell < 0)
      {
         printf("dlp2_1pts : LP unbounded\n\a");
         abort();
      }                                /* step 3 : update x, Bidx, Binv */

      for(i=0; i<na; i++) x[i] -= sigj * Binv[k][i];
      top = a[ell][0] * Binv[k][0];
      for(i=1; i<na; i++) top += a[ell][i] * Binv[k][i];
      if(fabs(top) < eps)
      {
         printf("dlp2_1pts : Base matrix is singular\a\n");
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
      if(ell >= TstPt && !PtIn[ell])          /* save Bidx, x, Binv */
      {
         PtIn[ell] = 1;
         IT_Add2(ItLp,ell,na,Bidx,x,Binv);
      }
   }
}

void dlp1_1pts
 ( int ma, int na, double **a, int nVar, double *c, int TstPt,
   int *Bidx, double *x, double **Binv, int *PtIn, IT_LP *ItLp )
{
   int ibrk,i,j,k,ell;
   double eps=1.0e-6,vmax,top,bot,sigj,dtmp;

   while(1)
   {
      vmax = -DBL_MAX;        /* step 1 : search for outgoing constraint */
      for(i=0; i<na; i++)
      {
         dtmp = c[0]*Binv[i][0];
         for(j=1; j<na; j++) dtmp += c[j]*Binv[i][j];
         if(dtmp > vmax)
         {
            vmax = dtmp;
            k = i;
         }
      }
      if (vmax < eps) return;             /* leave with optimal solution */
      sigj = DBL_MAX;         /* step 2 : search for incoming constraint */
      ell = -1;
      for(i=0; i<ma; i++)
      {
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
         if(bot > -eps) continue;
         top = -a[i][nVar];
         for(j=0; j<na; j++) top += a[i][j] * x[j];
         top /= bot;
         if(top > sigj) continue;
         sigj = top;
         ell = i;
      }
      if(ell < 0)
      {
         printf("dlp1_1pts.cpp: LP unbounded\n\a");
         abort();
      }                                  /* step 3: update x, Bidx, Binv */
      for(i=0; i<na; i++ ) x[i] -= sigj * Binv[k][i];
      top = a[ell][0] * Binv[k][0];
      for(i=1; i<na; i++) top += a[ell][i] * Binv[k][i];
      if(fabs(top) < eps)
      {
         printf("dlp1_1pts : Base matrix is singular\a\n");
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
      if(ell >= TstPt && !PtIn[ell])              /* save Bidx, x, Binv */
      {
         PtIn[ell] = 1;
         IT_Add2(ItLp,ell,na,Bidx,x,Binv);
      }
   }
}

void Sort ( int n, int *a )
{
   int j,itmp,i;

   if(verbose_sort>0)
   {
      printf("array (n=%d) before sort :",n);
      for(i=0; i<n; i++) printf(" %d",a[i]);
      printf("\n");
   }
   for(i=1; i<n; i++)
   {
      itmp = a[i];
      for(j = i; j>0 && itmp<a[j-1]; j--) a[j] = a[j-1];
      a[j] = itmp;
   }
   if(verbose_sort>0)
   {
      printf("array after sort :");
      for(i=0; i<n; i++) printf(" %d",a[i]);
      printf("\n");
   }
}
