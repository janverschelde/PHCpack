/* The file dbl_tabs_host.cpp defines the functions specified in
 * the file dbl_tabs_host.h. */

#include <cstdlib>
#include <ctime>
#include "dbl_factorizations.h"
#include "dbl_tabs_host.h"

void CPU_dbl_backsubs ( int dim, double **U, double *b, double *x )
{
   CPU_dbl_factors_backward(dim,U,b,x);
}

void CPU_cmplx_backsubs
 ( int dim, double **Ure, double **Uim, double *bre, double *bim,
   double *xre, double *xim )
{
   CPU_cmplx_factors_backward(dim,Ure,Uim,bre,bim,xre,xim);
}

void CPU_dbl_upper_inverse
 ( int dim, double **U, double **invU, double *lapsec )
{
   double *col = new double[dim];
   double *rhs = new double[dim];

   clock_t start = clock();

   for(int i=0; i<dim; i++) rhs[i] = 0.0;

   for(int j=0; j<dim; j++) // compute j-th column of the inverse
   {
      rhs[j] = 1.0;
      CPU_dbl_backsubs(dim,U,rhs,col);
      for(int i=0; i<dim; i++) invU[i][j] = col[i];
      rhs[j] = 0.0;
   }
   clock_t end = clock();
   *lapsec = double(end - start)/CLOCKS_PER_SEC;

   free(rhs); free(col);
}

void CPU_cmplx_upper_inverse
 ( int dim, double **Ure, double **Uim, double **invUre, double **invUim,
   double *lapsec )
{
   double *colre = new double[dim];
   double *colim = new double[dim];
   double *rhsre = new double[dim];
   double *rhsim = new double[dim];

   clock_t start = clock();

   for(int i=0; i<dim; i++)
   {
      rhsre[i] = 0.0;
      rhsim[i] = 0.0;
   }
   for(int j=0; j<dim; j++) // compute j-th column of the inverse
   {
      rhsre[j] = 1.0;
      rhsim[j] = 0.0;

      CPU_cmplx_backsubs(dim,Ure,Uim,rhsre,rhsim,colre,colim);

      for(int i=0; i<dim; i++)
      {
         invUre[i][j] = colre[i];
         invUim[i][j] = colim[i];
      }
      rhsre[j] = 0.0;
      rhsim[j] = 0.0;
   }
   clock_t end = clock();
   *lapsec = double(end - start)/CLOCKS_PER_SEC;

   free(rhsre); free(colre);
   free(rhsim); free(colim);
}

void CPU_dbl_matmatmul ( int dim, double **A, double **F )
{
   double **result = new double*[dim];
   for(int i=0; i<dim; i++) result[i] = new double[dim];

   for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++)
      {
         result[i][j] = 0.0;
         for(int k=0; k<dim; k++)
            result[i][j] = result[i][j] + F[i][k]*A[k][j];
      }

   for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++) A[i][j] = result[i][j];

   for(int i=0; i<dim; i++) free(result[i]);
   free(result);
}

void CPU_dbl_upper_tiled_solver
 ( int dim, int szt, int nbt, double **U, double *b, double *x,
   double *lapsec )
{
   double **T = new double*[szt];
   double **invT = new double*[szt];
   for(int i=0; i<szt; i++)
   {
      T[i] = new double[szt];
      invT[i] = new double[szt];
   }
   double timelapsed;

   clock_t start = clock();

   int idx = (nbt-1)*szt;
   for(int i=0; i<szt; i++)
      for(int j=0; j<szt; j++) T[i][j] = U[idx+i][idx+j];

   CPU_dbl_upper_inverse(szt,T,invT,&timelapsed);

   for(int i=0; i<szt; i++)
      for(int j=0; j<szt; j++) U[idx+i][idx+j] = invT[i][j];

   for(int i=0; i<szt; i++)
   {
      x[idx+i] = 0.0;
      for(int j=0; j<szt; j++) x[idx+i] = x[idx+i] + invT[i][j]*b[idx+j];
   }
   double *wb = new double[szt];    // work space for b
   double **wT = new double*[szt];  // work space for a tile
   for(int i=0; i<szt; i++)
      wT[i] = new double[szt];

   double prod;

   for(int k=nbt-1; k>0; k--)  // update with solution tile k
   {
      idx = idx - szt; // idx is start index of diagonal tile

      for(int i=0; i<szt; i++)
         for(int j=0; j<szt; j++) T[i][j] = U[idx+i][idx+j];

      CPU_dbl_upper_inverse(szt,T,invT,&timelapsed); // invert diagonal tile

      for(int i=0; i<szt; i++)
         for(int j=0; j<szt; j++) U[idx+i][idx+j] = invT[i][j];

      for(int L=0; L<k; L++)   // update wb as many times as k
      {
         int rowidx = L*szt;

         for(int i=0; i<szt; i++) // load the work space
         {
            wb[i] = b[rowidx+i];
            for(int j=0; j<szt; j++) wT[i][j] = U[rowidx+i][idx+szt+j];
         }
         for(int i=0; i<szt; i++) // update wb
         {
            prod = 0.0;
            for(int j=0; j<szt; j++) prod = prod + wT[i][j]*x[idx+szt+j];
            wb[i] = wb[i] - prod;
         }
         for(int i=0; i<szt; i++) b[rowidx+i] = wb[i]; // for next update
      }
      for(int i=0; i<szt; i++)   // wb = invT*b
      {
         prod = 0.0;
         for(int j=0; j<szt; j++) prod = prod + invT[i][j]*b[idx+j];
         wb[i] = prod;
      }
      for(int i=0; i<szt; i++) x[idx+i] = wb[i];
   }
   clock_t end = clock();
   *lapsec = double(end - start)/CLOCKS_PER_SEC;

   for(int i=0; i<szt; i++)
   {
      free(T[i]); free(invT[i]); free(wT[i]);
   }
   free(T); free(invT); free(wb); free(wT);
}

void CPU_cmplx_upper_tiled_solver
 ( int dim, int szt, int nbt, double **Ure, double **Uim,
   double *bre, double *bim, double *xre, double *xim, double *lapsec )
{
   double **Tre = new double*[szt];
   double **Tim = new double*[szt];
   double **invTre = new double*[szt];
   double **invTim = new double*[szt];

   for(int i=0; i<szt; i++)
   {
      Tre[i] = new double[szt];
      Tim[i] = new double[szt];
      invTre[i] = new double[szt];
      invTim[i] = new double[szt];
   }
   double timelapsed;

   clock_t start = clock();

   int idx = (nbt-1)*szt;
   for(int i=0; i<szt; i++)
      for(int j=0; j<szt; j++)
      {
         Tre[i][j] = Ure[idx+i][idx+j];
         Tim[i][j] = Uim[idx+i][idx+j];
      }

   CPU_cmplx_upper_inverse(szt,Tre,Tim,invTre,invTim,&timelapsed);

   for(int i=0; i<szt; i++)
      for(int j=0; j<szt; j++)
      {
         Ure[idx+i][idx+j] = invTre[i][j];
         Uim[idx+i][idx+j] = invTim[i][j];
      }

   double accre,accim;

   for(int i=0; i<szt; i++)
   {
      xre[idx+i] = 0.0;
      xim[idx+i] = 0.0;
      for(int j=0; j<szt; j++) // x[idx+i] = x[idx+i] + invT[i][j]*b[idx+j];
      {
         accre = invTre[i][j]*bre[idx+j] - invTim[i][j]*bim[idx+j];
         accim = invTim[i][j]*bre[idx+j] + invTre[i][j]*bim[idx+j];
         xre[idx+i] = xre[idx+i] + accre;
         xim[idx+i] = xim[idx+i] + accim;
      }
   }
   double *wbre = new double[szt];    // work space for b
   double *wbim = new double[szt];    // work space for b
   double **wTre = new double*[szt];  // work space for a tile
   double **wTim = new double*[szt];  // work space for a tile

   for(int i=0; i<szt; i++)
   {
      wTre[i] = new double[szt];
      wTim[i] = new double[szt];
   }
   double prodre,prodim;

   for(int k=nbt-1; k>0; k--)  // update with solution tile k
   {
      idx = idx - szt; // idx is start index of diagonal tile

      for(int i=0; i<szt; i++)
         for(int j=0; j<szt; j++)
         {
            Tre[i][j] = Ure[idx+i][idx+j];
            Tim[i][j] = Uim[idx+i][idx+j];
         }

      // invert diagonal tile
      CPU_cmplx_upper_inverse(szt,Tre,Tim,invTre,invTim,&timelapsed);

      for(int i=0; i<szt; i++)
         for(int j=0; j<szt; j++)
         {
            Ure[idx+i][idx+j] = invTre[i][j];
            Uim[idx+i][idx+j] = invTim[i][j];
         }

      for(int L=0; L<k; L++)   // update wb as many times as k
      {
         int rowidx = L*szt;

         for(int i=0; i<szt; i++) // load the work space
         {
            wbre[i] = bre[rowidx+i];
            wbim[i] = bim[rowidx+i];
            for(int j=0; j<szt; j++)
            {
               wTre[i][j] = Ure[rowidx+i][idx+szt+j];
               wTim[i][j] = Uim[rowidx+i][idx+szt+j];
            }
         }
         for(int i=0; i<szt; i++) // update wb
         {
            prodre = 0.0;
            prodim = 0.0;
            for(int j=0; j<szt; j++) // prod = prod + wT[i][j]*x[idx+szt+j];
            {
               accre = wTre[i][j]*xre[idx+szt+j] - wTim[i][j]*xim[idx+szt+j];
               accim = wTim[i][j]*xre[idx+szt+j] + wTre[i][j]*xim[idx+szt+j];
               prodre = prodre + accre;
               prodim = prodim + accim;
            }
            wbre[i] = wbre[i] - prodre;
            wbim[i] = wbim[i] - prodim;
         }
         for(int i=0; i<szt; i++)
         {
            bre[rowidx+i] = wbre[i]; // for next update
            bim[rowidx+i] = wbim[i];
         }
      }
      for(int i=0; i<szt; i++)   // wb = invT*b
      {
         prodre = 0.0;
         prodim = 0.0;
         for(int j=0; j<szt; j++) // prod = prod + invT[i][j]*b[idx+j];
         {
            accre = invTre[i][j]*bre[idx+j] - invTim[i][j]*bim[idx+j];
            accim = invTim[i][j]*bre[idx+j] + invTre[i][j]*bim[idx+j];
            prodre = prodre + accre;
            prodim = prodim + accim;
         }
         wbre[i] = prodre;
         wbim[i] = prodim;
      }
      for(int i=0; i<szt; i++)
      {
         xre[idx+i] = wbre[i];
         xim[idx+i] = wbim[i];
      }
   }
   clock_t end = clock();
   *lapsec = double(end - start)/CLOCKS_PER_SEC;

   for(int i=0; i<szt; i++)
   {
      free(Tre[i]); free(invTre[i]); free(wTre[i]);
      free(Tim[i]); free(invTim[i]); free(wTim[i]);
   }
   free(Tre); free(invTre); free(wbre); free(wTre);
   free(Tim); free(invTim); free(wbim); free(wTim);
}
