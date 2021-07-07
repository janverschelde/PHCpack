/* The file dbl_factorizations.cpp defines the functions specified in
 * the file dbl_factorizations.h. */

#include <cmath>
#include "dbl_factorizations.h"

void CPU_dbl_factors_matmatmul
 ( int rows, int dim, int cols, double **A, double **B, double **C )
{
   for(int i=0; i<rows; i++)
      for(int j=0; j<cols; j++)
      {
         C[i][j] = 0.0;
         for(int k=0; k<dim; k++)
            C[i][j] = C[i][j] + A[i][k]*B[k][j];
      }
}

void CPU_dbl_factors_forward ( int dim, double **L, double *b, double *x )
{
   x[0] = b[0];
   for(int i=1; i<dim; i++)
   {
      x[i] = b[i];
      for(int j=0; j<i; j++) x[i] = x[i] - L[i][j]*x[j];
   }
}

void CPU_cmplx_factors_forward
 ( int dim, double **Lre, double **Lim, double *bre, double *bim,
   double *xre, double *xim )
{
   double accre,accim;

   xre[0] = bre[0];
   xim[0] = bim[0];
   for(int i=1; i<dim; i++)
   {
      xre[i] = bre[i];
      xim[i] = bim[i];
      for(int j=0; j<i; j++) // x[i] = x[i] - L[i][j]*x[j];
      {
         accre = Lre[i][j]*xre[j] - Lim[i][j]*xim[j];
         accim = Lim[i][j]*xre[j] + Lre[i][j]*xim[j];
         xre[i] = xre[i] - accre;
         xim[i] = xim[i] - accim;
      }
   }
}

void CPU_dbl_factors_backward ( int dim, double **U, double *b, double *x )
{
   for(int i=dim-1; i>=0; i--)
   {
      x[i] = b[i];
      for(int j=dim-1; j>i; j--) x[i] = x[i] - U[i][j]*x[j];
      x[i] = x[i]/U[i][i];
   }
}

void CPU_cmplx_factors_backward
 ( int dim, double **Ure, double **Uim, double *bre, double *bim,
   double *xre, double *xim )
{
   double accre,accim,det,zre,zim;

   for(int i=dim-1; i>=0; i--)
   {
      xre[i] = bre[i];
      xim[i] = bim[i];
      for(int j=dim-1; j>i; j--) // x[i] = x[i] - U[i][j]*x[j];
      {
         accre = Ure[i][j]*xre[j] - Uim[i][j]*xim[j];
         accim = Uim[i][j]*xre[j] + Ure[i][j]*xim[j];
         xre[i] = xre[i] - accre;
         xim[i] = xim[i] - accim;
      }
      // x[i] = x[i]/U[i][i];
      det = Ure[i][i]*Ure[i][i] + Uim[i][i]*Uim[i][i];
      accre = Ure[i][i]/det;   // accre is the real part of 1/U[i][i]
      accim = -Uim[i][i]/det;  // accim is imaginary part of 1/U[i][i]
      zre = xre[i]*accre - xim[i]*accim;
      zim = xim[i]*accre + xre[i]*accim;
      xre[i] = zre;
      xim[i] = zim;
   }
}

void CPU_dbl_factors_lufac ( int dim, double **A, int *pivots )
{
   double valmax,valtmp;
   int idxmax,idxtmp;

   for(int j=0; j<dim; j++) pivots[j] = j;
   for(int j=0; j<dim; j++)
   {
      valmax = fabs(A[j][j]); idxmax = j;     // find the pivot
      for(int i=j+1; i< dim; i++)
      {
         valtmp = fabs(A[i][j]);
         if(valtmp > valmax)
         {
            valmax = valtmp; idxmax = i;
         }
      }
      if(idxmax != j)                        // swap rows
      {
         for(int k=0; k<dim; k++)
         {
            valtmp = A[idxmax][k];
            A[idxmax][k] = A[j][k];
            A[j][k] = valtmp;
         }
         idxtmp = pivots[idxmax];
         pivots[idxmax] = pivots[j];
         pivots[j] = idxtmp;
      }
      for(int i=j+1; i<dim; i++)             // reduce
      {
         A[i][j] = A[i][j]/A[j][j];
         for(int k=j+1; k<dim; k++) 
            A[i][k] = A[i][k] - A[i][j]*A[j][k];
      }
   }
}

void CPU_cmplx_factors_lufac
 ( int dim, double **Are, double **Aim, int *pivots )
{
   double valmax,valtmp,accre,accim,det,zre,zim;
   int idxmax,idxtmp;

   for(int j=0; j<dim; j++) pivots[j] = j;
   for(int j=0; j<dim; j++)
   {
      valmax = fabs(Are[j][j]) + fabs(Aim[j][j]); // find the pivot
      idxmax = j;
      for(int i=j+1; i< dim; i++)
      {
         valtmp = fabs(Are[i][j]) + fabs(Aim[i][j]);
         if(valtmp > valmax)
         {
            valmax = valtmp; idxmax = i;
         }
      }
      if(idxmax != j)                        // swap rows
      {
         for(int k=0; k<dim; k++)
         {
            valtmp = Are[idxmax][k];
            Are[idxmax][k] = Are[j][k];
            Are[j][k] = valtmp;
            valtmp = Aim[idxmax][k];
            Aim[idxmax][k] = Aim[j][k];
            Aim[j][k] = valtmp;
         }
         idxtmp = pivots[idxmax];
         pivots[idxmax] = pivots[j];
         pivots[j] = idxtmp;
      }
      for(int i=j+1; i<dim; i++)             // reduce
      {
         // A[i][j] = A[i][j]/A[j][j];
         det = Are[j][j]*Are[j][j] + Aim[j][j]*Aim[j][j];
         accre = Are[j][j]/det;
         accim = -Aim[j][j]/det;
         zre = Are[i][j]*accre - Aim[i][j]*accim;
         zim = Aim[i][j]*accre + Are[i][j]*accim;
         Are[i][j] = zre;
         Aim[i][j] = zim;
         for(int k=j+1; k<dim; k++) // A[i][k] = A[i][k] - A[i][j]*A[j][k];
         {
            accre = Are[i][j]*Are[j][k] - Aim[i][j]*Aim[j][k];
            accim = Aim[i][j]*Are[j][k] + Are[i][j]*Aim[j][k];
            Are[i][k] = Are[i][k] - accre;
            Aim[i][k] = Aim[i][k] - accim;
         }
      }
   }
}

void CPU_dbl_factors_lusolve
 ( int dim, double **A, int *pivots, double *b, double *x )
{
   CPU_dbl_factors_lufac(dim,A,pivots);
   for(int i=0; i<dim; i++) x[i] = b[pivots[i]];
   CPU_dbl_factors_forward(dim,A,x,b);
   CPU_dbl_factors_backward(dim,A,b,x);
}

void CPU_cmplx_factors_lusolve
 ( int dim, double **Are, double **Aim, int *pivots,
   double *bre, double *bim, double *xre, double *xim )
{
   CPU_cmplx_factors_lufac(dim,Are,Aim,pivots);
   for(int i=0; i<dim; i++)
   {
      xre[i] = bre[pivots[i]];
      xim[i] = bim[pivots[i]];
   }
   CPU_cmplx_factors_forward(dim,Are,Aim,xre,xim,bre,bim);
   CPU_cmplx_factors_backward(dim,Are,Aim,bre,bim,xre,xim);
}
