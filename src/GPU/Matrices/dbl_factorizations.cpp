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

void CPU_dbl_factors_backward ( int dim, double **U, double *b, double *x )
{
   for(int i=dim-1; i>=0; i--)
   {
      x[i] = b[i];
      for(int j=dim-1; j>i; j--) x[i] = x[i] - U[i][j]*x[j];
      x[i] = x[i]/U[i][i];
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

void CPU_dbl_factors_lusolve
 ( int dim, double **A, int *pivots, double *b, double *x )
{
   CPU_dbl_factors_lufac(dim,A,pivots);
   for(int i=0; i<dim; i++) x[i] = b[pivots[i]];
   CPU_dbl_factors_forward(dim,A,x,b);
   CPU_dbl_factors_backward(dim,A,b,x);
}
