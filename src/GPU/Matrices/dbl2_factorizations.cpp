/* The file dbl2_factorizations.cpp defines functions specified in
 * the file dbl2_factorizations.h. */

#include <cmath>
#include "double_double_functions.h"
#include "dbl2_factorizations.h"

void CPU_dbl2_factors_matmatmul
 ( int rows, int dim, int cols, double **Ahi, double **Alo,
   double **Bhi, double **Blo, double **Chi, double **Clo )
{
   double acchi,acclo;

   for(int i=0; i<rows; i++)
      for(int j=0; j<cols; j++)
      {
         Chi[i][j] = 0.0;
         Clo[i][j] = 0.0;
         for(int k=0; k<dim; k++) // C[i][j] = C[i][j] + A[i][k]*B[k][j];
         {
            ddf_mul(Ahi[i][k],Alo[i][k],Bhi[k][j],Blo[k][j],&acchi,&acclo);
            ddf_inc(&Chi[i][j],&Clo[i][j],acchi,acclo);
         }
      }
}

void CPU_dbl2_factors_forward
 ( int dim, double **Lhi, double **Llo, double *bhi, double *blo,
   double *xhi, double *xlo )
{
   double acchi,acclo;

   xhi[0] = bhi[0];
   xlo[0] = blo[0];
   for(int i=1; i<dim; i++)
   {
      xhi[i] = bhi[i];
      xlo[i] = blo[i];
      for(int j=0; j<i; j++) // x[i] = x[i] - L[i][j]*x[j];
      {
         ddf_mul(Lhi[i][j],Llo[i][j],xhi[j],xlo[j],&acchi,&acclo);
         ddf_dec(&xhi[i],&xlo[i],acchi,acclo);
      }
   }
}

void CPU_dbl2_factors_backward
 ( int dim, double **Uhi, double **Ulo, double *bhi, double *blo,
   double *xhi, double *xlo )
{
   double acchi,acclo;

   for(int i=dim-1; i>=0; i--)
   {
      xhi[i] = bhi[i];
      xlo[i] = blo[i];
      for(int j=dim-1; j>i; j--) // x[i] = x[i] - U[i][j]*x[j];
      {
         ddf_mul(Uhi[i][j],Ulo[i][j],xhi[j],xlo[j],&acchi,&acclo);
         ddf_dec(&xhi[i],&xlo[i],acchi,acclo);
      }
      // x[i] = x[i]/U[i][i];
      ddf_div(xhi[i],xlo[i],Uhi[i][i],Ulo[i][i],&acchi,&acclo);
      xhi[i] = acchi; xlo[i] = acclo;
   }
}

void CPU_dbl2_factors_lufac
 ( int dim, double **Ahi, double **Alo, int *pivots )
{
   double valmax,valtmp,acchi,acclo;
   int idxmax,idxtmp;

   for(int j=0; j<dim; j++) pivots[j] = j;
   for(int j=0; j<dim; j++)
   {
      valmax = fabs(Ahi[j][j]); idxmax = j;     // find the pivot
      for(int i=j+1; i< dim; i++)
      {
         valtmp = fabs(Ahi[i][j]);
         if(valtmp > valmax)
         {
            valmax = valtmp; idxmax = i;
         }
      }
      if(idxmax != j)                        // swap rows
      {
         for(int k=0; k<dim; k++)
         {
            valtmp = Ahi[idxmax][k];
            Ahi[idxmax][k] = Ahi[j][k];
            Ahi[j][k] = valtmp;
            valtmp = Alo[idxmax][k];
            Alo[idxmax][k] = Alo[j][k];
            Alo[j][k] = valtmp;
         }
         idxtmp = pivots[idxmax];
         pivots[idxmax] = pivots[j];
         pivots[j] = idxtmp;
      }
      for(int i=j+1; i<dim; i++)             // reduce
      {
         // A[i][j] = A[i][j]/A[j][j];
         ddf_div(Ahi[i][j],Alo[i][j],Ahi[j][j],Alo[j][j],&acchi,&acclo);
         Ahi[i][j] = acchi;
         Alo[i][j] = acclo;
         for(int k=j+1; k<dim; k++) // A[i][k] = A[i][k] - A[i][j]*A[j][k];
         {
            ddf_mul(Ahi[i][j],Alo[i][j],Ahi[j][k],Alo[j][k],&acchi,&acclo);
            ddf_dec(&Ahi[i][k],&Alo[i][k],acchi,acclo);
         }
      }
   }
}

void CPU_dbl2_factors_lusolve
 ( int dim, double **Ahi, double **Alo, int *pivots,
   double *bhi, double *blo, double *xhi, double *xlo )
{
   CPU_dbl2_factors_lufac(dim,Ahi,Alo,pivots);
   for(int i=0; i<dim; i++) 
   {
      xhi[i] = bhi[pivots[i]];
      xlo[i] = blo[pivots[i]];
   }
   CPU_dbl2_factors_forward(dim,Ahi,Alo,xhi,xlo,bhi,blo);
   CPU_dbl2_factors_backward(dim,Ahi,Alo,bhi,blo,xhi,xlo);
}
