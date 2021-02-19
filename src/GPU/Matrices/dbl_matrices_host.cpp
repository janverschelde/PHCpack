/* The file dbl_matrices_host.cpp defines the function specified
 * in dbl_matrices_host.h. */

#include <cstdlib>
#include <cmath>
#include "dbl_convolutions_host.h"

void real_inner_product
 ( int dim, int deg, double **x, double **y, double *z )
{
   double *prod = new double[deg+1];

   for(int i=0; i<=deg; i++) z[i] = 0.0;

   for(int k=0; k<dim; k++)
   {
      CPU_dbl_product(deg,x[k],y[k],prod);
      for(int i=0; i<=deg; i++)
         z[i] = z[i] + prod[i];
   }
   free(prod);
}

void cmplx_inner_product
 ( int dim, int deg,
   double **xre, double **xim, double **yre, double **yim,
   double *zre, double *zim )
{
   double *prodre = new double[deg+1];
   double *prodim = new double[deg+1];

   for(int i=0; i<=deg; i++)
   {
      zre[i] = 0.0; zim[i] = 0.0;
   }
   for(int k=0; k<dim; k++)
   {
      CPU_cmplx_product(deg,xre[k],xim[k],yre[k],yim[k],prodre,prodim);
      for(int i=0; i<=deg; i++)
      {
         zre[i] = zre[i] + prodre[i];
         zim[i] = zim[i] + prodim[i];
      }
   }
   free(prodre); free(prodim);
}

void real_matrix_vector_product
 ( int rows, int cols, int deg, double ***A, double **x, double **y )
{
   for(int k=0; k<rows; k++)
      real_inner_product(cols,deg,A[k],x,y[k]);
}

void cmplx_matrix_vector_product
 ( int rows, int cols, int deg, double ***Are, double ***Aim,
   double **xre, double **xim, double **yre, double **yim )
{
   for(int k=0; k<rows; k++)
      cmplx_inner_product(cols,deg,Are[k],Aim[k],xre,xim,yre[k],yim[k]);
}

void real_upper_solver
 ( int dim, int deg, double ***U, double **b, double **x  )
{
   double *prod = new double[deg+1];
   double *work = new double[deg+1];

   for(int i=dim-1; i>=0; i--)
   {
      for(int k=0; k<=deg; k++) prod[k] = b[i][k];

      for(int j=i+1; j<dim; j++)
      {
         CPU_dbl_product(deg,U[i][j],x[j],work);

         for(int k=0; k<=deg; k++) prod[k] = prod[k] - work[k];
      }
      CPU_dbl_inverse(deg,U[i][i],work);
      CPU_dbl_product(deg,work,prod,x[i]);
   }
   free(prod); free(work);
}

void cmplx_upper_solver
 ( int dim, int deg, double ***Ure, double ***Uim,
   double **bre, double **bim, double **xre, double **xim )
{
   double *prodre = new double[deg+1];
   double *prodim = new double[deg+1];
   double *workre = new double[deg+1];
   double *workim = new double[deg+1];

   for(int i=dim-1; i>=0; i--)
   {
      for(int k=0; k<=deg; k++)
      {
         prodre[k] = bre[i][k];
         prodim[k] = bim[i][k];
      }
      for(int j=i+1; j<dim; j++)
      {
         CPU_cmplx_product
            (deg,Ure[i][j],Uim[i][j],xre[j],xim[j],workre,workim);

         for(int k=0; k<=deg; k++)
         {
            prodre[k] = prodre[k] - workre[k];
            prodim[k] = prodim[k] - workim[k];
         }
      }
      CPU_cmplx_inverse(deg,Ure[i][i],Uim[i][i],workre,workim);
      CPU_cmplx_product(deg,workre,workim,prodre,prodim,xre[i],xim[i]);
   }
   free(prodre); free(workre);
   free(prodim); free(workim);
}

void real_lufac ( int dim, int deg, double ***A, int *pivots )
{
   double valmax,valtmp;
   int idxmax,idxtmp;
   double *work = new double[deg+1];
   double *prod = new double[deg+1];

   for(int j=0; j<dim; j++)
   {
      valmax = fabs(A[j][j][0]); idxmax = j;
      for(int i=j+1; i<dim; i++)
      {                                // find the pivot
         valtmp = fabs(A[i][j][0]);
         if(valtmp > valmax) 
         {
            valmax = valtmp; idxmax = i;
         }
      }
      if(idxmax != j)                  // swap rows
      {
         for(int k=0; k<dim; k++)
         {
            for(int i=0; i<=deg; i++)
            {
               valtmp = A[idxmax][k][i];
               A[idxmax][k][i] = A[j][k][i];
               A[j][k][i] = valtmp;
            }
         }
         idxtmp = pivots[idxmax];
         pivots[idxmax] = pivots[j];
         pivots[j] = idxtmp;
      }
      for(int i=j+1; i<dim; i++)
      {                                        // A[i][j] = A[i][j]/A[j][j]
         CPU_dbl_inverse(deg,A[j][j],work);  
         CPU_dbl_product(deg,A[i][j],work,prod);
         for(int k=0; k<=deg; k++) A[i][j][k] = prod[k];

         for(int k=j+1; k<dim; k++)  // A[i][k] = A[i][k] - A[i][j]*A[j][k]
         {
            CPU_dbl_product(deg,A[i][j],A[j][k],prod);
            for(int d=0; d<=deg; d++) A[i][k][d] = A[i][k][d] - prod[d];
         }
      }
   }
}
