/* The file dbl_matrices_host.cpp defines the function specified
 * in dbl_matrices_host.h. */

#include <cstdlib>
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
