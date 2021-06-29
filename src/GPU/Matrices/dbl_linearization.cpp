// The file dbl_linearization.cpp defines the functions
// specified in dbl_linearization.h.

#include "dbl_linearization.h"

void dbl_linear_series_vector
 ( int dim, int deg, double **v, double **w )
{
   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++) w[j][i] = v[i][j];
}

void cmplx_linear_series_vector
 ( int dim, int deg, double **vre, double **vim,
   double **wre, double **wim )
{
   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++)
      {
         wre[j][i] = vre[i][j];
         wim[j][i] = vim[i][j];
      }
}

void dbl_linear_series_matrix
 ( int nrows, int ncols, int deg, double ***A, double ***B )
{
   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
         for(int k=0; k<=deg; k++) B[k][i][j] = A[i][j][k];
}

void cmplx_linear_series_matrix
 ( int nrows, int ncols, int deg, double ***Are, double ***Aim,
   double ***Bre, double ***Bim )
{
   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
         for(int k=0; k<=deg; k++)
         {
            Bre[k][i][j] = Are[i][j][k];
            Bim[k][i][j] = Aim[i][j][k];
         }
}
