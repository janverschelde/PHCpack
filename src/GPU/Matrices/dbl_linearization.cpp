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
