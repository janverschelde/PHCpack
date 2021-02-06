// The file random_matrices.cpp defines the functions with prototypes
// in random_matrices.h.

#include "random_series.h"
#include "random_matrices.h"

void random_dbl_series_vectors
 ( int dim, int deg, double *x, double **plux, double **minx )
{
   for(int k=0; k<dim; k++)
      random_dbl_exponentials(deg,&x[k],plux[k],minx[k]);
}
