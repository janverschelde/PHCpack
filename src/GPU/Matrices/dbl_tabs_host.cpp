/* The file dbl_tabs_host.cpp defines the functions specified in
 * the file dbl_tabs_host.h. */

#include <cstdlib>
#include "dbl_matrices_host.h"
#include "dbl_tabs_host.h"

void CPU_dbl_backsubs ( int dim, double **U, double *b, double *x )
{
   CPU_dbl_upper_lead_solver(dim,U,b,x);
}

void CPU_dbl_upper_inverse ( int dim, double **U, double **invU )
{
   double *col = new double[dim];
   double *rhs = new double[dim];

   for(int i=0; i<dim; i++) rhs[i] = 0.0;

   for(int j=0; j<dim; j++) // compute j-th column of the inverse
   {
      rhs[j] = 1.0;
      CPU_dbl_backsubs(dim,U,rhs,col);
      for(int i=0; i<dim; i++) invU[i][j] = col[i];
      rhs[j] = 0.0;
   }
   free(rhs); free(col);
}
