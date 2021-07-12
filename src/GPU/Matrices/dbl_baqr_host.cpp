/* The file dbl_baqr_host.cpp defines the functions specified in
 * the file dbl_baqr_host.h. */

#include <cstdlib>
#include "dbl_factorizations.h"
#include "dbl_baqr_host.h"

void CPU_dbl_blocked_houseqr
 ( int nrows, int ncols, int szt, int nbt,
   double **A, double **Q, double **R )
{
   double *x = new double[nrows];
   double *v = new double[nrows];
   double beta;

   for(int i=0; i<nrows; i++)   // Q = I, R = A
   {
      for(int j=0; j<nrows; j++) Q[i][j] = 0.0;
      Q[i][i] = 1.0;
      for(int j=0; j<ncols; j++) R[i][j] = A[i][j];
   }
   for(int k=0; k<ncols; k++)
   {
      for(int i=k; i<nrows; i++) x[i-k] = R[i][k];
      CPU_dbl_factors_house(nrows-k,x,v,&beta);
   }
   free(x); free(v);
}
