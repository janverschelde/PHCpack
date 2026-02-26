/* Definitions of functions for double double matrix matrix multiplication. */

#include <iostream>
#include "double_double_functions.h"
#include "random2_matrices.h"

using namespace std;

void random_dd_matrices
 ( int nrows, int ncols, int dim,
   double **Ahi, double **Alo, double **Bhi, double **Blo,
   double **Chi, double **Clo, int vrblvl )
{
   if(vrblvl > 0)
      cout << "-> in random_dd_matrices, nrows : " << nrows
           << ", ncols : " << ncols << ", dim : " << dim << endl;

   for(int i=0; i<nrows; i++)
   {
      Chi[i] = new double[ncols];
      Clo[i] = new double[ncols];
      Ahi[i] = new double[dim];
      Alo[i] = new double[dim];
   }
   random_dbl2_matrix(nrows, dim, Ahi, Alo);

   for(int i=0; i<dim; i++)
   {
      Bhi[i] = new double[ncols];
      Blo[i] = new double[ncols];
   }
   random_dbl2_matrix(dim, ncols, Bhi, Blo);

   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         Chi[i][j] = 0.0;
         Clo[i][j] = 0.0;
      }
}

void double_double_matmatmul
 ( int nrows, int ncols, int dim,
   double **Ahi, double **Alo, double **Bhi, double **Blo,
   double **Chi, double **Clo )
{
   double acchi,acclo;

   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         Chi[i][j] = 0.0; Clo[i][j] = 0.0;

         for(int k=0; k<dim; k++)
         {
            ddf_mul(Ahi[i][k], Alo[i][k], Bhi[k][j], Blo[k][j],
                    &acchi, &acclo);
            ddf_inc(&Chi[i][j], &Clo[i][j], acchi, acclo);
         }
      }
}

void transpose_dd_matrix
 ( int nrows, int ncols,
   double **Ahi, double **Alo, double **Thi, double **Tlo )
{
   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         Thi[j][i] = Ahi[i][j];
         Tlo[j][i] = Alo[i][j];
      }
}

void double_double_transposed_mm
 ( int nrows, int ncols, int dim,
   double **Ahi, double **Alo, double **BThi, double **BTlo,
   double **Chi, double **Clo )
{
   double acchi,acclo;

   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         Chi[i][j] = 0.0; Clo[i][j] = 0.0;

         for(int k=0; k<dim; k++)
         {
            ddf_mul(Ahi[i][k], Alo[i][k], BThi[j][k], BTlo[j][k],
                    &acchi, &acclo);
            ddf_inc(&Chi[i][j], &Clo[i][j], acchi, acclo);
         }
      }
}

void flopcount_dd_matmatmul
 ( int nrows, int ncols, int dim,
   long long int *add, long long int *mul )
{
   const int nbops = nrows*ncols*dim;

   *add = 20*nbops; // one addition requires 20 double operations
   *mul = 23*nbops; // one multiplication requires 23 double operations
}
