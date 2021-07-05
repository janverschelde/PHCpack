// The file random2_matrices.cpp defines the functions specified in
// the file random2_matrices.h.

#include "random2_vectors.h"
#include "random2_matrices.h"

void random_dbl2_upper_matrix
 ( int rows, int cols, double **Ahi, double **Alo )
{
   for(int i=0; i<rows; i++)
   {
      for(int j=0; j<i; j++)
      {
         Ahi[i][j] = 0.0;
         Alo[i][j] = 0.0;
      }
      for(int j=i; j<cols; j++)
         random_double_double(&Ahi[i][j],&Alo[i][j]);
   }
}

void random_dbl2_matrix
 ( int rows, int cols, double **Ahi, double **Alo )
{
   for(int i=0; i<rows; i++)
   {
      for(int j=0; j<cols; j++)
         random_double_double(&Ahi[i][j],&Alo[i][j]);
   }
}
