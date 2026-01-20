/* Collection of functions for matrix matrix multiplication of doubles. */

void double_indexed_matrix_multiplication
 ( int nrows, int ncols, int dim, double **A, double **B, double **C )
{
   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         C[i][j] = 0.0;

         for(int k=0; k<dim; k++) C[i][j] += A[i][k]*B[k][j];
      }
}
