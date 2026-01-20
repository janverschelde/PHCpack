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

void single_indexed_matrix_multiplication
 ( int nrows, int ncols, int dim, double *A, double *B, double *C )
{
   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         double temp = 0.0;

         for(int k=0; k<dim; k++) temp += A[i*dim+k] * B[j*dim+k];

         C[i*ncols+j] += temp;
      }
}

void transpose_rows_columns
 ( int nrows, int ncols, double **A, double **T )
{
   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++) T[j][i] = A[i][j];
}

void double2single_row_major
 ( int nrows, int ncols, double **A, double *B )
{
   for(int i=0, k=0; i<nrows; i++)
      for(int j=0; j<ncols; j++) B[k++] = A[i][j];
}

void double2single_column_major
 ( int nrows, int ncols, double **A, double *B )
{
   for(int i=0, k=0; i<ncols; i++)
      for(int j=0; j<nrows; j++) B[k++] = A[i][j];
}

void single2double_row_major
 ( int nrows, int ncols, double *A, double **B )
{
   for(int i=0, k=0; i<nrows; i++)
      for(int j=0; j<ncols; j++) B[i][j] = A[k++];
}

void double2single_column_major
 ( int nrows, int ncols, double *A, double **B )
{
   for(int i=0, k=0; i<ncols; i++)
      for(int j=0; j<nrows; j++) B[i][j] = A[k++];
}
