/* Prototypes of functions for double double matrix matrix multiplication. */

#ifndef __DDMM_HOST_H__
#define __DDMM_HOST_H__

void random_dd_matrices
 ( int nrows, int ncols, int dim,
   double **Ahi, double **Alo, double **Bhi, double **Blo,
   double **Chi, double **Clo, int vrblvl=0 );
/*
 * Allocates and generates two random double double matrices.
 *
 * ON ENTRY :
 *   nrows    number of rows of A and the product C;
 *   ncols    number of columns of B and the product C;
 *   dim      number of columns of A and rows of B;
 *   Ahi      nrows rows allocated of the high doubles of the matrix A;
 *   Alo      nrows rows allocated of the low doubles of the matrix A;
 *   Bhi      dim rows allocated of the high doubles of the matrix B;
 *   Blo      dim rows allocated of the low doubles of the matrix B;
 *   Chi      dim rows allocated of the high doubles of the matrix C;
 *   Clo      dim rows allocated of the low doubles of the matrix C;
 *   vrblvl   is the verbose level.
 *
 * ON RETURN :
 *   Ahi      high doubles of an nrows-by-dim random double double matrix A;
 *   Alo      low doubles of an nrows-by-dim random double double matrix A;
 *   Bhi      high doubles of a dim-by-ncols random double double matrix B;
 *   Blo      low doubles of a dim-by-ncols random double double matrix B;
 *   Chi      allocated columns of high doubles of an nrows-by-ncols matrix;
 *   Clo      allocated columns of low doubles of an nrows-by-ncols matrix. */

void double_double_matmatmul
 ( int nrows, int ncols, int dim,
   double **Ahi, double **Alo, double **Bhi, double **Blo,
   double **Chi, double **Clo );
/*
 * Given in (Ahi, Alo) is an nrows-by-dim matrix A, and
 * given in (Bhi, Blo) is an dim-by-ncols matrix B,
 * returns in (Chi, Clo) the matrix matrix multiplication of A with B.
 * All matrices are row major which is not good for performance.  */

void transpose_dd_matrix
 ( int nrows, int ncols,
   double **Ahi, double **Alo, double **Thi, double **Tlo );
/*
 * Given an nrows-by-ncols double double matrix A (Ahi, Alo),
 * returns an ncols-by-nrows double double matrix T (Thi, Tlo),
 * which is the transpose of A. */

void double_double_transposed_mm
 ( int nrows, int ncols, int dim,
   double **Ahi, double **Alo, double **BThi, double **BTlo,
   double **Chi, double **Clo );
/*
 * Given in (Ahi, Alo) is an nrows-by-dim matrix A, and
 * given in (BThi, BTlo) is an ncols-by-dim matrix BT,
 * which stores the transpose of the matrix B,
 * returns in (Chi, Clo) the matrix matrix multiplication of A with B. */

void flopcount_dd_matmatmul
 ( int nrows, int ncols, int dim,
   long long int *add, long long int *mul );
/*
 * Returns in add the number of floating-point additions
 * and in mul the number of floating-point multiplications to
 * multiply an nrows-by-dim with a dim-by-ncols double double matrix. */

#endif
