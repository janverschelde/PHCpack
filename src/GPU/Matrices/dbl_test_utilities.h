// The file dbl_test_utilities.h specifies common test functions on
// vectors and matrices in double precision.

#ifndef __dbl_test_utilities_h__
#define __dbl_test_utilities_h__

double dbl_Difference_Sum ( int n, double *x, double *y );
/*
 * DESCRIPTION :
 *   Returns the sum of the absolute value of the differences
 *   between the entries of the n-dimensional vectors x and y. */

double cmplx_Difference_Sum
 ( int n, double *xre, double *xim, double *yre, double *yim );
/*
 * DESCRIPTION :
 *   Returns the sum of the absolute value of the differences
 *   between the entries of the n-dimensional vectors x and y,
 *   given by separated real and imaginary vectors. */

double dbl_Column_Sum ( int dim, int col, double **A );
/*
 * DESCRIPTION :
 *   Returns the sum of the absolute values of all dim elements
 *   of the matrix A in the column with index col. */

double cmplx_Column_Sum ( int dim, int col, double **Are, double **Aim );
/*
 * DESCRIPTION :
 *   Returns the sum of the absolute values of all dim elements
 *   of the matrix A in the column with index col. */

double dbl_Max_Column_Sum ( int dim, double **A );
/*
 * DESCRIPTION :
 *   Returns the maximum column sum of the elements
 *   in the matrix A of dimension dim. */

double cmplx_Max_Column_Sum ( int dim, double **Are, double **Aim );
/*
 * DESCRIPTION :
 *   Returns the maximum column sum of the elements
 *   in the matrix A of dimension dim. */

double dbl_condition ( int dim, double **A, double **invA );
/*
 * DESCRIPTION :
 *   Returns the condition number of A using the 1-norm
 *   on the matrix of dimension dim and its inverse invA. */

double cmplx_condition
 ( int dim, double **Are, double **Aim, double **invAre, double **invAim );
/*
 * DESCRIPTION :
 *   Returns the condition number of A using the 1-norm
 *   on the matrix of dimension dim and its inverse invA. */

double dbl_Matrix_Difference_Sum ( int n, double **A, double **B );
/*
 * DESCRIPTION :
 *   Returns the sum of the absolute value of the differences
 *   between two n-by-n matrices A and B. */

double cmplx_Matrix_Difference_Sum
 ( int n, double **Are, double **Aim, double **Bre, double **Bim );
/*
 * DESCRIPTION :
 *   Returns the sum of the absolute value of the differences
 *   between two n-by-n matrices A and B. */

double dbl_Diagonal_Difference_Sum
 ( int nbt, int szt, double **A, double **B );
/*
 * DESCRIPTION :
 *   Returns the sum of the absolute value of all differences 
 *   of the diagonal tiles of the matrices A and B.
 *   The number of tiles equals nbt and the size of each tile is szt. */

double cmplx_Diagonal_Difference_Sum
 ( int nbt, int szt, double **Are, double **Aim,
   double **Bre, double **Bim );
/*
 * DESCRIPTION :
 *   Returns the sum of the absolute value of all differences 
 *   of the diagonal tiles of the matrices A and B.
 *   The number of tiles equals nbt and the size of each tile is szt. */

void dbl_random_upper_factor ( int dim, double **A );
/*
 * DESCRIPTION :
 *   Returns the upper triangular factor of the LU factorization
 *   with row pivoting on a random matrix of dimension dim.
 *   This yields a much better conditioned upper triangular matrix
 *   than the direct generation of a random upper triangular matrix. */

void cmplx_random_upper_factor ( int dim, double **Are, double **Aim );
/*
 * DESCRIPTION :
 *   Returns the upper triangular factor of the LU factorization
 *   with row pivoting on a random matrix of dimension dim.
 *   This yields a much better conditioned upper triangular matrix
 *   than the direct generation of a random upper triangular matrix. */

#endif
