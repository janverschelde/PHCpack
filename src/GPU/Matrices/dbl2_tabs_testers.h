// The file dbl2_tabs_testers.h specifies test functions on the
// tiled accelerated back substitution in double double precision.

#ifndef __dbl2_tabs_testers_h__
#define __dbl2_tabs_testers_h__

double dbl2_Difference_Sum
 ( int n, double *xhi, double *xlo, double *yhi, double *ylo );
/*
 * DESCRIPTION :
 *   Returns the sum of the absolute value of the differences
 *   between the entries of the n-dimensional vectors x and y. */

double dbl2_Column_Sum ( int dim, int col, double **Ahi, double **Alo );
/*
 * DESCRIPTION :
 *   Returns the sum of the absolute values of all dim elements
 *   of the matrix A in the column with index col. */

double dbl2_Max_Column_Sum ( int dim, double **Ahi, double **Alo );
/*
 * DESCRIPTION :
 *   Returns the maximum column sum of the elements
 *   in the matrix A of dimension dim. */

double dbl2_condition
 ( int dim, double **Ahi, double **Alo, double **invAhi, double **invAlo );
/*
 * DESCRIPTION :
 *   Returns the condition number of A using the 1-norm
 *   on the matrix of dimension dim and its inverse invA. */

double dbl2_Matrix_Difference_Sum
 ( int n, double **Ahi, double **Alo, double **Bhi, double **Blo );
/*
 * DESCRIPTION :
 *   Returns the sum of the absolute value of the differences
 *   between two n-by-n matrices A and B. */

double dbl2_Diagonal_Difference_Sum
 ( int nbt, int szt, double **Ahi, double **Alo,
   double **Bhi, double **Blo );
/*
 * DESCRIPTION :
 *   Returns the sum of the absolute value of all differences 
 *   of the diagonal tiles of the matrices A and B.
 *   The number of tiles equals nbt and the size of each tile is szt. */

void dbl2_random_upper_factor ( int dim, double **Ahi, double **Alo );
/*
 * DESCRIPTION :
 *   Returns the upper triangular factor of the LU factorization
 *   with row pivoting on a random matrix of dimension dim.
 *   This yields a much better conditioned upper triangular matrix
 *   than the direct generation of a random upper triangular matrix. */

void test_real2_upper_inverse ( void );
/*
 * DESCRIPTION :
 *   Generates a random real upper triangular matrix
 *   to test the computation of its inverse. */

void test_real2_upper_tiling ( void );
/*
 * DESCRIPTION :
 *   Prompts for the size of each tile and the number of tiles
 *   and applies the tiled back substitution to a random system. */

#endif
