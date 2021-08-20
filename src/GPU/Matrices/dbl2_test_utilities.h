// The file dbl2_test_utilities.h specifies common test functions on
// vectors and matrices in double double precision.

#ifndef __dbl2_test_utilities_h__
#define __dbl2_test_utilities_h__

double dbl2_Difference_Sum
 ( int n, double *xhi, double *xlo, double *yhi, double *ylo );
/*
 * DESCRIPTION :
 *   Returns the sum of the absolute value of the differences
 *   between the entries of the n-dimensional vectors x and y. */

double cmplx2_Difference_Sum
 ( int n, double *xrehi, double *xrelo, double *ximhi, double *ximlo,
          double *yrehi, double *yrelo, double *yimhi, double *yimlo );
/*
 * DESCRIPTION :
 *   Returns the sum of the absolute value of the differences
 *   between the entries of the n-dimensional vectors x and y. */

double dbl2_Column_Sum ( int dim, int col, double **Ahi, double **Alo );
/*
 * DESCRIPTION :
 *   Returns the sum of the absolute values of all dim elements
 *   of the matrix A in the column with index col. */

double cmplx2_Column_Sum
 ( int dim, int col,
   double **Arehi, double **Arelo, double **Aimhi, double **Aimlo );
/*
 * DESCRIPTION :
 *   Returns the sum of the absolute values of all dim elements
 *   of the matrix A in the column with index col. */

double dbl2_Max_Column_Sum ( int dim, double **Ahi, double **Alo );
/*
 * DESCRIPTION :
 *   Returns the maximum column sum of the elements
 *   in the matrix A of dimension dim. */

double cmplx2_Max_Column_Sum
 ( int dim, double **Arehi, double **Arelo, double **Aimhi, double **Aimlo );
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

double cmplx2_condition
 ( int dim, double **Arehi, double **Arelo, double **Aimhi, double **Aimlo,
   double **invArehi, double **invArelo,
   double **invAimhi, double **invAimlo );
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

double cmplx2_Matrix_Difference_Sum
 ( int n, double **Arehi, double **Arelo, double **Aimhi, double **Aimlo,
   double **Brehi, double **Brelo, double **Bimhi, double **Bimlo );
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

double cmplx2_Diagonal_Difference_Sum
 ( int nbt, int szt,
   double **Arehi, double **Arelo, double **Aimhi, double **Aimlo,
   double **Brehi, double **Brelo, double **Bimhi, double **Bimlo );
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

void cmplx2_random_upper_factor
 ( int dim, double **Arehi, double **Arelo, double **Aimhi, double **Aimlo );
/*
 * DESCRIPTION :
 *   Returns the upper triangular factor of the LU factorization
 *   with row pivoting on a random matrix of dimension dim.
 *   This yields a much better conditioned upper triangular matrix
 *   than the direct generation of a random upper triangular matrix. */

#endif
