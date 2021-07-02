// The file dbl_tabs_testers.h specifies test functions on
// tiled accelerated back substitution in double precision.

#ifndef __dbl_tabs_testers_h__
#define __dbl_tabs_testers_h__

double dbl_Difference_Sum ( int n, double *x, double *y );
/*
 * DESCRIPTION :
 *   Returns the sum of the absolute value of the differences
 *   between the entries of the n-dimensional vectors x and y. */

double dbl_Column_Sum ( int dim, int col, double **A );
/*
 * DESCRIPTION :
 *   Returns the sum of the absolute values of all dim elements
 *   of the matrix A in the column with index col. */

double dbl_Max_Column_Sum ( int dim, double **A );
/*
 * DESCRIPTION :
 *   Returns the maximum column sum of the elements
 *   in the matrix A of dimension dim. */

double dbl_condition ( int dim, double **A, double **invA );
/*
 * DESCRIPTION :
 *   Returns the condition number of A using the 1-norm
 *   on the matrix of dimension dim and its inverse invA. */

void test_real_upper_inverse ( void );
/*
 * DESCRIPTION :
 *   Generates a random real upper triangular matrix
 *   to test the computation of its inverse. */

void test_real_upper_tiling ( void );
/*
 * DESCRIPTION :
 *   Prompts for the size of each tile and the number of tiles
 *   and applies the tiled back substitution to a random system. */

#endif
