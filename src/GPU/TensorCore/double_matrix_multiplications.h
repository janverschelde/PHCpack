/* Defines functions for matrix matrix multiplication of doubles. */

#ifndef __DOUBLE_MATRIX_MULTIPLICATION_H__
#define __DOUBLE_MATRIX_MULTIPLICATION_H__

void double_indexed_matrix_multiplication
 ( int nrows, int ncols, int dim, double **A, double **B, double **C );
/*
 * Given in A is an nrows-by-dim matrix and
 * given in B is an dim-by-ncols matrix,
 * returns in C the product of A with B,
 * using the common double indexed notation,
 * working with matrices as arrays of rows. */

#endif
