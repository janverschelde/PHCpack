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

void single_indexed_matrix_multiplication
 ( int nrows, int ncols, int dim, double *A, double *B, double *C );
/*
 * Multiplies the matrices A and B into C on the host.
 * All matrices are stored in one dimensional arrays.
 * Matrices A and C are row major, B is column major.
 * A is nrows-by-dim, B is dim-by-ncols, C is nrows-by-ncols. */

void transpose_rows_columns
 ( int nrows, int ncols, double **A, double **T );
/*
 * Given in A is an nrows-by-ncols matrix,
 * swaps rows with columns into the ncols-by-nrows matrix T
 * so T is the column major version of A. */

void double2single_row_major
 ( int nrows, int ncols, double **A, double *B );
/*
 * Maps the nrows-by-ncols matrix A,
 * given as an array of arrays, of nrows arrays of size ncols,
 * into the single indexed array B,
 * which has space for nrows*ncols doubles. */

void double2single_column_major
 ( int nrows, int ncols, double **A, double *B );
/*
 * Maps the ncols-by-nrows matrix A,
 * given as an array of arrays, of ncols arrays of size nrows,
 * into the single indexed array B,
 * which has space for nrows*ncols doubles. */

void single2double_row_major
 ( int nrows, int ncols, double *A, double **B );
/*
 * Maps the single indexed array A of size nrows*ncols
 * into the double indexed nrows-by-ncols matrix B,
 * stored as an array of arrays, of nrows arrays of size ncols. */

void double2single_column_major
 ( int nrows, int ncols, double *A, double **B );
/*
 * Maps the single indexed array A of size nrows*ncols,
 * into the double indexed ncols-by-nrows matrix B,
 * stored as an array of arrays, of ncols arrows of size nrows. */

#endif
