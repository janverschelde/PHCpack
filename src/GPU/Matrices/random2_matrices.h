// The file random2_matrices.h specifies functions define tests to work
// with matrices and vectors in double double precision.

#ifndef __random2_matrices_h__
#define __random2_matrices_h__

void random_dbl2_upper_matrix
 ( int rows, int cols, double **Ahi, double **Alo );
/*
 * DESCRIPTION :
 *   Returns in (Ahi, Alo) a real upper triangular matrix
 *   of dimensions rows and cols,
 *   of randomly generated double doubles.
 *
 * ON ENTRY :
 *   rows     the number of rows in Ahi and Alo;
 *   cols     the number of columns in Ahi and Alo;
 *   Ahi      space for a matrix of doubles, of dimensions rows and cols.
 *   Alo      space for a matrix of doubles, of dimensions rows and cols.
 *
 * ON RETURN :
 *   Ahi      the high doubles of an upper triangular matrix
 *            of randomly generated doubles;
 *   Alo      the low doubles of an upper triangular matrix
 *            of randomly generated doubles. */

void random_dbl2_matrix
 ( int rows, int cols, double **Ahi, double **Alo );
/*
 * DESCRIPTION :
 *   Returns in (Ahi, Alo) a real matrix of dimensions rows and cols,
 *   of randomly generated double doubles.
 *
 * ON ENTRY :
 *   rows     the number of rows in Ahi and Alo;
 *   cols     the number of columns in Ahi and Alo;
 *   Ahi      space for a matrix of doubles, of dimensions rows and cols.
 *   Alo      space for a matrix of doubles, of dimensions rows and cols.
 *
 * ON RETURN :
 *   Ahi      the high doubles of a matrix of randomly generated doubles;
 *   Alo      the low doubles of a matrix of randomly generated doubles. */

#endif
