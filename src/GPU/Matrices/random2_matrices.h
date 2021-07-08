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
 *   Ahi      space for a matrix, of dimensions rows and cols.
 *   Alo      space for a matrix, of dimensions rows and cols.
 *
 * ON RETURN :
 *   Ahi      the high doubles of an upper triangular matrix
 *            of randomly generated doubles;
 *   Alo      the low doubles of an upper triangular matrix
 *            of randomly generated doubles. */

void random_cmplx2_upper_matrix
 ( int rows, int cols,
   double **Arehi, double **Arelo, double **Aimhi, double **Aimlo );
/*
 * DESCRIPTION :
 *   Returns in (Arehi, Arelo, Aimhi, Aimlo) a complex upper triangular
 *   matrix of dimensions rows and cols,
 *   of randomly generated double doubles.
 *
 * ON ENTRY :
 *   rows     the number of rows in Arehi, Arelo, Aimhi, and Aimlo;
 *   cols     the number of columns in Arehi, Arelo, Aimhi, and Aimlo;
 *   Arehi    space for a matrix, of dimensions rows and cols;
 *   Arelo    space for a matrix, of dimensions rows and cols;
 *   Aimhi    space for a matrix, of dimensions rows and cols;
 *   Aimlo    space for a matrix, of dimensions rows and cols.
 *
 * ON RETURN :
 *   Arehi    high doubles of the real parts of a random matrix;
 *   Arelo    low doubles of the real parts of a random matrix;
 *   Aimhi    high doubles of the imaginary parts of a random matrix;
 *   Aimlo    low doubles of the imaginary parts of a random matrix. */

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
 *   Ahi      space for a matrix, of dimensions rows and cols;
 *   Alo      space for a matrix, of dimensions rows and cols.
 *
 * ON RETURN :
 *   Ahi      the high doubles of a matrix of randomly generated doubles;
 *   Alo      the low doubles of a matrix of randomly generated doubles. */

void random_cmplx2_matrix
 ( int rows, int cols,
   double **Arehi, double **Arelo, double **Aimhi, double **Aimlo );
/*
 * DESCRIPTION :
 *   Returns in (Arehi, Arelo, Aimhi, Aimlo) a complex matrix 
 *   of dimensions rows and cols, of randomly generated double doubles.
 *
 * ON ENTRY :
 *   rows     the number of rows in Arehi, Arelo, Aimhi, and Aimlo;
 *   cols     the number of columns in Arehi, Arelo, Aimhi, and Aimlo;
 *   Arehi    space for a matrix, of dimensions rows and cols;
 *   Arelo    space for a matrix, of dimensions rows and cols;
 *   Aimhi    space for a matrix, of dimensions rows and cols;
 *   Aimlo    space for a matrix, of dimensions rows and cols.
 *
 * ON RETURN :
 *   Arehi    high doubles of the real parts of a random matrix;
 *   Arelo    low doubles of the real parts of a random matrix;
 *   Aimhi    high doubles of the imaginary parts of a random matrix;
 *   Aimlo    low doubles of the imaginary parts of a random matrix. */

#endif
