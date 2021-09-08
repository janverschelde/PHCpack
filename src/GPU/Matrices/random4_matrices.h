// The file random4_matrices.h specifies functions define tests to work
// with matrices and vectors in quad double precision.

#ifndef __random4_matrices_h__
#define __random4_matrices_h__

void random_dbl4_upper_matrix
 ( int rows, int cols,
   double **Ahihi, double **Alohi, double **Ahilo, double **Alolo );
/*
 * DESCRIPTION :
 *   Returns in (Ahihi, Alohi, Ahilo, Alolo) a real upper triangular matrix
 *   of dimensions rows and cols, of randomly generated quad doubles.
 *
 * ON ENTRY :
 *   rows     the number of rows of all matrices;
 *   cols     the number of columns of all matrices;
 *   Ahihi    space for a matrix, of dimensions rows and cols;
 *   Alohi    space for a matrix, of dimensions rows and cols;
 *   Ahilo    space for a matrix, of dimensions rows and cols;
 *   Alolo    space for a matrix, of dimensions rows and cols.
 *
 * ON RETURN :
 *   Ahihi    the highest doubles of an upper triangular matrix
 *            of randomly generated quad doubles;
 *   Alohi    the second highest doubles of an upper triangular matrix
 *            of randomly generated quad doubles;
 *   Ahilo    the second lowes doubles of an upper triangular matrix
 *            of randomly generated quad doubles;
 *   Alolo    thees low doubles of an upper triangular matrix
 *            of randomly generated quad doubles. */

void random_cmplx4_upper_matrix
 ( int rows, int cols,
   double **Arehihi, double **Arelohi, double **Arehilo, double **Arelolo,
   double **Aimhihi, double **Aimlohi, double **Aimhilo, double **Aimlolo );
/*
 * DESCRIPTION :
 *   Returns in (Arehihi, Arelohi, Arehilo, Arelolo, Aimhihi, Aimlohi,
 *   Aimhilo, Aimlolo) a complex upper triangular matrix of dimensions
 *   rows and cols, of randomly generated quad doubles.
 *
 * ON ENTRY :
 *   rows     the number of rows of all matrices;
 *   cols     the number of columns of all matrices;
 *   Arehihi  space for a matrix, of dimensions rows and cols;
 *   Arelohi  space for a matrix, of dimensions rows and cols;
 *   Arehilo  space for a matrix, of dimensions rows and cols;
 *   Arelolo  space for a matrix, of dimensions rows and cols;
 *   Aimhihi  space for a matrix, of dimensions rows and cols;
 *   Aimlohi  space for a matrix, of dimensions rows and cols;
 *   Aimhilo  space for a matrix, of dimensions rows and cols;
 *   Aimlolo  space for a matrix, of dimensions rows and cols.
 *
 * ON RETURN :
 *   Arehihi  highest doubles of the real parts of a random matrix;
 *   Arelohi  second highest doubles of the real parts of a random matrix;
 *   Arehilo  second low doubles of the real parts of a random matrix;
 *   Arelolo  lowest doubles of the real parts of a random matrix;
 *   Aimhihi  highest doubles of the imaginary parts of a random matrix;
 *   Aimlohi  second highest doubles of the imaginary parts;
 *   Aimhilo  second lowest doubles of the imaginary parts;
 *   Aimlolo  lowest doubles of the imaginary parts of a random matrix. */

void random_dbl4_matrix
 ( int rows, int cols,
   double **Ahihi, double **Alohi, double **Ahilo, double **Alolo );
/*
 * DESCRIPTION :
 *   Returns in (Ahihi, Alohi, Ahilo, Alolo) a real matrix
 *   of dimensions rows and cols, of randomly generated quad doubles.
 *
 * ON ENTRY :
 *   rows     the number of rows of all matrices;
 *   cols     the number of columns of all matrices;
 *   Ahihi    space for a matrix, of dimensions rows and cols;
 *   Alohi    space for a matrix, of dimensions rows and cols;
 *   Ahilo    space for a matrix, of dimensions rows and cols;
 *   Alolo    space for a matrix, of dimensions rows and cols.
 *
 * ON RETURN :
 *   Ahihi    highest doubles of a random matrix;
 *   Alohi    second highest doubles of a random matrix;
 *   Alohi    second lowest doubles of a random matrix;
 *   Alolo    lowest doubles of a random matrix. */

void random_cmplx4_matrix
 ( int rows, int cols,
   double **Arehihi, double **Arelohi, double **Arehilo, double **Arelolo,
   double **Aimhihi, double **Aimlohi, double **Aimhilo, double **Aimlolo );
/*
 * DESCRIPTION :
 *   Returns in (Arehihi, Arelohi, Arehilo, Arelolo, Aimhihi, Aimlohi,
 *   Aimhilo, Aimlolo) a complex matrix of dimensions rows and cols,
 *   of randomly generated quad doubles.
 *
 * ON ENTRY :
 *   rows     the number of rows of all matrices;
 *   cols     the number of columns of all matrices;
 *   Arehihi  space for a matrix, of dimensions rows and cols;
 *   Arelohi  space for a matrix, of dimensions rows and cols;
 *   Arehilo  space for a matrix, of dimensions rows and cols;
 *   Arelolo  space for a matrix, of dimensions rows and cols;
 *   Aimhihi  space for a matrix, of dimensions rows and cols;
 *   Aimlohi  space for a matrix, of dimensions rows and cols;
 *   Aimhilo  space for a matrix, of dimensions rows and cols;
 *   Aimlolo  space for a matrix, of dimensions rows and cols.
 *
 * ON RETURN :
 *   Arehihi  highest doubles of the real parts of a random matrix;
 *   Arelohi  second highest doubles of the real parts of a random matrix;
 *   Arehilo  second low doubles of the real parts of a random matrix;
 *   Arelolo  lowest doubles of the real parts of a random matrix;
 *   Aimhihi  highest doubles of the imaginary parts of a random matrix;
 *   Aimlohi  second highest doubles of the imaginary parts;
 *   Aimhilo  second lowest doubles of the imaginary parts;
 *   Aimlolo  lowest doubles of the imaginary parts of a random matrix. */

#endif
