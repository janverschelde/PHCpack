// The file random16_matrices.h specifies functions define tests to work
// with matrices and vectors in hexa double precision.

#ifndef __random16_matrices_h__
#define __random16_matrices_h__

void random_dbl16_matrix
 ( int rows, int cols,
   double **Ahihihihi, double **Alohihihi,
   double **Ahilohihi, double **Alolohihi,
   double **Ahihihilo, double **Alohihilo,
   double **Ahilohilo, double **Alolohilo,
   double **Ahihilohi, double **Alohilohi,
   double **Ahilolohi, double **Alololohi,
   double **Ahihilolo, double **Alohilolo,
   double **Ahilololo, double **Alolololo );
/*
 * DESCRIPTION :
 *   Returns sixteen matrices of dimensions rows and cols,
 *   to represent a real matrix of hexa doubles.
 *
 * ON ENTRY :
 *   rows        the number of rows of all matrices;
 *   cols        the number of columns of all matrices;
 *   Ahihihihi   has space for a matrix, of dimensions rows and cols;
 *   Alohihihi   has space for a matrix, of dimensions rows and cols;
 *   Ahilohihi   has space for a matrix, of dimensions rows and cols;
 *   Alolohihi   has space for a matrix, of dimensions rows and cols;
 *   Ahihilohi   has space for a matrix, of dimensions rows and cols;
 *   Alohilohi   has space for a matrix, of dimensions rows and cols;
 *   Ahilolohi   has space for a matrix, of dimensions rows and cols;
 *   Alololohi   has space for a matrix, of dimensions rows and cols;
 *   Ahihihilo   has space for a matrix, of dimensions rows and cols;
 *   Alohihilo   has space for a matrix, of dimensions rows and cols;
 *   Ahilohilo   has space for a matrix, of dimensions rows and cols;
 *   Alolohilo   has space for a matrix, of dimensions rows and cols;
 *   Ahihilolo   has space for a matrix, of dimensions rows and cols;
 *   Alohilolo   has space for a matrix, of dimensions rows and cols;
 *   Ahilololo   has space for a matrix, of dimensions rows and cols;
 *   Alolololo   has space for a matrix, of dimensions rows and cols.
 *
 * ON RETURN :
 *   Ahihihihi   are the highest doubles of a matrix A
 *               of randomly generated octo doubles;
 *   Alohihihi   are the second highest doubles of A;
 *   Ahilohihi   are the third highest doubles of A;
 *   Alolohihi   are the fourth highest doubles of A;
 *   Ahihilohi   are the fifth highest doubles of A;
 *   Alohilohi   are the sixth highest doubles of A;
 *   Ahilolohi   are the seventh highest doubles of A;
 *   Alololohi   are the eighth highest doubles of A;
 *   Ahihihilo   are the eighth lowest doubles of A;
 *   Alohihilo   are the seventh lowest doubles of A;
 *   Ahilohilo   are the sixth lowest doubles of A;
 *   Alolohilo   are the fifth lowest doubles of A;
 *   Ahihilolo   are the fourth lowest doubles of A;
 *   Alohilolo   are the third lowest doubles of A;
 *   Ahilololo   are the second lowest doubles of A;
 *   Alolololo   are the lowest doubles of A. */

#endif
