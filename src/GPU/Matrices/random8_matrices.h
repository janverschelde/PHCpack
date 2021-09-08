// The file random8_matrices.h specifies functions define tests to work
// with matrices and vectors in octo double precision.

#ifndef __random8_matrices_h__
#define __random8_matrices_h__

void random_dbl8_upper_matrix
 ( int rows, int cols,
   double **Ahihihi, double **Alohihi, double **Ahilohi, double **Alolohi,
   double **Ahihilo, double **Alohilo, double **Ahilolo, double **Alololo );
/*
 * DESCRIPTION :
 *   Returns eight upper triangular matrices of dimensions rows and cols,
 *   to represent a real random upper triangular matrix of octo doubles.
 *
 * ON ENTRY :
 *   rows      the number of rows of all matrices;
 *   cols      the number of columns of all matrices;
 *   Ahihihi   has space for a matrix, of dimensions rows and cols;
 *   Alohihi   has space for a matrix, of dimensions rows and cols;
 *   Ahilohi   has space for a matrix, of dimensions rows and cols;
 *   Alolohi   has space for a matrix, of dimensions rows and cols;
 *   Ahihilo   has space for a matrix, of dimensions rows and cols;
 *   Alohilo   has space for a matrix, of dimensions rows and cols;
 *   Ahilolo   has space for a matrix, of dimensions rows and cols;
 *   Alololo   has space for a matrix, of dimensions rows and cols.
 *
 * ON RETURN :
 *   Ahihihi   are the highest doubles of an upper triangular matrix A
 *             of randomly generated octo doubles;
 *   Alohihi   are the second highest doubles of A;
 *   Ahilohi   are the third highest doubles of A;
 *   Alolohi   are the fourth highest doubles of A;
 *   Ahihilo   are the fourth lowest doubles of A;
 *   Alohilo   are the third lowest doubles of A;
 *   Ahilolo   are the second lowest doubles of A;
 *   Alololo   are the lowest doubles of A. */

void random_cmplx8_upper_matrix
 ( int rows, int cols,
   double **Arehihihi, double **Arelohihi,
   double **Arehilohi, double **Arelolohi,
   double **Arehihilo, double **Arelohilo,
   double **Arehilolo, double **Arelololo,
   double **Aimhihihi, double **Aimlohihi,
   double **Aimhilohi, double **Aimlolohi,
   double **Aimhihilo, double **Aimlohilo,
   double **Aimhilolo, double **Aimlololo );
/*
 * DESCRIPTION :
 *   Returns 16 upper triangular matrices of dimensions rows and cols,
 *   to represent a complex random upper triangular matrix of octo doubles.
 *
 * ON ENTRY :
 *   rows      the number of rows of all matrices;
 *   cols      the number of columns of all matrices;
 *   Arehihihi has space for a matrix, of dimensions rows and cols;
 *   Arelohihi has space for a matrix, of dimensions rows and cols;
 *   Arehilohi has space for a matrix, of dimensions rows and cols;
 *   Arelolohi has space for a matrix, of dimensions rows and cols;
 *   Arehihilo has space for a matrix, of dimensions rows and cols;
 *   Arelohilo has space for a matrix, of dimensions rows and cols;
 *   Arehilolo has space for a matrix, of dimensions rows and cols;
 *   Arelololo has space for a matrix, of dimensions rows and cols;
 *   Aimhihihi has space for a matrix, of dimensions rows and cols;
 *   Aimlohihi has space for a matrix, of dimensions rows and cols;
 *   Aimhihihi has space for a matrix, of dimensions rows and cols;
 *   Aimlohihi has space for a matrix, of dimensions rows and cols;
 *   Aimhilolo has space for a matrix, of dimensions rows and cols;
 *   Aimlololo has space for a matrix, of dimensions rows and cols;
 *   Aimhilolo has space for a matrix, of dimensions rows and cols;
 *   Aimlololo has space for a matrix, of dimensions rows and cols.
 *
 * ON RETURN :
 *   Arehihihi has the highest doubles of the real parts 
 *             of a random upper triangular matrix A;
 *   Arelohihi has the second highest doubles of the real parts of A;
 *   Arehilohi has the third highest doubles of the real parts of A;
 *   Arelolohi has the fourth highest doubles of the real parts of A;
 *   Arehihilo has the fourth lowest doubles of the real parts of A;
 *   Arelohilo has the third lowest doubles of the real parts of A;
 *   Arehilolo has the second lowest doubles of the real parts of A;
 *   Arelololo has the lowest doubles of the real parts of A;
 *   Aimhihihi has the highest doubles of the imaginary parts of A;
 *   Aimlohihi has the second highest doubles of the imaginary parts of A;
 *   Aimhilohi has the third highest doubles of the imaginary parts of A;
 *   Aimlolohi has the fourth highest doubles of the imaginary parts of A;
 *   Aimhihilo has the fourth lowest doubles of the imaginary parts of A;
 *   Aimlohilo has the third lowest doubles of the imaginary parts of A;
 *   Aimhilolo has the second lowest doubles of the imaginary parts of A;
 *   Aimlololo has the lowest doubles of the imaginary parts of A. */

void random_dbl8_matrix
 ( int rows, int cols,
   double **Ahihihi, double **Alohihi, double **Ahilohi, double **Alolohi,
   double **Ahihilo, double **Alohilo, double **Ahilolo, double **Alololo );
/*
 * DESCRIPTION :
 *   Returns eight upper triangular matrices of dimensions rows and cols,
 *   to represent a real random triangular matrix of octo doubles.
 *
 * ON ENTRY :
 *   rows      the number of rows of all matrices;
 *   cols      the number of columns of all matrices;
 *   Ahihihi   has space for a matrix, of dimensions rows and cols;
 *   Alohihi   has space for a matrix, of dimensions rows and cols;
 *   Ahilohi   has space for a matrix, of dimensions rows and cols;
 *   Alolohi   has space for a matrix, of dimensions rows and cols;
 *   Ahihilo   has space for a matrix, of dimensions rows and cols;
 *   Alohilo   has space for a matrix, of dimensions rows and cols;
 *   Ahilolo   has space for a matrix, of dimensions rows and cols;
 *   Alololo   has space for a matrix, of dimensions rows and cols.
 *
 * ON RETURN :
 *   Ahihihi   are the highest doubles of a matrix A
 *             of randomly generated octo doubles;
 *   Alohihi   are the second highest doubles of A;
 *   Ahilohi   are the third highest doubles of A;
 *   Alolohi   are the fourth highest doubles of A;
 *   Ahihilo   are the fourth lowest doubles of A;
 *   Alohilo   are the third lowest doubles of A;
 *   Ahilolo   are the second lowest doubles of A;
 *   Alololo   are the lowest doubles of A. */

void random_cmplx8_matrix
 ( int rows, int cols,
   double **Arehihihi, double **Arelohihi,
   double **Arehilohi, double **Arelolohi,
   double **Arehihilo, double **Arelohilo,
   double **Arehilolo, double **Arelololo,
   double **Aimhihihi, double **Aimlohihi,
   double **Aimhilohi, double **Aimlolohi,
   double **Aimhihilo, double **Aimlohilo,
   double **Aimhilolo, double **Aimlololo );
/*
 * DESCRIPTION :
 *   Returns 16 upper triangular matrices of dimensions rows and cols,
 *   to represent a complex random matrix of octo doubles.
 *
 * ON ENTRY :
 *   rows      the number of rows of all matrices;
 *   cols      the number of columns of all matrices;
 *   Arehihihi has space for a matrix, of dimensions rows and cols;
 *   Arelohihi has space for a matrix, of dimensions rows and cols;
 *   Arehilohi has space for a matrix, of dimensions rows and cols;
 *   Arelolohi has space for a matrix, of dimensions rows and cols;
 *   Arehihilo has space for a matrix, of dimensions rows and cols;
 *   Arelohilo has space for a matrix, of dimensions rows and cols;
 *   Arehilolo has space for a matrix, of dimensions rows and cols;
 *   Arelololo has space for a matrix, of dimensions rows and cols;
 *   Aimhihihi has space for a matrix, of dimensions rows and cols;
 *   Aimlohihi has space for a matrix, of dimensions rows and cols;
 *   Aimhihihi has space for a matrix, of dimensions rows and cols;
 *   Aimlohihi has space for a matrix, of dimensions rows and cols;
 *   Aimhilolo has space for a matrix, of dimensions rows and cols;
 *   Aimlololo has space for a matrix, of dimensions rows and cols;
 *   Aimhilolo has space for a matrix, of dimensions rows and cols;
 *   Aimlololo has space for a matrix, of dimensions rows and cols.
 *
 * ON RETURN :
 *   Arehihihi has the highest doubles of the real parts 
 *             of a random matrix A;
 *   Arelohihi has the second highest doubles of the real parts of A;
 *   Arehilohi has the third highest doubles of the real parts of A;
 *   Arelolohi has the fourth highest doubles of the real parts of A;
 *   Arehihilo has the fourth lowest doubles of the real parts of A;
 *   Arelohilo has the third lowest doubles of the real parts of A;
 *   Arehilolo has the second lowest doubles of the real parts of A;
 *   Arelololo has the lowest doubles of the real parts of A;
 *   Aimhihihi has the highest doubles of the imaginary parts of A;
 *   Aimlohihi has the second highest doubles of the imaginary parts of A;
 *   Aimhilohi has the third highest doubles of the imaginary parts of A;
 *   Aimlolohi has the fourth highest doubles of the imaginary parts of A;
 *   Aimhihilo has the fourth lowest doubles of the imaginary parts of A;
 *   Aimlohilo has the third lowest doubles of the imaginary parts of A;
 *   Aimhilolo has the second lowest doubles of the imaginary parts of A;
 *   Aimlololo has the lowest doubles of the imaginary parts of A. */

#endif
