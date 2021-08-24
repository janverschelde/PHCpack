// The file dbl4_test_utilities.h specifies common test functions on
// vectors and matrices in quad double precision.

#ifndef __dbl4_test_utilities_h__
#define __dbl4_test_utilities_h__

double dbl4_Difference_Sum
 ( int n, double *xhihi, double *xlohi, double *xhilo, double *xlolo,
   double *yhihi, double *ylohi, double *yhilo, double *ylolo );
/*
 * DESCRIPTION :
 *   Returns the sum of the absolute value of the differences
 *   between the entries of the n-dimensional vectors x and y. */

double cmplx4_Difference_Sum
 ( int n, double *xrehihi, double *xrelohi, double *xrehilo, double *xrelolo,
   double *ximhihi, double *ximlohi, double *ximhilo, double *ximlolo,
   double *yrehihi, double *yrelohi, double *yrehilo, double *yrelolo,
   double *yimhihi, double *yimlohi, double *yimhilo, double *yimlolo );
/*
 * DESCRIPTION :
 *   Returns the sum of the absolute value of the differences
 *   between the entries of the n-dimensional vectors x and y. */

double dbl4_Column_Sum
 ( int dim, int col,
   double **Ahihi, double **Alohi, double **Ahilo, double **Alolo );
/*
 * DESCRIPTION :
 *   Returns the sum of the absolute values of all dim elements
 *   of the matrix A in the column with index col. */

double cmplx4_Column_Sum
 ( int dim, int col,
   double **Arehihi, double **Arelohi, double **Arehilo, double **Arelolo,
   double **Aimhihi, double **Aimlohi, double **Aimhilo, double **Aimlolo );
/*
 * DESCRIPTION :
 *   Returns the sum of the absolute values of all dim elements
 *   of the matrix A in the column with index col. */

double dbl4_Max_Column_Sum
 ( int dim, double **Ahihi, double **Alohi, double **Ahilo, double **Alolo );
/*
 * DESCRIPTION :
 *   Returns the maximum column sum of the elements
 *   in the matrix A of dimension dim. */

double cmplx4_Max_Column_Sum
 ( int dim,
   double **Arehihi, double **Arelohi, double **Arehilo, double **Arelolo,
   double **Aimhihi, double **Aimlohi, double **Aimhilo, double **Aimlolo );
/*
 * DESCRIPTION :
 *   Returns the maximum column sum of the elements
 *   in the matrix A of dimension dim. */

double dbl4_condition
 ( int dim, double **Ahihi, double **Alohi, double **Ahilo, double **Alolo,
   double **invAhihi, double **invAlohi,
   double **invAhilo, double **invAlolo );
/*
 * DESCRIPTION :
 *   Returns the condition number of A using the 1-norm
 *   on the matrix of dimension dim and its inverse invA. */

double cmplx4_condition
 ( int dim,
   double **Arehihi, double **Arelohi, double **Arehilo, double **Arelolo,
   double **Aimhihi, double **Aimlohi, double **Aimhilo, double **Aimlolo,
   double **invArehihi, double **invArelohi,
   double **invArehilo, double **invArelolo,
   double **invAimhihi, double **invAimlohi,
   double **invAimhilo, double **invAimlolo );
/*
 * DESCRIPTION :
 *   Returns the condition number of A using the 1-norm
 *   on the matrix of dimension dim and its inverse invA. */

double dbl4_Matrix_Difference_Sum
 ( int n, double **Ahihi, double **Alohi, double **Ahilo, double **Alolo,
   double **Bhihi, double **Blohi, double **Bhilo, double **Blolo );
/*
 * DESCRIPTION :
 *   Returns the sum of the absolute value of the differences
 *   between two n-by-n matrices A and B. */

double cmplx4_Matrix_Difference_Sum
 ( int n,
   double **Arehihi, double **Arelohi, double **Arehilo, double **Arelolo,
   double **Aimhihi, double **Aimlohi, double **Aimhilo, double **Aimlolo,
   double **Brehihi, double **Brelohi, double **Brehilo, double **Brelolo,
   double **Bimhihi, double **Bimlohi, double **Bimhilo, double **Bimlolo );
/*
 * DESCRIPTION :
 *   Returns the sum of the absolute value of the differences
 *   between two n-by-n matrices A and B. */

double dbl4_Diagonal_Difference_Sum
 ( int nbt, int szt,
   double **Ahihi, double **Alohi, double **Ahilo, double **Alolo,
   double **Bhihi, double **Blohi, double **Bhilo, double **Blolo );
/*
 * DESCRIPTION :
 *   Returns the sum of the absolute value of all differences 
 *   of the diagonal tiles of the matrices A and B.
 *   The number of tiles equals nbt and the size of each tile is szt. */

double cmplx4_Diagonal_Difference_Sum
 ( int nbt, int szt,
   double **Arehihi, double **Arelohi, double **Arehilo, double **Arelolo,
   double **Aimhihi, double **Aimlohi, double **Aimhilo, double **Aimlolo,
   double **Brehihi, double **Brelohi, double **Brehilo, double **Brelolo,
   double **Bimhihi, double **Bimlohi, double **Bimhilo, double **Bimlolo );
/*
 * DESCRIPTION :
 *   Returns the sum of the absolute value of all differences 
 *   of the diagonal tiles of the matrices A and B.
 *   The number of tiles equals nbt and the size of each tile is szt. */

void dbl4_random_upper_factor
 ( int dim, double **Ahihi, double **Alohi, double **Ahilo, double **Alolo );
/*
 * DESCRIPTION :
 *   Returns the upper triangular factor of the LU factorization
 *   with row pivoting on a random matrix of dimension dim.
 *   This yields a much better conditioned upper triangular matrix
 *   than the direct generation of a random upper triangular matrix. */

void cmplx4_random_upper_factor
 ( int dim,
   double **Arehihi, double **Arelohi, double **Arehilo, double **Arelolo,
   double **Aimhihi, double **Aimlohi, double **Aimhilo, double **Aimlolo );
/*
 * DESCRIPTION :
 *   Returns the upper triangular factor of the LU factorization
 *   with row pivoting on a random matrix of dimension dim.
 *   This yields a much better conditioned upper triangular matrix
 *   than the direct generation of a random upper triangular matrix. */

#endif
