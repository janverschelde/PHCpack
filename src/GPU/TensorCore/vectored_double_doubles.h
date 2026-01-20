/* Collection of functions for vectored double double arithmetic. */

#ifndef __VECTORED_DOUBLE_DOUBLES_H__
#define __VECTORED_DOUBLE_DOUBLES_H__

void quarter_double_double
 ( double xhi, double xlo,
   double *xhi0, double *xhi1, double *xhi2, double *xhi3,
   double *xlo0, double *xlo1, double *xlo2, double *xlo3 );
/*
 * Quarters the high part xhi and the low part xlo of a double double,
 * resulting in (xhi0, xhi1, xhi2, xhi3) and in (xlo0, xlo1, xlo2, xlo3). */

void quarter_dd_vector
 ( int dim, double *xhi, double *xlo,
   double *xhi0, double *xhi1, double *xhi2, double *xhi3,
   double *xlo0, double *xlo1, double *xlo2, double *xlo3 );
/*
 * Given a vector of size dim in xhi and xlo,
 * quarters the high part xhi and the low part xlo of a double double,
 * resulting in (xhi0, xhi1, xhi2, xhi3) and in (xlo0, xlo1, xlo2, xlo3). */

void quarter_dd_matrix
 ( int nrows, int ncols, double **Ahi, double **Alo,
   double **Ahi0, double **Ahi1, double **Ahi2, double **Ahi3,
   double **Alo0, double **Alo1, double **Alo2, double **Alo3 );
/*
 * Given a matrix of nrows rows and ncols columns in (Ahi, Alo),
 * quarters the matrix into 8 matrices of the low and high parts. */

void to_double_double
 ( double xhi0, double xhi1, double xhi2, double xhi3,
   double xlo0, double xlo1, double xlo2, double xlo3,
   double *xhi, double *xlo );
/*
 * Given the quarters of the high part in (xhi0, xhi1, xhi2, xhi3)
 * and the quarters of the low part in (xlo0, xlo1, xlo2, xlo3),
 * returns in xhi and xlo the high and low parts of a double double,
 * using double double arithmetic. */

void to_double_double_matrix
 ( int nrows, int ncols,
   double **Ahi0, double **Ahi1, double **Ahi2, double **Ahi3,
   double **Alo0, double **Alo1, double **Alo2, double **Alo3,
   double **Ahi, double **Alo );
/*
 * Given the quarters of an nrows-by-ncols matrix,
 * returns the high and low parts of the double doubles in the matix. */

void dd_write_vector ( int dim, double *xhi, double *xlo );
/*
 * Writes the double double vector of size dim with high and low parts
 * respectively in xhi and xlo. */

void double_double_product
 ( int dim, double *xhi, double *xlo, double *yhi, double *ylo,
   double *prdhi, double *prdlo );
/*
 * Makes the product of two double double vectors of size dim,
 * given in (xhi, xlo) and (yhi, ylo), with result in (prdhi, prdlo). */

void double_double_matmatmul
 ( int nrows, int ncols, int dim,
   double **Ahi, double **Alo, double **Bhi, double **Blo,
   double **Chi, double **Clo );
/*
 * Given in (Ahi, Alo) is an nrows-by-dim matrix A, and
 * given in (Bhi, Blo) is an dim-by-ncols matrix B,
 * returns in (Chi, Clo) the matrix matrix multiplication of A with B. */

void vectored_dd_product
 ( int dim,
   double *x0, double *x1, double *x2, double *x3,
   double *x4, double *x5, double *x6, double *x7,
   double *y0, double *y1, double *y2, double *y3,
   double *y4, double *y5, double *y6, double *y7,
   double *s0, double *s1, double *s2, double *s3,
   double *s4, double *s5, double *s6, double *s7 );
/*
 * Makes the vectored product of x and y, with the sums of the product
 * in s0, s1, etc ... */

void transpose_dd_quarters
 ( int nrows, int ncols,
   double **A0, double **A1, double **A2, double **A3,
   double **A4, double **A5, double **A6, double **A7,
   double **T0, double **T1, double **T2, double **T3,
   double **T4, double **T5, double **T6, double **T7 );
/*
 * Returns in T0, T1, ... the transpose of A0, A1, ...
 * where A is nrows-by-ncols, T is ncols-by-nrows */

void vectored_dd_matmatmul
 ( int nrows, int ncols, int dim,
   double **A0, double **A1, double **A2, double **A3,
   double **A4, double **A5, double **A6, double **A7,
   double **B0, double **B1, double **B2, double **B3,
   double **B4, double **B5, double **B6, double **B7,
   double **C0, double **C1, double **C2, double **C3,
   double **C4, double **C5, double **C6, double **C7 );
/*
 * Makes the vectored product of the matrix A and B,
 * given by their quarters in A0, A1, .., B0, B1, ...,
 * resulting in the quarters in the nrows-by-ncols matrix C.
 * The number of columns of A and te number of rows in B is dim,
 * but the matrix B is column major, while A and C are row major. */

void dd_convolute_quarters
 ( int nrows, int ncols,
   double **A0, double **A1, double **A2, double **A3,
   double **A4, double **A5, double **A6, double **A7, double **cA );
/*
 * Given in A0, A1, ... the quarters of the parts of 
 * an nrows-by-ncols double double matrix, returns in cA
 * the convoluted matrix where each element in the original matrix A
 * is replaced by an 8-by-8 convolution matrix.
 * Therefore, cA is an 8*nrows-by-8*ncols matrix. */

void dd_stack_quarters
 ( int nrows, int ncols,
   double **A0, double **A1, double **A2, double **A3,
   double **A4, double **A5, double **A6, double **A7, double **sA );
/*
 * Given in A0, A1, ... the quarters of the parts of 
 * an nrows-by-ncols double double matrix, returns in sA
 * the stacked matrix where each element in the original matrix A
 * is replaced by an 8-by-1 column.
 * Therefore, cA is an 8*nrows-by-ncols matrix. */

void extract_dd_quarters
 ( int nrows, int ncols, double **qC,
   double **D0, double **D1, double **D2, double **D3,
   double **D4, double **D5, double **D6, double **D7 );
/*
 * Given is in qC the quartered project as an 8*nrows-by-ncols matrix,
 * extracts the quarters into the nrows-by-ncols matrices D0, D1, ... */

#endif
