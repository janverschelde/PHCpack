/* Collection of functions for vectored double double arithmetic. */

#ifndef __VECTORED_DOUBLE_DOUBLES_H__
#define __VECTORED_DOUBLE_DOUBLES_H__

bool is_dd_quarter_balanced
 ( double xhi0, double xhi1, double xhi2, double xhi3,
   double xlo0, double xlo1, double xlo2, double xlo3, int vrblvl=0 );
/*
 * Given the eight quarters of a double double,
 * returns true if they balanced, false otherwise.
 * If vrblvl > 0, then extra output is written. */

bool is_dd_12split_balanced
 ( double xhi0, double xhi1, double xhi2, double xhi3,
   double xlo0, double xlo1, double xlo2, double xlo3,
   double xlo4, double xlo5, double xlo6, double xlo7, int vrblvl=0 );
/*
 * Given the 12 parts of a double double,
 * returns true if they balanced, false otherwise.
 * If vrblvl > 0, then extra output is written. */

void dd_balance_quarters
 ( double *xhi0, double *xhi1, double *xhi2, double *xhi3,
   double *xlo0, double *xlo1, double *xlo2, double *xlo3, int vrblvl=0 );
/*
 * If unbalanced, balances the eight quarters of a double double.
 * If vrblvl > 0, then extra output is written. */

void dd_12split_balance
 ( double *xhi0, double *xhi1, double *xhi2, double *xhi3,
   double *xlo0, double *xlo1, double *xlo2, double *xlo3,
   double *xlo4, double *xlo5, double *xlo6, double *xlo7, int vrblvl=0 );
/*
 * If unbalanced, balances the 12 parts of a double double.
 * If vrblvl > 0, then extra output is written. */

void make_dd_exponent_zero ( double *xhi, double *xlo, int vrblvl=0 );
/*
 * Multiplies the double double so that the leading double has
 * exponent zero.  Writes extra output if vrblvl > 0. */

void quarter_double_double
 ( double xhi, double xlo,
   double *xhi0, double *xhi1, double *xhi2, double *xhi3,
   double *xlo0, double *xlo1, double *xlo2, double *xlo3, int vrblvl=0 );
/*
 * Quarters the high part xhi and the low part xlo of a double double,
 * resulting in (xhi0, xhi1, xhi2, xhi3) and in (xlo0, xlo1, xlo2, xlo3).
 * The quarters are balanced. */

void split_double_double
 ( double xhi, double xlo,
   double *xhi0, double *xhi1, double *xhi2, double *xhi3,
   double *xlo0, double *xlo1, double *xlo2, double *xlo3,
   double *xlo4, double *xlo5, double *xlo6, double *xlo7, int vrblvl=0 );
/*
 * Quarters the high part xhi and octo splits the low part xlo,
 * resulting in (xhi0, xhi1, xhi2, xhi3) and in (xlo0, xlo1, .., xlo4). */

void quarter_dd_vector
 ( int dim, double *xhi, double *xlo,
   double *xhi0, double *xhi1, double *xhi2, double *xhi3,
   double *xlo0, double *xlo1, double *xlo2, double *xlo3, int vrblvl=0 );
/*
 * Given a vector of size dim in xhi and xlo,
 * quarters the high part xhi and the low part xlo of a double double,
 * resulting in (xhi0, xhi1, xhi2, xhi3) and in (xlo0, xlo1, xlo2, xlo3). */

void split_dd_vector
 ( int dim, double *xhi, double *xlo,
   double *xhi0, double *xhi1, double *xhi2, double *xhi3,
   double *xlo0, double *xlo1, double *xlo2, double *xlo3,
   double *xlo4, double *xlo5, double *xlo6, double *xlo7, int vrblvl=0 );
/*
 * Given a vector of size dim in xhi and xlo,
 * quarters the high part xhi and splits the low part xlo of a double double,
 * resulting in (xhi0, xhi1, xhi2, xhi3) and in (xlo0, xlo1, .., xlo7). */

void quarter_dd_matrix
 ( int nrows, int ncols, double **Ahi, double **Alo,
   double **Ahi0, double **Ahi1, double **Ahi2, double **Ahi3,
   double **Alo0, double **Alo1, double **Alo2, double **Alo3, int vrblvl=0 );
/*
 * Given a matrix of nrows rows and ncols columns in (Ahi, Alo),
 * quarters the matrix into 8 matrices of the low and high parts. */

void split_dd_matrix
 ( int nrows, int ncols, double **Ahi, double **Alo,
   double **Ahi0, double **Ahi1, double **Ahi2, double **Ahi3,
   double **Alo0, double **Alo1, double **Alo2, double **Alo3,
   double **Alo4, double **Alo5, double **Alo6, double **Alo7, int vrblvl=0 );
/*
 * Given a matrix of nrows rows and ncols columns in (Ahi, Alo),
 * quarters the high parts of the matrix and octo splits the low parts
 * into 12 matrices of the splitted parts. */

void to_double_double8sum
 ( double xhi0, double xhi1, double xhi2, double xhi3,
   double xlo0, double xlo1, double xlo2, double xlo3,
   double *xhi, double *xlo );
/*
 * Given the quarters of the high part in (xhi0, xhi1, xhi2, xhi3)
 * and the quarters of the low part in (xlo0, xlo1, xlo2, xlo3),
 * returns in xhi and xlo the high and low parts of a double double,
 * using double double arithmetic. */

void to_double_double12sum
 ( double xhi0, double xhi1, double xhi2, double xhi3,
   double xlo0a, double xlo1a, double xlo2a, double xlo3a,
   double xlo0b, double xlo1b, double xlo2b, double xlo3b,
   double *xhi, double *xlo );
/*
 * Similar to to_double_double_8sum, but the quarters of the low part are
 * in (xlo0a, xlo1a, xlo2a, xl03a) and (xlo0b, xlo1b, xlo2b, xlo3b).
 * All quarters are added in double double arithmetic to make (xhi, xlo). */

void to_double_double8sum_matrix
 ( int nrows, int ncols,
   double **Ahi0, double **Ahi1, double **Ahi2, double **Ahi3,
   double **Alo0, double **Alo1, double **Alo2, double **Alo3,
   double **Ahi, double **Alo );
/*
 * Given the quarters of an nrows-by-ncols matrix,
 * returns the high and low parts of the double doubles in the matrix. */

void to_double_double12sum_matrix
 ( int nrows, int ncols,
   double **Ahi0, double **Ahi1, double **Ahi2, double **Ahi3,
   double **Alo0a, double **Alo1a, double **Alo2a, double **Alo3a,
   double **Alo0b, double **Alo1b, double **Alo2b, double **Alo3b,
   double **Ahi, double **Alo );
/*
 * Given the quarters of an nrows-by-ncols matrix,
 * returns the high and low parts of the double doubles in the matrix.
 * Similar to the 8-sum version, the quarters of the low parts are
 * given in separate arrays. */

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

void recursive_dd_product
 ( int dim, double *xhi, double *xlo, double *yhi, double *ylo,
   double *prdhi, double *prdlo );
/*
 * Uses recursive summatiom for better accuracy. */

void double_double_matmatmul
 ( int nrows, int ncols, int dim,
   double **Ahi, double **Alo, double **Bhi, double **Blo,
   double **Chi, double **Clo );
/*
 * Given in (Ahi, Alo) is an nrows-by-dim matrix A, and
 * given in (Bhi, Blo) is an dim-by-ncols matrix B,
 * returns in (Chi, Clo) the matrix matrix multiplication of A with B. */

void transpose_dd_matrix
 ( int nrows, int ncols,
   double **Ahi, double **Alo, double **Thi, double **Tlo );
/*
 * Given an nrows-by-ncols double double matrix A (Ahi, Alo),
 * returns an ncols-by-nrols double double matrix T (Thi, Tlo),
 * which is the transpose of A. */

void recursive_dd_matmatmul
 ( int nrows, int ncols, int dim,
   double **Ahi, double **Alo, double **Bhi, double **Blo,
   double **Chi, double **Clo );
/*
 * Given in (Ahi, Alo) is a row major nrows-by-dim matrix A, and
 * given in (Bhi, Blo) is a column major dim-by-ncols matrix B,
 * returns in (Chi, Clo) the matrix matrix multiplication of A with B,
 * using the recursive summation of the inner product,
 * for better accuracy, for which B must be column major. */

void vectored_dd_product8sum
 ( int dim,
   double *x0, double *x1, double *x2, double *x3,
   double *x4, double *x5, double *x6, double *x7,
   double *y0, double *y1, double *y2, double *y3,
   double *y4, double *y5, double *y6, double *y7,
   double *s0, double *s1, double *s2, double *s3,
   double *s4, double *s5, double *s6, double *s7 );
/*
 * Makes the vectored product of x and y, with the sums of the product
 * returned in s0, s1, .., s7. */

void vectored_dd_product12sum
 ( int dim,
   double *x0, double *x1, double *x2, double *x3,
   double *x4, double *x5, double *x6, double *x7,
   double *y0, double *y1, double *y2, double *y3,
   double *y4, double *y5, double *y6, double *y7,
   double *s0, double *s1, double *s2, double *s3,
   double *s4a, double *s5a, double *s6a, double *s7a,
   double *s4b, double *s5b, double *s6b, double *s7b );
/*
 * Similar to vectored_dd_product8sum, but the last eight doubles
 * of vectored_dd_product8sum are computed in two sums each,
 * that is the s4 of vectored_dd_product8sum is s4a + s4b
 * of this vectored_dd_product12sum, as is also s5, s6, s7. */

void vectored_dd_product
 ( int dim,
   double *x0, double *x1, double *x2, double *x3,
   double *x4, double *x5, double *x6, double *x7,
   double *x8, double *x9, double *x10, double *x11,
   double *y0, double *y1, double *y2, double *y3,
   double *y4, double *y5, double *y6, double *y7,
   double *y8, double *y9, double *y10, double *y11,
   double *s0, double *s1, double *s2, double *s3,
   double *s4, double *s5, double *s6, double *s7,
   double *s8, double *s9, double *s10, double *s11 );
/*
 * Makes the vectored product of x and y,
 * given by the quarters of their high parts
 * and the eight doubles of their low parts,
 * returning the sums of the product in s0, s1, .., s11 */

void transpose_dd_quarters
 ( int nrows, int ncols,
   double **A0, double **A1, double **A2, double **A3,
   double **A4, double **A5, double **A6, double **A7,
   double **T0, double **T1, double **T2, double **T3,
   double **T4, double **T5, double **T6, double **T7 );
/*
 * Returns in T0, T1, .., T7 the transpose of A0, A1, .., A7,
 * where A is nrows-by-ncols, T is ncols-by-nrows */

void transpose_dd_splits
 ( int nrows, int ncols,
   double **A0, double **A1, double **A2, double **A3,
   double **A4, double **A5, double **A6, double **A7,
   double **A8, double **A9, double **A10, double **A11,
   double **T0, double **T1, double **T2, double **T3,
   double **T4, double **T5, double **T6, double **T7,
   double **T8, double **T9, double **T10, double **T11 );
/*
 * Returns in T0, T1, .., T11 the transpose of A0, A1, .., A11,
 * where A is nrows-by-ncols, T is ncols-by-nrows */

void vectored_dd_matmatmul8sum
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

void vectored_dd_matmatmul12sum
 ( int nrows, int ncols, int dim,
   double **A0, double **A1, double **A2, double **A3,
   double **A4, double **A5, double **A6, double **A7,
   double **B0, double **B1, double **B2, double **B3,
   double **B4, double **B5, double **B6, double **B7,
   double **C0, double **C1, double **C2, double **C3,
   double **C4a, double **C5a, double **C6a, double **C7a,
   double **C4b, double **C5b, double **C6b, double **C7b );
/*
 * Similar to the 8-sum vectored double double matrix product,
 * except that the quarters for the low parts are summed in two halves. */

void vectored_dd_matmatmul
 ( int nrows, int ncols, int dim,
   double **A0, double **A1, double **A2, double **A3,
   double **A4, double **A5, double **A6, double **A7,
   double **A8, double **A9, double **A10, double **A11,
   double **B0, double **B1, double **B2, double **B3,
   double **B4, double **B5, double **B6, double **B7,
   double **B8, double **B9, double **B10, double **B11,
   double **C0, double **C1, double **C2, double **C3,
   double **C4, double **C5, double **C6, double **C7,
   double **C8, double **C9, double **C10, double **C11 );
/*
 * Makes the vectored product of the matrix A and B,
 * given by their parts in A0, A1, .., B0, B1, ...,
 * resulting in the 12 parts in the nrows-by-ncols matrix C.
 * The number of columns of A and te number of rows in B is dim,
 * but the matrix B is column major, while A and C are row major. */

void dd_convolute_quarters
 ( int nrows, int ncols,
   double **A0, double **A1, double **A2, double **A3,
   double **A4, double **A5, double **A6, double **A7, double **cA );
/*
 * Given in A0, A1, ..., A7, the quarters of the parts of 
 * an nrows-by-ncols double double matrix, returns in cA
 * the convoluted matrix where each element in the original matrix A
 * is replaced by an 8-by-8 convolution matrix.
 * Therefore, cA is an 8*nrows-by-8*ncols matrix. */

void dd_convolute_12splits
 ( int nrows, int ncols,
   double **A0, double **A1, double **A2, double **A3,
   double **A4, double **A5, double **A6, double **A7,
   double **A8, double **A9, double **A10, double **A11, double **cA );
/*
 * Given in A0, A1, ..., A11, the parts of a 12-splitted
 * nrows-by-ncols double double matrix, returns in cA
 * the convoluted matrix where each element in the original matrix A
 * is replaced by an 12-by-12 convolution matrix.
 * Therefore, cA is an 12*nrows-by-12*ncols matrix. */

void dd_stack_quarters
 ( int nrows, int ncols,
   double **A0, double **A1, double **A2, double **A3,
   double **A4, double **A5, double **A6, double **A7, double **sA );
/*
 * Given in A0, A1, ..., A7, the quarters of the parts of 
 * an nrows-by-ncols double double matrix, returns in sA
 * the stacked matrix where each element in the original matrix A
 * is replaced by an 8-by-1 column.
 * Therefore, cA is an 8*nrows-by-ncols matrix. */

void dd_stack_12splits
 ( int nrows, int ncols,
   double **A0, double **A1, double **A2, double **A3,
   double **A4, double **A5, double **A6, double **A7,
   double **A8, double **A9, double **A10, double **A11, double **sA );
/*
 * Given in A0, A1, ..., A11, the parts of a 12-splitted
 * nrows-by-ncols double double matrix, returns in sA
 * the stacked matrix where each element in the original matrix A
 * is replaced by an 12-by-1 column.
 * Therefore, cA is an 12*nrows-by-ncols matrix. */

void extract_dd_quarters
 ( int nrows, int ncols, double **qC,
   double **D0, double **D1, double **D2, double **D3,
   double **D4, double **D5, double **D6, double **D7 );
/*
 * Given is in qC the quartered project as an 8*nrows-by-ncols matrix,
 * extracts the quarters into the nrows-by-ncols matrices D0, D1, ..., D7. */

void extract_dd_12splits
 ( int nrows, int ncols, double **qC,
   double **D0, double **D1, double **D2, double **D3,
   double **D4, double **D5, double **D6, double **D7,
   double **D8, double **D9, double **D10, double **D11 );
/*
 * Given is in qC the quartered project as an 12*nrows-by-ncols matrix,
 * extracts the parts into the nrows-by-ncols matrices D0, D1, ..., D11. */

#endif
