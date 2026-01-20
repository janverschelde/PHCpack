/* Collection of functions for vectored quad double arithmetic. */

#ifndef __VECTORED_QUAD_DOUBLES_H__
#define __VECTORED_QUAD_DOUBLES_H__

void quarter_quad_double
 ( double xhihi, double xlohi, double xhilo, double xlolo,
   double *xhihi0, double *xhihi1, double *xhihi2, double *xhihi3,
   double *xlohi0, double *xlohi1, double *xlohi2, double *xlohi3,
   double *xhilo0, double *xhilo1, double *xhilo2, double *xhilo3,
   double *xlolo0, double *xlolo1, double *xlolo2, double *xlolo3 );
/*
 * Quarters the parts of a quad double (xhihi, xlohi, xhilo, xlolo)
 * resulting in (xhihi0, xhihi1, xhihi2, xhihi3),
 * (xlohi0, xlohi1, xlohi2, xlohi3), (xhilo0, xhilo1, xhilo2, xhilo3),
 * and (xlolo0, xlolo1, xlolo2, xlolo3). */

void quarter_qd_vector
 ( int dim, double *xhihi, double *xlohi, double *xhilo, double *xlolo,
   double *xhihi0, double *xhihi1, double *xhihi2, double *xhihi3,
   double *xlohi0, double *xlohi1, double *xlohi2, double *xlohi3,
   double *xhilo0, double *xhilo1, double *xhilo2, double *xhilo3,
   double *xlolo0, double *xlolo1, double *xlolo2, double *xlolo3 );
/*
 * Given a vector of size dim in (xhihi, xlohi, xhilo, xloxlo), quarters
 * the parts of a quad double, resulting in 16 vectors of size dim. */

void quarter_qd_matrix
 ( int nrows, int ncols,
   double **Ahihi, double **Alohi, double **Ahilo, double **Alolo,
   double **Ahihi0, double **Ahihi1, double **Ahihi2, double **Ahihi3,
   double **Alohi0, double **Alohi1, double **Alohi2, double **Alohi3,
   double **Ahilo0, double **Ahilo1, double **Ahilo2, double **Ahilo3,
   double **Alolo0, double **Alolo1, double **Alolo2, double **Alolo3 );
/*
 * Given a matrix of nrows rows and ncols columns 
 * in (Ahihi, Alohi, Ahilo, Alolo),
 * quarters the matrix into 16 matrices of the parts. */

void to_quad_double
 ( double xhihi0, double xhihi1, double xhihi2, double xhihi3,
   double xlohi0, double xlohi1, double xlohi2, double xlohi3,
   double xhilo0, double xhilo1, double xhilo2, double xhilo3,
   double xlolo0, double xlolo1, double xlolo2, double xlolo3,
   double *xhihi, double *xlohi, double *xhilo, double *xlolo );
/*
 * Given the quarters of the parts of a quad double,
 * returns in (xhihi, xlohi, xhilo, xlolo) the parts of a quad double,
 * using quad double arithmetic. */

void to_quad_double_matrix
 ( int nrows, int ncols,
   double **Ahihi0, double **Ahihi1, double **Ahihi2, double **Ahihi3,
   double **Alohi0, double **Alohi1, double **Alohi2, double **Alohi3,
   double **Ahilo0, double **Ahilo1, double **Ahilo2, double **Ahilo3,
   double **Alolo0, double **Alolo1, double **Alolo2, double **Alolo3,
   double **Ahihi, double **Alohi, double **Ahilo, double **Alolo );
/*
 * Given the quarters of an nrows-by-ncols matrix,
 * returns the high and low parts of the quad doubles in the matix. */

void qd_write_vector
 ( int dim, double *xhihi, double *xlohi, double *xhilo, double *xlolo );
/*
 * Writes the quad double vector of size dim with parts in
 * xhihi, xlohi, xhilo, and xlolo. */

void quad_double_product
 ( int dim, double *xhihi, double *xlohi, double *xhilo, double *xlolo,
            double *yhihi, double *ylohi, double *yhilo, double *ylolo,
   double *prdhihi, double *prdlohi, double *prdhilo, double *prdlolo );
/*
 * Makes the product of two quad double vectors of size dim,
 * given in (xhihi, xlohi, xhilo, xlolo) and (yhihi, ylohi, yhilo, ylolo),
 * with result in (prdhihi, prdlohi, prdhilo, prdlolo). */

void quad_double_matmatmul
 ( int nrows, int ncols, int dim,
   double **Ahihi, double **Alohi, double **Ahilo, double **Alolo,
   double **Bhihi, double **Blohi, double **Bhilo, double **Blolo,
   double **Chihi, double **Clohi, double **Chilo, double **Clolo );
/*
 * Given in (Ahihi, Alohi, Ahilo, Alolo) is an nrows-by-dim matrix A, and
 * given in (Bhihi, Blohi, Bhilo, Blolo) is an dim-by-ncols matrix B,
 * returns in (Chihi, Clohi, Chilo, Clolo) the matrix matrix multiplication
 * of A with B. */

void vectored_qd_product
 ( int dim,
   double *x0, double *x1, double *x2, double *x3,
   double *x4, double *x5, double *x6, double *x7,
   double *x8, double *x9, double *xA, double *xB,
   double *xC, double *xD, double *xE, double *xF,
   double *y0, double *y1, double *y2, double *y3,
   double *y4, double *y5, double *y6, double *y7,
   double *y8, double *y9, double *yA, double *yB,
   double *yC, double *yD, double *yE, double *yF,
   double *s0, double *s1, double *s2, double *s3,
   double *s4, double *s5, double *s6, double *s7,
   double *s8, double *s9, double *sA, double *sB,
   double *sC, double *sD, double *sE, double *sF );
/*
 * Makes the vectored product of x and y, with the sums of the product
 * in s0, s1, etc ... */

void transpose_qd_quarters
 ( int nrows, int ncols,
   double **A0, double **A1, double **A2, double **A3,
   double **A4, double **A5, double **A6, double **A7,
   double **A8, double **A9, double **A10, double **A11,
   double **A12, double **A13, double **A14, double **A15,
   double **T0, double **T1, double **T2, double **T3,
   double **T4, double **T5, double **T6, double **T7,
   double **T8, double **T9, double **T10, double **T11,
   double **T12, double **T13, double **T14, double **T15 );
/*
 * Returns in T0, T1, ... the transpose of A0, A1, ...
 * where A is nrows-by-ncols, T is ncols-by-nrows */

void vectored_qd_matmatmul
 ( int nrows, int ncols, int dim,
   double **A0, double **A1, double **A2, double **A3,
   double **A4, double **A5, double **A6, double **A7,
   double **A8, double **A9, double **A10, double **A11,
   double **A12, double **A13, double **A14, double **A15,
   double **B0, double **B1, double **B2, double **B3,
   double **B4, double **B5, double **B6, double **B7,
   double **B8, double **B9, double **B10, double **B11,
   double **B12, double **B13, double **B14, double **B15,
   double **C0, double **C1, double **C2, double **C3,
   double **C4, double **C5, double **C6, double **C7,
   double **C8, double **C9, double **C10, double **C11,
   double **C12, double **C13, double **C14, double **C15 );
/*
 * Makes the vectored product of the matrix A and B,
 * given by their quarters in A0, A1, .., B0, B1, ...,
 * resulting in the quarters in the nrows-by-ncols matrix C.
 * The number of columns of A and the number of rows in B is dim,
 * but the matrix B is column major, while A and C are row major. */

void qd_convolute_quarters
 ( int nrows, int ncols,
   double **A0, double **A1, double **A2, double **A3,
   double **A4, double **A5, double **A6, double **A7,
   double **A8, double **A9, double **A10, double **A11,
   double **A12, double **A13, double **A14, double **A15, double **cA );
/*
 * Given in A0, A1, ... the quarters of the parts of 
 * an nrows-by-ncols quad double matrix, returns in cA
 * the convoluted matrix where each element in the original matrix A
 * is replaced by a 16-by-16 convolution matrix.
 * Therefore, cA is a 16*nrows-by-16*ncols matrix. */

void qd_stack_quarters
 ( int nrows, int ncols,
   double **A0, double **A1, double **A2, double **A3,
   double **A4, double **A5, double **A6, double **A7,
   double **A8, double **A9, double **A10, double **A11,
   double **A12, double **A13, double **A14, double **A15, double **sA );
/*
 * Given in A0, A1, ... the quarters of the parts of 
 * an nrows-by-ncols quad double matrix, returns in sA
 * the stacked matrix where each element in the original matrix A
 * is replaced by a 16-by-1 column.
 * Therefore, cA is a 16*nrows-by-ncols matrix. */

void extract_qd_quarters
 ( int nrows, int ncols, double **qC,
   double **D0, double **D1, double **D2, double **D3,
   double **D4, double **D5, double **D6, double **D7,
   double **D8, double **D9, double **D10, double **D11,
   double **D12, double **D13, double **D14, double **D15 );
/*
 * Given is in qC the quartered project as a 16*nrows-by-ncols matrix,
 * extracts the quarters into the nrows-by-ncols matrices D0, D1, ... */

#endif
