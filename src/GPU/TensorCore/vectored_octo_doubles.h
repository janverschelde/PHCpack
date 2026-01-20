/* Collection of functions for vectored octo double arithmetic. */

#ifndef __VECTORED_OCTO_DOUBLES_H__
#define __VECTORED_OCTO_DOUBLES_H__

void quarter_octo_double
 ( double xhihihi, double xlohihi, double xhilohi, double xlolohi,
   double xhihilo, double xlohilo, double xhilolo, double xlololo,
   double *xhihihi0, double *xhihihi1, double *xhihihi2, double *xhihihi3,
   double *xlohihi0, double *xlohihi1, double *xlohihi2, double *xlohihi3,
   double *xhilohi0, double *xhilohi1, double *xhilohi2, double *xhilohi3,
   double *xlolohi0, double *xlolohi1, double *xlolohi2, double *xlolohi3,
   double *xhihilo0, double *xhihilo1, double *xhihilo2, double *xhihilo3,
   double *xlohilo0, double *xlohilo1, double *xlohilo2, double *xlohilo3,
   double *xhilolo0, double *xhilolo1, double *xhilolo2, double *xhilolo3,
   double *xlololo0, double *xlololo1, double *xlololo2, double *xlololo3 );
/*
 * Quarters the parts of an octo double
 * (xhihihi, xlohihi, xhilohi, xlolohi, xhihilo, xlohilo, xhilolo, xlololo),
 * resulting in (xhihihi0, xhihihi1, xhihihi2, xhihihi3),
 * (xlohihi0, xlohihi1, xlohihi2, xlohihi3),
 * (xhilohi0, xhilohi1, xhilohi2, xhilohi3),
 * (xlolohi0, xlolohi1, xlolohi2, xlolohi3),
 * (xhihilo0, xhihilo1, xhihilo2, xhihilo3),
 * (xlohilo0, xlohilo1, xlohilo2, xlohilo3),
 * (xhilolo0, xhilolo1, xhilolo2, xhilolo3), and
 * (xlololo0, xlololo1, xlololo2, xlololo3). */

void quarter_od_vector
 ( int dim,
   double *xhihihi, double *xlohihi, double *xhilohi, double *xlolohi,
   double *xhihilo, double *xlohilo, double *xhilolo, double *xlololo,
   double *xhihihi0, double *xhihihi1, double *xhihihi2, double *xhihihi3,
   double *xlohihi0, double *xlohihi1, double *xlohihi2, double *xlohihi3,
   double *xhilohi0, double *xhilohi1, double *xhilohi2, double *xhilohi3,
   double *xlolohi0, double *xlolohi1, double *xlolohi2, double *xlolohi3,
   double *xhihilo0, double *xhihilo1, double *xhihilo2, double *xhihilo3,
   double *xlohilo0, double *xlohilo1, double *xlohilo2, double *xlohilo3,
   double *xhilolo0, double *xhilolo1, double *xhilolo2, double *xhilolo3,
   double *xlololo0, double *xlololo1, double *xlololo2, double *xlololo3 );
/*
 * Given a vector of size dim in 
 * (xhihihi, xlohihi, xhilohi, xloxlohi, xhihilo, xlohilo, xhilolo, xloxlolo),
 * quarters the parts of an octo double,
 * resulting in 32 vectors of size dim. */

void quarter_od_matrix
 ( int nrows, int ncols,
   double **Ahihihi, double **Alohihi, double **Ahilohi, double **Alolohi,
   double **Ahihilo, double **Alohilo, double **Ahilolo, double **Alololo,
   double **Ahihihi0, double **Ahihihi1, double **Ahihihi2, double **Ahihihi3,
   double **Alohihi0, double **Alohihi1, double **Alohihi2, double **Alohihi3,
   double **Ahilohi0, double **Ahilohi1, double **Ahilohi2, double **Ahilohi3,
   double **Alolohi0, double **Alolohi1, double **Alolohi2, double **Alolohi3,
   double **Ahihilo0, double **Ahihilo1, double **Ahihilo2, double **Ahihilo3,
   double **Alohilo0, double **Alohilo1, double **Alohilo2, double **Alohilo3,
   double **Ahilolo0, double **Ahilolo1, double **Ahilolo2, double **Ahilolo3,
   double **Alololo0, double **Alololo1, double **Alololo2, double **Alololo3 );
/*
 * Given a matrix of nrows rows and ncols columns in
 * (Ahihihi, Alohihi, Ahilohi, Alolohi, Ahihilo, Alohilo, Ahilolo, Alololo),
 * quarters the matrix into 32 matrices of the parts. */

void to_octo_double
 ( double xhihihi0, double xhihihi1, double xhihihi2, double xhihihi3,
   double xlohihi0, double xlohihi1, double xlohihi2, double xlohihi3,
   double xhilohi0, double xhilohi1, double xhilohi2, double xhilohi3,
   double xlolohi0, double xlolohi1, double xlolohi2, double xlolohi3,
   double xhihilo0, double xhihilo1, double xhihilo2, double xhihilo3,
   double xlohilo0, double xlohilo1, double xlohilo2, double xlohilo3,
   double xhilolo0, double xhilolo1, double xhilolo2, double xhilolo3,
   double xlololo0, double xlololo1, double xlololo2, double xlololo3,
   double *xhihihi, double *xlohihi, double *xhilohi, double *xlolohi,
   double *xhihilo, double *xlohilo, double *xhilolo, double *xlololo );
/*
 * Given the quarters of the parts of a octo double, returns in
 * (xhihihi, xlohihi, xhilohi, xlolohi, xhihilo, xlohilo, xhilolo, xlololo) 
 * the parts of an octo double, using octo double arithmetic. */

void od_write_vector
 ( int dim,
   double *xhihihi, double *xlohihi, double *xhilohi, double *xlolohi,
   double *xhihilo, double *xlohilo, double *xhilolo, double *xlololo );
/*
 * Writes the octo double vector of size dim with parts in
 * xhihihi, xlohihi, xhilohi, xlolohi, xhihilo, xlohilo, xhilolo, xlololo. */

void octo_double_product
 ( int dim,
   double *xhihihi, double *xlohihi, double *xhilohi, double *xlolohi,
   double *xhihilo, double *xlohilo, double *xhilolo, double *xlololo,
   double *yhihihi, double *ylohihi, double *yhilohi, double *ylolohi,
   double *yhihilo, double *ylohilo, double *yhilolo, double *ylololo,
   double *phihihi, double *plohihi, double *philohi, double *plolohi,
   double *phihilo, double *plohilo, double *philolo, double *plololo );
/*
 * Makes the product of two octo double vectors of size dim, given in
 * (xhihihi, xlohihi, xhilohi, xlolohi, xhihilo, xlohilo, xhilolo, xlololo) 
 * and
 * (yhihihi, ylohihi, yhilohi, ylolohi, yhihilo, ylohilo, yhilolo, ylololo),
 * resulting in
 * (phihihi, plohihi, philohi, plolohi, phihilo, plohilo, philolo, plololo). */

void octo_double_matmatmul
 ( int nrows, int ncols, int dim,
   double **Ahihihi, double **Alohihi, double **Ahilohi, double **Alolohi,
   double **Ahihilo, double **Alohilo, double **Ahilolo, double **Alololo,
   double **Bhihihi, double **Blohihi, double **Bhilohi, double **Blolohi,
   double **Bhihilo, double **Blohilo, double **Bhilolo, double **Blololo,
   double **Chihihi, double **Clohihi, double **Chilohi, double **Clolohi,
   double **Chihilo, double **Clohilo, double **Chilolo, double **Clololo );
/*
 * Given in
 * (Ahihihi, Alohihi, Ahilohi, Alolohi, Ahihilo, Alohilo, Ahilolo, Alololo)
 * is an nrows-by-dim matrix A, and given in 
 * (Bhihihi, Blohihi, Bhilohi, Blolohi, Bhihilo, Blohilo, Bhilolo, Blololo)
 * is an dim-by-ncols matrix B, returns in
 * (Chihihi, Clohihi, Chilohi, Clolohi, Chihilo, Clohilo, Chilolo, Clololo)
 * the matrix matrix multiplication of A with B. */

void vectored_od_product
 ( int dim,
   double *x0, double *x1, double *x2, double *x3,
   double *x4, double *x5, double *x6, double *x7,
   double *x8, double *x9, double *x10, double *x11,
   double *x12, double *x13, double *x14, double *x15,
   double *x16, double *x17, double *x18, double *x19,
   double *x20, double *x21, double *x22, double *x23,
   double *x24, double *x25, double *x26, double *x27,
   double *x28, double *x29, double *x30, double *x31,
   double *y0, double *y1, double *y2, double *y3,
   double *y4, double *y5, double *y6, double *y7,
   double *y8, double *y9, double *y10, double *y11,
   double *y12, double *y13, double *y14, double *y15,
   double *y16, double *y17, double *y18, double *y19,
   double *y20, double *y21, double *y22, double *y23,
   double *y24, double *y25, double *y26, double *y27,
   double *y28, double *y29, double *y30, double *y31,
   double *s0, double *s1, double *s2, double *s3,
   double *s4, double *s5, double *s6, double *s7,
   double *s8, double *s9, double *s10, double *s11,
   double *s12, double *s13, double *s14, double *s15,
   double *s16, double *s17, double *s18, double *s19,
   double *s20, double *s21, double *s22, double *s23,
   double *s24, double *s25, double *s26, double *s27,
   double *s28, double *s29, double *s30, double *s31 );
/*
 * Makes the vectored product of x and y, with the sums of the product
 * in s0, s1, etc ... */

void transpose_od_quarters
 ( int nrows, int ncols,
   double **A0, double **A1, double **A2, double **A3,
   double **A4, double **A5, double **A6, double **A7,
   double **A8, double **A9, double **A10, double **A11,
   double **A12, double **A13, double **A14, double **A15,
   double **A16, double **A17, double **A18, double **A19,
   double **A20, double **A21, double **A22, double **A23,
   double **A24, double **A25, double **A26, double **A27,
   double **A28, double **A29, double **A30, double **A31,
   double **T0, double **T1, double **T2, double **T3,
   double **T4, double **T5, double **T6, double **T7,
   double **T8, double **T9, double **T10, double **T11,
   double **T12, double **T13, double **T14, double **T15,
   double **T16, double **T17, double **T18, double **T19,
   double **T20, double **T21, double **T22, double **T23,
   double **T24, double **T25, double **T26, double **T27,
   double **T28, double **T29, double **T30, double **T31 );
/*
 * Returns in T0, T1, ... the transpose of A0, A1, ...
 * where A is nrows-by-ncols, T is ncols-by-nrows */

void vectored_od_matmatmul
 ( int nrows, int ncols, int dim,
   double **A0, double **A1, double **A2, double **A3,
   double **A4, double **A5, double **A6, double **A7,
   double **A8, double **A9, double **A10, double **A11,
   double **A12, double **A13, double **A14, double **A15,
   double **A16, double **A17, double **A18, double **A19,
   double **A20, double **A21, double **A22, double **A23,
   double **A24, double **A25, double **A26, double **A27,
   double **A28, double **A29, double **A30, double **A31,
   double **B0, double **B1, double **B2, double **B3,
   double **B4, double **B5, double **B6, double **B7,
   double **B8, double **B9, double **B10, double **B11,
   double **B12, double **B13, double **B14, double **B15,
   double **B16, double **B17, double **B18, double **B19,
   double **B20, double **B21, double **B22, double **B23,
   double **B24, double **B25, double **B26, double **B27,
   double **B28, double **B29, double **B30, double **B31,
   double **C0, double **C1, double **C2, double **C3,
   double **C4, double **C5, double **C6, double **C7,
   double **C8, double **C9, double **C10, double **C11,
   double **C12, double **C13, double **C14, double **C15,
   double **C16, double **C17, double **C18, double **C19,
   double **C20, double **C21, double **C22, double **C23,
   double **C24, double **C25, double **C26, double **C27,
   double **C28, double **C29, double **C30, double **C31 );
/*
 * Makes the vectored product of the matrix A and B,
 * given by their quarters in A0, A1, .., B0, B1, ...,
 * resulting in the quarters in the nrows-by-ncols matrix C.
 * The number of columns of A and the number of rows in B is dim,
 * but the matrix B is column major, while A and C are row major. */

#endif
