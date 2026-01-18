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

#endif
