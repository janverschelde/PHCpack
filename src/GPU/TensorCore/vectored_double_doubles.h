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

void to_double_double
 ( double xhi0, double xhi1, double xhi2, double xhi3,
   double xlo0, double xlo1, double xlo2, double xlo3,
   double *xhi, double *xlo );
/*
 * Given the quarters of the high part in (xhi0, xhi1, xhi2, xhi3)
 * and the quarters of the low part in (xlo0, xlo1, xlo2, xlo3),
 * returns in xhi and xlo the high and low parts of a double double,
 * using double double arithmetic. */

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

#endif
