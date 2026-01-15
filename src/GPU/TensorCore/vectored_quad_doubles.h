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

#endif
