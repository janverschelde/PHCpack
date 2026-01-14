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

#endif
