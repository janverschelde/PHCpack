// The file random_series.h specifies functions to generate random series
// in double precision.

#ifndef __random_series_h__
#define __random_series_h__

void dbl_exponential ( int deg, double x, double *s );
/*
 * DESCRIPTION :
 *   Returns power series for exp(x) truncated at degree deg.
 *
 * ON ENTRY :
 *   deg      degree to truncate the series;
 *   x        some double;
 *   s        space for deg+1 doubles.
 *
 * ON RETURN :
 *   s        power series of exp(x) truncated to degree deg. */

void dbl_exponentials ( int deg, double x, double *plux, double *minx );
/*
 * DESCRIPTION :
 *   Returns power series for exp(x) and exp(-x),
 *   with real coefficients in double precision.
 *
 * ON ENTRY :
 *   deg      degree to truncate the series;
 *   x        some double;
 *   plux     space for deg+1 doubles for the exp(+x) series;
 *   minx     space for deg+1 doubles for the exp(-x) series.
 *
 * ON RETURN :
 *   plux     power series of exp(+x) truncated to degree deg;
 *   minx     power series of exp(-x) truncated to degree deg. */

void random_dbl_exponential ( int deg, double *x, double *s );
/*
 * DESCRIPTION :
 *   Returns a power series of exp(x) truncated at degree deg
 *   for a randomly generated double x. */

void random_dbl_exponentials
 ( int deg, double *x, double *plux, double *minx );
/*
 * DESCRIPTION :
 *   Returns power series for exp(x) and exp(-x),
 *   truncated to degree deg for a random double x.
 *   Parameters are the same as dbl_exponentials,
 *   except that x is a return parameter. */

void cmplx_exponentials
 ( int deg, double xre, double xim,
   double *pluxre, double *pluxim, double *minxre, double *minxim );
/*
 * DESCRIPTION :
 *   Returns power series for exp(x) and exp(-x),
 *   with complex coefficients in double precision.
 *
 * ON ENTRY :
 *   deg      degree to truncate the series;
 *   xre      real part of some complex number;
 *   xim      imaginary part of some complex number;
 *   pluxre   space for deg+1 doubles for the real parts
 *            of the exp(+x) series;
 *   pluxim   space for deg+1 doubles for the imaginary parts
 *            of the exp(+x) series;
 *   minxre   space for deg+1 doubles for the real parts
 *            of the exp(-x) series;
 *   minxim   space for deg+1 doubles for the imaginary parts
 *            of the exp(-x) series.
 *
 * ON RETURN :
 *   pluxre   real parts of the power series of exp(+x)
 *            truncated to degree deg;
 *   pluxim   imaginary parts of the power series of exp(+x)
 *            truncated to degree deg;
 *   minxre   real parts of the power series of exp(-x)
 *            truncated to degree deg;
 *   minxim   imaginary parts of the power series of exp(-x)
 *            truncated to degree deg. */

void random_cmplx_exponentials
 ( int deg, double *xre, double *xim,
   double *pluxre, double *pluxim, double *minxre, double *minxim );
/*
 * DESCRIPTION :
 *   Returns power series for exp(x) and exp(-x),
 *   truncated to degree deg for a random complex double x.
 *   Parameters are the same as dbl_exponentials,
 *   except that xre and xim are return parameters. */

#endif
