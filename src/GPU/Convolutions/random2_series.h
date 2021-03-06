// The file random2_series.h specifies functions to generate random series
// in double double precision.

#ifndef __random2_series_h__
#define __random2_series_h__

void dbl2_exponentials
 ( int deg, double xhi, double xlo, 
   double *pluxhi, double *pluxlo, double *minxhi, double *minxlo );
/*
 * DESCRIPTION :
 *   Returns power series for exp(x) and exp(-x),
 *   with real coefficients in double double precision.
 *
 * ON ENTRY :
 *   deg      degree to truncate the series;
 *   xhi      high part of some double double;
 *   xlo      low part of some double double;
 *   pluxhi   space for deg+1 doubles for high parts of the exp(+x) series;
 *   pluxlo   space for deg+1 doubles for low parts of the exp(+x) series;
 *   minxhi   space for deg+1 doubles for high parts of the exp(-x) series;
 *   minxlo   space for deg+1 doubles for low parts of the exp(-x) series.
 *
 * ON RETURN :
 *   pluxhi   high parts of the series of exp(+x), truncated to degree deg;
 *   pluxlo   low parts of the series of exp(+x), truncated to degree deg;
 *   minxhi   high parts of the series of exp(-x), truncated to degree deg;
 *   minxlo   low parts of the series of exp(-x), truncated to degree deg. */

void random_dbl2_exponentials
 ( int deg, double *xhi, double *xlo,
   double *pluxhi, double *pluxlo, double *minxhi, double *minxlo );
/*
 * DESCRIPTION :
 *   Returns power series for exp(x) and exp(-x),
 *   truncated to degree deg for a random double double x.
 *   Parameters are the same as dbl2_exponentials,
 *   except that xhi and xlo are return parameters. */

void cmplx2_exponentials
 ( int deg, double xrehi, double xrelo, double ximhi, double ximlo,
   double *pluxrehi, double *pluxrelo, double *pluximhi, double *pluximlo,
   double *minxrehi, double *minxrelo, double *minximhi, double *minximlo );
/*
 * DESCRIPTION :
 *   Returns power series for exp(x) and exp(-x),
 *   with complex coefficients in double double precision.
 *
 * ON ENTRY :
 *   deg       degree to truncate the series;
 *   xrehi     high double of the real part of some complex number;
 *   xrelo     low double of the real part of some complex number;
 *   ximhi     high double of the imaginary part of some complex number;
 *   ximlo     low double of the imaginary part of some complex number;
 *   pluxrehi  space for deg+1 doubles for the high doubles of the real parts
 *             of the exp(+x) series;
 *   pluxrelo  space for deg+1 doubles for the low doubles of the real parts
 *             of the exp(+x) series;
 *   pluximhi  space for deg+1 doubles for the high doubles of the imaginary
 *             parts of the exp(+x) series;
 *   pluximlo  space for deg+1 doubles for the low doubles of the imaginary
 *             parts of the exp(+x) series;
 *   minxrehi  space for deg+1 doubles for the high doubles of the real parts
 *             of the exp(-x) series;
 *   minxrelo  space for deg+1 doubles for the low doubles of the real parts
 *             of the exp(-x) series;
 *   minximhi  space for deg+1 doubles for the high doubles of the imaginary
 *             parts of the exp(-x) series.
 *   minximlo  space for deg+1 doubles for the low doubles of the imaginary
 *             parts of the exp(-x) series.
 *
 * ON RETURN :
 *   pluxrehi  high doubles of the real parts of the power series of exp(+x)
 *             truncated to degree deg;
 *   pluxrelo  low doubles of the real parts of the power series of exp(+x)
 *             truncated to degree deg;
 *   pluximhi  high doubles of the imaginary parts of the power series 
 *             of exp(+x) truncated to degree deg;
 *   pluximlo  low doubles of the imaginary parts of the power series 
 *             of exp(+x) truncated to degree deg;
 *   minxrehi  high doubles of the real parts of the power series of exp(-x)
 *             truncated to degree deg;
 *   minxrelo  low doubles of the real parts of the power series of exp(-x)
 *             truncated to degree deg;
 *   minximhi  high doubles of the imaginary parts of the power series
 *             of exp(-x) truncated to degree deg;
 *   minximlo  low doubles of the imaginary parts of the power series
 *             of exp(-x) truncated to degree deg. */

void random_cmplx2_exponentials
 ( int deg, double *xrehi, double *xrelo, double *ximhi, double *ximlo,
   double *pluxrehi, double *pluxrelo, double *pluximhi, double *pluximlo,
   double *minxrehi, double *minxrelo, double *minximhi, double *minximlo );
/*
 * DESCRIPTION :
 *   Returns power series for exp(x) and exp(-x),
 *   truncated to degree deg for a random complex double double x.
 *   Parameters are the same as cmplx2_exponentials,
 *   except that xrehi, xrelo, ximhi, and ximlo are return parameters. */

#endif
