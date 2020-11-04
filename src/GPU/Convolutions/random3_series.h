// The file random3_series.h specifies functions to generate random series
// in triple double precision.

#ifndef __random3_series_h__
#define __random3_series_h__

void dbl3_exponentials
 ( int deg, double xhi, double xmi, double xlo, 
   double *pluxhi, double *pluxmi, double *pluxlo,
   double *minxhi, double *minxmi, double *minxlo );
/*
 * DESCRIPTION :
 *   Returns power series for exp(x) and exp(-x),
 *   with real coefficients in triple double precision.
 *
 * ON ENTRY :
 *   deg      degree to truncate the series;
 *   xhi      high part of some triple double;
 *   xmi      middle part of some triple double;
 *   xlo      low part of some triple double;
 *   pluxhi   space for deg+1 doubles for high parts of the exp(+x) series;
 *   pluxmi   space for deg+1 doubles for middle parts of the exp(+x) series;
 *   pluxlo   space for deg+1 doubles for low parts of the exp(+x) series;
 *   minxhi   space for deg+1 doubles for high parts of the exp(-x) series;
 *   minxmi   space for deg+1 doubles for middle parts of the exp(-x) series;
 *   minxlo   space for deg+1 doubles for low parts of the exp(-x) series.
 *
 * ON RETURN :
 *   pluxhi   high parts of the series of exp(+x), truncated to degree deg;
 *   pluxmi   middle parts of the series of exp(+x), truncated to degree deg;
 *   pluxlo   low parts of the series of exp(+x), truncated to degree deg;
 *   minxhi   high parts of the series of exp(-x), truncated to degree deg;
 *   minxmi   middle parts of the series of exp(-x), truncated to degree deg;
 *   minxlo   low parts of the series of exp(-x), truncated to degree deg. */

void random_dbl3_exponentials
 ( int deg, double *xhi, double *xmi, double *xlo,
   double *pluxhi, double *pluxmi, double *pluxlo,
   double *minxhi, double *minxmi, double *minxlo );
/*
 * DESCRIPTION :
 *   Returns power series for exp(x) and exp(-x),
 *   truncated to degree deg for a random triple double x.
 *   Parameters are the same as dbl3_exponentials,
 *   except that xhi, xmi, and xlo are return parameters. */

void cmplx3_exponentials
 ( int deg, double xrehi, double xremi, double xrelo,
            double ximhi, double ximmi, double ximlo,
   double *pluxrehi, double *pluxremi, double *pluxrelo,
   double *pluximhi, double *pluximmi, double *pluximlo,
   double *minxrehi, double *minxremi, double *minxrelo,
   double *minximhi, double *minximmi, double *minximlo );
/*
 * DESCRIPTION :
 *   Returns power series for exp(x) and exp(-x),
 *   with complex coefficients in triple double precision.
 *
 * ON ENTRY :
 *   deg       degree to truncate the series;
 *   xrehi     high double of the real part of some complex number;
 *   xremi     middle double of the real part of some complex number;
 *   xrelo     low double of the real part of some complex number;
 *   ximhi     high double of the imaginary part of some complex number;
 *   ximmi     middle double of the imaginary part of some complex number;
 *   ximlo     low double of the imaginary part of some complex number;
 *   pluxrehi  space for deg+1 doubles for the high doubles of the real parts
 *             of the exp(+x) series;
 *   pluxremi  space for deg+1 doubles for the middle doubles of the real
 *             parts of the exp(+x) series;
 *   pluxrelo  space for deg+1 doubles for the low doubles of the real parts
 *             of the exp(+x) series;
 *   pluximhi  space for deg+1 doubles for the high doubles of the imaginary
 *             parts of the exp(+x) series;
 *   pluximmi  space for deg+1 doubles for the middle doubles of the imaginary
 *             parts of the exp(+x) series;
 *   pluximlo  space for deg+1 doubles for the low doubles of the imaginary
 *             parts of the exp(+x) series;
 *   minxrehi  space for deg+1 doubles for the high doubles of the real parts
 *             of the exp(-x) series;
 *   minxremi  space for deg+1 doubles for the middle doubles of the real
 *             parts of the exp(-x) series;
 *   minxrelo  space for deg+1 doubles for the low doubles of the real parts
 *             of the exp(-x) series;
 *   minximhi  space for deg+1 doubles for the high doubles of the imaginary
 *             parts of the exp(-x) series.
 *   minximmi  space for deg+1 doubles for the middle doubles of the imaginary
 *             parts of the exp(-x) series.
 *   minximlo  space for deg+1 doubles for the low doubles of the imaginary
 *             parts of the exp(-x) series.
 *
 * ON RETURN :
 *   pluxrehi  high doubles of the real parts of the power series of exp(+x)
 *             truncated to degree deg;
 *   pluxremi  middle doubles of the real parts of the power series of exp(+x)
 *             truncated to degree deg;
 *   pluxrelo  low doubles of the real parts of the power series of exp(+x)
 *             truncated to degree deg;
 *   pluximhi  high doubles of the imaginary parts of the power series 
 *             of exp(+x) truncated to degree deg;
 *   pluximmi  middle doubles of the imaginary parts of the power series 
 *             of exp(+x) truncated to degree deg;
 *   pluximlo  low doubles of the imaginary parts of the power series 
 *             of exp(+x) truncated to degree deg;
 *   minxrehi  high doubles of the real parts of the power series of exp(-x)
 *             truncated to degree deg;
 *   minxremi  middle doubles of the real parts of the power series of exp(-x)
 *             truncated to degree deg;
 *   minxrelo  low doubles of the real parts of the power series of exp(-x)
 *             truncated to degree deg;
 *   minximhi  high doubles of the imaginary parts of the power series
 *             of exp(-x) truncated to degree deg;
 *   minximmi  middle doubles of the imaginary parts of the power series
 *             of exp(-x) truncated to degree deg;
 *   minximlo  low doubles of the imaginary parts of the power series
 *             of exp(-x) truncated to degree deg. */

void random_cmplx3_exponentials
 ( int deg, double *xrehi, double *xremi, double *xrelo,
            double *ximhi, double *ximmi, double *ximlo,
   double *pluxrehi, double *pluxremi, double *pluxrelo,
   double *pluximhi, double *pluximmi, double *pluximlo,
   double *minxrehi, double *minxremi, double *minxrelo,
   double *minximhi, double *minximmi, double *minximlo );
/*
 * DESCRIPTION :
 *   Returns power series for exp(x) and exp(-x),
 *   truncated to degree deg for a random complex triple double x.
 *   Parameters are the same as cmplx3_exponentials,
 *   except that xrehi, xremi, xrelo, ximhi, ximmi, and ximlo
 *   are return parameters. */

#endif
