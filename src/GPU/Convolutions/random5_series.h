// The file random5_series.h specifies functions to generate random series
// in penta double precision.

#ifndef __random5_series_h__
#define __random5_series_h__

void dbl5_exponentials
 ( int deg, double xtb, double xix, double xmi, double xrg, double xpk, 
   double *pluxtb, double *pluxix, double *pluxmi, double *pluxrg,
   double *pluxpk, double *minxtb, double *minxix, double *minxmi,
   double *minxrg, double *minxpk );
/*
 * DESCRIPTION :
 *   Returns power series for exp(x) and exp(-x),
 *   with real coefficients in penta double precision.
 *
 * ON ENTRY :
 *   deg      degree to truncate the series;
 *   xtb      highest part of some penta double;
 *   xix      second highest part of some penta double;
 *   xmi      middle part of some penta double;
 *   xrg      second lowest part of some penta double;
 *   xpk      lowest part of some penta double;
 *   pluxtb   space for deg+1 doubles for the highest doubles
 *            of the exp(+x) series;
 *   pluxix   space for deg+1 doubles for the second highest doubles
 *            of the exp(+x) series;
 *   pluxmi   space for deg+1 doubles for the middle doubles
 *            of the exp(+x) series;
 *   pluxrg   space for deg+1 doubles for the second lowest doubles
 *            of the exp(+x) series;
 *   pluxpk   space for deg+1 doubles for the lowest doubles
 *            of the exp(+x) series;
 *   minxtb   space for deg+1 doubles for the highest doubles
 *            of the exp(-x) series;
 *   minxix   space for deg+1 doubles for the second highest doubles
 *            of the exp(-x) series;
 *   minxmi   space for deg+1 doubles for the middle doubles
 *            of the exp(-x) series;
 *   minxrg   space for deg+1 doubles for the second lowest doubles
 *            of the exp(-x) series;
 *   minxpk   space for deg+1 doubles for the lowest doubles
 *            of the exp(-x) series.
 *
 * ON RETURN :
 *   pluxtb   highest doubles of the series of exp(+x);
 *   pluxix   second highest doubles of the series of exp(+x);
 *   pluxmi   middle doubles of the series of exp(+x);
 *   pluxrg   second lowest doubles of the series of exp(+x);
 *   pluxpk   lowest doubles of the series of exp(+x);
 *   minxtb   highest doubles of the series of exp(-x);
 *   minxix   second highest doubles of the series of exp(-x);
 *   minxmi   middle doubles of the series of exp(-x);
 *   minxrg   second lowest doubles of the series of exp(-x);
 *   minxpk   lowest doubles of the series of exp(-x). */

void random_dbl5_exponentials
 ( int deg, double *xtb, double *xix, double *xmi, double *xrg, double *xpk,
   double *pluxtb, double *pluxix, double *pluxmi, double *pluxrg,
   double *pluxpk, double *minxtb, double *minxix, double *minxmi,
   double *minxrg, double *minxpk );
/*
 * DESCRIPTION :
 *   Returns power series for exp(x) and exp(-x),
 *   truncated to degree deg for a random penta double x.
 *   Parameters are the same as dbl5_exponentials,
 *   except that xtb, xix, xmi, xrg, and xpk are return parameters. */

void cmplx5_exponentials
 ( int deg,
   double xretb, double xreix, double xremi, double xrerg, double xrepk,
   double ximtb, double ximix, double ximmi, double ximrg, double ximpk,
   double *pluxretb, double *pluxreix, double *pluxremi, double *pluxrerg,
   double *pluxrepk, double *pluximtb, double *pluximix, double *pluximmi,
   double *pluximrg, double *pluximpk,
   double *minxretb, double *minxreix, double *minxremi, double *minxrerg,
   double *minxrepk, double *minximtb, double *minximix, double *minximmi,
   double *minximrg, double *minximpk );
/*
 * DESCRIPTION :
 *   Returns power series for exp(x) and exp(-x),
 *   with complex coefficients in penta double precision.
 *
 * ON ENTRY :
 *   deg       degree to truncate the series;
 *   xretb     high double of the real part of some complex number;
 *   xremi     middle double of the real part of some complex number;
 *   xrelo     low double of the real part of some complex number;
 *   ximtb     high double of the imaginary part of some complex number;
 *   ximmi     middle double of the imaginary part of some complex number;
 *   ximpk     low double of the imaginary part of some complex number;
 *   pluxretb  space for deg+1 doubles for the high doubles of the real parts
 *             of the exp(+x) series;
 *   pluxremi  space for deg+1 doubles for the middle doubles of the real
 *             parts of the exp(+x) series;
 *   pluxrelo  space for deg+1 doubles for the low doubles of the real parts
 *             of the exp(+x) series;
 *   pluximtb  space for deg+1 doubles for the high doubles of the imaginary
 *             parts of the exp(+x) series;
 *   pluximmi  space for deg+1 doubles for the middle doubles of the imaginary
 *             parts of the exp(+x) series;
 *   pluximpk  space for deg+1 doubles for the low doubles of the imaginary
 *             parts of the exp(+x) series;
 *   minxretb  space for deg+1 doubles for the high doubles of the real parts
 *             of the exp(-x) series;
 *   minxremi  space for deg+1 doubles for the middle doubles of the real
 *             parts of the exp(-x) series;
 *   minxrelo  space for deg+1 doubles for the low doubles of the real parts
 *             of the exp(-x) series;
 *   minximtb  space for deg+1 doubles for the high doubles of the imaginary
 *             parts of the exp(-x) series.
 *   minximmi  space for deg+1 doubles for the middle doubles of the imaginary
 *             parts of the exp(-x) series.
 *   minximpk  space for deg+1 doubles for the low doubles of the imaginary
 *             parts of the exp(-x) series.
 *
 * ON RETURN :
 *   pluxretb  high doubles of the real parts of the power series of exp(+x)
 *             truncated to degree deg;
 *   pluxremi  middle doubles of the real parts of the power series of exp(+x)
 *             truncated to degree deg;
 *   pluxrelo  low doubles of the real parts of the power series of exp(+x)
 *             truncated to degree deg;
 *   pluximtb  high doubles of the imaginary parts of the power series 
 *             of exp(+x) truncated to degree deg;
 *   pluximmi  middle doubles of the imaginary parts of the power series 
 *             of exp(+x) truncated to degree deg;
 *   pluximpk  low doubles of the imaginary parts of the power series 
 *             of exp(+x) truncated to degree deg;
 *   minxretb  high doubles of the real parts of the power series of exp(-x)
 *             truncated to degree deg;
 *   minxremi  middle doubles of the real parts of the power series of exp(-x)
 *             truncated to degree deg;
 *   minxrelo  low doubles of the real parts of the power series of exp(-x)
 *             truncated to degree deg;
 *   minximtb  high doubles of the imaginary parts of the power series
 *             of exp(-x) truncated to degree deg;
 *   minximmi  middle doubles of the imaginary parts of the power series
 *             of exp(-x) truncated to degree deg;
 *   minximpk  low doubles of the imaginary parts of the power series
 *             of exp(-x) truncated to degree deg. */

void random_cmplx5_exponentials
 ( int deg,
   double *xretb, double *xreix, double *xremi, double *xrerg, double *xrepk,
   double *ximtb, double *ximix, double *ximmi, double *ximrg, double *ximpk,
   double *pluxretb, double *pluxreix, double *pluxremi, double *pluxrerg,
   double *pluxrepk, double *pluximtb, double *pluximix, double *pluximmi,
   double *pluximrg, double *pluximpk,
   double *minxretb, double *minxreix, double *minxremi, double *minxrerg,
   double *minxrepk, double *minximtb, double *minximix, double *minximmi,
   double *minximrg, double *minximpk );
/*
 * DESCRIPTION :
 *   Returns power series for exp(x) and exp(-x),
 *   truncated to degree deg for a random complex penta double x.
 *   Parameters are the same as cmplx5_exponentials,
 *   except that xretb, xreix, xremi, xrerg, xrepk, ximtb, ximix, ximmi,
 *   ximrg, and ximpk are return parameters. */

#endif
