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
 *   xretb     highest double of the real part of a complex number;
 *   xreix     second highest double of the real part of a complex number;
 *   xremi     middle double of the real part of a complex number;
 *   xrerg     second lowest double of the real part of a complex number;
 *   xrelo     lowest double of the real part of a complex number;
 *   ximtb     highest double of the imaginary part of a complex number;
 *   ximix     second highest double of the imaginary part of a complex number;
 *   ximmi     middle double of the imaginary part of a complex number;
 *   ximrg     second lowest double of the imaginary part of a complex number;
 *   ximpk     lowest double of the imaginary part of a complex number;
 *   pluxretb  space for deg+1 doubles for the highest doubles of the real
 *             parts of the exp(+x) series;
 *   pluxreix  space for deg+1 doubles for the second highest doubles of
 *             the real parts of the exp(+x) series;
 *   pluxremi  space for deg+1 doubles for the middle doubles of the real
 *             parts of the exp(+x) series;
 *   pluxrerg  space for deg+1 doubles for the second lowest doubles of the
 *             real parts of the exp(+x) series;
 *   pluxrelo  space for deg+1 doubles for the lowest doubles of the real
 *             parts of the exp(+x) series;
 *   pluximtb  space for deg+1 doubles for the highest doubles of the
 *             imaginary parts of the exp(+x) series;
 *   pluximix  space for deg+1 doubles for the second highest doubles of the
 *             imaginary parts of the exp(+x) series;
 *   pluximmi  space for deg+1 doubles for the middle doubles of the imaginary
 *             parts of the exp(+x) series;
 *   pluximrg  space for deg+1 doubles for the second lowest doubles of the
 *             imaginary parts of the exp(+x) series;
 *   pluximpk  space for deg+1 doubles for the lowest doubles of the imaginary
 *             parts of the exp(+x) series;
 *   minxretb  space for deg+1 doubles for the highest doubles of the real
 *             parts of the exp(-x) series;
 *   minxreix  space for deg+1 doubles for the second highest doubles of the
 *             real parts of the exp(-x) series;
 *   minxremi  space for deg+1 doubles for the middle doubles of the real
 *             parts of the exp(-x) series;
 *   minxrerg  space for deg+1 doubles for the second lowest doubles of the
 *             real parts of the exp(-x) series;
 *   minxrelo  space for deg+1 doubles for the lowest doubles of the real
 *             parts of the exp(-x) series;
 *   minximtb  space for deg+1 doubles for the highest doubles of the
 *             imaginary parts of the exp(-x) series.
 *   minximix  space for deg+1 doubles for the second highest doubles of the
 *             imaginary parts of the exp(-x) series.
 *   minximmi  space for deg+1 doubles for the middle doubles of the imaginary
 *             parts of the exp(-x) series.
 *   minximrg  space for deg+1 doubles for the second lowest doubles of the
 *             imaginary parts of the exp(-x) series.
 *   minximpk  space for deg+1 doubles for the lowest doubles of the imaginary
 *             parts of the exp(-x) series.
 *
 * ON RETURN :
 *   pluxretb  highest doubles of the real parts of the series of exp(+x)
 *             truncated to degree deg;
 *   pluxreix  second highest doubles of the real parts of the series of
 *             exp(+x) truncated to degree deg;
 *   pluxremi  middle doubles of the real parts of the series of exp(+x)
 *             truncated to degree deg;
 *   pluxrerg  second lowest doubles of the real parts of the series of
 *             exp(+x) truncated to degree deg;
 *   pluxrepk  lowest doubles of the real parts of the series of exp(+x)
 *             truncated to degree deg;
 *   pluximtb  highest doubles of the imaginary parts of the series 
 *             of exp(+x) truncated to degree deg;
 *   pluximix  second highest doubles of the imaginary parts of the series 
 *             of exp(+x) truncated to degree deg;
 *   pluximmi  middle doubles of the imaginary parts of the series 
 *             of exp(+x) truncated to degree deg;
 *   pluximrg  second lowest doubles of the imaginary parts of the series 
 *             of exp(+x) truncated to degree deg;
 *   pluximpk  lowest doubles of the imaginary parts of the series 
 *             of exp(+x) truncated to degree deg;
 *   minxretb  highest doubles of the real parts of the series of exp(-x)
 *             truncated to degree deg;
 *   minxreix  second highest doubles of the real parts of the series of
 *             exp(-x) truncated to degree deg;
 *   minxremi  middle doubles of the real parts of the series of exp(-x)
 *             truncated to degree deg;
 *   minxrerg  second lowest doubles of the real parts of the series of
 *             exp(-x) truncated to degree deg;
 *   minxrepk  lowest doubles of the real parts of the series of exp(-x)
 *             truncated to degree deg;
 *   minximtb  highest doubles of the imaginary parts of the series
 *             of exp(-x) truncated to degree deg;
 *   minximix  second highest doubles of the imaginary parts of the series
 *             of exp(-x) truncated to degree deg;
 *   minximmi  middle doubles of the imaginary parts of the series
 *             of exp(-x) truncated to degree deg;
 *   minximpk  second lowest doubles of the imaginary parts of the series
 *             of exp(-x) truncated to degree deg;
 *   minximpk  lowest doubles of the imaginary parts of the series
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
