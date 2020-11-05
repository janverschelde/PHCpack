// The file random10_series.h specifies functions to generate random series
// in deca double precision.

#ifndef __random10_series_h__
#define __random10_series_h__

void dbl10_exponentials
 ( int deg,
   double xrtb, double xrix, double xrmi, double xrrg, double xrpk, 
   double xltb, double xlix, double xlmi, double xlrg, double xlpk, 
   double *pluxrtb, double *pluxrix, double *pluxrmi, double *pluxrrg,
   double *pluxrpk, double *pluxltb, double *pluxlix, double *pluxlmi,
   double *pluxlrg, double *pluxlpk,
   double *minxrtb, double *minxrix, double *minxrmi, double *minxrrg,
   double *minxrpk, double *minxltb, double *minxlix, double *minxlmi,
   double *minxlrg, double *minxlpk );
/*
 * DESCRIPTION :
 *   Returns power series for exp(x) and exp(-x),
 *   with real coefficients in deca double precision.
 *
 * ON ENTRY :
 *   deg      degree to truncate the series;
 *   xrtb      highest part of some deca double;
 *   xrix      second highest part of some deca double;
 *   xrmi      third highest part of some deca double;
 *   xrrg      fourth highest part of some deca double;
 *   xrpk      fifth highest part of some deca double;
 *   xltb      fifth lowest part of some deca double;
 *   xlix      fourth lowest part of some deca double;
 *   xlmi      third lowest part of some deca double;
 *   xlrg      second lowest part of some deca double;
 *   xlpk      lowest part of some deca double;
 *   pluxrtb   space for deg+1 doubles for the highest doubles
 *             of the exp(+x) series;
 *   pluxrix   space for deg+1 doubles for the second highest doubles
 *             of the exp(+x) series;
 *   pluxrmi   space for deg+1 doubles for the third highest doubles
 *             of the exp(+x) series;
 *   pluxrrg   space for deg+1 doubles for the fourth highest doubles
 *             of the exp(+x) series;
 *   pluxrpk   space for deg+1 doubles for the fifth highest doubles
 *             of the exp(+x) series;
 *   pluxltb   space for deg+1 doubles for the fifth lowest doubles
 *             of the exp(+x) series;
 *   pluxlix   space for deg+1 doubles for the fourth lowest doubles
 *             of the exp(+x) series;
 *   pluxlmi   space for deg+1 doubles for the third lowest doubles
 *             of the exp(+x) series;
 *   pluxlrg   space for deg+1 doubles for the second lowest doubles
 *             of the exp(+x) series;
 *   pluxlpk   space for deg+1 doubles for the lowest doubles
 *             of the exp(+x) series;
 *   minxrtb   space for deg+1 doubles for the highest doubles
 *             of the exp(-x) series;
 *   minxrix   space for deg+1 doubles for the second highest doubles
 *             of the exp(-x) series;
 *   minxrmi   space for deg+1 doubles for the third highest doubles
 *             of the exp(-x) series;
 *   minxrrg   space for deg+1 doubles for the fourth highest doubles
 *             of the exp(-x) series;
 *   minxrpk   space for deg+1 doubles for the fifth highest doubles
 *             of the exp(-x) series.
 *   minxltb   space for deg+1 doubles for the fifth lowest doubles
 *             of the exp(-x) series;
 *   minxlix   space for deg+1 doubles for the fourth lowest doubles
 *             of the exp(-x) series;
 *   minxlmi   space for deg+1 doubles for the third lowest doubles
 *             of the exp(-x) series;
 *   minxlrg   space for deg+1 doubles for the second lowest doubles
 *             of the exp(-x) series;
 *   minxlpk   space for deg+1 doubles for the lowest doubles
 *             of the exp(-x) series.
 *
 * ON RETURN :
 *   pluxrtb   highest doubles of the series of exp(+x);
 *   pluxrix   second highest doubles of the series of exp(+x);
 *   pluxrmi   third highest doubles of the series of exp(+x);
 *   pluxrrg   fourth highest doubles of the series of exp(+x);
 *   pluxrpk   fifth highest doubles of the series of exp(+x);
 *   pluxltb   fifth lowest doubles of the series of exp(+x);
 *   pluxlix   fourth lowest doubles of the series of exp(+x);
 *   pluxlmi   third lowest doubles of the series of exp(+x);
 *   pluxlrg   second lowest doubles of the series of exp(+x);
 *   pluxlpk   lowest doubles of the series of exp(+x);
 *   minxrtb   highest doubles of the series of exp(-x);
 *   minxrix   second highest doubles of the series of exp(-x);
 *   minxrmi   third highest doubles of the series of exp(-x);
 *   minxrrg   fourth highest doubles of the series of exp(-x);
 *   minxrpk   fifth highest doubles of the series of exp(-x);
 *   minxltb   fifth lowest doubles of the series of exp(-x);
 *   minxlix   fourth lowest doubles of the series of exp(-x);
 *   minxlmi   third lowest doubles of the series of exp(-x);
 *   minxlrg   second lowest doubles of the series of exp(-x);
 *   minxlpk   lowest doubles of the series of exp(-x). */

void random_dbl10_exponentials
 ( int deg,
   double *xrtb, double *xrix, double *xrmi, double *xrrg, double *xrpk,
   double *xltb, double *xlix, double *xlmi, double *xlrg, double *xlpk,
   double *pluxrtb, double *pluxrix, double *pluxrmi, double *pluxrrg,
   double *pluxrpk, double *pluxltb, double *pluxlix, double *pluxlmi,
   double *pluxlrg, double *pluxlpk,
   double *minxrtb, double *minxrix, double *minxrmi, double *minxrrg,
   double *minxrpk, double *minxltb, double *minxlix, double *minxlmi,
   double *minxlrg, double *minxlpk );
/*
 * DESCRIPTION :
 *   Returns power series for exp(x) and exp(-x),
 *   truncated to degree deg for a random deca double x.
 *   Parameters are the same as dbl10_exponentials, except that
 *   xrtb, xrix, xrmi, xrrg, xrpk, xltb, xlix, xlmi, xlrg, and xlpk
 *   are return parameters. */

void cmplx10_exponentials
 ( int deg,
   double xrertb, double xrerix, double xrermi, double xrerrg, double xrerpk,
   double xreltb, double xrelix, double xrelmi, double xrelrg, double xrelpk,
   double ximrtb, double ximrix, double ximrmi, double ximrrg, double ximrpk,
   double ximltb, double ximlix, double ximlmi, double ximlrg, double ximlpk,
   double *pluxrertb, double *pluxrerix, double *pluxrermi, double *pluxrerrg,
   double *pluxrerpk, double *pluxreltb, double *pluxrelix, double *pluxrelmi,
   double *pluxrelrg, double *pluxrelpk,
   double *pluximrtb, double *pluximrix, double *pluximrmi, double *pluximrrg,
   double *pluximrpk, double *pluximltb, double *pluximlix, double *pluximlmi,
   double *pluximlrg, double *pluximlpk,
   double *minxrertb, double *minxrerix, double *minxrermi, double *minxrerrg,
   double *minxrerpk, double *minxreltb, double *minxrelix, double *minxrelmi,
   double *minxrelrg, double *minxrelpk,
   double *minximrtb, double *minximrix, double *minximrmi, double *minximrrg,
   double *minximrpk, double *minximltb, double *minximlix, double *minximlmi,
   double *minximlrg, double *minximlpk );
/*
 * DESCRIPTION :
 *   Returns power series for exp(x) and exp(-x),
 *   with complex coefficients in deca double precision.
 *
 * ON ENTRY :
 *   deg       degree to truncate the series;
 *   xrertb    highest double of the real part of a complex number;
 *   xrerix    second highest double of the real part of a complex number;
 *   xrermi    third highest double of the real part of a complex number;
 *   xrerrg    fourth highest double of the real part of a complex number;
 *   xrerpk    fifth highest double of the real part of a complex number;
 *   xreltb    fifth lowest double of the real part of a complex number;
 *   xrelix    fourth lowest double of the real part of a complex number;
 *   xrelmi    third lowest double of the real part of a complex number;
 *   xrelrg    second lowest double of the real part of a complex number;
 *   xrelpk    lowest double of the real part of a complex number;
 *   ximrtb    highest double of the imaginary part of a complex number;
 *   ximrix    second highest double of the imaginary part of a complex number;
 *   ximrmi    third highest double of the imaginary part of a complex number;
 *   ximrrg    fourth highest double of the imaginary part of a complex number;
 *   ximrpk    fifth highest double of the imaginary part of a complex number;
 *   ximltb    fifth lowest double of the imaginary part of a complex number;
 *   ximlix    fourth lowest double of the imaginary part of a complex number;
 *   ximlmi    third lowest double of the imaginary part of a complex number;
 *   ximlrg    second lowest double of the imaginary part of a complex number;
 *   ximlpk    lowest double of the imaginary part of a complex number;
 *   pluxrertb holds space for deg+1 doubles for the highest doubles
 *             of the real parts of the exp(+x) series;
 *   pluxrerix holds space for deg+1 doubles for the second highest doubles
 *             of the real parts of the exp(+x) series;
 *   pluxrermi holds space for deg+1 doubles for the third highest doubles
 *             of the real parts of the exp(+x) series;
 *   pluxrerrg holds space for deg+1 doubles for the fourth highest doubles
 *             of the real parts of the exp(+x) series;
 *   pluxrerpk holds space for deg+1 doubles for the fifth highest doubles
 *             of the real parts of the exp(+x) series;
 *   pluxreltb holds space for deg+1 doubles for the fifth lowest doubles
 *             of the real parts of the exp(+x) series;
 *   pluxrelix holds space for deg+1 doubles for the fourth lowest doubles
 *             of the real parts of the exp(+x) series;
 *   pluxrelmi holds space for deg+1 doubles for the third lowest doubles
 *             of the real parts of the exp(+x) series;
 *   pluxrelrg holds space for deg+1 doubles for the second lowest doubles
 *             of the real parts of the exp(+x) series;
 *   pluxrelpk holds space for deg+1 doubles for the lowest doubles
 *             of the real parts of the exp(+x) series;
 *   pluximrtb holds space for deg+1 doubles for the highest doubles
 *             of the imaginary parts of the exp(+x) series;
 *   pluximrix holds space for deg+1 doubles for the second highest doubles
 *             of the imaginary parts of the exp(+x) series;
 *   pluximrmi holds space for deg+1 doubles for the third highest doubles
 *             of the imaginary parts of the exp(+x) series;
 *   pluximrrg holds space for deg+1 doubles for the fourth highest doubles
 *             of the imaginary parts of the exp(+x) series;
 *   pluximrpk holds space for deg+1 doubles for the fifth highest doubles
 *             of the imaginary parts of the exp(+x) series;
 *   pluximltb holds space for deg+1 doubles for the fifth lowest doubles
 *             of the imaginary parts of the exp(+x) series;
 *   pluximlix holds space for deg+1 doubles for the fourth lowest doubles
 *             of the imaginary parts of the exp(+x) series;
 *   pluximlmi holds space for deg+1 doubles for the third lowest doubles
 *             of the imaginary parts of the exp(+x) series;
 *   pluximlrg holds space for deg+1 doubles for the second lowest doubles
 *             of the imaginary parts of the exp(+x) series;
 *   pluximlpk holds space for deg+1 doubles for the lowest doubles
 *             of the imaginary parts of the exp(+x) series;
 *   minxrertb holds space for deg+1 doubles for the highest doubles
 *             of the real parts of the exp(-x) series;
 *   minxrerix holds space for deg+1 doubles for the second highest doubles
 *             of the real parts of the exp(-x) series;
 *   minxrermi holds space for deg+1 doubles for the third highest doubles
 *             of the real parts of the exp(-x) series;
 *   minxrerrg holds space for deg+1 doubles for the fourth highest doubles
 *             of the real parts of the exp(-x) series;
 *   minxrerpk holds space for deg+1 doubles for the fifth highest doubles
 *             of the real parts of the exp(-x) series;
 *   minximltb holds space for deg+1 doubles for the fifth lowest doubles
 *             of the imaginary parts of the exp(-x) series.
 *   minximlix holds space for deg+1 doubles for the fourth lowest doubles
 *             of the imaginary parts of the exp(-x) series.
 *   minximlmi holds space for deg+1 doubles for the third lowest doubles
 *             of the imaginary parts of the exp(-x) series.
 *   minximlrg holds space for deg+1 doubles for the second lowest doubles
 *             of the imaginary parts of the exp(-x) series.
 *   minximlpk holds space for deg+1 doubles for the lowest doubles
 *             of the imaginary parts of the exp(-x) series.
 *
 * ON RETURN :
 *   pluxrertb are the highest doubles of the real parts of the series
 *             of exp(+x) truncated to degree deg;
 *   pluxrerix are the second highest doubles of the real parts of the series
 *             of exp(+x) truncated to degree deg;
 *   pluxrermi are the third highest doubles of the real parts of the series
 *             of exp(+x) truncated to degree deg;
 *   pluxrerrg are the fourth highest doubles of the real parts of the series
 *             of exp(+x) truncated to degree deg;
 *   pluxrerpk are the fifth highest doubles of the real parts of the series
 *             of exp(+x) truncated to degree deg;
 *   pluxreltb are the fifth lowest doubles of the real parts of the series
 *             of exp(+x) truncated to degree deg;
 *   pluxrelix are the fourth lowest doubles of the real parts of the series
 *             of exp(+x) truncated to degree deg;
 *   pluxrelmi are the third lowest doubles of the real parts of the series
 *             of exp(+x) truncated to degree deg;
 *   pluxrelrg are the second lowest doubles of the real parts of the series
 *             of exp(+x) truncated to degree deg;
 *   pluxrelpk are the lowest doubles of the real parts of the series
 *             of exp(+x) truncated to degree deg;
 *   pluximrtb are the highest doubles of the imaginary parts of the series 
 *             of exp(+x) truncated to degree deg;
 *   pluximrix are the second highest doubles of the imaginary parts
 *             of the series of exp(+x) truncated to degree deg;
 *   pluximrmi are the third highest doubles of the imaginary parts
 *             of the series of exp(+x) truncated to degree deg;
 *   pluximrrg are the fourth highest doubles of the imaginary parts
 *             of the series of exp(+x) truncated to degree deg;
 *   pluximrpk are the fifth highest doubles of the imaginary parts
 *             of the series  of exp(+x) truncated to degree deg;
 *   pluximltb are the fifth lowest doubles of the imaginary parts
 *             of the series of exp(+x) truncated to degree deg;
 *   pluximlix are the fourth lowest doubles of the imaginary parts
 *             of the series of exp(+x) truncated to degree deg;
 *   pluximlmi are the third lowest doubles of the imaginary parts
 *             of the series of exp(+x) truncated to degree deg;
 *   pluximlrg are the second lowest doubles of the imaginary parts
 *             of the series of exp(+x) truncated to degree deg;
 *   pluximlpk are the lowest doubles of the imaginary parts of the series 
 *             of exp(+x) truncated to degree deg;
 *   minxrertb are the highest doubles of the real parts of the series
 *             of exp(-x) truncated to degree deg;
 *   minxrerix are the second highest doubles of the real parts of the series
 *             of exp(-x) truncated to degree deg;
 *   minxrermi are the third highest doubles of the real parts of the series
 *             of exp(-x) truncated to degree deg;
 *   minxrerrg are the fourth highest doubles of the real parts of the series
 *             of exp(-x) truncated to degree deg;
 *   minxrerpk are the fifth highest doubles of the real parts of the series
 *             of exp(-x) truncated to degree deg;
 *   minxreltb are the fifth lowest doubles of the real parts of the series
 *             of exp(-x) truncated to degree deg;
 *   minxrelix are the fourth lowest doubles of the real parts of the series
 *             of exp(-x) truncated to degree deg;
 *   minxrelmi are the third lowest doubles of the real parts of the series
 *             of exp(-x) truncated to degree deg;
 *   minxrelrg are the second lowest doubles of the real parts of the series
 *             of exp(-x) truncated to degree deg;
 *   minxrelpk are the lowest doubles of the real parts of the series
 *             of exp(-x) truncated to degree deg;
 *   minximrtb are the highest doubles of the imaginary parts of the series
 *             of exp(-x) truncated to degree deg;
 *   minximrix are the second highest doubles of the imaginary parts
 *             of the series of exp(-x) truncated to degree deg;
 *   minximrmi are the third highest doubles of the imaginary parts
 *             of the series of exp(-x) truncated to degree deg;
 *   minximrrg are the fourth highest doubles of the imaginary parts
 *             of the series of exp(-x) truncated to degree deg;
 *   minximrpk are the fifth highest doubles of the imaginary parts
 *             of the series of exp(-x) truncated to degree deg;
 *   minximltb are the fifth lowest doubles of the imaginary parts
 *             of the series of exp(-x) truncated to degree deg;
 *   minximlix are the fourth lowest doubles of the imaginary parts
 *             of the series of exp(-x) truncated to degree deg;
 *   minximlmi are the third lowest doubles of the imaginary parts
 *             of the series of exp(-x) truncated to degree deg;
 *   minximlrg are the second lowest doubles of the imaginary parts
 *             of the series of exp(-x) truncated to degree deg;
 *   minximlpk are the lowest doubles of the imaginary parts of the series
 *             of exp(-x) truncated to degree deg. */

void random_cmplx10_exponentials
 ( int deg,
   double *xrertb, double *xrerix, double *xrermi, double *xrerrg,
   double *xrerpk, double *xreltb, double *xrelix, double *xrelmi,
   double *xrelrg, double *xrelpk,
   double *ximrtb, double *ximrix, double *ximrmi, double *ximrrg,
   double *ximrpk, double *ximltb, double *ximlix, double *ximlmi,
   double *ximlrg, double *ximlpk,
   double *pluxrertb, double *pluxrerix, double *pluxrermi, double *pluxrerrg,
   double *pluxrerpk, double *pluxreltb, double *pluxrelix, double *pluxrelmi,
   double *pluxrelrg, double *pluxrelpk,
   double *pluximrtb, double *pluximrix, double *pluximrmi, double *pluximrrg,
   double *pluximrpk, double *pluximltb, double *pluximlix, double *pluximlmi,
   double *pluximlrg, double *pluximlpk,
   double *minxrertb, double *minxrerix, double *minxrermi, double *minxrerrg,
   double *minxrerpk, double *minxreltb, double *minxrelix, double *minxrelmi,
   double *minxrelrg, double *minxrelpk,
   double *minximrtb, double *minximrix, double *minximrmi, double *minximrrg,
   double *minximrpk, double *minximltb, double *minximlix, double *minximlmi,
   double *minximlrg, double *minximlpk );
/*
 * DESCRIPTION :
 *   Returns power series for exp(x) and exp(-x),
 *   truncated to degree deg for a random complex deca double x.
 *   Parameters are the same as cmplx10_exponentials, except that
 *   xrertb, xrerix, xrermi, xrerrg, xrerpk, xreltb, xrelix, xrelmi,
 *   xrelrg, xrelpk, ximrtb, ximrix, ximrmi, ximrrg, ximrpk, ximltb,
 *   ximlix, ximlmi, ximlrg, and ximlpk are return parameters. */

#endif
