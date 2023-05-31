// The file random8_series.h specifies functions to generate random series
// in octo double precision.

#ifndef __random8_series_h__
#define __random8_series_h__

void dbl8_exponential
 ( int deg, double xhihihi, double xlohihi, double xhilohi, double xlolohi,
            double xhihilo, double xlohilo, double xhilolo, double xlololo,
   double *shihihi, double *slohihi, double *shilohi, double *slolohi,
   double *shihilo, double *slohilo, double *shilolo, double *slololo );
/*
 * DESCRIPTION :
 *   Returns power series truncated at degree deg,
 *   following the expansion of exp(x).
 *
 * ON ENTRY :
 *   deg      degree to truncate the series;
 *   xhihihi  highest double of x;
 *   xlohihi  second highest double of x;
 *   xhilohi  third highest double of x;
 *   xlolohi  fourth highest double of x;
 *   xhihilo  fourth lowest double of x;
 *   xlohilo  third lowest double of x;
 *   xhilolo  second lowest double of x;
 *   xlololo  lowest double of x;
 *   shihihi  space for deg+1 doubles;
 *   slohihi  space for deg+1 doubles;
 *   shilohi  space for deg+1 doubles;
 *   slolohi  space for deg+1 doubles;
 *   shihilo  space for deg+1 doubles;
 *   slohilo  space for deg+1 doubles;
 *   shilolo  space for deg+1 doubles;
 *   slololo  space for deg+1 doubles.
 *
 * ON RETURN :
 *   shihihi  highest doubles of the coefficients of series
 *            following exp(x), truncated to degree deg;
 *   slohihi  second highest doubles of the coefficients of series
 *            following exp(x), truncated to degree deg;
 *   shilohi  third highest doubles of the coefficients of series
 *            following exp(x), truncated to degree deg;
 *   slolohi  fourth highest doubles of the coefficients of series
 *            following exp(x), truncated to degree deg;
 *   shihilo  fourth lowest doubles of the coefficients of series
 *            following exp(x), truncated to degree deg;
 *   slohilo  third lowest doubles of the coefficients of series
 *            following exp(x), truncated to degree deg;
 *   shilolo  second lowest doubles of the coefficients of series
 *            following exp(x), truncated to degree deg;
 *   slololo  lowest doubles of the coefficients of series
 *            following exp(x), truncated to degree deg. */

void dbl8_exponentials
 ( int deg, double xhihihi, double xlohihi, double xhilohi, double xlolohi, 
            double xhihilo, double xlohilo, double xhilolo, double xlololo, 
   double *pluxhihihi, double *pluxlohihi, double *pluxhilohi,
   double *pluxlolohi, double *pluxhihilo, double *pluxlohilo,
   double *pluxhilolo, double *pluxlololo,
   double *minxhihihi, double *minxlohihi, double *minxhilohi,
   double *minxlolohi, double *minxhihilo, double *minxlohilo,
   double *minxhilolo, double *minxlololo );
/*
 * DESCRIPTION :
 *   Returns power series for exp(x) and exp(-x),
 *   with real coefficients in octo double precision.
 *
 * ON ENTRY :
 *   deg          degree to truncate the series;
 *   xhihihi      highest part of some octo double;
 *   xlohihi      second highest part of some octo double;
 *   xhilohi      third highest part of some octo double;
 *   xlolohi      fourth highest part of some octo double;
 *   xhihilo      fourth lowest part of some octo double;
 *   xlohilo      third lowest part of some octo double;
 *   xhilolo      second lowest part of some octo double;
 *   xlololo      lowest part of some octo double;
 *   pluxhihihi   space for deg+1 doubles for the highest doubles
 *                of the exp(+x) series;
 *   pluxlohihi   space for deg+1 doubles for the second highest doubles
 *                of the exp(+x) series;
 *   pluxhilohi   space for deg+1 doubles for the third highest doubles
 *                of the exp(+x) series;
 *   pluxlolohi   space for deg+1 doubles for the fourth highest doubles
 *                of the exp(+x) series;
 *   pluxhihilo   space for deg+1 doubles for the fourth lowest doubles
 *                of the exp(+x) series;
 *   pluxlohilo   space for deg+1 doubles for the third lowest doubles
 *                of the exp(+x) series;
 *   pluxhilolo   space for deg+1 doubles for the second lowest doubles
 *                of the exp(+x) series;
 *   pluxlololo   space for deg+1 doubles for the lowest doubles
 *                of the exp(+x) series;
 *   minxhihihi   space for deg+1 doubles for the highest doubles
 *                of the exp(-x) series;
 *   minxlohihi   space for deg+1 doubles for the second highest doubles
 *                of the exp(-x) series;
 *   minxhilohi   space for deg+1 doubles for the third highest doubles
 *                of the exp(-x) series;
 *   minxlolohi   space for deg+1 doubles for the fourth highest doubles
 *                of the exp(-x) series;
 *   minxhihilo   space for deg+1 doubles for the fourth lowest doubles
 *                of the exp(-x) series.
 *   minxlohilo   space for deg+1 doubles for the third lowest doubles
 *                of the exp(-x) series.
 *   minxhilolo   space for deg+1 doubles for the second lowest doubles
 *                of the exp(-x) series.
 *   minxlololo   space for deg+1 doubles for the lowest doubles
 *                of the exp(-x) series.
 *
 * ON RETURN :
 *   pluxhihihi   highest doubles of the series of exp(+x);
 *   pluxlohihi   second highest doubles of the series of exp(+x);
 *   pluxhilohi   third highest doubles of the series of exp(+x);
 *   pluxlolohi   fourth highest doubles of the series of exp(+x);
 *   pluxhihilo   fourth lowest doubles of the series of exp(+x);
 *   pluxlohilo   third lowest doubles of the series of exp(+x);
 *   pluxhilolo   second lowest doubles of the series of exp(+x);
 *   pluxlololo   lowest doubles of the series of exp(+x);
 *   minxhihihi   highest doubles of the series of exp(-x);
 *   minxlohihi   second highest of the series of exp(-x);
 *   minxhilohi   third highest of the series of exp(-x);
 *   minxlolohi   fourth highest of the series of exp(-x);
 *   minxhihilo   fourth lowest doubles of the series of exp(-x);
 *   minxlohilo   third lowest doubles of the series of exp(-x);
 *   minxhilolo   second lowest doubles of the series of exp(-x);
 *   minxlololo   lowest doubles of the series of exp(-x). */

void random_dbl8_exponential
 ( int deg,
   double *xhihihi, double *xlohihi, double *xhilohi, double *xlolohi,
   double *xhihilo, double *xlohilo, double *xhilolo, double *xlololo,
   double *shihihi, double *slohihi, double *shilohi, double *slolohi,
   double *shihilo, double *slohilo, double *shilolo, double *slololo );
/*
 * DESCRIPTION :
 *   Returns a power series truncated at degree deg
 *   for a randomly generated octo double x,
 *   following the expansion of exp(x). */

void random_dbl8_exponentials
 ( int deg,
   double *xhihihi, double *xlohihi, double *xhilohi, double *xlolohi,
   double *xhihilo, double *xlohilo, double *xhilolo, double *xlololo,
   double *pluxhihihi, double *pluxlohihi, double *pluxhilohi,
   double *pluxlolohi, double *pluxhihilo, double *pluxlohilo,
   double *pluxhilolo, double *pluxlololo,
   double *minxhihihi, double *minxlohihi, double *minxhilohi,
   double *minxlolohi, double *minxhihilo, double *minxlohilo,
   double *minxhilolo, double *minxlololo );
/*
 * DESCRIPTION :
 *   Returns power series for exp(x) and exp(-x),
 *   truncated to degree deg for a random octo double x.
 *   Parameters are the same as dbl8_exponentials,
 *   except that xhihihi, xlohihi, xhilohi, xlolohi,
 *   xhihilo, xlohilo, xhilolo, and xlololo are return parameters. */

void cmplx8_exponential
 ( int deg,
   double xrehihihi, double xrelohihi, double xrehilohi, double xrelolohi, 
   double xrehihilo, double xrelohilo, double xrehilolo, double xrelololo, 
   double ximhihihi, double ximlohihi, double ximhilohi, double ximlolohi,
   double ximhihilo, double ximlohilo, double ximhilolo, double ximlololo,
   double *srehihihi, double *srelohihi, double *srehilohi, double *srelolohi,
   double *srehihilo, double *srelohilo, double *srehilolo, double *srelololo,
   double *simhihihi, double *simlohihi, double *simhilohi, double *simlolohi,
   double *simhihilo, double *simlohilo, double *simhilolo, double *simlololo );
/*
 * DESCRIPTION :
 *   Returns power series truncated at degree deg,
 *   following the expansion of exp(x).
 *
 * ON ENTRY :
 *   deg      degree to truncate the series;
 *   xrehihihi is the highest double of the real part of x;
 *   xrelohihi is the second highest double of the real part of x;
 *   xrehilohi is the third highest double of the real part of x;
 *   xrelolohi is the fourth highest double of the real part of x;
 *   xrehihilo is the fourth lowest double of the real part of x;
 *   xrelohilo is the third lowest double of the real part of x;
 *   xrehilolo is the second lowest double of the real part of x;
 *   xrelololo is the lowest double of the real part of x;
 *   ximhihihi is the highest double of the imaginary part of x;
 *   ximlohihi is the second highest double of the imaginary part of x;
 *   ximhilohi is the third highest double of the imaginary part of x;
 *   ximlolohi is the fourth highest double of the imaginary part of x;
 *   ximhihilo is the fourth lowest double of the imaginary part of x;
 *   ximlohilo is the third lowest double of the imaginary part of x;
 *   ximhilolo is the second lowest double of the imaginary part of x;
 *   ximlololo is the lowest double of the imaginary part of x;
 *   srehihihi has space for deg+1 doubles;
 *   srelohihi has space for deg+1 doubles;
 *   srehilohi has space for deg+1 doubles;
 *   srelolohi has space for deg+1 doubles;
 *   srehihilo has space for deg+1 doubles;
 *   srelohilo has space for deg+1 doubles;
 *   srehilolo has space for deg+1 doubles;
 *   srelololo has space for deg+1 doubles;
 *   simhihihi has space for deg+1 doubles;
 *   simlohihi has space for deg+1 doubles;
 *   simhilohi has space for deg+1 doubles;
 *   simlolohi has space for deg+1 doubles.
 *   simhihilo has space for deg+1 doubles;
 *   simlohilo has space for deg+1 doubles;
 *   simhilolo has space for deg+1 doubles;
 *   simlololo has space for deg+1 doubles.
 *
 * ON RETURN :
 *   srehihihi are the highest doubles of the real parts of series
 *            following exp(x), truncated to degree deg;
 *   srelohihi are the second highest doubles of the real parts of series
 *            following exp(x), truncated to degree deg;
 *   srehilohi are the third highest doubles of the real parts of series
 *            following exp(x), truncated to degree deg;
 *   srelolohi are the fourth highest doubles of the real parts of series
 *            following exp(x), truncated to degree deg;
 *   srehihilo are the fourth lowest doubles of the real parts of series
 *            following exp(x), truncated to degree deg;
 *   srelohilo are the third lowest doubles of the real parts of series
 *            following exp(x), truncated to degree deg;
 *   srehilolo are the second lowest doubles of the real parts of series
 *            following exp(x), truncated to degree deg;
 *   srelololo are the lowest doubles of the real parts of series
 *            following exp(x), truncated to degree deg;
 *   simhihihi are the highest doubles of the imaginary parts of series
 *            following exp(x), truncated to degree deg;
 *   simlohihi are the second highest doubles of the imaginary parts of series
 *            following exp(x), truncated to degree deg;
 *   simhilohi are the third highest doubles of the imaginary parts of series
 *            following exp(x), truncated to degree deg;
 *   simlolohi are the fourth highest doubles of the imaginary parts of series
 *            following exp(x), truncated to degree deg;
 *   simhihilo are the fourth lowest doubles of the imaginary parts of series
 *            following exp(x), truncated to degree deg;
 *   simlohilo are the third lowest doubles of the imaginary parts of series
 *            following exp(x), truncated to degree deg;
 *   simhilolo are the second lowest doubles of the imaginary parts of series
 *            following exp(x), truncated to degree deg;
 *   simlololo are the lowest doubles of the imaginary parts of series
 *            following exp(x), truncated to degree deg. */

void cmplx8_exponentials
 ( int deg,
   double xrehihihi, double xrelohihi, double xrehilohi, double xrelolohi,
   double xrehihilo, double xrelohilo, double xrehilolo, double xrelololo,
   double ximhihihi, double ximlohihi, double ximhilohi, double ximlolohi,
   double ximhihilo, double ximlohilo, double ximhilolo, double ximlololo,
   double *pluxrehihihi, double *pluxrelohihi,
   double *pluxrehilohi, double *pluxrelolohi,
   double *pluxrehihilo, double *pluxrelohilo,
   double *pluxrehilolo, double *pluxrelololo,
   double *pluximhihihi, double *pluximlohihi,
   double *pluximhilohi, double *pluximlolohi,
   double *pluximhihilo, double *pluximlohilo,
   double *pluximhilolo, double *pluximlololo,
   double *minxrehihihi, double *minxrelohihi,
   double *minxrehilohi, double *minxrelolohi,
   double *minxrehihilo, double *minxrelohilo,
   double *minxrehilolo, double *minxrelololo,
   double *minximhihihi, double *minximlohihi,
   double *minximhilohi, double *minximlolohi,
   double *minximhihilo, double *minximlohilo,
   double *minximhilolo, double *minximlololo );
/*
 * DESCRIPTION :
 *   Returns power series for exp(x) and exp(-x),
 *   with complex coefficients in octo double precision.
 *
 * ON ENTRY :
 *   deg           degree to truncate the series;
 *   xrehihihi     highest double of the real part of a complex number;
 *   xrelohihi     second highest double of the real part of a complex number;
 *   xrehilohi     third highest double of the real part of a complex number;
 *   xrelolohi     fourth highest double of the real part of a complex number;
 *   xrehihilo     second lowest double of the real part of a complex number;
 *   xrelohilo     second lowest double of the real part of a complex number;
 *   xrehilolo     second lowest double of the real part of a complex number;
 *   xrelololo     lowest double of the real part of a complex number;
 *   ximhihihi     highest double of the imaginary part of a complex number;
 *   ximlohihi     second highest double of the imaginary part
 *                 of a complex number;
 *   ximhilohi     third highest double of the imaginary part
 *                 of a complex number;
 *   ximlolohi     fourth highest double of the imaginary part
 *                 of a complex number;
 *   ximhihilo     fourth lowest double of the imaginary part
 *                 of a complex number;
 *   ximlohilo     third lowest double of the imaginary part
 *                 of a complex number;
 *   ximhilolo     second lowest double of the imaginary part
 *                 of a complex number;
 *   ximlololo     lowest double of the imaginary part of a complex number;
 *   pluxrehihihi  space for deg+1 doubles for the highest doubles
 *                 of the real parts of the exp(+x) series;
 *   pluxrelohihi  space for deg+1 doubles for the second highest doubles
 *                 of the real parts of the exp(+x) series;
 *   pluxrehilohi  space for deg+1 doubles for the third highest doubles
 *                 of the real parts of the exp(+x) series;
 *   pluxrelolohi  space for deg+1 doubles for the fourth highest doubles
 *                 of the real parts of the exp(+x) series;
 *   pluxrehihilo  space for deg+1 doubles for the fourth lowest doubles
 *                 of the real parts of the exp(+x) series;
 *   pluxrelohilo  space for deg+1 doubles for the third lowest doubles
 *                 of the real parts of the exp(+x) series;
 *   pluxrehilolo  space for deg+1 doubles for the second lowest doubles
 *                 of the real parts of the exp(+x) series;
 *   pluxrelololo  space for deg+1 doubles for the lowest doubles
 *                 of the real parts of the exp(+x) series;
 *   pluximhihihi  space for deg+1 doubles for the highest doubles
 *                 of the imaginary parts of the exp(+x) series;
 *   pluximlohihi  space for deg+1 doubles for the second highest doubles
 *                 of the imaginary parts of the exp(+x) series;
 *   pluximhilohi  space for deg+1 doubles for the third highest doubles
 *                 of the imaginary parts of the exp(+x) series;
 *   pluximlolohi  space for deg+1 doubles for the fourth highest doubles
 *                 of the imaginary parts of the exp(+x) series;
 *   pluximhihilo  space for deg+1 doubles for the fourth lowest doubles
 *                 of the imaginary parts of the exp(+x) series;
 *   pluximlohilo  space for deg+1 doubles for the third lowest doubles
 *                 of the imaginary parts of the exp(+x) series;
 *   pluximhilolo  space for deg+1 doubles for the second lowest doubles
 *                 of the imaginary parts of the exp(+x) series;
 *   pluximlololo  space for deg+1 doubles for the lowest doubles
 *                 of the imaginary parts of the exp(+x) series;
 *   minxrehihihi  space for deg+1 doubles for the highest doubles
 *                 of the real parts of the exp(-x) series;
 *   minxrelohihi  space for deg+1 doubles for the second highest doubles
 *                 of the real parts of the exp(-x) series;
 *   minxrehilohi  space for deg+1 doubles for the third highest doubles
 *                 of the real parts of the exp(-x) series;
 *   minxrelolohi  space for deg+1 doubles for the fourth highest doubles
 *                 of the real parts of the exp(-x) series;
 *   minxrehihilo  space for deg+1 doubles for the fourth lowest doubles
 *                 of the real parts of the exp(-x) series;
 *   minxrelohilo  space for deg+1 doubles for the third lowest doubles
 *                 of the real parts of the exp(-x) series;
 *   minxrehilolo  space for deg+1 doubles for the second lowest doubles
 *                 of the real parts of the exp(-x) series;
 *   minxrelololo  space for deg+1 doubles for the lowest doubles
 *                 of the real parts of the exp(-x) series;
 *   minximhihihi  space for deg+1 doubles for the highest doubles
 *                 of the imaginary parts of the exp(-x) series.
 *   minximlohihi  space for deg+1 doubles for the second highest doubles
 *                 of the imaginary parts of the exp(-x) series.
 *   minximhilohi  space for deg+1 doubles for the third highest doubles
 *                 of the imaginary parts of the exp(-x) series.
 *   minximlolohi  space for deg+1 doubles for the fourth highest doubles
 *                 of the imaginary parts of the exp(-x) series.
 *   minximhihilo  space for deg+1 doubles for the fourth lowest doubles
 *                 of the imaginary parts of the exp(-x) series.
 *   minximlohilo  space for deg+1 doubles for the third lowest doubles
 *                 of the imaginary parts of the exp(-x) series.
 *   minximhilolo  space for deg+1 doubles for the second lowest doubles
 *                 of the imaginary parts of the exp(-x) series.
 *   minximlololo  space for deg+1 doubles for the lowest doubles
 *                 of the imaginary parts of the exp(-x) series.
 *
 * ON RETURN :
 *   pluxrehihihi  highest doubles of the real parts of the exp(+x) series
 *                 truncated to degree deg;
 *   pluxrelohihi  second highest doubles of the real parts of the exp(+x)
 *                 series truncated to degree deg;
 *   pluxrehilohi  third highest doubles of the real parts of the exp(+x)
 *                 series truncated to degree deg;
 *   pluxrelolohi  fourth highest doubles of the real parts of the exp(+x)
 *                 series truncated to degree deg;
 *   pluxrehihilo  fourth lowest doubles of the real parts of the exp(+x)
 *                 series truncated to degree deg;
 *   pluxrelohilo  third lowest doubles of the real parts of the exp(+x)
 *                 series truncated to degree deg;
 *   pluxrehilolo  second lowest doubles of the real parts of the exp(+x)
 *                 series truncated to degree deg;
 *   pluxrelololo  lowest doubles of the real parts of the exp(+x) series
 *                 truncated to degree deg;
 *   pluximhihihi  highest doubles of the imaginary parts of the exp(+x)
 *                 series truncated to degree deg;
 *   pluximlohihi  second highest doubles of the imaginary parts of the
 *                 exp(+x) series truncated to degree deg;
 *   pluximhilohi  third highest doubles of the imaginary parts of the exp(+x)
 *                 series truncated to degree deg;
 *   pluximlolohi  fourth highest doubles of the imaginary parts of the
 *                 exp(+x) series truncated to degree deg;
 *   pluximhihilo  fourth lowest doubles of the imaginary parts of the exp(+x)
 *                 series truncated to degree deg;
 *   pluximlohilo  third lowest doubles of the imaginary parts of the exp(+x)
 *                 series truncated to degree deg;
 *   pluximhilolo  second lowest doubles of the imaginary parts of the exp(+x)
 *                 series truncated to degree deg;
 *   pluximlololo  lowest doubles of the imaginary parts of the exp(+x) series 
 *                 truncated to degree deg;
 *   minxrehihihi  highest doubles of the real parts of the exp(-x) series
 *                 truncated to degree deg;
 *   minxrelohihi  second highest doubles of the real parts of the exp(-x)
 *                 series truncated to degree deg;
 *   minxrehilohi  third highest doubles of the real parts of the exp(-x)
 *                 series truncated to degree deg;
 *   minxrelolohi  fourth highest doubles of the real parts of the exp(-x)
 *                 series truncated to degree deg;
 *   minxrehihilo  fourth lowest doubles of the real parts of the exp(-x)
 *                 series truncated to degree deg;
 *   minxrelohilo  third lowest doubles of the real parts of the exp(-x)
 *                 series truncated to degree deg;
 *   minxrehilolo  second lowest doubles of the real parts of the exp(-x)
 *                 series truncated to degree deg;
 *   minxrelololo  lowest doubles of the real parts of the exp(-x) series
 *                 truncated to degree deg;
 *   minximhihihi  highest doubles of the imaginary parts of the exp(-x) 
 *                 series truncated to degree deg;
 *   minximlohihi  second highest doubles of the imaginary parts of the
 *                 exp(-x) series truncated to degree deg;
 *   minximhilohi  third highest doubles of the imaginary parts of the
 *                 exp(-x) series truncated to degree deg;
 *   minximlolohi  fourth highest doubles of the imaginary parts of the
 *                 exp(-x) series truncated to degree deg;
 *   minximhihilo  fourth lowest doubles of the imaginary parts of the exp(-x)
 *                 series truncated to degree deg;
 *   minximlohilo  third lowest doubles of the imaginary parts of the exp(-x)
 *                 series truncated to degree deg;
 *   minximhilolo  second lowest doubles of the imaginary parts of the exp(-x)
 *                 series truncated to degree deg;
 *   minximlololo  lowest doubles of the imaginary parts of the exp(-x) series
 *                 truncated to degree deg. */

void random_cmplx8_exponential
 ( int deg,
   double *xrehihihi, double *xrelohihi, double *xrehilohi, double *xrelolohi,
   double *xrehihilo, double *xrelohilo, double *xrehilolo, double *xrelololo,
   double *ximhihihi, double *ximlohihi, double *ximhilohi, double *ximlolohi,
   double *ximhihilo, double *ximlohilo, double *ximhilolo, double *ximlololo,
   double *srehihihi, double *srelohihi, double *srehilohi, double *srelolohi,
   double *srehihilo, double *srelohilo, double *srehilolo, double *srelololo,
   double *simhihihi, double *simlohihi, double *simhilohi, double *simlolohi,
   double *simhihilo, double *simlohilo, double *simhilolo, double *simlololo );
/*
 * DESCRIPTION :
 *   Returns power series following the expansion of exp(x),
 *   for a complex x, with randomly generated real and imaginary parts
 *   returned in xre and xim, truncated at degree deg. */

void random_cmplx8_exponentials
 ( int deg, 
   double *xrehihihi, double *xrelohihi, double *xrehilohi, double *xrelolohi,
   double *xrehihilo, double *xrelohilo, double *xrehilolo, double *xrelololo,
   double *ximhihihi, double *ximlohihi, double *ximhilohi, double *ximlolohi,
   double *ximhihilo, double *ximlohilo, double *ximhilolo, double *ximlololo,
   double *pluxrehihihi, double *pluxrelohihi,
   double *pluxrehilohi, double *pluxrelolohi,
   double *pluxrehihilo, double *pluxrelohilo,
   double *pluxrehilolo, double *pluxrelololo,
   double *pluximhihihi, double *pluximlohihi,
   double *pluximhilohi, double *pluximlolohi,
   double *pluximhihilo, double *pluximlohilo,
   double *pluximhilolo, double *pluximlololo,
   double *minxrehihihi, double *minxrelohihi,
   double *minxrehilohi, double *minxrelolohi,
   double *minxrehihilo, double *minxrelohilo,
   double *minxrehilolo, double *minxrelololo,
   double *minximhihihi, double *minximlohihi,
   double *minximhilohi, double *minximlolohi,
   double *minximhihilo, double *minximlohilo,
   double *minximhilolo, double *minximlololo );
/*
 * DESCRIPTION :
 *   Returns power series for exp(x) and exp(-x),
 *   truncated to degree deg for a random complex octo double x.
 *   Parameters are the same as cmplx8_exponentials, except that
 *   xrehihihi, xrelohihi, xrehilohi, xrelolohi, xrehihilo, xrelohilo,
 *   xrehilolo, xrelololo, ximhihihi, ximlohihi, ximhilohi, ximlolohi,
 *   ximhihilo, ximlohilo, ximhilolo, ximlololo are return parameters. */

#endif
