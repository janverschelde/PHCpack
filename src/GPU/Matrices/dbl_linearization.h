// The file dbl_linearization.h specifies functions to linearize
// vectors and matrices of series in double_precision.

#ifndef __dbl_linearization_h__
#define __dbl_linearization_h__

void dbl_linear_series_vector
 ( int dim, int deg, double **v, double **w );
/*
 * DESCRIPTION :
 *   Linearization turns a vector of power series
 *   into a series with vectors as coefficients, for real data.
 *
 * ON ENTRY :
 *   dim      dimension of the vectors;
 *   deg      truncation degree of the series;
 *   v        v[k] is a power series truncated at degree deg,
 *            for k ranging from 0 to dim-1;
 *   w        space for deg+1 arrays of dim doubles.
 * 
 * ON RETURN :
 *   w        w[k] is the coefficient vector of the k-th term of the
 *            series, for k ranging from 0 to deg. */

void cmplx_linear_series_vector
 ( int dim, int deg, double **vre, double **vim,
   double **wre, double **wim );
/*
 * DESCRIPTION :
 *   Linearization turns a vector of power series
 *   into a series with vectors as coefficients, for complex data.
 *
 * ON ENTRY :
 *   dim      dimension of the vectors;
 *   deg      truncation degree of the series;
 *   vre      vre[k] stores the real part of a power series
 *            truncated at degree deg, for k ranging from 0 to dim-1;
 *   vim      vre[k] stores the imaginary part of a power series
 *            truncated at degree deg, for k ranging from 0 to dim-1;
 *   wre      space for deg+1 arrays of dim doubles;
 *   wim      space for deg+1 arrays of dim doubles.
 * 
 * ON RETURN :
 *   wre      wre[k] is the real part of the coefficient vector
 *            of the k-th term of the series, for k ranging from 0 to deg;
 *   wim      wim[k] is the imaginary part of the coefficient vector
 *            of the k-th term of the series, for k ranging from 0 to deg. */

#endif
