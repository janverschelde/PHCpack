// The file random_matrices.h specifies functions define tests to work
// with matrices and vectors of series in double precision.

#ifndef __random_matrices_h__
#define __random_matrices_h__

void random_dbl_series_vector ( int dim, int deg, double *x, double **v );
/*
 * DESCRIPTION :
 *   Returns in v a random vector of dimension dim,
 *   of series truncated at degree deg.
 *
 * ON ENTRY :
 *   dim      dimension of the vectors;
 *   deg      truncation degree of the series;
 *   x        space for dim doubles;
 *   v        space for dim arrays of deg+1 doubles.
 * 
 * ON RETURN :
 *   x        dim random doubles;
 *   v        v[k] is a power series of exp(x[k]) at degree deg,
 *            for k ranging from 0 to dim-1. */

void random_dbl_series_vectors
 ( int dim, int deg, double *x, double **plux, double **minx );
/*
 * DESCRIPTION :
 *   Returns in plux and minx a pair of two real vectors with series 
 *   for exp(x[k]) and exp(-x[k]), for k ranging from 0 to dim-1.
 *
 * ON ENTRY :
 *   dim      dimension of the vectors;
 *   deg      truncation degree of the series;
 *   x        space for dim doubles;
 *   plux     space for dim arrays of deg+1 doubles;
 *   minx     space for dim arrays of deg+1 doubles.
 *
 * ON RETURN :
 *   x        dim random doubles;
 *   plux     plux[k] is power series of exp(+x[k]), for k from 0 to dim-1,
 *            truncated to degree deg;
 *   minx     plux[k] is power series of exp(-x[k]), for k from 0 to dim-1,
 *            truncated to degree deg. */

void random_cmplx_series_vectors
 ( int dim, int deg, double *xre, double *xim,
   double **pluxre, double **pluxim, double **minxre, double **minxim );
/*
 * DESCRIPTION :
 *   Returns in plux and minx a pair of two complex vectors with series 
 *   for exp(x[k]) and exp(-x[k]), for k ranging from 0 to dim-1.
 *
 * ON ENTRY :
 *   dim      dimension of the vectors;
 *   deg      truncation degree of the series;
 *   xre      space for dim doubles;
 *   xim      space for dim doubles;
 *   pluxre   space for dim arrays of deg+1 doubles;
 *   pluxim   space for dim arrays of deg+1 doubles;
 *   minxre   space for dim arrays of deg+1 doubles;
 *   minxim   space for dim arrays of deg+1 doubles.
 *
 * ON RETURN :
 *   xre      real parts of dim random complex numbers,
 *   xim      imaginary parts of dim random complex numbers;
 *   pluxre   pluxre[k] is the real part of the series of exp(+x[k]),
 *            for k from 0 to dim-1, truncated to degree deg;
 *   pluxim   pluxim[k] is the imaginary part of the series of exp(+x[k]),
 *            for k from 0 to dim-1, truncated to degree deg;
 *   minxre   pluxre[k] is the real part of the series of exp(-x[k]),
 *            for k from 0 to dim-1, truncated to degree deg;
 *   minxim   pluxim[k] is the imaginary part of the series of exp(-x[k]),
 *            for k from 0 to dim-1, truncated to degree deg. */

void random_dbl_series_matrix
 ( int rows, int cols, int deg, double **x, double ***A );
/*
 * DESCRIPTION :
 *   Returns in A and matrix of dimensions rows and cols,
 *   of power series for exp(x), truncated at degree deg.
 *
 * ON ENTRY :
 *   rows     the number of rows in the matrices x and A;
 *   cols     the number of columns in the matrices x and A;
 *   deg      truncation degree of the series;
 *   x        space for a matrix of doubles, of dimensions rows and cols;
 *   A        space for a matrix of series, of dimensions rows and cols.
 *
 * ON RETURN :
 *   x        a matrix of randomly generated doubles;
 *   A        a matrix of randomly generated power series. */

#endif
