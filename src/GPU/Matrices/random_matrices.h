// The file random_matrices.h specifies functions define tests to work
// with matrices and vectors of series in double precision.

#ifndef __random_matrices_h__
#define __random_matrices_h__

void random_dbl_series_vector
 ( int dim, int deg, double *x, double **v, bool expform=true );
/*
 * DESCRIPTION :
 *   Returns in v a random real vector of dimension dim,
 *   of series truncated at degree deg.
 *
 * ON ENTRY :
 *   dim      dimension of the vectors;
 *   deg      truncation degree of the series;
 *   x        space for dim doubles;
 *   v        space for dim arrays of deg+1 doubles;
 *   expform  if true, then the exp() expansion is applied,
 *            otherwise, the log(1+x) expansion is applied.
 * 
 * ON RETURN :
 *   x        dim random doubles;
 *   v        v[k] is a power series of exp(x[k]) at degree deg,
 *            for k ranging from 0 to dim-1. */

void random_cmplx_series_vector
 ( int dim, int deg, double *xre, double *xim,
   double **vre, double **vim, bool expform=true );
/*
 * DESCRIPTION :
 *   Returns in v a random complex vector of dimension dim,
 *   of series truncated at degree deg.
 *
 * ON ENTRY :
 *   dim      dimension of the vectors;
 *   deg      truncation degree of the series;
 *   xre      space for dim doubles;
 *   xim      space for dim doubles;
 *   vre      space for dim arrays of deg+1 doubles;
 *   vim      space for dim arrays of deg+1 doubles;
 *   expform  if true, then the exp() expansion is applied,
 *            otherwise, the log(1+x) expansion is applied.
 * 
 * ON RETURN :
 *   xre      dim random doubles;
 *   xim      dim random doubles;
 *   vre      vre[k] is the real part of a power series truncated 
 *            at degree deg, for k ranging from 0 to dim-1;
 *   vim      vre[k] is the imaginary part of a power series truncated 
 *            at degree deg, for k ranging from 0 to dim-1. */

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
 ( int rows, int cols, int deg, double **x, double ***A, bool expform=true );
/*
 * DESCRIPTION :
 *   Returns in A a real matrix of dimensions rows and cols,
 *   of power series for exp(x), truncated at degree deg.
 *
 * ON ENTRY :
 *   rows     the number of rows in the matrices x and A;
 *   cols     the number of columns in the matrices x and A;
 *   deg      truncation degree of the series;
 *   x        space for a matrix of doubles, of dimensions rows and cols;
 *   A        space for a matrix of series, of dimensions rows and cols;
 *   expform  if true, then the exp() expansion is applied,
 *            otherwise, the log(1+x) expansion is applied.
 *
 * ON RETURN :
 *   x        a matrix of randomly generated doubles;
 *   A        a matrix of randomly generated power series. */

void random_dbl_upper_matrix ( int rows, int cols, double **A );
/*
 * DESCRIPTION :
 *   Returns in A a real upper triangular matrix of dimensions rows and cols,
 *   of randomly generated doubles.
 *
 * ON ENTRY :
 *   rows     the number of rows in A;
 *   cols     the number of columns in A;
 *   A        space for a matrix of doubles, of dimensions rows and cols.
 *
 * ON RETURN :
 *   A        a matrix of randomly generated doubles. */

void random_dbl_matrix ( int rows, int cols, double **A );
/*
 * DESCRIPTION :
 *   Returns in A a real matrix of dimensions rows and cols,
 *   of randomly generated doubles.
 *
 * ON ENTRY :
 *   rows     the number of rows in A;
 *   cols     the number of columns in A;
 *   A        space for a matrix of doubles, of dimensions rows and cols.
 *
 * ON RETURN :
 *   A        a matrix of randomly generated doubles. */

void random_cmplx_upper_matrix
 ( int rows, int cols, double **Are, double **Aim );
/*
 * DESCRIPTION :
 *   Returns in A a complex upper triangular matrix of dimensions 
 *   rows and cols, of randomly generated doubles.
 *
 * ON ENTRY :
 *   rows     the number of rows in A;
 *   cols     the number of columns in A;
 *   Are      space for the real parts of a matrix of doubles,
 *            of dimensions rows and cols;
 *   Aim      space for the imaginary parts of a matrix of doubles,
 *            of dimensions rows and cols.
 *
 * ON RETURN :
 *   Are      real parts of a matrix of randomly generated doubles;
 *   Aim      imaginary parts of a matrix of randomly generated doubles. */

void random_cmplx_matrix
 ( int rows, int cols, double **Are, double **Aim );
/*
 * DESCRIPTION :
 *   Returns in A a complex matrix of dimensions rows and cols,
 *   of randomly generated doubles.
 *
 * ON ENTRY :
 *   rows     the number of rows in A;
 *   cols     the number of columns in A;
 *   Are      space for the real parts of a matrix of doubles,
 *            of dimensions rows and cols;
 *   Aim      space for the imaginary parts of a matrix of doubles,
 *            of dimensions rows and cols.
 *
 * ON RETURN :
 *   Are      real parts of a matrix of randomly generated doubles;
 *   Aim      imaginary parts of a matrix of randomly generated doubles. */

void random_dbl_upper_series_matrix
 ( int rows, int cols, int deg, double **x, double ***A, bool expform=true );
/*
 * DESCRIPTION :
 *   Returns in A a real upper triangular matrix
 *   of dimensions rows and cols, of power series for exp(x),
 *   truncated at degree deg.
 *
 * ON ENTRY :
 *   rows     the number of rows in the matrices x and A;
 *   cols     the number of columns in the matrices x and A;
 *   deg      truncation degree of the series;
 *   x        space for a matrix of doubles, of dimensions rows and cols;
 *   A        space for a matrix of series, of dimensions rows and cols;
 *   expform  if true, then the exp() expansion is applied,
 *            otherwise, the log(1+x) expansion is applied.
 *
 * ON RETURN :
 *   x        upper triangular matrix of randomly generated doubles,
 *            undefined for lower triangular elements, below the diagonal;
 *   A        upper triangular matrix of randomly generated power series;
 *            elements in the lower triangular part, below the diagonal,
 *            are equal to zero. */

void random_dbl_lower_series_matrix
 ( int rows, int cols, int deg, double **x, double ***A, bool expform=true );
/*
 * DESCRIPTION :
 *   Returns in A a real lower triangular matrix
 *   of dimensions rows and cols, of power series for exp(x),
 *   truncated at degree deg.
 *   The elements on the diagonal are all equal to one.
 *
 * ON ENTRY :
 *   rows     the number of rows in the matrices x and A;
 *   cols     the number of columns in the matrices x and A;
 *   deg      truncation degree of the series;
 *   x        space for a matrix of doubles, of dimensions rows and cols;
 *   A        space for a matrix of series, of dimensions rows and cols;
 *   expform  if true, then the exp() expansion is applied,
 *            otherwise, the log(1+x) expansion is applied.
 *
 * ON RETURN :
 *   x        lower triangular matrix of randomly generated doubles,
 *            undefined for upper triangular elements, above the diagonal,
 *            with zeros on the diagonal;
 *   A        lower triangular matrix of randomly generated power series;
 *            elements in the upper triangular part, above the diagonal,
 *            are equal to zero. */

void random_cmplx_series_matrix
 ( int rows, int cols, int deg, double **xre, double **xim,
   double ***Are, double ***Aim, bool expform=true );
/*
 * DESCRIPTION :
 *   Returns in A a complex matrix of dimensions rows and cols,
 *   of power series for exp(x), truncated at degree deg.
 *
 * ON ENTRY :
 *   rows     the number of rows in the matrices x and A;
 *   cols     the number of columns in the matrices x and A;
 *   deg      truncation degree of the series;
 *   xre      space for a matrix of doubles, of dimensions rows and cols;
 *   xim      space for a matrix of doubles, of dimensions rows and cols;
 *   Are      space for a matrix of series, of dimensions rows and cols;
 *   Aim      space for a matrix of series, of dimensions rows and cols;
 *   expform  if true, then the exp() expansion is applied,
 *            otherwise, the log(1+x) expansion is applied.
 *
 * ON RETURN :
 *   xre      real parts of a matrix of random complex numbers;
 *   xim      imaginary parts of a matrix of random complex numbers;
 *   Are      real parts of a matrix of random power series;
 *   Aim      imaginary parts of a matrix of random power series. */

void random_cmplx_upper_series_matrix
 ( int rows, int cols, int deg, double **xre, double **xim,
   double ***Are, double ***Aim, bool expform=true );
/*
 * DESCRIPTION :
 *   Returns in A a complex upper triangular matrix
 *   of dimensions rows and cols, of power series for exp(x),
 *   for complex x, truncated at degree deg.
 *
 * ON ENTRY :
 *   rows     the number of rows in the matrices x and A;
 *   cols     the number of columns in the matrices x and A;
 *   deg      truncation degree of the series;
 *   xre      space for a matrix of doubles, of dimensions rows and cols;
 *   xim      space for a matrix of doubles, of dimensions rows and cols;
 *   Are      space for a matrix of series, of dimensions rows and cols;
 *   Aim      space for a matrix of series, of dimensions rows and cols;
 *   expform  if true, then the exp() expansion is applied,
 *            otherwise, the log(1+x) expansion is applied.
 *
 * ON RETURN :
 *   xre      real parts of a matrix of random complex numbers,
 *            undefined for lower triangular elements, below the diagonal;
 *   xim      imaginary parts of a matrix of random complex numbers,
 *            undefined for lower triangular elements, below the diagonal;
 *   Are      real parts of an upper triangular matrix of random power series,
 *            elements in the lower triangular part, below the diagonal,
 *            are equal to zero;
 *   Aim      imaginary parts of an upper triangular matrix of random power
 *            series, elements in the lower triangular part, below the
 *            diagonal, are equal to zero. */

void random_cmplx_lower_series_matrix
 ( int rows, int cols, int deg, double **xre, double **xim,
   double ***Are, double ***Aim, bool expform=true );
/*
 * DESCRIPTION :
 *   Returns in A a complex lower triangular matrix
 *   of dimensions rows and cols, of power series for exp(x),
 *   for complex x, truncated at degree deg.
 *   The elements on the diagonal are all equal to one.
 *
 * ON ENTRY :
 *   rows     the number of rows in the matrices x and A;
 *   cols     the number of columns in the matrices x and A;
 *   deg      truncation degree of the series;
 *   xre      space for a matrix of doubles, of dimensions rows and cols;
 *   xim      space for a matrix of doubles, of dimensions rows and cols;
 *   Are      space for a matrix of series, of dimensions rows and cols;
 *   Aim      space for a matrix of series, of dimensions rows and cols;
 *   expform  if true, then the exp() expansion is applied,
 *            otherwise, the log(1+x) expansion is applied.
 *
 * ON RETURN :
 *   xre      real parts of a matrix of random complex numbers,
 *            undefined for upper triangular elements, above the diagonal;
 *   xim      imaginary parts of a matrix of random complex numbers,
 *            undefined for upper triangular elements, above the diagonal;
 *   Are      real parts of an lower triangular matrix of random power series,
 *            elements in the upper triangular part, above the diagonal,
 *            are equal to zero;
 *   Aim      imaginary parts of an lower triangular matrix of random power
 *            series, elements in the upper triangular part, above the
 *            diagonal, are equal to zero. */

#endif
