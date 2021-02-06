// The file random_matrices.h specifies functions define tests to work
// with matrices and vectors of series in double precision.

#ifndef __random_matrices_h__
#define __random_matrices_h__

void random_dbl_series_vectors
 ( int dim, int deg, double *x, double **plux, double **minx );
/*
 * DESCRIPTION :
 *   Returns a pair of two vectors with series for exp(x[i]) and exp(-x[i])
 *   for i ranging from 0 to dim-1.
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

#endif
