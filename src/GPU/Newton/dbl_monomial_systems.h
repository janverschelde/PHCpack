// The file dbl_monomial_systems.h specifies functions to generate
// monomial systems with prescribed solutions in double precision.

#ifndef __dbl_monomial_systems_h__
#define __dbl_monomial_systems_h__

void make_complex_exponentials
 ( int dim, int deg, double *angles, double **sre, double **sim );
/*
 * DESCRIPTION :
 *   Returns dim expansions of exp(x) truncated at degree deg,
 *   for random complex values of x on the complex unit circle.
 *
 * ON ENTRY :
 *   dim      number of power series;
 *   deg      truncation degree;
 *   angles   has space for dim angles;
 *   sre      space for dim arrays of size deg+1;
 *   sim      space for dim arrays of size deg+1.
 *
 * ON RETURN :
 *   angles   contains the angles of the random complex numbers,
 *            cos(angles[i]) is sre[i][1] and sin(angles[i]) is sim[i][1];
 *   sre      real parts of the coefficients of the dim power series;
 *   sim      imaginary parts of the coefficients of the dim power series. */

void evaluate_complex_monomials
 ( int dim, int deg, int **rowsA,
   double **sre, double **sim, double **rhsre, double **rhsim );
/*
 * DESCRIPTION :
 *   Evaluates the monomials defined in the rows of a matrix at
 *   complex series to make the right hand side of a monomial system.
 *
 * ON ENTRY :
 *   dim      dimension of the monomial system;
 *   deg      truncation degree of the series;
 *   rowsA    the rows of A have the exponents of the monomials;
 *   sre      real parts of the series;
 *   sim      imaginary parts of the series;
 *   rhsre    space for dim arrays of size deg+1;
 *   rhsim    space for dim arrays of size deg+1.
 *
 * ON RETURN :
 *   rhsre    real parts of the evaluated monomials;
 *   rhsim    imaginary parts of the evaluated monomials. */

#endif
