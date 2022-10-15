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

#endif
