// The file dbl2_monomial_systems.h specifies functions to generate
// monomial systems with prescribed solutions in double double precision.

#ifndef __dbl2_monomial_systems_h__
#define __dbl2_monomial_systems_h__

void make_complex2_exponentials
 ( int dim, int deg,
   double **srehi, double **srelo, double **simhi, double **simlo );
/*
 * DESCRIPTION :
 *   Returns dim expansions of exp(x) truncated at degree deg,
 *   for random complex values of x on the complex unit circle.
 *
 * ON ENTRY :
 *   dim      number of power series;
 *   deg      truncation degree;
 *   srehi    space for dim arrays of size deg+1;
 *   srelo    space for dim arrays of size deg+1;
 *   simhi    space for dim arrays of size deg+1;
 *   simlo    space for dim arrays of size deg+1.
 *
 * ON RETURN :
 *   srehi    high doubles of the real parts of the dim power series;
 *   srelo    low doubles of the real parts of the dim power series;
 *   simhi    high doubles of the imaginary parts of the dim power series;
 *   simlo    low doubles of the imaginary parts of the dim power series. */

void evaluate_complex2_monomials
 ( int dim, int deg, int **rowsA,
   double **srehi, double **srelo, double **simhi, double **simlo,
   double **rhsrehi, double **rhsrelo, double **rhsimhi, double **rhsimlo );
/*
 * DESCRIPTION :
 *   Evaluates the monomials defined in the rows of a matrix at
 *   complex series to make the right hand side of a monomial system.
 *
 * ON ENTRY :
 *   dim      dimension of the monomial system;
 *   deg      truncation degree of the series;
 *   rowsA    the rows of A have the exponents of the monomials;
 *   srehi    high doubles of the real parts of the series;
 *   srelo    low doubles of the real parts of the series;
 *   simhi    high doubles of the imaginary parts of the series;
 *   simlo    low doubles of the imaginary parts of the series;
 *   rhsrehi  space for dim arrays of size deg+1;
 *   rhsrelo  space for dim arrays of size deg+1;
 *   rhsimhi  space for dim arrays of size deg+1;
 *   rhsimlo  space for dim arrays of size deg+1.
 *
 * ON RETURN :
 *   rhsrehi  high doubles of the real parts of the evaluated monomials;
 *   rhsrelo  low doubles of the real parts of the evaluated monomials;
 *   rhsimhi  high doubles of the imaginary parts of the evaluated monomials;
 *   rhsimlo  low doubles of the imaginary parts of the evaluated monomials. */

#endif
