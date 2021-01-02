// The file random4_polynomials.h specifies functions to setup tests for the 
// evaluation and differentiation of a polynomial in quad double precision.

#ifndef __random4_polynomials_h__
#define __random4_polynomials_h__

bool make_real4_polynomial
 ( int dim, int nbr, int pwr, int deg, int *nvr, int **idx, int **exp,
   double *csthihi, double *cstlohi, double *csthilo, double *cstlolo,
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo );
/*
 * DESCRIPTION :
 *   Makes a polynomial in several variables, with power series coefficients,
 *   generating random exponents and real coefficients.
 *   Writes an error message and returns true if for one k, nvr[k] > dim.
 *
 * ON ENTRY :
 *   dim       dimension, total number of variables;
 *   nbr       number of monomials, excluding the constant term,
 *             nbr is the number of monomials that have at least one variable;
 *   pwr       largest power of a variable;
 *   deg       degree of the power series coefficient;
 *   nvr       nbr integers with counts for the number of variables,
 *             nvr[k] has the number of variables in monomial k;
 *   exp       space allocated for nbr arrays of nvr[k] integers;
 *   csthihi   space allocated for deg+1 doubles;
 *   cstlohi   space allocated for deg+1 doubles;
 *   csthilo   space allocated for deg+1 doubles;
 *   cstlolo   space allocated for deg+1 doubles;
 *   cffhihi   space allocated for nbr arrays of deg+1 doubles;
 *   cfflohi   space allocated for nbr arrays of deg+1 doubles;
 *   cffhilo   space allocated for nbr arrays of deg+1 doubles;
 *   cfflolo   space allocated for nbr arrays of deg+1 doubles.
 *
 * ON RETURN :
 *   idx       the participating variables in each monomial,
 *             idx[k] is an array of nvr[k] integers,
 *             idx[k][i] is the index of variable i in monomial k;
 *   exp       the exponents of the variables in each monomial,
 *             exp[k] is an array of nvr[k] integers,
 *             exp[k][i] is the exponent of variable i in monomial k;
 *   csthihi   deg+1 highest coefficients of the constant monomial;
 *   cstlohi   deg+1 second highest coefficients of the constant monomial;
 *   csthilo   deg+1 second lowest coefficients of the constant monomial;
 *   cstlolo   deg+1 lowest coefficients of the constant monomial;
 *   cffhihi   the highest parts of the coefficient series for each monomial,
 *             cffhihi[k] has the deg+1 highest coefficients of monomial k;
 *   cfflohi   second highest parts of the coefficient series, cfflohi[k]
 *             has the deg+1 second highest coefficients of monomial k;
 *   cffhilo   second lowest parts of the coefficient series, cffhilo[k]
 *             has the deg+1 second lowest coefficients of monomial k;
 *   cfflolo   the lowest parts of the coefficient series for each monomial,
 *             cfflolo[k] has the deg+1 lowest coefficients of monomial k. */

void make_real4_products
 ( int dim, int nbr, int nva, int deg, int **idx,
   double *csthihi, double *cstlohi, double *csthilo, double *cstlolo,
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo );
/*
 * DESCRIPTION :
 *   Returns the sum of all products of size nbr out of dimension dim,
 *   with random power series coefficients truncated to degree deg.
 *
 * REQUIRED : dim > nbr.
 *
 * ON ENTRY :
 *   dim       dimension, total number of variables;
 *   nbr       number of monomials, excluding the constant;
 *   nva       number of variables in each monomial;
 *   deg       truncation degree of the power series;
 *   idx       space for nbr index vectors;
 *   csthihi   space allocated for deg+1 doubles;
 *   cstlohi   space allocated for deg+1 doubles;
 *   csthilo   space allocated for deg+1 doubles;
 *   cstlolo   space allocated for deg+1 doubles;
 *   cffhihi   space allocated for nbr arrays of deg+1 doubles;
 *   cfflohi   space allocated for nbr arrays of deg+1 doubles;
 *   cffhilo   space allocated for nbr arrays of deg+1 doubles;
 *   cfflolo   space allocated for nbr arrays of deg+1 doubles.
 *
 * ON RETURN :
 *   idx       the participating variables in each monomial,
 *             idx[k] is an array of nva integers,
 *             idx[k][i] is the index of variable i in monomial k;
 *   csthihi   deg+1 highest coefficients of the constant monomial;
 *   cstlohi   deg+1 second highest coefficients of the constant monomial;
 *   csthilo   deg+1 second lowest coefficients of the constant monomial;
 *   cstlolo   deg+1 lowest coefficients of the constant monomial;
 *   cffhihi   the highest parts of the coefficient series for each monomial,
 *             cffhihi[k] has the deg+1 highest coefficients of monomial k;
 *   cfflohi   second highest parts of the coefficient series, cfflohi[k]
 *             has the deg+1 second highest coefficients of monomial k;
 *   cffhilo   second lowest parts of the coefficient series, cffhilo[k]
 *             has the deg+1 second lowest coefficients of monomial k;
 *   cfflolo   the lowest parts of the coefficient series for each monomial,
 *             cfflolo[k] has the deg+1 lowest coefficients of monomial k. */

void make_real4_cyclic
 ( int dim, int nva, int deg, int **idx,
   double *csthihi, double *cstlohi, double *csthilo, double *cstlolo,
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo );
/*
 * DESCRIPTION :
 *   Returns the cyclic polynomial with nva variables in dimension dim,
 *   with random power series coefficients truncated to degree deg.
 *
 * REQUIRED : dim > nbr.
 *
 * ON ENTRY :
 *   dim       dimension, total number of variables;
 *   nva       number of variables in each monomial;
 *   deg       truncation degree of the power series;
 *   idx       space for nva index vectors;
 *   csthihi   space allocated for deg+1 doubles;
 *   cstlohi   space allocated for deg+1 doubles;
 *   csthilo   space allocated for deg+1 doubles;
 *   cstlolo   space allocated for deg+1 doubles;
 *   cffhihi   space allocated for nbr arrays of deg+1 doubles;
 *   cfflohi   space allocated for nbr arrays of deg+1 doubles;
 *   cffhilo   space allocated for nbr arrays of deg+1 doubles;
 *   cfflolo   space allocated for nbr arrays of deg+1 doubles.
 *
 * ON RETURN :
 *   idx       the participating variables in each monomial,
 *             idx[k] is an array of nva integers,
 *             idx[k][i] is the index of variable i in monomial k;
 *   csthihi   deg+1 highest coefficients of the constant monomial;
 *   cstlohi   deg+1 second highest coefficients of the constant monomial;
 *   csthilo   deg+1 second lowest coefficients of the constant monomial;
 *   cstlolo   deg+1 lowest coefficients of the constant monomial;
 *   cffhihi   the highest parts of the coefficient series for each monomial,
 *             cffhihi[k] has the deg+1 highest coefficients of monomial k;
 *   cfflohi   second highest parts of the coefficient series, cfflohi[k]
 *             has the deg+1 second highest coefficients of monomial k;
 *   cffhilo   second lowest parts of the coefficient series, cffhilo[k]
 *             has the deg+1 second lowest coefficients of monomial k;
 *   cfflolo   the lowest parts of the coefficient series for each monomial,
 *             cfflolo[k] has the deg+1 lowest coefficients of monomial k. */

#endif
