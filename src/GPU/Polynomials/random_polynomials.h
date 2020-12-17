// The file random_polynomials.h specifies functions to setup tests for the 
// evaluation and differentiation of a polynomial in double precision.

#ifndef __random_polynomials_h__
#define __random_polynomials_h__

void make_supports ( int dim, int nbr, int *nvr );
/*
 * DESCRIPTION :
 *   Generates the number of variables that will appear in each monomial.
 *
 * ON ENTRY :
 *   dim     total number of variables;
 *   nbr     number of monomials;
 *   nvr     space allocated for nbr integers.
 *
 * ON RETURN :
 *   nvr     nvr[k] has the number of variables in monomial k,
 *           nvr[k] ranges between 1 and dim. */

bool make_real_polynomial
 ( int dim, int nbr, int pwr, int deg, int *nvr, int **idx, int **exp,
   double *cst, double **cff );
/*
 * DESCRIPTION :
 *   Makes a polynomial in several variables, with power series coefficients,
 *   generating random exponents and real coefficients.
 *   Writes an error message and returns true if for one k, nvr[k] > dim.
 *
 * ON ENTRY :
 *   dim     dimension, total number of variables;
 *   nbr     number of monomials, excluding the constant term,
 *           nbr is the number of monomials that have at least one variable;
 *   pwr     largest power of a variable;
 *   deg     degree of the power series coefficient;
 *   nvr     nbr integers, with nvr[k] the number of variables in monomial k;
 *   exp     space allocated for nbr arrays of nvr[k] integers;
 *   cst     space allocated for deg+1 doubles;
 *   cff     space allocated for nbr arrays of deg+1 doubles.
 *
 * ON RETURN :
 *   idx     the participating variables in each monomial,
 *           idx[k] is an array of nvr[k] integers,
 *           idx[k][i] is the index of variable i in monomial k;
 *   exp     the exponents of the variables in each monomial,
 *           exp[k] is an array of nvr[k] integers,
 *           exp[k][i] is the exponent of variable i in monomial k;
 *   cst     deg+1 coefficients of the constant monomial;
 *   cff     the coefficient series for each monomial,
 *           cff[k] has the deg+1 coefficients of monomial k. */

#endif
