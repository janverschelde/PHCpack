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

int duplicate_index ( int dim, int nbr, int *nvr, int **idx, int monidx );
/*
 * DESCRIPTION :
 *   Returns the index < monidx for which idx[index][i] == idx[monidx][i]
 *   for all i in 0..nvr[idx]-1, where nvr[idx] == nvr[monidx],
 *   Returns -1 if no such index exists and there are no duplicate
 *   support vectors.
 *
 * REQUIRED : monidx < nbr.
 *
 * ON ENTRY :
 *   dim     dimension, total number of variables;
 *   nbr     number of monomials, excluding the constant term,
 *           nbr is the number of monomials that have at least one variable;
 *   nvr     nbr integers, with nvr[k] the number of variables in monomial k;
 *   idx     the participating variables in each monomial,
 *           idx[k] is an array of nvr[k] integers,
 *           idx[k][i] is the index of variable i in monomial k;
 *   monidx  index of a monomial, monidx < nbr.
 *
 * ON RETURN :
 *   -1      if no duplicate support vectors, or the index
 *           of the first duplicate support vector. */

bool duplicate_supports
 ( int dim, int nbr, int *nvr, int **idx, bool verbose );
/*
 * DESCRIPTION :
 *   Returns true if the polynomial with supports in idx has duplicate
 *   support vectors, returns false otherwise.
 *
 * ON ENTRY :
 *   dim     dimension, total number of variables;
 *   nbr     number of monomials, excluding the constant term,
 *           nbr is the number of monomials that have at least one variable;
 *   nvr     nbr integers, with nvr[k] the number of variables in monomial k;
 *   idx     the participating variables in each monomial,
 *           idx[k] is an array of nvr[k] integers,
 *           idx[k][i] is the index of variable i in monomial k;
 *   verbose indicates if all indices with duplicates are printed,
 *           if false, then there will be no output. */

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
