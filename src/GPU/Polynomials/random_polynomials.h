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

void read_supports ( int dim, int nbr, int *nvr );
/*
 * DESCRIPTION :
 *   Prompts for nbr number of variables in each monomial.
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

bool make_complex_polynomial
 ( int dim, int nbr, int pwr, int deg, int *nvr, int **idx, int **exp,
   double *cstre, double *cstim, double **cffre, double **cffim );
/*
 * DESCRIPTION :
 *   Makes a polynomial in several variables, with power series coefficients,
 *   generating random exponents and complex coefficients.
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
 *   cstre   space allocated for deg+1 doubles;
 *   cstim   space allocated for deg+1 doubles;
 *   cffre   space allocated for nbr arrays of deg+1 doubles.
 *   cffim   space allocated for nbr arrays of deg+1 doubles.
 *
 * ON RETURN :
 *   idx     the participating variables in each monomial,
 *           idx[k] is an array of nvr[k] integers,
 *           idx[k][i] is the index of variable i in monomial k;
 *   exp     the exponents of the variables in each monomial,
 *           exp[k] is an array of nvr[k] integers,
 *           exp[k][i] is the exponent of variable i in monomial k;
 *   cstre   deg+1 coefficients of the real parts of the constant;
 *   cstim   deg+1 coefficients of the imaginary parts of the constant;
 *   cffre   the real parts of the coefficient series for each monomial,
 *           cffre[k] has the deg+1 real parts of monomial k;
 *   cffim   the imaginary parts of the coefficient series for each monomial,
 *           cffim[k] has the deg+1 imaginary parts of monomial k. */

int products_count ( int dim, int nbr );
/*
 * DESCRIPTION :
 *   Returns the count of all choices of nbr distinct numbers
 *   in the range 0 to dim-1. */

void make_product_exponents
 ( int lvl, int dim, int nva, int *accu, int *moncnt, int **idx );
/*
 * DESCRIPTION :
 *   Generates all exponent vectors that define the products.
 *
 * REQUIRED : nva < dim.
 *
 * ON ENTRY :
 *   lvl     current level in the recursion, is zero initially;
 *   dim     dimension, total number of variables;
 *   nva     number of variables in each exponent vector;
 *   accu    the accumulator has space for nva integers;
 *   moncnt  current number of monomials, is zero initially;
 *   idx     has space for all exponent vectors.
 *
 * ON RETURN :
 *   moncnt  equals the number of exponent vectors;
 *   accu    the last values in the accumulator;
 *   idx     the participating variables in each monomial,
 *           idx[k] is an array of nvr[k] integers,
 *           idx[k][i] is the index of variable i in monomial k. */

void insert_sort ( int dim, int *data );
/*
 * DESCRIPTION :
 *   Applies insertion sort to the dim numbers in data. */

void make_cyclic_exponents ( int dim, int nva, int **idx );
/*
 * DESCRIPTION :
 *   Makes all exponents in a cyclic polynomial.
 *
 * REQUIRED : nva < dim.
 *
 * ON ENTRY :
 *   dim     dimension, total number of variables;
 *   nva     number of variables in each monomial;
 *   idx     has space for all nva exponent vectors.
 *
 * ON RETURN :
 *   idx     the participating variables in each monomial,
 *           idx[k] is an array of nvr[k] integers,
 *           idx[k][i] is the index of variable i in monomial k. */

void make_real_products
 ( int dim, int nbr, int nva, int deg, int **idx, double *cst, double **cff );
/*
 * DESCRIPTION :
 *   Returns the sum of all products of size nbr out of dimension dim,
 *   with random real power series coefficients truncated to degree deg.
 *
 * REQUIRED : dim > nbr.
 *
 * ON ENTRY :
 *   dim     dimension, total number of variables;
 *   nbr     number of monomials, excluding the constant;
 *   nva     number of variables in each monomial;
 *   deg     truncation degree of the power series;
 *   idx     space for nbr index vectors;
 *   cst     space for deg + 1 doubles;
 *   cff     space for nbr power series coefficients.
 *
 * ON RETURN :
 *   idx     the participating variables in each monomial,
 *           idx[k] is an array of nva integers,
 *           idx[k][i] is the index of variable i in monomial k;
 *   cst     deg+1 coefficients of the constant monomial;
 *   cff     the coefficient series for each monomial,
 *           cff[k] has the deg+1 coefficients of monomial k. */

void make_complex_products
 ( int dim, int nbr, int nva, int deg, int **idx,
   double *cstre, double *cstim, double **cffre, double **cffim );
/*
 * DESCRIPTION :
 *   Returns the sum of all products of size nbr out of dimension dim,
 *   with random complex power series coefficients truncated to degree deg.
 *
 * REQUIRED : dim > nbr.
 *
 * ON ENTRY :
 *   dim     dimension, total number of variables;
 *   nbr     number of monomials, excluding the constant;
 *   nva     number of variables in each monomial;
 *   deg     truncation degree of the power series;
 *   idx     space for nbr index vectors;
 *   cstre   space for deg + 1 doubles;
 *   cstim   space for deg + 1 doubles;
 *   cffre   space for nbr power series coefficients;
 *   cffim   space for nbr power series coefficients.
 *
 * ON RETURN :
 *   idx     the participating variables in each monomial,
 *           idx[k] is an array of nva integers,
 *           idx[k][i] is the index of variable i in monomial k;
 *   cstre   deg+1 coefficients of the real parts of the constant;
 *   cstim   deg+1 coefficients of the imaginary parts of the constant;
 *   cffre   the real parts of the coefficient series for each monomial,
 *           cffre[k] has the deg+1 real parts of monomial k;
 *   cffim   the imaginary parts of the coefficient series for each monomial,
 *           cffim[k] has the deg+1 imaginary parts of monomial k. */

void make_real_cyclic
 ( int dim, int nva, int deg, int **idx, double *cst, double **cff );
/*
 * DESCRIPTION :
 *   Returns the cyclic polynomial with nva variables in dimension dim,
 *   with random real power series coefficients truncated to degree deg.
 *
 * REQUIRED : dim > nbr.
 *
 * ON ENTRY :
 *   dim     dimension, total number of variables;
 *   nva     number of variables in each monomial;
 *   deg     truncation degree of the power series;
 *   idx     space for nva index vectors;
 *   cst     space for deg + 1 doubles;
 *   cff     space for nva power series coefficients.
 *
 * ON RETURN :
 *   idx     the participating variables in each monomial,
 *           idx[k] is an array of nva integers,
 *           idx[k][i] is the index of variable i in monomial k;
 *   cst     deg+1 coefficients of the constant monomial;
 *   cff     the coefficient series for each monomial,
 *           cff[k] has the deg+1 coefficients of monomial k. */

void make_complex_cyclic
 ( int dim, int nva, int deg, int **idx,
   double *cstre, double *cstim, double **cffre, double **cffim );
/*
 * DESCRIPTION :
 *   Returns the cyclic polynomial with nva variables in dimension dim,
 *   with random complex power series coefficients truncated to degree deg.
 *
 * REQUIRED : dim > nbr.
 *
 * ON ENTRY :
 *   dim     dimension, total number of variables;
 *   nva     number of variables in each monomial;
 *   deg     truncation degree of the power series;
 *   idx     space for nva index vectors;
 *   cstre   space for deg + 1 doubles;
 *   cstim   space for deg + 1 doubles;
 *   cffre   space for nbr power series coefficients;
 *   cffim   space for nbr power series coefficients.
 *
 * ON RETURN :
 *   idx     the participating variables in each monomial,
 *           idx[k] is an array of nva integers,
 *           idx[k][i] is the index of variable i in monomial k;
 *   cstre   deg+1 coefficients of the real parts of the constant;
 *   cstim   deg+1 coefficients of the imaginary parts of the constant;
 *   cffre   the real parts of the coefficient series for each monomial,
 *           cffre[k] has the deg+1 real parts of monomial k;
 *   cffim   the imaginary parts of the coefficient series for each monomial,
 *           cffim[k] has the deg+1 imaginary parts of monomial k. */

#endif
