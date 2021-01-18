// The file random2_polynomials.h specifies functions to setup tests for the 
// evaluation and differentiation of a polynomial in double double precision.

#ifndef __random2_polynomials_h__
#define __random2_polynomials_h__

bool make_real2_polynomial
 ( int dim, int nbr, int pwr, int deg, int *nvr, int **idx, int **exp,
   double *csthi, double *cstlo, double **cffhi, double **cfflo );
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
 *   csthi   space allocated for deg+1 doubles;
 *   cstlo   space allocated for deg+1 doubles;
 *   cffhi   space allocated for nbr arrays of deg+1 doubles;
 *   cfflo   space allocated for nbr arrays of deg+1 doubles.
 *
 * ON RETURN :
 *   idx     the participating variables in each monomial,
 *           idx[k] is an array of nvr[k] integers,
 *           idx[k][i] is the index of variable i in monomial k;
 *   exp     the exponents of the variables in each monomial,
 *           exp[k] is an array of nvr[k] integers,
 *           exp[k][i] is the exponent of variable i in monomial k;
 *   csthi   deg+1 high coefficients of the constant monomial;
 *   cstlo   deg+1 low coefficients of the constant monomial;
 *   cffhi   the high parts of the coefficient series for each monomial,
 *           cffhi[k] has the deg+1 high coefficients of monomial k;
 *   cfflo   the low parts of the coefficient series for each monomial,
 *           cfflo[k] has the deg+1 low coefficients of monomial k. */

void random_double_double_complex
 ( double *rehi, double *relo, double *imhi, double *imlo );
/*
 * DESCRIPTION :
 *   Generates the high and low double of a random cosine c
 *   and computes the corresponding sine as sqrt(1-c*c),
 *   so the randomly generated complex number has modulus one.
 *
 * ON RETURN :
 *   rehi     the high double of the real part,
 *   relo     the low double of the real part,
 *   imhi     the high double of the imaginary part,
 *   imlo     the low double of the imaginary part
 *            of the randomly generated complex number. */

bool make_complex2_polynomial
 ( int dim, int nbr, int pwr, int deg, int *nvr, int **idx, int **exp,
   double *cstrehi, double *cstrelo, double *cstimhi, double *cstimlo,
   double **cffrehi, double **cffrelo, double **cffimhi, double **cffimlo );
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
 *   cstrehi has space allocated for deg+1 doubles;
 *   cstrelo has space allocated for deg+1 doubles;
 *   cstimhi has space allocated for deg+1 doubles;
 *   cstimlo has space allocated for deg+1 doubles;
 *   cffrehi has space allocated for nbr arrays of deg+1 doubles;
 *   cffrelo has space allocated for nbr arrays of deg+1 doubles;
 *   cffimhi has space allocated for nbr arrays of deg+1 doubles;
 *   cffimlo has space allocated for nbr arrays of deg+1 doubles.
 *
 * ON RETURN :
 *   idx     the participating variables in each monomial,
 *           idx[k] is an array of nvr[k] integers,
 *           idx[k][i] is the index of variable i in monomial k;
 *   exp     the exponents of the variables in each monomial,
 *           exp[k] is an array of nvr[k] integers,
 *           exp[k][i] is the exponent of variable i in monomial k;
 *   cstrehi has the deg+1 high coefficients of the real parts
 *           of the constant monomial;
 *   cstrelo has the deg+1 low coefficients of the real parts
 *           of the constant monomial;
 *   cstimhi has the deg+1 high coefficients of the imaginary parts
 *           of the constant monomial;
 *   cstimlo has the deg+1 low coefficients of the imaginary parts
 *           of the constant monomial;
 *   cffrehi has the high doubles of the real parts of the coefficients,
 *           cffrehi[k] has deg+1 high coefficients of monomial k;
 *   cffrelo has the low doubles of the real parts of the coefficients,
 *           cffrelo[k] has deg+1 low coefficients of monomial k;
 *   cffimhi has the high doubles of the imaginary parts of the coefficients,
 *           cffimhi[k] has deg+1 high coefficients of monomial k;
 *   cffimlo has the low doubles of the imaginary parts of the coefficient,
 *           cffimlo[k] has the deg+1 low coefficients of monomial k. */

void make_real2_products
 ( int dim, int nbr, int nva, int deg, int **idx,
   double *csthi, double *cstlo, double **cffhi, double **cfflo );
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
 *   csthi   space allocated for deg+1 doubles;
 *   cstlo   space allocated for deg+1 doubles;
 *   cffhi   space allocated for nbr arrays of deg+1 doubles;
 *   cfflo   space allocated for nbr arrays of deg+1 doubles.
 *
 * ON RETURN :
 *   idx     the participating variables in each monomial,
 *           idx[k] is an array of nva integers,
 *           idx[k][i] is the index of variable i in monomial k;
 *   csthi   deg+1 high coefficients of the constant monomial;
 *   cstlo   deg+1 low coefficients of the constant monomial;
 *   cffhi   the high parts of the coefficient series for each monomial,
 *           cffhi[k] has the deg+1 high coefficients of monomial k;
 *   cfflo   the low parts of the coefficient series for each monomial,
 *           cfflo[k] has the deg+1 high coefficients of monomial k. */

void make_complex2_products
 ( int dim, int nbr, int nva, int deg, int **idx,
   double *cstrehi, double *cstrelo, double *cstimhi, double *cstimlo,
   double **cffrehi, double **cffrelo, double **cffimhi, double **cffimlo );
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
 *   cstrehi has space allocated for deg+1 doubles;
 *   cstrelo has space allocated for deg+1 doubles;
 *   cstimhi has space allocated for deg+1 doubles;
 *   cstimlo has space allocated for deg+1 doubles;
 *   cffrehi has space allocated for nbr arrays of deg+1 doubles;
 *   cffrelo has space allocated for nbr arrays of deg+1 doubles;
 *   cffimhi has space allocated for nbr arrays of deg+1 doubles;
 *   cffimlo has space allocated for nbr arrays of deg+1 doubles.
 *
 * ON RETURN :
 *   idx     the participating variables in each monomial,
 *           idx[k] is an array of nva integers,
 *           idx[k][i] is the index of variable i in monomial k;
 *   cstrehi has the deg+1 high coefficients of the real parts
 *           of the constant monomial;
 *   cstrelo has the deg+1 low coefficients of the real parts
 *           of the constant monomial;
 *   cstimhi has the deg+1 high coefficients of the imaginary parts
 *           of the constant monomial;
 *   cstimlo has the deg+1 low coefficients of the imaginary parts
 *           of the constant monomial;
 *   cffrehi has the high doubles of the real parts of the coefficients,
 *           cffrehi[k] has deg+1 high coefficients of monomial k;
 *   cffrelo has the low doubles of the real parts of the coefficients,
 *           cffrelo[k] has deg+1 low coefficients of monomial k;
 *   cffimhi has the high doubles of the imaginary parts of the coefficients,
 *           cffimhi[k] has deg+1 high coefficients of monomial k;
 *   cffimlo has the low doubles of the imaginary parts of the coefficient,
 *           cffimlo[k] has the deg+1 low coefficients of monomial k. */

void make_real2_cyclic
 ( int dim, int nva, int deg, int **idx,
   double *csthi, double *cstlo, double **cffhi, double **cfflo );
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
 *   csthi   space allocated for deg+1 doubles;
 *   cstlo   space allocated for deg+1 doubles;
 *   cffhi   space allocated for nbr arrays of deg+1 doubles;
 *   cfflo   space allocated for nbr arrays of deg+1 doubles.
 *
 * ON RETURN :
 *   idx     the participating variables in each monomial,
 *           idx[k] is an array of nva integers,
 *           idx[k][i] is the index of variable i in monomial k;
 *   csthi   deg+1 high coefficients of the constant monomial;
 *   cstlo   deg+1 low coefficients of the constant monomial;
 *   cffhi   the high parts of the coefficient series for each monomial,
 *           cffhi[k] has the deg+1 high coefficients of monomial k;
 *   cfflo   the low parts of the coefficient series for each monomial,
 *           cfflo[k] has the deg+1 high coefficients of monomial k. */

void make_complex2_cyclic
 ( int dim, int nva, int deg, int **idx,
   double *cstrehi, double *cstrelo, double *cstimhi, double *cstimlo,
   double **cffrehi, double **cffrelo, double **cffimhi, double **cffimlo );
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
 *   cstrehi has space allocated for deg+1 doubles;
 *   cstrelo has space allocated for deg+1 doubles;
 *   cstimhi has space allocated for deg+1 doubles;
 *   cstimlo has space allocated for deg+1 doubles;
 *   cffrehi has space allocated for nbr arrays of deg+1 doubles;
 *   cffrelo has space allocated for nbr arrays of deg+1 doubles;
 *   cffimhi has space allocated for nbr arrays of deg+1 doubles;
 *   cffimlo has space allocated for nbr arrays of deg+1 doubles.
 *
 * ON RETURN :
 *   idx     the participating variables in each monomial,
 *           idx[k] is an array of nva integers,
 *           idx[k][i] is the index of variable i in monomial k;
 *   cstrehi has the deg+1 high coefficients of the real parts
 *           of the constant monomial;
 *   cstrelo has the deg+1 low coefficients of the real parts
 *           of the constant monomial;
 *   cstimhi has the deg+1 high coefficients of the imaginary parts
 *           of the constant monomial;
 *   cstimlo has the deg+1 low coefficients of the imaginary parts
 *           of the constant monomial;
 *   cffrehi has the high doubles of the real parts of the coefficients,
 *           cffrehi[k] has deg+1 high coefficients of monomial k;
 *   cffrelo has the low doubles of the real parts of the coefficients,
 *           cffrelo[k] has deg+1 low coefficients of monomial k;
 *   cffimhi has the high doubles of the imaginary parts of the coefficients,
 *           cffimhi[k] has deg+1 high coefficients of monomial k;
 *   cffimlo has the low doubles of the imaginary parts of the coefficient,
 *           cffimlo[k] has the deg+1 low coefficients of monomial k. */

#endif
