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

void random_quad_double_complex
 ( double *rehihi, double *relohi, double *rehilo, double *relolo,
   double *imhihi, double *imlohi, double *imhilo, double *imlolo );
/*
 * DESCRIPTION :
 *   Generates the high, middle, and low double of a random cosine c
 *   and computes the corresponding sine as sqrt(1-c*c),
 *   so the randomly generated complex number has modulus one.
 *
 * ON RETURN :
 *   rehihi   the highest double of the real part,
 *   relohi   the second highest double of the real part,
 *   rehilo   the second lowest double of the real part,
 *   relolo   the lowest double of the real part,
 *   imhihi   the highest double of the imaginary part,
 *   imlohi   the second highest double of the imaginary part,
 *   imhilo   the second lowest double of the imaginary part,
 *   imlolo   the lowest double of the imaginary part
 *            of the randomly generated complex number. */

bool make_complex4_polynomial
 ( int dim, int nbr, int pwr, int deg, int *nvr, int **idx, int **exp,
   double *cstrehihi, double *cstrelohi, double *cstrehilo, double *cstrelolo,
   double *cstimhihi, double *cstimlohi, double *cstimhilo, double *cstimlolo,
   double **cffrehihi, double **cffrelohi,
   double **cffrehilo, double **cffrelolo,
   double **cffimhihi, double **cffimlohi,
   double **cffimhilo, double **cffimlolo );
/*
 * DESCRIPTION :
 *   Makes a polynomial in several variables, with power series coefficients,
 *   generating random exponents and complex coefficients.
 *   Writes an error message and returns true if for one k, nvr[k] > dim.
 *
 * ON ENTRY :
 *   dim       dimension, total number of variables;
 *   nbr       number of monomials, excluding the constant term,
 *             nbr is the number of monomials that have at least one variable;
 *   pwr       largest power of a variable;
 *   deg       degree of the power series coefficient;
 *   nvr       nbr integers, with nvr[k] the number of variables in monomial k;
 *   exp       space allocated for nbr arrays of nvr[k] integers;
 *   cstrehihi has space allocated for deg+1 doubles;
 *   cstrelohi has space allocated for deg+1 doubles;
 *   cstrehilo has space allocated for deg+1 doubles;
 *   cstrelolo has space allocated for deg+1 doubles;
 *   cstimhihi has space allocated for deg+1 doubles;
 *   cstimlohi has space allocated for deg+1 doubles;
 *   cstimhilo has space allocated for deg+1 doubles;
 *   cstimlolo has space allocated for deg+1 doubles;
 *   cffrehihi has space allocated for nbr arrays of deg+1 doubles;
 *   cffrelohi has space allocated for nbr arrays of deg+1 doubles;
 *   cffrehilo has space allocated for nbr arrays of deg+1 doubles;
 *   cffrelolo has space allocated for nbr arrays of deg+1 doubles;
 *   cffimhihi has space allocated for nbr arrays of deg+1 doubles;
 *   cffimlohi has space allocated for nbr arrays of deg+1 doubles;
 *   cffimhilo has space allocated for nbr arrays of deg+1 doubles.
 *   cffimlolo has space allocated for nbr arrays of deg+1 doubles.
 *
 * ON RETURN :
 *   idx       the participating variables in each monomial,
 *             idx[k] is an array of nvr[k] integers,
 *             idx[k][i] is the index of variable i in monomial k;
 *   exp       the exponents of the variables in each monomial,
 *             exp[k] is an array of nvr[k] integers,
 *             exp[k][i] is the exponent of variable i in monomial k;
 *   cstrehihi has the deg+1 highest coefficients of the real parts
 *             of the constant monomial;
 *   cstrelohi has the deg+1 second highest coefficients of the real parts
 *             of the constant monomial;
 *   cstrelolo has the deg+1 second lowest coefficients of the real parts
 *             of the constant monomial;
 *   cstrelolo has the deg+1 lowest coefficients of the real parts
 *             of the constant monomial;
 *   cstimhihi has the deg+1 highest coefficients of the imaginary parts
 *             of the constant monomial;
 *   cstimlohi has the deg+1 second highest coefficients of the imaginary parts
 *             of the constant monomial;
 *   cstimhilo has the deg+1 second lowest coefficients of the imaginary parts
 *             of the constant monomial;
 *   cstimlolo has the deg+1 lowest coefficients of the imaginary parts
 *             of the constant monomial;
 *   cffrehihi has the highest doubles of the real parts of the coefficients,
 *             cffrehihi[k] has deg+1 highest coefficients of monomial k;
 *   cffrelohi has the second highest doubles of the real parts
 *             of the coefficients, cffrelohi[k] has deg+1 second highest
 *             coefficients of monomial k;
 *   cffrehilo has the second lowest doubles of the real parts
 *             of the coefficients, cffrehilo[k] has deg+1 second lowest
 *             coefficients of monomial k;
 *   cffrelolo has the lowest doubles of the real parts of the coefficients,
 *             cffrelolo[k] has deg+1 lowest coefficients of monomial k;
 *   cffimhihi has the highest doubles of the imaginary parts
 *             of the coefficients, cffimhihi[k] has deg+1 highest
 *             coefficients of monomial k;
 *   cffimlohi has the second highest doubles of the imaginary parts
 *             of the coefficients, cffimlohi[k] has deg+1 second highest
 *             coefficients of monomial k;
 *   cffimhilo has the second lowest doubles of the imaginary parts
 *             of the coefficient, cffimlolo[k] has the deg+1 second lowest
 *             coefficients of monomial k;
 *   cffimlolo has the lowest doubles of the imaginary parts
 *             of the coefficient, cffimlolo[k] has the deg+1 lowest
 *             coefficients of monomial k. */

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

void make_complex4_products
 ( int dim, int nbr, int nva, int deg, int **idx,
   double *cstrehihi, double *cstrelohi, double *cstrehilo, double *cstrelolo,
   double *cstimhihi, double *cstimlohi, double *cstimhilo, double *cstimlolo,
   double **cffrehihi, double **cffrelohi,
   double **cffrehilo, double **cffrelolo,
   double **cffimhihi, double **cffimlohi,
   double **cffimhilo, double **cffimlolo );
/*
 * DESCRIPTION :
 *   Returns the sum of all products of size nbr out of dimension dim,
 *   with random complex power series coefficients truncated to degree deg.
 *
 * REQUIRED : dim > nbr.
 *
 * ON ENTRY :
 *   dim       dimension, total number of variables;
 *   nbr       number of monomials, excluding the constant term,
 *             nbr is the number of monomials that have at least one variable;
 *   pwr       largest power of a variable;
 *   deg       degree of the power series coefficient;
 *   nvr       nbr integers, with nvr[k] the number of variables in monomial k;
 *   exp       space allocated for nbr arrays of nvr[k] integers;
 *   cstrehihi has space allocated for deg+1 doubles;
 *   cstrelohi has space allocated for deg+1 doubles;
 *   cstrehilo has space allocated for deg+1 doubles;
 *   cstrelolo has space allocated for deg+1 doubles;
 *   cstimhihi has space allocated for deg+1 doubles;
 *   cstimlohi has space allocated for deg+1 doubles;
 *   cstimhilo has space allocated for deg+1 doubles;
 *   cstimlolo has space allocated for deg+1 doubles;
 *   cffrehihi has space allocated for nbr arrays of deg+1 doubles;
 *   cffrelohi has space allocated for nbr arrays of deg+1 doubles;
 *   cffrehilo has space allocated for nbr arrays of deg+1 doubles;
 *   cffrelolo has space allocated for nbr arrays of deg+1 doubles;
 *   cffimhihi has space allocated for nbr arrays of deg+1 doubles;
 *   cffimlohi has space allocated for nbr arrays of deg+1 doubles;
 *   cffimhilo has space allocated for nbr arrays of deg+1 doubles.
 *   cffimlolo has space allocated for nbr arrays of deg+1 doubles.
 *
 * ON RETURN :
 *   idx       the participating variables in each monomial,
 *             idx[k] is an array of nvr[k] integers,
 *             idx[k][i] is the index of variable i in monomial k;
 *   exp       the exponents of the variables in each monomial,
 *             exp[k] is an array of nvr[k] integers,
 *             exp[k][i] is the exponent of variable i in monomial k;
 *   cstrehihi has the deg+1 highest coefficients of the real parts
 *             of the constant monomial;
 *   cstrelohi has the deg+1 second highest coefficients of the real parts
 *             of the constant monomial;
 *   cstrelolo has the deg+1 second lowest coefficients of the real parts
 *             of the constant monomial;
 *   cstrelolo has the deg+1 lowest coefficients of the real parts
 *             of the constant monomial;
 *   cstimhihi has the deg+1 highest coefficients of the imaginary parts
 *             of the constant monomial;
 *   cstimlohi has the deg+1 second highest coefficients of the imaginary parts
 *             of the constant monomial;
 *   cstimhilo has the deg+1 second lowest coefficients of the imaginary parts
 *             of the constant monomial;
 *   cstimlolo has the deg+1 lowest coefficients of the imaginary parts
 *             of the constant monomial;
 *   cffrehihi has the highest doubles of the real parts of the coefficients,
 *             cffrehihi[k] has deg+1 highest coefficients of monomial k;
 *   cffrelohi has the second highest doubles of the real parts
 *             of the coefficients, cffrelohi[k] has deg+1 second highest
 *             coefficients of monomial k;
 *   cffrehilo has the second lowest doubles of the real parts
 *             of the coefficients, cffrehilo[k] has deg+1 second lowest
 *             coefficients of monomial k;
 *   cffrelolo has the lowest doubles of the real parts of the coefficients,
 *             cffrelolo[k] has deg+1 lowest coefficients of monomial k;
 *   cffimhihi has the highest doubles of the imaginary parts
 *             of the coefficients, cffimhihi[k] has deg+1 highest
 *             coefficients of monomial k;
 *   cffimlohi has the second highest doubles of the imaginary parts
 *             of the coefficients, cffimlohi[k] has deg+1 second highest
 *             coefficients of monomial k;
 *   cffimhilo has the second lowest doubles of the imaginary parts
 *             of the coefficient, cffimlolo[k] has the deg+1 second lowest
 *             coefficients of monomial k;
 *   cffimlolo has the lowest doubles of the imaginary parts
 *             of the coefficient, cffimlolo[k] has the deg+1 lowest
 *             coefficients of monomial k. */

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

void make_complex4_cyclic
 ( int dim, int nva, int deg, int **idx,
   double *cstrehihi, double *cstrelohi, double *cstrehilo, double *cstrelolo,
   double *cstimhihi, double *cstimlohi, double *cstimhilo, double *cstimlolo,
   double **cffrehihi, double **cffrelohi,
   double **cffrehilo, double **cffrelolo,
   double **cffimhihi, double **cffimlohi,
   double **cffimhilo, double **cffimlolo );
/*
 * DESCRIPTION :
 *   Returns the cyclic polynomial with nva variables in dimension dim,
 *   with random complex power series coefficients truncated to degree deg.
 *
 * REQUIRED : dim > nbr.
 *
 * ON ENTRY :
 *   dim       dimension, total number of variables;
 *   nbr       number of monomials, excluding the constant term,
 *             nbr is the number of monomials that have at least one variable;
 *   pwr       largest power of a variable;
 *   deg       degree of the power series coefficient;
 *   nvr       nbr integers, with nvr[k] the number of variables in monomial k;
 *   exp       space allocated for nbr arrays of nvr[k] integers;
 *   cstrehihi has space allocated for deg+1 doubles;
 *   cstrelohi has space allocated for deg+1 doubles;
 *   cstrehilo has space allocated for deg+1 doubles;
 *   cstrelolo has space allocated for deg+1 doubles;
 *   cstimhihi has space allocated for deg+1 doubles;
 *   cstimlohi has space allocated for deg+1 doubles;
 *   cstimhilo has space allocated for deg+1 doubles;
 *   cstimlolo has space allocated for deg+1 doubles;
 *   cffrehihi has space allocated for nbr arrays of deg+1 doubles;
 *   cffrelohi has space allocated for nbr arrays of deg+1 doubles;
 *   cffrehilo has space allocated for nbr arrays of deg+1 doubles;
 *   cffrelolo has space allocated for nbr arrays of deg+1 doubles;
 *   cffimhihi has space allocated for nbr arrays of deg+1 doubles;
 *   cffimlohi has space allocated for nbr arrays of deg+1 doubles;
 *   cffimhilo has space allocated for nbr arrays of deg+1 doubles.
 *   cffimlolo has space allocated for nbr arrays of deg+1 doubles.
 *
 * ON RETURN :
 *   idx       the participating variables in each monomial,
 *             idx[k] is an array of nvr[k] integers,
 *             idx[k][i] is the index of variable i in monomial k;
 *   exp       the exponents of the variables in each monomial,
 *             exp[k] is an array of nvr[k] integers,
 *             exp[k][i] is the exponent of variable i in monomial k;
 *   cstrehihi has the deg+1 highest coefficients of the real parts
 *             of the constant monomial;
 *   cstrelohi has the deg+1 second highest coefficients of the real parts
 *             of the constant monomial;
 *   cstrelolo has the deg+1 second lowest coefficients of the real parts
 *             of the constant monomial;
 *   cstrelolo has the deg+1 lowest coefficients of the real parts
 *             of the constant monomial;
 *   cstimhihi has the deg+1 highest coefficients of the imaginary parts
 *             of the constant monomial;
 *   cstimlohi has the deg+1 second highest coefficients of the imaginary parts
 *             of the constant monomial;
 *   cstimhilo has the deg+1 second lowest coefficients of the imaginary parts
 *             of the constant monomial;
 *   cstimlolo has the deg+1 lowest coefficients of the imaginary parts
 *             of the constant monomial;
 *   cffrehihi has the highest doubles of the real parts of the coefficients,
 *             cffrehihi[k] has deg+1 highest coefficients of monomial k;
 *   cffrelohi has the second highest doubles of the real parts
 *             of the coefficients, cffrelohi[k] has deg+1 second highest
 *             coefficients of monomial k;
 *   cffrehilo has the second lowest doubles of the real parts
 *             of the coefficients, cffrehilo[k] has deg+1 second lowest
 *             coefficients of monomial k;
 *   cffrelolo has the lowest doubles of the real parts of the coefficients,
 *             cffrelolo[k] has deg+1 lowest coefficients of monomial k;
 *   cffimhihi has the highest doubles of the imaginary parts
 *             of the coefficients, cffimhihi[k] has deg+1 highest
 *             coefficients of monomial k;
 *   cffimlohi has the second highest doubles of the imaginary parts
 *             of the coefficients, cffimlohi[k] has deg+1 second highest
 *             coefficients of monomial k;
 *   cffimhilo has the second lowest doubles of the imaginary parts
 *             of the coefficient, cffimlolo[k] has the deg+1 second lowest
 *             coefficients of monomial k;
 *   cffimlolo has the lowest doubles of the imaginary parts
 *             of the coefficient, cffimlolo[k] has the deg+1 lowest
 *             coefficients of monomial k. */

#endif
