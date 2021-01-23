// The file random3_polynomials.h specifies functions to setup tests for the 
// evaluation and differentiation of a polynomial in triple double precision.

#ifndef __random3_polynomials_h__
#define __random3_polynomials_h__

bool make_real3_polynomial
 ( int dim, int nbr, int pwr, int deg, int *nvr, int **idx, int **exp,
   double *csthi, double *cstmi, double *cstlo,
   double **cffhi, double **cffmi, double **cfflo );
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
 *   cstmi   space allocated for deg+1 doubles;
 *   cstlo   space allocated for deg+1 doubles;
 *   cffhi   space allocated for nbr arrays of deg+1 doubles;
 *   cffmi   space allocated for nbr arrays of deg+1 doubles;
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
 *   cstmi   deg+1 middle coefficients of the constant monomial;
 *   cstlo   deg+1 low coefficients of the constant monomial;
 *   cffhi   the high parts of the coefficient series for each monomial,
 *           cffhi[k] has the deg+1 high coefficients of monomial k;
 *   cffmi   the middle parts of the coefficient series for each monomial,
 *           cffhi[k] has the deg+1 middle coefficients of monomial k;
 *   cfflo   the low parts of the coefficient series for each monomial,
 *           cfflo[k] has the deg+1 low coefficients of monomial k. */

void random_triple_double_complex
 ( double *rehi, double *remi, double *relo,
   double *imhi, double *immi, double *imlo );
/*
 * DESCRIPTION :
 *   Generates the high, middle, and low double of a random cosine c
 *   and computes the corresponding sine as sqrt(1-c*c),
 *   so the randomly generated complex number has modulus one.
 *
 * ON RETURN :
 *   rehi     the high double of the real part,
 *   remi     the middle double of the real part,
 *   relo     the low double of the real part,
 *   imhi     the high double of the imaginary part,
 *   immi     the middle double of the imaginary part,
 *   imlo     the low double of the imaginary part
 *            of the randomly generated complex number. */

bool make_complex3_polynomial
 ( int dim, int nbr, int pwr, int deg, int *nvr, int **idx, int **exp,
   double *cstrehi, double *cstremi, double *cstrelo,
   double *cstimhi, double *cstimmi, double *cstimlo,
   double **cffrehi, double **cffremi, double **cffrelo,
   double **cffimhi, double **cffimmi, double **cffimlo );
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
 *   cstremi has space allocated for deg+1 doubles;
 *   cstrelo has space allocated for deg+1 doubles;
 *   cstimhi has space allocated for deg+1 doubles;
 *   cstimmi has space allocated for deg+1 doubles;
 *   cstimlo has space allocated for deg+1 doubles;
 *   cffrehi has space allocated for nbr arrays of deg+1 doubles;
 *   cffremi has space allocated for nbr arrays of deg+1 doubles;
 *   cffrelo has space allocated for nbr arrays of deg+1 doubles;
 *   cffimhi has space allocated for nbr arrays of deg+1 doubles;
 *   cffimmi has space allocated for nbr arrays of deg+1 doubles;
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
 *   cstremi has the deg+1 middle coefficients of the real parts
 *           of the constant monomial;
 *   cstrelo has the deg+1 low coefficients of the real parts
 *           of the constant monomial;
 *   cstimhi has the deg+1 high coefficients of the imaginary parts
 *           of the constant monomial;
 *   cstimmi has the deg+1 middle coefficients of the imaginary parts
 *           of the constant monomial;
 *   cstimlo has the deg+1 low coefficients of the imaginary parts
 *           of the constant monomial;
 *   cffrehi has the high doubles of the real parts of the coefficients,
 *           cffrehi[k] has deg+1 high coefficients of monomial k;
 *   cffremi has the middle doubles of the real parts of the coefficients,
 *           cffremi[k] has deg+1 high coefficients of monomial k;
 *   cffrelo has the low doubles of the real parts of the coefficients,
 *           cffrelo[k] has deg+1 low coefficients of monomial k;
 *   cffimhi has the high doubles of the imaginary parts of the coefficients,
 *           cffimhi[k] has deg+1 high coefficients of monomial k;
 *   cffimmi has the middle doubles of the imaginary parts of the coefficients,
 *           cffimmi[k] has deg+1 high coefficients of monomial k;
 *   cffimlo has the low doubles of the imaginary parts of the coefficient,
 *           cffimlo[k] has the deg+1 low coefficients of monomial k. */

void make_real3_products
 ( int dim, int nbr, int nva, int deg, int **idx,
   double *csthi, double *cstmi, double *cstlo,
   double **cffhi, double **cffmi, double **cfflo );
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
 *   cstmi   space allocated for deg+1 doubles;
 *   cstlo   space allocated for deg+1 doubles;
 *   cffhi   space allocated for nbr arrays of deg+1 doubles;
 *   cffmi   space allocated for nbr arrays of deg+1 doubles;
 *   cfflo   space allocated for nbr arrays of deg+1 doubles.
 *
 * ON RETURN :
 *   idx     the participating variables in each monomial,
 *           idx[k] is an array of nva integers,
 *           idx[k][i] is the index of variable i in monomial k;
 *   csthi   deg+1 high coefficients of the constant monomial;
 *   cstmi   deg+1 middle coefficients of the constant monomial;
 *   cstlo   deg+1 low coefficients of the constant monomial;
 *   cffhi   the high parts of the coefficient series for each monomial,
 *           cffhi[k] has the deg+1 high coefficients of monomial k;
 *   cffmi   the middle parts of the coefficient series for each monomial,
 *           cffmi[k] has the deg+1 middle coefficients of monomial k;
 *   cfflo   the low parts of the coefficient series for each monomial,
 *           cfflo[k] has the deg+1 low coefficients of monomial k. */

void make_complex3_products
 ( int dim, int nbr, int nva, int deg, int **idx,
   double *cstrehi, double *cstremi, double *cstrelo,
   double *cstimhi, double *cstimmi, double *cstimlo,
   double **cffrehi, double **cffremi, double **cffrelo,
   double **cffimhi, double **cffimmi, double **cffimlo );
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
 *   cstremi has space allocated for deg+1 doubles;
 *   cstrelo has space allocated for deg+1 doubles;
 *   cstimhi has space allocated for deg+1 doubles;
 *   cstimmi has space allocated for deg+1 doubles;
 *   cstimlo has space allocated for deg+1 doubles;
 *   cffrehi has space allocated for nbr arrays of deg+1 doubles;
 *   cffremi has space allocated for nbr arrays of deg+1 doubles;
 *   cffrelo has space allocated for nbr arrays of deg+1 doubles;
 *   cffimhi has space allocated for nbr arrays of deg+1 doubles;
 *   cffimmi has space allocated for nbr arrays of deg+1 doubles;
 *   cffimlo has space allocated for nbr arrays of deg+1 doubles.
 *
 * ON RETURN :
 *   idx     the participating variables in each monomial,
 *           idx[k] is an array of nva integers,
 *           idx[k][i] is the index of variable i in monomial k;
 *   cstrehi has the deg+1 high coefficients of the real parts
 *           of the constant monomial;
 *   cstremi has the deg+1 middle coefficients of the real parts
 *           of the constant monomial;
 *   cstrelo has the deg+1 low coefficients of the real parts
 *           of the constant monomial;
 *   cstimhi has the deg+1 high coefficients of the imaginary parts
 *           of the constant monomial;
 *   cstimmi has the deg+1 middle coefficients of the imaginary parts
 *           of the constant monomial;
 *   cstimlo has the deg+1 low coefficients of the imaginary parts
 *           of the constant monomial;
 *   cffrehi has the high doubles of the real parts of the coefficients,
 *           cffrehi[k] has deg+1 high coefficients of monomial k;
 *   cffremi has the middle doubles of the real parts of the coefficients,
 *           cffremi[k] has deg+1 high coefficients of monomial k;
 *   cffrelo has the low doubles of the real parts of the coefficients,
 *           cffrelo[k] has deg+1 low coefficients of monomial k;
 *   cffimhi has the high doubles of the imaginary parts of the coefficients,
 *           cffimhi[k] has deg+1 high coefficients of monomial k;
 *   cffimmi has the middle doubles of the imaginary parts of the coefficients,
 *           cffimmi[k] has deg+1 high coefficients of monomial k;
 *   cffimlo has the low doubles of the imaginary parts of the coefficient,
 *           cffimlo[k] has the deg+1 low coefficients of monomial k. */

void make_real3_cyclic
 ( int dim, int nva, int deg, int **idx,
   double *csthi, double *cstmi, double *cstlo,
   double **cffhi, double **cffmi, double **cfflo );
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
 *   cstmi   space allocated for deg+1 doubles;
 *   cstlo   space allocated for deg+1 doubles;
 *   cffhi   space allocated for nbr arrays of deg+1 doubles;
 *   cffmi   space allocated for nbr arrays of deg+1 doubles;
 *   cfflo   space allocated for nbr arrays of deg+1 doubles.
 *
 * ON RETURN :
 *   idx     the participating variables in each monomial,
 *           idx[k] is an array of nva integers,
 *           idx[k][i] is the index of variable i in monomial k;
 *   csthi   deg+1 high coefficients of the constant monomial;
 *   cstmi   deg+1 middle coefficients of the constant monomial;
 *   cstlo   deg+1 low coefficients of the constant monomial;
 *   cffhi   the high parts of the coefficient series for each monomial,
 *           cffhi[k] has the deg+1 high coefficients of monomial k;
 *   cffmi   the middle parts of the coefficient series for each monomial,
 *           cffmi[k] has the deg+1 middle coefficients of monomial k;
 *   cfflo   the low parts of the coefficient series for each monomial,
 *           cfflo[k] has the deg+1 low coefficients of monomial k. */

void make_complex3_cyclic
 ( int dim, int nva, int deg, int **idx,
   double *cstrehi, double *cstremi, double *cstrelo,
   double *cstimhi, double *cstimmi, double *cstimlo,
   double **cffrehi, double **cffremi, double **cffrelo,
   double **cffimhi, double **cffimmi, double **cffimlo );
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
 *   cstremi has space allocated for deg+1 doubles;
 *   cstrelo has space allocated for deg+1 doubles;
 *   cstimhi has space allocated for deg+1 doubles;
 *   cstimmi has space allocated for deg+1 doubles;
 *   cstimlo has space allocated for deg+1 doubles;
 *   cffrehi has space allocated for nbr arrays of deg+1 doubles;
 *   cffremi has space allocated for nbr arrays of deg+1 doubles;
 *   cffrelo has space allocated for nbr arrays of deg+1 doubles;
 *   cffimhi has space allocated for nbr arrays of deg+1 doubles;
 *   cffimmi has space allocated for nbr arrays of deg+1 doubles;
 *   cffimlo has space allocated for nbr arrays of deg+1 doubles.
 *
 * ON RETURN :
 *   idx     the participating variables in each monomial,
 *           idx[k] is an array of nva integers,
 *           idx[k][i] is the index of variable i in monomial k;
 *   cstrehi has the deg+1 high coefficients of the real parts
 *           of the constant monomial;
 *   cstremi has the deg+1 middle coefficients of the real parts
 *           of the constant monomial;
 *   cstrelo has the deg+1 low coefficients of the real parts
 *           of the constant monomial;
 *   cstimhi has the deg+1 high coefficients of the imaginary parts
 *           of the constant monomial;
 *   cstimmi has the deg+1 middle coefficients of the imaginary parts
 *           of the constant monomial;
 *   cstimlo has the deg+1 low coefficients of the imaginary parts
 *           of the constant monomial;
 *   cffrehi has the high doubles of the real parts of the coefficients,
 *           cffrehi[k] has deg+1 high coefficients of monomial k;
 *   cffremi has the middle doubles of the real parts of the coefficients,
 *           cffremi[k] has deg+1 high coefficients of monomial k;
 *   cffrelo has the low doubles of the real parts of the coefficients,
 *           cffrelo[k] has deg+1 low coefficients of monomial k;
 *   cffimhi has the high doubles of the imaginary parts of the coefficients,
 *           cffimhi[k] has deg+1 high coefficients of monomial k;
 *   cffimmi has the middle doubles of the imaginary parts of the coefficients,
 *           cffimmi[k] has deg+1 high coefficients of monomial k;
 *   cffimlo has the low doubles of the imaginary parts of the coefficient,
 *           cffimlo[k] has the deg+1 low coefficients of monomial k. */

#endif
