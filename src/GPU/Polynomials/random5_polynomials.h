// The file random5_polynomials.h specifies functions to setup tests for the 
// evaluation and differentiation of a polynomial in penta double precision.

#ifndef __random5_polynomials_h__
#define __random5_polynomials_h__

bool make_real5_polynomial
 ( int dim, int nbr, int pwr, int deg, int *nvr, int **idx, int **exp,
   double *csttb, double *cstix, double *cstmi, double *cstrg, double *cstpk,
   double **cfftb, double **cffix, double **cffmi, double **cffrg,
   double **cffpk );
/*
 * DESCRIPTION :
 *   Makes a polynomial in several variables, with power series coefficients,
 *   generating random exponents and real coefficients.
 *   Writes an error message and returns true if for one k, nvr[k] > dim.
 *
 * ON ENTRY :
 *   dim      dimension, total number of variables;
 *   nbr      number of monomials, excluding the constant term,
 *            nbr is the number of monomials that have at least one variable;
 *   pwr      largest power of a variable;
 *   deg      degree of the power series coefficient;
 *   nvr      nbr integers with counts for the number of variables,
 *            nvr[k] has the number of variables in monomial k;
 *   exp      space allocated for nbr arrays of nvr[k] integers;
 *   csttb    space allocated for deg+1 doubles;
 *   cstix    space allocated for deg+1 doubles;
 *   cstmi    space allocated for deg+1 doubles;
 *   cstrg    space allocated for deg+1 doubles;
 *   cstpk    space allocated for deg+1 doubles;
 *   cfftb    space allocated for nbr arrays of deg+1 doubles;
 *   cffix    space allocated for nbr arrays of deg+1 doubles;
 *   cffmi    space allocated for nbr arrays of deg+1 doubles;
 *   cffrg    space allocated for nbr arrays of deg+1 doubles;
 *   cffpk    space allocated for nbr arrays of deg+1 doubles.
 *
 * ON RETURN :
 *   idx      the participating variables in each monomial,
 *            idx[k] is an array of nvr[k] integers,
 *            idx[k][i] is the index of variable i in monomial k;
 *   exp      the exponents of the variables in each monomial,
 *            exp[k] is an array of nvr[k] integers,
 *            exp[k][i] is the exponent of variable i in monomial k;
 *   csttb    deg+1 highest coefficients of the constant monomial;
 *   cstix    deg+1 second highest coefficients of the constant monomial;
 *   cstmi    deg+1 middle coefficients of the constant monomial;
 *   cstrg    deg+1 second lowest coefficients of the constant monomial;
 *   cstpk    deg+1 lowest coefficients of the constant monomial;
 *   cfftb    the highest parts of the coefficient series for each monomial,
 *            cfftb[k] has the deg+1 highest coefficients of monomial k;
 *   cffix    second highest parts of the coefficient series, cffix[k]
 *            has the deg+1 second highest coefficients of monomial k;
 *   cffmi    middle parts of the coefficient series,
 *            cffix[k] has the deg+1 middle coefficients of monomial k;
 *   cffrg    second lowest parts of the coefficient series, cffrg[k]
 *            has the deg+1 second lowest coefficients of monomial k;
 *   cffpk    the lowest parts of the coefficient series for each monomial,
 *            cffpk[k] has the deg+1 lowest coefficients of monomial k. */

void random_penta_double_complex
 ( double *retb, double *reix, double *remi, double *rerg, double *repk,
   double *imtb, double *imix, double *immi, double *imrg, double *impk );
/*
 * DESCRIPTION :
 *   Generates the high, middle, and low double of a random cosine c
 *   and computes the corresponding sine as sqrt(1-c*c),
 *   so the randomly generated complex number has modulus one.
 *
 * ON RETURN :
 *   retb   the highest double of the real part,
 *   reix   the second highest double of the real part,
 *   remi   the middle double of the real part,
 *   rerg   the second lowest double of the real part,
 *   repk   the lowest double of the real part,
 *   imtb   the highest double of the imaginary part,
 *   imix   the second highest double of the imaginary part,
 *   immi   the middle double of the imaginary part,
 *   imrg   the second lowest double of the imaginary part,
 *   impk   the lowest double of the imaginary part
 *            of the randomly generated complex number. */

bool make_complex5_polynomial
 ( int dim, int nbr, int pwr, int deg, int *nvr, int **idx, int **exp,
   double *cstretb, double *cstreix, double *cstremi,
   double *cstrerg, double *cstrepk,
   double *cstimtb, double *cstimix, double *cstimmi,
   double *cstimrg, double *cstimpk,
   double **cffretb, double **cffreix, double **cffremi,
   double **cffrerg, double **cffrepk,
   double **cffimtb, double **cffimix, double **cffimmi,
   double **cffimrg, double **cffimpk );
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
 *   cstretb   has space allocated for deg+1 doubles;
 *   cstreix   has space allocated for deg+1 doubles;
 *   cstremi   has space allocated for deg+1 doubles;
 *   cstrerg   has space allocated for deg+1 doubles;
 *   cstrepk   has space allocated for deg+1 doubles;
 *   cstimtb   has space allocated for deg+1 doubles;
 *   cstimix   has space allocated for deg+1 doubles;
 *   cstimmi   has space allocated for deg+1 doubles;
 *   cstimrg   has space allocated for deg+1 doubles;
 *   cstimpk   has space allocated for deg+1 doubles;
 *   cffretb   has space allocated for nbr arrays of deg+1 doubles;
 *   cffreix   has space allocated for nbr arrays of deg+1 doubles;
 *   cffremi   has space allocated for nbr arrays of deg+1 doubles;
 *   cffrerg   has space allocated for nbr arrays of deg+1 doubles;
 *   cffrepk   has space allocated for nbr arrays of deg+1 doubles;
 *   cffimtb   has space allocated for nbr arrays of deg+1 doubles;
 *   cffimix   has space allocated for nbr arrays of deg+1 doubles;
 *   cffimmi   has space allocated for nbr arrays of deg+1 doubles;
 *   cffimrg   has space allocated for nbr arrays of deg+1 doubles.
 *   cffimpk   has space allocated for nbr arrays of deg+1 doubles.
 *
 * ON RETURN :
 *   idx       the participating variables in each monomial,
 *             idx[k] is an array of nvr[k] integers,
 *             idx[k][i] is the index of variable i in monomial k;
 *   exp       the exponents of the variables in each monomial,
 *             exp[k] is an array of nvr[k] integers,
 *             exp[k][i] is the exponent of variable i in monomial k;
 *   cstretb   has the deg+1 highest coefficients of the real parts
 *             of the constant monomial;
 *   cstreix   has the deg+1 second highest coefficients of the real parts
 *             of the constant monomial;
 *   cstremi   has the deg+1 middle highest coefficients of the real parts
 *             of the constant monomial;
 *   cstrerg   has the deg+1 second lowest coefficients of the real parts
 *             of the constant monomial;
 *   cstrepk   has the deg+1 lowest coefficients of the real parts
 *             of the constant monomial;
 *   cstimtb   has the deg+1 highest coefficients of the imaginary parts
 *             of the constant monomial;
 *   cstimix   has the deg+1 second highest coefficients of the imaginary parts
 *             of the constant monomial;
 *   cstimmi   has the deg+1 middle coefficients of the imaginary parts
 *             of the constant monomial;
 *   cstimrg   has the deg+1 second lowest coefficients of the imaginary parts
 *             of the constant monomial;
 *   cstimpk   has the deg+1 lowest coefficients of the imaginary parts
 *             of the constant monomial;
 *   cffretb   has the highest doubles of the real parts of the coefficients,
 *             cffretb[k] has deg+1 highest coefficients of monomial k;
 *   cffreix   has the second highest doubles of the real parts
 *             of the coefficients, cffreix[k] has deg+1 second highest
 *             coefficients of monomial k;
 *   cffremi   has the middle doubles of the real parts
 *             of the coefficients, cffreix[k] has deg+1 second highest
 *             coefficients of monomial k;
 *   cffrerg   has the second lowest doubles of the real parts
 *             of the coefficients, cffrerg[k] has deg+1 second lowest
 *             coefficients of monomial k;
 *   cffrepk   has the lowest doubles of the real parts of the coefficients,
 *             cffrepk[k] has deg+1 lowest coefficients of monomial k;
 *   cffimtb   has the highest doubles of the imaginary parts
 *             of the coefficients, cffimtb[k] has deg+1 highest
 *             coefficients of monomial k;
 *   cffimix   has the second highest doubles of the imaginary parts
 *             of the coefficients, cffimix[k] has deg+1 second highest
 *             coefficients of monomial k;
 *   cffimmi   has the middle doubles of the imaginary parts
 *             of the coefficients, cffimix[k] has deg+1 second highest
 *             coefficients of monomial k;
 *   cffimrg   has the second lowest doubles of the imaginary parts
 *             of the coefficient, cffimpk[k] has the deg+1 second lowest
 *             coefficients of monomial k;
 *   cffimpk   has the lowest doubles of the imaginary parts
 *             of the coefficient, cffimpk[k] has the deg+1 lowest
 *             coefficients of monomial k. */

void make_real5_products
 ( int dim, int nbr, int nva, int deg, int **idx,
   double *csttb, double *cstix, double *cstmi, double *cstrg, double *cstpk,
   double **cfftb, double **cffix, double **cffmi, double **cffrg,
   double **cffpk );
/*
 * DESCRIPTION :
 *   Returns the sum of all products of size nbr out of dimension dim,
 *   with random power series coefficients truncated to degree deg.
 *
 * REQUIRED : dim > nbr.
 *
 * ON ENTRY :
 *   dim      dimension, total number of variables;
 *   nbr      number of monomials, excluding the constant;
 *   nva      number of variables in each monomial;
 *   deg      truncation degree of the power series;
 *   idx      space for nbr index vectors;
 *   csttb    space allocated for deg+1 doubles;
 *   cstix    space allocated for deg+1 doubles;
 *   cstmi    space allocated for deg+1 doubles;
 *   cstrg    space allocated for deg+1 doubles;
 *   cstpk    space allocated for deg+1 doubles;
 *   cfftb    space allocated for nbr arrays of deg+1 doubles;
 *   cffix    space allocated for nbr arrays of deg+1 doubles;
 *   cffmi    space allocated for nbr arrays of deg+1 doubles;
 *   cffrg    space allocated for nbr arrays of deg+1 doubles;
 *   cffpk    space allocated for nbr arrays of deg+1 doubles.
 *
 * ON RETURN :
 *   idx      the participating variables in each monomial,
 *            idx[k] is an array of nva integers,
 *            idx[k][i] is the index of variable i in monomial k;
 *   csttb    deg+1 highest coefficients of the constant monomial;
 *   cstix    deg+1 second highest coefficients of the constant monomial;
 *   cstmi    deg+1 middle coefficients of the constant monomial;
 *   cstrg    deg+1 second lowest coefficients of the constant monomial;
 *   cstpk    deg+1 lowest coefficients of the constant monomial;
 *   cfftb    the highest parts of the coefficient series for each monomial,
 *            cfftb[k] has the deg+1 highest coefficients of monomial k;
 *   cffix    second highest parts of the coefficient series, cffix[k]
 *            has the deg+1 second highest coefficients of monomial k;
 *   cffmi    middle parts of the coefficient series,
 *            cffmi[k] has the deg+1 middle coefficients of monomial k;
 *   cffrg    second lowest parts of the coefficient series, cffrg[k]
 *            has the deg+1 second lowest coefficients of monomial k;
 *   cffpk    the lowest parts of the coefficient series for each monomial,
 *            cffpk[k] has the deg+1 lowest coefficients of monomial k. */

void make_complex5_products
 ( int dim, int nbr, int nva, int deg, int **idx,
   double *cstretb, double *cstreix, double *cstremi,
   double *cstrerg, double *cstrepk,
   double *cstimtb, double *cstimix, double *cstimmi,
   double *cstimrg, double *cstimpk,
   double **cffretb, double **cffreix, double **cffremi,
   double **cffrerg, double **cffrepk,
   double **cffimtb, double **cffimix, double **cffimmi,
   double **cffimrg, double **cffimpk );
/*
 * DESCRIPTION :
 *   Returns the sum of all products of size nbr out of dimension dim,
 *   with random complex power series coefficients truncated to degree deg.
 *
 * REQUIRED : dim > nbr.
 *
 * ON ENTRY :
 *   dim      dimension, total number of variables;
 *   nbr      number of monomials, excluding the constant term,
 *            nbr is the number of monomials that have at least one variable;
 *   pwr      largest power of a variable;
 *   deg      degree of the power series coefficient;
 *   nvr      nbr integers, with nvr[k] the number of variables in monomial k;
 *   exp      space allocated for nbr arrays of nvr[k] integers;
 *   csttb    space allocated for deg+1 doubles;
 *   cstix    space allocated for deg+1 doubles;
 *   cstmi    space allocated for deg+1 doubles;
 *   cstrg    space allocated for deg+1 doubles;
 *   cstpk    space allocated for deg+1 doubles;
 *   cfftb    space allocated for nbr arrays of deg+1 doubles;
 *   cffix    space allocated for nbr arrays of deg+1 doubles;
 *   cffmi    space allocated for nbr arrays of deg+1 doubles;
 *   cffrg    space allocated for nbr arrays of deg+1 doubles;
 *   cffpk    space allocated for nbr arrays of deg+1 doubles.
 *
 * ON RETURN :
 *   idx      the participating variables in each monomial,
 *            idx[k] is an array of nva integers,
 *            idx[k][i] is the index of variable i in monomial k;
 *   csttb    deg+1 highest coefficients of the constant monomial;
 *   cstix    deg+1 second highest coefficients of the constant monomial;
 *   cstmi    deg+1 middle coefficients of the constant monomial;
 *   cstrg    deg+1 second lowest coefficients of the constant monomial;
 *   cstpk    deg+1 lowest coefficients of the constant monomial;
 *   cfftb    the highest parts of the coefficient series for each monomial,
 *            cfftb[k] has the deg+1 highest coefficients of monomial k;
 *   cffix    second highest parts of the coefficient series, cffix[k]
 *            has the deg+1 second highest coefficients of monomial k;
 *   cffmi    middle parts of the coefficient series,
 *            cffmi[k] has the deg+1 middle coefficients of monomial k;
 *   cffrg    second lowest parts of the coefficient series, cffrg[k]
 *            has the deg+1 second lowest coefficients of monomial k;
 *   cffpk    the lowest parts of the coefficient series for each monomial,
 *            cffpk[k] has the deg+1 lowest coefficients of monomial k. */

void make_real5_cyclic
 ( int dim, int nva, int deg, int **idx,
   double *csttb, double *cstix, double *cstmi, double *cstrg, double *cstpk,
   double **cfftb, double **cffix, double **cffmi, double **cffrg,
   double **cffpk );
/*
 * DESCRIPTION :
 *   Returns the cyclic polynomial with nva variables in dimension dim,
 *   with random power series coefficients truncated to degree deg.
 *
 * REQUIRED : dim > nbr.
 *
 * ON ENTRY :
 *   dim      dimension, total number of variables;
 *   nva      number of variables in each monomial;
 *   deg      truncation degree of the power series;
 *   idx      space for nva index vectors;
 *   csttb    space allocated for deg+1 doubles;
 *   cstix    space allocated for deg+1 doubles;
 *   cstmi    space allocated for deg+1 doubles;
 *   cstrg    space allocated for deg+1 doubles;
 *   cstpk    space allocated for deg+1 doubles;
 *   cfftb    space allocated for nbr arrays of deg+1 doubles;
 *   cffix    space allocated for nbr arrays of deg+1 doubles;
 *   cffmi    space allocated for nbr arrays of deg+1 doubles;
 *   cffrg    space allocated for nbr arrays of deg+1 doubles;
 *   cffpk    space allocated for nbr arrays of deg+1 doubles.
 *
 * ON RETURN :
 *   idx      the participating variables in each monomial,
 *            idx[k] is an array of nva integers,
 *            idx[k][i] is the index of variable i in monomial k;
 *   csttb    deg+1 highest coefficients of the constant monomial;
 *   cstix    deg+1 second highest coefficients of the constant monomial;
 *   cstrg    deg+1 second lowest coefficients of the constant monomial;
 *   cstpk    deg+1 lowest coefficients of the constant monomial;
 *   cfftb    the highest parts of the coefficient series for each monomial,
 *            cfftb[k] has the deg+1 highest coefficients of monomial k;
 *   cffix    second highest parts of the coefficient series, cffix[k]
 *            has the deg+1 second highest coefficients of monomial k;
 *   cffmi    middle parts of the coefficient series,
 *            cffix[k] has the deg+1 middle coefficients of monomial k;
 *   cffrg    second lowest parts of the coefficient series, cffrg[k]
 *            has the deg+1 second lowest coefficients of monomial k;
 *   cffpk    the lowest parts of the coefficient series for each monomial,
 *            cffpk[k] has the deg+1 lowest coefficients of monomial k. */

void make_complex5_cyclic
 ( int dim, int nva, int deg, int **idx,
   double *cstretb, double *cstreix, double *cstremi,
   double *cstrerg, double *cstrepk,
   double *cstimtb, double *cstimix, double *cstimmi,
   double *cstimrg, double *cstimpk,
   double **cffretb, double **cffreix, double **cffremi, 
   double **cffrerg, double **cffrepk,
   double **cffimtb, double **cffimix, double **cffimmi,
   double **cffimrg, double **cffimpk );
/*
 * DESCRIPTION :
 *   Returns the cyclic polynomial with nva variables in dimension dim,
 *   with random complex power series coefficients truncated to degree deg.
 *
 * REQUIRED : dim > nbr.
 *
 * ON ENTRY :
 *   dim      dimension, total number of variables;
 *   nva      number of variables in each monomial;
 *   deg      truncation degree of the power series;
 *   idx      space for nva index vectors;
 *   csttb    space allocated for deg+1 doubles;
 *   cstix    space allocated for deg+1 doubles;
 *   cstmi    space allocated for deg+1 doubles;
 *   cstrg    space allocated for deg+1 doubles;
 *   cstpk    space allocated for deg+1 doubles;
 *   cfftb    space allocated for nbr arrays of deg+1 doubles;
 *   cffix    space allocated for nbr arrays of deg+1 doubles;
 *   cffmi    space allocated for nbr arrays of deg+1 doubles;
 *   cffrg    space allocated for nbr arrays of deg+1 doubles;
 *   cffpk    space allocated for nbr arrays of deg+1 doubles.
 *
 * ON RETURN :
 *   idx      the participating variables in each monomial,
 *            idx[k] is an array of nva integers,
 *            idx[k][i] is the index of variable i in monomial k;
 *   csttb    deg+1 highest coefficients of the constant monomial;
 *   cstix    deg+1 second highest coefficients of the constant monomial;
 *   cstrg    deg+1 second lowest coefficients of the constant monomial;
 *   cstpk    deg+1 lowest coefficients of the constant monomial;
 *   cfftb    the highest parts of the coefficient series for each monomial,
 *            cfftb[k] has the deg+1 highest coefficients of monomial k;
 *   cffix    second highest parts of the coefficient series, cffix[k]
 *            has the deg+1 second highest coefficients of monomial k;
 *   cffmi    middle parts of the coefficient series,
 *            cffix[k] has the deg+1 middle coefficients of monomial k;
 *   cffrg    second lowest parts of the coefficient series, cffrg[k]
 *            has the deg+1 second lowest coefficients of monomial k;
 *   cffpk    the lowest parts of the coefficient series for each monomial,
 *            cffpk[k] has the deg+1 lowest coefficients of monomial k. */

#endif
