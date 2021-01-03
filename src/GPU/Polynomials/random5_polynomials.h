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

#endif
