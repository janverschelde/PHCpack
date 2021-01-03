// The file random10_polynomials.h specifies functions to setup tests for the 
// evaluation and differentiation of a polynomial in deca double precision.

#ifndef __random10_polynomials_h__
#define __random10_polynomials_h__

bool make_real10_polynomial
 ( int dim, int nbr, int pwr, int deg, int *nvr, int **idx, int **exp,
   double *cstrtb, double *cstrix, double *cstrmi, double *cstrrg,
   double *cstrpk, double *cstltb, double *cstlix, double *cstlmi,
   double *cstlrg, double *cstlpk,
   double **cffrtb, double **cffrix, double **cffrmi, double **cffrrg,
   double **cffrpk, double **cffltb, double **cfflix, double **cfflmi,
   double **cfflrg, double **cfflpk );
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
 *   cstrtb   space allocated for deg+1 doubles;
 *   cstrix   space allocated for deg+1 doubles;
 *   cstrmi   space allocated for deg+1 doubles;
 *   cstrrg   space allocated for deg+1 doubles;
 *   cstrpk   space allocated for deg+1 doubles;
 *   cstltb   space allocated for deg+1 doubles;
 *   cstlix   space allocated for deg+1 doubles;
 *   cstlmi   space allocated for deg+1 doubles;
 *   cstlrg   space allocated for deg+1 doubles;
 *   cstlpk   space allocated for deg+1 doubles;
 *   cffrtb   space allocated for nbr arrays of deg+1 doubles;
 *   cffrix   space allocated for nbr arrays of deg+1 doubles;
 *   cffrmi   space allocated for nbr arrays of deg+1 doubles;
 *   cffrrg   space allocated for nbr arrays of deg+1 doubles;
 *   cffrpk   space allocated for nbr arrays of deg+1 doubles;
 *   cffltb   space allocated for nbr arrays of deg+1 doubles;
 *   cfflix   space allocated for nbr arrays of deg+1 doubles;
 *   cfflmi   space allocated for nbr arrays of deg+1 doubles;
 *   cfflrg   space allocated for nbr arrays of deg+1 doubles;
 *   cfflpk   space allocated for nbr arrays of deg+1 doubles.
 *
 * ON RETURN :
 *   idx      the participating variables in each monomial,
 *            idx[k] is an array of nvr[k] integers,
 *            idx[k][i] is the index of variable i in monomial k;
 *   exp      the exponents of the variables in each monomial,
 *            exp[k] is an array of nvr[k] integers,
 *            exp[k][i] is the exponent of variable i in monomial k;
 *   cstrtb   deg+1 highest coefficients of the constant monomial;
 *   cstrix   deg+1 second highest coefficients of the constant monomial;
 *   cstrmi   deg+1 third highest coefficients of the constant monomial;
 *   cstrrg   deg+1 fourth highest coefficients of the constant monomial;
 *   cstrpk   deg+1 fifth highest coefficients of the constant monomial;
 *   cstltb   deg+1 fifth lowest coefficients of the constant monomial;
 *   cstlix   deg+1 fourth lowest coefficients of the constant monomial;
 *   cstlmi   deg+1 third lowest coefficients of the constant monomial;
 *   cstlrg   deg+1 second lowest coefficients of the constant monomial;
 *   cstlpk   deg+1 lowest coefficients of the constant monomial;
 *   cffrtb   the highest parts of the coefficient series for each monomial,
 *            cffrtb[k] has the deg+1 highest coefficients of monomial k;
 *   cffrix   second highest parts of the coefficient series, cffrix[k]
 *            has the deg+1 second highest coefficients of monomial k;
 *   cffrmi   third highest parts of the coefficient series, cffrmi[k]
 *            has the deg+1 third highest coefficients of monomial k;
 *   cffrrg   fourth highest parts of the coefficient series, cffrrg[k]
 *            has the deg+1 fourth highest coefficients of monomial k;
 *   cffrpk   fifth highest parts of the coefficient series, cffrpk[k]
 *            has the deg+1 fifth highest coefficients of monomial k;
 *   cffltb   fifth lowest parts of the coefficient series, cffltb[k]
 *            has the deg+1 fifth lowest coefficients of monomial k;
 *   cfflix   fourth lowest parts of the coefficient series, cfflix[k]
 *            has the deg+1 fourth lowest coefficients of monomial k;
 *   cfflmi   third lowest parts of the coefficient series, cfflmi[k]
 *            has the deg+1 third lowest coefficients of monomial k;
 *   cfflrg   second lowest parts of the coefficient series, cfflrg[k]
 *            has the deg+1 second lowest coefficients of monomial k;
 *   cfflpk   the lowest parts of the coefficient series for each monomial,
 *            cfflpk[k] has the deg+1 lowest coefficients of monomial k. */

void make_real10_products
 ( int dim, int nbr, int nva, int deg, int **idx,
   double *cstrtb, double *cstrix, double *cstrmi, double *cstrrg,
   double *cstrpk, double *cstltb, double *cstlix, double *cstlmi,
   double *cstlrg, double *cstlpk,
   double **cffrtb, double **cffrix, double **cffrmi, double **cffrrg,
   double **cffrpk, double **cffltb, double **cfflix, double **cfflmi,
   double **cfflrg, double **cfflpk );
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
 *   cstrtb   space allocated for deg+1 doubles;
 *   cstrix   space allocated for deg+1 doubles;
 *   cstrmi   space allocated for deg+1 doubles;
 *   cstrrg   space allocated for deg+1 doubles;
 *   cstrpk   space allocated for deg+1 doubles;
 *   cstltb   space allocated for deg+1 doubles;
 *   cstlix   space allocated for deg+1 doubles;
 *   cstlmi   space allocated for deg+1 doubles;
 *   cstlrg   space allocated for deg+1 doubles;
 *   cstlpk   space allocated for deg+1 doubles;
 *   cffrtb   space allocated for nbr arrays of deg+1 doubles;
 *   cffrix   space allocated for nbr arrays of deg+1 doubles;
 *   cffrmi   space allocated for nbr arrays of deg+1 doubles;
 *   cffrrg   space allocated for nbr arrays of deg+1 doubles;
 *   cffrpk   space allocated for nbr arrays of deg+1 doubles;
 *   cffltb   space allocated for nbr arrays of deg+1 doubles;
 *   cfflix   space allocated for nbr arrays of deg+1 doubles;
 *   cfflmi   space allocated for nbr arrays of deg+1 doubles;
 *   cfflrg   space allocated for nbr arrays of deg+1 doubles;
 *   cfflpk   space allocated for nbr arrays of deg+1 doubles.
 *
 * ON RETURN :
 *   idx      the participating variables in each monomial,
 *            idx[k] is an array of nva integers,
 *            idx[k][i] is the index of variable i in monomial k;
 *   cstrtb   deg+1 highest coefficients of the constant monomial;
 *   cstrix   deg+1 second highest coefficients of the constant monomial;
 *   cstrmi   deg+1 third highest coefficients of the constant monomial;
 *   cstrrg   deg+1 fourth highest coefficients of the constant monomial;
 *   cstrpk   deg+1 fifth highest coefficients of the constant monomial;
 *   cstltb   deg+1 fifth lowest coefficients of the constant monomial;
 *   cstlix   deg+1 fourth lowest coefficients of the constant monomial;
 *   cstlmi   deg+1 third lowest coefficients of the constant monomial;
 *   cstlrg   deg+1 second lowest coefficients of the constant monomial;
 *   cstlpk   deg+1 lowest coefficients of the constant monomial;
 *   cffrtb   the highest parts of the coefficient series for each monomial,
 *            cffrtb[k] has the deg+1 highest coefficients of monomial k;
 *   cffrix   second highest parts of the coefficient series, cffrix[k]
 *            has the deg+1 second highest coefficients of monomial k;
 *   cffrmi   third highest parts of the coefficient series, cffrmi[k]
 *            has the deg+1 third highest coefficients of monomial k;
 *   cffrrg   fourth highest parts of the coefficient series, cffrrg[k]
 *            has the deg+1 fourth highest coefficients of monomial k;
 *   cffrpk   fifth highest parts of the coefficient series, cffrpk[k]
 *            has the deg+1 fifth highest coefficients of monomial k;
 *   cffltb   fifth lowest parts of the coefficient series, cffltb[k]
 *            has the deg+1 fifth lowest coefficients of monomial k;
 *   cfflix   fourth lowest parts of the coefficient series, cfflix[k]
 *            has the deg+1 fourth lowest coefficients of monomial k;
 *   cfflmi   third lowest parts of the coefficient series, cfflmi[k]
 *            has the deg+1 third lowest coefficients of monomial k;
 *   cfflrg   second lowest parts of the coefficient series, cfflrg[k]
 *            has the deg+1 second lowest coefficients of monomial k;
 *   cfflpk   the lowest parts of the coefficient series for each monomial,
 *            cfflpk[k] has the deg+1 lowest coefficients of monomial k. */

void make_real10_cyclic
 ( int dim, int nva, int deg, int **idx,
   double *cstrtb, double *cstrix, double *cstrmi, double *cstrrg,
   double *cstrpk, double *cstltb, double *cstlix, double *cstlmi,
   double *cstlrg, double *cstlpk,
   double **cffrtb, double **cffrix, double **cffrmi, double **cffrrg,
   double **cffrpk, double **cffltb, double **cfflix, double **cfflmi,
   double **cfflrg, double **cfflpk );
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
 *   cstrtb   space allocated for deg+1 doubles;
 *   cstrix   space allocated for deg+1 doubles;
 *   cstrmi   space allocated for deg+1 doubles;
 *   cstrrg   space allocated for deg+1 doubles;
 *   cstrpk   space allocated for deg+1 doubles;
 *   cstltb   space allocated for deg+1 doubles;
 *   cstlix   space allocated for deg+1 doubles;
 *   cstlmi   space allocated for deg+1 doubles;
 *   cstlrg   space allocated for deg+1 doubles;
 *   cstlpk   space allocated for deg+1 doubles;
 *   cffrtb   space allocated for nbr arrays of deg+1 doubles;
 *   cffrix   space allocated for nbr arrays of deg+1 doubles;
 *   cffrmi   space allocated for nbr arrays of deg+1 doubles;
 *   cffrrg   space allocated for nbr arrays of deg+1 doubles;
 *   cffrpk   space allocated for nbr arrays of deg+1 doubles;
 *   cffltb   space allocated for nbr arrays of deg+1 doubles;
 *   cfflix   space allocated for nbr arrays of deg+1 doubles;
 *   cfflmi   space allocated for nbr arrays of deg+1 doubles;
 *   cfflrg   space allocated for nbr arrays of deg+1 doubles;
 *   cfflpk   space allocated for nbr arrays of deg+1 doubles.
 *
 * ON RETURN :
 *   idx      the participating variables in each monomial,
 *            idx[k] is an array of nva integers,
 *            idx[k][i] is the index of variable i in monomial k;
 *   cstrtb   deg+1 highest coefficients of the constant monomial;
 *   cstrix   deg+1 second highest coefficients of the constant monomial;
 *   cstrmi   deg+1 third highest coefficients of the constant monomial;
 *   cstrrg   deg+1 fourth highest coefficients of the constant monomial;
 *   cstrpk   deg+1 fifth highest coefficients of the constant monomial;
 *   cstltb   deg+1 fifth lowest coefficients of the constant monomial;
 *   cstlix   deg+1 fourth lowest coefficients of the constant monomial;
 *   cstlmi   deg+1 third lowest coefficients of the constant monomial;
 *   cstlrg   deg+1 second lowest coefficients of the constant monomial;
 *   cstlpk   deg+1 lowest coefficients of the constant monomial;
 *   cffrtb   the highest parts of the coefficient series for each monomial,
 *            cffrtb[k] has the deg+1 highest coefficients of monomial k;
 *   cffrix   second highest parts of the coefficient series, cffrix[k]
 *            has the deg+1 second highest coefficients of monomial k;
 *   cffrmi   third highest parts of the coefficient series, cffrmi[k]
 *            has the deg+1 third highest coefficients of monomial k;
 *   cffrrg   fourth highest parts of the coefficient series, cffrrg[k]
 *            has the deg+1 fourth highest coefficients of monomial k;
 *   cffrpk   fifth highest parts of the coefficient series, cffrpk[k]
 *            has the deg+1 fifth highest coefficients of monomial k;
 *   cffltb   fifth lowest parts of the coefficient series, cffltb[k]
 *            has the deg+1 fifth lowest coefficients of monomial k;
 *   cfflix   fourth lowest parts of the coefficient series, cfflix[k]
 *            has the deg+1 fourth lowest coefficients of monomial k;
 *   cfflmi   third lowest parts of the coefficient series, cfflmi[k]
 *            has the deg+1 third lowest coefficients of monomial k;
 *   cfflrg   second lowest parts of the coefficient series, cfflrg[k]
 *            has the deg+1 second lowest coefficients of monomial k;
 *   cfflpk   the lowest parts of the coefficient series for each monomial,
 *            cfflpk[k] has the deg+1 lowest coefficients of monomial k. */

#endif
