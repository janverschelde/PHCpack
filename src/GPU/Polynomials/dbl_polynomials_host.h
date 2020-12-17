/* The file dbl_polynomials_host.h specifies functions to evaluate and
 * differentiate a polynomial at power series truncated to the same degree,
 * in double precision. */

#ifndef __dbl_polynomials_host_h__
#define __dbl_polynomials_host_h__

void CPU_dbl_poly_speel
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double **cff, double **input, double **output,
   double **forward, double **backward, double **cross );
/*
 * DESCRIPTION :
 *   Runs the reverse mode of algorithmic differentiation
 *   of a polynomial at power series truncated to the same degree,
 *   for real coefficients in double precision.
 *
 * REQUIRED : nbr > 1.
 *
 * ON ENTRY :
 *   dim      total number of variables;
 *   nbr      number of monomials, including the constant term,
 *            if nbr == 1, then there is no constant term,
 *            only one monomial multiplied with a coefficient series;
 *   deg      truncation degree of the series;
 *   nvr      nvr[k] holds the number of variables in monomial k;
 *   idx      idx[k] has as many indices as the value of nvr[k],
 *            idx[k][i] defines the place of the i-th variable,
 *            with input values in input[idx[k][i]];
 *   cff      cff[k] has deg+1 doubles for the coefficient series
 *            of monomial k;
 *   input    contains the coefficients of the power series
 *            for all variables in the polynomial;
 *   output   space allocated for the value and all derivatives;
 *   forward  contains work space for all nvr forward products,
 *            forward[k] contains space for deg+1 doubles;
 *   backward contains work space for all nvr-2 backward products;
 *            backward[k] contains space for deg+1 doubles;
 *   cross    contains work space for all nvr-2 cross products;
 *            cross[k] contains space for deg+1 doubles.
 *
 * ON RETURN :
 *   output   contains derivatives and the value of the polynomial,
 *            output[k], for k from 0 to dim-1, contains the
 *            derivative with respect to the variable k;
 *            output[dim] contains the value of the polynomial. */

void CPU_dbl_poly_evaldiff
 ( int dim, int nbr, int deg, int *nvr, int **idx, 
   double *cst, double **cff, double **input, double **output );
/*
 * DESCRIPTION :
 *   Allocates work space memory to store the forward, backward, and
 *   cross products in the evaluation and differentiation of a polynomial.
 *
 * ON ENTRY :
 *   dim      total number of variables;
 *   nbr      number of monomials, including the constant term,
 *            if nbr == 1, then there is no constant term,
 *            only one monomial multiplied with a coefficient series;
 *   deg      truncation degree of the series;
 *   nvr      nvr[k] holds the number of variables in monomial k;
 *   idx      idx[k] has as many indices as the value of nvr[k],
 *            idx[k][i] defines the place of the i-th variable,
 *            with input values in input[idx[k][i]];
 *   cst      deg+1 doubles for the constant coefficient series;
 *   cff      cff[k] has deg+1 doubles for the coefficient series
 *            of monomial k;
 *   input    contains the coefficients of the power series
 *            for all variables in the polynomial;
 *   output   space allocated for the value and all derivatives.
 *
 * ON RETURN :
 *   output   contains derivatives and the value of the polynomial,
 *            output[k], for k from 0 to dim-1, contains the
 *            derivative with respect to the variable k;
 *            output[dim] contains the value of the polynomial. */

#endif
