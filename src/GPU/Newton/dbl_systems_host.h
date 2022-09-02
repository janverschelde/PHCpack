// The file dbl_systems_host.h specifies functions to evaluate and
// differentiate monomials with common factors in double precision.

#ifndef __dbl_systems_host_h__
#define __dbl_systems_host_h__

void CPU_dbl_evaluate_monomials
 ( int dim, int deg, int *nvr, int **idx, int **exp, int *nbrfac,
   int **expfac, double **cff, double *acc, double **input,
   double ***output, int vrblvl );
/*
 * DESCRIPTION :
 *   Evaluates monomials at power series.
 *
 * ON ENTRY :
 *   dim       number of monomials;
 *   deg       degree of the power series;
 *   nvr       nvr[i] is the number of variables in the i-th monomial;
 *   idx       idx[i] are the indices of the variables in monomial i;
 *   exp       exp[i] are the exponents of the variables in monomial i;
 *   nbrfac    nbrfac[i] are the number of exponents > 1 in monomial i;
 *   expfac    expfac[i] are the exponents in the i-th polynomial
 *             that are larger than one, minus one in the factor,
 *             if exp[i][k] > 1, then expfac[i][k] = exp[i][k] - 1;
 *   cff       coefficients of the monomials;
 *   acc       space to accumulate one power series of degree deg;
 *   input     coefficients of the power series of degree deg,
 *             for dim variables;
 *   output    space for the output;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   cff       contains the evaluated common factors;
 *   output    evaluated and differentiated monomials in the system,
 *             output[i][dim] is the value of the input series
 *             at the i-th monomial, and for k in range 0..nvr[i]-1,
 *             output[i][idx[i]] is the derivative w.r.t. idx[k]. */

void dbl_linearize_evaldiff_output
 ( int dim, int degp1, int *nvr, int **idx, double ***output,
   double **funval, double **rhs, double ***jacval, int vrblvl );
/*
 * DESCRIPTION :
 *   Linearizes the output of the evaluation and differentiation
 *   of the monomials.
 *
 * ON ENTRY :
 *   dim       number of monomials;
 *   degp1     degree plus one;
 *   nvr       number of variables that occur in each monomial;
 *   idx       for each monomials the indices of each variable;
 *   output    output of the evaluation and differentiation
 *             of the monomials in the system, for the i-th monomial:
 *             output[i][dim] is the power series value, 
 *             output[i][idx[k]] is the derivative w.r.t. idx[k],
 *             for k in range 0..nvr[i]-1;
 *   funval    space allocated for dim power series;
 *   rhs       space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   jacval    space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   funval    collects the output[i][dim], the evaluated series;
 *   rhs       the linearized right hand side are the function values
 *             subtracted by 1 and added by t;
 *   jacval    a series with matrices as coefficients,
 *             the leading coefficient is the Jacobian matrix. */

#endif
