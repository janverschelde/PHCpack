// The file dbl2_systems_host.h specifies functions to evaluate and
// differentiate monomials with common factors in double double precision.

#ifndef __dbl2_systems_host_h__
#define __dbl2_systems_host_h__

void CPU_dbl2_evaluate_monomials
 ( int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffhi, double **cfflo, double *acchi, double *acclo,
   double **inputhi, double **inputlo,
   double ***outputhi, double ***outputlo, int vrblvl );
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
 *   cffhi     high doubles of the coefficients of the monomials;
 *   cfflo     low doubles of the coefficients of the monomials;
 *   acchi     space to accumulate one power series of degree deg;
 *   acclo     space to accumulate one power series of degree deg;
 *   inputhi   high doubles of the coefficients of the power series
 *             of degree deg, for dim variables;
 *   inputhi   low doubles of the coefficients of the power series
 *             of degree deg, for dim variables;
 *   outputhi  space for the high doubles of the output;
 *   outputlo  space for the low doubles of the output;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   cffhi     high doubles of the evaluated common factors;
 *   cfflo     low doubles of the evaluated common factors;
 *   outputhi  high doubles of the evaluated and differentiated monomials,
 *   outputlo  low doubles of the evaluated and differentiated monomials,
 *             output[i][dim] is the value of the input series
 *             at the i-th monomial, and for k in range 0..nvr[i]-1,
 *             output[i][idx[i]] is the derivative w.r.t. idx[k]. */

void dbl2_linearize_evaldiff_output
 ( int dim, int degp1, int *nvr, int **idx,
   double ***outputhi, double ***outputlo,
   double **funvalhi, double **funvallo, 
   double **rhshi, double **rhslo, double ***jacvalhi, double ***jacvallo,
   int vrblvl );
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
 *   outputhi  high part of the output of the evaluation and differentiation
 *             of the monomials in the system, for the i-th monomial:
 *             output[i][dim] is the power series value, 
 *             output[i][idx[k]] is the derivative w.r.t. idx[k],
 *             for k in range 0..nvr[i]-1;
 *   outputlo  low parts of the output of the evaluation and differentiation;
 *   funvalhi  space allocated for dim power series;
 *   funvallo  space allocated for dim power series;
 *   rhshi     space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   rhslo     space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   jacvalhi  space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   jacvallo  space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   funvalhi  collects the high parts of output[i][dim], the evaluated series;
 *   funvallo  collects the low parts of output[i][dim], the evaluated series;
 *   rhshi     the linearized right hand side are the high doubles of the
 *             function values subtracted by 1 and added by t;
 *   rhslo     the linearized right hand side are the low doubles of the
 *             function values subtracted by 1 and added by t;
 *   jacvalhi  a series with matrices as coefficients, of high doubles,
 *             the leading coefficient is the Jacobian matrix;
 *   jacvallo  a series with matrices as coefficients, of low doubles,
 *             the leading coefficient is the Jacobian matrix. */

#endif
