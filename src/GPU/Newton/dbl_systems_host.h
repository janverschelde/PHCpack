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
 *   Evaluates monomials at power series, on real data.
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

void CPU_cmplx_evaluate_monomials
 ( int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffre, double **cffim, double *accre, double *accim,
   double **inputre, double **inputim, double ***outputre, double ***outputim,
   int vrblvl );
/*
 * DESCRIPTION :
 *   Evaluates monomials at power series, on complex data.
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
 *   cffre     real parts of the coefficients of the monomials;
 *   cffim     imaginary parts of the coefficients of the monomials;
 *   accre     space to accumulate one power series of degree deg;
 *   accim     space to accumulate one power series of degree deg;
 *   inputre   real parts of coefficients of the series of degree deg,
 *             for dim variables;
 *   inputim   imaginary parts of coefficients of the series of degree deg,
 *             for dim variables;
 *   outputre  space for the output;
 *   outputim  space for the output;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   cffre     real parts of the evaluated common factors;
 *   cffim     imaginary parts of the evaluated common factors;
 *   outputre  real parts of evaluated and differentiated monomials,
 *             outputre[i][dim] is the real part of value of the input
 *             at the i-th monomial, and for k in range 0..nvr[i]-1,
 *             outputre[i][idx[i]] is the derivative w.r.t. idx[k];
 *   outputim  imaginary parts of evaluated and differentiated monomials,
 *             outputim[i][dim] is the real part of value of the input
 *             at the i-th monomial, and for k in range 0..nvr[i]-1,
 *             outputim[i][idx[i]] is the derivative w.r.t. idx[k]. */

void dbl_linearize_evaldiff_output
 ( int dim, int degp1, int *nvr, int **idx, double damper,
   double ***output, double **funval, double **rhs, double ***jacval,
   int vrblvl );
/*
 * DESCRIPTION :
 *   Linearizes the output of the evaluation and differentiation
 *   of the monomials, on real data.
 *
 * ON ENTRY :
 *   dim       number of monomials;
 *   degp1     degree plus one;
 *   nvr       number of variables that occur in each monomial;
 *   idx       for each monomials the indices of each variable;
 *   damper    the positive damping coefficient for t,
 *             if 1.0, then no damping, if > 1.0, then overdamping;
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

void cmplx_linearize_evaldiff_output
 ( int dim, int degp1, int *nvr, int **idx, double damper,
   double ***outputre, double ***outputim,
   double **funvalre, double **funvalim,
   double **rhsre, double **rhsim, double ***jacvalre, double ***jacvalim,
   int vrblvl );
/*
 * DESCRIPTION :
 *   Linearizes the output of the evaluation and differentiation
 *   of the monomials, on complex data.
 *
 * ON ENTRY :
 *   dim       number of monomials;
 *   degp1     degree plus one;
 *   nvr       number of variables that occur in each monomial;
 *   idx       for each monomials the indices of each variable;
 *   damper    the positive damping coefficient for t,
 *             if 1.0, then no damping, if > 1.0, then overdamping;
 *   outputre  real parts of the output of evaluated and differentiated
 *             system, for the i-th monomial:
 *             outputre[i][dim] is the power series value, 
 *             outputre[i][idx[k]] is the derivative w.r.t. idx[k],
 *             for k in range 0..nvr[i]-1;
 *   outputim  imaginary parts of the output of evaluated and
 *             differentiated system, for the i-th monomial:
 *             outputim[i][dim] is the power series value, 
 *             outputim[i][idx[k]] is the derivative w.r.t. idx[k],
 *             for k in range 0..nvr[i]-1;
 *   funvalre  space allocated for dim power series;
 *   funvalim  space allocated for dim power series;
 *   rhsre     space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   rhsim     space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   jacvalre  space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   jacvalim  space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   funvalre  real parts of the outputre[i][dim];
 *   funvalim  imaginary parts of the outputre[i][dim];
 *   rhsre     real parts of the linearized right hand side,
 *             subtracted by 1 and added by t;
 *   rhsim     imaginary parts of the linearized right hand side,
 *             subtracted by 1 and added by t;
 *   jacvalre  are the real parts of the matrix series;
 *   jacvalim  are the imaginary parts of the matrix series. */

#endif
