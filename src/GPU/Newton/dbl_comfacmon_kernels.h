// The file dbl_comfacmon_kernels.h specifies functions to define
// memory transfers ad kernel launches to evaluate and differentiate
// monomials with common factors in double precision.

#ifndef __dbl_comfacmon_kernels_h__
#define __dbl_comfacmon_kernels_h__

void GPU_dbl_evaluate_monomials
 ( int dim, int deg, int szt, int nbt,
   int *nvr, int **idx, int **exp, int *nbrfac,
   int **expfac, double **cff, double *acc, double **input,
   double ***output, int vrblvl );
/*
 * DESCRIPTION :
 *   Evaluates monomials at power series.
 *
 * ON ENTRY :
 *   dim       number of monomials;
 *   deg       degree of the power series;
 *   szt       size of each block of threads;
 *   nbt       number of thread blocks;
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

#endif
