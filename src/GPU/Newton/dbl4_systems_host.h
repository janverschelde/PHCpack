// The file dbl4_systems_host.h specifies functions to evaluate and
// differentiate monomials with common factors in quad double precision.

#ifndef __dbl4_systems_host_h__
#define __dbl4_systems_host_h__

void CPU_dbl4_evaluate_monomials
 ( int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo,
   double *acchihi, double *acclohi, double *acchilo, double *acclolo,
   double **inputhihi, double **inputlohi,
   double **inputhilo, double **inputlolo,
   double ***outputhihi, double ***outputlohi,
   double ***outputhilo, double ***outputlolo, int vrblvl );
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
 *   cffhihi   highest doubles of the coefficients of the monomials;
 *   cfflohi   2nd highest doubles of the coefficients of the monomials;
 *   cffhilo   2nd lowest doubles of the coefficients of the monomials;
 *   cfflolo   lowest doubles of the coefficients of the monomials;
 *   acchihi   space to accumulate one power series of degree deg;
 *   acclohi   space to accumulate one power series of degree deg;
 *   acchilo   space to accumulate one power series of degree deg;
 *   acclolo   space to accumulate one power series of degree deg;
 *   inputhihi are the highest doubles of the coefficients
 *             of the power series of degree deg, for dim variables;
 *   inputlohi are the second highest doubles of the coefficients
 *             of the power series of degree deg, for dim variables;
 *   inputhilo are the second lowest doubles of the coefficients
 *             of the power series of degree deg, for dim variables;
 *   inputlolo are the lowest doubles of the coefficients
 *             of the power series of degree deg, for dim variables;
 *   outputhihi has space for the highest doubles of the output;
 *   outputlohi has space for the second highest doubles of the output;
 *   outputhilo has space for the second lowest doubles of the output;
 *   outputlolo has space for the lowest doubles of the output;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   cffhihi   highest doubles of the evaluated common factors;
 *   cfflohi   second highest doubles of the evaluated common factors;
 *   cffhilo   second lowest doubles of the evaluated common factors;
 *   cfflolo   lowest doubles of the evaluated common factors;
 *   outputhihi are the highest doubles of the evaluated
 *             and differentiated monomials,
 *   outputlohi are the second highest doubles of the evaluated
 *             and differentiated monomials,
 *   outputhilo are the second lowest doubles of the evaluated
 *             and differentiated monomials,
 *   outputlolo are the lowest doubles of the evaluated
 *             and differentiated monomials,
 *             output[i][dim] is the value of the input series
 *             at the i-th monomial, and for k in range 0..nvr[i]-1,
 *             output[i][idx[i]] is the derivative w.r.t. idx[k]. */

void dbl4_linearize_evaldiff_output
 ( int dim, int degp1, int *nvr, int **idx,
   double ***outputhihi, double ***outputlohi,
   double ***outputhilo, double ***outputlolo,
   double **funvalhihi, double **funvallohi, 
   double **funvalhilo, double **funvallolo, 
   double **rhshihi, double **rhslohi, double **rhshilo, double **rhslolo,
   double ***jacvalhihi, double ***jacvallohi,
   double ***jacvalhilo, double ***jacvallolo, int vrblvl );
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
 *   outputhihi are the highest doubles of the output
 *             of the evaluation and differentiation
 *             of the monomials in the system, for the i-th monomial:
 *             output[i][dim] is the power series value, 
 *             output[i][idx[k]] is the derivative w.r.t. idx[k],
 *             for k in range 0..nvr[i]-1;
 *   outputlohi are the second highest doubles
 *             of the output of the evaluation and differentiation;
 *   outputhilo are the second lowest doubles
 *             of the output of the evaluation and differentiation;
 *   outputlolo are the lowest doubles
 *             of the output of the evaluation and differentiation;
 *   funvalhihi has space allocated for dim power series;
 *   funvallohi has space allocated for dim power series;
 *   funvalhilo has space allocated for dim power series;
 *   funvallolo has space allocated for dim power series;
 *   rhshihi   space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   rhslohi   space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   rhslohi   space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   rhslolo   space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   jacvalhihi has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   jacvallohi has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   jacvalhilo has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   jacvallolo has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   funvalhihi has the highest doubles of output[i][dim];
 *   funvallohi has the second highest doubles of output[i][dim];
 *   funvallohi has the second lowest doubles of output[i][dim];
 *   funvallolo has the lowest doubles of output[i][dim];
 *   rhshihi   the linearized right hand side are the highest doubles
 *             of the function values subtracted by 1 and added by t;
 *   rhslohi   the linearized right hand side are the 2nd highest doubles
 *             of the function values subtracted by 1 and added by t;
 *   rhshilo   the linearized right hand side are the 2nd lowest doubles
 *             of the function values subtracted by 1 and added by t;
 *   rhslolo   the linearized right hand side are the lowest doubles
 *             of the function values subtracted by 1 and added by t;
 *   jacvalhihi is a series with matrices as coefficients, of highest doubles,
 *             the leading coefficient is the Jacobian matrix;
 *   jacvallohi is a series with matrices as coefficients, of second highest
 *             doubles, the leading coefficient is the Jacobian matrix;
 *   jacvalhilo is a series with matrices as coefficients, of second lowest
 *             doubles, the leading coefficient is the Jacobian matrix
 *   jacvallolo is a series with matrices as coefficients, of lowest doubles,
 *             the leading coefficient is the Jacobian matrix. */

#endif
