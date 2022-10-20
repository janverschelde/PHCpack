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

void CPU_cmplx2_evaluate_monomials
 ( int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffrehi, double **cffrelo, double **cffimhi, double **cffimlo,
   double *accrehi, double *accrelo, double *accimhi, double *accimlo,
   double **inputrehi, double **inputrelo,
   double **inputimhi, double **inputimlo, 
   double ***outputrehi, double ***outputrelo, 
   double ***outputimhi, double ***outputimlo, int vrblvl );
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
 *   cffrehi   high doubles of real parts of monomial coefficients;
 *   cffrelo   low doubles of real parts of monomial coefficients;
 *   cffimhi   high doubles of imaginary parts of monomial coefficients;
 *   cffimlo   low doubles of imaginary parts of monomial coefficients;
 *   accrehi   space to accumulate one power series of degree deg;
 *   accrelo   space to accumulate one power series of degree deg;
 *   accimhi   space to accumulate one power series of degree deg;
 *   accimlo   space to accumulate one power series of degree deg;
 *   inputrehi are the high doubles of the real parts of coefficients
 *             of the series of degree deg, for dim variables;
 *   inputrelo are the low doubles of the real parts of coefficients
 *             of the series of degree deg, for dim variables;
 *   inputimhi are the high doubles of the imaginary parts of coefficients
 *             of the series of degree deg, for dim variables;
 *   inputimlo are the low doubles of the imaginary parts of coefficients
 *             of the series of degree deg, for dim variables;
 *   outputrehi has space for the output;
 *   outputrelo has space for the output;
 *   outputimhi has space for the output;
 *   outputimlo has space for the output;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   cffrehi   high doubles of real parts of the evaluated common factors;
 *   cffrelo   low doubles of real parts of the evaluated common factors;
 *   cffimhi   high doubles of imaginary parts of the evaluated common factors;
 *   cffimlo   low doubles of imaginary parts of the evaluated common factors;
 *   outputrehi are the high doubles of the real parts of evaluated
 *             and differentiated monomials, outputrehi[i][dim] is the high
 *             double of the real part of value of the input
 *             at the i-th monomial, and for k in range 0..nvr[i]-1,
 *             outputrehi[i][idx[k]] is the derivative w.r.t. idx[k];
 *   outputrelo are the low doubles of the real parts of evaluated
 *             and differentiated monomials, outputrelo[i][dim] is the low 
 *             double of the real part of value of the input
 *             at the i-th monomial, and for k in range 0..nvr[i]-1,
 *             outputrelo[i][idx[k]] is the derivative w.r.t. idx[k];
 *   outputimhi are the high doubles of the imaginary parts of evaluated
 *             and differentiated monomials, outputimhi[i][dim] is the high
 *             double of the imaginary part of value of the input
 *             at the i-th monomial, and for k in range 0..nvr[i]-1,
 *             outputimhi[i][idx[k]] is the derivative w.r.t. idx[k];
 *   outputimlo are the low doubles of the imaginary parts of evaluated
 *             and differentiated monomials, outputimlo[i][dim] is the low
 *             double of the imaginary part of value of the input
 *             at the i-th monomial, and for k in range 0..nvr[i]-1,
 *             outputimlo[i][idx[k]] is the derivative w.r.t. idx[k]. */

void dbl2_linearize_evaldiff_output
 ( int dim, int degp1, int *nvr, int **idx,
   double **mbhi, double **mblo, double damper,
   double ***outputhi, double ***outputlo,
   double **funvalhi, double **funvallo, 
   double **rhshi, double **rhslo, double ***jacvalhi, double ***jacvallo,
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
 *   mbhi      hight doubles of the right hand side of monomial system;
 *   mblo      low doubles of the right hand side of monomial system;
 *   damper    the positive damping coefficient for t,
 *             if 1.0, then no damping, if > 1.0, then overdamping;
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

void cmplx2_linearize_evaldiff_output
 ( int dim, int degp1, int *nvr, int **idx,
   double **mbrehi, double **mbrelo, double **mbimhi, double **mbimlo,
   double damper,
   double ***outputrehi, double ***outputrelo,
   double ***outputimhi, double ***outputimlo,
   double **funvalrehi, double **funvalrelo,
   double **funvalimhi, double **funvalimlo,
   double **rhsrehi, double **rhsrelo, double **rhsimhi, double **rhsimlo,
   double ***jacvalrehi, double ***jacvalrelo,
   double ***jacvalimhi, double ***jacvalimlo, int vrblvl );
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
 *   mbrehi    high real parts of the right hand side of monomial system;
 *   mbrelo    low real parts of the right hand side of monomial system;
 *   mbimhi    high imaginary parts of the right hand side series;
 *   mbimlo    low imaginary parts of the right hand side series;
 *   damper    the positive damping coefficient for t,
 *             if 1.0, then no damping, if > 1.0, then overdamping;
 *   outputrehi are the high doubles of the real parts of the output
 *             of evaluated and differentiated system, for the i-th monomial:
 *             outputrehi[i][dim] is the power series value, 
 *             outputrehi[i][idx[k]] is the derivative w.r.t. idx[k],
 *             for k in range 0..nvr[i]-1;
 *   outputrelo are the low doubles of the real parts of the output
 *             of evaluated and differentiated system, for the i-th monomial:
 *             outputrelo[i][dim] is the power series value, 
 *             outputrelo[i][idx[k]] is the derivative w.r.t. idx[k],
 *             for k in range 0..nvr[i]-1;
 *   outputimhi are the high doubles of the imaginary parts of the output
 *             of evaluated and differentiated system, for the i-th monomial:
 *             outputimhi[i][dim] is the power series value, 
 *             outputimhi[i][idx[k]] is the derivative w.r.t. idx[k],
 *             for k in range 0..nvr[i]-1;
 *   outputimlo are the low doubles of the imaginary parts of the output
 *             of evaluated and differentiated system, for the i-th monomial:
 *             outputimlo[i][dim] is the power series value, 
 *             outputimlo[i][idx[k]] is the derivative w.r.t. idx[k],
 *             for k in range 0..nvr[i]-1;
 *   funvalrehi has space allocated for dim power series;
 *   funvalrelo has space allocated for dim power series;
 *   funvalimhi has space allocated for dim power series;
 *   funvalimlo has space allocated for dim power series;
 *   rhsrehi   space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   rhsrelo   space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   rhsimhi   space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   rhsimlo   space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   jacvalrehi has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   jacvalrelo has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   jacvalimhi has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   jacvalimlo has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   funvalrehi are the high doubles of the real parts in output[i][dim];
 *   funvalrelo are the low doubles of the real parts in output[i][dim];
 *   funvalimhi are the high doubles of the imaginary parts 
 *             of the output[i][dim];
 *   funvalimlo are the low doubles of the imaginary parts
 *             of the output[i][dim];
 *   rhsrehi   high doubles of real parts of the linearized right hand side,
 *             subtracted by 1 and added by t;
 *   rhsrelo   low doubles of real parts of the linearized right hand side,
 *             subtracted by 1 and added by t;
 *   rhsimhi   high doubles of imaginary parts of the linearized right hand
 *             side, subtracted by 1 and added by t;
 *   rhsimlo   low doubles of imaginary parts of the linearized right hand
 *             side, subtracted by 1 and added by t;
 *   jacvalrehi are the high doubles of the real parts of the matrix series;
 *   jacvalrelo are the low doubles of the real parts of the matrix series;
 *   jacvalimhi are the high doubles of imaginary parts of the matrix series;
 *   jacvalimlo are the low doubles imaginary parts of the matrix series. */

#endif
