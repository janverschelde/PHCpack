// The file dbl2_newton_testers.h specifies test function for Newton's method
// on series in double double precision.

#ifndef __dbl2_newton_testers_h__
#define __dbl2_newton_testers_h__

void dbl2_unit_series_vector
 ( int dim, int deg, double **cffhi, double **cfflo );
/*
 * DESCRIPTION :
 *   Given space in cffhi, cfflo for the high and low doubles
 *   for the vector of dim power series,
 *   of series truncated after degree deg,
 *   initializes the coefficients in cff to one as leading coefficients,
 *   and zero for all other coefficients. */

void dbl2_evaluate_monomials
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

void dbl2_update_series
 ( int dim, int degp1, double **xhi, double **xlo,
   double **dxhi, double **dxlo, int vrblvl );
/*
 * DESCRIPTION :
 *   Adds the series in dx to x.
 *
 * ON ENTRY :
 *   dim       number of series in x and dx;
 *   degp1     degree plus one of the series;
 *   xhi       high doubles of the series to updated, not linearized;
 *   xlo       low doubles of the series to updated, not linearized;
 *   dxhi      linearized high doubles of the update of the series;
 *   dxlo      linearized low doubles of the update of the series;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   xhi       high doubles of the series x updated with dx;
 *   xlo       low doubles of the series x updated with dx. */

void dbl2_newton_step
 ( int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffhi, double **cfflo, double *acchi, double *acclo,
   double **inputhi, double **inputlo,
   double ***outputhi, double ***outputlo,
   double **funvalhi, double **funvallo,
   double ***jacvalhi, double ***jacvallo,
   double **rhshi, double **rhslo, double **solhi, double **sollo,
   double **workmathi, double **workmatlo,
   double *workvechi, double *workveclo,
   double **workrhshi, double **workrhslo,
   double **resvechi, double **resveclo,
   double *resmaxhi, double *resmaxlo, int *ipvt, int vrblvl );
/*
 * DESCRIPTION :
 *   Does one step with Newton's method to update a power series.
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
 *   input     low doubles of the coefficients of the power series
 *             of degree deg, for dim variables;
 *   outputhi  space for the evaluated and differentiated monomials;
 *   outputlo  space for the evaluated and differentiated monomials;
 *   funvalhi  space for the evaluated power series;
 *   funvallo  space for the evaluated power series;
 *   jacvalhi  space for deg+1 matrices of dimension dim;
 *   jacvallo  space for deg+1 matrices of dimension dim;
 *   rhshi     space for deg+1 vectors of dimension dim;
 *   rhslo     space for deg+1 vectors of dimension dim;
 *   solhi     space for deg+1 vectors of dimension dim;
 *   sollo     space for deg+1 vectors of dimension dim;
 *   wrkmathi  work space allocated for a matrix of dimension dim;
 *   wrkmatlo  work space allocated for a matrix of dimension dim;
 *   wrkvechi  work space allocated for a vector of dimension dim;
 *   wrkveclo  work space allocated for a vector of dimension dim;
 *   resvechi  space for deg+1 vectors of dimension dim;
 *   resveclo  space for deg+1 vectors of dimension dim;
 *   ipvt      space allocated for dim pivots;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   funvalhi  high doubles of the output[i][dim], the evaluated series;
 *   funvallo  low doubles of the output[i][dim], the evaluated series;
 *   jacvalhi  high doubles of a series with matrices as coefficients,
 *             the leading coefficient is the Jacobian matrix.
 *   jacvallo  low doubles of a series with matrices as coefficients,
 *             the leading coefficient is the Jacobian matrix.
 *   rhshi     high doubles of the linearized right hand side are the 
 *             function values subtracted by 1 and added by t;
 *   rhslo     low doubles of the linearized right hand side are the 
 *             function values subtracted by 1 and added by t;
 *   wrkmathi  high doubles of the LU factorization of the Jacobian matrix;
 *   wrkmatlo  low doubles of the LU factorization of the Jacobian matrix;
 *   resvechi  high doubles of the residual vectors;
 *   resveclo  low doubles of the residual vectors;
 *   resmaxhi  high double of the maximum element of the residual vectors;
 *   resmaxlo  low double of the maximum element of the residual vectors;
 *   ipvt      pivots used on the LU factorization of the lead matrix. */

#endif
