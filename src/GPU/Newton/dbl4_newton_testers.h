// The file dbl4_newton_testers.h specifies test function for Newton's method
// on series in quad double precision.

#ifndef __dbl4_newton_testers_h__
#define __dbl4_newton_testers_h__

void dbl4_unit_series_vector
 ( int dim, int deg,
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo );
/*
 * DESCRIPTION :
 *   Given space in cffhihi, cfflohi, cffhilo, cfflolo for the doubles
 *   for the vector of dim power series in quad double precision,
 *   of series truncated after degree deg,
 *   initializes the coefficients in cff to one as leading coefficients,
 *   and zero for all other coefficients. */

void dbl4_evaluate_monomials
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

void dbl4_update_series
 ( int dim, int degp1,
   double **xhihi, double **xlohi, double **xhilo, double **xlolo,
   double **dxhihi, double **dxlohi, double **dxhilo, double **dxlolo,
   int vrblvl );
/*
 * DESCRIPTION :
 *   Adds the series in dx to x.
 *
 * ON ENTRY :
 *   dim       number of series in x and dx;
 *   degp1     degree plus one of the series;
 *   xhihi     highest doubles of the series to updated, not linearized;
 *   xlohi     2nd highest doubles of the series to updated, not linearized;
 *   xlohi     2nd lowest doubles of the series to updated, not linearized;
 *   xlolo     lowest doubles of the series to updated, not linearized;
 *   dxhihi    linearized highest doubles of the update of the series;
 *   dxlohi    linearized 2nd highest doubles of the update of the series;
 *   dxhilo    linearized 2nd lowest doubles of the update of the series;
 *   dxlolo    linearized lowest doubles of the update of the series;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   xhihi     highest doubles of the series x updated with dx;
 *   xlohi     second highest doubles of the series x updated with dx;
 *   xhilo     second lowest doubles of the series x updated with dx;
 *   xlolo     lowest doubles of the series x updated with dx. */

void dbl4_newton_step
 ( int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo,
   double *acchihi, double *acclohi, double *acchilo, double *acclolo,
   double **inputhihi, double **inputlohi,
   double **inputhilo, double **inputlolo,
   double ***outputhihi, double ***outputlohi,
   double ***outputhilo, double ***outputlolo,
   double **funvalhihi, double **funvallohi,
   double **funvalhilo, double **funvallolo,
   double ***jacvalhihi, double ***jacvallohi,
   double ***jacvalhilo, double ***jacvallolo,
   double **rhshihi, double **rhslohi, double **rhshilo, double **rhslolo,
   double **solhihi, double **sollohi, double **solhilo, double **sollolo,
   double **workmathihi, double **workmatlohi,
   double **workmathilo, double **workmatlolo,
   double *workvechihi, double *workveclohi,
   double *workvechilo, double *workveclolo,
   double **workrhshihi, double **workrhslohi,
   double **workrhshilo, double **workrhslolo,
   double **resvechihi, double **resveclohi,
   double **resvechilo, double **resveclolo,
   double *resmaxhihi, double *resmaxlohi,
   double *resmaxhilo, double *resmaxlolo, int *ipvt, int vrblvl );
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
 *   cffhihi   highest doubles of the coefficients of the monomials;
 *   cfflohi   2nd highest doubles of the coefficients of the monomials;
 *   cffhilo   2nd lowest doubles of the coefficients of the monomials;
 *   cfflolo   lowest doubles of the coefficients of the monomials;
 *   acchihi   space to accumulate one power series of degree deg;
 *   acclohi   space to accumulate one power series of degree deg;
 *   acchilo   space to accumulate one power series of degree deg;
 *   acclolo   space to accumulate one power series of degree deg;
 *   inputhihi are the highest doubles of the coefficients of the power
 *             series of degree deg, for dim variables;
 *   inputlohi are the second highest doubles of the coefficients of the
 *             power series of degree deg, for dim variables;
 *   inputhilo are the second lowest doubles of the coefficients of the
 *             power series of degree deg, for dim variables;
 *   inputlolo are the lowest doubles of the coefficients of the power
 *             series of degree deg, for dim variables;
 *   outputhihi has space for the evaluated and differentiated monomials;
 *   outputlohi has space for the evaluated and differentiated monomials;
 *   outputhilo has space for the evaluated and differentiated monomials;
 *   outputlolo has space for the evaluated and differentiated monomials;
 *   funvalhihi has space for the evaluated power series;
 *   funvallohi has space for the evaluated power series;
 *   funvalhilo has space for the evaluated power series;
 *   funvallolo has space for the evaluated power series;
 *   jacvalhihi has space for deg+1 matrices of dimension dim;
 *   jacvallohi has space for deg+1 matrices of dimension dim;
 *   jacvalhilo has space for deg+1 matrices of dimension dim;
 *   jacvallolo has space for deg+1 matrices of dimension dim;
 *   rhshihi   space for deg+1 vectors of dimension dim;
 *   rhslohi   space for deg+1 vectors of dimension dim;
 *   rhshilo   space for deg+1 vectors of dimension dim;
 *   rhslolo   space for deg+1 vectors of dimension dim;
 *   solhihi   space for deg+1 vectors of dimension dim;
 *   sollohi   space for deg+1 vectors of dimension dim;
 *   solhilo   space for deg+1 vectors of dimension dim;
 *   sollolo   space for deg+1 vectors of dimension dim;
 *   wrkmathihi has work space allocated for a matrix of dimension dim;
 *   wrkmatlohi has work space allocated for a matrix of dimension dim;
 *   wrkmathilo has work space allocated for a matrix of dimension dim;
 *   wrkmatlolo has work space allocated for a matrix of dimension dim;
 *   wrkvechihi has work space allocated for a vector of dimension dim;
 *   wrkveclohi has work space allocated for a vector of dimension dim;
 *   wrkvechilo has work space allocated for a vector of dimension dim;
 *   wrkveclolo has work space allocated for a vector of dimension dim;
 *   resvechihi has space for deg+1 vectors of dimension dim;
 *   resveclohi has space for deg+1 vectors of dimension dim;
 *   resvechilo has space for deg+1 vectors of dimension dim;
 *   resveclolo has space for deg+1 vectors of dimension dim;
 *   ipvt      space allocated for dim pivots;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   funvalhihi is the highest doubles of the output[i][dim];
 *   funvallohi is the second highest doubles of the output[i][dim];
 *   funvalhilo is the second lowest doubles of the output[i][dim];
 *   funvallolo is the lowest doubles of the output[i][dim];
 *   jacvalhihi are the highest doubles of a matrix series,
 *             the leading coefficient is the Jacobian matrix.
 *   jacvallohi are the second highest doubles of a matrix series,
 *             the leading coefficient is the Jacobian matrix.
 *   jacvalhilo are the second lowest doubles of a matrix series,
 *             the leading coefficient is the Jacobian matrix.
 *   jacvallolo are the lowest doubles of a matrix series,
 *             the leading coefficient is the Jacobian matrix.
 *   rhshihi   highest doubles of the linearized right hand side
 *             are the function values subtracted by 1 and added by t;
 *   rhslohi   second highest doubles of the linearized right hand side
 *             are the function values subtracted by 1 and added by t;
 *   rhshilo   second lowest doubles of the linearized right hand side
 *             are the function values subtracted by 1 and added by t;
 *   rhslolo   lowest doubles ofthe linearized right hand side
 *             are the function values subtracted by 1 and added by t;
 *   wrkmathihi are the highest doubles of the LU factorization
 *             of the Jacobian matrix;
 *   wrkmatlohi are the second highest doubles of the LU factorization
 *             of the Jacobian matrix;
 *   wrkmathilo are the second lowest doubles of the LU factorization
 *             of the Jacobian matrix;
 *   wrkmatlolo are the lowest doubles of the LU factorization
 *             of the Jacobian matrix;
 *   resvechihi are the highest doubles of the residual vectors;
 *   resveclohi are the second highest doubles of the residual vectors;
 *   resvechilo are the second lowest doubles of the residual vectors;
 *   resveclolo are the lowest doubles of the residual vectors;
 *   resmaxhihi is the highest double of the maximum element
 *             of the residual vectors;
 *   resmaxlohi is the second highest double of the maximum element
 *             of the residual vectors;
 *   resmaxhilo is the second lowest double of the maximum element
 *             of the residual vectors;
 *   resmaxlolo is the lowest double of the maximum element
 *             of the residual vectors;
 *   ipvt      pivots used on the LU factorization of the lead matrix. */

#endif
