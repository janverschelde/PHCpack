// The file dbl8_newton_testers.h specifies test function for Newton's method
// on series in octo double precision.

#ifndef __dbl8_newton_testers_h__
#define __dbl8_newton_testers_h__

void dbl8_unit_series_vector
 ( int dim, int deg,
   double **cffhihihi, double **cfflohihi,
   double **cffhilohi, double **cfflolohi,
   double **cffhihilo, double **cfflohilo,
   double **cffhilolo, double **cfflololo );
/*
 * DESCRIPTION :
 *   Given space in cff* for the doubles
 *   for the vector of dim power series in octo double precision,
 *   of series truncated after degree deg,
 *   initializes the coefficients in cff to one as leading coefficients,
 *   and zero for all other coefficients. */

void dbl8_evaluate_monomials
 ( int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffhihihi, double **cfflohihi,
   double **cffhilohi, double **cfflolohi,
   double **cffhihilo, double **cfflohilo,
   double **cffhilolo, double **cfflololo,
   double *acchihihi, double *acclohihi, double *acchilohi, double *acclolohi,
   double *acchihilo, double *acclohilo, double *acchilolo, double *acclololo,
   double **inputhihihi, double **inputlohihi,
   double **inputhilohi, double **inputlolohi,
   double **inputhihilo, double **inputlohilo,
   double **inputhilolo, double **inputlololo,
   double ***outputhihihi, double ***outputlohihi,
   double ***outputhilohi, double ***outputlolohi,
   double ***outputhihilo, double ***outputlohilo,
   double ***outputhilolo, double ***outputlololo, int vrblvl );
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
 *   cffhihihi are the highest doubles of the monomial coefficients;
 *   cfflohihi are the 2nd highest doubles of the monomial coefficients;
 *   cffhilohi are the 3rd highest doubles of the monomial coefficients;
 *   cfflolohi are the 4th highest doubles of the monomial coefficients;
 *   cffhihilo are the 4th lowest doubles of the monomial coefficients;
 *   cfflohilo are the 3rd lowest doubles of the monomial coefficients;
 *   cffhilolo are the 2nd lowest doubles of the monomial coefficients;
 *   cfflololo are the lowest doubles of the monomial coefficients;
 *   acchihihi has space to accumulate one power series of degree deg;
 *   acclohihi has space to accumulate one power series of degree deg;
 *   acchilohi has space to accumulate one power series of degree deg;
 *   acclolohi has space to accumulate one power series of degree deg;
 *   acchihilo has space to accumulate one power series of degree deg;
 *   acclohilo has space to accumulate one power series of degree deg;
 *   acchilolo has space to accumulate one power series of degree deg;
 *   acclololo has space to accumulate one power series of degree deg;
 *   inputhihihi are the highest doubles of the coefficients
 *             of the power series of degree deg, for dim variables;
 *   inputlohihi are the second highest doubles of the coefficients
 *             of the power series of degree deg, for dim variables;
 *   inputhilohi are the third highest doubles of the coefficients
 *             of the power series of degree deg, for dim variables;
 *   inputlolohi are the fourth highest doubles of the coefficients
 *             of the power series of degree deg, for dim variables;
 *   inputhihilo are the fourth lowest doubles of the coefficients
 *             of the power series of degree deg, for dim variables;
 *   inputlohilo are the third lowest doubles of the coefficients
 *             of the power series of degree deg, for dim variables;
 *   inputhilolo are the second lowest doubles of the coefficients
 *             of the power series of degree deg, for dim variables;
 *   inputlololo are the lowest doubles of the coefficients
 *             of the power series of degree deg, for dim variables;
 *   outputhihihi has space for the highest output doubles;
 *   outputlohihi has space for the second highest output doubles;
 *   outputhilohi has space for the third highest output doubles;
 *   outputlolohi has space for the fourth highest output doubles;
 *   outputhihilo has space for the fourth lowest output doubles;
 *   outputlohilo has space for the third lowest output doubles;
 *   outputhilolo has space for the second lowest output doubles;
 *   outputlololo has space for the lowest output doubles;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   cffhihihi are the highest doubles of the evaluated factors;
 *   cfflohihi are the second highest doubles of the evaluated factors;
 *   cffhilohi are the third highest doubles of the evaluated factors;
 *   cfflolohi are the fourth highest doubles of the evaluated factors;
 *   cffhihilo are the fourth lowest doubles of the evaluated factors;
 *   cfflohilo are the third lowest doubles of the evaluated factors;
 *   cffhilolo are the second lowest doubles of the evaluated factors;
 *   cfflololo are the lowest doubles of the evaluated factors;
 *   outputhihihi are the highest doubles of the evaluated
 *             and differentiated monomials,
 *   outputlohihi are the second highest doubles of the evaluated
 *             and differentiated monomials,
 *   outputhilohi are the third highest doubles of the evaluated
 *             and differentiated monomials,
 *   outputlolohi are the fourth highest doubles of the evaluated
 *             and differentiated monomials,
 *   outputhihilo are the fourth lowest doubles of the evaluated
 *             and differentiated monomials,
 *   outputlohilo are the third lowest doubles of the evaluated
 *             and differentiated monomials,
 *   outputhilolo are the second lowest doubles of the evaluated
 *             and differentiated monomials,
 *   outputlololo are the lowest doubles of the evaluated
 *             and differentiated monomials,
 *             output[i][dim] is the value of the input series
 *             at the i-th monomial, and for k in range 0..nvr[i]-1,
 *             output[i][idx[i]] is the derivative w.r.t. idx[k]. */

void dbl8_linearize_evaldiff_output
 ( int dim, int degp1, int *nvr, int **idx,
   double ***outputhihihi, double ***outputlohihi,
   double ***outputhilohi, double ***outputlolohi,
   double ***outputhihilo, double ***outputlohilo,
   double ***outputhilolo, double ***outputlololo,
   double **funvalhihihi, double **funvallohihi, 
   double **funvalhilohi, double **funvallolohi,
   double **funvalhihilo, double **funvallohilo, 
   double **funvalhilolo, double **funvallololo, 
   double **rhshihihi, double **rhslohihi,
   double **rhshilohi, double **rhslolohi,
   double **rhshihilo, double **rhslohilo,
   double **rhshilolo, double **rhslololo,
   double ***jacvalhihihi, double ***jacvallohihi,
   double ***jacvalhilohi, double ***jacvallolohi,
   double ***jacvalhihilo, double ***jacvallohilo,
   double ***jacvalhilolo, double ***jacvallololo, int vrblvl );
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
 *   outputhihihi are the highest doubles of the output
 *             of the evaluation and differentiation
 *             of the monomials in the system, for the i-th monomial:
 *             output[i][dim] is the power series value, 
 *             output[i][idx[k]] is the derivative w.r.t. idx[k],
 *             for k in range 0..nvr[i]-1;
 *   outputlohihi are the second highest doubles
 *             of the output of the evaluation and differentiation;
 *   outputhilohi are the third highest doubles
 *             of the output of the evaluation and differentiation;
 *   outputlolohi are the fourth highest doubles
 *             of the output of the evaluation and differentiation;
 *   outputhihilo are the fourth lowest doubles
 *             of the output of the evaluation and differentiation;
 *   outputlohilo are the third lowest doubles
 *             of the output of the evaluation and differentiation;
 *   outputhilolo are the second lowest doubles
 *             of the output of the evaluation and differentiation;
 *   outputlololo are the lowest doubles
 *             of the output of the evaluation and differentiation;
 *   funvalhihihi has space allocated for dim power series;
 *   funvallohihi has space allocated for dim power series;
 *   funvalhilohi has space allocated for dim power series;
 *   funvallolohi has space allocated for dim power series;
 *   funvalhihilo has space allocated for dim power series;
 *   funvallohilo has space allocated for dim power series;
 *   funvalhilolo has space allocated for dim power series;
 *   funvallololo has space allocated for dim power series;
 *   rhshihihi has space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   rhslohihi has space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   rhslohihi has space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   rhslolohih has space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   rhshihilo has space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   rhslohilo has space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   rhslohilo has space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   rhslololo has space allocated for linearized power series,
 *             where degp1 is the leading dimension;
 *   jacvalhihihi has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   jacvallohihi has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   jacvalhilohi has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   jacvallolohi has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   jacvalhihilo has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   jacvallohilo has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   jacvalhilolo has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   jacvallololo has space allocated for the series of degp1 matrices,
 *             all matrices have dimension dim;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   funvalhihihi has the highest doubles of output[i][dim];
 *   funvallohihi has the second highest doubles of output[i][dim];
 *   funvallohihi has the third highest doubles of output[i][dim];
 *   funvallolohi has the fourth highest doubles of output[i][dim];
 *   funvalhihilo has the fourth lowest doubles of output[i][dim];
 *   funvallohilo has the third lowest doubles of output[i][dim];
 *   funvallohilo has the second lowest doubles of output[i][dim];
 *   funvallololo has the lowest doubles of output[i][dim];
 *   rhshihihi are the linearized right hand side are the highest doubles
 *             of the function values subtracted by 1 and added by t;
 *   rhslohihi are the linearized right hand side are the 2nd highest doubles
 *             of the function values subtracted by 1 and added by t;
 *   rhshilohi are the linearized right hand side are the 3rd highest doubles
 *             of the function values subtracted by 1 and added by t;
 *   rhslolohi are the linearized right hand side are the 4th highest doubles
 *             of the function values subtracted by 1 and added by t;
 *   rhshihilo are the linearized right hand side are the 4th lowest doubles
 *             of the function values subtracted by 1 and added by t;
 *   rhslohilo are the linearized right hand side are the 3th lowest doubles
 *             of the function values subtracted by 1 and added by t;
 *   rhshilolo are the linearized right hand side are the 2nd lowest doubles
 *             of the function values subtracted by 1 and added by t;
 *   rhslololo are the linearized right hand side are the lowest doubles
 *             of the function values subtracted by 1 and added by t;
 *   jacvalhihihi is a series with matrices as coefficients, of highest
 *             doubles, the leading coefficient is the Jacobian matrix;
 *   jacvallohihi is a series with matrices as coefficients, of second highest
 *             doubles, the leading coefficient is the Jacobian matrix;
 *   jacvalhilohi is a series with matrices as coefficients, of third highest
 *             doubles, the leading coefficient is the Jacobian matrix
 *   jacvallolohi is a series with matrices as coefficients, of fourth highest
 *             doubles, the leading coefficient is the Jacobian matrix;
 *   jacvalhihilo is a series with matrices as coefficients, of fourth lowest
 *             doubles, the leading coefficient is the Jacobian matrix;
 *   jacvallohilo is a series with matrices as coefficients, of third lowest
 *             doubles, the leading coefficient is the Jacobian matrix;
 *   jacvalhilolo is a series with matrices as coefficients, of second lowest
 *             doubles, the leading coefficient is the Jacobian matrix
 *   jacvallololo is a series with matrices as coefficients, of lowest
 *             doubles, the leading coefficient is the Jacobian matrix. */

void dbl8_update_series
 ( int dim, int degp1,
   double **xhihihi, double **xlohihi, double **xhilohi, double **xlolohi,
   double **xhihilo, double **xlohilo, double **xhilolo, double **xlololo,
   double **dxhihihi, double **dxlohihi, double **dxhilohi, double **dxlolohi,
   double **dxhihilo, double **dxlohilo, double **dxhilolo, double **dxlololo,
   int vrblvl );
/*
 * DESCRIPTION :
 *   Adds the series in dx to x.
 *
 * ON ENTRY :
 *   dim       number of series in x and dx;
 *   degp1     degree plus one of the series;
 *   xhihihi   highest doubles of the series to updated, not linearized;
 *   xlohihi   2nd highest doubles of the series to updated, not linearized;
 *   xlohihi   3rd highest doubles of the series to updated, not linearized;
 *   xlolohi   4th highest doubles of the series to updated, not linearized;
 *   xhihilo   4th lowest doubles of the series to updated, not linearized;
 *   xlohilo   3rd lowest doubles of the series to updated, not linearized;
 *   xlohilo   2nd lowest doubles of the series to updated, not linearized;
 *   xlololo   lowest doubles of the series to updated, not linearized;
 *   dxhihihi  linearized highest doubles of the update of the series;
 *   dxlohihi  linearized 2nd highest doubles of the update of the series;
 *   dxhilohi  linearized 3rd highest doubles of the update of the series;
 *   dxlolohi  linearized 4th highest doubles of the update of the series;
 *   dxhihilo  linearized 4th lowest doubles of the update of the series;
 *   dxlohilo  linearized 3rd lowest doubles of the update of the series;
 *   dxhilolo  linearized 2nd lowest doubles of the update of the series;
 *   dxlololo  linearized lowest doubles of the update of the series;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   xhihihi   highest doubles of the series x updated with dx;
 *   xlohihi   second highest doubles of the series x updated with dx;
 *   xhilohi   third highest doubles of the series x updated with dx;
 *   xlolohi   fourth highest doubles of the series x updated with dx;
 *   xhihilo   fourth lowest doubles of the series x updated with dx;
 *   xlohilo   third lowest doubles of the series x updated with dx;
 *   xhilolo   second lowest doubles of the series x updated with dx;
 *   xlololo   lowest doubles of the series x updated with dx. */

void dbl8_newton_step
 ( int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffhihihi, double **cfflohihi,
   double **cffhilohi, double **cfflolohi,
   double **cffhihilo, double **cfflohilo,
   double **cffhilolo, double **cfflololo,
   double *acchihihi, double *acclohihi,
   double *acchilohi, double *acclolohi,
   double *acchihilo, double *acclohilo,
   double *acchilolo, double *acclololo,
   double **inputhihihi, double **inputlohihi,
   double **inputhilohi, double **inputlolohi,
   double **inputhihilo, double **inputlohilo,
   double **inputhilolo, double **inputlololo,
   double ***outputhihihi, double ***outputlohihi,
   double ***outputhilohi, double ***outputlolohi,
   double ***outputhihilo, double ***outputlohilo,
   double ***outputhilolo, double ***outputlololo,
   double **funvalhihihi, double **funvallohihi,
   double **funvalhilohi, double **funvallolohi,
   double **funvalhihilo, double **funvallohilo,
   double **funvalhilolo, double **funvallololo,
   double ***jacvalhihihi, double ***jacvallohihi,
   double ***jacvalhilohi, double ***jacvallolohi,
   double ***jacvalhihilo, double ***jacvallohilo,
   double ***jacvalhilolo, double ***jacvallololo,
   double **rhshihihi, double **rhslohihi,
   double **rhshilohi, double **rhslolohi,
   double **rhshihilo, double **rhslohilo,
   double **rhshilolo, double **rhslololo,
   double **solhihihi, double **sollohihi,
   double **solhilohi, double **sollolohi,
   double **solhihilo, double **sollohilo,
   double **solhilolo, double **sollololo,
   double **workmathihihi, double **workmatlohihi,
   double **workmathilohi, double **workmatlolohi,
   double **workmathihilo, double **workmatlohilo,
   double **workmathilolo, double **workmatlololo,
   double *workvechihihi, double *workveclohihi,
   double *workvechilohi, double *workveclolohi,
   double *workvechihilo, double *workveclohilo,
   double *workvechilolo, double *workveclololo,
   double **workrhshihihi, double **workrhslohihi,
   double **workrhshilohi, double **workrhslolohi,
   double **workrhshihilo, double **workrhslohilo,
   double **workrhshilolo, double **workrhslololo,
   double **resvechihihi, double **resveclohihi,
   double **resvechilohi, double **resveclolohi,
   double **resvechihilo, double **resveclohilo,
   double **resvechilolo, double **resveclololo,
   double *resmaxhihihi, double *resmaxlohihi,
   double *resmaxhilohi, double *resmaxlolohi,
   double *resmaxhihilo, double *resmaxlohilo,
   double *resmaxhilolo, double *resmaxlololo, int *ipvt, int vrblvl );
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
 *   cffhihihi has the highest doubles of the monomial coefficients;
 *   cfflohihi has the 2nd highest doubles of the monomial coefficients;
 *   cffhilohi has the 3rd highest doubles of the monomial coefficients;
 *   cfflolohi has the 4th highest doubles of the monomial coefficients;
 *   cffhihilo has the 4th lowest doubles of the monomial coefficients;
 *   cfflohilo has the 3rd lowest doubles of the monomial coefficients;
 *   cffhilolo has the 2nd lowest doubles of the monomial coefficients;
 *   cfflololo has the lowest doubles of the monomial coefficients;
 *   acchihihi has space to accumulate one power series of degree deg;
 *   acclohihi has space to accumulate one power series of degree deg;
 *   acchilohi has space to accumulate one power series of degree deg;
 *   acclolohi has space to accumulate one power series of degree deg;
 *   acchihilo has space to accumulate one power series of degree deg;
 *   acclohilo has space to accumulate one power series of degree deg;
 *   acchilolo has space to accumulate one power series of degree deg;
 *   acclololo has space to accumulate one power series of degree deg;
 *   inputhihihi are the highest doubles of the coefficients of the power
 *             series of degree deg, for dim variables;
 *   inputlohihi are the second highest doubles of the coefficients of the
 *             power series of degree deg, for dim variables;
 *   inputhilohi are the third highest doubles of the coefficients of the
 *             power series of degree deg, for dim variables;
 *   inputlolohi are the fourth highest doubles of the coefficients of the
 *             power series of degree deg, for dim variables;
 *   inputhihilo are the fourth lowest doubles of the coefficients of the
 *             power series of degree deg, for dim variables;
 *   inputlohilo are the third lowest doubles of the coefficients of the
 *             power series of degree deg, for dim variables;
 *   inputhilolo are the second lowest doubles of the coefficients of the
 *             power series of degree deg, for dim variables;
 *   inputlololo are the lowest doubles of the coefficients of the power
 *             series of degree deg, for dim variables;
 *   outputhihihi has space for the evaluated and differentiated monomials;
 *   outputlohihi has space for the evaluated and differentiated monomials;
 *   outputhilohi has space for the evaluated and differentiated monomials;
 *   outputlolohi has space for the evaluated and differentiated monomials;
 *   outputhihilo has space for the evaluated and differentiated monomials;
 *   outputlohilo has space for the evaluated and differentiated monomials;
 *   outputhilolo has space for the evaluated and differentiated monomials;
 *   outputlololo has space for the evaluated and differentiated monomials;
 *   funvalhihihi has space for the evaluated power series;
 *   funvallohihi has space for the evaluated power series;
 *   funvalhilohi has space for the evaluated power series;
 *   funvallolohi has space for the evaluated power series;
 *   funvalhihilo has space for the evaluated power series;
 *   funvallohilo has space for the evaluated power series;
 *   funvalhilolo has space for the evaluated power series;
 *   funvallololo has space for the evaluated power series;
 *   jacvalhihihi has space for deg+1 matrices of dimension dim;
 *   jacvallohihi has space for deg+1 matrices of dimension dim;
 *   jacvalhilohi has space for deg+1 matrices of dimension dim;
 *   jacvallolohi has space for deg+1 matrices of dimension dim;
 *   jacvalhihilo has space for deg+1 matrices of dimension dim;
 *   jacvallohilo has space for deg+1 matrices of dimension dim;
 *   jacvalhilolo has space for deg+1 matrices of dimension dim;
 *   jacvallololo has space for deg+1 matrices of dimension dim;
 *   rhshihihi has space for deg+1 vectors of dimension dim;
 *   rhslohihi has space for deg+1 vectors of dimension dim;
 *   rhshilohi has space for deg+1 vectors of dimension dim;
 *   rhslolohi has space for deg+1 vectors of dimension dim;
 *   rhshihilo has space for deg+1 vectors of dimension dim;
 *   rhslohilo has space for deg+1 vectors of dimension dim;
 *   rhshilolo has pace for deg+1 vectors of dimension dim;
 *   rhslololo has space for deg+1 vectors of dimension dim;
 *   solhihihi has space for deg+1 vectors of dimension dim;
 *   sollohihi has space for deg+1 vectors of dimension dim;
 *   solhilohi has space for deg+1 vectors of dimension dim;
 *   sollolohi has space for deg+1 vectors of dimension dim;
 *   solhihilo has space for deg+1 vectors of dimension dim;
 *   sollohilo has space for deg+1 vectors of dimension dim;
 *   solhilolo has space for deg+1 vectors of dimension dim;
 *   sollololo has space for deg+1 vectors of dimension dim;
 *   wrkmathihihi has work space allocated for a matrix of dimension dim;
 *   wrkmatlohihi has work space allocated for a matrix of dimension dim;
 *   wrkmathilohi has work space allocated for a matrix of dimension dim;
 *   wrkmatlolohi has work space allocated for a matrix of dimension dim;
 *   wrkmathihilo has work space allocated for a matrix of dimension dim;
 *   wrkmatlohilo has work space allocated for a matrix of dimension dim;
 *   wrkmathilolo has work space allocated for a matrix of dimension dim;
 *   wrkmatlololo has work space allocated for a matrix of dimension dim;
 *   wrkvechihihi has work space allocated for a vector of dimension dim;
 *   wrkveclohihi has work space allocated for a vector of dimension dim;
 *   wrkvechilohi has work space allocated for a vector of dimension dim;
 *   wrkveclolohi has work space allocated for a vector of dimension dim;
 *   wrkvechihilo has work space allocated for a vector of dimension dim;
 *   wrkveclohilo has work space allocated for a vector of dimension dim;
 *   wrkvechilolo has work space allocated for a vector of dimension dim;
 *   wrkveclololo has work space allocated for a vector of dimension dim;
 *   resvechihihi has space for deg+1 vectors of dimension dim;
 *   resveclohihi has space for deg+1 vectors of dimension dim;
 *   resvechilohi has space for deg+1 vectors of dimension dim;
 *   resveclolohi has space for deg+1 vectors of dimension dim;
 *   resvechihilo has space for deg+1 vectors of dimension dim;
 *   resveclohilo has space for deg+1 vectors of dimension dim;
 *   resvechilolo has space for deg+1 vectors of dimension dim;
 *   resveclololo has space for deg+1 vectors of dimension dim;
 *   ipvt      space allocated for dim pivots;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   funvalhihihi are the highest doubles of the output[i][dim];
 *   funvallohihi are the second highest doubles of the output[i][dim];
 *   funvalhilohi are the third highest doubles of the output[i][dim];
 *   funvallolohi are the fourth highest doubles of the output[i][dim];
 *   funvalhihilo are the fourth lowest doubles of the output[i][dim];
 *   funvallohilo are the third lowest doubles of the output[i][dim];
 *   funvalhilolo are the second lowest doubles of the output[i][dim];
 *   funvallololo are the lowest doubles of the output[i][dim];
 *   jacvalhihihi are the highest doubles of a matrix series,
 *             the leading coefficient is the Jacobian matrix.
 *   jacvallohihi are the second highest doubles of a matrix series,
 *             the leading coefficient is the Jacobian matrix.
 *   jacvalhilohi are the third highest doubles of a matrix series,
 *             the leading coefficient is the Jacobian matrix.
 *   jacvallolohi are the fourth highest doubles of a matrix series,
 *             the leading coefficient is the Jacobian matrix.
 *   jacvalhihilo are the fourth lowest doubles of a matrix series,
 *             the leading coefficient is the Jacobian matrix.
 *   jacvallohilo are the third lowest doubles of a matrix series,
 *             the leading coefficient is the Jacobian matrix.
 *   jacvalhilolo are the second lowest doubles of a matrix series,
 *             the leading coefficient is the Jacobian matrix.
 *   jacvallololo are the lowest doubles of a matrix series,
 *             the leading coefficient is the Jacobian matrix.
 *   rhshihihi are the highest doubles of the linearized right hand side,
 *             the function values subtracted by 1 and added by t;
 *   rhslohihi are the second highest doubles of the linearized right hand
 *             side, the function values subtracted by 1 and added by t;
 *   rhshilohi are the third highest doubles of the linearized right hand
 *             side, are the function values subtracted by 1 and added by t;
 *   rhslolohi are the fourth highest doubles ofthe linearized right hand
 *             side, are the function values subtracted by 1 and added by t;
 *   rhshihilo are the fourth lowest doubles of the linearized right hand
 *             side, the function values subtracted by 1 and added by t;
 *   rhslohilo are the third lowest doubles of the linearized right hand
 *             side, the function values subtracted by 1 and added by t;
 *   rhshilolo are the second lowest doubles of the linearized right hand
 *             side, the function values subtracted by 1 and added by t;
 *   rhslololo are the lowest doubles ofthe linearized right hand
 *             side, the function values subtracted by 1 and added by t;
 *   wrkmathihihi are the highest doubles of the LU factorization
 *             of the Jacobian matrix;
 *   wrkmatlohihi are the second highest doubles of the LU factorization
 *             of the Jacobian matrix;
 *   wrkmathilohi are the third highest doubles of the LU factorization
 *             of the Jacobian matrix;
 *   wrkmatlolohi are the fourth highest doubles of the LU factorization
 *             of the Jacobian matrix;
 *   wrkmathihilo are the fourth lowest doubles of the LU factorization
 *             of the Jacobian matrix;
 *   wrkmatlohilo are the third lowest doubles of the LU factorization
 *             of the Jacobian matrix;
 *   wrkmathilolo are the second lowest doubles of the LU factorization
 *             of the Jacobian matrix;
 *   wrkmatlololo are the lowest doubles of the LU factorization
 *             of the Jacobian matrix;
 *   resvechihihi are the highest doubles of the residual vectors;
 *   resveclohihi are the second highest doubles of the residual vectors;
 *   resvechilohi are the third highest doubles of the residual vectors;
 *   resveclolohi are the fourth highest doubles of the residual vectors;
 *   resvechihilo are the fourth lowest doubles of the residual vectors;
 *   resveclohilo are the third lowest doubles of the residual vectors;
 *   resvechilolo are the second lowest doubles of the residual vectors;
 *   resveclololo are the lowest doubles of the residual vectors;
 *   resmaxhihihi is the highest double of the maximum element
 *             of the residual vectors;
 *   resmaxlohihi is the second highest double of the maximum element
 *             of the residual vectors;
 *   resmaxhilohi is the third highest double of the maximum element
 *             of the residual vectors;
 *   resmaxlolohi is the fourth highest double of the maximum element
 *             of the residual vectors;
 *   resmaxhihilo is the fourth lowest double of the maximum element
 *             of the residual vectors;
 *   resmaxlohilo is the third lowest double of the maximum element
 *             of the residual vectors;
 *   resmaxhilolo is the second lowest double of the maximum element
 *             of the residual vectors;
 *   resmaxlololo is the lowest double of the maximum element
 *             of the residual vectors;
 *   ipvt      pivots used on the LU factorization of the lead matrix. */

#endif
