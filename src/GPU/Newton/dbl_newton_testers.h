// The file dbl_newton_testers.h specifies test function for Newton's method
// on series in double precision.

#ifndef __dbl_newton_testers_h__
#define __dbl_newton_testers_h__

void prompt_newton_setup
 ( int *seed, int *szt, int*nbt, int *dim, int *deg, int *size, int *posvals,
   int *vrblvl, int *mode, int *nbritr, int *nbsteps );
/*
 * DESCRIPTION :
 *   Prompts for the parameters to test Newton's method.
 *
 * ON RETURN :
 *   seed      the seed for the random number generator (0 for time);
 *   szt       size of one tile;
 *   nbt       number of tiles, szt*nbt equals the dimension;
 *   dim       the dimension is the number of monomials
 *             and the maximum number of variables in each monomial;
 *   deg       degree of the series;
 *   size      size of the numbers;
 *   posvals   positive exponents (1 if yes);
 *   vrblvl    verbose level (0 if silent);
 *   mode      execution mode, 0, 1, or 2;
 *   nbritr    number of unimodular multiplications in the making of
 *             the exponent matrix;
 *   nbsteps   the number of Newton steps. */

void dbl_unit_series_vector ( int dim, int deg, double **cff );
/*
 * DESCRIPTION :
 *   Given space in cff for the vector of dim power series,
 *   of series truncated after degree deg,
 *   initializes the coefficients in cff to one as leading coefficients,
 *   and zero for all other coefficients. */

void dbl_evaluate_monomials
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

void dbl_update_series
 ( int dim, int degp1, double **x, double **dx, int vrblvl );
/*
 * DESCRIPTION :
 *   Adds the series in dx to x.
 *
 * ON ENTRY :
 *   dim       number of series in x and dx;
 *   degp1     degree plus one of the series;
 *   x         series to updated, not linearized;
 *   dx        linearized update of the series;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   x         the series x updated with dx. */

void dbl_newton_lustep
 ( int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cff, double *acc, double **input, double ***output,
   double **funval, double ***jacval, double **rhs, double **sol,
   double **workmat, double *workvec, double **workrhs, double **resvec,
   double *resmax, int *ipvt, int vrblvl );
/*
 * DESCRIPTION :
 *   Does one step with Newton's method to update a power series,
 *   using LU factorization to solve linear systems.
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
 *   output    space for the evaluated and differentiated monomials;
 *   funval    space for the evaluated power series;
 *   jacval    space for deg+1 matrices of dimension dim;
 *   rhs       space for deg+1 vectors of dimension dim;
 *   sol       space for deg+1 vectors of dimension dim;
 *   wrkmat    work space allocated for a matrix of dimension dim;
 *   wrkvec    work space allocated for a vector of dimension dim;
 *   resvec    space for deg+1 vectors of dimension dim;
 *   ipvt      space allocated for dim pivots;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   funval    collects the output[i][dim], the evaluated series;
 *   jacval    a series with matrices as coefficients,
 *             the leading coefficient is the Jacobian matrix.
 *   rhs       the linearized right hand side are the function values
 *             subtracted by 1 and added by t;
 *   wrkmat    has the LU factorization of the Jacobian matrix;
 *   resvec    residual vectors;
 *   resmax    the maximum element of the residual vectors;
 *   ipvt      pivots used on the LU factorization of the lead matrix. */

void dbl_newton_qrstep
 ( int szt, int nbt, int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cff, double *acc,
   double **input_h, double **input_d, double ***output,
   double **funval, double ***jacval, double **rhs,
   double **urhs_h, double **urhs_d, double **sol_h, double **sol_d,
   double **Q_h, double **Q_d, double **R_h, double **R_d,
   double **workmat, double *workvec, double **resvec, double *resmax,
   int vrblvl, int mode );
/*
 * DESCRIPTION :
 *   Does one step with Newton's method to update a power series,
 *   using QR factorization to solve linear systems.
 *
 * REQUIRED : szt*nbt = dim for GPU computing.
 *
 * ON ENTRY :
 *   szt       size of each tile and block;
 *   nbt       number of tiles and number of blocks;
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
 *   input_h   coefficients of the power series of degree deg,
 *             for dim variables, computed on host;
 *   input_d   space for power series computed on device;
 *   output    space for the evaluated and differentiated monomials;
 *   funval    space for the evaluated power series;
 *   jacval    space for deg+1 matrices of dimension dim;
 *   rhs       space for deg+1 vectors of dimension dim;
 *   urhs_h    space for updated right hand side vectors computed by host;
 *   urhs_d    space for updated right hand side vectors computed by device; 
 *   sol_h     space for deg+1 vectors of dimension dim;
 *   sol_d     space for deg+1 vectors of dimension dim;
 *   Q_h       space allocated for the Q computed by the host;
 *   Q_d       space allocated for the Q computed by the device;
 *   R_h       space allocated for the R computed by the host;
 *   R_d       space allocated for the R computed by the device;
 *   wrkmat    work space allocated for a matrix of dimension dim;
 *   wrkvec    work space allocated for a vector of dimension dim;
 *   resvec    space for deg+1 vectors of dimension dim;
 *   vrblvl    is the verbose level;
 *   mode      execution mode, 0 (GPU only), 1 (CPU only) or 2 (GPU+CPU).
 *
 * ON RETURN :
 *   input_h   power series computed on host (depending on mode);
 *   input_d   power series computed on device (depending on mode);
 *   funval    collects the output[i][dim], the evaluated series;
 *   jacval    a series with matrices as coefficients,
 *             the leading coefficient is the Jacobian matrix.
 *   rhs       the linearized right hand side are the function values
 *             subtracted by 1 and added by t;
 *   urhs_h    right hand side vector updated by the host;
 *   urhs_d    right hand side vector updated by the device;
 *   sol_h     solution computed by the host;
 *   sol_d     solution computed by the device;
 *   Q_h       Q of the QR factorization computed by the host;
 *   Q_d       Q of the QR factorization computed by the device;
 *   R_h       R of the QR factorization computed by the host;
 *   R_d       R of the QR factorization computed by the device;
 *   wrkmat    has a copy of the Jacobian matrix;
 *   resvec    residual vectors;
 *   resmax    the maximum element of the residual vectors. */

#endif
