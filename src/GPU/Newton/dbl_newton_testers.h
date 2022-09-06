// The file dbl_newton_testers.h specifies test functions for Newton's method
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

void cmplx_unit_series_vector
 ( int dim, int deg, double **cffre, double **cffim );
/*
 * DESCRIPTION :
 *   Given space allocated in cffre and cffim,
 *   returns the values for a complex unit series. */

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

void cmplx_update_series
 ( int dim, int degp1,
   double **xre, double **xim, double **dxre, double **dxim,
   int vrblvl );
/*
 * DESCRIPTION :
 *   Adds the series in dx to x.
 *
 * ON ENTRY :
 *   dim       number of series in x and dx;
 *   degp1     degree plus one of the series;
 *   xre       real parts of series to updated, not linearized;
 *   xim       imaginary parts of series to updated, not linearized;
 *   dxre      real parts of the linearized update of the series;
 *   dxim      imaginary parts of the linearized update of the series;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   xre       real parts of the series x updated with dx;
 *   xim       imaginary parts of the series x updated with dx. */

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
   double **input_h, double **input_d, double ***output_h, double ***output_d,
   double **funval_h, double **funval_d,
   double ***jacval_h, double ***jacval_d, double **rhs_h, double **rhs_d,
   double **urhs_h, double **urhs_d, double **sol_h, double **sol_d,
   double **Q_h, double **Q_d, double **R_h, double **R_d,
   double **workmat, double *workvec, double **resvec, double *resmax,
   int vrblvl, int mode );
/*
 * DESCRIPTION :
 *   Does one step with Newton's method to update a power series,
 *   using QR factorization to solve linear systems, on real data.
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
 *   output_h  space for the evaluated and differentiated monomials,
 *             computed on the host;
 *   output_d  space for the evaluated and differentiated monomials,
 *             computed on the device;
 *   funval_h  space for the evaluated power series computed by host;
 *   funval_d  space for the evaluated power series computed by device;
 *   jacval_h  space for deg+1 matrices of dimension dim on host;
 *   jacval_d  space for deg+1 matrices of dimension dim on device;
 *   rhs_h     space for deg+1 vectors of dimension dim on host;
 *   rhs_d     space for deg+1 vectors of dimension dim on device;
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
 *   output_h  evaluated power series computed on host (depending on mode);
 *   output_d  evaluated power series computed on device (depending on mode);
 *   funval_h  collects the output[i][dim], the evaluated series on host;
 *   funval_d  collects the output[i][dim], the evaluated series on device;
 *   jacval_h  a series with matrices as coefficients, computed by host,
 *             the leading coefficient is the Jacobian matrix;
 *   jacval_d  a series with matrices as coefficients, computed by device,
 *             the leading coefficient is the Jacobian matrix;
 *   rhs_h     the linearized right hand side are the function values
 *             subtracted by 1 and added by t, computed by host;
 *   rhs_d     the linearized right hand side are the function values
 *             subtracted by 1 and added by t, computed by device;
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

void cmplx_newton_qrstep
 ( int szt, int nbt, int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffre, double **cffim, double *accre, double *accim,
   double **inputre_h, double **inputim_h,
   double **inputre_d, double **inputim_d,
   double ***outputre_h, double ***outputim_h,
   double ***outputre_d, double ***outputim_d,
   double **funvalre_h, double **funvalim_h,
   double **funvalre_d, double **funvalim_d,
   double ***jacvalre_h, double ***jacvalim_h,
   double ***jacvalre_d, double ***jacvalim_d,
   double **rhsre_h, double **rhsim_h, double **rhsre_d, double **rhsim_d,
   double **urhsre_h, double **urhsim_h, double **urhsre_d, double **urhsim_d,
   double **solre_h, double **solim_h, double **solre_d, double **solim_d,
   double **Qre_h, double **Qim_h, double **Qre_d, double **Qim_d, 
   double **Rre_h, double **Rim_h, double **Rre_d, double **Rim_d,
   double **workmatre, double **workmatim,
   double *workvecre, double *workvecim,
   double **resvecre, double **resvecim, double *resmax,
   int vrblvl, int mode );
/*
 * DESCRIPTION :
 *   Does one step with Newton's method to update a power series,
 *   using QR factorization to solve linear systems, on real data.
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
 *   cffre     real parts of the coefficients of the monomials;
 *   cffim     imaginary parts of the coefficients of the monomials;
 *   accre     space to accumulate one power series of degree deg;
 *   accim     space to accumulate one power series of degree deg;
 *   inputre_h are the real parts coefficients of the series of degree deg,
 *             for dim variables, computed on host;
 *   inputim_h are the imaginary parts coefficients of the series
 *             of degree deg for dim variables, computed on host;
 *   inputre_d has space for power series computed on device;
 *   inputim_d has space for power series computed on device;
 *   outputre_h has space for the evaluated and differentiated monomials,
 *             computed on the host;
 *   outputim_h has space for the evaluated and differentiated monomials,
 *             computed on the host;
 *   outputre_d has space for the evaluated and differentiated monomials,
 *             computed on the device;
 *   outputim_d has space for the evaluated and differentiated monomials,
 *             computed on the device;
 *   funvalre_h has space for the evaluated power series computed by host;
 *   funvalim_h has space for the evaluated power series computed by host;
 *   funvalre_d has space for the evaluated power series computed by device;
 *   funvalim_d has space for the evaluated power series computed by device;
 *   jacvalre_h has space for deg+1 matrices of dimension dim on host;
 *   jacvalim_h has space for deg+1 matrices of dimension dim on host;
 *   jacvalre_d has space for deg+1 matrices of dimension dim on device;
 *   jacvalim_d has space for deg+1 matrices of dimension dim on device;
 *   rhsre_h   space for deg+1 vectors of dimension dim on host;
 *   rhsim_h   space for deg+1 vectors of dimension dim on host;
 *   rhsre_d   space for deg+1 vectors of dimension dim on device;
 *   rhsim_d   space for deg+1 vectors of dimension dim on device;
 *   urhsre_h  space for updated right hand side vectors computed by host;
 *   urhsim_h  space for updated right hand side vectors computed by host;
 *   urhsre_d  space for updated right hand side vectors computed by device; 
 *   urhsim_d  space for updated right hand side vectors computed by device; 
 *   solre_h   space for deg+1 vectors of dimension dim;
 *   solim_h   space for deg+1 vectors of dimension dim;
 *   solre_d   space for deg+1 vectors of dimension dim;
 *   solim_d   space for deg+1 vectors of dimension dim;
 *   Qre_h     space allocated for the Q computed by the host;
 *   Qim_h     space allocated for the Q computed by the host;
 *   Qre_d     space allocated for the Q computed by the device;
 *   Qim_d     space allocated for the Q computed by the device;
 *   Rre_h     space allocated for the R computed by the host;
 *   Rim_h     space allocated for the R computed by the host;
 *   Rre_d     space allocated for the R computed by the device;
 *   Rim_d     space allocated for the R computed by the device;
 *   wrkmatre  work space allocated for a matrix of dimension dim;
 *   wrkmatim  work space allocated for a matrix of dimension dim;
 *   wrkvecre  work space allocated for a vector of dimension dim;
 *   wrkvecim  work space allocated for a vector of dimension dim;
 *   resvecre  space for deg+1 vectors of dimension dim;
 *   resvecim  space for deg+1 vectors of dimension dim;
 *   vrblvl    is the verbose level;
 *   mode      execution mode, 0 (GPU only), 1 (CPU only) or 2 (GPU+CPU).
 *
 * ON RETURN :
 *   inputre_h has real parts of series, on host (depending on mode);
 *   inputim_h has imaginary parts of series, on host (depending on mode);
 *   inputre_d has real parts of series, on device (depending on mode);
 *   inputim_d has imaginary parts of series, on device (depending on mode);
 *   outputre_h is the real part of output, on host (depending on mode);
 *   outputim_h is the imag part of output, on host (depending on mode);
 *   outputre_d is the real part of the output, on device (depending on mode);
 *   outputim_d is the imag part of the output, on device (depending on mode);
 *   funvalre_h is outputre[i][dim], on host;
 *   funvalim_h is outputim[i][dim], on host;
 *   funvalre_d is outputre[i][dim], on device;
 *   funvalim_d is outputim[i][dim], on device;
 *   jacvalre_h are the real parts of a matrix series, on host;
 *   jacvalim_h are the real parts of a matrix series, on host;
 *   jacvalre_d are the real parts of a matrix series, on device;
 *   jacvalim_d are the real parts of a matrix series, on device;
 *   rhsre_h   real parts of the linearized right hand side,
 *             subtracted by 1 and added by t, computed by host;
 *   rhsim_h   imaginary parts of the linearized right hand side,
 *             subtracted by 1 and added by t, computed by host;
 *   rhsre_d   real parts of the linearized right hand side,
 *             subtracted by 1 and added by t, computed by device;
 *   rhsim_d   imaginary parts of the linearized right hand side,
 *             subtracted by 1 and added by t, computed by device;
 *   urhsre_h  real parts of right hand side updated by the host;
 *   urhsim_h  imaginary parts of right hand side updated by the host;
 *   urhsre_d  real parts of right hand side updated by the device;
 *   urhsim_d  imaginary parts of right hand side updated by the device;
 *   solre_h   real parts of solution computed by the host;
 *   solim_h   imaginary parts of solution computed by the host;
 *   solre_d   real parts of solution computed by the device;
 *   solim_d   imagnary parts of solution computed by the device;
 *   Qre_h     real Q of the QR factorization computed by the host;
 *   Qim_h     imaginary Q of the QR factorization computed by the host;
 *   Qre_d     real Q of the QR factorization computed by the device;
 *   Qim_d     imaginary Q of the QR factorization computed by the device;
 *   Rre_h     real R of the QR factorization computed by the host;
 *   Rim_h     imaginary R of the QR factorization computed by the host;
 *   Rre_d     real R of the QR factorization computed by the device;
 *   Rim_d     imaginary R of the QR factorization computed by the device;
 *   wrkmat    has a copy of the Jacobian matrix;
 *   resvecre  real parts of the residual vectors;
 *   resvecim  imaginary parts of the residual vectors;
 *   resmax    the maximum element of the residual vectors. */

int test_dbl_real_newton
 ( int szt, int nbt, int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   int nbsteps, int mode, int vrblvl );
/*
 * DESCRIPTION :
 *   Runs Newton on a monomial system with real double arithmetic.
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
 *   nbsteps   the number of Newton steps;
 *   mode      the mode of execution, 0 for GPU only, 1 for CPU only,
 *             2 for CPU and GPU;
 *   vrblvl    is the verbose level. */

int test_cmplx_real_newton
 ( int szt, int nbt, int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   int nbsteps, int mode, int vrblvl );
/*
 * DESCRIPTION :
 *   Runs Newton on a monomial system with complex double arithmetic.
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
 *   nbsteps   the number of Newton steps;
 *   mode      the mode of execution, 0 for GPU only, 1 for CPU only,
 *             2 for CPU and GPU;
 *   vrblvl    is the verbose level. */

#endif
