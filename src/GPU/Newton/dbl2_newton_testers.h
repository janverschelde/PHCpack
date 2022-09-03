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

void dbl2_newton_lustep
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

void dbl2_newton_qrstep
 ( int szt, int nbt, int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffhi, double **cfflo, double *acchi, double *acclo,
   double **inputhi_h, double **inputlo_h,
   double **inputhi_d, double **inputlo_d,
   double ***outputhi_h, double ***outputlo_h,
   double ***outputhi_d, double ***outputlo_d,
   double **funvalhi_h, double **funvallo_h,
   double **funvalhi_d, double **funvallo_d,
   double ***jacvalhi_h, double ***jacvallo_h,
   double ***jacvalhi_d, double ***jacvallo_d,
   double **rhshi_h, double **rhslo_h, double **rhshi_d, double **rhslo_d,
   double **urhshi_h, double **urhslo_h, double **urhshi_d, double **urhslo_d,
   double **solhi_h, double **sollo_h, double **solhi_d, double **sollo_d,
   double **Qhi_h, double **Qlo_h, double **Qhi_d, double **Qlo_d,
   double **Rhi_h, double **Rlo_h, double **Rhi_d, double **Rlo_d,
   double **workmathi, double **workmatlo,
   double *workvechi, double *workveclo,
   double **resvechi, double **resveclo, double *resmaxhi, double *resmaxlo,
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
 *   cffhi     high doubles of the coefficients of the monomials;
 *   cfflo     low doubles of the coefficients of the monomials;
 *   acchi     space for high doubles of one power series of degree deg;
 *   acclo     space for low doubles of one power series of degree deg;
 *   inputhi_h has high doubles of coefficients of the series of degree deg,
 *             for dim variables, computed on host;
 *   inputlo_h has low doubles of  coefficients of the series of degree deg,
 *             for dim variables, computed on host;
 *   inputhi_d has space for power series computed on device;
 *   inputlo_d has space for power series computed on device;
 *   outputhi_h has space for the high doubles of the evaluated and 
 *             differentiated monomials, computed on the host;
 *   outputlo_h has space for the low doubles of the evaluated and 
 *             differentiated monomials, computed on the host;
 *   outputhi_d has space for the high doubles of the evaluated and
 *             differentiated monomials, computed on the device;
 *   outputlo_d has space for the low doubles of the evaluated and
 *             differentiated monomials, computed on the device;
 *   funvalhi_h has space for the evaluated power series computed by host;
 *   funvallo_h has space for the evaluated power series computed by host;
 *   funvalhi_d has space for the evaluated power series computed by device;
 *   funvallo_d has space for the evaluated power series computed by device;
 *   jacvalhi_h has space for deg+1 matrices of dimension dim on host;
 *   jacvallo_h has space for deg+1 matrices of dimension dim on host;
 *   jacvalhi_d has space for deg+1 matrices of dimension dim on device;
 *   jacvallo_d has space for deg+1 matrices of dimension dim on device;
 *   rhshi_h   space for deg+1 vectors of dimension dim on host;
 *   rhslo_h   space for deg+1 vectors of dimension dim on host;
 *   rhshi_d   space for deg+1 vectors of dimension dim on device;
 *   rhslo_d   space for deg+1 vectors of dimension dim on device;
 *   urhshi_h  space for updated right hand side vectors computed by host;
 *   urhslo_h  space for updated right hand side vectors computed by host;
 *   urhshi_d  space for updated right hand side vectors computed by device; 
 *   urhslo_d  space for updated right hand side vectors computed by device; 
 *   solhi_h   space for deg+1 vectors of dimension dim;
 *   sollo_h   space for deg+1 vectors of dimension dim;
 *   solhi_d   space for deg+1 vectors of dimension dim;
 *   sollo_d   space for deg+1 vectors of dimension dim;
 *   Qhi_h     space allocated for the Q computed by the host;
 *   Qlo_h     space allocated for the Q computed by the host;
 *   Qhi_d     space allocated for the Q computed by the device;
 *   Qlo_d     space allocated for the Q computed by the device;
 *   Rhi_h     space allocated for the R computed by the host;
 *   Rlo_h     space allocated for the R computed by the host;
 *   Rhi_d     space allocated for the R computed by the device;
 *   Rlo_d     space allocated for the R computed by the device;
 *   wrkmathi  has work space allocated for a matrix of dimension dim;
 *   wrkmatlo  has work space allocated for a matrix of dimension dim;
 *   wrkvechi  has work space allocated for a vector of dimension dim;
 *   wrkveclo  has work space allocated for a vector of dimension dim;
 *   resvechi  has space for deg+1 vectors of dimension dim;
 *   resveclo  has space for deg+1 vectors of dimension dim;
 *   vrblvl    is the verbose level;
 *   mode      execution mode, 0 (GPU only), 1 (CPU only) or 2 (GPU+CPU).
 *
 * ON RETURN :
 *   inputhi_h has high doubles of series computed on host (mode 1 or 2);
 *   inputlo_h has low doubles of series computed on host (mode 1 or 2);
 *   inputhi_d has high doubles of series, computed by device (mode 0 or 2);
 *   inputlo_d has low doubles of series, computed by device (mode 0 or 2);
 *   outputhi_h has high doubles of evaluated series computed on host;
 *   outputlo_h has low doubles of evaluated series computed on host;
 *   outputhi_d has high doubles of evaluated series computed on device;
 *   outputlo_d has low doubles of evaluated series computed on device;
 *   funvalhi_h has high doubles of output[i][dim], function values on host;
 *   funvallo_h has low doubles of output[i][dim], function values on host;
 *   funvalhi_d has high doubles of output[i][dim], function values on device;
 *   funvallo_d has high doubles of output[i][dim], function values on device;
 *   jacvalhi_h has high doubles of matrices, computed by host,
 *             the leading coefficient is the Jacobian matrix;
 *   jacvallo_h has low doubles of matrices, computed by host,
 *             the leading coefficient is the Jacobian matrix;
 *   jacvalhi_d has high doubles of matrices, computed by device,
 *             the leading coefficient is the Jacobian matrix;
 *   jacvallo_d has high doubles of matrices, computed by device,
 *             the leading coefficient is the Jacobian matrix;
 *   rhshi_h   high doubles of the right hand side are linearized values
 *             subtracted by 1 and added by t, computed by host;
 *   rhslo_h   low doubles of the right hand side are linearized values
 *             subtracted by 1 and added by t, computed by host;
 *   rhshi_d   high doubles of the right hand side are linearized values
 *             subtracted by 1 and added by t, computed by device;
 *   rhslo_d   low doubles of the right hand side are linearized values
 *             subtracted by 1 and added by t, computed by device;
 *   urhshi_h  high doubles of right hand side vector updated by the host;
 *   urhslo_h  low doubles of right hand side vector updated by the host;
 *   urhshi_d  high doubles of right hand side vector updated by the device;
 *   urhslo_d  low doubles of right hand side vector updated by the device;
 *   solhi_h   high doubles of solution computed by the host;
 *   sollo_h   low doubles of solution computed by the host;
 *   solhi_d   high doubles of solution computed by the device;
 *   sollo_d   low doubles of solution computed by the device;
 *   Qhi_h     high doubles of Q of the QR computed by the host;
 *   Qlo_h     low doubles of Q of the QR computed by the host;
 *   Qhi_d     high doubles of Q of the QR computed by the device;
 *   Qlo_d     low doubles of Q of the QR computed by the device;
 *   Rhi_h     high doubles of R of the QR computed by the host;
 *   Rlo_h     low doubles of R of the QR computed by the host;
 *   Rhi_d     high doubles of R of the QR computed by the device;
 *   Rlo_d     low doubles of R of the QR computed by the device;
 *   wrkmathi  has a copy of the high doubles of the Jacobian matrix;
 *   wrkmatlo  has a copy of the low doubles of the Jacobian matrix;
 *   resvechi  high doubles of the residual vectors;
 *   resveclo  low doubles of the residual vectors;
 *   resmaxhi  high double of the maximum element of the residual vectors;
 *   resmaxlo  low double of the maximum element of the residual vectors. */

#endif
