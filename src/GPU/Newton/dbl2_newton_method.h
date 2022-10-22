// The file dbl2_newton_method.h specifies Newton's method
// on series in double double precision on real numbers.

#ifndef __dbl2_newton_method_h__
#define __dbl2_newton_method_h__

void dbl2_newton_qrstep
 ( int szt, int nbt, int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **mbhi, double **mblo, double dpr,
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
   bool *noqr_h, bool *noqr_d,
   int *upidx_h, int *bsidx_h, int *upidx_d, int *bsidx_d,
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
 *   mbhi      high doubles of the right hand side vector of series;
 *   mblo      low doubles of the right hand side vector of series;
 *   dpr       damper multiplier for t, should be in (0.0, 1.0];
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
 *   noqr_h    flag if true, then no qr on host;
 *   noqr_d    flag if true, then no qr on device;
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
 *   resmaxlo  low double of the maximum element of the residual vectors;
 *   noqr_h    updated flag if ||dx_0|| is zero for the first time on host;
 *   noqr_d    updated flag if ||dx_0|| is zero for the first time on device;
 *   upidx_h   counts the number of updates skipped by host;
 *   bsidx_h   counts the number of backsubstitutions skipped by host;
 *   upidx_d   counts the number of updates skipped by device;
 *   bsidx_d   counts the number of backsubstitutions skipped by device. */

int test_dbl2_real_newton
 ( int szt, int nbt, int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac, int **rowsA,
   double dpr, int nbsteps, int mode, int vrblvl );
/*
 * DESCRIPTION :
 *   Runs Newton on a monomial system with real double double arithmetic.
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
 *   rowsA     rows of the exponents of the dim monomials;
 *   dpr       damper multiplier for t, should be in (0.0, 1.0];
 *   nbsteps   the number of Newton steps;
 *   mode      the mode of execution, 0 for GPU only, 1 for CPU only,
 *             2 for CPU and GPU;
 *   vrblvl    is the verbose level. */

#endif
