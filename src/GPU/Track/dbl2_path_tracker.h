// The file dbl2_path_tracker.h specifies a path tracker
// on series in double double precision on real numbers.

#ifndef __dbl2_path_tracker_h__
#define __dbl2_path_tracker_h__

int dbl2_run_newton
 ( int szt, int nbt, int dim, int deg, int nbrcol, int nbsteps,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac,
   double **rhshi, double **rhslo,
   double ***cffhi, double ***cfflo, double **acchi, double **acclo,
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
   double *workvechi, double *workveclo,
   double **resvechi, double **resveclo, double *resmaxhi, double *resmaxlo,
   int vrblvl, int mode );
/*
 * DESCRIPTION :
 *   Runs Newton's method on power series on real data
 *   in double double precision.
 *
 * REQUIRED : szt*nbt = dim for GPU computing.
 *
 * ON ENTRY :
 *   szt       size of each tile and block;
 *   nbt       number of tiles and number of blocks;
 *   dim       number of monomials;
 *   deg       degree of the power series;
 *   nbrcol    is the number of columns, if 1, then the system is monomial,
 *             otherwise nbrcol columns are given on input;
 *   nbsteps   maximum number of Newton steps;
 *   nvr       nvr[i][j] is the number of variables in the j-th monomial
 *             of the i-th column;
 *   idx       idx[i][j] are the indices of the variables in monomial j
 *             of the i-th column;
 *   exp       exp[i] are the exponents of the variables in monomial i;
 *   nbrfac    nbrfac[i] are the number of exponents > 1 in monomial i;
 *   expfac    expfac[i] are the exponents in the i-th polynomial
 *             that are larger than one, minus one in the factor,
 *             if exp[i][k] > 1, then expfac[i][k] = exp[i][k] - 1;
 *   rhshi     high doubles of the right hand side vector of series;
 *   rhslo     low doubles of the right hand side vector of series;
 *   cffhi     high double coefficients of the monomials in each column;
 *   cfflo     low double coefficients of the monomials in each column;
 *   acchi     space to accumulate dim+1 power series of degree deg;
 *   acclo     space to accumulate dim+1 power series of degree deg;
 *   inputhi_h are the high double coefficients of the series of degree deg,
 *             for dim variables, computed on host;
 *   inputlo_h are the low double coefficients of the series of degree deg,
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
 *   resvechi  high doubles of the residual vectors;
 *   resveclo  low doubles of the residual vectors;
 *   resmaxhi  high double of the maximum element of the residual vectors;
 *   resmaxlo  low double of the maximum element of the residual vectors. */

int test_dbl2_real_track
 ( int szt, int nbt, int dim, int deg, int nbrcol,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac, int **rowsA,
   int nbsteps, int mode, int vrblvl );
/*
 * DESCRIPTION :
 *   Tracks a path on a system with real double double arithmetic.
 *
 * ON ENTRY :
 *   szt       size of each tile and block;
 *   nbt       number of tiles and number of blocks;
 *   dim       number of monomials;
 *   deg       degree of the power series;
 *   nbrcol    number of columns, if 1, then the system is monomial,
 *             otherwise nbrcol columns are expected;
 *   nvr       nvr[i][j] is the number of variables in the j-th monomial
 *             of the i-th column;
 *   idx       idx[i][j] are the indices of the variables in monomial j
 *             of the i-th column;
 *   exp       exp[i] are the exponents of the variables in monomial i;
 *   nbrfac    nbrfac[i] are the number of exponents > 1 in monomial i;
 *   expfac    expfac[i] are the exponents in the i-th polynomial
 *             that are larger than one, minus one in the factor,
 *             if exp[i][k] > 1, then expfac[i][k] = exp[i][k] - 1;
 *   rowsA     rows of the exponents of the dim monomials;
 *   nbsteps   the number of Newton steps;
 *   mode      the mode of execution, 0 for GPU only, 1 for CPU only,
 *             2 for CPU and GPU;
 *   vrblvl    is the verbose level. */

#endif
