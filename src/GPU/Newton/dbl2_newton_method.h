// The file dbl2_newton_method.h specifies Newton's method
// on series in double double precision on real numbers.

#ifndef __dbl2_newton_method_h__
#define __dbl2_newton_method_h__

int dbl2_errors_funjacrhs
 ( int dim, int deg,
   double **funvalhi_h, double **funvallo_h,
   double **funvalhi_d, double **funvallo_d,
   double ***jacvalhi_h, double ***jacvallo_h,
   double ***jacvalhi_d, double ***jacvallo_d,
   double **rhshi_h, double **rhslo_h, double **rhshi_d, double **rhslo_d,
   int vrblvl );
/*
 * DESCRIPTION :
 *   Compares results on host with those on the device
 *   and returns the number of errors larger than 1.0e-20.
 *
 * ON ENTRY :
 *   dim       number of equations and variables;
 *   deg       degree at which the series are truncated;
 *   funvalhi_h are the highest function values on the host;
 *   funvallo_h are the lowest function values on the host;
 *   funvalhi_d are the highest function values on the device;
 *   funvallo_d are the lowest function values on the device;
 *   jacvalhi_h are the highest doubles of the Jacobian matrix on the host;
 *   jacvallo_h are the lowest doubles of the Jacobian matrix on the host;
 *   jacvalhi_d are the highest doubles of the Jacobian matrix on the device;
 *   jacvallo_d are the lowest doubles of the Jacobian matrix on the device;
 *   rhshi_h   highest doubles of the right hand side vector on the host;
 *   rhslo_h   lowest doubles of the right hand side vector on the host;
 *   rhshi_d   highest doubles of the right hand side vector on the device;
 *   rhslo_d   lowest doubles of the right hand side vector on the device;
 *   vrblvl    verbose level. */

int dbl2_errors_inurhsQRsol
 ( int dim, int deg,
   double **inputhi_h, double **inputlo_h,
   double **inputhi_d, double **inputlo_d,
   double **Qhi_h, double **Qlo_h, double **Qhi_d, double **Qlo_d,
   double **Rhi_h, double **Rlo_h, double **Rhi_d, double **Rlo_d,
   double **urhshi_h, double **urhslo_h, double **urhshi_d, double **urhslo_d,
   double **solhi_h, double **sollo_h, double **solhi_d, double **sollo_d,
   int vrblvl );
/*
 * DESCRIPTION :
 *   Compares results on host with those on the device
 *   and returns the number of errors larger than 1.0e-8.
 *
 * ON ENTRY :
 *   dim       number of equations and variables;
 *   deg       degree at which the series are truncated;
 *   inputhi_h highest doubles of updated series on the host;
 *   inputlo_h lowest doubles of updated series on the host;
 *   inputhi_d highest doubles of updated series on the device;
 *   inputlo_d lowest doubles of updated series on the device;
 *   Qhi_h     highest doubles of Q on the host;
 *   Qlo_h     lowest doubles of Q on the host;
 *   Qhi_d     highest doubles of Q on the device;
 *   Qlo_d     lowest doubles of Q on the device;
 *   Rhi_h     highest doubles of R on the host;
 *   Rlo_h     lowest doubles of R on the host;
 *   Rhi_d     highest doubles of R on the device;
 *   Rlo_d     lowest doubles of R on the device;
 *   urhshi_h  highest doubles of updated right hand side on the host;
 *   urhslo_h  lowest doubles of updated right hand side on the host;
 *   urhshi_d  highest doubles of updated right hand side on the device;
 *   urhslo_d  lowest doubles of updated right hand side on the device;
 *   solhi_h   highest doubles of update to the solution on the host,
 *   sollo_h   lowest doubles of update to the solution on the host,
 *   solhi_d   highest doubles of update to the solution on the device;
 *   sollo_d   lowest doubles of update to the solution on the device;
 *   vrblvl    is the verbose level. */

int dbl2_update_newton_qrstep
 ( int szt, int nbt, int dim, int deg,
   int *tailidx_h, int *tailidx_d,
   double **inputhi_h, double **inputlo_h,
   double **inputhi_d, double **inputlo_d,
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
   bool *zeroQ_h, bool *noqr_h, bool *zeroQ_d, bool *noqr_d,
   int *upidx_h, int *bsidx_h, int *upidx_d, int *bsidx_d,
   double *totqrlapsedms, double *totqtblapsedms,
   double *totbslapsedms, double *totupdlapsedms, double *totreslapsedms,
   int vrblvl, int mode );
/*
 * DESCRIPTION :
 *   Given the function values and the matrix series of the Jacobian,
 *   computes the update to the solution series.
 *
 * REQUIRED : szt*nbt = dim for GPU acceleration, when mode is 2.
 *
 * ON ENTRY :
 *   szt       size of each tile and block;
 *   nbt       number of tiles and number of blocks;
 *   dim       number of monomials;
 *   deg       degree of the power series;
 *   tailidx_h is the start index of the update of the tail on the host;
 *   tailidx_d is the start index of the update of the tail on the device;
 *   inputhi_h has high doubles of coefficients of the series of degree deg,
 *             for dim variables, computed on host;
 *   inputlo_h has low doubles of  coefficients of the series of degree deg,
 *             for dim variables, computed on host;
 *   inputhi_d has space for power series computed on device;
 *   inputlo_d has space for power series computed on device;
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
 *   zeroQ_h   if true, then Q is zero and Q must be computed on host;
 *   noqr_h    flag if true, then no qr on host;
 *   zeroQ_d   if true, then Q is zero and Q must be computed on device;
 *   noqr_d    flag if true, then no qr on device;
 *   totqrlapsedms accumulates the milliseconds spent on the Householder QR;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs;
 *   totbslapsedms accumulates the milliseconds spent on back substitutions;
 *   totupdlapsedms accumulates the milliseconds spent on updates;
 *   totreslapsedms accumulates the milliseconds spent on residuals;
 *   vrblvl    is the verbose level;
 *   mode      execution mode, 0 (GPU only), 1 (CPU only) or 2 (GPU+CPU).
 *
 * ON RETURN :
 *   tailidx_h is the updated value for tailidx_h;
 *   tailidx_d is the updated value for tailidx_d;
 *   inputhi_h has high doubles of series computed on host (mode 1 or 2);
 *   inputlo_h has low doubles of series computed on host (mode 1 or 2);
 *   inputhi_d has high doubles of series, computed by device (mode 0 or 2);
 *   inputlo_d has low doubles of series, computed by device (mode 0 or 2);
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
 *   resmaxlo  low double of the maximum element of the residual vectors;
 *   zeroQ_h   false if Q was computed on host;
 *   noqr_h    updated flag if ||dx_0|| is zero for the first time on host;
 *   zeroQ_d   false if Q was computed on device;
 *   noqr_d    updated flag if ||dx_0|| is zero for the first time on device;
 *   upidx_h   counts the number of updates skipped by host;
 *   bsidx_h   counts the number of backsubstitutions skipped by host;
 *   upidx_d   counts the number of updates skipped by device;
 *   bsidx_d   counts the number of backsubstitutions skipped by device;
 *   totqrlapsedms accumulates the milliseconds spent on the Householder QR;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs;
 *   totbslapsedms accumulates the milliseconds spent on back substitutions;
 *   totupdlapsedms accumulates the milliseconds spent on updates;
 *   totreslapsedms accumulates the milliseconds spent on residuals. */

int dbl2_column_newton_qrstep
 ( int szt, int nbt, int dim, int deg, int nbrcol,
   int *tailidx_h, int *tailidx_d,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac,
   double **mbhi, double **mblo, double dpr,
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
   bool* zeroQ_h, bool *noqr_h, bool* zeroQ_d, bool *noqr_d,
   int *upidx_h, int *bsidx_h, int *upidx_d, int *bsidx_d,
   double *totcnvlapsedms, double *totqrlapsedms, double *totqtblapsedms,
   double *totbslapsedms, double *totupdlapsedms, double *totreslapsedms,
   int vrblvl, int mode );
/*
 * DESCRIPTION :
 *   Does one step with Newton's method to update a power series,
 *   using QR factorization to solve linear systems,
 *   on one or more columns of monomials.
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
 *   tailidx_h is the start index of the update of the tail on the host;
 *   tailidx_d is the start index of the update of the tail on the device;
 *   nvr       nvr[i][j] is the number of variables in the j-th monomial
 *             of the i-th column;
 *   idx       idx[i][j] are the indices of the variables in monomial j
 *             of the i-th column;
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
 *   acchi     space for high doubles of dim+1 power series of degree deg;
 *   acclo     space for low doubles of dim+1 power series of degree deg;
 *   inputhi_h has high doubles of coefficients of the series of degree deg,
 *             for dim variables, computed on host;
 *   inputlo_h has low doubles of coefficients of the series of degree deg,
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
 *   zeroQ_h   if true, then Q is zero and Q must be computed on host;
 *   noqr_h    flag if true, then no qr on host;
 *   zeroQ_d   if true, then Q is zero and Q must be computed on device;
 *   noqr_d    flag if true, then no qr on device;
 *   totcnvlapsedms accumulates the milliseconds spent on the convolutions;
 *   totqrlapsedms accumulates the milliseconds spent on the Householder QR;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs;
 *   totbslapsedms accumulates the milliseconds spent on back substitutions;
 *   totupdlapsedms accumulates the milliseconds spent on updates;
 *   totreslapsedms accumulates the milliseconds spent on residuals;
 *   vrblvl    is the verbose level;
 *   mode      execution mode, 0 (GPU only), 1 (CPU only) or 2 (GPU+CPU).
 *
 * ON RETURN :
 *   tailidx_h is the new value for tailidx_h;
 *   tailidx_d is the new value for tailidx_d;
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
 *   resmaxlo  low double of the maximum element of the residual vectors;
 *   zeroQ_h   false if Q was computed on host;
 *   noqr_h    updated flag if ||dx_0|| is zero for the first time on host;
 *   zeroQ_d   false if Q was computed on device;
 *   noqr_d    updated flag if ||dx_0|| is zero for the first time on device;
 *   upidx_h   counts the number of updates skipped by host;
 *   bsidx_h   counts the number of backsubstitutions skipped by host;
 *   upidx_d   counts the number of updates skipped by device;
 *   bsidx_d   counts the number of backsubstitutions skipped by device;
 *   totcnvlapsedms accumulates the milliseconds spent on the convolutions;
 *   totqrlapsedms accumulates the milliseconds spent on the Householder QR;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs;
 *   totbslapsedms accumulates the milliseconds spent on back substitutions;
 *   totupdlapsedms accumulates the milliseconds spent on updates;
 *   totreslapsedms accumulates the milliseconds spent on residuals. */

int dbl2_row_newton_qrstep
 ( int szt, int nbt, int dim, int deg, int *tailidx_h, int *tailidx_d,
   int *nbr, int **nvr, int ***idx,
   double **csthi, double **cstlo,
   double ***cffhi, double ***cfflo, double dpr,
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
   bool *zeroQ_h, bool *noqr_h, bool *zeroQ_d, bool *noqr_d,
   int *upidx_h, int *bsidx_h, int *upidx_d, int *bsidx_d,
   double *totcnvlapsedms, double *totqrlapsedms, double *totqtblapsedms,
   double *totbslapsedms, double *totupdlapsedms, double *totreslapsedms,
   int vrblvl, int mode );
/*
 * DESCRIPTION :
 *   Does one step with Newton's method to update a power series,
 *   using QR factorization to solve linear systems, on real data,
 *   on an indexed polynomial system.
 *
 * REQUIRED : szt*nbt = dim for GPU acceleration, when mode is 2.
 *
 * ON ENTRY :
 *   szt       size of each tile and block;
 *   nbt       number of tiles and number of blocks;
 *   dim       number of monomials;
 *   deg       degree of the power series;
 *   tailidx_h is the start index of the update of the tail on the host;
 *   tailidx_d is the start index of the update of the tail on the device;
 *   nbr       nbr[i] is the number of monomials in the i-th polynomial;
 *   nvr       nvr[i][j] is the number of variables in the j-th monomial
 *             of the i-th polynomial;
 *   idx       idx[i][j] are the indices of the variables in monomial j
 *             of the i-th polynomial;
 *   csthi     high doubles of the constant coefficients of the system;
 *   cstlo     low doubles of the constant coefficients of the system;
 *   cffhi     high doubles of the coefficients in the system;
 *   cfflo     low doubles of the coefficients in the system;
 *   dpr       damper multiplier for t, should be in (0.0, 1.0];
 *   inputhi_h has high doubles of coefficients of the series of degree deg,
 *             for dim variables, computed on host;
 *   inputlo_h has low doubles of coefficients of the series of degree deg,
 *             for dim variables, computed on host;
 *   inputhi_d has space for power series computed on device;
 *   inputlo_d has space for power series computed on device;
 *   outputhi_h has space for the high doubles of the evaluated
 *             and differentiated monomials, computed on the host;
 *   outputlo_h has space for the low doubles of the evaluated
 *             and differentiated monomials, computed on the host;
 *   outputhi_d has space for the high doubles of the evaluated
 *             and differentiated monomials, computed on the device;
 *   outputlo_d has space for the low doubles of the evaluated
 *             and differentiated monomials, computed on the device;
 *   funvalhi_h has space for series computed by host;
 *   funvallo_h has space for series computed by host;
 *   funvalhi_d has space for series computed by device;
 *   funvallo_d has space for series computed by device;
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
 *   wrkvechi  work space allocated for a vector of dimension dim;
 *   wrkveclo  work space allocated for a vector of dimension dim;
 *   resvechi  space for deg+1 vectors of dimension dim;
 *   resveclo  space for deg+1 vectors of dimension dim;
 *   zeroQ_h   if true, then Q is zero and Q must be computed on host;
 *   noqr_h    flag if true, then no qr on host;
 *   zeroQ_d   if true, then Q is zero and Q must be computed on device;
 *   noqr_d    flag if true, then no qr on device;
 *   totcnvlapsedms accumulates the milliseconds spent on the convolutions;
 *   totqrlapsedms accumulates the milliseconds spent on the Householder QR;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs;
 *   totbslapsedms accumulates the milliseconds spent on back substitutions;
 *   totupdlapsedms accumulates the milliseconds spent on updates;
 *   totreslapsedms accumulates the milliseconds spent on residuals;
 *   vrblvl    is the verbose level;
 *   mode      execution mode, 0 (GPU only), 1 (CPU only) or 2 (GPU+CPU).
 *
 * ON RETURN :
 *   tailidx_h is the updated value for tailidx_h;
 *   tailidx_d is the updated value for tailidx_d;
 *   inputhi_h are high doubles computed on host (depending on mode);
 *   inputlo_h are low doubles computed on host (depending on mode);
 *   inputhi_d are high doubles computed on device (depending on mode);
 *   inputlo_d are low doubles computed on device (depending on mode);
 *   outputhi_h are high doubles computed on host (depending on mode);
 *   outputlo_h are low doubles computed on host (depending on mode);
 *   outputhi_d are high doubles computed on device (depending on mode);
 *   outputlo_d are low doubles computed on device (depending on mode);
 *   funvalhi_h are high doubles of function values on host;
 *   funvallo_h are low doubles of function values on host;
 *   funvalhi_d are high doubles of function values on device;
 *   funvallo_d are low doubles of function values on device;
 *   jacvalhi_h are the high doubles of Jacobian computed on host;
 *   jacvallo_h are the high doubles of Jacobian computed on host;
 *   jacvalhi_d are the high doubles of Jacobian computed on device;
 *   jacvallo_d are the high doubles of Jacobian computed on device;
 *   rhshi_h   high doubles of linearized right hand side on host;
 *   rhslo_h   low doubles of linearized right hand side on host;
 *   rhshi_d   high doubles of linearized right hand side on device;
 *   rhslo_d   low doubles of linearized right hand side on device;
 *   urhshi_h  high doubles of right hand side vector updated by the host;
 *   urhslo_h  low doubles of right hand side vector updated by the host;
 *   urhshi_d  high doubles of right hand side vector updated by the device;
 *   urhslo_d  low doubles of right hand side vector updated by the device;
 *   solhi_h   high doubles of the solution computed by the host;
 *   sollo_h   low doubles of the solution computed by the host;
 *   solhi_d   high doubles of the solution computed by the device;
 *   sollo_d   low doubles of the solution computed by the device;
 *   Qhi_h     high doubles of the Q computed by the host;
 *   Qlo_h     low doubles of the Q computed by the host;
 *   Qhi_d     high doubles of the Q computed by the device;
 *   Qlo_d     low doubles of the Q computed by the device;
 *   Rhi_h     high doubles of the R computed by the host;
 *   Rlo_h     low doubles of the R computed by the host;
 *   Rhi_d     high doubles of the R computed by the device;
 *   Rlo_d     low doubles of the R computed by the device;
 *   resvechi  high doubles of the residual vectors;
 *   resveclo  low doubles of the residual vectors;
 *   resmaxhi  high double of the maximum element of the residual vectors;
 *   resmaxlo  low double of the maximum element of the residual vectors;
 *   zeroQ_h   false if Q was computed on host;
 *   noqr_h    updated flag if ||dx_0|| is zero for the first time on host;
 *   zeroQ_d   false if Q was computed on device;
 *   noqr_d    updated flag if ||dx_0|| is zero for the first time on device;
 *   upidx_h   counts the number of updates skipped by host;
 *   bsidx_h   counts the number of backsubstitutions skipped by host;
 *   upidx_d   counts the number of updates skipped by device;
 *   bsidx_d   counts the number of backsubstitutions skipped by device;
 *   totcnvlapsedms accumulates the milliseconds spent on the convolutions;
 *   totqrlapsedms accumulates the milliseconds spent on the Householder QR;
 *   totqtblapsedms accumulates the milliseconds spent on Q times rhs;
 *   totbslapsedms accumulates the milliseconds spent on back substitutions;
 *   totupdlapsedms accumulates the milliseconds spent on updates;
 *   totreslapsedms accumulates the milliseconds spent on residuals. */

int dbl2_allocate_inoutfunjac
 ( int dim, int deg, int mode,
   double **inputhi_h, double **inputlo_h,
   double **inputhi_d, double **inputlo_d,
   double ***outputhi_h, double ***outputlo_h, 
   double ***outputhi_d, double ***outputlo_d,
   double **funvalhi_h, double **funvallo_h,
   double **funvalhi_d, double **funvallo_d,
   double ***jacvalhi_h, double ***jacvallo_h,
   double ***jacvalhi_d, double ***jacvallo_d );
/*
 * DESCRIPTION :
 *   Allocates work space memory for input, output,
 *   the function values and the value of the Jacobian matrix.
 *
 * ON ENTRY :
 *   dim        dimension of the system;
 *   deg        degree at which the series are truncated;
 *   mode       0 (GPU), 1 (CPU), or 2 (GPU+CPU).
 *
 * ON RETURN :
 *   inputhi_h  high doubles of the input on the host;
 *   inputlo_h  low doubles of the input on the host;
 *   inputhi_d  high doubles of the input on the device;
 *   inputlo_d  low doubles of the input on the device;
 *   outputhi_h are the high doubles of the output on the host;
 *   outputlo_h are the low doubles of the output on the host;
 *   outputhi_d are the high doubles of the output on the device;
 *   outputlo_d are the low doubles of the output on the device;
 *   funvalhi_h are the high doubles of the function values on the host;
 *   funvallo_h are the low doubles of the function values on the host;
 *   funvalhi_d are the high doubles of the function values on the device;
 *   funvallo_d are the low doubles of the function values on the device;
 *   jacvalhi_h are the high doubles of the Jacobian on the host;
 *   jacvallo_h are the low doubles of the Jacobian on the host;
 *   jacvalhi_d are the high doubles of the Jacobian on the device;
 *   jacvallo_d are the low doubles of the Jacobian on the device. */

int dbl2_allocate_rhsqrsol
 ( int dim, int deg, int mode,
   double **rhshi_h, double **rhslo_h, double **rhshi_d, double **rhslo_d,
   double **urhshi_h, double **urhslo_h, double **urhshi_d, double **urhslo_d,
   double **Qhi_h, double **Qlo_h, double **Qhi_d, double **Qlo_d, 
   double **Rhi_h, double **Rlo_h, double **Rhi_d, double **Rlo_d,
   double **solhi_h, double **sollo_h, double **solhi_d, double **sollo_d );
/*
 * DESCRIPTION :
 *   Allocates work space memory for the linearized power series system.
 *
 * ON ENTRY :
 *   dim        dimension of the system;
 *   deg        degree at which the series are truncated;
 *   mode       0 (GPU), 1 (CPU), or 2 (GPU+CPU).
 *
 * ON RETURN :
 *   rhshi_h    high doubles ofright-hand side on the host;
 *   rhslo_h    low doubles of right-hand side on the host;
 *   rhshi_d    high doubles of right-hand side on the device;
 *   rhslo_d    low doubles of right-hand side on the device;
 *   urhshi_h   high doubles of updated right-hand side on the host;
 *   urhslo_h   low doubles of updated right-hand side on the host;
 *   urhshi_d   high doubles of updated right-hand side on the device;
 *   urhslo_d   low doubles of updated right-hand side on the device;
 *   Qhi_h      high doubles of Q on the host;
 *   Qlo_h      low doubles of Q on the host;
 *   Qhi_d      high doubles of Q on the device;
 *   Qlo_d      low doubles of Q on the device;
 *   Rhi_h      high doubles of R on the host;
 *   Rlo_h      low doubles of R on the host;
 *   Rhi_d      high doubles of R on the device;
 *   Rlo_d      low doubles of R on the device;
 *   solhi_h    high doubles of update to the solution on the host;
 *   sollo_h    low doubles of update to the solution on the host;
 *   solhi_h    high doubles of update to the solution of the device;
 *   sollo_h    low doubles of update to the solution of the device. */

void dbl2_start_setup
 ( int dim, int deg, double **testsolhi, double **testsollo,
   double **inputhi_h, double **inputlo_h,
   double **inputhi_d, double **inputlo_d, int mode, int vrblvl );
/*
 * DESCRIPTION :
 *   Given the test solution, defines the start vector.
 *
 * ON ENTRY :
 *   dim        dimension of the solution vectors;
 *   deg        degree at which the series are truncated;
 *   testsolhi  high doubles of the test solution;
 *   testsollo  low doubles of the test solution;
 *   inputhi_h  allocated on host if mode is 1 or 2;
 *   inputlo_h  allocated on host if mode is 1 or 2;
 *   inputhi_d  allocated on device if mode is 0 or 2;
 *   inputlo_d  allocated on device if mode is 0 or 2;
 *   mode       0 (GPU), 1 (CPU), or 2 (GPU+CPU);
 *   vrblvl     is the verbose level.
 *
 * ON RETURN :
 *   inputhi_h  high doubles of start vector for host if mode is 1 or 2;
 *   inputlo_h  low doubles of start vector for host if mode is 1 or 2;
 *   inputhi_d  high doubles of start vector for device if mode is 0 or 2;
 *   inputlo_d  low doubles of start vector for device if mode is 0 or 2. */

void dbl2_column_setup
 ( int dim, int deg, int nbrcol,
   int **nvr, int ***idx, int **rowsA, double ***cffhi, double ***cfflo,
   double **testsolhi, double **testsollo,
   double **mbrhshi, double **mbrhslo,
   double **inputhi_h, double **inputlo_h,
   double **inputhi_d, double **inputlo_d, int mode, int vrblvl );
/*
 * DESCRIPTION :
 *   Defines the test solution and start vectors to run Newton's method
 *   on one or more columns of monomials.
 *
 * ON ENTRY :
 *   dim        dimension of the solution vectors;
 *   deg        degree at which the series are truncated;
 *   nbrcol     number of columns, if 1, then the system is monomial,
 *              otherwise nbrcol columns are expected;
 *   nvr        nvr[i][j] is the number of variables in the j-th monomial
 *              of the i-th column;
 *   idx        idx[i][j] are the indices of the variables in monomial j
 *              of the i-th column;
 *   rowsA      exponents for monomials if only one column;
 *   cffhi      high doubles of coefficients, if more than one column;
 *   cfflo      low doubles of coefficients, if more than one column;
 *   testsolhi  space for dim pointers;
 *   testsollo  space for dim pointers;
 *   mbrhshi    space for dim pointers;
 *   mbrhslo    space for dim pointers;
 *   inputhi_h  allocated on host if mode is 1 or 2;
 *   inputlo_h  allocated on host if mode is 1 or 2;
 *   inputhi_d  allocated on device if mode is 0 or 2;
 *   inputlo_d  allocated on device if mode is 0 or 2;
 *   mode       0 (GPU), 1 (CPU), or 2 (GPU+CPU);
 *   vrblvl     is the verbose level.
 *
 * ON RETURN :
 *   testsolhi  high doubles of test solution;
 *   testsollo  low doubles of test solution;
 *   mbrhshi    high doubles of right hand side vector for the test solution;
 *   mbrhslo    low doubles of right hand side vector for the test solution;
 *   inputhi_h  high doubles of start vector for host if mode is 1 or 2;
 *   inputlo_h  low doubles of start vector for host if mode is 1 or 2;
 *   inputhi_d  high doubles of start vector for device if mode is 0 or 2;
 *   inputlo_d  low doubles of start vector for device if mode is 0 or 2. */

void dbl2_row_setup
 ( int dim, int deg, int *nbr, int **nvr, int ***idx,
   double **csthi, double **cstlo, double ***cffhi, double ***cfflo, 
   double **testsolhi, double **testsollo, 
   double **inputhi_h, double **inputlo_h,
   double **inputhi_d, double **inputlo_d,
   double ***outputhi_h, double ***outputlo_h,
   double ***outputhi_d, double ***outputlo_d, int mode, int vrblvl );
/*
 * DESCRIPTION :
 *   Defines the test solution and start vectors to run Newton's method
 *   on one or more columns of monomials.
 *
 * ON ENTRY :
 *   dim        dimension of the solution vectors;
 *   deg        degree at which the series are truncated;
 *   nbr        nbr[i] is the number of monomials in the i-th polynomial;
 *   nvr        nvr[i][j] is the number of variables in the j-th monomial
 *              of the i-th polynomial;
 *   idx        idx[i][j] are the indices of the variables in monomial j
 *              of the i-th polynomial;
 *   csthi      high doubles of constants of the polynomials;
 *   cstlo      low doubles of constants of the polynomials;
 *   cffhi      high doubles of coefficients of the monomials;
 *   cfflo      low doubles of coefficients of the monomials;
 *   testsolhi  space for dim pointers;
 *   testsollo  space for dim pointers;
 *   inputhi_h  allocated on host if mode is 1 or 2;
 *   inputlo_h  allocated on host if mode is 1 or 2;
 *   inputhi_d  allocated on device if mode is 0 or 2;
 *   inputlo_d  allocated on device if mode is 0 or 2;
 *   outpuhit_h allocated on host if mode is 1 or 2;
 *   outpulot_h allocated on host if mode is 1 or 2;
 *   outputhi_d allocated on device if mode is 0 or 2;
 *   outputlo_d allocated on device if mode is 0 or 2;
 *   mode       0 (GPU), 1 (CPU), or 2 (GPU+CPU);
 *   vrblvl     is the verbose level.
 *
 * ON RETURN :
 *   csthi      high doubles of constant, adjusted for the test solution;
 *   cstlo      low doubles of constant, adjusted for the test solution;
 *   testsolhi  high doubles of the test solution;
 *   testsollo  low doubles of the test solution;
 *   inputhi_h  high doubles of start vector for host if mode is 1 or 2;
 *   inputlo_h  low doubles of start vector for host if mode is 1 or 2;
 *   inputhi_d  high doubles of start vector for device if mode is 0 or 2;
 *   inputlo_d  low doubles of start vector for device if mode is 0 or 2;
 *   outputhi_h are the high doubles of evaluated test solution if mode is 1;
 *   outputlo_h are the low doubles of evaluated test solution if mode is 1;
 *   outputhi_d are the high doubles of evaluated test solution,
 *              if mode is 0 or 2;
 *   outputlo_d are the low doubles of evaluated test solution,
 *              if mode is 0 or 2. */

int dbl2_error_testsol
 ( int dim, int deg, int mode, double **testsolhi, double **testsollo,
   double **inputhi_h, double **inputlo_h,
   double **inputhi_d, double **inputlo_d );
/*
 * DESCRIPTION :
 *   Compares the computed solution against the test solution.
 *
 * ON ENTRY :
 *   dim        dimension of the solution vectors;
 *   deg        degree at which the series are truncated;
 *   mode       0 (GPU), 1 (CPU), or 2 (GPU+CPU).
 *   testsolhi  high doubles of the test solution;
 *   testsollo  low doubles of the test solution;
 *   inputhi_h  high doubles on host if mode is 1 or 2;
 *   inputlo_h  low doubles on host if mode is 1 or 2;
 *   inputhi_d  high doubles on device if mode is 0 or 2;
 *   inputlo_d  low doubles on device if mode is 0 or 2. */

int test_dbl2_column_newton
 ( int szt, int nbt, int dim, int deg, int nbrcol,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac, int **rowsA,
   double dpr, int nbsteps, int mode, int vrblvl );
/*
 * DESCRIPTION :
 *   Runs Newton with real double double arithmetic,
 *   on one or more columns of monomials.
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
 *   dpr       damper multiplier for t, should be in (0.0, 1.0];
 *   nbsteps   the number of Newton steps;
 *   mode      the mode of execution, 0 for GPU only, 1 for CPU only,
 *             2 for CPU and GPU;
 *   vrblvl    is the verbose level. */

int test_dbl2_row_newton
 ( int szt, int nbt, int dim, int deg, int *nbr, int **nvr, int ***idx,
   double dpr, int nbsteps, int mode, int vrblvl );
/*
 * DESCRIPTION :
 *   Runs Newton with real double double arithmetic,
 *   on an indexed polynomial system.
 *
 * ON ENTRY :
 *   szt       size of each tile and block;
 *   nbt       number of tiles and number of blocks;
 *   dim       number of equations and variables in the system;
 *   deg       degree of the power series;
 *   nbr       nbr[i] is the number of monomials in the i-th polynomial;
 *   nvr       nvr[i][j] is the number of variables in the j-th monomial
 *             of the i-th polynomial;
 *   idx       idx[i][j] are the indices of the variables in monomial j
 *             of the i-th polynomial;
 *   dpr       damper multiplier for t, should be in (0.0, 1.0];
 *   nbsteps   the number of Newton steps;
 *   mode      the mode of execution, 0 for GPU only, 1 for CPU only,
 *             2 for CPU and GPU;
 *   vrblvl    is the verbose level. */

#endif
