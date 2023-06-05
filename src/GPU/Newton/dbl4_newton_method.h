// The file dbl4_newton_method.h specifies Newton's method
// on series in quad double precision on real numbers.

#ifndef __dbl4_newton_method_h__
#define __dbl4_newton_method_h__

int dbl4_errors_funjacrhs
 ( int dim, int deg,
   double **funvalhihi_h, double **funvallohi_h,
   double **funvalhilo_h, double **funvallolo_h,
   double **funvalhihi_d, double **funvallohi_d,
   double **funvalhilo_d, double **funvallolo_d,
   double ***jacvalhihi_h, double ***jacvallohi_h,
   double ***jacvalhilo_h, double ***jacvallolo_h,
   double ***jacvalhihi_d, double ***jacvallohi_d,
   double ***jacvalhilo_d, double ***jacvallolo_d,
   double **rhshihi_h, double **rhslohi_h,
   double **rhshilo_h, double **rhslolo_h,
   double **rhshihi_d, double **rhslohi_d,
   double **rhshilo_d, double **rhslolo_d, int vrblvl );
/*
 * DESCRIPTION :
 *   Compares results on host with those on the device
 *   and returns the number of errors larger than 1.0e-50.
 *
 * ON ENTRY :
 *   dim       number of equations and variables;
 *   deg       degree at which the series are truncated;
 *   funvalhihi_h are the highest function values on the host;
 *   funvallohi_h are the second highest function values on the host;
 *   funvalhilo_h are the second lowest function values on the host;
 *   funvallolo_h are the lowest function values on the host;
 *   funvalhihi_d are the highest function values on the device;
 *   funvallohi_d are the second highest function values on the device;
 *   funvalhilo_d are the second lowest function values on the device;
 *   funvallolo_d are the lowest function values on the device;
 *   jacvalhihi_h are the highest doubles of the Jacobian on the host;
 *   jacvallohi_h are the second highest doubles of the Jacobian on the host;
 *   jacvalhilo_h are the second lowest doubles of the Jacobian on the host;
 *   jacvallolo_h are the lowest doubles of the Jacobian on the host;
 *   jacvalhihi_d are the highest doubles of the Jacobian on the device;
 *   jacvallohi_d are the second highest doubles of the Jacobian on the device;
 *   jacvalhilo_d are the second lowest doubles of the Jacobian on the device;
 *   jacvallolo_d are the lowest doubles of the Jacobian on the device;
 *   rhshihi_h are the highest doubles of the right hand side on the host;
 *   rhslohi_h are the 2nd highest doubles of the right hand side on the host;
 *   rhshilo_h are the 2nd lowest doubles of the right hand side on the host;
 *   rhslolo_h are the lowest doubles of the right hand side on the host;
 *   rhshihi_d are the highest doubles of the right hand side on the device;
 *   rhslohi_d are the 2nd highest doubles of the right hand side on the device;
 *   rhshilo_d are the 2nd lowest doubles of the right hand side on the device;
 *   rhslolo_d are the lowest doubles of the right hand side on the device;
 *   vrblvl    is the verbose level. */

int dbl4_errors_inurhsQRsol
 ( int dim, int deg,
   double **inputhihi_h, double **inputlohi_h,
   double **inputhilo_h, double **inputlolo_h,
   double **inputhihi_d, double **inputlohi_d,
   double **inputhilo_d, double **inputlolo_d,
   double **Qhihi_h, double **Qlohi_h, double **Qhilo_h, double **Qlolo_h,
   double **Qhihi_d, double **Qlohi_d, double **Qhilo_d, double **Qlolo_d,
   double **Rhihi_h, double **Rlohi_h, double **Rhilo_h, double **Rlolo_h,
   double **Rhihi_d, double **Rlohi_d, double **Rhilo_d, double **Rlolo_d,
   double **urhshihi_h, double **urhslohi_h,
   double **urhshilo_h, double **urhslolo_h,
   double **urhshihi_d, double **urhslohi_d,
   double **urhshilo_d, double **urhslolo_d,
   double **solhihi_h, double **sollohi_h,
   double **solhilo_h, double **sollolo_h,
   double **solhihi_d, double **sollohi_d,
   double **solhilo_d, double **sollolo_d, int vrblvl );
/*
 * DESCRIPTION :
 *   Compares results on host with those on the device
 *   and returns the number of errors larger than 1.0e-8.
 *
 * ON ENTRY :
 *   dim       number of equations and variables;
 *   deg       degree at which the series are truncated;
 *   inputhihi_h are the highest doubles of updated series on the host;
 *   inputlohi_h are the 2nd highest doubles of updated series on the host;
 *   inputhilo_h are the 2nd lowest doubles of updated series on the host;
 *   inputlolo_h are the lowest doubles of updated series on the host;
 *   inputhihi_d are the highest doubles of updated series on the device;
 *   inputlohi_d are the 2nd highest doubles of updated series on the device;
 *   inputhilo_d are the 2nd lowest doubles of updated series on the device;
 *   inputlolo_d are the lowest doubles of updated series on the device;
 *   Qhihi_h   highest doubles of Q on the host;
 *   Qlohi_h   second highest doubles of Q on the host;
 *   Qhilo_h   second lowest doubles of Q on the host;
 *   Qlolo_h   lowest doubles of Q on the host;
 *   Qhihi_d   highest doubles of Q on the device;
 *   Qlohi_d   second  highest doubles of Q on the device;
 *   Qhilo_d   second lowest doubles of Q on the device;
 *   Qlolo_d   lowest doubles of Q on the device;
 *   Rhihi_h   highest doubles of R on the host;
 *   Rlohi_h   second highest doubles of R on the host;
 *   Rhilo_h   second lowest doubles of R on the host;
 *   Rlolo_h   lowest doubles of R on the host;
 *   Rhihi_d   highest doubles of R on the device;
 *   Rlohi_d   second highest doubles of R on the device;
 *   Rhilo_d   second lowest doubles of R on the device;
 *   Rlolo_d   lowest doubles of R on the device;
 *   urhshihi_h are highest doubles of updated right hand side on host;
 *   urhslohi_h are 2nd highest doubles of updated right hand side on host;
 *   urhshilo_h are 2nd lowest doubles of updated right hand side on host;
 *   urhslolo_h are lowest doubles of updated right hand side on host;
 *   urhshihi_d are highest doubles of updated right hand side on device;
 *   urhslohi_d are 2nd highest doubles of updated right hand side on device;
 *   urhshilo_d are 2nd lowest doubles of updated right hand side on device;
 *   urhslolo_d are lowest doubles of updated right hand side on device;
 *   solhihi_h  are highest doubles of update to the solution on host,
 *   sollohi_h  are 2nd highest doubles of update to the solution on host,
 *   solhilo_h  are 2nd lowest doubles of update to the solution on host,
 *   sollolo_h  are lowest doubles of update to the solution on host,
 *   solhihi_d  are highest doubles of update to the solution on device;
 *   sollohi_d  are 2nd highest doubles of update to the solution on device;
 *   solhilo_d  are 2nd lowest doubles of update to the solution on device;
 *   sollolo_d  are lowest doubles of update to the solution on device;
 *   vrblvl    is the verbose level. */

int dbl4_update_newton_qrstep
 ( int szt, int nbt, int dim, int deg,
   int *tailidx_h, int *tailidx_d,
   double **inputhihi_h, double **inputlohi_h,
   double **inputhilo_h, double **inputlolo_h,
   double **inputhihi_d, double **inputlohi_d,
   double **inputhilo_d, double **inputlolo_d,
   double **funvalhihi_h, double **funvallohi_h,
   double **funvalhilo_h, double **funvallolo_h,
   double **funvalhihi_d, double **funvallohi_d,
   double **funvalhilo_d, double **funvallolo_d,
   double ***jacvalhihi_h, double ***jacvallohi_h,
   double ***jacvalhilo_h, double ***jacvallolo_h,
   double ***jacvalhihi_d, double ***jacvallohi_d,
   double ***jacvalhilo_d, double ***jacvallolo_d,
   double **rhshihi_h, double **rhslohi_h,
   double **rhshilo_h, double **rhslolo_h,
   double **rhshihi_d, double **rhslohi_d,
   double **rhshilo_d, double **rhslolo_d,
   double **urhshihi_h, double **urhslohi_h,
   double **urhshilo_h, double **urhslolo_h,
   double **urhshihi_d, double **urhslohi_d,
   double **urhshilo_d, double **urhslolo_d,
   double **solhihi_h, double **sollohi_h,
   double **solhilo_h, double **sollolo_h,
   double **solhihi_d, double **sollohi_d,
   double **solhilo_d, double **sollolo_d,
   double **Qhihi_h, double **Qlohi_h, double **Qhilo_h, double **Qlolo_h,
   double **Qhihi_d, double **Qlohi_d, double **Qhilo_d, double **Qlolo_d,
   double **Rhihi_h, double **Rlohi_h, double **Rhilo_h, double **Rlolo_h,
   double **Rhihi_d, double **Rlohi_d, double **Rhilo_d, double **Rlolo_d,
   double *workvechihi, double *workveclohi,
   double *workvechilo, double *workveclolo,
   double **resvechihi, double **resveclohi, 
   double **resvechilo, double **resveclolo, 
   double *resmaxhihi, double *resmaxlohi,
   double *resmaxhilo, double *resmaxlolo,
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
 *   inputhihi_h has highest doubles of coefficients of the series,
 *             of degree deg, for dim variables, computed on host;
 *   inputlohi_h has second highest doubles of coefficients of the series,
 *             of degree deg, for dim variables, computed on host;
 *   inputhilo_h has second lowest doubles of coefficients of the series,
 *             of degree deg, for dim variables, computed on host;
 *   inputlolo_h has lowest doubles of coefficients of the series,
 *             of degree deg, for dim variables, computed on host;
 *   inputhihi_d has space for series computed on device;
 *   inputlohi_d has space for series computed on device;
 *   inputhilo_d has space for series computed on device;
 *   inputlolo_d has space for series computed on device;
 *   funvalhihi_h has space for the evaluated power series computed by host;
 *   funvallohi_h has space for the evaluated power series computed by host;
 *   funvalhilo_h has space for the evaluated power series computed by host;
 *   funvallolo_h has space for the evaluated power series computed by host;
 *   funvalhihi_d has space for the evaluated power series computed by device;
 *   funvallohi_d has space for the evaluated power series computed by device;
 *   funvalhilo_d has space for the evaluated power series computed by device;
 *   funvallolo_d has space for the evaluated power series computed by device;
 *   jacvalhihi_h has space for deg+1 matrices of dimension dim on host;
 *   jacvallohi_h has space for deg+1 matrices of dimension dim on host;
 *   jacvalhilo_h has space for deg+1 matrices of dimension dim on host;
 *   jacvallolo_h has space for deg+1 matrices of dimension dim on host;
 *   jacvalhihi_d has space for deg+1 matrices of dimension dim on device;
 *   jacvallohi_d has space for deg+1 matrices of dimension dim on device;
 *   jacvalhilo_d has space for deg+1 matrices of dimension dim on device;
 *   jacvallolo_d has space for deg+1 matrices of dimension dim on device;
 *   rhshihi_h   space for deg+1 vectors of dimension dim on host;
 *   rhslohi_h   space for deg+1 vectors of dimension dim on host;
 *   rhshilo_h   space for deg+1 vectors of dimension dim on host;
 *   rhslolo_h   space for deg+1 vectors of dimension dim on host;
 *   rhshihi_d   space for deg+1 vectors of dimension dim on device;
 *   rhslohi_d   space for deg+1 vectors of dimension dim on device;
 *   rhshilo_d   space for deg+1 vectors of dimension dim on device;
 *   rhslolo_d   space for deg+1 vectors of dimension dim on device;
 *   urhshihi_h  space for updated right hand side vectors computed by host;
 *   urhslohi_h  space for updated right hand side vectors computed by host;
 *   urhshilo_h  space for updated right hand side vectors computed by host;
 *   urhslolo_h  space for updated right hand side vectors computed by host;
 *   urhshihi_d  space for updated right hand side vectors computed by device; 
 *   urhslohi_d  space for updated right hand side vectors computed by device; 
 *   urhshilo_d  space for updated right hand side vectors computed by device; 
 *   urhslolo_d  space for updated right hand side vectors computed by device; 
 *   solhihi_h   space for deg+1 vectors of dimension dim;
 *   sollohi_h   space for deg+1 vectors of dimension dim;
 *   solhilo_h   space for deg+1 vectors of dimension dim;
 *   sollolo_h   space for deg+1 vectors of dimension dim;
 *   solhihi_d   space for deg+1 vectors of dimension dim;
 *   sollohi_d   space for deg+1 vectors of dimension dim;
 *   solhilo_d   space for deg+1 vectors of dimension dim;
 *   sollolo_d   space for deg+1 vectors of dimension dim;
 *   Qhihi_h     space allocated for the Q computed by the host;
 *   Qlohi_h     space allocated for the Q computed by the host;
 *   Qhilo_h     space allocated for the Q computed by the host;
 *   Qlolo_h     space allocated for the Q computed by the host;
 *   Qhihi_d     space allocated for the Q computed by the device;
 *   Qlohi_d     space allocated for the Q computed by the device;
 *   Qhilo_d     space allocated for the Q computed by the device;
 *   Qlolo_d     space allocated for the Q computed by the device;
 *   Rhihi_h     space allocated for the R computed by the host;
 *   Rlohi_h     space allocated for the R computed by the host;
 *   Rhilo_h     space allocated for the R computed by the host;
 *   Rlolo_h     space allocated for the R computed by the host;
 *   Rhihi_d     space allocated for the R computed by the device;
 *   Rlohi_d     space allocated for the R computed by the device;
 *   Rhilo_d     space allocated for the R computed by the device;
 *   Rlolo_d     space allocated for the R computed by the device;
 *   wrkvechihi  has work space allocated for a vector of dimension dim;
 *   wrkveclohi  has work space allocated for a vector of dimension dim;
 *   wrkvechilo  has work space allocated for a vector of dimension dim;
 *   wrkveclolo  has work space allocated for a vector of dimension dim;
 *   resvechihi  has space for deg+1 vectors of dimension dim;
 *   resveclohi  has space for deg+1 vectors of dimension dim;
 *   resvechilo  has space for deg+1 vectors of dimension dim;
 *   resveclolo  has space for deg+1 vectors of dimension dim;
 *   zeroQ_h     if true, then Q is zero and Q must be computed on host;
 *   noqr_h      flag if true, then no qr on host;
 *   zeroQ_d     if true, then Q is zero and Q must be computed on device;
 *   noqr_d      flag if true, then no qr on device;
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
 *   inputhihi_h has highest doubles of series computed on host (mode 1, 2);
 *   inputlohi_h has second highest doubles of series computed on host;
 *   inputhilo_h has second lowest doubles of series computed on host;
 *   inputlolo_h has lowest doubles of series computed on host;
 *   inputhihi_d has highest doubles of series, computed by device (mode 0, 2);
 *   inputlohi_d has second highest doubles of series, computed by device;
 *   inputhilo_d has second lowest doubles of series, computed by device;
 *   inputlolo_d has lowest doubles of series, computed by device;
 *   outputhihi_h has highest doubles of evaluated series computed on host;
 *   outputlohi_h has second highest doubles of series computed on host;
 *   outputhilo_h has second lowest doubles of series computed on host;
 *   outputlolo_h has lowest doubles of series computed on host;
 *   outputhihi_d has highest doubles of series computed on device;
 *   outputlohi_d has second highest doubles of series computed on device;
 *   outputhilo_d has second lowest doubles of series computed on device;
 *   outputlolo_d has lowest doubles of series computed on device;
 *   funvalhihi_h has highest doubles of output[i][dim], values on host;
 *   funvallohi_h has second highest doubles of output[i][dim], values on host;
 *   funvalhilo_h has second lowest doubles of output[i][dim], values on host;
 *   funvallolo_h has lowest doubles of output[i][dim], values on host;
 *   funvalhihi_d has highest doubles of output[i][dim], values on device;
 *   funvallohi_d has second highest doubles of output[i][dim], on device;
 *   funvalhilo_d has second lowest doubles of output[i][dim], on device;
 *   funvallolo_d has lowest doubles of output[i][dim], values on device;
 *   jacvalhihi_h has highest doubles of matrices, computed by host,
 *             the leading coefficient is the Jacobian matrix;
 *   jacvallohi_h has second highest doubles of matrices, computed by host,
 *             the leading coefficient is the Jacobian matrix;
 *   jacvalhilo_h has second lowest doubles of matrices, computed by host,
 *             the leading coefficient is the Jacobian matrix;
 *   jacvallolo_h has lowest doubles of matrices, computed by host,
 *             the leading coefficient is the Jacobian matrix;
 *   jacvalhihi_d has highest doubles of matrices, computed by device,
 *             the leading coefficient is the Jacobian matrix;
 *   jacvallohi_d has second highest doubles of matrices, computed by device,
 *             the leading coefficient is the Jacobian matrix;
 *   jacvalhilo_d has second lowest doubles of matrices, computed by device,
 *             the leading coefficient is the Jacobian matrix;
 *   jacvallolo_d has lowest doubles of matrices, computed by device,
 *             the leading coefficient is the Jacobian matrix;
 *   rhshihi_h has highest doubles of the right hand side, linearized,
 *             subtracted by 1 and added by t, computed by host;
 *   rhslohi_h has second highest doubles of the right hand side, linearized,
 *             subtracted by 1 and added by t, computed by host;
 *   rhshilo_h has second lowest doubles of the right hand side, linearized,
 *             subtracted by 1 and added by t, computed by host;
 *   rhslolo_h has lowest doubles of the right hand side, linearized,
 *             subtracted by 1 and added by t, computed by host;
 *   rhshihi_d has highest doubles of the right hand side, linearized,
 *             subtracted by 1 and added by t, computed by device;
 *   rhslohi_d has second highest doubles of the right hand side, linearized,
 *             subtracted by 1 and added by t, computed by device;
 *   rhshilo_d has second lowest doubles of the right hand side, linearized,
 *             subtracted by 1 and added by t, computed by device;
 *   rhslolo_d has lowest doubles of the right hand side, linearized,
 *             subtracted by 1 and added by t, computed by device;
 *   urhshihi_h has highest doubles of right hand side, on host;
 *   urhslohi_h has second highest doubles of right hand side, on host;
 *   urhshilo_h has second lowest doubles of right hand side, on host;
 *   urhslolo_h has lowest doubles of right hand side, on host;
 *   urhshihi_d has highest doubles of right hand side, on device;
 *   urhslohi_d has second highest doubles of right hand side, on device;
 *   urhshilo_d has second lowest doubles of right hand side, on device;
 *   urhslolo_d has lowest doubles of right hand side, on device;
 *   solhihi_h has highest doubles of solution computed by the host;
 *   sollohi_h has second highest doubles of solution computed by the host;
 *   solhilo_h has second lowest doubles of solution computed by the host;
 *   sollolo_h has lowest doubles of solution computed by the host;
 *   solhihi_d has highest doubles of solution computed by the device;
 *   sollohi_d has second highest doubles of solution computed by the device;
 *   solhilo_d has second lowest doubles of solution computed by the device;
 *   sollolo_d has lowest doubles of solution computed by the device;
 *   Qhihi_h   highest doubles of Q of the QR computed by the host;
 *   Qlohi_h   second highest doubles of Q of the QR computed by the host;
 *   Qhilo_h   second lowest doubles of Q of the QR computed by the host;
 *   Qlolo_h   lowest doubles of Q of the QR computed by the host;
 *   Qhihi_d   highest doubles of Q of the QR computed by the device;
 *   Qlohi_d   second highest doubles of Q of the QR computed by the device;
 *   Qhilo_d   second lowest doubles of Q of the QR computed by the device;
 *   Qlolo_d   lowest doubles of Q of the QR computed by the device;
 *   Rhihi_h   highest doubles of R of the QR computed by the host;
 *   Rlohi_h   second highest doubles of R of the QR computed by the host;
 *   Rhilo_h   second lowest doubles of R of the QR computed by the host;
 *   Rlolo_h   lowest doubles of R of the QR computed by the host;
 *   Rhihi_d   highest doubles of R of the QR computed by the device;
 *   Rlohi_d   second highest doubles of R of the QR computed by the device;
 *   Rhilo_d   second lowest doubles of R of the QR computed by the device;
 *   Rlolo_d   lowest doubles of R of the QR computed by the device;
 *   wrkmathihi has a copy of the highest doubles of the Jacobian matrix;
 *   wrkmatlohi has a copy of the second highest doubles of the Jacobian;
 *   wrkmathilo has a copy of the second lowest doubles of the Jacobian;
 *   wrkmatlolo has a copy of the lowest doubles of the Jacobian matrix;
 *   resvechihi has highest doubles of the residual vectors;
 *   resveclohi has second highest doubles of the residual vectors;
 *   resvechilo has second lowest doubles of the residual vectors;
 *   resveclolo has lowest doubles of the residual vectors;
 *   resmaxhihi is highest double of the maximum of residual vectors;
 *   resmaxlohi is second highest double of the maximum of residual vectors;
 *   resmaxhilo is second lowest double of the maximum of residual vectors;
 *   resmaxlolo is lowest double of the maximum of residual vectors;
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

int dbl4_column_newton_qrstep
 ( int szt, int nbt, int dim, int deg, int nbrcol,
   int *tailidx_h, int *tailidx_d,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac, 
   double **mbhihi, double **mblohi, double **mbhilo, double **mblolo,
   double dpr,
   double ***cffhihi, double ***cfflohi, double ***cffhilo, double ***cfflolo,
   double **acchihi, double **acclohi, double **acchilo, double **acclolo,
   double **inputhihi_h, double **inputlohi_h,
   double **inputhilo_h, double **inputlolo_h,
   double **inputhihi_d, double **inputlohi_d,
   double **inputhilo_d, double **inputlolo_d,
   double ***outputhihi_h, double ***outputlohi_h,
   double ***outputhilo_h, double ***outputlolo_h,
   double ***outputhihi_d, double ***outputlohi_d,
   double ***outputhilo_d, double ***outputlolo_d,
   double **funvalhihi_h, double **funvallohi_h,
   double **funvalhilo_h, double **funvallolo_h,
   double **funvalhihi_d, double **funvallohi_d,
   double **funvalhilo_d, double **funvallolo_d,
   double ***jacvalhihi_h, double ***jacvallohi_h,
   double ***jacvalhilo_h, double ***jacvallolo_h,
   double ***jacvalhihi_d, double ***jacvallohi_d,
   double ***jacvalhilo_d, double ***jacvallolo_d,
   double **rhshihi_h, double **rhslohi_h,
   double **rhshilo_h, double **rhslolo_h,
   double **rhshihi_d, double **rhslohi_d,
   double **rhshilo_d, double **rhslolo_d,
   double **urhshihi_h, double **urhslohi_h,
   double **urhshilo_h, double **urhslolo_h,
   double **urhshihi_d, double **urhslohi_d,
   double **urhshilo_d, double **urhslolo_d,
   double **solhihi_h, double **sollohi_h,
   double **solhilo_h, double **sollolo_h,
   double **solhihi_d, double **sollohi_d,
   double **solhilo_d, double **sollolo_d,
   double **Qhihi_h, double **Qlohi_h, double **Qhilo_h, double **Qlolo_h,
   double **Qhihi_d, double **Qlohi_d, double **Qhilo_d, double **Qlolo_d,
   double **Rhihi_h, double **Rlohi_h, double **Rhilo_h, double **Rlolo_h,
   double **Rhihi_d, double **Rlohi_d, double **Rhilo_d, double **Rlolo_d,
   double *workvechihi, double *workveclohi,
   double *workvechilo, double *workveclolo,
   double **resvechihi, double **resveclohi, 
   double **resvechilo, double **resveclolo, 
   double *resmaxhihi, double *resmaxlohi,
   double *resmaxhilo, double *resmaxlolo,
   bool *zeroQ_h, bool *noqr_h, bool *zeroQ_d, bool *noqr_d,
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
 *   nbrcol    number of columns, if 1, then the system is monomial,
 *             otherwise nbrcol columns are expected;
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
 *   mbhihi    highest doubles of the right hand side of monomial system;
 *   mblohi    second highest doubles of the right hand side;
 *   mbhilo    second lowest doubles of the right hand side;
 *   mblolo    lowest doubles of the right hand side;
 *   dpr       damper multiplier for t, should be in (0.0, 1.0];
 *   cffhihi   highest doubles of the coefficients of the monomials;
 *   cfflohi   second highest doubles of the coefficients of the monomials;
 *   cffhilo   second lowest doubles of the coefficients of the monomials;
 *   cfflolo   lowest doubles of the coefficients of the monomials;
 *   acchihi   space for highest doubles of one series of degree deg;
 *   acclohi   space for second highest doubles of one series of degree deg;
 *   acchilo   space for second lowest doubles of one series of degree deg;
 *   acclolo   space for lowest doubles of one series of degree deg;
 *   inputhihi_h has highest doubles of coefficients of the series,
 *             of degree deg, for dim variables, computed on host;
 *   inputlohi_h has second highest doubles of coefficients of the series,
 *             of degree deg, for dim variables, computed on host;
 *   inputhilo_h has second lowest doubles of coefficients of the series,
 *             of degree deg, for dim variables, computed on host;
 *   inputlolo_h has lowest doubles of coefficients of the series,
 *             of degree deg, for dim variables, computed on host;
 *   inputhihi_d has space for series computed on device;
 *   inputlohi_d has space for series computed on device;
 *   inputhilo_d has space for series computed on device;
 *   inputlolo_d has space for series computed on device;
 *   outputhihi_h has space for the highest doubles of the evaluated
 *             and differentiated monomials, computed on the host;
 *   outputlohi_h has space for the second highest doubles of the evaluated 
 *             and differentiated monomials, computed on the host;
 *   outputhilo_h has space for the second lowest doubles of the evaluated
 *             and differentiated monomials, computed on the host;
 *   outputlolo_h has space for the lowest doubles of the evaluated 
 *             and differentiated monomials, computed on the host;
 *   outputhihi_d has space for the highest doubles of the evaluated and
 *             differentiated monomials, computed on the device;
 *   outputlohi_d has space for the second highest doubles of the evaluated
 *             and differentiated monomials, computed on the device;
 *   outputhilo_d has space for the second lowest doubles of the evaluated
 *             and differentiated monomials, computed on the device;
 *   outputlolo_d has space for the lowest doubles of the evaluated and
 *             differentiated monomials, computed on the device;
 *   funvalhihi_h has space for the evaluated power series computed by host;
 *   funvallohi_h has space for the evaluated power series computed by host;
 *   funvalhilo_h has space for the evaluated power series computed by host;
 *   funvallolo_h has space for the evaluated power series computed by host;
 *   funvalhihi_d has space for the evaluated power series computed by device;
 *   funvallohi_d has space for the evaluated power series computed by device;
 *   funvalhilo_d has space for the evaluated power series computed by device;
 *   funvallolo_d has space for the evaluated power series computed by device;
 *   jacvalhihi_h has space for deg+1 matrices of dimension dim on host;
 *   jacvallohi_h has space for deg+1 matrices of dimension dim on host;
 *   jacvalhilo_h has space for deg+1 matrices of dimension dim on host;
 *   jacvallolo_h has space for deg+1 matrices of dimension dim on host;
 *   jacvalhihi_d has space for deg+1 matrices of dimension dim on device;
 *   jacvallohi_d has space for deg+1 matrices of dimension dim on device;
 *   jacvalhilo_d has space for deg+1 matrices of dimension dim on device;
 *   jacvallolo_d has space for deg+1 matrices of dimension dim on device;
 *   rhshihi_h   space for deg+1 vectors of dimension dim on host;
 *   rhslohi_h   space for deg+1 vectors of dimension dim on host;
 *   rhshilo_h   space for deg+1 vectors of dimension dim on host;
 *   rhslolo_h   space for deg+1 vectors of dimension dim on host;
 *   rhshihi_d   space for deg+1 vectors of dimension dim on device;
 *   rhslohi_d   space for deg+1 vectors of dimension dim on device;
 *   rhshilo_d   space for deg+1 vectors of dimension dim on device;
 *   rhslolo_d   space for deg+1 vectors of dimension dim on device;
 *   urhshihi_h  space for updated right hand side vectors computed by host;
 *   urhslohi_h  space for updated right hand side vectors computed by host;
 *   urhshilo_h  space for updated right hand side vectors computed by host;
 *   urhslolo_h  space for updated right hand side vectors computed by host;
 *   urhshihi_d  space for updated right hand side vectors computed by device; 
 *   urhslohi_d  space for updated right hand side vectors computed by device; 
 *   urhshilo_d  space for updated right hand side vectors computed by device; 
 *   urhslolo_d  space for updated right hand side vectors computed by device; 
 *   solhihi_h   space for deg+1 vectors of dimension dim;
 *   sollohi_h   space for deg+1 vectors of dimension dim;
 *   solhilo_h   space for deg+1 vectors of dimension dim;
 *   sollolo_h   space for deg+1 vectors of dimension dim;
 *   solhihi_d   space for deg+1 vectors of dimension dim;
 *   sollohi_d   space for deg+1 vectors of dimension dim;
 *   solhilo_d   space for deg+1 vectors of dimension dim;
 *   sollolo_d   space for deg+1 vectors of dimension dim;
 *   Qhihi_h     space allocated for the Q computed by the host;
 *   Qlohi_h     space allocated for the Q computed by the host;
 *   Qhilo_h     space allocated for the Q computed by the host;
 *   Qlolo_h     space allocated for the Q computed by the host;
 *   Qhihi_d     space allocated for the Q computed by the device;
 *   Qlohi_d     space allocated for the Q computed by the device;
 *   Qhilo_d     space allocated for the Q computed by the device;
 *   Qlolo_d     space allocated for the Q computed by the device;
 *   Rhihi_h     space allocated for the R computed by the host;
 *   Rlohi_h     space allocated for the R computed by the host;
 *   Rhilo_h     space allocated for the R computed by the host;
 *   Rlolo_h     space allocated for the R computed by the host;
 *   Rhihi_d     space allocated for the R computed by the device;
 *   Rlohi_d     space allocated for the R computed by the device;
 *   Rhilo_d     space allocated for the R computed by the device;
 *   Rlolo_d     space allocated for the R computed by the device;
 *   wrkmathihi  has work space allocated for a matrix of dimension dim;
 *   wrkmatlohi  has work space allocated for a matrix of dimension dim;
 *   wrkmathilo  has work space allocated for a matrix of dimension dim;
 *   wrkmatlolo  has work space allocated for a matrix of dimension dim;
 *   wrkvechihi  has work space allocated for a vector of dimension dim;
 *   wrkveclohi  has work space allocated for a vector of dimension dim;
 *   wrkvechilo  has work space allocated for a vector of dimension dim;
 *   wrkveclolo  has work space allocated for a vector of dimension dim;
 *   resvechihi  has space for deg+1 vectors of dimension dim;
 *   resveclohi  has space for deg+1 vectors of dimension dim;
 *   resvechilo  has space for deg+1 vectors of dimension dim;
 *   resveclolo  has space for deg+1 vectors of dimension dim;
 *   zeroQ_h     if true, then Q is zero and Q must be computed on host;
 *   noqr_h      flag if true, then no qr on host;
 *   zeroQ_d     if true, then Q is zero and Q must be computed on device;
 *   noqr_d      flag if true, then no qr on device;
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
 *   inputhihi_h has highest doubles of series computed on host (mode 1, 2);
 *   inputlohi_h has second highest doubles of series computed on host;
 *   inputhilo_h has second lowest doubles of series computed on host;
 *   inputlolo_h has lowest doubles of series computed on host;
 *   inputhihi_d has highest doubles of series, computed by device (mode 0, 2);
 *   inputlohi_d has second highest doubles of series, computed by device;
 *   inputhilo_d has second lowest doubles of series, computed by device;
 *   inputlolo_d has lowest doubles of series, computed by device;
 *   outputhihi_h has highest doubles of evaluated series computed on host;
 *   outputlohi_h has second highest doubles of series computed on host;
 *   outputhilo_h has second lowest doubles of series computed on host;
 *   outputlolo_h has lowest doubles of series computed on host;
 *   outputhihi_d has highest doubles of series computed on device;
 *   outputlohi_d has second highest doubles of series computed on device;
 *   outputhilo_d has second lowest doubles of series computed on device;
 *   outputlolo_d has lowest doubles of series computed on device;
 *   funvalhihi_h has highest doubles of output[i][dim], values on host;
 *   funvallohi_h has second highest doubles of output[i][dim], values on host;
 *   funvalhilo_h has second lowest doubles of output[i][dim], values on host;
 *   funvallolo_h has lowest doubles of output[i][dim], values on host;
 *   funvalhihi_d has highest doubles of output[i][dim], values on device;
 *   funvallohi_d has second highest doubles of output[i][dim], on device;
 *   funvalhilo_d has second lowest doubles of output[i][dim], on device;
 *   funvallolo_d has lowest doubles of output[i][dim], values on device;
 *   jacvalhihi_h has highest doubles of matrices, computed by host,
 *             the leading coefficient is the Jacobian matrix;
 *   jacvallohi_h has second highest doubles of matrices, computed by host,
 *             the leading coefficient is the Jacobian matrix;
 *   jacvalhilo_h has second lowest doubles of matrices, computed by host,
 *             the leading coefficient is the Jacobian matrix;
 *   jacvallolo_h has lowest doubles of matrices, computed by host,
 *             the leading coefficient is the Jacobian matrix;
 *   jacvalhihi_d has highest doubles of matrices, computed by device,
 *             the leading coefficient is the Jacobian matrix;
 *   jacvallohi_d has second highest doubles of matrices, computed by device,
 *             the leading coefficient is the Jacobian matrix;
 *   jacvalhilo_d has second lowest doubles of matrices, computed by device,
 *             the leading coefficient is the Jacobian matrix;
 *   jacvallolo_d has lowest doubles of matrices, computed by device,
 *             the leading coefficient is the Jacobian matrix;
 *   rhshihi_h has highest doubles of the right hand side, linearized,
 *             subtracted by 1 and added by t, computed by host;
 *   rhslohi_h has second highest doubles of the right hand side, linearized,
 *             subtracted by 1 and added by t, computed by host;
 *   rhshilo_h has second lowest doubles of the right hand side, linearized,
 *             subtracted by 1 and added by t, computed by host;
 *   rhslolo_h has lowest doubles of the right hand side, linearized,
 *             subtracted by 1 and added by t, computed by host;
 *   rhshihi_d has highest doubles of the right hand side, linearized,
 *             subtracted by 1 and added by t, computed by device;
 *   rhslohi_d has second highest doubles of the right hand side, linearized,
 *             subtracted by 1 and added by t, computed by device;
 *   rhshilo_d has second lowest doubles of the right hand side, linearized,
 *             subtracted by 1 and added by t, computed by device;
 *   rhslolo_d has lowest doubles of the right hand side, linearized,
 *             subtracted by 1 and added by t, computed by device;
 *   urhshihi_h has highest doubles of right hand side, on host;
 *   urhslohi_h has second highest doubles of right hand side, on host;
 *   urhshilo_h has second lowest doubles of right hand side, on host;
 *   urhslolo_h has lowest doubles of right hand side, on host;
 *   urhshihi_d has highest doubles of right hand side, on device;
 *   urhslohi_d has second highest doubles of right hand side, on device;
 *   urhshilo_d has second lowest doubles of right hand side, on device;
 *   urhslolo_d has lowest doubles of right hand side, on device;
 *   solhihi_h has highest doubles of solution computed by the host;
 *   sollohi_h has second highest doubles of solution computed by the host;
 *   solhilo_h has second lowest doubles of solution computed by the host;
 *   sollolo_h has lowest doubles of solution computed by the host;
 *   solhihi_d has highest doubles of solution computed by the device;
 *   sollohi_d has second highest doubles of solution computed by the device;
 *   solhilo_d has second lowest doubles of solution computed by the device;
 *   sollolo_d has lowest doubles of solution computed by the device;
 *   Qhihi_h   highest doubles of Q of the QR computed by the host;
 *   Qlohi_h   second highest doubles of Q of the QR computed by the host;
 *   Qhilo_h   second lowest doubles of Q of the QR computed by the host;
 *   Qlolo_h   lowest doubles of Q of the QR computed by the host;
 *   Qhihi_d   highest doubles of Q of the QR computed by the device;
 *   Qlohi_d   second highest doubles of Q of the QR computed by the device;
 *   Qhilo_d   second lowest doubles of Q of the QR computed by the device;
 *   Qlolo_d   lowest doubles of Q of the QR computed by the device;
 *   Rhihi_h   highest doubles of R of the QR computed by the host;
 *   Rlohi_h   second highest doubles of R of the QR computed by the host;
 *   Rhilo_h   second lowest doubles of R of the QR computed by the host;
 *   Rlolo_h   lowest doubles of R of the QR computed by the host;
 *   Rhihi_d   highest doubles of R of the QR computed by the device;
 *   Rlohi_d   second highest doubles of R of the QR computed by the device;
 *   Rhilo_d   second lowest doubles of R of the QR computed by the device;
 *   Rlolo_d   lowest doubles of R of the QR computed by the device;
 *   wrkmathihi has a copy of the highest doubles of the Jacobian matrix;
 *   wrkmatlohi has a copy of the second highest doubles of the Jacobian;
 *   wrkmathilo has a copy of the second lowest doubles of the Jacobian;
 *   wrkmatlolo has a copy of the lowest doubles of the Jacobian matrix;
 *   resvechihi has highest doubles of the residual vectors;
 *   resveclohi has second highest doubles of the residual vectors;
 *   resvechilo has second lowest doubles of the residual vectors;
 *   resveclolo has lowest doubles of the residual vectors;
 *   resmaxhihi is highest double of the maximum of residual vectors;
 *   resmaxlohi is second highest double of the maximum of residual vectors;
 *   resmaxhilo is second lowest double of the maximum of residual vectors;
 *   resmaxlolo is lowest double of the maximum of residual vectors;
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

int dbl4_row_newton_qrstep
 ( int szt, int nbt, int dim, int deg, int *tailidx_h, int *tailidx_d,
   int *nbr, int **nvr, int ***idx,
   double **csthihi, double **cstlohi, double **csthilo, double **cstlolo,
   double ***cffhihi, double ***cfflohi, double ***cffhilo, double ***cfflolo,
   double dpr,
   double **inputhihi_h, double **inputlohi_h,
   double **inputhilo_h, double **inputlolo_h,
   double **inputhihi_d, double **inputlohi_d,
   double **inputhilo_d, double **inputlolo_d,
   double ***outputhihi_h, double ***outputlohi_h,
   double ***outputhilo_h, double ***outputlolo_h,
   double ***outputhihi_d, double ***outputlohi_d,
   double ***outputhilo_d, double ***outputlolo_d,
   double **funvalhihi_h, double **funvallohi_h,
   double **funvalhilo_h, double **funvallolo_h,
   double **funvalhihi_d, double **funvallohi_d,
   double **funvalhilo_d, double **funvallolo_d,
   double ***jacvalhihi_h, double ***jacvallohi_h,
   double ***jacvalhilo_h, double ***jacvallolo_h,
   double ***jacvalhihi_d, double ***jacvallohi_d,
   double ***jacvalhilo_d, double ***jacvallolo_d,
   double **rhshihi_h, double **rhslohi_h,
   double **rhshilo_h, double **rhslolo_h,
   double **rhshihi_d, double **rhslohi_d,
   double **rhshilo_d, double **rhslolo_d,
   double **urhshihi_h, double **urhslohi_h,
   double **urhshilo_h, double **urhslolo_h,
   double **urhshihi_d, double **urhslohi_d,
   double **urhshilo_d, double **urhslolo_d,
   double **solhihi_h, double **sollohi_h,
   double **solhilo_h, double **sollolo_h,
   double **solhihi_d, double **sollohi_d,
   double **solhilo_d, double **sollolo_d,
   double **Qhihi_h, double **Qlohi_h, double **Qhilo_h, double **Qlolo_h,
   double **Qhihi_d, double **Qlohi_d, double **Qhilo_d, double **Qlolo_d,
   double **Rhihi_h, double **Rlohi_h, double **Rhilo_h, double **Rlolo_h,
   double **Rhihi_d, double **Rlohi_d, double **Rhilo_d, double **Rlolo_d,
   double *workvechihi, double *workveclohi,
   double *workvechilo, double *workveclolo,
   double **resvechihi, double **resveclohi, 
   double **resvechilo, double **resveclolo, 
   double *resmaxhihi, double *resmaxlohi,
   double *resmaxhilo, double *resmaxlolo,
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
 *   nbr       nbr[i] is the number of terms in the i-th polynomial;
 *   nvr       nvr[i][j] is the number of variables in the j-th monomial
 *             of the i-th polynomial;
 *   idx       idx[i][j] are the indices of the variables in monomial j
 *             of the i-th polynomial;
 *   csthihi   highest doubles of the constants;
 *   cstlohi   second highest doubles of the constants;
 *   csthilo   second lowest doubles of the constants;
 *   cstlolo   lowest doubles of the constants;
 *   cffhihi   highest doubles of the coefficients of the monomials;
 *   cfflohi   second highest doubles of the coefficients of the monomials;
 *   cffhilo   second lowest doubles of the coefficients of the monomials;
 *   cfflolo   lowest doubles of the coefficients of the monomials;
 *   dpr       damper multiplier for t, should be in (0.0, 1.0];
 *   inputhihi_h has highest doubles of coefficients of the series,
 *             of degree deg, for dim variables, computed on host;
 *   inputlohi_h has second highest doubles of coefficients of the series,
 *             of degree deg, for dim variables, computed on host;
 *   inputhilo_h has second lowest doubles of coefficients of the series,
 *             of degree deg, for dim variables, computed on host;
 *   inputlolo_h has lowest doubles of coefficients of the series,
 *             of degree deg, for dim variables, computed on host;
 *   inputhihi_d has space for series computed on device;
 *   inputlohi_d has space for series computed on device;
 *   inputhilo_d has space for series computed on device;
 *   inputlolo_d has space for series computed on device;
 *   outputhihi_h has space for the highest doubles of the evaluated
 *             and differentiated monomials, computed on the host;
 *   outputlohi_h has space for the second highest doubles of the evaluated 
 *             and differentiated monomials, computed on the host;
 *   outputhilo_h has space for the second lowest doubles of the evaluated
 *             and differentiated monomials, computed on the host;
 *   outputlolo_h has space for the lowest doubles of the evaluated 
 *             and differentiated monomials, computed on the host;
 *   outputhihi_d has space for the highest doubles of the evaluated and
 *             differentiated monomials, computed on the device;
 *   outputlohi_d has space for the second highest doubles of the evaluated
 *             and differentiated monomials, computed on the device;
 *   outputhilo_d has space for the second lowest doubles of the evaluated
 *             and differentiated monomials, computed on the device;
 *   outputlolo_d has space for the lowest doubles of the evaluated and
 *             differentiated monomials, computed on the device;
 *   funvalhihi_h has space for the evaluated power series computed by host;
 *   funvallohi_h has space for the evaluated power series computed by host;
 *   funvalhilo_h has space for the evaluated power series computed by host;
 *   funvallolo_h has space for the evaluated power series computed by host;
 *   funvalhihi_d has space for the evaluated power series computed by device;
 *   funvallohi_d has space for the evaluated power series computed by device;
 *   funvalhilo_d has space for the evaluated power series computed by device;
 *   funvallolo_d has space for the evaluated power series computed by device;
 *   jacvalhihi_h has space for deg+1 matrices of dimension dim on host;
 *   jacvallohi_h has space for deg+1 matrices of dimension dim on host;
 *   jacvalhilo_h has space for deg+1 matrices of dimension dim on host;
 *   jacvallolo_h has space for deg+1 matrices of dimension dim on host;
 *   jacvalhihi_d has space for deg+1 matrices of dimension dim on device;
 *   jacvallohi_d has space for deg+1 matrices of dimension dim on device;
 *   jacvalhilo_d has space for deg+1 matrices of dimension dim on device;
 *   jacvallolo_d has space for deg+1 matrices of dimension dim on device;
 *   rhshihi_h   space for deg+1 vectors of dimension dim on host;
 *   rhslohi_h   space for deg+1 vectors of dimension dim on host;
 *   rhshilo_h   space for deg+1 vectors of dimension dim on host;
 *   rhslolo_h   space for deg+1 vectors of dimension dim on host;
 *   rhshihi_d   space for deg+1 vectors of dimension dim on device;
 *   rhslohi_d   space for deg+1 vectors of dimension dim on device;
 *   rhshilo_d   space for deg+1 vectors of dimension dim on device;
 *   rhslolo_d   space for deg+1 vectors of dimension dim on device;
 *   urhshihi_h  space for updated right hand side vectors computed by host;
 *   urhslohi_h  space for updated right hand side vectors computed by host;
 *   urhshilo_h  space for updated right hand side vectors computed by host;
 *   urhslolo_h  space for updated right hand side vectors computed by host;
 *   urhshihi_d  space for updated right hand side vectors computed by device; 
 *   urhslohi_d  space for updated right hand side vectors computed by device; 
 *   urhshilo_d  space for updated right hand side vectors computed by device; 
 *   urhslolo_d  space for updated right hand side vectors computed by device; 
 *   solhihi_h   space for deg+1 vectors of dimension dim;
 *   sollohi_h   space for deg+1 vectors of dimension dim;
 *   solhilo_h   space for deg+1 vectors of dimension dim;
 *   sollolo_h   space for deg+1 vectors of dimension dim;
 *   solhihi_d   space for deg+1 vectors of dimension dim;
 *   sollohi_d   space for deg+1 vectors of dimension dim;
 *   solhilo_d   space for deg+1 vectors of dimension dim;
 *   sollolo_d   space for deg+1 vectors of dimension dim;
 *   Qhihi_h     space allocated for the Q computed by the host;
 *   Qlohi_h     space allocated for the Q computed by the host;
 *   Qhilo_h     space allocated for the Q computed by the host;
 *   Qlolo_h     space allocated for the Q computed by the host;
 *   Qhihi_d     space allocated for the Q computed by the device;
 *   Qlohi_d     space allocated for the Q computed by the device;
 *   Qhilo_d     space allocated for the Q computed by the device;
 *   Qlolo_d     space allocated for the Q computed by the device;
 *   Rhihi_h     space allocated for the R computed by the host;
 *   Rlohi_h     space allocated for the R computed by the host;
 *   Rhilo_h     space allocated for the R computed by the host;
 *   Rlolo_h     space allocated for the R computed by the host;
 *   Rhihi_d     space allocated for the R computed by the device;
 *   Rlohi_d     space allocated for the R computed by the device;
 *   Rhilo_d     space allocated for the R computed by the device;
 *   Rlolo_d     space allocated for the R computed by the device;
 *   wrkmathihi  has work space allocated for a matrix of dimension dim;
 *   wrkmatlohi  has work space allocated for a matrix of dimension dim;
 *   wrkmathilo  has work space allocated for a matrix of dimension dim;
 *   wrkmatlolo  has work space allocated for a matrix of dimension dim;
 *   wrkvechihi  has work space allocated for a vector of dimension dim;
 *   wrkveclohi  has work space allocated for a vector of dimension dim;
 *   wrkvechilo  has work space allocated for a vector of dimension dim;
 *   wrkveclolo  has work space allocated for a vector of dimension dim;
 *   resvechihi  has space for deg+1 vectors of dimension dim;
 *   resveclohi  has space for deg+1 vectors of dimension dim;
 *   resvechilo  has space for deg+1 vectors of dimension dim;
 *   resveclolo  has space for deg+1 vectors of dimension dim;
 *   zeroQ_h     if true, then Q is zero and Q must be computed on host;
 *   noqr_h      flag if true, then no qr on host;
 *   zeroQ_d     if true, then Q is zero and Q must be computed on device;
 *   noqr_d      flag if true, then no qr on device;
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
 *   inputhihi_h has highest doubles of series computed on host (mode 1, 2);
 *   inputlohi_h has second highest doubles of series computed on host;
 *   inputhilo_h has second lowest doubles of series computed on host;
 *   inputlolo_h has lowest doubles of series computed on host;
 *   inputhihi_d has highest doubles of series, computed by device (mode 0, 2);
 *   inputlohi_d has second highest doubles of series, computed by device;
 *   inputhilo_d has second lowest doubles of series, computed by device;
 *   inputlolo_d has lowest doubles of series, computed by device;
 *   outputhihi_h has highest doubles of evaluated series computed on host;
 *   outputlohi_h has second highest doubles of series computed on host;
 *   outputhilo_h has second lowest doubles of series computed on host;
 *   outputlolo_h has lowest doubles of series computed on host;
 *   outputhihi_d has highest doubles of series computed on device;
 *   outputlohi_d has second highest doubles of series computed on device;
 *   outputhilo_d has second lowest doubles of series computed on device;
 *   outputlolo_d has lowest doubles of series computed on device;
 *   funvalhihi_h has highest doubles of output[i][dim], values on host;
 *   funvallohi_h has second highest doubles of output[i][dim], values on host;
 *   funvalhilo_h has second lowest doubles of output[i][dim], values on host;
 *   funvallolo_h has lowest doubles of output[i][dim], values on host;
 *   funvalhihi_d has highest doubles of output[i][dim], values on device;
 *   funvallohi_d has second highest doubles of output[i][dim], on device;
 *   funvalhilo_d has second lowest doubles of output[i][dim], on device;
 *   funvallolo_d has lowest doubles of output[i][dim], values on device;
 *   jacvalhihi_h has highest doubles of matrices, computed by host,
 *             the leading coefficient is the Jacobian matrix;
 *   jacvallohi_h has second highest doubles of matrices, computed by host,
 *             the leading coefficient is the Jacobian matrix;
 *   jacvalhilo_h has second lowest doubles of matrices, computed by host,
 *             the leading coefficient is the Jacobian matrix;
 *   jacvallolo_h has lowest doubles of matrices, computed by host,
 *             the leading coefficient is the Jacobian matrix;
 *   jacvalhihi_d has highest doubles of matrices, computed by device,
 *             the leading coefficient is the Jacobian matrix;
 *   jacvallohi_d has second highest doubles of matrices, computed by device,
 *             the leading coefficient is the Jacobian matrix;
 *   jacvalhilo_d has second lowest doubles of matrices, computed by device,
 *             the leading coefficient is the Jacobian matrix;
 *   jacvallolo_d has lowest doubles of matrices, computed by device,
 *             the leading coefficient is the Jacobian matrix;
 *   rhshihi_h has highest doubles of the right hand side, linearized,
 *             subtracted by 1 and added by t, computed by host;
 *   rhslohi_h has second highest doubles of the right hand side, linearized,
 *             subtracted by 1 and added by t, computed by host;
 *   rhshilo_h has second lowest doubles of the right hand side, linearized,
 *             subtracted by 1 and added by t, computed by host;
 *   rhslolo_h has lowest doubles of the right hand side, linearized,
 *             subtracted by 1 and added by t, computed by host;
 *   rhshihi_d has highest doubles of the right hand side, linearized,
 *             subtracted by 1 and added by t, computed by device;
 *   rhslohi_d has second highest doubles of the right hand side, linearized,
 *             subtracted by 1 and added by t, computed by device;
 *   rhshilo_d has second lowest doubles of the right hand side, linearized,
 *             subtracted by 1 and added by t, computed by device;
 *   rhslolo_d has lowest doubles of the right hand side, linearized,
 *             subtracted by 1 and added by t, computed by device;
 *   urhshihi_h has highest doubles of right hand side, on host;
 *   urhslohi_h has second highest doubles of right hand side, on host;
 *   urhshilo_h has second lowest doubles of right hand side, on host;
 *   urhslolo_h has lowest doubles of right hand side, on host;
 *   urhshihi_d has highest doubles of right hand side, on device;
 *   urhslohi_d has second highest doubles of right hand side, on device;
 *   urhshilo_d has second lowest doubles of right hand side, on device;
 *   urhslolo_d has lowest doubles of right hand side, on device;
 *   solhihi_h has highest doubles of solution computed by the host;
 *   sollohi_h has second highest doubles of solution computed by the host;
 *   solhilo_h has second lowest doubles of solution computed by the host;
 *   sollolo_h has lowest doubles of solution computed by the host;
 *   solhihi_d has highest doubles of solution computed by the device;
 *   sollohi_d has second highest doubles of solution computed by the device;
 *   solhilo_d has second lowest doubles of solution computed by the device;
 *   sollolo_d has lowest doubles of solution computed by the device;
 *   Qhihi_h   highest doubles of Q of the QR computed by the host;
 *   Qlohi_h   second highest doubles of Q of the QR computed by the host;
 *   Qhilo_h   second lowest doubles of Q of the QR computed by the host;
 *   Qlolo_h   lowest doubles of Q of the QR computed by the host;
 *   Qhihi_d   highest doubles of Q of the QR computed by the device;
 *   Qlohi_d   second highest doubles of Q of the QR computed by the device;
 *   Qhilo_d   second lowest doubles of Q of the QR computed by the device;
 *   Qlolo_d   lowest doubles of Q of the QR computed by the device;
 *   Rhihi_h   highest doubles of R of the QR computed by the host;
 *   Rlohi_h   second highest doubles of R of the QR computed by the host;
 *   Rhilo_h   second lowest doubles of R of the QR computed by the host;
 *   Rlolo_h   lowest doubles of R of the QR computed by the host;
 *   Rhihi_d   highest doubles of R of the QR computed by the device;
 *   Rlohi_d   second highest doubles of R of the QR computed by the device;
 *   Rhilo_d   second lowest doubles of R of the QR computed by the device;
 *   Rlolo_d   lowest doubles of R of the QR computed by the device;
 *   wrkmathihi has a copy of the highest doubles of the Jacobian matrix;
 *   wrkmatlohi has a copy of the second highest doubles of the Jacobian;
 *   wrkmathilo has a copy of the second lowest doubles of the Jacobian;
 *   wrkmatlolo has a copy of the lowest doubles of the Jacobian matrix;
 *   resvechihi has highest doubles of the residual vectors;
 *   resveclohi has second highest doubles of the residual vectors;
 *   resvechilo has second lowest doubles of the residual vectors;
 *   resveclolo has lowest doubles of the residual vectors;
 *   resmaxhihi is highest double of the maximum of residual vectors;
 *   resmaxlohi is second highest double of the maximum of residual vectors;
 *   resmaxhilo is second lowest double of the maximum of residual vectors;
 *   resmaxlolo is lowest double of the maximum of residual vectors;
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

int dbl4_allocate_inoutfunjac
 ( int dim, int deg, int mode,
   double **inputhihi_h, double **inputlohi_h,
   double **inputhilo_h, double **inputlolo_h,
   double **inputhihi_d, double **inputlohi_d,
   double **inputhilo_d, double **inputlolo_d,
   double ***outputhihi_h, double ***outputlohi_h, 
   double ***outputhilo_h, double ***outputlolo_h, 
   double ***outputhihi_d, double ***outputlohi_d,
   double ***outputhilo_d, double ***outputlolo_d,
   double **funvalhihi_h, double **funvallohi_h,
   double **funvalhilo_h, double **funvallolo_h,
   double **funvalhihi_d, double **funvallohi_d,
   double **funvalhilo_d, double **funvallolo_d,
   double ***jacvalhihi_h, double ***jacvallohi_h,
   double ***jacvalhilo_h, double ***jacvallolo_h,
   double ***jacvalhihi_d, double ***jacvallohi_d,
   double ***jacvalhilo_d, double ***jacvallolo_d );
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
 *   inputhihi_h are the highest doubles of the input on the host;
 *   inputlohi_h are the 2nd highest doubles of the input on the host;
 *   inputhilo_h are the 2nd lowest doubles of the input on the host;
 *   inputlolo_h are the lowest doubles of the input on the host;
 *   inputhihi_d are the highest doubles of the input on the device;
 *   inputlohi_d are the 2nd highest doubles of the input on the device;
 *   inputhilo_d are the 2nd lowest doubles of the input on the device;
 *   inputlolo_d are the lowest doubles of the input on the device;
 *   outputhihi_h are the highest doubles of the output on the host;
 *   outputlohi_h are the 2nd highest doubles of the output on the host;
 *   outputhilo_h are the 2nd lowest doubles of the output on the host;
 *   outputlolo_h are the lowest doubles of the output on the host;
 *   outputhihi_d are the highest doubles of the output on the device;
 *   outputlohi_d are the 2nd highest doubles of the output on the device;
 *   outputhilo_d are the 2nd lowest doubles of the output on the device;
 *   outputlolo_d are the lowest doubles of the output on the device;
 *   funvalhihi_h are the highest doubles of the function values on host;
 *   funvallohi_h are the 2nd highest doubles of the function values on host;
 *   funvalhilo_h are the 2nd lowest doubles of the function values on host;
 *   funvallolo_h are the lowest doubles of the function values on host;
 *   funvalhihi_d are the highest doubles of the function values on device;
 *   funvallohi_d are the 2nd highest doubles of the function values on device;
 *   funvalhilo_d are the 2nd lowest doubles of the function values on device;
 *   funvallolo_d are the lowest doubles of the function values on device;
 *   jacvalhihi_h are the highest doubles of the Jacobian on the host;
 *   jacvallohi_h are the 2nd highest doubles of the Jacobian on the host;
 *   jacvalhilo_h are the 2nd lowest doubles of the Jacobian on the host;
 *   jacvallolo_h are the lowest doubles of the Jacobian on the host;
 *   jacvalhihi_d are the highest doubles of the Jacobian on the device;
 *   jacvallohi_d are the 2nd highest doubles of the Jacobian on the device;
 *   jacvalhilo_d are the 2nd lowest doubles of the Jacobian on the device;
 *   jacvalhilo_d are the lowest doubles of the Jacobian on the device. */

int dbl4_allocate_rhsqrsol
 ( int dim, int deg, int mode,
   double **rhshihi_h, double **rhslohi_h,
   double **rhshilo_h, double **rhslolo_h,
   double **rhshihi_d, double **rhslohi_d,
   double **rhshilo_d, double **rhslolo_d,
   double **urhshihi_h, double **urhslohi_h,
   double **urhshilo_h, double **urhslolo_h,
   double **urhshihi_d, double **urhslohi_d,
   double **urhshilo_d, double **urhslolo_d,
   double **Qhihi_h, double **Qlohi_h, double **Qhilo_h, double **Qlolo_h,
   double **Qhihi_d, double **Qlohi_d, double **Qhilo_d, double **Qlolo_d,
   double **Rhihi_h, double **Rlohi_h, double **Rhilo_h, double **Rlolo_h,
   double **Rhihi_d, double **Rlohi_d, double **Rhilo_d, double **Rlolo_d,
   double **solhihi_h, double **sollohi_h,
   double **solhilo_h, double **sollolo_h,
   double **solhihi_d, double **sollohi_d,
   double **solhilo_d, double **sollolo_d );
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
 *   rhshihi_h  highest doubles of right-hand side on host;
 *   rhslohi_h  second highest doubles of right-hand side on host;
 *   rhshilo_h  second lowest doubles of right-hand side on host;
 *   rhslolo_h  lowest doubles of right-hand side on host;
 *   rhshihi_d  highest doubles of right-hand side on device;
 *   rhslohi_d  second highest doubles of right-hand side on device;
 *   rhshilo_d  second lowest doubles of right-hand side on device;
 *   rhslolo_d  lowest doubles of right-hand side on device;
 *   urhshihi_h are highest doubles of updated right-hand side on host;
 *   urhslohi_h are 2nd highest doubles of updated right-hand side on host;
 *   urhshilo_h are 2nd lowest doubles of updated right-hand side on host;
 *   urhslolo_h are lowest doubles of updated right-hand side on host;
 *   urhshihi_d are highest doubles of updated right hand on device;
 *   urhslohi_d are 2nd highest doubles of updated right hand on device;
 *   urhshilo_d are 2nd lowest doubles of updated right hand on device;
 *   urhslolo_d are  lowest doubles of updated right hand on device;
 *   Qhihi_h    highest doubles of Q of the QR computed by the host;
 *   Qlohi_h    second highest doubles of Q of the QR computed by the host;
 *   Qhilo_h    second lowest doubles of Q of the QR computed by the host;
 *   Qlolo_h    lowest doubles of Q of the QR computed by the host;
 *   Qhihi_d    highest doubles of Q of the QR computed by the device;
 *   Qlohi_d    second highest doubles of Q of the QR computed by the device;
 *   Qhilo_d    second lowest doubles of Q of the QR computed by the device;
 *   Qlolo_d    lowest doubles of Q of the QR computed by the device;
 *   Rhihi_h    highest doubles of R of the QR computed by the host;
 *   Rlohi_h    second highest doubles of R of the QR computed by the host;
 *   Rhilo_h    second lowest doubles of R of the QR computed by the host;
 *   Rlolo_h    lowest doubles of R of the QR computed by the host;
 *   Rhihi_d    highest doubles of R of the QR computed by the device;
 *   Rlohi_d    second highest doubles of R of the QR computed by the device;
 *   Rhilo_d    second lowest doubles of R of the QR computed by the device;
 *   Rlolo_d    lowest doubles of R of the QR computed by the device;
 *   solhihi_h has highest doubles of solution computed by the host;
 *   sollohi_h has second highest doubles of solution computed by the host;
 *   solhilo_h has second lowest doubles of solution computed by the host;
 *   sollolo_h has lowest doubles of solution computed by the host;
 *   solhihi_d has highest doubles of solution computed by the device;
 *   sollohi_d has second highest doubles of solution computed by the device;
 *   solhilo_d has second lowest doubles of solution computed by the device;
 *   sollolo_d has lowest doubles of solution computed by the device. */

void dbl4_start_setup
 ( int dim, int deg,
   double **testsolhihi, double **testsollohi,
   double **testsolhilo, double **testsollolo,
   double **inputhihi_h, double **inputlohi_h,
   double **inputhilo_h, double **inputlolo_h,
   double **inputhihi_d, double **inputlohi_d,
   double **inputhilo_d, double **inputlolo_d, int mode, int vrblvl );
/*
 * DESCRIPTION :
 *   Given the test solution, defines the start vector.
 *
 * ON ENTRY :
 *   dim        dimension of the solution vectors;
 *   deg        degree at which the series are truncated;
 *   testsolhihi are the highest doubles of the test solution;
 *   testsollohi are the highest doubles of the test solution;
 *   testsolhilo are the lowest doubles of the test solution;
 *   testsollolo are the lowest doubles of the test solution;
 *   inputhihi_h is allocated on host if mode is 1 or 2;
 *   inputlohi_h is allocated on host if mode is 1 or 2;
 *   inputhilo_h is allocated on host if mode is 1 or 2;
 *   inputlolo_h is allocated on host if mode is 1 or 2;
 *   inputhihi_d is allocated on device if mode is 0 or 2;
 *   inputlohi_d is allocated on device if mode is 0 or 2;
 *   inputhilo_d is allocated on device if mode is 0 or 2;
 *   inputlolo_d is allocated on device if mode is 0 or 2;
 *   mode       0 (GPU), 1 (CPU), or 2 (GPU+CPU);
 *   vrblvl     is the verbose level.
 *
 * ON RETURN :
 *   inputhihi_h are the highest doubles of start vector
 *              on host if mode is 1 or 2;
 *   inputlohi_h are the second highest doubles of start vector
 *              on host if mode is 1 or 2;
 *   inputhilo_h are the second lowest doubles of start vector
 *              on host if mode is 1 or 2;
 *   inputlolo_h are the lowest doubles of start vector
 *              on host if mode is 1 or 2;
 *   inputhihi_d are the highest doubles of start vector
 *              on device if mode is 0 or 2;
 *   inputlohi_d are the second highest doubles of start vector
 *              on device if mode is 0 or 2;
 *   inputhilo_d are the second lowest doubles of start vector
 *              on device if mode is 0 or 2;
 *   inputlolo_d are the lowest doubles of start vector
 *              on device if mode is 0 or 2. */

void dbl4_column_setup
 ( int dim, int deg, int nbrcol, int **nvr, int ***idx, int **rowsA,
   double ***cffhihi, double ***cfflohi,
   double ***cffhilo, double ***cfflolo,
   double **testsolhihi, double **testsollohi,
   double **testsolhilo, double **testsollolo,
   double **mbrhshihi, double **mbrhslohi,
   double **mbrhshilo, double **mbrhslolo,
   double **inputhihi_h, double **inputlohi_h,
   double **inputhilo_h, double **inputlolo_h,
   double **inputhihi_d, double **inputlohi_d,
   double **inputhilo_d, double **inputlolo_d, int mode, int vrblvl );
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
 *   testsolhihi has space for dim pointers;
 *   testsollohi has space for dim pointers;
 *   testsolhilo has space for dim pointers;
 *   testsollolo has space for dim pointers;
 *   mbrhshihi  space for dim pointers;
 *   mbrhslohi  space for dim pointers;
 *   mbrhshilo  space for dim pointers;
 *   mbrhslolo  space for dim pointers;
 *   inputhihi_h is allocated on host if mode is 1 or 2;
 *   inputlohi_h is allocated on host if mode is 1 or 2;
 *   inputhilo_h is allocated on host if mode is 1 or 2;
 *   inputlolo_h is allocated on host if mode is 1 or 2;
 *   inputhihi_d is allocated on device if mode is 0 or 2;
 *   inputlohi_d is allocated on device if mode is 0 or 2;
 *   inputhilo_d is allocated on device if mode is 0 or 2;
 *   inputlolo_d is allocated on device if mode is 0 or 2;
 *   mode       0 (GPU), 1 (CPU), or 2 (GPU+CPU);
 *   vrblvl     is the verbose level.
 *
 * ON RETURN :
 *   testsolhihi are the highest doubles of test solution;
 *   testsollohi are the second highest doubles of test solution;
 *   testsolhilo are the second lowest doubles of test solution;
 *   testsollolo are the lowest doubles of test solution;
 *   mbrhshihi  highest doubles of right hand side vector
 *              for the test solution;
 *   mbrhslohi  second highest doubles of right hand side vector
 *              for the test solution;
 *   mbrhshilo  second lowest doubles of right hand side vector
 *              for the test solution;
 *   mbrhslolo  lowest doubles of right hand side vector
 *              for the test solution;
 *   inputhihi_h are the highest doubles of start vector
 *              on host if mode is 1 or 2;
 *   inputlohi_h are the second highest doubles of start vector
 *              on host if mode is 1 or 2;
 *   inputhilo_h are the second lowest doubles of start vector
 *              on host if mode is 1 or 2;
 *   inputlolo_h are the lowest doubles of start vector
 *              on host if mode is 1 or 2;
 *   inputhihi_d are the highest doubles of start vector
 *              on device if mode is 0 or 2;
 *   inputlohi_d are the second highest doubles of start vector
 *              on device if mode is 0 or 2;
 *   inputhilo_d are the second lowest doubles of start vector
 *              on device if mode is 0 or 2;
 *   inputlolo_d are the lowest doubles of start vector
 *              on device if mode is 0 or 2. */

void dbl4_row_setup
 ( int dim, int deg, int *nbr, int **nvr, int ***idx,
   double **csthihi, double **cstlohi,
   double **csthilo, double **cstlolo,
   double ***cffhihi, double ***cfflohi,
   double ***cffhilo, double ***cfflolo,
   double **testsolhihi, double **testsollohi,
   double **testsolhilo, double **testsollolo,
   double **inputhihi_h, double **inputlohi_h,
   double **inputhilo_h, double **inputlolo_h,
   double **inputhihi_d, double **inputlohi_d,
   double **inputhilo_d, double **inputlolo_d,
   double ***outputhihi_h, double ***outputlohi_h,
   double ***outputhilo_h, double ***outputlolo_h,
   double ***outputhihi_d, double ***outputlohi_d,
   double ***outputhilo_d, double ***outputlolo_d, int mode, int vrblvl );
/*
 * DESCRIPTION :
 *   Defines the test solution and start vectors to run Newton's method
 *   on an indexed polynomial system.
 *
 * ON ENTRY :
 *   dim        dimension of the solution vectors;
 *   deg        degree at which the series are truncated;
 *   nbr        nbr[i] is the number of monomials in the i-th polynomial;
 *   nvr        nvr[i][j] is the number of variables in the j-th monomial
 *              of the i-th polynomial;
 *   idx        idx[i][j] are the indices of the variables in monomial j
 *              of the i-th polynomial;
 *   csthihi    highest doubles of constants;
 *   cstlohi    second highest doubles of constants;
 *   csthilo    second lowest doubles of constants;
 *   cstlolo    lowest doubles of constants;
 *   cffhihi    highest doubles of coefficients;
 *   cfflohi    second highest doubles of coefficients;
 *   cffhilo    second lowest doubles of coefficients;
 *   cfflolo    lowest doubles of coefficients;
 *   testsolhihi has space for dim pointers;
 *   testsollohi has space for dim pointers;
 *   testsolhilo has space for dim pointers;
 *   testsollolo has space for dim pointers;
 *   inputhihi_h is allocated on host if mode is 1 or 2;
 *   inputlohi_h is allocated on host if mode is 1 or 2;
 *   inputhilo_h is allocated on host if mode is 1 or 2;
 *   inputlolo_h is allocated on host if mode is 1 or 2;
 *   inputhihi_d is allocated on device if mode is 0 or 2;
 *   inputlohi_d is allocated on device if mode is 0 or 2;
 *   inputhilo_d is allocated on device if mode is 0 or 2;
 *   inputlolo_d is allocated on device if mode is 0 or 2;
 *   outputhihi_h is allocated on host if mode is 1 or 2;
 *   outputlohi_h is allocated on host if mode is 1 or 2;
 *   outputhilo_h is allocated on host if mode is 1 or 2;
 *   outputlolo_h is allocated on host if mode is 1 or 2;
 *   outputhihi_d is allocated on device if mode is 0 or 2;
 *   outputlohi_d is allocated on device if mode is 0 or 2;
 *   outputhilo_d is allocated on device if mode is 0 or 2;
 *   outputlolo_d is allocated on device if mode is 0 or 2;
 *   mode       0 (GPU), 1 (CPU), or 2 (GPU+CPU);
 *   vrblvl     is the verbose level.
 *
 * ON RETURN :
 *   csthihi    highest doubles of constant, 
 *              adjusted for the test solution;
 *   cstlohi    second highest doubles of constant, 
 *              adjusted for the test solution;
 *   csthilo    second lowest doubles of constant, 
 *              adjusted for the test solution;
 *   cstlolo    lowest doubles of constant, 
 *              adjusted for the test solution;
 *   testsolhihi are the highest doubles of test solution;
 *   testsollohi are the second highest doubles of test solution;
 *   testsolhilo are the second lowest doubles of test solution;
 *   testsollolo are the lowest doubles of test solution;
 *   inputhihi_h are the highest doubles of start vector
 *              on host if mode is 1 or 2;
 *   inputlohi_h are the second highest doubles of start vector
 *              on host if mode is 1 or 2;
 *   inputhilo_h are the second lowest doubles of start vector
 *              on host if mode is 1 or 2;
 *   inputlolo_h are the lowest doubles of start vector
 *              on host if mode is 1 or 2;
 *   inputhihi_d are the highest doubles of start vector
 *              on device if mode is 0 or 2;
 *   inputlohi_d are the second highest doubles of start vector
 *              on device if mode is 0 or 2;
 *   inputhilo_d are the second lowest doubles of start vector
 *              on device if mode is 0 or 2;
 *   inputlolo_d are the lowest doubles of start vector
 *              on device if mode is 0 or 2;
 *   outputhihi_h are the highest doubles 
 *              of evaluated test solution if mode is 1;
 *   outputlohi_h are the second highest doubles 
 *              of evaluated test solution if mode is 1;
 *   outputhilo_h are the second lowest doubles
 *              of evaluated test solution if mode is 1;
 *   outputlolo_h are the lowest doubles
 *              of evaluated test solution if mode is 1;
 *   outputhihi_d are the highest doubles
 *              of evaluated test solution, if mode is 0 or 2;
 *   outputlohi_d are the second highest doubles
 *              of evaluated test solution, if mode is 0 or 2;
 *   outputhilo_d are the second lowest doubles
 *              of evaluated test solution, if mode is 0 or 2;
 *   outputlolo_d are the lowest doubles
 *              of evaluated test solution, if mode is 0 or 2. */

int dbl4_error_testsol
 ( int dim, int deg, int mode,
   double **testsolhihi, double **testsollohi,
   double **testsolhilo, double **testsollolo,
   double **inputhihi_h, double **inputlohi_h,
   double **inputhilo_h, double **inputlolo_h,
   double **inputhihi_d, double **inputlohi_d,
   double **inputhilo_d, double **inputlolo_d );
/*
 * DESCRIPTION :
 *   Compares the computed solution against the test solution.
 *
 * ON ENTRY :
 *   dim        dimension of the solution vectors;
 *   deg        degree at which the series are truncated;
 *   mode       0 (GPU), 1 (CPU), or 2 (GPU+CPU).
 *   testsolhihi are the highest doubles of the test solution;
 *   testsollohi are the second highest doubles of the test solution;
 *   testsolhilo are the second lowest doubles of the test solution;
 *   testsollolo are the lowest doubles of the test solution;
 *   inputhihi_h are the highest doubles on host if mode is 1 or 2;
 *   inputlohi_h are the second highest doubles on host if mode is 1 or 2;
 *   inputhilo_h are the second lowest doubles on host if mode is 1 or 2;
 *   inputlolo_h are the lowest doubles on host if mode is 1 or 2;
 *   inputhihi_d are the highest doubles on device if mode is 0 or 2;
 *   inputlohi_d are the second highest doubles on device if mode is 0 or 2;
 *   inputhilo_d are the second lowest doubles on device if mode is 0 or 2;
 *   inputlolo_d are the lowest doubles on device if mode is 0 or 2. */

int test_dbl4_column_newton
 ( int szt, int nbt, int dim, int deg, int nbrcol,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac, int **rowsA,
   double dpr, int nbsteps, int mode, int vrblvl );
/*
 * DESCRIPTION :
 *   Runs Newton with real quad double arithmetic,
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

int test_dbl4_row_newton
 ( int szt, int nbt, int dim, int deg, int *nbr, int **nvr, int ***idx,
   double dpr, int nbsteps, int mode, int vrblvl );
/*
 * DESCRIPTION :
 *   Runs Newton with real quad double arithmetic,
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
