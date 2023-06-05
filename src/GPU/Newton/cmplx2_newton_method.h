// The file cmplx2_newton_method.h specifies Newton's method
// on series in double precision on complex numbers.

#ifndef __cmplx2_newton_method_h__
#define __cmplx2_newton_method_h__

int cmplx2_errors_funjacrhs
 ( int dim, int deg,
   double **funvalrehi_h, double **funvalrelo_h,
   double **funvalimhi_h, double **funvalimlo_h,
   double **funvalrehi_d, double **funvalrelo_d,
   double **funvalimhi_d, double **funvalimlo_d,
   double ***jacvalrehi_h, double ***jacvalrelo_h,
   double ***jacvalimhi_h, double ***jacvalimlo_h,
   double ***jacvalrehi_d, double ***jacvalrelo_d,
   double ***jacvalimhi_d, double ***jacvalimlo_d,
   double **rhsrehi_h, double **rhsrelo_h,
   double **rhsimhi_h, double **rhsimlo_h,
   double **rhsrehi_d, double **rhsrelo_d,
   double **rhsimhi_d, double **rhsimlo_d, int vrblvl );
/*
 * DESCRIPTION :
 *   Compares results on host with those on the device
 *   and returns the number of errors larger than 1.0e-20.
 *
 * ON ENTRY :
 *   dim       number of equations and variables;
 *   deg       degree at which the series are truncated;
 *   funvalrehi_h are the high doubles of the real parts
 *             of the function values on the host;
 *   funvalrelo_h are the low doubles of the real parts
 *             of the function values on the host;
 *   funvalimhi_h are the high doubles of the imaginary parts
 *             of the function values on the host;
 *   funvalimlo_h are the low doubles of the imaginary parts
 *             of the function values on the host;
 *   funvalrehi_d are the high doubles of the real parts
 *             of the function values on the device;
 *   funvalrelo_d are the low doubles of the real parts
 *             of the function values on the device;
 *   funvalimhi_d are the high doubles of the imaginary parts
 *             of the function values on the device;
 *   funvalimlo_d are the low doubles of the imaginary parts
 *             of the function values on the device;
 *   jacvalrehi_h are the high doubles of the real parts
 *             of the Jacobian matrix on the host;
 *   jacvalrelo_h are the low doubles of the real parts
 *             of the Jacobian matrix on the host;
 *   jacvalimhi_h are the high doubles of the imaginary parts
 *             of the Jacobian matrix on the host;
 *   jacvalimlo_h are the low doubles of the imaginary parts
 *             of the Jacobian matrix on the host;
 *   jacvalrehi_d are the high doubles of the real parts
 *             of the Jacobian matrix on the device;
 *   jacvalrelo_d are the low doubles of the real parts
 *             of the Jacobian matrix on the device;
 *   jacvalimhi_d are the high doubles of the imaginary parts
 *             of the Jacobian matrix on the device;
 *   jacvalimlo_d are the low doubles of the imaginary parts
 *             of the Jacobian matrix on the device;
 *   rhsrehi_h are the high doubles of real parts
 *             of the right hand side vector on the host;
 *   rhsrelo_h are the low doubles of real parts
 *             of the right hand side vector on the host;
 *   rhsimhi_h are the high doubles of imaginary parts
 *             of the right hand side vector on the host;
 *   rhsimlo_h are the low doubles of imaginary parts
 *             of the right hand side vector on the host;
 *   rhsrehi_d are the high doubles of real parts
 *             of the right hand side vector on the device;
 *   rhsrelo_d are the low doubles of real parts
 *             of the right hand side vector on the device;
 *   rhsimhi_d are the high doubles of imaginary parts
 *             of the right hand side vector on the device;
 *   rhsimlo_d are the low doubles of imaginary parts
 *             of the right hand side vector on the device;
 *   vrblvl    is the verbose level. */

int cmplx2_errors_inurhsQRsol
 ( int dim, int deg,
   double **inputrehi_h, double **inputrelo_h,
   double **inputimhi_h, double **inputimlo_h,
   double **inputrehi_d, double **inputrelo_d,
   double **inputimhi_d, double **inputimlo_d,
   double **Qrehi_h, double **Qrelo_h, double **Qimhi_h, double **Qimlo_h, 
   double **Qrehi_d, double **Qrelo_d, double **Qimhi_d, double **Qimlo_d,
   double **Rrehi_h, double **Rrelo_h, double **Rimhi_h, double **Rimlo_h, 
   double **Rrehi_d, double **Rrelo_d, double **Rimhi_d, double **Rimlo_d,
   double **urhsrehi_h, double **urhsrelo_h, 
   double **urhsimhi_h, double **urhsimlo_h,
   double **urhsrehi_d, double **urhsrelo_d, 
   double **urhsimhi_d, double **urhsimlo_d,
   double **solrehi_h, double **solrelo_h,
   double **solimhi_h, double **solimlo_h,
   double **solrehi_d, double **solrelo_d,
   double **solimhi_d, double **solimlo_d, int vrblvl );
/*
 * DESCRIPTION :
 *   Compares results on host with those on the device
 *   and returns the number of errors larger than 1.0e-20.
 *
 * ON ENTRY :
 *   dim       number of equations and variables;
 *   deg       degree at which the series are truncated;
 *   inputrehi_h are the high doubles of the real parts
 *             of the updated series on the host;
 *   inputrelo_h are the low doubles of the real parts
 *             of the updated series on the host;
 *   inputimhi_h are the high doubles of the imaginary parts
 *             of the updated series on the host;
 *   inputimlo_h are the low doubles of the imaginary parts
 *             of the updated series on the host;
 *   inputrehi_d are the high doubles of the real parts
 *             of the updated series on the device;
 *   inputrelo_d are the low doubles of the real parts
 *             of the updated series on the device;
 *   inputimhi_d are the high doubles of the imaginary parts
 *             of the updated series on the device;
 *   inputimlo_d are the low doubles of the imaginary parts
 *             of the updated series on the device;
 *   Qrehi_h   high doubles of the real parts of Q on the host;
 *   Qrelo_h   low doubles of the real parts of Q on the host;
 *   Qimhi_h   high doubles of the imaginary parts of Q on the host;
 *   Qimlo_h   low doubles of the imaginary parts of Q on the host;
 *   Qrehi_d   high doubles of the real parts of Q on the device;
 *   Qrelo_d   low doubles of the real parts of Q on the device;
 *   Qimhi_d   high doubles of the imaginary parts of Q on the device;
 *   Qimlo_d   low doubles of the imaginary parts of Q on the device;
 *   Rrehi_h   high doubles of the real parts of R on the host;
 *   Rrelo_h   low doubles of the real parts of R on the host;
 *   Rimhi_h   high doubles of the imaginary parts of R on the host;
 *   Rimlo_h   low doubles of the imaginary parts of R on the host;
 *   Rrehi_d   high doubles of the real parts of R on the device;
 *   Rrelo_d   low doubles of the real parts of R on the device;
 *   Rimhi_d   high doubles of the imaginary parts of R on the device;
 *   Rimlo_d   low doubles of the imaginary parts of R on the device;
 *   urhsrehi_h are the high doubles of the real parts
 *             of the updated right hand side on the host;
 *   urhsrelo_h are the low doubles of the real parts
 *             of the updated right hand side on the host;
 *   urhsimhi_h are the high doubles of the imaginary parts
 *             of the updated right hand side on the host;
 *   urhsimlo_h are the low doubles of the imaginary parts
 *             of the updated right hand side on the host;
 *   urhsrehi_d are the high doubles of the real parts
 *             of the updated right hand side on the device;
 *   urhsrelo_d are the low doubles of the real parts
 *             of the updated right hand side on the device;
 *   urhsimhi_d are the high doubles of the imaginary parts
 *             of the updated right hand side on the device;
 *   urhsimlo_d are the low doubles of the imaginary parts
 *             of the updated right hand side on the device;
 *   solrehi_h are the high doubles of the real parts
 *             of the update to the solution on the host,
 *   solrelo_h are the lowh doubles of the real parts
 *             of the update to the solution on the host,
 *   solimhi_h are the high doubles of the imaginary parts
 *             of the update to the solution on the host,
 *   solimlo_h are the low doubles of the imaginary parts
 *             of the update to the solution on the host,
 *   solrehi_d are the high doubles of the real parts
 *             of the update to the solution on the device;
 *   solrelo_d are the low doubles of the real parts
 *             of the update to the solution on the device;
 *   solimhi_d are the high doubles of the imaginary parts
 *             of the update to the solution on the device;
 *   solimlo_d are the low doubles of the imaginary parts
 *             of the update to the solution on the device;
 *   vrblvl    is the verbose level. */

int cmplx2_update_newton_qrstep
 ( int szt, int nbt, int dim, int deg,
   int *tailidx_h, int *tailidx_d,
   double **inputrehi_h, double **inputrelo_h,
   double **inputimhi_h, double **inputimlo_h,
   double **inputrehi_d, double **inputrelo_d,
   double **inputimhi_d, double **inputimlo_d,
   double **funvalrehi_h, double **funvalrelo_h,
   double **funvalimhi_h, double **funvalimlo_h,
   double **funvalrehi_d, double **funvalrelo_d,
   double **funvalimhi_d, double **funvalimlo_d,
   double ***jacvalrehi_h, double ***jacvalrelo_h,
   double ***jacvalimhi_h, double ***jacvalimlo_h,
   double ***jacvalrehi_d, double ***jacvalrelo_d,
   double ***jacvalimhi_d, double ***jacvalimlo_d,
   double **rhsrehi_h, double **rhsrelo_h,
   double **rhsimhi_h, double **rhsimlo_h,
   double **rhsrehi_d, double **rhsrelo_d, 
   double **rhsimhi_d, double **rhsimlo_d,
   double **urhsrehi_h, double **urhsrelo_h,
   double **urhsimhi_h, double **urhsimlo_h,
   double **urhsrehi_d, double **urhsrelo_d,
   double **urhsimhi_d, double **urhsimlo_d,
   double **solrehi_h, double **solrelo_h,
   double **solimhi_h, double **solimlo_h, 
   double **solrehi_d, double **solrelo_d, 
   double **solimhi_d, double **solimlo_d,
   double **Qrehi_h, double **Qrelo_h, double **Qimhi_h, double **Qimlo_h,
   double **Qrehi_d, double **Qrelo_d, double **Qimhi_d, double **Qimlo_d, 
   double **Rrehi_h, double **Rrelo_h, double **Rimhi_h, double **Rimlo_h, 
   double **Rrehi_d, double **Rrelo_d, double **Rimhi_d, double **Rimlo_d,
   double *workvecrehi, double *workvecrelo,
   double *workvecimhi, double *workvecimlo,
   double **resvecrehi, double **resvecrelo,
   double **resvecimhi, double **resvecimlo,
   double *resmaxhi, double *resmaxlo,
   bool* zeroQ_h, bool *noqr_h, bool* zeroQ_d, bool *noqr_d,
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
 *   inputrehi_h are the high doubles of the real parts coefficients of series
 *             of degree deg, for dim variables, computed on host;
 *   inputrelo_h are the low doubles of the real parts coefficients of series
 *             of degree deg, for dim variables, computed on host;
 *   inputimhi_h are the high doubles of the imaginary parts coefficients of
 *             the series of degree deg for dim variables, computed on host;
 *   inputimlo_h are the low doubles of the imaginary parts coefficients of
 *             the series of degree deg for dim variables, computed on host;
 *   inputrehi_d has space for power series computed on device;
 *   inputrelo_d has space for power series computed on device;
 *   inputimhi_d has space for power series computed on device;
 *   inputimlo_d has space for power series computed on device;
 *   funvalrehi_h has space for the evaluated series computed by host;
 *   funvalrelo_h has space for the evaluated series computed by host;
 *   funvalimhi_h has space for the evaluated series computed by host;
 *   funvalimlo_h has space for the evaluated series computed by host;
 *   funvalrehi_d has space for the evaluated series computed by device;
 *   funvalrelo_d has space for the evaluated series computed by device;
 *   funvalimhi_d has space for the evaluated series computed by device;
 *   funvalimlo_d has space for the evaluated series computed by device;
 *   jacvalrehi_h has space for deg+1 matrices of dimension dim on host;
 *   jacvalrelo_h has space for deg+1 matrices of dimension dim on host;
 *   jacvalimhi_h has space for deg+1 matrices of dimension dim on host;
 *   jacvalimlo_h has space for deg+1 matrices of dimension dim on host;
 *   jacvalrehi_d has space for deg+1 matrices of dimension dim on device;
 *   jacvalrelo_d has space for deg+1 matrices of dimension dim on device;
 *   jacvalimhi_d has space for deg+1 matrices of dimension dim on device;
 *   jacvalimlo_d has space for deg+1 matrices of dimension dim on device;
 *   rhsrehi_h has space for deg+1 vectors of dimension dim on host;
 *   rhsrelo_h has space for deg+1 vectors of dimension dim on host;
 *   rhsimhi_h has space for deg+1 vectors of dimension dim on host;
 *   rhsimlo_h has space for deg+1 vectors of dimension dim on host;
 *   rhsrehi_d has space for deg+1 vectors of dimension dim on device;
 *   rhsrelo_d has space for deg+1 vectors of dimension dim on device;
 *   rhsimhi_d has space for deg+1 vectors of dimension dim on device;
 *   rhsimlo_d has space for deg+1 vectors of dimension dim on device;
 *   urhsrehi_h has space for updated right hand side computed by host;
 *   urhsrelo_h has space for updated right hand side computed by host;
 *   urhsimhi_h has space for updated right hand side computed by host;
 *   urhsimlo_h has space for updated right hand side computed by host;
 *   urhsrehi_d has space for updated right hand side computed by device; 
 *   urhsrelo_d has space for updated right hand side computed by device; 
 *   urhsimhi_d has space for updated right hand side computed by device; 
 *   urhsimlo_d has space for updated right hand side computed by device; 
 *   solrehi_h has space for deg+1 vectors of dimension dim;
 *   solrelo_h has space for deg+1 vectors of dimension dim;
 *   solimhi_h has space for deg+1 vectors of dimension dim;
 *   solimlo_h has space for deg+1 vectors of dimension dim;
 *   solrehi_d has space for deg+1 vectors of dimension dim;
 *   solrelo_d has space for deg+1 vectors of dimension dim;
 *   solimhi_d has space for deg+1 vectors of dimension dim;
 *   solimlo_d has space for deg+1 vectors of dimension dim;
 *   Qrehi_h   space allocated for the Q computed by the host;
 *   Qrelo_h   space allocated for the Q computed by the host;
 *   Qimhi_h   space allocated for the Q computed by the host;
 *   Qimlo_h   space allocated for the Q computed by the host;
 *   Qrehi_d   space allocated for the Q computed by the device;
 *   Qrelo_d   space allocated for the Q computed by the device;
 *   Qimhi_d   space allocated for the Q computed by the device;
 *   Qimlo_d   space allocated for the Q computed by the device;
 *   Rrehi_h   space allocated for the R computed by the host;
 *   Rrelo_h   space allocated for the R computed by the host;
 *   Rimhi_h   space allocated for the R computed by the host;
 *   Rimlo_h   space allocated for the R computed by the host;
 *   Rrehi_d   space allocated for the R computed by the device;
 *   Rrelo_d   space allocated for the R computed by the device;
 *   Rimhi_d   space allocated for the R computed by the device;
 *   Rimlo_d   space allocated for the R computed by the device;
 *   workvecrehi is work space allocated for a vector of dimension dim;
 *   workvecrelo is work space allocated for a vector of dimension dim;
 *   workvecimhi is work space allocated for a vector of dimension dim;
 *   workvecimlo is work space allocated for a vector of dimension dim;
 *   resvecrehi has space for deg+1 vectors of dimension dim;
 *   resvecrelo has space for deg+1 vectors of dimension dim;
 *   resvecimhi has space for deg+1 vectors of dimension dim;
 *   resvecimlo has space for deg+1 vectors of dimension dim;
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
 *   tailidx_h is the new value for tailidx_h;
 *   tailidx_d is the new value for tailidx_d;
 *   inputrehi_h has the high doubles of the real parts of series, on host;
 *   inputrelo_h has the low doubles of the real parts of series, on host;
 *   inputimhi_h has the high doubles of the imag parts of series, on host;
 *   inputimlo_h has the low doubles of the imag parts of series, on host;
 *   inputrehi_d has the high doubles of the real parts of series, on device;
 *   inputrelo_d has the low doubles of the real parts of series, on device;
 *   inputimhi_d has the high doubles of the imag parts of series, on device;
 *   inputimlo_d has the low doubles of the imag parts of series, on device;
 *   funvalrehi_h is outputrehi[i][dim], on host;
 *   funvalrelo_h is outputrelo[i][dim], on host;
 *   funvalimhi_h is outputimhi[i][dim], on host;
 *   funvalimlo_h is outputimlo[i][dim], on host;
 *   funvalrehi_d is outputrehi[i][dim], on device;
 *   funvalrelo_d is outputrelo[i][dim], on device;
 *   funvalimhi_d is outputimhi[i][dim], on device;
 *   funvalimlo_d is outputimlo[i][dim], on device;
 *   jacvalrehi_h are the high doubles of the real parts of a matrix series,
 *             on host;
 *   jacvalrelo_h are the low doubles of the real parts of a matrix series,
 *             on host;
 *   jacvalimhi_h are the high doubles of the real parts of a matrix series,
 *             on host;
 *   jacvalimlo_h are the low doubles of the real parts of a matrix series,
 *             on host;
 *   jacvalrehi_d are the high doubles of the real parts of a matrix series,
 *             on device;
 *   jacvalrelo_d are the low doubles of the real parts of a matrix series,
 *             on device;
 *   jacvalimhi_d are the high doubles of the real parts of a matrix series,
 *             on device;
 *   jacvalimlo_d are the low doubles of the real parts of a matrix series,
 *             on device;
 *   rhsrehi_h are high doubles of real parts of the linearized right hand
 *             side, subtracted by 1 and added by t, computed by host;
 *   rhsrelo_h are low doubles of real parts of the linearized right hand,
 *             side, subtracted by 1 and added by t, computed by host;
 *   rhsimhi_h are high doubles of imag parts of the linearized right hand,
 *             side, subtracted by 1 and added by t, computed by host;
 *   rhsimlo_h are low doubles of imag parts of the linearized right hand,
 *             side, subtracted by 1 and added by t, computed by host;
 *   rhsrehi_d are high doubles of real parts of the linearized right hand 
 *             side, subtracted by 1 and added by t, computed by device;
 *   rhsrelo_d are low doubles of real parts of the linearized right hand,
 *             side, subtracted by 1 and added by t, computed by device;
 *   rhsimhi_d are high doubles of imag parts of the linearized right hand,
 *             side, subtracted by 1 and added by t, computed by device;
 *   rhsimlo_d are low doubles of imag parts of the linearized right hand
 *             side, subtracted by 1 and added by t, computed by device;
 *   urhsrehi_h are the high doubles of the real parts of right hand side
 *             updated by the host;
 *   urhsrelo_h are the high doubles of the real parts of right hand side
 *             updated by the host;
 *   urhsimhi_h are the high doubles of the imaginary parts of right hand side
 *             updated by the host;
 *   urhsimlo_h are the high doubles of the imaginary parts of right hand side
 *             updated by the host;
 *   urhsrehi_d are the high doubles of the real parts of right hand side
 *             updated by the device;
 *   urhsrelo_d are the high doubles of the real parts of right hand side
 *             updated by the device;
 *   urhsimhi_d are the high doubles of the imaginary parts of right hand side
 *             updated by the device;
 *   urhsimlo_d are the high doubles of the imaginary parts of right hand side
 *             updated by the device;
 *   solrehi_h are the high doubles of real parts of solution on host;
 *   solrelo_h are the low doubles of real parts of solution on host;
 *   solimhi_h are the high doubles of imag parts of solution on host;
 *   solimlo_h are the low doubles of imag parts of solution on host;
 *   solrehi_d are the high doubles of real parts of solution on device;
 *   solrelo_d are the low doubles of real parts of solution on device;
 *   solimhi_d are the high doubles of imag parts of solution on device;
 *   solimlo_d are the low doubles of imag parts of solution on device;
 *   Qrehi_h   high doubles of real Q of the QR on host;
 *   Qrelo_h   low doubles of real Q of the QR on host;
 *   Qimhi_h   high doubles of imaginary Q of the QR on host;
 *   Qimlo_h   low doubles of imaginary Q of the QR on host;
 *   Qrehi_d   high doubles of real Q of the QR on device;
 *   Qrelo_d   low doubles of real Q of the QR on device;
 *   Qimhi_d   high doubles of imaginary Q of the QR on device;
 *   Qimlo_d   low doubles of imaginary Q of the QR on device;
 *   Rrehi_h   high doubles of real R of the QR on host;
 *   Rrelo_h   low doubles of real R of the QR on host;
 *   Rimhi_h   high doubles of imaginary R of the QR on host;
 *   Rimlo_h   low doubles of imaginary R of the QR on host;
 *   Rrehi_d   high doubles of real R of the QR on device;
 *   Rrelo_d   low doubles of real R of the QR on device;
 *   Rimhi_d   high doubles of imaginary R of the QR on device;
 *   Rimlo_d   low doubles of imaginary R of the QR on device;
 *   resvecrehi are high doubles of the real parts of the residual vectors;
 *   resvecrelo are high doubles of the real parts of the residual vectors;
 *   resvecimhi are high doubles of the imag parts of the residual vectors;
 *   resvecimlo are high doubles of the imag parts of the residual vectors;
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
 *   totqrlapsedms accumulates the milliseconds spent on Q times rhs;
 *   totbslapsedms accumulates the milliseconds spent on back substitutions;
 *   totupdlapsedms accumulates the milliseconds spent on updates;
 *   totreslapsedms accumulates the milliseconds spent on residuals. */

int cmplx2_column_newton_qrstep
 ( int szt, int nbt, int dim, int deg, int nbrcol,
   int *tailidx_h, int *tailidx_d,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac,
   double **mbrehi, double **mbrelo, double **mbimhi, double **mbimlo,
   double dpr,
   double ***cffrehi, double ***cffrelo, double ***cffimhi, double ***cffimlo,
   double **accrehi, double **accrelo, double **accimhi, double **accimlo,
   double **inputrehi_h, double **inputrelo_h,
   double **inputimhi_h, double **inputimlo_h,
   double **inputrehi_d, double **inputrelo_d,
   double **inputimhi_d, double **inputimlo_d,
   double ***outputrehi_h, double ***outputrelo_h,
   double ***outputimhi_h, double ***outputimlo_h,
   double ***outputrehi_d, double ***outputrelo_d,
   double ***outputimhi_d, double ***outputimlo_d,
   double **funvalrehi_h, double **funvalrelo_h,
   double **funvalimhi_h, double **funvalimlo_h,
   double **funvalrehi_d, double **funvalrelo_d,
   double **funvalimhi_d, double **funvalimlo_d,
   double ***jacvalrehi_h, double ***jacvalrelo_h,
   double ***jacvalimhi_h, double ***jacvalimlo_h,
   double ***jacvalrehi_d, double ***jacvalrelo_d,
   double ***jacvalimhi_d, double ***jacvalimlo_d,
   double **rhsrehi_h, double **rhsrelo_h,
   double **rhsimhi_h, double **rhsimlo_h,
   double **rhsrehi_d, double **rhsrelo_d, 
   double **rhsimhi_d, double **rhsimlo_d,
   double **urhsrehi_h, double **urhsrelo_h,
   double **urhsimhi_h, double **urhsimlo_h,
   double **urhsrehi_d, double **urhsrelo_d,
   double **urhsimhi_d, double **urhsimlo_d,
   double **solrehi_h, double **solrelo_h,
   double **solimhi_h, double **solimlo_h, 
   double **solrehi_d, double **solrelo_d, 
   double **solimhi_d, double **solimlo_d,
   double **Qrehi_h, double **Qrelo_h, double **Qimhi_h, double **Qimlo_h,
   double **Qrehi_d, double **Qrelo_d, double **Qimhi_d, double **Qimlo_d, 
   double **Rrehi_h, double **Rrelo_h, double **Rimhi_h, double **Rimlo_h, 
   double **Rrehi_d, double **Rrelo_d, double **Rimhi_d, double **Rimlo_d,
   double *workvecrehi, double *workvecrelo,
   double *workvecimhi, double *workvecimlo,
   double **resvecrehi, double **resvecrelo,
   double **resvecimhi, double **resvecimlo,
   double *resmaxhi, double *resmaxlo,
   bool* zeroQ_h, bool *noqr_h, bool* zeroQ_d, bool *noqr_d,
   int *upidx_h, int *bsidx_h, int *upidx_d, int *bsidx_d,
   double *totcnvlapsedms, double *totqrlapsedms, double *totqtblapsedms,
   double *totbslapsedms, double *totupdlapsedms, double *totreslapsedms,
   int vrblvl, int mode );
/*
 * DESCRIPTION :
 *   Does one step with Newton's method to update a power series,
 *   using QR factorization to solve linear systems, on real data,
 *   on one or more columns of monomials.
 *
 * REQUIRED : szt*nbt = dim for GPU computing.
 *
 * ON ENTRY :
 *   szt       size of each tile and block;
 *   nbt       number of tiles and number of blocks;
 *   dim       number of monomials;
 *   deg       degree of the power series;
 *   nbrcol    is the number of columns;
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
 *   mbrehi    high real parts of the right hand side vector of series;
 *   mbrelo    low real parts of the right hand side vector of series;
 *   mbimhi    high imaginary parts of the right hand side vector of series;
 *   mbimlo    low imaginary parts of the right hand side vector of series;
 *   dpr       damper multiplier for t, should be in (0.0, 1.0];
 *   cffrehi   high doubles of the real parts of the coefficients
 *             of the monomials;
 *   cffrelo   low doubles of the real parts of the coefficients
 *             of the monomials;
 *   cffimhi   high doubles of the imaginary parts of the coefficients
 *             of the monomials;
 *   cffimlo   low doubles of the imaginary parts of the coefficients
 *             of the monomials;
 *   accrehi   space to accumulate one power series of degree deg;
 *   accrelo   space to accumulate one power series of degree deg;
 *   accimhi   space to accumulate one power series of degree deg;
 *   accimlo   space to accumulate one power series of degree deg;
 *   inputrehi_h are the high doubles of the real parts coefficients of series
 *             of degree deg, for dim variables, computed on host;
 *   inputrelo_h are the low doubles of the real parts coefficients of series
 *             of degree deg, for dim variables, computed on host;
 *   inputimhi_h are the high doubles of the imaginary parts coefficients of
 *             the series of degree deg for dim variables, computed on host;
 *   inputimlo_h are the low doubles of the imaginary parts coefficients of
 *             the series of degree deg for dim variables, computed on host;
 *   inputrehi_d has space for power series computed on device;
 *   inputrelo_d has space for power series computed on device;
 *   inputimhi_d has space for power series computed on device;
 *   inputimlo_d has space for power series computed on device;
 *   outputrehi_h has space for the evaluated and differentiated monomials,
 *             computed on the host;
 *   outputrelo_h has space for the evaluated and differentiated monomials,
 *             computed on the host;
 *   outputimhi_h has space for the evaluated and differentiated monomials,
 *             computed on the host;
 *   outputimlo_h has space for the evaluated and differentiated monomials,
 *             computed on the host;
 *   outputrehi_d has space for the evaluated and differentiated monomials,
 *             computed on the device;
 *   outputrelo_d has space for the evaluated and differentiated monomials,
 *             computed on the device;
 *   outputimhi_d has space for the evaluated and differentiated monomials,
 *             computed on the device;
 *   outputimlo_d has space for the evaluated and differentiated monomials,
 *             computed on the device;
 *   funvalrehi_h has space for the evaluated series computed by host;
 *   funvalrelo_h has space for the evaluated series computed by host;
 *   funvalimhi_h has space for the evaluated series computed by host;
 *   funvalimlo_h has space for the evaluated series computed by host;
 *   funvalrehi_d has space for the evaluated series computed by device;
 *   funvalrelo_d has space for the evaluated series computed by device;
 *   funvalimhi_d has space for the evaluated series computed by device;
 *   funvalimlo_d has space for the evaluated series computed by device;
 *   jacvalrehi_h has space for deg+1 matrices of dimension dim on host;
 *   jacvalrelo_h has space for deg+1 matrices of dimension dim on host;
 *   jacvalimhi_h has space for deg+1 matrices of dimension dim on host;
 *   jacvalimlo_h has space for deg+1 matrices of dimension dim on host;
 *   jacvalrehi_d has space for deg+1 matrices of dimension dim on device;
 *   jacvalrelo_d has space for deg+1 matrices of dimension dim on device;
 *   jacvalimhi_d has space for deg+1 matrices of dimension dim on device;
 *   jacvalimlo_d has space for deg+1 matrices of dimension dim on device;
 *   rhsrehi_h has space for deg+1 vectors of dimension dim on host;
 *   rhsrelo_h has space for deg+1 vectors of dimension dim on host;
 *   rhsimhi_h has space for deg+1 vectors of dimension dim on host;
 *   rhsimlo_h has space for deg+1 vectors of dimension dim on host;
 *   rhsrehi_d has space for deg+1 vectors of dimension dim on device;
 *   rhsrelo_d has space for deg+1 vectors of dimension dim on device;
 *   rhsimhi_d has space for deg+1 vectors of dimension dim on device;
 *   rhsimlo_d has space for deg+1 vectors of dimension dim on device;
 *   urhsrehi_h has space for updated right hand side computed by host;
 *   urhsrelo_h has space for updated right hand side computed by host;
 *   urhsimhi_h has space for updated right hand side computed by host;
 *   urhsimlo_h has space for updated right hand side computed by host;
 *   urhsrehi_d has space for updated right hand side computed by device; 
 *   urhsrelo_d has space for updated right hand side computed by device; 
 *   urhsimhi_d has space for updated right hand side computed by device; 
 *   urhsimlo_d has space for updated right hand side computed by device; 
 *   solrehi_h has space for deg+1 vectors of dimension dim;
 *   solrelo_h has space for deg+1 vectors of dimension dim;
 *   solimhi_h has space for deg+1 vectors of dimension dim;
 *   solimlo_h has space for deg+1 vectors of dimension dim;
 *   solrehi_d has space for deg+1 vectors of dimension dim;
 *   solrelo_d has space for deg+1 vectors of dimension dim;
 *   solimhi_d has space for deg+1 vectors of dimension dim;
 *   solimlo_d has space for deg+1 vectors of dimension dim;
 *   Qrehi_h   space allocated for the Q computed by the host;
 *   Qrelo_h   space allocated for the Q computed by the host;
 *   Qimhi_h   space allocated for the Q computed by the host;
 *   Qimlo_h   space allocated for the Q computed by the host;
 *   Qrehi_d   space allocated for the Q computed by the device;
 *   Qrelo_d   space allocated for the Q computed by the device;
 *   Qimhi_d   space allocated for the Q computed by the device;
 *   Qimlo_d   space allocated for the Q computed by the device;
 *   Rrehi_h   space allocated for the R computed by the host;
 *   Rrelo_h   space allocated for the R computed by the host;
 *   Rimhi_h   space allocated for the R computed by the host;
 *   Rimlo_h   space allocated for the R computed by the host;
 *   Rrehi_d   space allocated for the R computed by the device;
 *   Rrelo_d   space allocated for the R computed by the device;
 *   Rimhi_d   space allocated for the R computed by the device;
 *   Rimlo_d   space allocated for the R computed by the device;
 *   workvecrehi is work space allocated for a vector of dimension dim;
 *   workvecrelo is work space allocated for a vector of dimension dim;
 *   workvecimhi is work space allocated for a vector of dimension dim;
 *   workvecimlo is work space allocated for a vector of dimension dim;
 *   resvecrehi has space for deg+1 vectors of dimension dim;
 *   resvecrelo has space for deg+1 vectors of dimension dim;
 *   resvecimhi has space for deg+1 vectors of dimension dim;
 *   resvecimlo has space for deg+1 vectors of dimension dim;
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
 *   inputrehi_h has the high doubles of the real parts of series, on host;
 *   inputrelo_h has the low doubles of the real parts of series, on host;
 *   inputimhi_h has the high doubles of the imag parts of series, on host;
 *   inputimlo_h has the low doubles of the imag parts of series, on host;
 *   inputrehi_d has the high doubles of the real parts of series, on device;
 *   inputrelo_d has the low doubles of the real parts of series, on device;
 *   inputimhi_d has the high doubles of the imag parts of series, on device;
 *   inputimlo_d has the low doubles of the imag parts of series, on device;
 *   outputrehi_h has the high doubles of the real part of output, on host;
 *   outputrelo_h has the low doubles of the real part of output, on host;
 *   outputimhi_h has the high doubles of the imag part of output, on host;
 *   outputimlo_h has the low doubles of the imag part of output, on host;
 *   outputrehi_d has the high doubles of the real part of output, on device;
 *   outputrelo_d has the low doubles of the real part of output, on device;
 *   outputimhi_d has the high doubles of the imag part of output, on device;
 *   outputimlo_d has the low doubles of the imag part of output, on device;
 *   funvalrehi_h is outputrehi[i][dim], on host;
 *   funvalrelo_h is outputrelo[i][dim], on host;
 *   funvalimhi_h is outputimhi[i][dim], on host;
 *   funvalimlo_h is outputimlo[i][dim], on host;
 *   funvalrehi_d is outputrehi[i][dim], on device;
 *   funvalrelo_d is outputrelo[i][dim], on device;
 *   funvalimhi_d is outputimhi[i][dim], on device;
 *   funvalimlo_d is outputimlo[i][dim], on device;
 *   jacvalrehi_h are the high doubles of the real parts of a matrix series,
 *             on host;
 *   jacvalrelo_h are the low doubles of the real parts of a matrix series,
 *             on host;
 *   jacvalimhi_h are the high doubles of the real parts of a matrix series,
 *             on host;
 *   jacvalimlo_h are the low doubles of the real parts of a matrix series,
 *             on host;
 *   jacvalrehi_d are the high doubles of the real parts of a matrix series,
 *             on device;
 *   jacvalrelo_d are the low doubles of the real parts of a matrix series,
 *             on device;
 *   jacvalimhi_d are the high doubles of the real parts of a matrix series,
 *             on device;
 *   jacvalimlo_d are the low doubles of the real parts of a matrix series,
 *             on device;
 *   rhsrehi_h are high doubles of real parts of the linearized right hand
 *             side, subtracted by 1 and added by t, computed by host;
 *   rhsrelo_h are low doubles of real parts of the linearized right hand,
 *             side, subtracted by 1 and added by t, computed by host;
 *   rhsimhi_h are high doubles of imag parts of the linearized right hand,
 *             side, subtracted by 1 and added by t, computed by host;
 *   rhsimlo_h are low doubles of imag parts of the linearized right hand,
 *             side, subtracted by 1 and added by t, computed by host;
 *   rhsrehi_d are high doubles of real parts of the linearized right hand 
 *             side, subtracted by 1 and added by t, computed by device;
 *   rhsrelo_d are low doubles of real parts of the linearized right hand,
 *             side, subtracted by 1 and added by t, computed by device;
 *   rhsimhi_d are high doubles of imag parts of the linearized right hand,
 *             side, subtracted by 1 and added by t, computed by device;
 *   rhsimlo_d are low doubles of imag parts of the linearized right hand
 *             side, subtracted by 1 and added by t, computed by device;
 *   urhsrehi_h are the high doubles of the real parts of right hand side
 *             updated by the host;
 *   urhsrelo_h are the high doubles of the real parts of right hand side
 *             updated by the host;
 *   urhsimhi_h are the high doubles of the imaginary parts of right hand side
 *             updated by the host;
 *   urhsimlo_h are the high doubles of the imaginary parts of right hand side
 *             updated by the host;
 *   urhsrehi_d are the high doubles of the real parts of right hand side
 *             updated by the device;
 *   urhsrelo_d are the high doubles of the real parts of right hand side
 *             updated by the device;
 *   urhsimhi_d are the high doubles of the imaginary parts of right hand side
 *             updated by the device;
 *   urhsimlo_d are the high doubles of the imaginary parts of right hand side
 *             updated by the device;
 *   solrehi_h are the high doubles of real parts of solution on host;
 *   solrelo_h are the low doubles of real parts of solution on host;
 *   solimhi_h are the high doubles of imag parts of solution on host;
 *   solimlo_h are the low doubles of imag parts of solution on host;
 *   solrehi_d are the high doubles of real parts of solution on device;
 *   solrelo_d are the low doubles of real parts of solution on device;
 *   solimhi_d are the high doubles of imag parts of solution on device;
 *   solimlo_d are the low doubles of imag parts of solution on device;
 *   Qrehi_h   high doubles of real Q of the QR on host;
 *   Qrelo_h   low doubles of real Q of the QR on host;
 *   Qimhi_h   high doubles of imaginary Q of the QR on host;
 *   Qimlo_h   low doubles of imaginary Q of the QR on host;
 *   Qrehi_d   high doubles of real Q of the QR on device;
 *   Qrelo_d   low doubles of real Q of the QR on device;
 *   Qimhi_d   high doubles of imaginary Q of the QR on device;
 *   Qimlo_d   low doubles of imaginary Q of the QR on device;
 *   Rrehi_h   high doubles of real R of the QR on host;
 *   Rrelo_h   low doubles of real R of the QR on host;
 *   Rimhi_h   high doubles of imaginary R of the QR on host;
 *   Rimlo_h   low doubles of imaginary R of the QR on host;
 *   Rrehi_d   high doubles of real R of the QR on device;
 *   Rrelo_d   low doubles of real R of the QR on device;
 *   Rimhi_d   high doubles of imaginary R of the QR on device;
 *   Rimlo_d   low doubles of imaginary R of the QR on device;
 *   resvecrehi are high doubles of the real parts of the residual vectors;
 *   resvecrelo are high doubles of the real parts of the residual vectors;
 *   resvecimhi are high doubles of the imag parts of the residual vectors;
 *   resvecimlo are high doubles of the imag parts of the residual vectors;
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
 *   totqrlapsedms accumulates the milliseconds spent on Q times rhs;
 *   totbslapsedms accumulates the milliseconds spent on back substitutions;
 *   totupdlapsedms accumulates the milliseconds spent on updates;
 *   totreslapsedms accumulates the milliseconds spent on residuals. */

int cmplx2_row_newton_qrstep
 ( int szt, int nbt, int dim, int deg, int *tailidx_h, int *tailidx_d,
   int *nbr, int **nvr, int ***idx,
   double **cstrehi, double **cstrelo, double **cstimhi, double **cstimlo,
   double ***cffrehi, double ***cffrelo, double ***cffimhi, double ***cffimlo,
   double dpr,
   double **inputrehi_h, double **inputrelo_h,
   double **inputimhi_h, double **inputimlo_h,
   double **inputrehi_d, double **inputrelo_d,
   double **inputimhi_d, double **inputimlo_d,
   double ***outputrehi_h, double ***outputrelo_h,
   double ***outputimhi_h, double ***outputimlo_h,
   double ***outputrehi_d, double ***outputrelo_d,
   double ***outputimhi_d, double ***outputimlo_d,
   double **funvalrehi_h, double **funvalrelo_h,
   double **funvalimhi_h, double **funvalimlo_h,
   double **funvalrehi_d, double **funvalrelo_d,
   double **funvalimhi_d, double **funvalimlo_d,
   double ***jacvalrehi_h, double ***jacvalrelo_h,
   double ***jacvalimhi_h, double ***jacvalimlo_h,
   double ***jacvalrehi_d, double ***jacvalrelo_d,
   double ***jacvalimhi_d, double ***jacvalimlo_d,
   double **rhsrehi_h, double **rhsrelo_h,
   double **rhsimhi_h, double **rhsimlo_h,
   double **rhsrehi_d, double **rhsrelo_d, 
   double **rhsimhi_d, double **rhsimlo_d,
   double **urhsrehi_h, double **urhsrelo_h,
   double **urhsimhi_h, double **urhsimlo_h,
   double **urhsrehi_d, double **urhsrelo_d,
   double **urhsimhi_d, double **urhsimlo_d,
   double **solrehi_h, double **solrelo_h,
   double **solimhi_h, double **solimlo_h, 
   double **solrehi_d, double **solrelo_d, 
   double **solimhi_d, double **solimlo_d,
   double **Qrehi_h, double **Qrelo_h, double **Qimhi_h, double **Qimlo_h,
   double **Qrehi_d, double **Qrelo_d, double **Qimhi_d, double **Qimlo_d, 
   double **Rrehi_h, double **Rrelo_h, double **Rimhi_h, double **Rimlo_h, 
   double **Rrehi_d, double **Rrelo_d, double **Rimhi_d, double **Rimlo_d,
   double *workvecrehi, double *workvecrelo,
   double *workvecimhi, double *workvecimlo,
   double **resvecrehi, double **resvecrelo,
   double **resvecimhi, double **resvecimlo,
   double *resmaxhi, double *resmaxlo,
   bool* zeroQ_h, bool *noqr_h, bool* zeroQ_d, bool *noqr_d,
   int *upidx_h, int *bsidx_h, int *upidx_d, int *bsidx_d,
   double *totcnvlapsedms, double *totqrlapsedms, double *totqtblapsedms,
   double *totbslapsedms, double *totupdlapsedms, double *totreslapsedms,
   int vrblvl, int mode );
/*
 * DESCRIPTION :
 *   Does one step with Newton's method to update a power series,
 *   using QR factorization to solve linear systems, on real data,
 *   on one or more columns of monomials.
 *
 * REQUIRED : szt*nbt = dim for GPU computing.
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
 *   cstrehi   high doubles of the real parts of the constants;
 *   cstrelo   low doubles of the real parts of the constants;
 *   cstimhi   high doubles of the imaginary parts of the constants;
 *   cstimlo   low doubles of the imaginary parts of the constants;
 *   dpr       damper multiplier for t, should be in (0.0, 1.0];
 *   cffrehi   high doubles of the real parts of the coefficients
 *             of the monomials;
 *   cffrelo   low doubles of the real parts of the coefficients
 *             of the monomials;
 *   cffimhi   high doubles of the imaginary parts of the coefficients
 *             of the monomials;
 *   cffimlo   low doubles of the imaginary parts of the coefficients
 *             of the monomials;
 *   inputrehi_h are the high doubles of the real parts coefficients of series
 *             of degree deg, for dim variables, computed on host;
 *   inputrelo_h are the low doubles of the real parts coefficients of series
 *             of degree deg, for dim variables, computed on host;
 *   inputimhi_h are the high doubles of the imaginary parts coefficients of
 *             the series of degree deg for dim variables, computed on host;
 *   inputimlo_h are the low doubles of the imaginary parts coefficients of
 *             the series of degree deg for dim variables, computed on host;
 *   inputrehi_d has space for power series computed on device;
 *   inputrelo_d has space for power series computed on device;
 *   inputimhi_d has space for power series computed on device;
 *   inputimlo_d has space for power series computed on device;
 *   outputrehi_h has space for the evaluated and differentiated monomials,
 *             computed on the host;
 *   outputrelo_h has space for the evaluated and differentiated monomials,
 *             computed on the host;
 *   outputimhi_h has space for the evaluated and differentiated monomials,
 *             computed on the host;
 *   outputimlo_h has space for the evaluated and differentiated monomials,
 *             computed on the host;
 *   outputrehi_d has space for the evaluated and differentiated monomials,
 *             computed on the device;
 *   outputrelo_d has space for the evaluated and differentiated monomials,
 *             computed on the device;
 *   outputimhi_d has space for the evaluated and differentiated monomials,
 *             computed on the device;
 *   outputimlo_d has space for the evaluated and differentiated monomials,
 *             computed on the device;
 *   funvalrehi_h has space for the evaluated series computed by host;
 *   funvalrelo_h has space for the evaluated series computed by host;
 *   funvalimhi_h has space for the evaluated series computed by host;
 *   funvalimlo_h has space for the evaluated series computed by host;
 *   funvalrehi_d has space for the evaluated series computed by device;
 *   funvalrelo_d has space for the evaluated series computed by device;
 *   funvalimhi_d has space for the evaluated series computed by device;
 *   funvalimlo_d has space for the evaluated series computed by device;
 *   jacvalrehi_h has space for deg+1 matrices of dimension dim on host;
 *   jacvalrelo_h has space for deg+1 matrices of dimension dim on host;
 *   jacvalimhi_h has space for deg+1 matrices of dimension dim on host;
 *   jacvalimlo_h has space for deg+1 matrices of dimension dim on host;
 *   jacvalrehi_d has space for deg+1 matrices of dimension dim on device;
 *   jacvalrelo_d has space for deg+1 matrices of dimension dim on device;
 *   jacvalimhi_d has space for deg+1 matrices of dimension dim on device;
 *   jacvalimlo_d has space for deg+1 matrices of dimension dim on device;
 *   rhsrehi_h has space for deg+1 vectors of dimension dim on host;
 *   rhsrelo_h has space for deg+1 vectors of dimension dim on host;
 *   rhsimhi_h has space for deg+1 vectors of dimension dim on host;
 *   rhsimlo_h has space for deg+1 vectors of dimension dim on host;
 *   rhsrehi_d has space for deg+1 vectors of dimension dim on device;
 *   rhsrelo_d has space for deg+1 vectors of dimension dim on device;
 *   rhsimhi_d has space for deg+1 vectors of dimension dim on device;
 *   rhsimlo_d has space for deg+1 vectors of dimension dim on device;
 *   urhsrehi_h has space for updated right hand side computed by host;
 *   urhsrelo_h has space for updated right hand side computed by host;
 *   urhsimhi_h has space for updated right hand side computed by host;
 *   urhsimlo_h has space for updated right hand side computed by host;
 *   urhsrehi_d has space for updated right hand side computed by device; 
 *   urhsrelo_d has space for updated right hand side computed by device; 
 *   urhsimhi_d has space for updated right hand side computed by device; 
 *   urhsimlo_d has space for updated right hand side computed by device; 
 *   solrehi_h has space for deg+1 vectors of dimension dim;
 *   solrelo_h has space for deg+1 vectors of dimension dim;
 *   solimhi_h has space for deg+1 vectors of dimension dim;
 *   solimlo_h has space for deg+1 vectors of dimension dim;
 *   solrehi_d has space for deg+1 vectors of dimension dim;
 *   solrelo_d has space for deg+1 vectors of dimension dim;
 *   solimhi_d has space for deg+1 vectors of dimension dim;
 *   solimlo_d has space for deg+1 vectors of dimension dim;
 *   Qrehi_h   space allocated for the Q computed by the host;
 *   Qrelo_h   space allocated for the Q computed by the host;
 *   Qimhi_h   space allocated for the Q computed by the host;
 *   Qimlo_h   space allocated for the Q computed by the host;
 *   Qrehi_d   space allocated for the Q computed by the device;
 *   Qrelo_d   space allocated for the Q computed by the device;
 *   Qimhi_d   space allocated for the Q computed by the device;
 *   Qimlo_d   space allocated for the Q computed by the device;
 *   Rrehi_h   space allocated for the R computed by the host;
 *   Rrelo_h   space allocated for the R computed by the host;
 *   Rimhi_h   space allocated for the R computed by the host;
 *   Rimlo_h   space allocated for the R computed by the host;
 *   Rrehi_d   space allocated for the R computed by the device;
 *   Rrelo_d   space allocated for the R computed by the device;
 *   Rimhi_d   space allocated for the R computed by the device;
 *   Rimlo_d   space allocated for the R computed by the device;
 *   workvecrehi is work space allocated for a vector of dimension dim;
 *   workvecrelo is work space allocated for a vector of dimension dim;
 *   workvecimhi is work space allocated for a vector of dimension dim;
 *   workvecimlo is work space allocated for a vector of dimension dim;
 *   resvecrehi has space for deg+1 vectors of dimension dim;
 *   resvecrelo has space for deg+1 vectors of dimension dim;
 *   resvecimhi has space for deg+1 vectors of dimension dim;
 *   resvecimlo has space for deg+1 vectors of dimension dim;
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
 *   inputrehi_h has the high doubles of the real parts of series, on host;
 *   inputrelo_h has the low doubles of the real parts of series, on host;
 *   inputimhi_h has the high doubles of the imag parts of series, on host;
 *   inputimlo_h has the low doubles of the imag parts of series, on host;
 *   inputrehi_d has the high doubles of the real parts of series, on device;
 *   inputrelo_d has the low doubles of the real parts of series, on device;
 *   inputimhi_d has the high doubles of the imag parts of series, on device;
 *   inputimlo_d has the low doubles of the imag parts of series, on device;
 *   outputrehi_h has the high doubles of the real part of output, on host;
 *   outputrelo_h has the low doubles of the real part of output, on host;
 *   outputimhi_h has the high doubles of the imag part of output, on host;
 *   outputimlo_h has the low doubles of the imag part of output, on host;
 *   outputrehi_d has the high doubles of the real part of output, on device;
 *   outputrelo_d has the low doubles of the real part of output, on device;
 *   outputimhi_d has the high doubles of the imag part of output, on device;
 *   outputimlo_d has the low doubles of the imag part of output, on device;
 *   funvalrehi_h is outputrehi[i][dim], on host;
 *   funvalrelo_h is outputrelo[i][dim], on host;
 *   funvalimhi_h is outputimhi[i][dim], on host;
 *   funvalimlo_h is outputimlo[i][dim], on host;
 *   funvalrehi_d is outputrehi[i][dim], on device;
 *   funvalrelo_d is outputrelo[i][dim], on device;
 *   funvalimhi_d is outputimhi[i][dim], on device;
 *   funvalimlo_d is outputimlo[i][dim], on device;
 *   jacvalrehi_h are the high doubles of the real parts of a matrix series,
 *             on host;
 *   jacvalrelo_h are the low doubles of the real parts of a matrix series,
 *             on host;
 *   jacvalimhi_h are the high doubles of the real parts of a matrix series,
 *             on host;
 *   jacvalimlo_h are the low doubles of the real parts of a matrix series,
 *             on host;
 *   jacvalrehi_d are the high doubles of the real parts of a matrix series,
 *             on device;
 *   jacvalrelo_d are the low doubles of the real parts of a matrix series,
 *             on device;
 *   jacvalimhi_d are the high doubles of the real parts of a matrix series,
 *             on device;
 *   jacvalimlo_d are the low doubles of the real parts of a matrix series,
 *             on device;
 *   rhsrehi_h are high doubles of real parts of the linearized right hand
 *             side, subtracted by 1 and added by t, computed by host;
 *   rhsrelo_h are low doubles of real parts of the linearized right hand,
 *             side, subtracted by 1 and added by t, computed by host;
 *   rhsimhi_h are high doubles of imag parts of the linearized right hand,
 *             side, subtracted by 1 and added by t, computed by host;
 *   rhsimlo_h are low doubles of imag parts of the linearized right hand,
 *             side, subtracted by 1 and added by t, computed by host;
 *   rhsrehi_d are high doubles of real parts of the linearized right hand 
 *             side, subtracted by 1 and added by t, computed by device;
 *   rhsrelo_d are low doubles of real parts of the linearized right hand,
 *             side, subtracted by 1 and added by t, computed by device;
 *   rhsimhi_d are high doubles of imag parts of the linearized right hand,
 *             side, subtracted by 1 and added by t, computed by device;
 *   rhsimlo_d are low doubles of imag parts of the linearized right hand
 *             side, subtracted by 1 and added by t, computed by device;
 *   urhsrehi_h are the high doubles of the real parts of right hand side
 *             updated by the host;
 *   urhsrelo_h are the high doubles of the real parts of right hand side
 *             updated by the host;
 *   urhsimhi_h are the high doubles of the imaginary parts of right hand side
 *             updated by the host;
 *   urhsimlo_h are the high doubles of the imaginary parts of right hand side
 *             updated by the host;
 *   urhsrehi_d are the high doubles of the real parts of right hand side
 *             updated by the device;
 *   urhsrelo_d are the high doubles of the real parts of right hand side
 *             updated by the device;
 *   urhsimhi_d are the high doubles of the imaginary parts of right hand side
 *             updated by the device;
 *   urhsimlo_d are the high doubles of the imaginary parts of right hand side
 *             updated by the device;
 *   solrehi_h are the high doubles of real parts of solution on host;
 *   solrelo_h are the low doubles of real parts of solution on host;
 *   solimhi_h are the high doubles of imag parts of solution on host;
 *   solimlo_h are the low doubles of imag parts of solution on host;
 *   solrehi_d are the high doubles of real parts of solution on device;
 *   solrelo_d are the low doubles of real parts of solution on device;
 *   solimhi_d are the high doubles of imag parts of solution on device;
 *   solimlo_d are the low doubles of imag parts of solution on device;
 *   Qrehi_h   high doubles of real Q of the QR on host;
 *   Qrelo_h   low doubles of real Q of the QR on host;
 *   Qimhi_h   high doubles of imaginary Q of the QR on host;
 *   Qimlo_h   low doubles of imaginary Q of the QR on host;
 *   Qrehi_d   high doubles of real Q of the QR on device;
 *   Qrelo_d   low doubles of real Q of the QR on device;
 *   Qimhi_d   high doubles of imaginary Q of the QR on device;
 *   Qimlo_d   low doubles of imaginary Q of the QR on device;
 *   Rrehi_h   high doubles of real R of the QR on host;
 *   Rrelo_h   low doubles of real R of the QR on host;
 *   Rimhi_h   high doubles of imaginary R of the QR on host;
 *   Rimlo_h   low doubles of imaginary R of the QR on host;
 *   Rrehi_d   high doubles of real R of the QR on device;
 *   Rrelo_d   low doubles of real R of the QR on device;
 *   Rimhi_d   high doubles of imaginary R of the QR on device;
 *   Rimlo_d   low doubles of imaginary R of the QR on device;
 *   resvecrehi are high doubles of the real parts of the residual vectors;
 *   resvecrelo are high doubles of the real parts of the residual vectors;
 *   resvecimhi are high doubles of the imag parts of the residual vectors;
 *   resvecimlo are high doubles of the imag parts of the residual vectors;
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
 *   totqrlapsedms accumulates the milliseconds spent on Q times rhs;
 *   totbslapsedms accumulates the milliseconds spent on back substitutions;
 *   totupdlapsedms accumulates the milliseconds spent on updates;
 *   totreslapsedms accumulates the milliseconds spent on residuals. */

int cmplx2_allocate_inoutfunjac
 ( int dim, int deg, int mode,
   double **inputrehi_h, double **inputrelo_h,
   double **inputimhi_h, double **inputimlo_h,
   double **inputrehi_d, double **inputrelo_d,
   double **inputimhi_d, double **inputimlo_d,
   double ***outputrehi_h, double ***outputrelo_h,
   double ***outputimhi_h, double ***outputimlo_h,
   double ***outputrehi_d, double ***outputrelo_d,
   double ***outputimhi_d, double ***outputimlo_d,
   double **funvalrehi_h, double **funvalrelo_h,
   double **funvalimhi_h, double **funvalimlo_h,
   double **funvalrehi_d, double **funvalrelo_d,
   double **funvalimhi_d, double **funvalimlo_d,
   double ***jacvalrehi_h, double ***jacvalrelo_h,
   double ***jacvalimhi_h, double ***jacvalimlo_h,
   double ***jacvalrehi_d, double ***jacvalrelo_d,
   double ***jacvalimhi_d, double ***jacvalimlo_d );
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
 *   inputrehi_h are the high doubles of the real parts
 *              of the input on the host;
 *   inputrelo_h are the low doubles of the real parts
 *              of the input on the host;
 *   inputimhi_h are the high doubles of the imaginary parts
 *              of the input on the host;
 *   inputimlo_h are the low doubles of the imaginary parts
 *              of the input on the host;
 *   inputrehi_d are the high doubles of the real parts
 *              of the input on the device;
 *   inputrelo_d are the low doubles of the real parts
 *              of the input on the device;
 *   inputimhi_d are the high doubles of the imaginary parts
 *              of the input on the device;
 *   inputimlo_d are the low doubles of the imaginary parts
 *              of the input on the device;
 *   outputrehi_h are the high doubles of the real parts
 *              of the output on the host;
 *   outputrelo_h are the low doubles of the real parts
 *              of the output on the host;
 *   outputimhi_h are the high doubles of the imaginary parts
 *              of the output on the host;
 *   outputimlo_h are the low doubles of the imaginary parts
 *              of the output on the host;
 *   outputrehi_d are the high doubles of the real parts
 *              of the output on the device;
 *   outputrelo_d are the low doubles of the real parts
 *              of the output on the device;
 *   outputimhi_d are the high doubles of the imaginary parts 
 *              of the output on the device;
 *   outputimlo_d are the low doubles of the imaginary parts 
 *              of the output on the device;
 *   funvalrehi_h are the high doubles of the real parts
 *              of the function values on the host;
 *   funvalrelo_h are the low doubles of the real parts
 *              of the function values on the host;
 *   funvalimhi_h are the high doubles of the imaginary parts
 *              of the function values on the host;
 *   funvalimlo_h are the low doubles of the imaginary parts
 *              of the function values on the host;
 *   funvalrehi_d are the high doubles of the real parts
 *              of the function values on the device;
 *   funvalrelo_d are the low doubles of the real parts
 *              of the function values on the device;
 *   funvalimhi_d are the high doubles of the imaginary parts
 *              of the function values on the device;
 *   funvalimlo_d are the low doubles of the imaginary parts
 *              of the function values on the device;
 *   jacvalrehi_h are the high doubles of the real parts
 *              of the Jacobian matrix on the host;
 *   jacvalrelo_h are the low doubles of the real parts
 *              of the Jacobian matrix on the host;
 *   jacvalimhi_h are the high doubles of the imaginary parts
 *              of the Jacobian matrix on the host;
 *   jacvalimlo_h are the low doubles of the imaginary parts
 *              of the Jacobian matrix on the host;
 *   jacvalrehi_d are the high doubles of the real parts
 *              of the Jacobian matrix on the device;
 *   jacvalrelo_d are the low doubles of the real parts
 *              of the Jacobian matrix on the device;
 *   jacvalimhi_d are the high doubles of the imaginary parts
 *              of the Jacobian matrix on the device;
 *   jacvalimlo_d are the low doubles of the imaginary parts
 *              of the Jacobian matrix on the device. */

int cmplx2_allocate_rhsqrsol
 ( int dim, int deg, int mode,
   double **rhsrehi_h, double **rhsrelo_h,
   double **rhsimhi_h, double **rhsimlo_h,
   double **rhsrehi_d, double **rhsrelo_d,
   double **rhsimhi_d, double **rhsimlo_d,
   double **urhsrehi_h, double **urhsrelo_h, 
   double **urhsimhi_h, double **urhsimlo_h,
   double **urhsrehi_d, double **urhsrelo_d, 
   double **urhsimhi_d, double **urhsimlo_d,
   double **Qrehi_h, double **Qrelo_h, double **Qimhi_h, double **Qimlo_h, 
   double **Qrehi_d, double **Qrelo_d, double **Qimhi_d, double **Qimlo_d,
   double **Rrehi_h, double **Rrelo_h, double **Rimhi_h, double **Rimlo_h, 
   double **Rrehi_d, double **Rrelo_d, double **Rimhi_d, double **Rimlo_d,
   double **solrehi_h, double **solrelo_h,
   double **solimhi_h, double **solimlo_h,
   double **solrehi_d, double **solrelo_d,
   double **solimhi_d, double **solimlo_d );
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
 *   rhsrehi_h are the high doubles of real parts
 *              of the right hand side vector on the host;
 *   rhsrelo_h are the low doubles of real parts
 *              of the right hand side vector on the host;
 *   rhsimhi_h are the high doubles of imaginary parts
 *              of the right hand side vector on the host;
 *   rhsimlo_h are the low doubles of imaginary parts
 *              of the right hand side vector on the host;
 *   rhsrehi_d are the high doubles of real parts
 *              of the right hand side vector on the device;
 *   rhsrelo_d are the low doubles of real parts
 *              of the right hand side vector on the device;
 *   rhsimhi_d are the high doubles of imaginary parts
 *              of the right hand side vector on the device;
 *   rhsimlo_d are the low doubles of imaginary parts
 *              of the right hand side vector on the device;
 *   urhsrehi_h are the high doubles of the real parts
 *              of the updated right hand side on the host;
 *   urhsrelo_h are the low doubles of the real parts
 *              of the updated right hand side on the host;
 *   urhsimhi_h are the high doubles of the imaginary parts
 *              of the updated right hand side on the host;
 *   urhsimlo_h are the low doubles of the imaginary parts
 *              of the updated right hand side on the host;
 *   urhsrehi_d are the high doubles of the real parts
 *              of the updated right hand side on the device;
 *   urhsrelo_d are the low doubles of the real parts
 *              of the updated right hand side on the device;
 *   urhsimhi_d are the high doubles of the imaginary parts
 *              of the updated right hand side on the device;
 *   urhsimlo_d are the low doubles of the imaginary parts
 *              of the updated right hand side on the device;
 *   Qrehi_h    high doubles of the real parts of Q on the host;
 *   Qrelo_h    low doubles of the real parts of Q on the host;
 *   Qimhi_h    high doubles of the imaginary parts of Q on the host;
 *   Qimlo_h    low doubles of the imaginary parts of Q on the host;
 *   Qrehi_d    high doubles of the real parts of Q on the device;
 *   Qrelo_d    low doubles of the real parts of Q on the device;
 *   Qimhi_d    high doubles of the imaginary parts of Q on the device;
 *   Qimlo_d    low doubles of the imaginary parts of Q on the device;
 *   Rrehi_h    high doubles of the real parts of R on the host;
 *   Rrelo_h    low doubles of the real parts of R on the host;
 *   Rimhi_h    high doubles of the imaginary parts of R on the host;
 *   Rimlo_h    low doubles of the imaginary parts of R on the host;
 *   Rrehi_d    high doubles of the real parts of R on the device;
 *   Rrelo_d    low doubles of the real parts of R on the device;
 *   Rimhi_d    high doubles of the imaginary parts of R on the device;
 *   Rimlo_d    low doubles of the imaginary parts of R on the device;
 *   solrehi_h are the high doubles of the real parts
 *              of the update to the solution on the host,
 *   solrelo_h are the lowh doubles of the real parts
 *              of the update to the solution on the host,
 *   solimhi_h are the high doubles of the imaginary parts
 *              of the update to the solution on the host,
 *   solimlo_h are the low doubles of the imaginary parts
 *              of the update to the solution on the host,
 *   solrehi_d are the high doubles of the real parts
 *              of the update to the solution on the device;
 *   solrelo_d are the low doubles of the real parts
 *              of the update to the solution on the device;
 *   solimhi_d are the high doubles of the imaginary parts
 *              of the update to the solution on the device;
 *   solimlo_d are the low doubles of the imaginary parts
 *              of the update to the solution on the device. */

void cmplx2_start_setup
 ( int dim, int deg,
   double **testsolrehi, double **testsolrelo,
   double **testsolimhi, double **testsolimlo,
   double **inputrehi_h, double **inputrelo_h,
   double **inputimhi_h, double **inputimlo_h,
   double **inputrehi_d, double **inputrelo_d,
   double **inputimhi_d, double **inputimlo_d, int mode, int vrblvl );
/*
 * DESCRIPTION :
 *   Given the test solution, defines the start vector.
 *
 * ON ENTRY :
 *   dim        dimension of the solution vectors;
 *   deg        degree at which the series are truncated;
 *   testsolrehi has the high doubles of the real parts of the test solution;
 *   testsolrelo has the low doubles of the real parts of the test solution;
 *   testsolimhi has the high doubles of the imag parts of the test solution;
 *   testsolimlo has the high doubles of the imag parts of the test solution;
 *   inputrehi_h is allocated on host if mode is 1 or 2;
 *   inputrelo_h is allocated on host if mode is 1 or 2;
 *   inputimhi_h is allocated on host if mode is 1 or 2;
 *   inputimlo_h is allocated on host if mode is 1 or 2;
 *   inputrehi_d is allocated on device if mode is 0 or 2;
 *   inputrelo_d is allocated on device if mode is 0 or 2;
 *   inputimhi_d is allocated on device if mode is 0 or 2;
 *   inputimlo_d is allocated on device if mode is 0 or 2;
 *   mode       0 (GPU), 1 (CPU), or 2 (GPU+CPU);
 *   vrblvl     is the verbose level.
 *
 * ON RETURN :
 *   inputrehi_h are the high doubles of the real parts of start vector
 *              on host if mode is 1 or 2;
 *   inputrelo_h are the low doubles of the real parts of start vector
 *              on host if mode is 1 or 2;
 *   inputimhi_h are the high doubles of the imaginary parts of start vector
 *              on host if mode is 1 or 2;
 *   inputimlo_h are the low doubles of the imaginary parts of start vector
 *              on host if mode is 1 or 2;
 *   inputrehi_d are the high doubles of the real parts of start vector
 *              on device if mode is 0 or 2;
 *   inputrelo_d are the low doubles of the real parts of start vector
 *              on device if mode is 0 or 2;
 *   inputimhi_d are the high doubles of the imaginary parts start vector
 *              on device if mode is 0 or 2.
 *   inputimlo_d are the low doubles of the imaginary parts start vector
 *              on device if mode is 0 or 2. */

void cmplx2_column_setup
 ( int dim, int deg, int nbrcol,
   int **nvr, int ***idx, int **rowsA,
   double ***cffrehi, double ***cffrelo,
   double ***cffimhi, double ***cffimlo,
   double **testsolrehi, double **testsolrelo,
   double **testsolimhi, double **testsolimlo,
   double **mbrhsrehi, double **mbrhsrelo,
   double **mbrhsimhi, double **mbrhsimlo,
   double **inputrehi_h, double **inputrelo_h,
   double **inputimhi_h, double **inputimlo_h,
   double **inputrehi_d, double **inputrelo_d,
   double **inputimhi_d, double **inputimlo_d, int mode, int vrblvl );
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
 *   cffrehi    high doubles of real parts of coefficients,
 *              if more than one column;
 *   cffrelo    low doubles of real parts of coefficients,
 *              if more than one column;
 *   cffimhi    high doubles of imaginary parts of coefficients,
 *              if more than one column;
 *   cffimlo    low doubles of imaginary parts of coefficients,
 *              if more than one column;
 *   testsolrehi has space for dim pointers;
 *   testsolrelo has space for dim pointers;
 *   testsolimhi has space for dim pointers;
 *   testsolimlo has space for dim pointers;
 *   mbrhsrehi  has space for dim pointers;
 *   mbrhsrelo  has space for dim pointers;
 *   mbrhsimhi  has space for dim pointers;
 *   mbrhsimlo  has space for dim pointers;
 *   inputrehi_h allocated on host if mode is 1 or 2;
 *   inputrelo_h allocated on host if mode is 1 or 2;
 *   inputimhi_h allocated on host if mode is 1 or 2;
 *   inputimlo_h allocated on host if mode is 1 or 2;
 *   inputrehi_d allocated on device if mode is 0 or 2;
 *   inputrelo_d allocated on device if mode is 0 or 2;
 *   inputimhi_d allocated on device if mode is 0 or 2;
 *   inputimlo_d allocated on device if mode is 0 or 2;
 *   mode       0 (GPU), 1 (CPU), or 2 (GPU+CPU);
 *   vrblvl     is the verbose level.
 *
 * ON RETURN :
 *   testsolrehi has the high doubles of the real parts of the test solution;
 *   testsolrelo has the low doubles of the real parts of the test solution;
 *   testsolimhi has the high doubles of the imag parts of the test solution;
 *   testsolimlo has the low doubles of the imag parts of the test solution;
 *   mbrhsrehi  has the high doubles of the real parts
 *              of right hand side vector for the test solution;
 *   mbrhsrelo  has the low doubles of the real parts
 *              of right hand side vector for the test solution;
 *   mbrhsimhi  has the high doubles of the imag parts
 *              of right hand side vector for the test solution;
 *   mbrhsimlo  has the low doubles of the imag parts
 *              of right hand side vector for the test solution;
 *   inputrehi_h has the high doubles of the real parts
 *              of the start vector for host if mode is 1 or 2;
 *   inputrelo_h has the low doubles of the real parts
 *              of the start vector for host if mode is 1 or 2;
 *   inputimhi_h has the high doubles of the imaginary parts
 *              of the start vector for host if mode is 1 or 2;
 *   inputimlo_h has the low doubles of the imaginary parts
 *              of the start vector for host if mode is 1 or 2;
 *   inputrehi_d has the high doubles of the real parts
 *              of the start vector for device if mode is 0 or 2;
 *   inputrelo_d has the low doubles of the real parts
 *              of the start vector for device if mode is 0 or 2;
 *   inputimhi_d has the high doubles of the imaginary parts
 *              of the start vector for device if mode is 0 or 2;
 *   inputimlo_d has the low doubles of the imaginary parts
 *              of the start vector for device if mode is 0 or 2. */

void cmplx2_row_setup
 ( int dim, int deg, int *nbr, int **nvr, int ***idx,
   double **cstrehi, double **cstrelo, double **cstimhi, double **cstimlo,
   double ***cffrehi, double ***cffrelo,
   double ***cffimhi, double ***cffimlo,
   double **testsolrehi, double **testsolrelo,
   double **testsolimhi, double **testsolimlo,
   double **inputrehi_h, double **inputrelo_h,
   double **inputimhi_h, double **inputimlo_h,
   double **inputrehi_d, double **inputrelo_d,
   double **inputimhi_d, double **inputimlo_d,
   double ***outputrehi_h, double ***outputrelo_h,
   double ***outputimhi_h, double ***outputimlo_h,
   double ***outputrehi_d, double ***outputrelo_d,
   double ***outputimhi_d, double ***outputimlo_d, int mode, int vrblvl );
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
 *   cstrehi   high doubles of the real parts of the constants;
 *   cstrelo   low doubles of the real parts of the constants;
 *   cstimhi   high doubles of the imaginary parts of the constants;
 *   cstimlo   low doubles of the imaginary parts of the constants;
 *   cffrehi    high doubles of real parts of coefficients,
 *              if more than one column;
 *   cffrelo    low doubles of real parts of coefficients,
 *              if more than one column;
 *   cffimhi    high doubles of imaginary parts of coefficients,
 *              if more than one column;
 *   cffimlo    low doubles of imaginary parts of coefficients,
 *              if more than one column;
 *   testsolrehi has space for dim pointers;
 *   testsolrelo has space for dim pointers;
 *   testsolimhi has space for dim pointers;
 *   testsolimlo has space for dim pointers;
 *   inputrehi_h allocated on host if mode is 1 or 2;
 *   inputrelo_h allocated on host if mode is 1 or 2;
 *   inputimhi_h allocated on host if mode is 1 or 2;
 *   inputimlo_h allocated on host if mode is 1 or 2;
 *   inputrehi_d allocated on device if mode is 0 or 2;
 *   inputrelo_d allocated on device if mode is 0 or 2;
 *   inputimhi_d allocated on device if mode is 0 or 2;
 *   inputimlo_d allocated on device if mode is 0 or 2;
 *   mode       0 (GPU), 1 (CPU), or 2 (GPU+CPU);
 *   vrblvl     is the verbose level.
 *
 * ON RETURN :
 *   cstrehi    high doubles of the real parts of the constant,
 *              adjusted for the test solution;
 *   cstrelo    low doubles of the real parts of the constant,
 *              adjusted for the test solution;
 *   cstimhi    high doubles of the imag parts of the constant,
 *              adjusted for the test solution;
 *   cstimlo    low doubles of the imag parts of the constant,
 *              adjusted for the test solution;
 *   testsolrehi has the high doubles of the real parts
 *              of the test solution;
 *   testsolrehi has the low doubles of the real parts
 *              of the test solution;
 *   testsolimhi has the high doubles of the imaginary parts
 *              of the test solution;
 *   testsolimhi has the low doubles of the imaginary parts
 *              of the test solution;
 *   inputrehi_h has the high doubles of the real parts
 *              of the start vector for host if mode is 1 or 2;
 *   inputrelo_h has the low doubles of the real parts
 *              of the start vector for host if mode is 1 or 2;
 *   inputimhi_h has the high doubles of the imaginary parts
 *              of the start vector for host if mode is 1 or 2;
 *   inputimlo_h has the low doubles of the imaginary parts
 *              of the start vector for host if mode is 1 or 2;
 *   inputrehi_d has the high doubles of the real parts
 *              of the start vector for device if mode is 0 or 2;
 *   inputrelo_d has the low doubles of the real parts
 *              of the start vector for device if mode is 0 or 2;
 *   inputimhi_d has the high doubles of the imaginary parts
 *              of the start vector for device if mode is 0 or 2;
 *   inputimlo_d has the low doubles of the imaginary parts
 *              of the start vector for device if mode is 0 or 2;
 *   outputrehi_h are high doubles of the real parts
 *              of the evaluated test solution if mode is 1;
 *   outputrelo_h are low doubles of the real parts
 *              of the evaluated test solution if mode is 1;
 *   outputimhi_h are the high doubles of the imaginary parts
 *              of the evaluated test solution if mode is 1;
 *   outputimlo_h are the low doubles of the imaginary parts
 *              of the evaluated test solution if mode is 1;
 *   outputrehi_d are the high doubles of the real parts
 *              of the evaluated test solution, if mode is 0 or 2;
 *   outputrelo_d are the low doubles of the real parts
 *              of the evaluated test solution, if mode is 0 or 2;
 *   outputimhi_d are the high doubles of the imaginary parts
 *              of the evaluated test solution, if mode is 0 or 2;
 *   outputimlo_d are the low doubles of the imaginary parts
 *              of the evaluated test solution, if mode is 0 or 2. */

int cmplx2_error_testsol
 ( int dim, int deg, int mode,
   double **testsolrehi, double **testsolrelo,
   double **testsolimhi, double **testsolimlo,
   double **inputrehi_h, double **inputrelo_h,
   double **inputimhi_h, double **inputimlo_h,
   double **inputrehi_d, double **inputrelo_d,
   double **inputimhi_d, double **inputimlo_d );
/*
 * DESCRIPTION :
 *   Compares the computed solution against the test solution.
 *
 * ON ENTRY :
 *   dim        dimension of the solution vectors;
 *   deg        degree at which the series are truncated;
 *   mode       0 (GPU), 1 (CPU), or 2 (GPU+CPU).
 *   testsolrehi has the high doubles of the real parts of the test solution;
 *   testsolrelo has the low doubles of the real parts of the test solution;
 *   testsolimhi has the high doubles of the imag parts of the test solution;
 *   testsolimlo has the low doubles of the imag parts of the test solution;
 *   inputrehi_h has the high doubles of the real parts
 *              of the solution on the host if mode is 1 or 2;
 *   inputrelo_h has the low doubles of the real parts
 *              of the solution on the host if mode is 1 or 2;
 *   inputimhi_h has the high doubles of the imaginary parts
 *              of the solution on the host if mode is 1 or 2;
 *   inputimlo_h has the low doubles of the imaginary parts
 *              of the solution on the host if mode is 1 or 2;
 *   inputrehi_d has the high doubles of the real parts
 *              of the solution on the device if mode is 0 or 2;
 *   inputrelo_d has the low doubles of the real parts
 *              of the solution on the device if mode is 0 or 2;
 *   inputimhi_d has the high doubles of the imaginary parts
 *              of the solution on the device if mode is 0 or 2;
 *   inputimlo_d has the low doubles of the imaginary parts
 *              of the solution on the evice if mode is 0 or 2. */

int test_cmplx2_column_newton
 ( int szt, int nbt, int dim, int deg, int nbrcol,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac, int **rowsA,
   double dpr, int nbsteps, int mode, int vrblvl );
/*
 * DESCRIPTION :
 *   Runs Newton with complex double double arithmetic,
 *   on one or more columns of monomials.
 *
 * ON ENTRY :
 *   szt       size of each tile and block;
 *   nbt       number of tiles and number of blocks;
 *   dim       number of monomials;
 *   deg       degree of the power series;
 *   nbrcol    is the number of columns;
 *   nvr       nvr[i][j] is the number of variables in the j-th monomial
 *             of the i-th column;
 *   idx       idx[i][j] are the indices of the variables in monomial j
 *             of the i-th column;
 *   exp       exp[i] are the exponents of the variables in monomial i;
 *   nbrfac    nbrfac[i] are the number of exponents > 1 in monomial i;
 *   expfac    expfac[i] are the exponents in the i-th polynomial
 *             that are larger than one, minus one in the factor,
 *             if exp[i][k] > 1, then expfac[i][k] = exp[i][k] - 1;
 *   rowsA     rows of exponents of the dim monomials;
 *   dpr       damper multiplier for t, should be in (0.0, 1.0];
 *   nbsteps   the number of Newton steps;
 *   mode      the mode of execution, 0 for GPU only, 1 for CPU only,
 *             2 for CPU and GPU;
 *   vrblvl    is the verbose level. */

int test_cmplx2_row_newton
 ( int szt, int nbt, int dim, int deg, int *nbr, int **nvr, int ***idx,
   double dpr, int nbsteps, int mode, int vrblvl );
/*
 * DESCRIPTION :
 *   Runs Newton with complex double double arithmetic,
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
