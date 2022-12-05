// The file cmplx_newton_method.h specifies Newton's method
// on series in double precision on complex numbers.

#ifndef __cmplx_newton_method_h__
#define __cmplx_newton_method_h__

void cmplx_newton_qrstep
 ( int szt, int nbt, int dim, int deg, int nbrcol,
   int *tailidx_h, int *tailidx_d,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac,
   double **mbre, double **mbim, double dpr,
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
   double *workvecre, double *workvecim,
   double **resvecre, double **resvecim, double *resmax,
   bool *noqr_h, bool *noqr_d,
   int *upidx_h, int *bsidx_h, int *upidx_d, int *bsidx_d,
   int vrblvl, int mode );
/*
 * DESCRIPTION :
 *   Does one step with Newton's method to update a power series,
 *   using QR factorization to solve linear systems, on complex data.
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
 *   mbre      real parts of the right hand side vector of series
 *   mbim      imaginary parts of the right hand side vector of series;
 *   dpr       damper multiplier for t, should be in (0.0, 1.0];
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
 *   wrkvecre  work space allocated for a vector of dimension dim;
 *   wrkvecim  work space allocated for a vector of dimension dim;
 *   resvecre  space for deg+1 vectors of dimension dim;
 *   resvecim  space for deg+1 vectors of dimension dim;
 *   noqr_h    flag if true, then no qr on host;
 *   noqr_d    flag if true, then no qr on device;
 *   vrblvl    is the verbose level;
 *   mode      execution mode, 0 (GPU only), 1 (CPU only) or 2 (GPU+CPU).
 *
 * ON RETURN :
 *   tailidx_d is the new value for tailidx_d;
 *   tailidx_h is the new value for tailidx_h;
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
 *   resvecre  real parts of the residual vectors;
 *   resvecim  imaginary parts of the residual vectors;
 *   resmax    the maximum element of the residual vectors;
 *   noqr_h    updated flag if ||dx_0|| is zero for the first time on host;
 *   noqr_d    updated flag if ||dx_0|| is zero for the first time on device;
 *   upidx_h   counts the number of updates skipped by host;
 *   bsidx_h   counts the number of backsubstitutions skipped by host;
 *   upidx_d   counts the number of updates skipped by device;
 *   bsidx_d   counts the number of backsubstitutions skipped by device. */

int test_dbl_complex_newton
 ( int szt, int nbt, int dim, int deg, int nbrcol,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac, int **rowsA,
   double dpr, int nbsteps, int mode, int vrblvl );
/*
 * DESCRIPTION :
 *   Runs Newton on a monomial system with complex double arithmetic.
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
 *   rowsA     rows of the exponents of the dim monomials;
 *   dpr       damper multiplier for t, should be in (0.0, 1.0];
 *   nbsteps   the number of Newton steps;
 *   mode      the mode of execution, 0 for GPU only, 1 for CPU only,
 *             2 for CPU and GPU;
 *   vrblvl    is the verbose level. */

#endif
