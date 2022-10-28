// The file cmplx2_newton_method.h specifies Newton's method
// on series in double precision on complex numbers.

#ifndef __cmplx2_newton_method_h__
#define __cmplx2_newton_method_h__

void cmplx2_newton_qrstep
 ( int szt, int nbt, int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **mbrehi, double **mbrelo, double **mbimhi, double **mbimlo,
   double dpr,
   double **cffrehi, double **cffrelo, double **cffimhi, double **cffimlo,
   double *accrehi, double *accrelo, double *accimhi, double *accimlo,
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
   bool *noqr_h, bool *noqr_d,
   int *upidx_h, int *bsidx_h, int *upidx_d, int *bsidx_d,
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
 *   noqr_h    flag if true, then no qr on host;
 *   noqr_d    flag if true, then no qr on device;
 *   vrblvl    is the verbose level;
 *   mode      execution mode, 0 (GPU only), 1 (CPU only) or 2 (GPU+CPU).
 *
 * ON RETURN :
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
 *   noqr_h    updated flag if ||dx_0|| is zero for the first time on host;
 *   noqr_d    updated flag if ||dx_0|| is zero for the first time on device;
 *   upidx_h   counts the number of updates skipped by host;
 *   bsidx_h   counts the number of backsubstitutions skipped by host;
 *   upidx_d   counts the number of updates skipped by device;
 *   bsidx_d   counts the number of backsubstitutions skipped by device. */

int test_dbl2_complex_newton
 ( int szt, int nbt, int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac, int **rowsA,
   double dpr, int nbsteps, int mode, int vrblvl );
/*
 * DESCRIPTION :
 *   Runs Newton on a monomial system with complex double double arithmetic.
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
 *   rowsA     rows of exponents of the dim monomials;
 *   dpr       damper multiplier for t, should be in (0.0, 1.0];
 *   nbsteps   the number of Newton steps;
 *   mode      the mode of execution, 0 for GPU only, 1 for CPU only,
 *             2 for CPU and GPU;
 *   vrblvl    is the verbose level. */

#endif
