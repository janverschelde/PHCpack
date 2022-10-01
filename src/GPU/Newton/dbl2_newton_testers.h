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

void cmplx2_unit_series_vector
 ( int dim, int deg,
   double **cffrehi, double **cffrelo, double **cffimhi, double **cffimlo );
/*
 * DESCRIPTION :
 *   Given space allocated in the coefficient arrays,
 *   returns the values for a complex unit series. */

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

void cmplx2_update_series
 ( int dim, int degp1,
   double **xrehi, double **xrelo, double **ximhi, double **ximlo,
   double **dxrehi, double **dxrelo, double **dximhi, double **dximlo,
   int vrblvl );
/*
 * DESCRIPTION :
 *   Adds the series in dx to x.
 *
 * ON ENTRY :
 *   dim       number of series in x and dx;
 *   degp1     degree plus one of the series;
 *   xrehi     high doubles of the real parts of series, not linearized;
 *   xrelo     low doubles of the real parts of series, not linearized;
 *   ximhi     high doubles of the imaginary parts of series, not linearized;
 *   ximlo     low doubles of the imaginary parts of series, not linearized;
 *   dxrehi    high doubles of the real parts of the linearized update;
 *   dxrelo    low doubles of the real parts of the linearized update;
 *   dximhi    high doubles of the imaginary parts of the linearized update;
 *   dximlo    low doubles of the imaginary parts of the linearized update;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   xrehi     high doubles of the real parts of the updated series;
 *   xrelo     low doubles of the real parts of the updated series;
 *   ximhi     high doubles fo the imaginary parts of the updated series;
 *   ximlo     low doubles fo the imaginary parts of the updated series. */

double dbl2_error3sum
 ( int dim1, int dim2, int dim3,
   double ***datahi_h, double ***datalo_h,
   double ***datahi_d, double ***datalo_d,
   std::string banner, int vrblvl );
/*
 * DESCRIPTION :
 *   Returns the sum of the absolute values of the differences
 *   between host and device real data, in three dimensions.
 *
 * ON ENTRY :
 *   dim1      first dimension;
 *   dim2      second dimension;
 *   dim3      third dimension;
 *   datahi_h  high doubles of data computed on the host;
 *   datalo_h  low doubles of data computed on the host;
 *   datahi_d  high doubles of data computed on the device;
 *   datalo_d  low doubles of data computed on the device;
 *   banner    string for printing if verbose level vrblvl > 1;
 *   vrblvl    is the verbose level. */

double cmplx2_error3sum
 ( int dim1, int dim2, int dim3,
   double ***datarehi_h, double ***datarelo_h,
   double ***dataimhi_h, double ***dataimlo_h,
   double ***datarehi_d, double ***datarelo_d,
   double ***dataimhi_d, double ***dataimlo_d,
   std::string banner, int vrblvl );
/*
 * DESCRIPTION :
 *   Returns the sum of the absolute values of the differences
 *   between host and device complex data, in three dimensions.
 *
 * ON ENTRY :
 *   dim1      first dimension;
 *   dim2      second dimension;
 *   dim3      third dimension;
 *   datarehi_h are high doubles of real parts of data, on host;
 *   datarelo_h are low doubles of real parts of data, on host;
 *   dataimhi_h are high doubles of imaginary parts of data, on host;
 *   dataimlo_h are low doubles of imaginary parts of data, on host;
 *   datarehi_d are high doubles of real parts of data, on device;
 *   datarelo_d are low doubles of real parts of data, on device;
 *   dataimhi_d are high doubles of imaginary parts of data, on device;
 *   dataimlo_d are low doubles of imaginary parts of data, on device;
 *   banner    string for printing if verbose level vrblvl > 1;
 *   vrblvl    is the verbose level. */

double dbl2_error2sum
 ( int nrows, int ncols,
   double **datahi_h, double **datalo_h,
   double **datahi_d, double **datalo_d,
   std::string banner, int vrblvl );
/*
 * DESCRIPTION :
 *   Returns the sum of the absolute values of the differences
 *   between host and device real data, in two dimensions.
 *
 * ON ENTRY :
 *   nrows     number of rows of the matrices;
 *   ncols     number of columns of the matrices;
 *   datahi_h  high doubles of data computed on the host;
 *   datalo_h  low doubles of data computed on the host;
 *   datahi_d  high doubles of data computed on the device;
 *   datalo_d  low doubles of data computed on the device;
 *   banner    string for printing if verbose level vrblvl > 1;
 *   vrblvl    is the verbose level. */

double cmplx2_error2sum
 ( int nrows, int ncols,
   double **datarehi_h, double **datarelo_h,
   double **dataimhi_h, double **dataimlo_h,
   double **datarehi_d, double **datarelo_d,
   double **dataimhi_d, double **dataimlo_d,
   std::string banner, int vrblvl );
/*
 * DESCRIPTION :
 *   Returns the sum of the absolute values of the differences
 *   between host and device complex data, in two dimensions.
 *
 * ON ENTRY :
 *   nrows     number of rows of the matrices;
 *   ncols     number of columns of the matrices;
 *   datarehi_h are high doubles of real parts of data, on host;
 *   datarelo_h are low doubles of real parts of data, on host;
 *   dataimhi_h are high doubles of imaginary parts of data, on host;
 *   dataimlo_h are low doubles of imaginary parts of data, on host;
 *   datarehi_d are high doubles of real parts of data, on device;
 *   datarelo_d are low doubles of real parts of data, on device;
 *   dataimhi_d are high doubles of imaginary parts of data, on device;
 *   dataimlo_d are low doubles of imaginary parts of data, on device;
 *   banner    string for printing if verbose level vrblvl > 1;
 *   vrblvl    is the verbose level. */

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

void cmplx2_newton_qrstep
 ( int szt, int nbt, int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
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
   double **workmatrehi, double **workmatrelo,
   double **workmatimhi, double **workmatimlo,
   double *workvecrehi, double *workvecrelo,
   double *workvecimhi, double *workvecimlo,
   double **resvecrehi, double **resvecrelo,
   double **resvecimhi, double **resvecimlo,
   double *resmaxhi, double *resmaxlo, int vrblvl, int mode );
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
 *   wrkmatrehi is work space allocated for a matrix of dimension dim;
 *   wrkmatrelo is work space allocated for a matrix of dimension dim;
 *   wrkmatimhi is work space allocated for a matrix of dimension dim;
 *   wrkmatimlo is work space allocated for a matrix of dimension dim;
 *   wrkvecrehi is work space allocated for a vector of dimension dim;
 *   wrkvecrelo is work space allocated for a vector of dimension dim;
 *   wrkvecimhi is work space allocated for a vector of dimension dim;
 *   wrkvecimlo is work space allocated for a vector of dimension dim;
 *   resvecrehi has space for deg+1 vectors of dimension dim;
 *   resvecrelo has space for deg+1 vectors of dimension dim;
 *   resvecimhi has space for deg+1 vectors of dimension dim;
 *   resvecimlo has space for deg+1 vectors of dimension dim;
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
 *   wrkmat    has a copy of the Jacobian matrix;
 *   resvecrehi are high doubles of the real parts of the residual vectors;
 *   resvecrelo are high doubles of the real parts of the residual vectors;
 *   resvecimhi are high doubles of the imag parts of the residual vectors;
 *   resvecimlo are high doubles of the imag parts of the residual vectors;
 *   resmaxhi  high double of the maximum element of the residual vectors;
 *   resmaxlo  low double of the maximum element of the residual vectors. */

int test_dbl2_real_newton
 ( int szt, int nbt, int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   int nbsteps, int mode, int vrblvl );
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
 *   nbsteps   the number of Newton steps;
 *   mode      the mode of execution, 0 for GPU only, 1 for CPU only,
 *             2 for CPU and GPU;
 *   vrblvl    is the verbose level. */

int test_dbl2_complex_newton
 ( int szt, int nbt, int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   int nbsteps, int mode, int vrblvl );
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
 *   nbsteps   the number of Newton steps;
 *   mode      the mode of execution, 0 for GPU only, 1 for CPU only,
 *             2 for CPU and GPU;
 *   vrblvl    is the verbose level. */

#endif
