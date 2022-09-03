// The file dbl4_newton_testers.h specifies test function for Newton's method
// on series in quad double precision.

#ifndef __dbl4_newton_testers_h__
#define __dbl4_newton_testers_h__

void dbl4_unit_series_vector
 ( int dim, int deg,
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo );
/*
 * DESCRIPTION :
 *   Given space in cffhihi, cfflohi, cffhilo, cfflolo for the doubles
 *   for the vector of dim power series in quad double precision,
 *   of series truncated after degree deg,
 *   initializes the coefficients in cff to one as leading coefficients,
 *   and zero for all other coefficients. */

void dbl4_update_series
 ( int dim, int degp1,
   double **xhihi, double **xlohi, double **xhilo, double **xlolo,
   double **dxhihi, double **dxlohi, double **dxhilo, double **dxlolo,
   int vrblvl );
/*
 * DESCRIPTION :
 *   Adds the series in dx to x.
 *
 * ON ENTRY :
 *   dim       number of series in x and dx;
 *   degp1     degree plus one of the series;
 *   xhihi     highest doubles of the series to updated, not linearized;
 *   xlohi     2nd highest doubles of the series to updated, not linearized;
 *   xlohi     2nd lowest doubles of the series to updated, not linearized;
 *   xlolo     lowest doubles of the series to updated, not linearized;
 *   dxhihi    linearized highest doubles of the update of the series;
 *   dxlohi    linearized 2nd highest doubles of the update of the series;
 *   dxhilo    linearized 2nd lowest doubles of the update of the series;
 *   dxlolo    linearized lowest doubles of the update of the series;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   xhihi     highest doubles of the series x updated with dx;
 *   xlohi     second highest doubles of the series x updated with dx;
 *   xhilo     second lowest doubles of the series x updated with dx;
 *   xlolo     lowest doubles of the series x updated with dx. */

void dbl4_newton_lustep
 ( int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo,
   double *acchihi, double *acclohi, double *acchilo, double *acclolo,
   double **inputhihi, double **inputlohi,
   double **inputhilo, double **inputlolo,
   double ***outputhihi, double ***outputlohi,
   double ***outputhilo, double ***outputlolo,
   double **funvalhihi, double **funvallohi,
   double **funvalhilo, double **funvallolo,
   double ***jacvalhihi, double ***jacvallohi,
   double ***jacvalhilo, double ***jacvallolo,
   double **rhshihi, double **rhslohi, double **rhshilo, double **rhslolo,
   double **solhihi, double **sollohi, double **solhilo, double **sollolo,
   double **workmathihi, double **workmatlohi,
   double **workmathilo, double **workmatlolo,
   double *workvechihi, double *workveclohi,
   double *workvechilo, double *workveclolo,
   double **workrhshihi, double **workrhslohi,
   double **workrhshilo, double **workrhslolo,
   double **resvechihi, double **resveclohi,
   double **resvechilo, double **resveclolo,
   double *resmaxhihi, double *resmaxlohi,
   double *resmaxhilo, double *resmaxlolo, int *ipvt, int vrblvl );
/*
 * DESCRIPTION :
 *   Does one step with Newton's method to update a power series.
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
 *   cffhihi   highest doubles of the coefficients of the monomials;
 *   cfflohi   2nd highest doubles of the coefficients of the monomials;
 *   cffhilo   2nd lowest doubles of the coefficients of the monomials;
 *   cfflolo   lowest doubles of the coefficients of the monomials;
 *   acchihi   space to accumulate one power series of degree deg;
 *   acclohi   space to accumulate one power series of degree deg;
 *   acchilo   space to accumulate one power series of degree deg;
 *   acclolo   space to accumulate one power series of degree deg;
 *   inputhihi are the highest doubles of the coefficients of the power
 *             series of degree deg, for dim variables;
 *   inputlohi are the second highest doubles of the coefficients of the
 *             power series of degree deg, for dim variables;
 *   inputhilo are the second lowest doubles of the coefficients of the
 *             power series of degree deg, for dim variables;
 *   inputlolo are the lowest doubles of the coefficients of the power
 *             series of degree deg, for dim variables;
 *   outputhihi has space for the evaluated and differentiated monomials;
 *   outputlohi has space for the evaluated and differentiated monomials;
 *   outputhilo has space for the evaluated and differentiated monomials;
 *   outputlolo has space for the evaluated and differentiated monomials;
 *   funvalhihi has space for the evaluated power series;
 *   funvallohi has space for the evaluated power series;
 *   funvalhilo has space for the evaluated power series;
 *   funvallolo has space for the evaluated power series;
 *   jacvalhihi has space for deg+1 matrices of dimension dim;
 *   jacvallohi has space for deg+1 matrices of dimension dim;
 *   jacvalhilo has space for deg+1 matrices of dimension dim;
 *   jacvallolo has space for deg+1 matrices of dimension dim;
 *   rhshihi   space for deg+1 vectors of dimension dim;
 *   rhslohi   space for deg+1 vectors of dimension dim;
 *   rhshilo   space for deg+1 vectors of dimension dim;
 *   rhslolo   space for deg+1 vectors of dimension dim;
 *   solhihi   space for deg+1 vectors of dimension dim;
 *   sollohi   space for deg+1 vectors of dimension dim;
 *   solhilo   space for deg+1 vectors of dimension dim;
 *   sollolo   space for deg+1 vectors of dimension dim;
 *   wrkmathihi has work space allocated for a matrix of dimension dim;
 *   wrkmatlohi has work space allocated for a matrix of dimension dim;
 *   wrkmathilo has work space allocated for a matrix of dimension dim;
 *   wrkmatlolo has work space allocated for a matrix of dimension dim;
 *   wrkvechihi has work space allocated for a vector of dimension dim;
 *   wrkveclohi has work space allocated for a vector of dimension dim;
 *   wrkvechilo has work space allocated for a vector of dimension dim;
 *   wrkveclolo has work space allocated for a vector of dimension dim;
 *   resvechihi has space for deg+1 vectors of dimension dim;
 *   resveclohi has space for deg+1 vectors of dimension dim;
 *   resvechilo has space for deg+1 vectors of dimension dim;
 *   resveclolo has space for deg+1 vectors of dimension dim;
 *   ipvt      space allocated for dim pivots;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   funvalhihi is the highest doubles of the output[i][dim];
 *   funvallohi is the second highest doubles of the output[i][dim];
 *   funvalhilo is the second lowest doubles of the output[i][dim];
 *   funvallolo is the lowest doubles of the output[i][dim];
 *   jacvalhihi are the highest doubles of a matrix series,
 *             the leading coefficient is the Jacobian matrix.
 *   jacvallohi are the second highest doubles of a matrix series,
 *             the leading coefficient is the Jacobian matrix.
 *   jacvalhilo are the second lowest doubles of a matrix series,
 *             the leading coefficient is the Jacobian matrix.
 *   jacvallolo are the lowest doubles of a matrix series,
 *             the leading coefficient is the Jacobian matrix.
 *   rhshihi   highest doubles of the linearized right hand side
 *             are the function values subtracted by 1 and added by t;
 *   rhslohi   second highest doubles of the linearized right hand side
 *             are the function values subtracted by 1 and added by t;
 *   rhshilo   second lowest doubles of the linearized right hand side
 *             are the function values subtracted by 1 and added by t;
 *   rhslolo   lowest doubles ofthe linearized right hand side
 *             are the function values subtracted by 1 and added by t;
 *   wrkmathihi are the highest doubles of the LU factorization
 *             of the Jacobian matrix;
 *   wrkmatlohi are the second highest doubles of the LU factorization
 *             of the Jacobian matrix;
 *   wrkmathilo are the second lowest doubles of the LU factorization
 *             of the Jacobian matrix;
 *   wrkmatlolo are the lowest doubles of the LU factorization
 *             of the Jacobian matrix;
 *   resvechihi are the highest doubles of the residual vectors;
 *   resveclohi are the second highest doubles of the residual vectors;
 *   resvechilo are the second lowest doubles of the residual vectors;
 *   resveclolo are the lowest doubles of the residual vectors;
 *   resmaxhihi is the highest double of the maximum element
 *             of the residual vectors;
 *   resmaxlohi is the second highest double of the maximum element
 *             of the residual vectors;
 *   resmaxhilo is the second lowest double of the maximum element
 *             of the residual vectors;
 *   resmaxlolo is the lowest double of the maximum element
 *             of the residual vectors;
 *   ipvt      pivots used on the LU factorization of the lead matrix. */

void dbl4_newton_qrstep
 ( int szt, int nbt, int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo,
   double *acchihi, double *acclohi, double *acchilo, double *acclolo,
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
   double **workmathihi, double **workmatlohi,
   double **workmathilo, double **workmatlolo,
   double *workvechihi, double *workveclohi,
   double *workvechilo, double *workveclolo,
   double **resvechihi, double **resveclohi, 
   double **resvechilo, double **resveclolo, 
   double *resmaxhihi, double *resmaxlohi,
   double *resmaxhilo, double *resmaxlolo, int vrblvl, int mode );
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
 *   vrblvl    is the verbose level;
 *   mode      execution mode, 0 (GPU only), 1 (CPU only) or 2 (GPU+CPU).
 *
 * ON RETURN :
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
 *   resmaxlolo is lowest double of the maximum of residual vectors. */

#endif
