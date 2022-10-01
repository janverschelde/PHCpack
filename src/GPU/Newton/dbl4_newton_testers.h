// The file dbl4_newton_testers.h specifies test function for Newton's method
// on series in quad double precision.

#ifndef __dbl4_newton_testers_h__
#define __dbl4_newton_testers_h__

void dbl4_start_series_vector
 ( int dim, int deg,
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo );
/*
 * DESCRIPTION :
 *   Given space in cffhihi, cfflohi, cffhilo, cfflolo for the doubles
 *   for the vector of dim power series in quad double precision,
 *   of series truncated after degree deg,
 *   sets the coefficients of the start series. */

void cmplx4_start_series_vector
 ( int dim, int deg,
   double **cffrehihi, double **cffrelohi,
   double **cffrehilo, double **cffrelolo,
   double **cffimhihi, double **cffimlohi,
   double **cffimhilo, double **cffimlolo );
/*
 * DESCRIPTION :
 *   Given space allocated in the coefficient arrays,
 *   sets the coefficients of the start series. */

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

void cmplx4_unit_series_vector
 ( int dim, int deg,
   double **cffrehihi, double **cffrelohi,
   double **cffrehilo, double **cffrelolo,
   double **cffimhihi, double **cffimlohi,
   double **cffimhilo, double **cffimlolo );
/*
 * DESCRIPTION :
 *   Given space allocated in the coefficient arrays,
 *   returns the values for a complex unit series. */

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

void cmplx4_update_series
 ( int dim, int degp1,
   double **xrehihi, double **xrelohi, double **xrehilo, double **xrelolo,
   double **ximhihi, double **ximlohi, double **ximhilo, double **ximlolo,
   double **dxrehihi, double **dxrelohi, double **dxrehilo, double **dxrelolo,
   double **dximhihi, double **dximlohi, double **dximhilo, double **dximlolo,
   int vrblvl );
/*
 * DESCRIPTION :
 *   Adds the series in dx to x.
 *
 * ON ENTRY :
 *   dim       number of series in x and dx;
 *   degp1     degree plus one of the series;
 *   xrehihi   highest doubles of the real parts of series x, not linearized;
 *   xrelohi   second highest doubles of the real parts of series x;
 *   xrehilo   second lowest doubles of the real parts of series x;
 *   xrelolo   lowest doubles of the real parts of series x;
 *   ximhihi   highest doubles of the imaginary parts of series x;
 *   ximlohi   second highest doubles of the imaginary parts of series x;
 *   ximhilo   second lowest doubles of the imaginary parts of series x;
 *   ximlolo   lowest doubles of the imaginary parts of series x;
 *   dxrehihi  highest doubles of the real parts of the linearized update dx;
 *   dxrelohi  second highest doubles of the real parts of dx;
 *   dxrehilo  second lowest doubles of the real parts of dx;
 *   dxrelolo  lowest doubles of the real parts of dx;
 *   dximhihi  highest doubles of the imaginary parts of dx;
 *   dximlohi  second highest doubles of the imaginary parts of dx;
 *   dximhilo  second lowest doubles of the imaginary parts of the dx;
 *   dximlolo  lowest doubles of the imaginary parts of the dx;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   xrehihi   highest doubles of the real parts of the updated series x;
 *   xrelohi   second highest doubles of the real parts of x;
 *   xrehilo   second lowest doubles of the real parts of x;
 *   xrelolo   lowest doubles of the real parts of x;
 *   ximhihi   highest doubles fo the imaginary parts of x;
 *   ximlohi   second highest doubles fo the imaginary parts of x;
 *   ximhilo   second lowest doubles fo the imaginary parts of x;
 *   ximlolo   lowest doubles fo the imaginary parts of x. */

double dbl4_error3sum
 ( int dim1, int dim2, int dim3,
   double ***datahihi_h, double ***datalohi_h,
   double ***datahilo_h, double ***datalolo_h,
   double ***datahihi_d, double ***datalohi_d,
   double ***datahilo_d, double ***datalolo_d,
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
 *   datahihi_h are highest doubles of data computed on the host;
 *   datalohi_h are 2nd highest doubles of data computed on the host;
 *   datahilo_h are 2nd lowest doubles of data computed on the host;
 *   datalolo_h are lowest doubles of data computed on the host;
 *   datahihi_d are highest doubles of data computed on the device;
 *   datalohi_d are 2nd highest doubles of data computed on the device;
 *   datahilo_d are 2nd lowest doubles of data computed on the device;
 *   datalolo_d are lowest doubles of data computed on the device;
 *   banner    string for printing if verbose level vrblvl > 1;
 *   vrblvl    is the verbose level. */

double cmplx4_error3sum
 ( int dim1, int dim2, int dim3,
   double ***datarehihi_h, double ***datarelohi_h,
   double ***datarehilo_h, double ***datarelolo_h,
   double ***dataimhihi_h, double ***dataimlohi_h,
   double ***dataimhilo_h, double ***dataimlolo_h,
   double ***datarehihi_d, double ***datarelohi_d,
   double ***datarehilo_d, double ***datarelolo_d,
   double ***dataimhihi_d, double ***dataimlohi_d,
   double ***dataimhilo_d, double ***dataimlolo_d,
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
 *   datarehihi_h are highest doubles of real parts of data, on host;
 *   datarelohi_h are 2nd highest doubles of real parts of data, on host;
 *   datarehilo_h are 2nd lowest doubles of real parts of data, on host;
 *   datarelolo_h are lowest doubles of real parts of data, on host;
 *   dataimhihi_h are highest doubles of imaginary parts of data, on host;
 *   dataimlohi_h are 2nd highest doubles of imaginary parts of data, on host;
 *   dataimhilo_h are 2nd lowest doubles of imaginary parts of data, on host;
 *   dataimlolo_h are lowest doubles of imaginary parts of data, on host;
 *   datarehihi_d are highest doubles of real parts of data, on device;
 *   datarelohi_d are 2nd highest doubles of real parts of data, on device;
 *   datarehilo_d are 2nd lowest doubles of real parts of data, on device;
 *   datarelolo_d are lowest doubles of real parts of data, on device;
 *   dataimhihi_d are highest doubles of imaginary parts of data, on device;
 *   dataimlohi_d are 2nd highest doubles of imag parts of data, on device;
 *   dataimhilo_d are 2nd lowest doubles of imag parts of data, on device;
 *   dataimlolo_d are lowest doubles of imaginary parts of data, on device;
 *   banner    string for printing if verbose level vrblvl > 1;
 *   vrblvl    is the verbose level. */

double dbl4_error2sum
 ( int nrows, int ncols,
   double **datahihi_h, double **datalohi_h,
   double **datahilo_h, double **datalolo_h,
   double **datahihi_d, double **datalohi_d,
   double **datahilo_d, double **datalolo_d,
   std::string banner, int vrblvl );
/*
 * DESCRIPTION :
 *   Returns the sum of the absolute values of the differences
 *   between host and device real data, in two dimensions.
 *
 * ON ENTRY :
 *   nrows     number of rows of the matrices;
 *   ncols     number of columns of the matrices;
 *   datahihi_h are highest doubles of data computed on the host;
 *   datalohi_h are 2nd highest doubles of data computed on the host;
 *   datahilo_h are 2nd lowest doubles of data computed on the host;
 *   datalolo_h are lowest doubles of data computed on the host;
 *   datahihi_d are highest doubles of data computed on the device;
 *   datalohi_d are 2nd highest doubles of data computed on the device;
 *   datahilo_d are 2nd lowest doubles of data computed on the device;
 *   datalolo_d are lowest doubles of data computed on the device;
 *   banner    string for printing if verbose level vrblvl > 1;
 *   vrblvl    is the verbose level. */

double cmplx4_error2sum
 ( int nrows, int ncols,
   double **datarehihi_h, double **datarelohi_h,
   double **datarehilo_h, double **datarelolo_h,
   double **dataimhihi_h, double **dataimlohi_h,
   double **dataimhilo_h, double **dataimlolo_h,
   double **datarehihi_d, double **datarelohi_d,
   double **datarehilo_d, double **datarelolo_d,
   double **dataimhihi_d, double **dataimlohi_d,
   double **dataimhilo_d, double **dataimlolo_d,
   std::string banner, int vrblvl );
/*
 * DESCRIPTION :
 *   Returns the sum of the absolute values of the differences
 *   between host and device complex data, in two dimensions.
 *
 * ON ENTRY :
 *   nrows     number of rows of the matrices;
 *   ncols     number of columns of the matrices;
 *   datarehihi_h are highest doubles of real parts of data, on host;
 *   datarelohi_h are 2nd highest doubles of real parts of data, on host;
 *   datarehilo_h are 2nd lowest doubles of real parts of data, on host;
 *   datarelolo_h are lowest doubles of real parts of data, on host;
 *   dataimhihi_h are highest doubles of imaginary parts of data, on host;
 *   dataimlohi_h are 2nd highest doubles of imaginary parts of data, on host;
 *   dataimhilo_h are 2nd lowest doubles of imaginary parts of data, on host;
 *   dataimlolo_h are lowest doubles of imaginary parts of data, on host;
 *   datarehihi_d are highest doubles of real parts of data, on device;
 *   datarelohi_d are 2nd highest doubles of real parts of data, on device;
 *   datarehilo_d are 2nd lowest doubles of real parts of data, on device;
 *   datarelolo_d are lowest doubles of real parts of data, on device;
 *   dataimhihi_d are highest doubles of imaginary parts of data, on device;
 *   dataimlohi_d are 2nd highest doubles of imag parts of data, on device;
 *   dataimhilo_d are 2nd lowest doubles of imag parts of data, on device;
 *   dataimlolo_d are lowest doubles of imaginary parts of data, on device;
 *   banner    string for printing if verbose level vrblvl > 1;
 *   vrblvl    is the verbose level. */

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

void cmplx4_newton_qrstep
 ( int szt, int nbt, int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffrehihi, double **cffrelohi,
   double **cffrehilo, double **cffrelolo,
   double **cffimhihi, double **cffimlohi,
   double **cffimhilo, double **cffimlolo,
   double *accrehihi, double *accrelohi,
   double *accrehilo, double *accrelolo,
   double *accimhihi, double *accimlohi,
   double *accimhilo, double *accimlolo,
   double **inputrehihi_h, double **inputrelohi_h,
   double **inputrehilo_h, double **inputrelolo_h,
   double **inputimhihi_h, double **inputimlohi_h,
   double **inputimhilo_h, double **inputimlolo_h,
   double **inputrehihi_d, double **inputrelohi_d,
   double **inputrehilo_d, double **inputrelolo_d,
   double **inputimhihi_d, double **inputimlohi_d,
   double **inputimhilo_d, double **inputimlolo_d,
   double ***outputrehihi_h, double ***outputrelohi_h,
   double ***outputrehilo_h, double ***outputrelolo_h,
   double ***outputimhihi_h, double ***outputimlohi_h,
   double ***outputimhilo_h, double ***outputimlolo_h,
   double ***outputrehihi_d, double ***outputrelohi_d,
   double ***outputrehilo_d, double ***outputrelolo_d,
   double ***outputimhihi_d, double ***outputimlohi_d,
   double ***outputimhilo_d, double ***outputimlolo_d,
   double **funvalrehihi_h, double **funvalrelohi_h,
   double **funvalrehilo_h, double **funvalrelolo_h,
   double **funvalimhihi_h, double **funvalimlohi_h,
   double **funvalimhilo_h, double **funvalimlolo_h,
   double **funvalrehihi_d, double **funvalrelohi_d,
   double **funvalrehilo_d, double **funvalrelolo_d,
   double **funvalimhihi_d, double **funvalimlohi_d,
   double **funvalimhilo_d, double **funvalimlolo_d,
   double ***jacvalrehihi_h, double ***jacvalrelohi_h,
   double ***jacvalrehilo_h, double ***jacvalrelolo_h,
   double ***jacvalimhihi_h, double ***jacvalimlohi_h,
   double ***jacvalimhilo_h, double ***jacvalimlolo_h,
   double ***jacvalrehihi_d, double ***jacvalrelohi_d,
   double ***jacvalrehilo_d, double ***jacvalrelolo_d,
   double ***jacvalimhihi_d, double ***jacvalimlohi_d,
   double ***jacvalimhilo_d, double ***jacvalimlolo_d,
   double **rhsrehihi_h, double **rhsrelohi_h,
   double **rhsrehilo_h, double **rhsrelolo_h,
   double **rhsimhihi_h, double **rhsimlohi_h,
   double **rhsimhilo_h, double **rhsimlolo_h,
   double **rhsrehihi_d, double **rhsrelohi_d, 
   double **rhsrehilo_d, double **rhsrelolo_d, 
   double **rhsimhihi_d, double **rhsimlohi_d,
   double **rhsimhilo_d, double **rhsimlolo_d,
   double **urhsrehihi_h, double **urhsrelohi_h,
   double **urhsrehilo_h, double **urhsrelolo_h,
   double **urhsimhihi_h, double **urhsimlohi_h,
   double **urhsimhilo_h, double **urhsimlolo_h,
   double **urhsrehihi_d, double **urhsrelohi_d,
   double **urhsrehilo_d, double **urhsrelolo_d,
   double **urhsimhihi_d, double **urhsimlohi_d,
   double **urhsimhilo_d, double **urhsimlolo_d,
   double **solrehihi_h, double **solrelohi_h,
   double **solrehilo_h, double **solrelolo_h,
   double **solimhihi_h, double **solimlohi_h, 
   double **solimhilo_h, double **solimlolo_h, 
   double **solrehihi_d, double **solrelohi_d, 
   double **solrehilo_d, double **solrelolo_d, 
   double **solimhihi_d, double **solimlohi_d,
   double **solimhilo_d, double **solimlolo_d,
   double **Qrehihi_h, double **Qrelohi_h,
   double **Qrehilo_h, double **Qrelolo_h,
   double **Qimhihi_h, double **Qimlohi_h,
   double **Qimhilo_h, double **Qimlolo_h,
   double **Qrehihi_d, double **Qrelohi_d,
   double **Qrehilo_d, double **Qrelolo_d,
   double **Qimhihi_d, double **Qimlohi_d, 
   double **Qimhilo_d, double **Qimlolo_d, 
   double **Rrehihi_h, double **Rrelohi_h,
   double **Rrehilo_h, double **Rrelolo_h,
   double **Rimhihi_h, double **Rimlohi_h, 
   double **Rimhilo_h, double **Rimlolo_h, 
   double **Rrehihi_d, double **Rrelohi_d,
   double **Rrehilo_d, double **Rrelolo_d,
   double **Rimhihi_d, double **Rimlohi_d,
   double **Rimhilo_d, double **Rimlolo_d,
   double **workmatrehihi, double **workmatrelohi,
   double **workmatrehilo, double **workmatrelolo,
   double **workmatimhihi, double **workmatimlohi,
   double **workmatimhilo, double **workmatimlolo,
   double *workvecrehihi, double *workvecrelohi,
   double *workvecrehilo, double *workvecrelolo,
   double *workvecimhihi, double *workvecimlohi,
   double *workvecimhilo, double *workvecimlolo,
   double **resvecrehihi, double **resvecrelohi,
   double **resvecrehilo, double **resvecrelolo,
   double **resvecimhihi, double **resvecimlohi,
   double **resvecimhilo, double **resvecimlolo,
   double *resmaxhihi, double *resmaxlohi,
   double *resmaxhilo, double *resmaxlolo, int vrblvl, int mode );
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
 *   nvr       nvr[i] is the number of variables in the i-th monomial;
 *   idx       idx[i] are the indices of the variables in monomial i;
 *   exp       exp[i] are the exponents of the variables in monomial i;
 *   nbrfac    nbrfac[i] are the number of exponents > 1 in monomial i;
 *   expfac    expfac[i] are the exponents in the i-th polynomial
 *             that are larger than one, minus one in the factor,
 *             if exp[i][k] > 1, then expfac[i][k] = exp[i][k] - 1;
 *   cffrehihi are the highest doubles of the real parts of
 *             the coefficients of the monomials;
 *   cffrelohi are the second highest doubles of the real parts of
 *             the coefficients of the monomials;
 *   cffrehilo are the second lowest doubles of the real parts of
 *             the coefficients of the monomials;
 *   cffrelolo are the lowest doubles of the real parts of
 *             the coefficients of the monomials;
 *   cffimhihi are the highest doubles of the imaginary parts of
 *             the coefficients of the monomials;
 *   cffimlohi are the second highest doubles of the imaginary parts of
 *             the coefficients of the monomials;
 *   cffimhilo are the second lowest doubles of the imaginary parts of
 *             the coefficients of the monomials;
 *   cffimlolo are the lowest doubles of the imaginary parts of
 *             the coefficients of the monomials;
 *   accrehihi has space to accumulate one series of degree deg;
 *   accrelohi has space to accumulate one series of degree deg;
 *   accrehilo has space to accumulate one series of degree deg;
 *   accrelolo has space to accumulate one series of degree deg;
 *   accimhihi has space to accumulate one series of degree deg;
 *   accimlohi has space to accumulate one series of degree deg;
 *   accimhilo has space to accumulate one series of degree deg;
 *   accimlolo has space to accumulate one series of degree deg;
 *   inputrehihi_h are the highest doubles of the real parts of series
 *             of degree deg, for dim variables, computed on host;
 *   inputrelohi_h are the second highest doubles of the real parts of series
 *             of degree deg, for dim variables, computed on host;
 *   inputrehilo_h are the second lowest doubles of the real parts of series
 *             of degree deg, for dim variables, computed on host;
 *   inputrelolo_h are the lowest doubles of the real parts of series
 *             of degree deg, for dim variables, computed on host;
 *   inputimhihi_h are the highest doubles of the imaginary parts of
 *             the series of degree deg for dim variables, computed on host;
 *   inputimlohi_h are the second highest doubles of the imaginary parts of
 *             the series of degree deg for dim variables, computed on host;
 *   inputimhilo_h are the second lowest doubles of the imaginary parts of
 *             the series of degree deg for dim variables, computed on host;
 *   inputimlolo_h are the lowest doubles of the imaginary parts of
 *             the series of degree deg for dim variables, computed on host;
 *   inputrehihi_d has space for series computed on device;
 *   inputrelohi_d has space for series computed on device;
 *   inputrehilo_d has space for series computed on device;
 *   inputrelolo_d has space for series computed on device;
 *   inputimhihi_d has space for series computed on device;
 *   inputimlohi_d has space for series computed on device;
 *   inputimhilo_d has space for series computed on device;
 *   inputimlolo_d has space for series computed on device;
 *   outputrehihi_h has space for the evaluated and differentiated monomials,
 *             computed on the host;
 *   outputrelohi_h has space for the evaluated and differentiated monomials,
 *             computed on the host;
 *   outputrehilo_h has space for the evaluated and differentiated monomials,
 *             computed on the host;
 *   outputrelolo_h has space for the evaluated and differentiated monomials,
 *             computed on the host;
 *   outputimhihi_h has space for the evaluated and differentiated monomials,
 *             computed on the host;
 *   outputimlohi_h has space for the evaluated and differentiated monomials,
 *             computed on the host;
 *   outputimhilo_h has space for the evaluated and differentiated monomials,
 *             computed on the host;
 *   outputimlolo_h has space for the evaluated and differentiated monomials,
 *             computed on the host;
 *   outputrehihi_d has space for the evaluated and differentiated monomials,
 *             computed on the device;
 *   outputrelohi_d has space for the evaluated and differentiated monomials,
 *             computed on the device;
 *   outputrehilo_d has space for the evaluated and differentiated monomials,
 *             computed on the device;
 *   outputrelolo_d has space for the evaluated and differentiated monomials,
 *             computed on the device;
 *   outputimhihi_d has space for the evaluated and differentiated monomials,
 *             computed on the device;
 *   outputimlohi_d has space for the evaluated and differentiated monomials,
 *             computed on the device;
 *   outputimhilo_d has space for the evaluated and differentiated monomials,
 *             computed on the device;
 *   outputimlolo_d has space for the evaluated and differentiated monomials,
 *             computed on the device;
 *   funvalrehihi_h has space for the evaluated series computed by host;
 *   funvalrelohi_h has space for the evaluated series computed by host;
 *   funvalrehilo_h has space for the evaluated series computed by host;
 *   funvalrelolo_h has space for the evaluated series computed by host;
 *   funvalimhihi_h has space for the evaluated series computed by host;
 *   funvalimlohi_h has space for the evaluated series computed by host;
 *   funvalimhilo_h has space for the evaluated series computed by host;
 *   funvalimlolo_h has space for the evaluated series computed by host;
 *   funvalrehihi_d has space for the evaluated series computed by device;
 *   funvalrelohi_d has space for the evaluated series computed by device;
 *   funvalrehilo_d has space for the evaluated series computed by device;
 *   funvalrelolo_d has space for the evaluated series computed by device;
 *   funvalimhihi_d has space for the evaluated series computed by device;
 *   funvalimlohi_d has space for the evaluated series computed by device;
 *   funvalimhilo_d has space for the evaluated series computed by device;
 *   funvalimlolo_d has space for the evaluated series computed by device;
 *   jacvalrehihi_h has space for deg+1 matrices of dimension dim on host;
 *   jacvalrelohi_h has space for deg+1 matrices of dimension dim on host;
 *   jacvalrehilo_h has space for deg+1 matrices of dimension dim on host;
 *   jacvalrelolo_h has space for deg+1 matrices of dimension dim on host;
 *   jacvalimhihi_h has space for deg+1 matrices of dimension dim on host;
 *   jacvalimlohi_h has space for deg+1 matrices of dimension dim on host;
 *   jacvalimhilo_h has space for deg+1 matrices of dimension dim on host;
 *   jacvalimlolo_h has space for deg+1 matrices of dimension dim on host;
 *   jacvalrehihi_d has space for deg+1 matrices of dimension dim on device;
 *   jacvalrelohi_d has space for deg+1 matrices of dimension dim on device;
 *   jacvalrehilo_d has space for deg+1 matrices of dimension dim on device;
 *   jacvalrelolo_d has space for deg+1 matrices of dimension dim on device;
 *   jacvalimhihi_d has space for deg+1 matrices of dimension dim on device;
 *   jacvalimlohi_d has space for deg+1 matrices of dimension dim on device;
 *   jacvalimhilo_d has space for deg+1 matrices of dimension dim on device;
 *   jacvalimhilo_d has space for deg+1 matrices of dimension dim on device;
 *   rhsrehihi_h has space for deg+1 vectors of dimension dim on host;
 *   rhsrelohi_h has space for deg+1 vectors of dimension dim on host;
 *   rhsrehilo_h has space for deg+1 vectors of dimension dim on host;
 *   rhsrelolo_h has space for deg+1 vectors of dimension dim on host;
 *   rhsimhihi_h has space for deg+1 vectors of dimension dim on host;
 *   rhsimlohi_h has space for deg+1 vectors of dimension dim on host;
 *   rhsimhilo_h has space for deg+1 vectors of dimension dim on host;
 *   rhsimlolo_h has space for deg+1 vectors of dimension dim on host;
 *   rhsrehihi_d has space for deg+1 vectors of dimension dim on device;
 *   rhsrelohi_d has space for deg+1 vectors of dimension dim on device;
 *   rhsrehilo_d has space for deg+1 vectors of dimension dim on device;
 *   rhsrelolo_d has space for deg+1 vectors of dimension dim on device;
 *   rhsimhihi_d has space for deg+1 vectors of dimension dim on device;
 *   rhsimlohi_d has space for deg+1 vectors of dimension dim on device;
 *   rhsimhilo_d has space for deg+1 vectors of dimension dim on device;
 *   rhsimlolo_d has space for deg+1 vectors of dimension dim on device;
 *   urhsrehihi_h has space for updated right hand side computed by host;
 *   urhsrelohi_h has space for updated right hand side computed by host;
 *   urhsrehilo_h has space for updated right hand side computed by host;
 *   urhsrelolo_h has space for updated right hand side computed by host;
 *   urhsimhihi_h has space for updated right hand side computed by host;
 *   urhsimlohi_h has space for updated right hand side computed by host;
 *   urhsimhilo_h has space for updated right hand side computed by host;
 *   urhsimlolo_h has space for updated right hand side computed by host;
 *   urhsrehihi_d has space for updated right hand side computed by device; 
 *   urhsrelohi_d has space for updated right hand side computed by device; 
 *   urhsrehilo_d has space for updated right hand side computed by device; 
 *   urhsrelolo_d has space for updated right hand side computed by device; 
 *   urhsimhihi_d has space for updated right hand side computed by device; 
 *   urhsimlohi_d has space for updated right hand side computed by device; 
 *   urhsimhilo_d has space for updated right hand side computed by device; 
 *   urhsimlolo_d has space for updated right hand side computed by device; 
 *   solrehihi_h has space for deg+1 vectors of dimension dim;
 *   solrelohi_h has space for deg+1 vectors of dimension dim;
 *   solrehilo_h has space for deg+1 vectors of dimension dim;
 *   solrelolo_h has space for deg+1 vectors of dimension dim;
 *   solimhihi_h has space for deg+1 vectors of dimension dim;
 *   solimlohi_h has space for deg+1 vectors of dimension dim;
 *   solimhilo_h has space for deg+1 vectors of dimension dim;
 *   solimlolo_h has space for deg+1 vectors of dimension dim;
 *   solrehihi_d has space for deg+1 vectors of dimension dim;
 *   solrelohi_d has space for deg+1 vectors of dimension dim;
 *   solrehilo_d has space for deg+1 vectors of dimension dim;
 *   solrelolo_d has space for deg+1 vectors of dimension dim;
 *   solimhihi_d has space for deg+1 vectors of dimension dim;
 *   solimlohi_d has space for deg+1 vectors of dimension dim;
 *   solimhilo_d has space for deg+1 vectors of dimension dim;
 *   solimlolo_d has space for deg+1 vectors of dimension dim;
 *   Qrehihi_h has space for the Q computed by the host;
 *   Qrelohi_h has space for the Q computed by the host;
 *   Qrehilo_h has space for the Q computed by the host;
 *   Qrelolo_h has space for the Q computed by the host;
 *   Qimhihi_h has space for the Q computed by the host;
 *   Qimlohi_h has space for the Q computed by the host;
 *   Qimhilo_h has space for the Q computed by the host;
 *   Qimlolo_h has space for the Q computed by the host;
 *   Qrehihi_d has space for the Q computed by the device;
 *   Qrelohi_d has space for the Q computed by the device;
 *   Qrehilo_d has space for the Q computed by the device;
 *   Qrelolo_d has space for the Q computed by the device;
 *   Qimhihi_d has space for the Q computed by the device;
 *   Qimlohi_d has space for the Q computed by the device;
 *   Qimhilo_d has space for the Q computed by the device;
 *   Qimlolo_d has space for the Q computed by the device;
 *   Rrehihi_h has space for the R computed by the host;
 *   Rrelohi_h has space for the R computed by the host;
 *   Rrehilo_h has space for the R computed by the host;
 *   Rrelolo_h has space for the R computed by the host;
 *   Rimhihi_h has space for the R computed by the host;
 *   Rimlohi_h has space for the R computed by the host;
 *   Rimhilo_h has space for the R computed by the host;
 *   Rimlolo_h has space for the R computed by the host;
 *   Rrehihi_d has space for the R computed by the device;
 *   Rrelohi_d has space for the R computed by the device;
 *   Rrehilo_d has space for the R computed by the device;
 *   Rrelolo_d has space for the R computed by the device;
 *   Rimhihi_d has space for the R computed by the device;
 *   Rimlohi_d has space for the R computed by the device;
 *   Rimhilo_d has space for the R computed by the device;
 *   Rimlolo_d has space for the R computed by the device;
 *   wrkmatrehihi is work space allocated for a matrix of dimension dim;
 *   wrkmatrelohi is work space allocated for a matrix of dimension dim;
 *   wrkmatrehilo is work space allocated for a matrix of dimension dim;
 *   wrkmatrelolo is work space allocated for a matrix of dimension dim;
 *   wrkmatimhihi is work space allocated for a matrix of dimension dim;
 *   wrkmatimlohi is work space allocated for a matrix of dimension dim;
 *   wrkmatimhilo is work space allocated for a matrix of dimension dim;
 *   wrkmatimlolo is work space allocated for a matrix of dimension dim;
 *   wrkvecrehihi is work space allocated for a vector of dimension dim;
 *   wrkvecrelohi is work space allocated for a vector of dimension dim;
 *   wrkvecrehilo is work space allocated for a vector of dimension dim;
 *   wrkvecrelolo is work space allocated for a vector of dimension dim;
 *   wrkvecimhihi is work space allocated for a vector of dimension dim;
 *   wrkvecimlohi is work space allocated for a vector of dimension dim;
 *   wrkvecimhilo is work space allocated for a vector of dimension dim;
 *   wrkvecimlolo is work space allocated for a vector of dimension dim;
 *   resvecrehihi has space for deg+1 vectors of dimension dim;
 *   resvecrelohi has space for deg+1 vectors of dimension dim;
 *   resvecrehilo has space for deg+1 vectors of dimension dim;
 *   resvecrelolo has space for deg+1 vectors of dimension dim;
 *   resvecimhihi has space for deg+1 vectors of dimension dim;
 *   resvecimlohi has space for deg+1 vectors of dimension dim;
 *   resvecimhilo has space for deg+1 vectors of dimension dim;
 *   resvecimlolo has space for deg+1 vectors of dimension dim;
 *   vrblvl    is the verbose level;
 *   mode      execution mode, 0 (GPU only), 1 (CPU only) or 2 (GPU+CPU).
 *
 * ON RETURN :
 *   inputrehihi_h has the highest doubles of the real parts of series,
 *             computed on host;
 *   inputrelohi_h has the second highest doubles of the real parts of series,
 *             computed on host;
 *   inputrehilo_h has the second lowest doubles of the real parts of series,
 *             computed on host;
 *   inputrelolo_h has the lowest doubles of the real parts of series,
 *             computed on host;
 *   inputimhihi_h has the highest doubles of the imag parts of series,
 *             computed on host;
 *   inputimlohi_h has the second highest doubles of the imag parts of series,
 *             computed on host;
 *   inputimhilo_h has the second lowest doubles of the imag parts of series,
 *             computed on host;
 *   inputimlolo_h has the lowest doubles of the imag parts of series,
 *             computed on host;
 *   inputrehihi_d has the highest doubles of the real parts of series,
 *             computed on device;
 *   inputrelohi_d has the second highest doubles of the real parts of series,
 *             computed on device;
 *   inputrehilo_d has the second lowest doubles of the real parts of series,
 *             computed on device;
 *   inputrelolo_d has the lowest doubles of the real parts of series,
 *             computed on device;
 *   inputimhihi_d has the highest doubles of the imag parts of series,
 *             computed on device;
 *   inputimlohi_d has the second highest doubles of the imag parts of series,
 *             computed on device;
 *   inputimhilo_d has the second lowest doubles of the imag parts of series,
 *             computed on device;
 *   inputimlolo_d has the lowest doubles of the imag parts of series,
 *             computed on device;
 *   outputrehihi_h has the highest doubles of the real part of output,
 *             computed on host;
 *   outputrelohi_h has the second highest doubles of the real part of output,
 *             computed on host;
 *   outputrehilo_h has the second lowest doubles of the real part of output,
 *             computed on host;
 *   outputrelolo_h has the lowest doubles of the real part of output,
 *             computed on host;
 *   outputimhihi_h has the highest doubles of the imag part of output,
 *             computed on host;
 *   outputimlohi_h has the second highest doubles of the imag part of output,
 *             computed on host;
 *   outputimhilo_h has the second lowest doubles of the imag part of output,
 *             computed on host;
 *   outputimlolo_h has the lowest doubles of the imag part of output,
 *             computed on host;
 *   outputrehihi_d has the highest doubles of the real part of output,
 *             computed on device;
 *   outputrelohi_d has the second highest doubles of the real part of output,
 *             computed on device;
 *   outputrehilo_d has the second lowest doubles of the real part of output,
 *             computed on device;
 *   outputrelolo_d has the lowest doubles of the real part of output,
 *             computed on device;
 *   outputimhihi_d has the highest doubles of the imag part of output,
 *             computed on device;
 *   outputimlohi_d has the second highest doubles of the imag part of output,
 *             computed on device;
 *   outputimhilo_d has the second lowest doubles of the imag part of output,
 *             computed on device;
 *   outputimlolo_d has the lowest doubles of the imag part of output,
 *             computed on device;
 *   funvalrehihi_h is outputrehihi[i][dim], on host;
 *   funvalrelohi_h is outputrelohi[i][dim], on host;
 *   funvalrehilo_h is outputrehilo[i][dim], on host;
 *   funvalrelolo_h is outputrelolo[i][dim], on host;
 *   funvalimhihi_h is outputimhihi[i][dim], on host;
 *   funvalimlohi_h is outputimlohi[i][dim], on host;
 *   funvalimhilo_h is outputimhilo[i][dim], on host;
 *   funvalimlolo_h is outputimlolo[i][dim], on host;
 *   funvalrehihi_d is outputrehihi[i][dim], on device;
 *   funvalrelohi_d is outputrelohi[i][dim], on device;
 *   funvalrehilo_d is outputrehilo[i][dim], on device;
 *   funvalrelolo_d is outputrelolo[i][dim], on device;
 *   funvalimhihi_d is outputimhihi[i][dim], on device;
 *   funvalimlohi_d is outputimlohi[i][dim], on device;
 *   funvalimhilo_d is outputimhilo[i][dim], on device;
 *   funvalimlolo_d is outputimlolo[i][dim], on device;
 *   jacvalrehihi_h are the highest doubles of the real parts
 *             of a matrix series, computed on host;
 *   jacvalrelohi_h are the second highest doubles of the real parts
 *             of a matrix series, computed on host;
 *   jacvalrehilo_h are the second lowest doubles of the real parts
 *             of a matrix series, computed on host;
 *   jacvalrelolo_h are the lowest doubles of the real parts
 *             of a matrix series, computed on host;
 *   jacvalimhihi_h are the highest doubles of the real parts
 *             of a matrix series, computed on host;
 *   jacvalimlohi_h are the second highest doubles of the real parts
 *             of a matrix series, computed on host;
 *   jacvalimhilo_h are the second lowest doubles of the real parts
 *             of a matrix series, computed on host;
 *   jacvalimlolo_h are the lowest doubles of the real parts
 *             of a matrix series, computed on host;
 *   jacvalrehihi_d are the highest doubles of the real parts
 *             of a matrix series, computed on device;
 *   jacvalrelohi_d are the second highest doubles of the real parts
 *             of a matrix series, computed on device;
 *   jacvalrehilo_d are the second lowest doubles of the real parts
 *             of a matrix series, computed on device;
 *   jacvalrelolo_d are the lowest doubles of the real parts
 *             of a matrix series, computed on device;
 *   jacvalimhihi_d are the highest doubles of the real parts
 *             of a matrix series, computed on device;
 *   jacvalimlohi_d are the second highest doubles of the real parts
 *             of a matrix series, computed on device;
 *   jacvalimhilo_d are the second lowest doubles of the real parts
 *             of a matrix series, computed on device;
 *   jacvalimlolo_d are the lowest doubles of the real parts
 *             of a matrix series, computed on device;
 *   rhsrehihi_h are highest doubles of real parts of the linearized
 *             right hand side, subtracted by 1 and added by t, on host;
 *   rhsrelohi_h are second highest doubles of real parts of the linearized
 *             right hand side, subtracted by 1 and added by t, on host;
 *   rhsrehilo_h are second lowest doubles of real parts of the linearized
 *             right hand side, subtracted by 1 and added by t, on host;
 *   rhsrelolo_h are lowest doubles of real parts of the linearized
 *             right hand side, subtracted by 1 and added by t, on host;
 *   rhsimhihi_h are highest doubles of imag parts of the linearized
 *             right hand side, subtracted by 1 and added by t, on host;
 *   rhsimlohi_h are second highest doubles of imag parts of the linearized
 *             right hand side, subtracted by 1 and added by t, on host;
 *   rhsimhilo_h are second lowest doubles of imag parts of the linearized
 *             right hand side, subtracted by 1 and added by t, on host;
 *   rhsimlolo_h are lowest doubles of imag parts of the linearized
 *             right hand side, subtracted by 1 and added by t, on host;
 *   rhsrehihi_d are highest doubles of real parts of the linearized 
 *             right hand side, subtracted by 1 and added by t, on device;
 *   rhsrelohi_d are second highest doubles of real parts of the linearized 
 *             right hand side, subtracted by 1 and added by t, on device;
 *   rhsrehilo_d are second lowest doubles of real parts of the linearized,
 *             right hand side, subtracted by 1 and added by t, on device;
 *   rhsrelolo_d are lowest doubles of real parts of the linearized,
 *             right hand side, subtracted by 1 and added by t, on device;
 *   rhsimhihi_d are highest doubles of imag parts of the linearized
 *             right hand side, subtracted by 1 and added by t, on device;
 *   rhsimlohi_d are second highest doubles of imag parts of the linearized
 *             right hand side, subtracted by 1 and added by t, on device;
 *   rhsimhilo_d are second lowest doubles of imag parts of the linearized
 *             right hand side, subtracted by 1 and added by t, on device;
 *   rhsimlolo_d are lowest doubles of imag parts of the linearized
 *             right hand side, subtracted by 1 and added by t, on device;
 *   urhsrehihi_h are the highest doubles of the real parts
 *             of the right hand side updated by the host;
 *   urhsrelohi_h are the second highest doubles of the real parts
 *             of the right hand side updated by the host;
 *   urhsrehilo_h are the second lowest doubles of the real parts
 *             of the right hand side updated by the host;
 *   urhsrelolo_h are the lowest doubles of the real parts
 *             of the right hand side updated by the host;
 *   urhsimhihi_h are the highest doubles of the imaginary parts
 *             of the right hand side updated by the host;
 *   urhsimlohi_h are the second highest doubles of the imaginary parts
 *             of the right hand side updated by the host;
 *   urhsimhilo_h are the second lowest doubles of the imaginary parts
 *             of the right hand side updated by the host;
 *   urhsimlolo_h are the lowest doubles of the imaginary parts
 *             of the right hand side updated by the host;
 *   urhsrehihi_d are the highest doubles of the real parts
 *             of the right hand side updated by the device;
 *   urhsrelohi_d are the second highest doubles of the real parts
 *             of the right hand side updated by the device;
 *   urhsrehilo_d are the second lowest doubles of the real parts
 *             of the right hand side updated by the device;
 *   urhsrelolo_d are the lowest doubles of the real parts
 *             of the right hand side updated by the device;
 *   urhsimhihi_d are the highest doubles of the imaginary parts
 *             of the right hand side updated by the device;
 *   urhsimlohi_d are the second highest doubles of the imaginary parts
 *             of the right hand side updated by the device;
 *   urhsimhilo_d are the second lowest doubles of the imaginary parts
 *             of the right hand side updated by the device;
 *   urhsimlolo_d are the lowest doubles of the imaginary parts
 *             of the right hand side updated by the device;
 *   solrehihi_h are the highest doubles of real parts of solution on host;
 *   solrelohi_h are the second highest doubles of real parts of sol on host;
 *   solrehilo_h are the second lowest doubles of real parts of sol on host;
 *   solrelolo_h are the lowest doubles of real parts of solution on host;
 *   solimhihi_h are the highest doubles of imag parts of solution on host;
 *   solimlohi_h are the second highest doubles of imag parts of sol on host;
 *   solimhilo_h are the second lowest doubles of imag parts of sol on host;
 *   solimlolo_h are the lowest doubles of imag parts of solution on host;
 *   solrehihi_d are the highest doubles of real parts of solution on device;
 *   solrelohi_d are the second highest doubles of real parts of sol on device;
 *   solrehilo_d are the second lowest doubles of real parts of sol on device;
 *   solrelolo_d are the lowest doubles of real parts of solution on device;
 *   solimhihi_d are the highest doubles of imag parts of solution on device;
 *   solimlohi_d are the seoncd highest doubles of imag parts of sol on device;
 *   solimhilo_d are the seoncd lowest doubles of imag parts of sol on device;
 *   solimlolo_d are the lowest doubles of imag parts of solution on device;
 *   Qrehihi_h are the highest doubles of real Q of the QR on host;
 *   Qrelohi_h are the second highest doubles of real Q of the QR on host;
 *   Qrehilo_h are the second lowest doubles of real Q of the QR on host;
 *   Qrelolo_h are the lowest doubles of real Q of the QR on host;
 *   Qimhihi_h are the highest doubles of imaginary Q of the QR on host;
 *   Qimlohi_h are the second highest doubles of imag Q of the QR on host;
 *   Qimhilo_h are the second lowest doubles of imag Q of the QR on host;
 *   Qimlolo_h are the lowest doubles of imaginary Q of the QR on host;
 *   Qrehihi_d are the highest doubles of real Q of the QR on device;
 *   Qrelohi_d are the second highest doubles of real Q of the QR on device;
 *   Qrehilo_d are the second lowest doubles of real Q of the QR on device;
 *   Qrelolo_d are the lowest doubles of real Q of the QR on device;
 *   Qimhihi_d are the highest doubles of imaginary Q of the QR on device;
 *   Qimlohi_d are the second highest doubles of imag Q of the QR on device;
 *   Qimhilo_d are the second lowest doubles of imag Q of the QR on device;
 *   Qimlolo_d are the lowest doubles of imaginary Q of the QR on device;
 *   Rrehihi_h are the highest doubles of real R of the QR on host;
 *   Rrelohi_h are the second highest doubles of real R of the QR on host;
 *   Rrehilo_h are the second lowest doubles of real R of the QR on host;
 *   Rrelolo_h are the lowest doubles of real R of the QR on host;
 *   Rimhihi_h are the highest doubles of imaginary R of the QR on host;
 *   Rimlohi_h are the second highest doubles of imag R of the QR on host;
 *   Rimhilo_h are the second lowest doubles of imag R of the QR on host;
 *   Rimlolo_h are the lowest doubles of imaginary R of the QR on host;
 *   Rrehihi_d are the highest doubles of real R of the QR on device;
 *   Rrelohi_d are the second highest doubles of real R of the QR on device;
 *   Rrehilo_d are the second lowest doubles of real R of the QR on device;
 *   Rrelolo_d are the lowest doubles of real R of the QR on device;
 *   Rimhihi_d are the highest doubles of imaginary R of the QR on device;
 *   Rimlohi_d are the second highest doubles of imag R of the QR on device;
 *   Rimhilo_d are the second lowest doubles of imag R of the QR on device;
 *   Rimlolo_d are the lowest doubles of imaginary R of the QR on device;
 *   wrkmat    has a copy of the Jacobian matrix;
 *   resvecrehi are highest doubles of the real parts of the residual vectors;
 *   resvecrelo are second highest doubles of the real parts of resvec;
 *   resvecrelo are second lowest doubles of the real parts of resvec;
 *   resvecrelo are lowest doubles of the real parts of resvec;
 *   resvecimhi are highest doubles of the imag parts of resvec;
 *   resvecimhi are second highest doubles of the imag parts of resvec;
 *   resvecimlo are second lowest doubles of the imag parts of resvec;
 *   resvecimlo are lowest doubles of the imag parts of resvec;
 *   resmaxhihi is the highest double of the max norm of the residual;
 *   resmaxlohi is the second highest double of the max norm of the residual;
 *   resmaxhilo is the second lowest double of the max norm of the residual;
 *   resmaxlolo is the lowest double of the max norm of the residual. */

int test_dbl4_real_newton
 ( int szt, int nbt, int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   int nbsteps, int mode, int vrblvl );
/*
 * DESCRIPTION :
 *   Runs Newton on a monomial system with real quad double arithmetic.
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

int test_dbl4_complex_newton
 ( int szt, int nbt, int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   int nbsteps, int mode, int vrblvl );
/*
 * DESCRIPTION :
 *   Runs Newton on a monomial system with complex quad double arithmetic.
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
