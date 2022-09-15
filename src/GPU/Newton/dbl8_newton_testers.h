// The file dbl8_newton_testers.h specifies test function for Newton's method
// on series in octo double precision.

#ifndef __dbl8_newton_testers_h__
#define __dbl8_newton_testers_h__

void dbl8_unit_series_vector
 ( int dim, int deg,
   double **cffhihihi, double **cfflohihi,
   double **cffhilohi, double **cfflolohi,
   double **cffhihilo, double **cfflohilo,
   double **cffhilolo, double **cfflololo );
/*
 * DESCRIPTION :
 *   Given space in cff* for the doubles
 *   for the vector of dim power series in octo double precision,
 *   of series truncated after degree deg,
 *   initializes the coefficients in cff to one as leading coefficients,
 *   and zero for all other coefficients. */

void cmplx8_unit_series_vector
 ( int dim, int deg,
   double **cffrehihihi, double **cffrelohihi,
   double **cffrehilohi, double **cffrelolohi,
   double **cffrehihilo, double **cffrelohilo,
   double **cffrehilolo, double **cffrelololo,
   double **cffimhihihi, double **cffimlohihi,
   double **cffimhilohi, double **cffimlolohi,
   double **cffimhihilo, double **cffimlohilo,
   double **cffimhilolo, double **cffimlololo );
/*
 * DESCRIPTION :
 *   Given space allocated in the coefficient arrays,
 *   returns the values for a complex unit series. */

void dbl8_update_series
 ( int dim, int degp1,
   double **xhihihi, double **xlohihi, double **xhilohi, double **xlolohi,
   double **xhihilo, double **xlohilo, double **xhilolo, double **xlololo,
   double **dxhihihi, double **dxlohihi, double **dxhilohi, double **dxlolohi,
   double **dxhihilo, double **dxlohilo, double **dxhilolo, double **dxlololo,
   int vrblvl );
/*
 * DESCRIPTION :
 *   Adds the series in dx to x.
 *
 * ON ENTRY :
 *   dim       number of series in x and dx;
 *   degp1     degree plus one of the series;
 *   xhihihi   highest doubles of the series to updated, not linearized;
 *   xlohihi   2nd highest doubles of the series to updated, not linearized;
 *   xlohihi   3rd highest doubles of the series to updated, not linearized;
 *   xlolohi   4th highest doubles of the series to updated, not linearized;
 *   xhihilo   4th lowest doubles of the series to updated, not linearized;
 *   xlohilo   3rd lowest doubles of the series to updated, not linearized;
 *   xlohilo   2nd lowest doubles of the series to updated, not linearized;
 *   xlololo   lowest doubles of the series to updated, not linearized;
 *   dxhihihi  linearized highest doubles of the update of the series;
 *   dxlohihi  linearized 2nd highest doubles of the update of the series;
 *   dxhilohi  linearized 3rd highest doubles of the update of the series;
 *   dxlolohi  linearized 4th highest doubles of the update of the series;
 *   dxhihilo  linearized 4th lowest doubles of the update of the series;
 *   dxlohilo  linearized 3rd lowest doubles of the update of the series;
 *   dxhilolo  linearized 2nd lowest doubles of the update of the series;
 *   dxlololo  linearized lowest doubles of the update of the series;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   xhihihi   highest doubles of the series x updated with dx;
 *   xlohihi   second highest doubles of the series x updated with dx;
 *   xhilohi   third highest doubles of the series x updated with dx;
 *   xlolohi   fourth highest doubles of the series x updated with dx;
 *   xhihilo   fourth lowest doubles of the series x updated with dx;
 *   xlohilo   third lowest doubles of the series x updated with dx;
 *   xhilolo   second lowest doubles of the series x updated with dx;
 *   xlololo   lowest doubles of the series x updated with dx. */

void cmplx8_update_series
 ( int dim, int degp1,
   double **xrehihihi, double **xrelohihi,
   double **xrehilohi, double **xrelolohi,
   double **xrehihilo, double **xrelohilo,
   double **xrehilolo, double **xrelololo,
   double **ximhihihi, double **ximlohihi,
   double **ximhilohi, double **ximlolohi,
   double **ximhihilo, double **ximlohilo,
   double **ximhilolo, double **ximlololo,
   double **dxrehihihi, double **dxrelohihi,
   double **dxrehilohi, double **dxrelolohi,
   double **dxrehihilo, double **dxrelohilo,
   double **dxrehilolo, double **dxrelololo,
   double **dximhihihi, double **dximlohihi,
   double **dximhilohi, double **dximlolohi,
   double **dximhihilo, double **dximlohilo,
   double **dximhilolo, double **dximlololo, int vrblvl );
/*
 * DESCRIPTION :
 *   Adds the series in dx to x.
 *
 * ON ENTRY :
 *   dim       number of series in x and dx;
 *   degp1     degree plus one of the series;
 *   xrehihihi are the highest doubles of the real parts of the series x,
 *             x is not linearized;
 *   xrelohihi are the second highest doubles of the real parts of x;
 *   xrehilohi are the third highest doubles of the real parts of x;
 *   xrelolohi are the fourth lowest doubles of the real parts of x;
 *   xrelololo are the lowest doubles of the real parts of x;
 *   ximhihihi are the highest doubles of the imaginary parts of x;
 *   ximlohihi are the second highest doubles of the imaginary parts of x;
 *   ximlohihi are the third highest doubles of the imaginary parts of x;
 *   ximhilohi are the fourth highest doubles of the imaginary parts of x;
 *   ximlolohi are the fourth lowest doubles of the imaginary parts of x;
 *   ximhihilo are the third lowest doubles of the imaginary parts of x;
 *   ximlohilo are the second lowest doubles of the imagi parts of x;
 *   ximlololo are the lowest doubles of the imaginary parts of x;
 *   dxrehihihi are the highest doubles of the real parts
 *             of the linearized update dx;
 *   dxrelohihi are the second highest doubles of the real parts of dx;
 *   dxrehilohi are the third highest doubles of the real parts of dx;
 *   dxrelolohi are the fourth highest doubles of the real parts of dx;
 *   dxrehihilo are the fourth lowest doubles of the real parts of dx;
 *   dxrelohilo are the third lowest doubles of the real parts of dx;
 *   dxrehilolo are the second lowest doubles of the real parts of dx;
 *   dxrelololo are the lowest doubles of the real parts of dx;
 *   dximhihihi are the highest doubles of the imaginary parts of dx;
 *   dximlohihi are the second highest doubles of the imaginary parts of dx;
 *   dximhilohi are the third highest doubles of the imaginary parts of dx;
 *   dximlolohi are the fourth highest doubles of the imaginary parts of dx;
 *   dximhihilo are the fourth lowest doubles of the imaginary parts of dx;
 *   dximlohilo are the third lowest doubles of the imaginary parts of dx;
 *   dximhilolo are the second lowest doubles of the imaginary parts of dx;
 *   dximlololo are the lowest doubles of the imaginary parts of the dx;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   xrehihihi are the highest doubles of the real parts of the updated x;
 *   xrelohihi are the second highest doubles of the real parts of x;
 *   xrehilohi are the third highest doubles of the real parts of x;
 *   xrelolohi are the fourth highest doubles of the real parts of x;
 *   xrehihilo are the fourth lowest doubles of the real parts of x;
 *   xrelohilo are the third lowest doubles of the real parts of x;
 *   xrehilolo are the second lowest doubles of the real parts of x;
 *   xrelololo are the lowest doubles of the real parts of x;
 *   ximhihihi are the highest doubles fo the imaginary parts of x;
 *   ximlohihi are the second highest doubles fo the imaginary parts of x;
 *   ximhilohi are the third highest doubles fo the imaginary parts of x;
 *   ximlolohi are the fourth highest doubles fo the imaginary parts of x;
 *   ximhihilo are the fourth lowest doubles fo the imaginary parts of x;
 *   ximlohilo are the third lowest doubles fo the imaginary parts of x;
 *   ximhilolo are the second lowest doubles fo the imaginary parts of x;
 *   ximlololo are the lowest doubles fo the imaginary parts of x. */

void dbl8_newton_lustep
 ( int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffhihihi, double **cfflohihi,
   double **cffhilohi, double **cfflolohi,
   double **cffhihilo, double **cfflohilo,
   double **cffhilolo, double **cfflololo,
   double *acchihihi, double *acclohihi,
   double *acchilohi, double *acclolohi,
   double *acchihilo, double *acclohilo,
   double *acchilolo, double *acclololo,
   double **inputhihihi, double **inputlohihi,
   double **inputhilohi, double **inputlolohi,
   double **inputhihilo, double **inputlohilo,
   double **inputhilolo, double **inputlololo,
   double ***outputhihihi, double ***outputlohihi,
   double ***outputhilohi, double ***outputlolohi,
   double ***outputhihilo, double ***outputlohilo,
   double ***outputhilolo, double ***outputlololo,
   double **funvalhihihi, double **funvallohihi,
   double **funvalhilohi, double **funvallolohi,
   double **funvalhihilo, double **funvallohilo,
   double **funvalhilolo, double **funvallololo,
   double ***jacvalhihihi, double ***jacvallohihi,
   double ***jacvalhilohi, double ***jacvallolohi,
   double ***jacvalhihilo, double ***jacvallohilo,
   double ***jacvalhilolo, double ***jacvallololo,
   double **rhshihihi, double **rhslohihi,
   double **rhshilohi, double **rhslolohi,
   double **rhshihilo, double **rhslohilo,
   double **rhshilolo, double **rhslololo,
   double **solhihihi, double **sollohihi,
   double **solhilohi, double **sollolohi,
   double **solhihilo, double **sollohilo,
   double **solhilolo, double **sollololo,
   double **workmathihihi, double **workmatlohihi,
   double **workmathilohi, double **workmatlolohi,
   double **workmathihilo, double **workmatlohilo,
   double **workmathilolo, double **workmatlololo,
   double *workvechihihi, double *workveclohihi,
   double *workvechilohi, double *workveclolohi,
   double *workvechihilo, double *workveclohilo,
   double *workvechilolo, double *workveclololo,
   double **workrhshihihi, double **workrhslohihi,
   double **workrhshilohi, double **workrhslolohi,
   double **workrhshihilo, double **workrhslohilo,
   double **workrhshilolo, double **workrhslololo,
   double **resvechihihi, double **resveclohihi,
   double **resvechilohi, double **resveclolohi,
   double **resvechihilo, double **resveclohilo,
   double **resvechilolo, double **resveclololo,
   double *resmaxhihihi, double *resmaxlohihi,
   double *resmaxhilohi, double *resmaxlolohi,
   double *resmaxhihilo, double *resmaxlohilo,
   double *resmaxhilolo, double *resmaxlololo, int *ipvt, int vrblvl );
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
 *   cffhihihi has the highest doubles of the monomial coefficients;
 *   cfflohihi has the 2nd highest doubles of the monomial coefficients;
 *   cffhilohi has the 3rd highest doubles of the monomial coefficients;
 *   cfflolohi has the 4th highest doubles of the monomial coefficients;
 *   cffhihilo has the 4th lowest doubles of the monomial coefficients;
 *   cfflohilo has the 3rd lowest doubles of the monomial coefficients;
 *   cffhilolo has the 2nd lowest doubles of the monomial coefficients;
 *   cfflololo has the lowest doubles of the monomial coefficients;
 *   acchihihi has space to accumulate one power series of degree deg;
 *   acclohihi has space to accumulate one power series of degree deg;
 *   acchilohi has space to accumulate one power series of degree deg;
 *   acclolohi has space to accumulate one power series of degree deg;
 *   acchihilo has space to accumulate one power series of degree deg;
 *   acclohilo has space to accumulate one power series of degree deg;
 *   acchilolo has space to accumulate one power series of degree deg;
 *   acclololo has space to accumulate one power series of degree deg;
 *   inputhihihi are the highest doubles of the coefficients of the power
 *             series of degree deg, for dim variables;
 *   inputlohihi are the second highest doubles of the coefficients of the
 *             power series of degree deg, for dim variables;
 *   inputhilohi are the third highest doubles of the coefficients of the
 *             power series of degree deg, for dim variables;
 *   inputlolohi are the fourth highest doubles of the coefficients of the
 *             power series of degree deg, for dim variables;
 *   inputhihilo are the fourth lowest doubles of the coefficients of the
 *             power series of degree deg, for dim variables;
 *   inputlohilo are the third lowest doubles of the coefficients of the
 *             power series of degree deg, for dim variables;
 *   inputhilolo are the second lowest doubles of the coefficients of the
 *             power series of degree deg, for dim variables;
 *   inputlololo are the lowest doubles of the coefficients of the power
 *             series of degree deg, for dim variables;
 *   outputhihihi has space for the evaluated and differentiated monomials;
 *   outputlohihi has space for the evaluated and differentiated monomials;
 *   outputhilohi has space for the evaluated and differentiated monomials;
 *   outputlolohi has space for the evaluated and differentiated monomials;
 *   outputhihilo has space for the evaluated and differentiated monomials;
 *   outputlohilo has space for the evaluated and differentiated monomials;
 *   outputhilolo has space for the evaluated and differentiated monomials;
 *   outputlololo has space for the evaluated and differentiated monomials;
 *   funvalhihihi has space for the evaluated power series;
 *   funvallohihi has space for the evaluated power series;
 *   funvalhilohi has space for the evaluated power series;
 *   funvallolohi has space for the evaluated power series;
 *   funvalhihilo has space for the evaluated power series;
 *   funvallohilo has space for the evaluated power series;
 *   funvalhilolo has space for the evaluated power series;
 *   funvallololo has space for the evaluated power series;
 *   jacvalhihihi has space for deg+1 matrices of dimension dim;
 *   jacvallohihi has space for deg+1 matrices of dimension dim;
 *   jacvalhilohi has space for deg+1 matrices of dimension dim;
 *   jacvallolohi has space for deg+1 matrices of dimension dim;
 *   jacvalhihilo has space for deg+1 matrices of dimension dim;
 *   jacvallohilo has space for deg+1 matrices of dimension dim;
 *   jacvalhilolo has space for deg+1 matrices of dimension dim;
 *   jacvallololo has space for deg+1 matrices of dimension dim;
 *   rhshihihi has space for deg+1 vectors of dimension dim;
 *   rhslohihi has space for deg+1 vectors of dimension dim;
 *   rhshilohi has space for deg+1 vectors of dimension dim;
 *   rhslolohi has space for deg+1 vectors of dimension dim;
 *   rhshihilo has space for deg+1 vectors of dimension dim;
 *   rhslohilo has space for deg+1 vectors of dimension dim;
 *   rhshilolo has pace for deg+1 vectors of dimension dim;
 *   rhslololo has space for deg+1 vectors of dimension dim;
 *   solhihihi has space for deg+1 vectors of dimension dim;
 *   sollohihi has space for deg+1 vectors of dimension dim;
 *   solhilohi has space for deg+1 vectors of dimension dim;
 *   sollolohi has space for deg+1 vectors of dimension dim;
 *   solhihilo has space for deg+1 vectors of dimension dim;
 *   sollohilo has space for deg+1 vectors of dimension dim;
 *   solhilolo has space for deg+1 vectors of dimension dim;
 *   sollololo has space for deg+1 vectors of dimension dim;
 *   wrkmathihihi has work space allocated for a matrix of dimension dim;
 *   wrkmatlohihi has work space allocated for a matrix of dimension dim;
 *   wrkmathilohi has work space allocated for a matrix of dimension dim;
 *   wrkmatlolohi has work space allocated for a matrix of dimension dim;
 *   wrkmathihilo has work space allocated for a matrix of dimension dim;
 *   wrkmatlohilo has work space allocated for a matrix of dimension dim;
 *   wrkmathilolo has work space allocated for a matrix of dimension dim;
 *   wrkmatlololo has work space allocated for a matrix of dimension dim;
 *   wrkvechihihi has work space allocated for a vector of dimension dim;
 *   wrkveclohihi has work space allocated for a vector of dimension dim;
 *   wrkvechilohi has work space allocated for a vector of dimension dim;
 *   wrkveclolohi has work space allocated for a vector of dimension dim;
 *   wrkvechihilo has work space allocated for a vector of dimension dim;
 *   wrkveclohilo has work space allocated for a vector of dimension dim;
 *   wrkvechilolo has work space allocated for a vector of dimension dim;
 *   wrkveclololo has work space allocated for a vector of dimension dim;
 *   resvechihihi has space for deg+1 vectors of dimension dim;
 *   resveclohihi has space for deg+1 vectors of dimension dim;
 *   resvechilohi has space for deg+1 vectors of dimension dim;
 *   resveclolohi has space for deg+1 vectors of dimension dim;
 *   resvechihilo has space for deg+1 vectors of dimension dim;
 *   resveclohilo has space for deg+1 vectors of dimension dim;
 *   resvechilolo has space for deg+1 vectors of dimension dim;
 *   resveclololo has space for deg+1 vectors of dimension dim;
 *   ipvt      space allocated for dim pivots;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   funvalhihihi are the highest doubles of the output[i][dim];
 *   funvallohihi are the second highest doubles of the output[i][dim];
 *   funvalhilohi are the third highest doubles of the output[i][dim];
 *   funvallolohi are the fourth highest doubles of the output[i][dim];
 *   funvalhihilo are the fourth lowest doubles of the output[i][dim];
 *   funvallohilo are the third lowest doubles of the output[i][dim];
 *   funvalhilolo are the second lowest doubles of the output[i][dim];
 *   funvallololo are the lowest doubles of the output[i][dim];
 *   jacvalhihihi are the highest doubles of a matrix series,
 *             the leading coefficient is the Jacobian matrix.
 *   jacvallohihi are the second highest doubles of a matrix series,
 *             the leading coefficient is the Jacobian matrix.
 *   jacvalhilohi are the third highest doubles of a matrix series,
 *             the leading coefficient is the Jacobian matrix.
 *   jacvallolohi are the fourth highest doubles of a matrix series,
 *             the leading coefficient is the Jacobian matrix.
 *   jacvalhihilo are the fourth lowest doubles of a matrix series,
 *             the leading coefficient is the Jacobian matrix.
 *   jacvallohilo are the third lowest doubles of a matrix series,
 *             the leading coefficient is the Jacobian matrix.
 *   jacvalhilolo are the second lowest doubles of a matrix series,
 *             the leading coefficient is the Jacobian matrix.
 *   jacvallololo are the lowest doubles of a matrix series,
 *             the leading coefficient is the Jacobian matrix.
 *   rhshihihi are the highest doubles of the linearized right hand side,
 *             the function values subtracted by 1 and added by t;
 *   rhslohihi are the second highest doubles of the linearized right hand
 *             side, the function values subtracted by 1 and added by t;
 *   rhshilohi are the third highest doubles of the linearized right hand
 *             side, are the function values subtracted by 1 and added by t;
 *   rhslolohi are the fourth highest doubles ofthe linearized right hand
 *             side, are the function values subtracted by 1 and added by t;
 *   rhshihilo are the fourth lowest doubles of the linearized right hand
 *             side, the function values subtracted by 1 and added by t;
 *   rhslohilo are the third lowest doubles of the linearized right hand
 *             side, the function values subtracted by 1 and added by t;
 *   rhshilolo are the second lowest doubles of the linearized right hand
 *             side, the function values subtracted by 1 and added by t;
 *   rhslololo are the lowest doubles ofthe linearized right hand
 *             side, the function values subtracted by 1 and added by t;
 *   wrkmathihihi are the highest doubles of the LU factorization
 *             of the Jacobian matrix;
 *   wrkmatlohihi are the second highest doubles of the LU factorization
 *             of the Jacobian matrix;
 *   wrkmathilohi are the third highest doubles of the LU factorization
 *             of the Jacobian matrix;
 *   wrkmatlolohi are the fourth highest doubles of the LU factorization
 *             of the Jacobian matrix;
 *   wrkmathihilo are the fourth lowest doubles of the LU factorization
 *             of the Jacobian matrix;
 *   wrkmatlohilo are the third lowest doubles of the LU factorization
 *             of the Jacobian matrix;
 *   wrkmathilolo are the second lowest doubles of the LU factorization
 *             of the Jacobian matrix;
 *   wrkmatlololo are the lowest doubles of the LU factorization
 *             of the Jacobian matrix;
 *   resvechihihi are the highest doubles of the residual vectors;
 *   resveclohihi are the second highest doubles of the residual vectors;
 *   resvechilohi are the third highest doubles of the residual vectors;
 *   resveclolohi are the fourth highest doubles of the residual vectors;
 *   resvechihilo are the fourth lowest doubles of the residual vectors;
 *   resveclohilo are the third lowest doubles of the residual vectors;
 *   resvechilolo are the second lowest doubles of the residual vectors;
 *   resveclololo are the lowest doubles of the residual vectors;
 *   resmaxhihihi is the highest double of the maximum element
 *             of the residual vectors;
 *   resmaxlohihi is the second highest double of the maximum element
 *             of the residual vectors;
 *   resmaxhilohi is the third highest double of the maximum element
 *             of the residual vectors;
 *   resmaxlolohi is the fourth highest double of the maximum element
 *             of the residual vectors;
 *   resmaxhihilo is the fourth lowest double of the maximum element
 *             of the residual vectors;
 *   resmaxlohilo is the third lowest double of the maximum element
 *             of the residual vectors;
 *   resmaxhilolo is the second lowest double of the maximum element
 *             of the residual vectors;
 *   resmaxlololo is the lowest double of the maximum element
 *             of the residual vectors;
 *   ipvt      pivots used on the LU factorization of the lead matrix. */

void dbl8_newton_qrstep
 ( int szt, int nbt, int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffhihihi, double **cfflohihi,
   double **cffhilohi, double **cfflolohi,
   double **cffhihilo, double **cfflohilo,
   double **cffhilolo, double **cfflololo,
   double *acchihihi, double *acclohihi, double *acchilohi, double *acclolohi,
   double *acchihilo, double *acclohilo, double *acchilolo, double *acclololo,
   double **inputhihihi_h, double **inputlohihi_h,
   double **inputhilohi_h, double **inputlolohi_h,
   double **inputhihilo_h, double **inputlohilo_h,
   double **inputhilolo_h, double **inputlololo_h,
   double **inputhihihi_d, double **inputlohihi_d,
   double **inputhilohi_d, double **inputlolohi_d,
   double **inputhihilo_d, double **inputlohilo_d,
   double **inputhilolo_d, double **inputlololo_d,
   double ***outputhihihi_h, double ***outputlohihi_h,
   double ***outputhilohi_h, double ***outputlolohi_h,
   double ***outputhihilo_h, double ***outputlohilo_h,
   double ***outputhilolo_h, double ***outputlololo_h,
   double ***outputhihihi_d, double ***outputlohihi_d,
   double ***outputhilohi_d, double ***outputlolohi_d,
   double ***outputhihilo_d, double ***outputlohilo_d,
   double ***outputhilolo_d, double ***outputlololo_d,
   double **funvalhihihi_h, double **funvallohihi_h,
   double **funvalhilohi_h, double **funvallolohi_h,
   double **funvalhihilo_h, double **funvallohilo_h,
   double **funvalhilolo_h, double **funvallololo_h,
   double **funvalhihihi_d, double **funvallohihi_d,
   double **funvalhilohi_d, double **funvallolohi_d,
   double **funvalhihilo_d, double **funvallohilo_d,
   double **funvalhilolo_d, double **funvallololo_d,
   double ***jacvalhihihi_h, double ***jacvallohihi_h,
   double ***jacvalhilohi_h, double ***jacvallolohi_h,
   double ***jacvalhihilo_h, double ***jacvallohilo_h,
   double ***jacvalhilolo_h, double ***jacvallololo_h,
   double ***jacvalhihihi_d, double ***jacvallohihi_d,
   double ***jacvalhilohi_d, double ***jacvallolohi_d,
   double ***jacvalhihilo_d, double ***jacvallohilo_d,
   double ***jacvalhilolo_d, double ***jacvallololo_d,
   double **rhshihihi_h, double **rhslohihi_h,
   double **rhshilohi_h, double **rhslolohi_h,
   double **rhshihilo_h, double **rhslohilo_h,
   double **rhshilolo_h, double **rhslololo_h,
   double **rhshihihi_d, double **rhslohihi_d,
   double **rhshilohi_d, double **rhslolohi_d,
   double **rhshihilo_d, double **rhslohilo_d,
   double **rhshilolo_d, double **rhslololo_d,
   double **urhshihihi_h, double **urhslohihi_h,
   double **urhshilohi_h, double **urhslolohi_h,
   double **urhshihilo_h, double **urhslohilo_h,
   double **urhshilolo_h, double **urhslololo_h,
   double **urhshihihi_d, double **urhslohihi_d,
   double **urhshilohi_d, double **urhslolohi_d,
   double **urhshihilo_d, double **urhslohilo_d,
   double **urhshilolo_d, double **urhslololo_d,
   double **solhihihi_h, double **sollohihi_h,
   double **solhilohi_h, double **sollolohi_h,
   double **solhihilo_h, double **sollohilo_h,
   double **solhilolo_h, double **sollololo_h,
   double **solhihihi_d, double **sollohihi_d,
   double **solhilohi_d, double **sollolohi_d,
   double **solhihilo_d, double **sollohilo_d,
   double **solhilolo_d, double **sollololo_d,
   double **Qhihihi_h, double **Qlohihi_h,
   double **Qhilohi_h, double **Qlolohi_h,
   double **Qhihilo_h, double **Qlohilo_h,
   double **Qhilolo_h, double **Qlololo_h,
   double **Qhihihi_d, double **Qlohihi_d,
   double **Qhilohi_d, double **Qlolohi_d,
   double **Qhihilo_d, double **Qlohilo_d,
   double **Qhilolo_d, double **Qlololo_d,
   double **Rhihihi_h, double **Rlohihi_h,
   double **Rhilohi_h, double **Rlolohi_h,
   double **Rhihilo_h, double **Rlohilo_h,
   double **Rhilolo_h, double **Rlololo_h,
   double **Rhihihi_d, double **Rlohihi_d,
   double **Rhilohi_d, double **Rlolohi_d,
   double **Rhihilo_d, double **Rlohilo_d,
   double **Rhilolo_d, double **Rlololo_d,
   double **workmathihihi, double **workmatlohihi,
   double **workmathilohi, double **workmatlolohi,
   double **workmathihilo, double **workmatlohilo,
   double **workmathilolo, double **workmatlololo,
   double *workvechihihi, double *workveclohihi,
   double *workvechilohi, double *workveclolohi,
   double *workvechihilo, double *workveclohilo,
   double *workvechilolo, double *workveclololo,
   double **resvechihihi, double **resveclohihi, 
   double **resvechilohi, double **resveclolohi, 
   double **resvechihilo, double **resveclohilo, 
   double **resvechilolo, double **resveclololo, 
   double *resmaxhihihi, double *resmaxlohihi,
   double *resmaxhilohi, double *resmaxlolohi,
   double *resmaxhihilo, double *resmaxlohilo,
   double *resmaxhilolo, double *resmaxlololo, int vrblvl, int mode );
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
 *   cffhihihi are the highest doubles of the coefficients of the monomials;
 *   cfflohihi are the second highest doubles of the coefficients;
 *   cffhilohi are the thirnd highest doubles of the coefficients;
 *   cfflolohi are the fourth highest doubles of the coefficients;
 *   cffhihilo are the fourth lowest doubles of the coefficients;
 *   cfflohilo are the third lowest doubles of the coefficients;
 *   cffhilolo are the second lowest doubles of the coefficients;
 *   cfflololo are the lowest doubles of the coefficients of the monomials;
 *   acchihihi has has space for highest doubles of one series of degree deg;
 *   acclohihi has space for second highest doubles of one series;
 *   acchilohi has space for third highest doubles of one series;
 *   acclolohi has space for fourth highest doubles of one series;
 *   acchihilo has space for fourth lowest doubles of one series;
 *   acclohilo has space for third lowest doubles of one series;
 *   acchilolo has space for second lowest doubles of one series;
 *   acclololo has space for lowest doubles of one series of degree deg;
 *   inputhihihi_h has the highest doubles of the input series,
 *             of degree deg, for dim variables, computed on host;
 *   inputlohihi_h has the second highest doubles of the input series,
 *             of degree deg, for dim variables, computed on host;
 *   inputhilohi_h has the third highest doubles of the input series,
 *             of degree deg, for dim variables, computed on host;
 *   inputlolohi_h has the fourth highest doubles of the input series,
 *             of degree deg, for dim variables, computed on host;
 *   inputhihilo_h has the fourth lowest doubles of the input series,
 *             of degree deg, for dim variables, computed on host;
 *   inputlohilo_h has the third lowest doubles of the input series,
 *             of degree deg, for dim variables, computed on host;
 *   inputhilolo_h has the second lowest doubles of the input series,
 *             of degree deg, for dim variables, computed on host;
 *   inputlololo_h has the lowest doubles of the input series,
 *             of degree deg, for dim variables, computed on host;
 *   inputhihihi_d has space for series computed on device;
 *   inputlohihi_d has space for series computed on device;
 *   inputhilohi_d has space for series computed on device;
 *   inputlolohi_d has space for series computed on device;
 *   inputhihilo_d has space for series computed on device;
 *   inputlohilo_d has space for series computed on device;
 *   inputhilolo_d has space for series computed on device;
 *   inputlololo_d has space for series computed on device;
 *   outputhihihi_h has space for the highest doubles of the evaluated
 *             and differentiated monomials, computed on the host;
 *   outputlohihi_h has space for the second highest doubles of the evaluated 
 *             and differentiated monomials, computed on the host;
 *   outputhilohi_h has space for the third highest doubles of the evaluated
 *             and differentiated monomials, computed on the host;
 *   outputlolohi_h has space for the fourth highest doubles of the evaluated 
 *             and differentiated monomials, computed on the host;
 *   outputhihilo_h has space for the fourth lowest doubles of the evaluated
 *             and differentiated monomials, computed on the host;
 *   outputlohilo_h has space for the third lowest doubles of the evaluated 
 *             and differentiated monomials, computed on the host;
 *   outputhilolo_h has space for the second lowest doubles of the evaluated
 *             and differentiated monomials, computed on the host;
 *   outputlololo_h has space for the lowest doubles of the evaluated 
 *             and differentiated monomials, computed on the host;
 *   outputhihihi_d has space for the output, computed on the device;
 *   outputlohihi_d has space for the output, computed on the device;
 *   outputhilohi_d has space for the output, computed on the device;
 *   outputlolohi_d has space for the output, computed on the device;
 *   outputhihilo_d has space for the output, computed on the device;
 *   outputlohilo_d has space for the output, computed on the device;
 *   outputhilolo_d has space for the output, computed on the device;
 *   outputlololo_d has space for the output, computed on the device;
 *   funvalhihihi_h has space for the evaluated series on the host;
 *   funvallohihi_h has space for the evaluated series on the host;
 *   funvalhilohi_h has space for the evaluated series on the host;
 *   funvallolohi_h has space for the evaluated series on the host;
 *   funvalhihilo_h has space for the evaluated series on the host;
 *   funvallohilo_h has space for the evaluated series on the host;
 *   funvalhilolo_h has space for the evaluated series on the host;
 *   funvallololo_h has space for the evaluated series on the host;
 *   funvalhihihi_d has space for the evaluated series on the device;
 *   funvallohihi_d has space for the evaluated series on the device;
 *   funvalhilohi_d has space for the evaluated series on the device;
 *   funvallolohi_d has space for the evaluated series on the device;
 *   funvalhihilo_d has space for the evaluated series on the device;
 *   funvallohilo_d has space for the evaluated series on the device;
 *   funvalhilolo_d has space for the evaluated series on the device;
 *   funvallololo_d has space for the evaluated series on the device;
 *   jacvalhihihi_h has space for deg+1 matrices of dimension dim on host;
 *   jacvallohihi_h has space for deg+1 matrices of dimension dim on host;
 *   jacvalhilohi_h has space for deg+1 matrices of dimension dim on host;
 *   jacvallolohi_h has space for deg+1 matrices of dimension dim on host;
 *   jacvalhihilo_h has space for deg+1 matrices of dimension dim on host;
 *   jacvallohilo_h has space for deg+1 matrices of dimension dim on host;
 *   jacvalhilolo_h has space for deg+1 matrices of dimension dim on host;
 *   jacvallololo_h has space for deg+1 matrices of dimension dim on host;
 *   jacvalhihihi_d has space for deg+1 matrices of dimension dim on device;
 *   jacvallohihi_d has space for deg+1 matrices of dimension dim on device;
 *   jacvalhilohi_d has space for deg+1 matrices of dimension dim on device;
 *   jacvallolohi_d has space for deg+1 matrices of dimension dim on device;
 *   jacvalhihilo_d has space for deg+1 matrices of dimension dim on device;
 *   jacvallohilo_d has space for deg+1 matrices of dimension dim on device;
 *   jacvalhilolo_d has space for deg+1 matrices of dimension dim on device;
 *   jacvallololo_d has space for deg+1 matrices of dimension dim on device;
 *   rhshihihi_h has space for deg+1 vectors of dimension dim on host;
 *   rhslohihi_h has space for deg+1 vectors of dimension dim on host;
 *   rhshilohi_h has space for deg+1 vectors of dimension dim on host;
 *   rhslolohi_h has space for deg+1 vectors of dimension dim on host;
 *   rhshihilo_h has space for deg+1 vectors of dimension dim on host;
 *   rhslohilo_h has space for deg+1 vectors of dimension dim on host;
 *   rhshilolo_h has space for deg+1 vectors of dimension dim on host;
 *   rhslololo_h has space for deg+1 vectors of dimension dim on host;
 *   rhshihihi_d has space for deg+1 vectors of dimension dim on device;
 *   rhslohihi_d has space for deg+1 vectors of dimension dim on device;
 *   rhshilohi_d has space for deg+1 vectors of dimension dim on device;
 *   rhslolohi_d has space for deg+1 vectors of dimension dim on device;
 *   rhshihilo_d has space for deg+1 vectors of dimension dim on device;
 *   rhslohilo_d has space for deg+1 vectors of dimension dim on device;
 *   rhshilolo_d has space for deg+1 vectors of dimension dim on device;
 *   rhslololo_d has space for deg+1 vectors of dimension dim on device;
 *   urhshihihi_h has space for updated right hand side vectors on host;
 *   urhslohihi_h has space for updated right hand side vectors on host;
 *   urhshilohi_h has space for updated right hand side vectors on host;
 *   urhslolohi_h has space for updated right hand side vectors on host;
 *   urhshihilo_h has space for updated right hand side vectors on host;
 *   urhslohilo_h has space for updated right hand side vectors on host;
 *   urhshilolo_h has space for updated right hand side vectors on host;
 *   urhslololo_h has space for updated right hand side vectors on host;
 *   urhshihihi_d has space for updated right hand side vectors on device; 
 *   urhslohihi_d has space for updated right hand side vectors on device; 
 *   urhshilohi_d has space for updated right hand side vectors on device; 
 *   urhslolohi_d has space for updated right hand side vectors on device; 
 *   urhshihilo_d has space for updated right hand side vectors on device; 
 *   urhslohilo_d has space for updated right hand side vectors on device; 
 *   urhshilolo_d has space for updated right hand side vectors on device; 
 *   urhslololo_d has space for updated right hand side vectors on device; 
 *   solhihihi_h has space for deg+1 vectors of dimension dim;
 *   sollohihi_h has space for deg+1 vectors of dimension dim;
 *   solhilohi_h has space for deg+1 vectors of dimension dim;
 *   sollolohi_h has space for deg+1 vectors of dimension dim;
 *   solhihilo_h has space for deg+1 vectors of dimension dim;
 *   sollohilo_h has space for deg+1 vectors of dimension dim;
 *   solhilolo_h has space for deg+1 vectors of dimension dim;
 *   sollololo_h has space for deg+1 vectors of dimension dim;
 *   solhihihi_d has space for deg+1 vectors of dimension dim;
 *   sollohihi_d has space for deg+1 vectors of dimension dim;
 *   solhilohi_d has space for deg+1 vectors of dimension dim;
 *   sollolohi_d has space for deg+1 vectors of dimension dim;
 *   solhihilo_d has space for deg+1 vectors of dimension dim;
 *   sollohilo_d has space for deg+1 vectors of dimension dim;
 *   solhilolo_d has space for deg+1 vectors of dimension dim;
 *   sollololo_d has space for deg+1 vectors of dimension dim;
 *   Qhihihi_h has space allocated for the Q computed by the host;
 *   Qlohihi_h has space allocated for the Q computed by the host;
 *   Qhilohi_h has space allocated for the Q computed by the host;
 *   Qlolohi_h has space allocated for the Q computed by the host;
 *   Qhihilo_h has space allocated for the Q computed by the host;
 *   Qlohilo_h has space allocated for the Q computed by the host;
 *   Qhilolo_h has space allocated for the Q computed by the host;
 *   Qlololo_h has space allocated for the Q computed by the host;
 *   Qhihihi_d has space allocated for the Q computed by the device;
 *   Qlohihi_d has space allocated for the Q computed by the device;
 *   Qhilohi_d has space allocated for the Q computed by the device;
 *   Qlolohi_d has space allocated for the Q computed by the device;
 *   Qhihilo_d has space allocated for the Q computed by the device;
 *   Qlohilo_d has space allocated for the Q computed by the device;
 *   Qhilolo_d has space allocated for the Q computed by the device;
 *   Qlololo_d has space allocated for the Q computed by the device;
 *   Rhihihi_h has space allocated for the R computed by the host;
 *   Rlohihi_h has space allocated for the R computed by the host;
 *   Rhilohi_h has space allocated for the R computed by the host;
 *   Rlolohi_h has space allocated for the R computed by the host;
 *   Rhihilo_h has space allocated for the R computed by the host;
 *   Rlohilo_h has space allocated for the R computed by the host;
 *   Rhilolo_h has space allocated for the R computed by the host;
 *   Rlololo_h has space allocated for the R computed by the host;
 *   Rhihihi_d has space allocated for the R computed by the device;
 *   Rlohihi_d has space allocated for the R computed by the device;
 *   Rhilohi_d has space allocated for the R computed by the device;
 *   Rlolohi_d has space allocated for the R computed by the device;
 *   Rhihilo_d has space allocated for the R computed by the device;
 *   Rlohilo_d has space allocated for the R computed by the device;
 *   Rhilolo_d has space allocated for the R computed by the device;
 *   Rlololo_d has space allocated for the R computed by the device;
 *   wrkmathihihi has work space allocated for a matrix of dimension dim;
 *   wrkmatlohihi has work space allocated for a matrix of dimension dim;
 *   wrkmathilohi has work space allocated for a matrix of dimension dim;
 *   wrkmatlolohi has work space allocated for a matrix of dimension dim;
 *   wrkmathihilo has work space allocated for a matrix of dimension dim;
 *   wrkmatlohilo has work space allocated for a matrix of dimension dim;
 *   wrkmathilolo has work space allocated for a matrix of dimension dim;
 *   wrkmatlololo has work space allocated for a matrix of dimension dim;
 *   wrkvechihihi has work space allocated for a vector of dimension dim;
 *   wrkveclohihi has work space allocated for a vector of dimension dim;
 *   wrkvechilohi has work space allocated for a vector of dimension dim;
 *   wrkveclolohi has work space allocated for a vector of dimension dim;
 *   wrkvechihilo has work space allocated for a vector of dimension dim;
 *   wrkveclohilo has work space allocated for a vector of dimension dim;
 *   wrkvechilolo has work space allocated for a vector of dimension dim;
 *   wrkveclololo has work space allocated for a vector of dimension dim;
 *   resvechihihi has space for deg+1 vectors of dimension dim;
 *   resveclohihi has space for deg+1 vectors of dimension dim;
 *   resvechilohi has space for deg+1 vectors of dimension dim;
 *   resveclolohi has space for deg+1 vectors of dimension dim;
 *   resvechihilo has space for deg+1 vectors of dimension dim;
 *   resveclohilo has space for deg+1 vectors of dimension dim;
 *   resvechilolo has space for deg+1 vectors of dimension dim;
 *   resveclololo has space for deg+1 vectors of dimension dim;
 *   vrblvl    is the verbose level;
 *   mode      execution mode, 0 (GPU only), 1 (CPU only) or 2 (GPU+CPU).
 *
 * ON RETURN :
 *   inputhihihi_h has the highest doubles of series, on host (mode 1, 2);
 *   inputlohihi_h has the second highest doubles of series computed on host;
 *   inputhilohi_h has the third highest doubles of series computed on host;
 *   inputlolohi_h has the fourth highest doubles of series computed on host;
 *   inputhihilo_h has the fourth lowest doubles of series, on host;
 *   inputlohilo_h has the third lowest doubles of series computed on host;
 *   inputhilolo_h has the second lowest doubles of series computed on host;
 *   inputlololo_h has the lowest doubles of series computed on host;
 *   inputhihihi_d has the highest doubles of series, on device (mode 0, 2);
 *   inputlohihi_d has the second highest doubles of series, on device;
 *   inputhilohi_d has the third highest doubles of series, on device;
 *   inputlolohi_d has the fourth highest doubles of series, on device;
 *   inputhihihi_d has the fourth lowest doubles of series, on device;
 *   inputlohihi_d has the third lowest doubles of series, on device;
 *   inputhilohi_d has the second lowest doubles of series, on device;
 *   inputlolohi_d has the lowest doubles of series, on device;
 *   outputhihi_h has the highest doubles of evaluated series, on host;
 *   outputlohi_h has the second highest doubles of series, on host;
 *   outputhilo_h has the third highest doubles of series, on host;
 *   outputlolo_h has the fourth highest doubles of series, on host;
 *   outputhihi_h has the fourth lowest doubles of series, on host;
 *   outputlohi_h has the third lowest doubles of series, on host;
 *   outputhilo_h has the second lowest doubles of series, on host;
 *   outputlolo_h has the lowest doubles of series, on host;
 *   outputhihi_d has the highest doubles of series, on device;
 *   outputlohi_d has the second highest doubles of series, on device;
 *   outputhilo_d has the third highest doubles of series, on device;
 *   outputlolo_d has the fourth highest doubles of series, on device;
 *   outputhihi_d has the fourth lowest doubles of series, on device;
 *   outputlohi_d has the third lowest doubles of series, on device;
 *   outputhilo_d has the second lowest doubles of series, on device;
 *   outputlolo_d has the lowest doubles of series, on device;
 *   funvalhihihi_h has the highest doubles of output[i][dim], on host;
 *   funvallohihi_h has the second highest doubles of output[i][dim], on host;
 *   funvalhilohi_h has the third highest doubles of output[i][dim], on host;
 *   funvallolohi_h has the fourth highest doubles of output[i][dim], on host;
 *   funvalhihilo_h has the fourth lowest doubles of output[i][dim], on host;
 *   funvallohilo_h has the third lowest doubles of output[i][dim], on host;
 *   funvalhilolo_h has the second lowest doubles of output[i][dim], on host;
 *   funvallololo_h has the lowest doubles of output[i][dim], on host;
 *   funvalhihihi_d has the highest doubles of output[i][dim], on device;
 *   funvallohihi_d has the 2nd highest doubles of output[i][dim], on device;
 *   funvalhilohi_d has the 3rd highest doubles of output[i][dim], on device;
 *   funvallolohi_d has the 4th highest doubles of output[i][dim], on device;
 *   funvalhihilo_d has the 4th lowest doubles of output[i][dim], on device;
 *   funvallohilo_d has the 3rd lowest doubles of output[i][dim], on device;
 *   funvalhilolo_d has the 2nd lowest doubles of output[i][dim], on device;
 *   funvallololo_d has the lowest doubles of output[i][dim], on device;
 *   jacvalhihihi_h has the highest doubles of matrices, on host,
 *             the leading coefficient is the Jacobian matrix;
 *   jacvallohihi_h has the second highest doubles of matrices, on host;
 *   jacvalhilohi_h has the third highest doubles of matrices, on host;
 *   jacvallolohi_h has the fourth highest doubles of matrices, on host;
 *   jacvalhihilo_h has the fourth lowest doubles of matrices, on host;
 *   jacvallohilo_h has the third lowest doubles of matrices, on host;
 *   jacvalhilolo_h has the second lowest doubles of matrices, on host;
 *   jacvallololo_h has the lowest doubles of matrices, on host;
 *   jacvalhihihi_d has the highest doubles of matrices, on device,
 *             the leading coefficient is the Jacobian matrix;
 *   jacvallohihi_d has the second highest doubles of matrices, on device;
 *   jacvalhilohi_d has the third highest doubles of matrices, on device;
 *   jacvallolohi_d has the fourth highest doubles of matrices, on device;
 *   jacvalhihilo_d has the fourth lowest doubles of matrices, on device;
 *   jacvallohilo_d has the third lowest doubles of matrices, on device;
 *   jacvalhilolo_d has the second lowest doubles of matrices, on device:
 *   jacvallololo_d has the lowest doubles of matrices, on device;
 *   rhshihihi_h has the highest doubles of the right hand side, linearized,
 *             subtracted by 1 and added by t, computed by host;
 *   rhslohihi_h has the 2nd highest doubles of the right hand side, on host;
 *   rhshilohi_h has the 3rd highest doubles of the right hand side, on host;
 *   rhslolohi_h has the 4th highest doubles of the right hand side, on host;
 *   rhshihilo_h has the 4th lowest doubles of the right hand side, on host;
 *   rhslohilo_h has the 3rd lowest doubles of the right hand side, on host;
 *   rhshilolo_h has the 2nd lowest doubles of the right hand side, on host;
 *   rhslololo_h has the lowest doubles of the right hand side, on host;
 *   rhshihi_d has the highest doubles of the right hand side, linearized,
 *             subtracted by 1 and added by t, computed by device;
 *   rhslohi_d has the 2nd highest doubles of the right hand side, on device;
 *   rhshilo_d has the 3rd highest doubles of the right hand side, on device;
 *   rhslolo_d has the 4th highest doubles of the right hand side, on device;
 *   rhshihi_d has the 4th lowest doubles of the right hand side, on device;
 *   rhslohi_d has the 3rd lowest doubles of the right hand side, on device;
 *   rhshilo_d has the 2nd lowest doubles of the right hand side, on device;
 *   rhslolo_d has the lowest doubles of the right hand side, on device;
 *   urhshihihi_h has the highest doubles of updated right hand side, on host;
 *   urhslohihi_h has the second highest doubles of right hand side, on host;
 *   urhshilohi_h has the third highest doubles of right hand side, on host;
 *   urhslolohi_h has the fourth highest doubles of right hand side, on host;
 *   urhshihilo_h has the fourth lowest doubles of right hand side, on host;
 *   urhslohilo_h has the third lowest doubles of right hand side, on host;
 *   urhshilolo_h has the second lowest doubles of right hand side, on host;
 *   urhslololo_h has the lowest doubles of updated right hand side, on host;
 *   urhshihihi_d has the highest doubles of right hand side, on device;
 *   urhslohihi_d has the second highest doubles of right hand side, on device;
 *   urhshilohi_d has the third highest doubles of right hand side, on device;
 *   urhslolohi_d has the fourth highest doubles of right hand side, on device;
 *   urhshihilo_d has the fourth lowest doubles of right hand side, on device;
 *   urhslohilo_d has the third lowest doubles of right hand side, on device;
 *   urhshilolo_d has the second lowest doubles of right hand side, on device;
 *   urhslololo_d has the lowest doubles of right hand side, on device;
 *   solhihihi_h has the highest doubles of solution, on the host;
 *   sollohihi_h has the second highest doubles of solution, on the host;
 *   solhilohi_h has the third highest doubles of solution, on the host;
 *   sollolohi_h has the fourth highest doubles of solution, on the host;
 *   solhihilo_h has the fourth lowest doubles of solution, on the host;
 *   sollohilo_h has the third lowest doubles of solution, on the host;
 *   solhilolo_h has the second lowest doubles of solution, on the host;
 *   sollololo_h has the lowest doubles of solution, on the host;
 *   solhihihi_d has the highest doubles of solution, on the device;
 *   sollohihi_d has the second highest doubles of solution, on the device;
 *   solhilohi_d has the third highest doubles of solution, on the device;
 *   sollolohi_d has the fourth highest doubles of solution, on the device;
 *   solhihilo_d has the fourth lowest doubles of solution, on the device;
 *   sollohilo_d has the third lowest doubles of solution, on the device;
 *   solhilolo_d has the second lowest doubles of solution, on the device;
 *   sollololo_d has the lowest doubles of solution, on the device;
 *   Qhihihi_h are the highest doubles of Q of the QR, on the host;
 *   Qlohihi_h are the second highest doubles of Q of the QR, on the host;
 *   Qhilohi_h are the third highest doubles of Q of the QR, on the host;
 *   Qlolohi_h are the fourth highest doubles of Q of the QR, on the host;
 *   Qhihilo_h are the fourth lowest doubles of Q of the QR, on the host;
 *   Qlohilo_h are the third lowest doubles of Q of the QR, on the host;
 *   Qhilolo_h are the second lowest doubles of Q of the QR, on the host;
 *   Qlololo_h are the lowest doubles of Q of the QR, on the host;
 *   Qhihihi_d are the highest doubles of Q of the QR, on the device;
 *   Qlohihi_d are the second highest doubles of Q of the QR, on the device;
 *   Qhilohi_d are the third highest doubles of Q of the QR, on the device;
 *   Qlolohi_d are the fourth highest doubles of Q of the QR, on the device;
 *   Qhihilo_d are the fourth lowest doubles of Q of the QR, on the device;
 *   Qlohilo_d are the third lowest doubles of Q of the QR, on the device;
 *   Qhilolo_d are the second lowest doubles of Q of the QR, on the device;
 *   Qlololo_d are the lowest doubles of Q of the QR, on the device;
 *   Rhihihi_h are the highest doubles of R of the QR, on the host;
 *   Rlohihi_h are the second highest doubles of R of the QR, on the host;
 *   Rhilohi_h are the third highest doubles of R of the QR, on the host;
 *   Rlolohi_h are the fourth highest doubles of R of the QR, on the host;
 *   Rhihilo_h are the fourth lowest doubles of R of the QR, on the host;
 *   Rlohilo_h are the third lowest doubles of R of the QR, on the host;
 *   Rhilolo_h are the second lowest doubles of R of the QR, on the host;
 *   Rlololo_h are the lowest doubles of R of the QR, on the host;
 *   Rhihihi_d are the highest doubles of R of the QR, on the device;
 *   Rlohihi_d are the second highest doubles of R of the QR, on the device;
 *   Rhilohi_d are the third highest doubles of R of the QR, on the device;
 *   Rlolohi_d are the fourth highest doubles of R of the QR, on the device;
 *   Rhihilo_d are the fourth lowest doubles of R of the QR, on the device;
 *   Rlohilo_d are the third lowhest doubles of R of the QR, on the device;
 *   Rhilolo_d are the second lowest doubles of R of the QR, on the device;
 *   Rlololo_d are the lowest doubles of R of the QR, on the device;
 *   wrkmathihihi has a copy of the highest doubles of the Jacobian;
 *   wrkmatlohihi has a copy of the second highest doubles of the Jacobian;
 *   wrkmathilohi has a copy of the third highest doubles of the Jacobian;
 *   wrkmatlolohi has a copy of the fourth highest doubles of the Jacobian;
 *   wrkmathihilo has a copy of the fourth lowest doubles of the Jacobian;
 *   wrkmatlohilo has a copy of the third lowest doubles of the Jacobian;
 *   wrkmathilolo has a copy of the second lowest doubles of the Jacobian;
 *   wrkmatlololo has a copy of the lowest doubles of the Jacobian matrix;
 *   resvechihihi has the highest doubles of the residual vectors;
 *   resveclohihi has the second highest doubles of the residual vectors;
 *   resvechilohi has the third highest doubles of the residual vectors;
 *   resveclolohi has the fourth highest doubles of the residual vectors;
 *   resvechihilo has the fourth lowest doubles of the residual vectors;
 *   resveclohilo has the third lowest doubles of the residual vectors;
 *   resvechilolo has the second lowest doubles of the residual vectors;
 *   resveclololo has the lowest doubles of the residual vectors;
 *   resmaxhihihi is the highest double of the maximum of residual vectors;
 *   resmaxlohihi is the second highest double of the maximum of residual;
 *   resmaxhilohi is the third highest double of the maximum of residual;
 *   resmaxlolohi is the fourth highest double of the maximum of residual;
 *   resmaxhihilo is the fourth lowest double of the maximum of residual;
 *   resmaxlohilo is the third lowest double of the maximum of residual;
 *   resmaxhilolo is the second lowest double of the maximum of residual;
 *   resmaxlololo is the lowest double of the maximum of residual. */

void cmplx8_newton_qrstep
 ( int szt, int nbt, int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffrehihihi, double **cffrelohihi,
   double **cffrehilohi, double **cffrelolohi,
   double **cffrehihilo, double **cffrelohilo,
   double **cffrehilolo, double **cffrelololo,
   double **cffimhihihi, double **cffimlohihi,
   double **cffimhilohi, double **cffimlolohi,
   double **cffimhihilo, double **cffimlohilo,
   double **cffimhilolo, double **cffimlololo,
   double *accrehihihi, double *accrelohihi,
   double *accrehilohi, double *accrelolohi,
   double *accrehihilo, double *accrelohilo,
   double *accrehilolo, double *accrelololo,
   double *accimhihihi, double *accimlohihi,
   double *accimhilohi, double *accimlolohi,
   double *accimhihilo, double *accimlohilo,
   double *accimhilolo, double *accimlololo,
   double **inputrehihihi_h, double **inputrelohihi_h,
   double **inputrehilohi_h, double **inputrelolohi_h,
   double **inputrehihilo_h, double **inputrelohilo_h,
   double **inputrehilolo_h, double **inputrelololo_h,
   double **inputimhihihi_h, double **inputimlohihi_h,
   double **inputimhilohi_h, double **inputimlolohi_h,
   double **inputimhihilo_h, double **inputimlohilo_h,
   double **inputimhilolo_h, double **inputimlololo_h,
   double **inputrehihihi_d, double **inputrelohihi_d,
   double **inputrehilohi_d, double **inputrelolohi_d,
   double **inputrehihilo_d, double **inputrelohilo_d,
   double **inputrehilolo_d, double **inputrelololo_d,
   double **inputimhihihi_d, double **inputimlohihi_d,
   double **inputimhilohi_d, double **inputimlolohi_d,
   double **inputimhihilo_d, double **inputimlohilo_d,
   double **inputimhilolo_d, double **inputimlololo_d,
   double ***outputrehihihi_h, double ***outputrelohihi_h,
   double ***outputrehilohi_h, double ***outputrelolohi_h,
   double ***outputrehihilo_h, double ***outputrelohilo_h,
   double ***outputrehilolo_h, double ***outputrelololo_h,
   double ***outputimhihihi_h, double ***outputimlohihi_h,
   double ***outputimhilohi_h, double ***outputimlolohi_h,
   double ***outputimhihilo_h, double ***outputimlohilo_h,
   double ***outputimhilolo_h, double ***outputimlololo_h,
   double ***outputrehihihi_d, double ***outputrelohihi_d,
   double ***outputrehilohi_d, double ***outputrelolohi_d,
   double ***outputrehihilo_d, double ***outputrelohilo_d,
   double ***outputrehilolo_d, double ***outputrelololo_d,
   double ***outputimhihihi_d, double ***outputimlohihi_d,
   double ***outputimhilohi_d, double ***outputimlolohi_d,
   double ***outputimhihilo_d, double ***outputimlohilo_d,
   double ***outputimhilolo_d, double ***outputimlololo_d,
   double **funvalrehihihi_h, double **funvalrelohihi_h,
   double **funvalrehilohi_h, double **funvalrelolohi_h,
   double **funvalrehihilo_h, double **funvalrelohilo_h,
   double **funvalrehilolo_h, double **funvalrelololo_h,
   double **funvalimhihihi_h, double **funvalimlohihi_h,
   double **funvalimhilohi_h, double **funvalimlolohi_h,
   double **funvalimhihilo_h, double **funvalimlohilo_h,
   double **funvalimhilolo_h, double **funvalimlololo_h,
   double **funvalrehihihi_d, double **funvalrelohihi_d,
   double **funvalrehilohi_d, double **funvalrelolohi_d,
   double **funvalrehihilo_d, double **funvalrelohilo_d,
   double **funvalrehilolo_d, double **funvalrelololo_d,
   double **funvalimhihihi_d, double **funvalimlohihi_d,
   double **funvalimhilohi_d, double **funvalimlolohi_d,
   double **funvalimhihilo_d, double **funvalimlohilo_d,
   double **funvalimhilolo_d, double **funvalimlololo_d,
   double ***jacvalrehihihi_h, double ***jacvalrelohihi_h,
   double ***jacvalrehilohi_h, double ***jacvalrelolohi_h,
   double ***jacvalrehihilo_h, double ***jacvalrelohilo_h,
   double ***jacvalrehilolo_h, double ***jacvalrelololo_h,
   double ***jacvalimhihihi_h, double ***jacvalimlohihi_h,
   double ***jacvalimhilohi_h, double ***jacvalimlolohi_h,
   double ***jacvalimhihilo_h, double ***jacvalimlohilo_h,
   double ***jacvalimhilolo_h, double ***jacvalimlololo_h,
   double ***jacvalrehihihi_d, double ***jacvalrelohihi_d,
   double ***jacvalrehilohi_d, double ***jacvalrelolohi_d,
   double ***jacvalrehihilo_d, double ***jacvalrelohilo_d,
   double ***jacvalrehilolo_d, double ***jacvalrelololo_d,
   double ***jacvalimhihihi_d, double ***jacvalimlohihi_d,
   double ***jacvalimhilohi_d, double ***jacvalimlolohi_d,
   double ***jacvalimhihilo_d, double ***jacvalimlohilo_d,
   double ***jacvalimhilolo_d, double ***jacvalimlololo_d,
   double **rhsrehihihi_h, double **rhsrelohihi_h,
   double **rhsrehilohi_h, double **rhsrelolohi_h,
   double **rhsrehihilo_h, double **rhsrelohilo_h,
   double **rhsrehilolo_h, double **rhsrelololo_h,
   double **rhsimhihihi_h, double **rhsimlohihi_h,
   double **rhsimhilohi_h, double **rhsimlolohi_h,
   double **rhsimhihilo_h, double **rhsimlohilo_h,
   double **rhsimhilolo_h, double **rhsimlololo_h,
   double **rhsrehihihi_d, double **rhsrelohihi_d, 
   double **rhsrehilohi_d, double **rhsrelolohi_d, 
   double **rhsrehihilo_d, double **rhsrelohilo_d, 
   double **rhsrehilolo_d, double **rhsrelololo_d, 
   double **rhsimhihihi_d, double **rhsimlohihi_d,
   double **rhsimhilohi_d, double **rhsimlolohi_d,
   double **rhsimhihilo_d, double **rhsimlohilo_d,
   double **rhsimhilolo_d, double **rhsimlololo_d,
   double **urhsrehihihi_h, double **urhsrelohihi_h,
   double **urhsrehilohi_h, double **urhsrelolohi_h,
   double **urhsrehihilo_h, double **urhsrelohilo_h,
   double **urhsrehilolo_h, double **urhsrelololo_h,
   double **urhsimhihihi_h, double **urhsimlohihi_h,
   double **urhsimhilohi_h, double **urhsimlolohi_h,
   double **urhsimhihilo_h, double **urhsimlohilo_h,
   double **urhsimhilolo_h, double **urhsimlololo_h,
   double **urhsrehihihi_d, double **urhsrelohihi_d,
   double **urhsrehilohi_d, double **urhsrelolohi_d,
   double **urhsrehihilo_d, double **urhsrelohilo_d,
   double **urhsrehilolo_d, double **urhsrelololo_d,
   double **urhsimhihihi_d, double **urhsimlohihi_d,
   double **urhsimhilohi_d, double **urhsimlolohi_d,
   double **urhsimhihilo_d, double **urhsimlohilo_d,
   double **urhsimhilolo_d, double **urhsimlololo_d,
   double **solrehihihi_h, double **solrelohihi_h,
   double **solrehilohi_h, double **solrelolohi_h,
   double **solrehihilo_h, double **solrelohilo_h,
   double **solrehilolo_h, double **solrelololo_h,
   double **solimhihihi_h, double **solimlohihi_h, 
   double **solimhilohi_h, double **solimlolohi_h, 
   double **solimhihilo_h, double **solimlohilo_h, 
   double **solimhilolo_h, double **solimlololo_h, 
   double **solrehihihi_d, double **solrelohihi_d, 
   double **solrehilohi_d, double **solrelolohi_d, 
   double **solrehihilo_d, double **solrelohilo_d, 
   double **solrehilolo_d, double **solrelololo_d, 
   double **solimhihihi_d, double **solimlohihi_d,
   double **solimhilohi_d, double **solimlolohi_d,
   double **solimhihilo_d, double **solimlohilo_d,
   double **solimhilolo_d, double **solimlololo_d,
   double **Qrehihihi_h, double **Qrelohihi_h,
   double **Qrehilohi_h, double **Qrelolohi_h,
   double **Qrehihilo_h, double **Qrelohilo_h,
   double **Qrehilolo_h, double **Qrelololo_h,
   double **Qimhihihi_h, double **Qimlohihi_h,
   double **Qimhilohi_h, double **Qimlolohi_h,
   double **Qimhihilo_h, double **Qimlohilo_h,
   double **Qimhilolo_h, double **Qimlololo_h,
   double **Qrehihihi_d, double **Qrelohihi_d,
   double **Qrehilohi_d, double **Qrelolohi_d,
   double **Qrehihilo_d, double **Qrelohilo_d,
   double **Qrehilolo_d, double **Qrelololo_d,
   double **Qimhihihi_d, double **Qimlohihi_d, 
   double **Qimhilohi_d, double **Qimlolohi_d, 
   double **Qimhihilo_d, double **Qimlohilo_d, 
   double **Qimhilolo_d, double **Qimlololo_d, 
   double **Rrehihihi_h, double **Rrelohihi_h,
   double **Rrehilohi_h, double **Rrelolohi_h,
   double **Rrehihilo_h, double **Rrelohilo_h,
   double **Rrehilolo_h, double **Rrelololo_h,
   double **Rimhihihi_h, double **Rimlohihi_h, 
   double **Rimhilohi_h, double **Rimlolohi_h, 
   double **Rimhihilo_h, double **Rimlohilo_h, 
   double **Rimhilolo_h, double **Rimlololo_h, 
   double **Rrehihihi_d, double **Rrelohihi_d,
   double **Rrehilohi_d, double **Rrelolohi_d,
   double **Rrehihilo_d, double **Rrelohilo_d,
   double **Rrehilolo_d, double **Rrelololo_d,
   double **Rimhihihi_d, double **Rimlohihi_d,
   double **Rimhilohi_d, double **Rimlolohi_d,
   double **Rimhihilo_d, double **Rimlohilo_d,
   double **Rimhilolo_d, double **Rimlololo_d,
   double **workmatrehihihi, double **workmatrelohihi,
   double **workmatrehilohi, double **workmatrelolohi,
   double **workmatrehihilo, double **workmatrelohilo,
   double **workmatrehilolo, double **workmatrelololo,
   double **workmatimhihihi, double **workmatimlohihi,
   double **workmatimhilohi, double **workmatimlolohi,
   double **workmatimhihilo, double **workmatimlohilo,
   double **workmatimhilolo, double **workmatimlololo,
   double *workvecrehihihi, double *workvecrelohihi,
   double *workvecrehilohi, double *workvecrelolohi,
   double *workvecrehihilo, double *workvecrelohilo,
   double *workvecrehilolo, double *workvecrelololo,
   double *workvecimhihihi, double *workvecimlohihi,
   double *workvecimhilohi, double *workvecimlolohi,
   double *workvecimhihilo, double *workvecimlohilo,
   double *workvecimhilolo, double *workvecimlololo,
   double **resvecrehihihi, double **resvecrelohihi,
   double **resvecrehilohi, double **resvecrelolohi,
   double **resvecrehihilo, double **resvecrelohilo,
   double **resvecrehilolo, double **resvecrelololo,
   double **resvecimhihihi, double **resvecimlohihi,
   double **resvecimhilohi, double **resvecimlolohi,
   double **resvecimhihilo, double **resvecimlohilo,
   double **resvecimhilolo, double **resvecimlololo,
   double *resmaxhihihi, double *resmaxlohihi,
   double *resmaxhilohi, double *resmaxlolohi,
   double *resmaxhihilo, double *resmaxlohilo,
   double *resmaxhilolo, double *resmaxlololo, int vrblvl, int mode );
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
 *   cffrehihihi are the highest doubles of the real parts of
 *             the coefficients of the monomials;
 *   cffrelohihi are the second highest doubles of the real parts of
 *             the coefficients of the monomials;
 *   cffrehilohi are the third highest doubles of the real parts of
 *             the coefficients of the monomials;
 *   cffrelolohi are the fourth highest doubles of the real parts of
 *             the coefficients of the monomials;
 *   cffrehihilo are the fourth lowest doubles of the real parts of
 *             the coefficients of the monomials;
 *   cffrelohilo are the third lowest doubles of the real parts of
 *             the coefficients of the monomials;
 *   cffrehilolo are the second lowest doubles of the real parts of
 *             the coefficients of the monomials;
 *   cffrelololo are the lowest doubles of the real parts of
 *             the coefficients of the monomials;
 *   cffimhihihi are the highest doubles of the imaginary parts of
 *             the coefficients of the monomials;
 *   cffimlohihi are the second highest doubles of the imaginary parts of
 *             the coefficients of the monomials;
 *   cffimhilohi are the third highest doubles of the imaginary parts of
 *             the coefficients of the monomials;
 *   cffimlolohi are the fourth highest doubles of the imaginary parts of
 *             the coefficients of the monomials;
 *   cffimhihilo are the fourth lowest doubles of the imaginary parts of
 *             the coefficients of the monomials;
 *   cffimlohilo are the third lowest doubles of the imaginary parts of
 *             the coefficients of the monomials;
 *   cffimhilolo are the second lowest doubles of the imaginary parts of
 *             the coefficients of the monomials;
 *   cffimlololo are the lowest doubles of the imaginary parts of
 *             the coefficients of the monomials;
 *   accrehihihi has space to accumulate one series of degree deg;
 *   accrelohihi has space to accumulate one series of degree deg;
 *   accrehilohi has space to accumulate one series of degree deg;
 *   accrelolohi has space to accumulate one series of degree deg;
 *   accrehihilo has space to accumulate one series of degree deg;
 *   accrelohilo has space to accumulate one series of degree deg;
 *   accrehilolo has space to accumulate one series of degree deg;
 *   accrelololo has space to accumulate one series of degree deg;
 *   accimhihihi has space to accumulate one series of degree deg;
 *   accimlohihi has space to accumulate one series of degree deg;
 *   accimhilohi has space to accumulate one series of degree deg;
 *   accimlolohi has space to accumulate one series of degree deg;
 *   accimhihilo has space to accumulate one series of degree deg;
 *   accimlohilo has space to accumulate one series of degree deg;
 *   accimhilolo has space to accumulate one series of degree deg;
 *   accimlololo has space to accumulate one series of degree deg;
 *   inputrehihihi_h are the highest doubles of the real parts
 *             of series of degree deg, for dim variables, computed on host;
 *   inputrelohihi_h are the second highest doubles of the real parts
 *             of series of degree deg, for dim variables, computed on host;
 *   inputrehilohi_h are the third highest doubles of the real parts
 *             of series of degree deg, for dim variables, computed on host;
 *   inputrelolohi_h are the fourth highest doubles of the real parts
 *             of series of degree deg, for dim variables, computed on host;
 *   inputrehihilo_h are the fourth lowest doubles of the real parts
 *             of series of degree deg, for dim variables, computed on host;
 *   inputrelohilo_h are the third lowest doubles of the real parts
 *             of series of degree deg, for dim variables, computed on host;
 *   inputrehilolo_h are the second lowest doubles of the real parts
 *             of series of degree deg, for dim variables, computed on host;
 *   inputrelolohi_h are the lowest doubles of the real parts
 *             of series of degree deg, for dim variables, computed on host;
 *   inputimhihihi_h are the highest doubles of the imaginary parts
 *             of series of degree deg for dim variables, computed on host;
 *   inputimlohihi_h are the second highest doubles of the imaginary parts
 *             of series of degree deg for dim variables, computed on host;
 *   inputimhilohi_h are the third highest doubles of the imaginary parts
 *             of series of degree deg for dim variables, computed on host;
 *   inputimlolohi_h are the fourth highest doubles of the imaginary parts
 *             of series of degree deg for dim variables, computed on host;
 *   inputimhihilo_h are the fourth lowest doubles of the imaginary parts
 *             of series of degree deg for dim variables, computed on host;
 *   inputimlohilo_h are the third lowest doubles of the imaginary parts
 *             of series of degree deg for dim variables, computed on host;
 *   inputimhilolo_h are the second lowest doubles of the imaginary parts
 *             of series of degree deg for dim variables, computed on host;
 *   inputimlolohi_h are the lowest doubles of the imaginary parts
 *             of series of degree deg for dim variables, computed on host;
 *   inputrehihihi_d has space for series computed on device;
 *   inputrelohihi_d has space for series computed on device;
 *   inputrehilohi_d has space for series computed on device;
 *   inputrelolohi_d has space for series computed on device;
 *   inputrehihilo_d has space for series computed on device;
 *   inputrelohilo_d has space for series computed on device;
 *   inputrehilolo_d has space for series computed on device;
 *   inputrelololo_d has space for series computed on device;
 *   inputimhihihi_d has space for series computed on device;
 *   inputimlohihi_d has space for series computed on device;
 *   inputimhilohi_d has space for series computed on device;
 *   inputimlolohi_d has space for series computed on device;
 *   inputimhihilo_d has space for series computed on device;
 *   inputimlohilo_d has space for series computed on device;
 *   inputimhilolo_d has space for series computed on device;
 *   inputimlololo_d has space for series computed on device;
 *   outputrehihihi_h has space for the evaluated and differentiated
 *             monomials, computed on the host;
 *   outputrelohihi_h has space for the evaluated and differentiated
 *             monomials, computed on the host;
 *   outputrehilohi_h has space for the evaluated and differentiated
 *             monomials, computed on the host;
 *   outputrelolohi_h has space for the evaluated and differentiated
 *             monomials, computed on the host;
 *   outputrehihilo_h has space for the evaluated and differentiated
 *             monomials, computed on the host;
 *   outputrelohilo_h has space for the evaluated and differentiated
 *             monomials, computed on the host;
 *   outputrehilolo_h has space for the evaluated and differentiated
 *             monomials, computed on the host;
 *   outputrelololo_h has space for the evaluated and differentiated
 *             monomials, computed on the host;
 *   outputimhihihi_h has space for the evaluated and differentiated
 *             monomials, computed on the host;
 *   outputimlohihi_h has space for the evaluated and differentiated
 *             monomials, computed on the host;
 *   outputimhilohi_h has space for the evaluated and differentiated
 *             monomials, computed on the host;
 *   outputimlolohi_h has space for the evaluated and differentiated
 *             monomials, computed on the host;
 *   outputimhihilo_h has space for the evaluated and differentiated
 *             monomials, computed on the host;
 *   outputimlohilo_h has space for the evaluated and differentiated
 *             monomials, computed on the host;
 *   outputimhilolo_h has space for the evaluated and differentiated
 *             monomials, computed on the host;
 *   outputimlololo_h has space for the evaluated and differentiated
 *             monomials, computed on the host;
 *   outputrehihihi_d has space for the evaluated and differentiated
 *             monomials, computed on the device;
 *   outputrelohihi_d has space for the evaluated and differentiated
 *             monomials, computed on the device;
 *   outputrehilohi_d has space for the evaluated and differentiated
 *             monomials, computed on the device;
 *   outputrelolohi_d has space for the evaluated and differentiated
 *             monomials, computed on the device;
 *   outputrehihilo_d has space for the evaluated and differentiated
 *             monomials, computed on the device;
 *   outputrelohilo_d has space for the evaluated and differentiated
 *             monomials, computed on the device;
 *   outputrehilolo_d has space for the evaluated and differentiated
 *             monomials, computed on the device;
 *   outputrelololo_d has space for the evaluated and differentiated
 *             monomials, computed on the device;
 *   outputimhihihi_d has space for the evaluated and differentiated
 *             monomials, computed on the device;
 *   outputimlohihi_d has space for the evaluated and differentiated
 *             monomials, computed on the device;
 *   outputimhilohi_d has space for the evaluated and differentiated
 *             monomials, computed on the device;
 *   outputimlolohi_d has space for the evaluated and differentiated
 *             monomials, computed on the device;
 *   outputimhihilo_d has space for the evaluated and differentiated
 *             monomials, computed on the device;
 *   outputimlohilo_d has space for the evaluated and differentiated
 *             monomials, computed on the device;
 *   outputimhilolo_d has space for the evaluated and differentiated
 *             monomials, computed on the device;
 *   outputimlololo_d has space for the evaluated and differentiated
 *             monomials, computed on the device;
 *   funvalrehihihi_h has space for the evaluated series on the host;
 *   funvalrelohihi_h has space for the evaluated series on the host;
 *   funvalrehilohi_h has space for the evaluated series on the host;
 *   funvalrelolohi_h has space for the evaluated series on the host;
 *   funvalrehihilo_h has space for the evaluated series on the host;
 *   funvalrelohilo_h has space for the evaluated series on the host;
 *   funvalrehilolo_h has space for the evaluated series on the host;
 *   funvalrelololo_h has space for the evaluated series on the host;
 *   funvalimhihihi_h has space for the evaluated series on the host;
 *   funvalimlohihi_h has space for the evaluated series on the host;
 *   funvalimhilohi_h has space for the evaluated series on the host;
 *   funvalimlolohi_h has space for the evaluated series on the host;
 *   funvalimhihilo_h has space for the evaluated series on the host;
 *   funvalimlohilo_h has space for the evaluated series on the host;
 *   funvalimhilolo_h has space for the evaluated series on the host;
 *   funvalimlololo_h has space for the evaluated series on the host;
 *   funvalrehihihi_d has space for the evaluated series on the device;
 *   funvalrelohihi_d has space for the evaluated series on the device;
 *   funvalrehilohi_d has space for the evaluated series on the device;
 *   funvalrelolohi_d has space for the evaluated series on the device;
 *   funvalrehihilo_d has space for the evaluated series on the device;
 *   funvalrelohilo_d has space for the evaluated series on the device;
 *   funvalrehilolo_d has space for the evaluated series on the device;
 *   funvalrelololo_d has space for the evaluated series on the device;
 *   funvalimhihihi_d has space for the evaluated series on the device;
 *   funvalimlohihi_d has space for the evaluated series on the device;
 *   funvalimhilohi_d has space for the evaluated series on the device;
 *   funvalimlolohi_d has space for the evaluated series on the device;
 *   funvalimhihilo_d has space for the evaluated series on the device;
 *   funvalimlohilo_d has space for the evaluated series on the device;
 *   funvalimhilolo_d has space for the evaluated series on the device;
 *   funvalimlololo_d has space for the evaluated series on the device;
 *   jacvalrehihihi_h has space for deg+1 matrices of dimension dim on host;
 *   jacvalrelohihi_h has space for deg+1 matrices of dimension dim on host;
 *   jacvalrehilohi_h has space for deg+1 matrices of dimension dim on host;
 *   jacvalrelolohi_h has space for deg+1 matrices of dimension dim on host;
 *   jacvalrehihilo_h has space for deg+1 matrices of dimension dim on host;
 *   jacvalrelohilo_h has space for deg+1 matrices of dimension dim on host;
 *   jacvalrehilolo_h has space for deg+1 matrices of dimension dim on host;
 *   jacvalrelololo_h has space for deg+1 matrices of dimension dim on host;
 *   jacvalimhihihi_h has space for deg+1 matrices of dimension dim on host;
 *   jacvalimlohihi_h has space for deg+1 matrices of dimension dim on host;
 *   jacvalimhilohi_h has space for deg+1 matrices of dimension dim on host;
 *   jacvalimlolohi_h has space for deg+1 matrices of dimension dim on host;
 *   jacvalimhihilo_h has space for deg+1 matrices of dimension dim on host;
 *   jacvalimlohilo_h has space for deg+1 matrices of dimension dim on host;
 *   jacvalimhilolo_h has space for deg+1 matrices of dimension dim on host;
 *   jacvalimlololo_h has space for deg+1 matrices of dimension dim on host;
 *   jacvalrehihihi_d has space for deg+1 matrices of dimension dim on device;
 *   jacvalrelohihi_d has space for deg+1 matrices of dimension dim on device;
 *   jacvalrehilohi_d has space for deg+1 matrices of dimension dim on device;
 *   jacvalrelolohi_d has space for deg+1 matrices of dimension dim on device;
 *   jacvalrehihilo_d has space for deg+1 matrices of dimension dim on device;
 *   jacvalrelohilo_d has space for deg+1 matrices of dimension dim on device;
 *   jacvalrehilolo_d has space for deg+1 matrices of dimension dim on device;
 *   jacvalrelololo_d has space for deg+1 matrices of dimension dim on device;
 *   jacvalimhihihi_d has space for deg+1 matrices of dimension dim on device;
 *   jacvalimlohihi_d has space for deg+1 matrices of dimension dim on device;
 *   jacvalimhilohi_d has space for deg+1 matrices of dimension dim on device;
 *   jacvalimlolohi_d has space for deg+1 matrices of dimension dim on device;
 *   jacvalimhihilo_d has space for deg+1 matrices of dimension dim on device;
 *   jacvalimlohilo_d has space for deg+1 matrices of dimension dim on device;
 *   jacvalimhilolo_d has space for deg+1 matrices of dimension dim on device;
 *   jacvalimlololo_d has space for deg+1 matrices of dimension dim on device;
 *   rhsrehihihi_h has space for deg+1 vectors of dimension dim on host;
 *   rhsrelohihi_h has space for deg+1 vectors of dimension dim on host;
 *   rhsrehilohi_h has space for deg+1 vectors of dimension dim on host;
 *   rhsrelolohi_h has space for deg+1 vectors of dimension dim on host;
 *   rhsrehihilo_h has space for deg+1 vectors of dimension dim on host;
 *   rhsrelohilo_h has space for deg+1 vectors of dimension dim on host;
 *   rhsrehilolo_h has space for deg+1 vectors of dimension dim on host;
 *   rhsrelololo_h has space for deg+1 vectors of dimension dim on host;
 *   rhsimhihihi_h has space for deg+1 vectors of dimension dim on host;
 *   rhsimlohihi_h has space for deg+1 vectors of dimension dim on host;
 *   rhsimhilohi_h has space for deg+1 vectors of dimension dim on host;
 *   rhsimlolohi_h has space for deg+1 vectors of dimension dim on host;
 *   rhsimhihilo_h has space for deg+1 vectors of dimension dim on host;
 *   rhsimlohilo_h has space for deg+1 vectors of dimension dim on host;
 *   rhsimhilolo_h has space for deg+1 vectors of dimension dim on host;
 *   rhsimlololo_h has space for deg+1 vectors of dimension dim on host;
 *   rhsrehihihi_d has space for deg+1 vectors of dimension dim on device;
 *   rhsrelohihi_d has space for deg+1 vectors of dimension dim on device;
 *   rhsrehilohi_d has space for deg+1 vectors of dimension dim on device;
 *   rhsrelolohi_d has space for deg+1 vectors of dimension dim on device;
 *   rhsrehihilo_d has space for deg+1 vectors of dimension dim on device;
 *   rhsrelohilo_d has space for deg+1 vectors of dimension dim on device;
 *   rhsrehilolo_d has space for deg+1 vectors of dimension dim on device;
 *   rhsrelololo_d has space for deg+1 vectors of dimension dim on device;
 *   rhsimhihihi_d has space for deg+1 vectors of dimension dim on device;
 *   rhsimlohihi_d has space for deg+1 vectors of dimension dim on device;
 *   rhsimhilohi_d has space for deg+1 vectors of dimension dim on device;
 *   rhsimlolohi_d has space for deg+1 vectors of dimension dim on device;
 *   rhsimhihilo_d has space for deg+1 vectors of dimension dim on device;
 *   rhsimlohilo_d has space for deg+1 vectors of dimension dim on device;
 *   rhsimhilolo_d has space for deg+1 vectors of dimension dim on device;
 *   rhsimlololo_d has space for deg+1 vectors of dimension dim on device;
 *   urhsrehihihi_h has space for updated right hand side on the host;
 *   urhsrelohihi_h has space for updated right hand side on the host;
 *   urhsrehilohi_h has space for updated right hand side on the host;
 *   urhsrelolohi_h has space for updated right hand side on the host;
 *   urhsrehihilo_h has space for updated right hand side on the host;
 *   urhsrelohilo_h has space for updated right hand side on the host;
 *   urhsrehilolo_h has space for updated right hand side on the host;
 *   urhsrelololo_h has space for updated right hand side on the host;
 *   urhsimhihihi_h has space for updated right hand side on the host;
 *   urhsimlohihi_h has space for updated right hand side on the host;
 *   urhsimhilohi_h has space for updated right hand side on the host;
 *   urhsimlolohi_h has space for updated right hand side on the host;
 *   urhsimhihilo_h has space for updated right hand side on the host;
 *   urhsimlohilo_h has space for updated right hand side on the host;
 *   urhsimhilolo_h has space for updated right hand side on the host;
 *   urhsimlololo_h has space for updated right hand side on the host;
 *   urhsrehihihi_d has space for updated right hand side on the device; 
 *   urhsrelohihi_d has space for updated right hand side on the device; 
 *   urhsrehilohi_d has space for updated right hand side on the device; 
 *   urhsrelolohi_d has space for updated right hand side on the device; 
 *   urhsrehihilo_d has space for updated right hand side on the device; 
 *   urhsrelohilo_d has space for updated right hand side on the device; 
 *   urhsrehilolo_d has space for updated right hand side on the device; 
 *   urhsrelololo_d has space for updated right hand side on the device; 
 *   urhsimhihihi_d has space for updated right hand side on the device; 
 *   urhsimlohihi_d has space for updated right hand side on the device; 
 *   urhsimhilohi_d has space for updated right hand side on the device; 
 *   urhsimlolohi_d has space for updated right hand side on the device; 
 *   urhsimhihilo_d has space for updated right hand side on the device; 
 *   urhsimlohilo_d has space for updated right hand side on the device; 
 *   urhsimhilolo_d has space for updated right hand side on the device; 
 *   urhsimlololo_d has space for updated right hand side on the device; 
 *   solrehihihi_h has space for deg+1 vectors of dimension dim;
 *   solrelohihi_h has space for deg+1 vectors of dimension dim;
 *   solrehilohi_h has space for deg+1 vectors of dimension dim;
 *   solrelolohi_h has space for deg+1 vectors of dimension dim;
 *   solrehihilo_h has space for deg+1 vectors of dimension dim;
 *   solrelohilo_h has space for deg+1 vectors of dimension dim;
 *   solrehilolo_h has space for deg+1 vectors of dimension dim;
 *   solrelololo_h has space for deg+1 vectors of dimension dim;
 *   solimhihihi_h has space for deg+1 vectors of dimension dim;
 *   solimlohihi_h has space for deg+1 vectors of dimension dim;
 *   solimhilohi_h has space for deg+1 vectors of dimension dim;
 *   solimlolohi_h has space for deg+1 vectors of dimension dim;
 *   solimhihilo_h has space for deg+1 vectors of dimension dim;
 *   solimlohilo_h has space for deg+1 vectors of dimension dim;
 *   solimhilolo_h has space for deg+1 vectors of dimension dim;
 *   solimlololo_h has space for deg+1 vectors of dimension dim;
 *   solrehihihi_d has space for deg+1 vectors of dimension dim;
 *   solrelohihi_d has space for deg+1 vectors of dimension dim;
 *   solrehilohi_d has space for deg+1 vectors of dimension dim;
 *   solrelolohi_d has space for deg+1 vectors of dimension dim;
 *   solrehihilo_d has space for deg+1 vectors of dimension dim;
 *   solrelohilo_d has space for deg+1 vectors of dimension dim;
 *   solrehilolo_d has space for deg+1 vectors of dimension dim;
 *   solrelololo_d has space for deg+1 vectors of dimension dim;
 *   solimhihihi_d has space for deg+1 vectors of dimension dim;
 *   solimlohihi_d has space for deg+1 vectors of dimension dim;
 *   solimhilohi_d has space for deg+1 vectors of dimension dim;
 *   solimlolohi_d has space for deg+1 vectors of dimension dim;
 *   solimhihilo_d has space for deg+1 vectors of dimension dim;
 *   solimlohilo_d has space for deg+1 vectors of dimension dim;
 *   solimhilolo_d has space for deg+1 vectors of dimension dim;
 *   solimlololo_d has space for deg+1 vectors of dimension dim;
 *   Qrehihihi_h has space for the Q computed by the host;
 *   Qrelohihi_h has space for the Q computed by the host;
 *   Qrehilohi_h has space for the Q computed by the host;
 *   Qrelolohi_h has space for the Q computed by the host;
 *   Qrehihilo_h has space for the Q computed by the host;
 *   Qrelohilo_h has space for the Q computed by the host;
 *   Qrehilolo_h has space for the Q computed by the host;
 *   Qrelololo_h has space for the Q computed by the host;
 *   Qimhihihi_h has space for the Q computed by the host;
 *   Qimlohihi_h has space for the Q computed by the host;
 *   Qimhilohi_h has space for the Q computed by the host;
 *   Qimlolohi_h has space for the Q computed by the host;
 *   Qimhihilo_h has space for the Q computed by the host;
 *   Qimlohilo_h has space for the Q computed by the host;
 *   Qimhilolo_h has space for the Q computed by the host;
 *   Qimlololo_h has space for the Q computed by the host;
 *   Qrehihihi_d has space for the Q computed by the device;
 *   Qrelohihi_d has space for the Q computed by the device;
 *   Qrehilohi_d has space for the Q computed by the device;
 *   Qrelolohi_d has space for the Q computed by the device;
 *   Qrehihilo_d has space for the Q computed by the device;
 *   Qrelohilo_d has space for the Q computed by the device;
 *   Qrehilolo_d has space for the Q computed by the device;
 *   Qrelololo_d has space for the Q computed by the device;
 *   Qimhihihi_d has space for the Q computed by the device;
 *   Qimlohihi_d has space for the Q computed by the device;
 *   Qimhilohi_d has space for the Q computed by the device;
 *   Qimlolohi_d has space for the Q computed by the device;
 *   Qimhihilo_d has space for the Q computed by the device;
 *   Qimlohilo_d has space for the Q computed by the device;
 *   Qimhilolo_d has space for the Q computed by the device;
 *   Qimlololo_d has space for the Q computed by the device;
 *   Rrehihihi_h has space for the R computed by the host;
 *   Rrelohihi_h has space for the R computed by the host;
 *   Rrehilohi_h has space for the R computed by the host;
 *   Rrelolohi_h has space for the R computed by the host;
 *   Rrehihilo_h has space for the R computed by the host;
 *   Rrelohilo_h has space for the R computed by the host;
 *   Rrehilolo_h has space for the R computed by the host;
 *   Rrelololo_h has space for the R computed by the host;
 *   Rimhihihi_h has space for the R computed by the host;
 *   Rimlohihi_h has space for the R computed by the host;
 *   Rimhilohi_h has space for the R computed by the host;
 *   Rimlolohi_h has space for the R computed by the host;
 *   Rimhihilo_h has space for the R computed by the host;
 *   Rimlohilo_h has space for the R computed by the host;
 *   Rimhilolo_h has space for the R computed by the host;
 *   Rimlololo_h has space for the R computed by the host;
 *   Rrehihihi_d has space for the R computed by the device;
 *   Rrelohihi_d has space for the R computed by the device;
 *   Rrehilohi_d has space for the R computed by the device;
 *   Rrelolohi_d has space for the R computed by the device;
 *   Rrehihilo_d has space for the R computed by the device;
 *   Rrelohilo_d has space for the R computed by the device;
 *   Rrehilolo_d has space for the R computed by the device;
 *   Rrelololo_d has space for the R computed by the device;
 *   Rimhihihi_d has space for the R computed by the device;
 *   Rimlohihi_d has space for the R computed by the device;
 *   Rimhilohi_d has space for the R computed by the device;
 *   Rimlolohi_d has space for the R computed by the device;
 *   Rimhihilo_d has space for the R computed by the device;
 *   Rimlohilo_d has space for the R computed by the device;
 *   Rimhilolo_d has space for the R computed by the device;
 *   Rimlololo_d has space for the R computed by the device;
 *   wrkmatrehihihi is work space allocated for a matrix of dimension dim;
 *   wrkmatrelohihi is work space allocated for a matrix of dimension dim;
 *   wrkmatrehilohi is work space allocated for a matrix of dimension dim;
 *   wrkmatrelolohi is work space allocated for a matrix of dimension dim;
 *   wrkmatrehihilo is work space allocated for a matrix of dimension dim;
 *   wrkmatrelohilo is work space allocated for a matrix of dimension dim;
 *   wrkmatrehilolo is work space allocated for a matrix of dimension dim;
 *   wrkmatrelololo is work space allocated for a matrix of dimension dim;
 *   wrkmatimhihihi is work space allocated for a matrix of dimension dim;
 *   wrkmatimlohihi is work space allocated for a matrix of dimension dim;
 *   wrkmatimhilohi is work space allocated for a matrix of dimension dim;
 *   wrkmatimlolohi is work space allocated for a matrix of dimension dim;
 *   wrkmatimhihilo is work space allocated for a matrix of dimension dim;
 *   wrkmatimlohilo is work space allocated for a matrix of dimension dim;
 *   wrkmatimhilolo is work space allocated for a matrix of dimension dim;
 *   wrkmatimlololo is work space allocated for a matrix of dimension dim;
 *   wrkvecrehihihi is work space allocated for a vector of dimension dim;
 *   wrkvecrelohihi is work space allocated for a vector of dimension dim;
 *   wrkvecrehilohi is work space allocated for a vector of dimension dim;
 *   wrkvecrelolohi is work space allocated for a vector of dimension dim;
 *   wrkvecrehihilo is work space allocated for a vector of dimension dim;
 *   wrkvecrelohilo is work space allocated for a vector of dimension dim;
 *   wrkvecrehilolo is work space allocated for a vector of dimension dim;
 *   wrkvecrelololo is work space allocated for a vector of dimension dim;
 *   wrkvecimhihihi is work space allocated for a vector of dimension dim;
 *   wrkvecimlohihi is work space allocated for a vector of dimension dim;
 *   wrkvecimhilohi is work space allocated for a vector of dimension dim;
 *   wrkvecimlolohi is work space allocated for a vector of dimension dim;
 *   wrkvecimhihilo is work space allocated for a vector of dimension dim;
 *   wrkvecimlohilo is work space allocated for a vector of dimension dim;
 *   wrkvecimhilolo is work space allocated for a vector of dimension dim;
 *   wrkvecimlololo is work space allocated for a vector of dimension dim;
 *   resvecrehihihi has space for deg+1 vectors of dimension dim;
 *   resvecrelohihi has space for deg+1 vectors of dimension dim;
 *   resvecrehilohi has space for deg+1 vectors of dimension dim;
 *   resvecrelolohi has space for deg+1 vectors of dimension dim;
 *   resvecrehihilo has space for deg+1 vectors of dimension dim;
 *   resvecrelohilo has space for deg+1 vectors of dimension dim;
 *   resvecrehilolo has space for deg+1 vectors of dimension dim;
 *   resvecrelololo has space for deg+1 vectors of dimension dim;
 *   resvecimhihihi has space for deg+1 vectors of dimension dim;
 *   resvecimlohihi has space for deg+1 vectors of dimension dim;
 *   resvecimhilohi has space for deg+1 vectors of dimension dim;
 *   resvecimlolohi has space for deg+1 vectors of dimension dim;
 *   resvecimhihilo has space for deg+1 vectors of dimension dim;
 *   resvecimlohilo has space for deg+1 vectors of dimension dim;
 *   resvecimhilolo has space for deg+1 vectors of dimension dim;
 *   resvecimlololo has space for deg+1 vectors of dimension dim;
 *   vrblvl    is the verbose level;
 *   mode      execution mode, 0 (GPU only), 1 (CPU only) or 2 (GPU+CPU).
 *
 * ON RETURN :
 *   inputrehihihi_h has the highest doubles of the real parts
 *             of series, computed on host;
 *   inputrelohihi_h has the second highest doubles of the real parts
 *             of series, computed on host;
 *   inputrehilohi_h has the third highest doubles of the real parts
 *             of series, computed on host;
 *   inputrelolohi_h has the fourth highest doubles of the real parts
 *             of series, computed on host;
 *   inputrehihilo_h has the fourth lowest doubles of the real parts
 *             of series, computed on host;
 *   inputrelohilo_h has the third lowest doubles of the real parts
 *             of series, computed on host;
 *   inputrehilolo_h has the second lowest doubles of the real parts
 *             of series, computed on host;
 *   inputrelololo_h has the lowest doubles of the real parts
 *             of series, computed on host;
 *   inputimhihihi_h has the highest doubles of the imaginary parts
 *             of series, computed on host;
 *   inputimlohihi_h has the second highest doubles of the imaginary parts
 *             of series, computed on host;
 *   inputimhilohi_h has the third highest doubles of the imaginary parts
 *             of series, computed on host;
 *   inputimlolohi_h has the fourth highest doubles of the imaginary parts
 *             of series, computed on host;
 *   inputimhihilo_h has the fourth lowest doubles of the imaginary parts
 *             of series, computed on host;
 *   inputimlohilo_h has the third lowest doubles of the imaginary parts
 *             of series, computed on host;
 *   inputimhilolo_h has the second lowest doubles of the imaginary parts
 *             of series, computed on host;
 *   inputimlololo_h has the lowest doubles of the imaginary parts
 *             of series, computed on host;
 *   inputrehihihi_d has the highest doubles of the real parts
 *             of series, computed on device;
 *   inputrelohihi_d has the second highest doubles of the real parts
 *             of series, computed on device;
 *   inputrehilohi_d has the third highest doubles of the real parts
 *             of series, computed on device;
 *   inputrelolohi_d has the fourth highest doubles of the real parts
 *             of series, computed on device;
 *   inputrehihilo_d has the fourth lowest doubles of the real parts
 *             of series, computed on device;
 *   inputrelohilo_d has the third lowest doubles of the real parts
 *             of series, computed on device;
 *   inputrehilolo_d has the second lowest doubles of the real parts
 *             of series, computed on device;
 *   inputrelololo_d has the lowest doubles of the real parts
 *             of series, computed on device;
 *   inputimhihihi_d has the highest doubles of the imaginary parts
 *             of series, computed on device;
 *   inputimlohihi_d has the second highest doubles of the imaginary parts
 *             of series, computed on device;
 *   inputimhilohi_d has the third highest doubles of the imaginary parts
 *             of series, computed on device;
 *   inputimlolohi_d has the fourth highest doubles of the imaginary parts
 *             of series, computed on device;
 *   inputimhihilo_d has the fourth lowest doubles of the imaginary parts
 *             of series, computed on device;
 *   inputimlohilo_d has the third lowest doubles of the imaginary parts
 *             of series, computed on device;
 *   inputimhilolo_d has the second lowest doubles of the imaginary parts
 *             of series, computed on device;
 *   inputimlololo_d has the lowest doubles of the imaginary parts
 *             of series, computed on device;
 *   outputrehihihi_h has the highest doubles of the real parts
 *             of the output, computed on host;
 *   outputrelohihi_h has the second highest doubles of the real parts
 *             of the output, computed on host;
 *   outputrehilohi_h has the third highest doubles of the real parts
 *             of the output, computed on host;
 *   outputrelolohi_h has the fourth highest doubles of the real parts
 *             of the output, computed on host;
 *   outputrehihilo_h has the fourth lowest doubles of the real parts
 *             of the output, computed on host;
 *   outputrelohilo_h has the third lowest doubles of the real parts
 *             of the output, computed on host;
 *   outputrehilolo_h has the second lowest doubles of the real parts
 *             of the output, computed on host;
 *   outputrelololo_h has the lowest doubles of the real parts
 *             of the output, computed on host;
 *   outputimhihihi_h has the highest doubles of the imaginary parts
 *             of the output, computed on host;
 *   outputimlohihi_h has the second highest doubles of the imaginary parts
 *             of the output, computed on host;
 *   outputimhilohi_h has the third highest doubles of the imaginary parts
 *             of the output, computed on host;
 *   outputimlolohi_h has the fourth highest doubles of the imaginary parts
 *             of the output, computed on host;
 *   outputimhihilo_h has the fourth lowest doubles of the imaginary parts
 *             of the output, computed on host;
 *   outputimlohilo_h has the third lowest doubles of the imaginary parts
 *             of the output, computed on host;
 *   outputimhilolo_h has the second lowest doubles of the imaginary parts
 *             of the output, computed on host;
 *   outputimlololo_h has the lowest doubles of the imaginary parts
 *             of the output, computed on host;
 *   outputrehihihi_d has the highest doubles of the real parts
 *             of the output, computed on device;
 *   outputrelohihi_d has the second highest doubles of the real parts
 *             of the output, computed on device;
 *   outputrehilohi_d has the third highest doubles of the real parts
 *             of the output, computed on device;
 *   outputrelolohi_d has the fourth highest doubles of the real parts
 *             of the output, computed on device;
 *   outputrehihilo_d has the fourth lowest doubles of the real parts
 *             of the output, computed on device;
 *   outputrelohilo_d has the third lowest doubles of the real parts
 *             of the output, computed on device;
 *   outputrehihilo_d has the second lowest doubles of the real parts
 *             of the output, computed on device;
 *   outputrelololo_d has the lowest doubles of the real parts
 *             of the output, computed on device;
 *   outputimhihihi_d has the highest doubles of the imaginary parts
 *             of the output, computed on device;
 *   outputimlohihi_d has the second highest doubles of the imaginary parts
 *             of the output, computed on device;
 *   outputimhilohi_d has the third highest doubles of the imaginary parts
 *             of the output, computed on device;
 *   outputimlolohi_d has the fourth highest doubles of the imaginary parts
 *             of the output, computed on device;
 *   outputimhihilo_d has the fourth lowest doubles of the imaginary parts
 *             of the output, computed on device;
 *   outputimlohilo_d has the third lowest doubles of the imaginary parts
 *             of the output, computed on device;
 *   outputimhilolo_d has the second lowest doubles of the imaginary parts
 *             of the output, computed on device;
 *   outputimlololo_d has the lowest doubles of the imaginary parts
 *             of the output, computed on device;
 *   funvalrehihihi_h is outputrehihi[i][dim], on host;
 *   funvalrelohihi_h is outputrelohi[i][dim], on host;
 *   funvalrehilohi_h is outputrehilo[i][dim], on host;
 *   funvalrelolohi_h is outputrelolo[i][dim], on host;
 *   funvalrehihilo_h is outputrehihi[i][dim], on host;
 *   funvalrelohilo_h is outputrelohi[i][dim], on host;
 *   funvalrehilolo_h is outputrehilo[i][dim], on host;
 *   funvalrelololo_h is outputrelolo[i][dim], on host;
 *   funvalimhihihi_h is outputimhihi[i][dim], on host;
 *   funvalimlohihi_h is outputimlohi[i][dim], on host;
 *   funvalimhilohi_h is outputimhilo[i][dim], on host;
 *   funvalimlolohi_h is outputimlolo[i][dim], on host;
 *   funvalimhihilo_h is outputimhihi[i][dim], on host;
 *   funvalimlohilo_h is outputimlohi[i][dim], on host;
 *   funvalimhilolo_h is outputimhilo[i][dim], on host;
 *   funvalimlololo_h is outputimlolo[i][dim], on host;
 *   funvalrehihihi_d is outputrehihi[i][dim], on device;
 *   funvalrelohihi_d is outputrelohi[i][dim], on device;
 *   funvalrehilohi_d is outputrehilo[i][dim], on device;
 *   funvalrelolohi_d is outputrelolo[i][dim], on device;
 *   funvalrehihilo_d is outputrehihi[i][dim], on device;
 *   funvalrelohilo_d is outputrelohi[i][dim], on device;
 *   funvalrehilolo_d is outputrehilo[i][dim], on device;
 *   funvalrelololo_d is outputrelolo[i][dim], on device;
 *   funvalimhihihi_d is outputimhihi[i][dim], on device;
 *   funvalimlohihi_d is outputimlohi[i][dim], on device;
 *   funvalimhilohi_d is outputimhilo[i][dim], on device;
 *   funvalimlolohi_d is outputimlolo[i][dim], on device;
 *   funvalimhihilo_d is outputimhihi[i][dim], on device;
 *   funvalimlohilo_d is outputimlohi[i][dim], on device;
 *   funvalimhilolo_d is outputimhilo[i][dim], on device;
 *   funvalimlololo_d is outputimlolo[i][dim], on device;
 *   jacvalrehihihi_h are the highest doubles of the real parts
 *             of a matrix series, computed on host;
 *   jacvalrelohihi_h are the second highest doubles of the real parts
 *             of a matrix series, computed on host;
 *   jacvalrehilohi_h are the third highest doubles of the real parts
 *             of a matrix series, computed on host;
 *   jacvalrelolohi_h are the fourth highest doubles of the real parts
 *             of a matrix series, computed on host;
 *   jacvalrehihilo_h are the fourth lowest doubles of the real parts
 *             of a matrix series, computed on host;
 *   jacvalrelohilo_h are the third lowest doubles of the real parts
 *             of a matrix series, computed on host;
 *   jacvalrehilolo_h are the second lowest doubles of the real parts
 *             of a matrix series, computed on host;
 *   jacvalrelololo_h are the lowest doubles of the real parts
 *             of a matrix series, computed on host;
 *   jacvalimhihihi_h are the highest doubles of the real parts
 *             of a matrix series, computed on host;
 *   jacvalimlohihi_h are the second highest doubles of the real parts
 *             of a matrix series, computed on host;
 *   jacvalimhilohi_h are the third highest doubles of the real parts
 *             of a matrix series, computed on host;
 *   jacvalimlolohi_h are the fourth highest doubles of the real parts
 *             of a matrix series, computed on host;
 *   jacvalimhihilo_h are the fourth lowest doubles of the real parts
 *             of a matrix series, computed on host;
 *   jacvalimlohilo_h are the third lowest doubles of the real parts
 *             of a matrix series, computed on host;
 *   jacvalimhilolo_h are the second lowest doubles of the real parts
 *             of a matrix series, computed on host;
 *   jacvalimlololo_h are the lowest doubles of the real parts
 *             of a matrix series, computed on host;
 *   jacvalrehihihi_d are the highest doubles of the real parts
 *             of a matrix series, computed on device;
 *   jacvalrelohihi_d are the second highest doubles of the real parts
 *             of a matrix series, computed on device;
 *   jacvalrehilohi_d are the third highest doubles of the real parts
 *             of a matrix series, computed on device;
 *   jacvalrelolohi_d are the fourth highest doubles of the real parts
 *             of a matrix series, computed on device;
 *   jacvalrehihilo_d are the fourth lowest doubles of the real parts
 *             of a matrix series, computed on device;
 *   jacvalrelohilo_d are the third lowest doubles of the real parts
 *             of a matrix series, computed on device;
 *   jacvalrehilolo_d are the second lowest doubles of the real parts
 *             of a matrix series, computed on device;
 *   jacvalrelololo_d are the lowest doubles of the real parts
 *             of a matrix series, computed on device;
 *   jacvalimhihihi_d are the highest doubles of the real parts
 *             of a matrix series, computed on device;
 *   jacvalimlohihi_d are the second highest doubles of the real parts
 *             of a matrix series, computed on device;
 *   jacvalimhilohi_d are the third highest doubles of the real parts
 *             of a matrix series, computed on device;
 *   jacvalimlolohi_d are the fourth highest doubles of the real parts
 *             of a matrix series, computed on device;
 *   jacvalimhihilo_d are the fourth lowest doubles of the real parts
 *             of a matrix series, computed on device;
 *   jacvalimlohilo_d are the third lowest doubles of the real parts
 *             of a matrix series, computed on device;
 *   jacvalimhilolo_d are the second lowest doubles of the real parts
 *             of a matrix series, computed on device;
 *   jacvalimlololo_d are the lowest doubles of the real parts
 *             of a matrix series, computed on device;
 *   rhsrehihihi_h are highest doubles of real parts of the linearized
 *             right hand side, subtracted by 1 and added by t, on host;
 *   rhsrelohihi_h are second highest doubles of real parts of the linearized
 *             right hand side, subtracted by 1 and added by t, on host;
 *   rhsrehilohi_h are third highest doubles of real parts of the linearized
 *             right hand side, subtracted by 1 and added by t, on host;
 *   rhsrelolohi_h are fourth highest doubles of real parts of the linearized
 *             right hand side, subtracted by 1 and added by t, on host;
 *   rhsrehihilo_h are fourth lowest doubles of real parts of the linearized
 *             right hand side, subtracted by 1 and added by t, on host;
 *   rhsrelohilo_h are third lowest doubles of real parts of the linearized
 *             right hand side, subtracted by 1 and added by t, on host;
 *   rhsrehilolo_h are second lowest doubles of real parts of the linearized
 *             right hand side, subtracted by 1 and added by t, on host;
 *   rhsrelololo_h are lowest doubles of real parts of the linearized
 *             right hand side, subtracted by 1 and added by t, on host;
 *   rhsimhihihi_h are highest doubles of imag parts of the linearized
 *             right hand side, subtracted by 1 and added by t, on host;
 *   rhsimlohihi_h are second highest doubles of imag parts of the linearized
 *             right hand side, subtracted by 1 and added by t, on host;
 *   rhsimhilohi_h are third highest doubles of imag parts of the linearized
 *             right hand side, subtracted by 1 and added by t, on host;
 *   rhsimlohihi_h are fourth highest doubles of imag parts of the linearized
 *             right hand side, subtracted by 1 and added by t, on host;
 *   rhsimhihilo_h are fourth lowest doubles of imag parts of the linearized
 *             right hand side, subtracted by 1 and added by t, on host;
 *   rhsimlohilo_h are third lowest doubles of imag parts of the linearized
 *             right hand side, subtracted by 1 and added by t, on host;
 *   rhsimhilolo_h are second lowest doubles of imag parts of the linearized
 *             right hand side, subtracted by 1 and added by t, on host;
 *   rhsimlololo_h are lowest doubles of imag parts of the linearized
 *             right hand side, subtracted by 1 and added by t, on host;
 *   rhsrehihihi_d are highest doubles of real parts of the linearized 
 *             right hand side, subtracted by 1 and added by t, on device;
 *   rhsrelohihi_d are second highest doubles of real parts of the linearized 
 *             right hand side, subtracted by 1 and added by t, on device;
 *   rhsrehilohi_d are third highest doubles of real parts of the linearized 
 *             right hand side, subtracted by 1 and added by t, on device;
 *   rhsrelolohi_d are fourth highest doubles of real parts of the linearized 
 *             right hand side, subtracted by 1 and added by t, on device;
 *   rhsrehihilo_d are fourth lowest doubles of real parts of the linearized,
 *             right hand side, subtracted by 1 and added by t, on device;
 *   rhsrelohilo_d are third lowest doubles of real parts of the linearized,
 *             right hand side, subtracted by 1 and added by t, on device;
 *   rhsrehilolo_d are second lowest doubles of real parts of the linearized,
 *             right hand side, subtracted by 1 and added by t, on device;
 *   rhsrelololo_d are lowest doubles of real parts of the linearized,
 *             right hand side, subtracted by 1 and added by t, on device;
 *   rhsimhihihi_d are highest doubles of imag parts of the linearized
 *             right hand side, subtracted by 1 and added by t, on device;
 *   rhsimlohihi_d are second highest doubles of imag parts of the linearized
 *             right hand side, subtracted by 1 and added by t, on device;
 *   rhsimhilohi_d are third highest doubles of imag parts of the linearized
 *             right hand side, subtracted by 1 and added by t, on device;
 *   rhsimlolohi_d are fourth highest doubles of imag parts of the linearized
 *             right hand side, subtracted by 1 and added by t, on device;
 *   rhsimhihilo_d are fourth lowest doubles of imag parts of the linearized
 *             right hand side, subtracted by 1 and added by t, on device;
 *   rhsimlohilo_d are third lowest doubles of imag parts of the linearized
 *             right hand side, subtracted by 1 and added by t, on device;
 *   rhsimhilolo_d are second lowest doubles of imag parts of the linearized
 *             right hand side, subtracted by 1 and added by t, on device;
 *   rhsimlololo_d are lowest doubles of imag parts of the linearized
 *             right hand side, subtracted by 1 and added by t, on device;
 *   urhsrehihihi_h are the highest doubles of the real parts
 *             of the right hand side updated by the host;
 *   urhsrelohihi_h are the second highest doubles of the real parts
 *             of the right hand side updated by the host;
 *   urhsrehilohi_h are the third highest doubles of the real parts
 *             of the right hand side updated by the host;
 *   urhsrelolohi_h are the fourth highest doubles of the real parts
 *             of the right hand side updated by the host;
 *   urhsrehihilo_h are the fourth lowest doubles of the real parts
 *             of the right hand side updated by the host;
 *   urhsrelohilo_h are the third lowest doubles of the real parts
 *             of the right hand side updated by the host;
 *   urhsrehilolo_h are the second lowest doubles of the real parts
 *             of the right hand side updated by the host;
 *   urhsrelololo_h are the lowest doubles of the real parts
 *             of the right hand side updated by the host;
 *   urhsimhihihi_h are the highest doubles of the imaginary parts
 *             of the right hand side updated by the host;
 *   urhsimlohihi_h are the second highest doubles of the imaginary parts
 *             of the right hand side updated by the host;
 *   urhsimhilohi_h are the third highest doubles of the imaginary parts
 *             of the right hand side updated by the host;
 *   urhsimlolohi_h are the fourth highest doubles of the imaginary parts
 *             of the right hand side updated by the host;
 *   urhsimhihilo_h are the fourth lowest doubles of the imaginary parts
 *             of the right hand side updated by the host;
 *   urhsimlohilo_h are the third lowest doubles of the imaginary parts
 *             of the right hand side updated by the host;
 *   urhsimhilolo_h are the second lowest doubles of the imaginary parts
 *             of the right hand side updated by the host;
 *   urhsimlololo_h are the lowest doubles of the imaginary parts
 *             of the right hand side updated by the host;
 *   urhsrehihihi_d are the highest doubles of the real parts
 *             of the right hand side updated by the device;
 *   urhsrelohi_d are the second highest doubles of the real parts
 *             of the right hand side updated by the device;
 *   urhsrelohi_d are the third highest doubles of the real parts
 *             of the right hand side updated by the device;
 *   urhsrelohi_d are the fourth highest doubles of the real parts
 *             of the right hand side updated by the device;
 *   urhsrehilo_d are the fourth lowest doubles of the real parts
 *             of the right hand side updated by the device;
 *   urhsrehilo_d are the third lowest doubles of the real parts
 *             of the right hand side updated by the device;
 *   urhsrehilo_d are the second lowest doubles of the real parts
 *             of the right hand side updated by the device;
 *   urhsrelololo_d are the lowest doubles of the real parts
 *             of the right hand side updated by the device;
 *   urhsimhihihi_d are the highest doubles of the imaginary parts
 *             of the right hand side updated by the device;
 *   urhsimlohihi_d are the second highest doubles of the imaginary parts
 *             of the right hand side updated by the device;
 *   urhsimhilohi_d are the third highest doubles of the imaginary parts
 *             of the right hand side updated by the device;
 *   urhsimlolohi_d are the fourth highest doubles of the imaginary parts
 *             of the right hand side updated by the device;
 *   urhsimhihilo_d are the fourth lowest doubles of the imaginary parts
 *             of the right hand side updated by the device;
 *   urhsimlohilo_d are the third lowest doubles of the imaginary parts
 *             of the right hand side updated by the device;
 *   urhsimhilolo_d are the second lowest doubles of the imaginary parts
 *             of the right hand side updated by the device;
 *   urhsimlololo_d are the lowest doubles of the imaginary parts
 *             of the right hand side updated by the device;
 *   solrehihihi_h are the highest doubles of real parts of solution on host;
 *   solrelohihi_h are the 2nd highest doubles of real parts of sol on host;
 *   solrehilohi_h are the 3rd highest doubles of real parts of sol on host;
 *   solrelolohi_h are the 4th highest doubles of real parts of sol on host;
 *   solrehihilo_h are the 4th lowest doubles of real parts of sol on host;
 *   solrelohilo_h are the 3rd lowest doubles of real parts of sol on host;
 *   solrehilolo_h are the 2nd lowest doubles of real parts of sol on host;
 *   solrelololo_h are the lowest doubles of real parts of solution on host;
 *   solimhihihi_h are the highest doubles of imag parts of solution on host;
 *   solimlohihi_h are the 2nd highest doubles of imag parts of sol on host;
 *   solimhilohi_h are the 3rd highest doubles of imag parts of sol on host;
 *   solimlolohi_h are the 4th highest doubles of imag parts of sol on host;
 *   solimhihilo_h are the 4th lowest doubles of imag parts of sol on host;
 *   solimlohilo_h are the 3rd lowest doubles of imag parts of sol on host;
 *   solimhilolo_h are the 2nd lowest doubles of imag parts of sol on host;
 *   solimlololo_h are the lowest doubles of imag parts of solution on host;
 *   solrehihihi_d are the highest doubles of real parts of solution on device;
 *   solrelohihi_d are the 2nd highest doubles of real parts of sol on device;
 *   solrehilohi_d are the 3rd highest doubles of real parts of sol on device;
 *   solrelolohi_d are the 4th highest doubles of real parts of sol on device;
 *   solrehihilo_d are the 4th lowest doubles of real parts of sol on device;
 *   solrelohilo_d are the 3rd lowest doubles of real parts of sol on device;
 *   solrehilolo_d are the 2nd lowest doubles of real parts of sol on device;
 *   solrelololo_d are the lowest doubles of real parts of solution on device;
 *   solimhihihi_d are the highest doubles of imag parts of solution on device;
 *   solimlohihi_d are the 2nd highest doubles of imag parts of sol on device;
 *   solimhilohi_d are the 3rd highest doubles of imag parts of sol on device;
 *   solimlolohi_d are the 4th highest doubles of imag parts of sol on device;
 *   solimhihilo_d are the 4th lowest doubles of imag parts of sol on device;
 *   solimlohilo_d are the 3rd lowest doubles of imag parts of sol on device;
 *   solimhilolo_d are the 2nd lowest doubles of imag parts of sol on device;
 *   solimlololo_d are the lowest doubles of imag parts of solution on device;
 *   Qrehihihi_h are the highest doubles of real Q of the QR on host;
 *   Qrelohihi_h are the second highest doubles of real Q of the QR on host;
 *   Qrehilohi_h are the third highest doubles of real Q of the QR on host;
 *   Qrelolohi_h are the fourth highest doubles of real Q of the QR on host;
 *   Qrehihilo_h are the fourth lowest doubles of real Q of the QR on host;
 *   Qrelohilo_h are the third lowest doubles of real Q of the QR on host;
 *   Qrehilolo_h are the second lowest doubles of real Q of the QR on host;
 *   Qrelololo_h are the lowest doubles of real Q of the QR on host;
 *   Qimhihihi_h are the highest doubles of imaginary Q of the QR on host;
 *   Qimlohihi_h are the second highest doubles of imag Q of the QR on host;
 *   Qimhilohi_h are the third highest doubles of imag Q of the QR on host;
 *   Qimlolohi_h are the fourth highest doubles of imag Q of the QR on host;
 *   Qimhihilo_h are the fourth lowest doubles of imag Q of the QR on host;
 *   Qimlohilo_h are the third lowest doubles of imag Q of the QR on host;
 *   Qimhilolo_h are the second lowest doubles of imag Q of the QR on host;
 *   Qimlololo_h are the lowest doubles of imaginary Q of the QR on host;
 *   Qrehihihi_d are the highest doubles of real Q of the QR on device;
 *   Qrelohihi_d are the second highest doubles of real Q of the QR on device;
 *   Qrehilohi_d are the third highest doubles of real Q of the QR on device;
 *   Qrelolohi_d are the fourth highest doubles of real Q of the QR on device;
 *   Qrehihilo_d are the fourth lowest doubles of real Q of the QR on device;
 *   Qrelohilo_d are the third lowest doubles of real Q of the QR on device;
 *   Qrehilolo_d are the second lowest doubles of real Q of the QR on device;
 *   Qrelololo_d are the lowest doubles of real Q of the QR on device;
 *   Qimhihihi_d are the highest doubles of imaginary Q of the QR on device;
 *   Qimlohihi_d are the second highest doubles of imag Q of the QR on device;
 *   Qimhilohi_d are the third highest doubles of imag Q of the QR on device;
 *   Qimlolohi_d are the fourth highest doubles of imag Q of the QR on device;
 *   Qimhihilo_d are the fourth lowest doubles of imag Q of the QR on device;
 *   Qimlohilo_d are the third lowest doubles of imag Q of the QR on device;
 *   Qimhilolo_d are the second lowest doubles of imag Q of the QR on device;
 *   Qimlololo_d are the lowest doubles of imaginary Q of the QR on device;
 *   Rrehihihi_h are the highest doubles of real R of the QR on host;
 *   Rrelohihi_h are the second highest doubles of real R of the QR on host;
 *   Rrehilohi_h are the third highest doubles of real R of the QR on host;
 *   Rrelolohi_h are the fourth highest doubles of real R of the QR on host;
 *   Rrehihilo_h are the fourth lowest doubles of real R of the QR on host;
 *   Rrelohilo_h are the third lowest doubles of real R of the QR on host;
 *   Rrehilolo_h are the second lowest doubles of real R of the QR on host;
 *   Rrelololo_h are the lowest doubles of real R of the QR on host;
 *   Rimhihihi_h are the highest doubles of imaginary R of the QR on host;
 *   Rimlohihi_h are the second highest doubles of imag R of the QR on host;
 *   Rimhilohi_h are the third highest doubles of imag R of the QR on host;
 *   Rimlolohi_h are the fourth highest doubles of imag R of the QR on host;
 *   Rimhihilo_h are the fourth lowest doubles of imag R of the QR on host;
 *   Rimlohilo_h are the third lowest doubles of imag R of the QR on host;
 *   Rimhilolo_h are the second lowest doubles of imag R of the QR on host;
 *   Rimlololo_h are the lowest doubles of imaginary R of the QR on host;
 *   Rrehihihi_d are the highest doubles of real R of the QR on device;
 *   Rrelohihi_d are the second highest doubles of real R of the QR on device;
 *   Rrehilohi_d are the third highest doubles of real R of the QR on device;
 *   Rrelolohi_d are the fourth highest doubles of real R of the QR on device;
 *   Rrehihilo_d are the fourth lowest doubles of real R of the QR on device;
 *   Rrelohilo_d are the third lowest doubles of real R of the QR on device;
 *   Rrehilolo_d are the second lowest doubles of real R of the QR on device;
 *   Rrelololo_d are the lowest doubles of real R of the QR on device;
 *   Rimhihihi_d are the highest doubles of imaginary R of the QR on device;
 *   Rimlohihi_d are the second highest doubles of imag R of the QR on device;
 *   Rimhilohi_d are the third highest doubles of imag R of the QR on device;
 *   Rimlolohi_d are the fourth highest doubles of imag R of the QR on device;
 *   Rimhihilo_d are the fourth lowest doubles of imag R of the QR on device;
 *   Rimlohilo_d are the third lowest doubles of imag R of the QR on device;
 *   Rimhilolo_d are the second lowest doubles of imag R of the QR on device;
 *   Rimlololo_d are the lowest doubles of imaginary R of the QR on device;
 *   wrkmat    has a copy of the Jacobian matrix;
 *   resvecrehihihi are the highest doubles of the real parts of resvec,
 *             which are the residual vectors;
 *   resvecrelohihi are second highest doubles of the real parts of resvec;
 *   resvecrehilohi are third highest doubles of the real parts of resvec;
 *   resvecrelolohi are fourth highest doubles of the real parts of resvec;
 *   resvecrehihilo are fourth lowest doubles of the real parts of resvec;
 *   resvecrelohilo are third lowest doubles of the real parts of resvec;
 *   resvecrehilolo are second lowest doubles of the real parts of resvec;
 *   resvecrelololo are lowest doubles of the real parts of resvec;
 *   resvecimhihihi are highest doubles of the imag parts of resvec;
 *   resvecimlohihi are second highest doubles of the imag parts of resvec;
 *   resvecimhilohi are third highest doubles of the imag parts of resvec;
 *   resvecimlolohi are fourth highest doubles of the imag parts of resvec;
 *   resvecimhihilo are fourth lowest doubles of the imag parts of resvec;
 *   resvecimlohilo are third lowest doubles of the imag parts of resvec;
 *   resvecimhilolo are second lowest doubles of the imag parts of resvec;
 *   resvecimlololo are lowest doubles of the imag parts of resvec;
 *   resmaxhihihi is the highest double of the residual max norm;
 *   resmaxlohihi is the second highest double of the residual max norm;
 *   resmaxhilohi is the third highest double of the residual max norm;
 *   resmaxlolohi is the fourth highest double of the residual max norm;
 *   resmaxhihilo is the fourth lowest double of the residual max norm;
 *   resmaxlohilo is the third lowest double of the residual max norm;
 *   resmaxhilolo is the second lowest double of the residual max norm;
 *   resmaxlololo is the lowest double of the residual max norm. */

int test_dbl8_real_newton
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

int test_dbl8_complex_newton
 ( int szt, int nbt, int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   int nbsteps, int mode, int vrblvl );
/*
 * DESCRIPTION :
 *   Runs Newton on a monomial system with complex octo double arithmetic.
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
