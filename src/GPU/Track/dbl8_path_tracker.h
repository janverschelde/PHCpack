// The file dbl8_path_tracker.h specifies a path tracker
// on series in octo double precision on real numbers.

#ifndef __dbl8_path_tracker_h__
#define __dbl8_path_tracker_h__

int dbl8_run_newton
 ( int szt, int nbt, int dim, int deg, int nbrcol, int nbsteps,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac,
   double **rhshihihi, double **rhslohihi,
   double **rhshilohi, double **rhslolohi,
   double **rhshihilo, double **rhslohilo,
   double **rhshilolo, double **rhslololo,
   double ***cffhihihi, double ***cfflohihi,
   double ***cffhilohi, double ***cfflolohi,
   double ***cffhihilo, double ***cfflohilo,
   double ***cffhilolo, double ***cfflololo,
   double **acchihihi, double **acclohihi,
   double **acchilohi, double **acclolohi,
   double **acchihilo, double **acclohilo,
   double **acchilolo, double **acclololo,
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
   double *resmaxhilolo, double *resmaxlololo,
   int vrblvl, int mode );
/*
 * DESCRIPTION :
 *   Runs Newton's method on power series on real data
 *   in octo double precision.
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
 *   rhshihihi are the highest doubles of the right hand side;
 *   rhslohihi are the second highest doubles of the right hand side;
 *   rhshilohi are the third highest doubles of the right hand side;
 *   rhslolohi are the fourth highest doubles of the right hand side;
 *   rhshihilo are the fourth lowest doubles of the right hand side;
 *   rhslohilo are the third lowest doubles of the right hand side;
 *   rhshilolo are the second lowest doubles of the right hand side;
 *   rhslololo are the lowest doubles of the right hand side;
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

int test_dbl8_real_track
 ( int szt, int nbt, int dim, int deg, int nbrcol,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac, int **rowsA,
   int nbsteps, int mode, int vrblvl );
/*
 * DESCRIPTION :
 *   Tracks a path on a system with real octo double arithmetic.
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
