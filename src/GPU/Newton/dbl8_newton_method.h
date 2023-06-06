// The file dbl8_newton_method.h specifies Newton's method
// on series in octo double precision on real numbers.

#ifndef __dbl8_newton_method_h__
#define __dbl8_newton_method_h__

int dbl8_errors_funjacrhs
 ( int dim, int deg,
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
   double **rhshilolo_d, double **rhslololo_d, int vrblvl );
/*
 * DESCRIPTION :
 *   Compares results on host with those on the device
 *   and returns the number of errors larger than 1.0e-100.
 *
 * ON ENTRY :
 *   dim       number of equations and variables;
 *   deg       degree at which the series are truncated;
 *   funvalhihihi_h are the highest function values on the host;
 *   funvallohihi_h are the second highest function values on the host;
 *   funvalhilohi_h are the third highest function values on the host;
 *   funvallolohi_h are the fourth highest function values on the host;
 *   funvalhihilo_h are the fourth lowest function values on the host;
 *   funvallohilo_h are the third lowest function values on the host;
 *   funvalhilolo_h are the second lowest function values on the host;
 *   funvallololo_h are the lowest function values on the host;
 *   funvalhihihi_d are the highest function values on the device;
 *   funvallohihi_d are the second highest function values on the device;
 *   funvalhilohi_d are the third highest function values on the device;
 *   funvallolohi_d are the fourth highest function values on the device;
 *   funvalhihilo_d are the fourth lowest function values on the device;
 *   funvallohilo_d are the third lowest function values on the device;
 *   funvalhilolo_d are the second lowest function values on the device;
 *   funvallololo_d are the lowest function values on the device;
 *   jacvalhihihi_h are the highest doubles of the Jacobian on the host;
 *   jacvallohihi_h are the second highest doubles of the Jacobian on the host;
 *   jacvalhilohi_h are the third highest doubles of the Jacobian on the host;
 *   jacvallolohi_h are the fourth highest doubles of the Jacobian on the host;
 *   jacvalhihilo_h are the fourth lowest doubles of the Jacobian on the host;
 *   jacvallohilo_h are the third lowest doubles of the Jacobian on the host;
 *   jacvalhilolo_h are the second lowest doubles of the Jacobian on the host;
 *   jacvallololo_h are the lowest doubles of the Jacobian on the host;
 *   jacvalhihihi_d are the highest doubles of the Jacobian on device;
 *   jacvallohihi_d are the second highest doubles of the Jacobian on device;
 *   jacvalhilohi_d are the third highest doubles of the Jacobian on device;
 *   jacvallolohi_d are the fourth highest doubles of the Jacobian on device;
 *   jacvalhihilo_d are the fourth lowest doubles of the Jacobian on device;
 *   jacvallohilo_d are the third lowest doubles of the Jacobian on device;
 *   jacvalhilolo_d are the second lowest doubles of the Jacobian on device;
 *   jacvallololo_d are the lowest doubles of the Jacobian on device;
 *   rhshihihi_h are the highest doubles of the right hand side on host;
 *   rhslohihi_h are the 2nd highest doubles of the right hand side on host;
 *   rhshilohi_h are the 3rd highest doubles of the right hand side on host;
 *   rhslolohi_h are the 4th highest doubles of the right hand side on host;
 *   rhshihilo_h are the 4th lowest doubles of the right hand side on host;
 *   rhslohilo_h are the 3rd lowest doubles of the right hand side on host;
 *   rhshilolo_h are the 2nd lowest doubles of the right hand side on host;
 *   rhslololo_h are the lowest doubles of the right hand side on host;
 *   rhshihihi_d are the highest doubles of the right hand side on device;
 *   rhslohihi_d are the 2nd highest doubles of the right hand side on device;
 *   rhshilohi_d are the 3rd highest doubles of the right hand side on device;
 *   rhslolhi_d are the 4th highest doubles of the right hand side on device;
 *   rhshihilo_d are the 4th lowest doubles of the right hand side on device;
 *   rhslohilo_d are the 3rd lowest doubles of the right hand side on device;
 *   rhshilolo_d are the 2nd lowest doubles of the right hand side on device;
 *   rhslololo_d are the lowest doubles of the right hand side on device;
 *   vrblvl    is the verbose level. */

int dbl8_errors_inurhsQRsol
 ( int dim, int deg,
   double **inputhihihi_h, double **inputlohihi_h,
   double **inputhilohi_h, double **inputlolohi_h,
   double **inputhihilo_h, double **inputlohilo_h,
   double **inputhilolo_h, double **inputlololo_h,
   double **inputhihihi_d, double **inputlohihi_d,
   double **inputhilohi_d, double **inputlolohi_d,
   double **inputhihilo_d, double **inputlohilo_d,
   double **inputhilolo_d, double **inputlololo_d,
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
   double **solhilolo_d, double **sollololo_d, int vrblvl );
/*
 * DESCRIPTION :
 *   Compares results on host with those on the device
 *   and returns the number of errors larger than 1.0e-8.
 *
 * ON ENTRY :
 *   dim       number of equations and variables;
 *   deg       degree at which the series are truncated;
 *   inputhihihi_h are the highest doubles of updated series on the host;
 *   inputlohihi_h are the 2nd highest doubles of updated series on the host;
 *   inputhilohi_h are the 3rd highest doubles of updated series on the host;
 *   inputlolohi_h are the 4th highest doubles of updated series on the host;
 *   inputhihilo_h are the 4th lowest doubles of updated series on the host;
 *   inputlohilo_h are the 3rd lowest doubles of updated series on the host;
 *   inputhilolo_h are the 2nd lowest doubles of updated series on the host;
 *   inputlololo_h are the lowest doubles of updated series on the host;
 *   inputhihihi_d are the highest doubles of updated series on device;
 *   inputlohihi_d are the 2nd highest doubles of updated series on device;
 *   inputhilohi_d are the 3rd highest doubles of updated series on device;
 *   inputlolohi_d are the 4th highest doubles of updated series on device;
 *   inputhihilo_d are the 4th lowest doubles of updated series on device;
 *   inputlohilo_d are the 3rd lowest doubles of updated series on device;
 *   inputhilolo_d are the 2nd lowest doubles of updated series on device;
 *   inputlololo_d are the lowest doubles of updated series on device;
 *   Qhihihi_h   highest doubles of Q on the host;
 *   Qlohihi_h   second highest doubles of Q on the host;
 *   Qhilohi_h   third highest doubles of Q on the host;
 *   Qlolohi_h   fourth highest doubles of Q on the host;
 *   Qhihilo_h   fourth lowest doubles of Q on the host;
 *   Qlohilo_h   third lowest doubles of Q on the host;
 *   Qhilolo_h   second lowest doubles of Q on the host;
 *   Qlololo_h   lowest doubles of Q on the host;
 *   Qhihihi_d   highest doubles of Q on the device;
 *   Qlohihi_d   second highest doubles of Q on the device;
 *   Qhilohi_d   third highest doubles of Q on the device;
 *   Qlolohi_d   fourth highest doubles of Q on the device;
 *   Qhihilo_d   fourth lowest doubles of Q on the device;
 *   Qlohilo_d   third lowest doubles of Q on the device;
 *   Qhilolo_d   second lowest doubles of Q on the device;
 *   Qlololo_d   lowest doubles of Q on the device;
 *   Rhihihi_h   highest doubles of R on the host;
 *   Rlohihi_h   second highest doubles of R on the host;
 *   Rhilohi_h   third highest doubles of R on the host;
 *   Rlolohi_h   fourth highest doubles of R on the host;
 *   Rhihilo_h   fourth lowest doubles of R on the host;
 *   Rlohilo_h   third lowest doubles of R on the host;
 *   Rhilolo_h   second lowest doubles of R on the host;
 *   Rlololo_h   lowest doubles of R on the host;
 *   Rhihihi_d   highest doubles of R on the device;
 *   Rlohihi_d   second highest doubles of R on the device;
 *   Rhilohi_d   third highest doubles of R on the device;
 *   Rlolohi_d   fourth highest doubles of R on the device;
 *   Rhihilo_d   fourth lowest doubles of R on the device;
 *   Rlohilo_d   third lowest doubles of R on the device;
 *   Rhilolo_d   second lowest doubles of R on the device;
 *   Rlololo_d   lowest doubles of R on the device;
 *   urhshihihi_h are highest doubles of updated right hand side on host;
 *   urhslohihi_h are 2nd highest doubles of updated right hand side on host;
 *   urhshilohi_h are 3rd highest doubles of updated right hand side on host;
 *   urhslolohi_h are 4th highest doubles of updated right hand side on host;
 *   urhshihilo_h are 4th lowest doubles of updated right hand side on host;
 *   urhslohilo_h are 3rd lowest doubles of updated right hand side on host;
 *   urhshillo_h are 2nd lowest doubles of updated right hand side on host;
 *   urhslololo_h are lowest doubles of updated right hand side on host;
 *   urhshihihi_d are highest doubles of updated right hand side on device;
 *   urhslohihi_d are 2nd highest doubles of updated right hand side on device;
 *   urhshilohi_d are 3rd highest doubles of updated right hand side on device;
 *   urhslolohi_d are 4th highest doubles of updated right hand side on device;
 *   urhshihilo_d are 4th lowest doubles of updated right hand side on device;
 *   urhslohilo_d are 3rd lowest doubles of updated right hand side on device;
 *   urhshilolo_d are 2nd lowest doubles of updated right hand side on device;
 *   urhslololo_d are lowest doubles of updated right hand side on device;
 *   solhihihi_h are highest doubles of update to the solution on host,
 *   sollohihi_h are 2nd highest doubles of update to the solution on host,
 *   solhilohi_h are 3rd highest doubles of update to the solution on host,
 *   sollolohi_h are 4th highest doubles of update to the solution on host,
 *   solhihilo_h are 4th lowest doubles of update to the solution on host,
 *   sollohilo_h are 3rd lowest doubles of update to the solution on host,
 *   solhilolo_h are 2nd lowest doubles of update to the solution on host,
 *   sollololo_h are lowest doubles of update to the solution on host,
 *   solhihihi_d are highest doubles of update to the solution on device;
 *   sollohihi_d are 2nd highest doubles of update to the solution on device;
 *   solhilohi_d are 3rd highest doubles of update to the solution on device;
 *   sollolohi_d are 4th highest doubles of update to the solution on device;
 *   solhihilo_d are 4th lowest doubles of update to the solution on device;
 *   sollohilo_d are 3rd lowest doubles of update to the solution on device;
 *   solhilolo_d are 2nd lowest doubles of update to the solution on device;
 *   sollololo_d are lowest doubles of update to the solution on device;
 *   vrblvl    is the verbose level. */

int dbl8_update_newton_qrstep
 ( int szt, int nbt, int dim, int deg,
   int *tailidx_h, int *tailidx_d,
   double **inputhihihi_h, double **inputlohihi_h,
   double **inputhilohi_h, double **inputlolohi_h,
   double **inputhihilo_h, double **inputlohilo_h,
   double **inputhilolo_h, double **inputlololo_h,
   double **inputhihihi_d, double **inputlohihi_d,
   double **inputhilohi_d, double **inputlolohi_d,
   double **inputhihilo_d, double **inputlohilo_d,
   double **inputhilolo_d, double **inputlololo_d,
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
 * REQUIRED : szt*nbt = dim for GPU computing.
 *
 * ON ENTRY :
 *   szt       size of each tile and block;
 *   nbt       number of tiles and number of blocks;
 *   dim       number of monomials;
 *   deg       degree of the power series;
 *   tailidx_h is the start index of the update of the tail on the host;
 *   tailidx_d is the start index of the update of the tail on the device;
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
 *   resmaxlololo is the lowest double of the maximum of residual;
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

int dbl8_column_newton_qrstep
 ( int szt, int nbt, int dim, int deg, int nbrcol,
   int *tailidx_h, int *tailidx_d,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac,
   double **mbhihihi, double **mblohihi, double **mbhilohi, double **mblolohi,
   double **mbhihilo, double **mblohilo, double **mbhilolo, double **mblololo,
   double dpr,
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
 *   mbhihihi  highest doubles of the right hand side of monomial system;
 *   mblohihi  second highest doubles of the right hand side;
 *   mbhilohi  third highest doubles of the right hand side;
 *   mblolohi  fourth highest doubles of the right hand side;
 *   mbhihilo  fourth lowest doubles of the right hand side;
 *   mblohilo  third lowest doubles of the right hand side;
 *   mbhilolo  second lowest doubles of the right hand side;
 *   mblololo  lowest doubles of the right hand side;
 *   dpr       damper multiplier for t, should be in (0.0, 1.0];
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
 *   resmaxlololo is the lowest double of the maximum of residual;
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

int dbl8_row_newton_qrstep
 ( int szt, int nbt, int dim, int deg, int *tailidx_h, int *tailidx_d,
   int *nbr, int **nvr, int ***idx, 
   double **csthihihi, double **cstlohihi,
   double **csthilohi, double **cstlolohi,
   double **csthihilo, double **cstlohilo,
   double **csthilolo, double **cstlololo,
   double ***cffhihihi, double ***cfflohihi,
   double ***cffhilohi, double ***cfflolohi,
   double ***cffhihilo, double ***cfflohilo,
   double ***cffhilolo, double ***cfflololo, double dpr,
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
   bool *zeroQ_h, bool *noqr_h, bool *zeroQ_d, bool *noqr_d,
   int *upidx_h, int *bsidx_h, int *upidx_d, int *bsidx_d,
   double *totcnvlapsedms, double *totqrlapsedms, double *totqtblapsedms,
   double *totbslapsedms, double *totupdlapsedms, double *totreslapsedms,
   int vrblvl, int mode );
/*
 * DESCRIPTION :
 *   Does one step with Newton's method to update a power series,
 *   using QR factorization to solve linear systems,
 *   on an indexed polynomial system.
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
 *   csthihihi are the highest doubles of the constants;
 *   cstlohihi are the second highest doubles of the constants;
 *   csthilohi are the third highest doubles of the constants;
 *   cstlolohi are the fourth highest doubles of the constants;
 *   csthihilo are the fourth lowest doubles of the constants;
 *   cstlohilo are the third lowest doubles of the constants;
 *   csthilolo are the second lowest doubles of the constants;
 *   cstlololo are the lowest doubles of the constants;
 *   cffhihihi are the highest doubles of the coefficients of the monomials;
 *   cfflohihi are the second highest doubles of the coefficients;
 *   cffhilohi are the thirnd highest doubles of the coefficients;
 *   cfflolohi are the fourth highest doubles of the coefficients;
 *   cffhihilo are the fourth lowest doubles of the coefficients;
 *   cfflohilo are the third lowest doubles of the coefficients;
 *   cffhilolo are the second lowest doubles of the coefficients;
 *   cfflololo are the lowest doubles of the coefficients of the monomials;
 *   dpr       damper multiplier for t, should be in (0.0, 1.0];
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
 *   resmaxlololo is the lowest double of the maximum of residual;
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

int dbl8_allocate_inoutfunjac
 ( int dim, int deg, int mode,
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
   double ***jacvalhilolo_d, double ***jacvallololo_d );
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
 *   inputhihihi_h are the highest doubles of the input on the host;
 *   inputlohihi_h are the 2nd highest doubles of the input on the host;
 *   inputhilohi_h are the 3rd highest doubles of the input on the host;
 *   inputlolohi_h are the 4th highest doubles of the input on the host;
 *   inputhihilo_h are the 4th lowest doubles of the input on the host;
 *   inputlohilo_h are the 3rd lowest doubles of the input on the host;
 *   inputhilolo_h are the 2nd lowest doubles of the input on the host;
 *   inputlololo_h are the lowest doubles of the input on the host;
 *   inputhihihi_d are the highest doubles of the input on the device;
 *   inputlohihi_d are the 2nd highest doubles of the input on the device;
 *   inputhilohi_d are the 3rd highest doubles of the input on the device;
 *   inputlolohi_d are the 4th highest doubles of the input on the device;
 *   inputhihilo_d are the 4th lowest doubles of the input on the device;
 *   inputlohilo_d are the 3rd lowest doubles of the input on the device;
 *   inputhilolo_d are the 2nd lowest doubles of the input on the device;
 *   inputlololo_d are the lowest doubles of the input on the device;
 *   outputhihihi_h are the highest doubles of the output on the host;
 *   outputlohihi_h are the 2nd highest doubles of the output on the host;
 *   outputhilohi_h are the 3rd highest doubles of the output on the host;
 *   outputlolohi_h are the 4th highest doubles of the output on the host;
 *   outputhihilo_h are the 4th lowest doubles of the output on the host;
 *   outputlohilo_h are the 3rd lowest doubles of the output on the host;
 *   outputhilolo_h are the 2nd lowest doubles of the output on the host;
 *   outputlololo_h are the lowest doubles of the output on the host;
 *   outputhihihi_d are the highest doubles of the output on the device;
 *   outputlohihi_d are the 2nd highest doubles of the output on the device;
 *   outputhilohi_d are the 3rd highest doubles of the output on the device;
 *   outputlolohi_d are the 4th highest doubles of the output on the device;
 *   outputhihilo_d are the 4th lowest doubles of the output on the device;
 *   outputlohilo_d are the 3rd lowest doubles of the output on the device;
 *   outputhilolo_d are the 2nd lowest doubles of the output on the device;
 *   outputlololo_d are the lowest doubles of the output on the device;
 *   funvalhihihi_h are the highest function values on the host;
 *   funvallohihi_h are the second highest function values on the host;
 *   funvalhilohi_h are the third highest function values on the host;
 *   funvallolohi_h are the fourth highest function values on the host;
 *   funvalhihilo_h are the fourth lowest function values on the host;
 *   funvallohilo_h are the third lowest function values on the host;
 *   funvalhilolo_h are the second lowest function values on the host;
 *   funvallololo_h are the lowest function values on the host;
 *   funvalhihihi_d are the highest function values on the device;
 *   funvallohihi_d are the second highest function values on the device;
 *   funvalhilohi_d are the third highest function values on the device;
 *   funvallolohi_d are the fourth highest function values on the device;
 *   funvalhihilo_d are the fourth lowest function values on the device;
 *   funvallohilo_d are the third lowest function values on the device;
 *   funvalhilolo_d are the second lowest function values on the device;
 *   funvallololo_d are the lowest function values on the device;
 *   jacvalhihihi_h are the highest doubles of the Jacobian on the host;
 *   jacvallohihi_h are the second highest doubles of the Jacobian on the host;
 *   jacvalhilohi_h are the third highest doubles of the Jacobian on the host;
 *   jacvallolohi_h are the fourth highest doubles of the Jacobian on the host;
 *   jacvalhihilo_h are the fourth lowest doubles of the Jacobian on the host;
 *   jacvallohilo_h are the third lowest doubles of the Jacobian on the host;
 *   jacvalhilolo_h are the second lowest doubles of the Jacobian on the host;
 *   jacvallololo_h are the lowest doubles of the Jacobian on the host;
 *   jacvalhihihi_d are the highest doubles of the Jacobian on device;
 *   jacvallohihi_d are the second highest doubles of the Jacobian on device;
 *   jacvalhilohi_d are the third highest doubles of the Jacobian on device;
 *   jacvallolohi_d are the fourth highest doubles of the Jacobian on device;
 *   jacvalhihilo_d are the fourth lowest doubles of the Jacobian on device;
 *   jacvallohilo_d are the third lowest doubles of the Jacobian on device;
 *   jacvalhilolo_d are the second lowest doubles of the Jacobian on device;
 *   jacvallololo_d are the lowest doubles of the Jacobian on device. */

int dbl8_allocate_rhsqrsol
 ( int dim, int deg, int mode,
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
   double **solhihihi_h, double **sollohihi_h,
   double **solhilohi_h, double **sollolohi_h,
   double **solhihilo_h, double **sollohilo_h,
   double **solhilolo_h, double **sollololo_h,
   double **solhihihi_d, double **sollohihi_d,
   double **solhilohi_d, double **sollolohi_d,
   double **solhihilo_d, double **sollohilo_d,
   double **solhilolo_d, double **sollololo_d );
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
 *   rhshihihi_h are the highest doubles of the right hand side on host;
 *   rhslohihi_h are the 2nd highest doubles of the right hand side on host;
 *   rhshilohi_h are the 3rd highest doubles of the right hand side on host;
 *   rhslolohi_h are the 4th highest doubles of the right hand side on host;
 *   rhshihilo_h are the 4th lowest doubles of the right hand side on host;
 *   rhslohilo_h are the 3rd lowest doubles of the right hand side on host;
 *   rhshilolo_h are the 2nd lowest doubles of the right hand side on host;
 *   rhslololo_h are the lowest doubles of the right hand side on host;
 *   rhshihihi_d are the highest doubles of the right hand side on device;
 *   rhslohihi_d are the 2nd highest doubles of the right hand side on device;
 *   rhshilohi_d are the 3rd highest doubles of the right hand side on device;
 *   rhslolhi_d are the 4th highest doubles of the right hand side on device;
 *   rhshihilo_d are the 4th lowest doubles of the right hand side on device;
 *   rhslohilo_d are the 3rd lowest doubles of the right hand side on device;
 *   rhshilolo_d are the 2nd lowest doubles of the right hand side on device;
 *   rhslololo_d are the lowest doubles of the right hand side on device;
 *   urhshihihi_h are highest doubles of updated right hand side on host;
 *   urhslohihi_h are 2nd highest doubles of updated right hand side on host;
 *   urhshilohi_h are 3rd highest doubles of updated right hand side on host;
 *   urhslolohi_h are 4th highest doubles of updated right hand side on host;
 *   urhshihilo_h are 4th lowest doubles of updated right hand side on host;
 *   urhslohilo_h are 3rd lowest doubles of updated right hand side on host;
 *   urhshillo_h are 2nd lowest doubles of updated right hand side on host;
 *   urhslololo_h are lowest doubles of updated right hand side on host;
 *   urhshihihi_d are highest doubles of updated right hand side on device;
 *   urhslohihi_d are 2nd highest doubles of updated right hand side on device;
 *   urhshilohi_d are 3rd highest doubles of updated right hand side on device;
 *   urhslolohi_d are 4th highest doubles of updated right hand side on device;
 *   urhshihilo_d are 4th lowest doubles of updated right hand side on device;
 *   urhslohilo_d are 3rd lowest doubles of updated right hand side on device;
 *   urhshilolo_d are 2nd lowest doubles of updated right hand side on device;
 *   urhslololo_d are lowest doubles of updated right hand side on device;
 *   Qhihihi_h   highest doubles of Q on the host;
 *   Qlohihi_h   second highest doubles of Q on the host;
 *   Qhilohi_h   third highest doubles of Q on the host;
 *   Qlolohi_h   fourth highest doubles of Q on the host;
 *   Qhihilo_h   fourth lowest doubles of Q on the host;
 *   Qlohilo_h   third lowest doubles of Q on the host;
 *   Qhilolo_h   second lowest doubles of Q on the host;
 *   Qlololo_h   lowest doubles of Q on the host;
 *   Qhihihi_d   highest doubles of Q on the device;
 *   Qlohihi_d   second highest doubles of Q on the device;
 *   Qhilohi_d   third highest doubles of Q on the device;
 *   Qlolohi_d   fourth highest doubles of Q on the device;
 *   Qhihilo_d   fourth lowest doubles of Q on the device;
 *   Qlohilo_d   third lowest doubles of Q on the device;
 *   Qhilolo_d   second lowest doubles of Q on the device;
 *   Qlololo_d   lowest doubles of Q on the device;
 *   Rhihihi_h   highest doubles of R on the host;
 *   Rlohihi_h   second highest doubles of R on the host;
 *   Rhilohi_h   third highest doubles of R on the host;
 *   Rlolohi_h   fourth highest doubles of R on the host;
 *   Rhihilo_h   fourth lowest doubles of R on the host;
 *   Rlohilo_h   third lowest doubles of R on the host;
 *   Rhilolo_h   second lowest doubles of R on the host;
 *   Rlololo_h   lowest doubles of R on the host;
 *   Rhihihi_d   highest doubles of R on the device;
 *   Rlohihi_d   second highest doubles of R on the device;
 *   Rhilohi_d   third highest doubles of R on the device;
 *   Rlolohi_d   fourth highest doubles of R on the device;
 *   Rhihilo_d   fourth lowest doubles of R on the device;
 *   Rlohilo_d   third lowest doubles of R on the device;
 *   Rhilolo_d   second lowest doubles of R on the device;
 *   Rlololo_d   lowest doubles of R on the device;
 *   solhihihi_h are highest doubles of update to the solution on host,
 *   sollohihi_h are 2nd highest doubles of update to the solution on host,
 *   solhilohi_h are 3rd highest doubles of update to the solution on host,
 *   sollolohi_h are 4th highest doubles of update to the solution on host,
 *   solhihilo_h are 4th lowest doubles of update to the solution on host,
 *   sollohilo_h are 3rd lowest doubles of update to the solution on host,
 *   solhilolo_h are 2nd lowest doubles of update to the solution on host,
 *   sollololo_h are lowest doubles of update to the solution on host,
 *   solhihihi_d are highest doubles of update to the solution on device;
 *   sollohihi_d are 2nd highest doubles of update to the solution on device;
 *   solhilohi_d are 3rd highest doubles of update to the solution on device;
 *   sollolohi_d are 4th highest doubles of update to the solution on device;
 *   solhihilo_d are 4th lowest doubles of update to the solution on device;
 *   sollohilo_d are 3rd lowest doubles of update to the solution on device;
 *   solhilolo_d are 2nd lowest doubles of update to the solution on device;
 *   sollololo_d are lowest doubles of update to the solution on device. */

void dbl8_start_setup
 ( int dim, int deg,
   double **testsolhihihi, double **testsollohihi,
   double **testsolhilohi, double **testsollolohi,
   double **testsolhihilo, double **testsollohilo,
   double **testsolhilolo, double **testsollololo,
   double **inputhihihi_h, double **inputlohihi_h,
   double **inputhilohi_h, double **inputlolohi_h,
   double **inputhihilo_h, double **inputlohilo_h,
   double **inputhilolo_h, double **inputlololo_h,
   double **inputhihihi_d, double **inputlohihi_d,
   double **inputhilohi_d, double **inputlolohi_d,
   double **inputhihilo_d, double **inputlohilo_d,
   double **inputhilolo_d, double **inputlololo_d, int mode, int vrblvl );
/*
 * DESCRIPTION :
 *   Given the test solution, defines the start vector.
 *
 * ON ENTRY :
 *   dim        dimension of the solution vectors;
 *   deg        degree at which the series are truncated;
 *   testsolhihihi are the highest doubles of the test solution;
 *   testsollohihi are the second highest doubles of the test solution;
 *   testsolhilohi are the third highest doubles of the test solution;
 *   testsollolohi are the fourth highest doubles of the test solution;
 *   testsolhihilo are the fourth lowest doubles of the test solution;
 *   testsollohilo are the third lowest doubles of the test solution;
 *   testsolhilolo are the second lowest doubles of the test solution;
 *   testsollololo are the lowest doubles of the test solution;
 *   inputhihihi_h is allocated on host if mode is 1 or 2;
 *   inputlohihi_h is allocated on host if mode is 1 or 2;
 *   inputhilohi_h is allocated on host if mode is 1 or 2;
 *   inputlolohi_h is allocated on host if mode is 1 or 2;
 *   inputhihilo_h is allocated on host if mode is 1 or 2;
 *   inputlohilo_h is allocated on host if mode is 1 or 2;
 *   inputhilolo_h is allocated on host if mode is 1 or 2;
 *   inputlololo_h is allocated on host if mode is 1 or 2;
 *   inputhihihi_d is allocated on device if mode is 0 or 2;
 *   inputlohihi_d is allocated on device if mode is 0 or 2;
 *   inputhilohi_d is allocated on device if mode is 0 or 2;
 *   inputlolohi_d is allocated on device if mode is 0 or 2;
 *   inputhihilo_d is allocated on device if mode is 0 or 2;
 *   inputlohilo_d is allocated on device if mode is 0 or 2;
 *   inputhilolo_d is allocated on device if mode is 0 or 2;
 *   inputlololo_d is allocated on device if mode is 0 or 2;
 *   mode       0 (GPU), 1 (CPU), or 2 (GPU+CPU);
 *   vrblvl     is the verbose level.
 *
 * ON RETURN :
 *   inputhihihi_h are the highest doubles of start vector
 *              on host if mode is 1 or 2;
 *   inputlohihi_h are the second highest doubles of start vector
 *              on host if mode is 1 or 2;
 *   inputhilohi_h are the third highest doubles of start vector
 *              on host if mode is 1 or 2;
 *   inputlolohi_h are the fourth highest doubles of start vector
 *              on host if mode is 1 or 2;
 *   inputhihilo_h are the fourth lowest doubles of start vector
 *              on host if mode is 1 or 2;
 *   inputlohilo_h are the third lowest doubles of start vector
 *              on host if mode is 1 or 2;
 *   inputhilolo_h are the second lowest doubles of start vector
 *              on host if mode is 1 or 2;
 *   inputlololo_h are the lowest doubles of start vector
 *              on host if mode is 1 or 2;
 *   inputhihihi_d are the highest doubles of start vector
 *              on device if mode is 0 or 2;
 *   inputlohihi_d are the second highest doubles of start vector
 *              on device if mode is 0 or 2;
 *   inputhilohi_d are the third highest doubles of start vector
 *              on device if mode is 0 or 2;
 *   inputlohihi_d are the fourth highest doubles of start vector
 *              on device if mode is 0 or 2;
 *   inputhihilo_d are the fourth lowest doubles of start vector
 *              on device if mode is 0 or 2;
 *   inputlohilo_d are the third lowest doubles of start vector
 *              on device if mode is 0 or 2;
 *   inputhilolo_d are the second lowest doubles of start vector
 *              on device if mode is 0 or 2;
 *   inputlololo_d are the lowest doubles of start vector
 *              on device if mode is 0 or 2. */

void dbl8_column_setup
 ( int dim, int deg, int nbrcol, int **nvr, int ***idx, int **rowsA,
   double ***cffhihihi, double ***cfflohihi,
   double ***cffhilohi, double ***cfflolohi,
   double ***cffhihilo, double ***cfflohilo,
   double ***cffhilolo, double ***cfflololo,
   double **testsolhihihi, double **testsollohihi,
   double **testsolhilohi, double **testsollolohi,
   double **testsolhihilo, double **testsollohilo,
   double **testsolhilolo, double **testsollololo,
   double **mbrhshihihi, double **mbrhslohihi,
   double **mbrhshilohi, double **mbrhslolohi,
   double **mbrhshihilo, double **mbrhslohilo,
   double **mbrhshilolo, double **mbrhslololo,
   double **inputhihihi_h, double **inputlohihi_h,
   double **inputhilohi_h, double **inputlolohi_h,
   double **inputhihilo_h, double **inputlohilo_h,
   double **inputhilolo_h, double **inputlololo_h,
   double **inputhihihi_d, double **inputlohihi_d,
   double **inputhilohi_d, double **inputlolohi_d,
   double **inputhihilo_d, double **inputlohilo_d,
   double **inputhilolo_d, double **inputlololo_d, int mode, int vrblvl );
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
 *   cffhihihi  highest doubles of coefficients, if more than one column;
 *   cfflohihi  2nd highest doubles of coefficients, if more than one column;
 *   cffhilohi  3nd highest doubles of coefficients, if more than one column;
 *   cfflolohi  4th highest doubles of coefficients, if more than one column;
 *   cffhihilo  4th lowest doubles of coefficients, if more than one column;
 *   cfflohilo  3rd lowest doubles of coefficients, if more than one column;
 *   cffhilolo  2nd lowest doubles of coefficients, if more than one column;
 *   cfflololo  lowest doubles of coefficients, if more than one column;
 *   testsolhihihi has space for dim pointers;
 *   testsollohihi has space for dim pointers;
 *   testsolhilohi has space for dim pointers;
 *   testsollolohi has space for dim pointers;
 *   testsolhihilo has space for dim pointers;
 *   testsollohilo has space for dim pointers;
 *   testsolhilolo has space for dim pointers;
 *   testsollololo has space for dim pointers;
 *   mbrhshihihi has space for dim pointers;
 *   mbrhslohihi has space for dim pointers;
 *   mbrhshilohi has space for dim pointers;
 *   mbrhslolohi has space for dim pointers;
 *   mbrhshihilo has space for dim pointers;
 *   mbrhslohilo has space for dim pointers;
 *   mbrhshilolo has space for dim pointers;
 *   mbrhslololo has space for dim pointers;
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
 *   testsolhihihi are the highest doubles of the test solution;
 *   testsollohihi are the second highest doubles of the test solution;
 *   testsolhilohi are the third highest doubles of the test solution;
 *   testsollolohi are the fourth highest doubles of the test solution;
 *   testsolhihilo are the fourth lowest doubles of the test solution;
 *   testsollohilo are the third lowest doubles of the test solution;
 *   testsolhilolo are the second lowest doubles of the test solution;
 *   testsollololo are the lowest doubles of the test solution;
 *   mbrhshihihi are the highest doubles of right hand side vector
 *              for the test solution;
 *   mbrhslohihi are the second highest doubles of right hand side vector
 *              for the test solution;
 *   mbrhshilohi are the third highest doubles of right hand side vector
 *              for the test solution;
 *   mbrhslolohi are the fourth highest doubles of right hand side vector
 *              for the test solution;
 *   mbrhshihilo are the fourth lowest doubles of right hand side vector
 *              for the test solution;
 *   mbrhslohilo are the third lowest doubles of right hand side vector
 *              for the test solution;
 *   mbrhshilolo are the second lowest doubles of right hand side vector
 *              for the test solution;
 *   mbrhslololo are the lowest doubles of right hand side vector
 *              for the test solution;
 *   inputhihihi_h are the highest doubles of start vector
 *              on host if mode is 1 or 2;
 *   inputlohihi_h are the second highest doubles of start vector
 *              on host if mode is 1 or 2;
 *   inputhilohi_h are the third highest doubles of start vector
 *              on host if mode is 1 or 2;
 *   inputlolohi_h are the fourth highest doubles of start vector
 *              on host if mode is 1 or 2;
 *   inputhihilo_h are the fourth lowest doubles of start vector
 *              on host if mode is 1 or 2;
 *   inputlohilo_h are the third lowest doubles of start vector
 *              on host if mode is 1 or 2;
 *   inputhilolo_h are the second lowest doubles of start vector
 *              on host if mode is 1 or 2;
 *   inputlololo_h are the lowest doubles of start vector
 *              on host if mode is 1 or 2;
 *   inputhihihi_d are the highest doubles of start vector
 *              on device if mode is 0 or 2;
 *   inputlohihi_d are the second highest doubles of start vector
 *              on device if mode is 0 or 2;
 *   inputhilohi_d are the third highest doubles of start vector
 *              on device if mode is 0 or 2;
 *   inputlohihi_d are the fourth highest doubles of start vector
 *              on device if mode is 0 or 2;
 *   inputhihilo_d are the fourth lowest doubles of start vector
 *              on device if mode is 0 or 2;
 *   inputlohilo_d are the third lowest doubles of start vector
 *              on device if mode is 0 or 2;
 *   inputhilolo_d are the second lowest doubles of start vector
 *              on device if mode is 0 or 2;
 *   inputlololo_d are the lowest doubles of start vector
 *              on device if mode is 0 or 2. */

void dbl8_row_setup
 ( int dim, int deg, int *nbr, int **nvr, int ***idx,
   double **csthihihi, double **cstlohihi,
   double **csthilohi, double **cstlolohi,
   double **csthihilo, double **cstlohilo,
   double **csthilolo, double **cstlololo,
   double ***cffhihihi, double ***cfflohihi,
   double ***cffhilohi, double ***cfflolohi,
   double ***cffhihilo, double ***cfflohilo,
   double ***cffhilolo, double ***cfflololo,
   double **testsolhihihi, double **testsollohihi,
   double **testsolhilohi, double **testsollolohi,
   double **testsolhihilo, double **testsollohilo,
   double **testsolhilolo, double **testsollololo,
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
   double ***outputhilolo_d, double ***outputlololo_d, int mode, int vrblvl );
/*
 * DESCRIPTION :
 *   Defines the test solution and start vectors to run Newton's method
 *   on one or more columns of monomials.
 *
 * ON ENTRY :
 *   dim        dimension of the solution vectors;
 *   deg        degree at which the series are truncated;
 *   nbr        nbr[i] is the number of monomials of the i-th polynomials;
 *   nvr        nvr[i][j] is the number of variables in the j-th monomial
 *              of the i-th polynomial;
 *   idx        idx[i][j] are the indices of the variables in monomial j
 *              of the i-th polynomial;
 *   csthihihi  are the highest doubles of the constants;
 *   cstlohihi  are the second highest doubles of the constants;
 *   csthilohi  are the third highest doubles of the constants;
 *   cstlolohi  are the fourth highest doubles of the constants;
 *   csthihilo  are the fourth lowest doubles of the constants;
 *   cstlohilo  are the third lowest doubles of the constants;
 *   csthilolo  are the second lowest doubles of the constants;
 *   cstlololo  are the lowest doubles of the constants;
 *   cffhihihi  highest doubles of coefficients, if more than one column;
 *   cfflohihi  2nd highest doubles of coefficients, if more than one column;
 *   cffhilohi  3nd highest doubles of coefficients, if more than one column;
 *   cfflolohi  4th highest doubles of coefficients, if more than one column;
 *   cffhihilo  4th lowest doubles of coefficients, if more than one column;
 *   cfflohilo  3rd lowest doubles of coefficients, if more than one column;
 *   cffhilolo  2nd lowest doubles of coefficients, if more than one column;
 *   cfflololo  lowest doubles of coefficients, if more than one column;
 *   testsolhihihi has space for dim pointers;
 *   testsollohihi has space for dim pointers;
 *   testsolhilohi has space for dim pointers;
 *   testsollolohi has space for dim pointers;
 *   testsolhihilo has space for dim pointers;
 *   testsollohilo has space for dim pointers;
 *   testsolhilolo has space for dim pointers;
 *   testsollololo has space for dim pointers;
 *   inputhihihi_h is allocated on host if mode is 1 or 2;
 *   inputlohihi_h is allocated on host if mode is 1 or 2;
 *   inputhilohi_h is allocated on host if mode is 1 or 2;
 *   inputlolohi_h is allocated on host if mode is 1 or 2;
 *   inputhihilo_h is allocated on host if mode is 1 or 2;
 *   inputlohilo_h is allocated on host if mode is 1 or 2;
 *   inputhilolo_h is allocated on host if mode is 1 or 2;
 *   inputlololo_h is allocated on host if mode is 1 or 2;
 *   inputhihihi_d is allocated on device if mode is 0 or 2;
 *   inputlohihi_d is allocated on device if mode is 0 or 2;
 *   inputhilohi_d is allocated on device if mode is 0 or 2;
 *   inputlolohi_d is allocated on device if mode is 0 or 2;
 *   inputhihilo_d is allocated on device if mode is 0 or 2;
 *   inputlohilo_d is allocated on device if mode is 0 or 2;
 *   inputhilolo_d is allocated on device if mode is 0 or 2;
 *   inputlololo_d is allocated on device if mode is 0 or 2;
 *   outputhihihi_h is allocated on host if mode is 1 or 2;
 *   outputlohihi_h is allocated on host if mode is 1 or 2;
 *   outputhilohi_h is allocated on host if mode is 1 or 2;
 *   outputlolohi_h is allocated on host if mode is 1 or 2;
 *   outputhihilo_h is allocated on host if mode is 1 or 2;
 *   outputlohilo_h is allocated on host if mode is 1 or 2;
 *   outputhilolo_h is allocated on host if mode is 1 or 2;
 *   outputlololo_h is allocated on host if mode is 1 or 2;
 *   outputhihihi_d is allocated on device if mode is 0 or 2;
 *   outputlohihi_d is allocated on device if mode is 0 or 2;
 *   outputhilohi_d is allocated on device if mode is 0 or 2;
 *   outputlolohi_d is allocated on device if mode is 0 or 2;
 *   outputhihilo_d is allocated on device if mode is 0 or 2;
 *   outputlohilo_d is allocated on device if mode is 0 or 2;
 *   outputhilolo_d is allocated on device if mode is 0 or 2;
 *   outputlololo_d is allocated on device if mode is 0 or 2;
 *   mode       0 (GPU), 1 (CPU), or 2 (GPU+CPU);
 *   vrblvl     is the verbose level.
 *
 * ON RETURN :
 *   testsolhihihi are the highest doubles of the test solution;
 *   testsollohihi are the second highest doubles of the test solution;
 *   testsolhilohi are the third highest doubles of the test solution;
 *   testsollolohi are the fourth highest doubles of the test solution;
 *   testsolhihilo are the fourth lowest doubles of the test solution;
 *   testsollohilo are the third lowest doubles of the test solution;
 *   testsolhilolo are the second lowest doubles of the test solution;
 *   testsollololo are the lowest doubles of the test solution;
 *   inputhihihi_h are the highest doubles of start vector
 *              on host if mode is 1 or 2;
 *   inputlohihi_h are the second highest doubles of start vector
 *              on host if mode is 1 or 2;
 *   inputhilohi_h are the third highest doubles of start vector
 *              on host if mode is 1 or 2;
 *   inputlolohi_h are the fourth highest doubles of start vector
 *              on host if mode is 1 or 2;
 *   inputhihilo_h are the fourth lowest doubles of start vector
 *              on host if mode is 1 or 2;
 *   inputlohilo_h are the third lowest doubles of start vector
 *              on host if mode is 1 or 2;
 *   inputhilolo_h are the second lowest doubles of start vector
 *              on host if mode is 1 or 2;
 *   inputlololo_h are the lowest doubles of start vector
 *              on host if mode is 1 or 2;
 *   inputhihihi_d are the highest doubles of start vector
 *              on device if mode is 0 or 2;
 *   inputlohihi_d are the second highest doubles of start vector
 *              on device if mode is 0 or 2;
 *   inputhilohi_d are the third highest doubles of start vector
 *              on device if mode is 0 or 2;
 *   inputlohihi_d are the fourth highest doubles of start vector
 *              on device if mode is 0 or 2;
 *   inputhihilo_d are the fourth lowest doubles of start vector
 *              on device if mode is 0 or 2;
 *   inputlohilo_d are the third lowest doubles of start vector
 *              on device if mode is 0 or 2;
 *   inputhilolo_d are the second lowest doubles of start vector
 *              on device if mode is 0 or 2;
 *   inputlololo_d are the lowest doubles of start vector
 *              on device if mode is 0 or 2;
 *   outputhihihi_h are the highest doubles 
 *              of evaluated test solution if mode is 1;
 *   outputlohihi_h are the second highest doubles 
 *              of evaluated test solution if mode is 1;
 *   outputhilohi_h are the third highest doubles 
 *              of evaluated test solution if mode is 1;
 *   outputlolohi_h are the fourth highest doubles 
 *              of evaluated test solution if mode is 1;
 *   outputhihilo_h are the fourth lowest doubles
 *              of evaluated test solution if mode is 1;
 *   outputlohilo_h are the third lowest doubles
 *              of evaluated test solution if mode is 1;
 *   outputhilolo_h are the second lowest doubles
 *              of evaluated test solution if mode is 1;
 *   outputlololo_h are the lowest doubles
 *              of evaluated test solution if mode is 1;
 *   outputhihihi_d are the highest doubles
 *              of evaluated test solution, if mode is 0 or 2;
 *   outputlohihi_d are the second highest doubles
 *              of evaluated test solution, if mode is 0 or 2;
 *   outputhilohi_d are the third highest doubles
 *              of evaluated test solution, if mode is 0 or 2;
 *   outputlolohi_d are the fourth highest doubles
 *              of evaluated test solution, if mode is 0 or 2;
 *   outputhihilo_d are the fourth lowest doubles
 *              of evaluated test solution, if mode is 0 or 2;
 *   outputlohilo_d are the third lowest doubles
 *              of evaluated test solution, if mode is 0 or 2;
 *   outputhilolo_d are the second lowest doubles
 *              of evaluated test solution, if mode is 0 or 2;
 *   outputlololo_d are the lowest doubles
 *              of evaluated test solution, if mode is 0 or 2. */

int dbl8_error_testsol
 ( int dim, int deg, int mode,
   double **testsolhihihi, double **testsollohihi,
   double **testsolhilohi, double **testsollolohi,
   double **testsolhihilo, double **testsollohilo,
   double **testsolhilolo, double **testsollololo,
   double **inputhihihi_h, double **inputlohihi_h,
   double **inputhilohi_h, double **inputlolohi_h,
   double **inputhihilo_h, double **inputlohilo_h,
   double **inputhilolo_h, double **inputlololo_h,
   double **inputhihihi_d, double **inputlohihi_d,
   double **inputhilohi_d, double **inputlolohi_d,
   double **inputhihilo_d, double **inputlohilo_d,
   double **inputhilolo_d, double **inputlololo_d );
/*
 * DESCRIPTION :
 *   Compares the computed solution against the test solution.
 *
 * ON ENTRY :
 *   dim        dimension of the solution vectors;
 *   deg        degree at which the series are truncated;
 *   mode       0 (GPU), 1 (CPU), or 2 (GPU+CPU).
 *   testsolhihihi are the highest doubles of the test solution;
 *   testsollohihi are the second highest doubles of the test solution;
 *   testsolhilohi are the third highest doubles of the test solution;
 *   testsollolohi are the fourth highest doubles of the test solution;
 *   testsolhihilo are the fourth lowest doubles of the test solution;
 *   testsollohilo are the third lowest doubles of the test solution;
 *   testsolhilolo are the second lowest doubles of the test solution;
 *   testsollololo are the lowest doubles of the test solution;
 *   inputhihihi_h are the highest doubles on host if mode is 1 or 2;
 *   inputlohihi_h are the second highest doubles on host if mode is 1 or 2;
 *   inputhilohi_h are the third highest doubles on host if mode is 1 or 2;
 *   inputlolohi_h are the fourth highest doubles on host if mode is 1 or 2;
 *   inputhihilo_h are the fourth lowest doubles on host if mode is 1 or 2;
 *   inputlohilo_h are the third lowest doubles on host if mode is 1 or 2;
 *   inputhilolo_h are the second lowest doubles on host if mode is 1 or 2;
 *   inputlololo_h are the lowest doubles on host if mode is 1 or 2;
 *   inputhihihi_d are the highest doubles on device if mode is 0 or 2;
 *   inputlohihi_d are the second highest doubles on device if mode is 0 or 2;
 *   inputhilohi_d are the third highest doubles on device if mode is 0 or 2;
 *   inputlolohi_d are the fourth highest doubles on device if mode is 0 or 2;
 *   inputhihilo_d are the fourth lowest doubles on device if mode is 0 or 2;
 *   inputlohilo_d are the third lowest doubles on device if mode is 0 or 2;
 *   inputhilolo_d are the second lowest doubles on device if mode is 0 or 2;
 *   inputlololo_d are the lowest doubles on device if mode is 0 or 2. */

int test_dbl8_column_newton
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

int test_dbl8_row_newton
 ( int szt, int nbt, int dim, int deg, int *nbr, int **nvr, int ***idx,
   double dpr, int nbsteps, int mode, int vrblvl );
/*
 * DESCRIPTION :
 *   Runs Newton with real double arithmetic,
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
