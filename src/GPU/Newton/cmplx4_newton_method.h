// The file cmplx4_newton_method.h specifies Newton's method
// on series in quad precision on complex numbers.

#ifndef __cmplx4_newton_method_h__
#define __cmplx4_newton_method_h__

int cmplx4_errors_funjacrhs
 ( int dim, int deg,
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
   double **rhsimhilo_d, double **rhsimlolo_d, int vrblvl );
/*
 * DESCRIPTION :
 *   Compares results on host with those on the device
 *   and returns the number of errors larger than 1.0e-100.
 *
 * ON ENTRY :
 *   dim       number of equations and variables;
 *   deg       degree at which the series are truncated;
 *   funvalrehihi_h are the highest doubles of the real parts
 *             of the function values on the host;
 *   funvalrelohi_h are the second highest doubles of the real parts
 *             of the function values on the host;
 *   funvalrehilo_h are the second lowest doubles of the real parts
 *             of the function values on the host;
 *   funvalrelolo_h are the lowest doubles of the real parts
 *             of the function values on the host;
 *   funvalimhihi_h are the highest doubles of the imaginary parts
 *             of the function values on the host;
 *   funvalimlohi_h are the second highest doubles of the imaginary parts
 *             of the function values on the host;
 *   funvalimhilo_h are the second lowest doubles of the imaginary parts
 *             of the function values on the host;
 *   funvalimlolo_h are the lowest doubles of the imaginary parts
 *             of the function values on the host;
 *   funvalrehihi_d are the highest doubles of the real parts
 *             of the function values on the device;
 *   funvalrelohi_d are the second highest doubles of the real parts
 *             of the function values on the device;
 *   funvalrehilo_d are the second lowest doubles of the real parts
 *             of the function values on the device;
 *   funvalrelolo_d are the lowest doubles of the real parts
 *             of the function values on the device;
 *   funvalimhihi_d are the highest doubles of the imaginary parts
 *             of the function values on the device;
 *   funvalimlohi_d are the second highest doubles of the imaginary parts
 *             of the function values on the device;
 *   funvalimhilo_d are the second lowest doubles of the imaginary parts
 *             of the function values on the device;
 *   funvalimlolo_d are the lowest doubles of the imaginary parts
 *             of the function values on the device;
 *   jacvalrehihi_h are the highest doubles of the real parts
 *             of the Jacobian matrix on the host;
 *   jacvalrelohi_h are the second highest doubles of the real parts
 *             of the Jacobian matrix on the host;
 *   jacvalrehilo_h are the second lowest doubles of the real parts
 *             of the Jacobian matrix on the host;
 *   jacvalrelolo_h are the lowest doubles of the real parts
 *             of the Jacobian matrix on the host;
 *   jacvalimhihi_h are the highest doubles of the imaginary parts
 *             of the Jacobian matrix on the host;
 *   jacvalimlohi_h are the second highest doubles of the imaginary parts
 *             of the Jacobian matrix on the host;
 *   jacvalimhilo_h are the second lowest doubles of the imaginary parts
 *             of the Jacobian matrix on the host;
 *   jacvalimlolo_h are the lowest doubles of the imaginary parts
 *             of the Jacobian matrix on the host;
 *   jacvalrehihi_d are the highest doubles of the real parts
 *             of the Jacobian matrix on the device;
 *   jacvalrelohi_d are the second highest doubles of the real parts
 *             of the Jacobian matrix on the device;
 *   jacvalrehilo_d are the second lowest doubles of the real parts
 *             of the Jacobian matrix on the device;
 *   jacvalrelolo_d are the lowest doubles of the real parts
 *             of the Jacobian matrix on the device;
 *   jacvalimhihi_d are the highest doubles of the imaginary parts
 *             of the Jacobian matrix on the device;
 *   jacvalimlohi_d are the second highest doubles of the imaginary parts
 *             of the Jacobian matrix on the device;
 *   jacvalimhilo_d are the second lowest doubles of the imaginary parts
 *             of the Jacobian matrix on the device;
 *   jacvalimlolo_d are the lowest doubles of the imaginary parts
 *             of the Jacobian matrix on the device;
 *   rhsrehihi_h are the highest doubles of real parts
 *             of the right hand side vector on the host;
 *   rhsrelohi_h are the second highest doubles of real parts
 *             of the right hand side vector on the host;
 *   rhsrehilo_h are the second lowest doubles of real parts
 *             of the right hand side vector on the host;
 *   rhsrelolo_h are the lowest doubles of real parts
 *             of the right hand side vector on the host;
 *   rhsimhihi_h are the highest doubles of imaginary parts
 *             of the right hand side vector on the host;
 *   rhsimlohi_h are the second highest doubles of imaginary parts
 *             of the right hand side vector on the host;
 *   rhsimhilo_h are the second lowest doubles of imaginary parts
 *             of the right hand side vector on the host;
 *   rhsimlolo_h are the lowest doubles of imaginary parts
 *             of the right hand side vector on the host;
 *   rhsrehihi_d are the highest doubles of real parts
 *             of the right hand side vector on the device;
 *   rhsrelohi_d are the second highest doubles of real parts
 *             of the right hand side vector on the device;
 *   rhsrehilo_d are the second lowest doubles of real parts
 *             of the right hand side vector on the device;
 *   rhsrelolo_d are the lowest doubles of real parts
 *             of the right hand side vector on the device;
 *   rhsimhihi_d are the highest doubles of imaginary parts
 *             of the right hand side vector on the device;
 *   rhsimlohi_d are the second highest doubles of imaginary parts
 *             of the right hand side vector on the device;
 *   rhsimhilo_d are the second lowest doubles of imaginary parts
 *             of the right hand side vector on the device;
 *   rhsimlolo_d are the lowest doubles of imaginary parts
 *             of the right hand side vector on the device;
 *   vrblvl    is the verbose level. */

int cmplx4_errors_inurhsQRsol
 ( int dim, int deg,
   double **inputrehihi_h, double **inputrelohi_h,
   double **inputrehilo_h, double **inputrelolo_h,
   double **inputimhihi_h, double **inputimlohi_h,
   double **inputimhilo_h, double **inputimlolo_h,
   double **inputrehihi_d, double **inputrelohi_d,
   double **inputrehilo_d, double **inputrelolo_d,
   double **inputimhihi_d, double **inputimlohi_d,
   double **inputimhilo_d, double **inputimlolo_d,
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
   double **solimhilo_d, double **solimlolo_d, int vrblvl );
/*
 * DESCRIPTION :
 *   Compares results on host with those on the device
 *   and returns the number of errors larger than 1.0e-20.
 *
 * ON ENTRY :
 *   dim       number of equations and variables;
 *   deg       degree at which the series are truncated;
 *   inputrehihi_h are the highest doubles of the real parts
 *             of the updated series on the host;
 *   inputrelohi_h are the second highest doubles of the real parts
 *             of the updated series on the host;
 *   inputrehilo_h are the second lowest doubles of the real parts
 *             of the updated series on the host;
 *   inputrelolo_h are the lowest doubles of the real parts
 *             of the updated series on the host;
 *   inputimhihi_h are the highest doubles of the imaginary parts
 *             of the updated series on the host;
 *   inputimlohi_h are the second highest doubles of the imaginary parts
 *             of the updated series on the host;
 *   inputimhilo_h are the second lowest doubles of the imaginary parts
 *             of the updated series on the host;
 *   inputimlolo_h are the lowest doubles of the imaginary parts
 *             of the updated series on the host;
 *   inputrehihi_d are the highest doubles of the real parts
 *             of the updated series on the device;
 *   inputrelohi_d are the second highest doubles of the real parts
 *             of the updated series on the device;
 *   inputrehilo_d are the second lowest doubles of the real parts
 *             of the updated series on the device;
 *   inputrelolo_d are the lowest doubles of the real parts
 *             of the updated series on the device;
 *   inputimhihi_d are the highest doubles of the imaginary parts
 *             of the updated series on the device;
 *   inputimhilo_d are the second highest doubles of the imaginary parts
 *             of the updated series on the device;
 *   inputimhilo_d are the second lowest doubles of the imaginary parts
 *             of the updated series on the device;
 *   inputimlolo_d are the lowest doubles of the imaginary parts
 *             of the updated series on the device;
 *   Qrehihi_h are the highest doubles of the real parts of Q on host;
 *   Qrelohi_h are the 2nd highest doubles of the real parts of Q on host;
 *   Qrehilo_h are the 2nd lowest doubles of the real parts of Q on host;
 *   Qrelolo_h are the lowest doubles of the real parts of Q on host;
 *   Qimhihi_h are the highest doubles of the imag parts of Q on host;
 *   Qimlohi_h are the 2nd highest doubles of the imag parts of Q on host;
 *   Qimhilo_h are the 2nd lowest doubles of the imag parts of Q on host;
 *   Qimlolo_h are the lowest doubles of the imag parts of Q on host;
 *   Qrehihi_d are the highest doubles of the real parts of Q on device;
 *   Qrelohi_d are the 2nd highest doubles of the real parts of Q on device;
 *   Qrehilo_d are the 2nd lowest doubles of the real parts of Q on device;
 *   Qrelolo_d are the lowest doubles of the real parts of Q on device;
 *   Qimhihi_d are the highest doubles of the imag parts of Q on device;
 *   Qimlohi_d are the 2nd highest doubles of the imag parts of Q on device;
 *   Qimhilo_d are the 2nd lowest doubles of the imag parts of Q on device;
 *   Qimlolo_d are the lowest doubles of the imag parts of Q on device;
 *   Rrehihi_h are the highest doubles of the real parts of R on host;
 *   Rrelohi_h are the 2nd highest doubles of the real parts of R on host;
 *   Rrehilo_h are the 2nd lowest doubles of the real parts of R on host;
 *   Rrelolo_h are the lowest doubles of the real parts of R on host;
 *   Rimhihi_h are the highest doubles of the imag parts of R on host;
 *   Rimlohi_h are the 2nd highest doubles of the imag parts of R on host;
 *   Rimhilo_h are the 2nd lowest doubles of the imag parts of R on host;
 *   Rimlolo_h are the lowest doubles of the imag parts of R on host;
 *   Rrehihi_d are the highest doubles of the real parts of R on device;
 *   Rrelohi_d are the 2nd highest doubles of the real parts of R on device;
 *   Rrehilo_d are the 2nd lowest doubles of the real parts of R on device;
 *   Rrelolo_d are the lowest doubles of the real parts of R on device;
 *   Rimhihi_d are the highest doubles of the imag parts of R on device;
 *   Rimlohi_d are the 2nd highest doubles of the imag parts of R on device;
 *   Rimhilo_d are the 2nd lowest doubles of the imag parts of R on device;
 *   Rimlolo_d are the lowest doubles of the imag parts of R on device;
 *   urhsrehihi_h are the highest doubles of the real parts
 *             of the updated right hand side on the host;
 *   urhsrelohi_h are the second highest doubles of the real parts
 *             of the updated right hand side on the host;
 *   urhsrehilo_h are the second lowest doubles of the real parts
 *             of the updated right hand side on the host;
 *   urhsrelolo_h are the lowest doubles of the real parts
 *             of the updated right hand side on the host;
 *   urhsimhihi_h are the highest doubles of the imaginary parts
 *             of the updated right hand side on the host;
 *   urhsimlohi_h are the second highest doubles of the imaginary parts
 *             of the updated right hand side on the host;
 *   urhsimhilo_h are the second lowest doubles of the imaginary parts
 *             of the updated right hand side on the host;
 *   urhsimlolo_h are the lowest doubles of the imaginary parts
 *             of the updated right hand side on the host;
 *   urhsrehihi_d are the highest doubles of the real parts
 *             of the updated right hand side on the device;
 *   urhsrelohi_d are the second highest doubles of the real parts
 *             of the updated right hand side on the device;
 *   urhsrehilo_d are the second lowest doubles of the real parts
 *             of the updated right hand side on the device;
 *   urhsrelolo_d are the lowest doubles of the real parts
 *             of the updated right hand side on the device;
 *   urhsimhihi_d are the highest doubles of the imaginary parts
 *             of the updated right hand side on the device;
 *   urhsimlohi_d are the second highest doubles of the imaginary parts
 *             of the updated right hand side on the device;
 *   urhsimhilo_d are the second lowest doubles of the imaginary parts
 *             of the updated right hand side on the device;
 *   urhsimlolo_d are the lowest doubles of the imaginary parts
 *             of the updated right hand side on the device;
 *   solrehihi_h are the highest doubles of the real parts
 *             of the update to the solution on the host,
 *   solrelohi_h are the second highest doubles of the real parts
 *             of the update to the solution on the host,
 *   solrehilo_h are the second lowest doubles of the real parts
 *             of the update to the solution on the host,
 *   solrelolo_h are the lowest doubles of the real parts
 *             of the update to the solution on the host,
 *   solimhihi_h are the highest doubles of the imaginary parts
 *             of the update to the solution on the host,
 *   solimlohi_h are the second highest doubles of the imaginary parts
 *             of the update to the solution on the host,
 *   solimhilo_h are the second lowest doubles of the imaginary parts
 *             of the update to the solution on the host,
 *   solimlolo_h are the lowest doubles of the imaginary parts
 *             of the update to the solution on the host,
 *   solrehihi_d are the highest doubles of the real parts
 *             of the update to the solution on the device;
 *   solrelohi_d are the second highest doubles of the real parts
 *             of the update to the solution on the device;
 *   solrehilo_d are the second lowest doubles of the real parts
 *             of the update to the solution on the device;
 *   solrelolo_d are the lowest doubles of the real parts
 *             of the update to the solution on the device;
 *   solimhihi_d are the highest doubles of the imaginary parts
 *             of the update to the solution on the device;
 *   solimlohi_d are the second highest doubles of the imaginary parts
 *             of the update to the solution on the device;
 *   solimhilo_d are the second lowest doubles of the imaginary parts
 *             of the update to the solution on the device;
 *   solimlolo_d are the lowest doubles of the imaginary parts
 *             of the update to the solution on the device;
 *   vrblvl    is the verbose level. */

int cmplx4_update_newton_qrstep
 ( int szt, int nbt, int dim, int deg,
   int *tailidx_h, int *tailidx_d,
   double **inputrehihi_h, double **inputrelohi_h,
   double **inputrehilo_h, double **inputrelolo_h,
   double **inputimhihi_h, double **inputimlohi_h,
   double **inputimhilo_h, double **inputimlolo_h,
   double **inputrehihi_d, double **inputrelohi_d,
   double **inputrehilo_d, double **inputrelolo_d,
   double **inputimhihi_d, double **inputimlohi_d,
   double **inputimhilo_d, double **inputimlolo_d,
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
   double *workvecrehihi, double *workvecrelohi,
   double *workvecrehilo, double *workvecrelolo,
   double *workvecimhihi, double *workvecimlohi,
   double *workvecimhilo, double *workvecimlolo,
   double **resvecrehihi, double **resvecrelohi,
   double **resvecrehilo, double **resvecrelolo,
   double **resvecimhihi, double **resvecimlohi,
   double **resvecimhilo, double **resvecimlolo,
   double *resmaxhihi, double *resmaxlohi,
   double *resmaxhilo, double *resmaxlolo,
   bool *zeroQ_h, bool *noqr_h, bool *zeroQ_d, bool *noqr_d,
   int *upidx_h, int *bsidx_h, int *upidx_d, int *bsidx_d,
   double *totqrlapsedms, double *totqtblapsedms,
   double *totbslapsedms, double *totupdlapsedms, double *totreslapsedms,
   int vrblvl, int mode );
/*
 * DESCRIPTION :
 *   Does one step with Newton's method to update a power series,
 *   using QR factorization to solve linear systems, on complex data,
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
 *   resmaxlolo is the lowest double of the max norm of the residual;
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

int cmplx4_column_newton_qrstep
 ( int szt, int nbt, int dim, int deg, int nbrcol,
   int *tailidx_h, int *tailidx_d,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac,
   double **mbrehihi, double **mbrelohi, double **mbrehilo, double **mbrelolo,
   double **mbimhihi, double **mbimlohi, double **mbimhilo, double **mbimlolo,
   double dpr,
   double ***cffrehihi, double ***cffrelohi,
   double ***cffrehilo, double ***cffrelolo,
   double ***cffimhihi, double ***cffimlohi,
   double ***cffimhilo, double ***cffimlolo,
   double **accrehihi, double **accrelohi,
   double **accrehilo, double **accrelolo,
   double **accimhihi, double **accimlohi,
   double **accimhilo, double **accimlolo,
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
   double *workvecrehihi, double *workvecrelohi,
   double *workvecrehilo, double *workvecrelolo,
   double *workvecimhihi, double *workvecimlohi,
   double *workvecimhilo, double *workvecimlolo,
   double **resvecrehihi, double **resvecrelohi,
   double **resvecrehilo, double **resvecrelolo,
   double **resvecimhihi, double **resvecimlohi,
   double **resvecimhilo, double **resvecimlolo,
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
 *   using QR factorization to solve linear systems, on complex data,
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
 *   mbrehihi  highest real parts of the right hand side series;
 *   mbrelohi  second highest real parts of the right hand side series;
 *   mbrehilo  second lowest real parts of the right hand side series;
 *   mbrelolo  lowest real parts of the right hand side series;
 *   mbimhihi  highest imaginary parts of the right hand side series;
 *   mbimlohi  second highest imaginary parts of the right hand side series;
 *   mbimhilo  second lowest imaginary parts of the right hand side series;
 *   mbimlolo  lowest imaginary parts of the right hand side series;
 *   dpr       damper multiplier for t, should be in (0.0, 1.0];
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
 *   resmaxlolo is the lowest double of the max norm of the residual;
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

int cmplx4_row_newton_qrstep
 ( int szt, int nbt, int dim, int deg, int *tailidx_h, int *tailidx_d,
   int *nbr, int **nvr, int ***idx,
   double **cstrehihi, double **cstrelohi,
   double **cstrehilo, double **cstrelolo,
   double **cstimhihi, double **cstimlohi,
   double **cstimhilo, double **cstimlolo,
   double ***cffrehihi, double ***cffrelohi,
   double ***cffrehilo, double ***cffrelolo,
   double ***cffimhihi, double ***cffimlohi,
   double ***cffimhilo, double ***cffimlolo,
   double dpr,
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
   double *workvecrehihi, double *workvecrelohi,
   double *workvecrehilo, double *workvecrelolo,
   double *workvecimhihi, double *workvecimlohi,
   double *workvecimhilo, double *workvecimlolo,
   double **resvecrehihi, double **resvecrelohi,
   double **resvecrehilo, double **resvecrelolo,
   double **resvecimhihi, double **resvecimlohi,
   double **resvecimhilo, double **resvecimlolo,
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
 *   using QR factorization to solve linear systems, on complex data,
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
 *   cstrehihi are the highest real parts of the constants;
 *   cstrelohi are the second highest real parts of the constants;
 *   cstrehilo are the second lowest real parts of the constants;
 *   cstrelolo are the lowest real parts of the constants;
 *   cstimhihi are the highest imaginary parts of the constants;
 *   cstimlohi are the second highest imaginary parts of the constants;
 *   cstimhilo are the second lowest imaginary parts of the constants;
 *   cstimlolo are the lowest imaginary parts of the constants;
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
 *   dpr       damper multiplier for t, should be in (0.0, 1.0];
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
 *   resmaxlolo is the lowest double of the max norm of the residual;
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

int cmplx4_allocate_inoutfunjac
 ( int dim, int deg, int mode,
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
   double ***jacvalimhilo_d, double ***jacvalimlolo_d );
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
 *   inputrehihi_h are the highest doubles of the real parts
 *              of the input on the host;
 *   inputrelohi_h are the second highest doubles of the real parts
 *              of the input on the host;
 *   inputrehilo_h are the second lowest doubles of the real parts
 *              of the input on the host;
 *   inputrelolo_h are the lowest doubles of the real parts
 *              of the input on the host;
 *   inputimhihi_h are the highest doubles of the imaginary parts
 *              of the input on the host;
 *   inputimlohi_h are the second highest doubles of the imaginary parts
 *              of the input on the host;
 *   inputimhilo_h are the second lowest doubles of the imaginary parts
 *              of the input on the host;
 *   inputimlolo_h are the lowest doubles of the imaginary parts
 *              of the input on the host;
 *   inputrehihi_d are the highest doubles of the real parts
 *              of the input on the device;
 *   inputrelohi_d are the second highest doubles of the real parts
 *              of the input on the device;
 *   inputrehilo_d are the second lowest doubles of the real parts
 *              of the input on the device;
 *   inputrelolo_d are the lowest doubles of the real parts
 *              of the input on the device;
 *   inputimhihi_d are the second highest doubles of the imaginary parts
 *              of the input on the device;
 *   inputimlohi_d are the second highest doubles of the imaginary parts
 *              of the input on the device;
 *   inputimlolo_d are the lowest doubles of the imaginary parts
 *              of the input on the device;
 *   inputimlolo_d are the lowest doubles of the imaginary parts
 *              of the input on the device;
 *   outputrehihi_h are the highest doubles of the real parts
 *              of the output on the host;
 *   outputrelohi_h are the second highest doubles of the real parts
 *              of the output on the host;
 *   outputrehilo_h are the second lowest doubles of the real parts
 *              of the output on the host;
 *   outputrelolo_h are the lowest doubles of the real parts
 *              of the output on the host;
 *   outputimhihi_h are the highest doubles of the imaginary parts
 *              of the output on the host;
 *   outputimlohi_h are the second highest doubles of the imaginary parts
 *              of the output on the host;
 *   outputimhilo_h are the second lowest doubles of the imaginary parts
 *              of the output on the host;
 *   outputimlolo_h are the lowest doubles of the imaginary parts
 *              of the output on the host;
 *   outputrehihi_d are the highest doubles of the real parts
 *              of the output on the device;
 *   outputrelohi_d are the second highest doubles of the real parts
 *              of the output on the device;
 *   outputrehilo_d are the second lowest doubles of the real parts
 *              of the output on the device;
 *   outputrelolo_d are the lowest doubles of the real parts
 *              of the output on the device;
 *   outputimhihi_d are the highest doubles of the imaginary parts 
 *              of the output on the device;
 *   outputimlohi_d are the second highest doubles of the imaginary parts 
 *              of the output on the device;
 *   outputimhilo_d are the second lowest doubles of the imaginary parts 
 *              of the output on the device;
 *   outputimlolo_d are the lowest doubles of the imaginary parts 
 *              of the output on the device;
 *   funvalrehihi_h are the highest doubles of the real parts
 *             of the function values on the host;
 *   funvalrelohi_h are the second highest doubles of the real parts
 *             of the function values on the host;
 *   funvalrehilo_h are the second lowest doubles of the real parts
 *             of the function values on the host;
 *   funvalrelolo_h are the lowest doubles of the real parts
 *             of the function values on the host;
 *   funvalimhihi_h are the highest doubles of the imaginary parts
 *             of the function values on the host;
 *   funvalimlohi_h are the second highest doubles of the imaginary parts
 *             of the function values on the host;
 *   funvalimhilo_h are the second lowest doubles of the imaginary parts
 *             of the function values on the host;
 *   funvalimlolo_h are the lowest doubles of the imaginary parts
 *             of the function values on the host;
 *   funvalrehihi_d are the highest doubles of the real parts
 *             of the function values on the device;
 *   funvalrelohi_d are the second highest doubles of the real parts
 *             of the function values on the device;
 *   funvalrehilo_d are the second lowest doubles of the real parts
 *             of the function values on the device;
 *   funvalrelolo_d are the lowest doubles of the real parts
 *             of the function values on the device;
 *   funvalimhihi_d are the highest doubles of the imaginary parts
 *             of the function values on the device;
 *   funvalimlohi_d are the second highest doubles of the imaginary parts
 *             of the function values on the device;
 *   funvalimhilo_d are the second lowest doubles of the imaginary parts
 *             of the function values on the device;
 *   funvalimlolo_d are the lowest doubles of the imaginary parts
 *             of the function values on the device;
 *   jacvalrehihi_h are the highest doubles of the real parts
 *             of the Jacobian matrix on the host;
 *   jacvalrelohi_h are the second highest doubles of the real parts
 *             of the Jacobian matrix on the host;
 *   jacvalrehilo_h are the second lowest doubles of the real parts
 *             of the Jacobian matrix on the host;
 *   jacvalrelolo_h are the lowest doubles of the real parts
 *             of the Jacobian matrix on the host;
 *   jacvalimhihi_h are the highest doubles of the imaginary parts
 *             of the Jacobian matrix on the host;
 *   jacvalimlohi_h are the second highest doubles of the imaginary parts
 *             of the Jacobian matrix on the host;
 *   jacvalimhilo_h are the second lowest doubles of the imaginary parts
 *             of the Jacobian matrix on the host;
 *   jacvalimlolo_h are the lowest doubles of the imaginary parts
 *             of the Jacobian matrix on the host;
 *   jacvalrehihi_d are the highest doubles of the real parts
 *             of the Jacobian matrix on the device;
 *   jacvalrelohi_d are the second highest doubles of the real parts
 *             of the Jacobian matrix on the device;
 *   jacvalrehilo_d are the second lowest doubles of the real parts
 *             of the Jacobian matrix on the device;
 *   jacvalrelolo_d are the lowest doubles of the real parts
 *             of the Jacobian matrix on the device;
 *   jacvalimhihi_d are the highest doubles of the imaginary parts
 *             of the Jacobian matrix on the device;
 *   jacvalimlohi_d are the second highest doubles of the imaginary parts
 *             of the Jacobian matrix on the device;
 *   jacvalimhilo_d are the second lowest doubles of the imaginary parts
 *             of the Jacobian matrix on the device;
 *   jacvalimlolo_d are the lowest doubles of the imaginary parts
 *             of the Jacobian matrix on the device. */

int cmplx4_allocate_rhsqrsol
 ( int dim, int deg, int mode,
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
   double **solrehihi_h, double **solrelohi_h,
   double **solrehilo_h, double **solrelolo_h,
   double **solimhihi_h, double **solimlohi_h,
   double **solimhilo_h, double **solimlolo_h,
   double **solrehihi_d, double **solrelohi_d,
   double **solrehilo_d, double **solrelolo_d,
   double **solimhihi_d, double **solimlohi_d,
   double **solimhilo_d, double **solimlolo_d );
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
 *   rhsrehihi_h are the highest doubles of real parts
 *             of the right hand side vector on the host;
 *   rhsrelohi_h are the second highest doubles of real parts
 *             of the right hand side vector on the host;
 *   rhsrehilo_h are the second lowest doubles of real parts
 *             of the right hand side vector on the host;
 *   rhsrelolo_h are the lowest doubles of real parts
 *             of the right hand side vector on the host;
 *   rhsimhihi_h are the highest doubles of imaginary parts
 *             of the right hand side vector on the host;
 *   rhsimlohi_h are the second highest doubles of imaginary parts
 *             of the right hand side vector on the host;
 *   rhsimhilo_h are the second lowest doubles of imaginary parts
 *             of the right hand side vector on the host;
 *   rhsimlolo_h are the lowest doubles of imaginary parts
 *             of the right hand side vector on the host;
 *   rhsrehihi_d are the highest doubles of real parts
 *             of the right hand side vector on the device;
 *   rhsrelohi_d are the second highest doubles of real parts
 *             of the right hand side vector on the device;
 *   rhsrehilo_d are the second lowest doubles of real parts
 *             of the right hand side vector on the device;
 *   rhsrelolo_d are the lowest doubles of real parts
 *             of the right hand side vector on the device;
 *   rhsimhihi_d are the highest doubles of imaginary parts
 *             of the right hand side vector on the device;
 *   rhsimlohi_d are the second highest doubles of imaginary parts
 *             of the right hand side vector on the device;
 *   rhsimhilo_d are the second lowest doubles of imaginary parts
 *             of the right hand side vector on the device;
 *   rhsimlolo_d are the lowest doubles of imaginary parts
 *             of the right hand side vector on the device;
 *   urhsrehihi_h are the highest doubles of the real parts
 *             of the updated right hand side on the host;
 *   urhsrelohi_h are the second highest doubles of the real parts
 *             of the updated right hand side on the host;
 *   urhsrehilo_h are the second lowest doubles of the real parts
 *             of the updated right hand side on the host;
 *   urhsrelolo_h are the lowest doubles of the real parts
 *             of the updated right hand side on the host;
 *   urhsimhihi_h are the highest doubles of the imaginary parts
 *             of the updated right hand side on the host;
 *   urhsimlohi_h are the second highest doubles of the imaginary parts
 *             of the updated right hand side on the host;
 *   urhsimhilo_h are the second lowest doubles of the imaginary parts
 *             of the updated right hand side on the host;
 *   urhsimlolo_h are the lowest doubles of the imaginary parts
 *             of the updated right hand side on the host;
 *   urhsrehihi_d are the highest doubles of the real parts
 *             of the updated right hand side on the device;
 *   urhsrelohi_d are the second highest doubles of the real parts
 *             of the updated right hand side on the device;
 *   urhsrehilo_d are the second lowest doubles of the real parts
 *             of the updated right hand side on the device;
 *   urhsrelolo_d are the lowest doubles of the real parts
 *             of the updated right hand side on the device;
 *   urhsimhihi_d are the highest doubles of the imaginary parts
 *             of the updated right hand side on the device;
 *   urhsimlohi_d are the second highest doubles of the imaginary parts
 *             of the updated right hand side on the device;
 *   urhsimhilo_d are the second lowest doubles of the imaginary parts
 *             of the updated right hand side on the device;
 *   urhsimlolo_d are the lowest doubles of the imaginary parts
 *             of the updated right hand side on the device;
 *   Qrehihi_h are the highest doubles of the real parts of Q on host;
 *   Qrelohi_h are the 2nd highest doubles of the real parts of Q on host;
 *   Qrehilo_h are the 2nd lowest doubles of the real parts of Q on host;
 *   Qrelolo_h are the lowest doubles of the real parts of Q on host;
 *   Qimhihi_h are the highest doubles of the imag parts of Q on host;
 *   Qimlohi_h are the 2nd highest doubles of the imag parts of Q on host;
 *   Qimhilo_h are the 2nd lowest doubles of the imag parts of Q on host;
 *   Qimlolo_h are the lowest doubles of the imag parts of Q on host;
 *   Qrehihi_d are the highest doubles of the real parts of Q on device;
 *   Qrelohi_d are the 2nd highest doubles of the real parts of Q on device;
 *   Qrehilo_d are the 2nd lowest doubles of the real parts of Q on device;
 *   Qrelolo_d are the lowest doubles of the real parts of Q on device;
 *   Qimhihi_d are the highest doubles of the imag parts of Q on device;
 *   Qimlohi_d are the 2nd highest doubles of the imag parts of Q on device;
 *   Qimhilo_d are the 2nd lowest doubles of the imag parts of Q on device;
 *   Qimlolo_d are the lowest doubles of the imag parts of Q on device;
 *   Rrehihi_h are the highest doubles of the real parts of R on host;
 *   Rrelohi_h are the 2nd highest doubles of the real parts of R on host;
 *   Rrehilo_h are the 2nd lowest doubles of the real parts of R on host;
 *   Rrelolo_h are the lowest doubles of the real parts of R on host;
 *   Rimhihi_h are the highest doubles of the imag parts of R on host;
 *   Rimlohi_h are the 2nd highest doubles of the imag parts of R on host;
 *   Rimhilo_h are the 2nd lowest doubles of the imag parts of R on host;
 *   Rimlolo_h are the lowest doubles of the imag parts of R on host;
 *   Rrehihi_d are the highest doubles of the real parts of R on device;
 *   Rrelohi_d are the 2nd highest doubles of the real parts of R on device;
 *   Rrehilo_d are the 2nd lowest doubles of the real parts of R on device;
 *   Rrelolo_d are the lowest doubles of the real parts of R on device;
 *   Rimhihi_d are the highest doubles of the imag parts of R on device;
 *   Rimlohi_d are the 2nd highest doubles of the imag parts of R on device;
 *   Rimhilo_d are the 2nd lowest doubles of the imag parts of R on device;
 *   Rimlolo_d are the lowest doubles of the imag parts of R on device;
 *   solrehihi_h are the highest doubles of the real parts
 *             of the update to the solution on the host,
 *   solrelohi_h are the second highest doubles of the real parts
 *             of the update to the solution on the host,
 *   solrehilo_h are the second lowest doubles of the real parts
 *             of the update to the solution on the host,
 *   solrelolo_h are the lowest doubles of the real parts
 *             of the update to the solution on the host,
 *   solimhihi_h are the highest doubles of the imaginary parts
 *             of the update to the solution on the host,
 *   solimlohi_h are the second highest doubles of the imaginary parts
 *             of the update to the solution on the host,
 *   solimhilo_h are the second lowest doubles of the imaginary parts
 *             of the update to the solution on the host,
 *   solimlolo_h are the lowest doubles of the imaginary parts
 *             of the update to the solution on the host,
 *   solrehihi_d are the highest doubles of the real parts
 *             of the update to the solution on the device;
 *   solrelohi_d are the second highest doubles of the real parts
 *             of the update to the solution on the device;
 *   solrehilo_d are the second lowest doubles of the real parts
 *             of the update to the solution on the device;
 *   solrelolo_d are the lowest doubles of the real parts
 *             of the update to the solution on the device;
 *   solimhihi_d are the highest doubles of the imaginary parts
 *             of the update to the solution on the device;
 *   solimlohi_d are the second highest doubles of the imaginary parts
 *             of the update to the solution on the device;
 *   solimhilo_d are the second lowest doubles of the imaginary parts
 *             of the update to the solution on the device;
 *   solimlolo_d are the lowest doubles of the imaginary parts
 *             of the update to the solution on the device. */

void cmplx4_start_setup
 ( int dim, int deg,
   double **testsolrehihi, double **testsolrelohi,
   double **testsolrehilo, double **testsolrelolo,
   double **testsolimhihi, double **testsolimlohi,
   double **testsolimhilo, double **testsolimlolo,
   double **inputrehihi_h, double **inputrelohi_h,
   double **inputrehilo_h, double **inputrelolo_h,
   double **inputimhihi_h, double **inputimlohi_h,
   double **inputimhilo_h, double **inputimlolo_h,
   double **inputrehihi_d, double **inputrelohi_d,
   double **inputrehilo_d, double **inputrelolo_d,
   double **inputimhihi_d, double **inputimlohi_d,
   double **inputimhilo_d, double **inputimlolo_d, int mode, int vrblvl );
/*
 * DESCRIPTION :
 *   Given the test solution, defines the start vector.
 *
 * ON ENTRY :
 *   dim        dimension of the solution vectors;
 *   deg        degree at which the series are truncated;
 *   testsolrehihi has the highest doubles of the real parts
 *              of the test solution;
 *   testsolrelohi has the second highest doubles of the real parts
 *              of the test solution;
 *   testsolrehilo has the second lowest doubles of the real parts
 *              of the test solution;
 *   testsolrelolo has the lowest doubles of the real parts
 *              of the test solution;
 *   testsolimhihi has the highest doubles of the imaginary parts
 *              of the test solution;
 *   testsolimlohi has the second highest doubles of the imaginary parts
 *              of the test solution;
 *   testsolimhilo has the second lowest doubles of the imaginary parts
 *              of the test solution;
 *   testsolimlolo has the lowest doubles of the imaginary parts
 *              of the test solution;
 *   inputrehihi_h is allocated on host if mode is 1 or 2;
 *   inputrelohi_h is allocated on host if mode is 1 or 2;
 *   inputrehilo_h is allocated on host if mode is 1 or 2;
 *   inputrelolo_h is allocated on host if mode is 1 or 2;
 *   inputimhihi_h is allocated on host if mode is 1 or 2;
 *   inputimlohi_h is allocated on host if mode is 1 or 2;
 *   inputimhilo_h is allocated on host if mode is 1 or 2;
 *   inputimlolo_h is allocated on host if mode is 1 or 2;
 *   inputrehihi_d is allocated on device if mode is 0 or 2;
 *   inputrelohi_d is allocated on device if mode is 0 or 2;
 *   inputrehilo_d is allocated on device if mode is 0 or 2;
 *   inputrelolo_d is allocated on device if mode is 0 or 2;
 *   inputimhihi_d is allocated on device if mode is 0 or 2;
 *   inputimlohi_d is allocated on device if mode is 0 or 2;
 *   inputimhilo_d is allocated on device if mode is 0 or 2;
 *   inputimlolo_d is allocated on device if mode is 0 or 2;
 *   mode       0 (GPU), 1 (CPU), or 2 (GPU+CPU);
 *   vrblvl     is the verbose level.
 *
 * ON RETURN :
 *   inputrehihi_h are the highest doubles of the real parts
 *              of the start vector on host if mode is 1 or 2;
 *   inputrelohi_h are the second highest doubles of the real parts
 *              of the start vector on host if mode is 1 or 2;
 *   inputrehilo_h are the second lowest doubles of the real parts
 *              of the start vector on host if mode is 1 or 2;
 *   inputrelolo_h are the lowest doubles of the real parts
 *              of the start vector on host if mode is 1 or 2;
 *   inputimhihi_h are the highest doubles of the imaginary parts
 *              of the start vector on host if mode is 1 or 2;
 *   inputimlohi_h are the second highest doubles of the imaginary parts
 *              of the start vector on host if mode is 1 or 2;
 *   inputimhilo_h are the second lowest doubles of the imaginary parts
 *              of the start vector on host if mode is 1 or 2;
 *   inputimlolo_h are the lowest doubles of the imaginary parts
 *              of the start vector on host if mode is 1 or 2;
 *   inputrehihi_d are the highest doubles of the real parts
 *              of the start vector on device if mode is 0 or 2;
 *   inputrelohi_d are the second highest doubles of the real parts
 *              of the start vector on device if mode is 0 or 2;
 *   inputrehilo_d are the second lowest doubles of the real parts
 *              of the start vector on device if mode is 0 or 2;
 *   inputrelolo_d are the lowest doubles of the real parts
 *              of the start vector on device if mode is 0 or 2;
 *   inputimhihi_d are the highest doubles of the imaginary parts
 *              of the start vector on device if mode is 0 or 2.
 *   inputimlohi_d are the second highest doubles of the imaginary parts
 *              of the start vector on device if mode is 0 or 2.
 *   inputimhilo_d are the second lowest doubles of the imaginary parts
 *              of the start vector on device if mode is 0 or 2;
 *   inputimlolo_d are the lowest doubles of the imaginary parts
 *              of the start vector on device if mode is 0 or 2. */

void cmplx4_column_setup
 ( int dim, int deg, int nbrcol,
   int **nvr, int ***idx, int **rowsA,
   double ***cffrehihi, double ***cffrelohi,
   double ***cffrehilo, double ***cffrelolo,
   double ***cffimhihi, double ***cffimlohi,
   double ***cffimhilo, double ***cffimlolo,
   double **testsolrehihi, double **testsolrelohi,
   double **testsolrehilo, double **testsolrelolo,
   double **testsolimhihi, double **testsolimlohi,
   double **testsolimhilo, double **testsolimlolo,
   double **mbrhsrehihi, double **mbrhsrelohi,
   double **mbrhsrehilo, double **mbrhsrelolo,
   double **mbrhsimhihi, double **mbrhsimlohi,
   double **mbrhsimhilo, double **mbrhsimlolo,
   double **inputrehihi_h, double **inputrelohi_h,
   double **inputrehilo_h, double **inputrelolo_h,
   double **inputimhihi_h, double **inputimlohi_h,
   double **inputimhilo_h, double **inputimlolo_h,
   double **inputrehihi_d, double **inputrelohi_d,
   double **inputrehilo_d, double **inputrelolo_d,
   double **inputimhihi_d, double **inputimlohi_d,
   double **inputimhilo_d, double **inputimlolo_d, int mode, int vrblvl );
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
 *   cffrehihi  are the highest doubles of the real parts of
 *              the coefficients, if more than one column;
 *   cffrelohi  are the second highest doubles of the real parts of
 *              the coefficients, if more than one column;
 *   cffrehilo  are the second lowest doubles of the real parts of
 *              the coefficients, if more than one column;
 *   cffrelolo  are the lowest doubles of the real parts of
 *              the coefficients, if more than one column;
 *   cffimhihi  are the highest doubles of the imaginary parts of
 *              the coefficients, if more than one column;
 *   cffimlohi  are the second highest doubles of the imaginary parts of
 *              the coefficients, if more than one column;
 *   cffimhilo  are the second lowest doubles of the imaginary parts of
 *              the coefficients, if more than one column;
 *   cffimlolo  are the lowest doubles of the imaginary parts of
 *              the coefficients, if more than one column;
 *   testsolrehihi has space for dim pointers;
 *   testsolrelohi has space for dim pointers;
 *   testsolrehilo has space for dim pointers;
 *   testsolrelolo has space for dim pointers;
 *   testsolimhihi has space for dim pointers;
 *   testsolimlohi has space for dim pointers;
 *   testsolimhilo has space for dim pointers;
 *   testsolimlolo has space for dim pointers;
 *   mbrhsrehihi has space for dim pointers;
 *   mbrhsrelohi has space for dim pointers;
 *   mbrhsrehilo has space for dim pointers;
 *   mbrhsrelolo has space for dim pointers;
 *   mbrhsimhihi has space for dim pointers;
 *   mbrhsimlohi has space for dim pointers;
 *   mbrhsimhilo has space for dim pointers;
 *   mbrhsimlolo has space for dim pointers;
 *   inputrehihi_h is allocated on host if mode is 1 or 2;
 *   inputrelohi_h is allocated on host if mode is 1 or 2;
 *   inputrehilo_h is allocated on host if mode is 1 or 2;
 *   inputrelolo_h is allocated on host if mode is 1 or 2;
 *   inputimhihi_h is allocated on host if mode is 1 or 2;
 *   inputimlohi_h is allocated on host if mode is 1 or 2;
 *   inputimhilo_h is allocated on host if mode is 1 or 2;
 *   inputimlolo_h is allocated on host if mode is 1 or 2;
 *   inputrehihi_d is allocated on device if mode is 0 or 2;
 *   inputrelohi_d is allocated on device if mode is 0 or 2;
 *   inputrehilo_d is allocated on device if mode is 0 or 2;
 *   inputrelolo_d is allocated on device if mode is 0 or 2;
 *   inputimhihi_d is allocated on device if mode is 0 or 2;
 *   inputimlohi_d is allocated on device if mode is 0 or 2;
 *   inputimhilo_d is allocated on device if mode is 0 or 2;
 *   inputimlolo_d is allocated on device if mode is 0 or 2;
 *   mode       0 (GPU), 1 (CPU), or 2 (GPU+CPU);
 *   vrblvl     is the verbose level.
 *
 * ON RETURN :
 *   testsolrehihi has the highest doubles of the real parts
 *              of the test solution;
 *   testsolrelohi has the second highest doubles of the real parts
 *              of the test solution;
 *   testsolrehilo has the second lowest doubles of the real parts
 *              of the test solution;
 *   testsolrelolo has the lowest doubles of the real parts
 *              of the test solution;
 *   testsolimhihi has the highest doubles of the imaginary parts
 *              of the test solution;
 *   testsolimlohi has the second highest doubles of the imaginary parts
 *              of the test solution;
 *   testsolimhilo has the second lowest doubles of the imaginary parts
 *              of the test solution;
 *   testsolimlolo has the lowest doubles of the imaginary parts
 *              of the test solution;
 *   mbrhsrelohi has the highest doubles of the real parts
 *              of right hand side vector for the test solution;
 *   mbrhsrehilo has the lowest doubles of the real parts
 *              of right hand side vector for the test solution;
 *   mbrhsrelolo has the lowest doubles of the real parts
 *              of right hand side vector for the test solution;
 *   mbrhsimhihi has the highest doubles of the imag parts
 *              of right hand side vector for the test solution;
 *   mbrhsimlohi has the highest doubles of the imag parts
 *              of right hand side vector for the test solution;
 *   mbrhsimhilo has the lowest doubles of the imag parts
 *              of right hand side vector for the test solution;
 *   mbrhsimlolo has the lowest doubles of the imag parts
 *              of right hand side vector for the test solution;
 *   inputrehihi_h are the highest doubles of the real parts
 *              of the start vector on host if mode is 1 or 2;
 *   inputrelohi_h are the second highest doubles of the real parts
 *              of the start vector on host if mode is 1 or 2;
 *   inputrehilo_h are the second lowest doubles of the real parts
 *              of the start vector on host if mode is 1 or 2;
 *   inputrelolo_h are the lowest doubles of the real parts
 *              of the start vector on host if mode is 1 or 2;
 *   inputimhihi_h are the highest doubles of the imaginary parts
 *              of the start vector on host if mode is 1 or 2;
 *   inputimlohi_h are the second highest doubles of the imaginary parts
 *              of the start vector on host if mode is 1 or 2;
 *   inputimhilo_h are the second lowest doubles of the imaginary parts
 *              of the start vector on host if mode is 1 or 2;
 *   inputimlolo_h are the lowest doubles of the imaginary parts
 *              of the start vector on host if mode is 1 or 2;
 *   inputrehihi_d are the highest doubles of the real parts
 *              of the start vector on device if mode is 0 or 2;
 *   inputrelohi_d are the second highest doubles of the real parts
 *              of the start vector on device if mode is 0 or 2;
 *   inputrehilo_d are the second lowest doubles of the real parts
 *              of the start vector on device if mode is 0 or 2;
 *   inputrelolo_d are the lowest doubles of the real parts
 *              of the start vector on device if mode is 0 or 2;
 *   inputimhihi_d are the highest doubles of the imaginary parts
 *              of the start vector on device if mode is 0 or 2.
 *   inputimlohi_d are the second highest doubles of the imaginary parts
 *              of the start vector on device if mode is 0 or 2.
 *   inputimhilo_d are the second lowest doubles of the imaginary parts
 *              of the start vector on device if mode is 0 or 2;
 *   inputimlolo_d are the lowest doubles of the imaginary parts
 *              of the start vector on device if mode is 0 or 2. */

void cmplx4_row_setup
 ( int dim, int deg, int *nbr, int **nvr, int ***idx,
   double **cstrehihi, double **cstrelohi, 
   double **cstrehilo, double **cstrelolo,
   double **cstimhihi, double **cstimlohi,
   double **cstimhilo, double **cstimlolo,
   double ***cffrehihi, double ***cffrelohi,
   double ***cffrehilo, double ***cffrelolo,
   double ***cffimhihi, double ***cffimlohi,
   double ***cffimhilo, double ***cffimlolo,
   double **testsolrehihi, double **testsolrelohi,
   double **testsolrehilo, double **testsolrelolo,
   double **testsolimhihi, double **testsolimlohi,
   double **testsolimhilo, double **testsolimlolo,
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
   double ***outputimhilo_d, double ***outputimlolo_d, int mode, int vrblvl );
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
 *   cstrehihi  are the highest real parts of the constants;
 *   cstrelohi  are the second highest real parts of the constants;
 *   cstrehilo  are the second lowest real parts of the constants;
 *   cstrelolo  are the lowest real parts of the constants;
 *   cstimhihi  are the highest imaginary parts of the constants;
 *   cstimlohi  are the second highest imaginary parts of the constants;
 *   cstimhilo  are the second lowest imaginary parts of the constants;
 *   cstimlolo  are the lowest imaginary parts of the constants;
 *   cffrehihi  are the highest doubles of the real parts of
 *              the coefficients, if more than one column;
 *   cffrelohi  are the second highest doubles of the real parts of
 *              the coefficients, if more than one column;
 *   cffrehilo  are the second lowest doubles of the real parts of
 *              the coefficients, if more than one column;
 *   cffrelolo  are the lowest doubles of the real parts of
 *              the coefficients, if more than one column;
 *   cffimhihi  are the highest doubles of the imaginary parts of
 *              the coefficients, if more than one column;
 *   cffimlohi  are the second highest doubles of the imaginary parts of
 *              the coefficients, if more than one column;
 *   cffimhilo  are the second lowest doubles of the imaginary parts of
 *              the coefficients, if more than one column;
 *   cffimlolo  are the lowest doubles of the imaginary parts of
 *              the coefficients, if more than one column;
 *   testsolrehihi has space for dim pointers;
 *   testsolrelohi has space for dim pointers;
 *   testsolrehilo has space for dim pointers;
 *   testsolrelolo has space for dim pointers;
 *   testsolimhihi has space for dim pointers;
 *   testsolimlohi has space for dim pointers;
 *   testsolimhilo has space for dim pointers;
 *   testsolimlolo has space for dim pointers;
 *   inputrehihi_h is allocated on host if mode is 1 or 2;
 *   inputrelohi_h is allocated on host if mode is 1 or 2;
 *   inputrehilo_h is allocated on host if mode is 1 or 2;
 *   inputrelolo_h is allocated on host if mode is 1 or 2;
 *   inputimhihi_h is allocated on host if mode is 1 or 2;
 *   inputimlohi_h is allocated on host if mode is 1 or 2;
 *   inputimhilo_h is allocated on host if mode is 1 or 2;
 *   inputimlolo_h is allocated on host if mode is 1 or 2;
 *   inputrehihi_d is allocated on device if mode is 0 or 2;
 *   inputrelohi_d is allocated on device if mode is 0 or 2;
 *   inputrehilo_d is allocated on device if mode is 0 or 2;
 *   inputrelolo_d is allocated on device if mode is 0 or 2;
 *   inputimhihi_d is allocated on device if mode is 0 or 2;
 *   inputimlohi_d is allocated on device if mode is 0 or 2;
 *   inputimhilo_d is allocated on device if mode is 0 or 2;
 *   inputimlolo_d is allocated on device if mode is 0 or 2;
 *   mode       0 (GPU), 1 (CPU), or 2 (GPU+CPU);
 *   vrblvl     is the verbose level.
 *
 * ON RETURN :
 *   cstrehihi  highest doubles of the real parts of the constant,
 *              adjusted for the test solution;
 *   cstrelohi  second highest doubles of the real parts of the constant,
 *              adjusted for the test solution;
 *   cstrehilo  second lowest doubles of the real parts of the constant,
 *              adjusted for the test solution;
 *   cstrelolo  lowest doubles of the real parts of the constant,
 *              adjusted for the test solution;
 *   cstimhihi  highest doubles of the imag parts of the constant,
 *              adjusted for the test solution;
 *   cstimlohi  second highest doubles of the imag parts of the constant,
 *              adjusted for the test solution;
 *   cstimhilo  second lowest doubles of the imag parts of the constant,
 *              adjusted for the test solution;
 *   cstimlolo  lowest doubles of the imag parts of the constant,
 *              adjusted for the test solution;
 *   testsolrehihi has the highest doubles of the real parts
 *              of the test solution;
 *   testsolrelohi has the second highest doubles of the real parts
 *              of the test solution;
 *   testsolrehilo has the second lowest doubles of the real parts
 *              of the test solution;
 *   testsolrelolo has the lowest doubles of the real parts
 *              of the test solution;
 *   testsolimhihi has the highest doubles of the imaginary parts
 *              of the test solution;
 *   testsolimlohi has the second highest doubles of the imaginary parts
 *              of the test solution;
 *   testsolimhilo has the second lowest doubles of the imaginary parts
 *              of the test solution;
 *   testsolimlolo has the lowest doubles of the imaginary parts
 *              of the test solution;
 *   inputrehihi_h has the highest doubles of the real parts
 *              of the start vector for host if mode is 1 or 2;
 *   inputrelohi_h has the second highest doubles of the real parts
 *              of the start vector for host if mode is 1 or 2;
 *   inputrehilo_h has the second lowest doubles of the real parts
 *              of the start vector for host if mode is 1 or 2;
 *   inputrelolo_h has the lowest doubles of the real parts
 *              of the start vector for host if mode is 1 or 2;
 *   inputimhihi_h has the highest doubles of the imaginary parts
 *              of the start vector for host if mode is 1 or 2;
 *   inputimlohi_h has the second highest doubles of the imaginary parts
 *              of the start vector for host if mode is 1 or 2;
 *   inputimhilo_h has the second lowest doubles of the imaginary parts
 *              of the start vector for host if mode is 1 or 2;
 *   inputimlolo_h has the lowest doubles of the imaginary parts
 *              of the start vector for host if mode is 1 or 2;
 *   inputrehihi_d has the highest doubles of the real parts
 *              of the start vector for device if mode is 0 or 2;
 *   inputrelohi_d has the second highest doubles of the real parts
 *              of the start vector for device if mode is 0 or 2;
 *   inputrehilo_d has the second lowest doubles of the real parts
 *              of the start vector for device if mode is 0 or 2;
 *   inputrelolo_d has the lowest doubles of the real parts
 *              of the start vector for device if mode is 0 or 2;
 *   inputimhihi_d has the highest doubles of the imaginary parts
 *              of the start vector for device if mode is 0 or 2;
 *   inputimlohi_d has the second highest doubles of the imaginary parts
 *              of the start vector for device if mode is 0 or 2;
 *   inputimhilo_d has the second lowest doubles of the imaginary parts
 *              of the start vector for device if mode is 0 or 2;
 *   inputimlolo_d has the lowest doubles of the imaginary parts
 *              of the start vector for device if mode is 0 or 2;
 *   outputrehihi_h are highest doubles of the real parts
 *              of the evaluated test solution if mode is 1;
 *   outputrelohi_h are second highest doubles of the real parts
 *              of the evaluated test solution if mode is 1;
 *   outputrehilo_h are second lowest doubles of the real parts
 *              of the evaluated test solution if mode is 1;
 *   outputrelolo_h are lowest doubles of the real parts
 *              of the evaluated test solution if mode is 1;
 *   outputimhihi_h are the highest doubles of the imaginary parts
 *              of the evaluated test solution if mode is 1;
 *   outputimlohi_h are the second highest doubles of the imaginary parts
 *              of the evaluated test solution if mode is 1;
 *   outputimhilo_h are the second lowest doubles of the imaginary parts
 *              of the evaluated test solution if mode is 1;
 *   outputimlolo_h are the lowest doubles of the imaginary parts
 *              of the evaluated test solution if mode is 1;
 *   outputrehihi_d are the highest doubles of the real parts
 *              of the evaluated test solution, if mode is 0 or 2;
 *   outputrelohi_d are the second highest doubles of the real parts
 *              of the evaluated test solution, if mode is 0 or 2;
 *   outputrehilo_d are the second lowest doubles of the real parts
 *              of the evaluated test solution, if mode is 0 or 2;
 *   outputrelolo_d are the lowest doubles of the real parts
 *              of the evaluated test solution, if mode is 0 or 2;
 *   outputimhihi_d are the highest doubles of the imaginary parts
 *              of the evaluated test solution, if mode is 0 or 2;
 *   outputimlohi_d are the second highest doubles of the imaginary parts
 *              of the evaluated test solution, if mode is 0 or 2;
 *   outputimhilo_d are the second lowest doubles of the imaginary parts
 *              of the evaluated test solution, if mode is 0 or 2;
 *   outputimlolo_d are the lowest doubles of the imaginary parts
 *              of the evaluated test solution, if mode is 0 or 2. */

int cmplx4_error_testsol
 ( int dim, int deg, int mode,
   double **testsolrehihi, double **testsolrelohi,
   double **testsolrehilo, double **testsolrelolo,
   double **testsolimhihi, double **testsolimlohi,
   double **testsolimhilo, double **testsolimlolo,
   double **inputrehihi_h, double **inputrelohi_h,
   double **inputrehilo_h, double **inputrelolo_h,
   double **inputimhihi_h, double **inputimlohi_h,
   double **inputimhilo_h, double **inputimlolo_h,
   double **inputrehihi_d, double **inputrelohi_d,
   double **inputrehilo_d, double **inputrelolo_d,
   double **inputimhihi_d, double **inputimlohi_d,
   double **inputimhilo_d, double **inputimlolo_d );
/*
 * DESCRIPTION :
 *   Compares the computed solution against the test solution.
 *
 * ON ENTRY :
 *   dim        dimension of the solution vectors;
 *   deg        degree at which the series are truncated;
 *   mode       0 (GPU), 1 (CPU), or 2 (GPU+CPU).
 *   testsolrehihi has the highest doubles of the real parts
 *              of the test solution;
 *   testsolrelohi has the second highest doubles of the real parts
 *              of the test solution;
 *   testsolrehilo has the second lowest doubles of the real parts
 *              of the test solution;
 *   testsolrelolo has the lowest doubles of the real parts
 *              of the test solution;
 *   testsolimhihi has the highest doubles of the imaginary parts
 *              of the test solution;
 *   testsolimlohi has the second highest doubles of the imaginary parts
 *              of the test solution;
 *   testsolimhilo has the second lowest doubles of the imaginary parts
 *              of the test solution;
 *   testsolimlolo has the lowest doubles of the imaginary parts
 *              of the test solution;
 *   inputrehihi_h has the highest doubles of the real parts
 *              of the solution on the host if mode is 1 or 2;
 *   inputrelohi_h has the second highest doubles of the real parts
 *              of the solution on the host if mode is 1 or 2;
 *   inputrehilo_h has the second lowest doubles of the real parts
 *              of the solution on the host if mode is 1 or 2;
 *   inputrelolo_h has the lowest doubles of the real parts
 *              of the solution on the host if mode is 1 or 2;
 *   inputimhihi_h has the highest doubles of the imaginary parts
 *              of the solution on the host if mode is 1 or 2;
 *   inputimlohi_h has the second highest doubles of the imaginary parts
 *              of the solution on the host if mode is 1 or 2;
 *   inputimhilo_h has the second lowest doubles of the imaginary parts
 *              of the solution on the host if mode is 1 or 2;
 *   inputimlolo_h has the lowest doubles of the imaginary parts
 *              of the solution on the host if mode is 1 or 2;
 *   inputrehihi_d has the highest doubles of the real parts
 *              of the solution on the device if mode is 0 or 2;
 *   inputrelohi_d has the second highest doubles of the real parts
 *              of the solution on the device if mode is 0 or 2;
 *   inputrehilo_d has the second lowest doubles of the real parts
 *              of the solution on the device if mode is 0 or 2;
 *   inputrelolo_d has the lowest doubles of the real parts
 *              of the solution on the device if mode is 0 or 2;
 *   inputimhihi_d has the highest doubles of the imaginary parts
 *              of the solution on the device if mode is 0 or 2;
 *   inputimlohi_d has the second highest doubles of the imaginary parts
 *              of the solution on the device if mode is 0 or 2;
 *   inputimhilo_d has the second lowest doubles of the imaginary parts
 *              of the solution on the evice if mode is 0 or 2;
 *   inputimlolo_d has the lowest doubles of the imaginary parts
 *              of the solution on the evice if mode is 0 or 2. */

int test_cmplx4_column_newton
 ( int szt, int nbt, int dim, int deg, int nbrcol,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac, int **rowsA,
   double dpr, int nbsteps, int mode, int vrblvl );
/*
 * DESCRIPTION :
 *   Runs Newton with complex quad double arithmetic,
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
 *   rowsA     rows of exponents of the dim monomials;
 *   dpr       damper multiplier for t, should be in (0.0, 1.0];
 *   nbsteps   the number of Newton steps;
 *   mode      the mode of execution, 0 for GPU only, 1 for CPU only,
 *             2 for CPU and GPU;
 *   vrblvl    is the verbose level. */

int test_cmplx4_row_newton
 ( int szt, int nbt, int dim, int deg, int *nbr, int **nvr, int ***idx,
   double dpr, int nbsteps, int mode, int vrblvl );
/*
 * DESCRIPTION :
 *   Runs Newton with complex quad double arithmetic,
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
