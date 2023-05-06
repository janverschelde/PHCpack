// The file cmplx4_path_tracker.h specifies a path tracker
// on series in quad double precision on complex numbers.

#ifndef __cmplx4_path_tracker_h__
#define __cmplx4_path_tracker_h__

int cmplx4_run_newton
 ( int szt, int nbt, int dim, int deg, int nbrcol, int nbsteps,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac,
   double **mbrehihi, double **mbrelohi, double **mbrehilo, double **mbrelolo,
   double **mbimhihi, double **mbimlohi, double **mbimhilo, double **mbimlolo,
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
   int vrblvl, int mode );
/*
 * DESCRIPTION :
 *   Runs Newton's method on power series on complex data
 *   in quad double precision.
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
 *   mbrehihi  highest real parts of the right hand side series;
 *   mbrelohi  second highest real parts of the right hand side series;
 *   mbrehilo  second lowest real parts of the right hand side series;
 *   mbrelolo  lowest real parts of the right hand side series;
 *   mbimhihi  highest imaginary parts of the right hand side series;
 *   mbimlohi  second highest imaginary parts of the right hand side series;
 *   mbimhilo  second lowest imaginary parts of the right hand side series;
 *   mbimlolo  lowest imaginary parts of the right hand side series;
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

int test_dbl4_complex_track
 ( int szt, int nbt, int dim, int deg, int nbrcol,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac, int **rowsA,
   int nbsteps, int mode, int vrblvl );
/*
 * DESCRIPTION :
 *   Tracks a path on a system with complex quad double arithmetic.
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
