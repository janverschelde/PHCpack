// The file cmplx8_path_tracker.h specifies a path tracker
// on series in octo double precision on complex numbers.

#ifndef __cmplx8_path_tracker_h__
#define __cmplx8_path_tracker_h__

int cmplx8_run_newton
 ( int szt, int nbt, int dim, int deg, int nbrcol, int nbsteps,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac,
   double **mbrehihihi, double **mbrelohihi,
   double **mbrehilohi, double **mbrelolohi,
   double **mbrehihilo, double **mbrelohilo,
   double **mbrehilolo, double **mbrelololo,
   double **mbimhihihi, double **mbimlohihi,
   double **mbimhilohi, double **mbimlolohi,
   double **mbimhihilo, double **mbimlohilo,
   double **mbimhilolo, double **mbimlololo,
   double ***cffrehihihi, double ***cffrelohihi,
   double ***cffrehilohi, double ***cffrelolohi,
   double ***cffrehihilo, double ***cffrelohilo,
   double ***cffrehilolo, double ***cffrelololo,
   double ***cffimhihihi, double ***cffimlohihi,
   double ***cffimhilohi, double ***cffimlolohi,
   double ***cffimhihilo, double ***cffimlohilo,
   double ***cffimhilolo, double ***cffimlololo,
   double **accrehihihi, double **accrelohihi,
   double **accrehilohi, double **accrelolohi,
   double **accrehihilo, double **accrelohilo,
   double **accrehilolo, double **accrelololo,
   double **accimhihihi, double **accimlohihi,
   double **accimhilohi, double **accimlolohi,
   double **accimhihilo, double **accimlohilo,
   double **accimhilolo, double **accimlololo,
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
   double *resmaxhilolo, double *resmaxlololo,
   int vrblvl, int mode );
/*
 * DESCRIPTION :
 *   Runs Newton's method on power series on complex data
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
 *   mbrehihihi are the highest real parts of the right hand side series;
 *   mbrelohihi are the 2nd highest real parts of the right hand side series;
 *   mbrehilohi are the 3rd highest real parts of the right hand side series;
 *   mbrelolohi are the 4th highest real parts of the right hand side series;
 *   mbrehihilo are the 4th lowest real parts of the right hand side series;
 *   mbrelohilo are the 3rd lowest real parts of the right hand side series;
 *   mbrehilolo are the 2nd lowest real parts of the right hand side series;
 *   mbrelololo are the lowest real parts of the right hand side series;
 *   mbimhihihi are the highest imaginary parts of the right hand side series;
 *   mbimlohihi are the 2nd highest imag parts of the right hand side series;
 *   mbimhilohi are the 3rd highest imag parts of the right hand side series;
 *   mbimlolohi are the 4th highest imag parts of the right hand side series;
 *   mbimhihilo are the 4th lowest imag parts of the right hand side series;
 *   mbimlohilo are the 3rd lowest imag parts of the right hand side series;
 *   mbimhilolo are the 2nd lowest imag parts of the right hand side series;
 *   mbimlololo are the lowest imaginary parts of the right hand side series;
 *   dpr       damper multiplier for t, should be in (0.0, 1.0];
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
 *   workvecrehihihi is work space allocated for a vector of dimension dim;
 *   workvecrelohihi is work space allocated for a vector of dimension dim;
 *   workvecrehilohi is work space allocated for a vector of dimension dim;
 *   workvecrelolohi is work space allocated for a vector of dimension dim;
 *   workvecrehihilo is work space allocated for a vector of dimension dim;
 *   workvecrelohilo is work space allocated for a vector of dimension dim;
 *   workvecrehilolo is work space allocated for a vector of dimension dim;
 *   workvecrelololo is work space allocated for a vector of dimension dim;
 *   workvecimhihihi is work space allocated for a vector of dimension dim;
 *   workvecimlohihi is work space allocated for a vector of dimension dim;
 *   workvecimhilohi is work space allocated for a vector of dimension dim;
 *   workvecimlolohi is work space allocated for a vector of dimension dim;
 *   workvecimhihilo is work space allocated for a vector of dimension dim;
 *   workvecimlohilo is work space allocated for a vector of dimension dim;
 *   workvecimhilolo is work space allocated for a vector of dimension dim;
 *   workvecimlololo is work space allocated for a vector of dimension dim;
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

int test_dbl8_complex_track
 ( int szt, int nbt, int dim, int deg, int nbrcol,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac, int **rowsA,
   int nbsteps, int mode, int vrblvl );
/*
 * DESCRIPTION :
 *   Tracks a path on a system with complex octo double arithmetic.
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
