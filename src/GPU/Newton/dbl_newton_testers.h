// The file dbl_newton_testers.h specifies test functions for Newton's method
// on series in double precision.

#ifndef __dbl_newton_testers_h__
#define __dbl_newton_testers_h__

void real_start_series_vector ( int dim, int deg, double *r0, double **cff );
/*
 * DESCRIPTION :
 *   Given space in cff for the vector of dim power series,
 *   of series truncated after degree deg,
 *   sets the coefficients cff of the start series.
 *   The constant term equals r0[j] with a small error. */

void cmplx_start_series_vector
 ( int dim, int deg, double *r0re, double *r0im,
   double **cffre, double **cffim );
/*
 * DESCRIPTION :
 *   Given space allocated in cffre and cffim,
 *   sets the coefficients of the start series.
 *   The constant term equals r0re[j] + i*r0im[j] with a small error,
 *   for j ranging from 0 to dim-1. */

void dbl_unit_series_vector ( int dim, int deg, double **cff );
/*
 * DESCRIPTION :
 *   Given space in cff for the vector of dim power series,
 *   of series truncated after degree deg,
 *   initializes the coefficients in cff to one as leading coefficients,
 *   and zero for all other coefficients. */

void cmplx_unit_series_vector
 ( int dim, int deg, double **cffre, double **cffim );
/*
 * DESCRIPTION :
 *   Given space allocated in cffre and cffim,
 *   returns the values for a complex unit series. */

void dbl_unit_series_vectors ( int nbr, int dim, int deg, double ***cff );
/*
 * DESCRIPTION :
 *   Given space in cff for nbr vectors of dim power series,
 *   of series truncated after degree deg,
 *   initializes the coefficients in cff to one as leading coefficients,
 *   and zero for all other coefficients. */

void cmplx_unit_series_vectors
 ( int nbr, int dim, int deg, double ***cffre, double ***cffim );
/*
 * DESCRIPTION :
 *   Given space allocated in cffre and cffim for nbr series,
 *   returns the values for a complex unit series. */

void dbl_update_series
 ( int dim, int degp1, int startidx, double **x, double **dx, int vrblvl );
/*
 * DESCRIPTION :
 *   Adds the series in dx to x.
 *
 * ON ENTRY :
 *   dim       number of series in x and dx;
 *   degp1     degree plus one of the series;
 *   startidx  index of the start of the update;
 *   x         series to updated, not linearized;
 *   dx        linearized update of the series;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   x         the series x updated with dx. */

void cmplx_update_series
 ( int dim, int degp1, int startidx,
   double **xre, double **xim, double **dxre, double **dxim, int vrblvl );
/*
 * DESCRIPTION :
 *   Adds the series in dx to x.
 *
 * ON ENTRY :
 *   dim       number of series in x and dx;
 *   degp1     degree plus one of the series;
 *   startidx  index of the start of the update;
 *   xre       real parts of series to updated, not linearized;
 *   xim       imaginary parts of series to updated, not linearized;
 *   dxre      real parts of the linearized update of the series;
 *   dxim      imaginary parts of the linearized update of the series;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   xre       real parts of the series x updated with dx;
 *   xim       imaginary parts of the series x updated with dx. */

double dbl_error3sum
 ( int dim1, int dim2, int dim3,
   double ***data_h, double ***data_d,
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
 *   data_h    data computed on the host;
 *   data_d    data computed on the device;
 *   banner    string for printing if verbose level vrblvl > 1;
 *   vrblvl    is the verbose level. */

double cmplx_error3sum
 ( int dim1, int dim2, int dim3,
   double ***datare_h, double ***dataim_h,
   double ***datare_d, double ***dataim_d,
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
 *   datare_h  real parts of data computed on the host;
 *   dataim_h  imaginary parts of data computed on the host;
 *   datare_d  real parts of data computed on the device;
 *   dataim_d  imaginary parts of data computed on the device;
 *   banner    string for printing if verbose level vrblvl > 1;
 *   vrblvl    is the verbose level. */

double dbl_error2sum
 ( int nrows, int ncols, double **data_h, double **data_d,
   std::string banner, int vrblvl );
/*
 * DESCRIPTION :
 *   Returns the sum of the absolute values of the differences
 *   between host and device real data, in two dimensions.
 *
 * ON ENTRY :
 *   nrows     number of rows of the matrices;
 *   ncols     number of columns of the matrices;
 *   data_h    data computed on the host;
 *   data_d    data computed on the device;
 *   banner    string for printing if verbose level vrblvl > 1;
 *   vrblvl    is the verbose level. */

double cmplx_error2sum
 ( int nrows, int ncols,
   double **datare_h, double **dataim_h,
   double **datare_d, double **dataim_d,
   std::string banner, int vrblvl );
/*
 * DESCRIPTION :
 *   Returns the sum of the absolute values of the differences
 *   between host and device complex data, in two dimensions.
 *
 * ON ENTRY :
 *   nrows     number of rows of the matrices;
 *   ncols     number of columns of the matrices;
 *   datare_h  real parts of data computed on the host;
 *   dataim_h  imaginary parts of data computed on the host;
 *   datare_d  real parts of data computed on the device;
 *   dataim_d  imaginary parts of data computed on the device;
 *   banner    string for printing if verbose level vrblvl > 1;
 *   vrblvl    is the verbose level. */

#endif
