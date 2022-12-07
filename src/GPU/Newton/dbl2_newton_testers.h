// The file dbl2_newton_testers.h specifies test function for Newton's method
// on series in double double precision.

#ifndef __dbl2_newton_testers_h__
#define __dbl2_newton_testers_h__

void real2_start_series_vector
 ( int dim, int deg, double *r0hi, double *r0lo,
   double **cffhi, double **cfflo );
/*
 * DESCRIPTION :
 *   Given space in cffhi, cfflo for the high and low doubles
 *   for the vector of dim power series,
 *   of series truncated after degree deg,
 *   sets the coefficients of the start series.
 *   The constant term equals r0[j] with a small error. */

void cmplx2_start_series_vector
 ( int dim, int deg,
   double *r0rehi, double *r0relo, double *r0imhi, double *r0imlo,
   double **cffrehi, double **cffrelo, double **cffimhi, double **cffimlo );
/*
 * DESCRIPTION :
 *   Given space allocated in the coefficient arrays,
 *   sets the coefficients of the start series.
 *   The constant terms of each series is set to (r0re, r0im),
 *   with some small error. */

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

void dbl2_unit_series_vectors
 ( int nbr, int dim, int deg, double ***cffhi, double ***cfflo );
/*
 * DESCRIPTION :
 *   Given space in cffhi, cfflo for nbr vectors of dim power series,
 *   of series truncated after degree deg,
 *   initializes the coefficients in cff to one as leading coefficients,
 *   and zero for all other coefficients. */

void cmplx2_unit_series_vectors
 ( int nbr, int dim, int deg,
   double ***cffrehi, double ***cffrelo,
   double ***cffimhi, double ***cffimlo );
/*
 * DESCRIPTION :
 *   Given space allocated in cffre and cffim for nbr series,
 *   returns the values for a complex unit series. */

void dbl2_update_series
 ( int dim, int degp1, int startidx, double **xhi, double **xlo,
   double **dxhi, double **dxlo, int vrblvl );
/*
 * DESCRIPTION :
 *   Adds the series in dx to x.
 *
 * ON ENTRY :
 *   dim       number of series in x and dx;
 *   degp1     degree plus one of the series;
 *   startidx  index of the start of the update;
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
 ( int dim, int degp1, int startidx,
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
 *   startidx  index of the start of the update;
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

#endif
