// The file dbl4_newton_testers.h specifies test function for Newton's method
// on series in quad double precision.

#ifndef __dbl4_newton_testers_h__
#define __dbl4_newton_testers_h__

void real4_start_series_vector
 ( int dim, int deg,
   double *r0hihi, double *r0lohi, double *r0hilo, double *r0lolo,
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo );
/*
 * DESCRIPTION :
 *   Given space in cffhihi, cfflohi, cffhilo, cfflolo for the doubles
 *   for the vector of dim power series in quad double precision,
 *   of series truncated after degree deg,
 *   sets the coefficients of the start series.
 *   The constant terms equals r0[j] with a small error. */

void cmplx4_start_series_vector
 ( int dim, int deg,
   double *r0rehihi, double *r0relohi, double *r0rehilo, double *r0relolo,
   double *r0imhihi, double *r0imlohi, double *r0imhilo, double *r0imlolo,
   double **cffrehihi, double **cffrelohi,
   double **cffrehilo, double **cffrelolo,
   double **cffimhihi, double **cffimlohi,
   double **cffimhilo, double **cffimlolo );
/*
 * DESCRIPTION :
 *   Given space allocated in the coefficient arrays,
 *   sets the coefficients of the start series.
 *   The constant term of each series is set the (r0re, r0im),
 *   with some small error. */

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

void dbl4_unit_series_vectors
 ( int nbr, int dim, int deg,
   double ***cffhihi, double ***cfflohi,
   double ***cffhilo, double ***cfflolo );
/*
 * DESCRIPTION :
 *   Given space for nbr vectors of dim power series,
 *   of series truncated after degree deg,
 *   initializes the coefficients in cff to one as leading coefficients,
 *   and zero for all other coefficients. */

void cmplx4_unit_series_vectors
 ( int nbr, int dim, int deg,
   double ***cffrehihi, double ***cffrelohi,
   double ***cffrehilo, double ***cffrelolo,
   double ***cffimhihi, double ***cffimlohi,
   double ***cffimhilo, double ***cffimlolo );
/*
 * DESCRIPTION :
 *   Given space allocated in cffre and cffim for nbr series,
 *   returns the values for a complex unit series. */

void dbl4_update_series
 ( int dim, int degp1, int startidx,
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
 *   startidx  index of the start of the update;
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
 ( int dim, int degp1, int startidx,
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
 *   startidx  index of the start of the update;
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

#endif
