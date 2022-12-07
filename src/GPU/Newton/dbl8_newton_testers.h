// The file dbl8_newton_testers.h specifies test function for Newton's method
// on series in octo double precision.

#ifndef __dbl8_newton_testers_h__
#define __dbl8_newton_testers_h__

void real8_start_series_vector
 ( int dim, int deg,
   double *r0hihihi, double *r0lohihi, double *r0hilohi, double *r0lolohi,
   double *r0hihilo, double *r0lohilo, double *r0hilolo, double *r0lololo,
   double **cffhihihi, double **cfflohihi,
   double **cffhilohi, double **cfflolohi,
   double **cffhihilo, double **cfflohilo,
   double **cffhilolo, double **cfflololo );
/*
 * DESCRIPTION :
 *   Given space in cff* for the doubles
 *   for the vector of dim power series in octo double precision,
 *   of series truncated after degree deg,
 *   sets the coefficients of the start series.
 *   The constant term equals r0[j] with a small error. */

void cmplx8_start_series_vector
 ( int dim, int deg,
   double *r0rehihihi, double *r0relohihi,
   double *r0rehilohi, double *r0relolohi,
   double *r0rehihilo, double *r0relohilo,
   double *r0rehilolo, double *r0relololo,
   double *r0imhihihi, double *r0imlohihi,
   double *r0imhilohi, double *r0imlolohi,
   double *r0imhihilo, double *r0imlohilo,
   double *r0imhilolo, double *r0imlololo,
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
 *   sets the coefficients of the start series.
 *   The constant term of each series is set the (r0re, r0im),
 *   with some small error. */

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
 *   sets the coefficients of the star series. */

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

void dbl8_unit_series_vectors
 ( int nbr, int dim, int deg,
   double ***cffhihihi, double ***cfflohihi,
   double ***cffhilohi, double ***cfflolohi,
   double ***cffhihilo, double ***cfflohilo,
   double ***cffhilolo, double ***cfflololo );
/*
 * DESCRIPTION :
 *   Given space for nbr vectors of dim power series,
 *   of series truncated after degree deg,
 *   initializes the coefficients in cff to one as leading coefficients,
 *   and zero for all other coefficients. */

void cmplx8_unit_series_vectors
 ( int nbr, int dim, int deg,
   double ***cffrehihihi, double ***cffrelohihi,
   double ***cffrehilohi, double ***cffrelolohi,
   double ***cffrehihilo, double ***cffrelohilo,
   double ***cffrehilolo, double ***cffrelololo,
   double ***cffimhihihi, double ***cffimlohihi,
   double ***cffimhilohi, double ***cffimlolohi,
   double ***cffimhihilo, double ***cffimlohilo,
   double ***cffimhilolo, double ***cffimlololo );
/*
 * DESCRIPTION :
 *   Given space allocated in cffre and cffim for nbr series,
 *   returns the values for a complex unit series. */

void dbl8_update_series
 ( int dim, int degp1, int startidx,
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
 *   startidx  index of the start of the update;
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
 ( int dim, int degp1, int startidx,
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
 *   startidx  index of the start of the update;
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

double dbl8_error3sum
 ( int dim1, int dim2, int dim3,
   double ***datahihihi_h, double ***datalohihi_h,
   double ***datahilohi_h, double ***datalolohi_h,
   double ***datahihilo_h, double ***datalohilo_h,
   double ***datahilolo_h, double ***datalololo_h,
   double ***datahihihi_d, double ***datalohihi_d,
   double ***datahilohi_d, double ***datalolohi_d,
   double ***datahihilo_d, double ***datalohilo_d,
   double ***datahilolo_d, double ***datalololo_d,
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
 *   datahihihi_h are highest doubles of data computed on the host;
 *   datalohihi_h are 2nd highest doubles of data computed on the host;
 *   datahilohi_h are 3rd highest doubles of data computed on the host;
 *   datalolohi_h are 4th highest doubles of data computed on the host;
 *   datahihilo_h are 4th lowest doubles of data computed on the host;
 *   datalohilo_h are 3rd lowest doubles of data computed on the host;
 *   datahilolo_h are 2nd lowest doubles of data computed on the host;
 *   datalololo_h are lowest doubles of data computed on the host;
 *   datahihihi_d are highest doubles of data computed on the device;
 *   datalohihi_d are 2nd highest doubles of data computed on the device;
 *   datahilohi_d are 3rd highest doubles of data computed on the device;
 *   datalolohi_d are 4th highest doubles of data computed on the device;
 *   datahihilo_d are 4th lowest doubles of data computed on the device;
 *   datalohilo_d are 3rd lowest doubles of data computed on the device;
 *   datahilolo_d are 2nd lowest doubles of data computed on the device;
 *   datalololo_d are lowest doubles of data computed on the device;
 *   banner    string for printing if verbose level vrblvl > 1;
 *   vrblvl    is the verbose level. */

double cmplx8_error3sum
 ( int dim1, int dim2, int dim3,
   double ***datarehihihi_h, double ***datarelohihi_h,
   double ***datarehilohi_h, double ***datarelolohi_h,
   double ***datarehihilo_h, double ***datarelohilo_h,
   double ***datarehilolo_h, double ***datarelololo_h,
   double ***dataimhihihi_h, double ***dataimlohihi_h,
   double ***dataimhilohi_h, double ***dataimlolohi_h,
   double ***dataimhihilo_h, double ***dataimlohilo_h,
   double ***dataimhilolo_h, double ***dataimlololo_h,
   double ***datarehihihi_d, double ***datarelohihi_d,
   double ***datarehilohi_d, double ***datarelolohi_d,
   double ***datarehihilo_d, double ***datarelohilo_d,
   double ***datarehilolo_d, double ***datarelololo_d,
   double ***dataimhihihi_d, double ***dataimlohihi_d,
   double ***dataimhilohi_d, double ***dataimlolohi_d,
   double ***dataimhihilo_d, double ***dataimlohilo_d,
   double ***dataimhilolo_d, double ***dataimlololo_d,
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
 *   datarehihihi_h are highest doubles of real parts of data, on host;
 *   datarelohihi_h are 2nd highest doubles of real parts of data, on host;
 *   datarehilohi_h are 3rd highest doubles of real parts of data, on host;
 *   datarelolohi_h are 4th highest doubles of real parts of data, on host;
 *   datarehihilo_h are 4th lowest doubles of real parts of data, on host;
 *   datarelohilo_h are 3rd lowest doubles of real parts of data, on host;
 *   datarehilolo_h are 2nd lowest doubles of real parts of data, on host;
 *   datarelololo_h are lowest doubles of real parts of data, on host;
 *   dataimhihihi_h are highest doubles of imaginary parts of data, on host;
 *   dataimlohihi_h are 2nd highest doubles of imag parts of data, on host;
 *   dataimhilohi_h are 3rd highest doubles of imag parts of data, on host;
 *   dataimlolohi_h are 4th highest doubles of imag parts of data, on host;
 *   dataimhihilo_h are 4th lowest doubles of imag parts of data, on host;
 *   dataimlohilo_h are 3rd lowest doubles of imag parts of data, on host;
 *   dataimhilolo_h are 2nd lowest doubles of imag parts of data, on host;
 *   dataimlololo_h are lowest doubles of imaginary parts of data, on host;
 *   datarehihihi_d are highest doubles of real parts of data, on device;
 *   datarelohihi_d are 2nd highest doubles of real parts of data, on device;
 *   datarehilohi_d are 3rd highest doubles of real parts of data, on device;
 *   datarelolohi_d are 4th highest doubles of real parts of data, on device;
 *   datarehihilo_d are 4th lowest doubles of real parts of data, on device;
 *   datarelohilo_d are 3rd lowest doubles of real parts of data, on device;
 *   datarehilolo_d are 2nd lowest doubles of real parts of data, on device;
 *   datarelololo_d are lowest doubles of real parts of data, on device;
 *   dataimhihihi_d are highest doubles of imaginary parts of data, on device;
 *   dataimlohihi_d are 2nd highest doubles of imag parts of data, on device;
 *   dataimhilohi_d are 3rd highest doubles of imag parts of data, on device;
 *   dataimlolohi_d are 4th highest doubles of imag parts of data, on device;
 *   dataimhihilo_d are 4th lowest doubles of imag parts of data, on device;
 *   dataimlohilo_d are 3rd lowest doubles of imag parts of data, on device;
 *   dataimhilolo_d are 2nd lowest doubles of imag parts of data, on device;
 *   dataimlololo_d are lowest doubles of imaginary parts of data, on device;
 *   banner    string for printing if verbose level vrblvl > 1;
 *   vrblvl    is the verbose level. */

double dbl8_error2sum
 ( int nrows, int ncols,
   double **datahihihi_h, double **datalohihi_h,
   double **datahilohi_h, double **datalolohi_h,
   double **datahihilo_h, double **datalohilo_h,
   double **datahilolo_h, double **datalololo_h,
   double **datahihihi_d, double **datalohihi_d,
   double **datahilohi_d, double **datalolohi_d,
   double **datahihilo_d, double **datalohilo_d,
   double **datahilolo_d, double **datalololo_d,
   std::string banner, int vrblvl );
/*
 * DESCRIPTION :
 *   Returns the sum of the absolute values of the differences
 *   between host and device real data, in two dimensions.
 *
 * ON ENTRY :
 *   nrows     number of rows of the matrices;
 *   ncols     number of columns of the matrices;
 *   datahihihi_h are highest doubles of data computed on the host;
 *   datalohihi_h are 2nd highest doubles of data computed on the host;
 *   datahilohi_h are 3rd highest doubles of data computed on the host;
 *   datalolohi_h are 4th highest doubles of data computed on the host;
 *   datahihilo_h are 4th lowest doubles of data computed on the host;
 *   datalohilo_h are 3rd lowest doubles of data computed on the host;
 *   datahilolo_h are 2nd lowest doubles of data computed on the host;
 *   datalololo_h are lowest doubles of data computed on the host;
 *   datahihihi_d are highest doubles of data computed on the device;
 *   datalohihi_d are 2nd highest doubles of data computed on the device;
 *   datahilohi_d are 3rd highest doubles of data computed on the device;
 *   datalolohi_d are 4th highest doubles of data computed on the device;
 *   datahihilo_d are 4th lowest doubles of data computed on the device;
 *   datalohilo_d are 3rd lowest doubles of data computed on the device;
 *   datahilolo_d are 2nd lowest doubles of data computed on the device;
 *   datalololo_d are lowest doubles of data computed on the device;
 *   banner    string for printing if verbose level vrblvl > 1;
 *   vrblvl    is the verbose level. */

double cmplx8_error2sum
 ( int nrows, int ncols,
   double **datarehihihi_h, double **datarelohihi_h,
   double **datarehilohi_h, double **datarelolohi_h,
   double **datarehihilo_h, double **datarelohilo_h,
   double **datarehilolo_h, double **datarelololo_h,
   double **dataimhihihi_h, double **dataimlohihi_h,
   double **dataimhilohi_h, double **dataimlolohi_h,
   double **dataimhihilo_h, double **dataimlohilo_h,
   double **dataimhilolo_h, double **dataimlololo_h,
   double **datarehihihi_d, double **datarelohihi_d,
   double **datarehilohi_d, double **datarelolohi_d,
   double **datarehihilo_d, double **datarelohilo_d,
   double **datarehilolo_d, double **datarelololo_d,
   double **dataimhihihi_d, double **dataimlohihi_d,
   double **dataimhilohi_d, double **dataimlolohi_d,
   double **dataimhihilo_d, double **dataimlohilo_d,
   double **dataimhilolo_d, double **dataimlololo_d,
   std::string banner, int vrblvl );
/*
 * DESCRIPTION :
 *   Returns the sum of the absolute values of the differences
 *   between host and device complex data, in two dimensions.
 *
 * ON ENTRY :
 *   nrows     number of rows of the matrices;
 *   ncols     number of columns of the matrices;
 *   datarehihihi_h are highest doubles of real parts of data, on host;
 *   datarelohihi_h are 2nd highest doubles of real parts of data, on host;
 *   datarehilohi_h are 3rd highest doubles of real parts of data, on host;
 *   datarelolohi_h are 4th highest doubles of real parts of data, on host;
 *   datarehihilo_h are 4th lowest doubles of real parts of data, on host;
 *   datarelohilo_h are 3rd lowest doubles of real parts of data, on host;
 *   datarehilolo_h are 2nd lowest doubles of real parts of data, on host;
 *   datarelololo_h are lowest doubles of real parts of data, on host;
 *   dataimhihihi_h are highest doubles of imaginary parts of data, on host;
 *   dataimlohihi_h are 2nd highest doubles of imag parts of data, on host;
 *   dataimhilohi_h are 3rd highest doubles of imag parts of data, on host;
 *   dataimlolohi_h are 4th highest doubles of imag parts of data, on host;
 *   dataimhihilo_h are 4th lowest doubles of imag parts of data, on host;
 *   dataimlohilo_h are 3rd lowest doubles of imag parts of data, on host;
 *   dataimhilolo_h are 2nd lowest doubles of imag parts of data, on host;
 *   dataimlololo_h are lowest doubles of imaginary parts of data, on host;
 *   datarehihihi_d are highest doubles of real parts of data, on device;
 *   datarelohihi_d are 2nd highest doubles of real parts of data, on device;
 *   datarehilohi_d are 3rd highest doubles of real parts of data, on device;
 *   datarelolohi_d are 4th highest doubles of real parts of data, on device;
 *   datarehihilo_d are 4th lowest doubles of real parts of data, on device;
 *   datarelohilo_d are 3rd lowest doubles of real parts of data, on device;
 *   datarehilolo_d are 2nd lowest doubles of real parts of data, on device;
 *   datarelololo_d are lowest doubles of real parts of data, on device;
 *   dataimhihihi_d are highest doubles of imaginary parts of data, on device;
 *   dataimlohihi_d are 2nd highest doubles of imag parts of data, on device;
 *   dataimhilohi_d are 3rd highest doubles of imag parts of data, on device;
 *   dataimlolohi_d are 4th highest doubles of imag parts of data, on device;
 *   dataimhihilo_d are 4th lowest doubles of imag parts of data, on device;
 *   dataimlohilo_d are 3rd lowest doubles of imag parts of data, on device;
 *   dataimhilolo_d are 2nd lowest doubles of imag parts of data, on device;
 *   dataimlololo_d are lowest doubles of imaginary parts of data, on device;
 *   banner    string for printing if verbose level vrblvl > 1;
 *   vrblvl    is the verbose level. */

#endif
