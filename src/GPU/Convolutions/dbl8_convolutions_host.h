/* The file dbl8_convolutions_host.h specifies functions for the product of 
 * two series, truncated to the same degree, in octo double precision. */

#ifndef __dbl8_convolutions_host_h__
#define __dbl8_convolutions_host_h__

void CPU_dbl8_product
 ( int deg,
   double *xhihihi, double *xlohihi, double *xhilohi, double *xlolohi,
   double *xhihilo, double *xlohilo, double *xhilolo, double *xlololo,
   double *yhihihi, double *ylohihi, double *yhilohi, double *ylolohi,
   double *yhihilo, double *ylohilo, double *yhilolo, double *ylololo,
   double *zhihihi, double *zlohihi, double *zhilohi, double *zlolohi,
   double *zhihilo, double *zlohilo, double *zhilolo, double *zlololo );
/*
 * DESCRIPTION :
 *   Computes the product of two series truncated to the same degree,
 *   for real coefficients in octo double precision.
 *
 * REQUIRED :
 *   All arrays have allocated space for deg+1 doubles,
 *   for range 0 up to deg, deg included.
 *
 * ON ENTRY :
 *   deg       truncation degree of the series x, y, and z;
 *   xhihihi   highest parts of the coefficients of x;
 *   xlohihi   second highest parts of the coefficients of x;
 *   xhilohi   third highest parts of the coefficients of x;
 *   xlolohi   fourth highest parts of the coefficients of x;
 *   xhihilo   fourth lowest parts of the coefficients of x;
 *   xlohilo   third lowest parts of the coefficients of x;
 *   xhilolo   second lowest parts of the coefficients of x;
 *   xlololo   lowest parts of the coefficients of x;
 *   yhihihi   highest parts of the coefficients of y;
 *   ylohihi   second highest parts of the coefficients of y;
 *   yhilohi   third highest parts of the coefficients of y;
 *   ylolohi   fourth highest parts of the coefficients of y;
 *   yhihilo   fourth lowest parts of the coefficients of y;
 *   ylohilo   third lowest parts of the coefficients of y;
 *   yhilolo   second lowest parts of the coefficients of y;
 *   ylololo   lowest parts of the coefficients of y.
 *
 * ON RETURN :
 *   zhihihi   highest parts of the coefficients of the product x*y;
 *   zlohihi   second highest parts of the coefficients of the product x*y;
 *   zhilohi   third highest parts of the coefficients of the product x*y;
 *   zlolohi   fourth highest parts of the coefficients of the product x*y;
 *   zhihilo   fourth lowest parts of the coefficients of the product x*y;
 *   zlohilo   third lowest parts of the coefficients of the product x*y;
 *   zhilolo   second lowest parts of the coefficients of the product x*y;
 *   zlololo   lowest parts of the coefficients of the product x*y. */

void CPU_cmplx8_product
 ( int deg,
   double *xrehihihi, double *xrelohihi, double *xrehilohi, double *xrelolohi,
   double *xrehihilo, double *xrelohilo, double *xrehilolo, double *xrelololo,
   double *ximhihihi, double *ximlohihi, double *ximhilohi, double *ximlolohi,
   double *ximhihilo, double *ximlohilo, double *ximhilolo, double *ximlololo,
   double *yrehihihi, double *yrelohihi, double *yrehilohi, double *yrelolohi,
   double *yrehihilo, double *yrelohilo, double *yrehilolo, double *yrelololo,
   double *yimhihihi, double *yimlohihi, double *yimhilohi, double *yimlolohi,
   double *yimhihilo, double *yimlohilo, double *yimhilolo, double *yimlololo,
   double *zrehihihi, double *zrelohihi, double *zrehilohi, double *zrelolohi,
   double *zrehihilo, double *zrelohilo, double *zrehilolo, double *zrelololo,
   double *zimhihihi, double *zimlohihi, double *zimhilohi, double *zimlolohi,
   double *zimhihilo, double *zimlohilo, double *zimhilolo, double *zimlololo );
/*
 * DESCRIPTION :
 *   Computes the product of two series truncated to the same degree,
 *   for complex coefficients in octo double precision.
 *
 * REQUIRED :
 *   All arrays have space allocated for deg+1 doubles,
 *   for range 0 up to deg, deg included.
 *
 * ON ENTRY :
 *   deg        truncation degree of the series x, y, and z;
 *   xrehihihi  highest real parts of the coefficients of x;
 *   xrelohihi  second highest real parts of the coefficients of x;
 *   xrehihihi  third highest real parts of the coefficients of x;
 *   xrelohihi  fourth third highest real parts of the coefficients of x;
 *   xrehihilo  fourth lowest real parts of the coefficients of x;
 *   xrelohilo  third lowest real parts of the coefficients of x;
 *   xrehilolo  second lowest real parts of the coefficients of x;
 *   xrelololo  lowest real parts of the coefficients of x;
 *   ximhihihi  highest imaginary parts of the coefficients of x;
 *   ximlohihi  second highest imaginary parts of the coefficients of x;
 *   ximhilohi  third highest imaginary parts of the coefficients of x;
 *   ximlolohi  fourth highest imaginary parts of the coefficients of x;
 *   ximhihilo  fourth lowest imaginary parts of the coefficients of x;
 *   ximlohilo  third lowest imaginary parts of the coefficients of x;
 *   ximhilolo  second lowest imaginary parts of the coefficients of x;
 *   ximlololo  lowest imaginary parts of the coefficients of x;
 *   yrehihihi  highest real parts of the coefficients of y;
 *   yrelohihi  second highest real parts of the coefficients of y;
 *   yrehilohi  third highest real parts of the coefficients of y;
 *   yrelolohi  fourth highest real parts of the coefficients of y;
 *   yrehihilo  fourth lowest real parts of the coefficients of y;
 *   yrelohilo  third lowest real parts of the coefficients of y;
 *   yrehilolo  second lowest real parts of the coefficients of y;
 *   yrelololo  lowest real parts of the coefficients of y;
 *   yimhihihi  highest imaginary parts of the coefficients of y;
 *   yimlohihi  second highest imaginary parts of the coefficients of y;
 *   yimhilohi  third highest imaginary parts of the coefficients of y;
 *   yimlolohi  fourth highest imaginary parts of the coefficients of y;
 *   yimhihilo  fourth lowest imaginary parts of the coefficients of y;
 *   yimlohilo  third lowest imaginary parts of the coefficients of y;
 *   yimhilolo  second lowest imaginary parts of the coefficients of y;
 *   yimlololo  lowest imaginary parts of the coefficients of y.
 *
 * ON RETURN :
 *   zrehihihi  highest real parts of the coefficients of x*y;
 *   zrelohihi  second highest real parts of the coefficients of x*y;
 *   zrehilohi  third highest real parts of the coefficients of x*y;
 *   zrelolohi  fourth highest real parts of the coefficients of x*y;
 *   zrehihilo  fourth lowest real parts of the coefficients of x*y;
 *   zrelohilo  third lowest real parts of the coefficients of x*y;
 *   zrehilolo  second lowest real parts of the coefficients of x*y;
 *   zrehilolo  lowest real parts of the coefficients of x*y;
 *   zimhihihi  highest imaginary parts of the coefficients of x*y;
 *   zimlohihi  second highest imaginary parts of the coefficients of x*y;
 *   zimhilohi  third highest imaginary parts of the coefficients of x*y;
 *   zimlolohi  fourth highest imaginary parts of the coefficients of x*y;
 *   zimhihilo  fourth lowest imaginary parts of the coefficients of x*y;
 *   zimlohilo  third imaginary parts of the coefficients of x*y;
 *   zimhilolo  second lowest imaginary parts of the coefficients of x*y;
 *   zimlololo  lowest imaginary parts of the coefficients of x*y. */

#endif
