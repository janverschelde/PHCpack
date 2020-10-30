/* The file dbl4_convolutions_host.h specifies functions for the product of 
 * two series, truncated to the same degree, in quad double precision. */

#ifndef __dbl4_convolutions_host_h__
#define __dbl4_convolutions_host_h__

void CPU_dbl4_product
 ( int deg, double *xhihi, double *xlohi, double *xhilo, double *xlolo,
            double *yhihi, double *ylohi, double *yhilo, double *ylolo,
            double *zhihi, double *zlohi, double *zhilo, double *zlolo );
/*
 * DESCRIPTION :
 *   Computes the product of two series truncated to the same degree,
 *   for real coefficients in quad double precision.
 *
 * REQUIRED :
 *   The arrays xhihi, xlohi, xhilo, xlolo, yhihi, ylohi, yhilo, ylolo,
 *   zhihi, zlohi, zhilo, and zlolo have allocated space for deg+1 doubles,
 *   for range 0 up to deg, deg included.
 *
 * ON ENTRY :
 *   deg      truncation degree of the series x, y, and z;
 *   xhihi    highest parts of the coefficients of the first series;
 *   xlohi    second highest parts of the coefficients of the first series;
 *   xhilo    second lowest parts of the coefficients of the first series;
 *   xlolo    lowest parts of the coefficients of the first series;
 *   yhihi    highest parts of the coefficients of the second series;
 *   ylohi    second highest parts of the coefficients of the second series;
 *   yhilo    second lowest parts of the coefficients of the second series;
 *   ylolo    lowest parts of the coefficients of the second series.
 *
 * ON RETURN :
 *   zhihi    highest parts of the coefficients of the product x*y;
 *   zlohi    second highest parts of the coefficients of the product x*y;
 *   zhilo    second lowest parts of the coefficients of the product x*y;
 *   zlolo    lowest parts of the coefficients of the product x*y. */

void CPU_cmplx4_product
 ( int deg,
   double *xrehihi, double *xrelohi, double *xrehilo, double *xrelolo,
   double *ximhihi, double *ximlohi, double *ximhilo, double *ximlolo,
   double *yrehihi, double *yrelohi, double *yrehilo, double *yrelolo,
   double *yimhihi, double *yimlohi, double *yimhilo, double *yimlolo,
   double *zrehihi, double *zrelohi, double *zrehilo, double *zrelolo,
   double *zimhihi, double *zimlohi, double *zimhilo, double *zimlolo );
/*
 * DESCRIPTION :
 *   Computes the product of two series truncated to the same degree,
 *   for complex coefficients in quad double precision.
 *
 * REQUIRED :
 *   All arrays xrehihi, xrelohi, xrehilo, xrelolo, ximhihi, ximlohi, ximhilo,
 *   ximlolo, yrehihi, yrelohi, yrehilo, yrelolo, yimhihi, yimlohi, yimhilo,
 *   yimlolo, zrehihi, zrelohi, zrehilo, zimhihi, zimlohi, zimhilo, and 
 *   zimlolo have allocated space  for deg+1 doubles, for range 0 up to deg,
 *   deg included.
 *
 * ON ENTRY :
 *   deg      truncation degree of the series x, y, and z;
 *   xrehihi  highest real parts of the coefficients of the first series;
 *   xrelohi  2nd highest real parts of the coefficients of the first series;
 *   xrehilo  2nd lowest real parts of the coefficients of the first series;
 *   xrelolo  lowest real parts of the coefficients of the first series;
 *   ximhihi  highest imaginary parts of the coefficients of the first series;
 *   ximlohi  2nd highest imag parts of the coefficients of the first series;
 *   ximhilo  2nd lowest imag parts of the coefficients of the first series;
 *   ximlolo  lowest imaginary parts of the coefficients of the first series;
 *   yrehihi  highest real parts of the coefficients of the second series;
 *   yrelohi  2nd highest real parts of the coefficients of the second series;
 *   yrehilo  2nd lowest real parts of the coefficients of the second series;
 *   yrelolo  lowest real parts of the coefficients of the second series;
 *   yimhihi  highest imag parts of the coefficients of the second series;
 *   yimlohi  2nd highest imag parts of the coefficients of the second series;
 *   yimhilo  2nd lowest imag parts of the coefficients of the second series;
 *   yimlolo  lowest imaginary parts of the coefficients of the second series.
 *
 * ON RETURN :
 *   zrehihi  highest real parts of the coefficients of the product x*y;
 *   zrelohi  2nd highest real parts of the coefficients of the product x*y;
 *   zrehilo  2nd lowest real parts of the coefficients of the product x*y;
 *   zrelolo  lowest real parts of the coefficients of the product x*y;
 *   zimhihi  highest imaginary parts of the coefficients of the product x*y;
 *   zimlohi  2nd highest imag parts of the coefficients of the product x*y;
 *   zimhilo  2nd lowest imag parts of the coefficients of the product x*y;
 *   zimlolo  lowest imaginary parts of the coefficients of the product x*y. */

#endif
