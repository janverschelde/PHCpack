/* The file dbl2_convolutions_host.h specifies functions for the product of 
 * two series, truncated to the same degree, in double double precision. */

#ifndef __dbl2_convolutions_host_h__
#define __dbl2_convolutions_host_h__

void CPU_dbl2_product
 ( int deg, double *xhi, double *xlo, double *yhi, double *ylo,
            double *zhi, double *zlo );
/*
 * DESCRIPTION :
 *   Computes the product of two series truncated to the same degree,
 *   for real coefficients in double double precision.
 *
 * REQUIRED :
 *   The arrays xhi, xlo, yhi, ylo, zhi, and zlo have allocated space 
 *   for deg+1 doubles, for range 0 up to deg, deg included.
 *
 * ON ENTRY :
 *   deg      truncation degree of the series x, y, and z;
 *   xhi      high parts of the coefficients of the first series;
 *   xlo      low parts of the coefficients of the first series;
 *   yhi      high parts of the coefficients of the second series;
 *   ylo      low parts of the coefficients of the second series.
 *
 * ON RETURN :
 *   zhi      high parts of the coefficients of the product x*y;
 *   zlo      low parts of the coefficients of the product x*y. */

void CPU_cmplx2_product
 ( int deg, double *xrehi, double *xrelo, double *ximhi, double *ximlo,
            double *yrehi, double *yrelo, double *yimhi, double *yimlo,
            double *zrehi, double *zrelo, double *zimhi, double *zimlo );
/*
 * DESCRIPTION :
 *   Computes the product of two series truncated to the same degree,
 *   for complex coefficients in double double precision.
 *
 * REQUIRED :
 *   All arrays xrehi, xrelo, ximhi, ximlo, yrehi, yrelo, yimhi, yimlo,
 *   zrehi, zrelo, zimhi, and zimlo have allocated space 
 *   for deg+1 doubles, for range 0 up to deg, deg included.
 *
 * ON ENTRY :
 *   deg      truncation degree of the series x, y, and z;
 *   xrehi    high real parts of the coefficients of the first series;
 *   xrelo    low real parts of the coefficients of the first series;
 *   ximhi    high imaginary parts of the coefficients of the first series;
 *   ximlo    low imaginary parts of the coefficients of the first series;
 *   yrehi    high real parts of the coefficients of the second series.
 *   yrelo    low real parts of the coefficients of the second series.
 *   yimhi    high maginary parts of the coefficients of the second series.
 *   yimlo    low imaginary parts of the coefficients of the second series.
 *
 * ON RETURN :
 *   zrehi    high real parts of the coefficients of the product x*y;
 *   zrelo    low real parts of the coefficients of the product x*y;
 *   zimhi    high imaginary parts of the coefficients of the product x*y;
 *   zimlo    low imaginary parts of the coefficients of the product x*y. */

#endif
