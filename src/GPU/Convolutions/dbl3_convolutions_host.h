/* The file dbl3_convolutions_host.h specifies functions for the product of 
 * two series, truncated to the same degree, in triple double precision. */

#ifndef __dbl3_convolutions_host_h__
#define __dbl3_convolutions_host_h__

void CPU_dbl3_product
 ( int deg, double *xhi, double *xmi, double *xlo,
            double *yhi, double *ymi, double *ylo,
            double *zhi, double *zmi, double *zlo );
/*
 * DESCRIPTION :
 *   Computes the product of two series truncated to the same degree,
 *   for real coefficients in triple double precision.
 *
 * REQUIRED :
 *   The arrays xhi, xmi, xlo, yhi, ymi, ylo, zhi, zmi, and zlo have 
 *   allocated space for deg+1 doubles, for range 0 up to deg, deg included.
 *
 * ON ENTRY :
 *   deg      truncation degree of the series x, y, and z;
 *   xhi      high parts of the coefficients of the first series;
 *   xmi      middle parts of the coefficients of the first series;
 *   xlo      low parts of the coefficients of the first series;
 *   yhi      high parts of the coefficients of the second series;
 *   ymi      middle parts of the coefficients of the second series;
 *   ylo      low parts of the coefficients of the second series.
 *
 * ON RETURN :
 *   zhi      high parts of the coefficients of the product x*y;
 *   zmi      middle parts of the coefficients of the product x*y;
 *   zlo      low parts of the coefficients of the product x*y. */

void CPU_cmplx3_product
 ( int deg, double *xrehi, double *xremi, double *xrelo,
            double *ximhi, double *ximmi, double *ximlo,
            double *yrehi, double *yremi, double *yrelo,
            double *yimhi, double *yimmi, double *yimlo,
            double *zrehi, double *zremi, double *zrelo,
            double *zimhi, double *zimmi, double *zimlo );
/*
 * DESCRIPTION :
 *   Computes the product of two series truncated to the same degree,
 *   for complex coefficients in triple double precision.
 *
 * REQUIRED :
 *   All arrays xrehi, xremi, xrelo, ximhi, ximmi, ximlo, yrehi, yremi,
 *   yrelo, yimhi, yimmi, yimlo, zrehi, zremi, zrelo, zimhi, zimmi, and zimlo
 *   have allocated space  for deg+1 doubles, for range 0 up to deg,
 *   deg included.
 *
 * ON ENTRY :
 *   deg      truncation degree of the series x, y, and z;
 *   xrehi    high real parts of the coefficients of the first series;
 *   xremi    middle real parts of the coefficients of the first series;
 *   xrelo    low real parts of the coefficients of the first series;
 *   ximhi    high imaginary parts of the coefficients of the first series;
 *   ximmi    middle imaginary parts of the coefficients of the first series;
 *   ximlo    low imaginary parts of the coefficients of the first series;
 *   yrehi    high real parts of the coefficients of the second series.
 *   yremi    middle real parts of the coefficients of the second series.
 *   yrelo    low real parts of the coefficients of the second series.
 *   yimhi    high imaginary parts of the coefficients of the second series.
 *   yimmi    middle imaginary parts of the coefficients of the second series.
 *   yimlo    low imaginary parts of the coefficients of the second series.
 *
 * ON RETURN :
 *   zrehi    high real parts of the coefficients of the product x*y;
 *   zremi    middle real parts of the coefficients of the product x*y;
 *   zrelo    low real parts of the coefficients of the product x*y;
 *   zimhi    high imaginary parts of the coefficients of the product x*y;
 *   zimmi    middle imaginary parts of the coefficients of the product x*y;
 *   zimlo    low imaginary parts of the coefficients of the product x*y. */

#endif
