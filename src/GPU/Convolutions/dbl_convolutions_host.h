/* The file dbl_convolutions_host.h specifies functions for the product of 
 * two series, truncated to the same degree, in double precision. */

#ifndef __dbl_convolutions_host_h__
#define __dbl_convolutions_host_h__

void CPU_dbl_product ( int deg, double *x, double *y, double *z );
/*
 * DESCRIPTION :
 *   Computes the product of two series truncated to the same degree,
 *   for real coefficients in double precision.
 *
 * REQUIRED :
 *   The arrays x, y, and z have allocated space for deg+1 doubles,
 *   for range 0 up to deg, deg included.
 *
 * ON ENTRY :
 *   deg      truncation degree of the series x, y, and z;
 *   x        coefficients of the first series;
 *   y        coefficients of the second series.
 *
 * ON RETURN :
 *   z        coefficients of the product x*y. */

void CPU_cmplx_product
 ( int deg, double *xre, double *xim, double *yre, double *yim,
            double *zre, double *zim );
/*
 * DESCRIPTION :
 *   Computes the product of two series truncated to the same degree,
 *   for complex coefficients in double precision.
 *
 * REQUIRED :
 *   All arrays xre, xim, yre, yim, zre, and zim have allocated space 
 *   for deg+1 doubles, for range 0 up to deg, deg included.
 *
 * ON ENTRY :
 *   deg      truncation degree of the series x, y, and z;
 *   xre      real parts of the coefficients of the first series;
 *   xim      imaginary parts of the coefficients of the first series;
 *   yre      real parts of the coefficients of the second series.
 *   yim      imaginary parts of the coefficients of the second series.
 *
 * ON RETURN :
 *   zre      real parts of the coefficients of the product x*y;
 *   zim      imaginary parts of the coefficients of the product x*y. */

void CPU_dbl_inverse ( int deg, double *x, double *y );
/*
 * DESCRIPTION :
 *   Computes the inverse of the series for real coefficients
 *   in double precision.
 *
 * REQUIRED : The lead coefficient x(0) of x is nonzero.
 *   The arrays x and y have allocated space for deg+1 doubles,
 *   for range 0 up to deg, deg included.
 *
 * ON ENTRY :
 *   deg      truncation degree of the series x and y;
 *   x        coefficients of a power series;
 *   y        space for deg+1 coefficients.
 *
 * ON RETURN :
 *   y        coefficients of the inverse of the series x. */

#endif
