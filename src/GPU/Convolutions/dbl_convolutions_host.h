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

void CPU_dbl_Laurent_product
 ( int deg, int xe, int ye, int *ze, double *x, double *y, double *z );
/*
 * DESCRIPTION :
 *   Computes the product of two Laurent series with the same precision,
 *   for real double coefficients.
 *
 * REQUIRED :
 *   The arrays x, y, and z have allocated space for deg+1 doubles,
 *   for range 0 up to deg, deg included.
 *
 * ON ENTRY :
 *   deg      truncation degree of the series x, y, and z;
 *   xe       the leading exponent of the first series x;
 *   ye       the leading exponent of the second series y;
 *   x        coefficients of the first series;
 *   y        coefficients of the second series.
 *
 * ON RETURN :
 *   ze       the leading exponent of the product of x with y;
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

void CPU_cmplx_Laurent_product
 ( int deg, int xe, int ye, int *ze, double *xre, double *xim,
   double *yre, double *yim, double *zre, double *zim );
/*
 * DESCRIPTION :
 *   Computes the product of two Laurent series with the same precision,
 *   for complex double coefficients. 
 *
 * REQUIRED :
 *   All arrays xre, xim, yre, yim, zre, and zim have allocated space 
 *   for deg+1 doubles, for range 0 up to deg, deg included.
 *
 * ON ENTRY :
 *   deg      truncation degree of the series x, y, and z;
 *   xe       the leading exponent of the first series x;
 *   ye       the leading exponent of the second series y;
 *   xre      real parts of the coefficients of the first series;
 *   xim      imaginary parts of the coefficients of the first series;
 *   yre      real parts of the coefficients of the second series.
 *   yim      imaginary parts of the coefficients of the second series.
 *
 * ON RETURN :
 *   ze       the leading exponent of the product of x with y;
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

void CPU_dbl_Laurent_inverse
 ( int deg, int xe, int *ye, double *x, double *y );
/*
 * DESCRIPTION :
 *   Computes the inverse of the Laurent series for real double coefficients.
 *
 * REQUIRED : The lead coefficient x(0) of x is nonzero.
 *   The arrays x and y have allocated space for deg+1 doubles,
 *   for range 0 up to deg, deg included.
 *
 * ON ENTRY :
 *   deg      truncation degree of the series x and y;
 *   xe       the leading exponent of the series x;
 *   x        coefficients of a power series;
 *   y        space for deg+1 coefficients.
 *
 * ON RETURN :
 *   ye       the leading exponent of the inverse of the series x;
 *   y        coefficients of the inverse of the series x. */

void CPU_cmplx_inverse
 ( int deg, double *xre, double *xim, double *yre, double *yim );
/*
 * DESCRIPTION :
 *   Computes the inverse of the series for complex coefficients
 *   in double precision.
 *
 * REQUIRED : The lead coefficient x(0) of x is nonzero.
 *   The arrays x and y have allocated space for deg+1 doubles,
 *   for range 0 up to deg, deg included.
 *
 * ON ENTRY :
 *   deg      truncation degree of the series x and y;
 *   xre      real parts of the coefficients of a power series;
 *   xim      imaginary parts of the coefficients of a power series;
 *   yre      space for deg+1 coefficients;
 *   yim      space for deg+1 coefficients.
 *
 * ON RETURN :
 *   yre      real parts of the coefficients of the inverse of x;
 *   yim      imaginary parts of the coefficients of the inverse of x. */

void CPU_cmplx_Laurent_inverse
 ( int deg, int xe, int *ye, double *xre, double *xim,
   double *yre, double *yim );
/*
 * DESCRIPTION :
 *   Computes the inverse of the Laurent series for complex coefficients
 *   in double precision.
 *
 * REQUIRED : The lead coefficient x(0) of x is nonzero.
 *   The arrays x and y have allocated space for deg+1 doubles,
 *   for range 0 up to deg, deg included.
 *
 * ON ENTRY :
 *   deg      truncation degree of the series x and y;
 *   xe       the leading coefficient of the series x;
 *   xre      real parts of the coefficients of a power series;
 *   xim      imaginary parts of the coefficients of a power series;
 *   yre      space for deg+1 coefficients;
 *   yim      space for deg+1 coefficients.
 *
 * ON RETURN :
 *   ye       the leading exponent of the inverse;
 *   yre      real parts of the coefficients of the inverse of x;
 *   yim      imaginary parts of the coefficients of the inverse of x. */

#endif
