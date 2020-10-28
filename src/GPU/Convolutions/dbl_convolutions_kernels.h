// The file dbl_convolutions_kernels.h defines prototypes for the kernels 
// to multiply two truncated power series in double precision.
// Defines the constant bound on the shared memory size of one block.

#ifndef __DBL_CONVOLUTIONS_KERNELS_H__
#define __DBL_CONVOLUTIONS_KERNELS_H__

#define d_shmemsize 512

/* The constant d_shmemsize is the bound on the shared memory size,
 * to compute the product of series with complex double coefficients.
 * For degree 511, we have 512 complex numbers, for x, y, and z,
 * real and imaginary parts, so 512*3*2*8 = 24576 bytes.
 * This constants bounds the degree of the power series. */

__global__ void dbl_convolute ( double *x, double *y, double *z, int dim );
/*
 * DESCRIPTION :
 *   Convolutes real vectors x and y to make z.
 *   All vectors have dimension dim. */

__global__ void cmplx_convolute
 ( double *xre, double *xim, double *yre, double *yim,
   double *zre, double *zim, int dim );
/*
 * DESCRIPTION :
 *   Convolutes complex vectors x and y to make z.
 *   The complex vectors are given as double vectors with real
 *   and imaginary parts.  All vectors have dimension dim.
 *
 * ON ENTRY :
 *   xre      real parts of the first vector x;
 *   xim      imaginary parts of the first vector x;
 *   yre      real parts of the second vector y;
 *   yim      imaginary parts of the second vector y.
 *
 * ON RETURN :
 *   zre      real parts of the product of x and y;
 *   zim      imaginary part of the product of x and y. */

void GPU_dbl_product
 ( double *x_h, double *y_h, double *z_h, int deg, int freq, int BS );
/*
 * DESCRIPTION :
 *   Computes the product of two power series with double coefficients.
 *
 * ON ENTRY :
 *   x_h      deg+1 coefficients of the first series;
 *   y_h      deg+1 coefficients of the second series;
 *   z_h      space allocated for deg+1 doubles;
 *   deg      degree of the truncated power series;
 *   freq     frequency for timing purposes;
 *   BS       block size, the number of threads in a block.
 *
 * ON RETURN :
 *   z_h      coefficients of the product of x with y. */

void GPU_cmplx_product
 ( double *xre_h, double *xim_h, double *yre_h, double *yim_h,
   double *zre_h, double *zim_h, int deg, int freq, int BS );
/*
 * DESCRIPTION :
 *   Computes the product of two power series
 *   with complex double coefficients.
 *
 * ON ENTRY :
 *   xre_h    deg+1 real parts of the coefficients of the first series;
 *   xim_h    deg+1 imag parts of the coefficients of the first series;
 *   yre_h    deg+1 real parts of the coefficients of the second series;
 *   yim_h    deg+1 imag parts of the coefficients of the second series;
 *   zre_h    space allocated for deg+1 doubles;
 *   zim_h    space allocated for deg+1 doubles;
 *   deg      degree of the truncated power series;
 *   freq     frequency for timing purposes;
 *   BS       block size, the number of threads in a block.
 *
 * ON RETURN :
 *   zre_h    real parts of the coefficients of the product of x with y;
 *   zim_h    imag parts of the coefficients of the product of x with y. */

#endif
