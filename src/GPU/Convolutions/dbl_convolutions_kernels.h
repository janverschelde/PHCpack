// The file dbl_convolutions_kernels.h defines prototypes for the kernels 
// to multiply two truncated power series in double precision.
// Defines the constant bound on the shared memory size of one block.

#ifndef __DBL_CONVOLUTIONS_KERNELS_H__
#define __DBL_CONVOLUTIONS_KERNELS_H__

// #define d_shmemsize 512
// Reset to smaller size for C2050
#define d_shmemsize 256

/* The constant d_shmemsize is the bound on the shared memory size,
 * to compute the product of series with complex double coefficients.
 * For degree 511, we have 512 complex numbers, for x, y, and z,
 * real and imaginary parts, so 512*3*2*8 = 24576 bytes.
 * This constant bounds the degree of the power series. */

__global__ void dbl_increment ( double *x, double *y, double *z, int dim );
/*
 * DESCRIPTION : z = x + y.
 *   Adds y to x to make z.  All arrays are of dimension dim. */

__global__ void dbl_decrement ( double *x, double *y, double *z, int dim );
/*
 * DESCRIPTION : z = x - y.
 *   Subtracts y from x to make z.  All arrays are of dimension dim. */

__global__ void dbl_convolute ( double *x, double *y, double *z, int dim );
/*
 * DESCRIPTION :
 *   Convolutes real vectors x and y to make z.
 *   All arrays have dimension dim. */

__global__ void dbl_padded_convolute
 ( double *x, double *y, double *z, int dim );
/*
 * DESCRIPTION :
 *   Convolutes real vectors x and y to make z.
 *   All arrays have dimension dim.
 *   Zeros are inserted in the shared memory for the second vector,
 *   so all threads perform the same number of operations
 *   and thread divergence is avoided. */

__global__ void cmplx_convolute
 ( double *xre, double *xim, double *yre, double *yim,
   double *zre, double *zim, int dim );
/*
 * DESCRIPTION :
 *   Convolutes complex vectors x and y to make z.
 *   The complex vectors are given as double vectors with real
 *   and imaginary parts.  All arrays have dimension dim.
 *
 * ON ENTRY :
 *   xre      real parts of the first vector x;
 *   xim      imaginary parts of the first vector x;
 *   yre      real parts of the second vector y;
 *   yim      imaginary parts of the second vector y;
 *   zre      dim doubles allocated for the real parts of the product;
 *   zim      dim doubles allocated for the imaginary parts of the product;
 *   dim      dimension of all arrays.
 *
 * ON RETURN :
 *   zre      real parts of the product of x and y;
 *   zim      imaginary parts of the product of x and y. */

__global__ void cmplx_padded_convolute
 ( double *xre, double *xim, double *yre, double *yim,
   double *zre, double *zim, int dim );
/*
 * DESCRIPTION :
 *   Convolutes complex vectors x and y to make z.
 *   All arrays have dimension dim.
 *   Zeros are inserted in the shared memory for the second vector,
 *   so all threads perform the same number of operations
 *   and thread divergence is avoided. */

__global__ void cmplx_looped_convolute
 ( double *xre, double *xim, double *yre, double *yim,
   double *zre, double *zim, int dim );
/*
 * DESCRIPTION :
 *   Does the same as cmplx_convolute, but in two loops,
 *   one for the real and one for the imaginary parts.
 *   The loop for the imaginary parts has the indices reversed,
 *   so every thread executes the same amount of operations. */

void GPU_dbl_product
 ( double *x_h, double *y_h, double *z_h, int deg, int freq, int BS,
   int mode );
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
 *   BS       block size, the number of threads in a block;
 *   padded   if 1, then the padded convolution is called.
 *
 * ON RETURN :
 *   z_h      coefficients of the product of x with y. */

void GPU_cmplx_product
 ( double *xre_h, double *xim_h, double *yre_h, double *yim_h,
   double *zre_h, double *zim_h, int deg, int freq, int BS, int looped );
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
 *   BS       block size, the number of threads in a block;
 *   mode     0 : the plain version,
 *            1 : looped convolute,
 *            2 : multiple looped convolute,
 *            3 : padded version.
 *
 * ON RETURN :
 *   zre_h    real parts of the coefficients of the product of x with y;
 *   zim_h    imag parts of the coefficients of the product of x with y. */

#endif
