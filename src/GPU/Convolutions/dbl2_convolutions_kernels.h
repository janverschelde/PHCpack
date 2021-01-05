// The file dbl2_convolutions_kernels.h defines prototypes for the kernels 
// to multiply two truncated power series in double double precision.
// Defines the constant bound on the shared memory size of one block.

#ifndef __DBL2_CONVOLUTIONS_KERNELS_H__
#define __DBL2_CONVOLUTIONS_KERNELS_H__

// #define dd_shmemsize 256
// redefined for C2050
#define dd_shmemsize 192

/* The constant dd_shmemsize is the bound on the shared memory size,
 * to compute the product of series with complex double double coefficients.
 * For degree 255, we have 256 complex numbers, for x, y, and z,
 * real and imaginary parts, high and low parts, so 256*2*3*2*8 = 24576 bytes.
 * This constant bounds the degree of the power series. */

__global__ void dbl2_increment
 ( double *xhi, double *xlo, double *yhi, double *ylo,
   double *zhi, double *zlo, int dim );
/*
 * DESCRIPTION : z = x + y.
 *   Add y to x to make z.  All arrays are of dimension dim. */

__global__ void dbl2_decrement
 ( double *xhi, double *xlo, double *yhi, double *ylo,
   double *zhi, double *zlo, int dim );
/*
 * DESCRIPTION : z = x + y.
 *  Subtracts y from x to make z.  All arrays are of dimension dim. */

__global__ void dbl2_convolute
 ( double *xhi, double *xlo, double *yhi, double *ylo,
   double *zhi, double *zlo, int dim );
/*
 * DESCRIPTION :
 *   Convolutes real vectors x and y to make z.
 *   All arrays have dimension dim.
 *
 * ON ENTRY :
 *   xhi      high parts of the first vector x;
 *   xlo      low parts of the first vector x;
 *   yhi      high parts of the second vector y;
 *   ylo      low parts of the second vector y;
 *   zhi      dim doubles allocated for the high parts of the product;
 *   zlo      dim doubles allocated for the low parts of the product;
 *   dim      dimension of all arrays.
 *
 * ON RETURN :
 *   zhi      high parts of the product of x with y;
 *   zlo      low parts of the product of x with y. */

__global__ void dbl2_padded_convolute
 ( double *xhi, double *xlo, double *yhi, double *ylo,
   double *zhi, double *zlo, int dim );
/*
 * DESCRIPTION :
 *   Convolutes real vectors x and y to make z.
 *   All arrays have dimension dim.
 *   Zeros are inserted in the shared memory for the second vector,
 *   so all threads perform the same number of operations
 *   and thread divergence is avoided. */

__global__ void cmplx2_convolute
 ( double *xrehi, double *xrelo, double *ximhi, double *ximlo,
   double *yrehi, double *yrelo, double *yimhi, double *yimlo,
   double *zrehi, double *zrelo, double *zimhi, double *zimlo, int dim );
/*
 * DESCRIPTION :
 *   Convolutes complex vectors x and y to make z.
 *   The complex vectors are given as double arrays with real
 *   and imaginary parts, of high and low parts.
 *   All arrays have dimension dim.
 *
 * ON ENTRY :
 *   xrehi    high parts of the real parts of the coefficients of x;
 *   ximlo    low parts of the imaginary parts of the coefficients of x
 *   yrehi    high parts of the real parts of the coefficients of y;
 *   yimlo    low parts of the imaginary parts of the coefficients of y.
 *   zrehi    dim doubles allocated for the high parts
 *            of the real parts of the coefficients of the product;
 *   zrelo    dim doubles allocated for the low parts
 *            of the real parts of the coefficients of the product;
 *   zimhi    dim doubles allocated for the high parts
 *            of the imaginary parts of the coefficients of the product;
 *   zimlo    dim doubles allocated for the low parts
 *            of the imaginary parts of the coefficients of the product;
 *   dim      dimension of all arrays.
 *
 * ON RETURN :
 *   zrehi    high parts of the real parts of the product of x and y;
 *   zrelo    low parts of the real parts of the product of x and y;
 *   zimhi    high parts of the imaginary parts of the product of x and y;
 *   zimlo    low parts of the imaginary parts of the product of x and y. */

__global__ void cmplx2_looped_convolute
 ( double *xrehi, double *xrelo, double *ximhi, double *ximlo,
   double *yrehi, double *yrelo, double *yimhi, double *yimlo,
   double *zrehi, double *zrelo, double *zimhi, double *zimlo, int dim );
/*
 * DESCRIPTION :
 *   Does the same as cmplx2_convolute, but in two loops,
 *   one for the real and one for the imaginary parts.
 *   The loop for the imaginary parts has the indices reversed,
 *   so every thread executes the same amount of operations. */

__global__ void cmplx2_padded_convolute
 ( double *xrehi, double *xrelo, double *ximhi, double *ximlo,
   double *yrehi, double *yrelo, double *yimhi, double *yimlo,
   double *zrehi, double *zrelo, double *zimhi, double *zimlo, int dim );
/*
 * DESCRIPTION :
 *   Convolutes complex vectors x and y to make z.
 *   The complex vectors are given as double arrays with real
 *   and imaginary parts, of high and low parts.
 *   All arrays have dimension dim.
 *   Zeros are inserted in the shared memory for the second vector,
 *   so all threads perform the same number of operations
 *   and thread divergence is avoided. */

void GPU_dbl2_product
 ( double *xhi_h, double *xlo_h, double *yhi_h, double *ylo_h,
   double *zhi_h, double *zlo_h, int deg, int freq, int BS, int padded );
/*
 * DESCRIPTION :
 *   Computes the product of two power series x and y to make z,
 *   with double double coefficients.
 *
 * ON ENTRY :
 *   xhi_h    deg+1 high parts of the coefficients of x;
 *   xlo_h    deg+1 low parts of the coefficients of x;
 *   yhi_h    deg+1 high parts of the coefficients of y;
 *   ylo_h    deg+1 low parts of the coefficients of y;
 *   zhi_h    space allocated for deg+1 doubles for the high parts of z;
 *   zlo_h    space allocated for deg+1 doubles for the low parts of z;
 *   deg      degree of the truncated power series;
 *   freq     frequency for timing purposes;
 *   BS       block size, the number of threads in a block;
 *   padded   if 1, then the padded convolute is applied.
 *
 * ON RETURN :
 *   zhi_h    high parts of the coefficients of the product of x with y;
 *   zlo_h    low parts of the coefficients of the product of x with y. */

void GPU_cmplx2_product
 ( double *xrehi_h, double *xrelo_h, double *ximhi_h, double *ximlo_h,
   double *yrehi_h, double *yrelo_h, double *yimhi_h, double *yimlo_h,
   double *zrehi_h, double *zrelo_h, double *zimhi_h, double *zimlo_h,
   int deg, int freq, int BS, int mode );
/*
 * DESCRIPTION :
 *   Computes the product of two power series x and y to make z,
 *   with complex double double coefficients.
 *
 * ON ENTRY :
 *   xrehi_h    deg+1 high parts of the real parts of x;
 *   xrelo_h    deg+1 low parts of the real parts of x;
 *   ximhi_h    deg+1 high parts of the imag parts of x;
 *   ximlo_h    deg+1 low parts of the imag parts of x;
 *   yrehi_h    deg+1 high parts of the real parts of y;
 *   yrelo_h    deg+1 low parts of the real parts of y;
 *   yimhi_h    deg+1 high parts of the imag parts of y;
 *   yimlo_h    deg+1 low parts of the imag parts of y;
 *   zrehi_h    space for deg+1 doubles for the high real parts of z;
 *   zrelo_h    space for deg+1 doubles for the low real parts of z;
 *   zimhi_h    space for deg+1 doubles for the high imag parts of z;
 *   zimlo_h    space for deg+1 doubles for the low imag parts of z;
 *   deg        degree of the truncated power series;
 *   freq       frequency for timing purposes;
 *   BS         block size, the number of threads in a block;
 *   mode       if 1, then the looped convolute is applied.
 *
 * ON RETURN :
 *   zrehi_h    high parts of the real parts of the coefficients of z;
 *   zrelo_h    low parts of the real parts of the coefficients of z;
 *   zimhi_h    high parts of the imag parts of the coefficients of z;
 *   zimlo_h    low parts of the imag parts of the coefficients of z. */

#endif
