// The file dbl3_convolutions_kernels.h defines prototypes for the kernels 
// to multiply two truncated power series in triple double precision.
// Defines the constant bound on the shared memory size of one block.

#ifndef __DBL3_CONVOLUTIONS_KERNELS_H__
#define __DBL3_CONVOLUTIONS_KERNELS_H__

#define td_shmemsize 192

/* The constant td_shmemsize is the bound on the shared memory size,
 * to compute the product of series with complex triple double coefficients.
 * For degree 191, we have 192 complex numbers, for x, y, and z,
 * real and imaginary parts, high and low parts, so 192*2*3*3*8 = 27648 bytes.
 * This constant bounds the degree of the power series. */

__global__ void dbl3_increment
 ( double *xhi, double *xmi, double *xlo,
   double *yhi, double *ymi, double *ylo,
   double *zhi, double *zmi, double *zlo, int dim );
/*
 * DESCRIPTION : z = x + y.
 *   Adds y to x to make z.  All arrays are of dimension dim. */

__global__ void dbl3_decrement
 ( double *xhi, double *xmi, double *xlo,
   double *yhi, double *ymi, double *ylo,
   double *zhi, double *zmi, double *zlo, int dim );
/*
 * DESCRIPTION : z = x - y.
 *   Subtracts y from x to make z.  All arrays are of dimension dim. */
 
__global__ void dbl3_convolute
 ( double *xhi, double *xmi, double *xlo,
   double *yhi, double *ymi, double *ylo,
   double *zhi, double *zmi, double *zlo, int dim );
/*
 * DESCRIPTION :
 *   Convolutes real vectors x and y to make z.
 *   All arrays have dimension dim.
 *
 * ON ENTRY :
 *   xhi      high parts of the first vector x;
 *   xmi      middle parts of the first vector x;
 *   xlo      low parts of the first vector x;
 *   yhi      high parts of the second vector y;
 *   ymi      middle parts of the second vector y;
 *   ylo      low parts of the second vector y;
 *   zhi      dim doubles allocated for the high parts of the product;
 *   zmi      dim doubles allocated for the middle parts of the product;
 *   zlo      dim doubles allocated for the low parts of the product;
 *   dim      dimension of all arrays.
 *
 * ON RETURN :
 *   zhi      high parts of the product of x with y;
 *   zmi      middle parts of the product of x with y;
 *   zlo      low parts of the product of x with y. */

__global__ void dbl3_padded_convolute
 ( double *xhi, double *xmi, double *xlo,
   double *yhi, double *ymi, double *ylo,
   double *zhi, double *zmi, double *zlo, int dim );
/*
 * DESCRIPTION :
 *   Convolutes real vectors x and y to make z.
 *   All arrays have dimension dim.
 *   Zeros are inserted in the shared memory for the second vector,
 *   so all threads perform the same number of operations
 *   and thread divergence is avoided. */

__global__ void cmplx3_convolute
 ( double *xrehi, double *xremi, double *xrelo,
   double *ximhi, double *ximmi, double *ximlo,
   double *yrehi, double *yremi, double *yrelo,
   double *yimhi, double *yimmi, double *yimlo,
   double *zrehi, double *zremi, double *zrelo,
   double *zimhi, double *zimmi, double *zimlo, int dim );
/*
 * DESCRIPTION :
 *   Convolutes complex vectors x and y to make z.
 *   The complex vectors are given as double arrays with real
 *   and imaginary parts, of high, middle, and low parts.
 *   All arrays have dimension dim.
 *
 * ON ENTRY :
 *   xrehi    high parts of the real parts of the coefficients of x;
 *   xremi    middle parts of the real parts of the coefficients of x;
 *   ximlo    low parts of the imaginary parts of the coefficients of x
 *   yrehi    high parts of the real parts of the coefficients of y;
 *   yremi    middle parts of the real parts of the coefficients of y;
 *   yimlo    low parts of the imaginary parts of the coefficients of y;
 *   zrehi    dim doubles allocated for the highest parts
 *            of the real parts of the coefficients of the product;
 *   zremi    dim doubles allocated for the middle parts
 *            of the real parts of the coefficients of the product;
 *   zrelo    dim doubles allocated for the lowest parts
 *            of the real parts of the coefficients of the product;
 *   zimhi    dim doubles allocated for the highest parts
 *            of the imaginary parts of the coefficients of the product;
 *   zimmi    dim doubles allocated for the middle parts
 *            of the imaginary parts of the coefficients of the product;
 *   zimlo    dim doubles allocated for the lowest parts
 *            of the imaginary parts of the coefficients of the product;
 *   dim      dimension of all arrays.
 *
 * ON RETURN :
 *   zrehi    high parts of the real parts of the product z;
 *   zremi    middle parts of the real parts of the product z;
 *   zrelo    low parts of the real parts of the product z;
 *   zimhi    high parts of the imaginary parts of the product z;
 *   zimmi    middle parts of the imaginary parts of the product z;
 *   zimlo    low parts of the imaginary parts of the product z. */

__global__ void cmplx3_looped_convolute
 ( double *xrehi, double *xremi, double *xrelo,
   double *ximhi, double *ximmi, double *ximlo,
   double *yrehi, double *yremi, double *yrelo,
   double *yimhi, double *yimmi, double *yimlo,
   double *zrehi, double *zremi, double *zrelo,
   double *zimhi, double *zimmi, double *zimlo, int dim );
/*
 * DESCRIPTION :
 *   Does the same as cmplx3_convolute, but in two loops,
 *   one for the real and one for the imaginary parts.
 *   The loop for the imaginary parts has the indices reversed,
 *   so every thread executes the same amount of operations. */

__global__ void cmplx3_padded_convolute
 ( double *xrehi, double *xremi, double *xrelo,
   double *ximhi, double *ximmi, double *ximlo,
   double *yrehi, double *yremi, double *yrelo,
   double *yimhi, double *yimmi, double *yimlo,
   double *zrehi, double *zremi, double *zrelo,
   double *zimhi, double *zimmi, double *zimlo, int dim );
/*
 * DESCRIPTION :
 *   Does the same as cmplx3_convolute, but with padding.
 *   Zeros are inserted in the shared memory for the second vector,
 *   so all threads perform the same number of operations
 *   and thread divergence is avoided. */

void GPU_dbl3_product
 ( double *xhi_h, double *xmi_h, double *xlo_h,
   double *yhi_h, double *ymi_h, double *ylo_h,
   double *zhi_h, double *zmi_h, double *zlo_h, int deg, int freq, int BS,
   int padded );
/*
 * DESCRIPTION :
 *   Computes the product of two power series x and y to make z,
 *   with triple double coefficients.
 *
 * ON ENTRY :
 *   xhi_h    deg+1 high parts of the coefficients of x;
 *   xmi_h    deg+1 middle parts of the coefficients of x;
 *   xlo_h    deg+1 low parts of the coefficients of x;
 *   yhi_h    deg+1 high parts of the coefficients of y;
 *   ymi_h    deg+1 middle parts of the coefficients of y;
 *   ylo_h    deg+1 low parts of the coefficients of y;
 *   zhi_h    space allocated for deg+1 doubles for the high parts of z;
 *   zmi_h    space allocated for deg+1 doubles for the middle parts of z;
 *   zlo_h    space allocated for deg+1 doubles for the low parts of z;
 *   deg      degree of the truncated power series;
 *   freq     frequency for timing purposes;
 *   BS       block size, the number of threads in a block;
 *   padded   if 1, then the padded version is applied.
 *
 * ON RETURN :
 *   zhi_h    high parts of the coefficients of the product of x with y;
 *   zmi_h    middle parts of the coefficients of the product of x with y;
 *   zlo_h    low parts of the coefficients of the product of x with y. */

void GPU_cmplx3_product
 ( double *xrehi_h, double *xremi_h, double *xrelo_h,
   double *ximhi_h, double *ximmi_h, double *ximlo_h,
   double *yrehi_h, double *yremi_h, double *yrelo_h,
   double *yimhi_h, double *yimmi_h, double *yimlo_h,
   double *zrehi_h, double *zremi_h, double *zrelo_h,
   double *zimhi_h, double *zimmi_h, double *zimlo_h,
   int deg, int freq, int BS, int mode );
/*
 * DESCRIPTION :
 *   Computes the product of two power series x and y to make z,
 *   with complex triple double coefficients.
 *
 * ON ENTRY :
 *   xrehi_h    deg+1 high parts of the real parts of x;
 *   xremi_h    deg+1 middle parts of the real parts of x;
 *   xrelo_h    deg+1 low parts of the real parts of x;
 *   ximhi_h    deg+1 high parts of the imag parts of x;
 *   ximmi_h    deg+1 middle parts of the imag parts of x;
 *   ximlo_h    deg+1 low parts of the imag parts of x;
 *   yrehi_h    deg+1 high parts of the real parts of y;
 *   yremi_h    deg+1 middle parts of the real parts of y;
 *   yrelo_h    deg+1 low parts of the real parts of y;
 *   yimhi_h    deg+1 high parts of the imag parts of y;
 *   yimmi_h    deg+1 middle parts of the imag parts of y;
 *   yimlo_h    deg+1 low parts of the imag parts of y;
 *   zrehi_h    space for deg+1 doubles for the high real parts of z;
 *   zremi_h    space for deg+1 doubles for the middle real parts of z;
 *   zrelo_h    space for deg+1 doubles for the low real parts of z;
 *   zimhi_h    space for deg+1 doubles for the high imag parts of z;
 *   zimmi_h    space for deg+1 doubles for the middle imag parts of z;
 *   zimlo_h    space for deg+1 doubles for the low imag parts of z;
 *   deg        degree of the truncated power series;
 *   freq       frequency for timing purposes;
 *   BS         block size, the number of threads in a block;
 *   mode       0 : plain version,
 *              1 : looped convolute,
 *              2 : multiple real convolutes
 *              3 : padded convolute.
 *
 * ON RETURN :
 *   zrehi_h    high parts of the real parts of the coefficients of z;
 *   zremi_h    middle parts of the real parts of the coefficients of z;
 *   zrelo_h    low parts of the real parts of the coefficients of z;
 *   zimhi_h    high parts of the imag parts of the coefficients of z;
 *   zimmi_h    middle parts of the imag parts of the coefficients of z;
 *   zimlo_h    low parts of the imag parts of the coefficients of z. */

#endif
