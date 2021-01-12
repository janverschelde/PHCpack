// The file dbl4_convolutions_kernels.h defines prototypes for the kernels 
// to multiply two truncated power series in quad double precision.
// Defines the constant bound on the shared memory size of one block.

#ifndef __DBL4_CONVOLUTIONS_KERNELS_H__
#define __DBL4_CONVOLUTIONS_KERNELS_H__

#define qd_shmemsize 128

/* The constant qd_shmemsize is the bound on the shared memory size,
 * to compute the product of series with complex quad double coefficients.
 * For degree 127, we have 128 complex numbers, for x, y, and z,
 * real and imaginary parts, highest, second highest, second lowest,
 * and lowest parts, so 128*2*3*4*8 = 24576 bytes.
 * This constant bounds the degree of the power series. */

__global__ void dbl4_increment
 ( double *xhihi, double *xlohi, double *xhilo, double *xlolo,
   double *yhihi, double *ylohi, double *yhilo, double *ylolo,
   double *zhihi, double *zlohi, double *zhilo, double *zlolo, int dim );
/*
 * DESCRIPTION : z = x + y.
 *   Adds y to x to make z.  All arrays have dimension dim. */

__global__ void dbl4_decrement
 ( double *xhihi, double *xlohi, double *xhilo, double *xlolo,
   double *yhihi, double *ylohi, double *yhilo, double *ylolo,
   double *zhihi, double *zlohi, double *zhilo, double *zlolo, int dim );
/*
 * DESCRIPTION : z = x + y.
 *   Subtracts y from x to make z.  All arrays have dimension dim. */

__global__ void dbl4_convolute
 ( double *xhihi, double *xlohi, double *xhilo, double *xlolo,
   double *yhihi, double *ylohi, double *yhilo, double *ylolo,
   double *zhihi, double *zlohi, double *zhilo, double *zlolo, int dim );
/*
 * DESCRIPTION :
 *   Convolutes real vectors x and y to make z.
 *   All arrays have dimension dim.
 *
 * ON ENTRY :
 *   xhihi    highest parts of the first vector x;
 *   xlohi    second highest parts of the first vector x;
 *   xhilo    second lowest parts of the first vector x;
 *   xlolo    lowest parts of the first vector x;
 *   yhihi    highest parts of the second vector y;
 *   ylohi    second highest parts of the second vector y;
 *   yhilo    second lowest parts of the second vector y;
 *   ylolo    lowest parts of the second vector y;
 *   zhihi    dim doubles allocated for the highest parts of z;
 *   zlohi    dim doubles allocated for the second highest parts of z;
 *   zhilo    dim doubles allocated for the second lowest parts of z;
 *   zlolo    dim doubles allocated for the lowest parts of z;
 *   dim      dimension of all arrays.
 *
 * ON RETURN :
 *   zhihi    highest parts of the product of x with y;
 *   zlohi    second highest parts of the product of x with y;
 *   zhilo    second lowest parts of the product of x with y;
 *   zlolo    lowest parts of the product of x with y. */

__global__ void dbl4_padded_convolute
 ( double *xhihi, double *xlohi, double *xhilo, double *xlolo,
   double *yhihi, double *ylohi, double *yhilo, double *ylolo,
   double *zhihi, double *zlohi, double *zhilo, double *zlolo, int dim );
/*
 * DESCRIPTION :
 *   Convolutes real vectors x and y to make z.
 *   All arrays have dimension dim.
 *   Zeros are inserted in the shared memory for the second vector,
 *   so all threads perform the same number of operations
 *   and thread divergence is avoided. */

__global__ void cmplx4_convolute
 ( double *xrehihi, double *xrelohi, double *xrehilo, double *xrelolo,
   double *ximhihi, double *ximlohi, double *ximhilo, double *ximlolo,
   double *yrehihi, double *yrelohi, double *yrehilo, double *yrelolo,
   double *yimhihi, double *yimlohi, double *yimhilo, double *yimlolo,
   double *zrehihi, double *zrelohi, double *zrehilo, double *zrelolo,
   double *zimhihi, double *zimlohi, double *zimhilo, double *zimlolo,
   int dim );
/*
 * DESCRIPTION :
 *   Convolutes complex vectors x and y to make z.
 *   The complex vectors are given as double arrays with real
 *   and imaginary parts, of highest, second highest, second lowest,
 *   and lowest parts.  All arrays have dimension dim.
 *
 * ON ENTRY :
 *   xrehihi    highest parts of the real parts of x;
 *   xrelohi    second highest parts of the real parts of x;
 *   xrehilo    second lowest parts of the real parts of x;
 *   xrelolo    lowest parts of the real parts of x;
 *   ximhihi    highest parts of the imaginary parts of x;
 *   ximlohi    second highest parts of the imaginary parts of x;
 *   ximhilo    second lowest parts of the imaginary parts of x;
 *   ximlolo    lowest parts of the imaginary parts of x;
 *   yrehihi    highest parts of the real parts of y;
 *   yrelohi    second highest parts of the real parts of y;
 *   yrehilo    second lowest parts of the real parts of y;
 *   yrelolo    lowest parts of the real parts of y;
 *   yimhihi    highest parts of the imaginary parts of y;
 *   yimlohi    second highest parts of the imaginary parts of y;
 *   yimhilo    second lowest parts of the imaginary parts of y;
 *   yimlolo    lowest parts of the imaginary parts of y;
 *   zrehihi    dim doubles allocated for the highest parts
 *              of the real parts of the coefficients of the product;
 *   zrelohi    dim doubles allocated for the second highest parts
 *              of the real parts of the coefficients of the product;
 *   zrehilo    dim doubles allocated for the second lowest parts
 *              of the real parts of the coefficients of the product;
 *   zrelolo    dim doubles allocated for the lowest parts
 *              of the real parts of the coefficients of the product;
 *   zimhihi    dim doubles allocated for the highest parts
 *              of the imaginary parts of the coefficients of the product;
 *   zimlohi    dim doubles allocated for the second highest parts
 *              of the imaginary parts of the coefficients of the product;
 *   zimhilo    dim doubles allocated for the second lowest parts
 *              of the imaginary parts of the coefficients of the product;
 *   zimlolo    dim doubles allocated for the lowest parts
 *              of the imaginary parts of the coefficients of the product;
 *   dim        dimension of all arrays.
 *
 * ON RETURN :
 *   zrehihi    highest parts of the real parts of the product z;
 *   zrelohi    second highest parts of the real parts of the product z;
 *   zrehilo    second lowest parts of the real parts of the product z;
 *   zrelolo    lowest parts of the real parts of the product z;
 *   zimhihi    highest parts of the imaginary parts of the product z;
 *   zimlohi    second highest parts of the imaginary parts of z;
 *   zimhilo    second lowest parts of the imaginary parts of z;
 *   zimlolo    lowest parts of the imaginary parts of the product z. */

__global__ void cmplx4_padded_convolute
 ( double *xrehihi, double *xrelohi, double *xrehilo, double *xrelolo,
   double *ximhihi, double *ximlohi, double *ximhilo, double *ximlolo,
   double *yrehihi, double *yrelohi, double *yrehilo, double *yrelolo,
   double *yimhihi, double *yimlohi, double *yimhilo, double *yimlolo,
   double *zrehihi, double *zrelohi, double *zrehilo, double *zrelolo,
   double *zimhihi, double *zimlohi, double *zimhilo, double *zimlolo,
   int dim );
/*
 * DESCRIPTION :
 *   Convolutes complex vectors x and y to make z.
 *   The complex vectors are given as double arrays with real
 *   and imaginary parts, of highest, second highest, second lowest,
 *   and lowest parts.  All arrays have dimension dim.
 *   Zeros are inserted in the shared memory for the second vector,
 *   so all threads perform the same number of operations
 *   and thread divergence is avoided. */

void GPU_dbl4_product
 ( double *xhihi_h, double *xlohi_h, double *xhilo_h, double *xlolo_h,
   double *yhihi_h, double *ylohi_h, double *yhilo_h, double *ylolo_h,
   double *zhihi_h, double *zlohi_h, double *zhilo_h, double *zlolo_h,
   int deg, int freq, int BS, int padded );
/*
 * DESCRIPTION :
 *   Computes the product of two power series x and y to make z,
 *   with quad double coefficients.
 *
 * ON ENTRY :
 *   xhihi_h    deg+1 highest parts of the coefficients of x;
 *   xlohi_h    deg+1 second highest parts of the coefficients of x;
 *   xhilo_h    deg+1 second lowest parts of the coefficients of x;
 *   xlolo_h    deg+1 lowest parts of the coefficients of x;
 *   yhihi_h    deg+1 highest parts of the coefficients of y;
 *   ylohi_h    deg+1 second highest low parts of the coefficients of y;
 *   yhilo_h    deg+1 second lowest parts of the coefficients of y;
 *   ylolo_h    deg+1 lowest parts of the coefficients of y;
 *   zhihi_h    space for deg+1 doubles for the highest parts of z;
 *   zlohi_h    space for deg+1 doubles for the second highest parts of z;
 *   zhilo_h    space for deg+1 doubles for the second lowest parts of z;
 *   zlolo_h    space for deg+1 doubles for the lowest parts of z;
 *   deg        degree of the truncated power series;
 *   freq       frequency for timing purposes;
 *   BS         block size, the number of threads in a block;
 *   padded     if 1, then the padded convolute is applied.
 *
 * ON RETURN :
 *   zhihi_h    highest parts of the coefficients of the product z;
 *   zlohi_h    second highest parts of the coefficients of the product z;
 *   zhilo_h    second lowest parts of the coefficients of the product z;
 *   zlolo_h    lowest parts of the coefficients of the product z. */

void GPU_cmplx4_product
 ( double *xrehihi_h, double *xrelohi_h, double *xrehilo_h, double *xrelolo_h,
   double *ximhihi_h, double *ximlohi_h, double *ximhilo_h, double *ximlolo_h,
   double *yrehihi_h, double *yrelohi_h, double *yrehilo_h, double *yrelolo_h,
   double *yimhihi_h, double *yimlohi_h, double *yimhilo_h, double *yimlolo_h,
   double *zrehihi_h, double *zrelohi_h, double *zrehilo_h, double *zrelolo_h,
   double *zimhihi_h, double *zimlohi_h, double *zimhilo_h, double *zimlolo_h,
   int deg, int freq, int BS, int mode );
/*
 * DESCRIPTION :
 *   Computes the product of two power series x and y to make z,
 *   with complex quad double coefficients.
 *
 * ON ENTRY :
 *   xrehihi_h    deg+1 highest parts of the real parts of x;
 *   xrelohi_h    deg+1 second highest parts of the real parts of x;
 *   xrehilo_h    deg+1 second lowest parts of the real parts of x;
 *   xrelolo_h    deg+1 lowest parts of the real parts of x;
 *   ximhihi_h    deg+1 highest parts of the imag parts of x;
 *   ximlohi_h    deg+1 second highest parts of the imag parts of x;
 *   ximhilo_h    deg+1 second lowest parts of the imag parts of x;
 *   ximlolo_h    deg+1 lowest parts of the imag parts of x;
 *   yrehihi_h    deg+1 highest parts of the real parts of y;
 *   yrelohi_h    deg+1 second highest parts of the real parts of y;
 *   yrehilo_h    deg+1 second lowest parts of the real parts of y;
 *   yrelolo_h    deg+1 lowest parts of the real parts of y;
 *   yimhihi_h    deg+1 highest parts of the imag parts of y;
 *   yimlohi_h    deg+1 second highest parts of the imag parts of y;
 *   yimhilo_h    deg+1 second lowest parts of the imag parts of y;
 *   yimlolo_h    deg+1 lowest parts of the imag parts of y;
 *   zrehihi_h    space for deg+1 doubles for the highest real parts of z;
 *   zrelohi_h    space for deg+1 doubles for the 2nd highest real parts of z;
 *   zrehilo_h    space for deg+1 doubles for the 2nd lowest real parts of z;
 *   zrelolo_h    space for deg+1 doubles for the lowest real parts of z;
 *   zimhihi_h    space for deg+1 doubles for the highest imag parts of z;
 *   zimlohi_h    space for deg+1 doubles for the 2nd highest imag parts of z;
 *   zimhilo_h    space for deg+1 doubles for the 2nd lowest imag parts of z;
 *   zimlolo_h    space for deg+1 doubles for the lowest imag parts of z;
 *   deg          degree of the truncated power series;
 *   freq         frequency for timing purposes;
 *   BS           block size, the number of threads in a block;
 *   mode         if 1, then multiple double kernels are launched;
 *                if 2, then one padded complex convolute is applied.
 *
 * ON RETURN :
 *   zrehihi_h    highest parts of the real parts of z;
 *   zrelohi_h    second highest parts of the real parts of z;
 *   zrehilo_h    second lowest parts of the real parts of z;
 *   zrelolo_h    lowest parts of the real parts of z;
 *   zimhihi_h    highest parts of the imag parts of z;
 *   zimlohi_h    second highest parts of the imag parts of z;
 *   zimhilo_h    second lowest parts of the imag parts of z;
 *   zimlolo_h    lowest parts of the imag parts of z. */

#endif
