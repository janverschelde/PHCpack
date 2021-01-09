// The file dbl8_convolutions_kernels.h defines prototypes for the kernels 
// to multiply two truncated power series in octo double precision.
// Defines the constant bound on the shared memory size of one block.

#ifndef __DBL8_CONVOLUTIONS_KERNELS_H__
#define __DBL8_CONVOLUTIONS_KERNELS_H__

#define od_shmemsize 64

/* The constant od_shmemsize is the bound on the shared memory size,
 * to compute the product of series with complex octo double coefficients.
 * For degree 63, we have 64 complex numbers, for x, y, and z,
 * real and imaginary parts, eight doubles for every octo double number,
 * so 64*2*3*8*8 = 24576 bytes.
 * This constant bounds the degree of the power series. */

__global__ void dbl8_increment
 ( double *xhihihi, double *xlohihi, double *xhilohi, double *xlolohi,
   double *xhihilo, double *xlohilo, double *xhilolo, double *xlololo,
   double *yhihihi, double *ylohihi, double *yhilohi, double *ylolohi,
   double *yhihilo, double *ylohilo, double *yhilolo, double *ylololo,
   double *zhihihi, double *zlohihi, double *zhilohi, double *zlolohi,
   double *zhihilo, double *zlohilo, double *zhilolo, double *zlololo,
   int dim );
/*
 * DESCRIPTION : z = x + y.
 *   Adds y to x to make z.  All arrays are of dimension dim. */

__global__ void dbl8_decrement
 ( double *xhihihi, double *xlohihi, double *xhilohi, double *xlolohi,
   double *xhihilo, double *xlohilo, double *xhilolo, double *xlololo,
   double *yhihihi, double *ylohihi, double *yhilohi, double *ylolohi,
   double *yhihilo, double *ylohilo, double *yhilolo, double *ylololo,
   double *zhihihi, double *zlohihi, double *zhilohi, double *zlolohi,
   double *zhihilo, double *zlohilo, double *zhilolo, double *zlololo,
   int dim );
/*
 * DESCRIPTION : z = x - y.
 *   Subtracts y to x to make z.  All arrays are of dimension dim. */

__global__ void dbl8_convolute
 ( double *xhihihi, double *xlohihi, double *xhilohi, double *xlolohi,
   double *xhihilo, double *xlohilo, double *xhilolo, double *xlololo,
   double *yhihihi, double *ylohihi, double *yhilohi, double *ylolohi,
   double *yhihilo, double *ylohilo, double *yhilolo, double *ylololo,
   double *zhihihi, double *zlohihi, double *zhilohi, double *zlolohi,
   double *zhihilo, double *zlohilo, double *zhilolo, double *zlololo,
   int dim );
/*
 * DESCRIPTION :
 *   Convolutes real vectors x and y to make z.
 *   All arrays have dimension dim.
 *
 * ON ENTRY :
 *   xhihihi    highest parts of the first vector x;
 *   xlohihi    second highest parts of the first vector x;
 *   xhilohi    third highest parts of the first vector x;
 *   xlolohi    fourth highest parts of the first vector x;
 *   xhihilo    fourth lowest parts of the first vector x;
 *   xlohilo    third lowest parts of the first vector x;
 *   xhilolo    second lowest parts of the first vector x;
 *   xlololo    lowest parts of the first vector x;
 *   yhihihi    highest parts of the second vector y;
 *   ylohihi    second highest parts of the second vector y;
 *   yhilohi    third highest parts of the second vector y;
 *   yholohi    fourth highest parts of the second vector y;
 *   yhihilo    fourth lowest parts of the second vector y;
 *   ylohilo    third lowest parts of the second vector y;
 *   yhilolo    second lowest parts of the second vector y;
 *   ylololo    lowest parts of the second vector y;
 *   zhihihi    dim doubles allocated for the highest parts of z;
 *   zlohihi    dim doubles allocated for the second highest parts of z;
 *   zhilohi    dim doubles allocated for the third highest parts of z;
 *   zlolohi    dim doubles allocated for the fourth highest parts of z;
 *   zhihilo    dim doubles allocated for the fourth lowest parts of z;
 *   zlohilo    dim doubles allocated for the third lowest parts of z;
 *   zhilolo    dim doubles allocated for the second lowest parts of z;
 *   zlololo    dim doubles allocated for the lowest parts of z;
 *   dim        dimension of all arrays.
 *
 * ON RETURN :
 *   zhihihi    highest parts of the product of x with y;
 *   zlohihi    second highest parts of the product of x with y;
 *   zhilohi    third highest parts of the product of x with y;
 *   zlolohi    fourth highest parts of the product of x with y;
 *   zhihilo    fourth lowest parts of the product of x with y;
 *   zlohilo    third lowest parts of the product of x with y;
 *   zhilolo    second lowest parts of the product of x with y;
 *   zlololo    lowest parts of the product of x with y. */

__global__ void dbl8_padded_convolute
 ( double *xhihihi, double *xlohihi, double *xhilohi, double *xlolohi,
   double *xhihilo, double *xlohilo, double *xhilolo, double *xlololo,
   double *yhihihi, double *ylohihi, double *yhilohi, double *ylolohi,
   double *yhihilo, double *ylohilo, double *yhilolo, double *ylololo,
   double *zhihihi, double *zlohihi, double *zhilohi, double *zlolohi,
   double *zhihilo, double *zlohilo, double *zhilolo, double *zlololo,
   int dim );
/*
 * DESCRIPTION :
 *   Convolutes real vectors x and y to make z.
 *   All arrays have dimension dim.
 *   Zeros are inserted in the shared memory for the second vector,
 *   so all threads perform the same number of operations
 *   and thread divergence is avoided. */

__global__ void cmplx8_convolute
 ( double *xrehihihi, double *xrelohihi, double *xrehilohi, double *xrelolohi,
   double *xrehihilo, double *xrelohilo, double *xrehilolo, double *xrelololo,
   double *ximhihihi, double *ximlohihi, double *ximhilohi, double *ximlolohi,
   double *ximhihilo, double *ximlohilo, double *ximhilolo, double *ximlololo,
   double *yrehihihi, double *yrelohihi, double *yrehilohi, double *yrelolohi,
   double *yrehihilo, double *yrelohilo, double *yrehilolo, double *yrelololo,
   double *yimhihihi, double *yimlohihi, double *yimhilohi, double *yimlolohi,
   double *yimhihilo, double *yimlohilo, double *yimhilolo, double *yimlololo,
   double *zrehihihi, double *zrelohihi, double *zrehilohi, double *zrelolohi,
   double *zrehihilo, double *zrelohilo, double *zrehilolo, double *zrelololo,
   double *zimhihihi, double *zimlohihi, double *zimhilohi, double *zimlolohi,
   double *zimhihilo, double *zimlohilo, double *zimhilolo, double *zimlololo,
   int dim );
/*
 * DESCRIPTION :
 *   Convolutes complex vectors x and y to make z.
 *   The complex vectors are given as double arrays with real
 *   and imaginary parts, as eight double arrays for all eight
 *   parts of the octo double numbers.
 *   All arrays have dimension dim.
 *
 * ON ENTRY :
 *   xrehihihi    highest parts of the real parts of x;
 *   xrelohihi    second highest parts of the real parts of x;
 *   xrehilohi    third highest parts of the real parts of x;
 *   xrelolohi    fourth highest parts of the real parts of x;
 *   xrehihilo    fourthllowest parts of the real parts of x;
 *   xrelohilo    third lowest parts of the real parts of x;
 *   xrehilolo    second lowest parts of the real parts of x;
 *   xrelololo    lowest parts of the real parts of x;
 *   ximhihihi    highest parts of the imaginary parts of x;
 *   ximlohihi    second highest parts of the imaginary parts of x;
 *   ximhilohi    third highest parts of the imaginary parts of x;
 *   ximlolohi    fourth highest parts of the imaginary parts of x;
 *   ximhihilo    fourth lowest parts of the imaginary parts of x;
 *   ximlohilo    third lowest parts of the imaginary parts of x;
 *   ximhilolo    second lowest parts of the imaginary parts of x;
 *   ximlololo    lowest parts of the imaginary parts of x;
 *   yrehihihi    highest parts of the real parts of y;
 *   yrelohihi    second highest parts of the real parts of y;
 *   yrehilohi    third highest parts of the real parts of y;
 *   yrelolohi    fourth highest parts of the real parts of y;
 *   yrehihilo    fourth lowest parts of the real parts of y;
 *   yrelohilo    third lowest parts of the real parts of y;
 *   yrehilolo    second lowest parts of the real parts of y;
 *   yrelololo    lowest parts of the real parts of y;
 *   yimhihihi    highest parts of the imaginary parts of y;
 *   yimlohihi    second highest parts of the imaginary parts of y;
 *   yimhilohi    third highest parts of the imaginary parts of y;
 *   yimlolohi    fourth highest parts of the imaginary parts of y;
 *   yimhihilo    fourth lowest parts of the imaginary parts of y;
 *   yimlohilo    third lowest parts of the imaginary parts of y;
 *   yimhilolo    second lowest parts of the imaginary parts of y;
 *   yimlololo    lowest parts of the imaginary parts of y;
 *   zrehihihi    dim doubles allocated for the highest parts
 *                of the real parts of the coefficients of the product;
 *   zrelohihi    dim doubles allocated for the second highest parts
 *                of the real parts of the coefficients of the product;
 *   zrehilohi    dim doubles allocated for the third highest parts
 *                of the real parts of the coefficients of the product;
 *   zrelolohi    dim doubles allocated for the fourth highest parts
 *                of the real parts of the coefficients of the product;
 *   zrehihilo    dim doubles allocated for the fourth lowest parts
 *                of the real parts of the coefficients of the product;
 *   zrelohilo    dim doubles allocated for the third lowest parts
 *                of the real parts of the coefficients of the product;
 *   zrehilolo    dim doubles allocated for the second lowest parts
 *                of the real parts of the coefficients of the product;
 *   zrelololo    dim doubles allocated for the lowest parts
 *                of the real parts of the coefficients of the product;
 *   zimhihihi    dim doubles allocated for the highest parts
 *                of the imaginary parts of the coefficients of the product;
 *   zimlohihi    dim doubles allocated for the second highest parts
 *                of the imaginary parts of the coefficients of the product;
 *   zimhilohi    dim doubles allocated for the third highest parts
 *                of the imaginary parts of the coefficients of the product;
 *   zimlolohi    dim doubles allocated for the fourth highest parts
 *                of the imaginary parts of the coefficients of the product;
 *   zimhihilo    dim doubles allocated for the fourth lowest parts
 *                of the imaginary parts of the coefficients of the product;
 *   zimlohilo    dim doubles allocated for the third lowest parts
 *                of the imaginary parts of the coefficients of the product;
 *   zimhilolo    dim doubles allocated for the second lowest parts
 *                of the imaginary parts of the coefficients of the product;
 *   zimlololo    dim doubles allocated for the lowest parts
 *                of the imaginary parts of the coefficients of the product;
 *   dim          dimension of all arrays.
 *
 * ON RETURN :
 *   zrehihihi    highest parts of the real parts of z;
 *   zrelohihi    second highest parts of the real parts of z;
 *   zrehilohi    third highest parts of the real parts of z;
 *   zrelolohi    fourth highest parts of the real parts of z;
 *   zrehihilo    fourth lowest parts of the real parts of z;
 *   zrelohilo    third lowest parts of the real parts of z;
 *   zrehilolo    second lowest parts of the real parts of z;
 *   zrelololo    lowest parts of the real parts of z;
 *   zimhihihi    highest parts of the imaginary parts of z;
 *   zimlohihi    second highest parts of the imaginary parts of z;
 *   zimhilohi    third highest parts of the imaginary parts of z;
 *   zimlolohi    fourth highest parts of the imaginary parts of z;
 *   zimhihilo    fourth lowest parts of the imaginary parts of z;
 *   zimlohilo    third lowest parts of the imaginary parts of z;
 *   zimhilolo    second lowest parts of the imaginary parts of z;
 *   zimlololo    lowest parts of the imaginary parts of z. */

__global__ void cmplx8_padded_convolute
 ( double *xrehihihi, double *xrelohihi, double *xrehilohi, double *xrelolohi,
   double *xrehihilo, double *xrelohilo, double *xrehilolo, double *xrelololo,
   double *ximhihihi, double *ximlohihi, double *ximhilohi, double *ximlolohi,
   double *ximhihilo, double *ximlohilo, double *ximhilolo, double *ximlololo,
   double *yrehihihi, double *yrelohihi, double *yrehilohi, double *yrelolohi,
   double *yrehihilo, double *yrelohilo, double *yrehilolo, double *yrelololo,
   double *yimhihihi, double *yimlohihi, double *yimhilohi, double *yimlolohi,
   double *yimhihilo, double *yimlohilo, double *yimhilolo, double *yimlololo,
   double *zrehihihi, double *zrelohihi, double *zrehilohi, double *zrelolohi,
   double *zrehihilo, double *zrelohilo, double *zrehilolo, double *zrelololo,
   double *zimhihihi, double *zimlohihi, double *zimhilohi, double *zimlolohi,
   double *zimhihilo, double *zimlohilo, double *zimhilolo, double *zimlololo,
   int dim );
/*
 * DESCRIPTION :
 *   Convolutes complex vectors x and y to make z.
 *   The complex vectors are given as double arrays with real
 *   and imaginary parts, as eight double arrays for all eight
 *   parts of the octo double numbers.  All arrays have dimension dim.
 *   Zeros are inserted in the shared memory for the second vector,
 *   so all threads perform the same number of operations
 *   and thread divergence is avoided. */

void GPU_dbl8_product
 ( double *xhihihi_h, double *xlohihi_h, double *xhilohi_h, double *xlolohi_h,
   double *xhihilo_h, double *xlohilo_h, double *xhilolo_h, double *xlololo_h,
   double *yhihihi_h, double *ylohihi_h, double *yhilohi_h, double *ylolohi_h,
   double *yhihilo_h, double *ylohilo_h, double *yhilolo_h, double *ylololo_h,
   double *zhihihi_h, double *zlohihi_h, double *zhilohi_h, double *zlolohi_h,
   double *zhihilo_h, double *zlohilo_h, double *zhilolo_h, double *zlololo_h,
   int deg, int freq, int BS, int padded );
/*
 * DESCRIPTION :
 *   Computes the product of two power series x and y to make z,
 *   with octo double coefficients.
 *
 * ON ENTRY :
 *   xhihihi_h    deg+1 highest parts of the coefficients of x;
 *   xlohihi_h    deg+1 second highest parts of the coefficients of x;
 *   xhilohi_h    deg+1 third highest parts of the coefficients of x;
 *   xlolohi_h    deg+1 fourth highest parts of the coefficients of x;
 *   xhihilo_h    deg+1 fourth lowest parts of the coefficients of x;
 *   xlohilo_h    deg+1 third lowest parts of the coefficients of x;
 *   xhilolo_h    deg+1 second lowest parts of the coefficients of x;
 *   xlololo_h    deg+1 lowest parts of the coefficients of x;
 *   yhihihi_h    deg+1 highest parts of the coefficients of y;
 *   ylohihi_h    deg+1 second highest parts of the coefficients of y;
 *   yhilohi_h    deg+1 third highest parts of the coefficients of y;
 *   ylolohi_h    deg+1 fourth highest parts of the coefficients of y;
 *   yhihilo_h    deg+1 fourth lowest parts of the coefficients of y;
 *   ylohilo_h    deg+1 third lowest parts of the coefficients of y;
 *   yhilolo_h    deg+1 second lowest parts of the coefficients of y;
 *   ylololo_h    deg+1 lowest parts of the coefficients of y;
 *   zhihihi_h    space for deg+1 doubles for the highest parts of z;
 *   zlohihi_h    space for deg+1 doubles for the second highest parts of z;
 *   zhilohi_h    space for deg+1 doubles for the third highest parts of z;
 *   zlolohi_h    space for deg+1 doubles for the fourth highest parts of z;
 *   zhihilo_h    space for deg+1 doubles for the fourth lowest parts of z;
 *   zlohilo_h    space for deg+1 doubles for the third lowest parts of z;
 *   zhilolo_h    space for deg+1 doubles for the second lowest parts of z;
 *   zlololo_h    space for deg+1 doubles for the lowest parts of z;
 *   deg          degree of the truncated power series;
 *   freq         frequency for timing purposes;
 *   BS           block size, the number of threads in a block;
 *   padded       if 1, then the padded convolution will run.
 *
 * ON RETURN :
 *   zhihihi_h    highest parts of the coefficients of the product z;
 *   zlohihi_h    second highest parts of the coefficients of the product z;
 *   zhilohi_h    third highest parts of the coefficients of the product z;
 *   zlolohi_h    fourth highest parts of the coefficients of the product z;
 *   zhihilo_h    fourth lowest parts of the coefficients of the product z;
 *   zlohilo_h    third lowest parts of the coefficients of the product z;
 *   zhilolo_h    second lowest parts of the coefficients of the product z;
 *   zlololo_h    lowest parts of the coefficients of the product z. */

void GPU_cmplx8_product
 ( double *xrehihihi_h, double *xrelohihi_h,
   double *xrehilohi_h, double *xrelolohi_h,
   double *xrehihilo_h, double *xrelohilo_h,
   double *xrehilolo_h, double *xrelololo_h,
   double *ximhihihi_h, double *ximlohihi_h,
   double *ximhilohi_h, double *ximlolohi_h,
   double *ximhihilo_h, double *ximlohilo_h,
   double *ximhilolo_h, double *ximlololo_h,
   double *yrehihihi_h, double *yrelohihi_h,
   double *yrehilohi_h, double *yrelolohi_h,
   double *yrehihilo_h, double *yrelohilo_h,
   double *yrehilolo_h, double *yrelololo_h,
   double *yimhihihi_h, double *yimlohihi_h,
   double *yimhilohi_h, double *yimlolohi_h,
   double *yimhihilo_h, double *yimlohilo_h,
   double *yimhilolo_h, double *yimlololo_h,
   double *zrehihihi_h, double *zrelohihi_h,
   double *zrehilohi_h, double *zrelolohi_h,
   double *zrehihilo_h, double *zrelohilo_h,
   double *zrehilolo_h, double *zrelololo_h,
   double *zimhihihi_h, double *zimlohihi_h,
   double *zimhilohi_h, double *zimlolohi_h,
   double *zimhihilo_h, double *zimlohilo_h,
   double *zimhilolo_h, double *zimlololo_h,
   int deg, int freq, int BS, int mode );
/*
 * DESCRIPTION :
 *   Computes the product of two power series x and y to make z,
 *   with complex octo double coefficients.
 *
 * ON ENTRY :
 *   xrehihihi_h    deg+1 highest parts of the real parts of x;
 *   xrelohihi_h    deg+1 second highest parts of the real parts of x;
 *   xrehilohi_h    deg+1 third highest parts of the real parts of x;
 *   xrelolohi_h    deg+1 fourth highest parts of the real parts of x;
 *   xrehihilo_h    deg+1 fourth lowest parts of the real parts of x;
 *   xrelohilo_h    deg+1 third lowest parts of the real parts of x;
 *   xrehilolo_h    deg+1 second lowest parts of the real parts of x;
 *   xrelololo_h    deg+1 lowest parts of the real parts of x;
 *   ximhihihi_h    deg+1 highest parts of the imag parts of x;
 *   ximlohihi_h    deg+1 second highest parts of the imag parts of x;
 *   ximhilohi_h    deg+1 third highest parts of the imag parts of x;
 *   ximlolohi_h    deg+1 fourth highest parts of the imag parts of x;
 *   ximhihilo_h    deg+1 fourth lowest parts of the imag parts of x;
 *   ximlohilo_h    deg+1 third lowest parts of the imag parts of x;
 *   ximhilolo_h    deg+1 second lowest parts of the imag parts of x;
 *   ximlololo_h    deg+1 lowest parts of the imag parts of x;
 *   yrehihihi_h    deg+1 highest parts of the real parts of y;
 *   yrelohihi_h    deg+1 second highest parts of the real parts of y;
 *   yrehilohi_h    deg+1 third highest parts of the real parts of y;
 *   yrelolohi_h    deg+1 fourth highest parts of the real parts of y;
 *   yrehihilo_h    deg+1 fourth lowest parts of the real parts of y;
 *   yrelohilo_h    deg+1 third lowest parts of the real parts of y;
 *   yrehilolo_h    deg+1 second lowest parts of the real parts of y;
 *   yrelololo_h    deg+1 lowest parts of the real parts of y;
 *   yimhihihi_h    deg+1 highest parts of the imag parts of y;
 *   yimlohihi_h    deg+1 second highest parts of the imag parts of y;
 *   yimhilohi_h    deg+1 third highest parts of the imag parts of y;
 *   yimlolohi_h    deg+1 fourth highest parts of the imag parts of y;
 *   yimhihilo_h    deg+1 fourth lowest parts of the imag parts of y;
 *   yimlohilo_h    deg+1 third lowest parts of the imag parts of y;
 *   yimhilolo_h    deg+1 second lowest parts of the imag parts of y;
 *   yimlololo_h    deg+1 lowest parts of the imag parts of y;
 *   deg            degree of the truncated power series;
 *   freq           frequency for timing purposes;
 *   BS             block size, the number of threads in a block;
 *   mode           0 : plain kernel,
 *                  1 : multiple real convolutions run,
 *                  2 : one padded complex convolution runs.
 *
 * ON RETURN :
 *   zrehihihi_h    highest parts of the real parts of z;
 *   zrelohihi_h    second highest parts of the real parts of z;
 *   zrehilohi_h    third highest parts of the real parts of z;
 *   zrelolohi_h    fourth highest parts of the real parts of z;
 *   zrehihilo_h    fourth lowest parts of the real parts of z;
 *   zrelohilo_h    third lowest parts of the real parts of z;
 *   zrehilolo_h    second lowest parts of the real parts of z;
 *   zrelololo_h    lowest parts of the real parts of z;
 *   zimhihihi_h    highest parts of the imag parts of z;
 *   zimlohihi_h    second highest parts of the imag parts of z;
 *   zimhilohi_h    third highest parts of the imag parts of z;
 *   zimlolohi_h    fourth highest parts of the imag parts of z;
 *   zimhihilo_h    fourth lowest parts of the imag parts of z;
 *   zimlohilo_h    third lowest parts of the imag parts of z;
 *   zimhilolo_h    second lowest parts of the imag parts of z;
 *   zimlololo_h    lowest parts of the imag parts of z. */

#endif
