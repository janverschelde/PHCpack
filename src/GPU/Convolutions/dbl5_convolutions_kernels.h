// The file dbl5_convolutions_kernels.h defines prototypes for the kernels 
// to multiply two truncated power series in penta double precision.
// Defines the constant bound on the shared memory size of one block.

#ifndef __DBL5_CONVOLUTIONS_KERNELS_H__
#define __DBL5_CONVOLUTIONS_KERNELS_H__

#define pd_shmemsize 128

/* The constant pd_shmemsize is the bound on the shared memory size,
 * to compute the product of series with complex penta double coefficients.
 * For degree 127, we have 128 complex numbers, for x, y, and z,
 * real and imaginary parts, five doubles, so 192*2*3*5*8 = 30720 bytes.
 * This constant bounds the degree of the power series. */

__global__ void dbl5_increment
 ( double *xtb, double *xix, double *xmi, double *xrg, double *xpk,
   double *ytb, double *yix, double *ymi, double *yrg, double *ypk,
   double *ztb, double *zix, double *zmi, double *zrg, double *zpk, int dim );
/*
 * DESCRIPTION : z = x + y.
 *   Adds y to x to make z.  All arrays are of dimension dim. */

__global__ void dbl5_decrement
 ( double *xtb, double *xix, double *xmi, double *xrg, double *xpk,
   double *ytb, double *yix, double *ymi, double *yrg, double *ypk,
   double *ztb, double *zix, double *zmi, double *zrg, double *zpk, int dim );
/*
 * DESCRIPTION : z = x + y.
 *   Subtracts y from x to make z.  All arrays are of dimension dim. */

__global__ void dbl5_convolute
 ( double *xtb, double *xix, double *xmi, double *xrg, double *xpk,
   double *ytb, double *yix, double *ymi, double *yrg, double *ypk,
   double *ztb, double *zix, double *zmi, double *zrg, double *zpk, int dim );
/*
 * DESCRIPTION :
 *   Convolutes real vectors x and y to make z.
 *   All arrays have dimension dim.
 *
 * ON ENTRY :
 *   xtb      highest parts of the first vector x;
 *   xix      second highest parts of the first vector x;
 *   xmi      middle parts of the first vector x;
 *   xrg      second lowest parts of the first vector x;
 *   xpk      lowest parts of the first vector x;
 *   ytb      highest parts of the second vector y;
 *   yix      second highest parts of the second vector y;
 *   ymi      middle parts of the second vector y;
 *   yrg      second lowest parts of the second vector y;
 *   ypk      lowest parts of the second vector y;
 *   ztb      dim doubles allocated for the highest parts of the product;
 *   zix      dim doubles allocated for the 2nd highest parts of the product;
 *   zmi      dim doubles allocated for the middle parts of the product;
 *   zrg      dim doubles allocated for the 2nd lowest parts of the product;
 *   zpk      dim doubles allocated for the lowest parts of the product;
 *   dim      dimension of all arrays.
 *
 * ON RETURN :
 *   ztb      highest parts of the product of x with y;
 *   zix      second highest parts of the product of x with y;
 *   zmi      middle parts of the product of x with y;
 *   zrg      second lowest parts of the product of x with y;
 *   zpk      lowest parts of the product of x with y. */

__global__ void dbl5_padded_convolute
 ( double *xtb, double *xix, double *xmi, double *xrg, double *xpk,
   double *ytb, double *yix, double *ymi, double *yrg, double *ypk,
   double *ztb, double *zix, double *zmi, double *zrg, double *zpk, int dim );
/*
 * DESCRIPTION :
 *   Convolutes real vectors x and y to make z.
 *   All arrays have dimension dim.
 *   Zeros are inserted in the shared memory for the vector y,
 *   so all threads perform the same number of operations. */

__global__ void cmplx5_convolute
 ( double *xretb, double *xreix, double *xremi, double *xrerg, double *xrepk,
   double *ximtb, double *ximix, double *ximmi, double *ximrg, double *ximpk,
   double *yretb, double *yreix, double *yremi, double *yrerg, double *yrepk,
   double *yimtb, double *yimix, double *yimmi, double *yimrg, double *yimpk,
   double *zretb, double *zreix, double *zremi, double *zrerg, double *zrepk,
   double *zimtb, double *zimix, double *zimmi, double *zimrg, double *zimpk,
   int dim );
/*
 * DESCRIPTION :
 *   Convolutes complex vectors x and y to make z.
 *   The complex vectors are given as double arrays with real
 *   and imaginary parts, of highest, second highest, middle,
 *   second lowest and lowest parts.  All arrays have dimension dim.
 *
 * ON ENTRY :
 *   xretb    highest parts of the real parts of x;
 *   xreix    second highest parts of the real parts of x;
 *   xremi    middle parts of the real parts of x;
 *   ximrg    second lowest parts of the imaginary parts of x;
 *   ximpk    lowest parts of the imaginary parts of x;
 *   yretb    highest parts of the real parts of y;
 *   yreix    second highest parts of the real parts of y;
 *   yremi    middle parts of the real parts of y;
 *   yimrg    second lowest parts of the imaginary parts of y;
 *   yimpk    lowest parts of the imaginary parts of y;
 *   zretb    dim doubles allocated for the highest parts
 *            of the real parts of the coefficients of the product;
 *   zreix    dim doubles allocated for the second highest parts
 *            of the real parts of the coefficients of the product;
 *   zremi    dim doubles allocated for the middle parts
 *            of the real parts of the coefficients of the product;
 *   zrerg    dim doubles allocated for the second lowest parts
 *            of the real parts of the coefficients of the product;
 *   zrepk    dim doubles allocated for the lowest parts
 *            of the real parts of the coefficients of the product;
 *   zimtb    dim doubles allocated for the highest parts
 *            of the imaginary parts of the coefficients of the product;
 *   zimix    dim doubles allocated for the second highest parts
 *            of the imaginary parts of the coefficients of the product;
 *   zimmi    dim doubles allocated for the middle parts
 *            of the imaginary parts of the coefficients of the product;
 *   zimrg    dim doubles allocated for the second lowest parts
 *            of the imaginary parts of the coefficients of the product;
 *   zimpk    dim doubles allocated for the lowest parts
 *            of the imaginary parts of the coefficients of the product;
 *   dim      dimension of all arrays.
 *
 * ON RETURN :
 *   zretb    highest parts of the real parts of the product z;
 *   zreix    second highest parts of the real parts of the product z;
 *   zremi    middle parts of the real parts of the product z;
 *   zrerg    second lowest parts of the real parts of the product z;
 *   zrepk    lowest parts of the real parts of the product z;
 *   zimtb    highest parts of the imaginary parts of the product z;
 *   zimix    second highest parts of the imaginary parts of the product z;
 *   zimmi    middle parts of the imaginary parts of the product z;
 *   zimrg    second lowest parts of the imaginary parts of the product z;
 *   zimpk    lowest parts of the imaginary parts of the product z. */

__global__ void cmplx5_padded_convolute
 ( double *xretb, double *xreix, double *xremi, double *xrerg, double *xrepk,
   double *ximtb, double *ximix, double *ximmi, double *ximrg, double *ximpk,
   double *yretb, double *yreix, double *yremi, double *yrerg, double *yrepk,
   double *yimtb, double *yimix, double *yimmi, double *yimrg, double *yimpk,
   double *zretb, double *zreix, double *zremi, double *zrerg, double *zrepk,
   double *zimtb, double *zimix, double *zimmi, double *zimrg, double *zimpk,
   int dim );
/*
 * DESCRIPTION :
 *   Convolutes complex vectors x and y to make z.
 *   The complex vectors are given as double arrays with real
 *   and imaginary parts, of highest, second highest, middle,
 *   second lowest and lowest parts.  All arrays have dimension dim.
 *   Zeros are inserted in the shared memory for the second vector,
 *   so all threads perform the same number of operations
 *   and thread divergence is avoided. */

void GPU_dbl5_product
 ( double *xtb_h, double *xix_h, double *xmi_h, double *xrg_h, double *xpk_h,
   double *ytb_h, double *yix_h, double *ymi_h, double *yrg_h, double *ypk_h,
   double *ztb_h, double *zix_h, double *zmi_h, double *zrg_h, double *zpk_h,
   int deg, int freq, int BS, int padded );
/*
 * DESCRIPTION :
 *   Computes the product of two power series x and y to make z,
 *   with penta double coefficients.
 *
 * ON ENTRY :
 *   xtb_h    deg+1 highest parts of the coefficients of x;
 *   xix_h    deg+1 second highest parts of the coefficients of x;
 *   xmi_h    deg+1 middle parts of the coefficients of x;
 *   xrg_h    deg+1 second lowest parts of the coefficients of x;
 *   xpk_h    deg+1 lowest parts of the coefficients of x;
 *   ytb_h    deg+1 highest parts of the coefficients of y;
 *   yix_h    deg+1 second highest parts of the coefficients of y;
 *   ymi_h    deg+1 middle parts of the coefficients of y;
 *   yrg_h    deg+1 second lowest parts of the coefficients of y;
 *   ypk_h    deg+1 lowest parts of the coefficients of y;
 *   ztb_h    space allocated for deg+1 doubles for the highest parts of z;
 *   zix_h    space allocated for deg+1 doubles 
 *            for the second highest parts of z;
 *   zmi_h    space allocated for deg+1 doubles for the middle parts of z;
 *   zrg_h    space allocated for deg+1 doubles 
 *            for the second lowest parts of z;
 *   zpk_h    space allocated for deg+1 doubles for the lowest parts of z;
 *   deg      degree of the truncated power series;
 *   freq     frequency for timing purposes;
 *   BS       block size, the number of threads in a block;
 *   padded   if 1, then the zeros are inserted for the vector y.
 *
 * ON RETURN :
 *   ztb_h    highest parts of the coefficients of the product z;
 *   zix_h    second highest parts of the coefficients of the product z;
 *   zmi_h    middle parts of the coefficients of the product of z;
 *   zrg_h    second lowest parts of the coefficients of the product of z;
 *   zpk_h    lowest parts of the coefficients of the product of z. */

void GPU_cmplx5_product
 ( double *xretb_h, double *xreix_h, double *xremi_h, double *xrerg_h,
   double *xrepk_h, double *ximtb_h, double *ximix_h, double *ximmi_h,
   double *ximrg_h, double *ximpk_h, double *yretb_h, double *yreix_h,
   double *yremi_h, double *yrerg_h, double *yrepk_h, double *yimtb_h,
   double *yimix_h, double *yimmi_h, double *yimrg_h, double *yimpk_h,
   double *zretb_h, double *zreix_h, double *zremi_h, double *zrerg_h,
   double *zrepk_h, double *zimtb_h, double *zimix_h, double *zimmi_h,
   double *zimrg_h, double *zimpk_h, int deg, int freq, int BS,
   int mode );
/*
 * DESCRIPTION :
 *   Computes the product of two power series x and y to make z,
 *   with complex penta double coefficients.
 *
 * ON ENTRY :
 *   xretb_h    deg+1 highest parts of the real parts of x;
 *   xreix_h    deg+1 second highest parts of the real parts of x;
 *   xremi_h    deg+1 middle parts of the real parts of x;
 *   xrerg_h    deg+1 second lowest parts of the real parts of x;
 *   xrepk_h    deg+1 lowest parts of the real parts of x;
 *   ximtb_h    deg+1 highest parts of the imag parts of x;
 *   ximix_h    deg+1 second highest parts of the imag parts of x;
 *   ximmi_h    deg+1 middle parts of the imag parts of x;
 *   ximix_h    deg+1 second lowest parts of the imag parts of x;
 *   ximpk_h    deg+1 lowest parts of the imag parts of x;
 *   yretb_h    deg+1 highest parts of the real parts of y;
 *   yreix_h    deg+1 second highest parts of the real parts of y;
 *   yremi_h    deg+1 middle parts of the real parts of y;
 *   yrerg_h    deg+1 second lowest parts of the real parts of y;
 *   yrepk_h    deg+1 lowest parts of the real parts of y;
 *   yimtb_h    deg+1 highest parts of the imag parts of y;
 *   yimix_h    deg+1 second highest parts of the imag parts of y;
 *   yimmi_h    deg+1 middle parts of the imag parts of y;
 *   yimrg_h    deg+1 second lowest parts of the imag parts of y;
 *   yimpk_h    deg+1 lowest parts of the imag parts of y;
 *   zretb_h    space for deg+1 doubles for the highest real parts of z;
 *   zreix_h    space for deg+1 doubles 
 *              for the second highest real parts of z;
 *   zremi_h    space for deg+1 doubles for the middle real parts of z;
 *   zrerg_h    space for deg+1 doubles for the second lowest real parts of z;
 *   zrepk_h    space for deg+1 doubles for the lowest real parts of z;
 *   zimtb_h    space for deg+1 doubles for the highest imag parts of z;
 *   zimix_h    space for deg+1 doubles 
 *              for the second highest imaginary parts of z;
 *   zimmi_h    space for deg+1 doubles for the middle imag parts of z;
 *   zimrg_h    space for deg+1 doubles
 *              for the second lowest imaginary parts of z;
 *   zimpk_h    space for deg+1 doubles for the lowest imag parts of z;
 *   deg        degree of the truncated power series;
 *   freq       frequency for timing purposes;
 *   BS         block size, the number of threads in a block;
 *   looped     if 1, then multiple double kernels are launched;
 *              if 2, then one padded complex convolute is applied.
 *
 * ON RETURN :
 *   zretb_h    highest parts of the real parts of z;
 *   zreix_h    second highest parts of the real parts of z;
 *   zremi_h    middle parts of the real parts of z;
 *   zrerg_h    second lowest parts of the real parts of z;
 *   zrepk_h    lowest parts of the real parts of z;
 *   zimtb_h    highest parts of the imag parts of z;
 *   zimix_h    second highest parts of the imag parts of z;
 *   zimmi_h    middle parts of the imag parts of z;
 *   zimrg_h    second lowest parts of the imag parts of z;
 *   zimpk_h    lowest parts of the imag parts of z. */

#endif
