// The file dbl10_convolutions_kernels.h defines prototypes for the kernels 
// to multiply two truncated power series in deca double precision.
// Defines the constant bound on the shared memory size of one block.

#ifndef __DBL10_CONVOLUTIONS_KERNELS_H__
#define __DBL10_CONVOLUTIONS_KERNELS_H__

#define da_shmemsize 64

/* The constant da_shmemsize is the bound on the shared memory size,
 * to compute the product of series with complex deca double coefficients.
 * For degree 63, we have 64 complex numbers, for x, y, and z,
 * real and imaginary parts, ten doubles, so 64*2*3*10*8 = 30720 bytes.
 * This constant bounds the degree of the power series. */

__global__ void dbl10_increment
 ( double *xrtb, double *xrix, double *xrmi, double *xrrg, double *xrpk,
   double *xltb, double *xlix, double *xlmi, double *xlrg, double *xlpk,
   double *yrtb, double *yrix, double *yrmi, double *yrrg, double *yrpk,
   double *yltb, double *ylix, double *ylmi, double *ylrg, double *ylpk,
   double *zrtb, double *zrix, double *zrmi, double *zrrg, double *zrpk,
   double *zltb, double *zlix, double *zlmi, double *zlrg, double *zlpk,
   int dim );
/*
 * DESCRIPTION : z = x + y.
 *   Adds y to x to make z.  All arrays are of dimension dim. */

__global__ void dbl10_decrement
 ( double *xrtb, double *xrix, double *xrmi, double *xrrg, double *xrpk,
   double *xltb, double *xlix, double *xlmi, double *xlrg, double *xlpk,
   double *yrtb, double *yrix, double *yrmi, double *yrrg, double *yrpk,
   double *yltb, double *ylix, double *ylmi, double *ylrg, double *ylpk,
   double *zrtb, double *zrix, double *zrmi, double *zrrg, double *zrpk,
   double *zltb, double *zlix, double *zlmi, double *zlrg, double *zlpk,
   int dim );
/*
 * DESCRIPTION : z = x + y.
 *   Subtracts y from x to make z.  All arrays are of dimension dim. */

__global__ void dbl10_convolute
 ( double *xrtb, double *xrix, double *xrmi, double *xrrg, double *xrpk,
   double *xltb, double *xlix, double *xlmi, double *xlrg, double *xlpk,
   double *yrtb, double *yrix, double *yrmi, double *yrrg, double *yrpk,
   double *yltb, double *ylix, double *ylmi, double *ylrg, double *ylpk,
   double *zrtb, double *zrix, double *zrmi, double *zrrg, double *zrpk,
   double *zltb, double *zlix, double *zlmi, double *zlrg, double *zlpk,
   int dim );
/*
 * DESCRIPTION :
 *   Convolutes real vectors x and y to make z.
 *   All arrays have dimension dim.
 *
 * ON ENTRY :
 *   xrtb     highest parts of the first vector x;
 *   xrix     second highest parts of the first vector x;
 *   xrmi     third highest parts of the first vector x;
 *   xrrg     fourth highest parts of the first vector x;
 *   xrpk     fifth highest parts of the first vector x;
 *   xltb     fifth lowest parts of the first vector x;
 *   xlix     fourth lowest parts of the first vector x;
 *   xlmi     third lowest parts of the first vector x;
 *   xlrg     second lowest parts of the first vector x;
 *   xlpk     lowest parts of the first vector x;
 *   yrtb     highest parts of the second vector y;
 *   yrix     second highest parts of the second vector y;
 *   yrmi     third highest parts of the second vector y;
 *   yrrg     fourth highest parts of the second vector y;
 *   yrpk     fifth highest parts of the second vector y;
 *   yltb     fifth lowest parts of the second vector y;
 *   ylix     fourth lowest parts of the second vector y;
 *   ylmi     third lowest parts of the second vector y;
 *   ylrg     second lowest parts of the second vector y;
 *   ylpk     lowest parts of the second vector y;
 *   zrtb     dim doubles allocated for the highest parts of z;
 *   zrix     dim doubles allocated for the second highest parts z;
 *   zrmi     dim doubles allocated for the third highest parts of z;
 *   zrrg     dim doubles allocated for the fourth highest parts of z;
 *   zrpk     dim doubles allocated for the fifth highest parts of z;
 *   zltb     dim doubles allocated for the fifth lowest parts of z;
 *   zlix     dim doubles allocated for the fourth lowest parts z;
 *   zlmi     dim doubles allocated for the third lowest parts of z;
 *   zlrg     dim doubles allocated for the second lowest parts of z;
 *   zlpk     dim doubles allocated for the lowest parts of z;
 *   dim      dimension of all arrays.
 *
 * ON RETURN :
 *   zrtb     highest parts of the product of x with y;
 *   zrix     second highest parts of the product of x with y;
 *   zrmi     third highest parts of the product of x with y;
 *   zrrg     fourth highest parts of the product of x with y;
 *   zrpk     fifth highest parts of the product of x with y;
 *   zltb     fifth lowest parts of the product of x with y;
 *   zlix     fourth lowest parts of the product of x with y;
 *   zlmi     third lowest parts of the product of x with y;
 *   zlrg     second lowest parts of the product of x with y;
 *   zlpk     lowest parts of the product of x with y. */

__global__ void dbl10_padded_convolute
 ( double *xrtb, double *xrix, double *xrmi, double *xrrg, double *xrpk,
   double *xltb, double *xlix, double *xlmi, double *xlrg, double *xlpk,
   double *yrtb, double *yrix, double *yrmi, double *yrrg, double *yrpk,
   double *yltb, double *ylix, double *ylmi, double *ylrg, double *ylpk,
   double *zrtb, double *zrix, double *zrmi, double *zrrg, double *zrpk,
   double *zltb, double *zlix, double *zlmi, double *zlrg, double *zlpk,
   int dim );
/*
 * DESCRIPTION :
 *   Convolutes real vectors x and y to make z.
 *   All arrays have dimension dim.
 *   Zeros are inserted in the shared memory for the second vector,
 *   so all threads perform the same number of operations
 *   and thread divergence is avoided. */

__global__ void cmplx10_convolute
 ( double *xrertb, double *xrerix, double *xrermi, double *xrerrg,
   double *xrerpk, double *xreltb, double *xrelix, double *xrelmi,
   double *xrelrg, double *xrelpk,
   double *ximrtb, double *ximrix, double *ximrmi, double *ximrrg,
   double *ximrpk, double *ximltb, double *ximlix, double *ximlmi,
   double *ximlrg, double *ximlpk,
   double *yrertb, double *yrerix, double *yrermi, double *yrerrg,
   double *yrerpk, double *yreltb, double *yrelix, double *yrelmi,
   double *yrelrg, double *yrelpk,
   double *yimrtb, double *yimrix, double *yimrmi, double *yimrrg,
   double *yimrpk, double *yimltb, double *yimlix, double *yimlmi,
   double *yimlrg, double *yimlpk,
   double *zrertb, double *zrerix, double *zrermi, double *zrerrg,
   double *zrerpk, double *zreltb, double *zrelix, double *zrelmi,
   double *zrelrg, double *zrelpk,
   double *zimrtb, double *zimrix, double *zimrmi, double *zimrrg,
   double *zimrpk, double *zimltb, double *zimlix, double *zimlmi,
   double *zimlrg, double *zimlpk, int dim );
/*
 * DESCRIPTION :
 *   Convolutes complex vectors x and y to make z.
 *   The complex vectors are given as double arrays with real
 *   and imaginary parts, of ten arrays for all parts of the deca doubles.
 *   All arrays have dimension dim.
 *
 * ON ENTRY :
 *   xrertb   highest parts of the real parts of x;
 *   xrerix   second highest parts of the real parts of x;
 *   xrermi   third highest parts of the real parts of x;
 *   ximrrg   fourth highest parts of the imaginary parts of x;
 *   ximrpk   fifth highest parts of the imaginary parts of x;
 *   xreltb   fifth lowest parts of the real parts of x;
 *   xrelix   fourth lowest parts of the real parts of x;
 *   xrelmi   third lowest parts of the real parts of x;
 *   ximlrg   second lowest parts of the imaginary parts of x;
 *   ximlpk   lowest parts of the imaginary parts of x;
 *   yrertb   highest parts of the real parts of y;
 *   yrerix   second highest parts of the real parts of y;
 *   yrermi   third highest parts of the real parts of y;
 *   yimrrg   fourth highest parts of the imaginary parts of y;
 *   yimrpk   fifth highest parts of the imaginary parts of y;
 *   yreltb   fifth lowest parts of the real parts of y;
 *   yrelix   fourth lowest parts of the real parts of y;
 *   yrelmi   third lowest parts of the real parts of y;
 *   yimlrg   second lowest parts of the imaginary parts of y;
 *   yimlpk   lowest parts of the imaginary parts of y;
 *   zrertb   dim doubles allocated for the highest parts
 *            of the real parts of the coefficients of the product;
 *   zrerix   dim doubles allocated for the second highest parts
 *            of the real parts of the coefficients of the product;
 *   zrermi   dim doubles allocated for the third highest parts
 *            of the real parts of the coefficients of the product;
 *   zrerrg   dim doubles allocated for the fourth highest parts
 *            of the real parts of the coefficients of the product;
 *   zrerpk   dim doubles allocated for the fifth highest parts
 *            of the real parts of the coefficients of the product;
 *   zreltb   dim doubles allocated for the fifth lowest parts
 *            of the real parts of the coefficients of the product;
 *   zrelix   dim doubles allocated for the fourth lowest parts
 *            of the real parts of the coefficients of the product;
 *   zrelmi   dim doubles allocated for the third lowest parts
 *            of the real parts of the coefficients of the product;
 *   zrelrg   dim doubles allocated for the second lowest parts
 *            of the real parts of the coefficients of the product;
 *   zrelpk   dim doubles allocated for the lowest parts
 *            of the real parts of the coefficients of the product;
 *   zimrtb   dim doubles allocated for the highest parts
 *            of the imaginary parts of the coefficients of the product;
 *   zimrix   dim doubles allocated for the second highest parts
 *            of the imaginary parts of the coefficients of the product;
 *   zimrmi   dim doubles allocated for the third highest parts
 *            of the imaginary parts of the coefficients of the product;
 *   zimrrg   dim doubles allocated for the fourth highest parts
 *            of the imaginary parts of the coefficients of the product;
 *   zimrpk   dim doubles allocated for the fifth highest parts
 *            of the imaginary parts of the coefficients of the product;
 *   zimltb   dim doubles allocated for the fifth lowest parts
 *            of the imaginary parts of the coefficients of the product;
 *   zimlix   dim doubles allocated for the fourth highest parts
 *            of the imaginary parts of the coefficients of the product;
 *   zimlmi   dim doubles allocated for the third lowest parts
 *            of the imaginary parts of the coefficients of the product;
 *   zimlrg   dim doubles allocated for the second lowest parts
 *            of the imaginary parts of the coefficients of the product;
 *   zimlpk   dim doubles allocated for the lowest parts
 *            of the imaginary parts of the coefficients of the product;
 *   dim      dimension of all arrays.
 *
 * ON RETURN :
 *   zrertb   highest parts of the real parts of the product z;
 *   zrerix   second highest parts of the real parts of the product z;
 *   zrermi   third highest parts of the real parts of the product z;
 *   zrerrg   fourth highest parts of the real parts of the product z;
 *   zrerpk   fifth highest parts of the real parts of the product z;
 *   zreltb   fifth lowest parts of the real parts of the product z;
 *   zrelix   fourth lowest parts of the real parts of the product z;
 *   zrelmi   third lowest parts of the real parts of the product z;
 *   zrelrg   second lowest parts of the real parts of the product z;
 *   zrerpk   lowest parts of the real parts of the product z;
 *   zimrtb   highest parts of the imaginary parts of the product z;
 *   zimrix   second highest parts of the imaginary parts of the product z;
 *   zimrmi   third highest parts of the imaginary parts of the product z;
 *   zimrrg   fourth highest parts of the imaginary parts of the product z;
 *   zimrpk   fifth highest parts of the imaginary parts of the product z;
 *   zimltb   fifth lowest parts of the imaginary parts of the product z;
 *   zimlix   fourth lowest parts of the imaginary parts of the product z;
 *   zimlmi   third lowest parts of the imaginary parts of the product z;
 *   zimlrg   second lowest parts of the imaginary parts of the product z;
 *   zimlpk   lowest parts of the imaginary parts of the product z. */

__global__ void cmplx10_padded_convolute
 ( double *xrertb, double *xrerix, double *xrermi, double *xrerrg,
   double *xrerpk, double *xreltb, double *xrelix, double *xrelmi,
   double *xrelrg, double *xrelpk,
   double *ximrtb, double *ximrix, double *ximrmi, double *ximrrg,
   double *ximrpk, double *ximltb, double *ximlix, double *ximlmi,
   double *ximlrg, double *ximlpk,
   double *yrertb, double *yrerix, double *yrermi, double *yrerrg,
   double *yrerpk, double *yreltb, double *yrelix, double *yrelmi,
   double *yrelrg, double *yrelpk,
   double *yimrtb, double *yimrix, double *yimrmi, double *yimrrg,
   double *yimrpk, double *yimltb, double *yimlix, double *yimlmi,
   double *yimlrg, double *yimlpk,
   double *zrertb, double *zrerix, double *zrermi, double *zrerrg,
   double *zrerpk, double *zreltb, double *zrelix, double *zrelmi,
   double *zrelrg, double *zrelpk,
   double *zimrtb, double *zimrix, double *zimrmi, double *zimrrg,
   double *zimrpk, double *zimltb, double *zimlix, double *zimlmi,
   double *zimlrg, double *zimlpk, int dim );
/*
 * DESCRIPTION :
 *   Convolutes complex vectors x and y to make z.
 *   The complex vectors are given as double arrays with real
 *   and imaginary parts, of ten arrays for all parts of the deca doubles.
 *   All arrays have dimension dim.
 *   Zeros are inserted in the shared memory for the second vector,
 *   so all threads perform the same number of operations
 *   and thread divergence is avoided. */

void GPU_dbl10_product
 ( double *xrtb_h, double *xrix_h, double *xrmi_h, double *xrrg_h,
   double *xrpk_h, double *xltb_h, double *xlix_h, double *xlmi_h,
   double *xlrg_h, double *xlpk_h,
   double *yrtb_h, double *yrix_h, double *yrmi_h, double *yrrg_h,
   double *yrpk_h, double *yltb_h, double *ylix_h, double *ylmi_h,
   double *ylrg_h, double *ylpk_h,
   double *zrtb_h, double *zrix_h, double *zrmi_h, double *zrrg_h,
   double *zrpk_h, double *zltb_h, double *zlix_h, double *zlmi_h,
   double *zlrg_h, double *zlpk_h, int deg, int freq, int BS, int padded );
/*
 * DESCRIPTION :
 *   Computes the product of two power series x and y to make z,
 *   with deca double coefficients.
 *
 * ON ENTRY :
 *   xrtb_h    deg+1 highest parts of the coefficients of x;
 *   xrix_h    deg+1 second highest parts of the coefficients of x;
 *   xrmi_h    deg+1 third highest parts of the coefficients of x;
 *   xrrg_h    deg+1 fourth highest parts of the coefficients of x;
 *   xrpk_h    deg+1 fifth highest parts of the coefficients of x;
 *   xltb_h    deg+1 fifth lowest parts of the coefficients of x;
 *   xlix_h    deg+1 fourth lowest parts of the coefficients of x;
 *   xlmi_h    deg+1 third lowest parts of the coefficients of x;
 *   xlrg_h    deg+1 second lowest parts of the coefficients of x;
 *   xlpk_h    deg+1 lowest parts of the coefficients of x;
 *   yltb_h    deg+1 highest parts of the coefficients of y;
 *   yrix_h    deg+1 second highest parts of the coefficients of y;
 *   yrmi_h    deg+1 third highest parts of the coefficients of y;
 *   yrrg_h    deg+1 fourth highest lowest parts of the coefficients of y;
 *   yrpk_h    deg+1 fifth highest parts of the coefficients of y;
 *   yltb_h    deg+1 fifth lowest parts of the coefficients of y;
 *   ylix_h    deg+1 fourth lowest parts of the coefficients of y;
 *   ylmi_h    deg+1 third lowest parts of the coefficients of y;
 *   ylrg_h    deg+1 second lowest parts of the coefficients of y;
 *   ylpk_h    deg+1 lowest parts of the coefficients of y;
 *   zrtb_h    space allocated for deg+1 doubles for the highest parts of z;
 *   zrix_h    space allocated for deg+1 doubles 
 *             for the second highest parts of z;
 *   zrmi_h    space allocated for deg+1 doubles
 *             for the third highest parts of z;
 *   zrrg_h    space allocated for deg+1 doubles 
 *             for the fourth highest parts of z;
 *   zrpk_h    space allocated for deg+1 doubles
 *             for the fifth highest parts of z;
 *   zltb_h    space allocated for deg+1 doubles
 *             for the fifth lowest parts of z;
 *   zlix_h    space allocated for deg+1 doubles 
 *             for the fourth lowest parts of z;
 *   zlmi_h    space allocated for deg+1 doubles
 *             for the third lowest parts of z;
 *   zlrg_h    space allocated for deg+1 doubles 
 *             for the second lowest parts of z;
 *   zlpk_h    space allocated for deg+1 doubles for the lowest parts of z;
 *   deg       degree of the truncated power series;
 *   freq      frequency for timing purposes;
 *   BS        block size, the number of threads in a block;
 *   padded    if 1, then the padded convolute is called.
 *
 * ON RETURN :
 *   zrtb_h    highest parts of the coefficients of the product z;
 *   zrix_h    second highest parts of the coefficients of the product z;
 *   zrmi_h    third highest parts of the coefficients of the product of z;
 *   zrrg_h    fourth highest parts of the coefficients of the product of z;
 *   zrpk_h    fifth highest parts of the coefficients of the product of z;
 *   zltb_h    fifth lowest parts of the coefficients of the product z;
 *   zlix_h    fourth lowest parts of the coefficients of the product z;
 *   zlmi_h    third lowest parts of the coefficients of the product of z;
 *   zlrg_h    second lowest parts of the coefficients of the product of z;
 *   zlpk_h    lowest parts of the coefficients of the product of z. */

void GPU_cmplx10_product
 ( double *xrertb_h, double *xrerix_h, double *xrermi_h, double *xrerrg_h,
   double *xrerpk_h, double *xreltb_h, double *xrelix_h, double *xrelmi_h,
   double *xrelrg_h, double *xrelpk_h,
   double *ximrtb_h, double *ximrix_h, double *ximrmi_h, double *ximrrg_h,
   double *ximrpk_h, double *ximltb_h, double *ximlix_h, double *ximlmi_h,
   double *ximlrg_h, double *ximlpk_h,
   double *yrertb_h, double *yrerix_h, double *yrermi_h, double *yrerrg_h,
   double *yrerpk_h, double *yreltb_h, double *yrelix_h, double *yrelmi_h,
   double *yrelrg_h, double *yrelpk_h,
   double *yimrtb_h, double *yimrix_h, double *yimrmi_h, double *yimrrg_h,
   double *yimrpk_h, double *yimltb_h, double *yimlix_h, double *yimlmi_h,
   double *yimlrg_h, double *yimlpk_h,
   double *zrertb_h, double *zrerix_h, double *zrermi_h, double *zrerrg_h,
   double *zrerpk_h, double *zreltb_h, double *zrelix_h, double *zrelmi_h,
   double *zrelrg_h, double *zrelpk_h,
   double *zimrtb_h, double *zimrix_h, double *zimrmi_h, double *zimrrg_h,
   double *zimrpk_h, double *zimltb_h, double *zimlix_h, double *zimlmi_h,
   double *zimlrg_h, double *zimlpk_h, int deg, int freq, int BS, int mode );
/*
 * DESCRIPTION :
 *   Computes the product of two power series x and y to make z,
 *   with complex deca double coefficients.
 *
 * ON ENTRY :
 *   xrertb_h   deg+1 highest parts of the real parts of x;
 *   xrerix_h   deg+1 second highest parts of the real parts of x;
 *   xrermi_h   deg+1 third highest parts of the real parts of x;
 *   xrerrg_h   deg+1 fourth highest parts of the real parts of x;
 *   xrerpk_h   deg+1 fifth highest parts of the real parts of x;
 *   xreltb_h   deg+1 fifth lowest parts of the real parts of x;
 *   xrelix_h   deg+1 fourth lowest parts of the real parts of x;
 *   xrelmi_h   deg+1 third lowest parts of the real parts of x;
 *   xrelrg_h   deg+1 second lowest parts of the real parts of x;
 *   xrelpk_h   deg+1 lowest parts of the real parts of x;
 *   ximrtb_h   deg+1 highest parts of the imag parts of x;
 *   ximrix_h   deg+1 second highest parts of the imag parts of x;
 *   ximrmi_h   deg+1 third highest parts of the imag parts of x;
 *   ximrix_h   deg+1 fourth highest parts of the imag parts of x;
 *   ximrpk_h   deg+1 fifth highest parts of the imag parts of x;
 *   ximltb_h   deg+1 fifth lowest parts of the imag parts of x;
 *   ximlix_h   deg+1 fourth lowest parts of the imag parts of x;
 *   ximlmi_h   deg+1 third lowest parts of the imag parts of x;
 *   ximlix_h   deg+1 second lowest parts of the imag parts of x;
 *   ximlpk_h   deg+1 lowest parts of the imag parts of x;
 *   yrertb_h   deg+1 highest parts of the real parts of y;
 *   yrerix_h   deg+1 second highest parts of the real parts of y;
 *   yrermi_h   deg+1 third highest parts of the real parts of y;
 *   yrerrg_h   deg+1 fourth highest parts of the real parts of y;
 *   yrerpk_h   deg+1 fifth highest parts of the real parts of y;
 *   yreltb_h   deg+1 fifth lowest parts of the real parts of y;
 *   yrelix_h   deg+1 fourth lowest parts of the real parts of y;
 *   yrelmi_h   deg+1 third lowest parts of the real parts of y;
 *   yrelrg_h   deg+1 second lowest parts of the real parts of y;
 *   yrelpk_h   deg+1 lowest parts of the real parts of y;
 *   yimrtb_h   deg+1 highest parts of the imag parts of y;
 *   yimrix_h   deg+1 second highest parts of the imag parts of y;
 *   yimrmi_h   deg+1 third highest parts of the imag parts of y;
 *   yimrrg_h   deg+1 fourth highest parts of the imag parts of y;
 *   yimrpk_h   deg+1 fifth highest parts of the imag parts of y;
 *   yimltb_h   deg+1 fifth lowest parts of the imag parts of y;
 *   yimlix_h   deg+1 fourth lowest parts of the imag parts of y;
 *   yimlmi_h   deg+1 third lowest parts of the imag parts of y;
 *   yimlrg_h   deg+1 second lowest parts of the imag parts of y;
 *   yimlpk_h   deg+1 lowest parts of the imag parts of y;
 *   zrertb_h   space for deg+1 doubles for the highest real parts of z;
 *   zrerix_h   space for deg+1 doubles 
 *              for the second highest real parts of z;
 *   zrermi_h   space for deg+1 doubles for the third highest real parts of z;
 *   zrerrg_h   space for deg+1 doubles
 *              for the fourth highest real parts of z;
 *   zrerpk_h   space for deg+1 doubles for the fifth highest real parts of z;
 *   zreltb_h   space for deg+1 doubles for the fifth lowest real parts of z;
 *   zrelix_h   space for deg+1 doubles for the fourth lowest real parts of z;
 *   zrelmi_h   space for deg+1 doubles for the third lowest real parts of z;
 *   zrelrg_h   space for deg+1 doubles for the second lowest real parts of z;
 *   zrelpk_h   space for deg+1 doubles for the lowest real parts of z;
 *   zimrtb_h   space for deg+1 doubles for the highest imag parts of z;
 *   zimrix_h   space for deg+1 doubles 
 *              for the second highest imaginary parts of z;
 *   zimrmi_h   space for deg+1 doubles for the third highest imag parts of z;
 *   zimrrg_h   space for deg+1 doubles
 *              for the fourth highest imaginary parts of z;
 *   zimrpk_h   space for deg+1 doubles for the fifth highest imag parts of z;
 *   zimltb_h   space for deg+1 doubles for the fifth lowest imag parts of z;
 *   zimlix_h   space for deg+1 doubles 
 *              for the fourth lowest imaginary parts of z;
 *   zimlmi_h   space for deg+1 doubles for the third lowest imag parts of z;
 *   zimlrg_h   space for deg+1 doubles
 *              for the second lowest imaginary parts of z;
 *   zimlpk_h   space for deg+1 doubles for the lowest imaginary parts of z;
 *   deg        degree of the truncated power series;
 *   freq       frequency for timing purposes;
 *   BS         block size, the number of threads in a block;
 *   mode       if 1, then multiple kernels will be launched;
 *              if 2, then one padded complex convolute runs.
 *
 * ON RETURN :
 *   zrertb_h   highest parts of the real parts of z;
 *   zrerix_h   second highest parts of the real parts of z;
 *   zrermi_h   third highest parts of the real parts of z;
 *   zrerrg_h   fourth highest parts of the real parts of z;
 *   zrerpk_h   fifth highest parts of the real parts of z;
 *   zreltb_h   fifth lowest parts of the real parts of z;
 *   zrelix_h   fourth lowest parts of the real parts of z;
 *   zrelmi_h   third lowest parts of the real parts of z;
 *   zrelrg_h   second lowest parts of the real parts of z;
 *   zrelpk_h   lowest parts of the real parts of z;
 *   zimrtb_h   highest parts of the imag parts of z;
 *   zimrix_h   second highest parts of the imag parts of z;
 *   zimrmi_h   third highest parts of the imag parts of z;
 *   zimrrg_h   fourth highest parts of the imag parts of z;
 *   zimrpk_h   fifth highest parts of the imag parts of z;
 *   zimltb_h   fifth lowest parts of the imag parts of z;
 *   zimlix_h   fourth lowest parts of the imag parts of z;
 *   zimlmi_h   third lowest parts of the imag parts of z;
 *   zimlrg_h   second lowest parts of the imag parts of z;
 *   zimlpk_h   lowest parts of the imag parts of z. */

#endif
