/* The file dbl_baqr_flopcounts.h specifies functions to count
 * the floating-point operations of the blocked accelerated QR,
 * in double precision. */

/* The names of the functions correspond to the names of the kernels
 * for which the operations are counted.
 * As the long long int type is limited to 64-bits,
 * the accuracy of the counts can only be valid for computations
 * that take no more than 10 seconds,
 * in case teraflop performance is attained. */

#ifndef __dbl_baqr_flopcounts_h__
#define __dbl_baqr_flopcounts_h__

void flopcount_dbl_small_house
 ( int dim, int dimLog2, long long int *add, long long int *mul,
   long long int *div, long long int *sqrtfun );
/*
 * DESCRIPTION :
 *   Returns the accumulated floating-point operation counts to compute
 *   the Householder vector of dimension dim+1, on real data.
 *
 * ON ENTRY :
 *   dim      the dimension of the vector must equal the block size;
 *   dimLog2  equals ceil(log2((double) dim), used in sum reduction;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications;
 *   div      current number of divisions;
 *   sqrtfun  current number of calls to sqrt().
 *
 * ON RETURN :
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications;
 *   div      accumulated number of divisions;
 *   sqrtfun  accumulated number of calls to sqrt(). */

void flopcount_cmplx_small_house
 ( int dim, int dimLog2, long long int *add, long long int *mul,
   long long int *div, long long int *sqrtfun );
/*
 * DESCRIPTION :
 *   Returns the accumulated floating-point operation counts to compute
 *   the Householder vector of dimension dim+1, on complex data.
 *
 * ON ENTRY :
 *   dim      the dimension of the vector must equal the block size;
 *   dimLog2  equals ceil(log2((double) dim), used in sum reduction;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications;
 *   div      current number of divisions;
 *   sqrtfun  current number of calls to sqrt().
 *
 * ON RETURN :
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications;
 *   div      accumulated number of divisions;
 *   sqrtfun  accumulated number of calls to sqrt(). */

void flopcount_dbl_large_sum_of_squares
 ( int nblocks, int szt, int sztLog2,
   long long int *add, long long int *mul );
/*
 * DESCRIPTION :
 *   Accumulates the number of floating-point operations to compute the
 *   sums of squares, using multiple blocks, on real data.
 *
 * ON ENTRY :
 *   nblocks   number of blocks;
 *   szt       number of threads in one block;
 *   sztLog2   the 2-log of the number of threads in one block;
 *   add       current number of additions and subtractions;
 *   mul       current number of multiplications.
 *
 * ON RETURN :
 *   add       accumulated number of additions and subtractions;
 *   mul       accumulated number of multiplications. */

void flopcount_cmplx_large_sum_of_squares
 ( int nblocks, int szt, int sztLog2,
   long long int *add, long long int *mul );
/*
 * DESCRIPTION :
 *   Accumulates the number of floating-point operations to compute the
 *   sums of squares, using multiple blocks, on complex data.
 *
 * ON ENTRY :
 *   nblocks   number of blocks;
 *   szt       number of threads in one block;
 *   sztLog2   the 2-log of the number of threads in one block;
 *   add       current number of additions and subtractions;
 *   mul       current number of multiplications.
 *
 * ON RETURN :
 *   add       accumulated number of additions and subtractions;
 *   mul       accumulated number of multiplications. */

void flopcount_dbl_small_leftRupdate
 ( int nrows, int ncols, int szt, int k,
   long long int *add, long long int *mul );
/*
 * DESCRIPTION :
 *   Accumulates the number of floating point operations to update
 *   one tile of the matrix with one block of threads, on real data.
 *
 * ON ENTRY :
 *   nrows    number of rows of R;
 *   ncols    number of columns of R;
 *   szt      size of each block;
 *   k        index of the current column;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications.
 *
 * ON RETURN :
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications. */

void flopcount_dbl_sum_accumulator
 ( int nbt, int nbtLog2, long long int *add );
/*
 * DESCRIPTION :
 *   Accumulates the number of additions to compute the subsums,
 *   with one block of threads.
 *
 * ON ENTRY :
 *   nbt       the number of threads;
 *   nbtLog2   the 2-log of the number of threads;
 *   add       current number of additions.
 *
 * ON RETURN :
 *   add       accumulated number of additions. */

void flopcount_dbl_normalize ( int nblocks, int szt, long long int *div );
/*
 * DESCRIPTION :
 *   Accumulates the number of divisions to divide a vector,
 *   using multiple blocks of threads, on real data.
 *
 * ON ENTRY :
 *   nblocks   the number of blocks;
 *   szt       the number of threads in one block;
 *   div       current number of divisions.
 *
 * ON RETURN :
 *   div       accumulated number of divisions. */

void flopcount_cmplx_normalize
 ( int nblocks, int szt, long long int *add, long long int *mul );
/*
 * DESCRIPTION :
 *   Accumulates the number of divisions to divide a vector,
 *   using multiple blocks of threads, on complex data.
 *
 * ON ENTRY :
 *   nblocks   the number of blocks;
 *   szt       the number of threads in one block;
 *   add       current number of additions;
 *   mul       current number of multiplications.
 *
 * ON RETURN :
 *   add       accumulated number of additions;
 *   mul       accumulated number of multiplications. */

void flopcount_cmplx_small_leftRupdate
 ( int nrows, int ncols, int szt, int k,
   long long int *add, long long int *mul );
/*
 * DESCRIPTION :
 *   Accumulates the number of floating point operations to update
 *   one tile of the matrix with one block of threads, on complex data.
 *
 * ON ENTRY :
 *   nrows    number of rows of R;
 *   ncols    number of columns of R;
 *   szt      size of each block;
 *   k        index of the current column;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications.
 *
 * ON RETURN :
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications. */

void flopcount_dbl_small_betaRTv 
 ( int nrows, int ncols, int szt, int k,
   long long int *add, long long int *mul );
/*
 * DESCRIPTION :
 *   Accumulates the number of floating-point operations to compute
 *   the vector beta*R^T*v, with one block of ncols - k threads,
 *   on real data.
 *
 * ON ENTRY :
 *   nrows    number of rows of R;
 *   ncols    number of columns of R;
 *   szt      size of each block;
 *   k        index of the current column;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications.
 *
 * ON RETURN :
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications. */

void flopcount_cmplx_small_betaRHv 
 ( int nrows, int ncols, int szt, int k,
   long long int *add, long long int *mul );
/*
 * DESCRIPTION :
 *   Accumulates the number of floating-point operations to compute
 *   the vector beta*R^T*v, with one block of ncols - k threads,
 *   on complex data.
 *
 * ON ENTRY :
 *   nrows    number of rows of R;
 *   ncols    number of columns of R;
 *   szt      size of each block;
 *   k        index of the current column;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications.
 *
 * ON RETURN :
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications. */

void flopcount_dbl_RTdotv ( int nrows, int szt, long long int *mul );
/*
 * DESCRIPTION :
 *   Accumulates the number of multiplications to compute RTdotv,
 *   using multiple blocks of threads, on real data.
 *
 * ON ENTRY :
 *   nrows    number of rows of the Householder vector;
 *   szt      number of threads in one block;
 *   mul      current number of multiplications.
 *
 * ON RETURN :
 *   mul      accumulated number of multiplications. */

void flopcount_cmplx_RHdotv
 ( int nrows, int szt, long long int *add, long long int *mul );
/*
 * DESCRIPTION :
 *   Accumulates the number of floating-point operations to compute RHdotv,
 *   using multiple blocks of threads, on complex data.
 *
 * ON ENTRY :
 *   nrows    number of rows of the Householder vector;
 *   szt      number of threads in one block;
 *   add      current number of additions;
 *   mul      current number of multiplications.
 *
 * ON RETURN :
 *   add      accumulated number of additions;
 *   mul      accumulated number of multiplications. */

void flopcount_dbl_sum_betaRTdotv
 ( int nrows, int dim, long long int *add, long long int *mul );
/*
 * DESCRIPTION :
 *   Accumulates the number of additions to compute beta*R^T*v,
 *   given RTdotv, with one block of dim threads, on real data.
 *
 * ON ENTRY :
 *   nrows    number of rows in RTdotv;
 *   dim      dimension of w = beta*R^T*v,
 *            equals the number of threads in the block;
 *   add      current number of additions;
 *   mul      current number of multiplications.
 *
 * ON RETURN :
 *   add      accumulated number of additions;
 *   mul      accumulated number of multiplications. */

void flopcount_cmplx_sum_betaRHdotv
 ( int nrows, int dim, long long int *add, long long int *mul );
/*
 * DESCRIPTION :
 *   Accumulates the number of additions to compute beta*R^H*v,
 *   given RHdotv, with one block of dim threads, on complex data.
 *
 * ON ENTRY :
 *   nrows    number of rows in RHdotv;
 *   dim      dimension of w = beta*R^H*v,
 *            equals the number of threads in the block;
 *   add      current number of additions;
 *   mul      current number of multiplications.
 *
 * ON RETURN :
 *   add      accumulated number of additions;
 *   mul      accumulated number of multiplications. */

void flopcount_dbl_medium_subvbetaRTv
 ( int nrows, int ncols, int szt, int k,
   long long int *add, long long int *mul );
/*
 * DESCRIPTION :
 *   Accumulates the number of floating-point operations to apply
 *   the Householder vector to update R, with multiple blocks,
 *   on real data.
 *
 * ON ENTRY :
 *   nrows    number of rows of R;
 *   ncols    number of columns of R;
 *   szt      size of each block;
 *   k        index of the current column;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications.
 *
 * ON RETURN :
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications. */

void flopcount_cmplx_medium_subvbetaRHv
 ( int nrows, int ncols, int szt, int k,
   long long int *add, long long int *mul );
/*
 * DESCRIPTION :
 *   Accumulates the number of floating-point operations to apply
 *   the Householder vector to update R, with multiple blocks,
 *   on complex data.
 *
 * ON ENTRY :
 *   nrows    number of rows of R;
 *   ncols    number of columns of R;
 *   szt      size of each block;
 *   k        index of the current column;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications.
 *
 * ON RETURN :
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications. */

void flopcount_dbl_VB_to_W
 ( int nrows, int ncols, long long int *add, long long int *mul );
/*
 * DESCRIPTION :
 *   Accumulates the number of floating-point operations to
 *   compute the W, with one block of threads, on real data.
 *
 * ON ENTRY :
 *   nrows    number of rows in the matrices V, Y, and W;
 *   ncols    equals the size of one tile, or equivalently,
 *            is the number of elements in B,
 *            and the number of columns in V, Y, and W;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications.
 *
 * ON RETURN :
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications. */

void flopcount_cmplx_VB_to_W
 ( int nrows, int ncols, long long int *add, long long int *mul );
/*
 * DESCRIPTION :
 *   Accumulates the number of floating-point operations to
 *   compute the W, with one block of threads, on complex data.
 *
 * ON ENTRY :
 *   nrows    number of rows in the matrices V, Y, and W;
 *   ncols    equals the size of one tile, or equivalently,
 *            is the number of elements in B,
 *            and the number of columns in V, Y, and W;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications.
 *
 * ON RETURN :
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications. */

void flopcount_dbl_beta_times_V ( int nrows, long long int *mul );
/*
 * DESCRIPTION :
 *   Accumulates the number of multiplications to multipy beta with the
 *   first Householder vector, with multiple blocks, on real data.
 *
 * ON ENTRY :
 *   nrows    number of rows in V, minus the inserted zeros;
 *   mul      current number of multiplications.
 *
 * ON RETURN :
 *   mul      accumulated number of multiplications. */

void flopcount_cmplx_beta_times_V ( int nrows, long long int *mul );
/*
 * DESCRIPTION :
 *   Accumulates the number of multiplications to multipy beta with the
 *   first Householder vector, with multiple blocks, on complex data.
 *
 * ON ENTRY :
 *   nrows    number of rows in V, minus the inserted zeros;
 *   mul      current number of multiplications.
 *
 * ON RETURN :
 *   mul      accumulated number of multiplications. */

void flopcount_dbl_initialize_WYT ( int dim, long long int *mul );
/*
 * DESCRIPTION :
 *   Accumulates the number of multiplications to initialize W*Y^T,
 *   on real data.
 *
 * ON ENTRY :
 *   dim      dimension of the matrix W*Y^T;
 *   mul      current number of multiplications.
 *
 * ON RETURN :
 *   mul      accumulated number of multiplications. */

void flopcount_cmplx_initialize_WYH
 ( int dim, long long int *add, long long int *mul );
/*
 * DESCRIPTION :
 *   Accumulates the number of floating-point operations to initialize W*Y^H,
 *   on complex data.
 *
 * ON ENTRY :
 *   dim      dimension of the matrix W*Y^H;
 *   add      current number of additions;
 *   mul      current number of multiplications.
 *
 * ON RETURN :
 *   add      accumulated number of additions;
 *   mul      accumulated number of multiplications. */

void flopcount_dbl_update_WYT
 ( int dim, long long int *add, long long int *mul );
/*
 * DESCRIPTION :
 *   Accumulates the number of floating-point operations to update W*Y^T,
 *   on real data.
 *
 * ON ENTRY :
 *   dim      dimension of the matrix W*Y^T;
 *   add      current number of additions;
 *   mul      current number of multiplications.
 *
 * ON RETURN :
 *   add      accumulated number of additions;
 *   mul      accumulated number of multiplications. */

void flopcount_cmplx_update_WYH
 ( int dim, long long int *add, long long int *mul );
/*
 * DESCRIPTION :
 *   Accumulates the number of floating-point operations to update W*Y^H,
 *   on complex data.
 *
 * ON ENTRY :
 *   dim      dimension of the matrix W*Y^H;
 *   add      current number of additions;
 *   mul      current number of multiplications.
 *
 * ON RETURN :
 *   add      accumulated number of additions;
 *   mul      accumulated number of multiplications. */

void flopcount_dbl_beta_next_W
 ( int nrows, long long int *add, long long int *mul );
/*
 * DESCRIPTION :
 *   Accumulates the number of floating-point operations to compute the
 *   next column of W, with multiple blocks, on real data.
 *
 * ON ENTRY :
 *   nrows    the number of rows and the total number of threads;
 *   add      current number of additions;
 *   mul      current number of multipications.
 *
 * ON RETURN :
 *   add      accumulated number of additions;
 *   mul      accumulated number of multiplications. */

void flopcount_cmplx_beta_next_W
 ( int nrows, long long int *add, long long int *mul );
/*
 * DESCRIPTION :
 *   Accumulates the number of floating-point operations to compute the
 *   next column of W, with multiple blocks, on complex data.
 *
 * ON ENTRY :
 *   nrows    the number of rows and the total number of threads;
 *   add      current number of additions;
 *   mul      current number of multipications.
 *
 * ON RETURN :
 *   add      accumulated number of additions;
 *   mul      accumulated number of multiplications. */

void flopcount_dbl_small_WYT
 ( int nrows, int szt, long long int *add, long long int *mul );
/*
 * DESCRIPTION :
 *   Accumulates the number of accumulated floating-point operations
 *   for the multiplication of W with Y^T with multiple blocks,
 *   on real data.
 *
 * ON ENTRY :
 *   nrows    number of rows of all matrices;
 *   szt      number of columns in W and Y,
 *            equals the number of threads in a block;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications.
 *
 * ON RETURN :
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications. */

void flopcount_cmplx_small_WYH
 ( int nrows, int szt, long long int *add, long long int *mul );
/*
 * DESCRIPTION :
 *   Accumulates the number of accumulated floating-point operations
 *   for the multiplication of W with Y^T with multiple blocks,
 *   on complex data.
 *
 * ON ENTRY :
 *   nrows    number of rows of all matrices;
 *   szt      number of columns in W and Y,
 *            equals the number of threads in a block;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications.
 *
 * ON RETURN :
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications. */

void flopcount_dbl_small_QWYT
 ( int dim, int rowdim, int szt, int coloff,
   long long int *add, long long int *mul );
/*
 * DESCRIPTION :
 *   Accumulates the number of floating-point operations to multiply Q
 *   with the WYT, with multiple blocks of threads, on real data.
 *
 * ON ENTRY :
 *   dim      number of rows and columns of the Q matrix;
 *   rowdim   number of rows and columns of the WYT matrix;
 *   szt      the number of threads in a block;
 *   coloff   offset for the column index in QWYT;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications.
 *
 * ON RETURN :
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications. */

void flopcount_cmplx_small_QWYH
 ( int dim, int rowdim, int szt, int coloff,
   long long int *add, long long int *mul );
/*
 * DESCRIPTION :
 *   Accumulates the number of floating-point operations to multiply Q
 *   with the WYH, with multiple blocks of threads, on complex data.
 *
 * ON ENTRY :
 *   dim      number of rows and columns of the Q matrix;
 *   rowdim   number of rows and columns of the WYT matrix;
 *   szt      the number of threads in a block;
 *   coloff   offset for the column index in QWYT;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications.
 *
 * ON RETURN :
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications. */

void flopcount_dbl_small_YWTC
 ( int rowdim, int coldim, long long int *add, long long int *mul );
/*
 * DESCRIPTION :
 *   Accumulates the number of floating-point operations to multiply
 *   YWT with C, with multiple blocks, on real data.
 *
 * ON ENTRY :
 *   rowdim   number of rows minus the row offset;
 *   coldim   number of columns minus the column offset;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications.
 *
 * ON RETURN :
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications. */

void flopcount_cmplx_small_YWHC
 ( int rowdim, int coldim, long long int *add, long long int *mul );
/*
 * DESCRIPTION :
 *   Accumulates the number of floating-point operations to multiply
 *   YWT with C, with multiple blocks, on complex data.
 *
 * ON ENTRY :
 *   rowdim   number of rows minus the row offset;
 *   coldim   number of columns minus the column offset;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications.
 *
 * ON RETURN :
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications. */

void flopcount_dbl_small_Qupdate ( int dim, int rowdim, long long int *add );
/*
 * DESCRIPTION :
 *   Accumulates the number of additions and subtractions to update Q,
 *   with multiple blocks of threads, on real data.
 *
 * ON ENTRY :
 *   dim      number of rows and columns of all matrices;
 *   rowdim   dimension minus the column offset;
 *   add      current number of additions and subtractions.
 *
 * ON RETURN :
 *   add      accumulated number of additions and subtractions. */

void flopcount_cmplx_small_Qupdate
 ( int dim, int rowdim, long long int *add );
/*
 * DESCRIPTION :
 *   Accumulates the number of additions and subtractions to update Q,
 *   with multiple blocks of threads, on complex data.
 *
 * ON ENTRY :
 *   dim      number of rows and columns of all matrices;
 *   rowdim   dimension minus the column offset;
 *   add      current number of additions and subtractions.
 *
 * ON RETURN :
 *   add      accumulated number of additions and subtractions. */

void flopcount_dbl_small_R_add_YWTC
 ( int nrows, int coldim, int szt, int rowoff, int coloff,
   long long int *add );
/*
 * DESCRIPTION :
 *   Accumulates the number of additions and subtractions to update R
 *   by adding YWTC, with multiple blocks, on real data.
 *
 * ON ENTRY :
 *   nrows    total number of rows in R and YWTC;
 *   coldim   number of columns minus the column offset;
 *   szt      the number of threads in a block;
 *   rowoff   offset for the row index in R;
 *   coloff   offset for the column index in R;
 *   add      current number of additions and subtractions;
 *
 * ON RETURN :
 *   add      accumulated number of additions and subtractions. */

void flopcount_cmplx_small_R_add_YWHC
 ( int nrows, int coldim, int szt, int rowoff, int coloff,
   long long int *add );
/*
 * DESCRIPTION :
 *   Accumulates the number of additions and subtractions to update R
 *   by adding YWHC, with multiple blocks, on real data.
 *
 * ON ENTRY :
 *   nrows    total number of rows in R and YWHC;
 *   coldim   number of columns minus the column offset;
 *   szt      the number of threads in a block;
 *   rowoff   offset for the row index in R;
 *   coloff   offset for the column index in R;
 *   add      current number of additions and subtractions;
 *
 * ON RETURN :
 *   add      accumulated number of additions and subtractions. */

#endif
