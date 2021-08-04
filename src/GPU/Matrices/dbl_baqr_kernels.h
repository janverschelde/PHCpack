/* The file dbl_baqr_kernels.h specifies functions for the
 * blocked accelerated QR in double precision. */

#ifndef __dbl_baqr_kernels_h__
#define __dbl_baqr_kernels_h__

#define d_shmemsize 1024

__global__ void dbl_small_house
 ( double *x0, double *x1, int dim, int dimLog2, double *v, double *beta );
/*
 * DESCRIPTION :
 *   Computes the Householder vector of a vector of dimension dim+1,
 *   with one block of dim threads, on real data.
 *
 * ON ENTRY :
 *   x0       first element of the vector x; 
 *   x1       array of dim doubles;
 *   dim      the dimension of the vector must equal the block size;
 *   dimLog2  equals ceil(log2((double) dim), used in sum reduction;
 *   v        space allocated for dim+1 doubles.
 *
 * ON RETURN :
 *   v        the Householder vector;
 *   beta     equals 2/(transpose(v)*v). */

void flopcount_dbl_small_house
 ( int dim, int dimLog2, long int *add, long int *mul, long int *div );
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
 *   div      current number of divisions.
 *
 * ON RETURN :
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications;
 *   div      accumulated number of divisions. */

__global__ void cmplx_small_house
 ( double *x0re, double *x0im, double *x1re, double *x1im,
   int dim, int dimLog2, double *vre, double *vim, double *beta );
/*
 * DESCRIPTION :
 *   Computes the Householder vector of a vector of dimension dim+1,
 *   with one block of dim threads, on complex data.
 *
 * ON ENTRY :
 *   x0re     real part of the first element of the vector x; 
 *   x0im     imaginary part of the first element of the vector x; 
 *   x1re     real parts of an array of dim doubles;
 *   x1im     imaginary parts of an array of dim doubles;
 *   dim      the dimension of the vector must equal the block size;
 *   dimLog2  equals ceil(log2((double) dim), used in sum reduction;
 *   vre      space allocated for dim+1 doubles.
 *   vim      space allocated for dim+1 doubles.
 *
 * ON RETURN :
 *   vre      real parts of the Householder vector;
 *   vim      imaginary parts of the Householder vector;
 *   beta     equals 2/(transpose(v)*v). */

__global__ void dbl_small_leftRupdate
 ( int nrows, int ncols, int szt, int k, double *R, double *v, double *beta );
/*
 * DESCRIPTION :
 *   Updates the matrix R starting at column k
 *   with the Householder vector in v and beta, on real data,
 *   with a total of ncols - k threads in one block.
 *
 * ON ENTRY :
 *   nrows    number of rows of R;
 *   ncols    number of columns of R;
 *   szt      size of each block;
 *   k        index of the current column;
 *   R        an nrows-by-ncols matrix, stored column wise;
 *   v        the Householder vector;
 *   beta     beta corresponding with v.
 *
 * ON RETURN :
 *   R        the updated matrix is trapezoidal. */

void flopcount_small_leftRupdate
 ( int nrows, int ncols, int szt, int k,
   long int *add, long int *mul, long int *div );
/*
 * DESCRIPTION :
 *   Accumulates the number of floating point operations to update
 *   one tile of the matrix with one block of threads.
 *
 * ON ENTRY :
 *   nrows    number of rows of R;
 *   ncols    number of columns of R;
 *   szt      size of each block;
 *   k        index of the current column;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications;
 *   div      current number of divisions.
 *
 * ON RETURN :
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications;
 *   div      accumulated number of divisions. */

__global__ void cmplx_small_leftRupdate
 ( int nrows, int ncols, int szt, int k,
   double *Rre, double *Rim, double *vre, double *vim, double *beta );
/*
 * DESCRIPTION :
 *   Updates the matrix R starting at column k
 *   with the Householder vector in v and beta, on complex data,
 *   with a total of ncols - k threads in one block.
 *
 * ON ENTRY :
 *   nrows    number of rows of R;
 *   ncols    number of columns of R;
 *   szt      size of each block;
 *   k        index of the current column;
 *   Rre      real parts of an nrows-by-ncols matrix, stored column wise.
 *   Rim      imaginary parts of an nrows-by-ncols matrix,
 *            stored column wise;
 *   vre      real parts of the Householder vector;
 *   vim      imaginary parts of the Householder vector;
 *   beta     the beat corresponding with v.
 *
 * ON RETURN :
 *   Rre      real parts of the updated matrix, which is trapezoidal;
 *   Rim      imaginary parts of the updated matrix. */

__global__ void dbl_small_betaRTv
 ( int nrows, int ncols, int szt, int k,
   double *R, double *v, double *beta, double *w );
/*
 * DESCRIPTION :
 *   Computes the vector w = beta*R^T*v, on real data,
 *   with one block of ncols - k threads.
 *
 * ON ENTRY :
 *   nrows    number of rows of R;
 *   ncols    number of columns of R;
 *   szt      size of each block;
 *   k        index of the current column;
 *   R        nrows-by-ncols matrix, stored column wise;
 *   v        the Householder vector;
 *   beta     the beta corresponding with v,
 *   w        space for ncols numbers;
 *
 * ON RETURN :
 *   w        the first ncols-k numbers define beta*R^T*v. */

void flopcount_dbl_small_betaRTv 
 ( int nrows, int ncols, int szt, int k,
   long int *add, long int *mul, long int *div );
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
 *   mul      current number of multiplications;
 *   div      current number of divisions.
 *
 * ON RETURN :
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications;
 *   div      accumulated number of divisions. */

__global__ void cmplx_small_betaRTv
 ( int nrows, int ncols, int szt, int k,
   double *Rre, double *Rim, double *vre, double *vim, double *beta,
   double *wre, double *wim );
/*
 * DESCRIPTION :
 *   Computes the vector w = beta*R^T*v, on complex data,
 *   with one block of ncols - k threads.
 *
 * ON ENTRY :
 *   nrows    number of rows of R;
 *   ncols    number of columns of R;
 *   szt      size of each block;
 *   k        index of the current column;
 *   Rre      real parts of an nrows-by-ncols matrix, stored column wise;
 *   Rim      imaginary parts of an  nrows-by-ncols matrix,
 *            stored column wise;
 *   vre      real parts of the Householder vector;
 *   vim      imaginary parts of the Householder vector;
 *   beta     the beta corresponding with v,
 *   wre      space for ncols numbers;
 *   wim      space for ncols numbers;
 *
 * ON RETURN :
 *   wre      the first ncols-k numbers define
 *            the real parts of beta*R^T*v;
 *   wim      the first ncols-k numbers define
 *            the imaginary parts of beta*R^T*v. */

__global__ void dbl_medium_betaRTv
 ( int nrows, int ncols, int szt, int k,
   double *R, double *v, double *beta, double *w );
/*
 * DESCRIPTION :
 *   Computes the vector w = beta*R^T*v, on real data.
 *
 * REQUIRED : nrows - k > szt as multiple blocks are used.
 *
 * ON ENTRY :
 *   nrows    number of rows of R;
 *   ncols    number of columns of R;
 *   szt      size of each block;
 *   k        index of the current column;
 *   R        nrows-by-ncols matrix, stored column wise;
 *   v        the Householder vector;
 *   beta     the beta corresponding with v,
 *            of length nrows+szt, the extra szt for padding;
 *   w        space for ncols numbers, with additional szt spots for padding.
 *
 * ON RETURN :
 *   w        the first ncols-k numbers define beta*R^T*v. */

__global__ void dbl_medium_subvbetaRTv
 ( int nrows, int ncols, int szt, int k,
   double *R, double *v, double *beta, double *w );
/*
 * DESCRIPTION :
 *   Applies the Householder vector to update R, on real data.
 *
 * REQUIRED : nrows - k > szt as multiple blocks are used.
 *
 * ON ENTRY :
 *   nrows    number of rows of R;
 *   ncols    number of columns of R;
 *   szt      size of each block;
 *   k        index of the current column;
 *   R        nrows-by-ncols matrix, stored column wise;
 *   v        the Householder vector;
 *   beta     the beta corresponding with v,
 *            of length nrows+szt, the extra szt for padding;
 *   w        the vector beta*R^T*v.
 *
 * ON RETURN :
 *   R        the updated matrix is trapezoidal. */

void flopcount_dbl_medium_subvbetaRTv
 ( int nrows, int ncols, int szt, int k,
   long int *add, long int *mul, long int *div );
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
 *   mul      current number of multiplications;
 *   div      current number of divisions.
 *
 * ON RETURN :
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications;
 *   div      accumulated number of divisions. */

__global__ void cmplx_medium_subvbetaRTv
 ( int nrows, int ncols, int szt, int k,
   double *Rre, double *Rim, double *vre, double *vim, double *beta,
   double *wre, double *wim );
/*
 * DESCRIPTION :
 *   Applies the Householder vector to update R, on complex data.
 *
 * REQUIRED : nrows - k > szt as multiple blocks are used.
 *
 * ON ENTRY :
 *   nrows    number of rows of R;
 *   ncols    number of columns of R;
 *   szt      size of each block;
 *   k        index of the current column;
 *   Rre      real parts of the nrows-by-ncols matrix, stored column wise;
 *   Rim      imaginary parts of the nrows-by-ncols matrix,
 *            stored column wise;
 *   vre      real parts of the Householder vector;
 *   vim      imaginary parts of the Householder vector;
 *   beta     the beta corresponding with v,
 *            of length nrows+szt, the extra szt for padding;
 *   wre      real parts of the vector beta*R^T*v;
 *   wim      imaginary parts of the vector beta*R^T*v.
 *
 * ON RETURN :
 *   Rre      real parts of the updated matrix, which is trapezoidal;
 *   Rim      imaginary parts of the updated matrix. */

__global__ void dbl_VB_to_W
 ( int nrows, int ncols, double *B, double *V, double *W );
/*
 * DESCRIPTION :
 *   Computes the W in the WY representation of the Householder
 *   transformations defined by V and B, on real data,
 *   with one block of nrows threads.
 *
 * ON ENTRY :
 *   nrows    number of rows in the matrices V, Y, and W;
 *   ncols    equals the size of one tile, or equivalently,
 *            is the number of elements in B,
 *            and the number of columns in V, Y, and W;
 *   szt      the number of columns in one tile;
 *   B        B[i] is the i-th beta computed by house;
 *   V        V[nrows*i] is the start of the i-th Householder vector,
 *            with i zeros inserted so V is trapezoidal;
 *   W        space for ncols columns with rows from 0 to nrows-1.
 *
 * ON RETURN :
 *   W        the W matrix in the WY representation. */

void flopcount_dbl_VB_to_W
 ( int nrows, int ncols, long int *add, long int *mul, long int* div );
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
 *   mul      current number of multiplications;
 *   div      current number of divisions.
 *
 * ON RETURN :
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications;
 *   div      accumulated number of divisions. */

__global__ void cmplx_VB_to_W
 ( int nrows, int ncols, double *B,
   double *Vre, double *Vim, double *Wre, double *Wim );
/*
 * DESCRIPTION :
 *   Computes the W in the WY representation of the Householder
 *   transformations defined by V and B, on complex data,
 *   with one block of nrows threads.
 *
 * ON ENTRY :
 *   nrows    number of rows in the matrices V, Y, and W;
 *   ncols    equals the size of one tile, or equivalently,
 *            is the number of elements in B,
 *            and the number of columns in V, Y, and W;
 *   B        B[i] is the i-th beta computed by house;
 *   Vre      Vre[nrows*i] is the start of the real parts of the i-th
 *            Householder vector, with i zeros inserted so V is trapezoidal;
 *   Vim      Vim[nrows*i] is the start of the imaginary parts of the i-th
 *            Householder vector, with i zeros inserted so V is trapezoidal;
 *   Wre      space for ncols columns with rows from 0 to nrows-1;
 *   Wim      space for ncols columns with rows from 0 to nrows-1.
 *
 * ON RETURN :
 *   Wre      real parts of the W matrix in the WY representation;
 *   Wim      imaginary parts of the W matrix in the WY representation. */

__global__ void dbl_small_WYT
 ( int nrows, int szt, double *W, double *Y, double *WYT );
/*
 * DESCRIPTION :
 *   Multiplies W with Y^T into the matrix WYT, on real data.
 *   Computes Y*W^T, swapping W with Y in the input arguments.
 *
 * ON ENTRY :
 *   nrows    number of rows of all matrices;
 *   szt      number of columns in W and Y,
 *            equals the number of threads in a block;
 *   W        the W matrix in the WY representation;
 *   Y        the columns of Y are the Householder vectors;
 *   WYT      space for an nrows-by-nrows matrix,
 *            with extra padding for when nrows > szt.
 *
 * ON RETURN :
 *   WYT      the product of W with Y^T. */

void flopcount_dbl_small_WYT
 ( int nrows, int szt, long int *add, long int *mul, long int *div );
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
 *   mul      current number of multiplications;
 *   div      current number of divisions.
 *
 * ON RETURN :
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications;
 *   div      accumulated number of divisions. */

__global__ void cmplx_small_WYT
 ( int nrows, int szt, double *Wre, double *Wim,
   double *Yre, double *Yim, double *WYTre, double *WYTim );
/*
 * DESCRIPTION :
 *   Multiplies W with Y^T into the matrix WYT, on complex data.
 *   Computes Y*W^T, swapping W with Y in the input arguments.
 *
 * ON ENTRY :
 *   nrows    number of rows of all matrices;
 *   szt      number of columns in W and Y,
 *            equals the number of threads in a block;
 *   Wre      real parts of the W matrix in the WY representation;
 *   Wim      imaginary parts of the W matrix in the WY representation;
 *   Yre      real parts of the columns of Y are the Householder vectors;
 *   Yim      imaginary parts of the columns of Y are the Householder vectors;
 *   WYTre    space for an nrows-by-nrows matrix,
 *            with extra padding for when nrows > szt;
 *   WYTim    space for an nrows-by-nrows matrix,
 *            with extra padding for when nrows > szt.
 *
 * ON RETURN :
 *   WYTre    real parts of the product of W with Y^T;
 *   WYTim    imaginary parts of the product of W with Y^T. */

__global__ void dbl_small_QWYT
 ( int dim, int rowdim, int szt, int coloff,
   double *Q, double *WYT, double *QWYT );
/*
 * DESCRIPTION :
 *   Multiplies Q with WYT into the matrix QWYT, on real data.
 *
 * ON ENTRY :
 *   dim      number of rows and columns of the Q matrix;
 *   rowdim   number of rows and columns of the WYT matrix;
 *   szt      the number of threads in a block;
 *   coloff   offset for the column index in QWYT;
 *   Q        the current orthogonal matrix;
 *   WYT      the product of W with Y^T;
 *   QWYT     space for a dim-by-dim matrix.
 *
 * ON RETURN :
 *   QWYT     the product of Q with QWYT. */

void flopcount_dbl_small_QWYT
 ( int dim, int rowdim, int szt, int coloff,
   long int *add, long int *mul, long int *div );
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
 *   mul      current number of multiplications;
 *   div      current number of divisions.
 *
 * ON RETURN :
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications;
 *   div      accumulated number of divisions. */

__global__ void cmplx_small_QWYT
 ( int dim, int rowdim, int szt, int coloff,
   double *Qre, double *Qim, double *WYTre, double *WYTim,
   double *QWYTre, double *QWYTim );
/*
 * DESCRIPTION :
 *   Multiplies Q with WYT into the matrix QWYT, on complex data.
 *
 * ON ENTRY :
 *   dim      number of rows and columns of the Q matrix;
 *   rowdim   number of rows and columns of the WYT matrix;
 *   szt      the number of threads in a block;
 *   coloff   offset for the column index in QWYT;
 *   Qre      real parts of the current orthogonal matrix;
 *   Qim      imaginary parts of the current orthogonal matrix;
 *   WYTre    real parts of the product of W with Y^T;
 *   WYTim    imaginary parts of the product of W with Y^T;
 *   QWYTre   space for a dim-by-dim matrix;
 *   QWYTim   space for a dim-by-dim matrix.
 *
 * ON RETURN :
 *   QWYTre   real parts of the product of Q with QWYT;
 *   QWYTim   imaginary parts of the product of Q with QWYT. */

__global__ void dbl_small_YWTC
 ( int nrows, int ncols, int rowdim, int coldim, int szt,
   int rowoff, int coloff, double *YWT, double *C, double *YWTC );
/*
 * DESCRIPTION :
 *   Multiplies YWT with C into the matrix YWTC, on real data.
 *
 * ON ENTRY :
 *   nrows    total number of rows, nrows = rowdim + rowoff;
 *   ncols    total number of colums, ncols = coldim + coloff;
 *   rowdim   number of rows minus the row offset;
 *   coldim   number of columns minus the column offset;
 *   szt      the number of threads in a block;
 *   rowoff   offset for the row index in C;
 *   coloff   offset for the column index in C;
 *   C        the current matrix to be reduced;
 *   YWT      the product of Y with W^T;
 *   YWTC     space for a rowdim-by-coldim matrix.
 *
 * ON RETURN :
 *   YWTC     the product of YWT with C. */

void flopcount_dbl_small_YWTC
 ( int rowdim, int coldim, long int *add, long int *mul, long int *div );
/*
 * DESCRIPTION :
 *   Accumulates the number of floating-point operations to multiply
 *   YWT with C, with multiple blocks, on real data.
 *
 * ON ENTRY :
 *   rowdim   number of rows minus the row offset;
 *   coldim   number of columns minus the column offset;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications;
 *   div      current number of divisions.
 *
 * ON RETURN :
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications;
 *   div      accumulated number of divisions. */

__global__ void cmplx_small_YWTC
 ( int nrows, int ncols, int rowdim, int coldim, int szt,
   int rowoff, int coloff, double *YWTre, double *YWTim,
   double *Cre, double *Cim, double *YWTCre, double *YWTCim );
/*
 * DESCRIPTION :
 *   Multiplies YWT with C into the matrix YWTC, on complex data.
 *
 * ON ENTRY :
 *   nrows    total number of rows, nrows = rowdim + rowoff;
 *   ncols    total number of colums, ncols = coldim + coloff;
 *   rowdim   number of rows minus the row offset;
 *   coldim   number of columns minus the column offset;
 *   szt      the number of threads in a block;
 *   rowoff   offset for the row index in C;
 *   coloff   offset for the column index in C;
 *   Cre      real parts of the current matrix to be reduced;
 *   Cim      imaginary parts of the current matrix to be reduced;
 *   YWT      real parts of the product of Y with W^T;
 *   YWT      imaginary parts of the product of Y with W^T;
 *   YWTCre   has space for a rowdim-by-coldim matrix;
 *   YWTCim   has space for a rowdim-by-coldim matrix.
 *
 * ON RETURN :
 *   YWTCre   real parts of the product of YWT with C;
 *   YWTCim   imaginary parts of the product of YWT with C. */

__global__ void dbl_small_Qupdate
 ( int dim, int rowdim, int szt, int coloff, double *Q, double *QWYT );
/*
 * DESCRIPTION :
 *   Updates Q by adding QWYT, on real data.
 *
 * ON ENTRY :
 *   dim      number of rows and columns of all matrices;
 *   rowdim   dimension minus the column offset;
 *   szt      the number of threads in a block;
 *   coloff   offset for the column index in QWYT;
 *   Q        the current orthogonal matrix;
 *   QWYT     space for a dim-by-dim matrix.
 *
 * ON RETURN :
 *   Q        Q + QWYT. */

void flopcount_dbl_small_Qupdate
 ( int dim, int rowdim, long int *add, long int *mul, long int *div );
/*
 * DESCRIPTION :
 *   Accumulates the number of floating-point operations to update Q,
 *   with multiple blocks of threads, on real data.
 *
 * ON ENTRY :
 *   dim      number of rows and columns of all matrices;
 *   rowdim   dimension minus the column offset;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications;
 *   div      current number of divisions.
 *
 * ON RETURN :
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications;
 *   div      accumulated number of divisions. */

__global__ void cmplx_small_Qupdate
 ( int dim, int rowdim, int szt, int coloff,
   double *Qre, double *Qim, double *QWYTre, double *QWYTim );
/*
 * DESCRIPTION :
 *   Updates Q by adding QWYT, on complex data.
 *
 * ON ENTRY :
 *   dim      number of rows and columns of all matrices;
 *   rowdim   dimension minus the column offset;
 *   szt      the number of threads in a block;
 *   coloff   offset for the column index in QWYT;
 *   Qre      real parts of the current orthogonal matrix;
 *   Qim      imaginary parts of the current orthogonal matrix;
 *   QWYTre   space for a dim-by-dim matrix;
 *   QWYTim   space for a dim-by-dim matrix.
 *
 * ON RETURN :
 *   Qre      real parts of Q + QWYT;
 *   Qim      imaginary parts Q + QWYT. */

__global__ void dbl_small_R_add_YWTC
 ( int nrows, int coldim, int szt, int rowoff, int coloff,
   double *R, double *YWTC );
/*
 * DESCRIPTION :
 *   Updates R by adding YWTC, on real data.
 *
 * ON ENTRY :
 *   nrows    total number of rows in R and YWTC;
 *   coldim   number of columns minus the column offset;
 *   szt      the number of threads in a block;
 *   rowoff   offset for the row index in R;
 *   coloff   offset for the column index in R;
 *   R        the current matrix to be reduced;
 *   YWTC     the product of YWT with C.
 *
 * ON RETURN :
 *   R        R + YWTC. */

void flopcount_dbl_small_R_and_YWTC
 ( int nrows, int coldim, int szt, int rowoff, int coloff,
   long int *add, long int *mul, long int *div );
/*
 * DESCRIPTION :
 *   Accumulates the number of floating-point operations to update R
 *   by adding YWTC, with multiple blocks, on real data.
 *
 * ON ENTRY :
 *   nrows    total number of rows in R and YWTC;
 *   coldim   number of columns minus the column offset;
 *   szt      the number of threads in a block;
 *   rowoff   offset for the row index in R;
 *   coloff   offset for the column index in R;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications;
 *   div      current number of divisions.
 *
 * ON RETURN :
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications;
 *   div      accumulated number of divisions. */

__global__ void cmplx_small_R_add_YWTC
 ( int nrows, int coldim, int szt, int rowoff, int coloff,
   double *Rre, double *Rim, double *YWTCre, double *YWTCim );
/*
 * DESCRIPTION :
 *   Updates R by adding YWTC, on complex data.
 *
 * ON ENTRY :
 *   nrows    total number of rows in R and YWTC;
 *   coldim   number of columns minus the column offset;
 *   szt      the number of threads in a block;
 *   rowoff   offset for the row index in R;
 *   coloff   offset for the column index in R;
 *   Rre      real parts of the current matrix to be reduced;
 *   Rim      imaginary parts of the current matrix to be reduced;
 *   YWTCre   real parts of the product of YWT with C.
 *   YWTCim   imaginary parts of the product of YWT with C.
 *
 * ON RETURN :
 *   Rre      real parts of R + YWTC;
 *   Rim      imaginary parts R + YWTC. */

void GPU_dbl_small_house
 ( int nrows, int ncols, int szt, int nbt,
   int colidx, int nrows1, int k, int L,
   double *A_h, double *A_d,
   double *v_h, double *V_d, double *beta_h, double *beta_d,
   double *lapms, long int *add, long int *mul, long int *div,
   bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernel to compute the Householder vector for small
 *   enough matrices for the vector to fit entirely in shared memory.
 *
 * ON ENTRY :
 *   nrows    number of rows in the matrix A;
 *   ncols    number of columns in the matrix A;
 *   szt      size of one tile;
 *   nbt      number of tiles, szt*nbt = ncols;
 *   colidx   global index of the current column;
 *   nrows1   number of threads in the block equals the number
 *            of elements computed in the Householder vector;
 *   L        local index of the column in the current tile;
 *   A_h      matrix on the host;
 *   A_d      matrix on the device;
 *   v_h      space for the current Householder vector;
 *   V_d      space for the Householder vectors on the device;
 *   beta_h   space for the betas if verbose;
 *   beta_d   space on the device for the betas;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications;
 *   div      current number of divisions;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   v_h      the next Householder vector on the host, if verbose;
 *   V_d      contains the next computed Householder vector;
 *   beta_h   updated vector of betas, if verbose;
 *   beta_d   the next beta constant;
 *   lapms    elapsed time spent by the kernel;
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications;
 *   div      accumulated number of divisions. */

void GPU_cmplx_small_house
 ( int nrows, int ncols, int szt, int nbt,
   int colidx, int nrows1, int k, int L,
   double *Are_h, double *Aim_h, double *Are_d, double *Aim_d,
   double *vre_h, double *vim_h, double *Vre_d, double *Vim_d,
   double *beta_h, double *beta_d, double *lapms, bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernel to compute the Householder vector for small
 *   enough matrices for the vector to fit entirely in shared memory.
 *
 * ON ENTRY :
 *   nrows    number of rows in the matrix A;
 *   ncols    number of columns in the matrix A;
 *   szt      size of one tile;
 *   nbt      number of tiles, szt*nbt = ncols;
 *   colidx   global index of the current column;
 *   nrows1   number of threads in the block equals the number
 *            of elements computed in the Householder vector;
 *   L        local index of the column in the current tile;
 *   Are_h    real parts of the matrix on the host;
 *   Aim_h    imaginary parts of the matrix on the host;
 *   Are_d    real parts of the matrix on the device;
 *   Aim_d    imaginary parts of the matrix on the device;
 *   vre_h    space for the real parts of the current Householder vector;
 *   vim_h    space for the imaginary parts of the current Householder vector;
 *   Vre_d    space for the real parts of the Householder vectors
 *            on the device;
 *   Vim_d    space for the imaginary parts of the Householder vectors
 *            on the device;
 *   beta_h   space for the betas if verbose;
 *   beta_d   space on the device for the betas;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   v_h      the next Householder vector on the host, if verbose;
 *   Vre_d    real parts of the next computed Householder vector;
 *   Vim_d    imaginary parts of the next computed Householder vector;
 *   beta_h   updated vector of betas, if verbose;
 *   beta_d   the next beta constant;
 *   lapms    elapsed time spent by the kernel. */

void GPU_dbl_small_leftRupdate
 ( int nrows, int ncols, int szt, int colidx, int k, int L,
   double *A_h, double *A_d, double *V_d, double *beta_h, double *beta_d,
   double *lapms, long int *add, long int *mul, long int *div,
   bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernel to update one tile.
 *   Wraps the timer and the print statements if verbose.
 *   If verbose, then the reduced matrix is returned on the host.
 *
 * ON ENTRY :
 *   nrows    number of rows in the matrix A;
 *   ncols    number of columns in the matrix A;
 *   szt      size of one tile;
 *   colidx   global index of the current column;
 *   k        index of the current tile;
 *   L        local index of the column in the current tile;
 *   A_h      matrix on the host;
 *   A_d      matrix on the device;
 *   V_d      space for the Householder vectors on the device;
 *   beta_h   space for the betas if verbose;
 *   beta_d   space on the device for the betas;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications;
 *   div      current number of divisions;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   V_d      contains the Householder vectors on the device;
 *   beta_h   vector of betas, if verbose;
 *   beta_d   the next beta constant;
 *   lapms    elapsed time spent by the kernel.
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications;
 *   div      accumulated number of divisions. */

void GPU_cmplx_small_leftRupdate
 ( int nrows, int ncols, int szt, int colidx, int k, int L,
   double *Are_h, double *Aim_h, double *Are_d, double *Aim_d,
   double *Vre_d, double *Vim_d, double *beta_h, double *beta_d,
   double *lapms, bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernel to update one tile.
 *   Wraps the timer and the print statements if verbose.
 *   If verbose, then the reduced matrix is returned on the host.
 *
 * ON ENTRY :
 *   nrows    number of rows in the matrix A;
 *   ncols    number of columns in the matrix A;
 *   szt      size of one tile;
 *   colidx   global index of the current column;
 *   k        index of the current tile;
 *   L        local index of the column in the current tile;
 *   Are_h    real parts of the matrix on the host;
 *   Aim_h    imaginary parts of the matrix on the host;
 *   Are_d    real parts of the matrix on the device;
 *   Aim_d    imaginary parts of the matrix on the device;
 *   Vre_d    space for the real parts of the Householder vectors
 *            on the device;
 *   Vim_d    space for the imaginary parts of the Householder vectors
 *            on the device;
 *   beta_h   space allocated for the betas if verbose;
 *   beta_d   space allocated on the device for the betas;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Vre_d    real parts of the Householder vectors on the device;
 *   Vim_d    imaginary parts of the Householder vectors on the device;
 *   beta_h   vector of betas, if verbose;
 *   beta_d   the next beta constant;
 *   lapms    elapsed time spent by the kernel. */

void GPU_dbl_medium_leftRupdate
 ( int nrows, int ncols, int szt, int colidx, int k, int L,
   double *A_h, double *A_d, double *V_d, double *beta_h, double *beta_d,
   double *w_h, double *w_d,
   double *lapms, long int *add, long int *mul, long int *div,
   bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernels to update one tile.
 *   Wraps the timer and the print statements if verbose.
 *   If verbose, then the reduced matrix is returned on the host.
 *
 * REQUIRED : nrows - colidx > szt.
 *
 * ON ENTRY :
 *   nrows    number of rows in the matrix A;
 *   ncols    number of columns in the matrix A;
 *   szt      size of one tile;
 *   colidx   global index of the current column;
 *   k        index of the current tile;
 *   L        local index of the column in the current tile;
 *   A_h      matrix on the host;
 *   A_d      matrix on the device;
 *   V_d      space for the Householder vectors on the device;
 *   beta_h   space for the betas if verbose;
 *   beta_d   space on the device for the betas;
 *   w_h      space for the beta*R^T*v on the host, if verbose;
 *   w_d      space for the beta*R^T*v plus szt padding;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications;
 *   div      current number of divisions;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   V_d      contains the Householder vectors on the device;
 *   beta_h   vector of betas, if verbose;
 *   beta_d   the next beta constant;
 *   w_h      stores the beta*R^T*v, if verbose;
 *   lapms    elapsed time spent by the kernel;
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications;
 *   div      accumulated number of divisions. */

void GPU_cmplx_medium_leftRupdate
 ( int nrows, int ncols, int szt, int colidx, int k, int L,
   double *Are_h, double *Aim_h, double *Are_d, double *Aim_d,
   double *Vre_d, double *Vim_d, double *beta_h, double *beta_d,
   double *wre_h, double *wim_h, double *wre_d, double *wim_d,
   double *lapms, bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernels to update one tile.
 *   Wraps the timer and the print statements if verbose.
 *   If verbose, then the reduced matrix is returned on the host.
 *
 * REQUIRED : nrows - colidx > szt.
 *
 * ON ENTRY :
 *   nrows    number of rows in the matrix A;
 *   ncols    number of columns in the matrix A;
 *   szt      size of one tile;
 *   colidx   global index of the current column;
 *   k        index of the current tile;
 *   L        local index of the column in the current tile;
 *   Are_h    real parts of the matrix on the host;
 *   Aim_h    imaginary parts of the matrix on the host;
 *   Are_d    real parts of the matrix on the device;
 *   Aim_d    imaginary parts of the matrix on the device;
 *   Vre_d    space for the real parts of the Householder vectors
 *            on the device;
 *   Vim_d    space for the imaginary parts of the Householder vectors
 *            on the device;
 *   beta_h   space allocated for the betas if verbose;
 *   beta_d   space allocated on the device for the betas;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Vre_d    real parts of the Householder vectors on the device;
 *   Vim_d    imaginary parts of the Householder vectors on the device;
 *   beta_h   vector of betas, if verbose;
 *   beta_d   the next beta constant;
 *   lapms    elapsed time spent by the kernel. */

void GPU_dbl_VB_to_W
 ( int nrows, int ncols, int szt,
   double *V_h, double *V_d, double *W_h, double *W_d,
   double *beta_h, double *beta_d, double *lapms,
   long int *add, long int *mul, long int *div, bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernel to compute the W in the WY representation.
 *   Wraps the timer and the print statements if verbose.
 *   If verbose, then the W matrix is returned on the host.
 *
 * ON ENTRY :
 *   nrows    number of rows in the matrices V, Y, and W;
 *   ncols    equals the size of one tile, or equivalently,
 *            is the number of elements in B,
 *            and the number of columns in V, Y, and W;
 *   V_h      space for the Householder vectors, if verbose;
 *   V_d      space for V on the device;
 *   W_h      space for the W matrix, if verbose;
 *   W_d      space for W on the device;
 *   beta_h   space for the betas if verbose;
 *   beta_d   space on the device for the betas;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications;
 *   div      current number of divisions;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   V_h      equals the Y matrix, if verbose;
 *   V_d      the Y matrix on the device;
 *   W_d      the W matrix in the WY representation, on the device;
 *   W_h      the W matrix in the WY representation, if verbose;
 *   lapms    elapsed time spent by the kernel;
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications;
 *   div      accumulated number of divisions. */

void GPU_cmplx_VB_to_W
 ( int nrows, int ncols, int szt,
   double *Vre_h, double *Vim_h, double *Vre_d, double *Vim_d,
   double *Wre_h, double *Wim_h, double *Wre_d, double *Wim_d,
   double *beta_h, double *beta_d, double *lapms, bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernel to compute the W in the WY representation.
 *   Wraps the timer and the print statements if verbose.
 *   If verbose, then the W matrix is returned on the host.
 *
 * ON ENTRY :
 *   nrows    number of rows in the matrices V, Y, and W;
 *   ncols    equals the size of one tile, or equivalently,
 *            is the number of elements in B,
 *            and the number of columns in V, Y, and W;
 *   Vre_h    space for the real parts of the Householder vectors,
 *            if verbose;
 *   Vim_h    space for the imaginary parts of the Householder vectors,
 *            if verbose;
 *   Vre_d    space for the real parts of V on the device;
 *   Vim_d    space for the imaginary parts of V on the device;
 *   Wre_h    space for the real parts of the W matrix, if verbose;
 *   Wim_h    space for the imaginary parts of the W matrix, if verbose;
 *   Wre_d    space for the real parts of W on the device;
 *   Wim_d    space for the imaginary parts of W on the device;
 *   beta_h   space for the betas if verbose;
 *   beta_d   space on the device for the betas;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Vre_h    equals the real parts of the Y matrix, if verbose;
 *   Vim_h    equals the imaginary parts of the Y matrix, if verbose;
 *   Vre_d    the real parts of the Y matrix on the device;
 *   Vim_d    the imaginary parts of the Y matrix on the device;
 *   Wre_d    the real parts of the W matrix in the WY representation,
 *            on the device;
 *   Wim_d    the imaginary parts of the W matrix in the WY representation,
 *            on the device;
 *   Wre_h    the real parts of the W matrix in the WY representation,
 *            if verbose;
 *   Wim_h    the imaginary parts of the W matrix in the WY representation,
 *            if verbose;
 *   lapms    elapsed time spent by the kernel. */

void GPU_dbl_small_WYT
 ( int nrows, int szt, double *W_d, double *Y_d, double *WYT_d,
   double *WYT_h, double *lapms, long int *add, long int *mul, long int *div,
   bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernel to compute W*Y^T.
 *   Wraps the timer and the print statements if verbose.
 *   If verbose, then the W*Y^T matrix is returned.
 *
 * ON ENTRY :
 *   nrows    number of rows in all matrices;
 *   szt      size of one tile and the number of threads in a block;
 *   W_d      the matrix W in the WY representation;
 *   Y_d      the matrix of Householder vectors;
 *   WYT_d    space for an nrows-by-nrows matrix on the device;
 *   WYT_h    space for an nrows-by-nrows matrix on the host, if verbose;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications;
 *   div      current number of divisions;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   WYT_d    the product W*Y^T on the device;
 *   WYT_h    the product W*Y^T, if verbose;
 *   lapms    elapsed time spent by the kernel;
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications;
 *   div      accumulated number of divisions. */

void GPU_cmplx_small_WYT
 ( int nrows, int szt, double *Wre_d, double *Wim_d,
   double *Yre_d, double *Yim_d, double *WYTre_d, double *WYTim_d,
   double *WYTre_h, double *WYTim_h, double *lapms, bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernel to compute W*Y^T.
 *   Wraps the timer and the print statements if verbose.
 *   If verbose, then the W*Y^T matrix is returned.
 *
 * ON ENTRY :
 *   nrows    number of rows in all matrices;
 *   szt      size of one tile and the number of threads in a block;
 *   Wre_d    the real parts of the matrix W in the WY representation;
 *   Wim_d    the imaginary parts of the matrix W in the WY representation;
 *   Yre_d    the matrix of the real parts of the Householder vectors;
 *   Yim_d    the matrix of the imaginary parts of the Householder vectors;
 *   WYTre_d  has space for an nrows-by-nrows matrix on the device;
 *   WYTim_d  has space for an nrows-by-nrows matrix on the device;
 *   WYTre_h  has space for an nrows-by-nrows matrix on the host, if verbose;
 *   WYTim_h  has space for an nrows-by-nrows matrix on the host, if verbose;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   WYTre_d  are the real parts of the product W*Y^T on the device;
 *   WYTim_d  are the imaginary parts of the product W*Y^T on the device;
 *   WYTre_h  are the real parts of the product W*Y^T, if verbose;
 *   WYTim_h  are the imaginary parts of the product W*Y^T, if verbose;
 *   lapms    elapsed time spent by the kernel. */

void GPU_dbl_small_YWT
 ( int nrows, int szt, int idx, double *Y_d, double *W_d, double *YWT_d,
   double *YWT_h, double *lapms, long int *add, long int *mul, long int *div,
   bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernel to compute W*Y^T.
 *   Wraps the timer and the print statements if verbose.
 *   If verbose, then the W*Y^T matrix is returned.
 *
 * ON ENTRY :
 *   nrows    number of rows in all matrices;
 *   szt      size of one tile and the number of threads in a block;
 *   idx      index of the current tile;
 *   Y_d      the matrix of Householder vectors;
 *   W_d      the matrix W in the WY representation;
 *   YWT_d    space for an nrows-by-nrows matrix on the device;
 *   YWT_h    space for an nrows-by-nrows matrix on the host, if verbose;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications;
 *   div      current number of divisions;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   YWT_d    the product Y*W^T on the device;
 *   YWT_h    the product Y*W^T, if verbose;
 *   lapms    elapsed time spent by the kernel;
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications;
 *   div      accumulated number of divisions. */

void GPU_cmplx_small_YWT
 ( int nrows, int szt, int idx,
   double *Yre_d, double *Yim_d, double *Wre_d, double *Wim_d,
   double *YWTre_d, double *YWTim_d, double *YWTre_h, double *YWTim_h,
   double *lapms, bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernel to compute W*Y^T.
 *   Wraps the timer and the print statements if verbose.
 *   If verbose, then the W*Y^T matrix is returned.
 *
 * ON ENTRY :
 *   nrows    number of rows in all matrices;
 *   szt      size of one tile and the number of threads in a block;
 *   idx      index of the current tile;
 *   Yre_d    real parts of the matrix of Householder vectors;
 *   Yim_d    imaginary parts of the matrix of Householder vectors;
 *   Wre_d    real parts of the matrix W in the WY representation;
 *   Wim_d    imaginary parts of the matrix W in the WY representation;
 *   YWTre_d  space for an nrows-by-nrows matrix on the device;
 *   YWTim_d  space for an nrows-by-nrows matrix on the device;
 *   YWTre_h  space for an nrows-by-nrows matrix on the host, if verbose;
 *   YWTim_h  space for an nrows-by-nrows matrix on the host, if verbose;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   YWTre_d  real parts of the product Y*W^T on the device;
 *   YWTim_d  imaginary parts of the product Y*W^T on the device;
 *   YWTre_h  real parts of the product Y*W^T, if verbose;
 *   YWTim_h  imaginary parts of the product Y*W^T, if verbose;
 *   lapms    elapsed time spent by the kernel. */

void GPU_dbl_small_QWYT
 ( int dim, int szt, int idx, double *Q_d, double *WYT_d, double *QWYT_d,
   double *QWYT_h, double *Q_h,
   double *lapms, long int *add, long int *mul, long int *div,
   bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernel to compute Q*WYT.
 *   Wraps the timer and the print statement if verbose.
 *   If verbose, then the Q*WYT matrix is returned.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in Q and WYT;
 *   szt      size of one tile and the number of threads in a block;
 *   idx      index of the current tile;
 *   Q_d      a dim-by-dim matrix, on the device;
 *   WYT_d    the product W*Y^T, on the device;
 *   QWYT_d   space for the product Q*WYT, on the device;
 *   QWYT_h   space for the product Q*WYT, on the host, if verbose;
 *   Q_h      if verbose, then used to print Q before the product;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications;
 *   div      current number of divisions;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   QWYT_d   the product Q*WYT on the device;
 *   QWYT_h   the product Q*WYT, if verbose;
 *   lapms    elapsed time spent by the kernel;
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications;
 *   div      accumulated number of divisions. */

void GPU_cmplx_small_QWYT
 ( int dim, int szt, int idx, double *Qre_d, double *Qim_d,
   double *WYTre_d, double *WYTim_d, double *QWYTre_d, double *QWYTim_d,
   double *QWYTre_h, double *QWYTim_h, double *Qre_h, double *Qim_h,
   double *lapms, bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernel to compute Q*WYT.
 *   Wraps the timer and the print statement if verbose.
 *   If verbose, then the Q*WYT matrix is returned.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in Q and WYT;
 *   szt      size of one tile and the number of threads in a block;
 *   idx      index of the current tile;
 *   Qre_d    real parts of a dim-by-dim matrix, on the device;
 *   Qim_d    imaginary parts of a dim-by-dim matrix, on the device;
 *   WYTre_d  are the real parts of the product W*Y^T, on the device;
 *   WYTim_d  are the imaginary parts of the product W*Y^T, on the device;
 *   QWYTre_d has space for the real parts for the product Q*WYT,
 *            on the device;
 *   QWYTim_d has space for the imaginary parts of the product Q*WYT,
 *            on the device;
 *   QWYTre_h has space for the real parts of the product Q*WYT,
 *            on the host, if verbose;
 *   QWYTim_h has space for the imaginary parts of the product Q*WYT,
 *            on the host, if verbose;
 *   Qre_h    if verbose, then used to print the real parts
 *            of the matrix Q before the product;
 *   Qim_h    if verbose, then used to print the imaginary parts
 *            of the matrix Q before the product;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   QWYTre_d are the real parts of the product Q*WYT on the device;
 *   QWYTim_d are the imaginary parts of the product Q*WYT on the device;
 *   QWYTre_h are the real parts of the product Q*WYT, if verbose;
 *   QWYTim_h are the imaginary parts of the product Q*WYT, if verbose;
 *   lapms    elapsed time spent by the kernel. */

void GPU_dbl_small_YWTC
 ( int nrows, int ncols, int szt, int idx,
   double *YWT_d, double *C_d, double *YWTC_d, double *YWTC_h,
   double *lapms, long int *add, long int *mul, long int *div,
   bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernel to compute YWT*C.
 *   Wraps the timer and the print statement if verbose.
 *   If verbose, then the YWT*C matrix is returned.
 *
 * ON ENTRY :
 *   nrows    number of rows of the matrices C, YWT, YWTC,
 *            and the number of columns of the matrix YWT;
 *   ncols    number of columns of the matrix C;
 *   szt      size of one tile and the number of threads in a block;
 *   idx      index of the current tile;
 *   YWT_d    the product Y*W^T, on the device;
 *   C_d      an nrows-by-ncols matrix, on the device;
 *   YWTC_d   space for the product YWT*C, on the device;
 *   YWTC_h   space for the product YWT*C, on the host, if verbose;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications;
 *   div      current number of divisions;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   YWTC_d   the product YWT*C on the device;
 *   YWTC_h   the product YWT*C, if verbose;
 *   lapms    elapsed time spent by the kernel;
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications;
 *   div      accumulated number of divisions. */

void GPU_cmplx_small_YWTC
 ( int nrows, int ncols, int szt, int idx,
   double *YWTre_d, double *YWTim_d, double *Cre_d, double *Cim_d,
   double *YWTCre_d, double *YWTCim_d, double *YWTCre_h, double *YWTCim_h,
   double *lapms, bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernel to compute YWT*C.
 *   Wraps the timer and the print statement if verbose.
 *   If verbose, then the YWT*C matrix is returned.
 *
 * ON ENTRY :
 *   nrows    number of rows of the matrices C, YWT, YWTC,
 *            and the number of columns of the matrix YWT;
 *   ncols    number of columns of the matrix C;
 *   szt      size of one tile and the number of threads in a block;
 *   idx      index of the current tile;
 *   YWTre_d  are the real parts of the product Y*W^T, on the device;
 *   YWTim_d  are the imaginary parts of the product Y*W^T, on the device;
 *   Cre_d    real parts of an nrows-by-ncols matrix, on the device;
 *   Cim_d    imaginary parts of an nrows-by-ncols matrix, on the device;
 *   YWTCre_d has space for the real parts of the product YWT*C,
 *            on the device;
 *   YWTCim_d has space for the imaginary parts of the product YWT*C,
 *            on the device;
 *   YWTCre_h has space for the real parts of the product YWT*C,
 *            on the host, if verbose;
 *   YWTCim_h has space for the imaginary parts of the product YWT*C,
 *            on the host, if verbose;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   YWTCre_d are the real parts of the product YWT*C on the device;
 *   YWTCim_d are the imaginary parts of the product YWT*C on the device;
 *   YWTCre_h are the real parts of the product YWT*C, if verbose;
 *   YWTCim_h are the imaginary parts of the product YWT*C, if verbose;
 *   lapms    elapsed time spent by the kernel. */

void GPU_dbl_small_Qupdate
 ( int dim, int rowdim, int szt, int idx, double *Q_d, double *QWYT_d,
   double *Q_h, double *lapms, long int *add, long int *mul, long int *div,
   bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernel to update Q as Q + QWYT.
 *   Wraps the timer and the print statement if verbose.
 *   If verbose, then the updated Q is returned.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in Q,
 *            number of rows in QWYT;
 *   rowdim   number of columns in QWYT;
 *   szt      size of one tile and the number of threads in a block;
 *   idx      index of the current tile;
 *   Q_d      a dim-by-dim matrix, on the device;
 *   QWYT_d   the product Q*W*Y^T, on the device;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications;
 *   div      current number of divisions;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Q_d      the updated Q on the device;
 *   Q_h      the updated Q, if verbose;
 *   lapms    elapsed time spent by the kernel;
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications;
 *   div      accumulated number of divisions. */

void GPU_cmplx_small_Qupdate
 ( int dim, int szt, int idx, double *Qre_d, double *Qim_d,
   double *QWYTre_d, double *QWYTim_d, double *Qre_h, double *Qim_h,
   double *lapms, bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernel to update Q as Q + QWYT.
 *   Wraps the timer and the print statement if verbose.
 *   If verbose, then the updated Q is returned.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in Q,
 *            number of rows in QWYT;
 *   rowdim   number of columns in QWYT;
 *   szt      size of one tile and the number of threads in a block;
 *   idx      index of the current tile;
 *   Qre_d    real parts of a dim-by-dim matrix, on the device;
 *   Qim_d    imaginary parts of a dim-by-dim matrix, on the device;
 *   QWYTre_d are the real parts of the product Q*W*Y^T, on the device;
 *   QWYTim_d are the imaginary parts of the product Q*W*Y^T, on the device;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Qre_d    real parts of the updated Q on the device;
 *   Qim_d    imaginary parts of the updated Q on the device;
 *   Qre_h    real parts of the updated Q, if verbose;
 *   Qim_h    imaginary parts of the updated Q, if verbose;
 *   lapms    elapsed time spent by the kernel. */

void GPU_dbl_small_R_add_YWTC
 ( int nrows, int ncols, int szt, int idx, double *R_d, double *YWTC_d,
   double *R_h, double *lapms, long int *add, long int *mul, long int *div,
   bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernel to update R as R + YWTC.
 *   Wraps the timer and the print statement if verbose.
 *   If verbose, then the updated R is returned.
 *
 * ON ENTRY :
 *   nrows    total number of rows in R and YWTC;
 *   ncols    total number of columns in R and YWTC;
 *   szt      size of one tile and the number of threads in a block;
 *   idx      index of the current tile;
 *   R_d      an nrows-by-ncols matrix, on the device;
 *   YWTC_d   the product Y*W^T*C, on the device;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications;
 *   div      current number of divisions;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   R_d      the updated R on the device;
 *   R_h      the updated R, if verbose;
 *   lapms    elapsed time spent by the kernel;
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications;
 *   div      accumulated number of divisions. */

void GPU_cmplx_small_R_add_YWTC
 ( int nrows, int ncols, int szt, int idx,
   double *Rre_d, double *Rim_d, double *YWTCre_d, double *YWTCim_d,
   double *Rre_h, double *Rim_h, double *lapms, bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernel to update R as R + YWTC.
 *   Wraps the timer and the print statement if verbose.
 *   If verbose, then the updated R is returned.
 *
 * ON ENTRY :
 *   nrows    total number of rows in R and YWTC;
 *   ncols    total number of columns in R and YWTC;
 *   szt      size of one tile and the number of threads in a block;
 *   idx      index of the current tile;
 *   Rre_d    real parts of an nrows-by-ncols matrix, on the device;
 *   Rim_d    imaginary parts of an nrows-by-ncols matrix, on the device;
 *   YWTCre_d are the real parts of the product Y*W^T*C, on the device;
 *   YWTCim_d are the imaginary parts of the product Y*W^T*C,
 *            on the device;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Rre_d    real parts of the updated R on the device;
 *   Rim_d    imaginary parts of the updated R on the device;
 *   Rre_h    real parts of the updated R, if verbose;
 *   Rim_h    imaginary parts of the updated R, if verbose;
 *   lapms    elapsed time spent by the kernel. */

void GPU_dbl_blocked_houseqr
 ( int nrows, int ncols, int szt, int nbt,
   double **A, double **Q, double **R,
   double *houselapms, double *tileRlapms, double *vb2Wlapms,
   double *WYTlapms, double *QWYTlapms, double *Qaddlapms,
   double *YWTlapms, double *YWTClapms, double *Raddlapms,
   double *walltimesec, long int *addcnt, long int *mulcnt, long int *divcnt,
   bool verbose=true );
/*
 * DESCRIPTION :
 *   Applies Householder transformations in a blocked manner
 *   to compute a QR decomposition of A, on real data.
 *
 * REQUIRED : nrows >= ncols.
 *
 * ON ENTRY :
 *   nrows    number of rows of A;
 *   ncols    number of columns of A;
 *   szt      size of each block;
 *   nbt      number of tiles, ncols = szt*nbt;
 *   A        an nrows-by-ncols matrix,
 *            stored as nrows arrays of ncols numbers;
 *   Q        space for an nrows-by-nrows matrix;
 *   R        space for an nrows-by-ncols matrix;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Q        an orthogonal matrix, transpose(Q)*A = R;
 *   R        the reduced upper triangular form of A;
 *   houselapms is the elapsed time spent by the kernel
 *            to compute the Householder vector and the beta;
 *   tileRlapms is the elapsed time spent by the kernel
 *            to reduce one tile;
 *   vb2Wlapms is the elapsed time spent by the kernel
 *            to compute the W representation;
 *   WYTlapms is the elapsed time spent by the kernel
 *            to compute the W*Y^T product;
 *   QWYTlapms is the elapsed time spent by the kernel
 *            to compute the Q*WYT product;
 *   Qaddlapms is the elapsed time spent by the kernel
 *            to compute Q by adding the Q*W*Y^T matrix;
 *   YWTlapms is the elapsed time spent by the kernel
 *            to compute the Y*W^T product;
 *   YWTClapms is the elapsed time spent by the kernel
 *            to compute the YWT*C product;
 *   Raddlapms is the elapsed time spent by the kernel
 *            to compute R by adding the Y*W^T*C matrix;
 *   walltimesec is the elapsed wall clock computation time. */

void GPU_cmplx_blocked_houseqr
 ( int nrows, int ncols, int szt, int nbt,
   double **Are, double **Aim, double **Qre, double **Qim,
   double **Rre, double **Rim,
   double *houselapms, double *tileRlapms, double *vb2Wlapms,
   double *WYTlapms, double *QWYTlapms, double *Qaddlapms,
   double *YWTlapms, double *YWTClapms, double *Raddlapms,
   double *walltimesec, bool verbose=true );
/*
 * DESCRIPTION :
 *   Applies Householder transformations in a blocked manner
 *   to compute a QR decomposition of A, on complex data.
 *
 * REQUIRED : nrows >= ncols.
 *
 * ON ENTRY :
 *   nrows    number of rows of A;
 *   ncols    number of columns of A;
 *   szt      size of each block;
 *   nbt      number of tiles, ncols = szt*nbt;
 *   Are      real parts of an nrows-by-ncols matrix,
 *            stored as nrows arrays of ncols numbers;
 *   Aim      imaginary parts of an nrows-by-ncols matrix,
 *            stored as nrows arrays of ncols numbers;
 *   Qre      space for an nrows-by-nrows matrix;
 *   Qim      space for an nrows-by-nrows matrix;
 *   Rre      space for an nrows-by-ncols matrix;
 *   Rim      space for an nrows-by-ncols matrix;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Qre      real parts of an orthogonal matrix, transpose(Q)*A = R;
 *   Qim      imaginary parts of the orthogonal matrix Q;
 *   Rre      real parts of the reduced upper triangular form of A;
 *   Rim      imaginary parts of the reduced upper triangular form of A;
 *   houselapms is the elapsed time spent by the kernel
 *            to compute the Householder vector and the beta;
 *   tileRlapms is the elapsed time spent by the kernel
 *            to reduce one tile;
 *   vb2Wlapms is the elapsed time spent by the kernel
 *            to compute the W representation;
 *   WYTlapms is the elapsed time spent by the kernel
 *            to compute the W*Y^T product;
 *   QWYTlapms is the elapsed time spent by the kernel
 *            to compute the Q*WYT product;
 *   Qaddlapms is the elapsed time spent by the kernel
 *            to compute Q by adding the Q*W*Y^T matrix;
 *   YWTlapms is the elapsed time spent by the kernel
 *            to compute the Y*W^T product;
 *   YWTClapms is the elapsed time spent by the kernel
 *            to compute the YWT*C product;
 *   Raddlapms is the elapsed time spent by the kernel
 *            to compute R by adding the Y*W^T*C matrix;
 *   walltimesec is the elapsed wall clock computation time. */

#endif
