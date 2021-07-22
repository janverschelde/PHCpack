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
 *   with one block of dim threads.
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

__global__ void dbl_factors_leftRupdate
 ( int nrows, int ncols, int szt, int k, double *R, double *v, double *beta );
/*
 * DESCRIPTION :
 *   Updates the matrix R starting at column k
 *   with the Householder vector in v and beta,
 *   with a total of ncols - k threads 
 *   and (ncols - k)/szt + 1 number of blocks.
 *
 * NOTE :
 *   This kernel can only work for one block.
 *   because thread blocks do not synchronize.
 *
 * ON ENTRY :
 *   nrows    number of rows of R;
 *   ncols    number of columns of R;
 *   szt      size of each block;
 *   k        index of the current column;
 *   R        an nrows-by-ncols matrix, stored column wise.
 *
 * ON RETURN :
 *   R        the updated matrix is trapezoidal. */

__global__ void dbl_small_leftRupdate
 ( int nrows, int ncols, int szt, int k, double *R, double *v, double *beta );
/*
 * DESCRIPTION :
 *   Updates the matrix R starting at column k
 *   with the Householder vector in v and beta,
 *   with a total of ncols - k threads in one block.
 *
 * ON ENTRY :
 *   nrows    number of rows of R;
 *   ncols    number of columns of R;
 *   szt      size of each block;
 *   k        index of the current column;
 *   R        an nrows-by-ncols matrix, stored column wise.
 *
 * ON RETURN :
 *   R        the updated matrix is trapezoidal. */

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
 *   B        B[i] is the i-th beta computed by house;
 *   V        V[nrows*i] is the start of the i-th Householder vector,
 *            with i zeros inserted so V is trapezoidal;
 *   W        space for ncols columns with rows from 0 to nrows-1.
 *
 * ON RETURN :
 *   W        the W matrix in the WY representation. */

__global__ void dbl_small_WYT
 ( int nrows, int szt, double *W, double *Y, double *WYT );
/*
 * DESCRIPTION :
 *   Multiplies W with Y^T into the matrix WYT.
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

__global__ void dbl_small_QWYT
 ( int dim, int szt, int coloff, double *Q, double *WYT, double *QWYT );
/*
 * DESCRIPTION :
 *   Multiplies Q with WYT into the matrix QWYT.
 *
 * ON ENTRY :
 *   dim      number of rows and columns of all matrices;
 *   szt      the number of threads in a block;
 *   coloff   offset for the column index in QWYT;
 *   Q        the current orthogonal matrix;
 *   WYT      the product of W with Y^T;
 *   QWYT     space for a dim-by-dim matrix.
 *
 * ON RETURN :
 *   QWYT     the product of Q with QWYT. */

__global__ void dbl_small_YWTC
 ( int nrows, int ncols, int rowdim, int coldim, int szt,
   int rowoff, int coloff, double *YWT, double *C, double *YWTC );
/*
 * DESCRIPTION :
 *   Multiplies YWT with C into the matrix YWTC.
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

__global__ void dbl_small_Qupdate
 ( int dim, int szt, int coloff, double *Q, double *QWYT );
/*
 * DESCRIPTION :
 *   Updates Q by adding QWYT.
 *
 * ON ENTRY :
 *   dim      number of rows and columns of all matrices;
 *   szt      the number of threads in a block;
 *   coloff   offset for the column index in QWYT;
 *   Q        the current orthogonal matrix;
 *   QWYT     space for a dim-by-dim matrix.
 *
 * ON RETURN :
 *   Q        Q + QWYT. */

void GPU_dbl_small_house
 ( int nrows, int ncols, int szt, int nbt,
   int colidx, int nrows1, int k, int L,
   double *x0_d, double *A_h, double *A_d,
   double *v_h, double *V_d, double *beta_h, double *beta_d,
   double *lapms, bool verbose=true );
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
 *   x0_d     space on the device for the first element of the vector
 *            in the current column;
 *   A_h      matrix on the host;
 *   A_d      matrix on the device;
 *   v_h      space for the current Householder vector;
 *   V_d      space allocated for the Householder vectors on the device;
 *   beta_h   space allocated for the betas if verbose;
 *   beta_d   space allocated on the device for the betas;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   v_h      the next Householder vector on the host, if verbose;
 *   V_d      contains the next computed Householder vector;
 *   beta_h   updated vector of betas, if verbose;
 *   beta_d   the next beta constant;
 *   lapms    elapsed time spent by the kernel. */

void GPU_dbl_small_leftRupdate
 ( int nrows, int ncols, int szt, int colidx, int L,
   double *A_h, double *A_d, double *V_d, double *beta_h, double *beta_d,
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
 *   L        local index of the column in the current tile;
 *   A_h      matrix on the host;
 *   A_d      matrix on the device;
 *   V_d      space allocated for the Householder vectors on the device;
 *   beta_h   space allocated for the betas if verbose;
 *   beta_d   space allocated on the device for the betas;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   V_d      contains the Householder vectors on the device;
 *   beta_h   vector of betas, if verbose;
 *   beta_d   the next beta constant;
 *   lapms    elapsed time spent by the kernel. */

void GPU_dbl_VB_to_W
 ( int nrows, int ncols, int szt,
   double *V_h, double *V_d, double *W_h, double *W_d,
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
 *   V_h      space allocated for the Householder vectors, if verbose;
 *   V_d      space allocated for V on the device;
 *   W_h      space allocated for the W matrix, if verbose;
 *   W_d      space allocated for W on the device;
 *   beta_h   space allocated for the betas if verbose;
 *   beta_d   space allocated on the device for the betas;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   V_h      equals the Y matrix, if verbose;
 *   V_d      the Y matrix on the device;
 *   W_d      the W matrix in the WY representation, on the device;
 *   W_h      the W matrix in the WY representation, if verbose;
 *   lapms    elapsed time spent by the kernel. */

void GPU_dbl_small_WYT
 ( int nrows, int szt, double *W_d, double *Y_d, double *WYT_d,
   double *WYT_h, double *lapms, bool verbose=true );
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
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   WYT_d    the product W*Y^T on the device;
 *   WYT_h    the product W*Y^T, if verbose;
 *   lapms    elapsed time spent by the kernel. */

void GPU_dbl_small_YWT
 ( int nrows, int szt, double *Y_d, double *W_d, double *YWT_d,
   double *YWT_h, double *lapms, bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernel to compute W*Y^T.
 *   Wraps the timer and the print statements if verbose.
 *   If verbose, then the W*Y^T matrix is returned.
 *
 * ON ENTRY :
 *   nrows    number of rows in all matrices;
 *   szt      size of one tile and the number of threads in a block;
 *   Y_d      the matrix of Householder vectors;
 *   W_d      the matrix W in the WY representation;
 *   YWT_d    space for an nrows-by-nrows matrix on the device;
 *   YWT_h    space for an nrows-by-nrows matrix on the host, if verbose;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   YWT_d    the product Y*W^T on the device;
 *   YWT_h    the product Y*W^T, if verbose;
 *   lapms    elapsed time spent by the kernel. */

void GPU_dbl_small_QWYT
 ( int dim, int szt, int idx, double *Q_d, double *WYT_d, double *QWYT_d,
   double *QWYT_h, double *lapms, bool verbose=true );
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
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   QWYT_d   the product Q*WYT on the device;
 *   QWYT_h   the product Q*WYT, if verbose;
 *   lapms    elapsed time spent by the kernel. */

void GPU_dbl_small_YWTC
 ( int nrows, int ncols, int szt, int idx,
   double *YWT_d, double *C_d, double *YWTC_d, double *YWTC_h,
   double *lapms, bool verbose );
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
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   YWTC_d   the product YWT*C on the device;
 *   YWTC_h   the product YWT*C, if verbose;
 *   lapms    elapsed time spent by the kernel. */

void GPU_dbl_small_Qupdate
 ( int dim, int szt, int idx, double *Q_d, double *QWYT_d,
   double *Q_h, double *lapms, bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernel to update Q as Q + QWYT.
 *   Wraps the timer and the print statement if verbose.
 *   If verbose, then the update Q is returned.
 *
 * ON ENTRY :
 *   dim      number of rows and column in Q and QWYT;
 *   szt      size of one tile and the number of threads in a block;
 *   idx      index of the current tile;
 *   Q_d      a dim-by-dim matrix, on the device;
 *   QWYT_d   the product Q*W*Y^T, on the device;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Q_d      the updated Q on the device;
 *   Q_h      the updated Q, if verbose;
 *   lapms    elapsed time spent by the kernel. */

void GPU_dbl_blocked_houseqr
 ( int nrows, int ncols, int szt, int nbt,
   double **A, double **Q, double **R,
   double *houselapms, double *tileRlapms, double *vb2Wlapms,
   double *WYTlapms, double *QWYTlapms, double *Qaddlapms,
   double *YWTlapms, double *YWTClapms, 
   double *walltimesec, bool verbose=true );
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
 *   walltimesec is the elapsed wall clock computation time. */

#endif
