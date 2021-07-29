/* The file dbl2_baqr_kernels.h specifies functions for the
 * blocked accelerated QR in double double precision. */

#ifndef __dbl2_baqr_kernels_h__
#define __dbl2_baqr_kernels_h__

#define d_shmemsize 1024

__global__ void dbl2_small_house
 ( double *x0hi, double *x0lo, double *x1hi, double *x1lo,
   int dim, int dimLog2,
   double *vhi, double *vlo, double *betahi, double *betalo );
/*
 * DESCRIPTION :
 *   Computes the Householder vector of a vector of dimension dim+1,
 *   with one block of dim threads, on real data.
 *
 * ON ENTRY :
 *   x0hi     high double of the first element of the vector x; 
 *   x0lo     low double of the first element of the vector x; 
 *   x1hi     array of dim doubles, with the high doubles of x;
 *   x1lo     array of dim doubles, with the low doubles of x;
 *   dim      the dimension of the vector must equal the block size;
 *   dimLog2  equals ceil(log2((double) dim), used in sum reduction;
 *   vhi      space allocated for dim+1 doubles;
 *   vlo      space allocated for dim+1 doubles.
 *
 * ON RETURN :
 *   vhi      high doubles of the Householder vector;
 *   vlo      low doubles of the Householder vector;
 *   betahi   the high double of 2/(transpose(v)*v);
 *   betalo   the low double of 2/(transpose(v)*v). */

__global__ void dbl2_small_leftRupdate
 ( int nrows, int ncols, int szt, int k, double *Rhi, double *Rlo,
   double *vhi, double *vlo, double *betahi, double *betalo );
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
 *   Rhi      high doubles of an nrows-by-ncols matrix, stored column wise;
 *   Rlo      low doubles of an nrows-by-ncols matrix, stored column wise;
 *   vhi      high doubles of the Householder vector;
 *   vlo      low doubles of the Householder vector;
 *   betahi   the high double of 2/(transpose(v)*v);
 *   betalo   the low double of 2/(transpose(v)*v).
 *
 * ON RETURN :
 *   Rhi      high doubles of the updated matrix, which is trapezoidal;
 *   Rlo      low doubles of the updated matrix. */

__global__ void dbl2_VB_to_W
 ( int nrows, int ncols, double *Bhi, double *Blo,
   double *Vhi, double *Vlo, double *Whi, double *Wlo );
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
 *   Bhi      Bhi[i] is the high double of the i-th beta computed by house;
 *   Blo      Blo[i] is the low doubles of the i-th beta computed by house;
 *   Vhi      Vhi[nrows*i] is the start of the high doubles of the i-th
 *            Householder vector, with i zeros inserted so V is trapezoidal;
 *   Vlo      Vlo[nrows*i] is the start of the low doubles of the i-th
 *            Householder vector, with i zeros inserted so V is trapezoidal;
 *   Whi      space for ncols columns with rows from 0 to nrows-1;
 *   Wlo      space for ncols columns with rows from 0 to nrows-1.
 *
 * ON RETURN :
 *   Whi      high doubles of the W matrix in the WY representation;
 *   Wlo      low doubles of the W matrix in the WY representation. */

__global__ void dbl2_small_WYT
 ( int nrows, int szt, double *Whi, double *Wlo, double *Yhi, double *Ylo,
   double *WYThi, double *WYTlo );
/*
 * DESCRIPTION :
 *   Multiplies W with Y^T into the matrix WYT, on real data.
 *   Computes Y*W^T, swapping W with Y in the input arguments.
 *
 * ON ENTRY :
 *   nrows    number of rows of all matrices;
 *   szt      number of columns in W and Y,
 *            equals the number of threads in a block;
 *   Whi      high doubles of the W matrix in the WY representation;
 *   Wlo      low doubles of the W matrix in the WY representation;
 *   Yhi      the columns of Y are high doubles of the Householder vectors;
 *   Ylo      the columns of Y are low doubles of the Householder vectors;
 *   WYThi    space for an nrows-by-nrows matrix,
 *            with extra padding for when nrows > szt;
 *   WYTlo    space for an nrows-by-nrows matrix,
 *            with extra padding for when nrows > szt.
 *
 * ON RETURN :
 *   WYThi    high doubles of the product of W with Y^T;
 *   WYTlo    low doubles of the product of W with Y^T. */

__global__ void dbl2_small_QWYT
 ( int dim, int rowdim, int szt, int coloff,
   double *Qhi, double *Qlo, double *WYThi, double *WYTlo,
   double *QWYThi, double *QWYTlo );
/*
 * DESCRIPTION :
 *   Multiplies Q with WYT into the matrix QWYT, on real data.
 *
 * ON ENTRY :
 *   dim      number of rows and columns of the Q matrix;
 *   rowdim   number of rows and columns of the WYT matrix;
 *   szt      the number of threads in a block;
 *   coloff   offset for the column index in QWYT;
 *   Qhi      high doubles of the current orthogonal matrix;
 *   Qlo      low doubles of the current orthogonal matrix;
 *   WYThi    high doubles of the product of W with Y^T;
 *   WYTlo    low doubles of the product of W with Y^T;
 *   QWYThi   space for a dim-by-dim matrix;
 *   QWYTlo   space for a dim-by-dim matrix.
 *
 * ON RETURN :
 *   QWYThi   high doubles of the product of Q with QWYT;
 *   QWYTlo   low doubles of the product of Q with QWYT. */

__global__ void dbl2_small_YWTC
 ( int nrows, int ncols, int rowdim, int coldim, int szt,
   int rowoff, int coloff, double *YWThi, double *YWTlo,
   double *Chi, double *Clo, double *YWTChi, double *YWTClo );
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
 *   Chi      high doubles of the current matrix to be reduced;
 *   Clo      low doubles of the current matrix to be reduced;
 *   YWThi    high doubles of the product of Y with W^T;
 *   YWTlo    low doubles of the product of Y with W^T;
 *   YWTChi   space for a rowdim-by-coldim matrix;
 *   YWTClo   space for a rowdim-by-coldim matrix.
 *
 * ON RETURN :
 *   YWTChi   high doubles of the product of YWT with C;
 *   YWTClo   low doubles of the product of YWT with C. */

__global__ void dbl2_small_Qupdate
 ( int dim, int rowdim, int szt, int coloff,
   double *Qhi, double *Qlo, double *QWYThi, double *QWYTlo );
/*
 * DESCRIPTION :
 *   Updates Q by adding QWYT, on real data.
 *
 * ON ENTRY :
 *   dim      number of rows and columns of all matrices;
 *   szt      the number of threads in a block;
 *   coloff   offset for the column index in QWYT;
 *   Qhi      high doubles of the current orthogonal matrix;
 *   Qlo      low doubles of the current orthogonal matrix;
 *   QWYThi   space for a dim-by-dim matrix;
 *   QWYTlo   space for a dim-by-dim matrix.
 *
 * ON RETURN :
 *   Qhi      high doubles of Q + QWYT;
 *   Qlo      low doubles of Q + QWYT. */

__global__ void dbl2_small_R_add_YWTC
 ( int nrows, int coldim, int szt, int rowoff, int coloff,
   double *Rhi, double *Rlo, double *YWTChi, double *YWTClo );
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
 *   Rhi      high doubles of the current matrix to be reduced;
 *   Rlo      low doubles of the current matrix to be reduced;
 *   YWTChi   high doubles of the product of YWT with C;
 *   YWTClo   low doubles of the product of YWT with C.
 *
 * ON RETURN :
 *   Rhi      high doubles of R + YWTC;
 *   Rlo      low doubles of R + YWTC. */

void GPU_dbl2_small_house
 ( int nrows, int ncols, int szt, int nbt,
   int colidx, int nrows1, int k, int L,
   double *Ahi_h, double *Alo_h, double *Ahi_d, double *Alo_d,
   double *vhi_h, double *vlo_h, double *Vhi_d, double *Vlo_d,
   double *betahi_h, double *betalo_h, double *betahi_d, double *betalo_d,
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
 *   Ahi_h    high doubles of the matrix on the host;
 *   Alo_h    low doubles of the matrix on the host;
 *   Ahi_d    high doubles of the matrix on the device;
 *   Alo_d    low doubles of the matrix on the device;
 *   vhi_h    space for the current Householder vector;
 *   vlo_h    space for the current Householder vector;
 *   Vhi_d    space for the Householder vectors on the device;
 *   Vlo_d    space for the Householder vectors on the device;
 *   betahi_h has space for the high doubles of the betas if verbose;
 *   betalo_h has space for the low doubles of the betas if verbose;
 *   betahi_d has space on the device for the high doubles of the betas;
 *   betalo_d has space on the device for the low doubles of the betas;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   vhi_h    high doubles of the next Householder vector on the host,
 *            if verbose;
 *   vlo_h    low doubles of the next Householder vector on the host,
 *            if verbose;
 *   Vhi_d    high doubles of the next computed Householder vector;
 *   Vlo_d    low doubles of the next computed Householder vector;
 *   betahi_h has the updated vector of the high doubles of the betas,
 *            if verbose;
 *   betalo_h has the updated vector of low doubles of the betas,
 *            if verbose;
 *   betahi_d has the high double of the next beta constant;
 *   betalo_d has the low double of the next beta constant;
 *   lapms    elapsed time spent by the kernel. */

void GPU_dbl2_small_leftRupdate
 ( int nrows, int ncols, int szt, int colidx, int k, int L,
   double *Ahi_h, double *Alo_h, double *Ahi_d, double *Alo_d,
   double *Vhi_d, double *Vlo_d,
   double *betahi_h, double *betalo_h, double *betahi_d, double *betalo_d,
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
 *   Ahi_h    high doubles of the matrix on the host;
 *   Alo_h    low doubles of the matrix on the host;
 *   Ahi_d    high doubles of the matrix on the device;
 *   Alo_d    low doubles of the matrix on the device;
 *   Vhi_d    space allocated for the Householder vectors on the device;
 *   Vlo_d    space allocated for the Householder vectors on the device;
 *   betahi_h has space allocated for the betas if verbose;
 *   betalo_h has space allocated for the betas if verbose;
 *   betahi_d has space allocated on the device for the betas;
 *   betalo_d has space allocated on the device for the betas;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Ahi_d    high doubles of the reduced matrix on the device;
 *   Alo_d    low doubles of the reduced matrix on the device;
 *   Ahi_h    high doubles of the reduced matrix on the host, if verbose;
 *   Alo_h    low doubles of the reduced matrix on the host, if verbose;
 *   lapms    elapsed time spent by the kernel. */

void GPU_dbl2_VB_to_W
 ( int nrows, int ncols, int szt,
   double *Vhi_h, double *Vlo_h, double *Vhi_d, double *Vlo_d,
   double *Whi_h, double *Wlo_h, double *Whi_d, double *Wlo_d,
   double *betahi_h, double *betalo_h, double *betahi_d, double *betalo_d,
   double *lapms, bool verbose=true );
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
 *   Vhi_h    space for the Householder vectors, if verbose;
 *   Vlo_h    space for the Householder vectors, if verbose;
 *   Vhi_d    space for V on the device;
 *   Vlo_d    space for V on the device;
 *   Whi_h    space for the W matrix, if verbose;
 *   Wlo_h    space for the W matrix, if verbose;
 *   Whi_d    space for W on the device;
 *   Wlo_d    space for W on the device;
 *   betahi_h has space for the betas if verbose;
 *   betalo_h has space for the betas if verbose;
 *   betahi_d has space on the device for the betas;
 *   betalo_d has space on the device for the betas;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Vhi_h    equals the high doubles of the Y matrix, if verbose;
 *   Vlo_h    equals the low doubles of the Y matrix, if verbose;
 *   Vhi_d    high doubles of the Y matrix on the device;
 *   Vlo_d    low doubles of the Y matrix on the device;
 *   Whi_d    high doubles of the W matrix in the WY representation,
 *            on the device;
 *   Wlo_d    low doubles of the W matrix in the WY representation,
 *            on the device;
 *   Whi_h    high doubles of the W matrix in the WY representation,
 *            if verbose;
 *   Wlo_h    low doubles of the W matrix in the WY representation,
 *            if verbose;
 *   lapms    elapsed time spent by the kernel. */

void GPU_dbl2_small_WYT
 ( int nrows, int szt,
   double *Whi_d, double *Wlo_d, double *Yhi_d, double *Ylo_d,
   double *WYThi_d, double *WYTlo_d, double *WYThi_h, double *WYTlo_h,
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
 *   Whi_d    high doubles of the matrix W in the WY representation;
 *   Wlo_d    low doubles of the matrix W in the WY representation;
 *   Yhi_d    high doubles of the matrix of Householder vectors;
 *   Ylo_d    low doubles of the matrix of Householder vectors;
 *   WYThi_d  space for an nrows-by-nrows matrix on the device;
 *   WYTlo_d  space for an nrows-by-nrows matrix on the device;
 *   WYThi_h  space for an nrows-by-nrows matrix on the host, if verbose;
 *   WYTlo_h  space for an nrows-by-nrows matrix on the host, if verbose;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   WYThi_d  high doubles of the product W*Y^T on the device;
 *   WYTlo_d  low doubles of the product W*Y^T on the device;
 *   WYThi_h  high doubles of the product W*Y^T, if verbose;
 *   WYTlo_h  low doubles of the product W*Y^T, if verbose;
 *   lapms    elapsed time spent by the kernel. */

void GPU_dbl2_small_YWT
 ( int nrows, int szt, int idx,
   double *Yhi_d, double *Ylo_d, double *Whi_d, double *Wlo_d,
   double *YWThi_d, double *YWTlo_d, double *YWThi_h, double *YWTlo_h,
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
 *   Yhi_d    high doubles of the matrix of Householder vectors;
 *   Ylo_d    low doubles of the matrix of Householder vectors;
 *   Whi_d    high doubles of the matrix W in the WY representation;
 *   Wlo_d    low doubles of the matrix W in the WY representation;
 *   YWThi_d  space for an nrows-by-nrows matrix on the device;
 *   YWTlo_d  space for an nrows-by-nrows matrix on the device;
 *   YWThi_h  space for an nrows-by-nrows matrix on the host, if verbose;
 *   YWTlo_h  space for an nrows-by-nrows matrix on the host, if verbose;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   YWThi_d  high doubles of the product Y*W^T on the device;
 *   YWTlo_d  low doubles of the product Y*W^T on the device;
 *   YWThi_h  high doubles of the product Y*W^T, if verbose;
 *   YWTlo_h  low doubles of the product Y*W^T, if verbose;
 *   lapms    elapsed time spent by the kernel. */

void GPU_dbl2_small_QWYT
 ( int dim, int szt, int idx, double *Qhi_d, double *Qlo_d,
   double *WYThi_d, double *WYTlo_d, double *QWYThi_d, double *QWYTlo_d,
   double *QWYThi_h, double *QWYTlo_h, double *Qhi_h, double *Qlo_h,
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
 *   Qhi_d    high doubles of a dim-by-dim matrix, on the device;
 *   Qlo_d    low doubles of a dim-by-dim matrix, on the device;
 *   WYThi_d  high doubles of the product W*Y^T, on the device;
 *   WYTlo_d  low doubles of the product W*Y^T, on the device;
 *   QWYThi_d has space for the product Q*WYT, on the device;
 *   QWYTlo_d has space for the product Q*WYT, on the device;
 *   QWYThi_h has space for the product Q*WYT, on the host, if verbose;
 *   QWYTlo_h has space for the product Q*WYT, on the host, if verbose;
 *   Qhi_h    if verbose, then used to print Q before the product;
 *   Qlo_h    if verbose, then used to print Q before the product;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   QWYThi_d has the high doubles of the product Q*WYT on the device;
 *   QWYTlo_d has the low doubles of the product Q*WYT on the device;
 *   QWYThi_h has the high doubles of the product Q*WYT, if verbose;
 *   QWYTlo_h has the low doubles of the product Q*WYT, if verbose;
 *   lapms    elapsed time spent by the kernel. */

void GPU_dbl2_small_YWTC
 ( int nrows, int ncols, int szt, int idx, double *YWThi_d, double *YWTlo_d,
   double *Chi_d, double *Clo_d, double *YWTChi_d, double *YWTClo_d,
   double *YWTChi_h, double *YWTClo_h, double *lapms, bool verbose=true );
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
 *   YWThi_d  high doubles of the product Y*W^T, on the device;
 *   YWTlo_d  low doubles of the product Y*W^T, on the device;
 *   Chi_d    high doubles of an nrows-by-ncols matrix, on the device;
 *   Clo_d    low doubles of an nrows-by-ncols matrix, on the device;
 *   YWTChi_d has space for the product YWT*C, on the device;
 *   YWTClo_d has space for the product YWT*C, on the device;
 *   YWTChi_h has space for the product YWT*C, on the host, if verbose;
 *   YWTClo_h has space for the product YWT*C, on the host, if verbose;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   YWTChi_d has the high doubles of the product YWT*C on the device;
 *   YWTClo_d has the low doubles of the product YWT*C on the device;
 *   YWTChi_h has the high doubles of the product YWT*C, if verbose;
 *   YWTClo_h has the low doubles of the product YWT*C, if verbose;
 *   lapms    elapsed time spent by the kernel. */

void GPU_dbl2_small_Qupdate
 ( int dim, int szt, int idx,
   double *Qhi_d, double *Qlo_d, double *QWYThi_d, double *QWYTlo_d,
   double *Qhi_h, double *Qlo_h, double *lapms, bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernel to update Q as Q + QWYT.
 *   Wraps the timer and the print statement if verbose.
 *   If verbose, then the updated Q is returned.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in Q,
 *            number of rows in QWYT;
 *   szt      size of one tile and the number of threads in a block;
 *   idx      index of the current tile;
 *   Qhi_d    a dim-by-dim matrix, on the device;
 *   Qlo_d    a dim-by-dim matrix, on the device;
 *   QWYThi_d has the high doubles of the product Q*W*Y^T, on the device;
 *   QWYTlo_d has the low doubles of the product Q*W*Y^T, on the device;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Qhi_d    high doubles of the updated Q on the device;
 *   Qlo_d    low doubles of the updated Q on the device;
 *   Qhi_h    high doubles of the updated Q, if verbose;
 *   Qlo_h    low doubles of the updated Q, if verbose;
 *   lapms    elapsed time spent by the kernel. */

void GPU_dbl2_small_R_add_YWTC
 ( int nrows, int ncols, int szt, int idx, double *Rhi_d, double *Rlo_d,
   double *YWTChi_d, double *YWTClo_d, double *Rhi_h, double *Rlo_h,
   double *lapms, bool verbose=true );
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
 *   Rhi_d    an nrows-by-ncols matrix, on the device;
 *   Rlo_d    an nrows-by-ncols matrix, on the device;
 *   YWTChi_d has the high doubles of the product Y*W^T*C, on the device;
 *   YWTClo_d has the low doubles of the product Y*W^T*C, on the device;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Rhi_d    high doubles of the updated R on the device;
 *   Rlo_d    low doubles of the updated R on the device;
 *   Rhi_h    high doubles of the updated R, if verbose;
 *   Rlo_h    low doubles of the updated R, if verbose;
 *   lapms    elapsed time spent by the kernel. */

void GPU_dbl2_blocked_houseqr
 ( int nrows, int ncols, int szt, int nbt,
   double **Ahi, double **Alo, double **Qhi, double **Qlo,
   double **Rhi, double **Rlo,
   double *houselapms, double *tileRlapms, double *vb2Wlapms,
   double *WYTlapms, double *QWYTlapms, double *Qaddlapms,
   double *YWTlapms, double *YWTClapms, double *Raddlapms,
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
 *   Ahi      high doubles of an nrows-by-ncols matrix,
 *            stored as nrows arrays of ncols numbers;
 *   Alo      low doubles of an nrows-by-ncols matrix,
 *            stored as nrows arrays of ncols numbers;
 *   Qhi      space for an nrows-by-nrows matrix;
 *   Qlo      space for an nrows-by-nrows matrix;
 *   Rhi      space for an nrows-by-ncols matrix;
 *   Rlo      space for an nrows-by-ncols matrix;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Qhi      high doubles of an orthogonal matrix, transpose(Q)*A = R;
 *   Qlo      low doubles of an orthogonal matrix, transpose(Q)*A = R;
 *   Rhi      high doubles of the reduced upper triangular form of A;
 *   Rlo      low doubles of the reduced upper triangular form of A;
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
