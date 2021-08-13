/* The file dbl_baqr_kernels.h specifies functions for the
 * blocked accelerated QR in double precision. */

#ifndef __dbl_baqr_kernels_h__
#define __dbl_baqr_kernels_h__

#define d_shmemsize 2048
#define cd_shmemsize 1024
// for complex problems, the d_shmemsize is divided in half

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

__global__ void cmplx_small_betaRHv
 ( int nrows, int ncols, int szt, int k,
   double *Rre, double *Rim, double *vre, double *vim, double *beta,
   double *wre, double *wim );
/*
 * DESCRIPTION :
 *   Computes the vector w = beta*R^H*v, on complex data,
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

__global__ void dbl_RTdotv
 ( int nrows, int szt, int colidx, int Roffset, int dim,
   double *R, double *v, double *RTdotv );
/*
 * DESCRIPTION :
 *   The elements of the matrix RTdotv are the elements of R^T,
 *   multiplied with the corresponding element of v.
 *   Multiple blocks of threads operate on real data.
 *
 * ON ENTRY :
 *   nrows    number of rows of R;
 *   szt      size of one tile and number of threads in one block;
 *   colidx   index of the current column in R;
 *   Roffset  offset in R for the first row to start;
 *   dim      number of columns in R;
 *   R        column wise stored matrix with number of rows equal to nrows;
 *   v        start of the first nonzero element of a Householder vector;
 *   RTdotv   space for a matrix of nrows-by-szt, plus some padding.
 *
 * ON RETURN :
 *   RTdotv   the element-by-element products of R^T with v,
 *            stored row by row. */

__global__ void cmplx_RHdotv
 ( int nrows, int szt, int colidx, int Roffset, int dim,
   double *Rre, double *Rim, double *vre, double *vim,
   double *RHdotvre, double *RHdotvim );
/*
 * DESCRIPTION :
 *   The elements of the matrix RHdotv are the elements of R^H,
 *   multiplied with the corresponding element of v.
 *   Multiple blocks of threads operate on complex data.
 *
 * ON ENTRY :
 *   nrows    number of rows of R;
 *   szt      size of one tile and number of threads in one block;
 *   colidx   index of the current column in R;
 *   Roffset  offset in R for the first row to start;
 *   dim      number of columns in R;
 *   Rre      real parts of a column wise stored matrix
 *            with number of rows equal to nrows;
 *   Rim      imaginary parts of a column wise stored matrix
 *            with number of rows equal to nrows;
 *   vre      start of the real parts of the first nonzero element
 *            of a Householder vector;
 *   vim      start of the imaginary parts of the first nonzero element
 *            of a Householder vector;
 *   RTdotvre has space for a matrix of nrows-by-szt, plus some padding;
 *   RTdotvim has space for a matrix of nrows-by-szt, plus some padding.
 *
 * ON RETURN :
 *   RTdotvre has the real parts of the element-by-element products
 *            of R^T with v, stored row by row;
 *   RTdotvim has the imaginary parts of the element-by-element products
 *            of R^T with v, stored row by row. */

__global__ void dbl_sum_betaRTdotv
 ( int nrows, double *beta, double *RTdotv, double *w );
/*
 * DESCRIPTION :
 *   Adds the rows in RTdotv to obtain w = beta*R^T*v,
 *   with one block of threads, on real data.
 *
 * ON ENTRY :
 *   nrows    number of rows in RTdotv;
 *   beta     the beta value corresponding to the Householder vector;
 *   RTdotv   contains all products of the elements of R dotted v;
 *   w        space for beta*R^T*v.
 *
 * ON RETURN :
 *   w        contains beta*R^T*v. */

__global__ void cmplx_sum_betaRHdotv
 ( int nrows, double *beta, double *RHdotvre, double *RHdotvim,
   double *wre, double *wim );
/*
 * DESCRIPTION :
 *   Adds the rows in RHdotv to obtain w = beta*R^H*v,
 *   with one block of threads, on complex data.
 *
 * ON ENTRY :
 *   nrows    number of rows in RHdotv;
 *   beta     the beta value corresponding to the Householder vector;
 *   RTdotv   contains all products of the elements of R dotted v;
 *   wre      space for the real parts of beta*R^H*v.
 *   wim      space for the imaginary parts of beta*R^H*v.
 *
 * ON RETURN :
 *   wre      the real parts of beta*R^H*v;
 *   wim      the imaginary parts of beta*R^H*v. */

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

__global__ void cmplx_medium_subvbetaRHv
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

__global__ void dbl_beta_times_V
 ( int nrows, int szt, double *B, double *V, double *W );
/*
 * DESCRIPTION :
 *   Computes the first vector in the W representation of the Householder
 *   transformations, multiplying B[0] with the first vector of V,
 *   and flipping the sign, with multiple blocks of threads, on real data.
 *
 * ON ENTRY :
 *   nrows    number of rows in V and W;
 *   szt      size of one tile and the number of threads in a block;
 *   B        B[i] is the i-th beta computed by house;
 *   V        V[nrows*i] is the start of the i-th Householder vector,
 *            with i zeros inserted so V is trapezoidal;
 *   W        space for the W matrix.
 *
 * ON RETURN :
 *   W        the first nrows numbers store the first vector
 *            of the W matrix in the WY representation. */

__global__ void cmplx_beta_times_V
 ( int nrows, int szt, double *B,
   double *Vre, double *Vim, double *Wre, double *Wim );
/*
 * DESCRIPTION :
 *   Computes the first vector in the W representation of the Householder
 *   transformations, multiplying B[0] with the first vector of V,
 *   and flipping the sign, with multiple blocks of threads, on real data.
 *
 * ON ENTRY :
 *   nrows    number of rows in V and W;
 *   szt      size of one tile and the number of threads in a block;
 *   B        B[i] is the i-th beta computed by house;
 *   Vre      Vre[nrows*i] is the start of the real parts of the i-th
 *            Householder vector, with i zeros inserted so V is trapezoidal;
 *   Vim      Vim[nrows*i] is the start of the imaginary parts of the i-th
 *            Householder vector, with i zeros inserted so V is trapezoidal;
 *   Wre      space for the real parts of the W matrix;
 *   Wim      space for the imaginary parts of the W matrix.
 *
 * ON RETURN :
 *   Wre      the first nrows numbers store the real parts of the first
 *            vector of the W matrix in the WY representation;
 *   Wim      the first nrows numbers store the imaginary parts of the first
 *            vector of the W matrix in the WY representation. */

__global__ void dbl_initialize_WYT
 ( int dim, int szt, double *V, double *W, double *WYT );
/*
 * DESCRIPTION :
 *   Initializes the matrix YWT with the product of the first dim products
 *   of W with Y elements with multiple blocks of szt threads,
 *   on real data.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in YWT;
 *   szt      number of threads in one block;
 *   V        the first dim numbers define the first Householder vector;
 *   W        the first dim numbers define the first column in the W matrix;
 *   WYT      space for a dim-by-dim matrix.
 *
 * ON RETURN :
 *   WYT      equals w*y^T, where y and w are the first columns of V and W. */

__global__ void cmplx_initialize_WYH
 ( int dim, int szt, double *Vre, double *Vim, double *Wre, double *Wim,
   double *YWHre, double *YWHim );
/*
 * DESCRIPTION :
 *   Initializes the matrix YWT with the product of the first dim products
 *   of W with Y elements with multiple blocks of szt threads,
 *   on complex data.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in YWT;
 *   szt      number of threads in one block;
 *   Vre      the first dim numbers define the real parts
 *            of the first Householder vector;
 *   Vim      the first dim numbers define the imaginary parts
 *            of the first Householder vector;
 *   Wre      the first dim numbers define the real parts
 *            of the first column in the W matrix;
 *   Wim      the first dim numbers define the real parts
 *            of the first column in the W matrix;
 *   WYHre    space for a dim-by-dim matrix;
 *   WYHim    space for a dim-by-dim matrix.
 *
 * ON RETURN :
 *   WYHre    equals the real parts of w*y^H,
 *            where y and w are the first columns of V and W;
 *   WYHim    equals the imaginary parts of y*w^H,
 *            where y and w are the first columns of V and W. */

__global__ void dbl_update_WYT
 ( int dim, int szt, double *V, double *W, double *WYT );
/*
 * DESCRIPTION :
 *   Updates the matrix YWT with the product of the first dim products
 *   of W with Y elements with multiple blocks of szt threads,
 *   on real data.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in YWT;
 *   szt      number of threads in one block;
 *   V        the first dim numbers define the next Householder vector;
 *   W        the first dim numbers define the next column in the W matrix;
 *   WYT      an dim-by-dim matrix with the current values for WYT.
 *
 * ON RETURN :
 *   WYT      the updated matrix W*Y^T. */

__global__ void cmplx_update_WYH
 ( int dim, int szt, double *Vre, double *Vim, double *Wre, double *Wim,
   double *WYHre, double *WYHim );
/*
 * DESCRIPTION :
 *   Updates the matrix YWT with the product of the first dim products
 *   of W with Y elements with multiple blocks of szt threads,
 *   on complex data.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in YWT;
 *   szt      number of threads in one block;
 *   Vre      the first dim numbers define the real parts
 *            of the first Householder vector;
 *   Vim      the first dim numbers define the imaginary parts
 *            of the first Householder vector;
 *   Wre      the first dim numbers define the real parts
 *            of the first column in the W matrix;
 *   Wim      the first dim numbers define the real parts
 *            of the first column in the W matrix;
 *   WYHre    real parts of the current values for WYH;
 *   WYHim    imaginary parts of the current values for WYH.
 *
 * ON RETURN :
 *   WYHre    the real parts of the updated W*Y*H;
 *   WYHim    the imaginary parts of the update W*Y^H. */

__global__ void dbl_beta_next_W
 ( int nrows, int szt, double *B, double *V, double *W, double *WYT );
/*
 * DECRIPTION :
 *   Computes the next column in the W matrix, with multiple blocks,
 *   on real data.
 *
 * ON ENTRY :
 *   nrows    number of rows in V, W, and the dimension of WYT;
 *   szt      number of threads in one block;
 *   V        the first nrows numbers define the next Householder vector;
 *   W        the first nrows numbers define the next column in the W matrix;
 *   WYT      an nrows-by-nrows matrix with the current values for WYT.
 *
 * ON RETURN :
 *   W        contains the values of the next column of W. */

__global__ void cmplx_beta_next_W
 ( int nrows, int szt, double *B, double *Vre, double *Vim,
   double *Wre, double *Wim, double *WYHre, double *WYHim );
/*
 * DECRIPTION :
 *   Computes the next column in the W matrix, with multiple blocks,
 *   on complex data.
 *
 * ON ENTRY :
 *   nrows    number of rows in V, W, and the dimension of WYH;
 *   szt      number of threads in one block;
 *   Vre      the first nrows numbers define the real parts
 *            of the first Householder vector;
 *   Vim      the first nrows numbers define the imaginary parts
 *            of the first Householder vector;
 *   Wre      the first nrows numbers define the real parts
 *            of the first column in the W matrix;
 *   Wim      the first nrows numbers define the real parts
 *            of the first column in the W matrix;
 *   WYHre    real parts of WYH;
 *   WYHim    imaginary parts of WYH.
 *
 * ON RETURN :
 *   Wre      contains the real parts of the next column of W;
 *   Wim      contains the imaginary parts of the next column of W. */

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

__global__ void cmplx_small_WYH
 ( int nrows, int szt, double *Wre, double *Wim,
   double *Yre, double *Yim, double *WYHre, double *WYHim );
/*
 * DESCRIPTION :
 *   Multiplies W with Y^H into the matrix WYH, on complex data.
 *   Computes Y*W^H, swapping W with Y in the input arguments.
 *
 * ON ENTRY :
 *   nrows    number of rows of all matrices;
 *   szt      number of columns in W and Y,
 *            equals the number of threads in a block;
 *   Wre      real parts of the W matrix in the WY representation;
 *   Wim      imaginary parts of the W matrix in the WY representation;
 *   Yre      real parts of the columns of Y are the Householder vectors;
 *   Yim      imaginary parts of the columns of Y are the Householder vectors;
 *   WYHre    space for an nrows-by-nrows matrix,
 *            with extra padding for when nrows > szt;
 *   WYHim    space for an nrows-by-nrows matrix,
 *            with extra padding for when nrows > szt.
 *
 * ON RETURN :
 *   WYHre    real parts of the product of W with Y^H;
 *   WYHim    imaginary parts of the product of W with Y^H. */

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

__global__ void cmplx_small_QWYH
 ( int dim, int rowdim, int szt, int coloff,
   double *Qre, double *Qim, double *WYHre, double *WYHim,
   double *QWYHre, double *QWYHim );
/*
 * DESCRIPTION :
 *   Multiplies Q with WYH into the matrix QWYH, on complex data.
 *
 * ON ENTRY :
 *   dim      number of rows and columns of the Q matrix;
 *   rowdim   number of rows and columns of the WYH matrix;
 *   szt      the number of threads in a block;
 *   coloff   offset for the column index in QWYH;
 *   Qre      real parts of the current orthogonal matrix;
 *   Qim      imaginary parts of the current orthogonal matrix;
 *   WYHre    real parts of the product of W with Y^H;
 *   WYHim    imaginary parts of the product of W with Y^H;
 *   QWYHre   space for a dim-by-dim matrix;
 *   QWYHim   space for a dim-by-dim matrix.
 *
 * ON RETURN :
 *   QWYHre   real parts of the product of Q with QWYH;
 *   QWYHim   imaginary parts of the product of Q with QWYH. */

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

__global__ void cmplx_small_YWHC
 ( int nrows, int ncols, int rowdim, int coldim, int szt,
   int rowoff, int coloff, double *YWHre, double *YWHim,
   double *Cre, double *Cim, double *YWHCre, double *YWHCim );
/*
 * DESCRIPTION :
 *   Multiplies YWH with C into the matrix YWHC, on complex data.
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
 *   YWH      real parts of the product of Y with W^H;
 *   YWH      imaginary parts of the product of Y with W^H;
 *   YWHCre   has space for a rowdim-by-coldim matrix;
 *   YWHCim   has space for a rowdim-by-coldim matrix.
 *
 * ON RETURN :
 *   YWHCre   real parts of the product of YWH with C;
 *   YWHCim   imaginary parts of the product of YWH with C. */

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

__global__ void cmplx_small_Qupdate
 ( int dim, int rowdim, int szt, int coloff,
   double *Qre, double *Qim, double *QWYHre, double *QWYHim );
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
 *   QWYHre   space for a dim-by-dim matrix;
 *   QWYHim   space for a dim-by-dim matrix.
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

__global__ void cmplx_small_R_add_YWHC
 ( int nrows, int coldim, int szt, int rowoff, int coloff,
   double *Rre, double *Rim, double *YWHCre, double *YWHCim );
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
 *   YWHCre   real parts of the product of YWH with C.
 *   YWHCim   imaginary parts of the product of YWH with C.
 *
 * ON RETURN :
 *   Rre      real parts of R + YWHC;
 *   Rim      imaginary parts R + YWHC. */

void GPU_dbl_small_house
 ( int nrows, int ncols, int szt, int nbt,
   int colidx, int nrows1, int k, int L,
   double *A_h, double *A_d,
   double *v_h, double *V_d, double *beta_h, double *beta_d,
   double *lapms, long long int *add, long long int *mul, long long int *div,
   long long int *sqrtfun, bool verbose=true );
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
 *   sqrtfun  current number of calls to sqrt();
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
 *   div      accumulated number of divisions;
 *   sqrtfun  accumulated number of calls to sqrt(). */

void GPU_cmplx_small_house
 ( int nrows, int ncols, int szt, int nbt,
   int colidx, int nrows1, int k, int L,
   double *Are_h, double *Aim_h, double *Are_d, double *Aim_d,
   double *vre_h, double *vim_h, double *Vre_d, double *Vim_d,
   double *beta_h, double *beta_d,
   double *lapms, long long int *add, long int *mul, long long int *div,
   long long int *sqrtfun, bool verbose=true );
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
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications;
 *   div      current number of divisions;
 *   sqrtfun  current number of calls to sqrt();
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   v_h      the next Householder vector on the host, if verbose;
 *   Vre_d    real parts of the next computed Householder vector;
 *   Vim_d    imaginary parts of the next computed Householder vector;
 *   beta_h   updated vector of betas, if verbose;
 *   beta_d   the next beta constant;
 *   lapms    elapsed time spent by the kernel;
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications;
 *   div      accumulated number of divisions;
 *   sqrtfun  accumulated number of calls to sqrt(). */

void GPU_dbl_small_leftRupdate
 ( int nrows, int ncols, int szt, int colidx, int k, int L,
   double *A_h, double *A_d, double *V_d, double *beta_h, double *beta_d,
   double *lapms, long long int *add, long long int *mul, bool verbose=true );
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
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   V_d      contains the Householder vectors on the device;
 *   beta_h   vector of betas, if verbose;
 *   beta_d   the next beta constant;
 *   lapms    elapsed time spent by the kernel.
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications. */

void GPU_cmplx_small_leftRupdate
 ( int nrows, int ncols, int szt, int colidx, int k, int L,
   double *Are_h, double *Aim_h, double *Are_d, double *Aim_d,
   double *Vre_d, double *Vim_d, double *beta_h, double *beta_d,
   double *lapms, long long int *add, long long int *mul, bool verbose=true );
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
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Vre_d    real parts of the Householder vectors on the device;
 *   Vim_d    imaginary parts of the Householder vectors on the device;
 *   beta_h   vector of betas, if verbose;
 *   beta_d   the next beta constant;
 *   lapms    elapsed time spent by the kernel;
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications;
 *   div      accumulated number of divisions. */

void GPU_dbl_medium_leftRupdate
 ( int nrows, int ncols, int szt, int colidx, int k, int L,
   double *A_h, double *A_d, double *V_d, double *beta_h, double *beta_d,
   double *RTdotv_h, double *RTdotv_d, double *w_h, double *w_d,
   double *RTvlapms, double *redlapms,
   long long int *add, long long int *mul, bool verbose=true );
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
 *   RTdotv_h has space for the componentwise product of R^T with v,
 *            if verbose;
 *   RTdotv_d has space for the componentwise product of R^T with v,
 *            on the device;
 *   w_h      space for the beta*R^T*v on the host, if verbose;
 *   w_d      space for the beta*R^T*v plus szt padding;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   V_d      contains the Householder vectors on the device;
 *   beta_h   vector of betas, if verbose;
 *   beta_d   the next beta constant;
 *   RTdotv_h stores the componentwise product of R^T with v, if verbose;
 *   RTdotv_d stores R^T*v, on the device;
 *   w_h      stores the beta*R^T*v, if verbose;
 *   w_d      stores the beta*R^T*v, on the device;
 *   RTvlapms is the elapsed time spent to compute beta*R^T*v;
 *   redlapms is the elapsed time spent to reduce one tile;
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications. */

void GPU_cmplx_medium_leftRupdate
 ( int nrows, int ncols, int szt, int colidx, int k, int L,
   double *Are_h, double *Aim_h, double *Are_d, double *Aim_d,
   double *Vre_d, double *Vim_d, double *beta_h, double *beta_d,
   double *RHdotvre_h, double *RHdotvim_h,
   double *RHdotvre_d, double *RHdotvim_d,
   double *wre_h, double *wim_h, double *wre_d, double *wim_d,
   double *RHvlapms, double *redlapms,
   long long int *add, long long int *mul, bool verbose=true );
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
 *   RHdotvre_h has space for the real parts of the componentwise products
 *            of R^H with v, if verbose;
 *   RHdotvim_h has space for the imaginary parts of the componentwise
 *            products of R^H with v, if verbose;
 *   RHdotvre_d has space for RHdotvre, on the device;
 *   RHdotvim_d has space for RHdotvim, on the device;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Vre_d    real parts of the Householder vectors on the device;
 *   Vim_d    imaginary parts of the Householder vectors on the device;
 *   beta_h   vector of betas, if verbose;
 *   beta_d   the next beta constant;
 *   RHdotvre_h stores the real parts of the componentwise products
 *            of R^H with v, if verbose;
 *   RHdotvim_h stores the imaginary parts of the componentwise
 *            products of R^H with v, if verbose;
 *   RHdotvre_d stores RHdotvre, on the device;
 *   RHdotvim_d stores RHdotvim, on the device;
 *   RHvlapms elapsed time spent on beta*R^H*v;
 *   redlapms elapsed time spent to reduce one tile;
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications. */

void GPU_dbl_VB_to_W
 ( int nrows, int ncols, int szt,
   double *V_h, double *V_d, double *W_h, double *W_d,
   double *beta_h, double *beta_d, double *lapms,
   long long int *add, long long int *mul, long long int *div,
   bool verbose=true );
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
   double *beta_h, double *beta_d, double *lapms,
   long long int *add, long long int *mul, long long int *div,
   bool verbose=true );
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
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications;
 *   div      current number of divisions;
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
 *   lapms    elapsed time spent by the kernel;
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications;
 *   div      accumulated number of divisions. */

void GPU_dbl_medium_VB_to_W
 ( int nrows, int ncols, int szt, int idx,
   double *V_h, double *V_d, double *W_h, double *W_d,
   double *WYT_h, double *YWT_d, double *beta_h, double *beta_d,
   double *lapms, long long int *add, long long int *mul, bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernel to compute the W in the WY representation,
 *   on real data.
 *   Wraps the timer and the print statements if verbose.
 *   If verbose, then the W matrix is returned on the host.
 *
 * ON ENTRY :
 *   nrows    number of rows in the matrices V, Y, and W;
 *   ncols    equals the size of one tile, or equivalently,
 *            is the number of elements in B,
 *            and the number of columns in V, Y, and W;
 *   szt      size of one tile and the number of threads in a block;
 *   idx      index of the current tile;
 *   V_h      the Householder vectors, if verbose;
 *   V_d      the Householder vectors on the device;
 *   W_h      space for the W matrix, if verbose;
 *   W_d      space for W on the device;
 *   WYT_h    space for the outer product of W with Y^T, if verbose;
 *   WYT_d    space for the outer product of W with Y^T, on the device;
 *   beta_h   space for the betas if verbose;
 *   beta_d   space on the device for the betas;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   V_h      the Y matrix, if verbose;
 *   V_d      the Y matrix on the device;
 *   W_d      the W matrix in the WY representation, on the device;
 *   W_h      the W matrix in the WY representation, if verbose;
 *   WYT_h    the outer product of W with Y^T, if verbose;
 *   WYT_d    the outer product of W with Y^T, on the device;
 *   lapms    elapsed time spent by the kernel;
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications. */

void GPU_cmplx_medium_VB_to_W
 ( int nrows, int ncols, int szt, int idx,
   double *Vre_h, double *Vim_h, double *Vre_d, double *Vim_d,
   double *Wre_h, double *Wim_h, double *Wre_d, double *Wim_d,
   double *WYHre_h, double *WYHim_h, double *WYHre_d, double *WYHim_d,
   double *beta_h, double *beta_d,
   double *lapms, long long int *add, long long int *mul, bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernel to compute the W in the WY representation,
 *   on real data.
 *   Wraps the timer and the print statements if verbose.
 *   If verbose, then the W matrix is returned on the host.
 *
 * ON ENTRY :
 *   nrows    number of rows in the matrices V, Y, and W;
 *   ncols    equals the size of one tile, or equivalently,
 *            is the number of elements in B,
 *            and the number of columns in V, Y, and W;
 *   szt      size of one tile and the number of threads in a block;
 *   idx      index of the current tile;
 *   Vre_h    space for the real parts of the Householder vectors,
 *            if verbose;
 *   Vim_h    space for the imaginary parts of the Householder vectors,
 *            if verbose;
 *   Vre_d    the real parts of the Householder vectors on the device;
 *   Vim_d    the imaginary parts of the Householder vectors,
 *            on the device;
 *   Wre_h    space for the real parts of the W matrix, if verbose;
 *   Wim_h    space for the imaginary parts of the W matrix, if verbose;
 *   Wre_d    space for the real parts of W on the device;
 *   Wim_d    space for the imaginary parts of W on the device;
 *   WYHre_h  has space for the real parts of the outer product
 *            of W with Y^H, if verbose;
 *   WYHim_h  has space for the imaginary parts of the outer product
 *            of W with Y^H, if verbose;
 *   WYHre_d  has space for the real parts of the outer product
 *            of W with Y^H, on the device;
 *   WYHim_d  has space for the imaginary parts of the outer product
 *            of W with Y^H, on the device;
 *   beta_h   space for the betas if verbose;
 *   beta_d   space on the device for the betas;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Vre_h    the real parts of the Y matrix, if verbose;
 *   Vim_h    the imaginary parts of the Y matrix, if verbose;
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
 *   WYHre_h  has the real parts of the outer product of W with Y^H,
 *            if verbose;
 *   WYHim_h  has the imaginary parts of the outer product of W with Y^H,
 *            if verbose;
 *   WYHre_d  has the real parts of the outer product of W with Y^H,
 *            on the device;
 *   WYHim_d  has the imaginary parts of the outer product of W with Y^H,
 *            on the device;
 *   lapms    elapsed time spent by the kernel;
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications. */

void GPU_dbl_small_WYT
 ( int nrows, int szt, double *W_d, double *Y_d, double *WYT_d,
   double *WYT_h, double *lapms, long long int *add, long long int *mul,
   long long int *div, bool verbose=true );
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

void GPU_cmplx_small_WYH
 ( int nrows, int szt, double *Wre_d, double *Wim_d,
   double *Yre_d, double *Yim_d, double *WYHre_d, double *WYHim_d,
   double *WYHre_h, double *WYHim_h, double *lapms,
   long long int *add, long long int *mul, long long int *div,
   bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernel to compute W*Y^H.
 *   Wraps the timer and the print statements if verbose.
 *   If verbose, then the W*Y^H matrix is returned.
 *
 * ON ENTRY :
 *   nrows    number of rows in all matrices;
 *   szt      size of one tile and the number of threads in a block;
 *   Wre_d    the real parts of the matrix W in the WY representation;
 *   Wim_d    the imaginary parts of the matrix W in the WY representation;
 *   Yre_d    the matrix of the real parts of the Householder vectors;
 *   Yim_d    the matrix of the imaginary parts of the Householder vectors;
 *   WYHre_d  has space for an nrows-by-nrows matrix on the device;
 *   WYHim_d  has space for an nrows-by-nrows matrix on the device;
 *   WYHre_h  has space for an nrows-by-nrows matrix on the host, if verbose;
 *   WYHim_h  has space for an nrows-by-nrows matrix on the host, if verbose;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications;
 *   div      current number of divisions;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   WYHre_d  are the real parts of the product W*Y^H on the device;
 *   WYHim_d  are the imaginary parts of the product W*Y^H on the device;
 *   WYHre_h  are the real parts of the product W*Y^H, if verbose;
 *   WYHim_h  are the imaginary parts of the product W*Y^H, if verbose;
 *   lapms    elapsed time spent by the kernel;
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications;
 *   div      accumulated number of divisions. */

void GPU_dbl_small_YWT
 ( int nrows, int szt, int idx, double *Y_d, double *W_d, double *YWT_d,
   double *YWT_h, double *lapms, long long int *add, long long int *mul,
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
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   YWT_d    the product Y*W^T on the device;
 *   YWT_h    the product Y*W^T, if verbose;
 *   lapms    elapsed time spent by the kernel;
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications. */

void GPU_cmplx_small_YWH
 ( int nrows, int szt, int idx,
   double *Yre_d, double *Yim_d, double *Wre_d, double *Wim_d,
   double *YWHre_d, double *YWHim_d, double *YWHre_h, double *YWHim_h,
   double *lapms, long long int *add, long long int *mul, bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernel to compute W*Y^H.
 *   Wraps the timer and the print statements if verbose.
 *   If verbose, then the W*Y^H matrix is returned.
 *
 * ON ENTRY :
 *   nrows    number of rows in all matrices;
 *   szt      size of one tile and the number of threads in a block;
 *   idx      index of the current tile;
 *   Yre_d    real parts of the matrix of Householder vectors;
 *   Yim_d    imaginary parts of the matrix of Householder vectors;
 *   Wre_d    real parts of the matrix W in the WY representation;
 *   Wim_d    imaginary parts of the matrix W in the WY representation;
 *   YWHre_d  space for an nrows-by-nrows matrix on the device;
 *   YWHim_d  space for an nrows-by-nrows matrix on the device;
 *   YWHre_h  space for an nrows-by-nrows matrix on the host, if verbose;
 *   YWHim_h  space for an nrows-by-nrows matrix on the host, if verbose;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   YWHre_d  real parts of the product Y*W^H on the device;
 *   YWHim_d  imaginary parts of the product Y*W^H on the device;
 *   YWHre_h  real parts of the product Y*W^H, if verbose;
 *   YWHim_h  imaginary parts of the product Y*W^H, if verbose;
 *   lapms    elapsed time spent by the kernel;
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications. */

void GPU_dbl_small_QWYT
 ( int dim, int szt, int idx, double *Q_d, double *WYT_d, double *QWYT_d,
   double *QWYT_h, double *Q_h,
   double *lapms, long long int *add, long long int *mul, bool verbose=true );
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
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   QWYT_d   the product Q*WYT on the device;
 *   QWYT_h   the product Q*WYT, if verbose;
 *   lapms    elapsed time spent by the kernel;
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications. */

void GPU_cmplx_small_QWYH
 ( int dim, int szt, int idx, double *Qre_d, double *Qim_d,
   double *WYHre_d, double *WYHim_d, double *QWYHre_d, double *QWYHim_d,
   double *QWYHre_h, double *QWYHim_h, double *Qre_h, double *Qim_h,
   double *lapms, long long int *add, long long int *mul, bool verbose=true );
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
 *   WYHre_d  are the real parts of the product W*Y^H, on the device;
 *   WYHim_d  are the imaginary parts of the product W*Y^H, on the device;
 *   QWYHre_d has space for the real parts for the product Q*WYH,
 *            on the device;
 *   QWYHim_d has space for the imaginary parts of the product Q*WYH,
 *            on the device;
 *   QWYHre_h has space for the real parts of the product Q*WYH,
 *            on the host, if verbose;
 *   QWYHim_h has space for the imaginary parts of the product Q*WYH,
 *            on the host, if verbose;
 *   Qre_h    if verbose, then used to print the real parts
 *            of the matrix Q before the product;
 *   Qim_h    if verbose, then used to print the imaginary parts
 *            of the matrix Q before the product;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   QWYHre_d are the real parts of the product Q*WYH on the device;
 *   QWYHim_d are the imaginary parts of the product Q*WYH on the device;
 *   QWYHre_h are the real parts of the product Q*WYH, if verbose;
 *   QWYHim_h are the imaginary parts of the product Q*WYH, if verbose;
 *   lapms    elapsed time spent by the kernel;
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications. */

void GPU_dbl_small_YWTC
 ( int nrows, int ncols, int szt, int idx,
   double *YWT_d, double *C_d, double *YWTC_d, double *YWTC_h,
   double *lapms, long long int *add, long long int *mul, bool verbose=true );
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
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   YWTC_d   the product YWT*C on the device;
 *   YWTC_h   the product YWT*C, if verbose;
 *   lapms    elapsed time spent by the kernel;
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications. */

void GPU_cmplx_small_YWHC
 ( int nrows, int ncols, int szt, int idx,
   double *YWHre_d, double *YWHim_d, double *Cre_d, double *Cim_d,
   double *YWHCre_d, double *YWHCim_d, double *YWHCre_h, double *YWHCim_h,
   double *lapms, long long int *add, long long int *mul, bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernel to compute YWT*C.
 *   Wraps the timer and the print statement if verbose.
 *   If verbose, then the YWT*C matrix is returned.
 *
 * ON ENTRY :
 *   nrows    number of rows of the matrices C, YWH, YWHC,
 *            and the number of columns of the matrix YWH;
 *   ncols    number of columns of the matrix C;
 *   szt      size of one tile and the number of threads in a block;
 *   idx      index of the current tile;
 *   YWHre_d  are the real parts of the product Y*W^H, on the device;
 *   YWHim_d  are the imaginary parts of the product Y*W^H, on the device;
 *   Cre_d    real parts of an nrows-by-ncols matrix, on the device;
 *   Cim_d    imaginary parts of an nrows-by-ncols matrix, on the device;
 *   YWHCre_d has space for the real parts of the product YWH*C,
 *            on the device;
 *   YWHCim_d has space for the imaginary parts of the product YWH*C,
 *            on the device;
 *   YWHCre_h has space for the real parts of the product YWH*C,
 *            on the host, if verbose;
 *   YWHCim_h has space for the imaginary parts of the product YWH*C,
 *            on the host, if verbose;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   YWHCre_d are the real parts of the product YWH*C on the device;
 *   YWHCim_d are the imaginary parts of the product YWH*C on the device;
 *   YWHCre_h are the real parts of the product YWH*C, if verbose;
 *   YWHCim_h are the imaginary parts of the product YWH*C, if verbose;
 *   lapms    elapsed time spent by the kernel;
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications. */

void GPU_dbl_small_Qupdate
 ( int dim, int rowdim, int szt, int idx, double *Q_d, double *QWYT_d,
   double *Q_h, double *lapms, long long int *add, bool verbose=true );
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
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Q_d      the updated Q on the device;
 *   Q_h      the updated Q, if verbose;
 *   lapms    elapsed time spent by the kernel;
 *   add      accumulated number of additions and subtractions. */

void GPU_cmplx_small_Qupdate
 ( int dim, int szt, int idx, double *Qre_d, double *Qim_d,
   double *QWYHre_d, double *QWYHim_d, double *Qre_h, double *Qim_h,
   double *lapms, long long int *add, bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernel to update Q as Q + QWYH.
 *   Wraps the timer and the print statement if verbose.
 *   If verbose, then the updated Q is returned.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in Q,
 *            number of rows in QWYH;
 *   rowdim   number of columns in QWYH;
 *   szt      size of one tile and the number of threads in a block;
 *   idx      index of the current tile;
 *   Qre_d    real parts of a dim-by-dim matrix, on the device;
 *   Qim_d    imaginary parts of a dim-by-dim matrix, on the device;
 *   QWYHre_d are the real parts of the product Q*W*Y^H, on the device;
 *   QWYHim_d are the imaginary parts of the product Q*W*Y^H, on the device;
 *   add      current number of additions and subtractions;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Qre_d    real parts of the updated Q on the device;
 *   Qim_d    imaginary parts of the updated Q on the device;
 *   Qre_h    real parts of the updated Q, if verbose;
 *   Qim_h    imaginary parts of the updated Q, if verbose;
 *   lapms    elapsed time spent by the kernel;
 *   add      accumulated number of additions and subtractions. */

void GPU_dbl_small_R_add_YWTC
 ( int nrows, int ncols, int szt, int idx, double *R_d, double *YWTC_d,
   double *R_h, double *lapms, long long int *add, bool verbose=true );
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
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   R_d      the updated R on the device;
 *   R_h      the updated R, if verbose;
 *   lapms    elapsed time spent by the kernel;
 *   add      accumulated number of additions and subtractions. */

void GPU_cmplx_small_R_add_YWHC
 ( int nrows, int ncols, int szt, int idx,
   double *Rre_d, double *Rim_d, double *YWHCre_d, double *YWHCim_d,
   double *Rre_h, double *Rim_h,
   double *lapms, long long int *add, bool verbose=true );
/*
 * DESCRIPTION :
 *   Calls the kernel to update R as R + YWHC.
 *   Wraps the timer and the print statement if verbose.
 *   If verbose, then the updated R is returned.
 *
 * ON ENTRY :
 *   nrows    total number of rows in R and YWHC;
 *   ncols    total number of columns in R and YWHC;
 *   szt      size of one tile and the number of threads in a block;
 *   idx      index of the current tile;
 *   Rre_d    real parts of an nrows-by-ncols matrix, on the device;
 *   Rim_d    imaginary parts of an nrows-by-ncols matrix, on the device;
 *   YWHCre_d are the real parts of the product Y*W^H*C, on the device;
 *   YWHCim_d are the imaginary parts of the product Y*W^H*C,
 *            on the device;
 *   add      current number of additions and subtractions;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Rre_d    real parts of the updated R on the device;
 *   Rim_d    imaginary parts of the updated R on the device;
 *   Rre_h    real parts of the updated R, if verbose;
 *   Rim_h    imaginary parts of the updated R, if verbose;
 *   lapms    elapsed time spent by the kernel;
 *   add      accumulated number of additions and subtractions. */

void GPU_dbl_blocked_houseqr
 ( int nrows, int ncols, int szt, int nbt,
   double **A, double **Q, double **R,
   double *houselapms, double *RTvlapms, double *tileRlapms,
   double *vb2Wlapms, double *WYTlapms, double *QWYTlapms, double *Qaddlapms,
   double *YWTlapms, double *YWTClapms, double *Raddlapms,
   double *walltimesec, long long int *addcnt, long long int *mulcnt,
   long long int *divcnt, long long int *sqrtcnt, bool verbose=true );
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
 *   RTvlapms is the elapsed time to compute beta*R^T*v;
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
 *   walltimesec is the elapsed wall clock computation time;
 *   addcnt   counts the number of additions and subtractions;
 *   mulcnt   counts the number of multiplications;
 *   divcnt   counts the number of divisions;
 *   sqrtcnt  counts the number of calls to sqrt(). */

void GPU_cmplx_blocked_houseqr
 ( int nrows, int ncols, int szt, int nbt,
   double **Are, double **Aim, double **Qre, double **Qim,
   double **Rre, double **Rim,
   double *houselapms, double *RHvlapms, double *tileRlapms,
   double *vb2Wlapms, double *WYHlapms, double *QWYHlapms, double *Qaddlapms,
   double *YWHlapms, double *YWHClapms, double *Raddlapms,
   double *walltimesec, long long int *addcnt, long long int *mulcnt,
   long long int *divcnt, long long int *sqrtcnt, bool verbose=true );
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
 *   RHvlapms is the elapsed time to compute beta*R^H*v;
 *   tileRlapms is the elapsed time spent by the kernel
 *            to reduce one tile;
 *   vb2Wlapms is the elapsed time spent by the kernel
 *            to compute the W representation;
 *   WYHlapms is the elapsed time spent by the kernel
 *            to compute the W*Y^H product;
 *   QWYHlapms is the elapsed time spent by the kernel
 *            to compute the Q*WYH product;
 *   Qaddlapms is the elapsed time spent by the kernel
 *            to compute Q by adding the Q*W*Y^H matrix;
 *   YWHlapms is the elapsed time spent by the kernel
 *            to compute the Y*W^H product;
 *   YWHClapms is the elapsed time spent by the kernel
 *            to compute the YWH*C product;
 *   Raddlapms is the elapsed time spent by the kernel
 *            to compute R by adding the Y*W^H*C matrix;
 *   walltimesec is the elapsed wall clock computation time;
 *   addcnt   counts the number of additions and subtractions;
 *   mulcnt   counts the number of multiplications;
 *   divcnt   counts the number of divisions;
 *   sqrtcnt  counts the number of calls to sqrt(). */

#endif
