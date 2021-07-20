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

void GPU_dbl_blocked_houseqr
 ( int nrows, int ncols, int szt, int nbt,
   double **A, double **Q, double **R,
   double *houselapms, double *tileRlapms, double *vb2Wlapms,
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
 *   walltimesec is the elapsed wall clock computation time. */

#endif
