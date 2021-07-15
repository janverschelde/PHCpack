/* The file dbl_baqr_kernels.h specifies functions for the
 * blocked accelerated QR in double precision. */

#ifndef __dbl_baqr_kernels_h__
#define __dbl_baqr_kernels_h__

#define d_shmemsize 1024

__global__ void dbl_small_house
 ( double x0, double *x1, int dim, int dimLog2, double *v, double *beta );
/*
 * DESCRIPTION :
 *   Computes the Householder vector of a vector of dimension dim+1.
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

void GPU_dbl_blocked_qr
 ( int nrows, int ncols, int szt, int nbt,
   double **A, double **Q, double **R,
   double *lapms, double *walltimesec );
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
 *   R        space for an nrows-by-ncols matrix.
 *
 * ON RETURN :
 *   Q        an orthogonal matrix, transpose(Q)*A = R;
 *   R        the reduced upper triangular form of A;
 *   lapms    elapsed time spent by the kernels;
 *   walltimesec is the elapsed wall clock computation time. */

#endif
