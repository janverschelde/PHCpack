/* The file dbl_baqr_host.h specifies functions on the host for the
 * blocked accelerated qr decomposition in double precision. */

#ifndef __dbl_baqr_host_h__
#define __dbl_baqr_host_h__

void CPU_dbl_blocked_VB_to_W
 ( int nrows, int ncols, double *B, double **V, double **W );
/*
 * DESCRIPTION :
 *   Computes the W in the WY representation of the Householder
 *   transformations defined by V and B.
 *
 * ON ENTRY :
 *   nrows    number of rows in the matrices V, Y, and W;
 *   ncols    equals the size of one tile, or equivalently,
 *            is the number of elements in B,
 *            and the number of columns in V, Y, and W;
 *   B        B[i] is the i-th beta computed by house;
 *   V        V[i] is the i-th Householder vector,
 *            with i zeros inserted so V is trapezoidal;
 *   W        space for ncols columns with rows from 0 to nrows-1.
 *
 * ON RETURN :
 *   W        the W matrix in the WY representation. */

void CPU_dbl_blocked_leftRupdate
 ( int nrows, int ncols, int szt, int idx,
   double **C, double **Y, double **W, bool verbose=true );
/*
 * DESCRIPTION :
 *   Applies the blocked Householder update to C,
 *   as C = C + (Y*W^T)*C.
 *
 * ON ENTRY :
 *   nrows    number of rows in the matrix C;
 *   ncols    number of columns in the matrix C;
 *   szt      size of each block;
 *   idx      index of the current block;
 *   C        an nrows-by-ncols matrix;
 *   Y        matrix of Householder vectors;
 *   W        the W matrix in the WY representation;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   C        update with the Householder matrix. */

void CPU_dbl_blocked_rightQupdate
 ( int dim, int szt, int idx, double **Q, double **Y, double **W,
   bool verbose=true );
/*
 * DESCRIPTION :
 *   Applies the Householder matrix to Q.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in the matrix Q;
 *   szt      size of each block;
 *   idx      index of the current block;
 *   Q        a dim-by-dim matrix;
 *   Y        matrix of Householder vectors;
 *   W        the W matrix in the WY representation;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Q        update with the Householder matrix. */

void CPU_dbl_blocked_houseqr
 ( int nrows, int ncols, int szt, int nbt,
   double **A, double **Q, double **R, bool verbose=true );
/*
 * DESCRIPTION :
 *   Applies Householder transformations in a blocked manner
 *   to compute a QR decomposition of A.
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
 *   R        the reduced upper triangular form of A. */

#endif
