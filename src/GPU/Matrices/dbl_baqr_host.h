/* The file dbl_baqr_host.h specifies functions on the host for the
 * blocked accelerated qr decomposition in double precision. */

#ifndef __dbl_baqr_host_h__
#define __dbl_baqr_host_h__

void CPU_dbl_blocked_VB_to_W
 ( int nrows, int ncols, double *B, double **V, double **W );
/*
 * DESCRIPTION :
 *   Computes the W in the WY representation of the Householder
 *   transformations defined by V and B, on real data.
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

void CPU_cmplx_blocked_VB_to_W
 ( int nrows, int ncols, double *B,
   double **Vre, double **Vim, double **Wre, double **Wim );
/*
 * DESCRIPTION :
 *   Computes the W in the WY representation of the Householder
 *   transformations defined by V and B, on complex data.
 *
 * ON ENTRY :
 *   nrows    number of rows in the matrices V, Y, and W;
 *   ncols    equals the size of one tile, or equivalently,
 *            is the number of elements in B,
 *            and the number of columns in V, Y, and W;
 *   B        B[i] is the i-th beta computed by house;
 *   Vre      real parts of the V[i] is the i-th Householder vector,
 *            with i zeros inserted so V is trapezoidal;
 *   Vim      imaginary parts of the V[i] is the i-th Householder vector,
 *            with i zeros inserted so V is trapezoidal;
 *   Wre      space for ncols columns with rows from 0 to nrows-1;
 *   Wim      space for ncols columns with rows from 0 to nrows-1.
 *
 * ON RETURN :
 *   Wre      real parts of the W matrix in the WY representation;
 *   Wim      imaginary parts of the W matrix in the WY representation. */

void CPU_dbl_blocked_leftRupdate
 ( int nrows, int ncols, int szt, int idx,
   double **C, double **Y, double **W, bool verbose=true );
/*
 * DESCRIPTION :
 *   Applies the blocked Householder update to C,
 *   as C = C + (Y*W^T)*C, on real data.
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

void CPU_cmplx_blocked_leftRupdate
 ( int nrows, int ncols, int szt, int idx,
   double **Cre, double **Cim, double **Yre, double **Yim,
   double **Wre, double **Wim, bool verbose=true );
/*
 * DESCRIPTION :
 *   Applies the blocked Householder update to C,
 *   as C = C + (Y*W^T)*C, on complex data.
 *
 * ON ENTRY :
 *   nrows    number of rows in the matrix C;
 *   ncols    number of columns in the matrix C;
 *   szt      size of each block;
 *   idx      index of the current block;
 *   Cre      real parts of an nrows-by-ncols matrix;
 *   Cim      imaginary parts of an nrows-by-ncols matrix;
 *   Yre      real parts of the matrix of Householder vectors;
 *   Yim      imaginary parts of the matrix of Householder vectors;
 *   Wre      real parts of the W matrix in the WY representation;
 *   Wim      imaginary parts of the W matrix in the WY representation;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Cre      real parts of the update with the Householder matrix;
 *   Cim      imaginary parts of the update with the Householder matrix. */

void CPU_dbl_blocked_rightQupdate
 ( int dim, int szt, int idx, double **Q, double **Y, double **W,
   bool verbose=true );
/*
 * DESCRIPTION :
 *   Applies the Householder matrix to Q, on real data.
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

void CPU_cmplx_blocked_rightQupdate
 ( int dim, int szt, int idx, double **Qre, double **Qim,
   double **Yre, double **Yim, double **Wre, double **Wim,
   bool verbose=true );
/*
 * DESCRIPTION :
 *   Applies the Householder matrix to Q, on complex data.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in the matrix Q;
 *   szt      size of each block;
 *   idx      index of the current block;
 *   Qre      real parts of a dim-by-dim matrix;
 *   Qim      imaginary parts of a dim-by-dim matrix;
 *   Yre      real parts of the matrix of Householder vectors;
 *   Yim      imaginary parts of the matrix of Householder vectors;
 *   Wre      real parts of the W matrix in the WY representation;
 *   Wim      imaginary parts of the W matrix in the WY representation;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Qre      real parts of the update with the Householder matrix;
 *   Qim      imaginary parts of the update with the Householder matrix. */

void CPU_dbl_blocked_houseqr
 ( int nrows, int ncols, int szt, int nbt,
   double **A, double **Q, double **R, double *lapsec,
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
 *   lapsec   elapsed time in seconds. */

void CPU_cmplx_blocked_houseqr
 ( int nrows, int ncols, int szt, int nbt,
   double **Are, double **Aim, double **Qre, double **Qim,
   double **Rre, double **Rim, double *lapsec,
   bool verbose=true );
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
 *   Qim      imaginary parts of an orthogonal matrix, transpose(Q)*A = R;
 *   Rre      real parts of the reduced upper triangular form of A;
 *   Rim      imaginary parts of the reduced upper triangular form of A;
 *   lapsec   elapsed time in seconds. */

#endif
