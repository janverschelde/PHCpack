/* The file dbl2_baqr_host.h specifies functions on the host for the
 * blocked accelerated qr decomposition in double double precision. */

#ifndef __dbl2_baqr_host_h__
#define __dbl2_baqr_host_h__

void CPU_dbl2_blocked_VB_to_W
 ( int nrows, int ncols, double *Bhi, double *Blo,
   double **Vhi, double **Vlo, double **Whi, double **Wlo );
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
 *   Bhi      Bhi[i] is the high double of the i-th beta computed by house;
 *   Blo      Blo[i] is the low double of the i-th beta computed by house;
 *   Vhi      Vhi[i] is the high doubles of the i-th Householder vector,
 *            with i zeros inserted so V is trapezoidal;
 *   Vlo      Vlo[i] is the low double of the i-th Householder vector,
 *            with i zeros inserted so V is trapezoidal;
 *   Whi      space for ncols columns with rows from 0 to nrows-1;
 *   Wlo      space for ncols columns with rows from 0 to nrows-1.
 *
 * ON RETURN :
 *   Whi      high doubles of the W matrix in the WY representation;
 *   Wlo      low doubles of the W matrix in the WY representation. */

void CPU_dbl2_blocked_leftRupdate
 ( int nrows, int ncols, int szt, int idx, double **Chi, double **Clo,
   double **Yhi, double **Ylo, double **Whi, double **Wlo );
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
 *   Chi      high doubles of an nrows-by-ncols matrix;
 *   Clo      low doubles of an nrows-by-ncols matrix;
 *   Yhi      high doubles of the matrix of Householder vectors;
 *   Ylo      low doubles of the matrix of Householder vectors;
 *   Whi      high doubles of the W matrix in the WY representation;
 *   Wlo      low doubles of the W matrix in the WY representation.
 *
 * ON RETURN :
 *   Chi      high doubles of the update with the Householder matrix;
 *   Clo      low doubles of the update with the Householder matrix. */

void CPU_dbl2_blocked_rightQupdate
 ( int dim, int szt, int idx, double **Qhi, double **Qlo,
   double **Yhi, double **Ylo, double **Whi, double **Wlo,
   bool verbose=true );
/*
 * DESCRIPTION :
 *   Applies the Householder matrix to Q, on real data.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in the matrix Q;
 *   szt      size of each block;
 *   idx      index of the current block;
 *   Qhi      high doubles of a dim-by-dim matrix;
 *   Qlo      low doubles of a dim-by-dim matrix;
 *   Yhi      high doubles of the matrix of Householder vectors;
 *   Ylo      low doubles of the matrix of Householder vectors;
 *   Whi      high doubles of the W matrix in the WY representation;
 *   Wlo      low doubles of the W matrix in the WY representation;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Qhi      high doubles of the update with the Householder matrix;
 *   Qlo      low doubles of the update with the Householder matrix. */

void CPU_dbl2_blocked_houseqr
 ( int nrows, int ncols, int szt, int nbt, double **Ahi, double **Alo,
   double **Qhi, double **Qlo, double **Rhi, double **Rlo,
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
 *   Rlo      low doubles of the reduced upper triangular form of A. */

#endif
