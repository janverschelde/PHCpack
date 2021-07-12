/* The file dbl_baqr_host.h specifies functions on the host for the
 * blocked accelerated qr decomposition in double precision. */

#ifndef __dbl_baqr_host_h__
#define __dbl_baqr_host_h__

void CPU_dbl_blocked_houseqr
 ( int nrows, int ncols, int szt, int nbt,
   double **A, double **Q, double **R );
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
 *   R        space for an nrows-by-ncols matrix.
 *
 * ON RETURN :
 *   Q        an orthogonal matrix, transpose(Q)*A = R;
 *   R        the reduced upper triangular form of A. */

#endif
