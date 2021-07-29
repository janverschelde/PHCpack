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

void CPU_cmplx2_blocked_VB_to_W
 ( int nrows, int ncols, double *Bhi, double *Blo,
   double **Vrehi, double **Vrelo, double **Vimhi, double **Vimlo,
   double **Wrehi, double **Wrelo, double **Wimhi, double **Wimlo );
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
 *   Bhi      B[i] is the high double of the i-th beta computed by house;
 *   Blo      B[i] is the low double of the i-th beta computed by house;
 *   Vrehi    high doubles of the real parts of the Householder vectors,
 *            with zeros inserted so V is trapezoidal;
 *   Vrelo    low doubles of the real parts of the Householder vectors,
 *            with zeros inserted so V is trapezoidal;
 *   Vimhi    high doubles of the imaginary parts of the Householder vectors,
 *            with zeros inserted so V is trapezoidal;
 *   Vimlo    low doubles of the imaginary parts of the Householder vectors,
 *            with zeros inserted so V is trapezoidal;
 *   Wrehi    space for ncols columns with rows from 0 to nrows-1;
 *   Wrelo    space for ncols columns with rows from 0 to nrows-1;
 *   Wimhi    space for ncols columns with rows from 0 to nrows-1;
 *   Wimlo    space for ncols columns with rows from 0 to nrows-1.
 *
 * ON RETURN :
 *   Wrehi    high doubles of the real parts of the W matrix 
 *            in the WY representation;
 *   Wrelo    low doubles of the real parts of the W matrix
 *            in the WY representation;
 *   Wimhi    high doubles of the imaginary parts of the W matrix
 *            in the WY representation;
 *   Wimlo    low doubles of the imaginary parts of the W matrix
 *            in the WY representation. */

void CPU_dbl2_blocked_leftRupdate
 ( int nrows, int ncols, int szt, int idx, double **Chi, double **Clo,
   double **Yhi, double **Ylo, double **Whi, double **Wlo,
   bool verbose=true );
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
 *   Wlo      low doubles of the W matrix in the WY representation;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Chi      high doubles of the update with the Householder matrix;
 *   Clo      low doubles of the update with the Householder matrix. */

void CPU_cmplx2_blocked_leftRupdate
 ( int nrows, int ncols, int szt, int idx,
   double **Crehi, double **Crelo, double **Cimhi, double **Cimlo,
   double **Yrehi, double **Yrelo, double **Yimhi, double **Yimlo,
   double **Wrehi, double **Wrelo, double **Wimhi, double **Wimlo,
   bool verbose=true );
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
 *   Crehi    high doubles of the real parts of an nrows-by-ncols matrix;
 *   Crelo    low doubles of the real parts of an nrows-by-ncols matrix;
 *   Cimhi    high doubles of the imaginary parts of an nrows-by-ncols matrix;
 *   Cimlo    low doubles of the imaginary parts of an nrows-by-ncols matrix;
 *   Yrehi    high doubles of the real parts of the Householder vectors;
 *   Yrelo    low doubles of the real parts of the Householder vectors;
 *   Yimhi    high doubles of the imaginary parts of the Householder vectors;
 *   Yimlo    low doubles of the imaginary parts of the Householder vectors;
 *   Wrehi    high doubles of the real parts of the W matrix
 *            in the WY representation;
 *   Wrelo    low doubles of the real parts of the W matrix
 *            in the WY representation;
 *   Wimhi    high doubles of the imaginary parts of the W matrix
 *            in the WY representation;
 *   Wimlo    low doubles of the imaginary parts of the W matrix
 *            in the WY representation;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Crehi    high doubles of the real parts of the update
 *            with the Householder matrix;
 *   Crelo    low doubles of the real parts of the update
 *            with the Householder matrix;
 *   Cimhi    high doubles of the imaginary parts of the update
 *            with the Householder matrix;
 *   Cimlo    low doubles of the imaginary parts of the update
 *            with the Householder matrix. */

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

void CPU_cmplx2_blocked_rightQupdate
 ( int dim, int szt, int idx,
   double **Qrehi, double **Qrelo, double **Qimhi, double **Qimlo,
   double **Yrehi, double **Yrelo, double **Yimhi, double **Yimlo,
   double **Wrehi, double **Wrelo, double **Wimhi, double **Wimlo,
   bool verbose=true );
/*
 * DESCRIPTION :
 *   Applies the Householder matrix to Q, on complex data.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in the matrix Q;
 *   szt      size of each block;
 *   idx      index of the current block;
 *   Qrehi    high doubles of the real parts of a dim-by-dim matrix;
 *   Qrelo    low doubles of the real parts of a dim-by-dim matrix;
 *   Qimhi    high doubles of the imaginary parts of a dim-by-dim matrix;
 *   Qimlo    lowh doubles of the imaginary parts of a dim-by-dim matrix;
 *   Yrehi    high doubles of the real parts of the Householder vectors;
 *   Yrelo    low doubles of the real parts of the Householder vectors;
 *   Yimhi    high doubles of the imaginary parts of the Householder vectors;
 *   Yimlo    low doubles of the imaginary parts of the Householder vectors;
 *   Wrehi    high doubles of the real parts of the W matrix
 *            in the WY representation;
 *   Wrelo    low doubles of the real parts of the W matrix
 *            in the WY representation;
 *   Wimhi    high doubles of the imaginary parts of the W matrix
 *            in the WY representation;
 *   Wimlo    low doubles of the imaginary parts of the W matrix
 *            in the WY representation;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Qrehi    high doubles of the real parts of the update
 *            with the Householder matrix;
 *   Qrelo    low doubles of the real parts of the update
 *            with the Householder matrix;
 *   Qimhi    high doubles of the imaginary parts of the update
 *            with the Householder matrix;
 *   Qimlo    low doubles of the imaginary parts of the update
 *            with the Householder matrix. */

void CPU_dbl2_blocked_houseqr
 ( int nrows, int ncols, int szt, int nbt, double **Ahi, double **Alo,
   double **Qhi, double **Qlo, double **Rhi, double **Rlo,
   double *lapsed, bool verbose=true );
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
 *   lapsec   elapsed time in seconds. */

void CPU_cmplx2_blocked_houseqr
 ( int nrows, int ncols, int szt, int nbt,
   double **Arehi, double **Arelo, double **Aimhi, double **Aimlo,
   double **Qrehi, double **Qrelo, double **Qimhi, double **Qimlo,
   double **Rrehi, double **Rrelo, double **Rimhi, double **Rimlo,
   double *lapsec, bool verbose=true );
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
 *   Arehi    high doubles of the real parts of an nrows-by-ncols matrix,
 *            stored as nrows arrays of ncols numbers;
 *   Arelo    low doubles of the real parts of an nrows-by-ncols matrix,
 *            stored as nrows arrays of ncols numbers;
 *   Aimhi    high doubles of the imaginary parts of an nrows-by-ncols matrix,
 *            stored as nrows arrays of ncols numbers;
 *   Aimlo    low doubles of the imaginary parts of an nrows-by-ncols matrix,
 *            stored as nrows arrays of ncols numbers;
 *   Qrehi    space for an nrows-by-nrows matrix;
 *   Qrelo    space for an nrows-by-nrows matrix;
 *   Qimhi    space for an nrows-by-nrows matrix;
 *   Qimlo    space for an nrows-by-nrows matrix;
 *   Rrehi    space for an nrows-by-ncols matrix;
 *   Rrelo    space for an nrows-by-ncols matrix;
 *   Rimhi    space for an nrows-by-ncols matrix;
 *   Rimlo    space for an nrows-by-ncols matrix;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Qrehi    high doubles of the real parts of an orthogonal matrix,
 *            transpose(Q)*A = R;
 *   Qrelo    low doubles of the real parts of an orthogonal matrix,
 *            transpose(Q)*A = R;
 *   Qimhi    high doubles of the imaginary parts of an orthogonal matrix,
 *            transpose(Q)*A = R;
 *   Qimlo    low doubles of the imaginary parts of an orthogonal matrix,
 *            transpose(Q)*A = R;
 *   Rrehi    high doubles of the real parts of the reduced 
 *            upper triangular form of A;
 *   Rrelo    low doubles of the real parts of the reduced
 *            upper triangular form of A;
 *   Rimhi    high doubles of the imaginary parts of the reduced
 *            upper triangular form of A;
 *   Rimlo    low doubles of the imaginary parts of the reduced
 *            upper triangular form of A;
 *   lapsec   elapsed time in seconds. */


#endif
