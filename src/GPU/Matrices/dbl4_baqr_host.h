/* The file dbl4_baqr_host.h specifies functions on the host for the
 * blocked accelerated qr decomposition in quad double precision. */

#ifndef __dbl4_baqr_host_h__
#define __dbl4_baqr_host_h__

void CPU_dbl4_blocked_VB_to_W
 ( int nrows, int ncols,
   double *Bhihi, double *Blohi, double *Bhilo, double *Blolo,
   double **Vhihi, double **Vlohi, double **Vhilo, double **Vlolo,
   double **Whihi, double **Wlohi, double **Whilo, double **Wlolo );
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
 *   Bhihi    highest doubles of the beta computed by house;
 *   Blohi    second highest doubles of the beta computed by house;
 *   Bhilo    second lowest doubles of the beta computed by house;
 *   Blolo    lowest doubles of the beta computed by house;
 *   Vhihi    highest doubles of the Householder vectors,
 *            with zeros inserted so V is trapezoidal;
 *   Vlohi    second highest doubles of the Householder vectors,
 *            with zeros inserted so V is trapezoidal;
 *   Vhilo    second lowest double of the Householder vectors,
 *            with zeros inserted so V is trapezoidal;
 *   Vlolo    lowest double of the Householder vectors,
 *            with zeros inserted so V is trapezoidal;
 *   Whihi    space for ncols columns with rows from 0 to nrows-1;
 *   Wlohi    space for ncols columns with rows from 0 to nrows-1;
 *   Whilo    space for ncols columns with rows from 0 to nrows-1;
 *   Wlolo    space for ncols columns with rows from 0 to nrows-1.
 *
 * ON RETURN :
 *   Whihi    highest doubles of W in the WY representation;
 *   Wlohi    second highest doubles of W;
 *   Whilo    second lowest doubles of W;
 *   Wlolo    lowest doubles of W. */

void CPU_cmplx4_blocked_VB_to_W
 ( int nrows, int ncols,
   double *Bhihi, double *Blohi, double *Bhilo, double *Blolo,
   double **Vrehihi, double **Vrelohi, double **Vrehilo, double **Vrelolo,
   double **Vimhihi, double **Vimlohi, double **Vimhilo, double **Vimlolo,
   double **Wrehihi, double **Wrelohi, double **Wrehilo, double **Wrelolo,
   double **Wimhihi, double **Wimlohi, double **Wimhilo, double **Wimlolo );
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
 *   Bhihi    highest doubles of the beta computed by house;
 *   Blohi    second highest doubles of the beta computed by house;
 *   Bhilo    second lowest doubles of the beta computed by house;
 *   Blolo    lowest doubles of the beta computed by house;
 *   Vrehihi  highest doubles of the real parts of the
 *            Householder vectors, with zeros inserted so V is trapezoidal;
 *   Vrelohi  second highest doubles of the real parts of the
 *            Householder vectors, with zeros inserted so V is trapezoidal;
 *   Vrehilo  second lowest doubles of the real parts of the
 *            Householder vectors, with zeros inserted so V is trapezoidal;
 *   Vrelolo  lowest doubles of the real parts of the
 *            Householder vectors, with zeros inserted so V is trapezoidal;
 *   Vimhihi  highest doubles of the imaginary parts of the
 *            Householder vectors, with zeros inserted so V is trapezoidal;
 *   Vimlohi  second highest doubles of the imaginary parts of the
 *            Householder vectors, with zeros inserted so V is trapezoidal;
 *   Vimhilo  second lowest doubles of the imaginary parts of the
 *            Householder vectors, with zeros inserted so V is trapezoidal;
 *   Vimlolo  lowest doubles of the imaginary parts of the
 *            Householder vectors, with zeros inserted so V is trapezoidal;
 *   Wrehihi  has space for ncols columns with rows from 0 to nrows-1;
 *   Wrelohi  has space for ncols columns with rows from 0 to nrows-1;
 *   Wrehilo  has space for ncols columns with rows from 0 to nrows-1;
 *   Wrelolo  has space for ncols columns with rows from 0 to nrows-1;
 *   Wimhihi  has space for ncols columns with rows from 0 to nrows-1;
 *   Wimlohi  has space for ncols columns with rows from 0 to nrows-1;
 *   Wimhilo  has space for ncols columns with rows from 0 to nrows-1;
 *   Wimlolo  has space for ncols columns with rows from 0 to nrows-1.
 *
 * ON RETURN :
 *   Wrehihi  highest doubles of the real parts of W;
 *   Wrelohi  second highest doubles of the real parts of W;
 *   Wrehilo  second lowest doubles of the real parts of W;
 *   Wrelolo  lowest doubles of the real parts of W;
 *   Wimhihi  highest doubles of the imaginary parts of W;
 *   Wimlohi  second highest doubles of the imaginary parts of W;
 *   Wimhilo  second lowest doubles of the imaginary parts of W;
 *   Wimlolo  lowest doubles of the imaginary parts of W. */

void CPU_dbl4_blocked_leftRupdate
 ( int nrows, int ncols, int szt, int idx,
   double **Chihi, double **Clohi, double **Chilo, double **Clolo,
   double **Yhihi, double **Ylohi, double **Yhilo, double **Ylolo,
   double **Whihi, double **Wlohi, double **Whilo, double **Wlolo,
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
 *   Chihi    highest doubles of an nrows-by-ncols matrix;
 *   Clohi    second highest doubles of an nrows-by-ncols matrix;
 *   Chilo    second lowest doubles of an nrows-by-ncols matrix;
 *   Clolo    lowest doubles of an nrows-by-ncols matrix;
 *   Yhihi    highest doubles of the matrix of Householder vectors;
 *   Ylohi    second highest doubles of the matrix of Householder vectors;
 *   Yhilo    second lowest doubles of the matrix of Householder vectors;
 *   Ylolo    lowest doubles of the matrix of Householder vectors;
 *   Whihi    highest doubles of W in the WY representation;
 *   Wlohi    second highest doubles of W in the WY representation;
 *   Whilo    second lowest doubles of W in the WY representation;
 *   Wlolo    lowest doubles of W in the WY representation;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Chihi    highest doubles of the update with the Householder matrix;
 *   Clohi    second highest doubles of the update;
 *   Chilo    second lowest doubles of the update;
 *   Clolo    lowest doubles of the update. */

void CPU_cmplx4_blocked_leftRupdate
 ( int nrows, int ncols, int szt, int idx,
   double **Crehihi, double **Crelohi, double **Crehilo, double **Crelolo,
   double **Cimhihi, double **Cimlohi, double **Cimhilo, double **Cimlolo,
   double **Yrehihi, double **Yrelohi, double **Yrehilo, double **Yrelolo,
   double **Yimhihi, double **Yimlohi, double **Yimhilo, double **Yimlolo,
   double **Wrehihi, double **Wrelohi, double **Wrehilo, double **Wrelolo,
   double **Wimhihi, double **Wimlohi, double **Wimhilo, double **Wimlolo,
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
 *   Crehihi  highest doubles of the real parts of C;
 *   Crelohi  second highest doubles of the real parts of C;
 *   Crehilo  second lowest doubles of the real parts of C;
 *   Crelolo  lowest doubles of the real parts of C;
 *   Cimhihi  highest doubles of the imaginary parts of C;
 *   Cimlohi  second highest doubles of the imaginary parts of C;
 *   Cimhilo  second lowest doubles of the imaginary parts of C;
 *   Cimlolo  lowest doubles of the imaginary parts of C;
 *   Yrehihi  highest doubles of the real parts of the Householder vectors Y;
 *   Yrelohi  second highest doubles of the real parts of Y;
 *   Yrehilo  second lowest doubles of the real parts of Y;
 *   Yrelolo  lowest doubles of the real parts of Y;
 *   Yimhihi  highest doubles of the imaginary parts of Y;
 *   Yimlohi  second highest doubles of the imaginary parts of Y;
 *   Yimhilo  second lowest doubles of the imaginary parts of Y;
 *   Yimlolo  lowest doubles of the imaginary parts of Y;
 *   Wrehihi  highest doubles of the real parts of W;
 *   Wrelohi  second highest doubles of the real parts of W;
 *   Wrehilo  second lowest doubles of the real parts of W;
 *   Wrelolo  lowest doubles of the real parts of W;
 *   Wimhihi  highest doubles of the imaginary parts of W;
 *   Wimlohi  second highest doubles of the imaginary parts of W;
 *   Wimhilo  second lowest doubles of the imaginary parts of W;
 *   Wimlolo  lowest doubles of the imaginary parts of W;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Crehihi  highest doubles of the real parts of the update;
 *   Crelohi  second highest doubles of the real parts of the update;
 *   Crehilo  second lowest doubles of the real parts of the update;
 *   Crelolo  lowest doubles of the real parts of the update;
 *   Cimhihi  highest doubles of the imaginary parts of the update;
 *   Cimlohi  second highest doubles of the imaginary parts of the update;
 *   Cimhilo  second lowest doubles of the imaginary parts of the update;
 *   Cimlolo  lowest doubles of the imaginary parts of the update. */

void CPU_dbl4_blocked_rightQupdate
 ( int dim, int szt, int idx,
   double **Qhihi, double **Qlohi, double **Qhilo, double **Qlolo,
   double **Yhihi, double **Ylohi, double **Yhilo, double **Ylolo,
   double **Whihi, double **Wlohi, double **Whilo, double **Wlolo,
   bool verbose=true );
/*
 * DESCRIPTION :
 *   Applies the Householder matrix to Q, on real data.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in the matrix Q;
 *   szt      size of each block;
 *   idx      index of the current block;
 *   Qhihi    highest doubles of Q;
 *   Qlohi    second highest doubles of Q;
 *   Qhilo    second lowest doubles of Q;
 *   Qlolo    lowest doubles of Q;
 *   Yhihi    highest doubles of the matrix of Householder vectors Y;
 *   Ylohi    second highest doubles of Y;
 *   Yhilo    second lowest doubles of Y;
 *   Ylolo    lowest doubles of Y;
 *   Whihi    highest doubles of the W matrix in the WY representation;
 *   Wlohi    second highest doubles of W;
 *   Whilo    second lowest doubles of W;
 *   Wlolo    lowest doubles of W;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Qhihi    highest doubles of the updated Q;
 *   Qlohi    second highest doubles of the updated Q;
 *   Qhilo    second lowest doubles of the updated Q;
 *   Qlolo    lowest doubles of the updated Q. */

void CPU_cmplx4_blocked_rightQupdate
 ( int dim, int szt, int idx,
   double **Qrehihi, double **Qrelohi, double **Qrehilo, double **Qrelolo,
   double **Qimhihi, double **Qimlohi, double **Qimhilo, double **Qimlolo,
   double **Yrehihi, double **Yrelohi, double **Yrehilo, double **Yrelolo,
   double **Yimhihi, double **Yimlohi, double **Yimhilo, double **Yimlolo,
   double **Wrehihi, double **Wrelohi, double **Wrehilo, double **Wrelolo,
   double **Wimhihi, double **Wimlohi, double **Wimhilo, double **Wimlolo,
   bool verbose=true );
/*
 * DESCRIPTION :
 *   Applies the Householder matrix to Q, on complex data.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in the matrix Q;
 *   szt      size of each block;
 *   idx      index of the current block;
 *   Qrehihi  highest doubles of the real parts of Q;
 *   Qrelohi  second highest doubles of the real parts of Q;
 *   Qrehilo  second lowest doubles of the real parts of Q;
 *   Qrelolo  lowest doubles of the real parts of Q;
 *   Qimhihi  highest doubles of the imaginary parts of Q;
 *   Qimlohi  second highest doubles of the imaginary parts of Q;
 *   Qimhilo  second lowest doubles of the imaginary parts of Q;
 *   Qimlolo  lowest doubles of the imaginary parts of Q;
 *   Yrehihi  highest doubles of the real parts of the Householder vectors Y;
 *   Yrelohi  second highest doubles of the real parts of Y;
 *   Yrehilo  second lowest doubles of the real parts of Y;
 *   Yrelolo  lowest doubles of the real parts of Y;
 *   Yimhihi  highest doubles of the imaginary parts of Y;
 *   Yimlohi  second highest doubles of the imaginary parts of Y;
 *   Yimhilo  second lowest doubles of the imaginary parts of Ys;
 *   Yimlolo  lowest doubles of the imaginary parts of Ys;
 *   Wrehihi  highest doubles of the real parts of W;
 *   Wrelohi  second highest doubles of the real parts of W;
 *   Wrehilo  second lowest doubles of the real parts of W;
 *   Wrelolo  lowest doubles of the real parts of W;
 *   Wimhihi  highest doubles of the imaginary parts of W;
 *   Wimlohi  second highest doubles of the imaginary parts of W;
 *   Wimhilo  second lowest doubles of the imaginary parts of W;
 *   Wimlolo  lowest doubles of the imaginary parts of W;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Qrehihi  highest doubles of the real parts of the updated Q;
 *   Qrelohi  second highest doubles of the real parts of the updated Q;
 *   Qrehilo  second lowest doubles of the real parts of the updated Q;
 *   Qrelolo  lowest doubles of the real parts of the updated Q;
 *   Qimhihi  highest doubles of the imaginary parts of the updated Q;
 *   Qimlohi  second highest doubles of the imaginary parts of the updated Q;
 *   Qimhilo  second lowest doubles of the imaginary parts of the updated Q;
 *   Qimlolo  lowest doubles of the imaginary parts of the updated Q. */

void CPU_dbl4_blocked_houseqr
 ( int nrows, int ncols, int szt, int nbt,
   double **Ahihi, double **Alohi, double **Ahilo, double **Alolo,
   double **Qhihi, double **Qlohi, double **Qhilo, double **Qlolo,
   double **Rhihi, double **Rlohi, double **Rhilo, double **Rlolo,
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
 *   Ahihi    highest doubles of A,
 *            stored as nrows arrays of ncols numbers;
 *   Alohi    second highest doubles of A;
 *   Ahilo    second lowest doubles of A;
 *   Alolo    lowest doubles of A;
 *   Qhihi    space for an nrows-by-nrows matrix;
 *   Qlohi    space for an nrows-by-nrows matrix;
 *   Qhilo    space for an nrows-by-nrows matrix;
 *   Qlolo    space for an nrows-by-nrows matrix;
 *   Rhihi    space for an nrows-by-ncols matrix;
 *   Rlohi    space for an nrows-by-ncols matrix;
 *   Rhilo    space for an nrows-by-ncols matrix;
 *   Rlolo    space for an nrows-by-ncols matrix;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Qhihi    highest doubles of Q, Q is orthogonal and transpose(Q)*A = R;
 *   Qlohi    second highest doubles of Q;
 *   Qhilo    second lowest doubles of Q;
 *   Qlolo    lowest doubles of Q;
 *   Rhihi    highest doubles of R, the reduced upper triangular form of A;
 *   Rlohi    second highest doubles of R;
 *   Rhilo    second lowest doubles of R;
 *   Rlolo    lowest doubles of R;
 *   lapsec   elapsed time in seconds. */

void CPU_cmplx4_blocked_houseqr
 ( int nrows, int ncols, int szt, int nbt,
   double **Arehihi, double **Arelohi, double **Arehilo, double **Arelolo,
   double **Aimhihi, double **Aimlohi, double **Aimhilo, double **Aimlolo,
   double **Qrehihi, double **Qrelohi, double **Qrehilo, double **Qrelolo,
   double **Qimhihi, double **Qimlohi, double **Qimhilo, double **Qimlolo,
   double **Rrehihi, double **Rrelohi, double **Rrehilo, double **Rrelolo,
   double **Rimhihi, double **Rimlohi, double **Rimhilo, double **Rimlolo,
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
 *   Arehihi  highest doubles of the real parts of A,
 *            stored as nrows arrays of ncols numbers;
 *   Arelohi  second highest doubles of the real parts of A;
 *   Arehilo  second lowest doubles of the real parts of A;
 *   Arelolo  lowest doubles of the real parts of A;
 *   Aimhihi  highest doubles of the imaginary parts of A;
 *   Aimlohi  second highest doubles of the imaginary parts of A;
 *   Aimhilo  second lowest doubles of the imaginary parts of A;
 *   Aimlolo  lowest doubles of the imaginary parts of A;
 *   Qrehihi  space for an nrows-by-nrows matrix;
 *   Qrelohi  space for an nrows-by-nrows matrix;
 *   Qrehilo  space for an nrows-by-nrows matrix;
 *   Qrelolo  space for an nrows-by-nrows matrix;
 *   Qimhihi  space for an nrows-by-nrows matrix;
 *   Qimlohi  space for an nrows-by-nrows matrix;
 *   Qimhilo  space for an nrows-by-nrows matrix;
 *   Qimlolo  space for an nrows-by-nrows matrix;
 *   Rrehihi  space for an nrows-by-ncols matrix;
 *   Rrelohi  space for an nrows-by-ncols matrix;
 *   Rrehilo  space for an nrows-by-ncols matrix;
 *   Rrelolo  space for an nrows-by-ncols matrix;
 *   Rimhihi  space for an nrows-by-ncols matrix;
 *   Rimlohi  space for an nrows-by-ncols matrix;
 *   Rimhilo  space for an nrows-by-ncols matrix;
 *   Rimlolo  space for an nrows-by-ncols matrix;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   Qrehihi  highest doubles of the real parts of Q, Q is orthogonal
 *            and transpose(Q)*A = R;
 *   Qrelohi  second highest doubles of the real parts of Q;
 *   Qrehilo  second lowest doubles of the real parts of Q;
 *   Qrelolo  lowest doubles of the real parts of Q;
 *   Qimhihi  highest doubles of the imaginary parts of Q;
 *   Qimlohi  second highest doubles of the imaginary parts of Q;
 *   Qimhilo  second lowest doubles of the imaginary parts of Q;
 *   Qimlolo  lowest doubles of the imaginary parts of Q;
 *   Rrehihi  highest doubles of the real parts of R,
 *            the reduced upper triangular form of A;
 *   Rrelohi  second highest doubles of the real parts of R;
 *   Rrehilo  second lowest doubles of the real parts of R;
 *   Rrelolo  lowest doubles of the real parts of R;
 *   Rimhihi  highest doubles of the imaginary parts of R;
 *   Rimlohi  second highest doubles of the imaginary parts of R;
 *   Rimhilo  second lowest doubles of the imaginary parts of R;
 *   Rimlolo  lowest doubles of the imaginary parts of R;
 *   lapsec   elapsed time in seconds. */

#endif
