/* The file dbl_factorizations.h specifies functions to factor matrices
 * in double precision. */

#ifndef __dbl_factorizations_h__
#define __dbl_factorizations_h__

void CPU_dbl_factors_matmatmul
 ( int rows, int dim, int cols, double **A, double **B, double **C );
/*
 * DESCRIPTION :
 *   Computes the product C of the matrix A with B on real data.
 *
 * ON ENTRY :
 *   rows     the number of rows in the matrices A and C;
 *   dim      the number of columns in A and rows in B;
 *   cols     the number of columns in the matrices B and C;
 *   A        matrix of dimensions rows and dim;
 *   B        matrix of dimensions dim and cols;
 *   C        space allocated for a rows-by-cols matrix.
 *
 * ON RETURN :
 *   C        product of A with B. */

void CPU_cmplx_factors_matmatmul
 ( int rows, int dim, int cols, double **Are, double **Aim,
   double **Bre, double **Bim, double **Cre, double **Cim );
/*
 * DESCRIPTION :
 *   Computes the product C of the matrix A with B on complex data.
 *
 * ON ENTRY :
 *   rows     the number of rows in the matrices A and C;
 *   dim      the number of columns in A and rows in B;
 *   cols     the number of columns in the matrices B and C;
 *   Are      real parts of a matrix of dimensions rows and dim;
 *   Aim      imaginary parts of a matrix of dimensions rows and dim;
 *   Bre      real parts of a matrix of dimensions dim and cols;
 *   Bim      imaginary parts of a  matrix of dimensions dim and cols;
 *   Cre      space allocated for a rows-by-cols matrix;
 *   Cim      space allocated for a rows-by-cols matrix.
 *
 * ON RETURN :
 *   Cre      real parts of the product of A with B;
 *   Cim      imaginary parts of product of A with B. */

void CPU_dbl_factors_forward ( int dim, double **L, double *b, double *x );
/*
 * DESCRIPTION :
 *   Solves the lower triangular system L*x = b with forward substitution,
 *   on real data.
 *
 * REQUIRED : the matrix L has ones on the diagonal.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in the matrix L;
 *   L        lower triangular matrix;
 *   b        right hand side vector;
 *   x        space for dim doubles.
 *
 * ON RETURN :
 *   x        the solution to L*x = b. */

void CPU_cmplx_factors_forward
 ( int dim, double **Lre, double **Lim, double *bre, double *bim,
   double *xre, double *xim );
/*
 * DESCRIPTION :
 *   Solves the lower triangular system L*x = b with forward substitution,
 *   on complex data.
 *
 * REQUIRED : the matrix L has ones on the diagonal.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in the matrix L;
 *   Lre      real parts of a lower triangular matrix;
 *   Lim      imaginary parts of a lower triangular matrix;
 *   bre      real parts of the right hand side vector;
 *   bim      imaginary parts of the right hand side vector;
 *   xre      space for dim doubles;
 *   xim      space for dim doubles.
 *
 * ON RETURN :
 *   xre      real parts of the solution to L*x = b;
 *   xim      imaginary parts the solution to L*x = b. */

void CPU_dbl_factors_backward ( int dim, double **U, double *b, double *x );
/*
 * DESCRIPTION :
 *   Solves the upper triangular system U*x = b with back substitution,
 *   on real data.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in the matrix U;
 *   U        upper triangular matrix;
 *   b        right hand side vector;
 *   x        space for dim doubles.
 *
 * ON RETURN :
 *   x        the solution to U*x = b. */

void CPU_cmplx_factors_backward
 ( int dim, double **Ure, double **Uim, double *bre, double *bim,
   double *xre, double *xim );
/*
 * DESCRIPTION :
 *   Solves the upper triangular system U*x = b with back substitution,
 *   on real data.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in the matrix U;
 *   Ure      real parts of an upper triangular matrix;
 *   Uim      imaginary parts of upper triangular matrix;
 *   bre      real parts of the right hand side vector;
 *   bim      imaginary parts of the right hand side vector;
 *   xre      space for dim doubles;
 *   xim      space for dim doubles.
 *
 * ON RETURN :
 *   xre      real parts of the solution to U*x = b;
 *   xim      imaginary parts of the solution to U*x = b. */

void CPU_dbl_factors_lufac ( int dim, double **A, int *pivots );
/*
 * DESCRIPTION :
 *   Does an inplace LU factorization with pivoting on the matrix A,
 *   on real data.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in the matrix A;
 *   A        matrix of dimension dim;
 *   pivots   space for dim pivots.
 *
 * ON RETURN :
 *   A        the lower triangular part of A contains the multipliers
 *            and the upper triangular part of A the row reduced A;
 *   pivots   are the pivots used. */

void CPU_cmplx_factors_lufac
 ( int dim, double **Are, double **Aim, int *pivots );
/*
 * DESCRIPTION :
 *   Does an inplace LU factorization with pivoting on the matrix A,
 *   on complex data.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in the matrix A;
 *   Are      real parts of a matrix of dimension dim;
 *   Aim      imaginary parts of a matrix of dimension dim;
 *   pivots   space for dim pivots.
 *
 * ON RETURN :
 *   Are      the lower triangular part of Are contains the real parts
 *            of the multipliers and the upper triangular part of Are
 *            the real parts of the row reduced A;
 *   Aim      the lower triangular part of Are contains the imaginary parts
 *            of the multipliers and the upper triangular part of Are
 *            the imaginary parts of the row reduced A;
 *   pivots   are the pivots used. */

void CPU_dbl_factors_lusolve
 ( int dim, double **A, int *pivots, double *b, double *x );
/*
 * DESCRIPTION :
 *   Does an inplace LU factorization with pivoting on the matrix A,
 *   to solve the system A*x = b.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in the matrix A;
 *   A        matrix of dimension dim;
 *   pivots   space for dim pivots;
 *   b        right hand side vector;
 *   x        space for dim doubles.
 *
 * ON RETURN :
 *   A        the lower triangular part of A contains the multipliers
 *            and the upper triangular part of A the row reduced A;
 *   b        used as work space;
 *   x        the solution to A*x = b. */

void CPU_cmplx_factors_lusolve
 ( int dim, double **Are, double **Aim, int *pivots,
   double *bre, double *bim, double *xre, double *xim );
/*
 * DESCRIPTION :
 *   Does an inplace LU factorization with pivoting on the matrix A,
 *   to solve the system A*x = b.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in the matrix A;
 *   Are      real parts of a matrix of dimension dim;
 *   Aim      imaginary parts of a  matrix of dimension dim;
 *   pivots   space for dim pivots;
 *   bre      real parts of the right hand side vector;
 *   bim      imaginary parts of the right hand side vector;
 *   xre      space for dim doubles;
 *   xim      space for dim doubles.
 *
 * ON RETURN :
 *   Are      the lower triangular part of Are contains the real parts
 *            of the multipliers and the upper triangular part of Are
 *            the real parts of the row reduced A;
 *   Aim      the lower triangular part of Are contains the imaginary parts
 *            of the multipliers and the upper triangular part of Are
 *            the imaginary parts of the row reduced A;
 *   bre      used as work space;
 *   bim      used as work space;
 *   xre      real parts of the solution to A*x = b;
 *   xim      imaginary parts of the solution to A*x = b. */

void CPU_dbl_factors_house ( int n, double *x, double *v, double *beta );
/*
 * DESCRIPTION :
 *   Computes the Householder vector of an n-dimensional vector x.
 *
 * ON ENTRY :
 *   n        dimension of the vector x;
 *   x        n doubles;
 *   v        space for n doubles.
 *
 * ON RETURN :
 *   v        the Householder vector;
 *   beta     equals 2/(transpose(v)*v). */

void CPU_cmplx_factors_house 
 ( int n, double *xre, double *xim, double *vre, double *vim, double *beta );
/*
 * DESCRIPTION :
 *   Computes the Householder vector of an n-dimensional vector x.
 *
 * ON ENTRY :
 *   n        dimension of the vector x;
 *   xre      real parts of the vector x;
 *   xim      imaginary parts of the vector x;
 *   vre      space for n doubles;
 *   vim      space for n doubles.
 *
 * ON RETURN :
 *   vre      real parts of the Householder vector;
 *   vim      imaginary parts of the Householder vector;
 *   beta     equals 2/(transpose(v)*v). */

void CPU_dbl_factors_leftRupdate
 ( int nrows, int ncols, int k, double **R, double *v, double beta,
   bool verbose=false );
/*
 * DESCRIPTION :
 *   Applies the Householder matrix to R.
 *
 * ON ENTRY :
 *   nrows    number of rows in the matrix R;
 *   ncols    number of columns in the matrix R;
 *   k        current column index in R;
 *   R        an nrows-by-ncols matrix;
 *   v        the Householder vector;
 *   beta     the beta computed by CPU_dbl_house;
 *   verbose  is the verbose flag.
 *
 * ON RETURN :
 *   R        update with the Householder matrix. */

void CPU_cmplx_factors_leftRupdate
 ( int nrows, int ncols, int k, double **Rre, double **Rim,
   double *vre, double *vim, double beta );
/*
 * DESCRIPTION :
 *   Applies the Householder matrix to R.
 *
 * ON ENTRY :
 *   nrows    number of rows in the matrix R;
 *   ncols    number of columns in the matrix R;
 *   k        current column index in R;
 *   Rre      real parts of an nrows-by-ncols matrix;
 *   Rim      imaginary parts of an nrows-by-ncols matrix;
 *   vre      real parts of the Householder vector;
 *   vim      imaginary parts of the Householder vector;
 *   beta     the beta computed by CPU_dbl_house.
 *
 * ON RETURN :
 *   Rre      real parts of the update with the Householder matrix;
 *   Rim      imaginary parts of the update with the Householder matrix. */

void CPU_dbl_factors_rightQupdate
 ( int n, int k, double **Q, double *v, double beta );
/*
 * DESCRIPTION :
 *   Applies the Householder matrix to Q.
 *
 * ON ENTRY :
 *   n        dimension of the matrix Q;
 *   k        current column index in Q;
 *   Q        an n-by-n matrix;
 *   v        the Householder vector;
 *   beta     the beta computed by CPU_dbl_house.
 *
 * ON RETURN :
 *   Q        update with the Householder matrix. */

void CPU_cmplx_factors_rightQupdate
 ( int n, int k, double **Qre, double **Qim,
   double *vre, double *vim, double beta );
/*
 * DESCRIPTION :
 *   Applies the Householder matrix to Q.
 *
 * ON ENTRY :
 *   n        dimension of the matrix Q;
 *   k        current column index in Q;
 *   Qre      real parts of an n-by-n matrix;
 *   Qim      imaginary parts of an n-by-n matrix;
 *   vre      real parts of the Householder vector;
 *   vim      imaginary parts the Householder vector;
 *   beta     the beta computed by CPU_dbl_house.
 *
 * ON RETURN :
 *   Qre      real parts of the update with the Householder matrix;
 *   Qim      imaginary parts of the update with the Householder matrix. */

void CPU_dbl_factors_houseqr
 ( int nrows, int ncols, double **A, double **Q, double **R );
/*
 * DESCRIPTION :
 *   Applies Householder matrices to compute a QR decomposition of A.
 *
 * REQUIRED : nrows >= ncols.
 *
 * ON ENTRY :
 *   nrows    number of rows of A;
 *   ncols    number of columns of A;
 *   A        an nrows-by-ncols matrix,
 *            stored as nrows arrays of ncols numbers;
 *   Q        space for an nrows-by-nrows matrix;
 *   R        space for an nrows-by-ncols matrix.
 *
 * ON RETURN :
 *   Q        an orthogonal matrix, transpose(Q)*A = R;
 *   R        the reduced upper triangular form of A. */

void CPU_cmplx_factors_houseqr
 ( int nrows, int ncols, double **Are, double **Aim,
   double **Qre, double **Qim, double **Rre, double **Rim );
/*
 * DESCRIPTION :
 *   Applies Householder matrices to compute a QR decomposition of A.
 *
 * REQUIRED : nrows >= ncols.
 *
 * ON ENTRY :
 *   nrows    number of rows of A;
 *   ncols    number of columns of A;
 *   Are      real parts of an nrows-by-ncols matrix,
 *            stored as nrows arrays of ncols numbers;
 *   Aim      imaginary parts of an nrows-by-ncols matrix,
 *            stored as nrows arrays of ncols numbers;
 *   Qre      space for an nrows-by-nrows matrix;
 *   Qim      space for an nrows-by-nrows matrix;
 *   Rre      space for an nrows-by-ncols matrix;
 *   Rim      space for an nrows-by-ncols matrix.
 *
 * ON RETURN :
 *   Qre      real parts of an orthogonal matrix, transpose(Q)*A = R;
 *   Qim      imaginary parts of an orthogonal matrix, transpose(Q)*A = R;
 *   Rre      real parts of the reduced upper triangular form of A;
 *   Rim      imaginary parts of the reduced upper triangular form of A. */

void CPU_dbl_factors_qrbs
 ( int nrows, int ncols, double **Q, double **R,
   double *rhs, double *sol, double *wrkvec ); 
/*
 * DESCRIPTION :
 *   Applies back substitution for the least squares solution
 *   with the QR factorization, on real data.
 *
 * REQUIRED : nrows >= ncols;
 *
 * ON ENTRY :
 *   nrows    number of rows of R, dimension of Q;
 *   ncols    number of columns of R;
 *   Q        the Q of the QR factorization;
 *   R        the R of the QR factorization;
 *   rhs      right hand side vector, of size nrows;
 *   sol      space for the solution, of size ncols;
 *   wrkvec   work space for nrows elements.
 *
 * ON RETURN :
 *   sol      the least squares solution;
 *   wrkvec   contains the product of Q transpose with rhs. */

void CPU_cmplx_factors_qrbs
 ( int nrows, int ncols,
   double **Qre, double **Qim, double **Rre, double **Rim,
   double *rhsre, double *rhsim, double *solre, double *solim,
   double *wrkvecre, double *wrkvecim );
/*
 * DESCRIPTION :
 *   Applies back substitution for the least squares solution
 *   with the QR factorization, on complex data.
 *
 * REQUIRED : nrows >= ncols;
 *
 * ON ENTRY :
 *   nrows    number of rows of R, dimension of Q;
 *   ncols    number of columns of R;
 *   Qre      the real parts of Q of the QR factorization;
 *   Qim      the imaginary parts of Q of the QR factorization;
 *   Rre      the real parts of R of the QR factorization;
 *   Rim      the imaginary parts of R of the QR factorization;
 *   rhsre    real parts of right hand side, of size nrows;
 *   rhsim    imaginary parts of right hand side, of size nrows;
 *   solre    space for the real parts of solution, of size ncols;
 *   solim    space for the imaginary parts of solution, of size ncols;
 *   wrkvecre is work space for nrows elements;
 *   wrkvecim is work space for nrows elements.
 *
 * ON RETURN :
 *   solre    real parts of the least squares solution;
 *   solim    imaginary parts of the least squares solution;
 *   wrkvecre has the real parts of product of Q^H with rhs;
 *   wrkvecim has the imaginary parts of product of Q^H with rhs. */

#endif
