/* The file dbl2_factorizations.h specifies functions to factor matrices
 * in double double precision. */

#ifndef __dbl2_factorizations_h__
#define __dbl2_factorizations_h__

void CPU_dbl2_factors_matmatmul
 ( int rows, int dim, int cols, double **Ahi, double **Alo,
   double **Bhi, double **Blo, double **Chi, double **Clo );
/*
 * DESCRIPTION :
 *   Computes the product C of the matrix A with B on real data.
 *
 * ON ENTRY :
 *   rows     the number of rows in the matrices A and C;
 *   dim      the number of columns in A and rows in B;
 *   cols     the number of columns in the matrices B and C;
 *   Ahi      high doubles of a matrix of dimensions rows and dim;
 *   Alo      low doubles of a matrix of dimensions rows and dim;
 *   Bhi      high doubles of a matrix of dimensions dim and cols;
 *   Blo      low doubles of a matrix of dimensions dim and cols;
 *   Chi      space allocated for a rows-by-cols matrix;
 *   Clo      space allocated for a rows-by-cols matrix.
 *
 * ON RETURN :
 *   Chi      high doubles of the product of A with B;
 *   Clo      low doubles of the product of A with B. */

void CPU_cmplx2_factors_matmatmul
 ( int rows, int dim, int cols,
   double **Arehi, double **Arelo, double **Aimhi, double **Aimlo,
   double **Brehi, double **Brelo, double **Bimhi, double **Bimlo,
   double **Crehi, double **Crelo, double **Cimhi, double **Cimlo );
/*
 * DESCRIPTION :
 *   Computes the product C of the matrix A with B on complex data.
 *
 * ON ENTRY :
 *   rows     the number of rows in the matrices A and C;
 *   dim      the number of columns in A and rows in B;
 *   cols     the number of columns in the matrices B and C;
 *   Are      high doubles of the real parts of A;
 *   Are      low doubles of the real parts of A;
 *   Aim      high doubles of the imaginary parts of A; 
 *   Aim      low doubles of the imaginary parts of A;
 *   Brehi    high doubles of the real parts of B;
 *   Brelo    low doubles of the real parts of B;
 *   Bimhi    high doubles of the imaginary parts of B;
 *   Bimlo    low doubles of the imaginary parts of B;
 *   Crehi    space allocated for a rows-by-cols matrix;
 *   Crelo    space allocated for a rows-by-cols matrix;
 *   Cimhi    space allocated for a rows-by-cols matrix.
 *   Cimlo    space allocated for a rows-by-cols matrix.
 *
 * ON RETURN :
 *   Crehi    high doubles of the real parts of the product of A with B;
 *   Crelo    low doubles of the real parts of the product of A with B;
 *   Cimhi    high doubles of the imaginary parts of product of A with B;
 *   Cimlo    low doubles of the imaginary parts of product of A with B. */

void CPU_dbl2_factors_forward
 ( int dim, double **Lhi, double **Llo, double *bhi, double *blo,
   double *xhi, double *xlo );
/*
 * DESCRIPTION :
 *   Solves the lower triangular system L*x = b with forward substitution,
 *   on real data.
 *
 * REQUIRED : the matrix L has ones on the diagonal.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in the matrix L;
 *   Lhi      high doubles of a lower triangular matrix;
 *   Llo      low doubles of a lower triangular matrix;
 *   bhi      high doubles of a right hand side vector;
 *   blo      low doubles of a right hand side vector;
 *   xhi      space for dim doubles;
 *   xlo      space for dim doubles.
 *
 * ON RETURN :
 *   xhi      high doubles of the solution to L*x = b;
 *   xlo      low doubles of the solution to L*x = b. */

void CPU_cmplx2_factors_forward
 ( int dim, double **Lrehi, double **Lrelo, double **Limhi, double **Limlo,
   double *brehi, double *brelo, double *bimhi, double *bimlo,
   double *xrehi, double *xrelo, double *ximhi, double *ximlo );
/*
 * DESCRIPTION :
 *   Solves the lower triangular system L*x = b with forward substitution,
 *   on real data.
 *
 * REQUIRED : the matrix L has ones on the diagonal.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in the matrix L;
 *   Lrehi    high doubles of the real parts of a lower triangular matrix;
 *   Lrelo    low doubles of the real parts of a lower triangular matrix;
 *   Limhi    high doubles of the imaginary parts of a lower triangular matrix;
 *   Limlo    low doubles of the imaginary parts of a lower triangular matrix;
 *   brehi    high doubles of the real parts of a right hand side vector;
 *   brelo    low doubles of the real parts of a right hand side vector;
 *   bimhi    high doubles of the imaginary parts of a right hand side vector;
 *   bimlo    low doubles of the imaginary parts of a right hand side vector;
 *   xrehi    space for dim doubles;
 *   xrelo    space for dim doubles;
 *   ximhi    space for dim doubles;
 *   ximlo    space for dim doubles.
 *
 * ON RETURN :
 *   xrehi    high doubles of the real parts of the solution to L*x = b;
 *   xrelo    low doubles of the real parts of the solution to L*x = b;
 *   ximhi    high doubles of the imaginary parts of the solution to L*x = b;
 *   ximlo    low doubles of the imaginary parts of the solution to L*x = b. */

void CPU_dbl2_factors_backward
 ( int dim, double **Uhi, double **Ulo, double *bhi, double *blo,
   double *xhi, double *xlo );
/*
 * DESCRIPTION :
 *   Solves the upper triangular system U*x = b with back substitution,
 *   on real data.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in the matrix U;
 *   Uhi      high doubles of an upper triangular matrix;
 *   Ulo      low doubles of an upper triangular matrix;
 *   bhi      high doubles of a right hand side vector;
 *   blo      low doubles of a right hand side vector;
 *   xhi      space for dim doubles;
 *   xlo      space for dim doubles.
 *
 * ON RETURN :
 *   xhi      high doubles of the solution to U*x = b;
 *   xlo      low doubles of the solution to U*x = b. */

void CPU_cmplx2_factors_backward
 ( int dim, double **Urehi, double **Urelo, double **Uimhi, double **Uimlo,
   double *brehi, double *brelo, double *bimhi, double *bimlo,
   double *xrehi, double *xrelo, double *ximhi, double *ximlo );
/*
 * DESCRIPTION :
 *   Solves the upper triangular system U*x = b with back substitution,
 *   on real data.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in the upper triangular matrix U;
 *   Urehi    high doubles of the real parts of U;
 *   Urelo    low doubles of the real parts of U;
 *   Uimhi    high doubles of the imaginary parts of U;
 *   Uimlo    low doubles of the imaginary parts of U;
 *   brehi    high doubles of the real parts of a right hand side vector;
 *   brelo    low doubles of the real parts of a right hand side vector;
 *   bimhi    high doubles of the imaginary parts of a right hand side vector;
 *   bimlo    low doubles of the imaginary parts of a right hand side vector;
 *   xrehi    space for dim doubles;
 *   xrelo    space for dim doubles;
 *   ximhi    space for dim doubles;
 *   ximlo    space for dim doubles.
 *
 * ON RETURN :
 *   xrehi    high doubles of the real parts of the solution to U*x = b;
 *   xrelo    low doubles of the imaginary parts of the solution to U*x = b;
 *   ximhi    high doubles of the real parts of the solution to U*x = b;
 *   ximlo    low doubles of the imaginary parts of the solution to U*x = b. */

void CPU_dbl2_factors_lufac
 ( int dim, double **Ahi, double **Alo, int *pivots );
/*
 * DESCRIPTION :
 *   Does an inplace LU factorization with pivoting on the matrix A,
 *   on real data.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in the matrix A;
 *   Ahi      high doubles of a matrix of dimension dim;
 *   Alo      low doubles of a matrix of dimension dim;
 *   pivots   space for dim pivots.
 *
 * ON RETURN :
 *   Ahi      the lower triangular part of Ahi contains the high doubles
 *            of the multipliers and the upper triangular part of Ahi
 *            the high doubles of the row reduced A;
 *   Alo      the lower triangular part of Ahi contains the low doubles
 *            of the multipliers and the upper triangular part of Ahi
 *            the low doubles of the row reduced A;
 *   pivots   are the pivots used. */

void CPU_cmplx2_factors_lufac
 ( int dim, double **Arehi, double **Arelo, double **Aimhi, double **Aimlo,
   int *pivots );
/*
 * DESCRIPTION :
 *   Does an inplace LU factorization with pivoting on the matrix A,
 *   on real data.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in the matrix A;
 *   Arehi    high doubles of the real parts of A;
 *   Arelo    low doubles of the real parts of A;
 *   Aimhi    high doubles of the imaginary parts of A;
 *   Aimlo    low doubles of the imaginary parts of A;
 *   pivots   space for dim pivots.
 *
 * ON RETURN :
 *   Arehi    high doubles of the real parts
 *            of the multipliers and the row reduced A;
 *   Arelo    low doubles of the real parts
 *            of the multipliers and the row reduced A;
 *   Aimhi    high doubles of the imaginary parts
 *            of the multipliers and the row reduced A;
 *   Aimlo    low doubles of the imaginary parts
 *            of the multipliers and the row reduced A;
 *   pivots   are the pivots used. */

void CPU_dbl2_factors_lusolve
 ( int dim, double **Ahi, double **Alo, int *pivots,
   double *bhi, double *blo, double *xhi, double *xlo );
/*
 * DESCRIPTION :
 *   Does an inplace LU factorization with pivoting on the matrix A,
 *   to solve the system A*x = b, on real data.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in the matrix A;
 *   Ahi      high doubles of a matrix of dimension dim;
 *   Alo      low doubles of a matrix of dimension dim;
 *   pivots   space for dim pivots;
 *   bhi      high doubles of the right hand side vector;
 *   blo      low doubles of the right hand side vector;
 *   xhi      space for dim doubles;
 *   xlo      space for dim doubles.
 *
 * ON RETURN :
 *   Ahi      the lower triangular part of Ahi contains the high doubles
 *            of the multipliers and the upper triangular part of Ahi
 *            the high doubles of the row reduced A;
 *   Alo      the lower triangular part of Ahi contains the low doubles
 *            of the multipliers and the upper triangular part of Ahi
 *            the low doubles of the row reduced A;
 *   bhi      used as work space;
 *   blo      used as work space;
 *   xhi      high doubles of the solution to A*x = b;
 *   xlo      low doubles of the solution to A*x = b. */

void CPU_cmplx2_factors_lusolve
 ( int dim, double **Arehi, double **Arelo, double **Aimhi, double **Aimlo,
   int *pivots, double *brehi, double *brelo, double *bimhi, double *bimlo,
   double *xrehi, double *xrelo, double *ximhi, double *ximlo );
/*
 * DESCRIPTION :
 *   Does an inplace LU factorization with pivoting on the matrix A,
 *   to solve the system A*x = b, on complex data.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in the matrix A;
 *   Arehi    high doubles of the real parts of A;
 *   Arelo    low doubles of the imaginary parts of A;
 *   Aimhi    high doubles of the real parts of A;
 *   Aimlo    low doubles of the imaginary parts of A;
 *   pivots   space for dim pivots;
 *   brehi    high doubles of the real parts of b;
 *   brelo    low doubles of the real parts of b;
 *   bimhi    high doubles of the imaginary parts of b;
 *   bimlo    low doubles of the imaginary parts of b;
 *   xrehi    space for dim doubles;
 *   xrelo    space for dim doubles;
 *   ximhi    space for dim doubles;
 *   ximlo    space for dim doubles.
 *
 * ON RETURN :
 *   Arehi    high doubles of the real parts
 *            of the multipliers and the row reduced A;
 *   Arelo    low doubles of the real parts
 *            of the multipliers and the row reduced A;
 *   Aimhi    high doubles of the imaginary parts
 *            of the multipliers and the row reduced A;
 *   Aimlo    low doubles of the imaginary parts
 *            of the multipliers and the row reduced A;
 *   brehi    used as work space;
 *   brelo    used as work space;
 *   bimhi    used as work space;
 *   bimlo    used as work space;
 *   xrehi    high doubles of the real parts of the solution x;
 *   xrelo    low doubles of the real parts of the solution x;
 *   ximhi    high doubles of the imaginary parts of the solution x;
 *   ximlo    low doubles of the imaginary parts of the solution x. */

void CPU_dbl2_factors_house
 ( int n, double *xhi, double *xlo, double *vhi, double *vlo,
   double *betahi, double *betalo );
/*
 * DESCRIPTION :
 *   Computes the Householder vector of an n-dimensional vector x.
 *
 * ON ENTRY :
 *   n        dimension of the vector x;
 *   xhi      the n high doubles of x;
 *   xlo      the n low doubles of x;
 *   vhi      space for n doubles;
 *   vhi      space for n doubles.
 *
 * ON RETURN :
 *   vhi      high doubles of the Householder vector;
 *   vlo      low doubles of the Householder vector;
 *   betahi   equals the high double of 2/(transpose(v)*v);
 *   betalo   equals the high double of 2/(transpose(v)*v). */

void CPU_cmplx2_factors_house 
( int n, double *xrehi, double *xrelo, double *ximhi, double *ximlo,
         double *vrehi, double *vrelo, double *vimhi, double *vimlo,
  double *betahi, double *betalo );
/*
 * DESCRIPTION :
 *   Computes the Householder vector of an n-dimensional vector x.
 *
 * ON ENTRY :
 *   n        dimension of the vector x;
 *   xrehi    high doubles of the real parts of the vector x;
 *   xrelo    low doubles of the real parts of the vector x;
 *   ximhi    high doubles of the imaginary parts of the vector x;
 *   ximlo    low doubles of the imaginary parts of the vector x;
 *   vrehi    space for n doubles;
 *   vrelo    space for n doubles;
 *   vimhi    space for n doubles;
 *   vimlo    space for n doubles.
 *
 * ON RETURN :
 *   vrehi    high doubles of the real parts of the Householder vector;
 *   vrelo    low doubles of the real parts of the Householder vector;
 *   vimhi    high doubles of the imaginary parts of the Householder vector;
 *   vimlo    low doubles of the imaginary parts of the Householder vector;
 *   betahi   high double of 2/(transpose(v)*v);
 *   betalo   low double of 2/(transpose(v)*v). */

void CPU_dbl2_factors_leftRupdate
 ( int nrows, int ncols, int k, double **Rhi, double **Rlo,
   double *vhi, double *vlo, double betahi, double betalo );
/*
 * DESCRIPTION :
 *   Applies the Householder matrix to R.
 *
 * ON ENTRY :
 *   nrows    number of rows in the matrix R;
 *   ncols    number of columns in the matrix R;
 *   k        current column index in R;
 *   Rhi      high doubles of an nrows-by-ncols matrix;
 *   Rlo      low doubles of an nrows-by-ncols matrix;
 *   vhi      high doubles of the Householder vector;
 *   vlo      low doubles of the Householder vector;
 *   betahi   the betahi computed by CPU_dbl_house;
 *   betalo   the betalo computed by CPU_dbl_house.
 *
 * ON RETURN :
 *   Rhi      high doubles of the update with the Householder matrix;
 *   Rlo      low doubles of the update with the Householder matrix. */

void CPU_cmplx2_factors_leftRupdate
 ( int nrows, int ncols, int k,
   double **Rrehi, double **Rrelo, double **Rimhi, double **Rimlo,
   double *vrehi, double *vrelo, double *vimhi, double *vimlo,
   double betahi, double betalo );
/*
 * DESCRIPTION :
 *   Applies the Householder matrix to R.
 *
 * ON ENTRY :
 *   nrows    number of rows in the matrix R;
 *   ncols    number of columns in the matrix R;
 *   k        current column index in R;
 *   Rrehi    high doubles of the real parts of an nrows-by-ncols matrix;
 *   Rrelo    low doubles of the real parts of an nrows-by-ncols matrix;
 *   Rimhi    high doubles of the imaginary parts of an nrows-by-ncols matrix;
 *   Rimlo    low doubles of the imaginary parts of an nrows-by-ncols matrix;
 *   vrehi    high doubles of the real parts of the Householder vector;
 *   vrelo    low doubles of the real parts of the Householder vector;
 *   vimhi    high doubles of the imaginary parts of the Householder vector;
 *   vimlo    low doubles of the imaginary parts of the Householder vector;
 *   betahi   high double of the the beta computed by CPU_dbl_house.
 *   betalo   low double of the the beta computed by CPU_dbl_house.
 *
 * ON RETURN :
 *   Rrehi    high doubles of the real parts of the update
 *            with the Householder matrix;
 *   Rrelo    low doubles of the real parts of the update
 *            with the Householder matrix;
 *   Rimhi    high doubles of the imaginary parts of the update
 *            with the Householder matrix;
 *   Rimlo    low doubles of the imaginary parts of the update
 *            with the Householder matrix. */

void CPU_dbl2_factors_rightQupdate
 ( int n, int k, double **Qhi, double **Qlo,
   double *vhi, double *vlo, double betahi, double betalo );
/*
 * DESCRIPTION :
 *   Applies the Householder matrix to Q.
 *
 * ON ENTRY :
 *   n        dimension of the matrix Q;
 *   k        current column index in Q;
 *   Qhi      high doubles of an n-by-n matrix;
 *   Qlo      low doubles of an n-by-n matrix;
 *   vhi      high doubles of the Householder vector;
 *   vlo      low doubles of the Householder vector;
 *   betahi   the betahi computed by CPU_dbl_house;
 *   betalo   the betalo computed by CPU_dbl_house.
 *
 * ON RETURN :
 *   Qhi      high doubles of the update with the Householder matrix;
 *   Qlo      low doubles of the update with the Householder matrix. */

void CPU_cmplx2_factors_rightQupdate
 ( int n, int k,
   double **Qrehi, double **Qrelo, double **Qimhi, double **Qimlo,
   double *vrehi, double *vrelo, double *vimhi, double *vimlo,
   double betahi, double betalo );
/*
 * DESCRIPTION :
 *   Applies the Householder matrix to Q.
 *
 * ON ENTRY :
 *   n        dimension of the matrix Q;
 *   k        current column index in Q;
 *   Qrehi    high doubles of the real parts of an n-by-n matrix;
 *   Qrelo    low doubles of the real parts of an n-by-n matrix;
 *   Qimhi    high doubles of the imaginary parts of an n-by-n matrix;
 *   Qimlo    low doubles of the imaginary parts of an n-by-n matrix;
 *   vrehi    high doubles of the real parts of the Householder vector;
 *   vrelo    low doubles of the real parts of the Householder vector;
 *   vimhi    high doubles of the imaginary parts the Householder vector;
 *   vimlo    low doubles of the imaginary parts the Householder vector;
 *   betahi   high double of the the beta computed by CPU_dbl_house;
 *   betalo   low double of the the beta computed by CPU_dbl_house.
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

void CPU_dbl2_factors_houseqr
 ( int nrows, int ncols, double **Ahi, double **Alo,
   double **Qhi, double **Qlo, double **Rhi, double **Rlo );
/*
 * DESCRIPTION :
 *   Applies Householder matrices to compute a QR decomposition of A.
 *
 * REQUIRED : nrows >= ncols.
 *
 * ON ENTRY :
 *   nrows    number of rows of A;
 *   ncols    number of columns of A;
 *   Ahi      high doubles of an nrows-by-ncols matrix,
 *            stored as nrows arrays of ncols numbers;
 *   Alo      low doubles of an nrows-by-ncols matrix,
 *            stored as nrows arrays of ncols numbers;
 *   Qhi      space for an nrows-by-nrows matrix;
 *   Qlo      space for an nrows-by-nrows matrix;
 *   Rhi      space for an nrows-by-ncols matrix;
 *   Rlo      space for an nrows-by-ncols matrix.
 *
 * ON RETURN :
 *   Qhi      high doubles of an orthogonal matrix, transpose(Q)*A = R;
 *   Qlo      low doubles of an orthogonal matrix, transpose(Q)*A = R;
 *   Rhi      high doubles of the reduced upper triangular form of A;
 *   Rlo      low doubles of the reduced upper triangular form of A. */

void CPU_cmplx2_factors_houseqr
 ( int nrows, int ncols,
   double **Arehi, double **Arelo, double **Aimhi, double **Aimlo,
   double **Qrehi, double **Qrelo, double **Qimhi, double **Qimlo,
   double **Rrehi, double **Rrelo, double **Rimhi, double **Rimlo );
/*
 * DESCRIPTION :
 *   Applies Householder matrices to compute a QR decomposition of A.
 *
 * REQUIRED : nrows >= ncols.
 *
 * ON ENTRY :
 *   nrows    number of rows of A;
 *   ncols    number of columns of A;
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
 *   Rimlo    space for an nrows-by-ncols matrix.
 *
 * ON RETURN :
 *   Qrehi    high doubles of the real parts of the orthogonal matrix Q
 *            that triangulates A: transpose(Q)*A = R;
 *   Qrelo    low doubles of the real parts of the orthogonal matrix Q;
 *   Qimhi    high doubles of the imaginary parts of the orthogonal matrix Q;
 *   Qimlo    low doubles of the imaginary parts of the orthogonal matrix Q;
 *   Rrehi    high doubles of the real parts of R,
 *            the reduced upper triangular form of A;
 *   Rrelo    low doubles of the real parts of R;
 *   Rimhi    high doubles of the imaginary parts of R;
 *   Rimlo    low doubles of the imaginary parts of R. */

void CPU_dbl2_factors_qrbs
 ( int nrows, int ncols,
   double **Qhi, double **Qlo, double **Rhi, double **Rlo,
   double *rhshi, double *rhslo, double *solhi, double *sollo,
   double *wrkvechi, double *wrkveclo ); 
/*
 * DESCRIPTION :
 *   Applies back substitution for the least squares solution
 *   with the QR factorization.
 *
 * REQUIRED : nrows >= ncols;
 *
 * ON ENTRY :
 *   nrows    number of rows of R, dimension of Q;
 *   ncols    number of columns of R;
 *   Qhi      high doubles of the Q of the QR factorization;
 *   Qlo      low doubles of the Q of the QR factorization;
 *   Rhi      high doubles of the R of the QR factorization;
 *   Rlo      low doubles of the R of the QR factorization;
 *   rhshi    high doubles of the right hand side vector, of size nrows;
 *   rhslo    low doubles of the right hand side vector, of size nrows;
 *   solhi    space for the high doubles of the solution, of size ncols;
 *   sollo    space for the low doubles of the solution, of size ncols;
 *   wrkvechi has work space for nrows elements.
 *   wrkveclo has work space for nrows elements.
 *
 * ON RETURN :
 *   solhi    high doubles of the least squares solution;
 *   sollo    low doubles of the least squares solution;
 *   wrkvechi contains the high doubles of the product of Q^T with rhs;
 *   wrkveclo contains the low doubles of the product of Q^T with rhs. */

void CPU_cmplx2_factors_qrbs
 ( int nrows, int ncols,
   double **Qrehi, double **Qrelo, double **Qimhi, double **Qimlo, 
   double **Rrehi, double **Rrelo, double **Rimhi, double **Rimlo,
   double *rhsrehi, double *rhsrelo, double *rhsimhi, double *rhsimlo,
   double *solrehi, double *solrelo, double *solimhi, double *solimlo,
   double *wrkvecrehi, double *wrkvecrelo,
   double *wrkvecimhi, double *wrkvecimlo );
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
 *   Qrehi    high doubles of the real parts of Q of the QR factorization;
 *   Qrelo    low doubles of the real parts of Q of the QR factorization;
 *   Qimhi    high doubles of the imaginary parts of Q of the QR;
 *   Qimlo    low doubles of the imaginary parts of Q of the QR;
 *   Rrehi    high doubles of the real parts of R of the QR;
 *   Rrelo    low doubles of the real parts of R of the QR;
 *   Rimhi    high doubles of the imaginary parts of R of the QR;
 *   Rimlo    low doubles of the imaginary parts of R of the QR;
 *   rhsrehi  high doubles of real parts of right hand side, of size nrows;
 *   rhsrelo  low doubles of real parts of right hand side, of size nrows;
 *   rhsimhi  high doubles of the imaginary parts of right hand side;
 *   rhsimlo  low doubles of the imaginary parts of right hand side;
 *   solrehi  space for the high doubles of the real parts of solution,
 *            of size ncols;
 *   solrelo  space for the low doubles of the real parts of solution;
 *   solimhi  space for the high doubles of the imaginary parts of solution;
 *   solimlo  space for the low doubles of the imaginary parts of solution;
 *   wrkvecrehi is work space for nrows elements;
 *   wrkvecrelo is work space for nrows elements;
 *   wrkvecimhi is work space for nrows elements;
 *   wrkvecimlo is work space for nrows elements.
 *
 * ON RETURN :
 *   solrehi  high doubles of the real parts of the least squares solution;
 *   solrelo  low doubles of the real parts of the least squares solution;
 *   solimhi  high doubles of the imaginary parts of the solution;
 *   solimlo  low doubles of the imaginary parts of the solution;
 *   wrkvecrehi has the high doubles of the real parts of Q^H*rhs;
 *   wrkvecrelo has the low doubles of the real parts of Q^H*rhs;
 *   wrkvecimhi has the high doubles of the imaginary parts of Q^H*rhs;
 *   wrkvecimlo has the low doubles of the imaginary parts of Q^H*rhs. */

#endif
