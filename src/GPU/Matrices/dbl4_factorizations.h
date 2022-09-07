/* The file dbl4_factorizations.h specifies functions to factor matrices
 * in quad double precision. */

#ifndef __dbl4_factorizations_h__
#define __dbl4_factorizations_h__

void CPU_dbl4_factors_matmatmul
 ( int rows, int dim, int cols,
   double **Ahihi, double **Alohi, double **Ahilo, double **Alolo,
   double **Bhihi, double **Blohi, double **Bhilo, double **Blolo,
   double **Chihi, double **Clohi, double **Chilo, double **Clolo );
/*
 * DESCRIPTION :
 *   Computes the product C of the matrix A with B on real data.
 *
 * ON ENTRY :
 *   rows     the number of rows in the matrices A and C;
 *   dim      the number of columns in A and rows in B;
 *   cols     the number of columns in the matrices B and C;
 *   Ahihi    highest doubles of a matrix of dimensions rows and dim;
 *   Alohi    second highest doubles of a matrix of dimensions rows and dim;
 *   Ahilo    second lowest doubles of a matrix of dimensions rows and dim;
 *   Alolo    lowest doubles of a matrix of dimensions rows and dim;
 *   Bhihi    highest doubles of a matrix of dimensions dim and cols;
 *   Blohi    second highest doubles of a matrix of dimensions dim and cols;
 *   Bhilo    second lowest doubles of a matrix of dimensions dim and cols;
 *   Blolo    lowest doubles of a matrix of dimensions dim and cols;
 *   Chihi    space allocated for a rows-by-cols matrix;
 *   Clohi    space allocated for a rows-by-cols matrix;
 *   Chilo    space allocated for a rows-by-cols matrix;
 *   Clolo    space allocated for a rows-by-cols matrix.
 *
 * ON RETURN :
 *   Chihi    highest doubles of the product of A with B;
 *   Clohi    second highest doubles of the product of A with B;
 *   Chilo    second lowest doubles of the product of A with B;
 *   Clolo    lowest doubles of the product of A with B. */

void CPU_cmplx4_factors_matmatmul
 ( int rows, int dim, int cols,
   double **Arehihi, double **Arelohi, double **Arehilo, double **Arelolo,
   double **Aimhihi, double **Aimlohi, double **Aimhilo, double **Aimlolo,
   double **Brehihi, double **Brelohi, double **Brehilo, double **Brelolo,
   double **Bimhihi, double **Bimlohi, double **Bimhilo, double **Bimlolo,
   double **Crehihi, double **Crelohi, double **Crehilo, double **Crelolo,
   double **Cimhihi, double **Cimlohi, double **Cimhilo, double **Cimlolo );
/*
 * DESCRIPTION :
 *   Computes the product C of the matrix A with B on complex data.
 *
 * ON ENTRY :
 *   rows     the number of rows in the matrices A and C;
 *   dim      the number of columns in A and rows in B;
 *   cols     the number of columns in the matrices B and C;
 *   Arehihi  are the highest doubles of the real parts of A;
 *   Arelohi  are the second highest doubles of the real parts of A;
 *   Arehilo  are the second lowest doubles of the real parts of A;
 *   Arelolo  are the lowest doubles of the real parts of A;
 *   Aimhihi  are the highest doubles of the imaginary parts of A; 
 *   Aimlohi  are the second highest doubles of the imaginary parts of A; 
 *   Aimhilo  are the second lowest doubles of the imaginary parts of A;
 *   Aimlolo  are the lowest doubles of the imaginary parts of A;
 *   Brehihi  are the highest doubles of the real parts of B;
 *   Brelohi  are the second highest doubles of the real parts of B;
 *   Brehilo  are the second lowest doubles of the real parts of B;
 *   Brelolo  are the lowest doubles of the real parts of B;
 *   Bimhihi  are the highest doubles of the imaginary parts of B;
 *   Bimlohi  are the second highest doubles of the imaginary parts of B;
 *   Bimhilo  are the second lowest doubles of the imaginary parts of B;
 *   Bimlolo  are the lowest doubles of the imaginary parts of B;
 *   Crehihi  has space allocated for a rows-by-cols matrix;
 *   Crelohi  has space allocated for a rows-by-cols matrix;
 *   Crehilo  has space allocated for a rows-by-cols matrix;
 *   Crelolo  has space allocated for a rows-by-cols matrix;
 *   Cimhihi  has space allocated for a rows-by-cols matrix;
 *   Cimlohi  has space allocated for a rows-by-cols matrix;
 *   Cimhilo  has space allocated for a rows-by-cols matrix;
 *   Cimlolo  has space allocated for a rows-by-cols matrix.
 *
 * ON RETURN :
 *   Crehihi  are the highest doubles of the real parts of C,
              the product of A with B;
 *   Crelohi  second highest doubles of the real parts of C;
 *   Crehilo  second lowest doubles of the real parts of C;
 *   Crelolo  lowest doubles of the real parts of C;
 *   Cimhihi  highest doubles of the imaginary parts of C;
 *   Cimlohi  second highest doubles of the imaginary parts of C;
 *   Cimhilo  second lowest doubles of the imaginary parts of C;
 *   Cimlolo  lowest doubles of the imaginary parts of C. */

void CPU_dbl4_factors_forward
 ( int dim,
   double **Lhihi, double **Llohi, double **Lhilo, double **Llolo,
   double *bhihi, double *blohi, double *bhilo, double *blolo,
   double *xhihi, double *xlohi, double *xhilo, double *xlolo );
/*
 * DESCRIPTION :
 *   Solves the lower triangular system L*x = b with forward substitution,
 *   on real data.
 *
 * REQUIRED : the matrix L has ones on the diagonal.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in the matrix L;
 *   Lhihi    highest doubles of a lower triangular matrix;
 *   Llohi    second highest doubles of a lower triangular matrix;
 *   Lhilo    second lowest doubles of a lower triangular matrix;
 *   Llolo    lowest doubles of a lower triangular matrix;
 *   bhihi    highest doubles of a right hand side vector;
 *   blohi    second highest doubles of a right hand side vector;
 *   bhilo    second lowest doubles of a right hand side vector;
 *   blolo    lowest doubles of a right hand side vector;
 *   xhihi    space for dim doubles;
 *   xlohi    space for dim doubles;
 *   xhilo    space for dim doubles;
 *   xlolo    space for dim doubles.
 *
 * ON RETURN :
 *   xhihi    highest doubles of the solution to L*x = b;
 *   xlohi    second highest doubles of the solution to L*x = b;
 *   xhilo    second lowest doubles of the solution to L*x = b;
 *   xlolo    lowest doubles of the solution to L*x = b. */

void CPU_cmplx4_factors_forward
 ( int dim,
   double **Lrehihi, double **Lrelohi, double **Lrehilo, double **Lrelolo,
   double **Limhihi, double **Limlohi, double **Limhilo, double **Limlolo,
   double *brehihi, double *brelohi, double *brehilo, double *brelolo,
   double *bimhihi, double *bimlohi, double *bimhilo, double *bimlolo,
   double *xrehihi, double *xrelohi, double *xrehilo, double *xrelolo,
   double *ximhihi, double *ximlohi, double *ximhilo, double *ximlolo );
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

void CPU_dbl4_factors_backward
 ( int dim,
   double **Uhihi, double **Ulohi, double **Uhilo, double **Ulolo,
   double *bhihi, double *blohi, double *bhilo, double *blolo,
   double *xhihi, double *xlohi, double *xhilo, double *xlolo );
/*
 * DESCRIPTION :
 *   Solves the upper triangular system U*x = b with back substitution,
 *   on real data.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in the matrix U;
 *   Uhihi    highest doubles of an upper triangular matrix;
 *   Ulohi    second highest doubles of an upper triangular matrix;
 *   Uhilo    second lowest doubles of an upper triangular matrix;
 *   Ulolo    lowest doubles of an upper triangular matrix;
 *   bhihi    highest doubles of a right hand side vector;
 *   blohi    second highest doubles of a right hand side vector;
 *   bhilo    second lowest doubles of a right hand side vector;
 *   blolo    lowest doubles of a right hand side vector;
 *   xhihi    space for dim doubles;
 *   xlohi    space for dim doubles;
 *   xhilo    space for dim doubles.
 *   xlolo    space for dim doubles.
 *
 * ON RETURN :
 *   xhihi    highest doubles of the solution to U*x = b;
 *   xlohi    second highest doubles of the solution to U*x = b;
 *   xhilo    second lowest doubles of the solution to U*x = b;
 *   xlolo    lowest doubles of the solution to U*x = b. */

void CPU_cmplx4_factors_backward
 ( int dim,
   double **Urehihi, double **Urelohi, double **Urehilo, double **Urelolo,
   double **Uimhihi, double **Uimlohi, double **Uimhilo, double **Uimlolo,
   double *brehihi, double *brelohi, double *brehilo, double *brelolo,
   double *bimhihi, double *bimlohi, double *bimhilo, double *bimlolo,
   double *xrehihi, double *xrelohi, double *xrehilo, double *xrelolo,
   double *ximhihi, double *ximlohi, double *ximhilo, double *ximlolo );
/*
 * DESCRIPTION :
 *   Solves the upper triangular system U*x = b with back substitution,
 *   on real data.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in the upper triangular matrix U;
 *   Urehihi  are the highest doubles of the real parts of U;
 *   Urelohi  are the second highest doubles of the real parts of U;
 *   Urehilo  are the second lowest doubles of the real parts of U;
 *   Urelolo  are the lowest doubles of the real parts of U;
 *   Uimhihi  are the highest doubles of the imaginary parts of U;
 *   Uimlohi  are the second highest doubles of the imaginary parts of U;
 *   Uimhilo  are the second lowest doubles of the imaginary parts of U;
 *   Uimlolo  are the lowest doubles of the imaginary parts of U;
 *   brehihi  are the highest doubles of the real parts of b;
 *   brelohi  are the second highest doubles of the real parts of b;
 *   brehilo  are the second lowest doubles of the real parts of b;
 *   brelolo  are the lowest doubles of the real parts of b;
 *   bimhihi  are the highest doubles of the imaginary parts of b;
 *   bimlohi  are the second highest doubles of the imaginary parts of b;
 *   bimhilo  are the second lowest doubles of the imaginary parts of b;
 *   bimlolo  are the lowest doubles of the imaginary parts of b;
 *   xrehihi  has space for dim doubles;
 *   xrelohi  has space for dim doubles;
 *   xrehilo  has space for dim doubles;
 *   xrelolo  has space for dim doubles;
 *   ximhihi  has space for dim doubles;
 *   ximlohi  has space for dim doubles;
 *   ximhilo  has space for dim doubles;
 *   ximlolo  has space for dim doubles.
 *
 * ON RETURN :
 *   xrehihi  are the highest doubles of the real parts of x;
 *   xrelohi  are the second highest doubles of the real parts of x;
 *   xrehilo  are the second lowest doubles of the imaginary parts of x;
 *   xrelolo  are the lowest doubles of the imaginary parts of x;
 *   ximhihi  are the highest doubles of the real parts of x;
 *   ximlohi  are the second highest doubles of the real parts of x;
 *   ximhilo  are the second lowest doubles of the imaginary parts of x;
 *   ximlolo  are the lowest doubles of the imaginary parts of x. */

void CPU_dbl4_factors_lufac
 ( int dim, double **Ahihi, double **Alohi, double **Ahilo, double **Alolo,
   int *pivots );
/*
 * DESCRIPTION :
 *   Does an inplace LU factorization with pivoting on the matrix A,
 *   on real data.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in the matrix A;
 *   Ahihi    highest doubles of a matrix of dimension dim;
 *   Alohi    second highest doubles of a matrix of dimension dim;
 *   Ahilo    second lowest doubles of a matrix of dimension dim;
 *   Alolo    lowest doubles of a matrix of dimension dim;
 *   pivots   space for dim pivots.
 *
 * ON RETURN :
 *   Ahihi    the lower triangular part of Ahihi contains the highest
 *            doubles of the multipliers and the upper triangular part
 *            of Ahihi has the highest doubles of the row reduced A;
 *   Alohi    the lower triangular part of Alohi contains the second highest
 *            doubles of the multipliers and the upper triangular part
 *            of Alohi has the second highest doubles of the row reduced A;
 *   Ahilo    the lower triangular part of Ahilo contains the second lowest
 *            doubles of the multipliers and the upper triangular part
 *            of Ahilo has the second lowest doubles of the row reduced A;
 *   Alolo    the lower triangular part of Alolo contains the lowest
 *            doubles of the multipliers and the upper triangular part
 *            of Alolo has the lowest doubles of the row reduced A;
 *   pivots   are the pivots used. */

void CPU_cmplx4_factors_lufac
 ( int dim,
   double **Arehihi, double **Arelohi, double **Arehilo, double **Arelolo,
   double **Aimhihi, double **Aimlohi, double **Aimhilo, double **Aimlolo,
   int *pivots );
/*
 * DESCRIPTION :
 *   Does an inplace LU factorization with pivoting on the matrix A,
 *   on real data.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in the matrix A;
 *   Arehihi  highest doubles of the real parts of A;
 *   Arelohi  second highest doubles of the real parts of A;
 *   Arehilo  second lowest doubles of the real parts of A;
 *   Arelolo  lowest doubles of the real parts of A;
 *   Aimhihi  highest doubles of the imaginary parts of A;
 *   Aimlohi  second highest doubles of the imaginary parts of A;
 *   Aimhilo  second lowest doubles of the imaginary parts of A;
 *   Aimlolo  lowest doubles of the imaginary parts of A;
 *   pivots   space for dim pivots.
 *
 * ON RETURN :
 *   Arehihi  highest doubles of the real parts
 *            of the multipliers and the row reduced A;
 *   Arelohi  second highest doubles of the real parts
 *            of the multipliers and the row reduced A;
 *   Arehilo  second lowest doubles of the real parts
 *            of the multipliers and the row reduced A;
 *   Arelolo  lowest doubles of the real parts
 *            of the multipliers and the row reduced A;
 *   Aimhihi  highest doubles of the imaginary parts
 *            of the multipliers and the row reduced A;
 *   Aimlohi  second highest doubles of the imaginary parts
 *            of the multipliers and the row reduced A;
 *   Aimhilo  second lowest doubles of the imaginary parts
 *            of the multipliers and the row reduced A;
 *   Aimlolo  lowest doubles of the imaginary parts
 *            of the multipliers and the row reduced A;
 *   pivots   are the pivots used. */

void CPU_dbl4_factors_lusolve
 ( int dim,
   double **Ahihi, double **Alohi, double **Ahilo, double **Alolo,
   int *pivots,
   double *bhihi, double *blohi, double *bhilo, double *blolo,
   double *xhihi, double *xlohi, double *xhilo, double *xlolo );
/*
 * DESCRIPTION :
 *   Does an inplace LU factorization with pivoting on the matrix A,
 *   to solve the system A*x = b, on real data.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in the matrix A;
 *   Ahihi    highest doubles of a matrix of dimension dim;
 *   Alohi    second highest doubles of a matrix of dimension dim;
 *   Ahilo    second lowest doubles of a matrix of dimension dim;
 *   Alolo    lowest doubles of a matrix of dimension dim;
 *   pivots   space for dim pivots;
 *   bhihi    highest doubles of the right hand side vector;
 *   blohi    second highest doubles of the right hand side vector;
 *   bhilo    second lowest doubles of the right hand side vector;
 *   blolo    lowest doubles of the right hand side vector;
 *   xhihi    space for dim doubles;
 *   xlohi    space for dim doubles;
 *   xhilo    space for dim doubles;
 *   xlolo    space for dim doubles.
 *
 * ON RETURN :
 *   Ahihi    the lower triangular part of Ahihi contains the highest
 *            doubles of the multipliers and the upper triangular part
 *            of Ahihi has the highest doubles of the row reduced A;
 *   Alohi    the lower triangular part of Alohi contains the second highest
 *            doubles of the multipliers and the upper triangular part
 *            of Alohi has the second highest doubles of the row reduced A;
 *   Ahilo    the lower triangular part of Ahilo contains the second lowest
 *            doubles of the multipliers and the upper triangular part
 *            of Ahilo has the second lowest doubles of the row reduced A;
 *   Alolo    the lower triangular part of Alolo contains the lowest
 *            doubles of the multipliers and the upper triangular part
 *            of Alolo has the lowest doubles of the row reduced A;
 *   bhihi    used as work space;
 *   blohi    used as work space;
 *   bhilo    used as work space;
 *   blolo    used as work space;
 *   xhihi    highest doubles of the solution to A*x = b;
 *   xlohi    second highest doubles of the solution to A*x = b;
 *   xhilo    second lowest doubles of the solution to A*x = b;
 *   xlolo    lowest doubles of the solution to A*x = b. */

void CPU_cmplx4_factors_lusolve
 ( int dim,
   double **Arehihi, double **Arelohi, double **Arehilo, double **Arelolo,
   double **Aimhihi, double **Aimlohi, double **Aimhilo, double **Aimlolo,
   int *pivots,
   double *brehihi, double *brelohi, double *brehilo, double *brelolo,
   double *bimhihi, double *bimlohi, double *bimhilo, double *bimlolo,
   double *xrehihi, double *xrelohi, double *xrehilo, double *xrelolo,
   double *ximhihi, double *ximlohi, double *ximhilo, double *ximlolo );
/*
 * DESCRIPTION :
 *   Does an inplace LU factorization with pivoting on the matrix A,
 *   to solve the system A*x = b, on complex data.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in the matrix A;
 *   Arehihi  highest doubles of the real parts of A;
 *   Arelohi  second highest doubles of the real parts of A;
 *   Arehilo  second lowest doubles of the imaginary parts of A;
 *   Arelolo  lowest doubles of the imaginary parts of A;
 *   Aimhihi  highest doubles of the real parts of A;
 *   Aimlohi  second highest doubles of the real parts of A;
 *   Aimhilo  second lowest doubles of the imaginary parts of A;
 *   Aimlolo  lowest doubles of the imaginary parts of A;
 *   pivots   space for dim pivots;
 *   brehihi  highest doubles of the real parts of b;
 *   brelohi  second highest doubles of the real parts of b;
 *   brehilo  second lowest doubles of the real parts of b;
 *   brelolo  lowest doubles of the real parts of b;
 *   bimhihi  highest doubles of the imaginary parts of b;
 *   bimlohi  second highest doubles of the imaginary parts of b;
 *   bimhilo  second lowest doubles of the imaginary parts of b;
 *   bimlolo  lowest doubles of the imaginary parts of b;
 *   xrehihi  space for dim doubles;
 *   xrelohi  space for dim doubles;
 *   xrehilo  space for dim doubles;
 *   xrelolo  space for dim doubles;
 *   ximhihi  space for dim doubles;
 *   ximlohi  space for dim doubles;
 *   ximhilo  space for dim doubles;
 *   ximlolo  space for dim doubles.
 *
 * ON RETURN :
 *   Arehihi  highest doubles of the real parts
 *            of the multipliers and the row reduced A;
 *   Arelohi  second highest doubles of the real parts
 *            of the multipliers and the row reduced A;
 *   Arehilo  second lowest doubles of the real parts
 *            of the multipliers and the row reduced A;
 *   Arelolo  lowest doubles of the real parts
 *            of the multipliers and the row reduced A;
 *   Aimhihi  highest doubles of the imaginary parts
 *            of the multipliers and the row reduced A;
 *   Aimlohi  second highest doubles of the imaginary parts
 *            of the multipliers and the row reduced A;
 *   Aimhilo  second lowest doubles of the imaginary parts
 *            of the multipliers and the row reduced A;
 *   Aimlolo  lowest doubles of the imaginary parts
 *            of the multipliers and the row reduced A;
 *   brehihi  used as work space;
 *   brelohi  used as work space;
 *   brehilo  used as work space;
 *   brelolo  used as work space;
 *   bimhihi  used as work space;
 *   bimlohi  used as work space;
 *   bimhilo  used as work space;
 *   bimlolo  used as work space;
 *   xrehihi  highest doubles of the real parts of the solution x;
 *   xrelohi  second highest doubles of the real parts of x;
 *   xrehilo  second lowest doubles of the real parts of x;
 *   xrelolo  lowest doubles of the real parts of x;
 *   ximhihi  highest doubles of the imaginary parts of x;
 *   ximlohi  second highest doubles of the imaginary parts of x;
 *   ximhilo  second lowest doubles of the imaginary parts of x;
 *   ximlolo  lowest doubles of the imaginary parts of x. */

void CPU_dbl4_factors_house
 ( int n,
   double *xhihi, double *xlohi, double *xhilo, double *xlolo,
   double *vhihi, double *vlohi, double *vhilo, double *vlolo,
   double *betahihi, double *betalohi, double *betahilo, double *betalolo );
/*
 * DESCRIPTION :
 *   Computes the Householder vector of an n-dimensional vector x.
 *
 * ON ENTRY :
 *   n        dimension of the vector x;
 *   xhihi    the n high doubles of x;
 *   xlohi    the n high doubles of x;
 *   xhilo    the n low doubles of x;
 *   xlolo    the n low doubles of x;
 *   vhihi    space for n doubles;
 *   vlohi    space for n doubles;
 *   vhilo    space for n doubles;
 *   vlolo    space for n doubles.
 *
 * ON RETURN :
 *   vhihi    highest doubles of the Householder vector;
 *   vlohi    second highest doubles of the Householder vector;
 *   vhilo    second lowest doubles of the Householder vector;
 *   vlolo    lowest doubles of the Householder vector;
 *   betahihi equals the highest double of 2/(transpose(v)*v);
 *   betalohi equals the second highest double of 2/(transpose(v)*v);
 *   betahilo equals the second lowest double of 2/(transpose(v)*v);
 *   betalolo equals the lowest double of 2/(transpose(v)*v). */

void CPU_cmplx4_factors_house 
( int n,
  double *xrehihi, double *xrelohi, double *xrehilo, double *xrelolo,
  double *ximhihi, double *ximlohi, double *ximhilo, double *ximlolo,
  double *vrehihi, double *vrelohi, double *vrehilo, double *vrelolo,
  double *vimhihi, double *vimlohi, double *vimhilo, double *vimlolo,
  double *betahihi, double *betalohi, double *betahilo, double *betalolo );
/*
 * DESCRIPTION :
 *   Computes the Householder vector of an n-dimensional vector x.
 *
 * ON ENTRY :
 *   n        dimension of the vector x;
 *   xrehihi  highest doubles of the real parts of the vector x;
 *   xrelohi  second highest doubles of the real parts of the vector x;
 *   xrehilo  second lowest doubles of the real parts of the vector x;
 *   xrelolo  lowest doubles of the real parts of the vector x;
 *   ximhihi  highest doubles of the imaginary parts of the vector x;
 *   ximlohi  second highest doubles of the imaginary parts of the vector x;
 *   ximhilo  second lowest doubles of the imaginary parts of the vector x;
 *   ximlolo  lowest doubles of the imaginary parts of the vector x;
 *   vrehihi  space for n doubles;
 *   vrelohi  space for n doubles;
 *   vrehilo  space for n doubles;
 *   vrelolo  space for n doubles;
 *   vimhihi  space for n doubles;
 *   vimlohi  space for n doubles;
 *   vimhilo  space for n doubles;
 *   vimlolo  space for n doubles.
 *
 * ON RETURN :
 *   vrehihi  highest doubles of the real parts of the Householder vector v;
 *   vrelohi  second highest doubles of the real parts of v;
 *   vrehilo  second lowest doubles of the real parts of v;
 *   vrelolo  lowest doubles of the real parts of v;
 *   vimhihi  highest doubles of the imaginary parts of v;
 *   vimlohi  second highest doubles of the imaginary parts of v;
 *   vimhilo  second lowest doubles of the imaginary parts of v;
 *   vimlolo  lowest doubles of the imaginary parts of v;
 *   betahihi is the highest double of 2/(transpose(v)*v);
 *   betalohi is the second highest double of 2/(transpose(v)*v);
 *   betahilo is the second lowest double of 2/(transpose(v)*v);
 *   betalolo is the lowest double of 2/(transpose(v)*v). */

void CPU_dbl4_factors_leftRupdate
 ( int nrows, int ncols, int k,
   double **Rhihi, double **Rlohi, double **Rhilo, double **Rlolo,
   double *vhihi, double *vlohi, double *vhilo, double *vlolo,
   double betahihi, double betalohi, double betahilo, double betalolo );
/*
 * DESCRIPTION :
 *   Applies the Householder matrix to R.
 *
 * ON ENTRY :
 *   nrows    number of rows in the matrix R;
 *   ncols    number of columns in the matrix R;
 *   k        current column index in R;
 *   Rhihi    highest doubles of an nrows-by-ncols matrix;
 *   Rlohi    second highest doubles of an nrows-by-ncols matrix;
 *   Rhilo    second lowest doubles of an nrows-by-ncols matrix;
 *   Rlolo    lowest doubles of an nrows-by-ncols matrix;
 *   vhihi    highest doubles of the Householder vector;
 *   vlohi    second highest doubles of the Householder vector;
 *   vhilo    second lowest doubles of the Householder vector;
 *   vlolo    lowest doubles of the Householder vector;
 *   betahihi is the betahihi computed by CPU_dbl4_factors_house;
 *   betalohi is the betalohi computed by CPU_dbl4_factors_house;
 *   betahilo is the betahilo computed by CPU_dbl4_factors_house;
 *   betalolo is the betalolo computed by CPU_dbl4_factors_house.
 *
 * ON RETURN :
 *   Rhihi    highest doubles of the update with the Householder matrix;
 *   Rlohi    second highest doubles of the update;
 *   Rhilo    second lowest doubles of the update;
 *   Rlolo    lowest doubles of the update with the Householder matrix. */

void CPU_cmplx4_factors_leftRupdate
 ( int nrows, int ncols, int k,
   double **Rrehihi, double **Rrelohi, double **Rrehilo, double **Rrelolo,
   double **Rimhihi, double **Rimlohi, double **Rimhilo, double **Rimlolo,
   double *vrehihi, double *vrelohi, double *vrehilo, double *vrelolo,
   double *vimhihi, double *vimlohi, double *vimhilo, double *vimlolo,
   double betahihi, double betalohi, double betahilo, double betalolo );
/*
 * DESCRIPTION :
 *   Applies the Householder matrix to R.
 *
 * ON ENTRY :
 *   nrows    number of rows in the matrix R;
 *   ncols    number of columns in the matrix R;
 *   k        current column index in R;
 *   Rrehihi  highest doubles of the real parts of an nrows-by-ncols matrix R;
 *   Rrelohi  second highest doubles of the real parts of R;
 *   Rrehilo  second lowest doubles of the real parts of R;
 *   Rrelolo  lowest doubles of the real parts of R;
 *   Rimhihi  highest doubles of the imaginary parts of R;
 *   Rimlohi  second highest doubles of the imaginary parts of R;
 *   Rimhilo  second lowest doubles of the imaginary parts of R;
 *   Rimlolo  lowest doubles of the imaginary parts of R;
 *   vrehihi  highest doubles of the real parts of the Householder vector v;
 *   vrelohi  second highest doubles of the real parts of v;
 *   vrehilo  second lowest doubles of the real parts of v;
 *   vrelolo  lowest doubles of the real parts of v;
 *   vimhihi  highest doubles of the imaginary parts of v;
 *   vimlohi  second highest doubles of the imaginary parts of v;
 *   vimhilo  second lowest doubles of the imaginary parts of v;
 *   vimlolo  lowest doubles of the imaginary parts of v;
 *   betahihi is the betahihi computed by CPU_cmplx4_factors_house;
 *   betalohi is the betalohi computed by CPU_cmplx4_factors_house;
 *   betahilo is the betahilo computed by CPU_cmplx4_factors_house;
 *   betalolo is the betalolo computed by CPU_cmplx4_factors_house.
 *
 * ON RETURN :
 *   Rrehihi  highest doubles of the real parts of the update
 *            with the Householder matrix;
 *   Rrelohi  second highest doubles of the real parts of the update;
 *   Rrehilo  second lowest doubles of the real parts of the update;
 *   Rrelolo  lowest doubles of the real parts of the update;
 *   Rimhihi  highest doubles of the imaginary parts of the update;
 *   Rimlohi  second highest doubles of the imaginary parts of the update;
 *   Rimhilo  second lowest doubles of the imaginary parts of the update;
 *   Rimlolo  lowest doubles of the imaginary parts of the update. */

void CPU_dbl4_factors_rightQupdate
 ( int n, int k,
   double **Qhihi, double **Qlohi, double **Qhilo, double **Qlolo,
   double *vhihi, double *vlohi, double *vhilo, double *vlolo,
   double betahihi, double betalohi, double betahilo, double betalolo );
/*
 * DESCRIPTION :
 *   Applies the Householder matrix to Q.
 *
 * ON ENTRY :
 *   n        dimension of the matrix Q;
 *   k        current column index in Q;
 *   Qhihi    highest doubles of an n-by-n matrix;
 *   Qlohi    highest doubles of an n-by-n matrix;
 *   Qhilo    lowest doubles of an n-by-n matrix;
 *   Qlolo    lowest doubles of an n-by-n matrix;
 *   vhihi    highest doubles of the Householder vector;
 *   vlohi    second highest doubles of the Householder vector;
 *   vhilo    second lowest doubles of the Householder vector;
 *   vlolo    lowest doubles of the Householder vector;
 *   betahihi is the betahihi computed by CPU_dbl4_factors_house;
 *   betalohi is the betalohi computed by CPU_dbl4_factors_house;
 *   betahilo is the betahilo computed by CPU_dbl4_factors_house;
 *   betalolo is the betalolo computed by CPU_dbl4_factors_house.
 *
 * ON RETURN :
 *   Qhihi    highest doubles of the update with the Householder matrix;
 *   Qlohi    second highest doubles of the update;
 *   Qhilo    second lowest doubles of the update;
 *   Qlolo    lowest doubles of the update with the Householder matrix. */

void CPU_cmplx4_factors_rightQupdate
 ( int n, int k,
   double **Qrehihi, double **Qrelohi, double **Qrehilo, double **Qrelolo,
   double **Qimhihi, double **Qimlohi, double **Qimhilo, double **Qimlolo,
   double *vrehihi, double *vrelohi, double *vrehilo, double *vrelolo,
   double *vimhihi, double *vimlohi, double *vimhilo, double *vimlolo,
   double betahihi, double betalohi, double betahilo, double betalolo );
/*
 * DESCRIPTION :
 *   Applies the Householder matrix to Q.
 *
 * ON ENTRY :
 *   n        dimension of the matrix Q;
 *   k        current column index in Q;
 *   Qrehihi  highest doubles of the real parts of an n-by-n matrix Q;
 *   Qrelohi  second highest doubles of the real parts of Q;
 *   Qrehilo  second lowest doubles of the real parts of Q;
 *   Qrelolo  lowest doubles of the real parts of Q;
 *   Qimhihi  highest doubles of the imaginary parts of Q;
 *   Qimlohi  second highest doubles of the imaginary parts of Q;
 *   Qimhilo  second lowest doubles of the imaginary parts of Q;
 *   Qimlolo  lowest doubles of the imaginary parts of Q;
 *   vrehihi  highest doubles of the real parts of the Householder vector v;
 *   vrelohi  second highest doubles of the real parts of v;
 *   vrehilo  second lowest doubles of the real parts of v;
 *   vrelolo  lowest doubles of the real parts of v;
 *   vimhihi  highest doubles of the imaginary parts of v;
 *   vimlohi  second highest doubles of the imaginary parts of v;
 *   vimhilo  second lowest doubles of the imaginary parts of v;
 *   vimlolo  lowest doubles of the imaginary parts of v;
 *   betahihi is the betahihi computed by CPU_cmplx4_factors_house;
 *   betalohi is the betalohi computed by CPU_cmplx4_factors_house;
 *   betahilo is the betahilo computed by CPU_cmplx4_factors_house;
 *   betalolo is the betalolo computed by CPU_cmplx4_factors_house.
 *
 * ON RETURN :
 *   Qrehihi  highest doubles of the real parts of the update
 *            with the Householder matrix;
 *   Qrelohi  second highest doubles of the real parts of the update;
 *   Qrehilo  second lowest doubles of the real parts of the update;
 *   Qrelolo  lowest doubles of the real parts of the update;
 *   Qimhihi  highest doubles of the imaginary parts of the update;
 *   Qimlohi  second highest doubles of the imaginary parts of the update;
 *   Qimhilo  second lowest doubles of the imaginary parts of the update;
 *   Qimlolo  lowest doubles of the imaginary parts of the update. */

void CPU_dbl4_factors_houseqr
 ( int nrows, int ncols,
   double **Ahihi, double **Alohi, double **Ahilo, double **Alolo,
   double **Qhihi, double **Qlohi, double **Qhilo, double **Qlolo,
   double **Rhihi, double **Rlohi, double **Rhilo, double **Rlolo );
/*
 * DESCRIPTION :
 *   Applies Householder matrices to compute a QR decomposition of A.
 *
 * REQUIRED : nrows >= ncols.
 *
 * ON ENTRY :
 *   nrows    number of rows of A;
 *   ncols    number of columns of A;
 *   Ahihi    highest doubles of an nrows-by-ncols matrix,
 *            stored as nrows arrays of ncols numbers;
 *   Alohi    second highest doubles of an nrows-by-ncols matrix,
 *            stored as nrows arrays of ncols numbers;
 *   Ahilo    second lowest doubles of an nrows-by-ncols matrix,
 *            stored as nrows arrays of ncols numbers;
 *   Alolo    lowest doubles of an nrows-by-ncols matrix,
 *            stored as nrows arrays of ncols numbers;
 *   Qhihi    space for an nrows-by-nrows matrix;
 *   Qlohi    space for an nrows-by-nrows matrix;
 *   Qhilo    space for an nrows-by-nrows matrix;
 *   Qlolo    space for an nrows-by-nrows matrix;
 *   Rhihi    space for an nrows-by-ncols matrix;
 *   Rlohi    space for an nrows-by-ncols matrix;
 *   Rhilo    space for an nrows-by-ncols matrix;
 *   Rlolo    space for an nrows-by-ncols matrix.
 *
 * ON RETURN :
 *   Qhihi    highest doubles of an orthogonal matrix, transpose(Q)*A = R;
 *   Qlohi    second highest doubles of an orthogonal matrix;
 *   Qhilo    second lowest doubles of an orthogonal matrix;
 *   Qlolo    lowest doubles of an orthogonal matrix, transpose(Q)*A = R;
 *   Rhihi    highest doubles of the reduced upper triangular form;
 *   Rlohi    second highest doubles of the reduced upper triangular form;
 *   Rhilo    second lowest doubles of the reduced upper triangular form;
 *   Rlolo    lowest doubles of the reduced upper triangular form. */

void CPU_cmplx4_factors_houseqr
 ( int nrows, int ncols,
   double **Arehihi, double **Arelohi, double **Arehilo, double **Arelolo,
   double **Aimhihi, double **Aimlohi, double **Aimhilo, double **Aimlolo,
   double **Qrehihi, double **Qrelohi, double **Qrehilo, double **Qrelolo,
   double **Qimhihi, double **Qimlohi, double **Qimhilo, double **Qimlolo,
   double **Rrehihi, double **Rrelohi, double **Rrehilo, double **Rrelolo,
   double **Rimhihi, double **Rimlohi, double **Rimhilo, double **Rimlolo );
/*
 * DESCRIPTION :
 *   Applies Householder matrices to compute a QR decomposition of A.
 *
 * REQUIRED : nrows >= ncols.
 *
 * ON ENTRY :
 *   nrows    number of rows of A;
 *   ncols    number of columns of A;
 *   Arehihi  highest doubles of the real parts of an nrows-by-ncols
 *            matrix A, stored as nrows arrays of ncols numbers;
 *   Arelohi  second highest doubles of the real parts of A;
 *   Arehilo  second lowest doubles of the real parts of :
 *   Arelolo  lowest doubles of the real parts of A;
 *   Aimhi    highest doubles of the imaginary parts of A;
 *   Aimhi    second highest doubles of the imaginary parts of A;
 *   Aimlo    second lowest doubles of the imaginary parts of A;
 *   Aimlo    lowest doubles of the imaginary parts of A;
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
 *   Rimlolo  space for an nrows-by-ncols matrix.
 *
 * ON RETURN :
 *   Qrehihi  highest doubles of the real parts
 *            of the orthogonal matrix Q, transpose(Q)*A = R;
 *   Qrelohi  second highest doubles of the real parts of Q;
 *   Qrehilo  second lowest doubles of the real parts of Q;
 *   Qrelolo  lowest doubles of the real parts of Q;
 *   Qimhihi  highest doubles of the imaginary parts of Q;
 *   Qimlohi  second highest doubles of the imaginary parts of Q;
 *   Qimhilo  second lowest doubles of the imaginary parts of Q;
 *   Qimlolo  lowest doubles of the imaginary parts of Q;
 *   Rrehi    highest doubles of the real parts of R,
 *            the reduced upper triangular form of A;
 *   Rrelohi  second highest doubles of the real parts of R,
 *   Rrehilo  second lowest doubles of the real parts of R;
 *   Rrelolo  lowest doubles of the real parts of R;
 *   Rimhihi  highest doubles of the imaginary parts of R;
 *   Rimlohi  second highest doubles of the imaginary parts of R;
 *   Rimhilo  second lowest doubles of the imaginary parts of R;
 *   Rimlolo  lowest doubles of the imaginary parts of R. */

void CPU_dbl4_factors_qrbs
 ( int nrows, int ncols,
   double **Qhihi, double **Qlohi, double **Qhilo, double **Qlolo,
   double **Rhihi, double **Rlohi, double **Rhilo, double **Rlolo,
   double *rhshihi, double *rhslohi, double *rhshilo, double *rhslolo,
   double *solhihi, double *sollohi, double *solhilo, double *sollolo,
   double *wrkvechihi, double *wrkveclohi,
   double *wrkvechilo, double *wrkveclolo ); 
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
 *   Qhihi    highest doubles of the Q of the QR factorization;
 *   Qlohi    second highest doubles of the Q of the QR factorization;
 *   Qhilo    second lowest doubles of the Q of the QR factorization;
 *   Qlolo    lowest doubles of the Q of the QR factorization;
 *   Rhihi    highest doubles of the R of the QR factorization;
 *   Rlohi    second highest doubles of the R of the QR factorization;
 *   Rhilo    second lowest doubles of the R of the QR factorization;
 *   Rlolo    lowest doubles of the R of the QR factorization;
 *   rhshihi  highest doubles of the right hand side vector, of size nrows;
 *   rhslohi  second highest doubles of the right hand side vector;
 *   rhshilo  second lowest doubles of the right hand side vector;
 *   rhslolo  lowest doubles of the right hand side vector;
 *   solhihi  space for the highest doubles of the solution, of size ncols;
 *   sollohi  space for the second highest doubles of the solution;
 *   solhilo  space for the second lowest doubles of the solution;
 *   sollolo  space for the lowest doubles of the solution;
 *   wrkvechihi has work space for nrows elements;
 *   wrkveclohi has work space for nrows elements;
 *   wrkvechilo has work space for nrows elements;
 *   wrkveclolo has work space for nrows elements.
 *
 * ON RETURN :
 *   solhihi  highest doubles of the least squares solution;
 *   sollohi  second highest doubles of the least squares solution;
 *   solhilo  second lowest doubles of the least squares solution;
 *   sollolo  lowest doubles of the least squares solution;
 *   wrkvechihi has the highest doubles of the product of Q^T with rhs;
 *   wrkveclohi has the second highest doubles of Q^T*rhs;
 *   wrkvechilo has the second lowest doubles of Q^T*rhs;
 *   wrkveclolo has the lowest doubles of Q^T*rhs. */

void CPU_cmplx4_factors_qrbs
 ( int nrows, int ncols,
   double **Qrehihi, double **Qrelohi, double **Qrehilo, double **Qrelolo,
   double **Qimhihi, double **Qimlohi, double **Qimhilo, double **Qimlolo, 
   double **Rrehihi, double **Rrelohi, double **Rrehilo, double **Rrelolo,
   double **Rimhihi, double **Rimlohi, double **Rimhilo, double **Rimlolo,
   double *rhsrehihi, double *rhsrelohi, double *rhsrehilo, double *rhsrelolo,
   double *rhsimhihi, double *rhsimlohi, double *rhsimhilo, double *rhsimlolo,
   double *solrehihi, double *solrelohi, double *solrehilo, double *solrelolo,
   double *solimhihi, double *solimlohi, double *solimhilo, double *solimlolo,
   double *wrkvecrehihi, double *wrkvecrelohi,
   double *wrkvecrehilo, double *wrkvecrelolo,
   double *wrkvecimhihi, double *wrkvecimlohi,
   double *wrkvecimhilo, double *wrkvecimlolo );
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
 *   Qrehihi  highest doubles of the real parts of Q of the QR factorization;
 *   Qrelohi  second highest doubles of the real parts of Q of the QR;
 *   Qrehilo  second lowest doubles of the real parts of Q of the QR;
 *   Qrelolo  lowest doubles of the real parts of Q of the QR factorization;
 *   Qimhihi  highest doubles of the imaginary parts of Q of the QR;
 *   Qimlohi  second highest doubles of the imaginary parts of Q of the QR;
 *   Qimhilo  second lowest doubles of the imaginary parts of Q of the QR;
 *   Qimlolo  lowest doubles of the imaginary parts of Q of the QR;
 *   Rrehihi  highest doubles of the real parts of R of the QR;
 *   Rrelohi  second highest doubles of the real parts of R of the QR;
 *   Rrehilo  second lowest doubles of the real parts of R of the QR;
 *   Rrelolo  lowest doubles of the real parts of R of the QR;
 *   Rimhihi  highest doubles of the imaginary parts of R of the QR;
 *   Rimlohi  second highest doubles of the imaginary parts of R of the QR;
 *   Rimhilo  second lowest doubles of the imaginary parts of R of the QR;
 *   Rimlolo  lowest doubles of the imaginary parts of R of the QR;
 *   rhsrehihi are the highest doubles of the real parts of right hand side,
 *            rhs has size nrows;
 *   rhsrelohi are the second highest doubles of the real parts of rhs;
 *   rhsrehilo are the second lowest doubles of the real parts of rhs;
 *   rhsrelolo are the lowest doubles of the real parts of rhs;
 *   rhsimhihi are the highest doubles of the imaginary parts of rhs;
 *   rhsimlohi are the second highest doubles of the imaginary parts of rhs;
 *   rhsimhilo are the lowest doubles of the imaginary parts of rhs;
 *   rhsimlolo are lowest doubles of the imaginary parts of rhs;
 *   solrehihi has space for the highest doubles of the real parts
 *            of the solution, of size ncols;
 *   solrelohi has space for the second highest doubles of the real parts
 *            of the solution, of size ncols;
 *   solrehilo has space for the second lowest doubles of the real parts
 *            of the solution;
 *   solrelolo has space for the lowest doubles of the real parts
 *            of the solution;
 *   solimhihi has space for the highest doubles of the imaginary parts
 *            of the solution;
 *   solimlohi has space for the second highest doubles of the imaginary parts
 *            of the solution;
 *   solimhilo has space for the second lowest doubles of the imaginary parts
 *            of the solution;
 *   solimlolo has space for the lowest doubles of the imaginary parts
 *            of the solution;
 *   wrkvecrehihi is work space for nrows elements;
 *   wrkvecrelohi is work space for nrows elements;
 *   wrkvecrehilo is work space for nrows elements;
 *   wrkvecrelolo is work space for nrows elements;
 *   wrkvecimhihi is work space for nrows elements;
 *   wrkvecimlohi is work space for nrows elements;
 *   wrkvecimhilo is work space for nrows elements;
 *   wrkvecimlolo is work space for nrows elements.
 *
 * ON RETURN :
 *   solrehi are the highest doubles of the real parts of the solution;
 *   solrehi are the second highest doubles of the real parts of the solution;
 *   solrelo are the second lowest doubles of the real parts of the solution;
 *   solrelo are the lowest doubles of the real parts of the solution;
 *   solimhi are the highest doubles of the imaginary parts of the solution;
 *   solimhi are the second highest doubles of the imaginary parts of sol;
 *   solimlo are the second lowest doubles of the imaginary parts of sol;
 *   solimlo are the lowest doubles of the imaginary parts of the solution;
 *   wrkvecrehihi has the highest doubles of the real parts of Q^H*rhs;
 *   wrkvecrelohi has the second highest doubles of the real parts of Q^H*rhs;
 *   wrkvecrehilo has the second lowest doubles of the real parts of Q^H*rhs;
 *   wrkvecrelolo has the lowest doubles of the real parts of Q^H*rhs;
 *   wrkvecimhihi has the highest doubles of the imaginary parts of Q^H*rhs;
 *   wrkvecimlohi has the second highest doubles of the imaginary parts
 *            of Q^H*rhs;
 *   wrkvecimhilo has the second lowest doubles of the imaginary parts 
 *            of Q^H*rhs;
 *   wrkvecimlolo has the lowest doubles of the imaginary parts of Q^H*rhs. */

#endif
