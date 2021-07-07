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

#endif