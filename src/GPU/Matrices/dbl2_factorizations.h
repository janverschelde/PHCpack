/* The file dbl2_factorizations.h specifies functions to factor matrices
 * in double double precision. */

#ifndef __dbl_factorizations_h__
#define __dbl_factorizations_h__

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

void CPU_dbl2_factors_lusolve
 ( int dim, double **Ahi, double **Alo, int *pivots,
   double *bhi, double *blo, double *xhi, double *xlo );
/*
 * DESCRIPTION :
 *   Does an inplace LU factorization with pivoting on the matrix A,
 *   to solve the system A*x = b.
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

#endif
