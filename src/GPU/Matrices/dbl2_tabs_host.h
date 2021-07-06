/* The file dbl2_tabs_host.h specifies functions on the host for the
 * tiled accelerated back substitution in double double precision. */

#ifndef __dbl2_tabs_host_h__
#define __dbl2_tabs_host_h__

void CPU_dbl2_backsubs
 ( int dim, double **Uhi, double **Ulo, double *bhi, double *blo,
   double *xhi, double *xlo );
/*
 * DESCRIPTION :
 *   Applies back substitution to solve an upper triangular system.
 *
 * ON ENTRY :
 *   dim      dimension of the upper triangular matrix U
 *            and the vectors b and x;
 *   Uhi      high doubles of an upper triangular matrix of dimension dim;
 *   Ulo      low doubles of an upper triangular matrix of dimension dim;
 *   bhi      high doubles of the right hand side of the linear system;
 *   blo      low doubles of the right hand side of the linear system;
 *   xhi      space allocated for dim doubles;
 *   xlo      space allocated for dim doubles.
 *
 * ON RETURN :
 *   xhi      high doubles of the solution to the system U*x = b;
 *   xlo      low doubles of the solution to the system U*x = b. */

void CPU_dbl2_upper_inverse
 ( int dim, double **Uhi, double **Ulo, double **invUhi, double **invUlo );
/*
 * DESCRIPTION :
 *   Computes the inverse of an upper triangular matrix.
 *
 * ON ENTRY :
 *   dim      dimension of the upper triangular matrix U;
 *   Uhi      high doubles of an upper triangular matrix of dimension dim;
 *   Ulo      low doubles of an upper triangular matrix of dimension dim;
 *   invUhi   space allocated for a matrix of dimension dim;
 *   invUlo   space allocated for a matrix of dimension dim.
 *
 * ON RETURN :
 *   invUhi   high doubles of the inverse of the matrix U;
 *   invUlo   low doubles of the inverse of the matrix U. */

void CPU_dbl2_matmatmul
 ( int dim, double **Ahi, double **Alo, double **Fhi, double **Flo );
/*
 * DESCRIPTION :
 *   Replaces A with the product F*A.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in A and F;
 *   Ahi      high doubles of a matrix of dimension dim;
 *   Alo      low doubles of a matrix of dimension dim;
 *   Fhi      high doubles of a matrix of dimension dim;
 *   Flo      low doubles of a matrix of dimension dim.
 *
 * ON RETURN :
 *   A        high doubles of the product of F with A;
 *   A        low doubles of the product of F with A. */

void CPU_dbl2_upper_tiled_solver
 ( int dim, int szt, int nbt, double **Uhi, double **Ulo,
   double *bhi, double *blo, double *xhi, double *xlo );
/*
 * DESCRIPTION :
 *   Solves an upper triangular system with a tiled algorithm.
 *
 * ON ENTRY :
 *   dim      dimension of the upper triangular matrix U;
 *   szt      size of each tile;
 *   nbt      number of tiles, dim = szt*nbt;
 *   Uhi      high doubles of an upper triangular matrix of dimension dim;
 *   Ulo      low doubles of an upper triangular matrix of dimension dim;
 *   bhi      high doubles of the right hand side of the linear system;
 *   blo      low doubles of the right hand side of the linear system;
 *   xhi      space allocated for dim doubles;
 *   xlo      space allocated for dim doubles.
 *
 * ON RETURN :
 *   xhi      high doubles of the solution to the system U*x = b;
 *   xlo      low doubles of the solution to the system U*x = b. */

#endif
