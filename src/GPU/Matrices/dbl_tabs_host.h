/* The file dbl_tabs_host.h specifies functions on the host for the
 * tiled accelerated back substitution in double precision. */

#ifndef __dbl_tabs_host_h__
#define __dbl_tabs_host_h__

void CPU_dbl_backsubs ( int dim, double **U, double *b, double *x );
/*
 * DESCRIPTION :
 *   Applies back substitution to solve an upper triangular system.
 *
 * ON ENTRY :
 *   dim      dimension of the upper triangular matrix U
 *            and the vectors b and x;
 *   U        an upper triangular matrix of dimension dim;
 *   b        the right hand side of the linear system;
 *   x        space allocated for dim doubles.
 *
 * ON RETURN :
 *   x        the solution to the system U*x = b. */

void CPU_dbl_upper_inverse ( int dim, double **U, double **invU );
/*
 * DESCRIPTION :
 *   Computes the inverse of an upper triangular matrix.
 *
 * ON ENTRY :
 *   dim      dimension of the upper triangular matrix U;
 *   U        an upper triangular matrix of dimension dim;
 *   invU     space allocated for a matrix of dimension dim.
 *
 * ON RETURN :
 *   invU     the inverse of the matrix U. */

void CPU_dbl_matmatmul ( int dim, double **A, double **F );
/*
 * DESCRIPTION :
 *   Replaces A with the product F*A.
 *
 * ON ENTRY :
 *   dim      number of rows and columns in A and F;
 *   A        a matrix of dimension dim;
 *   F        a matrix of dimension dim.
 *
 * ON RETURN :
 *   A        the product of F with A. */

void CPU_dbl_upper_tiled_solver
 ( int dim, int szt, int nbt, double **U, double *b, double *x );
/*
 * DESCRIPTION :
 *   Solves an upper triangular system with a tiled algorithm.
 *
 * ON ENTRY :
 *   dim      dimension of the upper triangular matrix U;
 *   szt      size of each tile;
 *   nbt      number of tiles, dim = szt*nbt;
 *   U        an upper triangular matrix of dimension dim;
 *   b        the right hand side of the linear system;
 *   x        space allocated for dim doubles.
 *
 * ON RETURN :
 *   x        the solution to the system U*x = b. */

#endif
