/* The file dbl_tabs_kernels.h specifies functions for the
 * tiled accelerated back substitution in double precision. */

#ifndef __dbl_tabs_kernels_h__
#define __dbl_tabs_kernels_h__

#define d_shmemsize 256

__global__ void dbl_small_invert_upper ( int dim, double *U, double *invU );
/*
 * DESCRIPTION :
 *   Computes the inverse of an upper triangular matrix.
 *   The U matrix is stored in columnwise fashion,
 *   as the row-by-row computation of the inverse invU
 *   applies a column-by-column load of U.
 *
 * REQUIRED : dim <= 16.
 *   Because the inverse is stored entirely in shared memory,
 *   the dimension dim is limited to 16 = 2^4, as 16^2 = 256,
 *   the upper limit on the shared memory, d_shmemsize.
 *
 * ON ENTRY :
 *   dim      dimension of the upper triangular matrix U;
 *   U        an upper triangular matrix stored column wise;
 *   invU     space allocated for a matrix of dimension dim.
 *
 * ON RETURN :
 *   invU     the inverse of the matrix U, stored row wise. */

void test_dbl_small_invert_upper ( int dim, double *U, double *invU );
/*
 * DESCRIPTION :
 *   Runs the same code as dbl_small_invert_upper,
 *   but in a serialized manner, one column after the other,
 *   with many print statements to verified the correctness.
 *   The parameters are the same as dbl_small_invert_upper. */

void GPU_dbl_upper_inverse ( int dim, double **U, double **invU );
/*
 * DESCRIPTION :
 *   Calls the kernel to invert the upper triangular matrix U.
 *   The matrices are stored in the conventional rowwise fashion.
 *
 * ON ENTRY :
 *   dim      dimension of the upper triangular matrix U;
 *   U        an upper triangular matrix of dimension dim;
 *   invU     space allocated for a matrix of dimension dim.
 *
 * ON RETURN :
 *   invU     the inverse of the matrix U. */

#endif
