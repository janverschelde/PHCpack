/* The file dbl_tabs_kernels.h specifies functions for the
 * tiled accelerated back substitution in double precision. */

#ifndef __dbl_tabs_kernels_h__
#define __dbl_tabs_kernels_h__

#define d_shmemsize 256

__global__ void dbl_invert_upper ( int dim, double *U, double *invU );
/*
 * DESCRIPTION :
 *   Computes the inverse of an upper triangular matrix.
 *   The matrices are linearized in columnwise fashion.
 *
 * ON ENTRY :
 *   dim      dimension of the upper triangular matrix U;
 *   U        an upper triangular matrix of dimension dim;
 *   invU     space allocated for a matrix of dimension dim.
 *
 * ON RETURN :
 *   invU     the inverse of the matrix U. */

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
