/* The file dbl_tabs_kernels.h specifies functions for the
 * tiled accelerated back substitution in double precision. */

#ifndef __dbl_tabs_kernels_h__
#define __dbl_tabs_kernels_h__

#define tabsd_shmemsize 256

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

__global__ void cmplx_small_invert_upper
 ( int dim, double *Ure, double *Uim, double *invUre, double *invUim );
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
 *   Ure      real parts of an upper triangular matrix stored column wise;
 *   Uim      imaginary parts of an upper triangular matrix;
 *   invUre   space allocated for a matrix of dimension dim;
 *   invUim   space allocated for a matrix of dimension dim.
 *
 * ON RETURN :
 *   invUre   real parts of the inverse of the matrix U, stored row wise;
 *   invUim   imaginary parts of the inverse of the matrix U,
 *            also stored row wise. */

void test_dbl_small_invert_upper ( int dim, double *U, double *invU );
/*
 * DESCRIPTION :
 *   Runs the same code as dbl_small_invert_upper,
 *   but in a serialized manner, one column after the other,
 *   with many print statements to verified the correctness.
 *   The parameters are the same as dbl_small_invert_upper. */

__global__ void dbl_medium_invert_upper ( int dim, double *U, double *invU );
/*
 * DESCRIPTION :
 *   Computes the inverse of an upper triangular matrix.
 *   The U matrix is stored in columnwise fashion,
 *   as the row-by-row computation of the inverse invU
 *   applies a column-by-column load of U.
 *
 * REQUIRED : dim <= 256.
 *   Because the columns of U are loaded entirely into shared memory
 *   and the rows of the inverses are computed first entirely in
 *   shared memory before storing, the dimension dim is limited 
 *   to 256, the upper limit on the shared memory, d_shmemsize.
 *
 * ON ENTRY :
 *   dim      dimension of the upper triangular matrix U;
 *   U        an upper triangular matrix stored column wise;
 *   invU     space allocated for a matrix of dimension dim.
 *
 * ON RETURN :
 *   invU     the inverse of the matrix U, stored row wise. */

__global__ void cmplx_medium_invert_upper
 ( int dim, double *Ure, double *Uim, double *invUre, double *invUim );
/*
 * DESCRIPTION :
 *   Computes the inverse of an upper triangular matrix.
 *   The U matrix is stored in columnwise fashion,
 *   as the row-by-row computation of the inverse invU
 *   applies a column-by-column load of U.
 *
 * REQUIRED : dim <= 256.
 *   Because the columns of U are loaded entirely into shared memory
 *   and the rows of the inverses are computed first entirely in
 *   shared memory before storing, the dimension dim is limited 
 *   to 256, the upper limit on the shared memory, d_shmemsize.
 *
 * ON ENTRY :
 *   dim      dimension of the upper triangular matrix U;
 *   Ure      real parts of an upper triangular matrix stored column wise;
 *   Uim      imaginary parts of an upper triangular matrix,
 *            also store column wise;
 *   invUre   space allocated for a matrix of dimension dim;
 *   invUim   space allocated for a matrix of dimension dim.
 *
 * ON RETURN :
 *   invUre   real parts of the inverse of the matrix U, stored row wise;
 *   invUim   imaginary parts of the inverse of the matrix U,
 *            also stored row wise. */

__global__ void  dbl_invert_tiles ( int dim, double *U, double *invU );
/*
 * DESCRIPTION :
 *   Replaces the columns of the tiles with the rows of the inverses.
 *   The number of blocks equals the number of tiles in U.
 *   The number of threads per block equals the dimension of each tile.
 *
 * REQUIRED : dim <= 256 = d_shmemsize.
 *
 * ON ENTRY :
 *   dim      the dimension of each tile;
 *   U        columns of all tiles on the diagonal 
 *            of an upper triangular matrix;
 *   invU     space allocated for the inverse of all tiles in U.
 *
 * ON RETURN :
 *   invU     rows of the inverse of the tiles in U. */

__global__ void  cmplx_invert_tiles
 ( int dim, double *Ure, double *Uim, double *invUre, double *invUim );
/*
 * DESCRIPTION :
 *   Replaces the columns of the tiles with the rows of the inverses.
 *   The number of blocks equals the number of tiles in U.
 *   The number of threads per block equals the dimension of each tile.
 *
 * REQUIRED : dim <= 256 = d_shmemsize.
 *
 * ON ENTRY :
 *   dim      the dimension of each tile;
 *   Ure      real parts of the columns of all tiles on the diagonal 
 *            of an upper triangular matrix;
 *   Uim      imaginary parts of the columns of all tiles on the diagonal 
 *            of an upper triangular matrix;
 *   invUre   space allocated for the real parts of
 *            the inverse of all tiles in U;
 *   invUim   space allocated for the imaginary parts of
 *            the inverse of all tiles in U.
 *
 * ON RETURN :
 *   invUre   real parts of rows of the inverse of the tiles in U;
 *   invU     imaginary parts of rows of the inverse of the tiles in U. */

__global__ void dbl_multiply_inverse
 ( int dim, int idx, double *invU, double *w );
/*
 * DESCRIPTION :
 *   Replaces b with the product of the inverse tile in U.
 *
 * ON ENTRY :
 *   dim      the dimension of each tile;
 *   idx      index of the diagonal tile;
 *   invU     contains the inverse of the diagonal tiles;
 *   w        right hand side vector.
 *
 * ON RETURN :
 *   w        product of the proper inverse of the diagonal tile
 *            defines by the index idx with the w on input. */

__global__ void cmplx_multiply_inverse
 ( int dim, int idx, double *invUre, double *invUim,
   double *wre, double *wim );
/*
 * DESCRIPTION :
 *   Replaces b with the product of the inverse tile in U.
 *
 * ON ENTRY :
 *   dim      the dimension of each tile;
 *   idx      index of the diagonal tile;
 *   invUre   real parts of the inverse of the diagonal tiles;
 *   invUim   imaginary parts of the inverse of the diagonal tiles;
 *   wre      real parts of the right hand side vector;
 *   wim      imaginary parts of the right hand side vector.
 *
 * ON RETURN :
 *   wre      real parts of the product of the inverse of the diagonal
 *            tile defined by the index idx with the w on input;
 *   wim      imaginary parts of the product of the inverse of the diagonal
 *            tile defined by the index idx with the w on input. */

__global__ void dbl_back_substitute
 ( int dim, int idx, double *U, double *w );
/*
 * DESCRIPTION :
 *   Updates the right hand side vector subtracting the solution
 *   defined by idx, multiplied with the corresponding rows in U.
 *
 * ON ENTRY :
 *   dim      dimension of each tile;
 *   idx      index of the solution tile in the multiplication;
 *   U        contains tiles to multiply the solution with;
 *   w        current right hand side vector.
 *
 * ON RETURN :
 *   w        updated right hand side vector. */

__global__ void cmplx_back_substitute
 ( int dim, int idx, double *Ure, double *Uim, double *wre, double *wim );
/*
 * DESCRIPTION :
 *   Updates the right hand side vector subtracting the solution
 *   defined by idx, multiplied with the corresponding rows in U.
 *
 * ON ENTRY :
 *   dim      dimension of each tile;
 *   idx      index of the solution tile in the multiplication;
 *   Ure      real parts of tiles to multiply the solution with;
 *   Uim      imaginary parts of tiles to multiply the solution with;
 *   wre      real parts of the current right hand side vector;
 *   wim      imaginary parts of the current right hand side vector.
 *
 * ON RETURN :
 *   wre      real parts of the updated right hand side vector;
 *   wim      imaginary parts of the updated right hand side vector. */

void GPU_dbl_upper_inverse
 ( int dim, double **U, double **invU, double *lapms, double *walltimesec );
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
 *   invU     the inverse of the matrix U;
 *   lapms    elapsed time spent by the kernels;
 *   walltimesec is the elapsed wall clock computation time. */

void GPU_cmplx_upper_inverse
 ( int dim, double **Ure, double **Uim, double **invUre, double **invUim,
   double *lapms, double *walltimesec );
/*
 * DESCRIPTION :
 *   Calls the kernel to invert the upper triangular matrix U.
 *   The matrices are stored in the conventional rowwise fashion.
 *
 * ON ENTRY :
 *   dim      dimension of the upper triangular matrix U;
 *   Ure      real parts of an upper triangular matrix of dimension dim;
 *   Uim      imaginary parts of an upper triangular matrix of dimension dim;
 *   invUre   space allocated for a matrix of dimension dim;
 *   invUim   space allocated for a matrix of dimension dim.
 *
 * ON RETURN :
 *   invUre   real parts of the inverse of the matrix U;
 *   invUim   imaginary parts of the inverse of the matrix U;
 *   lapms    elapsed time spent by the kernels;
 *   walltimesec is the elapsed wall clock computation time. */

void GPU_dbl_upper_tiled_solver
 ( int dim, int szt, int nbt, double **U, double *b, double *x,
   double *invlapms, double *mullapms, double *sublapms, double *totlapms,
   double *walltimesec,
   long long int *addcnt, long long int *mulcnt, long long int *divcnt );
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
 *   x        the solution to the system U*x = b;
 *   invlapms is the elapsed time spent by the kernel to invert a tile;
 *   mullapms is the elapsed time spent by the kernel to multiply
 *            with the inversed diagonal tile;
 *   sublapms is the elapsed time spent by the kernel for back substitution;
 *   totlapms is the total elapsed time spent by all kernels;
 *   walltimesec is the elapsed wall clock computation time;
 *   addcnt   counts the number of additions and subtractions;
 *   mulcnt   counts the number of multiplications;
 *   divcnt   counts the number of divisions. */

void GPU_cmplx_upper_tiled_solver
 ( int dim, int szt, int nbt, double **Ure, double **Uim,
   double *bre, double *bim, double *xre, double *xim,
   double *invlapms, double *mullapms, double *sublapms, double *totlapms,
   double *walltimesec,
   long long int *addcnt, long long int *mulcnt, long long int *divcnt );
/*
 * DESCRIPTION :
 *   Solves an upper triangular system with a tiled algorithm.
 *
 * ON ENTRY :
 *   dim      dimension of the upper triangular matrix U;
 *   szt      size of each tile;
 *   nbt      number of tiles, dim = szt*nbt;
 *   Ure      real parts of an upper triangular matrix of dimension dim;
 *   Uim      imaginary parts of an upper triangular matrix of dimension dim;
 *   bre      real parts of the right hand side of the linear system;
 *   bim      imaginary parts of the right hand side of the linear system;
 *   xre      space allocated for dim doubles;
 *   xim      space allocated for dim doubles.
 *
 * ON RETURN :
 *   xre      real parts of the solution to the system U*x = b;
 *   xim      imaginary parts of the solution to the system U*x = b;
 *   invlapms is the elapsed time spent by the kernel to invert a tile;
 *   mullapms is the elapsed time spent by the kernel to multiply
 *            with the inversed diagonal tile;
 *   sublapms is the elapsed time spent by the kernel for back substitution;
 *   totlapms is the total elapsed time spent by all kernels;
 *   walltimesec is the elapsed wall clock computation time;
 *   addcnt   counts the number of additions and subtractions;
 *   mulcnt   counts the number of multiplications;
 *   divcnt   counts the number of divisions. */

#endif
