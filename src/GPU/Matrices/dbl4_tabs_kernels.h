/* The file dbl4_tabs_kernels.h specifies functions for the
 * tiled accelerated back substitution in quad double precision. */

#ifndef __dbl4_tabs_kernels_h__
#define __dbl4_tabs_kernels_h__

#define tabsqd_shmemsize 256

__global__ void dbl4_small_invert_upper 
( int dim, double *Uhihi, double *Ulohi, double *Uhilo, double *Ulolo,
  double *invUhihi, double *invUlohi, double *invUhilo, double *invUlolo );
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
 *   the upper limit on the shared memory, tabsqd_shmemsize.
 *
 * ON ENTRY :
 *   dim      dimension of the upper triangular matrix U;
 *   Uhihi    highest doubles of U stored column wise;
 *   Ulohi    second highest doubles of U stored column wise;
 *   Uhilo    second lowest doubles of U stored column wise;
 *   Ulolo    lowest doubles of U stored column wise;
 *   invUhihi has space allocated for a matrix of dimension dim;
 *   invUlohi has space allocated for a matrix of dimension dim;
 *   invUhilo has space allocated for a matrix of dimension dim;
 *   invUlolo has space allocated for a matrix of dimension dim.
 *
 * ON RETURN :
 *   invUhihi has the highest doubles of the inverse of the matrix U,
 *            stored row wise;
 *   invUlohi has the second highest doubles of the inverse of the matrix U,
 *            stored row wise;
 *   invUhilo has the second lowest doubles the inverse of the matrix U,
 *            stored row wise;
 *   invUlolo are the lowest doubles the inverse of the matrix U,
 *            stored row wise. */

__global__ void cmplx4_small_invert_upper
 ( int dim,
   double *Urehihi, double *Urelohi, double *Urehilo, double *Urelolo,
   double *Uimhihi, double *Uimlohi, double *Uimhilo, double *Uimlolo,
   double *invUrehihi, double *invUrelohi,
   double *invUrehilo, double *invUrelolo,
   double *invUimhihi, double *invUimlohi,
   double *invUimhilo, double *invUimlolo );
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
 *   the upper limit on the shared memory, tabsqd_shmemsize.
 *
 * ON ENTRY :
 *   dim      dimension of the upper triangular matrix U;
 *   Urehihi  highest doubles of the real parts of U stored column wise;
 *   Urelohi  second highest doubles of the real parts of U;
 *   Urehilo  second lowest doubles of the real parts of U;
 *   Urelolo  lowest doubles of the real parts of U stored column wise;
 *   Uimhihi  highest doubles of the imaginary parts of U;
 *   Uimlohi  second highest doubles of the imaginary parts of U;
 *   Uimhilo  second lowest doubles of the imaginary parts of U;
 *   Uimlolo  lowest doubles of the imaginary parts of U;
 *   invUrehihi has space allocated for a matrix of dimension dim;
 *   invUrelohi has space allocated for a matrix of dimension dim;
 *   invUrehilo has space allocated for a matrix of dimension dim;
 *   invUrelolo has space allocated for a matrix of dimension dim;
 *   invUimhihi has space allocated for a matrix of dimension dim;
 *   invUimlohi has space allocated for a matrix of dimension dim;
 *   invUimhilo has space allocated for a matrix of dimension dim;
 *   invUimlolo has space allocated for a matrix of dimension dim.
 *
 * ON RETURN :
 *   invUrehihi has the highest doubles of the real parts of the inverse,
 *            stored row wise;
 *   invUrelohi has the second highest doubles of the real parts
 *            of the inverse, stored row wise;
 *   invUrehilo has the second lowest doubles of the real parts
 *            of the inverse, stored row wise;
 *   invUrelolo has the lowest doubles of the real parts of the inverse,
 *            stored row wise;
 *   invUimhihi has the highest doubles of the imaginary parts of the inverse,
 *            also stored row wise;
 *   invUimlohi has the second highest doubles of the imaginary parts
 *            of the inverse, also stored row wise;
 *   invUimhilo has the second lowest doubles of the imaginary parts
 *            of the inverse, also stored row wise.
 *   invUimlolo has the lowest doubles of the imaginary parts of the inverse,
 *            also stored row wise. */

__global__ void dbl4_medium_invert_upper
 ( int dim, double *Uhihi, double *Ulohi, double *Uhilo, double *Ulolo,
   double *invUhihi, double *invUlohi, double *invUhilo, double *invUlolo);
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
 *   to 256, the upper limit on the shared memory, tabsqd_shmemsize.
 *
 * ON ENTRY :
 *   dim      dimension of the upper triangular matrix U;
 *   Uhihi    highest doubles of U stored column wise;
 *   Ulohi    second highest doubles of U stored column wise;
 *   Uhilo    second lowest doubles of U stored column wise;
 *   Ulolo    lowest doubles of U stored column wise;
 *   invUhihi has space allocated for a matrix of dimension dim;
 *   invUlohi has space allocated for a matrix of dimension dim;
 *   invUhilo has space allocated for a matrix of dimension dim;
 *   invUlolo has space allocated for a matrix of dimension dim.
 *
 * ON RETURN :
 *   invUhihi has the highest doubles of the inverse of the matrix U,
 *            stored row wise;
 *   invUlohi has the second highest doubles of the inverse of the matrix U,
 *            stored row wise;
 *   invUhilo has the second lowest doubles the inverse of the matrix U,
 *            stored row wise;
 *   invUlolo are the lowest doubles the inverse of the matrix U,
 *            stored row wise. */

__global__ void cmplx4_medium_invert_upper
 ( int dim, 
   double *Urehihi, double *Urelohi, double *Urehilo, double *Urelolo,
   double *Uimhihi, double *Uimlohi, double *Uimhilo, double *Uimlolo,
   double *invUrehihi, double *invUrelohi,
   double *invUrehilo, double *invUrelolo,
   double *invUimhihi, double *invUimlohi,
   double *invUimhilo, double *invUimlolo );
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
 *   to 256, the upper limit on the shared memory, tabsqd_shmemsize.
 *
 * ON ENTRY :
 *   dim      dimension of the upper triangular matrix U;
 *   Urehihi  highest doubles of the real parts of U stored column wise;
 *   Urelohi  second highest doubles of the real parts of U;
 *   Urehilo  second lowest doubles of the real parts of U;
 *   Urelolo  lowest doubles of the real parts of U stored column wise;
 *   Uimhihi  highest doubles of the imaginary parts of U;
 *   Uimlohi  second highest doubles of the imaginary parts of U;
 *   Uimhilo  second lowest doubles of the imaginary parts of U;
 *   Uimlolo  lowest doubles of the imaginary parts of U;
 *   invUrehihi has space allocated for a matrix of dimension dim;
 *   invUrelohi has space allocated for a matrix of dimension dim;
 *   invUrehilo has space allocated for a matrix of dimension dim;
 *   invUrelolo has space allocated for a matrix of dimension dim;
 *   invUimhihi has space allocated for a matrix of dimension dim;
 *   invUimlohi has space allocated for a matrix of dimension dim;
 *   invUimhilo has space allocated for a matrix of dimension dim;
 *   invUimlolo has space allocated for a matrix of dimension dim.
 *
 * ON RETURN :
 *   invUrehihi has the highest doubles of the real parts of the inverse,
 *            stored row wise;
 *   invUrelohi has the second highest doubles of the real parts
 *            of the inverse, stored row wise;
 *   invUrehilo has the second lowest doubles of the real parts
 *            of the inverse, stored row wise;
 *   invUrelolo has the lowest doubles of the real parts of the inverse,
 *            stored row wise;
 *   invUimhihi has the highest doubles of the imaginary parts of the inverse,
 *            also stored row wise;
 *   invUimlohi has the second highest doubles of the imaginary parts
 *            of the inverse, also stored row wise;
 *   invUimhilo has the second lowest doubles of the imaginary parts
 *            of the inverse, also stored row wise.
 *   invUimlolo has the lowest doubles of the imaginary parts of the inverse,
 *            also stored row wise. */

__global__ void  dbl4_invert_tiles
 ( int dim, double *Uhihi, double *Ulohi, double *Uhilo, double *Ulolo,
   double *invUhihi, double *invUlohi, double *invUhilo, double *invUlolo);
/*
 * DESCRIPTION :
 *   Replaces the columns of the tiles with the rows of the inverses.
 *   The number of blocks equals the number of tiles in U.
 *   The number of threads per block equals the dimension of each tile.
 *
 * REQUIRED : dim <= 256 = tabsqd_shmemsize.
 *
 * ON ENTRY :
 *   dim      the dimension of each tile;
 *   Uhihi    highest doubles of columns of all tiles on the diagonal 
 *            of an upper triangular matrix;
 *   Ulohi    second highest doubles of columns of all tiles on the diagonal 
 *            of an upper triangular matrix;
 *   Uhilo    second lowest doubles of columns of all tiles on the diagonal 
 *            of an upper triangular matrix;
 *   Ulolo    lowest doubles of columns of all tiles on the diagonal 
 *            of an upper triangular matrix;
 *   invUhihi has space allocated for the inverse of all tiles in U;
 *   invUlohi has space allocated for the inverse of all tiles in U;
 *   invUhilo has space allocated for the inverse of all tiles in U;
 *   invUlolo has space allocated for the inverse of all tiles in U.
 *
 * ON RETURN :
 *   invUhihi are the highest doubles of the inverse of the tiles in U;
 *   invUlohi are the second highest doubles of the inverse of the tiles in U;
 *   invUhilo are the second lowest doubles of the inverse of the tiles in U;
 *   invUlolo are the lowest doubles of the inverse of the tiles in U. */

__global__ void  cmplx4_invert_tiles
 ( int dim, 
   double *Urehihi, double *Urelohi, double *Urehilo, double *Urelolo,
   double *Uimhihi, double *Uimlohi, double *Uimhilo, double *Uimlolo,
   double *invUrehihi, double *invUrelohi,
   double *invUrehilo, double *invUrelolo,
   double *invUimhihi, double *invUimlohi,
   double *invUimhilo, double *invUimlolo );
/*
 * DESCRIPTION :
 *   Replaces the columns of the tiles with the rows of the inverses.
 *   The number of blocks equals the number of tiles in U.
 *   The number of threads per block equals the dimension of each tile.
 *
 * REQUIRED : dim <= 256 = tabsqd_shmemsize.
 *
 * ON ENTRY :
 *   dim      the dimension of each tile;
 *   Urehihi  highest doubles of the real parts of the columns
 *            of all tiles on the diagonal of an upper triangular matrix;
 *   Urelohi  second highest doubles of the real parts of the columns
 *            of all tiles on the diagonal of an upper triangular matrix;
 *   Urehilo  second lowest doubles of the real parts of the columns
 *            of all tiles on the diagonal of an upper triangular matrix;
 *   Urelolo  lowest doubles of the real parts of the columns
 *            of all tiles on the diagonal of an upper triangular matrix;
 *   Uimhihi  highest doubles of the imaginary parts of the columns
 *            of all tiles on the diagonal of an upper triangular matrix;
 *   Uimlohi  second highest doubles of the imaginary parts of the columns
 *            of all tiles on the diagonal of an upper triangular matrix;
 *   Uimhilo  second lowest doubles of the imaginary parts of the columns
 *            of all tiles on the diagonal of an upper triangular matrix;
 *   Uimlolo  lowest doubles of the imaginary parts of the columns
 *            of all tiles on the diagonal of an upper triangular matrix;
 *   invUrehihi has space allocated for the highest doubles of
 *            the real parts of the inverse of all tiles in U;
 *   invUrelohi has space allocated for the second highest doubles of
 *            the real parts of the inverse of all tiles in U;
 *   invUrehilo has space allocated for the second lowest doubles of
 *            the real parts of the inverse of all tiles in U;
 *   invUrelolo has space allocated for the lowest doubles of
 *            the real parts of the inverse of all tiles in U;
 *   invUimhihi has space allocated for the highest doubles of
 *            the imaginary parts of the inverse of all tiles in U;
 *   invUimlohi has space allocated for the second highest doubles of
 *            the imaginary parts of the inverse of all tiles in U;
 *   invUimhilo has space allocated for the second lowest doubles of
 *            the imaginary parts of the inverse of all tiles in U.
 *   invUimlolo has space allocated for the lowest doubles of
 *            the imaginary parts of the inverse of all tiles in U.
 *
 * ON RETURN :
 *   invUrehihi has the highest doubles of the real parts
 *            of rows of the inverse tiles;
 *   invUrelohi has the second highest doubles of the real parts
 *            of rows of the inverse tiles;
 *   invUrehilo has the second lowest doubles of the real parts
 *            of rows of the inverse tiles;
 *   invUrelolo has the lowest doubles of the real parts
 *            of rows of the inverse tiles;
 *   invUimhihi has the highest doubles of the imaginary parts
 *            of rows of the inverse tiles;
 *   invUimlohi has the second highest doubles of the imaginary parts
 *            of rows of the inverse tiles;
 *   invUimhilo has the second lowest doubles of the imaginary parts
 *            of rows of the inverse tiles;
 *   invUimlolo has the lowest doubles of the imaginary parts
 *            of rows of the inverse tiles. */

__global__ void dbl4_multiply_inverse
 ( int dim, int idx,
   double *invUhihi, double *invUlohi, double *invUhilo, double *invUlolo,
   double *whihi, double *wlohi, double *whilo, double *wlolo );
/*
 * DESCRIPTION :
 *   Replaces b with the product of the inverse tile in U.
 *
 * ON ENTRY :
 *   dim      the dimension of each tile;
 *   idx      index of the diagonal tile;
 *   invUhihi are the highest doubles of the inverse of the diagonal tiles;
 *   invUlohi are the second highest doubles of the inverse diagonal tiles;
 *   invUhilo are the second lowest doubles of the inverse diagonal tiles;
 *   invUlolo are the  lowest doubles of the inverse of the diagonal tiles;
 *   whihi    highest doubles of the right hand side vector;
 *   wlohi    second highest doubles of the right hand side vector;
 *   whilo    second lowest doubles of the right hand side vector.
 *   wlolo    lowest doubles of the right hand side vector.
 *
 * ON RETURN :
 *   whihi    highest doubles of the product of the inverse of the
 *            diagonal tile defined by the index idx with the w on input;
 *   wlohi    second highest doubles of the product of the inverse of the
 *            diagonal tile defined by the index idx with the w on input;
 *   whilo    second lowest doubles of the product of the inverse of the
 *            diagonal tile defined by the index idx with the w on input;
 *   wlolo    lowest doubles of the product of the inverse of the
 *            diagonal tile defined by the index idx with the w on input. */

__global__ void cmplx4_multiply_inverse
 ( int dim, int idx,
   double *invUrehihi, double *invUrelohi,
   double *invUrehilo, double *invUrelolo,
   double *invUimhihi, double *invUimlohi,
   double *invUimhilo, double *invUimlolo,
   double *wrehihi, double *wrelohi, double *wrehilo, double *wrelolo,
   double *wimhihi, double *wimlohi, double *wimhilo, double *wimlolo );
/*
 * DESCRIPTION :
 *   Replaces b with the product of the inverse tile in U.
 *
 * ON ENTRY :
 *   dim      the dimension of each tile;
 *   idx      index of the diagonal tile;
 *   invUrehihi are the highest doubles of the real parts of
 *            the inverse of the diagonal tiles;
 *   invUrelohi are the second highest doubles of the real parts of
 *            the inverse of the diagonal tiles;
 *   invUrehilo are the second lowest doubles of the real parts of
 *            the inverse of the diagonal tiles;
 *   invUrelolo are the lowest doubles of the real parts of
 *            the inverse of the diagonal tiles;
 *   invUimhihi are the highest doubles of the imaginary parts of
 *            the inverse of the diagonal tiles;
 *   invUimlohi are the second highest doubles of the imaginary parts of
 *            the inverse of the diagonal tiles;
 *   invUimhilo are the second lowest doubles of the imaginary parts of
 *            the inverse of the diagonal tiles;
 *   invUimlolo are the lowest doubles of the imaginary parts of
 *            the inverse of the diagonal tiles;
 *   wrehihi  highest doubles of the real parts of the right hand side b;
 *   wrelohi  second highest doubles of the real parts of b;
 *   wrehilo  second lowest doubles of the real parts of b;
 *   wrelolo  lowest doubles of the real parts of b;
 *   wimhihi  highest doubles of the imaginary parts of b;
 *   wimlohi  second highest doubles of the imaginary parts of b;
 *   wimhilo  second lowest doubles of the imaginary parts of b.
 *   wimlolo  lowest doubles of the imaginary parts of b.
 *
 * ON RETURN :
 *   wrehihi  highest doubles of the real parts of the product of the
 *            inverse of the tile defined by the index idx with the input w;
 *   wrelohi  second highest doubles of the real parts of the product
 *            of the inverse of the tile defined by the index idx;
 *   wrehilo  second lowest doubles of the real parts of the product
 *            of the inverse of the tile defined by the index idx;
 *   wrelolo  lowest doubles of the real parts of the product
 *            of the inverse of the tile defined by the index idx;
 *   wimhihi  highest doubles of the imaginary parts of the product
 *            of the inverse of the tile defined by the index idx;
 *   wimlohi  second highest doubles of the imaginary parts of the product
 *            of the inverse of the tile defined by the index idx;
 *   wimhilo  second lowest doubles of the imaginary parts of the product
 *            of the inverse of the tile defined by the index idx;
 *   wimlolo  lowest doubles of the imaginary parts of the product
 *            of the inverse of the tile defined by the index idx. */

__global__ void dbl4_back_substitute
 ( int dim, int idx, 
   double *Uhihi, double *Ulohi, double *Uhilo, double *Ulolo, 
   double *whihi, double *wlohi, double *whilo, double *wlolo );
/*
 * DESCRIPTION :
 *   Updates the right hand side vector subtracting the solution
 *   defined by idx, multiplied with the corresponding rows in U.
 *
 * ON ENTRY :
 *   dim      dimension of each tile;
 *   idx      index of the solution tile in the multiplication;
 *   Uhihi    highest doubles of tiles to multiply the solution with;
 *   Ulohi    second highest doubles of tiles to multiply the solution with;
 *   Uhilo    second lowest doubles of tiles to multiply the solution with;
 *   Ulolo    lowest doubles of tiles to multiply the solution with;
 *   whihi    highest doubles of the current right hand side vector;
 *   wlohi    second highest doubles of the current right hand side vector;
 *   whilo    second lowest doubles of the current right hand side vector;
 *   wlolo    lowest doubles of the current right hand side vector.
 *
 * ON RETURN :
 *   whihi    highest doubles of the updated right hand side vector;
 *   wlohi    second highest doubles of the updated right hand side vector;
 *   whilo    second lowest doubles of the updated right hand side vector;
 *   wlolo    lowest doubles of the updated right hand side vector. */

__global__ void cmplx4_back_substitute
 ( int dim, int idx,
   double *Urehihi, double *Urelohi, double *Urehilo, double *Urelolo,
   double *Uimhihi, double *Uimlohi, double *Uimhilo, double *Uimlolo,
   double *wrehihi, double *wrelohi, double *wrehilo, double *wrelolo,
   double *wimhihi, double *wimlohi, double *wimhilo, double *wimlolo );
/*
 * DESCRIPTION :
 *   Updates the right hand side vector subtracting the solution
 *   defined by idx, multiplied with the corresponding rows in U.
 *
 * ON ENTRY :
 *   dim      dimension of each tile;
 *   idx      index of the solution tile in the multiplication;
 *   Urehihi  highest doubles of the real parts of tiles 
 *            to multiply the solution with;
 *   Urelohi  second highest doubles of the real parts of tiles 
 *            to multiply the solution with;
 *   Urehilo  second lowest doubles of the real parts of tiles
 *            to multiply the solution with;
 *   Urelolo  lowest doubles of the real parts of tiles
 *            to multiply the solution with;
 *   Uimhihi  highest doubles of the imaginary parts of tiles
 *            to multiply the solution with;
 *   Uimlohi  second highest doubles of the imaginary parts of tiles
 *            to multiply the solution with;
 *   Uimhilo  second lowest doubles of the imaginary parts of tiles
 *            to multiply the solution with;
 *   Uimlolo  lowest doubles of the imaginary parts of tiles
 *            to multiply the solution with;
 *   wrehihi  highest doubles of the real parts of 
 *            the current right hand side vector;
 *   wrelohi  second highest doubles of the real parts of 
 *            the current right hand side vector;
 *   wrehilo  second lowest doubles of the real parts of 
 *            the current right hand side vector;
 *   wrelolo  lowest doubles of the real parts of 
 *            the current right hand side vector;
 *   wimhihi  highest doubles of the imaginary parts
 *            of the current right hand side vector;
 *   wimlohi  second highest doubles of the imaginary parts
 *            of the current right hand side vector;
 *   wimhilo  second lowest doubles of the imaginary parts
 *            of the current right hand side vector;
 *   wimlolo  lowest doubles of the imaginary parts
 *            of the current right hand side vector.
 *
 * ON RETURN :
 *   wrehihi  highest doubles of the real parts
 *            of the updated right hand side vector;
 *   wrelohi  second highest doubles of the real parts
 *            of the updated right hand side vector;
 *   wrehilo  second lowest doubles of the real parts
 *            of the updated right hand side vector;
 *   wrelolo  lowest doubles of the real parts
 *            of the updated right hand side vector;
 *   wimhihi  highest doubles of the imaginary parts
 *            of the updated right hand side vector;
 *   wimlohi  second highest doubles of the imaginary parts
 *            of the updated right hand side vector;
 *   wimhilo  second lowest doubles of the imaginary parts
 *            of the updated right hand side vector;
 *   wimlolo  lowest doubles of the imaginary parts
 *            of the updated right hand side vector. */

void GPU_dbl4_upper_inverse
 ( int dim, double **Uhihi, double **Ulohi, double **Uhilo, double **Ulolo,
   double **invUhihi, double **invUlohi, double **invUhilo, double **invUlolo,
   double *lapms, double *walltimesec );
/*
 * DESCRIPTION :
 *   Calls the kernel to invert the upper triangular matrix U.
 *   The matrices are stored in the conventional rowwise fashion.
 *
 * ON ENTRY :
 *   dim      dimension of the upper triangular matrix U;
 *   Uhihi    highest doubles of U;
 *   Ulohi    second highest doubles of U;
 *   Uhilo    second lowest doubles of U;
 *   Ulolo    lowest doubles of U;
 *   invUhihi has space allocated for a matrix of dimension dim;
 *   invUlohi has space allocated for a matrix of dimension dim;
 *   invUhilo has space allocated for a matrix of dimension dim;
 *   invUlolo has space allocated for a matrix of dimension dim.
 *
 * ON RETURN :
 *   invUhihi has the highest doubles of the inverse of U;
 *   invUlohi has the second highest doubles of the inverse of U;
 *   invUhilo has the second lowest doubles of the inverse of U;
 *   invUlolo has the lowest doubles of the inverse of U;
 *   lapms    elapsed time spent by the kernels;
 *   walltimesec is the elapsed wall clock computation time. */

void GPU_cmplx4_upper_inverse
 ( int dim,
   double **Urehihi, double **Urelohi, double **Urehilo, double **Urelolo,
   double **Uimhihi, double **Uimlohi, double **Uimhilo, double **Uimlolo,
   double **invUrehihi, double **invUrelohi,
   double **invUrehilo, double **invUrelolo,
   double **invUimhihi, double **invUimlohi, 
   double **invUimhilo, double **invUimlolo, 
   double *lapms, double *walltimesec );
/*
 * DESCRIPTION :
 *   Calls the kernel to invert the upper triangular matrix U.
 *   The matrices are stored in the conventional row wise fashion.
 *
 * ON ENTRY :
 *   dim      dimension of the upper triangular matrix U;
 *   Urehihi  highest doubles of the real parts of U;
 *   Urelohi  second highest doubles of the real parts of U;
 *   Urehilo  second lowest doubles of the real parts of U;
 *   Urelolo  lowest doubles of the real parts of U;
 *   Uimhihi  highest doubles of the imaginary parts of U;
 *   Uimlohi  second highest doubles of the imaginary parts of U;
 *   Uimhilo  second lowest doubles of the imaginary parts of U;
 *   Uimlolo  lowest doubles of the imaginary parts of U;
 *   invUrehihi has space allocated for a matrix of dimension dim;
 *   invUrelohi has space allocated for a matrix of dimension dim;
 *   invUrehilo has space allocated for a matrix of dimension dim;
 *   invUrelolo has space allocated for a matrix of dimension dim;
 *   invUimhihi has space allocated for a matrix of dimension dim;
 *   invUimlohi has space allocated for a matrix of dimension dim;
 *   invUimhilo has space allocated for a matrix of dimension dim;
 *   invUimlolo has space allocated for a matrix of dimension dim.
 *
 * ON RETURN :
 *   invUrehihi has the highest doubles of the real parts of the inverse;
 *   invUrelohi has the second highest doubles of the real parts
 *            of the inverse;
 *   invUrehilo has the second lowest doubles of the real partls
 *            of the inverse;
 *   invUrelolo has the lowest doubles of the real parts of the inverse;
 *   invUimhihi has the highest doubles of the imaginary parts of the inverse;
 *   invUimlohi has the second highest doubles of the imaginary parts
 *            of the inverse;
 *   invUimhilo has the second lowest doubles of the imaginary parts
 *            of the inverse;
 *   invUimlolo has the lowest doubles of the imaginary parts of the inverse;
 *   lapms    elapsed time spent by the kernels;
 *   walltimesec is the elapsed wall clock computation time. */

void GPU_dbl4_upper_tiled_solver
 ( int dim, int szt, int nbt,
   double **Uhihi, double **Ulohi, double **Uhilo, double **Ulolo,
   double *bhihi, double *blohi, double *bhilo, double *blolo,
   double *xhihi, double *xlohi, double *xhilo, double *xlolo,
   double *invlapms, double *mullapms, double *sublapms, double *totlapms,
   double *walltimesec,
   long long int *addcnt, double *addover,
   long long int *mulcnt, double *mulover,
   long long int *divcnt, double *divover );
/*
 * DESCRIPTION :
 *   Solves an upper triangular system with a tiled algorithm.
 *
 * ON ENTRY :
 *   dim      dimension of the upper triangular matrix U;
 *   szt      size of each tile;
 *   nbt      number of tiles, dim = szt*nbt;
 *   Uhihi    highest doubles of U;
 *   Ulohi    second highest doubles of U;
 *   Uhilo    second lowest doubles of U;
 *   Ulolo    lowest doubles of U;
 *   bhihi    highest doubles of the right hand side b;
 *   blohi    second highest doubles of the right hand side b;
 *   bhilo    second lowest doubles of the right hand side b;
 *   blolo    lowest doubles of the right hand side b;
 *   xhihi    space allocated for dim doubles;
 *   xlohi    space allocated for dim doubles;
 *   xhilo    space allocated for dim doubles;
 *   xlolo    space allocated for dim doubles.
 *
 * ON RETURN :
 *   xhihi    highest doubles of the solution to U*x = b;
 *   xlohi    second highest doubles of the solution to U*x = b;
 *   xhilo    second lowest doubles of the solution to U*x = b;
 *   xlolo    lowest doubles of the solution to U*x = b;
 *   invlapms is the elapsed time spent by the kernel to invert a tile;
 *   mullapms is the elapsed time spent by the kernel to multiply
 *            with the inversed diagonal tile;
 *   sublapms is the elapsed time spent by the kernel for back substitution;
 *   totlapms is the total elapsed time spent by all kernels;
 *   walltimesec is the elapsed wall clock computation time;
 *   addcnt   counts the number of additions and subtractions;
 *   addover  overflowed number of additions and subtractions;
 *   mulcnt   counts the number of multiplications;
 *   mulover  overflowed the number of multiplications;
 *   divcnt   counts the number of divisions;
 *   divover  overflowed number of divisions. */

void GPU_cmplx4_upper_tiled_solver
 ( int dim, int szt, int nbt,
   double **Urehihi, double **Urelohi, double **Urehilo, double **Urelolo,
   double **Uimhihi, double **Uimlohi, double **Uimhilo, double **Uimlolo,
   double *brehihi, double *brelohi, double *brehilo, double *brelolo,
   double *bimhihi, double *bimlohi, double *bimhilo, double *bimlolo,
   double *xrehihi, double *xrelohi, double *xrehilo, double *xrelolo,
   double *ximhihi, double *ximlohi, double *ximhilo, double *ximlolo,
   double *invlapms, double *mullapms, double *sublapms, double *totlapms,
   double *walltimesec,
   long long int *addcnt, double *addover,
   long long int *mulcnt, double *mulover,
   long long int *divcnt, double *divover );
/*
 * DESCRIPTION :
 *   Solves an upper triangular system with a tiled algorithm.
 *
 * ON ENTRY :
 *   dim      dimension of the upper triangular matrix U;
 *   szt      size of each tile;
 *   nbt      number of tiles, dim = szt*nbt;
 *   Urehihi  highest doubles of the real parts of U;
 *   Urelohi  second highest doubles of the real parts of U;
 *   Urehilo  second lowest doubles of the real parts of U;
 *   Urelolo  lowest doubles of the real parts of U;
 *   Uimhihi  highest doubles of the imaginary parts of U;
 *   Uimlohi  second highest doubles of the imaginary parts of U;
 *   Uimhilo  second lowest doubles of the imaginary parts of U;
 *   Uimlolo  lowest doubles of the imaginary parts of U;
 *   brehihi  highest doubles of the real parts of the right hand side b;
 *   brelohi  second highest doubles of the real parts of b;
 *   brehilo  second lowest doubles of the real parts of b;
 *   brelolo  lowest doubles of the real parts of b;
 *   bimhihi  highest doubles of the imaginary parts of b;
 *   bimlohi  second highest doubles of the imaginary parts of b;
 *   bimhilo  second lowest doubles of the imaginary parts of b;
 *   bimlolo  lowest doubles of the imaginary parts of b;
 *   xrehihi  has space allocated for dim doubles;
 *   xrelohi  has space allocated for dim doubles;
 *   xrehilo  has space allocated for dim doubles;
 *   xrelolo  has space allocated for dim doubles;
 *   ximhihi  has space allocated for dim doubles;
 *   ximlohi  has space allocated for dim doubles;
 *   ximhilo  has space allocated for dim doubles;
 *   ximlolo  has space allocated for dim doubles.
 *
 * ON RETURN :
 *   xrehihi  highest doubles of the real parts of the solution;
 *   xrelohi  second highest doubles of the real parts of the solution;
 *   xrehilo  second lowest doubles of the real parts of the solution;
 *   xrelolo  lowest doubles of the real parts of the solution;
 *   ximhihi  highest doubles of the imaginary parts of the solution;
 *   ximlohi  second highest doubles of the imaginary parts of the solution;
 *   ximhilo  second lowest doubles of the imaginary parts of the solution;
 *   ximlolo  lowest doubles of the imaginary parts of the solution;
 *   invlapms is the elapsed time spent by the kernel to invert a tile;
 *   mullapms is the elapsed time spent by the kernel to multiply
 *            with the inversed diagonal tile;
 *   sublapms is the elapsed time spent by the kernel for back substitution;
 *   totlapms is the total elapsed time spent by all kernels;
 *   walltimesec is the elapsed wall clock computation time;
 *   addcnt   counts the number of additions and subtractions;
 *   addover  overflowed number of additions and subtractions;
 *   mulcnt   counts the number of multiplications;
 *   mulover  overflowed the number of multiplications;
 *   divcnt   counts the number of divisions;
 *   divover  overflowed the number of divisions. */

#endif
