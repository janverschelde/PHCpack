/* The file dbl8_tabs_kernels.h specifies functions for the
 * tiled accelerated back substitution in octo double precision. */

#ifndef __dbl8_tabs_kernels_h__
#define __dbl8_tabs_kernels_h__

#define tabsod_shmemsize 169

__global__ void dbl8_small_invert_upper 
( int dim,
  double *Uhihihi, double *Ulohihi, double *Uhilohi, double *Ulolohi,
  double *Uhihilo, double *Ulohilo, double *Uhilolo, double *Ulololo,
  double *invUhihihi, double *invUlohihi,
  double *invUhilohi, double *invUlolohi,
  double *invUhihilo, double *invUlohilo,
  double *invUhilolo, double *invUlololo );
/*
 * DESCRIPTION :
 *   Computes the inverse of an upper triangular matrix.
 *   The U matrix is stored in columnwise fashion,
 *   as the row-by-row computation of the inverse invU
 *   applies a column-by-column load of U.
 *
 * REQUIRED : dim <= 13
 *   Because the inverse is stored entirely in shared memory,
 *   the dimension dim is limited to 13, as 13^2 = 169,
 *   the upper limit on the shared memory, tabsod_shmemsize.
 *
 * ON ENTRY :
 *   dim      dimension of the upper triangular matrix U;
 *   Uhihihi  are the highest doubles of U stored column wise;
 *   Ulohihi  are the second highest doubles of U;
 *   Uhilohi  are the third highest doubles of U;
 *   Ulolohi  are the fourth highest doubles of U;
 *   Uhihilo  are the fourth lowest doubles of U;
 *   Ulohilo  are the third lowest doubles of U;
 *   Uhilolo  are the second lowest doubles of U;
 *   Ulololo  are the lowest doubles of U;
 *   invUhihihi has space allocated for a matrix of dimension dim;
 *   invUlohihi has space allocated for a matrix of dimension dim;
 *   invUhilohi has space allocated for a matrix of dimension dim;
 *   invUlolohi has space allocated for a matrix of dimension dim;
 *   invUhihilo has space allocated for a matrix of dimension dim;
 *   invUlohilo has space allocated for a matrix of dimension dim;
 *   invUhilolo has space allocated for a matrix of dimension dim;
 *   invUlololo has space allocated for a matrix of dimension dim.
 *
 * ON RETURN :
 *   invUhihihi has the highest doubles of the inverse of U, stored row wise;
 *   invUlohihi has the second highest doubles of the inverse of U;
 *   invUhilohi has the third highest doubles of the inverse of U;
 *   invUlolohi has the fourth highest doubles of the inverse of U;
 *   invUhihilo has the fourth lowest doubles the inverse of U;
 *   invUlohilo has the third lowest doubles the inverse of U;
 *   invUhilolo has the second lowest doubles the inverse of U;
 *   invUlololo are the lowest doubles the inverse of the matrix U. */

__global__ void cmplx8_small_invert_upper
 ( int dim,
   double *Urehihihi, double *Urelohihi, double *Urehilohi, double *Urelolohi,
   double *Urehihilo, double *Urelohilo, double *Urehilolo, double *Urelololo,
   double *Uimhihihi, double *Uimlohihi, double *Uimhilohi, double *Uimlolohi,
   double *Uimhihilo, double *Uimlohilo, double *Uimhilolo, double *Uimlololo,
   double *invUrehihihi, double *invUrelohihi,
   double *invUrehilohi, double *invUrelolohi,
   double *invUrehihilo, double *invUrelohilo,
   double *invUrehilolo, double *invUrelololo,
   double *invUimhihihi, double *invUimlohihi,
   double *invUimhilohi, double *invUimlolohi,
   double *invUimhihilo, double *invUimlohilo,
   double *invUimhilolo, double *invUimlololo );
/*
 * DESCRIPTION :
 *   Computes the inverse of an upper triangular matrix.
 *   The U matrix is stored in columnwise fashion,
 *   as the row-by-row computation of the inverse invU
 *   applies a column-by-column load of U.
 *
 * REQUIRED : dim <= 13.
 *   Because the inverse is stored entirely in shared memory,
 *   the dimension dim is limited to 13, as 13^2 = 169,
 *   the upper limit on the shared memory, tabsod_shmemsize.
 *
 * ON ENTRY :
 *   dim       dimension of the upper triangular matrix U,
 *             U is stored column wise;
 *   Urehihihi has the highest doubles of the real parts of U;
 *   Urelohihi has the second highest doubles of the real parts of U;
 *   Urehilohi has the third highest doubles of the real parts of U;
 *   Urelolohi has the fourth highest doubles of the real parts of U;
 *   Urehihilo has the fourth lowest doubles of the real parts of U;
 *   Urelohilo has the third lowest doubles of the real parts of U;
 *   Urehilolo has the second lowest doubles of the real parts of U;
 *   Urelololo has the lowest doubles of the real parts of U;
 *   Uimhihihi has the highest doubles of the imaginary parts of U;
 *   Uimlohihi has the second highest doubles of the imaginary parts of U;
 *   Uimhilohi has the third highest doubles of the imaginary parts of U;
 *   Uimlolohi has the fourth highest doubles of the imaginary parts of U;
 *   Uimhihilo has the fourth lowest doubles of the imaginary parts of U;
 *   Uimlohilo has the third lowest doubles of the imaginary parts of U;
 *   Uimhilolo has the second lowest doubles of the imaginary parts of U;
 *   Uimlololo has the lowest doubles of the imaginary parts of U;
 *   invUrehihihi has space allocated for a matrix of dimension dim;
 *   invUrelohihi has space allocated for a matrix of dimension dim;
 *   invUrehilohi has space allocated for a matrix of dimension dim;
 *   invUrelolohi has space allocated for a matrix of dimension dim;
 *   invUrehihilo has space allocated for a matrix of dimension dim;
 *   invUrelohilo has space allocated for a matrix of dimension dim;
 *   invUrehilolo has space allocated for a matrix of dimension dim;
 *   invUrelololo has space allocated for a matrix of dimension dim;
 *   invUimhihihi has space allocated for a matrix of dimension dim;
 *   invUimlohihi has space allocated for a matrix of dimension dim;
 *   invUimhilohi has space allocated for a matrix of dimension dim;
 *   invUimlolohi has space allocated for a matrix of dimension dim.
 *   invUimhihilo has space allocated for a matrix of dimension dim;
 *   invUimlohilo has space allocated for a matrix of dimension dim;
 *   invUimhilolo has space allocated for a matrix of dimension dim;
 *   invUimlololo has space allocated for a matrix of dimension dim.
 *
 * ON RETURN :
 *   invUrehihihi has the highest doubles of the real parts of the inverse,
 *            stored row wise;
 *   invUrelohihi has the second highest doubles of the real parts
 *            of the inverse, stored row wise;
 *   invUrehilohi has the third highest doubles of the real parts
 *            of the inverse, stored row wise;
 *   invUrelolohi has the fourth highest doubles of the real parts
 *            of the inverse, stored row wise;
 *   invUrehihilo has the fourth lowest doubles of the real parts
 *            of the inverse, stored row wise;
 *   invUrelohilo has the third lowest doubles of the real parts
 *            of the inverse, stored row wise;
 *   invUrehilolo has the second lowest doubles of the real parts
 *            of the inverse, stored row wise;
 *   invUrelololo has the lowest doubles of the real parts of the inverse,
 *            stored row wise;
 *   invUimhihihi has the highest doubles of the imaginary parts
 *            of the inverse, also stored row wise;
 *   invUimlohihi has the second highest doubles of the imaginary parts
 *            of the inverse, also stored row wise;
 *   invUimhilohi has the third highest doubles of the imaginary parts
 *            of the inverse, also stored row wise;
 *   invUimlolohi has the fourth highest doubles of the imaginary parts
 *            of the inverse, also stored row wise;
 *   invUimhihilo has the fourth lowest doubles of the imaginary parts
 *            of the inverse, also stored row wise.
 *   invUimlohilo has the third lowest doubles of the imaginary parts
 *            of the inverse, also stored row wise.
 *   invUimhilolo has the second lowest doubles of the imaginary parts
 *            of the inverse, also stored row wise.
 *   invUimlololo has the lowest doubles of the imaginary parts
 *            of the inverse, also stored row wise. */

__global__ void dbl8_medium_invert_upper
 ( int dim,
   double *Uhihihi, double *Ulohihi, double *Uhilohi, double *Ulolohi,
   double *Uhihilo, double *Ulohilo, double *Uhilolo, double *Ulololo,
   double *invUhihihi, double *invUlohihi,
   double *invUhilohi, double *invUlolohi,
   double *invUhihilo, double *invUlohilo,
   double *invUhilolo, double *invUlololo );
/*
 * DESCRIPTION :
 *   Computes the inverse of an upper triangular matrix.
 *   The U matrix is stored in columnwise fashion,
 *   as the row-by-row computation of the inverse invU
 *   applies a column-by-column load of U.
 *
 * REQUIRED : dim <= 169.
 *   Because the columns of U are loaded entirely into shared memory
 *   and the rows of the inverses are computed first entirely in
 *   shared memory before storing, the dimension dim is limited 
 *   to 169, the upper limit on the shared memory, tabsod_shmemsize.
 *
 * ON ENTRY :
 *   dim      dimension of the upper triangular matrix U;
 *   Uhihihi  are the highest doubles of U stored column wise;
 *   Ulohihi  are the second highest doubles of U;
 *   Uhilohi  are the third highest doubles of U;
 *   Ulolohi  are the fourth highest doubles of U;
 *   Uhihilo  are the fourth lowest doubles of U;
 *   Ulohilo  are the third lowest doubles of U;
 *   Uhilolo  are the second lowest doubles of U;
 *   Ulololo  are the lowest doubles of U;
 *   invUhihihi has space allocated for a matrix of dimension dim;
 *   invUlohihi has space allocated for a matrix of dimension dim;
 *   invUhilohi has space allocated for a matrix of dimension dim;
 *   invUlolohi has space allocated for a matrix of dimension dim;
 *   invUhihilo has space allocated for a matrix of dimension dim;
 *   invUlohilo has space allocated for a matrix of dimension dim;
 *   invUhilolo has space allocated for a matrix of dimension dim;
 *   invUlololo has space allocated for a matrix of dimension dim.
 *
 * ON RETURN :
 *   invUhihihi has the highest doubles of the inverse of U,
 *            stored row wise;
 *   invUlohihi has the second highest doubles of the inverse of U,
 *            stored row wise;
 *   invUhilohi has the third highest doubles of the inverse of U,
 *            stored row wise;
 *   invUlolohi has the fourth highest doubles of the inverse of U,
 *            stored row wise;
 *   invUhihilo has the fourth lowest doubles the inverse of U,
 *            stored row wise;
 *   invUlohilo has the third lowest doubles the inverse of U,
 *            stored row wise;
 *   invUhilolo has the second lowest doubles the inverse of U,
 *            stored row wise;
 *   invUlololo are the lowest doubles the inverse of U,
 *            stored row wise. */

__global__ void cmplx8_medium_invert_upper
 ( int dim, 
   double *Urehihihi, double *Urelohihi, double *Urehilohi, double *Urelolohi,
   double *Urehihilo, double *Urelohilo, double *Urehilolo, double *Urelololo,
   double *Uimhihihi, double *Uimlohihi, double *Uimhilohi, double *Uimlolohi,
   double *Uimhihilo, double *Uimlohilo, double *Uimhilolo, double *Uimlololo,
   double *invUrehihihi, double *invUrelohihi,
   double *invUrehilohi, double *invUrelolohi,
   double *invUrehihilo, double *invUrelohilo,
   double *invUrehilolo, double *invUrelololo,
   double *invUimhihihi, double *invUimlohihi,
   double *invUimhilohi, double *invUimlolohi,
   double *invUimhihilo, double *invUimlohilo,
   double *invUimhilolo, double *invUimlololo );
/*
 * DESCRIPTION :
 *   Computes the inverse of an upper triangular matrix.
 *   The U matrix is stored in columnwise fashion,
 *   as the row-by-row computation of the inverse invU
 *   applies a column-by-column load of U.
 *
 * REQUIRED : dim <= 169.
 *   Because the columns of U are loaded entirely into shared memory
 *   and the rows of the inverses are computed first entirely in
 *   shared memory before storing, the dimension dim is limited 
 *   to 169, the upper limit on the shared memory, tabsod_shmemsize.
 *
 * ON ENTRY :
 *   dim       dimension of the upper triangular matrix U;
 *   Urehihihi are the highest doubles of the real parts of U;
 *   Urelohihi are the second highest doubles of the real parts of U;
 *   Urehilohi are the third highest doubles of the real parts of U;
 *   Urelolohi are the fourth highest doubles of the real parts of U;
 *   Urehihilo are the fourth lowest doubles of the real parts of U;
 *   Urelohilo are the third lowest doubles of the real parts of U;
 *   Urehilolo are the second lowest doubles of the real parts of U;
 *   Urelololo are the lowest doubles of the real parts of U;
 *   Uimhihihi are the highest doubles of the imaginary parts of U;
 *   Uimlohihi are the second highest doubles of the imaginary parts of U;
 *   Uimhilohi are the third highest doubles of the imaginary parts of U;
 *   Uimlolohi are the fourth highest doubles of the imaginary parts of U;
 *   Uimhihilo are the fourth lowest doubles of the imaginary parts of U;
 *   Uimlohilo are the third lowest doubles of the imaginary parts of U;
 *   Uimhilolo are the second lowest doubles of the imaginary parts of U;
 *   Uimlololo are the lowest doubles of the imaginary parts of U;
 *   invUrehihihi has space allocated for a matrix of dimension dim;
 *   invUrelohihi has space allocated for a matrix of dimension dim;
 *   invUrehilohi has space allocated for a matrix of dimension dim;
 *   invUrelolohi has space allocated for a matrix of dimension dim;
 *   invUrehihilo has space allocated for a matrix of dimension dim;
 *   invUrelohilo has space allocated for a matrix of dimension dim;
 *   invUrehilolo has space allocated for a matrix of dimension dim;
 *   invUrelololo has space allocated for a matrix of dimension dim;
 *   invUimhihihi has space allocated for a matrix of dimension dim;
 *   invUimlohihi has space allocated for a matrix of dimension dim;
 *   invUimhilohi has space allocated for a matrix of dimension dim;
 *   invUimlolohi has space allocated for a matrix of dimension dim;
 *   invUimhihilo has space allocated for a matrix of dimension dim;
 *   invUimlohilo has space allocated for a matrix of dimension dim;
 *   invUimhilolo has space allocated for a matrix of dimension dim;
 *   invUimlololo has space allocated for a matrix of dimension dim.
 *
 * ON RETURN :
 *   invUrehihihi has the highest doubles of the real parts of the inverse,
 *            stored row wise;
 *   invUrelohihi has the second highest doubles of the real parts
 *            of the inverse, stored row wise;
 *   invUrehilohi has the third highest doubles of the real parts
 *            of the inverse, stored row wise;
 *   invUrelolohi has the fourth highest doubles of the real parts
 *            of the inverse, stored row wise;
 *   invUrehihilo has the fourth lowest doubles of the real parts
 *            of the inverse, stored row wise;
 *   invUrelohilo has the third lowest doubles of the real parts
 *            of the inverse, stored row wise;
 *   invUrehilolo has the second lowest doubles of the real parts
 *            of the inverse, stored row wise;
 *   invUrelololo has the lowest doubles of the real parts of the inverse,
 *            stored row wise;
 *   invUimhihihi has the highest doubles of the imaginary parts of the
 *            inverse, also stored row wise;
 *   invUimlohihi has the second highest doubles of the imaginary parts
 *            of the inverse, also stored row wise;
 *   invUimhilohi has the third highest doubles of the imaginary parts
 *            of the inverse, also stored row wise;
 *   invUimlolohi has the fourth highest doubles of the imaginary parts
 *            of the inverse, also stored row wise;
 *   invUimhihilo has the fourth lowest doubles of the imaginary parts
 *            of the inverse, also stored row wise.
 *   invUimlohilo has the third lowest doubles of the imaginary parts
 *            of the inverse, also stored row wise.
 *   invUimhilolo has the second lowest doubles of the imaginary parts
 *            of the inverse, also stored row wise.
 *   invUimlololo has the lowest doubles of the imaginary parts of the
 *            inverse, also stored row wise. */

__global__ void  dbl8_invert_tiles
 ( int dim,
   double *Uhihihi, double *Ulohihi, double *Uhilohi, double *Ulolohi,
   double *Uhihilo, double *Ulohilo, double *Uhilolo, double *Ulololo,
   double *invUhihihi, double *invUlohihi,
   double *invUhilohi, double *invUlolohi,
   double *invUhihilo, double *invUlohilo,
   double *invUhilolo, double *invUlololo );
/*
 * DESCRIPTION :
 *   Replaces the columns of the tiles with the rows of the inverses.
 *   The number of blocks equals the number of tiles in U.
 *   The number of threads per block equals the dimension of each tile.
 *
 * REQUIRED : dim <= 169 = tabsod_shmemsize.
 *
 * ON ENTRY :
 *   dim      the dimension of each tile;
 *   Uhihihi  are the highest doubles of columns of all tiles on the
 *            diagonal of an upper triangular matrix;
 *   Ulohihi  are the second highest doubles of columns of all tiles on the
 *            diagonal of an upper triangular matrix;
 *   Uhilohi  are the third highest doubles of columns of all tiles on the
 *            diagonal of an upper triangular matrix;
 *   Ulolohi  are the fourth highest doubles of columns of all tiles on the
 *            diagonal of an upper triangular matrix;
 *   Uhihilo  are the fourth lowest doubles of columns of all tiles on the
 *            diagonal of an upper triangular matrix;
 *   Ulohilo  are the third lowest doubles of columns of all tiles on the
 *            diagonal of an upper triangular matrix;
 *   Uhilolo  are the second lowest doubles of columns of all tiles on the
 *            diagonal of an upper triangular matrix;
 *   Ulololo  are the lowest doubles of columns of all tiles on the
 *            diagonal of an upper triangular matrix;
 *   invUhihihi has space allocated for the inverse of all tiles in U;
 *   invUlohihi has space allocated for the inverse of all tiles in U;
 *   invUhilohi has space allocated for the inverse of all tiles in U;
 *   invUlolohi has space allocated for the inverse of all tiles in U;
 *   invUhihilo has space allocated for the inverse of all tiles in U;
 *   invUlohilo has space allocated for the inverse of all tiles in U;
 *   invUhilolo has space allocated for the inverse of all tiles in U;
 *   invUlololo has space allocated for the inverse of all tiles in U.
 *
 * ON RETURN :
 *   invUhihihi are the highest doubles of the inverse of the tiles in U;
 *   invUlohihi are the second highest doubles
 *              of the inverse of the tiles in U;
 *   invUhilohi are the third highest doubles
 *              of the inverse of the tiles in U;
 *   invUlolohi are the fourth highest doubles
 *              of the inverse of the tiles in U;
 *   invUhihilo are the fourth lowest doubles
 *              of the inverse of the tiles in U;
 *   invUlohilo are the third lowest doubles
 *              of the inverse of the tiles in U;
 *   invUhilolo are the second lowest doubles
 *              of the inverse of the tiles in U;
 *   invUlololo are the lowest doubles of the inverse of the tiles in U. */

__global__ void cmplx8_invert_tiles
 ( int dim, 
   double *Urehihihi, double *Urelohihi, double *Urehilohi, double *Urelolohi,
   double *Urehihilo, double *Urelohilo, double *Urehilolo, double *Urelololo,
   double *Uimhihihi, double *Uimlohihi, double *Uimhilohi, double *Uimlolohi,
   double *Uimhihilo, double *Uimlohilo, double *Uimhilolo, double *Uimlololo,
   double *invUrehihihi, double *invUrelohihi,
   double *invUrehilohi, double *invUrelolohi,
   double *invUrehihilo, double *invUrelohilo,
   double *invUrehilolo, double *invUrelololo,
   double *invUimhihihi, double *invUimlohihi,
   double *invUimhilohi, double *invUimlolohi,
   double *invUimhihilo, double *invUimlohilo,
   double *invUimhilolo, double *invUimlololo );
/*
 * DESCRIPTION :
 *   Replaces the columns of the tiles with the rows of the inverses.
 *   The number of blocks equals the number of tiles in U.
 *   The number of threads per block equals the dimension of each tile.
 *
 * REQUIRED : dim <= 169 = tabsod_shmemsize.
 *
 * ON ENTRY :
 *   dim       the dimension of each tile;
 *   Urehihihi are the highest doubles of the real parts of the columns
 *             of all tiles on the diagonal of an upper triangular matrix;
 *   Urelohihi are the second highest doubles of the real parts of the columns
 *             of all tiles on the diagonal of an upper triangular matrix;
 *   Urehilohi are the third highest doubles of the real parts of the columns
 *             of all tiles on the diagonal of an upper triangular matrix;
 *   Urelolohi are the fourth highest doubles of the real parts of the columns
 *             of all tiles on the diagonal of an upper triangular matrix;
 *   Urehihilo are the fourth lowest doubles of the real parts of the columns
 *             of all tiles on the diagonal of an upper triangular matrix;
 *   Urelohilo are the third lowest doubles of the real parts of the columns
 *             of all tiles on the diagonal of an upper triangular matrix;
 *   Urehilolo are the second lowest doubles of the real parts of the columns
 *             of all tiles on the diagonal of an upper triangular matrix;
 *   Urelololo are the lowest doubles of the real parts of the columns
 *             of all tiles on the diagonal of an upper triangular matrix;
 *   Uimhihihi are the highest doubles of the imaginary parts of the columns
 *             of all tiles on the diagonal of an upper triangular matrix;
 *   Uimlohihi are the second highest doubles of the imaginary parts
 *             of all tiles on the diagonal of an upper triangular matrix;
 *   Uimhilohi are the third highest doubles of the imaginary parts
 *             of all tiles on the diagonal of an upper triangular matrix;
 *   Uimlohihi are the fourth highest doubles of the imaginary parts
 *             of all tiles on the diagonal of an upper triangular matrix;
 *   Uimhihilo are the fourth lowest doubles of the imaginary parts
 *             of all tiles on the diagonal of an upper triangular matrix;
 *   Uimlohilo are the third lowest doubles of the imaginary parts
 *             of all tiles on the diagonal of an upper triangular matrix;
 *   Uimhilolo are the second lowest doubles of the imaginary parts
 *             of all tiles on the diagonal of an upper triangular matrix;
 *   Uimlololo are the lowest doubles of the imaginary parts of the columns
 *             of all tiles on the diagonal of an upper triangular matrix;
 *   invUrehihihi has space allocated for the highest doubles of
 *             the real parts of the inverse of all tiles in U;
 *   invUrelohihi has space allocated for the second highest doubles of
 *             the real parts of the inverse of all tiles in U;
 *   invUrehilohi has space allocated for the third highest doubles of
 *             the real parts of the inverse of all tiles in U;
 *   invUrelolohi has space allocated for the fourth highest doubles of
 *             the real parts of the inverse of all tiles in U;
 *   invUrehihilo has space allocated for the fourth lowest doubles of
 *             the real parts of the inverse of all tiles in U;
 *   invUrelohilo has space allocated for the third lowest doubles of
 *             the real parts of the inverse of all tiles in U;
 *   invUrehilolo has space allocated for the second lowest doubles of
 *             the real parts of the inverse of all tiles in U;
 *   invUrelololo has space allocated for the lowest doubles of
 *             the real parts of the inverse of all tiles in U;
 *   invUimhihihi has space allocated for the highest doubles of
 *             the imaginary parts of the inverse of all tiles in U;
 *   invUimlohihi has space allocated for the second highest doubles of
 *             the imaginary parts of the inverse of all tiles in U;
 *   invUimhilohi has space allocated for the third highest doubles of
 *            the imaginary parts of the inverse of all tiles in U.
 *   invUimlolohi has space allocated for the fourth highest doubles of
 *            the imaginary parts of the inverse of all tiles in U.
 *   invUimhihilo has space allocated for the fourth lowest doubles of
 *             the imaginary parts of the inverse of all tiles in U;
 *   invUimlohilo has space allocated for the third lowest doubles of
 *             the imaginary parts of the inverse of all tiles in U;
 *   invUimhilolo has space allocated for the second lowest doubles of
 *             the imaginary parts of the inverse of all tiles in U.
 *   invUimlololo has space allocated for the lowest doubles of
 *             the imaginary parts of the inverse of all tiles in U.
 *
 * ON RETURN :
 *   invUrehihihi has the highest doubles of the real parts
 *             of rows of the inverse tiles;
 *   invUrelohihi has the second highest doubles of the real parts
 *             of rows of the inverse tiles;
 *   invUrehilohi has the third highest doubles of the real parts
 *             of rows of the inverse tiles;
 *   invUrelolohi has the fourth highest doubles of the real parts
 *             of rows of the inverse tiles;
 *   invUrehihilo has the fourth lowest doubles of the real parts
 *             of rows of the inverse tiles;
 *   invUrelohilo has the third lowest doubles of the real parts
 *             of rows of the inverse tiles;
 *   invUrehilolo has the second lowest doubles of the real parts
 *             of rows of the inverse tiles;
 *   invUrelololo has the lowest doubles of the real parts
 *             of rows of the inverse tiles;
 *   invUimhihihi has the highest doubles of the imaginary parts
 *             of rows of the inverse tiles;
 *   invUimlohihi has the second highest doubles of the imaginary parts
 *             of rows of the inverse tiles;
 *   invUimhilohi has the third highest doubles of the imaginary parts
 *             of rows of the inverse tiles;
 *   invUimlolohi has the fourth highest doubles of the imaginary parts
 *             of rows of the inverse tiles;
 *   invUimhihilo has the fourth lowest doubles of the imaginary parts
 *             of rows of the inverse tiles;
 *   invUimlohilo has the third lowest doubles of the imaginary parts
 *             of rows of the inverse tiles;
 *   invUimhilolo has the second lowest doubles of the imaginary parts
 *             of rows of the inverse tiles;
 *   invUimlololo has the lowest doubles of the imaginary parts
 *             of rows of the inverse tiles. */

__global__ void dbl8_multiply_inverse
 ( int dim, int idx,
   double *invUhihihi, double *invUlohihi,
   double *invUhilohi, double *invUlolohi,
   double *invUhihilo, double *invUlohilo,
   double *invUhilolo, double *invUlololo,
   double *whihihi, double *wlohihi, double *whilohi, double *wlolohi,
   double *whihilo, double *wlohilo, double *whilolo, double *wlololo );
/*
 * DESCRIPTION :
 *   Replaces b with the product of the inverse tile in U.
 *
 * ON ENTRY :
 *   dim      the dimension of each tile;
 *   idx      index of the diagonal tile;
 *   invUhihihi are the highest doubles of the inverse of the diagonal tiles;
 *   invUlohihi are the second highest doubles of the inverse diagonal tiles;
 *   invUhilohi are the third highest doubles of the inverse diagonal tiles;
 *   invUlolohi are the fourth highest doubles of the inverse diagonal tiles;
 *   invUhihilo are the fourth lowest doubles of the inverse diagonal tiles;
 *   invUlohilo are the third lowest doubles of the inverse diagonal tiles;
 *   invUhilolo are the second lowest doubles of the inverse diagonal tiles;
 *   invUlololo are the  lowest doubles of the inverse of the diagonal tiles;
 *   whihihi  are the highest doubles of the right hand side vector;
 *   wlohihi  are the second highest doubles of the right hand side vector;
 *   whilohi  are the third highest doubles of the right hand side vector;
 *   wlolohi  are the fourth highest doubles of the right hand side vector;
 *   whihilo  are the fourth lowest doubles of the right hand side vector.
 *   wlohilo  are the third lowest doubles of the right hand side vector.
 *   whilolo  are the second lowest doubles of the right hand side vector.
 *   wlololo  are the lowest doubles of the right hand side vector.
 *
 * ON RETURN :
 *   whihihi  are the highest doubles of the product of the inverse of the
 *            diagonal tile defined by the index idx with the w on input;
 *   wlohihi  are the second highest doubles of the product of the inverse of
 *            the diagonal tile defined by the index idx with the w on input;
 *   whilohi  are the third highest doubles of the product of the inverse of
 *            the diagonal tile defined by the index idx with the w on input;
 *   wlolohi  are the fourth lowest doubles of the product of the inverse of
 *            the diagonal tile defined by the index idx with the w on input;
 *   whihilo  are the fourth lowest doubles of the product of the inverse of
 *            the diagonal tile defined by the index idx with the w on input;
 *   wlohilo  are the third lowest doubles of the product of the inverse of
 *            the diagonal tile defined by the index idx with the w on input;
 *   whilolo  are the second lowest doubles of the product of the inverse of
 *            the diagonal tile defined by the index idx with the w on input;
 *   wlololo  are the lowest doubles of the product of the inverse of the
 *            diagonal tile defined by the index idx with the w on input. */

__global__ void cmplx8_multiply_inverse
 ( int dim, int idx,
   double *invUrehihihi, double *invUrelohihi,
   double *invUrehilohi, double *invUrelolohi,
   double *invUrehihilo, double *invUrelohilo,
   double *invUrehilolo, double *invUrelololo,
   double *invUimhihihi, double *invUimlohihi,
   double *invUimhilohi, double *invUimlolohi,
   double *invUimhihilo, double *invUimlohilo,
   double *invUimhilolo, double *invUimlololo,
   double *wrehihihi, double *wrelohihi,
   double *wrehilohi, double *wrelolohi,
   double *wrehihilo, double *wrelohilo,
   double *wrehilolo, double *wrelololo,
   double *wimhihihi, double *wimlohihi,
   double *wimhilohi, double *wimlolohi,
   double *wimhihilo, double *wimlohilo,
   double *wimhilolo, double *wimlololo );
/*
 * DESCRIPTION :
 *   Replaces b with the product of the inverse tile in U.
 *
 * ON ENTRY :
 *   dim      the dimension of each tile;
 *   idx      index of the diagonal tile;
 *   invUrehihihi are the highest doubles of the real parts of
 *            the inverse of the diagonal tiles;
 *   invUrelohihi are the second highest doubles of the real parts of
 *            the inverse of the diagonal tiles;
 *   invUrehilohi are the third highest doubles of the real parts of
 *            the inverse of the diagonal tiles;
 *   invUrelolohi are the fourth highest doubles of the real parts of
 *            the inverse of the diagonal tiles;
 *   invUrehihilo are the fourth lowest doubles of the real parts of
 *            the inverse of the diagonal tiles;
 *   invUrelohilo are the third lowest doubles of the real parts of
 *            the inverse of the diagonal tiles;
 *   invUrehilolo are the second lowest doubles of the real parts of
 *            the inverse of the diagonal tiles;
 *   invUrelololo are the lowest doubles of the real parts of
 *            the inverse of the diagonal tiles;
 *   invUimhihihi are the highest doubles of the imaginary parts of
 *            the inverse of the diagonal tiles;
 *   invUimlohihi are the second highest doubles of the imaginary parts of
 *            the inverse of the diagonal tiles;
 *   invUimhilohi are the third highest doubles of the imaginary parts of
 *            the inverse of the diagonal tiles;
 *   invUimlolohi are the fourth highest doubles of the imaginary parts of
 *            the inverse of the diagonal tiles;
 *   invUimhihilo are the fourth lowest doubles of the imaginary parts of
 *            the inverse of the diagonal tiles;
 *   invUimlohilo are the third lowest doubles of the imaginary parts of
 *            the inverse of the diagonal tiles;
 *   invUimhilolo are the second lowest doubles of the imaginary parts of
 *            the inverse of the diagonal tiles;
 *   invUimlololo are the lowest doubles of the imaginary parts of
 *            the inverse of the diagonal tiles;
 *   wrehihihi are the highest doubles of the real parts of b;
 *   wrelohihi are the second highest doubles of the real parts of b;
 *   wrehilohi are the third highest doubles of the real parts of b;
 *   wrelolohi are the fourth highest doubles of the real parts of b;
 *   wrehihilo are the fourth lowest doubles of the real parts of b;
 *   wrelohilo are the third lowest doubles of the real parts of b;
 *   wrehilolo are the second lowest doubles of the real parts of b;
 *   wrelololo are the lowest doubles of the real parts of b;
 *   wimhihihi are the highest doubles of the imaginary parts of b;
 *   wimlohihi are the second highest doubles of the imaginary parts of b;
 *   wimhilohi are the third highest doubles of the imaginary parts of b;
 *   wimlolohi are the fourth highest doubles of the imaginary parts of b;
 *   wimhihilo are the fourth lowest doubles of the imaginary parts of b.
 *   wimlohilo are the third lowest doubles of the imaginary parts of b.
 *   wimhilolo are the second lowest doubles of the imaginary parts of b.
 *   wimlololo are the lowest doubles of the imaginary parts of b.
 *
 * ON RETURN :
 *   wrehihihi are the highest doubles of the real parts of the product
 *            of the inverse of the tile defined by the index idx;
 *   wrelohihi are the second highest doubles of the real parts of the
 *            product of the inverse of the tile defined by the index idx;
 *   wrehilohi are the third highest doubles of the real parts of the
 *            product of the inverse of the tile defined by the index idx;
 *   wrelolohi are the fourth highest doubles of the real parts of the
 *            product of the inverse of the tile defined by the index idx;
 *   wrehihilo are the fourth lowest doubles of the real parts of the
 *            product of the inverse of the tile defined by the index idx;
 *   wrelohilo are the third lowest doubles of the real parts of the
 *            product of the inverse of the tile defined by the index idx;
 *   wrehilolo are the second lowest doubles of the real parts of the
 *            product of the inverse of the tile defined by the index idx;
 *   wrelololo are the lowest doubles of the real parts of the product
 *            of the inverse of the tile defined by the index idx;
 *   wimhihihi are the highest doubles of the imaginary parts of the product
 *            of the inverse of the tile defined by the index idx;
 *   wimlohihi are the second highest doubles of the imaginary parts of the
 *            product of the inverse of the tile defined by the index idx;
 *   wimhilohi are the third highest doubles of the imaginary parts of the
 *            product of the inverse of the tile defined by the index idx;
 *   wimlolohi are the fourth highest doubles of the imaginary parts of the
 *            product of the inverse of the tile defined by the index idx;
 *   wimhihilo are the fourth lowest doubles of the imaginary parts of the
 *            product of the inverse of the tile defined by the index idx;
 *   wimlohilo are the third lowest doubles of the imaginary parts of the
 *            product of the inverse of the tile defined by the index idx;
 *   wimhilolo are the second lowest doubles of the imaginary parts of the
 *            product of the inverse of the tile defined by the index idx;
 *   wimlololo are the lowest doubles of the imaginary parts of the product
 *            of the inverse of the tile defined by the index idx. */

__global__ void dbl8_back_substitute
 ( int dim, int idx, 
   double *Uhihihi, double *Ulohihi, double *Uhilohi, double *Ulolohi, 
   double *Uhihilo, double *Ulohilo, double *Uhilolo, double *Ulololo, 
   double *whihihi, double *wlohihi, double *whilohi, double *wlolohi,
   double *whihilo, double *wlohilo, double *whilolo, double *wlololo );
/*
 * DESCRIPTION :
 *   Updates the right hand side vector subtracting the solution
 *   defined by idx, multiplied with the corresponding rows in U.
 *
 * ON ENTRY :
 *   dim      dimension of each tile;
 *   idx      index of the solution tile in the multiplication;
 *   Uhihihi  are the highest doubles of tiles to multiply with;
 *   Ulohihi  are the second highest doubles of tiles to multiply with;
 *   Uhilohi  are the third highest doubles of tiles to multiply with;
 *   Ulolohi  are the fourth highest doubles of tiles to multiply with;
 *   Uhihilo  are the fourth lowest doubles of tiles to multiply with;
 *   Ulohilo  are the third lowest doubles of tiles to multiply with;
 *   Uhilolo  are the second lowest doubles of tiles to multiply with;
 *   Ulololo  are the lowest doubles of tiles to multiply with;
 *   whihihi  are the highest doubles of w;
 *   wlohihi  are the second highest doubles of w;
 *   whilohi  are the third highest doubles of w;
 *   wlolohi  are the fourth highest doubles of w;
 *   whihilo  are the fourth lowest doubles of w;
 *   wlohilo  are the third lowest doubles of w;
 *   whilolo  are the second lowest doubles of w;
 *   wlololo  are the lowest doubles of w.
 *
 * ON RETURN :
 *   whihihi  are the highest doubles of the updated w;
 *   wlohihi  are the second highest doubles of the updated w;
 *   whilohi  are the third highest doubles of the updated w;
 *   wlolohi  are the fourth highest doubles of the updated w;
 *   whihilo  are the fourth lowest doubles of the updated w;
 *   wlohilo  are the third lowest doubles of the updated w;
 *   whilolo  are the second lowest doubles of the updated w;
 *   wlololo  are the lowest doubles of the updated w. */

__global__ void cmplx8_back_substitute
 ( int dim, int idx,
   double *Urehihihi, double *Urelohihi,
   double *Urehilohi, double *Urelolohi,
   double *Urehihilo, double *Urelohilo,
   double *Urehilolo, double *Urelololo,
   double *Uimhihihi, double *Uimlohihi,
   double *Uimhilohi, double *Uimlolohi,
   double *Uimhihilo, double *Uimlohilo,
   double *Uimhilolo, double *Uimlololo,
   double *wrehihihi, double *wrelohihi,
   double *wrehilohi, double *wrelolohi,
   double *wrehihilo, double *wrelohilo,
   double *wrehilolo, double *wrelololo,
   double *wimhihihi, double *wimlohihi,
   double *wimhilohi, double *wimlolohi,
   double *wimhihilo, double *wimlohilo,
   double *wimhilolo, double *wimlololo );
/*
 * DESCRIPTION :
 *   Updates the right hand side vector subtracting the solution
 *   defined by idx, multiplied with the corresponding rows in U.
 *
 * ON ENTRY :
 *   dim      dimension of each tile;
 *   idx      index of the solution tile in the multiplication;
 *   Urehihihi are the highest doubles of the real parts of tiles 
 *            to multiply the solution with;
 *   Urelohihi are the second highest doubles of the real parts of tiles 
 *            to multiply the solution with;
 *   Urehilohi are the third highest doubles of the real parts of tiles 
 *            to multiply the solution with;
 *   Urelolohi are the fourth highest doubles of the real parts of tiles 
 *            to multiply the solution with;
 *   Urehihilo are the fourth lowest doubles of the real parts of tiles
 *            to multiply the solution with;
 *   Urelohilo are the third lowest doubles of the real parts of tiles
 *            to multiply the solution with;
 *   Urehilolo are the second lowest doubles of the real parts of tiles
 *            to multiply the solution with;
 *   Urelololo are the lowest doubles of the real parts of tiles
 *            to multiply the solution with;
 *   Uimhihihi are the highest doubles of the imaginary parts of tiles
 *            to multiply the solution with;
 *   Uimlohihi are the second highest doubles of the imaginary parts of tiles
 *            to multiply the solution with;
 *   Uimhilohi are the third highest doubles of the imaginary parts of tiles
 *            to multiply the solution with;
 *   Uimlolohi are the fourth highest doubles of the imaginary parts of tiles
 *            to multiply the solution with;
 *   Uimhihilo are the fourth lowest doubles of the imaginary parts of tiles
 *            to multiply the solution with;
 *   Uimlohilo are the third lowest doubles of the imaginary parts of tiles
 *            to multiply the solution with;
 *   Uimhilolo are the second lowest doubles of the imaginary parts of tiles
 *            to multiply the solution with;
 *   Uimlololo are the lowest doubles of the imaginary parts of tiles
 *            to multiply the solution with;
 *   wrehihihi are the highest doubles of the real parts of 
 *            the current right hand side vector;
 *   wrelohihi are the second highest doubles of the real parts of 
 *            the current right hand side vector;
 *   wrehilohi are the third highest doubles of the real parts of 
 *            the current right hand side vector;
 *   wrelolohi are the fourth highest doubles of the real parts of 
 *            the current right hand side vector;
 *   wrehihilo are the fourth lowest doubles of the real parts of 
 *            the current right hand side vector;
 *   wrelohilo are the third lowest doubles of the real parts of 
 *            the current right hand side vector;
 *   wrehilolo are the second lowest doubles of the real parts of 
 *            the current right hand side vector;
 *   wrelololo are the lowest doubles of the real parts of 
 *            the current right hand side vector;
 *   wimhihihi are the highest doubles of the imaginary parts
 *            of the current right hand side vector;
 *   wimlohihi are the second highest doubles of the imaginary parts
 *            of the current right hand side vector;
 *   wimhilohi are the third highest doubles of the imaginary parts
 *            of the current right hand side vector;
 *   wimlolohi are the fourth highest doubles of the imaginary parts
 *            of the current right hand side vector;
 *   wimhihilo are the fourth lowest doubles of the imaginary parts
 *            of the current right hand side vector;
 *   wimlohilo are the third lowest doubles of the imaginary parts
 *            of the current right hand side vector;
 *   wimhilolo are the second lowest doubles of the imaginary parts
 *            of the current right hand side vector;
 *   wimlololo are the lowest doubles of the imaginary parts
 *            of the current right hand side vector.
 *
 * ON RETURN :
 *   wrehihihi are the highest doubles of the real parts
 *            of the updated right hand side vector;
 *   wrelohihi are the second highest doubles of the real parts
 *            of the updated right hand side vector;
 *   wrehilohi are the third highest doubles of the real parts
 *            of the updated right hand side vector;
 *   wrelolohi are the fourth highest doubles of the real parts
 *            of the updated right hand side vector;
 *   wrehihilo are the fourth lowest doubles of the real parts
 *            of the updated right hand side vector;
 *   wrelohilo are the third lowest doubles of the real parts
 *            of the updated right hand side vector;
 *   wrehilolo are the second lowest doubles of the real parts
 *            of the updated right hand side vector;
 *   wrelololo are the lowest doubles of the real parts
 *            of the updated right hand side vector;
 *   wimhihihi are the highest doubles of the imaginary parts
 *            of the updated right hand side vector;
 *   wimlohihi are the second highest doubles of the imaginary parts
 *            of the updated right hand side vector;
 *   wimhilohi are the third highest doubles of the imaginary parts
 *            of the updated right hand side vector;
 *   wimlolohi are the fourth highest doubles of the imaginary parts
 *            of the updated right hand side vector;
 *   wimhihilo are the fourth lowest doubles of the imaginary parts
 *            of the updated right hand side vector;
 *   wimlohilo are the third lowest doubles of the imaginary parts
 *            of the updated right hand side vector;
 *   wimhilolo are the second lowest doubles of the imaginary parts
 *            of the updated right hand side vector;
 *   wimlololo are the lowest doubles of the imaginary parts
 *            of the updated right hand side vector. */

void GPU_dbl8_upper_inverse
 ( int dim,
   double **Uhihihi, double **Ulohihi, double **Uhilohi, double **Ulolohi,
   double **Uhihilo, double **Ulohilo, double **Uhilolo, double **Ulololo,
   double **invUhihihi, double **invUlohihi,
   double **invUhilohi, double **invUlolohi,
   double **invUhihilo, double **invUlohilo,
   double **invUhilolo, double **invUlololo,
   double *lapms, double *walltimesec );
/*
 * DESCRIPTION :
 *   Calls the kernel to invert the upper triangular matrix U.
 *   The matrices are stored in the conventional rowwise fashion.
 *
 * ON ENTRY :
 *   dim      dimension of the upper triangular matrix U;
 *   Uhihihi  has the highest doubles of U;
 *   Ulohihi  has the second highest doubles of U;
 *   Uhilohi  has the third highest doubles of U;
 *   Ulolohi  has the fourth highest doubles of U;
 *   Uhihilo  has the fourth lowest doubles of U;
 *   Ulohilo  has the third lowest doubles of U;
 *   Uhilolo  has the second lowest doubles of U;
 *   Ulololo  has the lowest doubles of U;
 *   invUhihi has space allocated for a matrix of dimension dim;
 *   invUlohi has space allocated for a matrix of dimension dim;
 *   invUhilo has space allocated for a matrix of dimension dim;
 *   invUlolo has space allocated for a matrix of dimension dim.
 *
 * ON RETURN :
 *   invUhihihi has the highest doubles of the inverse of U;
 *   invUlohihi has the second highest doubles of the inverse of U;
 *   invUhilohi has the third highest doubles of the inverse of U;
 *   invUlolohi has the fourth highest doubles of the inverse of U;
 *   invUhihilo has the fourth lowest doubles of the inverse of U;
 *   invUlohilo has the third lowest doubles of the inverse of U;
 *   invUhilolo has the second lowest doubles of the inverse of U;
 *   invUlololo has the lowest doubles of the inverse of U;
 *   lapms    elapsed time spent by the kernels;
 *   walltimesec is the elapsed wall clock computation time. */

void GPU_cmplx8_upper_inverse
 ( int dim,
   double **Urehihihi, double **Urelohihi,
   double **Urehilohi, double **Urelolohi,
   double **Urehihilo, double **Urelohilo,
   double **Urehilolo, double **Urelololo,
   double **Uimhihihi, double **Uimlohihi,
   double **Uimhilohi, double **Uimlolohi,
   double **Uimhihilo, double **Uimlohilo,
   double **Uimhilolo, double **Uimlololo,
   double **invUrehihihi, double **invUrelohihi,
   double **invUrehilohi, double **invUrelolohi,
   double **invUrehihilo, double **invUrelohilo,
   double **invUrehilolo, double **invUrelololo,
   double **invUimhihihi, double **invUimlohihi, 
   double **invUimhilohi, double **invUimlolohi, 
   double **invUimhihilo, double **invUimlohilo, 
   double **invUimhilolo, double **invUimlololo, 
   double *lapms, double *walltimesec );
/*
 * DESCRIPTION :
 *   Calls the kernel to invert the upper triangular matrix U.
 *   The matrices are stored in the conventional row wise fashion.
 *
 * ON ENTRY :
 *   dim       dimension of the upper triangular matrix U;
 *   Urehihihi are the highest doubles of the real parts of U;
 *   Urelohihi are the second highest doubles of the real parts of U;
 *   Urehilohi are the third highest doubles of the real parts of U;
 *   Urelolohi are the fourth highest doubles of the real parts of U;
 *   Urehihilo are the fourth lowest doubles of the real parts of U;
 *   Urelohilo are the third lowest doubles of the real parts of U;
 *   Urehilolo are the second lowest doubles of the real parts of U;
 *   Urelololo are the lowest doubles of the real parts of U;
 *   Uimhihihi are the highest doubles of the imaginary parts of U;
 *   Uimlohihi are the second highest doubles of the imaginary parts of U;
 *   Uimhilohi are the third highest doubles of the imaginary parts of U;
 *   Uimlolohi are the fourth highest doubles of the imaginary parts of U;
 *   Uimhihilo are the fourth lowest doubles of the imaginary parts of U;
 *   Uimlohilo are the third lowest doubles of the imaginary parts of U;
 *   Uimhilolo are the second lowest doubles of the imaginary parts of U;
 *   Uimlololo are the lowest doubles of the imaginary parts of U;
 *   invUrehihihi has space allocated for a matrix of dimension dim;
 *   invUrelohihi has space allocated for a matrix of dimension dim;
 *   invUrehilohi has space allocated for a matrix of dimension dim;
 *   invUrelolohi has space allocated for a matrix of dimension dim;
 *   invUrehihilo has space allocated for a matrix of dimension dim;
 *   invUrelohilo has space allocated for a matrix of dimension dim;
 *   invUrehilolo has space allocated for a matrix of dimension dim;
 *   invUrelololo has space allocated for a matrix of dimension dim;
 *   invUimhihihi has space allocated for a matrix of dimension dim;
 *   invUimlohihi has space allocated for a matrix of dimension dim;
 *   invUimhilohi has space allocated for a matrix of dimension dim;
 *   invUimlolohi has space allocated for a matrix of dimension dim;
 *   invUimhihilo has space allocated for a matrix of dimension dim;
 *   invUimlohilo has space allocated for a matrix of dimension dim;
 *   invUimhilolo has space allocated for a matrix of dimension dim;
 *   invUimlololo has space allocated for a matrix of dimension dim.
 *
 * ON RETURN :
 *   invUrehihihi has the highest doubles of the real parts of the inverse;
 *   invUrelohihi has the second highest doubles of the real parts
 *            of the inverse;
 *   invUrehilohi has the third highest doubles of the real parts
 *            of the inverse;
 *   invUrelolohi has the fourth highest doubles of the real parts
 *            of the inverse;
 *   invUrehihilo has the fourth lowest doubles of the real partls
 *            of the inverse;
 *   invUrelohilo has the third lowest doubles of the real partls
 *            of the inverse;
 *   invUrehilolo has the second lowest doubles of the real partls
 *            of the inverse;
 *   invUrelololo has the lowest doubles of the real parts of the inverse;
 *   invUimhihihi has the highest doubles of the imaginary parts
 *            of the inverse;
 *   invUimlohihi has the second highest doubles of the imaginary parts
 *            of the inverse;
 *   invUimlohihi has the third highest doubles of the imaginary parts
 *            of the inverse;
 *   invUimlolohi has the fourth highest doubles of the imaginary parts
 *            of the inverse;
 *   invUimhihilo has the fourth lowest doubles of the imaginary parts
 *            of the inverse;
 *   invUimlohilo has the third lowest doubles of the imaginary parts
 *            of the inverse;
 *   invUimhilolo has the second lowest doubles of the imaginary parts
 *            of the inverse;
 *   invUimlololo has the lowest doubles of the imaginary parts
 *            of the inverse;
 *   lapms    elapsed time spent by the kernels;
 *   walltimesec is the elapsed wall clock computation time. */

void GPU_dbl8_upper_tiled_solver
 ( int dim, int szt, int nbt,
   double **Uhihihi, double **Ulohihi, double **Uhilohi, double **Ulolohi,
   double **Uhihilo, double **Ulohilo, double **Uhilolo, double **Ulololo,
   double *bhihihi, double *blohihi, double *bhilohi, double *blolohi,
   double *bhihilo, double *blohilo, double *bhilolo, double *blololo,
   double *xhihihi, double *xlohihi, double *xhilohi, double *xlolohi,
   double *xhihilo, double *xlohilo, double *xhilolo, double *xlololo,
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
 *   Uhihihi  are the highest doubles of U;
 *   Ulohihi  are the second highest doubles of U;
 *   Uhilohi  are the third highest doubles of U;
 *   Ulolohi  are the fourth highest doubles of U;
 *   Uhihilo  are the fourth lowest doubles of U;
 *   Ulohilo  are the third lowest doubles of U;
 *   Uhilolo  are the second lowest doubles of U;
 *   Ulololo  are the lowest doubles of U;
 *   bhihihi  are the highest doubles of the right hand side b;
 *   blohihi  are the second highest doubles of the right hand side b;
 *   bhilohi  are the third highest doubles of the right hand side b;
 *   blolohi  are the fourth highest doubles of the right hand side b;
 *   bhihilo  are the fourth lowest doubles of the right hand side b;
 *   blohilo  are the third lowest doubles of the right hand side b;
 *   bhilolo  are the second lowest doubles of the right hand side b;
 *   blololo  are the lowest doubles of the right hand side b;
 *   xhihihi  has space allocated for dim doubles;
 *   xlohihi  has space allocated for dim doubles;
 *   xhilohi  has space allocated for dim doubles;
 *   xlolohi  has space allocated for dim doubles;
 *   xhihilo  has space allocated for dim doubles;
 *   xlohilo  has space allocated for dim doubles;
 *   xhilolo  has space allocated for dim doubles;
 *   xlololo  has space allocated for dim doubles.
 *
 * ON RETURN :
 *   xhihihi  are the highest doubles of the solution to U*x = b;
 *   xlohihi  are the second highest doubles of the solution to U*x = b;
 *   xhilohi  are the third highest doubles of the solution to U*x = b;
 *   xlolohi  are the fourth highest doubles of the solution to U*x = b;
 *   xhihilo  are the fourth lowest doubles of the solution to U*x = b;
 *   xlohilo  are the third lowest doubles of the solution to U*x = b;
 *   xhilolo  are the second lowest doubles of the solution to U*x = b;
 *   xlololo  are the lowest doubles of the solution to U*x = b;
 *   invlapms is the elapsed time spent by the kernel to invert a tile;
 *   mullapms is the elapsed time spent by the kernel to multiply
 *            with the inversed diagonal tile;
 *   sublapms is the elapsed time spent by the kernel for back substitution;
 *   totlapms is the total elapsed time spent by all kernels;
 *   walltimesec is the elapsed wall clock computation time;
 *   addcnt   counts the number of additions and subtractions;
 *   mulcnt   counts the number of multiplications;
 *   divcnt   counts the number of divisions. */

void GPU_cmplx8_upper_tiled_solver
 ( int dim, int szt, int nbt,
   double **Urehihihi, double **Urelohihi,
   double **Urehilohi, double **Urelolohi,
   double **Urehihilo, double **Urelohilo,
   double **Urehilolo, double **Urelololo,
   double **Uimhihihi, double **Uimlohihi,
   double **Uimhilohi, double **Uimlolohi,
   double **Uimhihilo, double **Uimlohilo,
   double **Uimhilolo, double **Uimlololo,
   double *brehihihi, double *brelohihi, double *brehilohi, double *brelolohi,
   double *brehihilo, double *brelohilo, double *brehilolo, double *brelololo,
   double *bimhihihi, double *bimlohihi, double *bimhilohi, double *bimlolohi,
   double *bimhihilo, double *bimlohilo, double *bimhilolo, double *bimlololo,
   double *xrehihihi, double *xrelohihi, double *xrehilohi, double *xrelolohi,
   double *xrehihilo, double *xrelohilo, double *xrehilolo, double *xrelololo,
   double *ximhihihi, double *ximlohihi, double *ximhilohi, double *ximlolohi,
   double *ximhihilo, double *ximlohilo, double *ximhilolo, double *ximlololo,
   double *invlapms, double *mullapms, double *sublapms, double *totlapms,
   double *walltimesec,
   long long int *addcnt, long long int *mulcnt, long long int *divcnt );
/*
 * DESCRIPTION :
 *   Solves an upper triangular system with a tiled algorithm.
 *
 * ON ENTRY :
 *   dim       dimension of the upper triangular matrix U;
 *   szt       size of each tile;
 *   nbt       number of tiles, dim = szt*nbt;
 *   Urehihihi are the highest doubles of the real parts of U;
 *   Urelohihi are the second highest doubles of the real parts of U;
 *   Urehilohi are the third highest doubles of the real parts of U;
 *   Urelolohi are the fourth highest doubles of the real parts of U;
 *   Urehihilo are the fourth lowest doubles of the real parts of U;
 *   Urelohilo are the third lowest doubles of the real parts of U;
 *   Urehilolo are the second lowest doubles of the real parts of U;
 *   Urelololo are the lowest doubles of the real parts of U;
 *   Uimhihihi are the highest doubles of the imaginary parts of U;
 *   Uimlohihi are the second highest doubles of the imaginary parts of U;
 *   Uimhilohi are the third highest doubles of the imaginary parts of U;
 *   Uimlolohi are the fourth highest doubles of the imaginary parts of U;
 *   Uimhihilo are the fourth lowest doubles of the imaginary parts of U;
 *   Uimlohilo are the third lowest doubles of the imaginary parts of U;
 *   Uimhilolo are the second lowest doubles of the imaginary parts of U;
 *   Uimlololo are the lowest doubles of the imaginary parts of U;
 *   brehihihi are the highest doubles of the real parts of b;
 *   brelohihi are the second highest doubles of the real parts of b;
 *   brehilohi are the third highest doubles of the real parts of b;
 *   brelolohi are the fourth highest doubles of the real parts of b;
 *   brehihilo are the fourth lowest doubles of the real parts of b;
 *   brelohilo are the third lowest doubles of the real parts of b;
 *   brehilolo are the second lowest doubles of the real parts of b;
 *   brelololo are the lowest doubles of the real parts of b;
 *   bimhihihi are the highest doubles of the imaginary parts of b;
 *   bimlohihi are the second highest doubles of the imaginary parts of b;
 *   bimhilohi are the third highest doubles of the imaginary parts of b;
 *   bimlolohi are the fourth highest doubles of the imaginary parts of b;
 *   bimhihilo are the fourth lowest doubles of the imaginary parts of b;
 *   bimlohilo are the third lowest doubles of the imaginary parts of b;
 *   bimhilolo are the second lowest doubles of the imaginary parts of b;
 *   bimlololo are the lowest doubles of the imaginary parts of b;
 *   xrehihihi has space allocated for dim doubles;
 *   xrelohihi has space allocated for dim doubles;
 *   xrehilohi has space allocated for dim doubles;
 *   xrelolohi has space allocated for dim doubles;
 *   ximhihilo has space allocated for dim doubles;
 *   ximlohilo has space allocated for dim doubles;
 *   ximhilolo has space allocated for dim doubles;
 *   ximlololo has space allocated for dim doubles;
 *   xrehihihi has space allocated for dim doubles;
 *   xrelohihi has space allocated for dim doubles;
 *   xrehilohi has space allocated for dim doubles;
 *   xrelolohi has space allocated for dim doubles;
 *   ximhihilo has space allocated for dim doubles;
 *   ximlohilo has space allocated for dim doubles;
 *   ximhilolo has space allocated for dim doubles;
 *   ximlololo has space allocated for dim doubles.
 *
 * ON RETURN :
 *   xrehihihi are the highest doubles of the real parts of the solution x;
 *   xrelohihi are the second highest doubles of the real parts of x;
 *   xrehilohi are the third highest doubles of the real parts of x;
 *   xrelolohi are the fourth highest doubles of the real parts of x;
 *   xrehihilo are the fourth lowest doubles of the real parts of x;
 *   xrelohilo are the third lowest doubles of the real parts of x;
 *   xrehilolo are the second lowest doubles of the real parts of x;
 *   xrelololo are the lowest doubles of the real parts of the x;
 *   ximhihihi are the highest doubles of the imaginary parts of x;
 *   ximlohihi are the second highest doubles of the imaginary parts of x;
 *   ximhilohi are the third highest doubles of the imaginary parts of x;
 *   ximlolohi are the fourth highest doubles of the imaginary parts of x;
 *   ximhihilo are the fourth lowest doubles of the imaginary parts of x;
 *   ximlohilo are the third lowest doubles of the imaginary parts of x;
 *   ximhilolo are the second lowest doubles of the imaginary parts of x;
 *   ximlololo are the lowest doubles of the imaginary parts of x; 
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
