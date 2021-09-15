/* The file dbl8_tabs_kernels.h specifies functions for the
 * tiled accelerated back substitution in octo double precision. */

#ifndef __dbl8_tabs_kernels_h__
#define __dbl8_tabs_kernels_h__

#define tabsod_shmemsize 128

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
 * REQUIRED : dim <= 16.
 *   Because the inverse is stored entirely in shared memory,
 *   the dimension dim is limited to 16 = 2^4, as 16^2 = 256,
 *   the upper limit on the shared memory, d_shmemsize.
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
 * REQUIRED : dim <= 16.
 *   Because the inverse is stored entirely in shared memory,
 *   the dimension dim is limited to 16 = 2^4, as 16^2 = 256,
 *   the upper limit on the shared memory, d_shmemsize.
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
 * REQUIRED : dim <= 256.
 *   Because the columns of U are loaded entirely into shared memory
 *   and the rows of the inverses are computed first entirely in
 *   shared memory before storing, the dimension dim is limited 
 *   to 256, the upper limit on the shared memory, dd_shmemsize.
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
 * REQUIRED : dim <= 256.
 *   Because the columns of U are loaded entirely into shared memory
 *   and the rows of the inverses are computed first entirely in
 *   shared memory before storing, the dimension dim is limited 
 *   to 256, the upper limit on the shared memory, od_shmemsize.
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

/*

__global__ void  dbl8_invert_tiles
 ( int dim, double *Uhihi, double *Ulohi, double *Uhilo, double *Ulolo,
   double *invUhihi, double *invUlohi, double *invUhilo, double *invUlolo);
*
 * DESCRIPTION :
 *   Replaces the columns of the tiles with the rows of the inverses.
 *   The number of blocks equals the number of tiles in U.
 *   The number of threads per block equals the dimension of each tile.
 *
 * REQUIRED : dim <= 256 = dd_shmemsize.
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
 *   invUhi   space allocated for the inverse of all tiles in U;
 *   invUhi   space allocated for the inverse of all tiles in U;
 *   invUlo   space allocated for the inverse of all tiles in U;
 *   invUlo   space allocated for the inverse of all tiles in U.
 *
 * ON RETURN :
 *   invUhihi are the highest doubles of the inverse of the tiles in U;
 *   invUlohi are the second highest doubles of the inverse of the tiles in U;
 *   invUhilo are the second lowest doubles of the inverse of the tiles in U;
 *   invUlolo are the lowest doubles of the inverse of the tiles in U. *

__global__ void cmplx8_invert_tiles
 ( int dim, 
   double *Urehihi, double *Urelohi, double *Urehilo, double *Urelolo,
   double *Uimhihi, double *Uimlohi, double *Uimhilo, double *Uimlolo,
   double *invUrehihi, double *invUrelohi,
   double *invUrehilo, double *invUrelolo,
   double *invUimhihi, double *invUimlohi,
   double *invUimhilo, double *invUimlolo );
*
 * DESCRIPTION :
 *   Replaces the columns of the tiles with the rows of the inverses.
 *   The number of blocks equals the number of tiles in U.
 *   The number of threads per block equals the dimension of each tile.
 *
 * REQUIRED : dim <= 256 = d_shmemsize.
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
 *            of rows of the inverse tiles. *

__global__ void dbl8_multiply_inverse
 ( int dim, int idx,
   double *invUhihi, double *invUlohi, double *invUhilo, double *invUlolo,
   double *whihi, double *wlohi, double *whilo, double *wlolo );
*
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
 *            diagonal tile defined by the index idx with the w on input. *

__global__ void cmplx8_multiply_inverse
 ( int dim, int idx,
   double *invUrehihi, double *invUrelohi,
   double *invUrehilo, double *invUrelolo,
   double *invUimhihi, double *invUimlohi,
   double *invUimhilo, double *invUimlolo,
   double *wrehihi, double *wrelohi, double *wrehilo, double *wrelolo,
   double *wimhihi, double *wimlohi, double *wimhilo, double *wimlolo );
*
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
 *   wrehi    highest doubles of the real parts of the product of the
 *            inverse of the tile defined by the index idx with the input w;
 *   wrehi    second highest doubles of the real parts of the product
 *            of the inverse of the tile defined by the index idx;
 *   wrelo    second lowest doubles of the real parts of the product
 *            of the inverse of the tile defined by the index idx;
 *   wrelo    lowest doubles of the real parts of the product
 *            of the inverse of the tile defined by the index idx;
 *   wimhi    highest doubles of the imaginary parts of the product
 *            of the inverse of the tile defined by the index idx;
 *   wimhi    second highest doubles of the imaginary parts of the product
 *            of the inverse of the tile defined by the index idx;
 *   wimlo    second lowest doubles of the imaginary parts of the product
 *            of the inverse of the tile defined by the index idx;
 *   wimlo    lowest doubles of the imaginary parts of the product
 *            of the inverse of the tile defined by the index idx. *

__global__ void dbl8_back_substitute
 ( int dim, int idx, 
   double *Uhihi, double *Ulohi, double *Uhilo, double *Ulolo, 
   double *whihi, double *wlohi, double *whilo, double *wlolo );
*
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
 *   wlolo    lowest doubles of the updated right hand side vector. *

__global__ void cmplx8_back_substitute
 ( int dim, int idx,
   double *Urehihi, double *Urelohi, double *Urehilo, double *Urelolo,
   double *Uimhihi, double *Uimlohi, double *Uimhilo, double *Uimlolo,
   double *wrehihi, double *wrelohi, double *wrehilo, double *wrelolo,
   double *wimhihi, double *wimlohi, double *wimhilo, double *wimlolo );
*
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

/*

void GPU_dbl8_upper_tiled_solver
 ( int dim, int szt, int nbt,
   double **Uhihi, double **Ulohi, double **Uhilo, double **Ulolo,
   double *bhihi, double *blohi, double *bhilo, double *blolo,
   double *xhihi, double *xlohi, double *xhilo, double *xlolo,
   double *invlapms, double *mullapms, double *sublapms, double *totlapms,
   double *walltimesec,
   long long int *addcnt, long long int *mulcnt, long long int *divcnt );
*
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
 *   mulcnt   counts the number of multiplications;
 *   divcnt   counts the number of divisions. *

void GPU_cmplx8_upper_tiled_solver
 ( int dim, int szt, int nbt,
   double **Urehihi, double **Urelohi, double **Urehilo, double **Urelolo,
   double **Uimhihi, double **Uimlohi, double **Uimhilo, double **Uimlolo,
   double *brehihi, double *brelohi, double *brehilo, double *brelolo,
   double *bimhihi, double *bimlohi, double *bimhilo, double *bimlolo,
   double *xrehihi, double *xrelohi, double *xrehilo, double *xrelolo,
   double *ximhihi, double *ximlohi, double *ximhilo, double *ximlolo,
   double *invlapms, double *mullapms, double *sublapms, double *totlapms,
   double *walltimesec,
   long long int *addcnt, long long int *mulcnt, long long int *divcnt );
*
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
 *   xrehi    highest doubles of the real parts of the solution;
 *   xrehi    second highest doubles of the real parts of the solution;
 *   xrelo    second lowest doubles of the real parts of the solution;
 *   xrelo    lowest doubles of the real parts of the solution;
 *   ximhi    highest doubles of the imaginary parts of the solution;
 *   ximhi    second highest doubles of the imaginary parts of the solution;
 *   ximlo    second lowest doubles of the imaginary parts of the solution;
 *   ximlo    lowest doubles of the imaginary parts of the solution;
 *   invlapms is the elapsed time spent by the kernel to invert a tile;
 *   mullapms is the elapsed time spent by the kernel to multiply
 *            with the inversed diagonal tile;
 *   sublapms is the elapsed time spent by the kernel for back substitution;
 *   totlapms is the total elapsed time spent by all kernels;
 *   walltimesec is the elapsed wall clock computation time;
 *   addcnt   counts the number of additions and subtractions;
 *   mulcnt   counts the number of multiplications;
 *   divcnt   counts the number of divisions. *

*/

#endif
