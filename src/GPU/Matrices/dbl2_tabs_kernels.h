/* The file dbl2_tabs_kernels.h specifies functions for the
 * tiled accelerated back substitution in double double precision. */

#ifndef __dbl2_tabs_kernels_h__
#define __dbl2_tabs_kernels_h__

#define tabsdd_shmemsize 256

__global__ void dbl2_small_invert_upper 
( int dim, double *Uhi, double *Ulo, double *invUhi, double *invUlo );
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
 *   Uhi      high doubles of an upper triangular matrix stored column wise;
 *   Ulo      low doubles of an upper triangular matrix stored column wise;
 *   invUhi   space allocated for a matrix of dimension dim;
 *   invUlo   space allocated for a matrix of dimension dim.
 *
 * ON RETURN :
 *   invUhi   high doubles of the inverse of the matrix U, stored row wise;
 *   invU     low doubles the inverse of the matrix U, stored row wise. */

__global__ void cmplx2_small_invert_upper
 ( int dim, double *Urehi, double *Urelo, double *Uimhi, double *Uimlo,
   double *invUrehi, double *invUrelo, double *invUimhi, double *invUimlo );
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
 *   Urehi    high doubles of the real parts of U stored column wise;
 *   Urelo    low doubles of the real parts of U stored column wise;
 *   Uimhi    high doubles of the imaginary parts of U;
 *   Uimlo    low doubles of the imaginary parts of U;
 *   invUrehi has space allocated for a matrix of dimension dim;
 *   invUrelo has space allocated for a matrix of dimension dim;
 *   invUimhi has space allocated for a matrix of dimension dim;
 *   invUimlo has space allocated for a matrix of dimension dim.
 *
 * ON RETURN :
 *   invUrehi has the high doubles of the real parts of the inverse of U,
 *            stored row wise;
 *   invUrelo has the low doubles of the real parts of the inverse of U,
 *            stored row wise;
 *   invUimhi has the high doubles of the imaginary parts of the inverse of U,
 *            also stored row wise;
 *   invUimlo has the low doubles of the imaginary parts of the inverse of U,
 *            also stored row wise. */

__global__ void dbl2_medium_invert_upper
 ( int dim, double *Uhi, double *Ulo, double *invUhi, double *invUlo);
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
 *   Uhi      high doubles of an upper triangular matrix stored column wise;
 *   Ulo      low doubles of an upper triangular matrix stored column wise;
 *   invUhi   space allocated for a matrix of dimension dim;
 *   invUlo   space allocated for a matrix of dimension dim.
 *
 * ON RETURN :
 *   invUhi   high doubles of the inverse of the matrix U, stored row wise;
 *   invUlo   low doubles of the inverse of the matrix U, stored row wise. */

__global__ void cmplx2_medium_invert_upper
 ( int dim, double *Urehi, double *Urelo, double *Uimhi, double *Uimlo,
   double *invUrehi, double *invUrelo, double *invUimhi, double *invUimlo );
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
 *   Urehi    high doubles of the real parts of U stored column wise;
 *   Urelo    low doubles of the real parts of U stored column wise;
 *   Uimhi    high doubles of the imaginary parts of U;
 *   Uimlo    low doubles of the imaginary parts of U;
 *   invUrehi has space allocated for a matrix of dimension dim;
 *   invUrelo has space allocated for a matrix of dimension dim;
 *   invUimhi has space allocated for a matrix of dimension dim;
 *   invUimlo has space allocated for a matrix of dimension dim.
 *
 * ON RETURN :
 *   invUrehi has the high doubles of the real parts of the inverse of U,
 *            stored row wise;
 *   invUrelo has the low doubles of the real parts of the inverse of U,
 *            stored row wise;
 *   invUimhi has the high doubles of the imaginary parts of the inverse of U,
 *            also stored row wise;
 *   invUimlo has the low doubles of the imaginary parts of the inverse of U,
 *            also stored row wise. */

__global__ void  dbl2_invert_tiles
 ( int dim, double *Uhi, double *Ulo, double *invUhi, double *invUlo );
/*
 * DESCRIPTION :
 *   Replaces the columns of the tiles with the rows of the inverses.
 *   The number of blocks equals the number of tiles in U.
 *   The number of threads per block equals the dimension of each tile.
 *
 * REQUIRED : dim <= 256 = dd_shmemsize.
 *
 * ON ENTRY :
 *   dim      the dimension of each tile;
 *   Uhi      high doubles of columns of all tiles on the diagonal 
 *            of an upper triangular matrix;
 *   Ulo      low doubles of columns of all tiles on the diagonal 
 *            of an upper triangular matrix;
 *   invUhi   space allocated for the inverse of all tiles in U;
 *   invUlo   space allocated for the inverse of all tiles in U.
 *
 * ON RETURN :
 *   invUhi   high doubles of the rows of the inverse of the tiles in U;
 *   invUlo   low doubles of the rows of the inverse of the tiles in U. */

__global__ void  cmplx2_invert_tiles
 ( int dim, double *Urehi, double *Urelo, double *Uimhi, double *Uimlo,
   double *invUrehi, double *invUrelo, double *invUimhi, double *invUimlo );
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
 *   Urehi    high doubles of the real parts of the columns of all tiles
 *            on the diagonal of an upper triangular matrix;
 *   Urelo    low doubles of the real parts of the columns of all tiles
 *            on the diagonal of an upper triangular matrix;
 *   Uimhi    high doubles of the imaginary parts of the columns of all
 *            tiles on the diagonal of an upper triangular matrix;
 *   Uimlo    low doubles of the imaginary parts of the columns of all
 *            tiles on the diagonal of an upper triangular matrix;
 *   invUrehi has space allocated for the high doubles of the real parts
 *            of the inverse of all tiles in U;
 *   invUrelo has space allocated for the low doubles of the real parts
 *            of the inverse of all tiles in U;
 *   invUimhi has space allocated for the high doubles of the imaginary
 *            parts of the inverse of all tiles in U;
 *   invUimlo has space allocated for the low doubles of the imaginary
 *            parts of the inverse of all tiles in U.
 *
 * ON RETURN :
 *   invUrehi has the high doubles of the real parts
 *            of rows of the inverse tiles;
 *   invUrelo has the low doubles of the real parts
 *            of rows of the inverse tiles;
 *   invUimhi has the high doubles of the imaginary parts
 *            of rows of the inverse tiles;
 *   invUimlo has the low doubles of the imaginary parts
 *            of rows of the inverse tiles. */

__global__ void dbl2_multiply_inverse
 ( int dim, int idx, double *invUhi, double *invUlo,
   double *whi, double *wlo );
/*
 * DESCRIPTION :
 *   Replaces b with the product of the inverse tile in U.
 *
 * ON ENTRY :
 *   dim      the dimension of each tile;
 *   idx      index of the diagonal tile;
 *   invUhi   high doubles of the inverse of the diagonal tiles;
 *   invUlo   low doubles of the inverse of the diagonal tiles;
 *   whi      high doubles of the right hand side vector;
 *   wlo      low doubles of the right hand side vector.
 *
 * ON RETURN :
 *   whi      high doubles of the product of the inverse of the diagonal
 *            tile defined by the index idx with the w on input;
 *   wlo      low doubles of the product of the inverse of the diagonal
 *            tile defined by the index idx with the w on input. */

__global__ void cmplx2_multiply_inverse
 ( int dim, int idx,
   double *invUrehi, double *invUrelo, double *invUimhi, double *invUimlo,
   double *wrehi, double *wrelo, double *wimhi, double *wimlo );
/*
 * DESCRIPTION :
 *   Replaces b with the product of the inverse tile in U.
 *
 * ON ENTRY :
 *   dim      the dimension of each tile;
 *   idx      index of the diagonal tile;
 *   invUrehi are the high doubles of the real parts of the inverse
 *            of the diagonal tiles;
 *   invUrelo are the low doubles of the real parts of the inverse
 *            of the diagonal tiles;
 *   invUimhi are the high doubles of the imaginary parts of the inverse
 *            of the diagonal tiles;
 *   invUimlo are the low doubles of the imaginary parts of the inverse
 *            of the diagonal tiles;
 *   wrehi    high doubles of the real parts of the right hand side;
 *   wrelo    low doubles of the real parts of the right hand side;
 *   wimhi    high doubles of the imaginary parts of the right hand side;
 *   wimlo    low doubles of the imaginary parts of the right hand side.
 *
 * ON RETURN :
 *   wrehi    high doubles of the real parts of the product of the inverse
 *            of the tile defined by the index idx with the input w;
 *   wrelo    low doubles of the real parts of the product of the inverse
 *            of the tile defined by the index idx with the input w;
 *   wimhi    high doubles of the imaginary parts of the product of the inverse
 *            of the tile defined by the index idx with the input w;
 *   wimlo    low doubles of the imaginary parts of the product of the inverse
 *            of the tile defined by the index idx with the input w. */

__global__ void dbl2_back_substitute
 ( int dim, int idx, double *Uhi, double *Ulo, double *whi, double *wlo );
/*
 * DESCRIPTION :
 *   Updates the right hand side vector subtracting the solution
 *   defined by idx, multiplied with the corresponding rows in U.
 *
 * ON ENTRY :
 *   dim      dimension of each tile;
 *   idx      index of the solution tile in the multiplication;
 *   Uhi      high doubles of tiles to multiply the solution with;
 *   Ulo      low doubles of tiles to multiply the solution with;
 *   whi      high doubles of the current right hand side vector;
 *   wlo      low doubles of the current right hand side vector.
 *
 * ON RETURN :
 *   whi      high doubles of the updated right hand side vector;
 *   wlo      low doubles of the updated right hand side vector. */

__global__ void cmplx2_back_substitute
 ( int dim, int idx,
   double *Urehi, double *Urelo, double *Uimhi, double *Uimlo,
   double *wrehi, double *wrelo, double *wimhi, double *wimlo );
/*
 * DESCRIPTION :
 *   Updates the right hand side vector subtracting the solution
 *   defined by idx, multiplied with the corresponding rows in U.
 *
 * ON ENTRY :
 *   dim      dimension of each tile;
 *   idx      index of the solution tile in the multiplication;
 *   Urehi    high doubles of the real parts of tiles 
 *            to multiply the solution with;
 *   Urelo    low doubles of the real parts of tiles
 *            to multiply the solution with;
 *   Uimhi    high doubles of the imaginary parts of tiles
 *            to multiply the solution with;
 *   Uimlo    low doubles of the imaginary parts of tiles
 *            to multiply the solution with;
 *   wrehi    high doubles of the real parts of 
 *            the current right hand side vector;
 *   wrelo    low doubles of the real parts of 
 *            the current right hand side vector;
 *   wimhi    high doubles of the imaginary parts
 *            of the current right hand side vector;
 *   wimlo    high doubles of the imaginary parts
 *            of the current right hand side vector.
 *
 * ON RETURN :
 *   wrehi    high doubles of the real parts
 *            of the updated right hand side vector;
 *   wrelo    low doubles of the real parts
 *            of the updated right hand side vector;
 *   wimhi    high doubles of the imaginary parts
 *            of the updated right hand side vector;
 *   wimlo    low doules of the imaginary parts
 *            of the updated right hand side vector. */

void GPU_dbl2_upper_inverse
 ( int dim, double **Uhi, double **Ulo, double **invUhi, double **invUlo,
   double *lapms, double *walltimesec );
/*
 * DESCRIPTION :
 *   Calls the kernel to invert the upper triangular matrix U.
 *   The matrices are stored in the conventional rowwise fashion.
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
 *   invUlo   low doubles of the inverse of the matrix U;
 *   lapms    elapsed time spent by the kernels;
 *   walltimesec is the elapsed wall clock computation time. */

void GPU_cmplx2_upper_inverse
 ( int dim, double **Urehi, double **Urelo, double **Uimhi, double **Uimlo,
   double **invUrehi, double **invUrelo,
   double **invUimhi, double **invUimlo, double *lapms, double *walltimesec );
/*
 * DESCRIPTION :
 *   Calls the kernel to invert the upper triangular matrix U.
 *   The matrices are stored in the conventional row wise fashion.
 *
 * ON ENTRY :
 *   dim      dimension of the upper triangular matrix U;
 *   Urehi    high doubles of the real parts of U;
 *   Urelo    low doubles of the real parts of U;
 *   Uimhi    high doubles of the imaginary parts of U;
 *   Uimlo    low doubles of the imaginary parts of U;
 *   invUrehi has space allocated for a matrix of dimension dim;
 *   invUrelo has space allocated for a matrix of dimension dim;
 *   invUimhi has space allocated for a matrix of dimension dim;
 *   invUimlo has space allocated for a matrix of dimension dim.
 *
 * ON RETURN :
 *   invUrehi has the high doubles of the real parts of the inverse;
 *   invUrelo has the low doubles of the real parts of the inverse;
 *   invUimhi has the high doubles of the imaginary parts of the inverse;
 *   invUimlo has the low doubles of the imaginary parts of the inverse;
 *   lapms    elapsed time spent by the kernels;
 *   walltimesec is the elapsed wall clock computation time. */

void GPU_dbl2_upper_tiled_solver
 ( int dim, int szt, int nbt, double **Uhi, double **Ulo,
   double *bhi, double *blo, double *xhi, double *xlo,
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
 *   Uhi      high doubles of an upper triangular matrix of dimension dim;
 *   Ulo      low doubles of an upper triangular matrix of dimension dim;
 *   bhi      high doubles of the right hand side of the linear system;
 *   blo      low doubles of the right hand side of the linear system;
 *   xhi      space allocated for dim doubles;
 *   xlo      space allocated for dim doubles.
 *
 * ON RETURN :
 *   xhi      high doubles of the solution to the system U*x = b;
 *   xlo      low doubles of the solution to the system U*x = b;
 *   invlapms is the elapsed time spent by the kernel to invert a tile;
 *   mullapms is the elapsed time spent by the kernel to multiply
 *            with the inversed diagonal tile;
 *   sublapms is the elapsed time spent by the kernel for back substitution;
 *   totlapms is the total elapsed time spent by all kernels;
 *   walltimesec is the elapsed wall clock computation time;
 *   addcnt   counts the number of additions and subtractions;
 *   mulcnt   counts the number of multiplications;
 *   divcnt   counts the number of divisions. */

void GPU_cmplx2_upper_tiled_solver
 ( int dim, int szt, int nbt,
   double **Urehi, double **Urelo, double **Uimhi, double **Uimlo,
   double *brehi, double *brelo, double *bimhi, double *bimlo,
   double *xrehi, double *xrelo, double *ximhi, double *ximlo,
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
 *   Urehi    high doubles of the real parts of U;
 *   Urelo    low doubles of the real parts of U;
 *   Uimhi    high doubles of the imaginary parts of U;
 *   Uimlo    low doubles of the imaginary parts of U;
 *   brehi    high doubles of the real parts of the right hand side;
 *   brelo    low doubles of the real parts of the right hand side;
 *   bimhi    high doubles of the imaginary parts of the right hand side;
 *   bimlo    low doubles of the imaginary parts of the right hand side;
 *   xrehi    space allocated for dim doubles;
 *   xrelo    space allocated for dim doubles;
 *   ximhi    space allocated for dim doubles;
 *   ximlo    space allocated for dim doubles.
 *
 * ON RETURN :
 *   xrehi    high doubles of the real parts of the solution;
 *   xrelo    low doubles of the real parts of the solution;
 *   ximhi    high doubles of the imaginary parts of the solution;
 *   ximlo    low doubles of the imaginary parts of the solution;
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
