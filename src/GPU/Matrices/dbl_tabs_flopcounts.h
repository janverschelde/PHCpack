/* The file dbl_tabs_flopcounts.h specifies functions to count
 * the floating-point operations of the tiled accelerated back substitution,
 * in double precision. */

/* The names of the functions correspond to the names of the kernels
 * for which the operations are counted.
 * As the long long int type is limited to 64-bits,
 * the accuracy of the counts can only be valid for computations
 * that end in about 10 seconds,
 * in case teraflop performance would be attained. */

#ifndef __dbl_tabs_flopcounts_h__
#define __dbl_tabs_flopcounts_h__

void flopcount_dbl_invert_tiles
 ( int nbt, int szt, long long int *add, long long int *mul,
   long long int *div );
/*
 * DESCRIPTION :
 *   Accumulates the number of floating-point operations to invert
 *   an upper triangular matrix of dimension szt with nbt*szt threads,
 *   on real data.
 *
 * ON ENTRY :
 *   nbt      the number of tiles equals the number of blocks;
 *   szt      equals the number of threads in one block;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications;
 *   mul      current number of divisions.
 *
 * ON RETURN :
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications;
 *   div      accumulated number of divisions. */

void flopcount_cmplx_invert_tiles
 ( int nbt, int szt, long long int *add, long long int *mul,
   long long int *div );
/*
 * DESCRIPTION :
 *   Accumulates the number of floating-point operations to invert
 *   an upper triangular matrix of dimension szt with nbt*szt threads,
 *   on complex data.
 *
 * ON ENTRY :
 *   nbt      the number of tiles equals the number of blocks;
 *   szt      equals the number of threads in one block;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications;
 *   mul      current number of divisions.
 *
 * ON RETURN :
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications;
 *   div      accumulated number of divisions. */

void flopcount_dbl_multiply_inverse
 ( int szt, long long int *add, long long int *mul );
/*
 * DESCRIPTION :
 *   Accumulates the number of floating-point operations to compute the
 *   product of a right hand side slice with an inverse diagonal tile,
 *   on real data.
 *
 * ON ENTRY :
 *   szt      number of threads in a block
 *            and the dimension of the tiles;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications.
 *
 * ON RETURN :
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications. */

void flopcount_cmplx_multiply_inverse
 ( int szt, long long int *add, long long int *mul );
/*
 * DESCRIPTION :
 *   Accumulates the number of floating-point operations to compute the
 *   product of a right hand side slice with an inverse diagonal tile,
 *   on complex data.
 *
 * ON ENTRY :
 *   szt      number of threads in a block
 *            and the dimension of the tiles;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications.
 *
 * ON RETURN :
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications. */

void flopcount_dbl_back_substitute
 ( int nblocks, int szt, long long int *add, long long int *mul );
/*
 * DESCRIPTION :
 *   Accumulates the number of floating-point operations to update
 *   the right hand side vector with back substitution, on real data.
 *
 * ON ENTRY :
 *   nblocks  the number of blocks;
 *   szt      number of threads in each block;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications.
 *
 * ON RETURN :
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications. */

void flopcount_cmplx_back_substitute
 ( int nblocks, int szt, long long int *add, long long int *mul );
/*
 * DESCRIPTION :
 *   Accumulates the number of floating-point operations to update
 *   the right hand side vector with back substitution, on complex data.
 *
 * ON ENTRY :
 *   nblocks  the number of blocks;
 *   szt      number of threads in each block;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications.
 *
 * ON RETURN :
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications. */

#endif
