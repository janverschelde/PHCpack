// The file dbl_bals_flopcounts.h defines the prototypes of functions to
// count the number of floating-point operations in the kernels to
// solve a linear system of power series, in double precision.

#ifndef __dbl_bals_kernels_h__
#define __dbl_bals_kernels_h__

void flopcount_dbl_bals_tail
 ( int dim, long long int *add, long long int *mul );
/*
 * DESCRIPTION :
 *   Returns the accumulated number of floating-point operations
 *   executed by dbl_bals_tail.
 *
 * ON ENTRY :
 *   dim      the dimension of the matrix and vectors;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications.
 *
 * ON RETURN :
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications. */

void flopcount_cmplx_bals_tail
 ( int dim, long long int *add, long long int *mul );
/*
 * DESCRIPTION :
 *   Returns the accumulated number of floating-point operations
 *   executed by cmplx_bals_tail.
 *
 * ON ENTRY :
 *   dim      the dimension of the matrix and vectors;
 *   add      current number of additions and subtractions;
 *   mul      current number of multiplications.
 *
 * ON RETURN :
 *   add      accumulated number of additions and subtractions;
 *   mul      accumulated number of multiplications. */

void prompt_flopbals_setup
 ( int *seed, int *dim, int *deg, int *szt, int *nbt, int *cdata );
/*
 * DESCRIPTION :
 *   Prompts for the parameters in the testing of the flops. 
 *
 * ON RETURN :
 *   seed     if 0, then random seed, otherwise fixed seed;
 *   dim      dimension of the linear system;
 *   deg      degree of truncation;
 *   szt      number of threads in a block;
 *   nbt      number of blocks, szt*nbt must equal dim;
 *   cdata    if 0 for real data, 1 for complex data. */

#endif
