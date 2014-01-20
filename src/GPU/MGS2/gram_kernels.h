// prototypes for the kernels to compute the Gram matrix
// of a sequence of complex vectors

#ifndef __KERNELS__
#define __KERNELS__

#include <iostream>
#include <cmath>
#include <gqd_type.h>
#include "DefineType.h"
#include "complex.h"
#include "vector_types.h"
#include "vector_functions.h"

__global__ void small_gram
 ( complex<T>* v, complex<T>* g, int pivot, int dim, int dimLog2 );
/*
 * DESCRIPTION :
 *   Kernel function to compute the conjugated complex inner product
 *   of the pivot column in v with as second vector the vector in v
 *   at column index equal to the block index, for small dimension.
 *
 * REQUIRED :
 *   The dimension dim equals the block size,
 *   which is the number of threads in the block.
 *
 * ON ENTRY :
 *   v         dim*dim complex numbers for dim vectors of dimension dim;
 *   g         space allocated for dim*dim complex numbers;
 *   pivot     index to the first vector in the conjugated inner product,
 *             the second vector is determined by the block index;
 *   dim       dimension of the vector must equal the block size;
 *   dimLog2   equals ceil(log2((double) dim), used in sum reduction.
 *
 * ON RETURN :
 *   g         for block index b: g[pivot*dim+b] contains the conjugated
 *             complex inner product of vector in column pivot with the
 *             vectors in column b of v. */

__global__ void large_gram
 ( complex<T>* v, complex<T>* g, int pivot, int dim,
   int rnd, int rndLog2, int BS, int BSLog2 );
/*
 * DESCRIPTION :
 *   Kernel function to compute the conjugated complex inner product
 *   of the pivot column in v with as second vector the vector in v 
 *   at column index equal to the block index, for larger dimensions.
 *
 * ON ENTRY :
 *   v         dim*dim complex numbers for dim vectors of dimension dim;
 *   g         space allocated for dim*dim complex numbers;
 *   pivot     index to the first vector in the conjugated inner product,
 *             the second vector is determined by the block index;
 *   rnd       the number of rounds or the multiplier for the dimension;
 *   rndLog2   the 2-logarithm of rnd for use in sum reduction;
 *   BS        the block size or the number of threads in the block;
 *   BSLog2    equals ceil(log2((double) BS), used in sum reduction.
 *
 * ON RETURN :
 *   g         for block index b: g[pivot*dim+b] contains the conjugated
 *             complex inner product of vector in column pivot with the
 *             vectors in column b of v. */

void GPU_gram ( complex<T>* v_h, complex<T> *g_h, int dim, int freq, int BS );
/*
 * DESCRIPTION :
 *   Allocates global memory on the card, transfers the data for the vectors
 *   into the allocated memory, computes the Gram matrix,
 *   and transfers the result from global memory of the card 
 *   to the memory of the host.
 *
 * ON ENTRY :
 *   v_h       contains dim*dim complex numbers for dim vectors of size dim;
 *   g_h       space allocated for dim*dim complex numbers;
 *   dim       dimension of the vectors in v_h;
 *   freq      frequency of the number of kernel launches (for timings);
 *   BS        block size, number of threads per block.
 *
 * ON RETURN :
 *   g_h       the Gram matrix of the vectors represented by v_h. */

#endif
