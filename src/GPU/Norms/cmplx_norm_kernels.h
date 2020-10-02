// Defines the prototypes and constants for the kernels
// to compute the 2-norm of a complex vector in double precision.

#ifndef __CMPLX_NORM_KERNELS_H__
#define __CMPLX_NORM_KERNELS_H__

/*
  The constant d_shmemsize determines the size of the vectors
  stored in shared memory.  As every thread works on one entry
  in the shared memory vectors, for double precision, this size
  is bounded by the number of threads in a block.
  The largest dimension for which the small normalization runs
  is thus the value of d_shmemsize.
 */

#define d_shmemsize 1024

/*
  The constant maxrounds determines the number of rounds
  in the normalization of medium sized vectors.
  The largest dimension for a medium size normalization
  is thus d_shemsize*maxrounds, for instance: 1024*32 = 32768.
 */

#define maxrounds 32

__global__ void small_normalize_vector
 ( double* vre, double* vim, int dim, int dimLog2, double* twonorm );
/*
 * DESCRIPTION :
 *   Kernel function to compute the norm of the complex vector,
 *   given by a vector of real parts and a vector of imaginary parts,
 *   for vectors of small dimension, in double precision.
 *
 * REQUIRED :
 *   The dimension dim equals the block size,
 *   which is the number of threads in the block.
 *
 * ON ENTRY :
 *   vre       real parts of a complex vector of dimension dim;
 *   vim       imaginary parts of a complex vector of dimension dim;
 *   dim       the dimension of the vector must equal the block size;
 *   dimLog2   equals ceil(log2((double) dim), used in sum reduction.
 *
 * ON RETURN :
 *   vre       real parts of the complex vector in same direction,
 *             of the original vector but with norm one;
 *   vim       imaginary parts of the complex vector in same direction,
 *             of the original vector but with norm one;
 *   twonorm   the 2-norm of the given vector. */

__global__ void medium_normalize_vector
 ( double* vre, double* vim, int dim, int rnd, int rndLog2,
   int BS, int BSLog2, double* twonorm );
/*
 * DESCRIPTION :
 *   Kernel function to compute the norm of the vector in v,
 *   for vectors of medium dimension, in double precision.
 *
 * REQUIRED :
 *   The dimension dim equals the block size BS times rnd.
 *   The block size BS is the number of threads in the block.
 *
 * ON ENTRY :
 *   vre       real parts of a complex vector of dimension dim;
 *   vim       imaginary parts of a complex vector of dimension dim;
 *   dim       dimension of the vector must equal rnd*BS;
 *   rnd       the number of rounds or the multiplier for the dimension;
 *   rndLog2   the 2-logarithm of rnd for use in sum reduction;
 *   BS        the block size or the number of threads in the block;
 *   BSLog2    equals ceil(log2((double) BS), used in sum reduction.
 *
 * ON RETURN :
 *   vre       real parts of the complex vector in same direction,
 *             of the original vector but with norm one;
 *   vim       imaginary parts of the complex vector in same direction,
 *             of the original vector but with norm one;
 *   twonorm   the 2-norm of the given vector. */

void GPU_norm
 ( double* vre_h, double* vim_h, int dim, int freq, int BS, double* twonorm );
/*
 * DESCRIPTION :
 *   Allocates global memory on the card, transfers the data for the vector
 *   into the allocated memory, computes the 2-norm, normalizes the vector,
 *   and transfers the result from global memory of the card 
 *   to the memory of the host.
 *
 * ON ENTRY :
 *   vre_h     real parts of a complex vector of dimension dim;
 *   vim_h     imaginary parts of a complex vector of dimension dim;
 *   dim       dimension of the given vector;
 *   freq      frequency of the number of kernel launches (for timings);
 *   BS        block size, number of threads per block.
 *
 * ON RETURN :
 *   vre_h     real parts of the the normalized vector on input;
 *   vim_h     imaginary parts of the the normalized vector on input;
 *   twonorm   the 2-norm of the given vector. */

#endif
