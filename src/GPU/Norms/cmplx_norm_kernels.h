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
  The largest dimension for which the large normalization runs
  is the square of the value of d_shmemsize because the algorithm
  sums two arrays of squares, the squares of the real and the
  squares of the imaginary part, thus very similar to case of
  the small normalization.  For example, for d_shmemsize = 1024,
  the largest dimension is thus 1024*1024 = 1048576.
 */

#define d_shmemsize 1024

/*
  The constant maxrounds determines the number of rounds
  in the normalization of medium sized vectors.
  The largest dimension for a medium size normalization
  is thus d_shemsize*maxrounds, for instance: 1024*512 = 524288.
 */

#define maxrounds 512

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

__global__ void large_sum_the_squares
 ( double* vre, double* vim, int dim, double* sums, int BS, int BSLog2 );
/*
 * DESCRIPTION :
 *   Computes the sums of the squares of the numbers in v,
 *   as needed in the 2-norm of the vector, with many blocks.
 *
 * REQUIRED :
 *   Space in sums is allocated for as many blocks in the launch.
 *
 * ON ENTRY :
 *   vre       real parts of a complex vector of dimension dim;
 *   vim       imaginary parts of a complex vector of dimension dim;
 *   dim       dimension of the vector must equal BS*nbsums;
 *   BS        the block size or the number of threads in the block;
 *   BSLog2    equals ceil(log2((double) BS), used in sum reduction.
 *
 * ON RETURN :
 *   sums      computed sums of squares of slices of v. */

__global__ void large_normalize_vector
 ( double* vre, double* vim, int dim, double* sums,
   int nbsums, int nbsumsLog2, int BS, double* twonorm );
/*
 * DESCRIPTION :
 *   Computes the norm of the vector v for vectors of large dimension,
 *   for given sums of squares, and normalizes the vector.
 *
 * REQUIRED :
 *   The dimension dim equals the block size BS times nbsums.
 *   The block size BS is the number of threads in the block.
 *   The number of blocks equals nbsums.
 *
 * ON ENTRY :
 *   vre       real parts of a complex vector of dimension dim;
 *   vim       imaginary parts of a complex vector of dimension dim;
 *   v         vector of doubles of dimension dim;
 *   dim       dimension of the vector must equal BS*nbsums;
 *   sums      computed sums of squares of slices of v;
 *   nbsums    the number of elements in sums equals
 *             the number of blocks in the kernel launch;
 *   nbsumsLog2 is ceil(log2((double) nbsums), used in sum reduction;
 *   BS        is the block size, the number of threads in a block.
 *
 * ON RETURN :
 *   vre       real parts of the complex vector in same direction,
 *             of the original vector but with norm one;
 *   vim       imaginary parts of the complex vector in same direction,
 *             of the original vector but with norm one;
 *   twonorm   the 2-norm of the original vector v. */

void GPU_norm
 ( double* vre_h, double* vim_h, int dim, int freq, int BS,
   double* twonorm, int blocked );
/*
 * DESCRIPTION :
 *   Allocates global memory on the card, transfers the data for the vector
 *   into the allocated memory, computes the 2-norm, normalizes the vector,
 *   and transfers the result from global memory of the card 
 *   to the memory of the host.
 *
 * REQUIRED :
 *   The dimension dim is a multiple of the block size BS.
 *
 * ON ENTRY :
 *   vre_h     real parts of a complex vector of dimension dim;
 *   vim_h     imaginary parts of a complex vector of dimension dim;
 *   dim       dimension of the given vector;
 *   freq      frequency of the number of kernel launches (for timings);
 *   BS        block size, number of threads per block;
 *   blocked   is 0 or 1, if 0, then only one block will be launched,
 *             if 1, then as many blocks as dim/BS will be lauched and
 *             dim must be a multiple of BS.
 *
 * ON RETURN :
 *   vre_h     real parts of the the normalized vector on input;
 *   vim_h     imaginary parts of the the normalized vector on input;
 *   twonorm   the 2-norm of the given vector. */

#endif
