// Defines the prototypes and constants for the kernels
// to compute the 2-norm of a complex vector in double double precision.

#ifndef __CMPLX2_NORM_KERNELS_H__
#define __CMPLX2_NORM_KERNELS_H__

/*
  The constant dd_shmemsize determines the size of the vectors
  stored in shared memory.  As every thread works on one entry
  in the shared memory vectors, for double double precision,
  this size is bounded by the number of threads in a block.
  The largest dimension for which the small normalization runs
  is thus the value of dd_shmemsize.
 */

#define dd_shmemsize 256

/*
  The constant maxrounds determines the number of rounds
  in the normalization of medium sized vectors.
  The largest dimension for a medium size normalization
  is thus d_shemsize*maxrounds, for instance: 1024*32 = 32768.
 */

#define maxrounds 32

__global__ void small_normalize_vector
 ( double* vrehi, double* vrelo, double* vimhi, double* vimlo,
   int dim, int dimLog2, double* normhi, double* normlo );
/*
 * DESCRIPTION :
 *   Kernel function to compute the norm of the complex vector,
 *   given by a vector of real parts and a vector of imaginary parts,
 *   for vectors of small dimension, in double double precision.
 *
 * REQUIRED :
 *   The dimension dim equals the block size,
 *   which is the number of threads in the block.
 *
 * ON ENTRY :
 *   vrehi     high real parts of a complex vector of dimension dim;
 *   vrelo     low real parts of a complex vector of dimension dim;
 *   vimhi     high imaginary parts of a complex vector of dimension dim;
 *   vimlo     low imaginary parts of a complex vector of dimension dim;
 *   dim       the dimension of the vector must equal the block size;
 *   dimLog2   equals ceil(log2((double) dim), used in sum reduction.
 *
 * ON RETURN :
 *   vrehi     high real parts of the complex vector in same direction,
 *             of the original vector but with norm one;
 *   vrelo     low real parts of the complex vector in same direction,
 *             of the original vector but with norm one;
 *   vimhi     high imaginary parts of the complex vector in same direction,
 *             of the original vector but with norm one;
 *   vimlo     low imaginary parts of the complex vector in same direction,
 *             of the original vector but with norm one;
 *   normhi    high part of the 2-norm of the given vector;
 *   normlo    low part of the 2-norm of the given vector. */

void GPU_norm
 ( double* vrehi_h, double* vrelo_h, double* vimhi_h, double* vimlo_h,
   int dim, int freq, int BS, double* normhi, double* normlo, int blocked );
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
 *   vrehi_h   high real parts of a complex vector of dimension dim;
 *   vrelo_h   low real parts of a complex vector of dimension dim;
 *   vimhi_h   high imaginary parts of a complex vector of dimension dim;
 *   vimlo_h   low imaginary parts of a complex vector of dimension dim;
 *   dim       dimension of the given vector;
 *   freq      frequency of the number of kernel launches (for timings);
 *   BS        block size, number of threads per block;
 *   blocked   is 0 or 1, if 0, then only one block will be launched,
 *             if 1, then as many blocks as dim/BS will be lauched and
 *             dim must be a multiple of BS.
 *
 * ON RETURN :
 *   vrehi_h   high real parts of the the normalized vector on input;
 *   vrelo_h   low real parts of the the normalized vector on input;
 *   vimhi_h   high imaginary parts of the the normalized vector on input;
 *   vimlo_h   low imaginary parts of the the normalized vector on input;
 *   normhi    high part the 2-norm of the given vector;
 *   normlo    low part the 2-norm of the given vector. */

#endif
