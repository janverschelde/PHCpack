// Defines the prototypes and constants for the kernels
// to compute the 2-norm of a complex vector in triple double precision.

#ifndef __CMPLX3_NORM_KERNELS_H__
#define __CMPLX3_NORM_KERNELS_H__

/*
  The constant td_shmemsize determines the size of the vectors
  stored in shared memory.  As every thread works on one entry
  in the shared memory vectors, for triple double precision,
  this size is bounded by the number of threads in a block.
  The largest dimension for which the small normalization runs
  is thus the value of td_shmemsize.
  The largest dimension for which the large normalization runs
  is the square of the value of td_shmemsize because the algorithm
  sums two arrays of squares, the squares of the real and the
  squares of the imaginary part, thus very similar to case of
  the small normalization.  For example, for td_shmemsize = 400,
  the largest dimension is thus 400*400 = 160000.
 */

#define td_shmemsize 400

/*
  The constant maxrounds determines the number of rounds
  in the normalization of medium sized vectors.
  The largest dimension for a medium size normalization
  is thus td_shemsize*maxrounds, for instance: 400*320 = 128000.
 */

#define maxrounds 320

__global__ void small_normalize_vector
 ( double *vrehi, double *vremi, double *vrelo,
   double *vimhi, double *vimmi, double *vimlo, int dim,
   int dimLog2, double *normhi, double *normmi, double *normlo );
/*
 * DESCRIPTION :
 *   Kernel function to compute the norm of the complex vector,
 *   given by a vector of real parts and a vector of imaginary parts,
 *   for vectors of small dimension, in triple double precision.
 *
 * REQUIRED :
 *   The dimension dim equals the block size,
 *   which is the number of threads in the block.
 *
 * ON ENTRY :
 *   vrehi     high real parts of a complex vector of dimension dim;
 *   vremi     middle real parts of a complex vector of dimension dim;
 *   vrelo     low real parts of a complex vector of dimension dim;
 *   vimhi     high imaginary parts of a complex vector of dimension dim;
 *   vimmi     middle imaginary parts of a complex vector of dimension dim;
 *   vimlo     low imaginary parts of a complex vector of dimension dim;
 *   dim       the dimension of the vector must equal the block size;
 *   dimLog2   equals ceil(log2((double) dim), used in sum reduction.
 *
 * ON RETURN :
 *   vrehi     high real parts of the complex vector in same direction,
 *             of the original vector but with norm one;
 *   vremi     middle real parts of the complex vector in same direction,
 *             of the original vector but with norm one;
 *   vrelo     low real parts of the complex vector in same direction,
 *             of the original vector but with norm one;
 *   vimhi     high imaginary parts of the complex vector in same direction,
 *             of the original vector but with norm one;
 *   vimmi     middle imaginary parts of the complex vector in same direction,
 *             of the original vector but with norm one;
 *   vimlo     low imaginary parts of the complex vector in same direction,
 *             of the original vector but with norm one;
 *   normhi    high part of the 2-norm of the given vector;
 *   normmi    middle part of the 2-norm of the given vector;
 *   normlo    low part of the 2-norm of the given vector. */

__global__ void medium_normalize_vector
 ( double *vrehi, double *vremi, double *vrelo,
   double *vimhi, double *vimmi, double *vimlo, int dim,
   int rnd, int rndLog2, int BS, int BSLog2,
   double *normhi, double *normmi, double *normlo );
/*
 * DESCRIPTION :
 *   Kernel function to compute the norm of the vector in v,
 *   for vectors of medium dimension, in triple double precision.
 *
 * REQUIRED :
 *   The dimension dim equals the block size BS times rnd.
 *   The block size BS is the number of threads in the block.
 *
 * ON ENTRY :
 *   vrehi     high real parts of a complex vector of dimension dim;
 *   vremi     middle real parts of a complex vector of dimension dim;
 *   vrelo     low real parts of a complex vector of dimension dim;
 *   vimhi     high imaginary parts of a complex vector of dimension dim;
 *   vimmi     middle imaginary parts of a complex vector of dimension dim;
 *   vimlo     low imaginary parts of a complex vector of dimension dim;
 *   dim       dimension of the vector must equal rnd*BS;
 *   rnd       the number of rounds or the multiplier for the dimension;
 *   rndLog2   the 2-logarithm of rnd for use in sum reduction;
 *   BS        the block size or the number of threads in the block;
 *   BSLog2    equals ceil(log2((double) BS), used in sum reduction.
 *
 * ON RETURN :
 *   vrehi     high real parts of the complex vector in same direction,
 *             of the original vector but with norm one;
 *   vremi     middle real parts of the complex vector in same direction,
 *             of the original vector but with norm one;
 *   vrelo     low real parts of the complex vector in same direction,
 *             of the original vector but with norm one;
 *   vimhi     high imaginary parts of the complex vector in same direction,
 *             of the original vector but with norm one;
 *   vimmi     middle imaginary parts of the complex vector in same direction,
 *             of the original vector but with norm one;
 *   vimlo     high imaginary parts of the complex vector in same direction,
 *             of the original vector but with norm one;
 *   normhi    high part of the 2-norm of the given vector;
 *   normmi    middle part of the 2-norm of the given vector;
 *   normlo    low part of the 2-norm of the given vector. */

__global__ void large_sum_the_squares
 ( double *vrehi, double *vremi, double *vrelo,
   double *vimhi, double *vimmi, double *vimlo, int dim,
   double *sumshi, double *sumsmi, double *sumslo, int BS, int BSLog2 );
/*
 * DESCRIPTION :
 *   Computes the sums of the squares of the numbers in v,
 *   as needed in the 2-norm of the vector, with many blocks.
 *
 * REQUIRED :
 *   Space in sums is allocated for as many blocks in the launch.
 *
 * ON ENTRY :
 *   vrehi     high real parts of a complex vector of dimension dim;
 *   vremi     middle real parts of a complex vector of dimension dim;
 *   vrelo     low real parts of a complex vector of dimension dim;
 *   vimhi     high imaginary parts of a complex vector of dimension dim;
 *   vimmi     middle imaginary parts of a complex vector of dimension dim;
 *   vimlo     low imaginary parts of a complex vector of dimension dim;
 *   dim       dimension of the vector must equal BS*nbsums;
 *   BS        the block size or the number of threads in the block;
 *   BSLog2    equals ceil(log2((double) BS), used in sum reduction.
 *
 * ON RETURN :
 *   sumshi    high parts of sums of squares of slices of v;
 *   sumsmi    middle parts of sums of squares of slices of v;
 *   sumslo    low parts of sums of squares of slices of v. */

__global__ void large_normalize_vector
 ( double *vrehi, double *vremi, double *vrelo,
   double *vimhi, double *vimmi, double *vimlo, int dim,
   double *sumshi, double *sumsmi, double *sumslo,
   int nbsums, int nbsumsLog2, int BS,
   double *normhi, double *normmi, double *normlo );
/*
 * DESCRIPTION :
 *   Computes the norm of the vector v for vectors of large dimension,
 *   for given sums of squares, and normalizes the vector.
 *
 * REQUIRED :
 *   The dimension dim equals the block size BS times nbsums.
 *   The block size BS is the number of threads in the block.
 *   The number of blocks equals nbsums.
 *   All arrays have size dim.
 *
 * ON ENTRY :
 *   vrehi     high real parts of a complex vector;
 *   vremi     middle real parts of a complex vector;
 *   vrelo     low real parts of a complex vector;
 *   vimhi     high imaginary parts of a complex vector;
 *   vimmi     middle imaginary parts of a complex vector;
 *   vimlo     low imaginary parts of a complex vector;
 *   dim       dimension of the vector must equal BS*nbsums;
 *   sumshi    high parts of sums of squares of slices of v;
 *   sumsmi    middle parts of sums of squares of slices of v;
 *   sumslo    low parts of sums of squares of slices of v;
 *   nbsums    the number of elements in sums equals
 *             the number of blocks in the kernel launch;
 *   nbsumsLog2 is ceil(log2((double) nbsums), used in sum reduction;
 *   BS        is the block size, the number of threads in a block.
 *
 * ON RETURN :
 *   vrehi     high real parts of the complex vector in same direction,
 *             of the original vector but with norm one;
 *   vremi     middle real parts of the complex vector in same direction,
 *             of the original vector but with norm one;
 *   vrelo     low real parts of the complex vector in same direction,
 *             of the original vector but with norm one;
 *   vimhi     high imaginary parts of the complex vector in same direction,
 *             of the original vector but with norm one;
 *   vimhi     middle imaginary parts of the complex vector in same direction,
 *             of the original vector but with norm one;
 *   vimlo     low imaginary parts of the complex vector in same direction,
 *             of the original vector but with norm one;
 *   normhi    high part of the 2-norm of the original vector v;
 *   normhi    middle part of the 2-norm of the original vector v;
 *   normlo    low part of the 2-norm of the original vector v. */

void GPU_norm
 ( double *vrehi_h, double *vremi_h, double *vrelo_h,
   double *vimhi_h, double *vimmi_h, double *vimlo_h,
   int dim, int freq, int BS,
   double *normhi, double *normmi, double *normlo, int blocked );
/*
 * DESCRIPTION :
 *   Allocates global memory on the card, transfers the data for the vector
 *   into the allocated memory, computes the 2-norm, normalizes the vector,
 *   and transfers the result from global memory of the card 
 *   to the memory of the host.
 *
 * REQUIRED :
 *   The dimension dim is a multiple of the block size BS.
 *   All arrays have size dim.
 *
 * ON ENTRY :
 *   vrehi_h   high real parts of a complex vector;
 *   vremi_h   middle real parts of a complex vector;
 *   vrelo_h   low real parts of a complex vector;
 *   vimhi_h   high imaginary parts of a complex vector;
 *   vimmi_h   middle imaginary parts of a complex vector;
 *   vimlo_h   low imaginary parts of a complex vector;
 *   dim       dimension of the given vector;
 *   freq      frequency of the number of kernel launches (for timings);
 *   BS        block size, number of threads per block;
 *   blocked   is 0 or 1, if 0, then only one block will be launched,
 *             if 1, then as many blocks as dim/BS will be lauched and
 *             dim must be a multiple of BS.
 *
 * ON RETURN :
 *   vrehi_h   high real parts of the normalized vector;
 *   vremi_h   middle real parts of the normalized vector;
 *   vrelo_h   low real parts of the normalized vector;
 *   vimhi_h   high imaginary parts of the normalized vector;
 *   vimmi_h   middle imaginary parts of the normalized vector;
 *   vimlo_h   low imaginary parts of the normalized vector;
 *   normhi    high part the 2-norm of the given vector;
 *   normmi    middle part the 2-norm of the given vector;
 *   normlo    low part the 2-norm of the given vector. */

#endif
