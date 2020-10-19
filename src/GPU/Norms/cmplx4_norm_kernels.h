// Defines the prototypes and constants for the kernels
// to compute the 2-norm of a complex vector in quad double precision.

#ifndef __CMPLX4_NORM_KERNELS_H__
#define __CMPLX4_NORM_KERNELS_H__

/*
  The constant qd_shmemsize determines the size of the vectors
  stored in shared memory.  As every thread works on one entry
  in the shared memory vectors, for quad double precision,
  thihis size is bounded by the number of threads in a block.
  The largest dimension for which the small normalization runs
  is thus the value of qd_shmemsize.
  The largest dimension for which the large normalization runs
  is the square of the value of qd_shmemsize because the algorithm
  sums two arrays of squares, the squares of the real and the
  squares of the imaginary part, thus very similar to case of
  the small normalization.  For example, for od_shmemsize = 320,
  the largest dimension is thus 320*320 = 102400.
 */

#define qd_shmemsize 320

/*
  The constant maxrounds determines the number of rounds
  in the normalization of medium sized vectors.
  The largest dimension for a medium size normalization
  is thus qd_shemsize*maxrounds, for instance: 320*96 = 30720.
 */

#define maxrounds 96

__global__ void small_normalize_vector
 ( double *vrehihi, double *vrelohi, double *vrehilo, double *vrelolo,
   double *vimhihi, double *vimlohi, double *vimhilo, double *vimlolo,
   int dim, int dimLog2,
   double *normhihi, double *normlohi, double *normhilo, double *normlolo );
/*
 * DESCRIPTION :
 *   Kernel function to compute the norm of the complex vector,
 *   given by a vector of real parts and a vector of imaginary parts,
 *   for vectors of small dimension, in quad double precision.
 *
 * REQUIRED :
 *   The dimension dim equals the block size,
 *   which is the number of threads in the block.
 *   Moreover, all arrays have size dim.
 *
 * ON ENTRY :
 *   vrehihi     highest real parts of a complex vector;
 *   vrelohi     second highest real parts of a complex vector;
 *   vrehilo     second lowest real parts of a complex vector;
 *   vrelolo     lowest real parts of a complex vector;
 *   vimhihi     highest imaginary parts of a complex vector;
 *   vimlohi     second highest imaginary parts of a complex vector;
 *   vimhilo     second lowest imaginary parts of a complex vector;
 *   vimlolo     lowest imaginary parts of a complex vector;
 *   dim         the dimension of the vector must equal the block size;
 *   dimLog2     equals ceil(log2((double) dim), used in sum reduction.
 *
 * ON RETURN :
 *   vrehihi     highest real parts of the normalized vector;
 *   vrelohi     second highest real parts of the normalized vector;
 *   vrehilo     second lowest real parts of the normalized vector;
 *   vrelolo     lowest real parts of the normalized vector;
 *   vimhihi     highest imaginary parts of the normalized vector;
 *   vimlohi     second highest imaginary parts of the normalized vector;
 *   vimhilo     second lowest imaginary parts of the normalized vector;
 *   vimlolo     lowest imaginary parts of the normalized vector;
 *   normhihi    highest part of the 2-norm of the given vector;
 *   normlohi    second highest part of the 2-norm of the given vector;
 *   normhilo    second lowest part of the 2-norm of the given vector;
 *   normlolo    lowest part of the 2-norm of the given vector. */

__global__ void medium_normalize_vector
 ( double *vrehihi, double *vrelohi, double *vrehilo, double *vrelolo,
   double *vimhihi, double *vimlohi, double *vimhilo, double *vimlolo,
   int dim, int rnd, int rndLog2, int BS, int BSLog2,
   double *normhihi, double *normlohi, double *normhilo, double *normlolo );
/*
 * DESCRIPTION :
 *   Kernel function to compute the norm of the vector in v,
 *   for vectors of medium dimension, in quad double precision.
 *
 * REQUIRED :
 *   The dimension dim equals the block size BS times rnd.
 *   The block size BS is the number of threads in the block.
 *   All arrays are of size dim.
 *
 * ON ENTRY :
 *   vrehihi     highest real parts of a complex vector;
 *   vrelohi     second highest real parts of a complex vector;
 *   vrehilo     second lowest real parts of a complex vector;
 *   vrelolo     lowest real parts of a complex vector;
 *   vimhihi     highest imaginary parts of a complex vector;
 *   vimlohi     second highest imaginary parts of a complex vector;
 *   vimhilo     second lowest imaginary parts of a complex vector;
 *   vimlolo     lowest imaginary parts of a complex vector;
 *   dim         dimension of the vector must equal rnd*BS;
 *   rnd         the number of rounds or the multiplier for the dimension;
 *   rndLog2     the 2-logarithm of rnd for use in sum reduction;
 *   BS          the block size or the number of threads in the block;
 *   BSLog2      equals ceil(log2((double) BS), used in sum reduction.
 *
 * ON RETURN :
 *   vrehihi     highest real parts of the normalized vector;
 *   vrelohi     second highest real parts of the normalized vector;
 *   vrehilo     second lowest real parts of the normalized vector;
 *   vrelolo     lowest real parts of the normalized vector;
 *   vimhihi     highest imaginary parts of the normalized vector;
 *   vimlohi     second highest imaginary parts of the normalized vector;
 *   vimhilo     second lowest imaginary parts of the normalized vector;
 *   vimlolo     lowest imaginary parts of the normalized vector;
 *   normhihi    highest part of the 2-norm of the given vector;
 *   normlohi    second highest part of the 2-norm of the given vector;
 *   normhilo    second lowest part of the 2-norm of the given vector;
 *   normlolo    lowest part of the 2-norm of the given vector. */

__global__ void large_sum_the_squares
 ( double *vrehihi, double *vrelohi, double *vrehilo, double *vrelolo,
   double *vimhihi, double *vimlohi, double *vimhilo, double *vimlolo,
   int dim,
   double *sumshihi, double *sumslohi, double *sumshilo, double *sumslolo,
   int BS, int BSLog2 );
/*
 * DESCRIPTION :
 *   Computes the sums of the squares of the numbers in v,
 *   as needed in the 2-norm of the vector, with many blocks.
 *
 * REQUIRED :
 *   Space in sums is allocated for as many blocks in the launch.
 *
 * ON ENTRY :
 *   vrehihi     highest real parts of a complex vector;
 *   vrelohi     second highest real parts of a complex vector;
 *   vrehilo     second lowest real parts of a complex vector;
 *   vrelolo     lowest real parts of a complex vector;
 *   vimhihi     highest imaginary parts of a complex vector;
 *   vimlohi     second highest imaginary parts of a complex vector;
 *   vimhilo     second lowest imaginary parts of a complex vector;
 *   vimlolo     lowest imaginary parts of a complex vector;
 *   dim         dimension of the vector must equal BS*nbsums;
 *   BS          the block size or the number of threads in the block;
 *   BSLog2      equals ceil(log2((double) BS), used in sum reduction.
 *
 * ON RETURN :
 *   sumshihi    highest parts of sums of squares of slices of v;
 *   sumslohi    second highest parts of sums of squares of slices of v;
 *   sumshilo    second lowest parts of sums of squares of slices of v;
 *   sumslolo    lowest parts of sums of squares of slices of v. */

__global__ void large_normalize_vector
 ( double *vrehihi, double *vrelohi, double *vrehilo, double *vrelolo,
   double *vimhihi, double *vimlohi, double *vimhilo, double *vimlolo,
   int dim,
   double *sumshihi, double *sumslohi, double *sumshilo, double *sumslolo,
   int nbsums, int nbsumsLog2, int BS,
   double *normhihi, double *normlohi, double *normhilo, double *normlolo );
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
 *   vrehihi     highest real parts of a complex vector;
 *   vrelohi     second highest real parts of a complex vector;
 *   vrehilo     second lowest real parts of a complex vector;
 *   vrelolo     lowest real parts of a complex vector;
 *   vimhihi     highest imaginary parts of a complex vector;
 *   vimlohi     second highest imaginary parts of a complex vector;
 *   vimhilo     second lowest imaginary parts of a complex vector;
 *   vimlolo     lowest imaginary parts of a complex vector;
 *   dim         dimension of the vector must equal BS*nbsums;
 *   sumshihi    highest parts of sums of squares of slices of v;
 *   sumslohi    second highest parts of sums of squares of slices of v;
 *   sumshilo    second lowest parts of sums of squares of slices of v;
 *   sumslolo    lowest parts of sums of squares of slices of v;
 *   nbsums      the number of elements in sums equals
 *               the number of blocks in the kernel launch;
 *   nbsumsLog2  is ceil(log2((double) nbsums), used in sum reduction;
 *   BS          is the block size, the number of threads in a block.
 *
 * ON RETURN :
 *   vrehihi     highest real parts of the normalized vector;
 *   vrelohi     second highest real parts of the normalized vector;
 *   vrehilo     second lowest real parts of the normalized vector;
 *   vrelolo     lowest real parts of the normalized vector;
 *   vimhihi     highest imaginary parts of the normalized vector;
 *   vimlohi     second highest imaginary parts of the normalized vector;
 *   vimhilo     second lowest imaginary parts of the normalized vector;
 *   vimlolo     lowest imaginary parts of the normalized vector;
 *   normhihi    highest part of the 2-norm of the given vector;
 *   normlohi    second highest part of the 2-norm of the given vector;
 *   normhilo    second lowest part of the 2-norm of the given vector;
 *   normlolo    lowest part of the 2-norm of the given vector. */

void GPU_norm
 ( double *vrehihi_h, double *vrelohi_h, double *vrehilo_h, double *vrelolo_h,
   double *vimhihi_h, double *vimlohi_h, double *vimhilo_h, double *vimlolo_h,
   int dim, int freq, int BS,
   double *normhihi, double *normlohi, double *normhilo, double *normlolo,
   int blocked );
/*
 * DESCRIPTION :
 *   Allocates global memory on the card, transfers the data for the vector
 *   into the allocated memory, computes the 2-norm, normalizes the vector,
 *   and transfers the result from global memory of the card 
 *   to the memory of the host.
 *
 * REQUIRED :
 *   The dimension dim is a multiple of the block size BS.
 *   All arrays are of dimension dim.
 *
 * ON ENTRY :
 *   vrehihi_h   highest real parts of a complex vector;
 *   vrelohi_h   second highest real parts of a complex vector;
 *   vrehilo_h   second lowest real parts of a complex vector;
 *   vrelolo_h   lowest real parts of a complex vector;
 *   vimhihi_h   highest imaginary parts of a complex vector;
 *   vimlohi_h   second highest imaginary parts of a complex vector;
 *   vimhilo_h   second lowest imaginary parts of a complex vector;
 *   vimlolo_h   lowest imaginary parts of a complex vector;
 *   dim         dimension of the given vector;
 *   freq        frequency of the number of kernel launches (for timings);
 *   BS          block size, number of threads per block;
 *   blocked     is 0 or 1, if 0, then only one block will be launched,
 *               if 1, then as many blocks as dim/BS will be lauched and
 *               dim must be a multiple of BS.
 *
 * ON RETURN :
 *   vrehihi_h   highest real parts of the normalized vector;
 *   vrelohi_h   second highest real parts of the normalized vector;
 *   vrehilo_h   second lowest real parts of the normalized vector;
 *   vrelolo_h   lowest real parts of the normalized vector;
 *   vimhihi_h   highest imaginary parts of the normalized vector;
 *   vimlohi_h   second highest imaginary parts of the normalized vector;
 *   vimhilo_h   second lowest imaginary parts of the normalized vector;
 *   vimlolo_h   lowest imaginary parts of the normalized vector;
 *   normhihi    highest part of the 2-norm of the given vector;
 *   normlohi    second highest part of the 2-norm of the given vector;
 *   normhilo    second lowest part of the 2-norm of the given vector;
 *   normlolo    lowest part of the 2-norm of the given vector. */

#endif
