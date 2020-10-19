// Defines the prototypes and constants for the kernels
// to compute the 2-norm of a complex vector in octo double precision.

#ifndef __CMPLX8_NORM_KERNELS_H__
#define __CMPLX8_NORM_KERNELS_H__

/*
  The constant od_shmemsize determines the size of the vectors
  stored in shared memory.  As every thread works on one entry
  in the shared memory vectors, for octo double precision,
  thihis size is bounded by the number of threads in a block.
  The largest dimension for which the small normalization runs
  is thus the value of od_shmemsize.
  The largest dimension for which the large normalization runs
  is the square of the value of od_shmemsize because the algorithm
  sums two arrays of squares, the squares of the real and the
  squares of the imaginary part, thus very similar to case of
  the small normalization.  For example, for od_shmemsize = 160,
  the largest dimension is thus 160*160 = 25600.
 */

#define od_shmemsize 160

/*
  The constant maxrounds determines the number of rounds
  in the normalization of medium sized vectors.
  The largest dimension for a medium size normalization
  is thus od_shemsize*maxrounds, for instance: 160*96 = 15360.
 */

#define maxrounds 96

__global__ void small_normalize_vector
 ( double *vrehihihi, double *vrelohihi, double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo, double *vrehilolo, double *vrelololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo, double *vimhilolo, double *vimlololo,
   int dim, int dimLog2, double *normhihihi, double *normlohihi,
   double *normhilohi, double *normlolohi, double *normhihilo, 
   double *normlohilo, double *normhilolo, double *normlololo );
/*
 * DESCRIPTION :
 *   Kernel function to compute the norm of the complex vector,
 *   given by a vector of real parts and a vector of imaginary parts,
 *   for vectors of small dimension, in octo double precision.
 *
 * REQUIRED :
 *   The dimension dim equals the block size,
 *   which is the number of threads in the block.
 *   Moreover, all arrays have size dim.
 *
 * ON ENTRY :
 *   vrehihihi     highest real parts of a complex vector;
 *   vrelohihi     second highest real parts of a complex vector;
 *   vrehilohi     third highest real parts of a complex vector;
 *   vrelolohi     fourth highest real parts of a complex vector;
 *   vrehihilo     fourth lowest real parts of a complex vector;
 *   vrelohilo     third lowest real parts of a complex vector;
 *   vrehilolo     second lowest real parts of a complex vector;
 *   vrelololo     lowest real parts of a complex vector;
 *   vimhihihi     highest imaginary parts of a complex vector;
 *   vimlohihi     second highest imaginary parts of a complex vector;
 *   vimhilohi     third highest imaginary parts of a complex vector;
 *   vimlolohi     fourth highest imaginary parts of a complex vector;
 *   vimhihilo     fourth lowest imaginary parts of a complex vector;
 *   vimlohilo     third lowest imaginary parts of a complex vector;
 *   vimhilolo     second lowest imaginary parts of a complex vector;
 *   vimlololo     lowest imaginary parts of a complex vector;
 *   dim           the dimension of the vector must equal the block size;
 *   dimLog2       equals ceil(log2((double) dim), used in sum reduction.
 *
 * ON RETURN :
 *   vrehihihi     highest real parts of the normalized vector;
 *   vrelohihi     second highest real parts of the normalized vector;
 *   vrehilohi     third highest real parts of the normalized vector;
 *   vrelolohi     fourth highest real parts of the normalized vector;
 *   vrehihilo     fourth lowest real parts of the normalized vector;
 *   vrelohilo     third highest real parts of the normalized vector;
 *   vrehilolo     second lowest real parts of the normalized vector;
 *   vrelololo     lowest real parts of the normalized vector;
 *   vimhihihi     highest imaginary parts of the normalized vector;
 *   vimlohihi     second highest imaginary parts of the normalized vector;
 *   vimhilohi     third highest imaginary parts of the normalized vector;
 *   vimlolohi     fourth highest imaginary parts of the normalized vector;
 *   vimhihilo     fourth lowest imaginary parts of the normalized vector;
 *   vimlohilo     third lowest imaginary parts of the normalized vector;
 *   vimhilolo     second lowest imaginary parts of the normalized vector;
 *   vimlololo     lowest imaginary parts of the normalized vector;
 *   normhihihi    highest part of the 2-norm of the given vector;
 *   normlohihi    second highest part of the 2-norm of the given vector;
 *   normhilohi    third highest part of the 2-norm of the given vector;
 *   normlolohi    fourth highest part of the 2-norm of the given vector;
 *   normhihilo    fourth lowest part of the 2-norm of the given vector;
 *   normlohilo    third lowest part of the 2-norm of the given vector;
 *   normhilolo    second lowest part of the 2-norm of the given vector;
 *   normlololo    lowest part of the 2-norm of the given vector. */

__global__ void medium_normalize_vector
 ( double *vrehihihi, double *vrelohihi, double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo, double *vrehilolo, double *vrelololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo, double *vimhilolo, double *vimlololo,
   int dim, int rnd, int rndLog2, int BS, int BSLog2,
   double *normhihihi, double *normlohihi, double *normhilohi,
   double *normlolohi, double *normhihilo, double *normlohilo,
   double *normhilolo, double *normlololo );
/*
 * DESCRIPTION :
 *   Kernel function to compute the norm of the vector in v,
 *   for vectors of medium dimension, in octo double precision.
 *
 * REQUIRED :
 *   The dimension dim equals the block size BS times rnd.
 *   The block size BS is the number of threads in the block.
 *   All arrays are of size dim.
 *
 * ON ENTRY :
 *   vrehihihi     highest real parts of a complex vector;
 *   vrelohihi     second highest real parts of a complex vector;
 *   vrehilohi     third highest real parts of a complex vector;
 *   vrelolohi     fourth highest real parts of a complex vector;
 *   vrehihilo     fourth lowest real parts of a complex vector;
 *   vrelohilo     third lowest real parts of a complex vector;
 *   vrehilolo     second lowest real parts of a complex vector;
 *   vrelololo     lowest real parts of a complex vector;
 *   vimhihihi     highest imaginary parts of a complex vector;
 *   vimlohihi     second highest imaginary parts of a complex vector;
 *   vimhilohi     third highest imaginary parts of a complex vector;
 *   vimlolohi     fourth highest imaginary parts of a complex vector;
 *   vimhihilo     fourth lowest imaginary parts of a complex vector;
 *   vimlohilo     third lowest imaginary parts of a complex vector;
 *   vimhilolo     second lowest imaginary parts of a complex vector;
 *   vimlololo     lowest imaginary parts of a complex vector;
 *   dim         dimension of the vector must equal rnd*BS;
 *   rnd         the number of rounds or the multiplier for the dimension;
 *   rndLog2     the 2-logarithm of rnd for use in sum reduction;
 *   BS          the block size or the number of threads in the block;
 *   BSLog2      equals ceil(log2((double) BS), used in sum reduction.
 *
 * ON RETURN :
 *   vrehihihi     highest real parts of the normalized vector;
 *   vrelohihi     second highest real parts of the normalized vector;
 *   vrehilohi     third highest real parts of the normalized vector;
 *   vrelolohi     fourth highest real parts of the normalized vector;
 *   vrehihilo     fourth lowest real parts of the normalized vector;
 *   vrelohilo     third highest real parts of the normalized vector;
 *   vrehilolo     second lowest real parts of the normalized vector;
 *   vrelololo     lowest real parts of the normalized vector;
 *   vimhihihi     highest imaginary parts of the normalized vector;
 *   vimlohihi     second highest imaginary parts of the normalized vector;
 *   vimhilohi     third highest imaginary parts of the normalized vector;
 *   vimlolohi     fourth highest imaginary parts of the normalized vector;
 *   vimhihilo     fourth lowest imaginary parts of the normalized vector;
 *   vimlohilo     third lowest imaginary parts of the normalized vector;
 *   vimhilolo     second lowest imaginary parts of the normalized vector;
 *   vimlololo     lowest imaginary parts of the normalized vector;
 *   normhihihi    highest part of the 2-norm of the given vector;
 *   normlohihi    second highest part of the 2-norm of the given vector;
 *   normhilohi    third highest part of the 2-norm of the given vector;
 *   normlolohi    fourth highest part of the 2-norm of the given vector;
 *   normhihilo    fourth lowest part of the 2-norm of the given vector;
 *   normlohilo    third lowest part of the 2-norm of the given vector;
 *   normhilolo    second lowest part of the 2-norm of the given vector;
 *   normlololo    lowest part of the 2-norm of the given vector. */

__global__ void large_sum_the_squares
 ( double *vrehihihi, double *vrelohihi, double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo, double *vrehilolo, double *vrelololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo, double *vimhilolo, double *vimlololo,
   int dim,
   double *sumshihihi, double *sumslohihi, double *sumshilohi,
   double *sumslolohi, double *sumshihilo, double *sumslohilo,
   double *sumshilolo, double *sumslololo, int BS, int BSLog2 );
/*
 * DESCRIPTION :
 *   Computes the sums of the squares of the numbers in v,
 *   as needed in the 2-norm of the vector, with many blocks.
 *
 * REQUIRED :
 *   Space in sums is allocated for as many blocks in the launch.
 *
 * ON ENTRY :
 *   vrehihihi     highest real parts of a complex vector;
 *   vrelohihi     second highest real parts of a complex vector;
 *   vrehilohi     third highest real parts of a complex vector;
 *   vrelolohi     fourth highest real parts of a complex vector;
 *   vrehihilo     fourth lowest real parts of a complex vector;
 *   vrelohilo     third lowest real parts of a complex vector;
 *   vrehilolo     second lowest real parts of a complex vector;
 *   vrelololo     lowest real parts of a complex vector;
 *   vimhihihi     highest imaginary parts of a complex vector;
 *   vimlohihi     second highest imaginary parts of a complex vector;
 *   vimhilohi     third highest imaginary parts of a complex vector;
 *   vimlolohi     fourth highest imaginary parts of a complex vector;
 *   vimhihilo     fourth lowest imaginary parts of a complex vector;
 *   vimlohilo     third lowest imaginary parts of a complex vector;
 *   vimhilolo     second lowest imaginary parts of a complex vector;
 *   vimlololo     lowest imaginary parts of a complex vector;
 *   dim           dimension of the vector must equal BS*nbsums;
 *   BS            the block size or the number of threads in the block;
 *   BSLog2        equals ceil(log2((double) BS), used in sum reduction.
 *
 * ON RETURN :
 *   sumshihihi    highest parts of sums of squares of slices of v;
 *   sumslohihi    second highest parts of sums of squares of slices of v;
 *   sumshilohi    third highest parts of sums of squares of slices of v;
 *   sumslolohi    fourth highest parts of sums of squares of slices of v;
 *   sumshihilo    fourth lowest parts of sums of squares of slices of v;
 *   sumslohilo    third lowest parts of sums of squares of slices of v;
 *   sumshilolo    second lowest parts of sums of squares of slices of v;
 *   sumslololo    lowest parts of sums of squares of slices of v. */

__global__ void large_normalize_vector
 ( double *vrehihihi, double *vrelohihi, double *vrehilohi, double *vrelolohi,
   double *vrehihilo, double *vrelohilo, double *vrehilolo, double *vrelololo,
   double *vimhihihi, double *vimlohihi, double *vimhilohi, double *vimlolohi,
   double *vimhihilo, double *vimlohilo, double *vimhilolo, double *vimlololo,
   int dim,
   double *sumshihihi, double *sumslohihi, double *sumshilohi,
   double *sumslolohi, double *sumshihilo, double *sumslohilo,
   double *sumshilolo, double *sumslololo,
   int nbsums, int nbsumsLog2, int BS, double *normhihihi, double *normlohihi,
   double *normhilohi, double *normlolohi, double *normhihilo,
   double *normlohilo, double *normhilolo, double *normlololo );
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
 *   vrehihihi     highest real parts of a complex vector;
 *   vrelohihi     second highest real parts of a complex vector;
 *   vrehilohi     third highest real parts of a complex vector;
 *   vrelolohi     fourth highest real parts of a complex vector;
 *   vrehihilo     fourth lowest real parts of a complex vector;
 *   vrelohilo     third lowest real parts of a complex vector;
 *   vrehilolo     second lowest real parts of a complex vector;
 *   vrelololo     lowest real parts of a complex vector;
 *   vimhihihi     highest imaginary parts of a complex vector;
 *   vimlohihi     second highest imaginary parts of a complex vector;
 *   vimhilohi     third highest imaginary parts of a complex vector;
 *   vimlolohi     fourth highest imaginary parts of a complex vector;
 *   vimhihilo     fourth lowest imaginary parts of a complex vector;
 *   vimlohilo     third lowest imaginary parts of a complex vector;
 *   vimhilolo     second lowest imaginary parts of a complex vector;
 *   vimlololo     lowest imaginary parts of a complex vector;
 *   dim           dimension of the vector must equal BS*nbsums;
 *   sumshihihi    highest parts of sums of squares of slices of v;
 *   sumslohihi    second highest parts of sums of squares of slices of v;
 *   sumshilohi    third highest parts of sums of squares of slices of v;
 *   sumslolohi    fourth highest parts of sums of squares of slices of v;
 *   sumshihilo    fourth lowest parts of sums of squares of slices of v;
 *   sumslohilo    third lowest parts of sums of squares of slices of v;
 *   sumshilolo    second lowest parts of sums of squares of slices of v;
 *   sumslololo    lowest parts of sums of squares of slices of v;
 *   nbsums        the number of elements in sums equals
 *                 the number of blocks in the kernel launch;
 *   nbsumsLog2    is ceil(log2((double) nbsums), used in sum reduction;
 *   BS            is the block size, the number of threads in a block.
 *
 * ON RETURN :
 *   vrehihihi     highest real parts of the normalized vector;
 *   vrelohihi     second highest real parts of the normalized vector;
 *   vrehilohi     third highest real parts of the normalized vector;
 *   vrelolohi     fourth highest real parts of the normalized vector;
 *   vrehihilo     fourth lowest real parts of the normalized vector;
 *   vrelohilo     third highest real parts of the normalized vector;
 *   vrehilolo     second lowest real parts of the normalized vector;
 *   vrelololo     lowest real parts of the normalized vector;
 *   vimhihihi     highest imaginary parts of the normalized vector;
 *   vimlohihi     second highest imaginary parts of the normalized vector;
 *   vimhilohi     third highest imaginary parts of the normalized vector;
 *   vimlolohi     fourth highest imaginary parts of the normalized vector;
 *   vimhihilo     fourth lowest imaginary parts of the normalized vector;
 *   vimlohilo     third lowest imaginary parts of the normalized vector;
 *   vimhilolo     second lowest imaginary parts of the normalized vector;
 *   vimlololo     lowest imaginary parts of the normalized vector;
 *   normhihihi    highest part of the 2-norm of the given vector;
 *   normlohihi    second highest part of the 2-norm of the given vector;
 *   normhilohi    third highest part of the 2-norm of the given vector;
 *   normlolohi    fourth highest part of the 2-norm of the given vector;
 *   normhihilo    fourth lowest part of the 2-norm of the given vector;
 *   normlohilo    third lowest part of the 2-norm of the given vector;
 *   normhilolo    second lowest part of the 2-norm of the given vector;
 *   normlololo    lowest part of the 2-norm of the given vector. */

void GPU_norm
 ( double *vrehihihi_h, double *vrelohihi_h, double *vrehilohi_h,
   double *vrelolohi_h, double *vrehihilo_h, double *vrelohilo_h,
   double *vrehilolo_h, double *vrelololo_h,
   double *vimhihihi_h, double *vimlohihi_h, double *vimhilohi_h,
   double *vimlolohi_h, double *vimhihilo_h, double *vimlohilo_h,
   double *vimhilolo_h, double *vimlololo_h,
   int dim, int freq, int BS, double *normhihihi, double *normlohihi,
   double *normhilohi, double *normlolohi, double *normhihilo,
   double *normlohilo, double *normhilolo, double *normlololo, int blocked );
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
 *   vrehihihi_h    highest real parts of a complex vector;
 *   vrelohihi_h    second highest real parts of a complex vector;
 *   vrehilohi_h    third highest real parts of a complex vector;
 *   vrelolohi_h    fourth highest real parts of a complex vector;
 *   vrehihilo_h    fourth lowest real parts of a complex vector;
 *   vrelohilo_h    third lowest real parts of a complex vector;
 *   vrehilolo_h    second lowest real parts of a complex vector;
 *   vrelololo_h    lowest real parts of a complex vector;
 *   vimhihihi_h    highest imaginary parts of a complex vector;
 *   vimlohihi_h    second highest imaginary parts of a complex vector;
 *   vimhilohi_h    third highest imaginary parts of a complex vector;
 *   vimlolohi_h    fourth highest imaginary parts of a complex vector;
 *   vimhihilo_h    fourth lowest imaginary parts of a complex vector;
 *   vimlohilo_h    third lowest imaginary parts of a complex vector;
 *   vimhilolo_h    second lowest imaginary parts of a complex vector;
 *   vimlololo_h    lowest imaginary parts of a complex vector;
 *   dim            dimension of the given vector;
 *   freq           frequency of the number of kernel launches
 *                  (for timings);
 *   BS             block size, number of threads per block;
 *   blocked        is 0 or 1, if 0, then only one block will be launched,
 *                  if 1, then as many blocks as dim/BS will be lauched
 *                  and dim must be a multiple of BS.
 *
 * ON RETURN :
 *   vrehihihi_h    highest real parts of the normalized vector;
 *   vrelohihi_h    second highest real parts of the normalized vector;
 *   vrehilohi_h    third highest real parts of the normalized vector;
 *   vrelolohi_h    fourth highest real parts of the normalized vector;
 *   vrehihilo_h    fourth lowest real parts of the normalized vector;
 *   vrelohilo_h    third highest real parts of the normalized vector;
 *   vrehilolo_h    second lowest real parts of the normalized vector;
 *   vrelololo_h    lowest real parts of the normalized vector;
 *   vimhihihi_h    highest imaginary parts of the normalized vector;
 *   vimlohihi_h    second highest imaginary parts of the normalized vector;
 *   vimhilohi_h    third highest imaginary parts of the normalized vector;
 *   vimlolohi_h    fourth highest imaginary parts of the normalized vector;
 *   vimhihilo_h    fourth lowest imaginary parts of the normalized vector;
 *   vimlohilo_h    third lowest imaginary parts of the normalized vector;
 *   vimhilolo_h    second lowest imaginary parts of the normalized vector;
 *   vimlololo_h    lowest imaginary parts of the normalized vector;
 *   normhihihi     highest part of the 2-norm of the given vector;
 *   normlohihi     second highest part of the 2-norm of the given vector;
 *   normhilohi     third highest part of the 2-norm of the given vector;
 *   normlolohi     fourth highest part of the 2-norm of the given vector;
 *   normhihilo     fourth lowest part of the 2-norm of the given vector;
 *   normlohilo     third lowest part of the 2-norm of the given vector;
 *   normhilolo     second lowest part of the 2-norm of the given vector;
 *   normlololo     lowest part of the 2-norm of the given vector. */

#endif
