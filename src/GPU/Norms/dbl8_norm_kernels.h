// Defines the constant and the prototypes for the kernels 
// to compute the 2-norm of a vector in octo double precision.

#ifndef __DBL8_NORM_KERNELS_H__
#define __DBL8_NORM_KERNELS_H__

/*
  The constant od_shmemsize determines the size of the vectors
  stored in shared memory.  As every thread works on one entry
  in the shared memory vectors, for octo double precision,
  this size is bounded by the number of threads in a block.
  The largest dimension for which the small normalization runs
  is thus the value of od_shmemsize.
  The largest dimension for which the large normalization runs
  is the square of the value of od_shmemsize because the algorithm
  sums two arrays of squares, the squares of the real and the
  squares of the imaginary part, thus very similar to case of
  the small normalization.  For example, for od_shmemsize = 320,
  the largest dimension is thus 320*320 = 102400.
 */

#define od_shmemsize 320

/*
  The constant maxrounds determines the number of rounds
  in the normalization of medium sized vectors.
  The largest dimension for a medium size normalization
  is thus od_shemsize*maxrounds, for instance: 320*96 = 30720.
 */

#define maxrounds 96

__global__ void small_normalize_vector
 ( double *vhihihi, double *vlohihi, double *vhilohi, double *vlolohi,
   double *vhihilo, double *vlohilo, double *vhilolo, double *vlololo,
   int dim, int dimLog2, double *normhihihi, double *normlohihi,
   double *normhilohi, double *normlolohi, double *normhihilo,
   double *normlohilo, double *normhilolo, double *normlololo );
/*
 * DESCRIPTION :
 *   Kernel function to compute the norm of the vector,
 *   given with eight arrays for all parts of an octo double,
 *   for vectors of small dimension, in octo double precision.
 *
 * REQUIRED :
 *   The dimension dim equals the block size,
 *   which is the number of threads in the block.
 *
 * ON ENTRY :
 *   vhihihi   highest parts of an octo double vector of dimension dim;
 *   vlohihi   2nd highest parts of an octo double vector of dimension dim;
 *   vhilohi   3rd highest parts of an octo double vector of dimension dim;
 *   vlolohi   4th highest parts of an octo double vector of dimension dim;
 *   vhihilo   4th lowest parts of an octo double vector of dimension dim;
 *   vlohilo   3rd lowest parts of an octo double vector of dimension dim;
 *   vhilolo   2nd lowest parts of an octo double vector of dimension dim;
 *   vlololo   lowest parts of an octo double vector of dimension dim;
 *   dim       the dimension of the vector must equal the block size;
 *   dimLog2   equals ceil(log2((double) dim), used in sum reduction.
 *
 * ON RETURN :
 *   vhihihi   highest parts of the octo double vector in the same direction 
 *             of the original vector but with norm one;
 *   vlohihi   second highest parts of the octo double vector in the same
 *             direction of the original vector but with norm one;
 *   vhilohi   third highest parts of the octo double vector in the same
 *             direction of the original vector but with norm one;
 *   vlolohi   fourth highest parts of the octo double vector in the same
 *             direction of the original vector but with norm one;
 *   vhihilo   fourth lowest parts of the octo double vector in the same
 *             direction of the original vector but with norm one;
 *   vhilolo   third lowest parts of the octo double vector in the same
 *             direction of the original vector but with norm one;
 *   vlohilo   second lowest parts of the octo double vector in the same
 *             direction of the original vector but with norm one;
 *   vlololo   lowest parts of the octo double vector in the same direction 
 *             of the original vector but with norm one;
 *   normhihihi is the highest part of the 2-norm of the original vector;
 *   normlohihi is the 2nd highest part of the 2-norm of the original vector;
 *   normhilohi is the 3rd highest part of the 2-norm of the original vector;
 *   normlolohi is the 4th highest part of the 2-norm of the original vector;
 *   normhihilo is the 4th lowest part of the 2-norm of the original vector;
 *   normlohilo is the 3rd lowest part of the 2-norm of the original vector;
 *   normhilolo is the 2nd lowest part of the 2-norm of the original vector;
 *   normlololo is the lowest part of the 2-norm of the original vector. */

__global__ void medium_normalize_vector
 ( double *vhihihi, double *vlohihi, double *vhilohi, double *vlolohi,
   double *vhihilo, double *vlohilo, double *vhilolo, double *vlololo,
   int dim, int rnd, int rndLog2, int BS, int BSLog2,
   double *normhihihi, double *normlohihi, double *normhilohi,
   double *normlolohi, double *normhihilo, double *normlohilo,
   double *normhilolo, double *normlololo );
/*
 * DESCRIPTION :
 *   Kernel function to compute the norm of the vector,
 *   given by eight arrays for all parts of the octo double numbers,
 *   for vectors of medium dimension, in octo double precision.
 *
 * REQUIRED :
 *   The dimension dim equals the block size BS times rnd.
 *   The block size BS is the number of threads in the block.
 *
 * ON ENTRY :
 *   vhihihi   highest parts of an octo double vector of dimension dim;
 *   vlohihi   2nd highest parts of an octo double vector of dimension dim;
 *   vhilohi   3rd highest parts of an octo double vector of dimension dim;
 *   vlolohi   4th highest parts of an octo double vector of dimension dim;
 *   vhihilo   4th lowest parts of an octo double vector of dimension dim;
 *   vlohilo   3rd lowest parts of an octo double vector of dimension dim;
 *   vhilolo   2nd lowest parts of an octo double vector of dimension dim;
 *   vlololo   lowest parts of an octo double vector of dimension dim;
 *   dim       dimension of the vector must equal rnd*BS;
 *   rnd       the number of rounds or the multiplier for the dimension;
 *   rndLog2   the 2-logarithm of rnd for use in sum reduction;
 *   BS        the block size or the number of threads in the block;
 *   BSLog2    equals ceil(log2((double) BS), used in sum reduction.
 *
 * ON RETURN :
 *   vhihihi   highest parts of the octo double vector in the same direction 
 *             of the original vector but with norm one;
 *   vlohihi   second highest parts of the octo double vector in the same
 *             direction of the original vector but with norm one;
 *   vhilohi   third highest parts of the octo double vector in the same
 *             direction of the original vector but with norm one;
 *   vlolohi   fourth highest parts of the octo double vector in the same
 *             direction of the original vector but with norm one;
 *   vhihilo   fourth lowest parts of the octo double vector in the same
 *             direction of the original vector but with norm one;
 *   vhilolo   third lowest parts of the octo double vector in the same
 *             direction of the original vector but with norm one;
 *   vlohilo   second lowest parts of the octo double vector in the same
 *             direction of the original vector but with norm one;
 *   vlololo   lowest parts of the octo double vector in the same direction 
 *             of the original vector but with norm one;
 *   normhihihi is the highest part of the 2-norm of the original vector;
 *   normlohihi is the 2nd highest part of the 2-norm of the original vector;
 *   normhilohi is the 3rd highest part of the 2-norm of the original vector;
 *   normlolohi is the 4th highest part of the 2-norm of the original vector;
 *   normhihilo is the 4th lowest part of the 2-norm of the original vector;
 *   normlohilo is the 3rd lowest part of the 2-norm of the original vector;
 *   normhilolo is the 2nd lowest part of the 2-norm of the original vector;
 *   normlololo is the lowest part of the 2-norm of the original vector. */

__global__ void large_sum_the_squares
 ( double *vhihihi, double *vlohihi, double *vhilohi, double *vlolohi,
   double *vhihilo, double *vlohilo, double *vhilolo, double *vlololo,
   int dim, double *sumshihihi, double *sumslohihi, double *sumshilohi,
   double *sumslolohi, double *sumshihilo, double *sumslohilo,
   double *sumshilolo, double *sumslololo, int BS, int BSLog2 );
/*
 * DESCRIPTION :
 *   Computes the sums of the squares of the numbers in a vector,
 *   given by eight arrays for all parts of the octo double numbers,
 *   as needed in the 2-norm of the vector, with many blocks.
 *
 * REQUIRED :
 *   Space in sums is allocated for as many blocks in the launch.
 *
 * ON ENTRY :
 *   vhihihi   highest parts of an octo double vector of dimension dim;
 *   vlohihi   2nd highest parts of an octo double vector of dimension dim;
 *   vhilohi   3rd highest parts of an octo double vector of dimension dim;
 *   vlolohi   4th highest parts of an octo double vector of dimension dim;
 *   vhihilo   4th lowest parts of an octo double vector of dimension dim;
 *   vlohilo   3rd lowest parts of an octo double vector of dimension dim;
 *   vhilolo   2nd lowest parts of an octo double vector of dimension dim;
 *   vlololo   lowest parts of an octo double vector of dimension dim;
 *   dim       dimension of the vector must equal BS*nbsums;
 *   BS        the block size or the number of threads in the block;
 *   BSLog2    equals ceil(log2((double) BS), used in sum reduction.
 *
 * ON RETURN :
 *   sumshihihi is the highest part of computed sums of squares of slices;
 *   sumslohihi is the 2nd highest part of computed sums of squares of slices;
 *   sumshilohi is the 3rd highest part of computed sums of squares of slices;
 *   sumslolohi is the 4th highest part of computed sums of squares of slices;
 *   sumshihilo is the 4th lowest part of computed sums of squares of slices;
 *   sumslohilo is the 3rd lowest part of computed sums of squares of slices;
 *   sumshilolo is the 2nd lowest part of computed sums of squares of slices;
 *   sumslololo is the lowest part of computed sums of squares of slices. */

__global__ void large_normalize_vector
 ( double *vhihihi, double *vlohihi, double *vhilohi, double *vlolohi,
   double *vhihilo, double *vlohilo, double *vhilolo, double *vlololo,
   int dim, double *sumshihihi, double *sumslohihi, double *sumshilohi,
   double *sumslolohi, double *sumshihilo, double *sumslohilo,
   double *sumshilolo, double *sumslololo, int nbsums, int nbsumsLog2,
   int BS, double *normhihihi, double *normlohihi, double *normhilohi,
   double *normlolohi, double *normhihilo, double *normlohilo,
   double *normhilolo, double *normlololo );
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
 *   vhihihi   highest parts of an octo double vector of dimension dim;
 *   vlohihi   2nd highest parts of an octo double vector of dimension dim;
 *   vhilohi   3rd highest parts of an octo double vector of dimension dim;
 *   vlolohi   4th highest parts of an octo double vector of dimension dim;
 *   vhihilo   4th lowest parts of an octo double vector of dimension dim;
 *   vlohilo   3rd lowest parts of an octo double vector of dimension dim;
 *   vhilolo   2nd lowest parts of an octo double vector of dimension dim;
 *   vlololo   lowest parts of an octo double vector of dimension dim;
 *   dim       dimension of the vector must equal BS*nbsums;
 *   sumshihihi is the highest part of computed sums of squares of slices;
 *   sumslohihi is the 2nd highest part of computed sums of squares of slices;
 *   sumshilohi is the 3rd highest part of computed sums of squares of slices;
 *   sumslolohi is the 4th highest part of computed sums of squares of slices;
 *   sumshihilo is the 4th lowest part of computed sums of squares of slices;
 *   sumslohilo is the 3rd lowest part of computed sums of squares of slices;
 *   sumshilolo is the 2nd lowest part of computed sums of squares of slices;
 *   sumslololo is the lowest part of computed sums of squares of slices;
 *   nbsums    the number of elements in sums equals
 *             the number of blocks in the kernel launch;
 *   nbsumsLog2 is ceil(log2((double) nbsums), used in sum reduction;
 *   BS        is the block size, the number of threads in a block.
 *
 * ON RETURN :
 *   vhihihi   highest parts of the octo double vector in the same direction 
 *             of the original vector but with norm one;
 *   vlohihi   second highest parts of the octo double vector in the same
 *             direction of the original vector but with norm one;
 *   vhilohi   third highest parts of the octo double vector in the same
 *             direction of the original vector but with norm one;
 *   vlolohi   fourth highest parts of the octo double vector in the same
 *             direction of the original vector but with norm one;
 *   vhihilo   fourth lowest parts of the octo double vector in the same
 *             direction of the original vector but with norm one;
 *   vhilolo   third lowest parts of the octo double vector in the same
 *             direction of the original vector but with norm one;
 *   vlohilo   second lowest parts of the octo double vector in the same
 *             direction of the original vector but with norm one;
 *   vlololo   lowest parts of the octo double vector in the same direction 
 *             of the original vector but with norm one;
 *   normhihihi is the highest part of the 2-norm of the original vector;
 *   normlohihi is the 2nd highest part of the 2-norm of the original vector;
 *   normhilohi is the 3rd highest part of the 2-norm of the original vector;
 *   normlolohi is the 4th highest part of the 2-norm of the original vector;
 *   normhihilo is the 4th lowest part of the 2-norm of the original vector;
 *   normlohilo is the 3rd lowest part of the 2-norm of the original vector;
 *   normhilolo is the 2nd lowest part of the 2-norm of the original vector;
 *   normlololo is the lowest part of the 2-norm of the original vector. */

void GPU_norm
 ( double *vhihihi_h, double *vlohihi_h, double *vhilohi_h, double *vlolohi_h,
   double *vhihilo_h, double *vlohilo_h, double *vhilolo_h, double *vlololo_h,
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
 *
 * ON ENTRY :
 *   vhihihi_h   highest parts of an octo double vector of dimension dim;
 *   vlohihi_h   2nd highest parts of an octo double vector of dimension dim;
 *   vhilohi_h   3rd highest parts of an octo double vector of dimension dim;
 *   vlolohi_h   4th highest parts of an octo double vector of dimension dim;
 *   vhihilo_h   4th lowest parts of an octo double vector of dimension dim;
 *   vlohilo_h   3rd lowest parts of an octo double vector of dimension dim;
 *   vhilolo_h   2nd lowest parts of an octo double vector of dimension dim;
 *   vlololo_h   lowest parts of an octo double vector of dimension dim;
 *   dim         dimension of the vector v_h;
 *   freq        frequency of the number of kernel launches (for timings);
 *   BS          block size, number of threads per block;
 *   blocked     is 0 or 1, if 0, then only one block will be launched,
 *               if 1, then as many blocks as dim/BS will be lauched and
 *               dim must be a multiple of BS.
 *
 * ON RETURN :
 *   vhihihi_h   highest parts of the normalized vector;
 *   vlohihi_h   second highest parts of the normalized vector;
 *   vhilohi_h   third highest parts of the normalized vector;
 *   vlolohi_h   fourth highest parts of the normalized vector;
 *   vhihilo_h   fourth lowest parts of the normalized vector;
 *   vlohilo_h   third lowest parts of the normalized vector;
 *   vhilolo_h   second lowest parts of the normalized vector;
 *   vlololo_h   lowest parts of the normalized vector;
 *   normhihihi  highest part of the 2-norm of the original vector;
 *   normlohihi  second highest part of the 2-norm of the original vector;
 *   normhilohi  third highest part of the 2-norm of the original vector;
 *   normlolohi  fourth highest part of the 2-norm of the original vector;
 *   normhihilo  fourth lowest part of the 2-norm of the original vector;
 *   normlohilo  third lowest part of the 2-norm of the original vector;
 *   normhilolo  second lowest part of the 2-norm of the original vector;
 *   normlololo  lowest part of the 2-norm of the original vector. */

#endif
