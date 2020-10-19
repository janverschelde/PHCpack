// Defines the constants and the prototypes for the kernels 
// to compute the 2-norm of a vector in quad double precision.

#ifndef __DBL4_NORM_KERNELS_H__
#define __DBL4_NORM_KERNELS_H__

/*
  The constant qd_shmemsize determines the size of the vectors
  stored in shared memory.  As every thread works on one entry
  in the shared memory vectors, for quad double precision,
  this size is bounded by the number of threads in a block.
  The largest dimension for which the small normalization runs
  is thus the value of qd_shmemsize.
  The largest dimension for which the large normalization runs
  is the square of the value of qd_shmemsize because the algorithm
  sums two arrays of squares, the squares of the real and the
  squares of the imaginary part, thus very similar to case of
  the small normalization.  For example, for od_shmemsize = 640
  the largest dimension is thus 640*640 = 409600.
 */

#define qd_shmemsize 640

/*
  The constant maxrounds determines the number of rounds
  in the normalization of medium sized vectors.
  The largest dimension for a medium size normalization
  is thus qd_shemsize*maxrounds, for instance: 640*192 = 122880.
 */

#define maxrounds 192

__global__ void small_normalize_vector
 ( double *vhihi, double *vlohi, double *vhilo, double *vlolo,
   int dim, int dimLog2,
   double *normhihi, double *normlohi, double *normhilo, double *normlolo );
/*
 * DESCRIPTION :
 *   Kernel function to compute the norm of the vector,
 *   given with highest parts in vhihi, second highest parts in vlohi,
 *   second lowest parts in vhilo, and lowest parts in vlolo,
 *   for vectors of small dimension, in quad double precision.
 *
 * REQUIRED :
 *   The dimension dim equals the block size,
 *   which is the number of threads in the block.
 *
 * ON ENTRY :
 *   vhihi     highest parts of a quad double vector of dimension dim;
 *   vlohi     2nd highest parts of a quad double vector of dimension dim;
 *   vhilo     2nd lowest parts of a quad double vector of dimension dim;
 *   vlolo     lowest parts of a quad double vector of dimension dim;
 *   dim       the dimension of the vector must equal the block size;
 *   dimLog2   equals ceil(log2((double) dim), used in sum reduction.
 *
 * ON RETURN :
 *   vhi       high parts of the quad double vector in the same direction 
 *             of the original vector but with norm one;
 *   vmi       middle parts of the quad double vector in the same direction 
 *             of the original vector but with norm one;
 *   vlo       low parts of the quad double vector in the same direction 
 *             of the original vector but with norm one;
 *   normhi    high part of the 2-norm of the original vector;
 *   normmi    middle part of the 2-norm of the original vector;
 *   normlo    low part of the 2-norm of the original vector. */

__global__ void medium_normalize_vector
 ( double *vhihi, double *vlohi, double *vhilo, double *vlolo,
   int dim, int rnd, int rndLog2, int BS, int BSLog2,
   double *normhihi, double *normlohi, double *normhilo, double *normlolo );
/*
 * DESCRIPTION :
 *   Kernel function to compute the norm of the vector,
 *   given with highest parts in vhihi, second highest parts in vlohi,
 *   second lowest parts in vhilo, and lowest parts in vlolo,
 *   for vectors of medium dimension, in quad double precision.
 *
 * REQUIRED :
 *   The dimension dim equals the block size BS times rnd.
 *   The block size BS is the number of threads in the block.
 *
 * ON ENTRY :
 *   vhihi     highest parts of a quad double vector of dimension dim;
 *   vlohi     2nd highest parts of a quad double vector of dimension dim;
 *   vhilo     2nd lowest parts of a quad double vector of dimension dim;
 *   vlolo     lowest parts of a quad double vector of dimension dim;
 *   dim       dimension of the vector must equal rnd*BS;
 *   rnd       the number of rounds or the multiplier for the dimension;
 *   rndLog2   the 2-logarithm of rnd for use in sum reduction;
 *   BS        the block size or the number of threads in the block;
 *   BSLog2    equals ceil(log2((double) BS), used in sum reduction.
 *
 * ON RETURN :
 *   vhihi     highest parts of the quad double vector in the same direction 
 *             of the original vector but with norm one;
 *   vlohi     second highest parts of the quad double vector in the same
 *             direction of the original vector but with norm one;
 *   vhilo     second lowest parts of the quad double vector in the same
 *             direction of the original vector but with norm one;
 *   vlolo     lowest parts of the quad double vector in the same direction 
 *             of the original vector but with norm one;
 *   normhihi  highest part of the 2-norm of the original vector;
 *   normlomi  second highest part of the 2-norm of the original vector;
 *   normhilo  second lowest part of the 2-norm of the original vector;
 *   normlolo  lowest part of the 2-norm of the original vector. */

__global__ void large_sum_the_squares
 ( double *vhihi, double *vlohi, double *vhilo, double *vlolo, int dim,
   double *sumshihi, double *sumslohi, double *sumshilo, double *sumslolo,
   int BS, int BSLog2 );
/*
 * DESCRIPTION :
 *   Computes the sums of the squares of the numbers in a vector,
 *   given with highest parts in vhihi, second highest parts in vlohi,
 *   second lowest parts in vhilo, and lowes parts in vlolo,
 *   as needed in the 2-norm of the vector, with many blocks.
 *
 * REQUIRED :
 *   Space in sums is allocated for as many blocks in the launch.
 *
 * ON ENTRY :
 *   vhihi     highest parts of a quad double vector of dimension dim;
 *   vlohi     second highest parts of a quad double vector of dimension dim;
 *   vhilo     second lowest parts of a quad double vector of dimension dim;
 *   vlolo     lowest parts of a quad double vector of dimension dim;
 *   dim       dimension of the vector must equal BS*nbsums;
 *   BS        the block size or the number of threads in the block;
 *   BSLog2    equals ceil(log2((double) BS), used in sum reduction.
 *
 * ON RETURN :
 *   sumshihi  highest parts of computed sums of squares of vector slices;
 *   sumslohi  2nd highest parts of computed sums of squares of vector slices;
 *   sumshilo  2nd lowest parts of computed sums of squares of vector slices;
 *   sumslolo  lowest parts of computed sums of squares of vector slices. */

__global__ void large_normalize_vector
 ( double *vhihi, double *vlohi, double *vhilo, double *vlolo, int dim,
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
 *   vhihi     highest parts of a quad double vector of dimension dim;
 *   vlohi     2nd highest parts of a quad double vector of dimension dim;
 *   vhilo     2nd lowest parts of a quad double vector of dimension dim;
 *   vlolo     lowest parts of a quad double vector of dimension dim;
 *   dim       dimension of the vector must equal BS*nbsums;
 *   sumshihi  highest parts of computed sums of squares of vector slices;
 *   sumslohi  2nd highest parts of computed sums of squares of vector slices;
 *   sumshilo  lowest parts of computed sums of squares of vector slices;
 *   sumslolo  lowest parts of computed sums of squares of vector slices;
 *   nbsums    the number of elements in sums equals
 *             the number of blocks in the kernel launch;
 *   nbsumsLog2 is ceil(log2((double) nbsums), used in sum reduction;
 *   BS        is the block size, the number of threads in a block.
 *
 * ON RETURN :
 *   vhihi     highest parts of the quad double vector in the same direction 
 *             of the original vector but with norm one;
 *   vlohi     second highest parts of the quad double vector in the same
 *             direction of the original vector but with norm one;
 *   vhilo     second lowest parts of the quad double vector in the same
 *             direction of the original vector but with norm one;
 *   vlolo     lowest parts of the quad double vector in the same direction 
 *             of the original vector but with norm one;
 *   normhihi  highest part of the 2-norm of the original vector;
 *   normlohi  second highest part of the 2-norm of the original vector;
 *   normhilo  second lowest part of the 2-norm of the original vector;
 *   normlolo  lowest part of the 2-norm of the original vector. */

void GPU_norm
 ( double *vhihi_h, double *vlohi_h, double *vhilo_h, double *vlolo_h,
   int dim, int freq, int BS,
   double *normhihi, double *normlohi, double *normhilo,double *normlolo,
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
 *
 * ON ENTRY :
 *   vhihi_h   highest parts of a quad double vector of dimension dim;
 *   vlohi_h   2nd highest parts of a quad double vector of dimension dim;
 *   vhilo_h   2nd lowest parts of a quad double vector of dimension dim;
 *   vlolo_h   lowest parts of a quad double vector of dimension dim;
 *   dim       dimension of the vector v_h;
 *   freq      frequency of the number of kernel launches (for timings);
 *   BS        block size, number of threads per block;
 *   blocked   is 0 or 1, if 0, then only one block will be launched,
 *             if 1, then as many blocks as dim/BS will be lauched and
 *             dim must be a multiple of BS.
 *
 * ON RETURN :
 *   vhihi_h   highest parts of the normalized vector;
 *   vlohi_h   second highest parts of the normalized vector;
 *   vhilo_h   second lowest parts of the normalized vector;
 *   vlolo_h   lowest parts of the normalized vector;
 *   normhihi  highest part of the 2-norm of the original vector;
 *   normlohi  second highest part of the 2-norm of the original vector;
 *   normhilo  second lowest part of the 2-norm of the original vector;
 *   normlolo  lowest part of the 2-norm of the original vector. */

#endif
