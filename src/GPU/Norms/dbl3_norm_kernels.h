// Defines the constants and the prototypes for the kernels 
// to compute the 2-norm of a vector in triple double precision.

#ifndef __DBL3_NORM_KERNELS_H__
#define __DBL3_NORM_KERNELS_H__

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
  the small normalization.  For example, for od_shmemsize = 800,
  the largest dimension is thus 800*800 = 640000.
 */

#define td_shmemsize 800

/*
  The constant maxrounds determines the number of rounds
  in the normalization of medium sized vectors.
  The largest dimension for a medium size normalization
  is thus td_shemsize*maxrounds, for instance: 800*320 = 25600.
 */

#define maxrounds 320

__global__ void small_normalize_vector
 ( double *vhi, double *vmi, double *vlo, int dim, int dimLog2,
   double *normhi, double *normmi, double *normlo );
/*
 * DESCRIPTION :
 *   Kernel function to compute the norm of the vector,
 *   given with high parts in vhi, middle parts in vmi, and low parts in vlo,
 *   for vectors of small dimension, in triple double precision.
 *
 * REQUIRED :
 *   The dimension dim equals the block size,
 *   which is the number of threads in the block.
 *
 * ON ENTRY :
 *   vhi       high parts of a triple double vector of dimension dim;
 *   vmi       middle parts of a triple double vector of dimension dim;
 *   vlo       low parts of a triple double vector of dimension dim;
 *   dim       the dimension of the vector must equal the block size;
 *   dimLog2   equals ceil(log2((double) dim), used in sum reduction.
 *
 * ON RETURN :
 *   vhi       high parts of the triple double vector in the same direction 
 *             of the original vector but with norm one;
 *   vmi       middle parts of the triple double vector in the same direction 
 *             of the original vector but with norm one;
 *   vlo       low parts of the triple double vector in the same direction 
 *             of the original vector but with norm one;
 *   normhi    high part of the 2-norm of the original vector;
 *   normmi    middle part of the 2-norm of the original vector;
 *   normlo    low part of the 2-norm of the original vector. */

__global__ void medium_normalize_vector
 ( double *vhi, double *vmi, double *vlo, int dim, int rnd, int rndLog2,
   int BS, int BSLog2, double *normhi, double *normmi, double *normlo );
/*
 * DESCRIPTION :
 *   Kernel function to compute the norm of the vector,
 *   given with high parts in vhi, middle parts in vmi, and low parts in vlo,
 *   for vectors of medium dimension, in triple double precision.
 *
 * REQUIRED :
 *   The dimension dim equals the block size BS times rnd.
 *   The block size BS is the number of threads in the block.
 *
 * ON ENTRY :
 *   vhi       high parts of a triple double vector of dimension dim;
 *   vmi       middle parts of a triple double vector of dimension dim;
 *   vlo       low parts of a triple double vector of dimension dim;
 *   dim       dimension of the vector must equal rnd*BS;
 *   rnd       the number of rounds or the multiplier for the dimension;
 *   rndLog2   the 2-logarithm of rnd for use in sum reduction;
 *   BS        the block size or the number of threads in the block;
 *   BSLog2    equals ceil(log2((double) BS), used in sum reduction.
 *
 * ON RETURN :
 *   vhi       high parts of the triple double vector in the same direction 
 *             of the original vector but with norm one;
 *   vmi       middle parts of the triple double vector in the same direction 
 *             of the original vector but with norm one;
 *   vlo       low parts of the triple double vector in the same direction 
 *             of the original vector but with norm one;
 *   normhi    high part of the 2-norm of the original vector;
 *   normmi    middle part of the 2-norm of the original vector;
 *   normlo    low part of the 2-norm of the original vector. */

__global__ void large_sum_the_squares
 ( double *vhi, double *vmi, double *vlo, int dim,
   double *sumshi, double *sumsmi, double *sumslo, int BS, int BSLog2 );
/*
 * DESCRIPTION :
 *   Computes the sums of the squares of the numbers in a vector,
 *   given with high parts in vhi, middle parts in vmi, and low parts in vlo,
 *   as needed in the 2-norm of the vector, with many blocks.
 *
 * REQUIRED :
 *   Space in sums is allocated for as many blocks in the launch.
 *
 * ON ENTRY :
 *   vhi       high parts of a triple double vector of dimension dim;
 *   vmi       middle parts of a triple double vector of dimension dim;
 *   vlo       low parts of a triple double vector of dimension dim;
 *   dim       dimension of the vector must equal BS*nbsums;
 *   BS        the block size or the number of threads in the block;
 *   BSLog2    equals ceil(log2((double) BS), used in sum reduction.
 *
 * ON RETURN :
 *   sumshi    high parts of computed sums of squares of vector slices;
 *   sumsmi    middle parts of computed sums of squares of vector slices;
 *   sumslo    low parts of computed sums of squares of vector slices. */

__global__ void large_normalize_vector
 ( double *vhi, double *vmi, double *vlo, int dim,
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
 *
 * ON ENTRY :
 *   vhi       high parts of a triple double vector of dimension dim;
 *   vmi       middle parts of a triple double vector of dimension dim;
 *   vlo       low parts of a triple double vector of dimension dim;
 *   dim       dimension of the vector must equal BS*nbsums;
 *   sumshi    high parts of computed sums of squares of vector slices;
 *   sumsmi    middle parts of computed sums of squares of vector slices;
 *   sumslo    low parts of computed sums of squares of vector slices;
 *   nbsums    the number of elements in sums equals
 *             the number of blocks in the kernel launch;
 *   nbsumsLog2 is ceil(log2((double) nbsums), used in sum reduction;
 *   BS        is the block size, the number of threads in a block.
 *
 * ON RETURN :
 *   vhi       high parts of the triple double vector in the same direction 
 *             of the original vector but with norm one;
 *   vmi       middle parts of the triple double vector in the same direction 
 *             of the original vector but with norm one;
 *   vlo       low parts of the triple double vector in the same direction 
 *             of the original vector but with norm one;
 *   normhi    high part of the 2-norm of the original vector;
 *   normmi    middle part of the 2-norm of the original vector;
 *   normlo    low part of the 2-norm of the original vector. */

void GPU_norm
 ( double *vhi_h, double *vmi_h, double *vlo_h, int dim, int freq, int BS,
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
 *
 * ON ENTRY :
 *   vhi_h     high parts of a triple double vector of dimension dim;
 *   vmi_h     middle parts of a triple double vector of dimension dim;
 *   vlo_h     low parts of a triple double vector of dimension dim;
 *   dim       dimension of the vector v_h;
 *   freq      frequency of the number of kernel launches (for timings);
 *   BS        block size, number of threads per block;
 *   blocked   is 0 or 1, if 0, then only one block will be launched,
 *             if 1, then as many blocks as dim/BS will be lauched and
 *             dim must be a multiple of BS.
 *
 * ON RETURN :
 *   vhi_h     high parts of the normalized vector;
 *   vmi_h     middle parts of the normalized vector;
 *   vlo_h     low parts of the normalized vector;
 *   normhi    high part of the 2-norm of the original vector;
 *   normmi    middle part of the 2-norm of the original vector;
 *   normlo    low part of the 2-norm of the original vector. */

#endif
