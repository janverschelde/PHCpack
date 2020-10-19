// Defines the constants and the prototypes for the kernels 
// to compute the 2-norm of a vector in penta double precision.

#ifndef __DBL5_NORM_KERNELS_H__
#define __DBL5_NORM_KERNELS_H__

/*
  The constant pd_shmemsize determines the size of the vectors
  stored in shared memory.  As every thread works on one entry
  in the shared memory vectors, for penta double precision,
  this size is bounded by the number of threads in a block.
  The largest dimension for which the small normalization runs
  is thus the value of pd_shmemsize.
  The largest dimension for which the large normalization runs
  is the square of the value of da_shmemsize because the algorithm
  sums two arrays of squares, the squares of the real and the
  squares of the imaginary part, thus very similar to case of
  the small normalization.  For example, for pd_shmemsize = 512,
  the largest dimension is thus 512*512 = 262144.
 */

#define pd_shmemsize 512

/*
  The constant maxrounds determines the number of rounds
  in the normalization of medium sized vectors.
  The largest dimension for a medium size normalization
  is thus pd_shmemsize*maxrounds, for instance: 512*192 = 98304.
 */

#define maxrounds 192

__global__ void small_normalize_vector
 ( double *vtb, double *vix, double *vmi, double *vrg, double *vpk, int dim,
   int dimLog2, double *normtb, double *normix, double *normmi,
   double *normrg, double *normpk );
/*
 * DESCRIPTION :
 *   Kernel function to compute the norm of the vector,
 *   given by five arrays with all parts of the penta double numbers,
 *   for vectors of small dimension.
 *
 * REQUIRED :
 *   The dimension dim equals the block size,
 *   which is the number of threads in the block.
 *
 * ON ENTRY :
 *   vtb       highest parts of a penta double vector of dimension dim;
 *   vix       second highest parts of a penta double vector of dimension dim;
 *   vmi       middle parts of a penta double vector of dimension dim;
 *   vrg       second lowest parts of a penta double vector of dimension dim;
 *   vpk       lowest parts of a penta double vector of dimension dim;
 *   dim       the dimension of the vector must equal the block size;
 *   dimLog2   equals ceil(log2((double) dim), used in sum reduction.
 *
 * ON RETURN :
 *   vtb       highest parts of the penta double vector in the same direction 
 *             of the original vector but with norm one;
 *   vix       second highest parts of the penta double vector in the same
 *             direction of the original vector but with norm one;
 *   vmi       middle parts of the penta double vector in the same direction 
 *             of the original vector but with norm one;
 *   vrg       second lowest parts of the penta double vector in the same
 *             direction of the original vector but with norm one;
 *   vpk       lowest parts of the penta double vector in the same direction 
 *             of the original vector but with norm one;
 *   normtb    highest part of the 2-norm of the original vector;
 *   normix    second highest part of the 2-norm of the original vector;
 *   normmi    middle part of the 2-norm of the original vector;
 *   normrg    second lowest part of the 2-norm of the original vector;
 *   normpk    lowest part of the 2-norm of the original vector. */

__global__ void medium_normalize_vector
 ( double *vtb, double *vix, double *vmi, double *vrg, double *vpk, int dim,
   int rnd, int rndLog2, int BS, int BSLog2, double *normtb, double *normix,
   double *normmi, double *normrg, double *normpk );
/*
 * DESCRIPTION :
 *   Kernel function to compute the norm of the vector,
 *   given by five arrays for all parts of the penta double numbers,
 *   for vectors of medium dimension.
 *
 * REQUIRED :
 *   The dimension dim equals the block size BS times rnd.
 *   The block size BS is the number of threads in the block.
 *
 * ON ENTRY :
 *   vtb       highest parts of a penta double vector of dimension dim;
 *   vix       second highest parts of a penta double vector of dimension dim;
 *   vmi       middle parts of a penta double vector of dimension dim;
 *   vrg       second lowest parts of a penta double vector of dimension dim;
 *   vpk       lowest parts of a penta double vector of dimension dim;
 *   dim       dimension of the vector must equal rnd*BS;
 *   rnd       the number of rounds or the multiplier for the dimension;
 *   rndLog2   the 2-logarithm of rnd for use in sum reduction;
 *   BS        the block size or the number of threads in the block;
 *   BSLog2    equals ceil(log2((double) BS), used in sum reduction.
 *
 * ON RETURN :
 *   vtb       highest parts of the penta double vector in the same direction 
 *             of the original vector but with norm one;
 *   vix       second highest parts of the penta double vector in the same
 *             direction of the original vector but with norm one;
 *   vmi       middle parts of the penta double vector in the same direction 
 *             of the original vector but with norm one;
 *   vrg       second lowest parts of the penta double vector in the same
 *             direction of the original vector but with norm one;
 *   vpk       lowest parts of the penta double vector in the same direction 
 *             of the original vector but with norm one;
 *   normtb    highest part of the 2-norm of the original vector;
 *   normix    second highest part of the 2-norm of the original vector;
 *   normmi    middle part of the 2-norm of the original vector;
 *   normrg    second lowest part of the 2-norm of the original vector;
 *   normpk    lowest part of the 2-norm of the original vector. */

__global__ void large_sum_the_squares
 ( double *vtb, double *vix, double *vmi, double *vrg, double *vpk, int dim,
   double *sumstb, double *sumsix, double *sumsmi, double *sumsrg,
   double *sumspk, int BS, int BSLog2 );
/*
 * DESCRIPTION :
 *   Computes the sums of the squares of the numbers in a vector,
 *   given by five arrays for all parts of the penta double numbers,
 *   as needed in the 2-norm of the vector, with many blocks.
 *
 * REQUIRED :
 *   Space in sums is allocated for as many blocks in the launch.
 *
 * ON ENTRY :
 *   vtb       highest parts of a penta double vector of dimension dim;
 *   vix       second highest parts of a penta double vector of dimension dim;
 *   vmi       middle parts of a penta double vector of dimension dim;
 *   vrg       second lowest parts of a penta double vector of dimension dim;
 *   vpk       lowest parts of a penta double vector of dimension dim;
 *   dim       dimension of the vector must equal BS*nbsums;
 *   BS        the block size or the number of threads in the block;
 *   BSLog2    equals ceil(log2((double) BS), used in sum reduction.
 *
 * ON RETURN :
 *   sumstb    highest parts of computed sums of squares of vector slices;
 *   sumsix    2nd highest parts of computed sums of squares of vector slices;
 *   sumsmi    middle parts of computed sums of squares of vector slices;
 *   sumsrg    2nd lowest parts of computed sums of squares of vector slices;
 *   sumspk    lowest parts of computed sums of squares of vector slices. */

__global__ void large_normalize_vector
 ( double *vtb, double *vix, double *vmi, double *vrg, double *vpk, int dim,
   double *sumstb, double *sumsix, double *sumsmi, double *sumsrg,
   double *sumspk, int nbsums, int nbsumsLog2, int BS, double *normtb,
   double *normix, double *normmi, double *normrg, double *normpk );
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
 *   vtb       highest parts of a penta double vector of dimension dim;
 *   vix       second highest parts of a penta double vector of dimension dim;
 *   vmi       middle parts of a penta double vector of dimension dim;
 *   vrg       second lowest parts of a penta double vector of dimension dim;
 *   vpk       lowest parts of a penta double vector of dimension dim;
 *   dim       dimension of the vector must equal BS*nbsums;
 *   sumstb    highest parts of computed sums of squares of vector slices;
 *   sumsix    2nd highest parts of computed sums of squares of vector slices;
 *   sumsmi    middle parts of computed sums of squares of vector slices;
 *   sumsrg    2nd low parts of computed sums of squares of vector slices;
 *   sumspk    lowest parts of computed sums of squares of vector slices;
 *   nbsums    the number of elements in sums equals
 *             the number of blocks in the kernel launch;
 *   nbsumsLog2 is ceil(log2((double) nbsums), used in sum reduction;
 *   BS        is the block size, the number of threads in a block.
 *
 * ON RETURN :
 *   vtb       highest parts of the penta double vector in the same direction 
 *             of the original vector but with norm one;
 *   vix       second highest parts of the penta double vector in the same
 *             direction of the original vector but with norm one;
 *   vmi       middle parts of the penta double vector in the same direction 
 *             of the original vector but with norm one;
 *   vrg       second lowest parts of the penta double vector in the same
 *             direction of the original vector but with norm one;
 *   vpk       lowest parts of the penta double vector in the same direction 
 *             of the original vector but with norm one;
 *   normtb    highest part of the 2-norm of the original vector;
 *   normix    second highest part of the 2-norm of the original vector;
 *   normmi    middle part of the 2-norm of the original vector;
 *   normrg    second lowest part of the 2-norm of the original vector;
 *   normpk    lowest part of the 2-norm of the original vector. */

void GPU_norm
 ( double *vtb_h, double *vix_h, double *vmi_h, double *vrg_h, double *vpk_h,
   int dim, int freq, int BS, double *normtb, double *normix, double *normmi,
   double *normrg, double *normpk, int blocked );
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
 *   vtb_h     highest parts of a penta double vector of dimension dim;
 *   vix_h     second highest parts of a penta double vector of dimension dim;
 *   vmi_h     middle parts of a penta double vector of dimension dim;
 *   vrg_h     second lowest parts of a penta double vector of dimension dim;
 *   vpk_h     lowest parts of a penta double vector of dimension dim;
 *   dim       dimension of the vector v_h;
 *   freq      frequency of the number of kernel launches (for timings);
 *   BS        block size, number of threads per block;
 *   blocked   is 0 or 1, if 0, then only one block will be launched,
 *             if 1, then as many blocks as dim/BS will be lauched and
 *             dim must be a multiple of BS.
 *
 * ON RETURN :
 *   vtb_h     highest parts of the normalized vector;
 *   vix_h     second highest parts of the normalized vector;
 *   vmi_h     middle parts of the normalized vector;
 *   vrg_h     second lowest parts of the normalized vector;
 *   vpk_h     lowest parts of the normalized vector;
 *   normtb    highest part of the 2-norm of the original vector;
 *   normix    second highest part of the 2-norm of the original vector;
 *   normmi    middle part of the 2-norm of the original vector;
 *   normrg    second lowest part of the 2-norm of the original vector;
 *   normpk    lowest part of the 2-norm of the original vector. */

#endif
