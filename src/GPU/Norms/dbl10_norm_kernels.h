// Defines the constants and the prototypes for the kernels 
// to compute the 2-norm of a vector in deca double precision.

#ifndef __DBL10_NORM_KERNELS_H__
#define __DBL10_NORM_KERNELS_H__

/*
  The constant da_shmemsize determines the size of the vectors
  stored in shared memory.  As every thread works on one entry
  in the shared memory vectors, for deca double precision,
  this size is bounded by the number of threads in a block.
  The largest dimension for which the small normalization runs
  is thus the value of da_shmemsize.
  The largest dimension for which the large normalization runs
  is the square of the value of da_shmemsize because the algorithm
  sums two arrays of squares, the squares of the real and the
  squares of the imaginary part, thus very similar to case of
  the small normalization.  For example, for da_shmemsize = 256,
  the largest dimension is thus 256*256 = 65536.
 */

#define da_shmemsize 256

/*
  The constant maxrounds determines the number of rounds
  in the normalization of medium sized vectors.
  The largest dimension for a medium size normalization
  is thus da_shmemsize*maxrounds, for instance: 256*96 = 24576.
 */

#define maxrounds 96

__global__ void small_normalize_vector
 ( double *vrtb, double *vrix, double *vrmi, double *vrrg, double *vrpk, 
   double *vltb, double *vlix, double *vlmi, double *vlrg, double *vlpk,
   int dim, int dimLog2, double *normrtb, double *normrix, double *normrmi,
   double *normrrg, double *normrpk, double *normltb, double *normlix,
   double *normlmi, double *normlrg, double *normlpk );
/*
 * DESCRIPTION :
 *   Kernel function to compute the norm of the vector,
 *   given by ten arrays with all parts of the deca double numbers,
 *   for vectors of small dimension.
 *
 * REQUIRED :
 *   The dimension dim equals the block size,
 *   which is the number of threads in the block.
 *
 * ON ENTRY :
 *   vrtb      highest parts of a deca double vector of size dim;
 *   vrix      second highest parts of a deca double vector of size dim;
 *   vrmi      third highest parts of a deca double vector of size dim;
 *   vrrg      fourth highest parts of a deca double vector of size dim;
 *   vrpk      fifth highest parts of a deca double vector of size dim;
 *   vltb      fifth lowest parts of a deca double vector of size dim;
 *   vlix      fourth lowest parts of a deca double vector of size dim;
 *   vlmi      third lowest parts of a deca double vector of size dim;
 *   vlrg      second lowest parts of a deca double vector of size dim;
 *   vlpk      lowest parts of a deca double vector of size dim;
 *   dim       the dimension of the vector must equal the block size;
 *   dimLog2   equals ceil(log2((double) dim), used in sum reduction.
 *
 * ON RETURN :
 *   vrtb      highest parts of the deca double vector in the same direction 
 *             of the original vector but with norm one;
 *   vrix      second highest parts of the deca double vector in the same
 *             direction of the original vector but with norm one;
 *   vrmi      third highest parts of the deca double vector in the same
 *             direction of the original vector but with norm one;
 *   vrrg      fourth highest parts of the deca double vector in the same
 *             direction of the original vector but with norm one;
 *   vrpk      fifth highest parts of the deca double vector in the same
 *             direction of the original vector but with norm one;
 *   vltb      fifth lowest parts of the deca double vector in the same
 *             direction of the original vector but with norm one;
 *   vlix      fourth lowest parts of the deca double vector in the same
 *             direction of the original vector but with norm one;
 *   vlmi      third lowest parts of the deca double vector in the same
 *             direction of the original vector but with norm one;
 *   vlrg      second lowest parts of the deca double vector in the same
 *             direction of the original vector but with norm one;
 *   vlpk      lowest parts of the deca double vector in the same direction 
 *             of the original vector but with norm one;
 *   normrtb   highest part of the 2-norm of the vector v;
 *   normrix   second highest part of the 2-norm of the vector v;
 *   normrmi   third highest part of the 2-norm of the vector v;
 *   normrrg   fourth highest lowest part of the 2-norm of the vector v;
 *   normrpk   fifth highest part of the 2-norm of the vector v;
 *   normltb   fifth lowest part of the 2-norm of the vector v;
 *   normlix   fourth lowest part of the 2-norm of the vector v;
 *   normlmi   third lowest part of the 2-norm of the vector v;
 *   normlrg   second lowest lowest part of the 2-norm of the vector v;
 *   normlpk   lowest part of the 2-norm of the vector v. */

__global__ void medium_normalize_vector
 ( double *vrtb, double *vrix, double *vrmi, double *vrrg, double *vrpk,
   double *vltb, double *vlix, double *vlmi, double *vlrg, double *vlpk,
   int dim, int rnd, int rndLog2, int BS, int BSLog2,
   double *normrtb, double *normrix, double *normrmi, double *normrrg,
   double *normrpk, double *normltb, double *normlix, double *normlmi,
   double *normlrg, double *normlpk );
/*
 * DESCRIPTION :
 *   Kernel function to compute the norm of the vector,
 *   given by five arrays for all parts of the deca double numbers,
 *   for vectors of medium dimension.
 *
 * REQUIRED :
 *   The dimension dim equals the block size BS times rnd.
 *   The block size BS is the number of threads in the block.
 *
 * ON ENTRY :
 *   vrtb      highest parts of a deca double vector of size dim;
 *   vrix      second highest parts of a deca double vector of size dim;
 *   vrmi      third highest parts of a deca double vector of size dim;
 *   vrrg      fourth highest parts of a deca double vector of size dim;
 *   vrpk      fifth highest parts of a deca double vector of size dim;
 *   vltb      fifth lowest parts of a deca double vector of size dim;
 *   vlix      fourth lowest parts of a deca double vector of size dim;
 *   vlmi      third lowest parts of a deca double vector of size dim;
 *   vlrg      second lowest parts of a deca double vector of size dim;
 *   vlpk      lowest parts of a deca double vector of size dim;
 *   dim       dimension of the vector must equal rnd*BS;
 *   rnd       the number of rounds or the multiplier for the dimension;
 *   rndLog2   the 2-logarithm of rnd for use in sum reduction;
 *   BS        the block size or the number of threads in the block;
 *   BSLog2    equals ceil(log2((double) BS), used in sum reduction.
 *
 * ON RETURN :
 *   vrtb      highest parts of the deca double vector in the same direction 
 *             of the original vector but with norm one;
 *   vrix      second highest parts of the deca double vector in the same
 *             direction of the original vector but with norm one;
 *   vrmi      third highest parts of the deca double vector in the same
 *             direction of the original vector but with norm one;
 *   vrrg      fourth highest parts of the deca double vector in the same
 *             direction of the original vector but with norm one;
 *   vrpk      fifth highest parts of the deca double vector in the same
 *             direction of the original vector but with norm one;
 *   vltb      fifth lowest parts of the deca double vector in the same
 *             direction of the original vector but with norm one;
 *   vlix      fourth lowest parts of the deca double vector in the same
 *             direction of the original vector but with norm one;
 *   vlmi      third lowest parts of the deca double vector in the same
 *             direction of the original vector but with norm one;
 *   vlrg      second lowest parts of the deca double vector in the same
 *             direction of the original vector but with norm one;
 *   vlpk      lowest parts of the deca double vector in the same direction 
 *             of the original vector but with norm one;
 *   normrtb   highest part of the 2-norm of the original vector v;
 *   normrix   second highest part of the 2-norm of the original vector v;
 *   normrmi   third highest part of the 2-norm of the original vector v;
 *   normrrg   fourth highest part of the 2-norm of the original vector v;
 *   normrpk   fifth highest part of the 2-norm of the original vector v;
 *   normltb   fifth lowest part of the 2-norm of the original vector v;
 *   normlix   fourth lowest part of the 2-norm of the original vector v;
 *   normlmi   third lowest part of the 2-norm of the original vector v;
 *   normlrg   second lowest part of the 2-norm of the original vector v;
 *   normlpk   lowest part of the 2-norm of the original vector v. */

__global__ void large_sum_the_squares
 ( double *vrtb, double *vrix, double *vrmi, double *vrrg, double *vrpk,
   double *vltb, double *vlix, double *vlmi, double *vlrg, double *vlpk,
   int dim, double *sumsrtb, double *sumsrix, double *sumsrmi,
   double *sumsrrg, double *sumsrpk, double *sumsltb, double *sumslix,
   double *sumslmi, double *sumslrg, double *sumslpk, int BS, int BSLog2 );
/*
 * DESCRIPTION :
 *   Computes the sums of the squares of the numbers in a vector,
 *   given by five arrays for all parts of the deca double numbers,
 *   as needed in the 2-norm of the vector, with many blocks.
 *
 * REQUIRED :
 *   Space in sums is allocated for as many blocks in the launch.
 *
 * ON ENTRY :
 *   vrtb      highest parts of a deca double vector of size dim;
 *   vrix      second highest parts of a deca double vector of size dim;
 *   vrmi      third highest parts of a deca double vector of size dim;
 *   vrrg      fourth highest parts of a deca double vector of size dim;
 *   vrpk      fifth highest parts of a deca double vector of size dim;
 *   vltb      fifth lowest parts of a deca double vector of size dim;
 *   vlix      fourth lowest parts of a deca double vector of size dim;
 *   vlmi      third lowest parts of a deca double vector of size dim;
 *   vlrg      second lowest parts of a deca double vector of size dim;
 *   vlpk      lowest parts of a deca double vector of size dim;
 *   dim       dimension of the vector must equal BS*nbsums;
 *   BS        the block size or the number of threads in the block;
 *   BSLog2    equals ceil(log2((double) BS), used in sum reduction.
 *
 * ON RETURN :
 *   sumsrtb   highest parts of computed sums of squares of vector slices;
 *   sumsrix   2nd highest parts of computed sums of squares of vector slices;
 *   sumsrmi   3rd highest parts of computed sums of squares of vector slices;
 *   sumsrrg   4th highest parts of computed sums of squares of vector slices;
 *   sumsrpk   5th highest parts of computed sums of squares of vector slices;
 *   sumsltb   5th lowest parts of computed sums of squares of vector slices;
 *   sumslix   4th lowest parts of computed sums of squares of vector slices;
 *   sumslmi   3rd lowest parts of computed sums of squares of vector slices;
 *   sumslrg   2nd lowest parts of computed sums of squares of vector slices;
 *   sumslpk   lowest parts of computed sums of squares of vector slices. */

__global__ void large_normalize_vector
 ( double *vrtb, double *vrix, double *vrmi, double *vrrg, double *vrpk,
   double *vltb, double *vlix, double *vlmi, double *vlrg, double *vlpk,
   int dim, double *sumsrtb, double *sumsrix, double *sumsrmi,
   double *sumsrrg, double *sumsrpk, double *sumsltb, double *sumslix,
   double *sumslmi, double *sumslrg, double *sumslpk,
   int nbsums, int nbsumsLog2, int BS, double *normrtb, double *normrix,
   double *normrmi, double *normrrg, double *normrpk, double *normltb,
   double *normlix, double *normlmi, double *normlrg, double *normlpk );
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
 *   vrtb      highest parts of a deca double vector of size dim;
 *   vrix      second highest parts of a deca double vector of size dim;
 *   vrmi      third highest parts of a deca double vector of size dim;
 *   vrrg      fourth highest parts of a deca double vector of size dim;
 *   vrpk      fifth highest parts of a deca double vector of size dim;
 *   vltb      fifth lowest parts of a deca double vector of size dim;
 *   vlix      fourth lowest parts of a deca double vector of size dim;
 *   vlmi      third lowest parts of a deca double vector of size dim;
 *   vlrg      second lowest parts of a deca double vector of size dim;
 *   vlpk      lowest parts of a deca double vector of size dim;
 *   dim       dimension of the vector must equal BS*nbsums;
 *   sumsrtb   highest parts of computed sums of squares of vector slices;
 *   sumsrix   2nd highest parts of computed sums of squares of vector slices;
 *   sumsrmi   3rd highest parts of computed sums of squares of vector slices;
 *   sumsrrg   4th highest parts of computed sums of squares of vector slices;
 *   sumsrpk   5th highest parts of computed sums of squares of vector slices;
 *   sumsltb   5th lowest parts of computed sums of squares of vector slices;
 *   sumslix   4th lowest parts of computed sums of squares of vector slices;
 *   sumslmi   3rd lowest parts of computed sums of squares of vector slices;
 *   sumslrg   2nd lowest parts of computed sums of squares of vector slices;
 *   sumslpk   lowest parts of computed sums of squares of vector slices;
 *   nbsums    the number of elements in sums equals
 *             the number of blocks in the kernel launch;
 *   nbsumsLog2 is ceil(log2((double) nbsums), used in sum reduction;
 *   BS        is the block size, the number of threads in a block.
 *
 * ON RETURN :
 *   vrtb      highest parts of the deca double vector in the same direction 
 *             of the original vector but with norm one;
 *   vrix      second highest parts of the deca double vector in the same
 *             direction of the original vector but with norm one;
 *   vrmi      third highest parts of the deca double vector in the same
 *             direction of the original vector but with norm one;
 *   vrrg      fourth highest parts of the deca double vector in the same
 *             direction of the original vector but with norm one;
 *   vrpk      fifth highest parts of the deca double vector in the same
 *             direction of the original vector but with norm one;
 *   vltb      fifth lowest parts of the deca double vector in the same
 *             direction of the original vector but with norm one;
 *   vlix      fourth lowest parts of the deca double vector in the same
 *             direction of the original vector but with norm one;
 *   vlmi      third lowest parts of the deca double vector in the same
 *             direction of the original vector but with norm one;
 *   vlrg      second lowest parts of the deca double vector in the same
 *             direction of the original vector but with norm one;
 *   vlpk      lowest parts of the deca double vector in the same direction 
 *             of the original vector but with norm one;
 *   normrtb   highest part of the 2-norm of the original vector v;
 *   normrix   second highest part of the 2-norm of the original vector v;
 *   normrmi   third highest part of the 2-norm of the original vector v;
 *   normrrg   fourth highest part of the 2-norm of the original vector v;
 *   normrpk   fifth highest part of the 2-norm of the original vector v;
 *   normltb   fifth lowest part of the 2-norm of the original vector v;
 *   normlix   fourth lowest part of the 2-norm of the original vector v;
 *   normlmi   third lowest part of the 2-norm of the original vector v;
 *   normlrg   second lowest part of the 2-norm of the original vector v;
 *   normlpk   lowest part of the 2-norm of the original vector v. */

void GPU_norm
 ( double *vrtb_h, double *vrix_h, double *vrmi_h, double *vrrg_h,
   double *vrpk_h, double *vltb_h, double *vlix_h, double *vlmi_h,
   double *vlrg_h, double *vlpk_h, int dim, int freq, int BS,
   double *normrtb, double *normrix, double *normrmi, double *normrrg,
   double *normrpk, double *normltb, double *normlix, double *normlmi,
   double *normlrg, double *normlpk, int blocked );
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
 *   vrtb_h    highest parts of a deca double vector of size dim;
 *   vrix_h    second highest parts of a deca double vector of size dim;
 *   vrmi_h    third highest parts of a deca double vector of size dim;
 *   vrrg_h    fourth highest parts of a deca double vector of size dim;
 *   vrpk_h    fifth highest parts of a deca double vector of size dim;
 *   vltb_h    fifth lowest parts of a deca double vector of size dim;
 *   vlix_h    fourth lowest parts of a deca double vector of size dim;
 *   vlmi_h    third lowest parts of a deca double vector of size dim;
 *   vlrg_h    second lowest parts of a deca double vector of size dim;
 *   vlpk_h    lowest parts of a deca double vector of size dim;
 *   dim       dimension of the vector v_h;
 *   freq      frequency of the number of kernel launches (for timings);
 *   BS        block size, number of threads per block;
 *   blocked   is 0 or 1, if 0, then only one block will be launched,
 *             if 1, then as many blocks as dim/BS will be lauched and
 *             dim must be a multiple of BS.
 *
 * ON RETURN :
 *   vrtb_h    highest parts of the normalized vector;
 *   vrix_h    second highest parts of the normalized vector;
 *   vrmi_h    third highest parts of the normalized vector;
 *   vrrg_h    fourth highest parts of the normalized vector;
 *   vrpk_h    fifth highest parts of the normalized vector;
 *   vltb_h    fifth lowest parts of the normalized vector;
 *   vlix_h    fourth lowest parts of the normalized vector;
 *   vlmi_h    third lowest parts of the normalized vector;
 *   vlrg_h    second lowest parts of the normalized vector;
 *   vlpk_h    lowest parts of the normalized vector;
 *   normrtb   highest part of the 2-norm of the original vector v;
 *   normrix   second highest part of the 2-norm of the original vector v;
 *   normrmi   third highest part of the 2-norm of the original vector v;
 *   normrrg   fourth highest part of the 2-norm of the original vector v;
 *   normrpk   fifth highest part of the 2-norm of the original vector v;
 *   normltb   fifth lowest part of the 2-norm of the original vector v;
 *   normlix   fourth lowest part of the 2-norm of the original vector v;
 *   normlmi   third lowest part of the 2-norm of the original vector v;
 *   normlrg   second lowest part of the 2-norm of the original vector v;
 *   normlpk   lowest part of the 2-norm of the original vector v. */

#endif
