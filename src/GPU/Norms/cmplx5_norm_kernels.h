// Defines the prototypes and constants for the kernels
// to compute the 2-norm of a complex vector in penta double precision.

#ifndef __CMPLX5_NORM_KERNELS_H__
#define __CMPLX5_NORM_KERNELS_H__

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
  the small normalization.  For example, for pd_shmemsize = 256,
  the largest dimension is thus 256*256 = 65536.
 */

#define pd_shmemsize 256

/*
  The constant maxrounds determines the number of rounds
  in the normalization of medium sized vectors.
  The largest dimension for a medium size normalization
  is thus pd_shmemsize*maxrounds, for instance: 256*192 = 49152.
 */

#define maxrounds 192

__global__ void small_normalize_vector
 ( double *vretb, double *vreix, double *vremi, double *vrerg, double *vrepk,
   double *vimtb, double *vimix, double *vimmi, double *vimrg, double *vimpk,
   int dim, int dimLog2,
   double *normtb, double *normix, double *normmi, double *normrg,
   double *normpk );
/*
 * DESCRIPTION :
 *   Kernel function to compute the norm of the complex vector,
 *   given by a vector of real parts and a vector of imaginary parts,
 *   for vectors of small dimension, in penta double precision.
 *
 * REQUIRED :
 *   The dimension dim equals the block size,
 *   which is the number of threads in the block.
 *   Moreover, all arrays have size dim.
 *
 * ON ENTRY :
 *   vretb     highest real parts of a complex vector;
 *   vreix     second highest real parts of a complex vector;
 *   vremix    middle real parts of a complex vector;
 *   vrerg     second lowest real parts of a complex vector;
 *   vrepk     lowest real parts of a complex vector;
 *   vimtb     highest imaginary parts of a complex vector;
 *   vimix     second highest imaginary parts of a complex vector;
 *   vimmi     middle imaginary parts of a complex vector;
 *   vimrg     second lowest imaginary parts of a complex vector;
 *   vimpk     lowest imaginary parts of a complex vector;
 *   dim       the dimension of the vector must equal the block size;
 *   dimLog2   equals ceil(log2((double) dim), used in sum reduction.
 *
 * ON RETURN :
 *   vretb     highest real parts of the normalized vector;
 *   vreix     second highest real parts of the normalized vector;
 *   vremi     middle real parts of the normalized vector;
 *   vrerg     second lowest real parts of the normalized vector;
 *   vrepk     lowest real parts of the normalized vector;
 *   vimtb     highest imaginary parts of the normalized vector;
 *   vimix     second highest imaginary parts of the normalized vector;
 *   vimmi     middle imaginary parts of the normalized vector;
 *   vimrg     second lowest imaginary parts of the normalized vector;
 *   vimpk     lowest imaginary parts of the normalized vector;
 *   normtb    highest part of the 2-norm of the given vector;
 *   normix    second highest part of the 2-norm of the given vector;
 *   normmi    middle part of the 2-norm of the given vector;
 *   normrg    second lowest part of the 2-norm of the given vector;
 *   normpk    lowest part of the 2-norm of the given vector. */

__global__ void medium_normalize_vector
 ( double *vretb, double *vreix, double *vremi, double *vrerg, double *vrepk,
   double *vimtb, double *vimix, double *vimmi, double *vimrg, double *vimpk,
   int dim, int rnd, int rndLog2, int BS, int BSLog2,
   double *normtb, double *normix, double *normmi, double *normrg,
   double *normpk );
/*
 * DESCRIPTION :
 *   Kernel function to compute the norm of the vector in v,
 *   for vectors of medium dimension, in penta double precision.
 *
 * REQUIRED :
 *   The dimension dim equals the block size BS times rnd.
 *   The block size BS is the number of threads in the block.
 *   All arrays are of size dim.
 *
 * ON ENTRY :
 *   vretb     highest real parts of a complex vector;
 *   vreix     second highest real parts of a complex vector;
 *   vremi     middle real parts of a complex vector;
 *   vrerg     second lowest real parts of a complex vector;
 *   vrepk     lowest real parts of a complex vector;
 *   vimtb     highest imaginary parts of a complex vector;
 *   vimix     second highest imaginary parts of a complex vector;
 *   vimmi     middle imaginary parts of a complex vector;
 *   vimrg     second lowest imaginary parts of a complex vector;
 *   vimpk     lowest imaginary parts of a complex vector;
 *   dim       dimension of the vector must equal rnd*BS;
 *   rnd       the number of rounds or the multiplier for the dimension;
 *   rndLog2   the 2-logarithm of rnd for use in sum reduction;
 *   BS        the block size or the number of threads in the block;
 *   BSLog2    equals ceil(log2((double) BS), used in sum reduction.
 *
 * ON RETURN :
 *   vretb     highest real parts of the normalized vector;
 *   vreix     second highest real parts of the normalized vector;
 *   vremi     middle real parts of the normalized vector;
 *   vrerg     second lowest real parts of the normalized vector;
 *   vrepk     lowest real parts of the normalized vector;
 *   vimtb     highest imaginary parts of the normalized vector;
 *   vimix     second highest imaginary parts of the normalized vector;
 *   vimmi     middle imaginary parts of the normalized vector;
 *   vimrg     second lowest imaginary parts of the normalized vector;
 *   vimpk     lowest imaginary parts of the normalized vector;
 *   normtb    highest part of the 2-norm of the given vector;
 *   normix    second highest part of the 2-norm of the given vector;
 *   normmi    middle part of the 2-norm of the given vector;
 *   normrg    second lowest part of the 2-norm of the given vector;
 *   normpk    lowest part of the 2-norm of the given vector. */

__global__ void large_sum_the_squares
 ( double *vretb, double *vreix, double *vremi, double *vrerg, double *vrepk,
   double *vimtb, double *vimix, double *vimmi, double *vimrg, double *vimpk,
   int dim,
   double *sumstb, double *sumsix, double *sumsmi, double *sumsrg,
   double *sumspk, int BS, int BSLog2 );
/*
 * DESCRIPTION :
 *   Computes the sums of the squares of the numbers in v,
 *   as needed in the 2-norm of the vector, with many blocks.
 *
 * REQUIRED :
 *   Space in sums is allocated for as many blocks in the launch.
 *
 * ON ENTRY :
 *   vretb     highest real parts of a complex vector;
 *   vreix     second highest real parts of a complex vector;
 *   vremi     middle real parts of a complex vector;
 *   vrerg     second lowest real parts of a complex vector;
 *   vrepk     lowest real parts of a complex vector;
 *   vimtb     highest imaginary parts of a complex vector;
 *   vimix     second highest imaginary parts of a complex vector;
 *   vimmi     middle imaginary parts of a complex vector;
 *   vimrg     second lowest imaginary parts of a complex vector;
 *   vimpk     lowest imaginary parts of a complex vector;
 *   dim       dimension of the vector must equal BS*nbsums;
 *   BS        the block size or the number of threads in the block;
 *   BSLog2    equals ceil(log2((double) BS), used in sum reduction.
 *
 * ON RETURN :
 *   sumstb    highest parts of sums of squares of slices of v;
 *   sumsix    second highest parts of sums of squares of slices of v;
 *   sumsmi    middle parts of sums of squares of slices of v;
 *   sumsrg    second lowest parts of sums of squares of slices of v;
 *   sumspk    lowest parts of sums of squares of slices of v. */

__global__ void large_normalize_vector
 ( double *vretb, double *vreix, double *vremi, double *vrerg, double *vrepk,
   double *vimtb, double *vimix, double *vimmi, double *vimrg, double *vimpk,
   int dim,
   double *sumstb, double *sumsix, double *sumsmi, double *sumsrg,
   double *sumspk, int nbsums, int nbsumsLog2, int BS,
   double *normtb, double *normix, double *normmi, double *normrg,
   double *normpk );
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
 *   vretb     highest real parts of a complex vector;
 *   vreix     second highest real parts of a complex vector;
 *   vremi     middle real parts of a complex vector;
 *   vrerg     second lowest real parts of a complex vector;
 *   vrepk     lowest real parts of a complex vector;
 *   vimtb     highest imaginary parts of a complex vector;
 *   vimix     second highest imaginary parts of a complex vector;
 *   vimmi     middle imaginary parts of a complex vector;
 *   vimrg     second lowest imaginary parts of a complex vector;
 *   vimpk     lowest imaginary parts of a complex vector;
 *   dim       dimension of the vector must equal BS*nbsums;
 *   sumstb    highest parts of sums of squares of slices of v;
 *   sumsix    second highest parts of sums of squares of slices of v;
 *   sumsmi    middle parts of sums of squares of slices of v;
 *   sumsrg    second lowest parts of sums of squares of slices of v;
 *   sumspk    lowest parts of sums of squares of slices of v;
 *   nbsums    the number of elements in sums equals
 *             the number of blocks in the kernel launch;
 *   nbsumsLog2 is ceil(log2((double) nbsums), used in sum reduction;
 *   BS        is the block size, the number of threads in a block.
 *
 * ON RETURN :
 *   vretb     highest real parts of the normalized vector;
 *   vreix     second highest real parts of the normalized vector;
 *   vremi     middle real parts of the normalized vector;
 *   vrerg     second lowest real parts of the normalized vector;
 *   vrepk     lowest real parts of the normalized vector;
 *   vimtb     highest imaginary parts of the normalized vector;
 *   vimix     second highest imaginary parts of the normalized vector;
 *   vimmi     middle imaginary parts of the normalized vector;
 *   vimrg     second lowest imaginary parts of the normalized vector;
 *   vimpk     lowest imaginary parts of the normalized vector;
 *   normtb    highest part of the 2-norm of the given vector;
 *   normix    second highest part of the 2-norm of the given vector;
 *   normmi    middle part of the 2-norm of the given vector;
 *   normrg    second lowest part of the 2-norm of the given vector;
 *   normpk    lowest part of the 2-norm of the given vector. */

void GPU_norm
 ( double *vretb_h, double *vreix_h, double *vremi_h, double *vrerg_h,
   double *vrepk_h,
   double *vimtb_h, double *vimix_h, double *vimmi_h, double *vimrg_h,
   double *vimpk_h,
   int dim, int freq, int BS,
   double *normtb, double *normix, double *normmi, double *normrg,
   double *normpk, int blocked );
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
 *   vretb_h   highest real parts of a complex vector;
 *   vreix_h   second highest real parts of a complex vector;
 *   vremi_h   middle real parts of a complex vector;
 *   vrerg_h   second lowest real parts of a complex vector;
 *   vrepk_h   lowest real parts of a complex vector;
 *   vimtb_h   highest imaginary parts of a complex vector;
 *   vimix_h   second highest imaginary parts of a complex vector;
 *   vimmi_h   middle imaginary parts of a complex vector;
 *   vimrg_h   second lowest imaginary parts of a complex vector;
 *   vimpk_h   lowest imaginary parts of a complex vector;
 *   dim       dimension of the given vector;
 *   freq      frequency of the number of kernel launches (for timings);
 *   BS        block size, number of threads per block;
 *   blocked   is 0 or 1, if 0, then only one block will be launched,
 *             if 1, then as many blocks as dim/BS will be lauched and
 *             dim must be a multiple of BS.
 *
 * ON RETURN :
 *   vretb_h   highest real parts of the normalized vector;
 *   vreix_h   second highest real parts of the normalized vector;
 *   vremi_h   middle real parts of the normalized vector;
 *   vrerg_h   second lowest real parts of the normalized vector;
 *   vrepk_h   lowest real parts of the normalized vector;
 *   vimtb_h   highest imaginary parts of the normalized vector;
 *   vimix_h   second highest imaginary parts of the normalized vector;
 *   vimmi_h   middle imaginary parts of the normalized vector;
 *   vimrg_h   second lowest imaginary parts of the normalized vector;
 *   vimpk_h   lowest imaginary parts of the normalized vector;
 *   normtb    highest part the 2-norm of the given vector;
 *   normix    second highest part of the 2-norm of the given vector;
 *   normmi    middle part of the 2-norm of the given vector;
 *   normrg    second lowest part of the 2-norm of the given vector;
 *   normpk    lowest part of the 2-norm of the given vector. */

#endif
