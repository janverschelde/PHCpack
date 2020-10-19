// Defines the prototypes and constants for the kernels
// to compute the 2-norm of a complex vector in deca double precision.

#ifndef __CMPLX10_NORM_KERNELS_H__
#define __CMPLX10_NORM_KERNELS_H__

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
  the small normalization.  For example, for da_shmemsize = 128,
  the largest dimension is thus 128*128 = 16384.
 */

#define da_shmemsize 128

/*
  The constant maxrounds determines the number of rounds
  in the normalization of medium sized vectors.
  The largest dimension for a medium size normalization
  is thus da_shmemsize*maxrounds, for instance: 128*96 = 12288.
 */

#define maxrounds 96

__global__ void small_normalize_vector
 ( double *vrertb, double *vrerix, double *vrermi, double *vrerrg,
   double *vrerpk, double *vreltb, double *vrelix, double *vrelmi,
   double *vrelrg, double *vrelpk,
   double *vimrtb, double *vimrix, double *vimrmi, double *vimrrg,
   double *vimrpk, double *vimltb, double *vimlix, double *vimlmi,
   double *vimlrg, double *vimlpk, int dim, int dimLog2,
   double *normrtb, double *normrix, double *normrmi, double *normrrg,
   double *normrpk, double *normltb, double *normlix, double *normlmi,
   double *normlrg, double *normlpk );
/*
 * DESCRIPTION :
 *   Kernel function to compute the norm of the complex vector,
 *   given by ten arrays of real parts and ten arrays of imaginary parts,
 *   for vectors of small dimension, in deca double precision.
 *
 * REQUIRED :
 *   The dimension dim equals the block size,
 *   which is the number of threads in the block.
 *   Moreover, all arrays have size dim.
 *
 * ON ENTRY :
 *   vrertb    highest real parts of a deca double vector;
 *   vrerix    second highest real parts of a deca double vector;
 *   vrermi    third highest real parts of a deca double vector;
 *   vrerrg    fourth highest real parts of a deca double vector;
 *   vrerpk    fifth highest real parts of a deca double vector;
 *   vreltb    fifth lowest real parts of a deca double vector;
 *   vrelix    fourth lowest real parts of a deca double vector;
 *   vrelmi    third lowest real parts of a deca double vector;
 *   vrelrg    second lowest real parts of a deca double vector;
 *   vrelpk    lowest real parts of a deca double vector;
 *   vimrtb    highest imaginary parts of a deca double vector;
 *   vimrix    second highest imaginary parts of a deca double vector;
 *   vimrmi    third highest imaginary parts of a deca double vector;
 *   vimrrg    fourth highest imaginary parts of a deca double vector;
 *   vimrpk    fifth highest imaginary parts of a deca double vector;
 *   vimltb    fifth lowest imaginary parts of a deca double vector;
 *   vimlix    fourth lowest imaginary parts of a deca double vector;
 *   vimlmi    third lowest imaginary parts of a deca double vector;
 *   vimlrg    second lowest imaginary parts of a deca double vector;
 *   vimlpk    lowest imaginary parts of a deca double vector;
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
 *   normrtb   highest part of the 2-norm of the complex vector;
 *   normrix   second highest part of the 2-norm;
 *   normrmi   third highest part of the 2-norm;
 *   normrrg   fourth highest part of the 2-norm;
 *   normrpk   fifth highest part of the 2-norm;
 *   normltb   fifth lowest part of the 2-norm;
 *   normlix   fourth lowest part of the 2-norm;
 *   normlmi   third lowest part of the 2-norm;
 *   normlrg   second lowest part of the 2-norm;
 *   normlpk   lowest part of the 2-norm. */

__global__ void medium_normalize_vector
 ( double *vrertb, double *vrerix, double *vrermi, double *vrerrg,
   double *vrerpk, double *vreltb, double *vrelix, double *vrelmi,
   double *vrelrg, double *vrelpk,
   double *vimrtb, double *vimrix, double *vimrmi, double *vimrrg,
   double *vimrpk, double *vimltb, double *vimlix, double *vimlmi,
   double *vimlrg, double *vimlpk,
   int dim, int rnd, int rndLog2, int BS, int BSLog2,
   double *normrtb, double *normrix, double *normrmi, double *normrrg,
   double *normrpk, double *normltb, double *normlix, double *normlmi,
   double *normlrg, double *normlpk );
/*
 * DESCRIPTION :
 *   Kernel function to compute the norm of the vector in v,
 *   for vectors of medium dimension, in deca double precision.
 *
 * REQUIRED :
 *   The dimension dim equals the block size BS times rnd.
 *   The block size BS is the number of threads in the block.
 *   All arrays are of size dim.
 *
 * ON ENTRY :
 *   vrertb    highest real parts of a deca double vector;
 *   vrerix    second highest real parts of a deca double vector;
 *   vrermi    third highest real parts of a deca double vector;
 *   vrerrg    fourth highest real parts of a deca double vector;
 *   vrerpk    fifth highest real parts of a deca double vector;
 *   vreltb    fifth lowest real parts of a deca double vector;
 *   vrelix    fourth lowest real parts of a deca double vector;
 *   vrelmi    third lowest real parts of a deca double vector;
 *   vrelrg    second lowest real parts of a deca double vector;
 *   vrelpk    lowest real parts of a deca double vector;
 *   vimrtb    highest imaginary parts of a deca double vector;
 *   vimrix    second highest imaginary parts of a deca double vector;
 *   vimrmi    third highest imaginary parts of a deca double vector;
 *   vimrrg    fourth highest imaginary parts of a deca double vector;
 *   vimrpk    fifth highest imaginary parts of a deca double vector;
 *   vimltb    fifth lowest imaginary parts of a deca double vector;
 *   vimlix    fourth lowest imaginary parts of a deca double vector;
 *   vimlmi    third lowest imaginary parts of a deca double vector;
 *   vimlrg    second lowest imaginary parts of a deca double vector;
 *   vimlpk    lowest imaginary parts of a deca double vector;
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
 *   normrtb   highest part of the 2-norm of the complex vector;
 *   normrix   second highest part of the 2-norm;
 *   normrmi   third highest part of the 2-norm;
 *   normrrg   fourth highest part of the 2-norm;
 *   normrpk   fifth highest part of the 2-norm;
 *   normltb   fifth lowest part of the 2-norm;
 *   normlix   fourth lowest part of the 2-norm;
 *   normlmi   third lowest part of the 2-norm;
 *   normlrg   second lowest part of the 2-norm;
 *   normlpk   lowest part of the 2-norm. */

__global__ void large_sum_the_squares
 ( double *vrertb, double *vrerix, double *vrermi, double *vrerrg,
   double *vrerpk, double *vreltb, double *vrelix, double *vrelmi,
   double *vrelrg, double *vrelpk,
   double *vimrtb, double *vimrix, double *vimrmi, double *vimrrg,
   double *vimrpk, double *vimltb, double *vimlix, double *vimlmi,
   double *vimlrg, double *vimlpk, int dim,
   double *sumsrtb, double *sumsrix, double *sumsrmi, double *sumsrrg,
   double *sumsrpk, double *sumsltb, double *sumslix, double *sumslmi,
   double *sumslrg, double *sumslpk, int BS, int BSLog2 );
/*
 * DESCRIPTION :
 *   Computes the sums of the squares of the numbers in v,
 *   as needed in the 2-norm of the vector, with many blocks.
 *
 * REQUIRED :
 *   Space in sums is allocated for as many blocks in the launch.
 *
 * ON ENTRY :
 *   vrertb    highest real parts of a deca double vector;
 *   vrerix    second highest real parts of a deca double vector;
 *   vrermi    third highest real parts of a deca double vector;
 *   vrerrg    fourth highest real parts of a deca double vector;
 *   vrerpk    fifth highest real parts of a deca double vector;
 *   vreltb    fifth lowest real parts of a deca double vector;
 *   vrelix    fourth lowest real parts of a deca double vector;
 *   vrelmi    third lowest real parts of a deca double vector;
 *   vrelrg    second lowest real parts of a deca double vector;
 *   vrelpk    lowest real parts of a deca double vector;
 *   vimrtb    highest imaginary parts of a deca double vector;
 *   vimrix    second highest imaginary parts of a deca double vector;
 *   vimrmi    third highest imaginary parts of a deca double vector;
 *   vimrrg    fourth highest imaginary parts of a deca double vector;
 *   vimrpk    fifth highest imaginary parts of a deca double vector;
 *   vimltb    fifth lowest imaginary parts of a deca double vector;
 *   vimlix    fourth lowest imaginary parts of a deca double vector;
 *   vimlmi    third lowest imaginary parts of a deca double vector;
 *   vimlrg    second lowest imaginary parts of a deca double vector;
 *   vimlpk    lowest imaginary parts of a deca double vector;
 *   dim       dimension of the vector must equal BS*nbsums;
 *   BS        the block size or the number of threads in the block;
 *   BSLog2    equals ceil(log2((double) BS), used in sum reduction.
 *
 * ON RETURN :
 *   sumsrtb    highest parts of sums of squares of slices of v;
 *   sumsrix    second highest parts of sums of squares of slices of v;
 *   sumsrmi    third highest parts of sums of squares of slices of v;
 *   sumsrrg    fourth highest parts of sums of squares of slices of v;
 *   sumsrpk    fifth highest parts of sums of squares of slices of v;
 *   sumsltb    fifth lowest parts of sums of squares of slices of v;
 *   sumslix    fourth lowest parts of sums of squares of slices of v;
 *   sumslmi    third lowest parts of sums of squares of slices of v;
 *   sumslrg    second lowest parts of sums of squares of slices of v;
 *   sumslpk    lowest parts of sums of squares of slices of v. */

__global__ void large_normalize_vector
 ( double *vrertb, double *vrerix, double *vrermi, double *vrerrg,
   double *vrerpk, double *vreltb, double *vrelix, double *vrelmi,
   double *vrelrg, double *vrelpk,
   double *vimrtb, double *vimrix, double *vimrmi, double *vimrrg,
   double *vimrpk, double *vimltb, double *vimlix, double *vimlmi,
   double *vimlrg, double *vimlpk, int dim,
   double *sumsrtb, double *sumsrix, double *sumsrmi, double *sumsrrg,
   double *sumsrpk, double *sumsltb, double *sumslix, double *sumslmi,
   double *sumslrg, double *sumslpk,
   int nbsums, int nbsumsLog2, int BS,
   double *normrtb, double *normrix, double *normrmi, double *normrrg,
   double *normrpk, double *normltb, double *normlix, double *normlmi,
   double *normlrg, double *normlpk );
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
 *   vrertb    highest real parts of a deca double vector;
 *   vrerix    second highest real parts of a deca double vector;
 *   vrermi    third highest real parts of a deca double vector;
 *   vrerrg    fourth highest real parts of a deca double vector;
 *   vrerpk    fifth highest real parts of a deca double vector;
 *   vreltb    fifth lowest real parts of a deca double vector;
 *   vrelix    fourth lowest real parts of a deca double vector;
 *   vrelmi    third lowest real parts of a deca double vector;
 *   vrelrg    second lowest real parts of a deca double vector;
 *   vrelpk    lowest real parts of a deca double vector;
 *   vimrtb    highest imaginary parts of a deca double vector;
 *   vimrix    second highest imaginary parts of a deca double vector;
 *   vimrmi    third highest imaginary parts of a deca double vector;
 *   vimrrg    fourth highest imaginary parts of a deca double vector;
 *   vimrpk    fifth highest imaginary parts of a deca double vector;
 *   vimltb    fifth lowest imaginary parts of a deca double vector;
 *   vimlix    fourth lowest imaginary parts of a deca double vector;
 *   vimlmi    third lowest imaginary parts of a deca double vector;
 *   vimlrg    second lowest imaginary parts of a deca double vector;
 *   vimlpk    lowest imaginary parts of a deca double vector;
 *   dim       dimension of the vector must equal BS*nbsums;
 *   sumsrtb   highest parts of sums of squares of slices of v;
 *   sumsrix   second highest parts of sums of squares of slices of v;
 *   sumsrmi   third highest parts of sums of squares of slices of v;
 *   sumsrrg   fourth highest parts of sums of squares of slices of v;
 *   sumsrpk   fifth highest parts of sums of squares of slices of v;
 *   sumsltb   fifth lowest parts of sums of squares of slices of v;
 *   sumslix   fourth lowest parts of sums of squares of slices of v;
 *   sumslmi   third lowest parts of sums of squares of slices of v;
 *   sumslrg   second lowest parts of sums of squares of slices of v;
 *   sumslpk   lowest parts of sums of squares of slices of v;
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
 *   normrtb   highest part of the 2-norm of the given vector;
 *   normrix   second highest part of the 2-norm;
 *   normrmi   third highest part of the 2-norm;
 *   normrrg   fourth highest part of the 2-norm;
 *   normrpk   fifth highest part of the 2-norm;
 *   normltb   fifth lowest part of the 2-norm;
 *   normlix   fourth lowest part of the 2-norm;
 *   normlmi   third lowest part of the 2-norm;
 *   normlrg   second lowest part of the 2-norm;
 *   normlpk   lowest part of the 2-norm. */

void GPU_norm
 ( double *vrertb_h, double *vrerix_h, double *vrermi_h, double *vrerrg_h,
   double *vrerpk_h, double *vreltb_h, double *vrelix_h, double *vrelmi_h,
   double *vrelrg_h, double *vrelpk_h,
   double *vimrtb_h, double *vimrix_h, double *vimrmi_h, double *vimrrg_h,
   double *vimrpk_h, double *vimltb_h, double *vimlix_h, double *vimlmi_h,
   double *vimlrg_h, double *vimlpk_h, int dim, int freq, int BS,
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
 *   All arrays are of dimension dim.
 *
 * ON ENTRY :
 *   vrertb_h    highest real parts of a deca double vector;
 *   vrerix_h    second highest real parts of a deca double vector;
 *   vrermi_h    third highest real parts of a deca double vector;
 *   vrerrg_h    fourth highest real parts of a deca double vector;
 *   vrerpk_h    fifth highest real parts of a deca double vector;
 *   vreltb_h    fifth lowest real parts of a deca double vector;
 *   vrelix_h    fourth lowest real parts of a deca double vector;
 *   vrelmi_h    third lowest real parts of a deca double vector;
 *   vrelrg_h    second lowest real parts of a deca double vector;
 *   vrelpk_h    lowest real parts of a deca double vector;
 *   vimrtb_h    highest imaginary parts of a deca double vector;
 *   vimrix_h    second highest imaginary parts of a deca double vector;
 *   vimrmi_h    third highest imaginary parts of a deca double vector;
 *   vimrrg_h    fourth highest imaginary parts of a deca double vector;
 *   vimrpk_h    fifth highest imaginary parts of a deca double vector;
 *   vimltb_h    fifth lowest imaginary parts of a deca double vector;
 *   vimlix_h    fourth lowest imaginary parts of a deca double vector;
 *   vimlmi_h    third lowest imaginary parts of a deca double vector;
 *   vimlrg_h    second lowest imaginary parts of a deca double vector;
 *   vimlpk_h    lowest imaginary parts of a deca double vector;
 *   dim         dimension of the given vector;
 *   freq        frequency of the number of kernel launches (for timings);
 *   BS          block size, number of threads per block;
 *   blocked     is 0 or 1, if 0, then only one block will be launched,
 *               if 1, then as many blocks as dim/BS will be lauched and
 *               dim must be a multiple of BS.
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
 *   normrtb   highest part of the 2-norm of the given vector;
 *   normrix   second highest part of the 2-norm;
 *   normrmi   third highest part of the 2-norm;
 *   normrrg   fourth highest part of the 2-norm;
 *   normrpk   fifth highest part of the 2-norm;
 *   normltb   fifth lowest part of the 2-norm;
 *   normlix   fourth lowest part of the 2-norm;
 *   normlmi   third lowest part of the 2-norm;
 *   normlrg   second lowest part of the 2-norm;
 *   normlpk   lowest part of the 2-norm. */

#endif
