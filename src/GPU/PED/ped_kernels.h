// Kernels to evaluate and differentiate polynomials in several variables,
// written for double, double double, and quad double precision.
// The definitions of the functions are in all_ped_kernels.cu.

#ifndef __ALL_PED_KERNELS_H__
#define __ALL_PED_KERNELS_H__

void GPU_evaldiff
 ( int BS, int dim, int NM, int NV, int deg, int r, int m, int ncoefs,
   char *pos, char *exp,
   complexD<double> *x_h, complexD<double> *c_h,
   complexD<double> *factors_h, complexD<double> *polvalues_h );

void GPU_evaldiff
 ( int BS, int dim, int NM, int NV, int deg, int r, int m, int ncoefs,
   char *pos, char *exp,
   complexD<gdd_real> *x_h, complexD<gdd_real> *c_h,
   complexD<gdd_real> *factors_h, complexD<gdd_real> *polvalues_h );

void GPU_evaldiff
 ( int BS, int dim, int NM, int NV, int deg, int r, int m, int ncoefs,
   char *pos, char *exp,
   complexD<gqd_real> *x_h, complexD<gqd_real> *c_h,
   complexD<gqd_real> *factors_h, complexD<gqd_real> *polvalues_h );
/*
 * The GPU is used to evaluate a system and its Jacobian matrix,
 * in double, double double, and quad double precision.
 *
 * ON ENTRY :
 *   BS           number of threads in a block, the block size;
 *   dim          dimension of the problem;
 *   NM           number of monomials;
 *   NV           number of variables;
 *   deg          highest degree of the variables;
 *   r            frequency of the runs;
 *   m            NM divided by dim;
 *   ncoefs       number of coefficients;
 *   pos          indicate the participating variables in the monomials;
 *   exp          the exponents of the participating variables; 
 *   c            coefficients of the monomials;
 *   x            point where to evaluate.
 *
 * ON RETURN :
 *   factors_h    the factors common to evaluation and differentiation;
 *   polvalues_h  are the values of polynomials and their derivatives. */

#endif
