// The file dbl5_sqrt_kernels.h defines the prototypes for the square root
// in penta double precision.

#ifndef __DBL5_SQRT_KERNELS_H__
#define __DBL5_SQRT_KERNELS_H__

__global__ void dbl5_sqrt
 ( double *tb, double *ix, double *mi, double *rg, double *pk,
   int max_steps );
/*
 * DESCRIPTION :
 *   Applies Newton's method to compute the square root of the
 *   penta double number given by the five parts,
 *   respectively in tb, ix, mi, rg, and pk.
 *   Executes as many steps as the value of max_steps.
 *   Returns the value of the square root. */

void GPU_dbl5_sqrt
 ( double *tb, double *ix, double *mi, double *rg, double *pk,
   int maxstp );
/*
 * DESCRIPTION :
 *   Calls Newton's method to compute the square root of a penta double.
 *
 * ON ENTRY :
 *   tb       highest part of a penta double;
 *   ix       second highest part of a penta double;
 *   mi       middle part of a penta double;
 *   rg       second lowest part of a penta double;
 *   pk       lowest part of a penta double;
 *   maxstp   the total number of iterations.
 *
 * ON RETURN :
 *   tb       highest part of the square root of a penta double;
 *   ix       second highest part of the square root of a penta double;
 *   mi       middle part of the square root of a penta double;
 *   rg       second lowest part of the square root of a penta double;
 *   pk       lowest part of the square root of a penta double. */

#endif
