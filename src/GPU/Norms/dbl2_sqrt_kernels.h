// The file dbl2_sqrt_kernels.h defines the prototypes for the square root
// in double double precision.

#ifndef __DBL2_SQRT_KERNELS_H__
#define __DBL2_SQRT_KERNELS_H__

__global__ void dbl2_sqrt
 ( double *hi, double *lo, int max_steps );
/*
 * DESCRIPTION :
 *   Applies Newton's method to compute the square root of the
 *   double double number given by the high and low parts,
 *   respectively in hi and lo.
 *   Executes as many steps as the value of max_steps.
 *   Returns in hi and lo the value of the square root. */

void GPU_dbl2_sqrt ( double *hi, double *lo, int maxstp );
/*
 * DESCRIPTION :
 *   Calls Newton's method to compute the square root of a double double.
 *
 * ON ENTRY :
 *   hi       high part of a double double;
 *   lo       low part of a double double;
 *   maxstp   the total number of iterations.
 *
 * ON RETURN :
 *   hi       high part of the square root of a double double;
 *   lo       low part of the square root of a double double. */

#endif
