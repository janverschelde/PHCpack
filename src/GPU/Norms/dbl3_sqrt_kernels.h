// The file dbl3_sqrt_kernels.h defines the prototypes for the square root
// in triple double precision.

#ifndef __DBL3_SQRT_KERNELS_H__
#define __DBL3_SQRT_KERNELS_H__

__global__ void dbl3_sqrt
 ( double *hi, double *mi, double *lo, int max_steps );
/*
 * DESCRIPTION :
 *   Applies Newton's method to compute the square root of the
 *   triple double number given by the high, middle and low parts,
 *   respectively in hi, mi, and lo.
 *   Executes as many steps as the value of max_steps.
 *   Returns in hi, mi, and lo the value of the square root. */

void GPU_dbl3_sqrt
 ( double *hi, double *mi, double *lo, int maxstp );
/*
 * DESCRIPTION :
 *   Calls Newton's method to compute the square root of a triple double.
 *
 * ON ENTRY :
 *   hi       high part of a triple double;
 *   mi       middle part of a triple double;
 *   lo       low part of a triple double;
 *   maxstp   the total number of iterations.
 *
 * ON RETURN :
 *   hi       high part of the square root of a triple double;
 *   mi       middle part of the square root of a triple double;
 *   lo       low part of the square root of a triple double. */

#endif
