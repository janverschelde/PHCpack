// The file dbl4_sqrt_kernels.h defines the prototypes for the square root
// in quad double precision.

#ifndef __DBL4_SQRT_KERNELS_H__
#define __DBL4_SQRT_KERNELS_H__

__global__ void dbl4_sqrt
 ( double *hihi, double *lohi, double *hilo, double *lolo, int max_steps );
/*
 * DESCRIPTION :
 *   Applies Newton's method to compute the square root of the
 *   quad double number given by the four parts,
 *   respectively in hihi, lohi, hilo, and lolo.
 *   Executes as many steps as the value of max_steps.
 *   Returns the value of the square root. */

void GPU_dbl4_sqrt
 ( double *hihi, double *lohi, double *hilo, double *lolo, int maxstp );
/*
 * DESCRIPTION :
 *   Calls Newton's method to compute the square root of a quad double.
 *
 * ON ENTRY :
 *   hihi      highest part of a quad double;
 *   lohi      second highest part of a quad double;
 *   hilo      second lowest part of a quad double;
 *   lolo      lowest part of a quad double;
 *   maxstp    the total number of iterations.
 *
 * ON RETURN :
 *   hihi      highest part of the square root of a quad double;
 *   lohi      second highest part of the square root of a quad double;
 *   hilo      second lowest part of the square root of a quad double;
 *   lolo      lowest part of the square root of a quad double. */

#endif
