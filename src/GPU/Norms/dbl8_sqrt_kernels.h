// The file dbl8_sqrt_kernels.h defines the prototypes for the square root
// in octo double precision.

#ifndef __DBL8_SQRT_KERNELS_H__
#define __DBL8_SQRT_KERNELS_H__

__global__ void dbl8_sqrt
 ( double *hihihi, double *lohihi, double *hilohi, double *lolohi,
   double *hihilo, double *lohilo, double *hilolo, double *lololo,
   int max_steps );
/*
 * DESCRIPTION :
 *   Applies Newton's method to compute the square root of the
 *   octo double number given by the eight parts, respectively in
 *   hihihi, lohihi, hilohi, lolohi, hihilo, lohilo, hilolo, and lololo.
 *   Executes as many steps as the value of max_steps.
 *   Returns the value of the square root. */

void GPU_dbl8_sqrt
 ( double *hihihi, double *lohihi, double *hilohi, double *lolohi,
   double *hihilo, double *lohilo, double *hilolo, double *lololo,
   int maxstp );
/*
 * DESCRIPTION :
 *   Calls Newton's method to compute the square root of a octo double.
 *
 * ON ENTRY :
 *   hihihi    highest part of a octo double;
 *   lohihi    second highest part of a octo double;
 *   hilohi    third highest part of a octo double;
 *   lohihi    fourth highest part of a octo double;
 *   hihilo    fourth lowest part of a octo double;
 *   lohilo    third lowest part of a octo double;
 *   hilolo    second lowest part of a octo double;
 *   lololo    lowest part of a octo double;
 *   maxstp    the total number of iterations.
 *
 * ON RETURN :
 *   hihihi    highest part of the square root of a octo double;
 *   lohihi    second highest part of the square root of a octo double;
 *   lohihi    third highest part of the square root of a octo double;
 *   lolohi    fourth highest part of the square root of a octo double;
 *   hihilo    fourth lowest part of the square root of a octo double;
 *   lohilo    third lowest part of the square root of a octo double;
 *   hilolo    second lowest part of the square root of a octo double;
 *   lololo    lowest part of the square root of a octo double. */

#endif
