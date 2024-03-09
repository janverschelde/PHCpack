// The file dbl16_sqrt_kernels.h defines the prototypes for the square root
// in hexa double precision.

#ifndef __DBL16_SQRT_KERNELS_H__
#define __DBL16_SQRT_KERNELS_H__

__global__ void dbl16_sqrt
 ( double *hihihihi, double *lohihihi, double *hilohihi, double *lolohihi,
   double *hihilohi, double *lohilohi, double *hilolohi, double *lololohi,
   double *hihihilo, double *lohihilo, double *hilohilo, double *lolohilo,
   double *hihilolo, double *lohilolo, double *hilololo, double *lolololo,
   int max_steps );
/*
 * DESCRIPTION :
 *   Applies Newton's method to compute the square root of the
 *   hexa double number given by the sixteen doubles.
 *   Executes as many steps as the value of max_steps.
 *   Returns the value of the square root. */

void GPU_dbl16_sqrt
 ( double *hihihihi, double *lohihihi, double *hilohihi, double *lolohihi,
   double *hihilohi, double *lohilohi, double *hilolohi, double *lololohi,
   double *hihihilo, double *lohihilo, double *hilohilo, double *lolohilo,
   double *hihilolo, double *lohilolo, double *hilololo, double *lolololo,
   int maxstp );
/*
 * DESCRIPTION :
 *   Calls Newton's method to compute the square root of a hexa double.
 *
 * ON ENTRY :
 *   hihihihi    highest part of a hexa double;
 *   lohihihi    second highest part of a hexa double;
 *   hilohihi    third highest part of a hexa double;
 *   lohihihi    fourth highest part of a hexa double;
 *   hihilohi    fifth highest part of a hexa double;
 *   lohilohi    sixth highest part of a hexa double;
 *   hilolohi    seventh highest part of a hexa double;
 *   lololohi    eighth highest part of a hexa double;
 *   hihihilo    eighth lowest part of a hexa double;
 *   lohihilo    seventh lowest part of a hexa double;
 *   hilohilo    sixth lowest part of a hexa double;
 *   lolohilo    fifth lowest part of a hexa double;
 *   hihilolo    fourth lowest part of a hexa double;
 *   lohilolo    third lowest part of a hexa double;
 *   hilololo    second lowest part of a hexa double;
 *   lolololo    lowest part of a hexa double;
 *   maxstp      the total number of iterations.
 *
 * ON RETURN :
 *   hihihihi    highest part of the square root of a hexa double;
 *   lohihihi    second highest part of the square root of a hexa double;
 *   lohihihi    third highest part of the square root of a hexa double;
 *   lolohihi    fourth highest part of the square root of a hexa double;
 *   hihilohi    fifth lowest part of the square root of a hexa double;
 *   lohilohi    sixth lowest part of the square root of a hexa double;
 *   hilolohi    seventh lowest part of the square root of a hexa double;
 *   lololohi    eighth lowest part of the square root of a hexa double;
 *   hihihilo    eighth lowest part of the square root of a hexa double;
 *   lohihilo    seventh lowest part of the square root of a hexa double;
 *   hilohilo    sixth lowest part of the square root of a hexa double;
 *   lolohilo    fifth lowest part of the square root of a hexa double;
 *   hihilolo    fourth lowest part of the square root of a hexa double;
 *   lohilolo    third lowest part of the square root of a hexa double;
 *   hilololo    second lowest part of the square root of a hexa double;
 *   lolololo    lowest part of the square root of a hexa double. */

#endif
