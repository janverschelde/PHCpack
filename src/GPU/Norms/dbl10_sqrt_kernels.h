// The file dbl10_sqrt_kernels.h defines the prototypes for the square root
// in deca double precision.

#ifndef __DBL10_SQRT_KERNELS_H__
#define __DBL10_SQRT_KERNELS_H__

__global__ void dbl10_sqrt
 ( double *rtb, double *rix, double *rmi, double *rrg, double *rpk,
   double *ltb, double *lix, double *lmi, double *lrg, double *lpk,
   int max_steps );
/*
 * DESCRIPTION :
 *   Applies Newton's method to compute the square root of the
 *   deca double number given by the ten parts, respectively in 
 *   rtb, rix, rmi, rrg, rpk, ltb, lix, lmi, lrg, and lpk.
 *   Executes as many steps as the value of max_steps.
 *   Returns the value of the square root. */

void GPU_dbl10_sqrt
 ( double *rtb, double *rix, double *rmi, double *rrg, double *rpk,
   double *ltb, double *lix, double *lmi, double *lrg, double *lpk,
   int maxstp );
/*
 * DESCRIPTION :
 *   Calls Newton's method to compute the square root of a deca double.
 *
 * ON ENTRY :
 *   rtb      highest part of a deca double;
 *   rix      second highest part of a deca double;
 *   rmi      third highest part of a deca double;
 *   rrg      fourth highest part of a deca double;
 *   rpk      fifth highest part of a deca double;
 *   ltb      fifth lowest part of a deca double;
 *   lix      fourth lowest part of a deca double;
 *   lmi      third lowest part of a deca double;
 *   lrg      second lowest part of a deca double;
 *   lpk      lowest part of a deca double;
 *   maxstp   the total number of iterations.
 *
 * ON RETURN :
 *   rtb      highest part of the square root of a deca double;
 *   rix      second highest part of the square root of a deca double;
 *   rmi      third highest part of the square root of a deca double;
 *   rrg      fourth highest part of the square root of a deca double;
 *   rpk      fifth highest part of the square root of a deca double;
 *   ltb      fifth lowest part of the square root of a deca double;
 *   lix      fourth lowest part of the square root of a deca double;
 *   lmi      third lowest part of the square root of a deca double;
 *   lrg      second lowest part of the square root of a deca double;
 *   lpk      lowest part of the square root of a deca double. */

#endif
