/* Utility to write the timings and the flop counts on the GPU
 * for the Householder QR in quad double precision. */

#ifndef __write_dbl4_qrtimeflops_h__
#define __write_dbl4_qrtimeflops_h__

void write_dbl4_qrtimeflops
 ( int ctype, int nrows, int ncols,
   double houselapsedms, double RTvlapsedms, double tileRlapsedms,
   double vb2Wlapsedms, double WYTlapsedms, double QWYTlapsedms,
   double Qaddlapsedms, double YWTlapsedms, double YWTClapsedms,
   double Raddlapsedms, double timelapsed,
   long long int addcnt, long long int mulcnt,
   long long int divcnt, long long int sqrtcnt );
/*
 * DESCRIPTION :
 *   Writes the timings, operational counts and computes the flops.
 *
 * ON ENTRY :
 *   ctype          0 if real, otherwise complex;
 *   nrows          number of rows in the matrix;
 *   ncols          number of columns in the matrix;
 *   houselapsedms  time in milliseconds by the Householder kernel;
 *   RTvlapsedms    time in milliseconds for beta*R^T*v;
 *   tileRlapsedms  time in milliseconds to reduce on tile;
 *   vb2Wlapsedms   time in milliseconds for the W matrix;
 *   WYTlapsedms    time in milliseconds for Y*W^T;
 *   QWYTlapsedms   time in milliseconds for Q*YWT;
 *   Qaddlapsedms   time in milliseconds for QWYT + Q;
 *   YWTClapsedms   time in milliseconds for YWT*C;
 *   Raddlapsedms   time in milliseconds for R + YWTC;
 *   timelapsed     total GPU time wall clock computation time in seconds;
 *   addcnt         number of additions and subtractions;
 *   mulcnt         number of multiplications;
 *   divcnt         number of divisions;
 *   sqrtcnt        number of square roots. */

#endif
