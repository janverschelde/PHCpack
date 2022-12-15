/* Utility to write the timings and the flop counts on the GPU
 * for the upper tiled back substitution solver in double precision. */

#ifndef __write_dbl_bstimeflops_h__
#define __write_dbl_bstimeflops_h__

void write_dbl_bstimeflops
 ( int sizetile, int numtiles, int ctype,
   double invlapsed, double mullapsed, double sublapsed, double elapsedms,
   double timelapsed, long int addcnt, long int mulcnt, long int divcnt );
/*
 * DESCRIPTION :
 *   Writes the timings, operational counts and computes the flops.
 *
 * ON ENTRY :
 *   sizetile     size of each tile;
 *   numtiles     number of tiles;
 *   ctype        0 if real, otherwise complex;
 *   invlapsed    time in milliseconds to invert diagonal tiles;
 *   mullapsed    time in milliseconds to multiply with inverted tiles;
 *   sublapsed    time in milliseconds for back substitution;
 *   elapsedms    total time in milliseconds spent by all kernels;
 *   timelapsed   total GPU time wall clock computation time in seconds;
 *   addcnt       number of additions and subtractions;
 *   mulcnt       number of multiplications;
 *   divcnt       number of divisions. */

#endif
