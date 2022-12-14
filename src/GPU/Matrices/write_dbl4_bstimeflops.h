/* Utility to write the timings and the flop counts on the GPU
 * for the upper tiled back substitution solver in quad double precision. */

#ifndef __write_dbl4_bstimeflops_h__
#define __write_dbl4_bstimeflops_h__

void write_dbl4_bstimeflops
 ( int sizetile, int numtiles, int ctype,
   double invlapsed, double mullapsed, double sublapsed, double elapsedms,
   double timelapsed, long int addcnt, double addover,
   long int mulcnt, double mulover, long int divcnt, double divover );
/*
 * DESCRIPTION :
 *   Writes timings, operational counts, arithmetical intensity and flops.
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
 *   addover      overflowed number of additions and subtractions;
 *   mulcnt       number of multiplications;
 *   mulover      overflowed number of multiplications;
 *   divcnt       number of divisions;
 *   divover      overflowed number of divisions. */

#endif
