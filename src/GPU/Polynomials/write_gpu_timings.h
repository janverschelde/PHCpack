/* Utility to write the timings on the GPU. */

#ifndef __write_GPU_timings_h__
#define __write_GPU_timings_h__

void write_GPU_timings
 ( double cnvlapms, double addlapms, double elapsedms, double walltimesec );
/*
 * DESCRIPTION :
 *   Writes the timings.
 *
 * ON ENTRY :
 *   cnvlapms is the elapsed time spent by all convolution kernels,
 *            expressed in milliseconds;
 *   addlapms is the elapsed time spent by all addition kernels,
 *            expressed in milliseconds;
 *   elapsedms is the elapsed time spent by all kernels,
 *            expressed in milliseconds;
 *   walltimesec is the elapsed wall clock time for all computations
 *            (excluding memory copies) in seconds. */

#endif
