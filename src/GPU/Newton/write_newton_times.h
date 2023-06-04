// The file write_newton_times.h specifies a function to write the
// times of running Newton's method.

#ifndef __write_newton_times_h__
#define __write_newton_times_h__

void write_newton_times
 ( int stepcnt, double walltimesec, double totcnvlapsedms,
   double totqrlapsedms, double totqtblapsedms, double totbslapsedms,
   double totupdlapsedms, double totreslapsedms );
/*
 * DESCRIPTION :
 *   Writes the times of running Newton's method.
 *
 * ON ENTRY :
 *   stepcnt    number of Newton steps done;
 *   walltimesec is the elapsed wall clock time in seconds;
 *   totcnvlapsedms is the convolutions time in milliseconds;
 *   totqrlapsedms is the qr decomposition time in milliseconds;
 *   totqtblapsedms is the q^T*b time in milliseconds;
 *   totbslapsedms is the back substitution time in milliseconds;
 *   totupdlapsedms is the update time in milliseconds;
 *   totreslapsedms is the residual time in milliseconds. */

#endif
