// C++ function called by Ada routine to accelerate the tracking of one path,
// this file gpuonepath_dd.h documents the prototypes of the functions.

#ifndef __GPUONEPATH_DD_H__
#define __GPUONEPATH_DD_H__

#include "poly.h"
#include "polysol.h"

int track
 ( int mode, int verbose, double regamma, double imgamma,
   PolySys& p, PolySys& q, PolySolSet& s );
/*
 * DESCRIPTION :
 *   Tracks one path defined by an artificial parameter homotopy,
 *   starting at a solution s of q and ending at a solution of p.
 *
 * ON ENTRY :
 *   mode      execution mode is 0 (CPU + GPU), 1 (CPU), or 2 (GPU);
 *   verbose   if > 0, then additional output is written to screen;
 *   regamma   real part of the gamma constant;
 *   imgamma   imaginary part of the gamma constant;
 *   p         target system in the homotopy;
 *   q         start system in the homotopy;
 *   s         a solution of the start system q. */

extern "C" int gpuonepath_dd
 ( int mode, int verbose, double regamma, double imgamma );
/*
 * DESCRIPTION :
 *   A C++ function to track one solution path,
 *   encapsulated as a C function for to be called from Ada.
 *
 * ON ENTRY :
 *   mode      execution mode equals 0, 1, or 2:
 *             if mode = 0, then both CPU and GPU will execute,
 *             if mode = 1, then only CPU runs Newton's method,
 *             if mode = 2, then only GPU runs Newton's method;
 *   verbose   if > 0, then additional output is written to screen;
 *   regamma   real part of the gamma constant;
 *   imgamma   imaginary part of the gamma constant. */

#endif
