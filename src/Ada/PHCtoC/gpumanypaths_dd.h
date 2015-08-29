/* This file gpumanypaths_dd.h collects the prototypes to track many paths
 * with GPU acceleration in double double precision. */

#ifndef __GPUMANYPATHS_DD_H__
#define __GPUMANYPATHS_DD_H__

#include "poly.h"
#include "polysol.h"

int manytrack
 ( int mode, int verbose, double regamma, double imgamma,
   PolySys& p, PolySys& q, PolySolSet& s );
/*
 * DESCRIPTION :
 *   Tracks many paths defined by an artificial parameter homotopy,
 *   starting at solutions in s of q and ending at solutions of p.
 *
 * ON ENTRY :
 *   mode      execution mode, 0 (GPU+CPU), 1 (CPU), or 2 (GPU);
 *   verbose   if > 0, then additional output is written to screen;
 *   regamma   real part of the gamma constant;
 *   imgamma   imaginary part of the gamma constant;
 *   p         target system in the homotopy;
 *   q         start system in the homotopy;
 *   s         solutions of the start system q. */

extern "C" int gpumanypaths_dd
 ( int mode, int verbose, double regamma, double imgamma );
/*
 * DESCRIPTION :
 *   A C++ function to track one solution path,
 *   encapsulated as a C function for to be called from Ada.
 *
 * ON ENTRY :
 *   mode      execution mode for CPU and GPU equals 0, 1, or 2:
 *             if mode = 0, then both CPU and GPU will execute,
 *             if mode = 1, then only CPU runs Newton's method,
 *             if mode = 2, then only GPU runs Newton's method.
 *   verbose   if > 0, then additional output is written to screen;
 *   regamma   real part of the gamma constant;
 *   imgamma   imaginary part of the gamma constant. */

#endif
