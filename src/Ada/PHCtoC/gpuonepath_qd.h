/* This file gpuonepath_qd.h collects the prototypes to track one path
 * with GPU acceleration in quad double precision. */

#ifndef __GPUONEPATH_QD_H__
#define __GPUONEPATH_QD_H__

#include "poly.h"
#include "polysol.h"

using namespace std;

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

extern "C" int gpuonepath_qd
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
