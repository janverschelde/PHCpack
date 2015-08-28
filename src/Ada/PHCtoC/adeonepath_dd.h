/* This file adeonepath_dd.h collects the prototypes to track one path
 * using algorithmic differentiation in double double precision */

#ifndef __ADEONEPATH_DD_H__
#define __ADEONEPATH_DD_H__

#include "poly.h"
#include "polysol.h"

int track
 ( int verbose, double regamma, double imgamma,
   PolySys& p, PolySys& q, PolySolSet& s );
/*
 * DESCRIPTION :
 *   Tracks one path defined by an artificial parameter homotopy,
 *   starting at a solution s of q and ending at a solution of p.
 *
 * ON ENTRY :
 *   verbose   if > 0, then additional output is written to screen;
 *   regamma   real part of the gamma constant;
 *   imgamma   imaginary part of the gamma constant;
 *   p         target system in the homotopy;
 *   q         start system in the homotopy;
 *   s         a solution of the start system q. */

extern "C" int adeonepath_dd ( int verbose, double regamma, double imgamma );
/*
 * DESCRIPTION :
 *   A C++ function to track one solution path,
 *   encapsulated as a C function for to be called from Ada.
 *
 * ON ENTRY :
 *   verbose   if > 0, then additional output is written to screen;
 *   regamma   real part of the gamma constant;
 *   imgamma   imaginary part of the gamma constant; */

#endif
