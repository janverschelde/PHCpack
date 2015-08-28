/* This file ademanypaths_dd.h collects the prototypes to track many paths
 * using algorithmic differentiation in double double precision. */

#ifndef __ADEMANYPATHS_DD_H__
#define __ADEMANYPATHS_DD_H__

#include "poly.h"
#include "polysol.h"

int manytrack
 ( int verbose, double regamma, double imgamma,
   PolySys& p, PolySys& q, PolySolSet& s );
/*
 * DESCRIPTION :
 *   Tracks many paths defined by an artificial parameter homotopy,
 *   starting at solutions in s of q and ending at solutions of p.
 *
 * ON ENTRY :
 *   verbose   if > 0, then additional output is written to screen;
 *   regamma   real part of the gamma constant;
 *   imgamma   imaginary part of the gamma constant;
 *   p         target system in the homotopy;
 *   q         start system in the homotopy;
 *   s         solutions of the start system q. */

extern "C" int ademanypaths_dd ( int verbose, double regamma, double imgamma );
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
