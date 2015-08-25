/* This file adeonepath_dd.h collects the prototypes to track one path
 * using algorithmic differentiation in double double precision */

#ifndef __ADEONEPATH_DD_H__
#define __ADEONEPATH_DD_H__

#include "poly.h"
#include "polysol.h"

int track ( int verbose, PolySys& p, PolySys& q, PolySolSet& s );
/*
 * DESCRIPTION :
 *   Tracks one path defined by an artificial parameter homotopy,
 *   starting at a solution s of q and ending at a solution of p.
 *   If verbose > 0, then additional output is written to screen. */

extern "C" int adeonepath_dd ( int verbose );
/*
 * DESCRIPTION :
 *   A C++ function to track one solution path,
 *   encapsulated as a C function for to be called from Ada.
 *   If verbose > 0, then additional output will be written. */

#endif
