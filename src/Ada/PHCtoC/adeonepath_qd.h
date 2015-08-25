/* This file adeonepath_qd.h collects the prototypes to track one path
 * using algorithmic differentiation in quad double precision */

#ifndef __ADEONEPATH_QD_H__
#define __ADEONEPATH_QD_H__

#include "poly.h"
#include "polysol.h"

int track ( int verbose, PolySys& p, PolySys& q, PolySolSet& s );
/*
 * DESCRIPTION :
 *   Tracks one path defined by an artificial parameter homotopy,
 *   starting at a solution s of q and ending at a solution of p.
 *   If verbose > 0, then additional output is written to screen. */

extern "C" int adeonepath_qd ( int verbose );
/*
 * DESCRIPTION :
 *   A C++ function to track one solution path,
 *   encapsulated as a C function for to be called from Ada.
 *   If verbose > 0, then additional output will be written. */

#endif
