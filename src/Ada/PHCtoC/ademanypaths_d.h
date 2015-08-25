/* This file ademanypaths_d.h collects the prototypes to track many paths
 * using algorithmic differentiation in double precision */

#ifndef __ADEMANYPATHS_D_H__
#define __ADEMANYPATHS_D_H__

#include "poly.h"
#include "polysol.h"

int manytrack ( int verbose, PolySys& p, PolySys& q, PolySolSet& s );
/*
 * DESCRIPTION :
 *   Tracks many paths defined by an artificial parameter homotopy,
 *   starting at solutions in s of q and ending at solutions of p.
 *   If verbose > 0, then additional output is written to screen. */

extern "C" int ademanypaths_d ( int verbose );
/*
 * DESCRIPTION :
 *   A C++ function to track one solution path,
 *   encapsulated as a C function for to be called from Ada.
 *   If verbose > 0, then additional output will be written. */

#endif
