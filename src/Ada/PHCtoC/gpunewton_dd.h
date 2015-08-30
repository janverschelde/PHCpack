/* This file gpunewton_dd.h collects the prototypes to run Newton's method
 * using algorithmic differentiation in double double precision on the GPU. */

#ifndef __GPUNEWTON_DD_H__
#define __GPUNEWTON_DD_H__

#include "poly.h"
#include "polysol.h"

int newton ( int mode, int verbose, PolySys& p, PolySolSet& s );
/*
 * DESCRIPTION :
 *   Applies Newton's method to the first solution in s,
 *   on the polynomial system p, and with execution mode
 *   equal to 0 (CPU+GPU), 1 (CPU), or 2 (GPU).
 *   If verbose > 0, then additional output is written to screen. */

extern "C" int gpunewton_dd ( int mode, int verbose );
/*
 * DESCRIPTION :
 *   A C++ function to accelerate Newton's method,
 *   encapsulated as a C function for to be called from Ada.
 *   If mode = 0, then both CPU and GPU will execute,
 *   if mode = 1, then only CPU runs Newton's method,
 *   if mode = 2, then only GPU runs Newton's method.
 *   If verbose > 0, then additional output will be written. */

#endif
