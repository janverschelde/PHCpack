/* This file gpunewton_d.h collects the prototypes to run Newton's method
 * using algorithmic differentiation in double precision on the GPU. */

#ifndef __GPUNEWTON_D_H__
#define __GPUNEWTON_D_H__

#include "poly.h"
#include "polysol.h"

int newton ( int mode, int verbose, PolySys& p, PolySolSet& s );
/*
 * DESCRIPTION :
 *   Applies Newton's method to the first solution in s,
 *   on the polynomial system p.
 *   The execution mode is either 0 (CPU+GPU), 1 (CPU), or 2 (GPU).
 *   If verbose > 0, then additional output is written to screen. */

extern "C" int gpunewton_d ( int mode, int verbose );
/*
 * DESCRIPTION :
 *   A C++ function to accelerate Newton's method,
 *   encapsulated as a C function.
 *   If mode = 0, then both CPU and GPU will execute,
 *   if mode = 1, then only CPU runs Newton's method,
 *   if mode = 2, then only GPU runs Newton's method.
 *   If verbose > 0, then additional output will be written. */

#endif
