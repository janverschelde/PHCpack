/* This file adenewton_qd.h collects the prototypes to run Newton's method
 * using algorithmic differentiation in quad double precision */

#ifndef __ADENEWTON_QD_H__
#define __ADENEWTON_QD_H__

#include "poly.h"
#include "polysol.h"

int newton ( int verbose, PolySys& p, PolySolSet& s );
/*
 * DESCRIPTION :
 *   Applies Newton's method to the first solution in s,
 *   on the polynomial system p.
 *   If verbose > 0, then additional output is written to screen. */

extern "C" int adenewton_qd ( int verbose );
/*
 * DESCRIPTION :
 *   A C++ function to accelerate Newton's method,
 *   encapsulated as a C function for to be called from Ada.
 *   If verbose > 0, then additional output will be written. */

#endif
