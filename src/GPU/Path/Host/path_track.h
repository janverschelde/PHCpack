#ifndef _PATH_TRACK_H_
#define _PATH_TRACK_H_

#include "poly.h"
#include "polysol.h"
#include "workspace_host.h"

template <class ComplexType, class RealType>
int manytrack
 ( int verbose, double regamma, double imgamma, int prec,
   PolySys<ComplexType,RealType>& p, PolySys<ComplexType,RealType>& q,
   PolySolSet<ComplexType,RealType>& s );
/*
 * DESCRIPTION :
 *
 * ON ENTRY :
 *   verbose   0 if silent, 1 for intermediate output;
 *   regamma   real part of the gamma constant;
 *   imgamma   imaginary part of the gamma constant;
 *   prec      precision, 16, 32, or 64;
 *   p         target system;
 *   q         start system;
 *   s         start solutions. */

#include "path_track.tpp"

#endif
