// The file path_track.h contains the prototypes to track many paths,
// with or without multithreading.

#ifndef __PATH_TRACK_H__
#define __PATH_TRACK_H__

#include "polysys.h"
#include "polysolset.h"
#include "workspace_host.h"
#include "jobqueue.h"
#include "worker.h"

template <class ComplexType, class RealType>
int manytrack
 ( int verbose, double regamma, double imgamma, Parameter pars, int prec,
   PolySys<ComplexType,RealType>& p, PolySys<ComplexType,RealType>& q,
   PolySolSet<ComplexType,RealType>& s );
/*
 * DESCRIPTION :
 *   Calls the path tracker to track all paths starting from the solutions
 *   in s, solutions of the start system q, to solve the system p.
 *
 * ON ENTRY :
 *   verbose   0 if silent, 1 for intermediate output;
 *   regamma   real part of the gamma constant;
 *   imgamma   imaginary part of the gamma constant;
 *   pars      settings for the parameters and tolerances;
 *   prec      precision, 16, 32, or 64;
 *   p         target system;
 *   q         start system;
 *   s         start solutions. */

template <class ComplexType, class RealType>
int crewmanytrack
 ( int crewsize,
   int verbose, double regamma, double imgamma, Parameter pars, int prec,
   PolySys<ComplexType,RealType>& p, PolySys<ComplexType,RealType>& q,
   PolySolSet<ComplexType,RealType>& s );
/*
 * DESCRIPTION :
 *   Launches a work crew to track all paths starting from the solutions
 *   in s, solutions of the start system q, to solve the system p.
 *
 * ON ENTRY :
 *   crewsize  number of workers in the crew;
 *   verbose   0 if silent, 1 for intermediate output;
 *   regamma   real part of the gamma constant;
 *   imgamma   imaginary part of the gamma constant;
 *   pars      settings for the parameters and tolerances;
 *   prec      precision, 16, 32, or 64;
 *   p         target system;
 *   q         start system;
 *   s         start solutions. */

#include "path_track.tpp"

#endif /* __PATH_TRACK_H__ */
