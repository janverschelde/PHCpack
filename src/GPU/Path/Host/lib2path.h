/* This file contains the prototypes for the development of the C++ and C
 * interface to the Ada code in PHCpack. */

#ifndef _LIB2PATH_H_
#define _LIB2PATH_H_

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <string.h>

#include <qd/qd_real.h>

#include "syscon.h"
#include "solcon.h"
#include "phcpack.h"
#include "jump_track.h"

#include "complexH.h"
#include "poly.h"
#include "polysol.h"

#ifdef compilewgpp
extern "C" void adainit( void );
extern "C" void adafinal( void );
#else
extern void adainit( void );
extern void adafinal( void );
#endif

void var_name ( char* x_string, int x_string_len, string*& x_names, int& dim );
/*
 * DESCRIPTION :
 *   Auxiliary function to read the names of the variables in a system. */

void lib2path_read_standard_sys
 ( int verbose, PolySys<complexH<double>,double>& sys );
void lib2path_read_dobldobl_sys
 ( int verbose, PolySys<complexH<dd_real>,dd_real>& sys );
void lib2path_read_quaddobl_sys
 ( int verbose, PolySys<complexH<qd_real>,qd_real>& sys );
/*
 * DESCRIPTION :
 *   Reads a polynomial system from the systems container
 *   into the data structure provided by the parameter sys.
 *   If verbose > 0, then extra output is written to screen. */

void lib2path_read_standard_sols
 ( PolySys<complexH<double>,double>& start_sys,
   PolySolSet<complexH<double>,double>& sols );
void lib2path_read_dobldobl_sols
 ( PolySys<complexH<dd_real>,dd_real>& start_sys,
   PolySolSet<complexH<dd_real>,dd_real>& sols );
void lib2path_read_quaddobl_sols
 ( PolySys<complexH<qd_real>,qd_real>& start_sys,
   PolySolSet<complexH<qd_real>,qd_real>& sols );
/* 
 * DESCRIPTION :
 *   Reads the corresponding solutions of the start system start_sys,
 *   into the provided data structure, respectively sols. 
 *   The dimension of the start_sys is used on input. */

void lib2path_write_standard_sols
 ( PolySolSet<complexH<double>,double>& sols );
void lib2path_write_dobldobl_sols
 ( PolySolSet<complexH<dd_real>,dd_real>& sols );
void lib2path_write_quaddobl_sols
 ( PolySolSet<complexH<qd_real>,qd_real>& sols );
/*
 * DESCRIPTION :
 *   Takes the solutions in sols and places them in the solutions container. */

void lib2path_read_standard_homotopy
 ( char* start_file, char* target_file,
   PolySys<complexH<double>,double>& start_sys,
   PolySys<complexH<double>,double>& target_sys,
   PolySolSet<complexH<double>,double>& sols );
void lib2path_read_dobldobl_homotopy
 ( char* start_file, char* target_file,
   PolySys<complexH<dd_real>,dd_real>& start_sys,
   PolySys<complexH<dd_real>,dd_real>& target_sys,
   PolySolSet<complexH<dd_real>,dd_real>& sols );
void lib2path_read_quaddobl_homotopy
 ( char* start_file, char* target_file,
   PolySys<complexH<qd_real>,qd_real>& start_sys,
   PolySys<complexH<qd_real>,qd_real>& target_sys,
   PolySolSet<complexH<qd_real>,qd_real>& sols );
/*
 * DESCRIPTION :
 *   Given in the strings the names of the files for start and target system,
 *   the C interface to PHCpack is used to parse the systems,
 *   and the solutions of the start system. */


/* This file ademanypaths.h collects the prototypes to track many paths
 * using algorithmic differentiation in double, double double,
 * and quad double precision. */

int standard_manytrack
 ( int verbose, double regamma, double imgamma,
   PolySys<complexH<double>,double>& p,
   PolySys<complexH<double>,double>& q,
   PolySolSet<complexH<double>,double>& s );
int dobldobl_manytrack
 ( int verbose, double regamma, double imgamma,
   PolySys<complexH<dd_real>,dd_real>& p,
   PolySys<complexH<dd_real>,dd_real>& q,
   PolySolSet<complexH<dd_real>,dd_real>& s );
int quaddobl_manytrack
 ( int verbose, double regamma, double imgamma,
   PolySys<complexH<qd_real>,qd_real>& p,
   PolySys<complexH<qd_real>,qd_real>& q,
   PolySolSet<complexH<qd_real>,qd_real>& s );
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

extern "C" int standard_ademanypaths
 ( int verbose, double regamma, double imgamma );
extern "C" int dobldobl_ademanypaths
 ( int verbose, double regamma, double imgamma );
extern "C" int quaddobl_ademanypaths
 ( int verbose, double regamma, double imgamma );
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
