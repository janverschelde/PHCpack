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

int standard_newton
 ( int verbose, PolySys<complexH<double>,double>& p,
   PolySolSet<complexH<double>,double>& s );
int dobldobl_newton
 ( int verbose, PolySys<complexH<dd_real>,dd_real>& p,
   PolySolSet<complexH<dd_real>,dd_real>& s );
int quaddobl_newton
 ( int verbose, PolySys<complexH<qd_real>,qd_real>& p,
   PolySolSet<complexH<qd_real>,qd_real>& s );
/*
 * DESCRIPTION :
 *   Applies Newton's method to the first solution in s,
 *   on the polynomial system p.
 *   If verbose > 0, then additional output is written to screen.
 *   Default values for the path parameters will be applied. */

int standard_newton_with_pars
 ( int verbose, Parameter pars,
   PolySys<complexH<double>,double>& p,
   PolySolSet<complexH<double>,double>& s );
int dobldobl_newton_with_pars
 ( int verbose, Parameter pars,
   PolySys<complexH<dd_real>,dd_real>& p,
   PolySolSet<complexH<dd_real>,dd_real>& s );
int quaddobl_newton_with_pars
 ( int verbose, Parameter pars,
   PolySys<complexH<qd_real>,qd_real>& p,
   PolySolSet<complexH<qd_real>,qd_real>& s );
/*
 * DESCRIPTION :
 *   Applies Newton's method to the first solution in s,
 *   on the polynomial system p.
 *   If verbose > 0, then additional output is written to screen.
 *   Tuned values for the path parameter can be given on input. */

int standard_onetrack
 ( int verbose, double regamma, double imgamma,
   PolySys<complexH<double>,double>& p,
   PolySys<complexH<double>,double>& q,
   PolySolSet<complexH<double>,double>& s );
int dobldobl_onetrack
 ( int verbose, double regamma, double imgamma,
   PolySys<complexH<dd_real>,dd_real>& p,
   PolySys<complexH<dd_real>,dd_real>& q,
   PolySolSet<complexH<dd_real>,dd_real>& s );
int quaddobl_onetrack
 ( int verbose, double regamma, double imgamma,
   PolySys<complexH<qd_real>,qd_real>& p,
   PolySys<complexH<qd_real>,qd_real>& q,
   PolySolSet<complexH<qd_real>,qd_real>& s );
/*
 * DESCRIPTION :
 *   Tracks one path defined by an artificial parameter homotopy,
 *   starting at a solution s of q and ending at a solution of p.
 *   Default values for the parameter are applied.
 *
 * ON ENTRY :
 *   verbose   if > 0, then additional output is written to screen;
 *   regamma   real part of the gamma constant;
 *   imgamma   imaginary part of the gamma constant;
 *   p         target system in the homotopy;
 *   q         start system in the homotopy;
 *   s         a solution of the start system q. */

int standard_onetrack_with_pars
 ( int verbose, double regamma, double imgamma, Parameter pars,
   PolySys<complexH<double>,double>& p,
   PolySys<complexH<double>,double>& q,
   PolySolSet<complexH<double>,double>& s );
int dobldobl_onetrack_with_pars
 ( int verbose, double regamma, double imgamma, Parameter pars,
   PolySys<complexH<dd_real>,dd_real>& p,
   PolySys<complexH<dd_real>,dd_real>& q,
   PolySolSet<complexH<dd_real>,dd_real>& s );
int quaddobl_onetrack_with_pars
 ( int verbose, double regamma, double imgamma, Parameter pars,
   PolySys<complexH<qd_real>,qd_real>& p,
   PolySys<complexH<qd_real>,qd_real>& q,
   PolySolSet<complexH<qd_real>,qd_real>& s );
/*
 * DESCRIPTION :
 *   Tracks one path defined by an artificial parameter homotopy,
 *   starting at a solution s of q and ending at a solution of p.
 *   Tuned values for the parameter can be given on input.
 *
 * ON ENTRY :
 *   verbose   if > 0, then additional output is written to screen;
 *   regamma   real part of the gamma constant;
 *   imgamma   imaginary part of the gamma constant;
 *   pars      values for the path parameters;
 *   p         target system in the homotopy;
 *   q         start system in the homotopy;
 *   s         a solution of the start system q. */

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
 *   Default values for the path parameters are applied.
 *
 * ON ENTRY :
 *   verbose   if > 0, then additional output is written to screen;
 *   regamma   real part of the gamma constant;
 *   imgamma   imaginary part of the gamma constant;
 *   p         target system in the homotopy;
 *   q         start system in the homotopy;
 *   s         solutions of the start system q. */

int standard_manytrack_with_pars
 ( int verbose, double regamma, double imgamma, Parameter pars,
   PolySys<complexH<double>,double>& p,
   PolySys<complexH<double>,double>& q,
   PolySolSet<complexH<double>,double>& s );
int dobldobl_manytrack_with_pars
 ( int verbose, double regamma, double imgamma, Parameter pars,
   PolySys<complexH<dd_real>,dd_real>& p,
   PolySys<complexH<dd_real>,dd_real>& q,
   PolySolSet<complexH<dd_real>,dd_real>& s );
int quaddobl_manytrack_with_pars
 ( int verbose, double regamma, double imgamma, Parameter pars,
   PolySys<complexH<qd_real>,qd_real>& p,
   PolySys<complexH<qd_real>,qd_real>& q,
   PolySolSet<complexH<qd_real>,qd_real>& s );
/*
 * DESCRIPTION :
 *   Tracks many paths defined by an artificial parameter homotopy,
 *   starting at solutions in s of q and ending at solutions of p.
 *   Tuned values for the parameters can be provided on input.
 *
 * ON ENTRY :
 *   verbose   if > 0, then additional output is written to screen;
 *   regamma   real part of the gamma constant;
 *   imgamma   imaginary part of the gamma constant;
 *   pars      values for the path parameters;
 *   p         target system in the homotopy;
 *   q         start system in the homotopy;
 *   s         solutions of the start system q. */

int standard_ade_newton ( int verbose );
int dobldobl_ade_newton ( int verbose );
int quaddobl_ade_newton ( int verbose );
/*
 * DESCRIPTION :
 *   Calls Newton's method with algorithmic differentiation.
 *   If verbose > 0, then additional output will be written.
 *   Default values for the path parameters are applied. */

int standard_ade_newton_with_pars ( int verbose, Parameter pars );
int dobldobl_ade_newton_with_pars ( int verbose, Parameter pars );
int quaddobl_ade_newton_with_pars ( int verbose, Parameter pars );
/*
 * DESCRIPTION :
 *   Calls Newton's method with algorithmic differentiation.
 *   If verbose > 0, then additional output will be written.
 *   Tuned values of the parameters can be given on input in pars. */

int standard_ade_onepath
 ( int verbose, double regamma, double imgamma );
int dobldobl_ade_onepath
 ( int verbose, double regamma, double imgamma );
int quaddobl_ade_onepath
 ( int verbose, double regamma, double imgamma );
/*
 * DESCRIPTION :
 *   A C++ function to track one solution path.
 *   Default values for the path parameters are applied.
 *
 * ON ENTRY :
 *   verbose   if > 0, then additional output is written to screen;
 *   regamma   real part of the gamma constant;
 *   imgamma   imaginary part of the gamma constant. */

int standard_ade_onepath_with_pars
 ( int verbose, double regamma, double imgamma, Parameter pars );
int dobldobl_ade_onepath_with_pars
 ( int verbose, double regamma, double imgamma, Parameter pars );
int quaddobl_ade_onepath_with_pars
 ( int verbose, double regamma, double imgamma, Parameter pars );
/*
 * DESCRIPTION :
 *   A C++ function to track one solution path.
 *   Tuned values for the path parameters can be given in pars.
 *
 * ON ENTRY :
 *   verbose   if > 0, then additional output is written to screen;
 *   regamma   real part of the gamma constant;
 *   imgamma   imaginary part of the gamma constant;
 *   pars      values for the path parameters. */

int standard_ade_manypaths
 ( int verbose, double regamma, double imgamma );
int dobldobl_ade_manypaths
 ( int verbose, double regamma, double imgamma );
int quaddobl_ade_manypaths
 ( int verbose, double regamma, double imgamma );
/*
 * DESCRIPTION :
 *   A C++ function to track one solution path,
 *   encapsulated as a C function for to be called from Ada.
 *   Default values for the path parameters are applied.
 *
 * ON ENTRY :
 *   verbose   if > 0, then additional output is written to screen;
 *   regamma   real part of the gamma constant;
 *   imgamma   imaginary part of the gamma constant. */

int standard_ade_manypaths_with_pars
 ( int verbose, double regamma, double imgamma, Parameter pars );
int dobldobl_ade_manypaths_with_pars
 ( int verbose, double regamma, double imgamma, Parameter pars );
int quaddobl_ade_manypaths_with_pars
 ( int verbose, double regamma, double imgamma, Parameter pars );
/*
 * DESCRIPTION :
 *   A C++ function to track one solution path,
 *   encapsulated as a C function for to be called from Ada.
 *   Tuned values for the path parameters can be given in pars.
 *
 * ON ENTRY :
 *   verbose   if > 0, then additional output is written to screen;
 *   regamma   real part of the gamma constant;
 *   imgamma   imaginary part of the gamma constant;
 *   pars      values for the path parameters. */

extern "C" int standard_adenewton ( int verbose );
extern "C" int dobldobl_adenewton ( int verbose );
extern "C" int quaddobl_adenewton ( int verbose );
/*
 * DESCRIPTION :
 *   Calls Newton's method with algorithmic differentiaton,
 *   encapsulated as a C function for to be called from Ada.
 *   Notice the lack of an underscore after the _ade in the functin names.
 *   If verbose > 0, then additional output will be written. */

extern "C" int standard_adeonepath
 ( int verbose, double regamma, double imgamma );
extern "C" int dobldobl_adeonepath
 ( int verbose, double regamma, double imgamma );
extern "C" int quaddobl_adeonepath
 ( int verbose, double regamma, double imgamma );
/*
 * DESCRIPTION :
 *   A C++ function to track one solution path,
 *   encapsulated as a C function for to be called from Ada.
 *   Notice the lack of an underscore after the _ade in the function names.
 *
 * ON ENTRY :
 *   verbose   if > 0, then additional output is written to screen;
 *   regamma   real part of the gamma constant;
 *   imgamma   imaginary part of the gamma constant; */

extern "C" int standard_ademanypaths
 ( int verbose, double regamma, double imgamma );
extern "C" int dobldobl_ademanypaths
 ( int verbose, double regamma, double imgamma );
extern "C" int quaddobl_ademanypaths
 ( int verbose, double regamma, double imgamma );
/*
 * DESCRIPTION :
 *   A C++ function to track many solution path,
 *   encapsulated as a C function for to be called from Ada.
 *   Notice the lack of an underscore after the _ade in the function names.
 *
 * ON ENTRY :
 *   verbose   if > 0, then additional output is written to screen;
 *   regamma   real part of the gamma constant;
 *   imgamma   imaginary part of the gamma constant; */

#endif
