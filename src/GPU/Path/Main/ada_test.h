/* This file contains the prototypes for the development of the C++ and C
 * interface to the Ada code in PHCpack. */

#ifndef ADA_TEST_H_
#define ADA_TEST_H_

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <string.h>

#include "syscon.h"
#include "solcon.h"
#include "phcpack.h"
#include "jump_track.h"

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

void ada_read_sys ( int verbose, PolySys& sys );
/*
 * DESCRIPTION :
 *   Reads a polynomial system from the systems container
 *   into the data structure provided by the parameter sys.
 *   If verbose > 0, then extra output is written to screen. */

void ada_read_sols ( PolySys& start_sys, PolySolSet& sols );
/* 
 * DESCRIPTION :
 *   Reads the corresponding solutions of the start system start_sys,
 *   into the provided data structure, respectively sols. 
 *   The dimension of the start_sys is used on input. */

void ada_write_sols ( PolySolSet& sols );
/*
 * DESCRIPTION :
 *   Takes the solutions in sols and places them in the solutions container. */

void ada_read_homotopy
 ( char* start_file, char* target_file,
   PolySys& start_sys, PolySys& target_sys, PolySolSet& sols );
/*
 * DESCRIPTION :
 *   Given in the strings the names of the files for start and target system,
 *   the C interface to PHCpack is used to parse the systems,
 *   and the solutions of the start system. */

#endif
