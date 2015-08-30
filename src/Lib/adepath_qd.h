/* The file adepath_qd.h contains prototypes to call Newton's method and
 * the path tracking methods using algorithmic differentiation in quad
 * double precision which use the Ada code for the data structures in PHCpack.
 * By default, compilation with gcc is assumed.
 * To compile with a C++ compiler such as g++, the flag compilewgpp must
 * be defined as "g++ -Dcompilewgpp=1." */

#ifndef __ADEPATH_DD_H__
#define __ADEPATH_DD_H__

int ade_newton_qd ( int verbose );
/*
 * DESCRIPTION :
 *   Runs Newton's method with algorithmic differentation in quad
 *   double precision on the data in the systems and solutions container.
 *
 * REQUIRED :
 *   The dobldobl systems container contains a valid polynomial system
 *   and the quaddobl solutions container holds a valid solution.
 *
 * ON ENTRY :
 *   verbose  0 if no intermediate output is wanted,
 *            1 if extra information should be written to screen.
 *
 * ON RETURN :
 *   fail     0 if all went well, and the corrected solution is in the
 *              solution container,
 *            if different from zero, then an error happened. */

int ade_onepath_qd ( int verbose, double regamma, double imgamma );
/*
 * DESCRIPTION :
 *   Tracks one solution path with algorithmic differentation in quad
 *   double precision on the data in the systems and solutions container.
 *
 * REQUIRED :
 *   The start and target systems have been defined
 *   and the quaddobl solutions container holds a valid solution.
 *
 * ON ENTRY :
 *   verbose  0 if no intermediate output is wanted,
 *            1 if extra information should be written to screen;
 *   regamma  real part of the random complex constant gamma;
 *   imgamma  imaginary part of the random complex constant gamma.
 *
 * ON RETURN :
 *   fail     0 if all went well, and the solution at the end of the path 
 *              is in the  solution container,
 *            if different from 0, then an error happened. */

int ade_manypaths_qd ( int verbose, double regamma, double imgamma );
/*
 * DESCRIPTION :
 *   Tracks many solution paths with algorithmic differentation in quad
 *   double precision on the data in the systems and solutions container.
 *
 * REQUIRED :
 *   The start and target systems have been defined
 *   and the quaddobl solutions container holds valid solutions.
 *
 * ON ENTRY :
 *   verbose  0 if no intermediate output is wanted,
 *            1 if extra information should be written to screen;
 *   regamma  real part of the random complex constant gamma;
 *   imgamma  imaginary part of the random complex constant gamma.
 *
 * ON RETURN :
 *   fail     0 if all went well, and the solutions at the end of paths
 *              are in the solution container,
 *            if different from 0, then an error happened. */

#endif
