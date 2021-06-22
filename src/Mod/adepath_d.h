/* The file adepath_d.h contains prototypes to call Newton's method and
 * the path tracking methods using algorithmic differentiation in double
 * precision wrapped so it can be used as a C library.
 * To compile with a C++ compiler such as g++, the flag compilewgpp must
 * be defined as "g++ -Dcompilewgpp=1." */

#ifndef __ADEPATH_D_H__
#define __ADEPATH_D_H__

int ade_newton_d ( int verbose );
/*
 * DESCRIPTION :
 *   Runs Newton's method with algorithmic differentation
 *   in double precision on the data in the systems and solutions container.
 *
 * REQUIRED :
 *   The standard systems container contains a valid polynomial system
 *   and the standard solutions container holds a valid solution.
 *
 * ON ENTRY :
 *   verbose  0 if no intermediate output is wanted,
 *            1 if extra information should be written to screen;
 *
 * ON RETURN :
 *   fail     0 if all went well, and the corrected solution is in the
 *              solution container,
 *            if different from zero, then an error happened. */

int ade_onepath_d ( int verbose, double regamma, double imgamma );
/*
 * DESCRIPTION :
 *   Tracks one solution path with algorithmic differentation
 *   in double precision on the data in the systems and solutions container.
 *
 * REQUIRED :
 *   The start and target systems have been defined
 *   and the standard solutions container holds a valid solution.
 *
 * ON ENTRY :
 *   verbose  0 if no intermediate output is wanted,
 *            1 if extra information should be written to screen;
 *
 * ON RETURN :
 *   fail     0 if all went well, and the solution at the end of the path 
 *              is in the  solution container,
 *            if different from 0, then an error happened. */

int ade_manypaths_d ( int verbose, double regamma, double imgamma );
/*
 * DESCRIPTION :
 *   Tracks many solution paths with algorithmic differentation
 *   in double precision on the data in the systems and solutions container.
 *
 * REQUIRED :
 *   The start and target systems have been defined
 *   and the standard solutions container holds valid solutions.
 *
 * ON ENTRY :
 *   verbose  0 if no intermediate output is wanted,
 *            1 if extra information should be written to screen;
 *   regamma  real part of the random gamma constant;
 *   imgamma  imaginary part of the random gamma constant.
 *
 * ON RETURN :
 *   fail     0 if all went well, and the solutions at the end of paths
 *              are in the solution container,
 *            if different from 0, then an error happened. */

#endif
