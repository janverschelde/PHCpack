/* The file gpupath_dd.h contains prototypes to call Newton's method and
 * the path tracking methods using algorithmic differentiation on the GPU
 * in double double precision wrapped so it can be used as a C library.
 * To compile with a C++ compiler such as g++, the flag compilewgpp must
 * be defined as "g++ -Dcompilewgpp=1." */

#ifndef __GPUPATH_DD_H__
#define __GPUPATH_DD_H__

int gpu_newton_dd ( int mode, int verbose );
/*
 * DESCRIPTION :
 *   Runs Newton's method with algorithmic differentation on CPU and GPU
 *   in double double precision on the data in the systems and solutions
 *   containers.
 *
 * REQUIRED :
 *   The dobldobl systems container contains a valid polynomial system
 *   and the dobldobl solutions container holds a valid solution.
 *
 * ON ENTRY :
 *   mode     execution mode equals 0, 1, or 2:
 *            if mode = 0, then both CPU and GPU will execute,
 *            if mode = 1, then only CPU runs Newton's method,
 *            if mode = 2, then only GPU runs Newton's method;
 *   verbose  0 if no intermediate output is wanted,
 *            1 if extra information should be written to screen.
 *
 * ON RETURN :
 *   fail     0 if all went well, and the corrected solution is in the
 *              solution container,
 *            if different from zero, then an error happened. */

int gpu_onepath_dd ( int mode, int verbose, double regamma, double imgamma );
/*
 * DESCRIPTION :
 *   Tracks one solution path with algorithmic differentation on CPU and GPU
 *   in double double precision on the data in the systems and solutions
 *   containers.
 *
 * REQUIRED :
 *   The start and target systems have been defined
 *   and the dobldobl solutions container holds a valid solution.
 *
 * ON ENTRY :
 *   mode     execution mode equals 0, 1, or 2:
 *            if mode = 0, then both CPU and GPU will execute,
 *            if mode = 1, then only CPU runs Newton's method,
 *            if mode = 2, then only GPU runs Newton's method;
 *   verbose  0 if no intermediate output is wanted,
 *            1 if extra information should be written to screen;
 *   regamma  real part of the random gamma constant;
 *   imgamma  imaginary part of the random gamma constant.
 *
 * ON RETURN :
 *   fail     0 if all went well, and the solution at the end of the path 
 *              is in the  solution container,
 *            if different from 0, then an error happened. */

int gpu_manypaths_dd ( int mode, int verbose, double regamma, double imgamma );
/*
 * DESCRIPTION :
 *   Tracks many solution paths with algorithmic differentation on CPU and GPU
 *   in double double precision on the data in the systems and solutions 
 *   containers.
 *
 * REQUIRED :
 *   The start and target systems have been defined
 *   and the dobldobl solutions container holds valid solutions.
 *
 * ON ENTRY :
 *   mode     execution mode equals 0, 1, or 2:
 *            if mode = 0, then both CPU and GPU will execute,
 *            if mode = 1, then only CPU runs Newton's method,
 *            if mode = 2, then only GPU runs Newton's method;
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
