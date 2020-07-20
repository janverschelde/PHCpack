/* The file multiplicity.h contains prototypes to compute the
 * multiplicity structure of an isolated solution of a polynomial system.
 * By default, compilation with gcc is assumed.
 * To compile with a C++ compiler such as g++, the flag compilewgpp must
 * be defined as "g++ -Dcompilewgpp=1." */

#ifndef __MULTIPLICITY_H__
#define __MULTIPLICITY_H__

#ifdef compilewgpp
extern "C" void adainit( void );
extern "C" int _ada_use_c2phc4c ( int task, int *a, int *b, double *c, int v );
extern "C" void adafinal( void );
#else
extern void adainit( void );
extern int _ada_use_c2phc4c ( int task, int *a, int *b, double *c, int v );
extern void adafinal( void );
#endif

int standard_multiplicity_structure
 ( int order, int verbose, double tol, int *multiplicity, int *hilbert );
/*
 * DESCRIPTION :
 *   Computes the multiplicity structure of the solution stored in
 *   the standard solutions container for the polynomial system stored
 *   in the standard systems container, in standard double precision.
 *
 * REQUIRED :
 *   hilbert points to space for order+1 integers.
 *
 * ON ENTRY :
 *   order    maximal order of differentiation;
 *   verbose  1 if output is requested, 0 otherwise;
 *   tol      tolerance on the numerical rank.
 *
 * ON RETURN :
 *   multiplicity is the computed multiplicity;
 *   hilbert contains the values of the Hilbert function. */

int dobldobl_multiplicity_structure
 ( int order, int verbose, double tol, int *multiplicity, int *hilbert );
/*
 * DESCRIPTION :
 *   Computes the multiplicity structure of the solution stored in
 *   the dobldobl solutions container for the polynomial system stored
 *   in the dobldobl systems container, in double double precision.
 *
 * REQUIRED :
 *   hilbert points to space for order+1 integers.
 *
 * ON ENTRY :
 *   order    maximal order of differentiation;
 *   verbose  1 if output is requested, 0 otherwise;
 *   tol      tolerance on the numerical rank.
 *
 * ON RETURN :
 *   multiplicity is the computed multiplicity;
 *   hilbert contains the values of the Hilbert function. */

int quaddobl_multiplicity_structure
 ( int order, int verbose, double tol, int *multiplicity, int *hilbert );
/*
 * DESCRIPTION :
 *   Computes the multiplicity structure of the solution stored in
 *   the quaddobl solutions container for the polynomial system stored
 *   in the quaddobl systems container, in quad double precision.
 *
 * REQUIRED :
 *   hilbert points to space for order+1 integers.
 *
 * ON ENTRY :
 *   order    maximal order of differentiation;
 *   verbose  1 if output is requested, 0 otherwise;
 *   tol      tolerance on the numerical rank.
 *
 * ON RETURN :
 *   multiplicity is the computed multiplicity;
 *   hilbert contains the values of the Hilbert function. */

#endif
