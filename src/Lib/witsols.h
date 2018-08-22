/* The file witsols.h contains prototypes of functions
 * for a numerical irreducible decomposition.
 * By default, compilation with gcc is assumed.
 * To compile with a C++ compiler such as g++, the flag compilewgpp must
 * be defined as "g++ -Dcompilewgpp=1." */

#ifndef __WITSOLS_H__
#define __WITSOLS_H__

#ifdef compilewgpp
extern "C" void adainit( void );
extern "C" int _ada_use_c2phc4c ( int task, int *a, int *b, double *c );
extern "C" void adafinal( void );
#else
extern void adainit( void );
extern int _ada_use_c2phc4c ( int task, int *a, int *b, double *c );
extern void adafinal( void );
#endif

/* the solvers */

int standard_polysys_solve
 ( int nbtasks, int topdim, int filter, int factor, int verbose );
/*
 * DESCRIPTION :
 *   Runs the cascades of homotopies on the polynomial system in
 *   the standard systems container.  Runs in standard double precision.
 *
 * ON ENTRY :
 *   nbtasks   equals the number of tasks for multitasking,
 *   topdim    the top dimension to start the homotopy cascades,
 *   filter    0 or 1 flag to filter the witness supersets, 
 *   factor    0 or 1 flag to factor the witness sets,
 *   verbose   for intermediate output. */

int standard_laursys_solve
 ( int nbtasks, int topdim, int filter, int factor, int verbose );
/*
 * DESCRIPTION :
 *   Runs the cascades of homotopies on the Laurent polynomial system in
 *   the standard systems container.  Runs in standard double precision.
 *
 * ON ENTRY :
 *   nbtasks   equals the number of tasks for multitasking,
 *   topdim    the top dimension to start the homotopy cascades,
 *   filter    0 or 1 flag to filter the witness supersets, 
 *   factor    0 or 1 flag to factor the witness sets,
 *   verbose   for intermediate output. */

int dobldobl_polysys_solve
 ( int nbtasks, int topdim, int filter, int factor, int verbose );
/*
 * DESCRIPTION :
 *   Runs the cascades of homotopies on the polynomial system in
 *   the dobldobl systems container.  Runs in double double precision.
 *
 * ON ENTRY :
 *   nbtasks   equals the number of tasks for multitasking,
 *   topdim    the top dimension to start the homotopy cascades,
 *   filter    0 or 1 flag to filter the witness supersets, 
 *   factor    0 or 1 flag to factor the witness sets,
 *   verbose   for intermediate output. */

int dobldobl_laursys_solve
 ( int nbtasks, int topdim, int filter, int factor, int verbose );
/*
 * DESCRIPTION :
 *   Runs the cascades of homotopies on the Laurent polynomial system in
 *   the dobldobl systems container.  Runs in double double precision.
 *
 * ON ENTRY :
 *   nbtasks   equals the number of tasks for multitasking,
 *   topdim    the top dimension to start the homotopy cascades,
 *   filter    0 or 1 flag to filter the witness supersets, 
 *   factor    0 or 1 flag to factor the witness sets,
 *   verbose   for intermediate output. */

int quaddobl_polysys_solve
 ( int nbtasks, int topdim, int filter, int factor, int verbose );
/*
 * DESCRIPTION :
 *   Runs the cascades of homotopies on the polynomial system in
 *   the quaddobl systems container.  Runs in quad double precision.
 *
 * ON ENTRY :
 *   nbtasks   equals the number of tasks for multitasking,
 *   topdim    the top dimension to start the homotopy cascades,
 *   filter    0 or 1 flag to filter the witness supersets, 
 *   factor    0 or 1 flag to factor the witness sets,
 *   verbose   for intermediate output. */

int quaddobl_laursys_solve
 ( int nbtasks, int topdim, int filter, int factor, int verbose );
/*
 * DESCRIPTION :
 *   Runs the cascades of homotopies on the Laurent polynomial system in
 *   the quaddobl systems container.  Runs in quad double precision.
 *
 * ON ENTRY :
 *   nbtasks   equals the number of tasks for multitasking,
 *   topdim    the top dimension to start the homotopy cascades,
 *   filter    0 or 1 flag to filter the witness supersets, 
 *   factor    0 or 1 flag to factor the witness sets,
 *   verbose   for intermediate output. */

#endif
