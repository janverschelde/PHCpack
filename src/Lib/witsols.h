/* The file witsols.h contains prototypes of functions
 * for a numerical irreducible decomposition.
 * By default, compilation with gcc is assumed.
 * To compile with a C++ compiler such as g++, the flag compilewgpp must
 * be defined as "g++ -Dcompilewgpp=1." */

#ifndef __WITSOLS_H__
#define __WITSOLS_H__

#ifdef compilewgpp
extern "C" void adainit( void );
extern "C" int _ada_use_c2phc4c ( int task, int *a, int *b, double *c, int v );
extern "C" void adafinal( void );
#else
extern void adainit( void );
extern int _ada_use_c2phc4c ( int task, int *a, int *b, double *c, int v );
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

/* extracting solution data */

int copy_standard_polysys_witset ( int dim );
/*
 * DESCRIPTION :
 *   Copies the witness set representation for a solution set
 *   of dimension dim into the systems and solutions container,
 *   in standard double precision.
 *
 * REQUIRED :
 *   1) standard_polysys_solve was executed successfully, and
 *   2) dim is in the range 0..topdim. */ 

int copy_standard_laursys_witset ( int dim );
/*
 * DESCRIPTION :
 *   Copies the witness set representation for a solution set
 *   of dimension dim into the Laurent systems and solutions container,
 *   in standard double precision.
 *
 * REQUIRED :
 *   1) standard_laursys_solve was executed successfully, and
 *   2) dim is in the range 0..topdim. */ 

int copy_dobldobl_polysys_witset ( int dim );
/*
 * DESCRIPTION :
 *   Copies the witness set representation for a solution set
 *   of dimension dim into the systems and solutions container,
 *   in double double precision.
 *
 * REQUIRED :
 *   1) dobldobl_polysys_solve was executed successfully, and
 *   2) dim is in the range 0..topdim. */ 

int copy_dobldobl_laursys_witset ( int dim );
/*
 * DESCRIPTION :
 *   Copies the witness set representation for a solution set
 *   of dimension dim into the Laurent systems and solutions container,
 *   in double double precision.
 *
 * REQUIRED :
 *   1) dobldobl_laursys_solve was executed successfully, and
 *   2) dim is in the range 0..topdim. */ 

int copy_quaddobl_polysys_witset ( int dim );
/*
 * DESCRIPTION :
 *   Copies the witness set representation for a solution set
 *   of dimension dim into the systems and solutions container,
 *   in quad double precision.
 *
 * REQUIRED :
 *   1) quaddobl_polysys_solve was executed successfully, and
 *   2) dim is in the range 0..topdim. */ 

int copy_quaddobl_laursys_witset ( int dim );
/*
 * DESCRIPTION :
 *   Copies the witness set representation for a solution set
 *   of dimension dim into the Laurent systems and solutions container,
 *   in quad double precision.
 *
 * REQUIRED :
 *   1) quaddobl_laursys_solve was executed successfully, and
 *   2) dim is in the range 0..topdim. */ 

int clear_standard_witsols ( void );
/*
 * DESCRIPTION :
 *   Clears the witness solutions in standard double precision. */

int clear_dobldobl_witsols ( void );
/*
 * DESCRIPTION :
 *   Clears the witness solutions in double double precision. */

int clear_quaddobl_witsols ( void );
/*
 * DESCRIPTION :
 *   Clears the witness solutions in quad double precision. */

#endif
