/* The file reducers.h contains prototypes to reduce systems.
 * By default, compilation with gcc is assumed.
 * To compile with a C++ compiler such as g++, the flag compilewgpp must
 * be defined as "g++ -Dcompilewgpp=1." */

#ifndef __REDUCERS_H__
#define __REDUCERS_H__

#ifdef compilewgpp
extern "C" void adainit( void );
extern "C" int _ada_use_c2phc ( int task, int *a, int *b, double *c );
extern "C" void adafinal( void );
#else
extern void adainit( void );
extern int _ada_use_c2phc ( int task, int *a, int *b, double *c );
extern void adafinal( void );
#endif

int standard_reduce_system ( int diag );
/*
 * DESCRIPTION :
 *   Applies row reduction to the coefficient matrix of the polynomial
 *   system in the standard systems container,
 *   with standard double precision arithmetic.
 *   The system in the standard systems container is replaced by
 *   the reduced system.
 *
 * REQUIRED :
 *   An equal number of equations as variables is assumed.
 *
 * ON ENTRY :
 *   diag     1 if the coefficient matrix needs to be diagonalized,
 *            0 if bringing in triangular form suffices. */

#endif
