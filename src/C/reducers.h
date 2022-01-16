/* The file reducers.h contains prototypes to reduce systems.
 * By default, compilation with gcc is assumed.
 * To compile with a C++ compiler such as g++, the flag compilewgpp must
 * be defined as "g++ -Dcompilewgpp=1." */

#ifndef __REDUCERS_H__
#define __REDUCERS_H__

#ifdef compilewgpp
extern "C" void adainit( void );
extern "C" int _ada_use_c2phc4c ( int task, int *a, int *b, double *c, int v );
extern "C" void adafinal( void );
#else
extern void adainit( void );
extern int _ada_use_c2phc4c ( int task, int *a, int *b, double *c, int v );
extern void adafinal( void );
#endif

int standard_row_reduce_system ( int diag );
/*
 * DESCRIPTION :
 *   Applies row reduction to the coefficient matrix of the polynomial
 *   system in the standard double systems container,
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

int dobldobl_row_reduce_system ( int diag );
/*
 * DESCRIPTION :
 *   Applies row reduction to the coefficient matrix of the polynomial
 *   system in the double double systems container,
 *   with double double precision arithmetic.
 *   The system in the double double systems container is replaced by
 *   the reduced system.
 *
 * REQUIRED :
 *   An equal number of equations as variables is assumed.
 *
 * ON ENTRY :
 *   diag     1 if the coefficient matrix needs to be diagonalized,
 *            0 if bringing in triangular form suffices. */

int quaddobl_row_reduce_system ( int diag );
/*
 * DESCRIPTION :
 *   Applies row reduction to the coefficient matrix of the polynomial
 *   system in the quad double systems container,
 *   with quad double precision arithmetic.
 *   The system in the quad double systems container is replaced by
 *   the reduced system.
 *
 * REQUIRED :
 *   An equal number of equations as variables is assumed.
 *
 * ON ENTRY :
 *   diag     1 if the coefficient matrix needs to be diagonalized,
 *            0 if bringing in triangular form suffices. */

int standard_nonlinear_reduce_system
 ( int eqmax, int spmax, int rpmax, int *eqcnt, int *spcnt, int *rpcnt );
/*
 * DESCRIPTION :
 *   Applies nonlinear reduction to the system in the standard container.
 *
 * REQUIRED :
 *   The system is square, as many equations as variables is assumed.
 *
 * ON ENTRY :
 *   eqmax    the maximal number of equal degree replacements;
 *   spmax    the maximal number of computed S-polynomials;
 *   rpmax    the maximal number of computed R-polynomials.
 *
 * ON RETURN :
 *   eqcnt    the number of equal degree replacements;
 *   spcnt    the number of computed S-polynomials;
 *   rpcnt    the number of computed R-polynomials. */

#endif
