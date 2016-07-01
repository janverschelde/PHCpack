/* The file series.h contains prototypes to the function to compute power
 * series solutions of polynomial systems by Newton's method.
 * By default, compilation with gcc is assumed.
 * To compile with a C++ compiler such as g++, the flag compilewgpp must
 * be defined as "g++ -Dcompilewgpp=1." */

#ifndef __SERIES_H__
#define __SERIES_H__

#ifdef compilewgpp
extern "C" void adainit( void );
extern "C" int _ada_use_c2phc ( int task, int *a, int *b, double *c );
extern "C" void adafinal( void );
#else
extern void adainit( void );
extern int _ada_use_c2phc ( int task, int *a, int *b, double *c );
extern void adafinal( void );
#endif

int standard_Newton_series ( int idx, int nbr, int verbose );
/*
 * DESCRIPTION :
 *   Given in the systems container a polynomial system with coefficients
 *   in standard double precision, and in the solutions container the
 *   leading coefficients of the power series, this function runs Newton's
 *   method to compute power series solutions of the system in the container,
 *   in standard double precision.
 *
 * ON ENTRY :
 *   idx     index of the series parameter;
 *   nbr     number of Newton steps to be done on each solution;
 *   verbose is 0 or 1 to indicate whether additional diagnostic output needs
 *           to be written to screen.
 *
 * ON RETURN :
 *   The failure code, which equals zero if all went well. */

int dobldobl_Newton_series ( int idx, int nbr, int verbose );
/*
 * DESCRIPTION :
 *   Given in the systems container a polynomial system with coefficients
 *   in standard double precision, and in the solutions container the
 *   leading coefficients of the power series, this function runs Newton's
 *   method to compute power series solutions of the system in the container,
 *   in double double precision.
 *
 * ON ENTRY :
 *   idx     index of the series parameter;
 *   nbr     number of Newton steps to be done on each solution;
 *   verbose is 0 or 1 to indicate whether additional diagnostic output needs
 *           to be written to screen.
 *
 * ON RETURN :
 *   The failure code, which equals zero if all went well. */

int quaddobl_Newton_series ( int idx, int nbr, int verbose );
/*
 * DESCRIPTION :
 *   Given in the systems container a polynomial system with coefficients
 *   in standard double precision, and in the solutions container the
 *   leading coefficients of the power series, this function runs Newton's
 *   method to compute power series solutions of the system in the container,
 *   in quad double precision.
 *
 * ON ENTRY :
 *   idx     index of the series parameter;
 *   nbr     number of Newton steps to be done on each solution;
 *   verbose is 0 or 1 to indicate whether additional diagnostic output needs
 *           to be written to screen.
 *
 * ON RETURN :
 *   The failure code, which equals zero if all went well. */

#endif
