/* The file sweep.h contains prototypes to the parameter homotopy continuation.
 * By default, compilation with gcc is assumed.
 * To compile with a C++ compiler such as g++, the flag compilewgpp must
 * be defined as "g++ -Dcompilewgpp=1." */

#ifndef __SWEEP_H__
#define __SWEEP_H__

#ifdef compilewgpp
extern "C" void adainit( void );
extern "C" int _ada_use_c2phc ( int task, int *a, int *b, double *c );
extern "C" void adafinal( void );
#else
extern void adainit( void );
extern int _ada_use_c2phc ( int task, int *a, int *b, double *c );
extern void adafinal( void );
#endif

int sweep_define_parameters_numerically ( int nq, int nv, int np, int *pars );
/*
 * DESCRIPTION :
 *   Defines the indices to the variables that serve as parameters
 *   numerically, that is: via integer indices.
 *
 * ON ENTRY :
 *   nq      the number of equations;
 *   nv      the number of variables, which includes the number of parameters;
 *   np      the number of parameters;
 *   pars    indices to the parameter variables, should be np in number. */

int sweep_get_number_of_equations ( int *n );
/*
 * DESCRIPTION :
 *   Returns in n the number of equations. */

int sweep_get_number_of_variables ( int *n );
/*
 * DESCRIPTION :
 *   Returns in n the number of variables. */

int sweep_get_number_of_parameters ( int *n );
/*
 * DESCRIPTION :
 *   Returns in n the number of parameters. */

int sweep_get_indices_numerically ( int *idx );
/*
 * DESCRIPTION :
 *   Returns in idx the indices of the variables that are parameters,
 *   if the argument contains sufficiently allocated space. */

int sweep_clear_definitions ( void );
/*
 * DESCRIPTION :
 *   Clears the definitions of the parameters. */

#endif
