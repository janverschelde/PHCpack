/* The file sweep.h contains prototypes to the parameter homotopy continuation.
 * By default, compilation with gcc is assumed.
 * To compile with a C++ compiler such as g++, the flag compilewgpp must
 * be defined as "g++ -Dcompilewgpp=1." */

#ifndef __SWEEP_H__
#define __SWEEP_H__

#ifdef compilewgpp
extern "C" void adainit( void );
extern "C" int _ada_use_c2phc4c ( int task, int *a, int *b, double *c, int v );
extern "C" void adafinal( void );
#else
extern void adainit( void );
extern int _ada_use_c2phc4c ( int task, int *a, int *b, double *c, int v );
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

int sweep_define_parameters_symbolically
 ( int nq, int nv, int np, int nc, char *pars );
/*
 * DESCRIPTION :
 *   Defines the indices to the variables that serve as parameters
 *   symbolically, that is, as names of variables.
 *   For this to work, the symbol table must be initialized.
 *
 * ON ENTRY :
 *   nq      the number of equations;
 *   nv      the number of variables, which includes the number of parameters;
 *   np      the number of parameters;
 *   nc      the number of characters in the string pars;
 *   pars    names of the parameter variables, should be np in number,
 *           the names are separated by one space each. */

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

int sweep_get_indices_symbolically ( int *nc, char *pars );
/*
 * DESCRIPTION :
 *   Returns in the string pars with nc characters the names
 *   of the parameters, each separated by one space. */

int sweep_clear_definitions ( void );
/*
 * DESCRIPTION :
 *   Clears the definitions of the parameters. */

int sweep_set_standard_start ( int n, double *c );
/*
 * DESCRIPTION :
 *   Sets the start values for the parameters in standard double precision,
 *   giving on input n doubles in c, with the consecutive real and imaginary
 *   parts for the start values of all parameters. */

int sweep_set_standard_target ( int n, double *c );
/*
 * DESCRIPTION :
 *   Sets the target values for the parameters in standard double precision,
 *   giving on input n doubles in c, with the consecutive real and imaginary
 *   parts for the target values of all parameters. */

int sweep_set_dobldobl_start ( int n, double *c );
/*
 * DESCRIPTION :
 *   Sets the start values for the parameters in double double precision,
 *   giving on input n doubles in c, with the consecutive real and imaginary
 *   parts for the start values of all parameters. */

int sweep_set_dobldobl_target ( int n, double *c );
/*
 * DESCRIPTION :
 *   Sets the target values for the parameters in double double precision,
 *   giving on input n doubles in c, with the consecutive real and imaginary
 *   parts for the target values of all parameters. */

int sweep_set_quaddobl_start ( int n, double *c );
/*
 * DESCRIPTION :
 *   Sets the start values for the parameters in quad double precision,
 *   giving on input n doubles in c, with the consecutive real and imaginary
 *   parts for the start values of all parameters. */

int sweep_set_quaddobl_target ( int n, double *c );
/*
 * DESCRIPTION :
 *   Sets the target values for the parameters in quad double precision,
 *   giving on input n doubles in c, with the consecutive real and imaginary
 *   parts for the target values of all parameters. */

int sweep_get_standard_start ( int n, double *c );
/*
 * DESCRIPTION :
 *   Gets the start values for the parameters in standard double precision,
 *   giving on input in n the number doubles that need to be returned in c.
 *   On return in c will be n doubles, for the consecutive real and imaginary
 *   parts for the start values of all parameters. */

int sweep_get_standard_target ( int n, double *c );
/*
 * DESCRIPTION :
 *   Gets the target values for the parameters in standard double precision,
 *   giving on input in n the number doubles that need to be returned in c.
 *   On return in c will be n doubles, for the consecutive real and imaginary
 *   parts for the target values of all parameters. */

int sweep_get_dobldobl_start ( int n, double *c );
/*
 * DESCRIPTION :
 *   Gets the start values for the parameters in double double precision,
 *   giving on input in n the number doubles that need to be returned in c.
 *   On return in c will be n doubles, for the consecutive real and imaginary
 *   parts for the start values of all parameters. */

int sweep_get_dobldobl_target ( int n, double *c );
/*
 * DESCRIPTION :
 *   Gets the target values for the parameters in double double precision,
 *   giving on input in n the number doubles that need to be returned in c.
 *   On return in c will be n doubles, for the consecutive real and imaginary
 *   parts for the target values of all parameters. */

int sweep_get_quaddobl_start ( int n, double *c );
/*
 * DESCRIPTION :
 *   Gets the start values for the parameters in quad double precision,
 *   giving on input in n the number doubles that need to be returned in c.
 *   On return in c will be n doubles, for the consecutive real and imaginary
 *   parts for the start values of all parameters. */

int sweep_get_quaddobl_target ( int n, double *c );
/*
 * DESCRIPTION :
 *   Gets the target values for the parameters in quad double precision,
 *   giving on input in n the number doubles that need to be returned in c.
 *   On return in c will be n doubles, for the consecutive real and imaginary
 *   parts for the target values of all parameters. */

int sweep_standard_complex_run
 ( int gchoice, double *regamma, double *imgamma );
/*
 * DESCRIPTION :
 *   Starts the trackers in a complex convex parameter homotopy,
 *   in standard double precision, where the indices to the parameters,
 *   start and target values are already defined.  Moreover, the containers
 *   of systems and solutions in standard double precision have been
 *   initialized with a parametric systems and start solutions.
 *   The input parameter gchoice is 0, 1, or 2, for respectively
 *   a randomly generated gamma (0), or no gamma (1), or a user given
 *   gamma with real and imaginary parts in regamma and imgamma. */

int sweep_dobldobl_complex_run
 ( int gchoice, double *regamma, double *imgamma );
/*
 * DESCRIPTION :
 *   Starts the trackers in a complex convex parameter homotopy,
 *   in double double precision, where the indices to the parameters,
 *   start and target values are already defined.  Moreover, the containers
 *   of systems and solutions in double double precision have been
 *   initialized with a parametric systems and start solutions.
 *   The input parameter gchoice is 0, 1, or 2, for respectively
 *   a randomly generated gamma (0), or no gamma (1), or a user given
 *   gamma with real and imaginary parts in regamma and imgamma.
 *   For gchoice == 2, regamma and imgamma are double doubles. */

int sweep_quaddobl_complex_run
 ( int gchoice, double *regamma, double *imgamma );
/*
 * DESCRIPTION :
 *   Starts the trackers in a complex convex parameter homotopy,
 *   in quad double precision, where the indices to the parameters,
 *   start and target values are already defined.  Moreover, the containers
 *   of systems and solutions in quad double precision have been
 *   initialized with a parametric systems and start solutions.
 *   The input parameter gchoice is 0, 1, or 2, for respectively
 *   a randomly generated gamma (0), or no gamma (1), or a user given
 *   gamma with real and imaginary parts in regamma and imgamma.
 *   For gchoice == 2, regamma and imgamma are quad doubles. */

int sweep_standard_real_run ( void );
/*
 * DESCRIPTION :
 *   Starts a sweep with a natural parameter in a family of n equations
 *   in n+1 variables, where the last variable is the artificial parameter s
 *   that moves the one natural parameter from a start to target value.
 *   The last equation is of the form (1-s)*(A - v[0]) + s*(A - v[1]),
 *   where A is the natural parameter, going from the start value v[0]
 *   to the target value v[1].
 *   This family must be stored in the systems container in standard double
 *   precision and the corresponding start solutions in the standard solutions
 *   container, where every solution has the value v[0] for the A variable.
 *   The sweep stops when s reaches the value v[1], or when a singularity
 *   is encountered on the path. */

int sweep_dobldobl_real_run ( void );
/*
 * DESCRIPTION :
 *   Starts a sweep with a natural parameter in a family of n equations
 *   in n+1 variables, where the last variable is the artificial parameter s
 *   that moves the one natural parameter from a start to target value.
 *   The last equation is of the form (1-s)*(A - v[0]) + s*(A - v[1]),
 *   where A is the natural parameter, going from the start value v[0]
 *   to the target value v[1].
 *   This family must be stored in the systems container in double double
 *   precision and the corresponding start solutions in the dobldobl solutions
 *   container, where every solution has the value v[0] for the A variable.
 *   The sweep stops when s reaches the value v[1], or when a singularity
 *   is encountered on the path. */

int sweep_quaddobl_real_run ( void );
/*
 * DESCRIPTION :
 *   Starts a sweep with a natural parameter in a family of n equations
 *   in n+1 variables, where the last variable is the artificial parameter s
 *   that moves the one natural parameter from a start to target value.
 *   The last equation is of the form (1-s)*(A - v[0]) + s*(A - v[1]),
 *   where A is the natural parameter, going from the start value v[0]
 *   to the target value v[1].
 *   This family must be stored in the systems container in quad double
 *   precision and the corresponding start solutions in the quaddobl solutions
 *   container, where every solution has the value v[0] for the A variable.
 *   The sweep stops when s reaches the value v[1], or when a singularity
 *   is encountered on the path. */

#endif
