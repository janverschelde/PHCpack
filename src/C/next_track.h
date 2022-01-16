/* The file next_track.h contains prototypes for the operations to
 * track a solution path with a generator, i.e.: a get_next() method.
 * By default, compilation with gcc is assumed.
 * To compile with a C++ compiler such as g++, the flag compilewgpp must
 * be defined as "g++ -Dcompilewgpp=1." */

#ifndef __NEXT_TRACK_H__
#define __NEXT_TRACK_H__

#ifdef compilewgpp
extern "C" void adainit( void );
extern "C" int _ada_use_c2phc4c ( int task, int *a, int *b, double *c, int v );
extern "C" void adafinal( void );
#else
extern void adainit( void );
extern int _ada_use_c2phc4c ( int task, int *a, int *b, double *c, int v );
extern void adafinal( void );
#endif

int initialize_standard_homotopy
 ( int fixed_gamma, double regamma, double imgamma );
/*
 * DESCRIPTION :
 *   Takes start and target system as stored in standard double precision
 *   in the PHCpack containers and initializes the homotopy for tracking
 *   in standard double complex arithmetic. 
 *   If fixed_gamma equals 1 (true), then gamma will be a fixed default,
 *   otherwise, a new random complex constant for gamma is generated,
 *   but only if regamma and imgamma are both zero.
 *   If fixed_gamma is not 1 and regamma or imgamma is nonzero,
 *   then the value of the complex number defined by regamma and imgamma
 *   will define the gamma constant in the homotopy. */

int initialize_dobldobl_homotopy
 ( int fixed_gamma, double regamma, double imgamma );
/*
 * DESCRIPTION :
 *   Takes start and target system as stored in double double precision
 *   in the PHCpack containers and initializes the homotopy for tracking
 *   in double double complex arithmetic.
 *   If fixed_gamma equals 1 (true), then gamma will be a fixed default,
 *   otherwise, a new random complex constant for gamma is generated,
 *   but only if regamma and imgamma are both zero.
 *   If fixed_gamma is not 1 and regamma or imgamma is nonzero,
 *   then the value of the complex number defined by regamma and imgamma
 *   will define the gamma constant in the homotopy. */

int initialize_quaddobl_homotopy
 ( int fixed_gamma, double regamma, double imgamma );
/*
 * DESCRIPTION :
 *   Takes start and target system as stored in quad double precision
 *   in the PHCpack containers and initializes the homotopy for tracking
 *   in quad double complex arithmetic.
 *   If fixed_gamma equals 1 (true), then gamma will be a fixed default,
 *   otherwise, a new random complex constant for gamma is generated,
 *   but only if regamma and imgamma are both zero.
 *   If fixed_gamma is not 1 and regamma or imgamma is nonzero,
 *   then the value of the complex number defined by regamma and imgamma
 *   will define the gamma constant in the homotopy. */

int initialize_multprec_homotopy ( int fixed_gamma, int decimals );
/*
 * DESCRIPTION :
 *   Takes start and target system as stored in multiprecision
 *   in the PHCpack containers and initializes the homotopy for tracking
 *   in multiprecision complex arithmetic with in the working precision
 *   as many decimal places as the value of decimals.
 *   If fixed_gamma equals 1 (true), then gamma will be a fixed default,
 *   otherwise, a new random complex constant for gamma is generated. */

int initialize_varbprec_homotopy 
 ( int fixed_gamma, int nc_target, char *target, int nc_start, char *start );
/*
 * DESCRIPTION :
 *   Initializes the variable precision homotopy with the target and
 *   start system stored in the strings.
 *
 * ON ENTRY :
 *   fixed_gamma   if 1, then a fixed value for the gamma constant is used,
 *                 if 0, a random value for gamma will be generated;
 *   nc_target     number of characters in the string target;
 *   target        string representation of the target system;
 *   nc_start      number of characters in the string start;
 *   start         string representation of the start system. */

int initialize_standard_solution ( int k );
/*
 * DESCRIPTION :
 *   Takes the k-th solution in the standard solution container
 *   and initializes the standard double path tracker with generator. */

int initialize_dobldobl_solution ( int k );
/*
 * DESCRIPTION :
 *   Takes the k-th solution in the double double solution container
 *   and initializes the double double path tracker with generator. */

int initialize_quaddobl_solution ( int k );
/*
 * DESCRIPTION :
 *   Takes the k-th solution in the quad double solution container
 *   and initializes the double double path tracker with generator. */

int initialize_multprec_solution ( int k );
/*
 * DESCRIPTION :
 *   Takes the k-th solution in the multiprecision solution container
 *   and initializes the multiprecision path tracker with generator. */

int initialize_varbprec_solution ( int nv, int nc, char *sol );
/*
 * DESCRIPTION :
 *   Uses the string representation of a solution to initialize the
 *   variable precision path tracker with.
 *
 * ON ENTRY :
 *   nv      the number of variables in the solution;
 *   nc      the number of characters in the string sol;
 *   sol     string representation of a solution. */

int next_standard_solution ( int k );
/*
 * DESCRIPTION :
 *   Applies one predictor-corrector step to the initialized path tracker
 *   in standard double precision and replaces the k-th solution in the 
 *   standard solution container with a new solution on the path.
 *
 * REQUIRED :
 *   The standard homotopy has been initialized and the standard path
 *   tracker was initialized with the k-th start solution.  */

int next_dobldobl_solution ( int k );
/*
 * DESCRIPTION :
 *   Applies one predictor-corrector step to the initialized path tracker
 *   in double double precision and replaces the k-th solution in the 
 *   double double solution container with a new solution on the path.
 *
 * REQUIRED :
 *   The double double homotopy has been initialized and the double double
 *   path tracker was initialized with the k-th start solution.  */

int next_quaddobl_solution ( int k );
/*
 * DESCRIPTION :
 *   Applies one predictor-corrector step to the initialized path tracker
 *   in quad double precision and replaces the k-th solution in the 
 *   quad double solution container with a new solution on the path.
 *
 * REQUIRED :
 *   The quad double homotopy has been initialized and the quad double
 *   path tracker was initialized with the k-th start solution.  */

int next_multprec_solution ( int k );
/*
 * DESCRIPTION :
 *   Applies one predictor-corrector step to the initialized path tracker
 *   in multiple precision and replaces the k-th solution in the 
 *   multiprecision solution container with a new solution on the path.
 *
 * REQUIRED :
 *   The multiprecision homotopy has been initialized and the multiprecision
 *   path tracker was initialized with the k-th start solution.  */

char *next_varbprec_solution
 ( int want, int maxprc, int maxitr, int verbose, int *nc, int *fail ); 
/*
 * DESCRIPTION :
 *   Returns the next point along the path computed with variable precision.
 *
 * REQUIRED :
 *   The variable precision path tracker is initialized with a homotopy
 *   and a solution.
 *
 * ON ENTRY :
 *   want     number of decimal places wanted as accurate in the solution;
 *   maxprc   maximum number of decimal places in the working precision;
 *   maxitr   maximum number of corrector steps;
 *   verbose  if 1, then extra output is written to screen,
 *            if 0, then the next step is done in silence.
 *
 * ON RETURN :
 *   nc       number of characters allocated in the returned string,
 *            which represents the next solution along a path;
 *   fail     if 0, then no failure occurred. */

int clear_standard_tracker ( void );
/*
 * DESCRIPTION :
 *   Deallocates and resets data for tracking paths with a generator
 *   in standard double precision complex arithmetic. */

int clear_dobldobl_tracker ( void );
/*
 * DESCRIPTION :
 *   Deallocates and resets data for tracking paths with a generator
 *   in double double precision complex arithmetic. */

int clear_quaddobl_tracker ( void );
/*
 * DESCRIPTION :
 *   Deallocates and resets data for tracking paths with a generator
 *   in quad double precision complex arithmetic. */

int clear_multprec_tracker ( void );
/*
 * DESCRIPTION :
 *   Deallocates and resets data for tracking paths with a generator
 *   in multiprecision complex arithmetic. */

int clear_varbprec_tracker ( void );
/*
 * DESCRIPTION :
 *   Deallocates and resets data for tracking paths with a generator
 *   in variable precision complex arithmetic. */

#endif
