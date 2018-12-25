/* The file padcon.h contains prototypes to the path trackers which apply
 * Pade approximants to predict solutions, padcon = Pade continuation.
 * By default, compilation with gcc is assumed.
 * To compile with a C++ compiler such as g++,
 * the flag compilewgpp must be defined as "g++ -Dcompilewgpp=1." */

#ifndef __PADCON_H__
#define __PADCON_H__

#ifdef compilewgpp
extern "C" void adainit( void );
extern "C" int _ada_use_c2phc4c ( int task, int *a, int *b, double *c );
extern "C" void adafinal( void );
#else
extern void adainit( void );
extern int _ada_use_c2phc4c ( int task, int *a, int *b, double *c );
extern void adafinal( void );
#endif

int padcon_set_default_parameters ( void );
/*
 * DESCRIPTION :
 *   Sets the default values of the homotopy continuation parameters. */

int padcon_clear_parameters ( void );
/*
 * DESCRIPTION :
 *   Deallocates the allocated space for the parameters. */

int padcon_get_homotopy_continuation_parameter ( int k, double *val );
/*
 * DESCRIPTION :
 *   Returns in val the value of the k-th continuation parameter,
 *   if k ranges between 1 and 12.
 *
 * ON ENTRY :
 *   k        an integer number between 1 and 12.
 *
 * ON RETURN 
 *   val      the value for the k-th homotopy continuation parameter. */

int padcon_set_homotopy_continuation_parameter ( int k, double *val );
/*
 * DESCRIPTION :
 *   Sets the value of the k-th continuation parameter to val,
 *   if k ranges between 1 and 12.
 *
 * ON ENTRY :
 *   k        an integer number between 1 and 12.
 *
 * ON RETURN 
 *   val      the value for the k-th homotopy continuation parameter. */

int padcon_standard_track ( int nbc, char* name, int verbose );
/*
 * DESCRIPTION :
 *   For the defined target, start system, and start solutions,
 *   launches the Pade continuation in standard double precision.
 *
 * ON ENTRY :
 *   nbc      number of characters of the output file name,
 *            equals 0 if no output will be written to file;
 *   name     defines the name of the output file,
 *            the file name has nbc characters;
 *   verbose  if > 0, then more information is written. */

int padcon_dobldobl_track ( int nbc, char* name, int verbose );
/*
 * DESCRIPTION :
 *   For the defined target, start system, and start solutions,
 *   launches the Pade continuation in double double precision.
 *
 * ON ENTRY :
 *   nbc      number of characters of the output file name,
 *            equals 0 if no output will be written to file;
 *   name     defines the name of the output file,
 *            the file name has nbc characters;
 *   verbose  if > 0, then more information is written. */

int padcon_quaddobl_track ( int nbc, char* name, int verbose );
/*
 * DESCRIPTION :
 *   For the defined target, start system, and start solutions,
 *   launches the Pade continuation in quad double precision.
 *
 * ON ENTRY :
 *   nbc      number of characters of the output file name,
 *            equals 0 if no output will be written to file;
 *   name     defines the name of the output file,
 *            the file name has nbc characters;
 *   verbose  if > 0, then more information is written. */

int padcon_standard_initialize_homotopy ( int verbose );
/*
 * DESCRIPTION :
 *   For the defined target and start system,
 *   initializes the homotopy in standard double precision,
 *   for the step-by-step Pade continuation.
 *   If verbose = 1, then extra output will be written. */

int padcon_dobldobl_initialize_homotopy ( int verbose );
/*
 * DESCRIPTION :
 *   For the defined target and start system,
 *   initializes the homotopy in double double precision,
 *   for the step-by-step Pade continuation.
 *   If verbose = 1, then extra output will be written. */

int padcon_quaddobl_initialize_homotopy ( int verbose );
/*
 * DESCRIPTION :
 *   For the defined target and start system,
 *   initializes the homotopy in quad double precision,
 *   for the step-by-step Pade continuation.
 *   If verbose = 1, then extra output will be written. */

int padcon_initialize_standard_solution ( int idx, int verbose );
/*
 * DESCRIPTION :
 *   Takes the solution with index idx in the solutions container in
 *   standard double precision and initializes the series-Pade tracker.
 *   If verbose = 1, then extra output will be written. */

int padcon_initialize_dobldobl_solution ( int idx, int verbose );
/*
 * DESCRIPTION :
 *   Takes the solution with index idx in the solutions container in
 *   double double precision and initializes the series-Pade tracker.
 *   If verbose = 1, then extra output will be written. */

int padcon_initialize_quaddobl_solution ( int idx, int verbose );
/*
 * DESCRIPTION :
 *   Takes the solution with index idx in the solutions container in
 *   quad double precision and initializes the series-Pade tracker.
 *   If verbose = 1, then extra output will be written. */

int padcon_standard_predict_correct ( int* fail, int verbose );
/*
 * DESCRIPTION :
 *   Executes one predict-correct step on the current solution and
 *   the defined homotopy in standard double precision.
 *   On return in fail is the failure code of the predict-correct step:
 *   if fail is zero, then the required accuracies were met,
 *   otherwise, either the predict or the correct step failed.
 *   If verbose = 1, then extra output will be written. */

int padcon_dobldobl_predict_correct ( int* fail, int verbose );
/*
 * DESCRIPTION :
 *   Executes one predict-correct step on the current solution and
 *   the defined homotopy in double double precision.
 *   On return in fail is the failure code of the predict-correct step:
 *   if fail is zero, then the required accuracies were met,
 *   otherwise, either the predict or the correct step failed.
 *   If verbose = 1, then extra output will be written. */

int padcon_quaddobl_predict_correct ( int* fail, int verbose );
/*
 * DESCRIPTION :
 *   Executes one predict-correct step on the current solution and
 *   the defined homotopy in quad double precision.
 *   On return in fail is the failure code of the predict-correct step:
 *   if fail is zero, then the required accuracies were met,
 *   otherwise, either the predict or the correct step failed.
 *   If verbose = 1, then extra output will be written. */

int padcon_get_standard_solution ( int idx, int verbose );
/*
 * DESCRIPTION :
 *   Retrieves the current solution and places it at position idx
 *   in the solutions container in standard double precision.
 *   If verbose = 1, then extra output will be written. */

int padcon_get_dobldobl_solution ( int idx, int verbose );
/*
 * DESCRIPTION :
 *   Retrieves the current solution and places it at position idx
 *   in the solutions container in double double precision.
 *   If verbose = 1, then extra output will be written. */

int padcon_get_quaddobl_solution ( int idx, int verbose );
/*
 * DESCRIPTION :
 *   Retrieves the current solution and places it at position idx
 *   in the solutions container in quad double precision.
 *   If verbose = 1, then extra output will be written. */

int padcon_clear_standard_data ( void );
/*
 * DESCRIPTION :
 *   Deallocates data for the series-Pade tracker in double precision. */

int padcon_clear_dobldobl_data ( void );
/*
 * DESCRIPTION :
 *   Deallocates data for the series-Pade tracker
 *   in double double precision. */

int padcon_clear_quaddobl_data ( void );
/*
 * DESCRIPTION :
 *   Deallocates data for the series-Pade tracker
 *   in quad double precision. */

#endif
