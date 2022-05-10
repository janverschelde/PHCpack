/* The file padcon.h contains prototypes to the path trackers which apply
 * Pade approximants to predict solutions, padcon = Pade continuation.
 * By default, compilation with gcc is assumed.
 * To compile with a C++ compiler such as g++,
 * the flag compilewgpp must be defined as "g++ -Dcompilewgpp=1." */

#ifndef __PADCON_H__
#define __PADCON_H__

#ifdef compilewgpp
extern "C" void adainit( void );
extern "C" int _ada_use_c2phc4c ( int task, int *a, int *b, double *c, int v );
extern "C" void adafinal( void );
#else
extern void adainit( void );
extern int _ada_use_c2phc4c ( int task, int *a, int *b, double *c, int v );
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

void padcon_write_homotopy_continuation_parameters ( void );
/*
 * DESCRIPTION :
 *   Writes the current values of the homotopy continuation parameters.
 *   As a pure C routine which makes no direct access to the Ada code,
 *   but it is too useful to be omitted from the padcon library. */

int padcon_write_homotopy_continuation_parameters_to_defined_output ( void );
/*
 * DESCRIPTION :
 *   Writes the current values of the homotopy continuation parameters
 *   to the defined output file, or to screen if there is no defined
 *   output file set by define_output_file* in phcpack. */

void padcon_tune_homotopy_continuation_parameters ( void );
/*
 * DESCRIPTION :
 *   Interactive loop to tune the homotopy continuation parameters.
 *   As a pure C routine it makes no direct access to the Ada code,
 *   but it is too useful to be omitted from the padcon library. */

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

int padcon_reset_homotopy_continuation_parameters ( int prc );
/*
 * DESCRIPTION :
 *   Resets the values of the homotopy continuation parameters
 *   in the step-by-step path trackers, using the value of prc,
 *   for double, double double, or quad double,
 *   depending whether prc is respectively zero, one, or two. */

int padcon_prompt_for_homogenization ( void );
/*
 * DESCRIPTION :
 *   Asks the user if a projective transformation must be applied.
 *   Return 1 for homogeneous coordinates, return 0 otherwise. */

int padcon_prompt_for_multi_homogenization ( int nvr );
/*
 * DESCRIPTION :
 *   Displays the menu for the multi-homogenization, for a number of
 *   variables equal to nvr.  On return is a positive integer in the range
 *   of 0 to nvr, where 0 standard for affine, 1 for 1-homogenization,
 *   and m for m-homogenization. */

void padcon_define_partition ( int m, int nvr, int *idz );
/*
 * DESCRIPTION :
 *   Interactively prompts the user for the partition of the variables.
 *
 * ON ENTRY :
 *   m        number of sets in the partition;
 *   nvr      the number of variables, which equals the number of symbols
 *            in the table and the length of the partition in idz;
 *   idz      space allocated for nvr integers.
 *
 * ON RETURN :
 *   idz      index representation of the partition,
 *            idz[k] is a number in the range 1 to m,
 *            defines the set the k-th variable in the partition. */

int padcon_add_Z0 ( void );
/*
 * DESCRIPTION :
 *   Adds the symbol 'Z0' to denote the 1-homogeneous variable. */

int padcon_add_symbols ( int m );
/*
 * DESCRIPTION :
 *   Augments the symbol table with Z1, Z2, .., Zm.
 *   Returns the failure code of the syscon function. */

void padcon_standard_projective_transformation ( void );
/*
 * DESCRIPTION :
 *   Replaces the start solutions in the solutions container by
 *   their 1-homogeneously transformed versions.
 *   Transforms as well the start and target systems
 *   in standard double precision.  Adds "Z0" to the symbol table. */

void padcon_dobldobl_projective_transformation ( void );
/*
 * DESCRIPTION :
 *   Replaces the start solutions in the solutions container by
 *   their 1-homogeneously transformed versions.
 *   Transforms as well the start and target systems
 *   in double double precision.  Adds "Z0" to the symbol table. */

void padcon_quaddobl_projective_transformation ( void );
/*
 * DESCRIPTION :
 *   Replaces the start solutions in the solutions container by
 *   their 1-homogeneously transformed versions.
 *   Transforms as well the start and target systems
 *   in quad double precision.  Adds "Z0" to the symbol table. */

void padcon_standard_multi_projective_transformation
 ( int n, int m, int *idz );
/*
 * DESCRIPTION :
 *   Replaces the start solutions in the solutions container by
 *   their m-homogeneously transformed versions.
 *   The index representation for the partition of the set of n variables
 *   is defined by the n integers in idz.
 *   Transforms as well the start and target systems
 *   in standard double precision.  Augments the symbol table. */

void padcon_dobldobl_multi_projective_transformation
 ( int n, int m, int *idz );
/*
 * DESCRIPTION :
 *   Replaces the start solutions in the solutions container by
 *   their m-homogeneously transformed versions.
 *   The index representation for the partition of the set of n variables
 *   is defined by the n integers in idz.
 *   Transforms as well the start and target systems
 *   in double double precision.  Augments the symbol table. */

void padcon_quaddobl_multi_projective_transformation
 ( int n, int m, int *idz );
/*
 * DESCRIPTION :
 *   Replaces the start solutions in the solutions container by
 *   their m-homogeneously transformed versions.
 *   The index representation for the partition of the set of n variables
 *   is defined by the n integers in idz.
 *   Transforms as well the start and target systems
 *   in quad double precision.  Augments the symbol table. */

int padcon_standard_track
 ( int nbc, char *name, int locfile, int verbose, int mhom,
   int nvr, int *idz );
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
 *   locfile  a flag to indicate if the file is local to the function,
 *            if 0, then the file is the defined output file, which is
 *            created but left open, so output may still be written to it,
 *            if 1, only this function will write to file;
 *   verbose  if > 0, then more information is written;
 *   mhom     if 0, then tracking happens in affine coordinates,
 *            if 1, then 1-homogeneous coordinates are applied,
 *            if m > 1, then m-homogenization happens,
 *            with the partition defined in idz;
 *   nvr      the number of variables, only if mhom > 1;
 *   idz      the index representation of the partition, only if mhom > 1,
 *            as an array of as many integers as the value of nvr. */

int padcon_dobldobl_track
 ( int nbc, char *name, int locfile, int verbose, int mhom,
   int nvr, int *idz );
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
 *   locfile  a flag to indicate if the file is local to the function,
 *            if 0, then the file is the defined output file, which is
 *            created but left open, so output may still be written to it,
 *   verbose  if > 0, then more information is written;
 *   mhom     if 0, then tracking happens in affine coordinates,
 *            if 1, then 1-homogeneous coordinates are applied,
 *            if m > 1, then m-homogenization happens,
 *            with the partition defined in idz;
 *   nvr      the number of variables, only if mhom > 1;
 *   idz      the index representation of the partition, only if mhom > 1,
 *            as an array of as many integers as the value of nvr. */

int padcon_quaddobl_track
 ( int nbc, char *name, int locfile, int verbose, int mhom,
   int nvr, int *idz );
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
 *   locfile  a flag to indicate if the file is local to the function,
 *            if 0, then the file is the defined output file, which is
 *            created but left open, so output may still be written to it,
 *   verbose  if > 0, then more information is written;
 *   mhom     if 0, then tracking happens in affine coordinates,
 *            if 1, then 1-homogeneous coordinates are applied,
 *            if m > 1, then m-homogenization happens,
 *            with the partition defined in idz;
 *   nvr      the number of variables, only if mhom > 1;
 *   idz      the index representation of the partition, only if mhom > 1,
 *            as an array of as many integers as the value of nvr. */

int padcon_standard_initialize_homotopy ( int verbose, int homo );
/*
 * DESCRIPTION :
 *   For the defined target and start system,
 *   initializes the homotopy in standard double precision,
 *   for the step-by-step Pade continuation.
 *   If verbose = 1, then extra output will be written.
 *   If homo = 0, then tracking happens in affine coordinates,
 *   if homo = 1, then homogeneous coordinate transformations are applied. */

int padcon_dobldobl_initialize_homotopy ( int verbose, int homo );
/*
 * DESCRIPTION :
 *   For the defined target and start system,
 *   initializes the homotopy in double double precision,
 *   for the step-by-step Pade continuation.
 *   If verbose = 1, then extra output will be written.
 *   If homo = 0, then tracking happens in affine coordinates,
 *   if homo = 1, then homogeneous coordinate transformations are applied. */

int padcon_quaddobl_initialize_homotopy ( int verbose, int homo );
/*
 * DESCRIPTION :
 *   For the defined target and start system,
 *   initializes the homotopy in quad double precision,
 *   for the step-by-step Pade continuation.
 *   If verbose = 1, then extra output will be written.
 *   If homo = 0, then tracking happens in affine coordinates,
 *   if homo = 1, then homogeneous coordinate transformations are applied. */

int padcon_standard_initialize_parameter_homotopy ( int idx, int verbose );
/*
 * DESCRIPTION :
 *   With the system, defined as target system, and the index idx
 *   for the continuation parameter, initializes the homotopy in
 *   standard double precision for the step-by-step Pade continuation.
 *   If verbose = 1, then extra output will be written. */

int padcon_dobldobl_initialize_parameter_homotopy ( int idx, int verbose );
/*
 * DESCRIPTION :
 *   With the system, defined as target system, and the index idx
 *   for the continuation parameter, initializes the homotopy in
 *   double double precision for the step-by-step Pade continuation.
 *   If verbose = 1, then extra output will be written. */

int padcon_quaddobl_initialize_parameter_homotopy ( int idx, int verbose );
/*
 * DESCRIPTION :
 *   With the system, defined as target system, and the index idx
 *   for the continuation parameter, initializes the homotopy in
 *   quad double precision for the step-by-step Pade continuation.
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

int padcon_get_standard_pole_radius ( double* frp );
/*
 * DESCRIPTION :
 *   Returns in frp the smallest pole radius computed 
 *   by the predictor in standard double precision. */

int padcon_get_dobldobl_pole_radius ( double* frp );
/*
 * DESCRIPTION :
 *   Returns in frp the smallest pole radius computed
 *   by the predictor in double double precision.
 *   The returned frp is the high part of the double double number. */

int padcon_get_quaddobl_pole_radius ( double* frp );
/*
 * DESCRIPTION :
 *   Returns in frp the smallest pole radius computed
 *   by the predictor in quad double precision.
 *   The returned frp is the highest part of the quad double number. */

int padcon_get_standard_closest_pole ( double* cre, double* cim );
/*
 * DESCRIPTION :
 *   Returns the real and imaginary parts of the closest pole
 *   respectively in cre and cim, computed by the predictor
 *   in standard double precision.
 *   Results are meaningful only if cre >= 0.0. */

int padcon_get_dobldobl_closest_pole ( double* cre, double* cim );
/*
 * DESCRIPTION :
 *   Returns the real and imaginary parts of the closest pole
 *   respectively in cre and cim, computed by the predictor
 *   in double double precision.
 *   The cre and cim are the high parts of the double doubles.
 *   Results are meaningful only if cre >= 0.0.  */

int padcon_get_quaddobl_closest_pole ( double* cre, double* cim );
/*
 * DESCRIPTION :
 *   Returns the real and imaginary parts of the closest pole
 *   respectively in cre and cim, computed by the predictor
 *   in quad double precision.
 *   The cre and cim are the highest parts of the quad doubles.
 *   Results are meaningful only if cre >= 0.0. */

int padcon_get_standard_t_value ( double *tval );
/*
 * DESCRIPTION :
 *   Returns in tval the current t value of the path tracker
 *   which runs in standard double precision. */

int padcon_get_dobldobl_t_value ( double *tval );
/*
 * DESCRIPTION :
 *   Returns in tval the current t value of the path tracker
 *   which runs in double double precision. */

int padcon_get_quaddobl_t_value ( double *tval );
/*
 * DESCRIPTION :
 *   Returns in tval the current t value of the path tracker
 *   which runs in quad double precision. */

int padcon_get_standard_step_size ( double *step );
/*
 * DESCRIPTION :
 *   Returns in step the current step size of the path tracker
 *   which runs in standard double precision. */

int padcon_get_dobldobl_step_size ( double *step );
/*
 * DESCRIPTION :
 *   Returns in step the current step size of the path tracker
 *   which runs in double double precision. */

int padcon_get_quaddobl_step_size ( double *step );
/*
 * DESCRIPTION :
 *   Returns in step the current step size of the path tracker
 *   which runs in quad double precision. */

int padcon_get_standard_series_step ( double *step );
/*
 * DESCRIPTION :
 *   Returns in step the current series step size of the path tracker
 *   which runs in standard double precision. */

int padcon_get_dobldobl_series_step ( double *step );
/*
 * DESCRIPTION :
 *   Returns in step the current series step size of the path tracker
 *   which runs in double double precision. */

int padcon_get_quaddobl_series_step ( double *step );
/*
 * DESCRIPTION :
 *   Returns in step the current series step size of the path tracker
 *   which runs in quad double precision. */

int padcon_get_standard_pole_step ( double *step );
/*
 * DESCRIPTION :
 *   Returns in step the current pole step size of the path tracker
 *   which runs in standard double precision. */

int padcon_get_dobldobl_pole_step ( double *step );
/*
 * DESCRIPTION :
 *   Returns in step the current pole step size of the path tracker
 *   which runs in double double precision. */

int padcon_get_quaddobl_pole_step ( double *step );
/*
 * DESCRIPTION :
 *   Returns in step the current pole step size of the path tracker
 *   which runs in quad double precision. */

int padcon_get_standard_estimated_distance ( double *eta );
/*
 * DESCRIPTION :
 *   Returns in eta the estimated distance to the closest solution by
 *   the path tracker which runs in standard double precision. */

int padcon_get_dobldobl_estimated_distance ( double *eta );
/*
 * DESCRIPTION :
 *   Returns in eta the estimated distance to the closest solution by 
 *   the path tracker which runs in double double precision. */

int padcon_get_quaddobl_estimated_distance ( double *eta );
/*
 * DESCRIPTION :
 *   Returns in eta the estimated distance to the closest solution by
 *   the path tracker which runs in quad double precision. */

int padcon_get_standard_hessian_step ( double *step );
/*
 * DESCRIPTION :
 *   Returns in step the current Hessian step size of the path tracker
 *   which runs in standard double precision. */

int padcon_get_dobldobl_hessian_step ( double *step );
/*
 * DESCRIPTION :
 *   Returns in step the current Hessian step size of the path tracker
 *   which runs in double double precision. */

int padcon_get_quaddobl_hessian_step ( double *step );
/*
 * DESCRIPTION :
 *   Returns in step the current Hessian step size of the path tracker
 *   which runs in quad double precision. */

int padcon_get_standard_series_coefficient
 ( int leadidx, int idx, int verbose, double* cre, double* cim );
/*
 * DESCRIPTION :
 *   Returns in cre and cim the real and imaginary parts of the series
 *   coefficient of component with leadidx at position idx, of the
 *   series computed by the predictor in double precision. */

int padcon_get_dobldobl_series_coefficient
 ( int leadidx, int idx, int verbose, double* cre, double* cim );
/*
 * DESCRIPTION :
 *   Returns in cre and cim the real and imaginary parts of the series
 *   coefficient of component with leadidx at position idx, of the
 *   series computed by the predictor in double double precision.
 *   The doubles are the highest parts of the double doubles. */

int padcon_get_quaddobl_series_coefficient
 ( int leadidx, int idx, int verbose, double* cre, double* cim );
/*
 * DESCRIPTION :
 *   Returns in cre and cim the real and imaginary parts of the series
 *   coefficient of component with leadidx at position idx, of the
 *   series computed by the predictor in quad double precision.
 *   The doubles are the highest parts of the quad doubles. */

int padcon_get_standard_numerator_coefficient
 ( int leadidx, int idx, int verbose, double* cre, double* cim );
/*
 * DESCRIPTION :
 *   Returns in cre and cim the real and imaginary parts of the 
 *   coefficient of the numerator of the Pade approximant,
 *   at the component with leadidx at position idx,
 *   computed by the predictor in double precision. */

int padcon_get_dobldobl_numerator_coefficient
 ( int leadidx, int idx, int verbose, double* cre, double* cim );
/*
 * DESCRIPTION :
 *   Returns in cre and cim the real and imaginary parts of the
 *   coefficient of the numerator of the Pade approximant,
 *   at the component with leadidx at position idx,
 *   computed by the predictor in double double precision.
 *   The doubles are the highest parts of the double doubles. */

int padcon_get_quaddobl_numerator_coefficient
 ( int leadidx, int idx, int verbose, double* cre, double* cim );
/*
 * DESCRIPTION :
 *   Returns in cre and cim the real and imaginary parts of the series
 *   coefficient of the numerator of the Pade approximant,
 *   at the component with leadidx at position idx,
 *   computed by the predictor in quad double precision.
 *   The doubles are the highest parts of the quad doubles. */

int padcon_get_standard_denominator_coefficient
 ( int leadidx, int idx, int verbose, double* cre, double* cim );
/*
 * DESCRIPTION :
 *   Returns in cre and cim the real and imaginary parts of the 
 *   coefficient of the denominator of the Pade approximant,
 *   at the component with leadidx at position idx,
 *   computed by the predictor in double precision. */

int padcon_get_dobldobl_denominator_coefficient
 ( int leadidx, int idx, int verbose, double* cre, double* cim );
/*
 * DESCRIPTION :
 *   Returns in cre and cim the real and imaginary parts of the
 *   coefficient of the denominator of the Pade approximant,
 *   at the component with leadidx at position idx,
 *   computed by the predictor in double double precision.
 *   The doubles are the highest parts of the double doubles. */

int padcon_get_quaddobl_denominator_coefficient
 ( int leadidx, int idx, int verbose, double* cre, double* cim );
/*
 * DESCRIPTION :
 *   Returns in cre and cim the real and imaginary parts of the series
 *   coefficient of the denominator of the Pade approximant,
 *   at the component with leadidx at position idx,
 *   computed by the predictor in quad double precision.
 *   The doubles are the highest parts of the quad doubles. */

int padcon_get_standard_pole
 ( int leadidx, int poleidx, int verbose, double* cre, double* cim );
/*
 * DESCRIPTION :
 *   Returns in cre and cim the real and imaginary parts of the pole
 *   Pade approximant with leadidx at position poleidx,
 *   computed by the predictor in double precision. */

int padcon_get_dobldobl_pole
 ( int leadidx, int poleidx, int verbose, double* cre, double* cim );
/*
 * DESCRIPTION :
 *   Returns in cre and cim the real and imaginary parts of the pole
 *   Pade approximant with leadidx at position poleidx,
 *   computed by the predictor in double double precision.
 *   The returned doubles are the highest parts of the double doubles. */

int padcon_get_quaddobl_pole
 ( int leadidx, int poleidx, int verbose, double* cre, double* cim );
/*
 * DESCRIPTION :
 *   Returns in cre and cim the real and imaginary parts of the pole
 *   Pade approximant with leadidx at position poleidx,
 *   computed by the predictor in quad double precision.
 *   The returned doubles are the highest parts of the quad doubles. */

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
