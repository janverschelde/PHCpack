/* The file jump_track.h contains prototypes for the operations to perform
 * path tracking with jumpstarting homotopies.
 * By default, compilation with gcc is assumed.
 * To compile with a C++ compiler such as g++, the flag compilewgpp must
 * be defined as "g++ -Dcompilewgpp=1." */

#ifndef __JUMP_TRACK_H__
#define __JUMP_TRACK_H__

#ifdef compilewgpp
extern "C" void adainit( void );
extern "C" int _ada_use_c2phc4c ( int task, int *a, int *b, double *c, int v );
extern "C" void adafinal( void );
#else
extern void adainit( void );
extern int _ada_use_c2phc4c ( int task, int *a, int *b, double *c, int v );
extern void adafinal( void );
#endif

int read_target_system_without_solutions ( void );
/*
 * DESCRIPTION :
 *   Prompts the user for a file name and reads the target system from file.
 *   Even if available on file, its solutions will not be read and stored. */

int read_dobldobl_target_system_without_solutions ( void );
/*
 * DESCRIPTION :
 *   Prompts the user for a file name and reads the target system from file.
 *   The coefficients are parsed in double double precision and the system
 *   is stored as a double double target system.
 *   Even if available on file, its solutions will not be read and stored. */

int read_quaddobl_target_system_without_solutions ( void );
/*
 * DESCRIPTION :
 *   Prompts the user for a file name and reads the target system from file.
 *   The coefficients are parsed in double double precision and the system
 *   is stored as a double double target system.
 *   Even if available on file, its solutions will not be read and stored. */

int read_named_target_without_solutions ( int n, char *s );
/*
 * DESCRITPION :
 *   Reads the target system system from file using the given name in s,
 *   a string of n characters.  Solutions will not be read. */

int read_start_system_without_solutions ( void );
/*
 * DESCRIPTION :
 *   Prompts the user for a file name and reads the start system from file.
 *   Even if start solutions are in that file, they will not be read. */

int read_named_start_without_solutions ( int n, char *s );
/*
 * DESCRITPION :
 *   Reads the start system system from file using the given name in s,
 *   a string of n characters.  Solutions will not be read. */

int read_named_linear_product_start_system ( int n, char *s );
/*
 * DESCRIPTION :
 *   Reads a linear-product start system from file using the name in s.
 *   The appropriate data structures in PHCpack are initialized. */

int silent_path_tracker
 ( int n, int *m, double *c,
   int *nbstep, int *nbfail, int *nbiter, int *nbsyst );
/*
 * DESCRIPTION :
 *   Tracks one path using a silent path tracker,
 *   in standard double precision.
 *
 * ON ENTRY :
 *   n       number of variables in the solution;
 *   m       multiplicity flag;
 *   c       coordinates of the solution.
 *
 * ON RETURN :
 *   m       value for the multiplicity;
 *   c       real and imaginary part of the continuation parameter,
 *           the real and imaginary parts of the solution coordinates,
 *           diagnostics: (err,rco,res) as the last 3 doubles;
 *   nbstep  #steps taken along the path;
 *   nbfail  #failures occurred along the path;
 *   nbiter  #corrector iterations done;
 *   nbsyst  #linear systems solved. */

int silent_dobldobl_path_tracker
 ( int n, int *m, double *c,
   int *nbstep, int *nbfail, int *nbiter, int *nbsyst );
/*
 * DESCRIPTION :
 *   Tracks one path using a silent path tracker,
 *   in double double precision.
 *
 * ON ENTRY :
 *   n       number of variables in the solution;
 *   m       multiplicity flag;
 *   c       coordinates of the solution.
 *
 * ON RETURN :
 *   m       value for the multiplicity;
 *   c       real and imaginary part of the continuation parameter,
 *           the real and imaginary parts of the solution coordinates,
 *           diagnostics: (err,rco,res) as the last 3 doubles;
 *   nbstep  #steps taken along the path;
 *   nbfail  #failures occurred along the path;
 *   nbiter  #corrector iterations done;
 *   nbsyst  #linear systems solved. */

int silent_quaddobl_path_tracker
 ( int n, int *m, double *c,
   int *nbstep, int *nbfail, int *nbiter, int *nbsyst );
/*
 * DESCRIPTION :
 *   Tracks one path using a silent path tracker,
 *   in quad double precision.
 *
 * ON ENTRY :
 *   n       number of variables in the solution;
 *   m       multiplicity flag;
 *   c       coordinates of the solution.
 *
 * ON RETURN :
 *   m       value for the multiplicity;
 *   c       real and imaginary part of the continuation parameter,
 *           the real and imaginary parts of the solution coordinates,
 *           diagnostics: (err,rco,res) as the last 3 doubles;
 *   nbstep  #steps taken along the path;
 *   nbfail  #failures occurred along the path;
 *   nbiter  #corrector iterations done;
 *   nbsyst  #linear systems solved. */

int reporting_path_tracker
 ( int n, int *m, double *c,
   int *nbstep, int *nbfail, int *nbiter, int *nbsyst );
/*
 * DESCRIPTION :
 *   Tracks one path, writing extra information while path tracking.
 *
 * ON ENTRY :
 *   n       number of variables in the solution;
 *   m       multiplicity flag;
 *   c       coordinates of the solution.
 *
 * ON RETURN :
 *   m       value for the multiplicity;
 *   c       real and imaginary part of the continuation parameter,
 *           the real and imaginary parts of the solution coordinates,
 *           diagnostics: (err,rco,res) as the last 3 doubles.
 *   nbstep  #steps taken along the path;
 *   nbfail  #failures occurred along the path;
 *   nbiter  #corrector iterations done;
 *   nbsyst  #linear systems solved. */

int reporting_dobldobl_path_tracker
 ( int n, int *m, double *c,
   int *nbstep, int *nbfail, int *nbiter, int *nbsyst );
/*
 * DESCRIPTION :
 *   Tracks one path, writing extra information while path tracking,
 *   in double double precision.
 *
 * ON ENTRY :
 *   n       number of variables in the solution;
 *   m       multiplicity flag;
 *   c       coordinates of the solution.
 *
 * ON RETURN :
 *   m       value for the multiplicity;
 *   c       real and imaginary part of the continuation parameter,
 *           the real and imaginary parts of the solution coordinates,
 *           diagnostics: (err,rco,res) as the last 3 doubles.
 *   nbstep  #steps taken along the path;
 *   nbfail  #failures occurred along the path;
 *   nbiter  #corrector iterations done;
 *   nbsyst  #linear systems solved. */

int reporting_quaddobl_path_tracker
 ( int n, int *m, double *c,
   int *nbstep, int *nbfail, int *nbiter, int *nbsyst );
/*
 * DESCRIPTION :
 *   Tracks one path, writing extra information while path tracking,
 *   in quad double precision.
 *
 * ON ENTRY :
 *   n       number of variables in the solution;
 *   m       multiplicity flag;
 *   c       coordinates of the solution.
 *
 * ON RETURN :
 *   m       value for the multiplicity;
 *   c       real and imaginary part of the continuation parameter,
 *           the real and imaginary parts of the solution coordinates,
 *           diagnostics: (err,rco,res) as the last 3 doubles.
 *   nbstep  #steps taken along the path;
 *   nbfail  #failures occurred along the path;
 *   nbiter  #corrector iterations done;
 *   nbsyst  #linear systems solved. */

int write_next_solution_with_diagnostics
 ( int *k, int n, int m, double *sol,
   int nbstep, int nbfail, int nbiter, int nbsyst );
/*
 * DESCRIPTION :
 *   The next solution is written to file, along with statistics
 *   of the path tracking process, and a preliminary diagnostic,
 *   for standard double precision.
 *
 * ON ENTRY :
 *   k        current number of solutions already written to file;
 *   n        dimension of the solution vector;
 *   m        multiplicity flag of the solution vector;
 *   sol      values of the continuation parameter, coordinates,
 *            and diagnostics of the solution, of dimension 2*n+5;
 *   nbstep   #steps taken along the path;
 *   nbfail   #failures occurred along the path;
 *   nbiter   #corrector iterations done;
 *   nbsyst   #linear systems solved. 
 *
 * ON RETURN :
 *   k         updated counter on number of solutions written to file. */

int write_next_dobldobl_solution_with_diagnostics
 ( int *k, int n, int m, double *sol,
   int nbstep, int nbfail, int nbiter, int nbsyst );
/*
 * DESCRIPTION :
 *   The next solution is written to file, along with statistics
 *   of the path tracking process, and a preliminary diagnostic,
 *   for double double precision.
 *
 * ON ENTRY :
 *   k        current number of solutions already written to file;
 *   n        dimension of the solution vector;
 *   m        multiplicity flag of the solution vector;
 *   sol      values of the continuation parameter, coordinates,
 *            and diagnostics of the solution, of dimension 2*n+5;
 *   nbstep   #steps taken along the path;
 *   nbfail   #failures occurred along the path;
 *   nbiter   #corrector iterations done;
 *   nbsyst   #linear systems solved. 
 *
 * ON RETURN :
 *   k         updated counter on number of solutions written to file. */

int write_next_quaddobl_solution_with_diagnostics
 ( int *k, int n, int m, double *sol,
   int nbstep, int nbfail, int nbiter, int nbsyst );
/*
 * DESCRIPTION :
 *   The next solution is written to file, along with statistics
 *   of the path tracking process, and a preliminary diagnostic,
 *   for quad double precision.
 *
 * ON ENTRY :
 *   k        current number of solutions already written to file;
 *   n        dimension of the solution vector;
 *   m        multiplicity flag of the solution vector;
 *   sol      values of the continuation parameter, coordinates,
 *            and diagnostics of the solution, of dimension 2*n+5;
 *   nbstep   #steps taken along the path;
 *   nbfail   #failures occurred along the path;
 *   nbiter   #corrector iterations done;
 *   nbsyst   #linear systems solved. 
 *
 * ON RETURN :
 *   k         updated counter on number of solutions written to file. */

int standard_crude_tracker ( int verbose );
/*
 * DESCRIPTION :
 *   A crude tracker appends the end point of a path directly to
 *   the solutions container, without refinement or postprocessing.
 *
 * REQUIRED :
 *   The target system, start system, and start solutions in standard
 *   double precision have been initialized in the containers.
 *
 * ON ENTRY :
 *   verbose  if 1, then the solution vectors are written to screen,
 *            otherwise, the crude tracker remains mute. */

int dobldobl_crude_tracker ( int verbose );
/*
 * DESCRIPTION :
 *   A crude tracker appends the end point of a path directly to
 *   the solutions container, without refinement or postprocessing.
 *
 * REQUIRED :
 *   The target system, start system, and start solutions in double
 *   double precision have been initialized in the containers.
 *
 * ON ENTRY :
 *   verbose  if 1, then the solution vectors are written to screen,
 *            otherwise, the crude tracker remains mute. */

int quaddobl_crude_tracker ( int verbose );
/*
 * DESCRIPTION :
 *   A crude tracker appends the end point of a path directly to
 *   the solutions container, without refinement or postprocessing.
 *
 * REQUIRED :
 *   The target system, start system, and start solutions in quad
 *   double precision have been initialized in the containers.
 *
 * ON ENTRY :
 *   verbose  if 1, then the solution vectors are written to screen,
 *            otherwise, the crude tracker remains mute. */

#endif
