/* This file contains the prototypes for the development of the C++ and C
 * interface to the Ada code in PHCpack. */

#ifndef _LIB2PATH_H_
#define _LIB2PATH_H_

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <string.h>

#include <qd/qd_real.h>

#include "syscon.h"
#include "solcon.h"
#include "phcpack.h"
#include "jump_track.h"

#include "complexH.h"
#include "poly.h"
#include "polysol.h"

#ifdef compilewgpp
extern "C" void adainit( void );
extern "C" void adafinal( void );
#else
extern void adainit( void );
extern void adafinal( void );
#endif

void var_name ( char* x_string, int x_string_len, string*& x_names, int& dim );
/*
 * DESCRIPTION :
 *   Auxiliary function to read the names of the variables in a system. */

void lib2path_read_standard_sys
 ( int verbose, PolySys<complexH<double>,double>& sys );
void lib2path_read_dobldobl_sys
 ( int verbose, PolySys<complexH<dd_real>,dd_real>& sys );
void lib2path_read_quaddobl_sys
 ( int verbose, PolySys<complexH<qd_real>,qd_real>& sys );
/*
 * DESCRIPTION :
 *   Reads a polynomial system from the systems container
 *   into the data structure provided by the parameter sys.
 *   If verbose > 0, then extra output is written to screen. */

void lib2path_read_standard_sols
 ( PolySys<complexH<double>,double>& start_sys,
   PolySolSet<complexH<double>,double>& sols );
void lib2path_read_dobldobl_sols
 ( PolySys<complexH<dd_real>,dd_real>& start_sys,
   PolySolSet<complexH<dd_real>,dd_real>& sols );
void lib2path_read_quaddobl_sols
 ( PolySys<complexH<qd_real>,qd_real>& start_sys,
   PolySolSet<complexH<qd_real>,qd_real>& sols );
/* 
 * DESCRIPTION :
 *   Reads the corresponding solutions of the start system start_sys,
 *   into the provided data structure, respectively sols. 
 *   The dimension of the start_sys is used on input. */

void lib2path_write_standard_sols
 ( PolySolSet<complexH<double>,double>& sols );
void lib2path_write_dobldobl_sols
 ( PolySolSet<complexH<dd_real>,dd_real>& sols );
void lib2path_write_quaddobl_sols
 ( PolySolSet<complexH<qd_real>,qd_real>& sols );
/*
 * DESCRIPTION :
 *   Takes the solutions in sols and places them in the solutions container. */

void lib2path_read_standard_homotopy
 ( char* start_file, char* target_file,
   PolySys<complexH<double>,double>& start_sys,
   PolySys<complexH<double>,double>& target_sys,
   PolySolSet<complexH<double>,double>& sols );
void lib2path_read_dobldobl_homotopy
 ( char* start_file, char* target_file,
   PolySys<complexH<dd_real>,dd_real>& start_sys,
   PolySys<complexH<dd_real>,dd_real>& target_sys,
   PolySolSet<complexH<dd_real>,dd_real>& sols );
void lib2path_read_quaddobl_homotopy
 ( char* start_file, char* target_file,
   PolySys<complexH<qd_real>,qd_real>& start_sys,
   PolySys<complexH<qd_real>,qd_real>& target_sys,
   PolySolSet<complexH<qd_real>,qd_real>& sols );
/*
 * DESCRIPTION :
 *   Given in the strings the names of the files for start and target system,
 *   the C interface to PHCpack is used to parse the systems,
 *   and the solutions of the start system. */

int standard_newton
 ( int verbose, PolySys<complexH<double>,double>& p,
   PolySolSet<complexH<double>,double>& s );
int dobldobl_newton
 ( int verbose, PolySys<complexH<dd_real>,dd_real>& p,
   PolySolSet<complexH<dd_real>,dd_real>& s );
int quaddobl_newton
 ( int verbose, PolySys<complexH<qd_real>,qd_real>& p,
   PolySolSet<complexH<qd_real>,qd_real>& s );
/*
 * DESCRIPTION :
 *   Applies Newton's method to the first solution in s,
 *   on the polynomial system p.
 *   If verbose > 0, then additional output is written to screen.
 *   Default values for the path parameters will be applied. */

int standard_newton_with_pars
 ( int verbose, Parameter pars,
   PolySys<complexH<double>,double>& p,
   PolySolSet<complexH<double>,double>& s );
int dobldobl_newton_with_pars
 ( int verbose, Parameter pars,
   PolySys<complexH<dd_real>,dd_real>& p,
   PolySolSet<complexH<dd_real>,dd_real>& s );
int quaddobl_newton_with_pars
 ( int verbose, Parameter pars,
   PolySys<complexH<qd_real>,qd_real>& p,
   PolySolSet<complexH<qd_real>,qd_real>& s );
/*
 * DESCRIPTION :
 *   Applies Newton's method to the first solution in s,
 *   on the polynomial system p.
 *   If verbose > 0, then additional output is written to screen.
 *   Tuned values for the path parameter can be given on input. */

int standard_onetrack
 ( int verbose, double regamma, double imgamma,
   PolySys<complexH<double>,double>& p,
   PolySys<complexH<double>,double>& q,
   PolySolSet<complexH<double>,double>& s );
int dobldobl_onetrack
 ( int verbose, double regamma, double imgamma,
   PolySys<complexH<dd_real>,dd_real>& p,
   PolySys<complexH<dd_real>,dd_real>& q,
   PolySolSet<complexH<dd_real>,dd_real>& s );
int quaddobl_onetrack
 ( int verbose, double regamma, double imgamma,
   PolySys<complexH<qd_real>,qd_real>& p,
   PolySys<complexH<qd_real>,qd_real>& q,
   PolySolSet<complexH<qd_real>,qd_real>& s );
/*
 * DESCRIPTION :
 *   Tracks one path defined by an artificial parameter homotopy,
 *   starting at a solution s of q and ending at a solution of p.
 *   Default values for the parameter are applied.
 *
 * ON ENTRY :
 *   verbose   if > 0, then additional output is written to screen;
 *   regamma   real part of the gamma constant;
 *   imgamma   imaginary part of the gamma constant;
 *   p         target system in the homotopy;
 *   q         start system in the homotopy;
 *   s         a solution of the start system q. */

int standard_onetrack_with_pars
 ( int verbose, double regamma, double imgamma, Parameter pars,
   PolySys<complexH<double>,double>& p,
   PolySys<complexH<double>,double>& q,
   PolySolSet<complexH<double>,double>& s );
int dobldobl_onetrack_with_pars
 ( int verbose, double regamma, double imgamma, Parameter pars,
   PolySys<complexH<dd_real>,dd_real>& p,
   PolySys<complexH<dd_real>,dd_real>& q,
   PolySolSet<complexH<dd_real>,dd_real>& s );
int quaddobl_onetrack_with_pars
 ( int verbose, double regamma, double imgamma, Parameter pars,
   PolySys<complexH<qd_real>,qd_real>& p,
   PolySys<complexH<qd_real>,qd_real>& q,
   PolySolSet<complexH<qd_real>,qd_real>& s );
/*
 * DESCRIPTION :
 *   Tracks one path defined by an artificial parameter homotopy,
 *   starting at a solution s of q and ending at a solution of p.
 *   Tuned values for the parameter can be given on input.
 *
 * ON ENTRY :
 *   verbose   if > 0, then additional output is written to screen;
 *   regamma   real part of the gamma constant;
 *   imgamma   imaginary part of the gamma constant;
 *   pars      values for the path parameters;
 *   p         target system in the homotopy;
 *   q         start system in the homotopy;
 *   s         a solution of the start system q. */

int standard_manytrack
 ( int verbose, double regamma, double imgamma,
   PolySys<complexH<double>,double>& p,
   PolySys<complexH<double>,double>& q,
   PolySolSet<complexH<double>,double>& s );
int dobldobl_manytrack
 ( int verbose, double regamma, double imgamma,
   PolySys<complexH<dd_real>,dd_real>& p,
   PolySys<complexH<dd_real>,dd_real>& q,
   PolySolSet<complexH<dd_real>,dd_real>& s );
int quaddobl_manytrack
 ( int verbose, double regamma, double imgamma,
   PolySys<complexH<qd_real>,qd_real>& p,
   PolySys<complexH<qd_real>,qd_real>& q,
   PolySolSet<complexH<qd_real>,qd_real>& s );
/*
 * DESCRIPTION :
 *   Tracks many paths defined by an artificial parameter homotopy,
 *   starting at solutions in s of q and ending at solutions of p.
 *   Default values for the path parameters are applied.
 *
 * ON ENTRY :
 *   verbose   if > 0, then additional output is written to screen;
 *   regamma   real part of the gamma constant;
 *   imgamma   imaginary part of the gamma constant;
 *   p         target system in the homotopy;
 *   q         start system in the homotopy;
 *   s         solutions of the start system q. */

int standard_manytrack_with_pars
 ( int verbose, double regamma, double imgamma, Parameter pars,
   PolySys<complexH<double>,double>& p,
   PolySys<complexH<double>,double>& q,
   PolySolSet<complexH<double>,double>& s );
int dobldobl_manytrack_with_pars
 ( int verbose, double regamma, double imgamma, Parameter pars,
   PolySys<complexH<dd_real>,dd_real>& p,
   PolySys<complexH<dd_real>,dd_real>& q,
   PolySolSet<complexH<dd_real>,dd_real>& s );
int quaddobl_manytrack_with_pars
 ( int verbose, double regamma, double imgamma, Parameter pars,
   PolySys<complexH<qd_real>,qd_real>& p,
   PolySys<complexH<qd_real>,qd_real>& q,
   PolySolSet<complexH<qd_real>,qd_real>& s );
/*
 * DESCRIPTION :
 *   Tracks many paths defined by an artificial parameter homotopy,
 *   starting at solutions in s of q and ending at solutions of p.
 *   Tuned values for the parameters can be provided on input.
 *
 * ON ENTRY :
 *   verbose   if > 0, then additional output is written to screen;
 *   regamma   real part of the gamma constant;
 *   imgamma   imaginary part of the gamma constant;
 *   pars      values for the path parameters;
 *   p         target system in the homotopy;
 *   q         start system in the homotopy;
 *   s         solutions of the start system q. */

int standard_ade_newton ( int verbose );
int dobldobl_ade_newton ( int verbose );
int quaddobl_ade_newton ( int verbose );
/*
 * DESCRIPTION :
 *   Calls Newton's method with algorithmic differentiation.
 *   If verbose > 0, then additional output will be written.
 *   Default values for the path parameters are applied. */

int standard_ade_newton_with_pars ( int verbose, Parameter pars );
int dobldobl_ade_newton_with_pars ( int verbose, Parameter pars );
int quaddobl_ade_newton_with_pars ( int verbose, Parameter pars );
/*
 * DESCRIPTION :
 *   Calls Newton's method with algorithmic differentiation.
 *   If verbose > 0, then additional output will be written.
 *   Tuned values of the parameters can be given on input in pars. */

int standard_ade_onepath
 ( int verbose, double regamma, double imgamma );
int dobldobl_ade_onepath
 ( int verbose, double regamma, double imgamma );
int quaddobl_ade_onepath
 ( int verbose, double regamma, double imgamma );
/*
 * DESCRIPTION :
 *   A C++ function to track one solution path.
 *   Default values for the path parameters are applied.
 *
 * ON ENTRY :
 *   verbose   if > 0, then additional output is written to screen;
 *   regamma   real part of the gamma constant;
 *   imgamma   imaginary part of the gamma constant. */

int standard_ade_onepath_with_pars
 ( int verbose, double regamma, double imgamma, Parameter pars );
int dobldobl_ade_onepath_with_pars
 ( int verbose, double regamma, double imgamma, Parameter pars );
int quaddobl_ade_onepath_with_pars
 ( int verbose, double regamma, double imgamma, Parameter pars );
/*
 * DESCRIPTION :
 *   A C++ function to track one solution path.
 *   Tuned values for the path parameters can be given in pars.
 *
 * ON ENTRY :
 *   verbose   if > 0, then additional output is written to screen;
 *   regamma   real part of the gamma constant;
 *   imgamma   imaginary part of the gamma constant;
 *   pars      values for the path parameters. */

int standard_ade_manypaths
 ( int verbose, double regamma, double imgamma );
int dobldobl_ade_manypaths
 ( int verbose, double regamma, double imgamma );
int quaddobl_ade_manypaths
 ( int verbose, double regamma, double imgamma );
/*
 * DESCRIPTION :
 *   A C++ function to track one solution path,
 *   encapsulated as a C function for to be called from Ada.
 *   Default values for the path parameters are applied.
 *
 * ON ENTRY :
 *   verbose   if > 0, then additional output is written to screen;
 *   regamma   real part of the gamma constant;
 *   imgamma   imaginary part of the gamma constant. */

int standard_ade_manypaths_with_pars
 ( int verbose, double regamma, double imgamma, Parameter pars );
int dobldobl_ade_manypaths_with_pars
 ( int verbose, double regamma, double imgamma, Parameter pars );
int quaddobl_ade_manypaths_with_pars
 ( int verbose, double regamma, double imgamma, Parameter pars );
/*
 * DESCRIPTION :
 *   A C++ function to track one solution path,
 *   encapsulated as a C function for to be called from Ada.
 *   Tuned values for the path parameters can be given in pars.
 *
 * ON ENTRY :
 *   verbose   if > 0, then additional output is written to screen;
 *   regamma   real part of the gamma constant;
 *   imgamma   imaginary part of the gamma constant;
 *   pars      values for the path parameters. */

extern "C" int standard_adenewton ( int verbose );
extern "C" int dobldobl_adenewton ( int verbose );
extern "C" int quaddobl_adenewton ( int verbose );
/*
 * DESCRIPTION :
 *   Calls Newton's method with algorithmic differentiation,
 *   encapsulated as a C function for to be called from Ada.
 *   Notice the lack of an underscore after the _ade in the function names.
 *   If verbose > 0, then additional output will be written.
 *   Default values for the path parameters are applied. */

extern "C" int standard_adenewton_with_pars
 ( int verbose, int max_step, int n_predictor,
   double step_increase, double step_decrease,
   double max_delta_t, double max_delta_t_end, double min_delta_t,
   double err_max_res, double err_max_delta_x, double err_max_first_delta_x,
   int max_it, double err_min_round_off,
   int max_it_refine, double err_min_round_off_refine );
extern "C" int dobldobl_adenewton_with_pars
 ( int verbose, int max_step, int n_predictor,
   double step_increase, double step_decrease,
   double max_delta_t, double max_delta_t_end, double min_delta_t,
   double err_max_res, double err_max_delta_x, double err_max_first_delta_x,
   int max_it, double err_min_round_off,
   int max_it_refine, double err_min_round_off_refine );
extern "C" int quaddobl_adenewton_with_pars
 ( int verbose, int max_step, int n_predictor,
   double step_increase, double step_decrease,
   double max_delta_t, double max_delta_t_end, double min_delta_t,
   double err_max_res, double err_max_delta_x, double err_max_first_delta_x,
   int max_it, double err_min_round_off,
   int max_it_refine, double err_min_round_off_refine );
/*
 * DESCRIPTION :
 *   Calls Newton's method with algorithmic differentiation,
 *   encapsulated as a C function for to be called from Ada.
 *   Notice the lack of an underscore after the _ade in the function names.
 *   If verbose > 0, then additional output will be written.
 *   All 14 values for the path parameters must be provided.
 *
 * ON ENTRY :
 *   verbose              if > 0, then additional output is written to screen;
 *   max_step             maximum number of steps along a path;
 *   n_predictor          number of points used in the predictor;
 *   step_increase        increase factor of the predictor;
 *   step_decrease        decrease factor of the precdictor;
 *   max_delta_t          maximum step size along a path;
 *   max_delta_t_end      maximum step size at the end of a path;
 *   min_delta_t          minimum step size;
 *   err_max_res          tolerance on the residual;
 *   err_max_delta_x      tolerance on the corrector update;
 *   err_max_first_delta_x is the tolerance on the first correction update;
 *   max_it               maximum number of iterations of the corrector;
 *   err_min_round_off    tolerance on the corrector;
 *   max_it_refine        number of steps in the Newton root refiner;
 *   err_min_round_off_refine is the tolerance for the final refinement. */

extern "C" int standard_adeonepath
 ( int verbose, double regamma, double imgamma );
extern "C" int dobldobl_adeonepath
 ( int verbose, double regamma, double imgamma );
extern "C" int quaddobl_adeonepath
 ( int verbose, double regamma, double imgamma );
/*
 * DESCRIPTION :
 *   A C++ function to track one solution path,
 *   encapsulated as a C function for to be called from Ada.
 *   Notice the lack of an underscore after the _ade in the function names.
 *   Default values for the path parameters are applied.
 *
 * ON ENTRY :
 *   verbose   if > 0, then additional output is written to screen;
 *   regamma   real part of the gamma constant;
 *   imgamma   imaginary part of the gamma constant. */

extern "C" int standard_adeonepath_with_pars
 ( int verbose, double regamma, double imgamma,
   int max_step, int n_predictor,
   double step_increase, double step_decrease,
   double max_delta_t, double max_delta_t_end, double min_delta_t,
   double err_max_res, double err_max_delta_x, double err_max_first_delta_x,
   int max_it, double err_min_round_off,
   int max_it_refine, double err_min_round_off_refine );
extern "C" int dobldobl_adeonepath_with_pars
 ( int verbose, double regamma, double imgamma,
   int max_step, int n_predictor,
   double step_increase, double step_decrease,
   double max_delta_t, double max_delta_t_end, double min_delta_t,
   double err_max_res, double err_max_delta_x, double err_max_first_delta_x,
   int max_it, double err_min_round_off,
   int max_it_refine, double err_min_round_off_refine );
extern "C" int quaddobl_adeonepath_with_pars
 ( int verbose, double regamma, double imgamma,
   int max_step, int n_predictor,
   double step_increase, double step_decrease,
   double max_delta_t, double max_delta_t_end, double min_delta_t,
   double err_max_res, double err_max_delta_x, double err_max_first_delta_x,
   int max_it, double err_min_round_off,
   int max_it_refine, double err_min_round_off_refine );
/*
 * DESCRIPTION :
 *   A C++ function to track one solution path,
 *   encapsulated as a C function for to be called from Ada.
 *   Notice the lack of an underscore after the _ade in the function names.
 *   Values for all 14 path parameters have to be provided.
 *
 * ON ENTRY :
 *   verbose              if > 0, then additional output is written to screen;
 *   regamma              real part of the gamma constant;
 *   imgamma              imaginary part of the gamma constant;
 *   max_step             maximum number of steps along a path;
 *   n_predictor          number of points used in the predictor;
 *   step_increase        increase factor of the predictor;
 *   step_decrease        decrease factor of the precdictor;
 *   max_delta_t          maximum step size along a path;
 *   max_delta_t_end      maximum step size at the end of a path;
 *   min_delta_t          minimum step size;
 *   err_max_res          tolerance on the residual;
 *   err_max_delta_x      tolerance on the corrector update;
 *   err_max_first_delta_x is the tolerance on the first correction update;
 *   max_it               maximum number of iterations of the corrector;
 *   err_min_round_off    tolerance on the corrector;
 *   max_it_refine        number of steps in the Newton root refiner;
 *   err_min_round_off_refine is the tolerance for the final refinement. */

extern "C" int standard_ademanypaths
 ( int verbose, double regamma, double imgamma );
extern "C" int dobldobl_ademanypaths
 ( int verbose, double regamma, double imgamma );
extern "C" int quaddobl_ademanypaths
 ( int verbose, double regamma, double imgamma );
/*
 * DESCRIPTION :
 *   A C++ function to track many solution paths,
 *   encapsulated as a C function for to be called from Ada.
 *   Notice the lack of an underscore after the _ade in the function names.
 *   Default values for the path parameters are applied.
 *
 * ON ENTRY :
 *   verbose              if > 0, then additional output is written to screen;
 *   regamma              real part of the gamma constant;
 *   imgamma              imaginary part of the gamma constant. */

extern "C" int standard_ademanypaths_with_pars
 ( int verbose, double regamma, double imgamma,
   int max_step, int n_predictor,
   double step_increase, double step_decrease,
   double max_delta_t, double max_delta_t_end, double min_delta_t,
   double err_max_res, double err_max_delta_x, double err_max_first_delta_x,
   int max_it, double err_min_round_off,
   int max_it_refine, double err_min_round_off_refine );
extern "C" int dobldobl_ademanypaths_with_pars
 ( int verbose, double regamma, double imgamma,
   int max_step, int n_predictor,
   double step_increase, double step_decrease,
   double max_delta_t, double max_delta_t_end, double min_delta_t,
   double err_max_res, double err_max_delta_x, double err_max_first_delta_x,
   int max_it, double err_min_round_off,
   int max_it_refine, double err_min_round_off_refine );
extern "C" int quaddobl_ademanypaths_with_pars
 ( int verbose, double regamma, double imgamma,
   int max_step, int n_predictor,
   double step_increase, double step_decrease,
   double max_delta_t, double max_delta_t_end, double min_delta_t,
   double err_max_res, double err_max_delta_x, double err_max_first_delta_x,
   int max_it, double err_min_round_off,
   int max_it_refine, double err_min_round_off_refine );
/*
 * DESCRIPTION :
 *   A C++ function to track many solution paths,
 *   encapsulated as a C function for to be called from Ada.
 *   Notice the lack of an underscore after the _ade in the function names.
 *   Values for all 14 path parameters have to be provided.
 *
 * ON ENTRY :
 *   verbose              if > 0, then additional output is written to screen;
 *   regamma              real part of the gamma constant;
 *   imgamma              imaginary part of the gamma constant;
 *   max_step             maximum number of steps along a path;
 *   n_predictor          number of points used in the predictor;
 *   step_increase        increase factor of the predictor;
 *   step_decrease        decrease factor of the precdictor;
 *   max_delta_t          maximum step size along a path;
 *   max_delta_t_end      maximum step size at the end of a path;
 *   min_delta_t          minimum step size;
 *   err_max_res          tolerance on the residual;
 *   err_max_delta_x      tolerance on the corrector update;
 *   err_max_first_delta_x is the tolerance on the first correction update;
 *   max_it               maximum number of iterations of the corrector;
 *   err_min_round_off    tolerance on the corrector;
 *   max_it_refine        number of steps in the Newton root refiner;
 *   err_min_round_off_refine is the tolerance for the final refinement. */

int set_path_parameter_value ( Parameter pars, int idx, double val );
/*
 * DESCRIPTION :
 *   Sets the value of the path parameter idx to the value val.
 *   Returns 0 if the value of idx is in the range 1 to 14.
 *   Returns -1 if the value of idx is out of the range 1 to 14.
 *
 * ON ENTRY :
 *   The value val sets the value corresponding to the index idx,
 *   where idx is the index to one of the 14 parameters:
 *    1 : maximum number of steps along a path;
 *    2 : number of points used in the predictor;
 *    3 : increase factor of the predictor;
 *    4 : decrease factor of the precdictor;
 *    5 : maximum step size along a path;
 *    6 : maximum step size at the end of a path;
 *    7 : minimum step size;
 *    8 : tolerance on the residual;
 *    9 : tolerance on the corrector update;
 *   10 : the tolerance on the first correction update;
 *   11 : maximum number of iterations of the corrector;
 *   12 : tolerance on the corrector;
 *   13 : number of steps in the Newton root refiner;
 *   14 : tolerance for the final refinement. */

int get_path_parameter_value ( Parameter pars, int idx, double* val );
/*
 * DESCRIPTION :
 *   Returns in val the value of the path parameter with index idx.
 *   Returns 0 if the value of idx is in the range 1 to 14.
 *   Returns -1 if the value of idx is out of the range 1 to 14.
 *
 * ON ENTRY :
 *   The index corresponding to the value should be one of the following
 *    1 : maximum number of steps along a path;
 *    2 : number of points used in the predictor;
 *    3 : increase factor of the predictor;
 *    4 : decrease factor of the precdictor;
 *    5 : maximum step size along a path;
 *    6 : maximum step size at the end of a path;
 *    7 : minimum step size;
 *    8 : tolerance on the residual;
 *    9 : tolerance on the corrector update;
 *   10 : the tolerance on the first correction update;
 *   11 : maximum number of iterations of the corrector;
 *   12 : tolerance on the corrector;
 *   13 : number of steps in the Newton root refiner;
 *   14 : tolerance for the final refinement.
 *
 * ON RETURN :
 *   val    the value with the corresponding index. */

int get_default_path_parameters
 ( int precision, int* max_step, int* n_predictor,
   double* step_increase, double* step_decrease,
   double* max_delta_t, double* max_delta_t_end, double* min_delta_t,
   double* err_max_res, double* err_max_delta_x, double* err_max_first_delta_x,
   int* max_it, double* err_min_round_off,
   int* max_it_refine, double* err_min_round_off_refine );
/*
 * DESCRIPTION :
 *   Returns all default values for the path parameters,
 *   corresponding to the value of the precision.
 *
 * ON ENTRY :
 *   precision is the number of decimal places in the working precision,
 *             should be 16, 32, or 64 for double, double double,
 *             or quad double precision respectively.
 *
 * ON RETURN :
 *   max_step             maximum number of steps along a path;
 *   n_predictor          number of points used in the predictor;
 *   step_increase        increase factor of the predictor;
 *   step_decrease        decrease factor of the precdictor;
 *   max_delta_t          maximum step size along a path;
 *   max_delta_t_end      maximum step size at the end of a path;
 *   min_delta_t          minimum step size;
 *   err_max_res          tolerance on the residual;
 *   err_max_delta_x      tolerance on the corrector update;
 *   err_max_first_delta_x is the tolerance on the first correction update;
 *   max_it               maximum number of iterations of the corrector;
 *   err_min_round_off    tolerance on the corrector;
 *   max_it_refine        number of steps in the Newton root refiner;
 *   err_min_round_off_refine is the tolerance for the final refinement. */

int standard_ademanypaths_with_parameters
 ( int verbose, double regamma, double imgamma,
   int max_step, int n_predictor,
   double step_increase, double step_decrease,
   double max_delta_t, double max_delta_t_end, double min_delta_t,
   double err_max_res, double err_max_delta_x, double err_max_first_delta_x,
   int max_it, double err_min_round_off,
   int max_it_refine, double err_min_round_off_refine );
int dobldobl_ademanypaths_with_parameters
 ( int verbose, double regamma, double imgamma,
   int max_step, int n_predictor,
   double step_increase, double step_decrease,
   double max_delta_t, double max_delta_t_end, double min_delta_t,
   double err_max_res, double err_max_delta_x, double err_max_first_delta_x,
   int max_it, double err_min_round_off,
   int max_it_refine, double err_min_round_off_refine );
int quaddobl_ademanypaths_with_parameters
 ( int verbose, double regamma, double imgamma,
   int max_step, int n_predictor,
   double step_increase, double step_decrease,
   double max_delta_t, double max_delta_t_end, double min_delta_t,
   double err_max_res, double err_max_delta_x, double err_max_first_delta_x,
   int max_it, double err_min_round_off,
   int max_it_refine, double err_min_round_off_refine );
/*
 * DESCRIPTION :
 *   Track many solution paths, to be used in the phcpy2c module.
 *   Notice the lack of an underscore after the _ade in the function names.
 *   Values for all 14 path parameters have to be provided.
 *
 * ON ENTRY :
 *   verbose              if > 0, then additional output is written to screen;
 *   regamma              real part of the gamma constant;
 *   imgamma              imaginary part of the gamma constant;
 *   max_step             maximum number of steps along a path;
 *   n_predictor          number of points used in the predictor;
 *   step_increase        increase factor of the predictor;
 *   step_decrease        decrease factor of the precdictor;
 *   max_delta_t          maximum step size along a path;
 *   max_delta_t_end      maximum step size at the end of a path;
 *   min_delta_t          minimum step size;
 *   err_max_res          tolerance on the residual;
 *   err_max_delta_x      tolerance on the corrector update;
 *   err_max_first_delta_x is the tolerance on the first correction update;
 *   max_it               maximum number of iterations of the corrector;
 *   err_min_round_off    tolerance on the corrector;
 *   max_it_refine        number of steps in the Newton root refiner;
 *   err_min_round_off_refine is the tolerance for the final refinement. */

#endif
