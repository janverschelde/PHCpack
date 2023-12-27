/* This file contains the prototypes for the py2c interface functions,
 * linked with the quad double version of the Path library. */

void initialize ( void );
/*
 * DESCRIPTION :
 *   Calls adainit(), initializing the interface to the Ada code,
 *   setting the initialized flag to one and the finalized flag to zero,
 *   if the initialized flag was zero.
 *   Nothing happens if the initialized flag equals one. */

void finalize ( void );
/*
 * DESCRIPTION :
 *   Calls adafinal(), finalizing the interface to the Ada code,
 *   setting the finalized flag to one and the initialized flag to zero,
 *   if the finalized flag was zero.
 *   Nothing happens if the finalized flag equals one. */

/* The wrapping of functions with prototypes in phcpack.h starts here. */

static PyObject *py2c_PHCpack_version_string
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the version string of PHCpack.
 *   The version string is 40 characters long. */

static PyObject *py2c_set_seed ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Takes the value of the integer given on input
 *   and sets the seed for the random number generators.
 *   This fixing of the seed enables reproducible runs. */

static PyObject *py2c_get_seed ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the current value of the seed.
 *   Using this value in py2c_set_seed will ensure that the
 *   results of previous runs can be reproduced. */

static PyObject *py2c_read_standard_target_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Prompts the user to enter a target system that will
 *   be parsed in standard double precision.
 *   The failure code is returned, which is zero if all went well. */

static PyObject *py2c_read_standard_target_system_from_file
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   The two input arguments are a number and a string:
 *   1) The number equals the number of characters in the string.
 *   2) The string given on input is the name of a file which contains
 *   a target system to be parsed in standard double precision.
 *   The failure code is returned, which is zero if all went well. */

static PyObject *py2c_read_standard_start_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Prompts the user to enter a start system that will
 *   be parsed in standard double precision.
 *   The failure code is returned, which is zero if all went well. */

static PyObject *py2c_read_standard_start_system_from_file
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   The two input arguments are a number and a string:
 *   1) The number equals the number of characters in the string.
 *   2) The string given on input is the name of a file which contains
 *   a start system to be parsed in standard double precision.
 *   The failure code is returned, which is zero if all went well. */

static PyObject *py2c_read_dobldobl_target_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Prompts the user to enter a target system that will
 *   be parsed in double double precision.
 *   The failure code is returned, which is zero if all went well. */

static PyObject *py2c_read_dobldobl_target_system_from_file
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   The two input arguments are a number and a string:
 *   1) The number equals the number of characters in the string.
 *   2) The string given on input is the name of a file which contains
 *   a target system to be parsed in double double precision.
 *   The failure code is returned, which is zero if all went well. */

static PyObject *py2c_read_dobldobl_start_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Prompts the user to enter a start system that will
 *   be parsed in double double precision.
 *   The failure code is returned, which is zero if all went well. */

static PyObject *py2c_read_dobldobl_start_system_from_file
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   The two input arguments are a number and a string:
 *   1) The number equals the number of characters in the string.
 *   2) The string given on input is the name of a file which contains
 *   a start system to be parsed in double double precision.
 *   The failure code is returned, which is zero if all went well. */

static PyObject *py2c_read_quaddobl_target_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Prompts the user to enter a target system that will
 *   be parsed in quad double precision.
 *   The failure code is returned, which is zero if all went well. */

static PyObject *py2c_read_quaddobl_target_system_from_file
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   The two input arguments are a number and a string:
 *   1) The number equals the number of characters in the string.
 *   2) The string given on input is the name of a file which contains
 *   a target system to be parsed in quad double precision.
 *   The failure code is returned, which is zero if all went well. */

static PyObject *py2c_read_quaddobl_start_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Prompts the user to enter a start system that will
 *   be parsed in quad double precision.
 *   The failure code is returned, which is zero if all went well. */

static PyObject *py2c_read_quaddobl_start_system_from_file
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   The two input arguments are a number and a string:
 *   1) The number equals the number of characters in the string.
 *   2) The string given on input is the name of a file which contains
 *   a start system to be parsed in quad double precision.
 *   The failure code is returned, which is zero if all went well. */

static PyObject *py2c_define_output_file ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Prompts the user to define the output file.
 *   On return is the failure code, which is zero if all went well. */

static PyObject *py2c_write_standard_target_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Writes the target system as stored in standard double precision
 *   to screen or to the defined output file. */

static PyObject *py2c_write_standard_start_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Writes the start system as stored in standard double precision
 *   to screen or to the defined output file. */

/* Copying systems from and to containers. */

static PyObject *py2c_copy_standard_target_system_to_container
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Copies the target system to the container for systems
 *   with coefficients in standard double precision. */

static PyObject *py2c_copy_dobldobl_target_system_to_container
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Copies the target system to the container for systems
 *   with coefficients in double double precision. */

static PyObject *py2c_copy_quaddobl_target_system_to_container
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Copies the target system to the container for systems
 *   with coefficients in quad double precision. */

static PyObject *py2c_copy_multprec_target_system_to_container
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Copies the target system to the container for systems
 *   with coefficients in arbitrary multiprecision. */

static PyObject *py2c_copy_standard_container_to_target_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Copies the system in the container for systems with coefficients 
 *   in standard double precision to the target system. */

static PyObject *py2c_copy_dobldobl_container_to_target_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Copies the system in the container for systems with coefficients
 *   in double double precision to the target system. */

static PyObject *py2c_copy_quaddobl_container_to_target_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Copies the system in the container for systems with coefficients 
 *   in quad double precision to the target system. */

static PyObject *py2c_copy_multprec_container_to_target_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Copies the system in the container for systems with coefficients
 *   in arbitrary multiprecision to the target system. */

static PyObject *py2c_copy_start_system_to_container
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Copies the start system to the container for systems
 *   with coefficients in standard double precision. */

static PyObject *py2c_copy_dobldobl_start_system_to_container
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Copies the start system to the container for systems
 *   with coefficients in double double precision. */

static PyObject *py2c_copy_quaddobl_start_system_to_container
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Copies the start system to the container for systems
 *   with coefficients in quad double precision. */

static PyObject *py2c_copy_multprec_start_system_to_container
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Copies the start system to the container for systems
 *   with coefficients in arbitrary multiprecision. */

static PyObject *py2c_copy_standard_container_to_start_system 
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Copies the system in the container for systems with coefficients
 *   in standard double precision to the start system. */

static PyObject *py2c_copy_dobldobl_container_to_start_system 
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Copies the system in the container for systems with coefficients
 *   in double double precision to the start system. */

static PyObject *py2c_copy_quaddobl_container_to_start_system 
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Copies the system in the container for systems with coefficients
 *   in quad double precision to the start system. */

static PyObject *py2c_copy_multprec_container_to_start_system 
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Copies the system in the container for systems with coefficients
 *   in arbitrary multiprecision to the start system. */

/* Creation of homotopy and the tracking of all paths. */

static PyObject *py2c_create_standard_homotopy
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Initializes the data for a homotopy in standard double precision.
 *   The failure code is returned, which is zero when all goes well. */

static PyObject *py2c_create_standard_homotopy_with_gamma
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Initializes the data for a homotopy in standard double precision.
 *   On input are two doubles: the real and imaginary part of the
 *   gamma constant.
 *   The failure code is returned, which is zero when all goes well. */

static PyObject *py2c_create_dobldobl_homotopy
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Initializes the data for a homotopy in double double precision.
 *   The failure code is returned, which is zero when all goes well. */

static PyObject *py2c_create_dobldobl_homotopy_with_gamma
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Initializes the data for a homotopy in double double precision.
 *   On input are two doubles: the real and imaginary part of the
 *   gamma constant.
 *   The failure code is returned, which is zero when all goes well. */

static PyObject *py2c_create_quaddobl_homotopy
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Initializes the data for a homotopy in quad double precision.
 *   The failure code is returned, which is zero when all goes well. */

static PyObject *py2c_create_quaddobl_homotopy_with_gamma
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Initializes the data for a homotopy in quad double precision.
 *   On input are two doubles: the real and imaginary part of the
 *   gamma constant.
 *   The failure code is returned, which is zero when all goes well. */

static PyObject *py2c_create_multprec_homotopy
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Initializes the data for a homotopy in arbitrary multiprecision.
 *   The failure code is returned, which is zero when all goes well. */

static PyObject *py2c_create_multprec_homotopy_with_gamma
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Initializes the data for a homotopy in arbitrary multiprecision.
 *   On input are two doubles: the real and imaginary part of the
 *   gamma constant.
 *   The failure code is returned, which is zero when all goes well. */

static PyObject *py2c_clear_standard_homotopy
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Deallocation of the homotopy stored in standard double precision.
 *   On return is the failure code, which equals zero if all is well. */

static PyObject *py2c_clear_dobldobl_homotopy
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Deallocation of the homotopy stored in double double precision.
 *   On return is the failure code, which equals zero if all is well. */

static PyObject *py2c_clear_quaddobl_homotopy
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Deallocation of the homotopy stored in quad double precision.
 *   On return is the failure code, which equals zero if all is well. */

static PyObject *py2c_clear_multprec_homotopy
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Deallocation of the homotopy stored in arbitrary multiprecision.
 *   On return is the failure code, which equals zero if all is well. */

static PyObject *py2c_write_start_solutions ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Writes the start solutions in standard double precision either to
 *   the screen (standard output) or to the defined output file.
 *   On return is the failure code, which is zero if all is well.  */

static PyObject *py2c_tune_continuation_parameters
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Interactive procedure to tune the continuation parameters. */

static PyObject *py2c_show_continuation_parameters
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Shows the current values of the continuation parameters. */

static PyObject *py2c_autotune_continuation_parameters
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Tunes the values of the continuation parameters.
 *   On input are two integers:
 *   1) the difficulty level of the solution paths; and
 *   2) the number of decimal places in the precision. */

static PyObject *py2c_determine_output_during_continuation
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Interactive procedure to determine the level of output
 *   during the path tracking. */

static PyObject *py2c_solve_by_standard_homotopy_continuation
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Tracks the paths defined by the homotopy in standard double precision.
 *   On input is one integer: the number of tasks for path tracking.
 *   If that input number is zero, then no multitasking is applied.
 *   On return is the failure code, which is zero when all went well. */

static PyObject *py2c_solve_by_dobldobl_homotopy_continuation
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Tracks the paths defined by the homotopy in double double precision.
 *   On input is one integer: the number of tasks for path tracking.
 *   If that input number is zero, then no multitasking is applied.
 *   On return is the failure code, which is zero when all went well. */

static PyObject *py2c_solve_by_quaddobl_homotopy_continuation
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Tracks the paths defined by the homotopy in quad double precision.
 *   On input is one integer: the number of tasks for path tracking.
 *   If that input number is zero, then no multitasking is applied.
 *   On return is the failure code, which is zero when all went well. */

static PyObject *py2c_solve_by_multprec_homotopy_continuation
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Tracks the paths defined by the homotopy in arbitrary multiprecision.
 *   On input is one integer: the number of decimal places in the precision.
 *   On return is the failure code, which is zero when all went well. */

/* copying solutions from and to containers */

static PyObject *py2c_copy_standard_target_solutions_to_container
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Copies the target solutions in standard double precision to the
 *   container for solutions in standard double precision. */

static PyObject *py2c_copy_dobldobl_target_solutions_to_container
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Copies the target solutions in double double precision to the
 *   container for solutions in double double precision. */

static PyObject *py2c_copy_quaddobl_target_solutions_to_container
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Copies the target solutions in quad double precision to the
 *   container for solutions in quad double precision. */

static PyObject *py2c_copy_multprec_target_solutions_to_container
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Copies the target solutions in arbitrary multiprecision to the
 *   container for solutions in arbitrary multiprecision. */

static PyObject *py2c_copy_standard_container_to_target_solutions
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Copies the solutions in standard double precision from the
 *   container to the target solutions in standard double precision. */

static PyObject *py2c_copy_dobldobl_container_to_target_solutions
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Copies the solutions in double double precision from the
 *   container to the target solutions in double double precision. */

static PyObject *py2c_copy_quaddobl_container_to_target_solutions
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Copies the solutions in quad double precision from the
 *   container to the target solutions in quad double precision. */

static PyObject *py2c_copy_multprec_container_to_target_solutions
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Copies the solutions in arbitrary multiprecision from the
 *   container to the target solutions in arbitrary multiprecision. */

static PyObject *py2c_copy_start_solutions_to_container
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Copies the start solutions in standard double precision to the
 *   container for solutions in standard double precision. */

static PyObject *py2c_copy_dobldobl_start_solutions_to_container
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Copies the start solutions in double double precision to the
 *   container for solutions in double double precision. */

static PyObject *py2c_copy_quaddobl_start_solutions_to_container
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Copies the start solutions in quad double precision to the
 *   container for solutions in quad double precision. */

static PyObject *py2c_copy_multprec_start_solutions_to_container
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Copies the start solutions in arbitrary multiprecision to the
 *   container for solutions in arbitrary multiprecision. */

static PyObject *py2c_copy_standard_container_to_start_solutions
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Copies the solutions in standard double precision from the
 *   container to the start solutions in standard double precision. */

static PyObject *py2c_copy_dobldobl_container_to_start_solutions
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Copies the solutions in double double precision from the
 *   container to the start solutions in double double precision. */

static PyObject *py2c_copy_quaddobl_container_to_start_solutions
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Copies the solutions in quad double precision from the
 *   container to the start solutions in quad double precision. */

static PyObject *py2c_copy_multprec_container_to_start_solutions
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Copies the solutions in arbitrary multiprecision from the
 *   container to the start solutions in arbitrary multiprecision. */

/* black box solver, mixed volume calculator, and Newton step */

static PyObject *py2c_solve_system ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Calls the blackbox solver on the system stored in the container for
 *   systems with coefficients in standard double precision.
 *   One integer is expected on input: the number of tasks.
 *   If that number is zero, then no multitasking is applied.
 *   On return, the container for solutions in standard double precision
 *   contains the solutions to the system in the standard systems container. */

static PyObject *py2c_solve_Laurent_system ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Calls the blackbox solver on the system stored in the container for
 *   Laurent systems with coefficients in standard double precision.
 *   Two integers are expected on input:
 *   1) a boolean flag silent: if 1, then no intermediate output about
 *   the root counts is printed, if 0, then the solver is verbose; and 
 *   2) the number of tasks: if 0, then no multitasking is applied,
 *   otherwise as many tasks as the number will run.
 *   On return, the container for solutions in standard double precision
 *   contains the solutions to the system in the standard Laurent systems
 *   container. */

static PyObject *py2c_mixed_volume ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Computes the mixed volume, and the stable mixed volume as well if
 *   the input parameter equals 1.  On return is the mixed volume, or
 *   a tuple with the mixed volume and the stable mixed volume. */

static PyObject *py2c_standard_deflate ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Applies deflation in standard double precision to the system and
 *   the solutions stored in the containers.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_dobldobl_deflate ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Applies deflation in double double precision to the system and
 *   the solutions stored in the containers.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_quaddobl_deflate ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Applies deflation in quad double precision to the system and
 *   the solutions stored in the containers.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_standard_Newton_step ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Applies one Newton step in standard double precision to the system in
 *   the standard systems container and to the solutions in the container.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_dobldobl_Newton_step ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Applies one Newton step in double double precision to the system in
 *   the dobldobl systems container and to the solutions in the container.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_quaddobl_Newton_step ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Applies one Newton step in quad double precision to the system in
 *   the quaddobl systems container and to the solutions in the container.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_multprec_Newton_step ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Applies one Newton step in arbitrary multiprecision to the system in
 *   the multprec systems container and to the solutions in the container.
 *   On input is an integer, the number of decimal places in the precision.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_standard_Newton_Laurent_step
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Applies one Newton step in standard double precision to the Laurent
 *   system in the standard Laurent systems container and to the solutions
 *   in the container.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_dobldobl_Newton_Laurent_step
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Applies one Newton step in double double precision to the Laurent
 *   system in the dobldobl Laurent systems container and to the solutions
 *   in the container.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_quaddobl_Newton_Laurent_step
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Applies one Newton step in quad double precision to the Laurent
 *   system in the quaddobl Laurent systems container and to the solutions
 *   in the container.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_multprec_Newton_Laurent_step
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Applies one Newton step in arbitrary multiprecision to the Laurent
 *   system in the multprec Laurent systems container and to the solutions
 *   in the container.
 *   On input is an integer: the number of decimal places in the precision.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_varbprec_Newton_Laurent_steps
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Applies Newton's method in variable precision.
 *   There are six input parameters:
 *   1) the dimension: the number of variables and equations;
 *   2) the accuracy, expressed as the correct number of decimal places;
 *   3) the maximum number of iterations in Newton's method;
 *   4) an upper bound on the number of decimal places in the precision;
 *   5) a string, with the representation of the polynomials in the system.
 *   On return is the failure code, which equals zero if all went well. */

/* The wrapping of functions with prototypes in unisolvers.h starts here. */

static PyObject *py2c_usolve_standard ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Applies the method of Weierstrass to compute all roots of a
 *   polynomial in one variable with standard double precision arithmetic.
 *   On input are two numbers:
 *   1) the maximum number of iterations in the method of Weierstrass; and
 *   2) the epsilon requirement on the accuracy of the roots.
 *   Before calling this function, the polynomial should be stored in
 *   the standard systems container.  After the call of this function,
 *   the standard solutions container contains the roots of the polynomial.
 *   On return is the number of iterations done by the solver. */

static PyObject *py2c_usolve_dobldobl ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Applies the method of Weierstrass to compute all roots of a
 *   polynomial in one variable with double double precision arithmetic.
 *   On input are two numbers: 
 *   1) the maximum number of iterations in the method of Weierstrass; and
 *   2) the epsilon requirement on the accuracy of the roots.
 *   Before calling this function, the polynomial should be stored in
 *   the dobldobl systems container.  After the call of this function,
 *   the dobldobl solutions container contains the roots of the polynomial.
 *   On return is the number of iterations done by the solver. */

static PyObject *py2c_usolve_quaddobl ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Applies the method of Weierstrass to compute all roots of a
 *   polynomial in one variable with quad double precision arithmetic.
 *   On input are two numbers:
 *   1) the maximum number of iterations in the method of Weierstrass; and
 *   2) the epsilon requirement on the accuracy of the roots.
 *   Before calling this function, the polynomial should be stored in
 *   the quaddobl systems container.  After the call of this function,
 *   the quaddobl solutions container contains the roots of the polynomial.
 *   On return is the number of iterations done by the solver. */

static PyObject *py2c_usolve_multprec ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Applies the method of Weierstrass to compute all roots of a
 *   polynomial in one variable with arbitrary multiprecision arithmetic.
 *   On input are three numbers: 
 *   1) the number of decimal places in the working precision; 
 *   2) the maximum number of iterations in the method of Weierstrass; and
 *   3) the epsilon requirement on the accuracy of the roots.
 *   Before calling this function, the polynomial should be stored in
 *   the multprec systems container.  After the call of this function,
 *   the multprec solutions container contains the roots of the polynomial.
 *   On return is the number of iterations done by the solver. */

/* The wrapping of functions with prototypes in giftwrappers.h starts here. */

static PyObject *py2c_giftwrap_planar ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Applies the giftwrapping algorithm to a planar point configuration.
 *   On input are an integer and a string:
 *   1) the number of points in the list;
 *   2) the string representation of a Python list of tuples.
 *   On return is the string representation of the vertex points,
 *   sorted so that each two consecutive points define an edge. */

static PyObject *py2c_giftwrap_convex_hull ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Applies the giftwrapping algorithm to a point configuration.
 *   On input are an integer and a string:
 *   1) the number of points in the list;
 *   2) the string representation of a Python list of tuples.
 *   When the function returns, the internal data structures
 *   to store the convex hull are defined.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_giftwrap_number_of_facets
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the number of facets of the given dimension.
 *   On input is an integer, the dimension of the facet. */

static PyObject *py2c_giftwrap_retrieve_facet
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the string representation of a facet.
 *   On input are two integer numbers:
 *   1) the dimension of the facet;
 *   2) the index of the facet. */

static PyObject *py2c_giftwrap_clear_3d_facets
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Deallocates the data for a convex hull in dimension three. */

static PyObject *py2c_giftwrap_clear_4d_facets
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Deallocates the data for a convex hull in dimension four. */

static PyObject *py2c_giftwrap_support_size
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the number of characters in the string representation of
 *   the support of the first Laurent polynomial in the container. */

static PyObject *py2c_giftwrap_support_string
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the string representation of the support
 *   of a Laurent polynomial. */

static PyObject *py2c_giftwrap_clear_support_string
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Deallocates the string representation of the support set
 *   that was stored internally by the call py2c_giftwrap_support_size. */

static PyObject *py2c_giftwrap_initial_form ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Replaces the system in the Laurent systems container by its initial form.
 *   There are three input parameters:
 *   1) the dimension, number of coordinates in the inner normal;
 *   2) the number of characters in the string representation for the normal;
 *   3) the string representation of the inner normal.
 *   On return is the failure code, which equals zero if all went well. */

/* wrapping functions in syscon.h starts from here */

static PyObject *py2c_syscon_read_standard_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Interactive procedure to read a polynomial system with coefficients
 *   in standard double precision.
 *   The system will be placed in the standard systems container.
 *   The failure code is returned, which equals zero if all went well. */

static PyObject *py2c_syscon_read_standard_Laurent_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Interactive procedure to read a Laurent polynomial system with
 *   coefficients in standard double precision.
 *   The system will be placed in the standard Laurent systems container.
 *   The failure code is returned, which equals zero if all went well. */

static PyObject *py2c_syscon_read_dobldobl_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Interactive procedure to read a polynomial system with coefficients
 *   in double double precision.
 *   The system will be placed in the dobldobl systems container.
 *   The failure code is returned, which equals zero if all went well. */

static PyObject *py2c_syscon_read_dobldobl_Laurent_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Interactive procedure to read a Laurent polynomial system with
 *   coefficients in double double precision.
 *   The system will be placed in the dobldobl Laurent systems container.
 *   The failure code is returned, which equals zero if all went well. */

static PyObject *py2c_syscon_read_quaddobl_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Interactive procedure to read a polynomial system with coefficients
 *   in quad double precision.
 *   The system will be placed in the quaddobl systems container.
 *   The failure code is returned, which equals zero if all went well. */

static PyObject *py2c_syscon_read_quaddobl_Laurent_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Interactive procedure to read a Laurent polynomial system with
 *   coefficients in quad double precision.
 *   The system will be placed in the quaddobl Laurent systems container.
 *   The failure code is returned, which equals zero if all went well. */

static PyObject *py2c_syscon_read_multprec_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Interactive procedure to read a polynomial system with coefficients
 *   in arbitrary multiprecision.  The one input parameter is an integer,
 *   the number of decimal places in the working precision.
 *   The system will be placed in the multprec systems container.
 *   The failure code is returned, which equals zero if all went well. */

static PyObject *py2c_syscon_read_multprec_Laurent_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Interactive procedure to read a Laurent polynomial system with
 *   coefficients in arbitrary multiprecision.  The one input parameter is
 *   an integer, the number of decimal places in the working precision.
 *   The system will be placed in the multprec Laurent systems container.
 *   The failure code is returned, which equals zero if all went well. */

static PyObject *py2c_syscon_random_system ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Places in the systems container a random polynomial system
 *   with coefficients in standard double precision.
 *   There are four integers as input parameters:
 *   1) n, the number of polynomials and variables;
 *   2) m, the number of monomials per equation;
 *   3) d, the largest degree of each monomial;
 *   4) c, the type of coefficient: 0 if on the complex unit circle,
 *   1, if all coefficients are one, 2, if all coefficients are
 *   random floats in [-1,+1]. */

static PyObject *py2c_syscon_write_standard_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Writes the polynomial system with standard double precision coefficients
 *   that is stored in the container. */

static PyObject *py2c_syscon_write_standard_Laurent_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Writes the Laurent polynomial system with standard double precision
 *   coefficients that is stored in the container. */

static PyObject *py2c_syscon_write_dobldobl_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Writes the polynomial system with double double precision coefficients
 *   that is stored in the container. */

static PyObject *py2c_syscon_write_dobldobl_Laurent_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Writes the Laurent polynomial system with double double precision
 *   coefficients that is stored in the container. */

static PyObject *py2c_syscon_write_quaddobl_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Writes the polynomial system with quad double precision coefficients
 *   that is stored in the container. */

static PyObject *py2c_syscon_write_quaddobl_Laurent_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Writes the Laurent polynomial system with quad double precision
 *   coefficients that is stored in the container. */

static PyObject *py2c_syscon_write_multprec_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Writes the polynomial system with arbitrary multiprecision coefficients
 *   that is stored in the container. */

static PyObject *py2c_syscon_write_multprec_Laurent_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Writes the Laurent polynomial system with arbitrary multiprecision
 *   coefficients that is stored in the container. */

static PyObject *py2c_syscon_clear_standard_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Deallocates the container for polynomial systems
 *   with coefficients in standard double precision. */

static PyObject *py2c_syscon_clear_standard_Laurent_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Deallocates the container for Laurent polynomial systems
 *   with coefficients in standard double precision. */

static PyObject *py2c_syscon_clear_dobldobl_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Deallocates the container for polynomial systems
 *   with coefficients in double double precision. */

static PyObject *py2c_syscon_clear_dobldobl_Laurent_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Deallocates the container for Laurent polynomial systems
 *   with coefficients in double double precision. */

static PyObject *py2c_syscon_clear_quaddobl_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Deallocates the container for polynomial systems
 *   with coefficients in quad double precision. */

static PyObject *py2c_syscon_clear_quaddobl_Laurent_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Deallocates the container for Laurent polynomial systems
 *   with coefficients in quad double precision. */

static PyObject *py2c_syscon_clear_multprec_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Deallocates the container for polynomial systems
 *   with coefficients in arbitrary multiprecision. */

static PyObject *py2c_syscon_clear_multprec_Laurent_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Deallocates the container for Laurent polynomial systems
 *   with coefficients in arbitrary multiprecision. */

static PyObject *py2c_syscon_number_of_symbols
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the number of symbols in the symbol table. */

static PyObject *py2c_syscon_write_symbols
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Writes the symbols in the symbol table to screen.
 *   Returns the failure code, which equals zero if all went well. */

static PyObject *py2c_syscon_string_of_symbols
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns a string that contains the symbols in the symbol table.
 *   The symbols are separate from each other by one space. */

static PyObject *py2c_syscon_remove_symbol_name
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Removes a symbol, given by name, from the symbol table.
 *   On input are two arguments:
 *   1) an integer, as the number of characters in the name;
 *   2) a string of characters with the name of the symbol.
 *   The failure code is returned, which equals zero when all went well. */

static PyObject *py2c_syscon_clear_symbol_table
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Clears the symbol table. */

static PyObject *py2c_syscon_number_of_standard_polynomials
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the number of polynomials with coefficients in standard
 *   double precision as stored in the systems container. */

static PyObject *py2c_syscon_number_of_dobldobl_polynomials
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the number of polynomials with coefficients in double
 *   double precision as stored in the systems container. */

static PyObject *py2c_syscon_number_of_quaddobl_polynomials
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the number of polynomials with coefficients in quad
 *   double precision as stored in the systems container. */

static PyObject *py2c_syscon_number_of_multprec_polynomials
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the number of polynomials with coefficients in arbitrary
 *   multiprecision as stored in the systems container. */

static PyObject *py2c_syscon_number_of_standard_Laurentials
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the number of Laurent polynomials with coefficients in
 *   standard double precision as stored in the systems container. */

static PyObject *py2c_syscon_number_of_dobldobl_Laurentials
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the number of Laurent polynomials with coefficients in
 *   double double precision as stored in the systems container. */

static PyObject *py2c_syscon_number_of_quaddobl_Laurentials
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the number of Laurent polynomials with coefficients in
 *   quad double precision as stored in the systems container. */

static PyObject *py2c_syscon_number_of_multprec_Laurentials
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the number of Laurent polynomials with coefficients in
 *   arbitrary multiprecision as stored in the systems container. */

static PyObject *py2c_syscon_initialize_number_of_standard_polynomials
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Initializes the container for polynomials with coefficients in
 *   standard double precision.  The input argument is an integer,
 *   the number of polynomials in the container.
 *   The failure code is returned, which equals zero if all went well. */

static PyObject *py2c_syscon_initialize_number_of_dobldobl_polynomials
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Initializes the container for polynomials with coefficients in
 *   double double precision.  The input argument is an integer,
 *   the number of polynomials in the container.
 *   The failure code is returned, which equals zero if all went well. */

static PyObject *py2c_syscon_initialize_number_of_quaddobl_polynomials
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Initializes the container for polynomials with coefficients in
 *   quad double precision.  The input argument is an integer,
 *   the number of polynomials in the container.
 *   The failure code is returned, which equals zero if all went well. */

static PyObject *py2c_syscon_initialize_number_of_multprec_polynomials
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Initializes the container for polynomials with coefficients in
 *   arbitrary multiprecision.  The input argument is an integer,
 *   the number of polynomials in the container.
 *   The failure code is returned, which equals zero if all went well. */

static PyObject *py2c_syscon_initialize_number_of_standard_Laurentials
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Initializes the container for Laurent polynomials with coefficients
 *   in standard double precision.  The input argument is an integer,
 *   the number of polynomials in the container.
 *   The failure code is returned, which equals zero if all went well. */

static PyObject *py2c_syscon_initialize_number_of_dobldobl_Laurentials
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Initializes the container for Laurent polynomials with coefficients
 *   in double double precision.  The input argument is an integer,
 *   the number of polynomials in the container.
 *   The failure code is returned, which equals zero if all went well. */

static PyObject *py2c_syscon_initialize_number_of_quaddobl_Laurentials
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Initializes the container for Laurent polynomials with coefficients
 *   in quad double precision.  The input argument is an integer,
 *   the number of polynomials in the container.
 *   The failure code is returned, which equals zero if all went well. */

static PyObject *py2c_syscon_initialize_number_of_multprec_Laurentials
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Initializes the container for Laurent polynomials with coefficients
 *   in arbitrary multiprecision.  The input argument is an integer,
 *   the number of polynomials in the container.
 *   The failure code is returned, which equals zero if all went well. */

static PyObject *py2c_syscon_degree_of_standard_polynomial
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the degree of the k-th polynomial in the container for
 *   polynomials with coefficients in standard double precision.
 *   The index k of the polynomial is the one input argument. */

static PyObject *py2c_syscon_degree_of_dobldobl_polynomial
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the degree of the k-th polynomial in the container for
 *   polynomials with coefficients in double double precision.
 *   The index k of the polynomial is the one input argument. */

static PyObject *py2c_syscon_degree_of_quaddobl_polynomial
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the degree of the k-th polynomial in the container for
 *   polynomials with coefficients in quad double precision.
 *   The index k of the polynomial is the one input argument. */

static PyObject *py2c_syscon_degree_of_multprec_polynomial
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the degree of the k-th polynomial in the container for
 *   polynomials with coefficients in arbitrary multiprecision.
 *   The index k of the polynomial is the one input argument. */

static PyObject *py2c_syscon_number_of_terms
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the number of terms in the k-th polynomial stored in the
 *   container for systems with coefficients in standard double precision.
 *   The input parameter k is the index of the polynomial k. */

static PyObject *py2c_syscon_number_of_Laurent_terms
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the number of terms in the k-th Laurent polynomial stored
 *   in the container for Laurent polynomials systems with coefficients
 *   in standard double precision.
 *   The input parameter k is the index of the polynomial k. */

static PyObject *py2c_syscon_retrieve_term ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Retrieves one term of a polynomial with coefficients in standard
 *   double precision, that is stored in the systems container.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_syscon_store_standard_polynomial
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Defines the k-th polynomial in the systems container for polynomials
 *   with coefficients in standard double precision.
 *   As a precondition for this function, the container must be initialized
 *   for sufficiently many polynomials, in any case >= k.
 *   There are four input parameters, three integers and one string:
 *   1) nc, the number of characters in the string p;
 *   2) n, the number of variables in the multivariate polynomial;
 *   3) k, the index of the polynomial in the system;
 *   4) p, a valid string representation for a polynomial.
 *   On return is the failure code, which equals zero when all went well. */

static PyObject *py2c_syscon_store_dobldobl_polynomial
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Defines the k-th polynomial in the systems container for polynomials
 *   with coefficients in double double precision.
 *   As a precondition for this function, the container must be initialized
 *   for sufficiently many polynomials, in any case >= k.
 *   There are four input parameters, three integers and one string:
 *   1) nc, the number of characters in the string p;
 *   2) n, the number of variables in the multivariate polynomial;
 *   3) k, the index of the polynomial in the system;
 *   4) p, a valid string representation for a polynomial.
 *   On return is the failure code, which equals zero when all went well. */

static PyObject *py2c_syscon_store_quaddobl_polynomial
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Defines the k-th polynomial in the systems container for polynomials
 *   with coefficients in quad double precision.
 *   As a precondition for this function, the container must be initialized
 *   for sufficiently many polynomials, in any case >= k.
 *   There are four input parameters, three integers and one string:
 *   1) nc, the number of characters in the string p;
 *   2) n, the number of variables in the multivariate polynomial;
 *   3) k, the index of the polynomial in the system;
 *   4) p, a valid string representation for a polynomial.
 *   On return is the failure code, which equals zero when all went well. */

static PyObject *py2c_syscon_store_multprec_polynomial
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Defines the k-th polynomial in the systems container for polynomials
 *   with coefficients in arbitrary multiprecision.
 *   As a precondition for this function, the container must be initialized
 *   for sufficiently many polynomials, in any case >= k.
 *   There are five input parameters, four integers and one string:
 *   1) nc, the number of characters in the string p;
 *   2) n, the number of variables in the multivariate polynomial;
 *   3) k, the index of the polynomial in the system;
 *   4) dp, the number of decimal places to parse the coefficients;
 *   5) p, a valid string representation for a polynomial.
 *   On return is the failure code, which equals zero when all went well. */

static PyObject *py2c_syscon_load_standard_polynomial
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the k-th polynomial in the systems container 
 *   with standard double complex coefficients as a string.
 *   The value for k is in the one integer parameter of this function. */

static PyObject *py2c_syscon_load_dobldobl_polynomial
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the k-th polynomial in the systems container 
 *   with double double complex coefficients as a string.
 *   The value for k is in the one integer parameter of this function. */

static PyObject *py2c_syscon_load_quaddobl_polynomial
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the k-th polynomial in the systems container 
 *   with quad double complex coefficients as a string.
 *   The value for k is in the one integer parameter of this function. */

static PyObject *py2c_syscon_load_multprec_polynomial
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the k-th polynomial in the systems container 
 *   with arbitrary multiprecision complex coefficients as a string.
 *   The value for k is in the one integer parameter of this function. */

static PyObject *py2c_syscon_store_standard_Laurential
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Defines the k-th polynomial in the systems container for Laurent
 *   polynomials with coefficients in standard double precision.
 *   As a precondition for this function, the container must be initialized
 *   for sufficiently many polynomials, in any case >= k.
 *   There are four input parameters, three integers and one string:
 *   1) nc, the number of characters in the string p;
 *   2) n, the number of variables in the multivariate polynomial;
 *   3) k, the index of the polynomial in the system;
 *   4) p, a valid string representation for a polynomial.
 *   On return is the failure code, which equals zero when all went well. */

static PyObject *py2c_syscon_store_dobldobl_Laurential
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Defines the k-th polynomial in the systems container for Laurent
 *   polynomials with coefficients in double double precision.
 *   As a precondition for this function, the container must be initialized
 *   for sufficiently many polynomials, in any case >= k.
 *   There are four input parameters, three integers and one string:
 *   1) nc, the number of characters in the string p;
 *   2) n, the number of variables in the multivariate polynomial;
 *   3) k, the index of the polynomial in the system;
 *   4) p, a valid string representation for a polynomial.
 *   On return is the failure code, which equals zero when all went well. */

static PyObject *py2c_syscon_store_quaddobl_Laurential
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Defines the k-th polynomial in the systems container for Laurent
 *   polynomials with coefficients in quad double precision.
 *   As a precondition for this function, the container must be initialized
 *   for sufficiently many polynomials, in any case >= k.
 *   There are four input parameters, three integers and one string:
 *   1) nc, the number of characters in the string p;
 *   2) n, the number of variables in the multivariate polynomial;
 *   3) k, the index of the polynomial in the system;
 *   4) p, a valid string representation for a polynomial.
 *   On return is the failure code, which equals zero when all went well. */

static PyObject *py2c_syscon_store_multprec_Laurential
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Defines the k-th polynomial in the systems container for Laurent
 *   polynomials with coefficients in arbitrary multiprecision.
 *   As a precondition for this function, the container must be initialized
 *   for sufficiently many polynomials, in any case >= k.
 *   There are five input parameters, four integers and one string:
 *   1) nc, the number of characters in the string p;
 *   2) n, the number of variables in the multivariate polynomial;
 *   3) k, the index of the polynomial in the system;
 *   4) dp, the number of decimal places to parse the coefficients;
 *   5) p, a valid string representation for a polynomial.
 *   On return is the failure code, which equals zero when all went well. */

static PyObject *py2c_syscon_load_standard_Laurential
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the k-th polynomial in the Laurent systems container 
 *   with standard double complex coefficients as a string.
 *   The value for k is in the one integer parameter of this function. */

static PyObject *py2c_syscon_load_dobldobl_Laurential
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the k-th polynomial in the Laurent systems container 
 *   with double double complex coefficients as a string.
 *   The value for k is in the one integer parameter of this function. */

static PyObject *py2c_syscon_load_quaddobl_Laurential
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the k-th polynomial in the Laurent systems container 
 *   with quad double complex coefficients as a string.
 *   The value for k is in the one integer parameter of this function. */

static PyObject *py2c_syscon_load_multprec_Laurential
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the k-th polynomial in the Laurent systems container 
 *   with arbitrary multiprecision complex coefficients as a string.
 *   The value for k is in the one integer parameter of this function. */

static PyObject *py2c_syscon_total_degree ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns in d the total degree of the system with coefficients in
 *   standard double precision, as stored in the container. */

static PyObject *py2c_syscon_standard_drop_variable_by_index
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Replaces the system in the standard double precision container 
 *   with the same system that has its k-th variable dropped.
 *   The index k of the variable is given as an input parameter.
 *   On return is the failure code, which equals zero if all went well.  */

static PyObject *py2c_syson_standard_drop_variable_by_name
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Replaces the system in the standard double precision container 
 *   with the same system that have that variable dropped
 *   corresponding to the name in the string s of nc characters long.
 *   The function has two input parameters, an integer and a string:
 *   1) nc, the number of characters in the string with the name;
 *   2) s, a string that holds the name of the variable.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_syscon_dobldobl_drop_variable_by_index
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Replaces the system in the double double precision container 
 *   with the same system that has its k-th variable dropped.
 *   The index k of the variable is given as an input parameter.
 *   On return is the failure code, which equals zero if all went well.  */

static PyObject *py2c_syscon_dobldobl_drop_variable_by_name
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Replaces the system in the double double precision container 
 *   with the same system that have that variable dropped
 *   corresponding to the name in the string s of nc characters long.
 *   The function has two input parameters, an integer and a string:
 *   1) nc, the number of characters in the string with the name;
 *   2) s, a string that holds the name of the variable.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_syscon_quaddobl_drop_variable_by_index
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Replaces the system in the quad double precision container 
 *   with the same system that has its k-th variable dropped.
 *   The index k of the variable is given as an input parameter.
 *   On return is the failure code, which equals zero if all went well.  */

static PyObject *py2c_syscon_quaddobl_drop_variable_by_name
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Replaces the system in the quad double precision container 
 *   with the same system that have that variable dropped
 *   corresponding to the name in the string s of nc characters long.
 *   The function has two input parameters, an integer and a string:
 *   1) nc, the number of characters in the string with the name;
 *   2) s, a string that holds the name of the variable.
 *   On return is the failure code, which equals zero if all went well. */

/* The wrapping of functions with prototypes in solcon.h starts here. */

static PyObject *py2c_solcon_read_standard_solutions
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Interactive function to read the solutions into the container,
 *   in standard double precision.
 *   Returns the failure code, which is zero when all went well. */

static PyObject *py2c_solcon_read_dobldobl_solutions
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Interactive function to read the solutions into the container,
 *   in double double precision.
 *   Returns the failure code, which is zero when all went well. */

static PyObject *py2c_solcon_read_quaddobl_solutions
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Interactive function to read the solutions into the container,
 *   in quad double precision.
 *   Returns the failure code, which is zero when all went well. */

static PyObject *py2c_solcon_read_multprec_solutions
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Interactive function to read the solutions into the container,
 *   in arbitrary multiprecision.
 *   Returns the failure code, which is zero when all went well. */

static PyObject *py2c_solcon_write_standard_solutions
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Writes the solutions in standard double precision to screen.
 *   Returns the failure code, which equals zero when all went well. */

static PyObject *py2c_solcon_write_dobldobl_solutions
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Writes the solutions in double double precision to screen.
 *   Returns the failure code, which equals zero when all went well. */

static PyObject *py2c_solcon_write_quaddobl_solutions
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Writes the solutions in quad double precision to screen.
 *   Returns the failure code, which equals zero when all went well. */

static PyObject *py2c_solcon_write_multprec_solutions
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Writes the solutions in arbitrary multiprecision to screen.
 *   Returns the failure code, which equals zero when all went well. */

static PyObject *py2c_solcon_clear_standard_solutions
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Deallocates the container for solutions in standard double precision.
 *   Returns the failure code, which equals zero when all went well. */

static PyObject *py2c_solcon_clear_dobldobl_solutions
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Deallocates the container for solutions in double double precision.
 *   Returns the failure code, which equals zero when all went well. */

static PyObject *py2c_solcon_clear_quaddobl_solutions
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Deallocates the container for solutions in quad double precision.
 *   Returns the failure code, which equals zero when all went well. */

static PyObject *py2c_solcon_clear_multprec_solutions
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Deallocates the container for solutions in arbitrary multiprecision.
 *   Returns the failure code, which equals zero when all went well. */

static PyObject *py2c_solcon_open_solution_input_file
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Prompts the user for the name of the input file for the solutions and
 *   opens the input file.  All subsequent reading happens from this input.
 *   Returns the failure code, which equals zero when all went well. */

static PyObject *py2c_solcon_length_standard_solution_string
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   On input is the index k to a solution in standard double precision,
 *   stored in the container.  On return is the length of the string 
 *   representation for that k-th solution in the container. */

static PyObject *py2c_solcon_length_dobldobl_solution_string
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   On input is the index k to a solution in double double precision,
 *   stored in the container.  On return is the length of the string 
 *   representation for that k-th solution in the container. */

static PyObject *py2c_solcon_length_quaddobl_solution_string
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   On input is the index k to a solution in quad double precision,
 *   stored in the container.  On return is the length of the string 
 *   representation for that k-th solution in the container. */

static PyObject *py2c_solcon_length_multprec_solution_string
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   On input is the index k to a solution in arbitrary multiprecision,
 *   stored in the container.  On return is the length of the string 
 *   representation for that k-th solution in the container. */

static PyObject *py2c_solcon_write_standard_solution_string
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the string representation for the k-th solution stored
 *   in standard double precision in the container.
 *   On input are two integers:
 *   1) the index to the solution; and
 *   2) the number of characters in the string representation
 *   for that solution. */

static PyObject *py2c_solcon_write_dobldobl_solution_string
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the string representation for the k-th solution stored
 *   in double double precision in the container.
 *   On input are two integers:
 *   1) the index to the solution; and
 *   2) the number of characters in the string representation
 *   for that solution. */

static PyObject *py2c_solcon_write_quaddobl_solution_string
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the string representation for the k-th solution stored
 *   in quad double precision in the container.
 *   On input are two integers:
 *   1) the index to the solution; and
 *   2) the number of characters in the string representation
 *   for that solution. */

static PyObject *py2c_solcon_write_multprec_solution_string
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the string representation for the k-th solution stored
 *   in arbitrary multiprecision in the container.
 *   On input are two integers:
 *   1) the index to the solution; and
 *   2) the number of characters in the string representation
 *   for that solution. */

static PyObject *py2c_solcon_retrieve_next_standard_initialize
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Resets the pointer to the current standard solution in the container
 *   to the first solution in the list.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_solcon_retrieve_next_dobldobl_initialize
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Resets the pointer to the current dobldobl solution in the container
 *   to the first solution in the list.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_solcon_retrieve_next_quaddobl_initialize
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Resets the pointer to the current quaddobl solution in the container
 *   to the first solution in the list.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_solcon_retrieve_next_multprec_initialize
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Resets the pointer to the current multprec solution in the container
 *   to the first solution in the list.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_solcon_move_current_standard_to_next
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Moves the pointer to the current solution in standard double precision
 *   to the next solution and returns the value of the cursor.
 *   If cursor on return is zero, then either the pointer was null
 *   or there is no next solution. */

static PyObject *py2c_solcon_move_current_dobldobl_to_next
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Moves the pointer to the current solution in double double precision
 *   to the next solution and returns the value of the cursor.
 *   If cursor on return is zero, then either the pointer was null
 *   or there is no next solution. */

static PyObject *py2c_solcon_move_current_quaddobl_to_next
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Moves the pointer to the current solution in double double precision
 *   to the next solution and returns the value of the cursor.
 *   If cursor on return is zero, then either the pointer was null
 *   or there is no next solution. */

static PyObject *py2c_solcon_move_current_multprec_to_next
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Moves the pointer to the current solution in arbitrary multiprecision
 *   to the next solution and returns the value of the cursor.
 *   If cursor on return is zero, then either the pointer was null
 *   or there is no next solution. */

static PyObject *py2c_solcon_length_current_standard_solution_string
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the number of characters in the string representation
 *   of the current standard double solution in the container,
 *   at the place indicated by the value of the cursor.
 *   If this value equals zero, then there is no current solution,
 *   and then the length on return equals zero. */

static PyObject *py2c_solcon_length_current_dobldobl_solution_string
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the number of characters in the string representation
 *   of the current double double solution in the container,
 *   at the place indicated by the value of the cursor.
 *   If this value equals zero, then there is no current solution,
 *   and then the length on return equals zero. */

static PyObject *py2c_solcon_length_current_quaddobl_solution_string
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the number of characters in the string representation
 *   of the current quad double solution in the container,
 *   at the place indicated by the value of the cursor.
 *   If this value equals zero, then there is no current solution,
 *   and then the length on return equals zero. */

static PyObject *py2c_solcon_length_current_multprec_solution_string
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the number of characters in the string representation
 *   of the current arbitrary multiprecision solution in the container,
 *   at the place indicated by the value of the cursor.
 *   If this value equals zero, then there is no current solution,
 *   and then the length on return equals zero. */

static PyObject *py2c_solcon_write_current_standard_solution_string
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Writes the current standard double solution in the solution container
 *   to the string s of n+1 characters.  The last character is \0.
 *   The value of n is given as the one input parameter to this function.
 *   On return is the string that contains the string representation of
 *   the current solution in standard double precision in the container. */

static PyObject *py2c_solcon_write_current_dobldobl_solution_string
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Writes the current double double solution in the solution container
 *   to the string s of n+1 characters.  The last character is \0.
 *   The value of n is given as the one input parameter to this function.
 *   On return is the string that contains the string representation of
 *   the current solution in double double precision in the container. */

static PyObject *py2c_solcon_write_current_quaddobl_solution_string
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Writes the current quad double solution in the solution container
 *   to the string s of n+1 characters.  The last character is \0.
 *   The value of n is given as the one input parameter to this function.
 *   On return is the string that contains the string representation of
 *   the current solution in quad double precision in the container. */

static PyObject *py2c_solcon_write_current_multprec_solution_string
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Writes the current multiprecision solution in the solution container
 *   to the string s of n+1 characters.  The last character is \0.
 *   The value of n is given as the one input parameter to this function.
 *   On return is the string that contains the string representation of
 *   the current solution in arbitrary multiprecision in the container. */

static PyObject *py2c_solcon_append_standard_solution_string
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Appends a solution in standard double precision to the list
 *   of solutions already stored in the container.
 *   There are three input parameters:
 *   1) the number of variables;
 *   2) the number of characters in the string;
 *   3) the string representing the solution to append to the list.
 *   Returns the failure code, which equals zero if all went well. */

static PyObject *py2c_solcon_append_dobldobl_solution_string
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Appends a solution in double double precision to the list
 *   of solutions already stored in the container.
 *   There are three input parameters:
 *   1) the number of variables;
 *   2) the number of characters in the string;
 *   3) the string representing the solution to append to the list.
 *   Returns the failure code, which equals zero if all went well. */

static PyObject *py2c_solcon_append_quaddobl_solution_string
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Appends a solution in quad double precision to the list
 *   of solutions already stored in the container.
 *   There are three input parameters:
 *   1) the number of variables;
 *   2) the number of characters in the string;
 *   3) the string representing the solution to append to the list.
 *   Returns the failure code, which equals zero if all went well. */

static PyObject *py2c_solcon_append_multprec_solution_string
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Appends a solution in arbitrary multiprecision to the list
 *   of solutions already stored in the container.
 *   There are three input parameters:
 *   1) the number of variables;
 *   2) the number of characters in the string;
 *   3) the string representing the solution to append to the list.
 *   Returns the failure code, which equals zero if all went well. */

static PyObject *py2c_solcon_number_of_standard_solutions
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the number of solutions in standard double precision,
 *   as stored in the container. */

static PyObject *py2c_solcon_number_of_dobldobl_solutions
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the number of solutions in double double precision,
 *   as stored in the container. */

static PyObject *py2c_solcon_number_of_quaddobl_solutions
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the number of solutions in quad double precision,
 *   as stored in the container. */

static PyObject *py2c_solcon_number_of_multprec_solutions
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the number of solutions in arbitrary multiprecision,
 *   as stored in the container. */

static PyObject *py2c_solcon_standard_drop_coordinate_by_index
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Replaces the solutions in the standard double precision container 
 *   with the same solutions that have their k-th coordinate dropped.
 *   There is one input parameter: the index k of the coordinate.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_solcon_standard_drop_coordinate_by_name
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Replaces the solutions in the standard double precision container 
 *   with the same solutions that have their coordinate dropped
 *   corresponding to the name in the string s of nc characters long.
 *   There are two input parameters, an integer and a string:
 *   1) nc, the number of characters in the string with the name;
 *   2) s, the string with the name of the variable.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_solcon_dobldobl_drop_coordinate_by_index
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Replaces the solutions in the double double precision container 
 *   with the same solutions that have their k-th coordinate dropped.
 *   There is one input parameter: the index k of the coordinate.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_solcon_dobldobl_drop_coordinate_by_name
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Replaces the solutions in the double double precision container 
 *   with the same solutions that have their coordinate dropped
 *   corresponding to the name in the string s of nc characters long.
 *   There are two input parameters, an integer and a string:
 *   1) nc, the number of characters in the string with the name;
 *   2) s, the string with the name of the variable.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_solcon_quaddobl_drop_coordinate_by_index
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Replaces the solutions in the quad double precision container 
 *   with the same solutions that have their k-th coordinate dropped.
 *   There is one input parameter: the index k of the coordinate.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_solcon_quaddobl_drop_coordinate_by_name
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Replaces the solutions in the quad double precision container 
 *   with the same solutions that have their coordinate dropped
 *   corresponding to the name in the string s of nc characters long.
 *   There are two input parameters, an integer and a string:
 *   1) nc, the number of characters in the string with the name;
 *   2) s, the string with the name of the variable.
 *   On return is the failure code, which equals zero if all went well. */

/* The wrapping of the functions in product.h starts here. */

static PyObject *py2c_product_supporting_set_structure 
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Builds a supporting set structure for the system stored in the
 *   container with coefficients in standard double precision. */

static PyObject *py2c_product_write_set_structure
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Writes the supporting set structure to screen. */

static PyObject *py2c_product_set_structure_string
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the string representation of the set structure. */

static PyObject *py2c_product_parse_set_structure
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Parses a given string into a set structure.
 *   On input are two parameters, one integer and one string:
 *   1) the number of characters in the given string; and
 *   2) the characters in the string.
 *   On return is the failure code, if zero, then the string
 *   has been parsed into a valid set structure. */

static PyObject *py2c_product_is_set_structure_supporting
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Checks whether the stored set structure is supporting
 *   for the system in the standard systems container.
 *   Returns an integer which represents true (1) or false (0). */

static PyObject *py2c_product_linear_product_root_count
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the linear-product root count, computed from
 *   the supporting set structure. */

static PyObject *py2c_product_random_linear_product_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Builds a random linear-product system based on the
 *   stored set structure.   On return is the failure code,
 *   which equals zero if all went well. */

static PyObject *py2c_product_solve_linear_product_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Computes all solutions to the linear-product system
 *   and stores the solutions in the container for solutions
 *   in standard double precision.  On return is the failure
 *   code, which equals zero if all went well. */

static PyObject *py2c_product_clear_set_structure
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Deallocates the set structure. */

static PyObject *py2c_product_m_homogeneous_Bezout_number
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   For the system in the standard systems container,
 *   a heuristic partition of the set of variables may
 *   lead to a Bezout number that is smaller than the total degree.
 *   On return is the m-homogeneous Bezout number for the
 *   string representation of the partition that is returned
 *   as the second argument in the tuple. */

static PyObject *py2c_product_m_partition_Bezout_number
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Given a partition of the set of variables, computes 
 *   the m-homogeneous Bezout number for the system in
 *   the standard systems container.
 *   On input are two arguments:
 *   1) the number of characters in the string (second argument); and
 *   2) the string representation for a partition of the variables.
 *   On return is the m-homogeneous Bezout number. */

static PyObject *py2c_product_m_homogeneous_start_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Given a partition of the set of variables, constructs
 *   an m-homogeneous Bezout number for the system in
 *   the standard systems container.
 *   On input are two arguments:
 *   1) the number of characters in the string (second argument); and
 *   2) the string representation for a partition of the variables.
 *   On return is the m-homogeneous Bezout number. */

/* The wrapping of functions with prototypes in celcon.h starts here. */

static PyObject *py2c_celcon_initialize_supports
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Initializes the cell container with the number of distinct supports,
 *   this number is given as the one input parameter.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_celcon_set_type_of_mixture
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Defines the type of mixture of the support sets.
 *   On input are two parameters, an integer and a string:
 *   1) the integer equals the number of distinct supports;
 *   2) the string is a string representation of a Python list of integers,
 *   there are as many integers as the value of the first parameter.
 *   Each integer is a positive number, equal to the number of occurrences
 *   of each support set. */

static PyObject *py2c_celcon_type_of_mixture
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the string representation of the type of mixture
 *   of the support sets.  This string is the string representation
 *   of a Python list of integers. */

static PyObject *py2c_celcon_append_lifted_point
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Appends a lifted point to the cells container.
 *   There are three input parameters:
 *   1) the dimension of the point;
 *   2) the index of the support to where to append to; and
 *   3) the string representation of the lifted point.
 *   Returns the failure code, which equals zero when all went well. */

static PyObject *py2c_celcon_retrieve_lifted_point
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns a string representation of a lifted point.
 *   On input are three integer numbers:
 *   1) the number of coordinates in the lifted point;
 *   2) the index to the support set; and
 *   3) the index to the point in that support set. */

static PyObject *py2c_celcon_mixed_volume_of_supports
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the mixed volume of the supports stored
 *   in the cell container. */

static PyObject *py2c_celcon_number_of_cells
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the number of cells in the cell container. */

static PyObject *py2c_celcon_standard_random_coefficient_system 
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Based on the lifted supports stored in the container,
 *   a random coefficient system with coefficients in standard double
 *   precision is stored in the cell container. */

static PyObject *py2c_celcon_dobldobl_random_coefficient_system 
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Based on the lifted supports stored in the container,
 *   a random coefficient system with coefficients in double double
 *   precision is stored in the cell container. */

static PyObject *py2c_celcon_quaddobl_random_coefficient_system 
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Based on the lifted supports stored in the container,
 *   a random coefficient system with coefficients in quad double
 *   precision is stored in the cell container. */

static PyObject *py2c_celcon_copy_into_standard_systems_container
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   The random coefficient system in standard double precision is copied 
 *   from the cell container to the container for systems with 
 *   coefficients in standard double precision. */

static PyObject *py2c_celcon_copy_into_dobldobl_systems_container
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   The random coefficient system in double double precision is copied 
 *   from the cell container to the container for systems with 
 *   coefficients in double double precision. */

static PyObject *py2c_celcon_copy_into_quaddobl_systems_container
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   The random coefficient system in quad double precision is copied 
 *   from the cell container to the container for systems with 
 *   coefficients in quad double precision. */

static PyObject *py2c_celcon_standard_polyhedral_homotopy
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Based on the lifting and the random coefficient system,
 *   the polyhedral homotopy to solve the random coefficient system 
 *   in standard double precision is constructed.
 *   This function also initializes the internal data structures to store
 *   the solutions of start and target systems.
 *   The lifted supports and the random coefficient system are defined.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_celcon_dobldobl_polyhedral_homotopy
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Based on the lifting and the random coefficient system,
 *   the polyhedral homotopy to solve the random coefficient system 
 *   in double double precision is constructed.
 *   This function also initializes the internal data structures to store
 *   the solutions of start and target systems.
 *   The lifted supports and the random coefficient system are defined.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_celcon_quaddobl_polyhedral_homotopy
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Based on the lifting and the random coefficient system,
 *   the polyhedral homotopy to solve the random coefficient system 
 *   in quad double precision is constructed.
 *   This function also initializes the internal data structures to store
 *   the solutions of start and target systems.
 *   The lifted supports and the random coefficient system are defined.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_celcon_solve_standard_start_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Solves the start system corresponding to the k-th mixed cell,
 *   using standard double precision arithmetic.
 *   The precondition for this function is that the creation of
 *   the polyhedral homotopy in standard double precision ended well.
 *   On return is the number of solution found, which must equal
 *   the mixed volume of the k-th mixed cell. */

static PyObject *py2c_celcon_solve_dobldobl_start_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Solves the start system corresponding to the k-th mixed cell,
 *   using double double precision arithmetic.
 *   The precondition for this function is that the creation of
 *   the polyhedral homotopy in double double precision ended well.
 *   On return is the number of solution found, which must equal
 *   the mixed volume of the k-th mixed cell. */

static PyObject *py2c_celcon_solve_quaddobl_start_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Solves the start system corresponding to the k-th mixed cell,
 *   using quad double precision arithmetic.
 *   The precondition for this function is that the creation of
 *   the polyhedral homotopy in quad double precision ended well.
 *   On return is the number of solution found, which must equal
 *   the mixed volume of the k-th mixed cell. */

static PyObject *py2c_celcon_track_standard_solution_path
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Tracks a solution path starting at the i-th solution of the k-th cell,
 *   using standard double precision arithmetic.
 *   The precondition for this function is that the start system defined
 *   by the k-th mixed cell is solved in standard double precision.
 *   There are three input parameters:
 *   1) k, the index to a mixed cell in the cell container;
 *   2) i, the index to a solution path defined by that mixed cell;
 *   3) otp, the level for intermediate output during path tracking.
 *   A target solution corresponding to the k-th cell is added on return. */

static PyObject *py2c_celcon_track_dobldobl_solution_path
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Tracks a solution path starting at the i-th solution of the k-th cell,
 *   using double double precision arithmetic.
 *   The precondition for this function is that the start system defined
 *   by the k-th mixed cell is solved in double double precision.
 *   There are three input parameters:
 *   1) k, the index to a mixed cell in the cell container;
 *   2) i, the index to a solution path defined by that mixed cell;
 *   3) otp, the level for intermediate output during path tracking.
 *   A target solution corresponding to the k-th cell is added on return. */

static PyObject *py2c_celcon_track_quaddobl_solution_path
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Tracks a solution path starting at the i-th solution of the k-th cell,
 *   using quad double precision arithmetic.
 *   The precondition for this function is that the start system defined
 *   by the k-th mixed cell is solved in quad double precision.
 *   There are three input parameters:
 *   1) k, the index to a mixed cell in the cell container;
 *   2) i, the index to a solution path defined by that mixed cell;
 *   3) otp, the level for intermediate output during path tracking.
 *   A target solution corresponding to the k-th cell is added on return. */

static PyObject *py2c_celcon_copy_target_standard_solution_to_container
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Copies the i-th target solution corresponding to the k-th mixed cell
 *   to the container for solutions in standard double precision.
 *   There are two input parameters for this function:
 *   1) k, the index to the mixed cell;
 *   2) i, the index to the i-th solution path defined by the cell.
 *   On return is the failure code, which equals zero when all went well. */

static PyObject *py2c_celcon_copy_target_dobldobl_solution_to_container
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Copies the i-th target solution corresponding to the k-th mixed cell
 *   to the container for solutions in double double precision.
 *   There are two input parameters for this function:
 *   1) k, the index to the mixed cell;
 *   2) i, the index to the i-th solution path defined by the cell.
 *   On return is the failure code, which equals zero when all went well. */

static PyObject *py2c_celcon_copy_target_quaddobl_solution_to_container
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Copies the i-th target solution corresponding to the k-th mixed cell
 *   to the container for solutions in quad double precision.
 *   There are two input parameters for this function:
 *   1) k, the index to the mixed cell;
 *   2) i, the index to the i-th solution path defined by the cell.
 *   On return is the failure code, which equals zero when all went well. */

static PyObject *py2c_celcon_permute_standard_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Permutes the systems in the container for polynomial and Laurent systems
 *   with standard double coefficients corresponding to the permutation
 *   used to compute the mixed-cell configuration.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_celcon_permute_dobldobl_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Permutes the systems in the container for polynomial and Laurent systems
 *   with double double coefficients corresponding to the permutation
 *   used to compute the mixed-cell configuration.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_celcon_permute_quaddobl_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Permutes the systems in the container for polynomial and Laurent systems
 *   with quad double coefficients corresponding to the permutation
 *   used to compute the mixed-cell configuration.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_celcon_clear_container
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Deallocates the data in the cell container. */

/* wrapping functions to scale polynomial systems and solutions */

static PyObject *py2c_scale_standard_system ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Applies scaling to the system in the standard systems container,
 *   with standard double precision arithmetic.  The system in the standard
 *   systems container is replaced by the scaled system.
 *   On entry is one integer, which should be either 0, 1, or 2:
 *   0 for only scaling of the equations,
 *   1 variable scaling without variability reduction,
 *   2 variable scaling with variability reduction.
 *   On return is a tuple with the scaling coefficients (if mode > 0)
 *   and the estimated inverse condition number of the scaling problem. */

static PyObject *py2c_scale_dobldobl_system ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Applies scaling to the system in the dobldobl systems container,
 *   with double double precision arithmetic.  The system in the dobldobl
 *   systems container is replaced by the scaled system.
 *   On entry is one integer, which should be either 0, 1, or 2:
 *   0 for only scaling of the equations,
 *   1 variable scaling without variability reduction,
 *   2 variable scaling with variability reduction.
 *   On return is a tuple with the scaling coefficients (if mode > 0)
 *   and the estimated inverse condition number of the scaling problem. */

static PyObject *py2c_scale_quaddobl_system ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Applies scaling to the system in the quaddobl systems container,
 *   with quad double precision arithmetic.  The system in the quaddobl
 *   systems container is replaced by the scaled system.
 *   On entry is one integer, which should be either 0, 1, or 2:
 *   0 for only scaling of the equations,
 *   1 variable scaling without variability reduction,
 *   2 variable scaling with variability reduction.
 *   On return is a tuple with the scaling coefficients (if mode > 0)
 *   and the estimated inverse condition number of the scaling problem. */

static PyObject *py2c_scale_standard_solutions 
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Replaces the solutions in the standard solutions container with
 *   the scaled solutions, scaled with standard double precision arithmetic,
 *   using the given scaling coefficients.
 *   On entry are two parameters: an integer and a string.
 *   The integer contains the number of elements in the list
 *   of scaling coefficients (doubles) stored in the string.
 *   The format of the string is the Python string representation
 *   of a list of doubles, i.e.: starting with '[' and ending with ']'. */

static PyObject *py2c_scale_dobldobl_solutions 
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Replaces the solutions in the dobldobl solutions container with
 *   the scaled solutions, scaled with double double precision arithmetic,
 *   using the given scaling coefficients.
 *   On entry are two parameters: an integer and a string.
 *   The integer contains the number of elements in the list
 *   of scaling coefficients (doubles) stored in the string.
 *   The format of the string is the Python string representation
 *   of a list of doubles, i.e.: starting with '[' and ending with ']'. */

static PyObject *py2c_scale_quaddobl_solutions 
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Replaces the solutions in the quaddobl solutions container with
 *   the scaled solutions, scaled with quad double precision arithmetic,
 *   using the given scaling coefficients.
 *   On entry are two parameters: an integer and a string.
 *   The integer contains the number of elements in the list
 *   of scaling coefficients (doubles) stored in a the string.
 *   The format of the string is the Python string representation
 *   of a list of doubles, i.e.: starting with '[' and ending with ']'. */

/* wrapping functions to manipulate algebraic sets */

static PyObject *py2c_embed_system ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Replaces the system with coefficients in standard double precision
 *   in the container with its embedding of dimension d.
 *   The dimension d is given as an integer parameter on input.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_standard_cascade_homotopy
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Creates a homotopy in standard double precision using the stored
 *   systems to go one level down the cascade, removing one slice.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_dobldobl_cascade_homotopy
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Creates a homotopy in double double precision using the stored
 *   systems to go one level down the cascade, removing one slice.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_quaddobl_cascade_homotopy
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Creates a homotopy in quad double precision using the stored
 *   systems to go one level down the cascade, removing one slice.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_factor_set_to_mute ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Sets the state of monodromy permutations to silent. */

static PyObject *py2c_factor_define_output_file_with_string
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Defines the output file for the factorization.
 *   On input are an integer and a string:
 *   1) the integer equals the number of characters in the string; and
 *   2) the string contains the name of a file.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_factor_assign_labels ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Assigns labels, replacing the multiplicity field of each solution
 *   in standard double precision stored in the container.
 *   On entry are two integers:
 *   1) n, the number of coordinates of the solutions;
 *   2) nbsols, the number of solutions in the container.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_factor_initialize_sampler
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Initializes the sampling machine with a witness set.
 *   On entry is the dimension or the number of hyperplanes
 *   to slide the positive dimensional solution set. */

static PyObject *py2c_factor_initialize_monodromy
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Initializes the internal data structures for n loops,
 *   to factor a k-dimensional solution component of degree d.
 *   There are three integers on input, in the following order:
 *   1) n, the number of loops;
 *   2) d, the degree of the solution set;
 *   3) k, the dimensional of the solution set.
 *   On return is the failure code, which equals zero when all went well. */

static PyObject *py2c_factor_store_solutions
 ( PyObject *self, PyObject *args );
/* 
 * DESCRIPTION :
 *   Stores the solutions in the container to the data for monodromy loops. */

static PyObject *py2c_factor_restore_solutions
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Restores the first initialized solutions from sampler to the container. */

static PyObject *py2c_factor_track_paths ( PyObject *self, PyObject *args );
/* 
 * DESCRIPTION :
 *   Tracks as many paths as defined by witness set.
 *   On return is the failure code, which is zero when all went well. */

static PyObject *py2c_factor_swap_slices ( PyObject *self, PyObject *args );
/* 
 * DESCRIPTION :
 *   Swaps the current slices with new slices and takes new solutions
 *   as start to turn back.
 *   On return is the failure code, which is zero when all went well. */

static PyObject *py2c_factor_new_slices ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Generates k random slides in n-space.
 *   The k and the n are the two input parameters.
 *   On return is the failure code, which is zero when all went well. */

static PyObject *py2c_factor_set_trace_slice
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Assigns the constant coefficient of the first slice.
 *   On entry is a flag to indicate if it was the first time or not.
 *   On return is the failure code, which is zero if all went well. */

static PyObject *py2c_factor_store_gammas ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Stores the gamma constants for the sampler in the monodromy loops.
 *   Generates as many random complex constants as the value on input.
 *   On return is the failure code, which is zero if all went well. */

static PyObject *py2c_factor_permutation_after_loop
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   For a solution set of degree d, computes the permutation using the
 *   solutions most recently stored, after a loop. 
 *   The number d is the input parameter of this function.
 *   On return is the string representation of the permutation. */

static PyObject *py2c_factor_update_decomposition
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Updates the decomposition with the given permutation of d elements.
 *   On entry are two integers and one string:
 *   1) d, the number of elements in the permutation;
 *   2) nc, the number of characters in the string;
 *   3) p, the string representation of the permutation.
 *   Returns one if the current decomposition is certified,
 *   otherwise returns zero. */

static PyObject *py2c_factor_number_of_components
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the number of irreducible factors in the current
 *   decomposition of the witness set. */

static PyObject *py2c_factor_witness_points_of_component
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns a string which represents an irreducible component.
 *   On entry are two integers:
 *   1) the sum of the degrees of all components;
 *   2) the index of the component. */

static PyObject *py2c_factor_trace_sum_difference
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the difference between the actual sum at the samples
 *   defined by the labels to the generic points in the factor,
 *   and the trace sum.
 *   On entry are three integer numbers and one string:
 *   1) d, the number of points in the witness set;
 *   2) k, the dimension of the solution set;
 *   3) nc, the number of characters in the string;
 *   4) ws, the string representing the labels of the witness set. */

static PyObject *py2c_witness_set_of_hypersurface
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Given in the string p of nc characters a polynomial in nv variables,
 *   terminated by a semicolon, the systems and solutions container on
 *   return contain a witness set for the hypersurface defined by p.
 *   On entry are two integers and one string, in the following order:
 *   1) nv, the number of variables of the polynomials;
 *   2) nc, the number of characters in the string p;
 *   3) p, string representation of a polynomial, terminates with ';'.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_create_diagonal_homotopy
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Creates a diagonal homotopy to intersect two solution sets of
 *   dimensions a and b respectively, where a >= b.
 *   The two input parameters are values for a and b.
 *   The systems stored as target and start system in the container
 *   define the witness sets for these two solution sets. */

static PyObject *py2c_start_diagonal_cascade_solutions
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Makes the start solutions to start the cascade homotopy to
 *   intersect two solution sets of dimensions a and b, where a >= b.
 *   The dimensions a and b are given as input parameters.
 *   The systems stored as target and start system in the container
 *   define the witness sets for these two solution sets.
 *   On return is the failure code, which equals zero when all went well. */

static PyObject *py2c_extrinsic_top_diagonal_dimension
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the dimension of the start and target system to
 *   start the extrinsic cascade to intersect two witness sets,
 *   respectively of dimensions a and b, with ambient dimensions
 *   respectively equal to n1 and n2.
 *   There are four integers as parameters on input: n1, n2, a and b. */

static PyObject *py2c_collapse_diagonal ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Eliminates the extrinsic diagonal for the system and solutions
 *   in the containers.  On input are two integers:
 *   1) k, the current number of slack variables in the embedding;
 *   2) d, the number of slack variables to add to the final embedding.
 *   The system in the container has its diagonal eliminated and is
 *   embedded with k+d slack variables.  The solutions corresponding
 *   to this system are in the solutions container.
 *   On return is the failure code, which equals zero if all went well. */

/* The wrapping of Pieri and Littlewood-Richardson homotopies,
 * with prototypes in schubert.h starts here. */

static PyObject *py2c_schubert_pieri_count
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the number of p-plane producing curves of degree q
 *   that meet m*p + q*(m+p) given general m-planes.
 *   On input are three integer numbers:
 *   1) m, the dimension of the input planes;
 *   2) p, the dimension of the output planes; and
 *   3) q, the degree of the curve that produces p-planes.
 *   The dimension of the ambient space of this Pieri problem is m+p. */

static PyObject *py2c_schubert_resolve_conditions
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Resolves a general Schubert intersection condition in n-space
 *   for k-planes subject to conditions defined by brackers.
 *   On return is the root count, the number of k-planes that satisfy
 *   the intersection conditions imposed by the brackets for general flags.
 *   On entry are five integers and one string:
 *   1) n, the ambient dimension, where the k-planes live;
 *   2) k, the dimension of the solution planes;
 *   3) c, the number of intersection conditions;
 *   4) nc, the number of characters in the string brackets;
 *   5) brackets is a string representation of c brackets, where the numbers
 *   in each bracket are separated by spaces;
 *   6) the flag verbose: when 0, no intermediate output is written,
 *   when 1, then the resolution is dispayed on screen. */

static PyObject *py2c_schubert_littlewood_richardson_homotopies
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Runs the Littlewood-Richardson homotopies to resolve a number of
 *   general Schubert intersection conditions on k-planes in n-space.
 *   The polynomial system that was solved is in the container for
 *   systems with coefficients in standard double precision and the
 *   corresponding solutions are in the standard solutions container.
 *   On entry are six integers and two strings, in the following order:
 *   1) n, the ambient dimension, where the k-planes live;
 *   2) k, the dimension of the solution planes;
 *   3) c,the number of intersection conditions;
 *   4) nc, the number of characters in the string brackets;
 *   5) brackets is a string representation of c brackets, where the numbers
 *   in each bracket are separated by spaces;
 *   6) the flag verbose: when 0, no intermediate output is written,
 *   when 1, then the resolution is dispayed on screen;
 *   7) nbchar, the number of characters in the string filename;
 *   8) filename is the name of the output file.
 *   The function returns a tuple of an integer and a string:
 *   0) r is the formal root count as the number of k-planes
 *   for conditions imposed by the brackets for general flags;
 *   1) flags, a string with the coefficients of the general flags. */

static PyObject *py2c_schubert_localization_poset
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the string representation of the localization poset for the
 *   Pieri root count for m, p, and q.  The input parameters are the
 *   integer values for m, p, and q:
 *   1) m, the dimension of the input planes;
 *   2) p, the dimension of the output planes;
 *   3) q, the degree of the curves that produce p-planes. */

static PyObject *py2c_schubert_pieri_homotopies
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Runs the Pieri homotopies for (m,p,q) dimensions on generic input data.
 *   On return the systems container for systems with coefficients in standard
 *   double precision contains the polynomial system solved and in the
 *   solutions in standard double precision are in the solutions container.
 *   On entry are four integers and two strings:
 *   1) m, the dimension of the input planes;
 *   2) p, the dimension of the output planes;
 *   3) q, the degree of the solution maps;
 *   4) nc, the number of characters in the string A;
 *   5) A, the string with m*p + q*(m+p) random complex input m-planes,
 *   where the real and imaginary parts are separated by a space;
 *   6) pts, the string with m*p + q*(m+p) random complex interpolation
 *   points, only needed if q > 0.
 *   The function returns the combinatorial Pieri root count,
 *   which should equal the number of solutions in the container. */

static PyObject *py2c_schubert_osculating_planes
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the string representation of n real m-planes in
 *   d-space osculating a rational normal curve
 *   at the n points in s, where n = m*p + q*(m+p) and d = m+p.
 *   On entry are four integers and one string:
 *   1) m, the dimension of the input planes;
 *   2) p, the dimension of the output planes;
 *   3) q, the degree of the solution maps;
 *   4) nc, the number of characters in the string pts; and
 *   5) pts, the string with m*p + q*(m+p) interpolation points. */

static PyObject *py2c_schubert_pieri_system
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Fills the container of systems with coefficients in standard
 *   double precision with a polynomial system that expresses the
 *   intersection conditions of a general Pieri problem.
 *   On input are five integers and one string:
 *   1) m, the dimension of the input planes;
 *   2) p, the dimension of the output planes;
 *   3) q, the degree of the solution maps;
 *   4) nc, the number of characters in the string A;
 *   5) A,  m*p + q*(m+p) random complex input m-planes, where
 *   the real and imaginary parts are separated by a space;
 *   6) a flag is_real: if == 1, then the coefficients of A are real,
 *   if == 0, then the coefficients of A are complex.
 *   Returns the failure code, which equals zero if all went well. */

/* The wrapping functions in mapcon.h starts here. */

static PyObject *py2c_mapcon_solve_system ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *  Solves the binomial system stored in the Laurent systems container.
 *  There is one input argument, either one or zero.
 *  If one, then only the pure top dimensional solutions are computed.
 *  If zero, then all solution sets are computed.
 *  Returns the failure code, which equals zero if all went well. */

static PyObject *py2c_mapcon_write_maps ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Writes the maps stored in the container to screen.
 *   Returns the failure code, which equals zero if all went well. */

static PyObject *py2c_mapcon_clear_maps ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Deallocates the maps stored in the container.
 *   Returns the failure code, which equals zero if all went well. */

static PyObject *py2c_mapcon_top_dimension ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the top dimension of the maps in the container. */

static PyObject *py2c_mapcon_number_of_maps ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the number of maps in the container. */

static PyObject *py2c_mapcon_degree_of_map ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Given the dimension and index of a map, given as two integers as
 *   input parameters, returns the degree of that map. */

static PyObject *py2c_mapcon_coefficients_of_map
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the coefficients of a monomial map stored in the container.
 *   On entry are three parameters:
 *   1) the dimension of the map;
 *   2) the index of the map in all maps of that dimension;
 *   3) the number of variables.
 *   On return is a Python list of complex doubles. */

static PyObject *py2c_mapcon_exponents_of_map
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Returns the exponents of a monomial map stored in the container.
 *   On entry are three parameters:
 *   1) the dimension of the map;
 *   2) the index of the map in all maps of that dimension;
 *   3) the number of variables.
 *   On return is a Python list of integers. */

/* The wrapping of functions with prototypes in next_track.h starts below. */

static PyObject *py2c_initialize_standard_homotopy
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Initializes the homotopy to track a path with a generator,
 *   using standard double precision arithmetic.
 *   There is one integer number on input to be considered as a boolean,
 *   as an indicator whether a fixed gamma constant will be used.
 *   Before calling this routine the target and start system must
 *   be copied over from the standard systems container. */

static PyObject *py2c_initialize_dobldobl_homotopy
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Initializes the homotopy to track a path with a generator,
 *   using double double precision arithmetic.
 *   There is one integer number on input to be considered as a boolean,
 *   as an indicator whether a fixed gamma constant will be used.
 *   Before calling this routine the target and start system must
 *   be copied over from the dobldobl systems container. */

static PyObject *py2c_initialize_quaddobl_homotopy
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Initializes the homotopy to track a path with a generator,
 *   using quad double precision arithmetic.
 *   There is one integer number on input to be considered as a boolean,
 *   as an indicator whether a fixed gamma constant will be used.
 *   Before calling this routine the target and start system must
 *   be copied over from the quaddobl systems container. */

static PyObject *py2c_initialize_multprec_homotopy
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Initializes the homotopy to track a path with a generator,
 *   using arbitrary multiprecision arithmetic.
 *   There is are two integer numbers on input:
 *   1) one to be considered as a boolean,
 *   as an indicator whether a fixed gamma constant will be used; and
 *   2) the number of decimal places in the working precision.
 *   Before calling this routine the target and start system must
 *   be copied over from the multprec systems container. */

static PyObject *py2c_initialize_varprec_homotopy
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Initializes the variable precision homotopy with the target and
 *   start system stored in the strings.
 *   On entry are three integers and two strings, in the following order:
 *   1) fixed_gamma is a flag: if 1, then a fixed value for the gamma constant
 *   is used, if 0, a random value for gamma will be generated;
 *   2) nc_target, the number of characters in the string target;
 *   3) target, the string representation of the target system;
 *   4) nc_start, the number of characters in the string start;
 *   5) start, the string representation of the start system.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_initialize_standard_solution
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Initializes the path tracker with a generator with a solution
 *   from the standard solutions container.  The index to the solution
 *   is given as an integer input parameter.  The counting of the
 *   indices starts at one, so the first solution has index one. */

static PyObject *py2c_initialize_dobldobl_solution
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Initializes the path tracker with a generator with a solution
 *   from the dobldobl solutions container.  The index to the solution
 *   is given as an integer input parameter.  The counting of the
 *   indices starts at one, so the first solution has index one. */

static PyObject *py2c_initialize_quaddobl_solution
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Initializes the path tracker with a generator with a solution
 *   from the quaddobl solutions container.  The index to the solution
 *   is given as an integer input parameter.  The counting of the
 *   indices starts at one, so the first solution has index one. */

static PyObject *py2c_initialize_multprec_solution
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Initializes the path tracker with a generator with a solution
 *   from the multprec solutions container.  The index to the solution
 *   is given as an integer input parameter.  The counting of the
 *   indices starts at one, so the first solution has index one. */

static PyObject *py2c_initialize_varbprec_solution
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Uses the string representation of a solution to initialize the
 *   variable precision path tracker with.
 *   There are three input parameters, two integers and one string:
 *   1) nv, the number of variables in the solution;
 *   2) nc, the number of characters in the string sol;
 *   3) sol, the string representation of a solution.
 *   On return is the failure code, which equals zero if all went well. */

static PyObject *py2c_next_standard_solution
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Computes the next point on the solution path with standard double
 *   precision for the given index.  This index is given as an input
 *   parameter.  The index to the solution path starts its count at one.
 *   The point itself is stored in the standard solutions container.
 *   The functions py2c_initialized_standard_tracker and
 *   py2c_initialize_standard_solution must have been executed earlier.
 *   The failcode is returned, which equals zero if all is well. */

static PyObject *py2c_next_dobldobl_solution
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Computes the next point on the solution path with double double
 *   precision for the given index.  This index is given as an input
 *   parameter.  The index to the solution path starts its count at one.
 *   The point itself is stored in the dobldobl solutions container.
 *   The functions py2c_initialized_dobldobl_tracker and
 *   py2c_initialize_dobldobl_solution must have been executed earlier.
 *   The failcode is returned, which equals zero if all is well. */

static PyObject *py2c_next_quaddobl_solution
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Computes the next point on the solution path with quad double
 *   precision for the given index.  This index is given as an input
 *   parameter.  The index to the solution path starts its count at one.
 *   The point itself is stored in the quaddobl solutions container.
 *   The functions py2c_initialized_quaddobl_tracker and
 *   py2c_initialize_quaddobl_solution must have been executed earlier.
 *   The failcode is returned, which equals zero if all is well. */

static PyObject *py2c_next_multprec_solution
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Computes the next point on the solution path with arbitrary
 *   multiprecision for the given index.  This index is given as an input
 *   parameter.  The index to the solution path starts its count at one.
 *   The point itself is stored in the multprec solutions container.
 *   The functions py2c_initialized_multprec_tracker and
 *   py2c_initialize_multprec_solution must have been executed earlier.
 *   The failcode is returned, which equals zero if all is well. */

static PyObject *py2c_next_varbprec_solution
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Computes the next point on a solution path in variable precision.
 *   There are four integer input parameters:
 *   1) the number of correct decimal places in the solution;
 *   2) an upper bound on the number of decimal places in the precision;
 *   3) the maximum number of Newton iterations;
 *   4) a flag zero or one to indicate the verbose level.
 *   On return is a tuple:
 *   0) the failure code, which equals zero if all went well; and
 *   1) the string representation of the next solution on the path. */

static PyObject *py2c_clear_standard_tracker
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Deallocates data used in the standard double precision tracker
 *   with a generator. */

static PyObject *py2c_clear_dobldobl_tracker
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Deallocates data used in the double double precision tracker
 *   with a generator. */

static PyObject *py2c_clear_quaddobl_tracker
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Deallocates data used in the quad double precision tracker
 *   with a generator. */

static PyObject *py2c_clear_multprec_tracker
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Deallocates data used in the arbitrary multiprecision tracker
 *   with a generator. */

static PyObject *py2c_clear_varbprec_tracker
 ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Deallocates data used in the variable precision tracker
 *   with a generator. */

/* The wrapping of Newton's method and path trackers with the evaluation
 * done by algorithmic differentiation in quad double precision by the 
 * functions in adepath_qd.h, starts here. */

static PyObject *py2c_ade_newton_qd ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Runs Newton's method with algorithmic differentation in quad
 *   double precision on the data in the systems and solutions container.
 *
 * REQUIRED :
 *   The quaddobl systems container contains a valid polynomial system
 *   and the quaddobl solutions container holds a valid solution.
 *
 * ON ENTRY :
 *   verbose  0 if no intermediate output is wanted,
 *            1 if extra information should be written to screen.
 *
 * ON RETURN :
 *   fail     0 if all went well, and the corrected solution is in the
 *              solution container,
 *            if different from zero, then an error happened. */

static PyObject *py2c_ade_onepath_qd ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Tracks one solution path with algorithmic differentation in quad
 *   double precision on the data in the systems and solutions container.
 *
 * REQUIRED :
 *   The start and target systems have been defined
 *   and the quaddobl solutions container holds a valid solution.
 *
 * ON ENTRY :
 *   verbose  0 if no intermediate output is wanted,
 *            1 if extra information should be written to screen;
 *   regamma  real part of the random gamma constant;
 *   imgamma  imaginary part of the random constant.
 *
 * ON RETURN :
 *   fail     0 if all went well, and the solution at the end of the path 
 *              is in the  solution container,
 *            if different from 0, then an error happened. */

static PyObject *py2c_ade_manypaths_qd ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Tracks many solution paths with algorithmic differentation in quad
 *   double precision on the data in the systems and solutions container.
 *
 * REQUIRED :
 *   The start and target systems have been defined
 *   and the quaddobl solutions container holds valid solutions.
 *
 * ON ENTRY :
 *   verbose  0 if no intermediate output is wanted,
 *            1 if extra information should be written to screen;
 *   regamma  real part of the random gamma constant;
 *   imgamma  imaginary part of the random constant.
 *
 * ON RETURN :
 *   fail     0 if all went well, and the solutions at the end of paths
 *              are in the solution container,
 *            if different from 0, then an error happened. */

/* The wrapping of Newton's method and path trackers with the evaluation
 * done by algorithmic differentiation in quad double precision on the GPU 
 * by the functions in gpupath_qd.h, starts here. */

static PyObject *py2c_gpu_newton_qd ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Runs Newton's method with algorithmic differentation on CPU and GPU
 *   in quad double precision on the data in the systems and solutions
 *   container.
 *
 * REQUIRED :
 *   The quaddobl systems container contains a valid polynomial system
 *   and the quaddobl solutions container holds a valid solution.
 *
 * ON ENTRY :
 *   mode     execution mode is 0 (CPU+GPU), 1 (CPU), or 2 (GPU);
 *   verbose  0 if no intermediate output is wanted,
 *            1 if extra information should be written to screen.
 *
 * ON RETURN :
 *   fail     0 if all went well, and the corrected solution is in the
 *              solution container,
 *            if different from zero, then an error happened. */

static PyObject *py2c_gpu_onepath_qd ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Tracks one solution path with algorithmic differentation on CPU and GPU
 *   in quad double precision on the data in the systems and solutions
 *   container.
 *
 * REQUIRED :
 *   The start and target systems have been defined
 *   and the quaddobl solutions container holds a valid solution.
 *
 * ON ENTRY :
 *   mode     execution mode is 0 (CPU+GPU), 1 (CPU), or 2 (GPU);
 *   verbose  0 if no intermediate output is wanted,
 *            1 if extra information should be written to screen;
 *   regamma  real part of the random gamma constant;
 *   imgamma  imaginary part of the random constant.
 *
 * ON RETURN :
 *   fail     0 if all went well, and the solution at the end of the path 
 *              is in the  solution container,
 *            if different from 0, then an error happened. */

static PyObject *py2c_gpu_manypaths_qd ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Tracks many solution paths with algorithmic differentation on CPU and GPU
 *   in quad double precision on the data in the systems and solutions
 *   container.
 *
 * REQUIRED :
 *   The start and target systems have been defined
 *   and the quaddobl solutions container holds valid solutions.
 *
 * ON ENTRY :
 *   mode     execution mode is 0 (CPU+GPU), 1 (CPU), or 2 (GPU);
 *   verbose  0 if no intermediate output is wanted,
 *            1 if extra information should be written to screen;
 *   regamma  real part of the random gamma constant;
 *   imgamma  imaginary part of the random constant.
 *
 * ON RETURN :
 *   fail     0 if all went well, and the solutions at the end of paths
 *              are in the solution container,
 *            if different from 0, then an error happened. */
