/* This file contains the prototypes for the py2c interface functions. */

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

static PyObject *py2c_write_target_system ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Writes the target system as stored in standard double precision
 *   to screen or to the defined output file. */

static PyObject *py2c_write_start_system ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Writes the target system as stored in standard double precision
 *   to screen or to the defined output file. */

/* Copying systems from and to containers. */

static PyObject *py2c_copy_target_system_to_container
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

static PyObject *py2c_copy_container_to_target_system
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

static PyObject *py2c_copy_container_to_start_system 
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

/* creation of homotopy and tracking all paths */

static PyObject *py2c_create_homotopy ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Initializes the data for a homotopy in standard double precision.
 *   The failure code is returned, which is zero when all goes well. */

static PyObject *py2c_create_homotopy_with_gamma
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

static PyObject *py2c_clear_homotopy ( PyObject *self, PyObject *args );
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

/* moving solutions from and to containers */

static PyObject *py2c_copy_target_solutions_to_container
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

static PyObject *py2c_copy_container_to_target_solutions
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

static PyObject *py2c_copy_container_to_start_solutions
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

static PyObject *py2c_syscon_read_system ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Interactive procedure to read a polynomial system with coefficients
 *   in standard double precision.
 *   The system will be placed in the standard systems container.
 *   The failure code is returned, which equals zero if all went well. */

static PyObject *py2c_syscon_read_Laurent_system
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

static PyObject *py2c_syscon_write_system ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Writes the polynomial system with standard double precision coefficients
 *   that is stored in the container. */

static PyObject *py2c_syscon_write_Laurent_system
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

static PyObject *py2c_syscon_clear_system ( PyObject *self, PyObject *args );
/*
 * DESCRIPTION :
 *   Deallocates the container for polynomial systems
 *   with coefficients in standard double precision. */

static PyObject *py2c_syscon_clear_Laurent_system
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
static PyObject *py2c_syscon_write_symbols
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_string_of_symbols
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_remove_symbol_name
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_clear_symbol_table
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_number_of_polynomials
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_number_of_dobldobl_polynomials
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_number_of_quaddobl_polynomials
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_number_of_multprec_polynomials
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_number_of_Laurentials
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_number_of_dobldobl_Laurentials
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_number_of_quaddobl_Laurentials
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_number_of_multprec_Laurentials
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_initialize_number
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_initialize_number_of_dobldobl_polynomials
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_initialize_number_of_quaddobl_polynomials
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_initialize_number_of_multprec_polynomials
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_initialize_number_of_Laurentials
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_initialize_number_of_dobldobl_Laurentials
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_initialize_number_of_quaddobl_Laurentials
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_initialize_number_of_multprec_Laurentials
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_degree_of_polynomial
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_degree_of_dobldobl_polynomial
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_degree_of_quaddobl_polynomial
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_degree_of_multprec_polynomial
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_number_of_terms ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_number_of_Laurent_terms
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_retrieve_term ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_store_polynomial
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_store_dobldobl_polynomial
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_store_quaddobl_polynomial
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_store_multprec_polynomial
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_load_polynomial
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_load_dobldobl_polynomial
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_load_quaddobl_polynomial
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_load_multprec_polynomial
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_store_Laurential
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_store_dobldobl_Laurential
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_store_quaddobl_Laurential
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_store_multprec_Laurential
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_load_standard_Laurential
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_load_dobldobl_Laurential
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_load_quaddobl_Laurential
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_load_multprec_Laurential
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_total_degree ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_standard_drop_variable_by_index
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syson_standard_drop_variable_by_name
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_dobldobl_drop_variable_by_index
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_dobldobl_drop_variable_by_name
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_quaddobl_drop_variable_by_index
 ( PyObject *self, PyObject *args );
static PyObject *py2c_syscon_quaddobl_drop_variable_by_name
 ( PyObject *self, PyObject *args );

/* wrapping functions in solcon.h starts from here */

static PyObject *py2c_solcon_read_solutions ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_read_dobldobl_solutions
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_read_quaddobl_solutions
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_read_multprec_solutions
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_write_solutions ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_write_dobldobl_solutions
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_write_quaddobl_solutions
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_write_multprec_solutions
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_clear_solutions ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_clear_dobldobl_solutions
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_clear_quaddobl_solutions
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_clear_multprec_solutions
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_open_solution_input_file
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_length_solution_string
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_length_dobldobl_solution_string
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_length_quaddobl_solution_string
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_length_multprec_solution_string
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_write_solution_string
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_write_dobldobl_solution_string
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_write_quaddobl_solution_string
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_write_multprec_solution_string
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_retrieve_next_standard_initialize
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_retrieve_next_dobldobl_initialize
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_retrieve_next_quaddobl_initialize
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_retrieve_next_multprec_initialize
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_move_current_standard_to_next
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_move_current_dobldobl_to_next
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_move_current_quaddobl_to_next
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_move_current_multprec_to_next
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_length_current_standard_solution_string
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_length_current_dobldobl_solution_string
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_length_current_quaddobl_solution_string
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_length_current_multprec_solution_string
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_write_current_standard_solution_string
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_write_current_dobldobl_solution_string
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_write_current_quaddobl_solution_string
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_write_current_multprec_solution_string
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_append_solution_string
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_append_dobldobl_solution_string
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_append_quaddobl_solution_string
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_append_multprec_solution_string
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_number_of_solutions
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_number_of_dobldobl_solutions
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_number_of_quaddobl_solutions
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_number_of_multprec_solutions
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_standard_drop_coordinate_by_index
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_standard_drop_coordinate_by_name
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_dobldobl_drop_coordinate_by_index
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_dobldobl_drop_coordinate_by_name
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_quaddobl_drop_coordinate_by_index
 ( PyObject *self, PyObject *args );
static PyObject *py2c_solcon_quaddobl_drop_coordinate_by_name
 ( PyObject *self, PyObject *args );

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
 *   Writes the supporting set structure to screen .*/

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
 *   stored set structure.   On returns is the failure code,
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
 *   Give a partition of the set of variables, computes 
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
 *   Give a partition of the set of variables, constructs
 *   an m-homogeneous Bezout number for the system in
 *   the standard systems container.
 *   On input are two arguments:
 *   1) the number of characters in the string (second argument); and
 *   2) the string representation for a partition of the variables.
 *   On return is the m-homogeneous Bezout number. */

/* wrapping functions in celcon.h starts here */

static PyObject *py2c_celcon_initialize_supports
 ( PyObject *self, PyObject *args );
static PyObject *py2c_celcon_set_type_of_mixture
 ( PyObject *self, PyObject *args );
static PyObject *py2c_celcon_type_of_mixture
 ( PyObject *self, PyObject *args );
static PyObject *py2c_celcon_append_lifted_point
 ( PyObject *self, PyObject *args );
static PyObject *py2c_celcon_retrieve_lifted_point
 ( PyObject *self, PyObject *args );
static PyObject *py2c_celcon_mixed_volume_of_supports
 ( PyObject *self, PyObject *args );
static PyObject *py2c_celcon_number_of_cells
 ( PyObject *self, PyObject *args );
static PyObject *py2c_celcon_create_random_coefficient_system 
 ( PyObject *self, PyObject *args );
static PyObject *py2c_celcon_dobldobl_random_coefficient_system 
 ( PyObject *self, PyObject *args );
static PyObject *py2c_celcon_quaddobl_random_coefficient_system 
 ( PyObject *self, PyObject *args );
static PyObject *py2c_celcon_copy_into_systems_container
 ( PyObject *self, PyObject *args );
static PyObject *py2c_celcon_copy_into_dobldobl_systems_container
 ( PyObject *self, PyObject *args );
static PyObject *py2c_celcon_copy_into_quaddobl_systems_container
 ( PyObject *self, PyObject *args );
static PyObject *py2c_celcon_create_polyhedral_homotopy
 ( PyObject *self, PyObject *args );
static PyObject *py2c_celcon_dobldobl_polyhedral_homotopy
 ( PyObject *self, PyObject *args );
static PyObject *py2c_celcon_quaddobl_polyhedral_homotopy
 ( PyObject *self, PyObject *args );
static PyObject *py2c_celcon_solve_start_system
 ( PyObject *self, PyObject *args );
static PyObject *py2c_celcon_solve_dobldobl_start_system
 ( PyObject *self, PyObject *args );
static PyObject *py2c_celcon_solve_quaddobl_start_system
 ( PyObject *self, PyObject *args );
static PyObject *py2c_celcon_track_solution_path
 ( PyObject *self, PyObject *args );
static PyObject *py2c_celcon_track_dobldobl_solution_path
 ( PyObject *self, PyObject *args );
static PyObject *py2c_celcon_track_quaddobl_solution_path
 ( PyObject *self, PyObject *args );
static PyObject *py2c_celcon_copy_target_solution_to_container
 ( PyObject *self, PyObject *args );
static PyObject *py2c_celcon_copy_target_dobldobl_solution_to_container
 ( PyObject *self, PyObject *args );
static PyObject *py2c_celcon_copy_target_quaddobl_solution_to_container
 ( PyObject *self, PyObject *args );
static PyObject *py2c_celcon_permute_system
 ( PyObject *self, PyObject *args );
static PyObject *py2c_celcon_permute_dobldobl_system
 ( PyObject *self, PyObject *args );
static PyObject *py2c_celcon_permute_quaddobl_system
 ( PyObject *self, PyObject *args );
static PyObject *py2c_celcon_clear_container
 ( PyObject *self, PyObject *args );

/* wrapping functions to manipulate algebraic sets */

static PyObject *py2c_embed_system ( PyObject *self, PyObject *args );
static PyObject *py2c_standard_cascade_homotopy
 ( PyObject *self, PyObject *args );
static PyObject *py2c_dobldobl_cascade_homotopy
 ( PyObject *self, PyObject *args );
static PyObject *py2c_quaddobl_cascade_homotopy
 ( PyObject *self, PyObject *args );
static PyObject *py2c_factor_set_to_mute ( PyObject *self, PyObject *args );
static PyObject *py2c_factor_define_output_file
 ( PyObject *self, PyObject *args );
static PyObject *py2c_factor_assign_labels ( PyObject *self, PyObject *args );
static PyObject *py2c_factor_initialize_sampler
 ( PyObject *self, PyObject *args );
static PyObject *py2c_factor_initialize_monodromy
 ( PyObject *self, PyObject *args );
static PyObject *py2c_factor_store_solutions
 ( PyObject *self, PyObject *args );
static PyObject *py2c_factor_restore_solutions
 ( PyObject *self, PyObject *args );
static PyObject *py2c_factor_track_paths ( PyObject *self, PyObject *args );
static PyObject *py2c_factor_swap_slices ( PyObject *self, PyObject *args );
static PyObject *py2c_factor_new_slices ( PyObject *self, PyObject *args );
static PyObject *py2c_factor_set_trace_slice
 ( PyObject *self, PyObject *args );
static PyObject *py2c_factor_store_gammas ( PyObject *self, PyObject *args );
static PyObject *py2c_factor_permutation_after_loop
 ( PyObject *self, PyObject *args );
static PyObject *py2c_factor_update_decomposition
 ( PyObject *self, PyObject *args );
static PyObject *py2c_factor_number_of_components
 ( PyObject *self, PyObject *args );
static PyObject *py2c_factor_witness_points_of_component
 ( PyObject *self, PyObject *args );
static PyObject *py2c_factor_trace_sum_difference
 ( PyObject *self, PyObject *args );
static PyObject *py2c_witness_set_of_hypersurface
 ( PyObject *self, PyObject *args );
static PyObject *py2c_create_diagonal_homotopy
 ( PyObject *self, PyObject *args );
static PyObject *py2c_start_diagonal_cascade_solutions
 ( PyObject *self, PyObject *args );
static PyObject *py2c_extrinsic_top_diagonal_dimension
 ( PyObject *self, PyObject *args );
static PyObject *py2c_collapse_diagonal ( PyObject *self, PyObject *args );

/* wrapping of Pieri and Littlewood-Richardson homotopies */

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
static PyObject *py2c_schubert_littlewood_richardson_homotopies
 ( PyObject *self, PyObject *args );
static PyObject *py2c_schubert_localization_poset
 ( PyObject *self, PyObject *args );
static PyObject *py2c_schubert_pieri_homotopies
 ( PyObject *self, PyObject *args );
static PyObject *py2c_schubert_osculating_planes
 ( PyObject *self, PyObject *args );
static PyObject *py2c_schubert_pieri_system
 ( PyObject *self, PyObject *args );

/* wrapping functions in mapcon.h starts from here */

static PyObject *py2c_mapcon_solve_system ( PyObject *self, PyObject *args );
static PyObject *py2c_mapcon_write_maps ( PyObject *self, PyObject *args );
static PyObject *py2c_mapcon_clear_maps ( PyObject *self, PyObject *args );
static PyObject *py2c_mapcon_top_dimension ( PyObject *self, PyObject *args );
static PyObject *py2c_mapcon_number_of_maps ( PyObject *self, PyObject *args );
static PyObject *py2c_mapcon_degree_of_map ( PyObject *self, PyObject *args );
static PyObject *py2c_mapcon_coefficients_of_map
 ( PyObject *self, PyObject *args );
static PyObject *py2c_mapcon_exponents_of_map
 ( PyObject *self, PyObject *args );

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
